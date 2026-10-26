/*******************************************************************************
 * CGoGN: Combinatorial and Geometric modeling with Generic N-dimensional Maps  *
 * Copyright (C), IGG Group, ICube, University of Strasbourg, France            *
 *                                                                              *
 * This library is free software; you can redistribute it and/or modify it      *
 * under the terms of the GNU Lesser General Public License as published by the *
 * Free Software Foundation; either version 2.1 of the License, or (at your     *
 * option) any later version.                                                   *
 *                                                                              *
 * This library is distributed in the hope that it will be useful, but WITHOUT *
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or        *
 * FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License *
 * for more details.                                                            *
 *                                                                              *
 * You should have received a copy of the GNU Lesser General Public License     *
 * along with this library; if not, write to the Free Software Foundation,      *
 * Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301 USA.           *
 *                                                                              *
 * Web site: http://cgogn.unistra.fr/                                           *
 * Contact information: cgogn@unistra.fr                                        *
 *                                                                              *
 *******************************************************************************/
#ifndef CGOGN_GEOMETRY_ALGOS_OFFSET_AXIS_NEURAL_NEURAL_FIELD_FORWARD_H_
#define CGOGN_GEOMETRY_ALGOS_OFFSET_AXIS_NEURAL_NEURAL_FIELD_FORWARD_H_

#include <cgogn/geometry/types/vector_traits.h>

#include <torch/script.h>
#include <torch/torch.h>

#include <cmath>
#include <iostream>
#include <memory>
#include <utility>
#include <vector>

namespace cgogn
{

namespace geometry
{

namespace offset_axis
{

class NeuralFieldForward
{
public:
	NeuralFieldForward(std::shared_ptr<torch::jit::Module> module, bool loaded, const torch::Device& device)
		: module_(module), loaded_(loaded), device_(device)
	{
	}

	bool is_loaded() const
	{
		return loaded_ && module_;
	}

	const torch::Device& device() const
	{
		return device_;
	}

	// Forward points in GPU tensor
	torch::Tensor forward_values_gpu(const torch::Tensor& points_gpu) const
	{
		if (!is_loaded())
			return torch::Tensor();
		try
		{
			torch::NoGradGuard no_grad;
			torch::Tensor x = points_gpu;
			if (x.device() != device_)
				x = x.to(device_);

			torch::IValue output = module_->forward({x});
			ParsedOutput parsed;
			OutputKind kind = OutputKind::Unknown;
			if (!parse_output(output, parsed, kind))
				return torch::Tensor();
			if (kind != OutputKind::Unknown)
				output_kind_ = kind;
			return parsed.values.detach();
		}
		catch (const c10::Error& e)
		{
			std::cerr << "NeuralField forward (values) failed: " << e.what() << std::endl;
			return torch::Tensor();
		}
	}

	// Forward points in GPU tensor and return UDF values + optional SDF values
	std::pair<torch::Tensor, torch::Tensor> forward_values_sdf_gpu(const torch::Tensor& points_gpu) const
	{
		if (!is_loaded())
			return {torch::Tensor(), torch::Tensor()};
		try
		{
			torch::NoGradGuard no_grad;
			torch::Tensor x = points_gpu;
			if (x.device() != device_)
				x = x.to(device_);

			torch::IValue output = module_->forward({x});
			ParsedOutput parsed;
			OutputKind kind = OutputKind::Unknown;
			if (!parse_output(output, parsed, kind))
				return {torch::Tensor(), torch::Tensor()};
			if (kind != OutputKind::Unknown)
				output_kind_ = kind;

			if (parsed.has_sdf)
				return {parsed.values.detach(), parsed.sdf.detach()};
			return {parsed.values.detach(), torch::Tensor()};
		}
		catch (const c10::Error& e)
		{
			std::cerr << "NeuralField forward (values+sdf) failed: " << e.what() << std::endl;
			return {torch::Tensor(), torch::Tensor()};
		}
	}

	// Forward points in GPU tensor with gradient; sdf_out (optional) receives the sdf output of MF models
	// from the same forward, and stays undefined for models without one
	std::pair<torch::Tensor, torch::Tensor> forward_values_grad_gpu(const torch::Tensor& points_gpu,
																	torch::Tensor* sdf_out = nullptr) const
	{
		if (sdf_out)
			*sdf_out = torch::Tensor();
		if (!is_loaded())
			return {torch::Tensor(), torch::Tensor()};
		try
		{
			torch::Tensor x = points_gpu;
			if (x.device() != device_)
				x = x.to(device_);

			if (output_kind_ == OutputKind::MFTuple)
				return forward_values_grad_gpu_predicted(x, sdf_out);

			return forward_values_grad_gpu_autograd(x, sdf_out);
		}
		catch (const c10::Error& e)
		{
			std::cerr << "NeuralField forward (values+grad) failed: " << e.what() << std::endl;
			return {torch::Tensor(), torch::Tensor()};
		}
	}

private:
	enum class OutputKind
	{
		Unknown,
		TensorScalar,
		MFTuple
	};

	struct ParsedOutput
	{
		torch::Tensor values;
		torch::Tensor pred_grad;
		bool has_pred_grad = false;
		torch::Tensor sdf;
		bool has_sdf = false;
	};

	std::pair<torch::Tensor, torch::Tensor> forward_values_grad_gpu_predicted(const torch::Tensor& points,
																			  torch::Tensor* sdf_out) const
	{
		torch::NoGradGuard no_grad;
		torch::IValue output = module_->forward({points});

		ParsedOutput parsed;
		OutputKind kind = OutputKind::Unknown;
		if (!parse_output(output, parsed, kind))
			return {torch::Tensor(), torch::Tensor()};
		if (kind != OutputKind::Unknown)
			output_kind_ = kind;
		if (!parsed.has_pred_grad)
			return {torch::Tensor(), torch::Tensor()};
		if (sdf_out && parsed.has_sdf)
			*sdf_out = parsed.sdf.detach();

		return {parsed.values.detach(), parsed.pred_grad.detach()};
	}

	std::pair<torch::Tensor, torch::Tensor> forward_values_grad_gpu_autograd(const torch::Tensor& points,
																			 torch::Tensor* sdf_out) const
	{
		torch::Tensor x = points.detach().requires_grad_(true);
		torch::IValue output = module_->forward({x});

		ParsedOutput parsed;
		OutputKind kind = OutputKind::Unknown;
		if (!parse_output(output, parsed, kind))
			return {torch::Tensor(), torch::Tensor()};
		if (kind != OutputKind::Unknown)
			output_kind_ = kind;
		if (sdf_out && parsed.has_sdf)
			*sdf_out = parsed.sdf.detach();

		if (parsed.has_pred_grad)
			return {parsed.values.detach(), parsed.pred_grad.detach()};

		torch::Tensor grad_outputs = torch::ones_like(parsed.values);
		std::vector<torch::Tensor> grads = torch::autograd::grad(
			/*outputs=*/{parsed.values},
			/*inputs=*/{x},
			/*grad_outputs=*/{grad_outputs},
			/*retain_graph=*/false,
			/*create_graph=*/false,
			/*allow_unused=*/false);

		torch::Tensor dx = grads[0];
		if (!validate_predicted_grad_tensor(dx, parsed.values.size(0)))
			return {torch::Tensor(), torch::Tensor()};

		return {parsed.values.detach(), dx.detach()};
	}

	bool parse_output(const torch::IValue& output, ParsedOutput& parsed, OutputKind& kind) const
	{
		parsed = ParsedOutput{};
		if (output.isTensor())
		{
			kind = OutputKind::TensorScalar;
			torch::Tensor values = output.toTensor();
			if (!normalize_scalar_tensor(values))
				return false;
			parsed.values = values;
			return true;
		}

		if (output.isTuple())
		{
			kind = OutputKind::MFTuple;
			auto tuple = output.toTuple();
			if (!tuple)
				return false;
			const auto& elements = tuple->elements();
			if (elements.size() != 3)
				return false;
			if (!elements[0].isTensor() || !elements[1].isTensor() || !elements[2].isTensor())
				return false;

			torch::Tensor sdf = elements[0].toTensor();
			torch::Tensor pgrad = elements[1].toTensor();
			torch::Tensor mf = elements[2].toTensor();

			if (!normalize_scalar_tensor(sdf) || !normalize_scalar_tensor(mf))
				return false;
			if (!validate_predicted_grad_tensor(pgrad, sdf.size(0)))
				return false;
			if (mf.size(0) != sdf.size(0))
				return false;

			parsed.sdf = sdf;
			parsed.has_sdf = true;
			// medial field (Deep Medial Field / QMDF): mf - |sdf| vanishes on the medial axis, not on the surface, so its
			// alpha level set is a thin layer around the medial axis and its minimum is not zero; filtering MF samples by
			// their distance to the input surface would remove most of them
			parsed.values = mf - torch::abs(sdf);
			parsed.pred_grad = pgrad;
			parsed.has_pred_grad = true;
			return true;
		}

		kind = OutputKind::Unknown;
		return false;
	}

	// (n, 1) -> (n): every value / sdf tensor returned by this class is 1-D
	bool normalize_scalar_tensor(torch::Tensor& values) const
	{
		if (values.dim() == 2 && values.size(1) == 1)
			values = values.squeeze(1);
		return values.dim() == 1;
	}

	bool validate_predicted_grad_tensor(const torch::Tensor& grad, int64_t expected_rows) const
	{
		return grad.defined() && grad.dim() == 2 && grad.size(0) == expected_rows && grad.size(1) == 3;
	}

	std::shared_ptr<torch::jit::Module> module_;
	bool loaded_;
	torch::Device device_;
	mutable OutputKind output_kind_ = OutputKind::Unknown;
};

} // namespace offset_axis

} // namespace geometry

} // namespace cgogn

#endif // CGOGN_GEOMETRY_ALGOS_OFFSET_AXIS_NEURAL_NEURAL_FIELD_FORWARD_H_
