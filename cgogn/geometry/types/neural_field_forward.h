#ifndef CGOGN_GEOMETRY_TYPES_NEURAL_FIELD_FORWARD_H_
#define CGOGN_GEOMETRY_TYPES_NEURAL_FIELD_FORWARD_H_

#include <cgogn/geometry/types/vector_traits.h>

#include <torch/script.h>
#include <torch/torch.h>

#include <cmath>
#include <iostream>
#include <limits>
#include <utility>
#include <vector>

namespace cgogn
{

namespace geometry
{

struct BatchUDFResult
{
	std::vector<Scalar> values; // size N
	std::vector<Vec3> gradients; // size N
	bool ok = false;
};

class NeuralFieldForward
{
public:
	using Result = BatchUDFResult;

	NeuralFieldForward(torch::jit::Module* module, bool loaded, const torch::Device& device)
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

	bool output_kind_known() const
	{
		return output_kind_ != OutputKind::Unknown;
	}

	bool is_mf_model() const
	{
		return output_kind_ == OutputKind::MFTuple;
	}

	//Forward a single point
	Scalar forward_point(const Vec3& query_point) const
	{
		return forward_point_impl(query_point, device_);
	}

	Scalar forward_point_cpu(const Vec3& query_point) const
	{
		return forward_point_impl(query_point, torch::kCPU);
	}

	Scalar forward_point_gpu(const Vec3& query_point) const
	{
		return forward_point_impl(query_point, device_);
	}

	//Forward a single point with gradient
	std::pair<Scalar, Vec3> forward_point_with_grad(const Vec3& query_point) const
	{
		return forward_point_with_grad_impl(query_point, device_);
	}

	std::pair<Scalar, Vec3> forward_point_with_grad_cpu(const Vec3& query_point) const
	{
		return forward_point_with_grad_impl(query_point, torch::kCPU);
	}

	std::pair<Scalar, Vec3> forward_point_with_grad_gpu(const Vec3& query_point) const
	{
		return forward_point_with_grad_impl(query_point, device_);
	}

	// Forward a batch of points
	Result forward_batch(const std::vector<Vec3>& query_points) const
	{
		return forward_batch_impl(query_points, device_, false);
	}

	Result forward_batch_cpu(const std::vector<Vec3>& query_points) const
	{
		return forward_batch_impl(query_points, torch::kCPU, false);
	}

	Result forward_batch_gpu(const std::vector<Vec3>& query_points) const
	{
		return forward_batch_impl(query_points, device_, false);
	}

	// Forward a batch of points with gradient
	Result forward_batch_with_grad(const std::vector<Vec3>& query_points) const
	{
		return forward_batch_impl(query_points, device_, true);
	}

	Result forward_batch_with_grad_cpu(const std::vector<Vec3>& query_points) const
	{
		return forward_batch_impl(query_points, torch::kCPU, true);
	}

	Result forward_batch_with_grad_gpu(const std::vector<Vec3>& query_points) const
	{
		return forward_batch_impl(query_points, device_, true);
	}

	// Forward points in GPU tensor
	torch::Tensor forward_values_gpu(const torch::Tensor& points_gpu) const
	{
		if (!is_loaded())
			return torch::Tensor();
		try
		{
			ModuleDeviceScope module_scope(module_, device_, device_);
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

	// Forward points in GPU tensor while preserving the autograd graph
	torch::Tensor forward_values_autograd_gpu(const torch::Tensor& points_gpu) const
	{
		if (!is_loaded())
			return torch::Tensor();
		try
		{
			ModuleDeviceScope module_scope(module_, device_, device_);
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
			return parsed.values;
		}
		catch (const c10::Error& e)
		{
			std::cerr << "NeuralField forward (values+autograd) failed: " << e.what() << std::endl;
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
			ModuleDeviceScope module_scope(module_, device_, device_);
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

	// Forward points in GPU tensor with gradient
	std::pair<torch::Tensor, torch::Tensor> forward_values_grad_gpu(const torch::Tensor& points_gpu) const
	{
		if (!is_loaded())
			return {torch::Tensor(), torch::Tensor()};
		try
		{
			ModuleDeviceScope module_scope(module_, device_, device_);
			torch::Tensor x = points_gpu;
			if (x.device() != device_)
				x = x.to(device_);

			if (output_kind_ == OutputKind::MFTuple)
				return forward_values_grad_gpu_predicted(x);

			return forward_values_grad_gpu_autograd(x);
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

	class ModuleDeviceScope
	{
	public:
		ModuleDeviceScope(torch::jit::Module* module, const torch::Device& current_device,
						  const torch::Device& target_device)
			: module_(module), current_device_(current_device), moved_(false)
		{
			if (module_ && current_device_ != target_device)
			{
				module_->to(target_device);
				moved_ = true;
			}
		}

		~ModuleDeviceScope()
		{
			if (module_ && moved_)
				module_->to(current_device_);
		}

	private:
		torch::jit::Module* module_;
		torch::Device current_device_;
		bool moved_;
	};

	Scalar forward_point_impl(const Vec3& query_point, const torch::Device& device) const
	{
		Result r = forward_batch_impl({query_point}, device, false);
		if (!r.ok || r.values.empty())
			return std::numeric_limits<Scalar>::max();
		return r.values[0];
	}

	std::pair<Scalar, Vec3> forward_point_with_grad_impl(const Vec3& query_point, const torch::Device& device) const
	{
		Result r = forward_batch_impl({query_point}, device, true);
		if (!r.ok || r.values.empty() || r.gradients.empty())
			return {std::numeric_limits<Scalar>::max(), Vec3(0, 0, 0)};
		return {r.values[0], r.gradients[0]};
	}

	Result forward_batch_impl(const std::vector<Vec3>& query_points, const torch::Device& device, bool with_grad) const
	{
		Result r;
		r.ok = false;
		const size_t N = query_points.size();
		if (!is_loaded() || N == 0)
			return r;
		try
		{
			ModuleDeviceScope module_scope(module_, device_, device);
			if (!with_grad)
			{
				torch::NoGradGuard no_grad;
				return forward_batch_impl_no_grad(query_points, device);
			}
			return forward_batch_impl_with_grad(query_points, device);
		}
		catch (const c10::Error& e)
		{
			std::cerr << "NeuralField batch forward failed: " << e.what() << std::endl;
			return r;
		}
	}

	Result forward_batch_impl_no_grad(const std::vector<Vec3>& query_points, const torch::Device& device) const
	{
		Result r;
		r.ok = false;
		const size_t N = query_points.size();

		auto cpu_opts = torch::TensorOptions().dtype(torch::kFloat32).device(torch::kCPU);
		if (device.is_cuda())
			cpu_opts = cpu_opts.pinned_memory(true);

		torch::Tensor points_cpu = torch::empty({static_cast<long>(N), 3}, cpu_opts);
		auto acc = points_cpu.accessor<float, 2>();
		for (size_t i = 0; i < N; ++i)
		{
			acc[i][0] = query_points[i].x();
			acc[i][1] = query_points[i].y();
			acc[i][2] = query_points[i].z();
		}

		torch::Tensor points = points_cpu.to(device);
		torch::IValue output = module_->forward({points});

		ParsedOutput parsed;
		OutputKind kind = OutputKind::Unknown;
		if (!parse_output(output, parsed, kind))
			return r;
		if (kind != OutputKind::Unknown)
			output_kind_ = kind;

		torch::Tensor values_cpu = parsed.values.detach().to(torch::kCPU).contiguous();
		if (!fill_udf_values_from_tensor(r, values_cpu, N))
			return r;

		r.ok = true;
		return r;
	}

	Result forward_batch_impl_with_grad(const std::vector<Vec3>& query_points, const torch::Device& device) const
	{
		Result r;
		r.ok = false;
		const size_t N = query_points.size();

		auto cpu_opts = torch::TensorOptions().dtype(torch::kFloat32).device(torch::kCPU);
		if (device.is_cuda())
			cpu_opts = cpu_opts.pinned_memory(true);

		torch::Tensor points_cpu = torch::empty({static_cast<long>(N), 3}, cpu_opts);
		auto acc = points_cpu.accessor<float, 2>();
		for (size_t i = 0; i < N; ++i)
		{
			acc[i][0] = query_points[i].x();
			acc[i][1] = query_points[i].y();
			acc[i][2] = query_points[i].z();
		}

		if (output_kind_ == OutputKind::MFTuple)
			return forward_batch_impl_with_grad_predicted(points_cpu, device, N);

		return forward_batch_impl_with_grad_autograd(points_cpu, device, N);
	}

	Result forward_batch_impl_with_grad_predicted(const torch::Tensor& points_cpu, const torch::Device& device,
											 size_t N) const
	{
		Result r;
		r.ok = false;

		torch::NoGradGuard no_grad;
		torch::Tensor points = points_cpu.to(device);
		torch::IValue output = module_->forward({points});

		ParsedOutput parsed;
		OutputKind kind = OutputKind::Unknown;
		if (!parse_output(output, parsed, kind))
			return r;
		if (kind != OutputKind::Unknown)
			output_kind_ = kind;
		if (!parsed.has_pred_grad)
			return r;

		torch::Tensor values_cpu = parsed.values.detach().to(torch::kCPU).contiguous();
		if (!fill_udf_values_from_tensor(r, values_cpu, N))
			return r;

		torch::Tensor grad_cpu = parsed.pred_grad.detach().to(torch::kCPU).contiguous();
		if (!fill_gradients_from_tensor(r, grad_cpu, N))
			return r;

		r.ok = true;
		return r;
	}

	Result forward_batch_impl_with_grad_autograd(const torch::Tensor& points_cpu, const torch::Device& device,
											 size_t N) const
	{
		Result r;
		r.ok = false;

		torch::Tensor points = points_cpu.to(device);
		points.set_requires_grad(true);
		torch::IValue output = module_->forward({points});

		ParsedOutput parsed;
		OutputKind kind = OutputKind::Unknown;
		if (!parse_output(output, parsed, kind))
			return r;
		output_kind_ = kind;

		if (parsed.has_pred_grad)
		{
			torch::Tensor values_cpu = parsed.values.detach().to(torch::kCPU).contiguous();
			if (!fill_udf_values_from_tensor(r, values_cpu, N))
				return r;

			torch::Tensor grad_cpu = parsed.pred_grad.detach().to(torch::kCPU).contiguous();
			if (!fill_gradients_from_tensor(r, grad_cpu, N))
				return r;

			r.ok = true;
			return r;
		}

		torch::Tensor grad_outputs = torch::ones_like(parsed.values);
		std::vector<torch::Tensor> grads = torch::autograd::grad(
			/*outputs=*/{parsed.values},
			/*inputs=*/{points},
			/*grad_outputs=*/{grad_outputs},
			/*retain_graph=*/false,
			/*create_graph=*/false,
			/*allow_unused=*/false);

		torch::Tensor grad = grads[0];
		if (!validate_predicted_grad_tensor(grad, parsed.values.size(0)))
			return r;

		torch::Tensor values_cpu = parsed.values.detach().to(torch::kCPU).contiguous();
		if (!fill_udf_values_from_tensor(r, values_cpu, N))
			return r;

		torch::Tensor grad_cpu = grad.detach().to(torch::kCPU).contiguous();
		if (!fill_gradients_from_tensor(r, grad_cpu, N))
			return r;

		r.ok = true;
		return r;
	}

	std::pair<torch::Tensor, torch::Tensor> forward_values_grad_gpu_predicted(const torch::Tensor& points) const
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

		return {parsed.values.detach(), parsed.pred_grad.detach()};
	}

	std::pair<torch::Tensor, torch::Tensor> forward_values_grad_gpu_autograd(const torch::Tensor& points) const
	{
		torch::Tensor x = points.detach().requires_grad_(true);
		torch::IValue output = module_->forward({x});

		ParsedOutput parsed;
		OutputKind kind = OutputKind::Unknown;
		if (!parse_output(output, parsed, kind))
			return {torch::Tensor(), torch::Tensor()};
		output_kind_ = kind;

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
			parsed.values = mf - torch::abs(sdf);
			parsed.pred_grad = pgrad;
			parsed.has_pred_grad = true;
			return true;
		}

		kind = OutputKind::Unknown;
		return false;
	}

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

	bool fill_udf_values_from_tensor(Result& r, const torch::Tensor& values_cpu, size_t N) const
	{
		if (values_cpu.dim() != 1 || values_cpu.size(0) != static_cast<long>(N))
			return false;
		auto out_acc = values_cpu.accessor<float, 1>();
		r.values.resize(N);
		r.gradients.clear();
		for (size_t i = 0; i < N; ++i)
			r.values[i] = static_cast<Scalar>(out_acc[(long)i]);
		return true;
	}

	bool fill_gradients_from_tensor(Result& r, const torch::Tensor& grad_cpu, size_t N) const
	{
		if (!validate_predicted_grad_tensor(grad_cpu, static_cast<long>(N)))
			return false;
		auto grad_acc = grad_cpu.accessor<float, 2>();
		r.gradients.resize(N);
		for (size_t i = 0; i < N; ++i)
		{
			r.gradients[i] = Vec3(static_cast<Scalar>(grad_acc[(long)i][0]), static_cast<Scalar>(grad_acc[(long)i][1]),
							  static_cast<Scalar>(grad_acc[(long)i][2]));
		}
		return true;
	}

	torch::jit::Module* module_;
	bool loaded_;
	torch::Device device_;
	mutable OutputKind output_kind_ = OutputKind::Unknown;
};

} // namespace geometry

} // namespace cgogn

#endif // CGOGN_GEOMETRY_TYPES_NEURAL_FIELD_FORWARD_H_
