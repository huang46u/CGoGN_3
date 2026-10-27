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
#ifndef CGOGN_GEOMETRY_ALGOS_OFFSET_AXIS_NEURAL_NEURAL_DISTANCE_FIELD_H_
#define CGOGN_GEOMETRY_ALGOS_OFFSET_AXIS_NEURAL_NEURAL_DISTANCE_FIELD_H_

#include <cgogn/core/functions/attributes.h>
#include <cgogn/core/functions/mesh_info.h>
#include <cgogn/geometry/algos/offset_axis/alpha_projection.h>
#include <cgogn/geometry/algos/offset_axis/distance_field.h>
#include <cgogn/geometry/algos/offset_axis/input_geometry.h>
#include <cgogn/geometry/algos/offset_axis/neural/neural_alpha_projection.h>
#include <cgogn/geometry/algos/offset_axis/neural/neural_field_forward.h>
#include <cgogn/geometry/algos/offset_axis/neural/neural_field_query.h>
#include <cgogn/geometry/algos/offset_axis/neural/neural_ray_traits.h>
#include <cgogn/geometry/algos/offset_axis/neural/ray_level_set_sampler_torch.h>
#include <cgogn/geometry/algos/offset_axis/neural_model_type.h>
#include <cgogn/geometry/algos/offset_axis/reconstruction_options.h>
#include <cgogn/geometry/functions/bounding_box.h>
#include <cgogn/geometry/types/vector_traits.h>

#include <torch/script.h>
#include <torch/torch.h>

#include <algorithm>
#include <cstdint>
#include <filesystem>
#include <memory>
#include <string>
#include <tuple>
#include <vector>

namespace cgogn
{

namespace geometry
{

namespace offset_axis
{

// true if LibTorch sees a CUDA device (neural fields then run on GPU 0)
inline bool neural_cuda_available()
{
	return torch::cuda::is_available();
}

// Distance field of a neural network (TorchScript module) loaded from a file. The network predicts the UDF or the MF
// value and its gradient; the samples are rays cast against the network, and the projections are Newton-like steps on
// the network. See distance_field.h for the interface.
template <typename SURFACE, typename POINTS>
class NeuralDistanceField
{
public:
	using PVertex = typename mesh_traits<POINTS>::Vertex;
	using SVertex = typename mesh_traits<SURFACE>::Vertex;
	using RayTraits = NeuralFieldRayTraits;
	static constexpr bool is_neural_type = true;

	NeuralDistanceField(InputGeometry<SURFACE, POINTS>& input, const ReconstructionOptions& options)
		: input_(input), options_(options)
	{
	}

	// loads the TorchScript module; the field keeps the model type (UDF or MF) given here
	bool load_model(const std::filesystem::path& path, NeuralModelType type)
	{
		try
		{
			neural_model_ = std::make_shared<torch::jit::Module>(torch::jit::load(path.string(), device_));
			neural_model_->eval();
			neural_model_loaded_ = true;
			model_type_ = type;
			neural_field_ = NeuralFieldForward(neural_model_, true, device_);
			load_error_.clear();
			return true;
		}
		catch (const c10::Error& e)
		{
			neural_model_loaded_ = false;
			neural_field_ = NeuralFieldForward(neural_model_, false, device_);
			load_error_ = e.what_without_backtrace();
			return false;
		}
	}

	// the reason of the last failed load_model (empty after a successful load)
	const std::string& load_error() const
	{
		return load_error_;
	}

	void reset()
	{
		neural_model_loaded_ = false;
		neural_field_ = NeuralFieldForward(neural_model_, false, device_);
	}

	bool loaded() const
	{
		return neural_model_loaded_;
	}

	NeuralModelType model_type() const
	{
		return model_type_;
	}

	const torch::Device& device() const
	{
		return device_;
	}

	NeuralFieldForward& forward()
	{
		return neural_field_;
	}

	NeuralProjectionWorkspace& workspace()
	{
		return projection_workspace_;
	}

	FieldPrepareResult prepare(const PrepareInputs<SURFACE, POINTS>& /*inputs*/)
	{
		if (!neural_model_loaded_)
			return FieldPrepareResult::model_not_loaded;
		return FieldPrepareResult::ok;
	}

	FieldPrepareResult sampling_box(const PrepareInputs<SURFACE, POINTS>& inputs, Vec3& bbox_min, Vec3& bbox_max)
	{
		if (model_type_ == NeuralModelType::udf)
		{
			bbox_min = Vec3(-0.5, -0.5, -0.5);
			bbox_max = Vec3(0.5, 0.5, 0.5);
		}
		if (inputs.surface && nb_cells<SVertex>(*inputs.surface) > 0)
		{
			auto position = get_attribute<Vec3, SVertex>(*inputs.surface, "position");
			if (position)
				std::tie(bbox_min, bbox_max) = geometry::bounding_box(*position);
		}
		else if (inputs.input_points && nb_cells<PVertex>(*inputs.input_points) > 0)
			std::tie(bbox_min, bbox_max) = geometry::bounding_box(*inputs.input_position);
		return FieldPrepareResult::ok;
	}

	RayTraits ray_traits(const PrepareInputs<SURFACE, POINTS>& /*inputs*/)
	{
		RayTraits traits;
		// the traits point to the member neural_field_, which outlives the sampling
		traits.field = &neural_field_;
		// the GPU level set sampler keeps its ray buffers here; they live on the device of the field
		ray_buffers_.device = device_;
		traits.ray_buffers = &ray_buffers_;
		return traits;
	}

	bool normals(const std::vector<Vec3>& points, std::vector<Vec3>& normals, std::vector<uint8_t>& keep)
	{
		if (!evaluate_udf_normals(neural_field_, points, static_cast<size_t>(options_.batch_size),
								  projection_workspace_, normals))
			return false;
		std::fill(keep.begin(), keep.end(), uint8_t(1));
		return true;
	}

	bool project(const std::vector<Vec3>& points, const AlphaProjectionParameters& pp,
				 AlphaProjectionResult& projection)
	{
		const bool filter_positive = model_type_ == NeuralModelType::mf;
		if (!offset_axis::project_points_to_alpha(neural_field_, points, pp, filter_positive, projection_workspace_,
											   projection))
			return false;
		return true;
	}

	bool distances(const std::vector<Vec3>& points, std::vector<Scalar>& values)
	{
		return evaluate_udf_values(forward(), points, static_cast<size_t>(options_.batch_size), values);
	}

	bool values_sdf(const std::vector<Vec3>& points, std::vector<Scalar>& values, std::vector<Scalar>& sdf)
	{
		return evaluate_udf_values_sdf(forward(), points, static_cast<size_t>(options_.batch_size), values, sdf);
	}

	// the normals of the samples are the network gradient
	static constexpr bool recompute_normals_from_input_after_filtering = false;

private:
	torch::Device device_ = torch::cuda::is_available() ? torch::Device(torch::kCUDA, 0) : torch::Device(torch::kCPU);
	std::shared_ptr<torch::jit::Module> neural_model_;
	bool neural_model_loaded_ = false;
	NeuralModelType model_type_ = NeuralModelType::udf;
	std::string load_error_;
	// one wrapper for the loaded model: it remembers the output layout found by its first forward, so MF models
	// do not go through autograd again on every call
	NeuralFieldForward neural_field_{neural_model_, false, device_};
	NeuralProjectionWorkspace projection_workspace_;
	TorchRayBuffers ray_buffers_;
	InputGeometry<SURFACE, POINTS>& input_;
	const ReconstructionOptions& options_;
};

} // namespace offset_axis

} // namespace geometry

} // namespace cgogn

#endif // CGOGN_GEOMETRY_ALGOS_OFFSET_AXIS_NEURAL_NEURAL_DISTANCE_FIELD_H_
