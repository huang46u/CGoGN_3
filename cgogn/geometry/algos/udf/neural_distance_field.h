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
#ifndef CGOGN_GEOMETRY_ALGOS_UDF_NEURAL_DISTANCE_FIELD_H_
#define CGOGN_GEOMETRY_ALGOS_UDF_NEURAL_DISTANCE_FIELD_H_

#include <cgogn/core/functions/attributes.h>
#include <cgogn/core/functions/mesh_info.h>
#include <cgogn/geometry/algos/udf/alpha_projection.h>
#include <cgogn/geometry/algos/udf/distance_field.h>
#include <cgogn/geometry/algos/udf/input_geometry.h>
#include <cgogn/geometry/algos/udf/neural_alpha_projection.h>
#include <cgogn/geometry/algos/udf/neural_field_query.h>
#include <cgogn/geometry/algos/udf/neural_model_type.h>
#include <cgogn/geometry/algos/udf/reconstruction_options.h>
#include <cgogn/geometry/functions/bounding_box.h>
#include <cgogn/geometry/types/neural_field_forward.h>
#include <cgogn/geometry/types/neural_ray_traits.h>
#include <cgogn/geometry/types/ray_level_set_sampler_torch.h>
#include <cgogn/geometry/types/vector_traits.h>

#include <torch/script.h>
#include <torch/torch.h>

#include <algorithm>
#include <cstdint>
#include <filesystem>
#include <memory>
#include <tuple>
#include <vector>

namespace cgogn::geometry
{

// Distance field of a neural network (TorchScript module) loaded from a file. The network predicts the UDF or the MF
// value and its gradient; the samples are rays cast against the network, and the projections are Newton-like steps on
// the network. See distance_field.h for the interface.
template <typename SURFACE, typename POINTS>
class NeuralDistanceField
{
public:
	using PVertex = typename mesh_traits<POINTS>::Vertex;
	using SVertex = typename mesh_traits<SURFACE>::Vertex;
	using RayConfig = RaySamplerTraits<RaySamplerNeural, SURFACE, POINTS>;
	using RayTraits = typename RayConfig::SamplerTraits;
	using RaySampler = TorchRayLevelSetSampler<RayTraits>;
	static constexpr bool is_neural_type = true;

	NeuralDistanceField(InputGeometry<SURFACE, POINTS>& input, const ReconstructionOptions& options)
		: input_(input), options_(options)
	{
	}

	// the options (neural_model_type) are written by the caller
	bool load_model(const std::filesystem::path& path, NeuralModelType /*type*/)
	{
		try
		{
			neural_model_ = torch::jit::load(path.string(), device_);
			neural_model_.eval();
			neural_model_loaded_ = true;
			neural_field_ = NeuralFieldForward(&neural_model_, true, device_);
			return true;
		}
		catch (const c10::Error&)
		{
			neural_model_loaded_ = false;
			neural_field_ = NeuralFieldForward(&neural_model_, false, device_);
			return false;
		}
	}

	void reset()
	{
		neural_model_loaded_ = false;
		neural_field_ = NeuralFieldForward(&neural_model_, false, device_);
	}

	bool loaded() const
	{
		return neural_model_loaded_;
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
		if (options_.neural_model_type == NeuralModelType::udf)
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
		// NeuralFieldRayTraits stores a pointer, so the copy lives in a member and outlives the traits
		sampling_field_ = neural_field_;
		return RayConfig::make(sampling_field_);
	}

	std::unique_ptr<RaySampler> make_ray_sampler(const typename RaySampler::Params& params) const
	{
		return std::make_unique<RaySampler>(params, device_);
	}

	void configure_ray_sampler(RaySampler& sampler, const typename RaySampler::Params& params) const
	{
		sampler.set_params(params);
		sampler.set_device(device_);
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
		const bool filter_positive = options_.neural_model_type == NeuralModelType::mf;
		if (!geometry::project_points_to_alpha(neural_field_, points, pp, filter_positive, projection_workspace_,
											   projection))
			return false;
		return true;
	}

	// topology scores: the UDF value of the network at the points
	bool distances(const std::vector<Vec3>& points, std::vector<Scalar>& values)
	{
		return evaluate_udf_values(forward(), points, static_cast<size_t>(options_.batch_size), values);
	}

	// MF value and sdf of the network at the points
	bool values_sdf(const std::vector<Vec3>& points, std::vector<Scalar>& values, std::vector<Scalar>& sdf)
	{
		return evaluate_udf_values_sdf(forward(), points, static_cast<size_t>(options_.batch_size), values, sdf);
	}

	// the normals of the samples are the network gradient: they are not recomputed from the input after the filter
	static constexpr bool recompute_normals_from_input_after_filtering = false;

private:
	// the neural state of the field: the device, the module, whether it is loaded and the wrapper of the module
	torch::Device device_ = torch::cuda::is_available() ? torch::Device(torch::kCUDA, 0) : torch::Device(torch::kCPU);
	torch::jit::Module neural_model_;
	bool neural_model_loaded_ = false;
	// one wrapper for the loaded model: it remembers the output layout found by its first forward, so MF models
	// do not go through autograd again on every call
	NeuralFieldForward neural_field_{&neural_model_, false, device_};
	NeuralProjectionWorkspace projection_workspace_;
	// copy of neural_field_ taken by ray_traits(): the ray traits keep a pointer to it during a sampling
	NeuralFieldForward sampling_field_{&neural_model_, false, device_};
	InputGeometry<SURFACE, POINTS>& input_;
	const ReconstructionOptions& options_;
};

} // namespace cgogn::geometry

#endif // CGOGN_GEOMETRY_ALGOS_UDF_NEURAL_DISTANCE_FIELD_H_
