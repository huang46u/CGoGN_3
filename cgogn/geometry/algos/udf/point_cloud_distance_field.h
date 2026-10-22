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
#ifndef CGOGN_GEOMETRY_ALGOS_UDF_POINT_CLOUD_DISTANCE_FIELD_H_
#define CGOGN_GEOMETRY_ALGOS_UDF_POINT_CLOUD_DISTANCE_FIELD_H_

#include <cgogn/core/functions/attributes.h>
#include <cgogn/core/functions/mesh_info.h>
#include <cgogn/geometry/algos/udf/alpha_projection.h>
#include <cgogn/geometry/algos/udf/distance_field.h>
#include <cgogn/geometry/algos/udf/input_geometry.h>
#include <cgogn/geometry/algos/udf/neural_model_type.h>
#include <cgogn/geometry/algos/udf/normal_utils.h>
#include <cgogn/geometry/algos/udf/parallel_for.h>
#include <cgogn/geometry/algos/udf/reconstruction_options.h>
#include <cgogn/geometry/algos/udf/sample_processing.h>
#include <cgogn/geometry/algos/udf/spatial_query.h>
#include <cgogn/geometry/functions/bounding_box.h>
#include <cgogn/geometry/types/ray_level_set_sampler.h>
#include <cgogn/geometry/types/ray_level_set_sampler_traits.h>
#include <cgogn/geometry/types/vector_traits.h>

#include <atomic>
#include <cstdint>
#include <filesystem>
#include <memory>
#include <tuple>
#include <vector>

namespace cgogn::geometry
{

// Distance field of an oriented point cloud input: the normals are estimated from the k nearest neighbours of the
// input points, and the projections are nearest-point queries on the k-d tree of the input. See distance_field.h for
// the interface.
template <typename SURFACE, typename POINTS>
class PointCloudDistanceField
{
public:
	using PVertex = typename mesh_traits<POINTS>::Vertex;
	using RayConfig = RaySamplerTraits<RaySamplerPointCloud, SURFACE, POINTS>;
	using RayTraits = typename RayConfig::SamplerTraits;
	using RaySampler = RayLevelSetSampler<RayTraits>;
	static constexpr bool is_neural_type = false;

	PointCloudDistanceField(InputGeometry<SURFACE, POINTS>& input, const ReconstructionOptions& options)
		: input_(input), options_(options)
	{
	}

	FieldPrepareResult prepare(const PrepareInputs<SURFACE, POINTS>& inputs)
	{
		if (!input_.input_kdtree() || !inputs.input_normal || !inputs.input_knn)
			return FieldPrepareResult::invalid_data;
		compute_input_normals(*inputs.input_points, *inputs.input_position, *inputs.input_normal, *inputs.input_knn,
							  *input_.input_kdtree(), input_.input_vertices(), options_.knn_k);
		return FieldPrepareResult::ok;
	}

	FieldPrepareResult sampling_box(const PrepareInputs<SURFACE, POINTS>& inputs, Vec3& bbox_min, Vec3& bbox_max)
	{
		if (inputs.input_points && nb_cells<PVertex>(*inputs.input_points) > 0)
			std::tie(bbox_min, bbox_max) = geometry::bounding_box(*inputs.input_position);
		return FieldPrepareResult::ok;
	}

	RayTraits ray_traits(const PrepareInputs<SURFACE, POINTS>& inputs)
	{
		return RayConfig::make(*inputs.input_points, inputs.input_position.get(), inputs.input_normal.get(),
							   inputs.input_knn.get(), input_.input_kdtree().get(), &input_.input_vertices());
	}

	std::unique_ptr<RaySampler> make_ray_sampler(const typename RaySampler::Params& params) const
	{
		return std::make_unique<RaySampler>(params);
	}

	void configure_ray_sampler(RaySampler& sampler, const typename RaySampler::Params& params) const
	{
		sampler.set_params(params);
	}

	bool normals(const std::vector<Vec3>& points, std::vector<Vec3>& normals, std::vector<uint8_t>& keep)
	{
		detail::parallel_for_index(points.size(), [&](std::size_t i) {
			AlphaProjectionInfo info;
			const bool query_ok = input_.query_point_cloud_projection_info(points[i], info);
			if (!query_ok || !info.valid || !points[i].allFinite() || !info.normal.allFinite())
				return;
			normals[i] = unit_normal_or_z(info.normal);
			keep[i] = uint8_t(1);
		});
		return true;
	}

	bool project(const std::vector<Vec3>& points, const AlphaProjectionParameters& pp,
				 AlphaProjectionResult& projection)
	{
		return geometry::project_points_to_alpha(
			points, pp,
			[&](const Vec3& p, AlphaProjectionInfo& info) { return input_.query_point_cloud_projection_info(p, info); },
			projection);
	}

	// topology scores: the unsigned distance from the points to the input points (nearest point on the k-d tree)
	bool distances(const std::vector<Vec3>& points, std::vector<Scalar>& values)
	{
		input_.build_input_kdtree();
		if (!input_.input_kdtree())
			return false;
		values.resize(points.size());
		std::atomic<bool> all_found{true};
		detail::parallel_for_index(points.size(), [&](std::size_t i) {
			uint32 point;
			Scalar distance;
			if (!query_point_cloud_nearest_point(*input_.input_kdtree(), points[i], point, distance))
			{
				all_found = false;
				return;
			}
			values[i] = distance;
		});
		if (!all_found)
			return false;
		return true;
	}

	// MF values are only defined for neural fields; this method is not called for the point cloud field
	bool values_sdf(const std::vector<Vec3>& /*points*/, std::vector<Scalar>& /*values*/, std::vector<Scalar>& /*sdf*/)
	{
		return false;
	}

	// the normals of the samples are recomputed from the input after the sampling filter
	static constexpr bool recompute_normals_from_input_after_filtering = true;

	// the point cloud field never holds a neural model (see NeuralDistanceField)
	bool load_model(const std::filesystem::path& /*path*/, NeuralModelType /*type*/)
	{
		return false;
	}
	void reset()
	{
	}
	bool loaded() const
	{
		return false;
	}

private:
	InputGeometry<SURFACE, POINTS>& input_;
	const ReconstructionOptions& options_;
};

} // namespace cgogn::geometry

#endif // CGOGN_GEOMETRY_ALGOS_UDF_POINT_CLOUD_DISTANCE_FIELD_H_
