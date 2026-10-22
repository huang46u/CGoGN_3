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
#ifndef CGOGN_GEOMETRY_ALGOS_UDF_SURFACE_DISTANCE_FIELD_H_
#define CGOGN_GEOMETRY_ALGOS_UDF_SURFACE_DISTANCE_FIELD_H_

#include <cgogn/core/functions/attributes.h>
#include <cgogn/core/functions/mesh_info.h>
#include <cgogn/geometry/algos/udf/alpha_projection.h>
#include <cgogn/geometry/algos/udf/distance_field.h>
#include <cgogn/geometry/algos/udf/input_geometry.h>
#include <cgogn/geometry/algos/udf/neural_model_type.h>
#include <cgogn/geometry/algos/udf/normal_utils.h>
#include <cgogn/geometry/algos/udf/parallel_for.h>
#include <cgogn/geometry/algos/udf/reconstruction_options.h>
#include <cgogn/geometry/algos/udf/spatial_query.h>
#include <cgogn/geometry/functions/bounding_box.h>
#include <cgogn/geometry/types/ray_level_set_sampler_traits.h>
#include <cgogn/geometry/types/vector_traits.h>

#include <atomic>
#include <cstdint>
#include <filesystem>
#include <tuple>
#include <vector>

namespace cgogn::geometry
{

// Distance field of a triangle surface input: rays are cast against the BVH of the surface, and the normals and the
// projections are closest-point queries on the surface. See distance_field.h for the interface.
template <typename SURFACE, typename POINTS>
class SurfaceDistanceField
{
public:
	using SVertex = typename mesh_traits<SURFACE>::Vertex;
	using RayConfig = RaySamplerTraits<RaySamplerSurface, SURFACE, POINTS>;
	using RayTraits = typename RayConfig::SamplerTraits;
	static constexpr bool is_neural_type = false;

	// the surface field does not use the options
	SurfaceDistanceField(InputGeometry<SURFACE, POINTS>& input, const ReconstructionOptions& /*options*/)
		: input_(input)
	{
	}

	FieldPrepareResult prepare(const PrepareInputs<SURFACE, POINTS>& /*inputs*/)
	{
		input_.build_surface_bvh();
		return FieldPrepareResult::ok;
	}

	FieldPrepareResult sampling_box(const PrepareInputs<SURFACE, POINTS>& inputs, Vec3& bbox_min, Vec3& bbox_max)
	{
		if (!inputs.surface || !input_.surface_bvh())
			return FieldPrepareResult::invalid_data;
		auto s_pos = get_attribute<Vec3, SVertex>(*inputs.surface, "position");
		std::tie(bbox_min, bbox_max) = geometry::bounding_box(*s_pos);
		return FieldPrepareResult::ok;
	}

	RayTraits ray_traits(const PrepareInputs<SURFACE, POINTS>& inputs)
	{
		auto s_pos = get_attribute<Vec3, SVertex>(*inputs.surface, "position");
		return RayConfig::make(*inputs.surface, s_pos.get(), input_.surface_bvh().get(), &input_.surface_faces());
	}

	torch::Device sampling_device() const
	{
		return torch::kCPU;
	}

	bool normals(const std::vector<Vec3>& points, std::vector<Vec3>& normals, std::vector<uint8_t>& keep)
	{
		detail::parallel_for_index(points.size(), [&](std::size_t i) {
			AlphaProjectionInfo info;
			const bool query_ok = input_.query_surface_projection_info(points[i], info);
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
			[&](const Vec3& p, AlphaProjectionInfo& info) { return input_.query_surface_projection_info(p, info); },
			projection);
	}

	// topology scores: the unsigned distance from the points to the surface (closest point on the BVH)
	bool distances(const std::vector<Vec3>& points, std::vector<Scalar>& values)
	{
		input_.build_surface_bvh();
		if (!input_.surface_bvh())
			return false;
		values.resize(points.size());
		std::atomic<bool> all_found{true};
		detail::parallel_for_index(points.size(), [&](std::size_t i) {
			uint32 primitive;
			Vec3 closest;
			Scalar distance;
			if (!query_surface_closest_point(*input_.surface_bvh(), points[i], primitive, closest, distance))
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

	// MF values are only defined for neural fields; this method is not called for the surface field
	bool values_sdf(const std::vector<Vec3>& /*points*/, std::vector<Scalar>& /*values*/, std::vector<Scalar>& /*sdf*/)
	{
		return false;
	}

	// the normals of the samples are recomputed from the input after the sampling filter
	static constexpr bool recompute_normals_from_input_after_filtering = true;

	// the surface field never holds a neural model (see NeuralDistanceField)
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
};

} // namespace cgogn::geometry

#endif // CGOGN_GEOMETRY_ALGOS_UDF_SURFACE_DISTANCE_FIELD_H_
