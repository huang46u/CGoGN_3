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
#ifndef CGOGN_GEOMETRY_ALGOS_OFFSET_AXIS_SURFACE_DISTANCE_FIELD_H_
#define CGOGN_GEOMETRY_ALGOS_OFFSET_AXIS_SURFACE_DISTANCE_FIELD_H_

#include <cgogn/core/functions/attributes.h>
#include <cgogn/core/functions/mesh_info.h>
#include <cgogn/geometry/algos/offset_axis/alpha_projection.h>
#include <cgogn/geometry/algos/offset_axis/cpu_distance_field.h>
#include <cgogn/geometry/algos/offset_axis/distance_field.h>
#include <cgogn/geometry/algos/offset_axis/input_geometry.h>
#include <cgogn/geometry/algos/offset_axis/reconstruction_options.h>
#include <cgogn/geometry/algos/offset_axis/spatial_query.h>
#include <cgogn/geometry/functions/bounding_box.h>
#include <cgogn/geometry/algos/offset_axis/ray_level_set_sampler_traits.h>
#include <cgogn/geometry/types/vector_traits.h>

#include <tuple>
#include <vector>

namespace cgogn
{

namespace geometry
{

namespace offset_axis
{

// Distance field of a triangle surface input: rays are cast against the BVH of the surface, and the normals and the
// projections are closest-point queries on the surface. See distance_field.h for the interface.
template <typename SURFACE, typename POINTS>
class SurfaceDistanceField : public CpuDistanceField<SurfaceMeshRayTraits<SURFACE>>
{
	using Base = CpuDistanceField<SurfaceMeshRayTraits<SURFACE>>;

public:
	using SVertex = typename mesh_traits<SURFACE>::Vertex;
	using typename Base::RayTraits;
	using typename Base::RaySampler;

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
		RayTraits traits;
		traits.mesh = inputs.surface;
		traits.position = s_pos.get();
		traits.bvh = input_.surface_bvh().get();
		traits.bvh_faces = &input_.surface_faces();
		return traits;
	}

	bool normals(const std::vector<Vec3>& points, std::vector<Vec3>& normals, std::vector<uint8_t>& keep)
	{
		return this->normals_from(points, normals, keep, [&](const Vec3& p, AlphaProjectionInfo& info) {
			return input_.query_surface_projection_info(p, info);
		});
	}

	bool project(const std::vector<Vec3>& points, const AlphaProjectionParameters& pp,
				 AlphaProjectionResult& projection)
	{
		return this->project_with(points, pp, projection, [&](const Vec3& p, AlphaProjectionInfo& info) {
			return input_.query_surface_projection_info(p, info);
		});
	}

	// topology scores: the unsigned distance from the points to the surface (closest point on the BVH)
	bool distances(const std::vector<Vec3>& points, std::vector<Scalar>& values)
	{
		input_.build_surface_bvh();
		if (!input_.surface_bvh())
			return false;
		return this->distances_from(points, values, [&](const Vec3& p, Scalar& distance) {
			uint32 primitive;
			Vec3 closest;
			return query_surface_closest_point(*input_.surface_bvh(), p, primitive, closest, distance);
		});
	}

private:
	InputGeometry<SURFACE, POINTS>& input_;
};

} // namespace offset_axis

} // namespace geometry

} // namespace cgogn

#endif // CGOGN_GEOMETRY_ALGOS_OFFSET_AXIS_SURFACE_DISTANCE_FIELD_H_
