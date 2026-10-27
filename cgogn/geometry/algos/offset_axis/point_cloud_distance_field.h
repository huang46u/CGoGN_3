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
#ifndef CGOGN_GEOMETRY_ALGOS_OFFSET_AXIS_POINT_CLOUD_DISTANCE_FIELD_H_
#define CGOGN_GEOMETRY_ALGOS_OFFSET_AXIS_POINT_CLOUD_DISTANCE_FIELD_H_

#include <cgogn/core/functions/attributes.h>
#include <cgogn/core/functions/mesh_info.h>
#include <cgogn/geometry/algos/offset_axis/alpha_projection.h>
#include <cgogn/geometry/algos/offset_axis/cpu_distance_field.h>
#include <cgogn/geometry/algos/offset_axis/distance_field.h>
#include <cgogn/geometry/algos/offset_axis/input_geometry.h>
#include <cgogn/geometry/algos/offset_axis/reconstruction_options.h>
#include <cgogn/geometry/algos/offset_axis/sample_processing.h>
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

// Distance field of an oriented point cloud input: the normals are estimated from the k nearest neighbours of the
// input points, and the projections are nearest-point queries on the k-d tree of the input. See distance_field.h for
// the interface.
template <typename SURFACE, typename POINTS>
class PointCloudDistanceField : public CpuDistanceField<PointCloudRayTraits<POINTS>>
{
	using Base = CpuDistanceField<PointCloudRayTraits<POINTS>>;

public:
	using PVertex = typename mesh_traits<POINTS>::Vertex;
	using typename Base::RayTraits;

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
		// prepare succeeded: the point cloud is not empty
		std::tie(bbox_min, bbox_max) = geometry::bounding_box(*inputs.input_position);
		return FieldPrepareResult::ok;
	}

	RayTraits ray_traits(const PrepareInputs<SURFACE, POINTS>& inputs)
	{
		RayTraits traits;
		traits.mesh = inputs.input_points;
		traits.position = inputs.input_position.get();
		traits.normal = inputs.input_normal.get();
		traits.knn = inputs.input_knn.get();
		traits.kdtree = input_.input_kdtree().get();
		traits.kdtree_vertices = &input_.input_vertices();
		return traits;
	}

	bool normals(const std::vector<Vec3>& points, std::vector<Vec3>& normals, std::vector<uint8_t>& keep)
	{
		return this->normals_from(points, normals, keep, [&](const Vec3& p, AlphaProjectionInfo& info) {
			return input_.query_point_cloud_projection_info(p, info);
		});
	}

	bool project(const std::vector<Vec3>& points, const AlphaProjectionParameters& pp,
				 AlphaProjectionResult& projection)
	{
		return this->project_with(points, pp, projection, [&](const Vec3& p, AlphaProjectionInfo& info) {
			return input_.query_point_cloud_projection_info(p, info);
		});
	}

	bool distances(const std::vector<Vec3>& points, std::vector<Scalar>& values)
	{
		input_.build_input_kdtree();
		if (!input_.input_kdtree())
			return false;
		return this->distances_from(points, values, [&](const Vec3& p, Scalar& distance) {
			uint32 point;
			return query_point_cloud_nearest_point(*input_.input_kdtree(), p, point, distance);
		});
	}

private:
	InputGeometry<SURFACE, POINTS>& input_;
	const ReconstructionOptions& options_;
};

} // namespace offset_axis

} // namespace geometry

} // namespace cgogn

#endif // CGOGN_GEOMETRY_ALGOS_OFFSET_AXIS_POINT_CLOUD_DISTANCE_FIELD_H_
