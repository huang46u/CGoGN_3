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
#ifndef CGOGN_GEOMETRY_ALGOS_OFFSET_AXIS_SAMPLING_FILTER_H_
#define CGOGN_GEOMETRY_ALGOS_OFFSET_AXIS_SAMPLING_FILTER_H_

#include <cgogn/core/functions/attributes.h>
#include <cgogn/core/functions/mesh_info.h>
#include <cgogn/core/types/mesh_traits.h>
#include <cgogn/geometry/algos/offset_axis/input_geometry.h>
#include <cgogn/geometry/algos/offset_axis/normal_utils.h>
#include <cgogn/geometry/types/vector_traits.h>

#include <utility>
#include <vector>

namespace cgogn
{

namespace geometry
{

namespace offset_axis
{

namespace constants
{

// Samples farther than alpha + filter_margin from the input are removed by the sampling filter.
inline constexpr Scalar filter_margin = Scalar(1e-2);

} // namespace constants

// Removes the samples farther than max_distance from the input: the input point cloud if its kd-tree is built,
// otherwise the input surface through its BVH (one of them must be built). Returns the number of removed samples.
template <typename SURFACE, typename POINTS>
uint32 filter_samples_far_from_input(const InputGeometry<SURFACE, POINTS>& input, Scalar max_distance,
									 POINTS& samples,
									 typename mesh_traits<POINTS>::template Attribute<Vec3>* sample_position,
									 typename mesh_traits<POINTS>::template Attribute<Vec3>* sample_normal)
{
	using PVertex = typename mesh_traits<POINTS>::Vertex;

	std::vector<Vec3> kept, kept_normals;
	kept.reserve(nb_cells<PVertex>(samples));
	foreach_cell(samples, [&](PVertex v) {
		const Vec3 p = (*sample_position)[index_of(samples, v)];
		bool within_distance = false;
		if (input.input_kdtree())
		{
			std::pair<uint32, Scalar> nearest;
			within_distance = input.input_kdtree()->find_nn(p, &nearest, max_distance);
		}
		else
		{
			std::pair<uint32, Vec3> closest;
			within_distance = input.surface_bvh()->closest_point(p, &closest, max_distance);
		}
		if (within_distance)
		{
			kept.push_back(p);
			kept_normals.push_back((*sample_normal)[index_of(samples, v)]);
		}
		return true;
	});
	const size_t before = nb_cells<PVertex>(samples);
	if (kept.size() == before)
		return 0;
	clear(samples);
	for (size_t i = 0; i < kept.size(); ++i)
	{
		PVertex v = add_vertex(samples);
		(*sample_position)[index_of(samples, v)] = kept[i];
		(*sample_normal)[index_of(samples, v)] = kept_normals[i];
	}
	return static_cast<uint32>(before - kept.size());
}

// Sets the normal of each sample from the input: the direction from the closest point of the input surface when the
// input point cloud has no kd-tree, otherwise the direction from the nearest input point (its normal when the sample
// is on it), otherwise +z.
template <typename SURFACE, typename POINTS>
void recompute_samples_normals_from_input(InputGeometry<SURFACE, POINTS>& input, const POINTS& samples,
										  const typename mesh_traits<POINTS>::template Attribute<Vec3>& sample_position,
										  typename mesh_traits<POINTS>::template Attribute<Vec3>& sample_normal)
{
	using PVertex = typename mesh_traits<POINTS>::Vertex;

	if (!input.input_kdtree() && input.data_.surface)
	{
		input.build_surface_bvh();
		if (input.surface_bvh())
		{
			foreach_cell(samples, [&](PVertex v) {
				const uint32 idx = index_of(samples, v);
				const Vec3& position = sample_position[idx];
				std::pair<uint32, Vec3> closest;
				Vec3 normal(0, 0, 1);
				if (input.surface_bvh()->closest_point(position, &closest))
					normal = unit_normal_or_z(position - closest.second);
				sample_normal[idx] = normal;
				return true;
			});
			return;
		}
	}
	if (input.input_kdtree() && input.data_.input_points)
	{
		const POINTS& input_points = *input.data_.input_points;
		const auto& input_position = *input.input_position();
		const auto& input_normal = *input.input_normal();
		foreach_cell(samples, [&](PVertex v) {
			const uint32 idx = index_of(samples, v);
			const Vec3& position = sample_position[idx];
			std::pair<uint32, Scalar> nearest;
			Vec3 normal(0, 0, 1);
			if (input.input_kdtree()->find_nn(position, &nearest))
			{
				const uint32 input_idx = index_of(input_points, input.input_vertices()[nearest.first]);
				normal = position - input_position[input_idx];
				if (normal.squaredNorm() < Scalar(1e-12))
					normal = input_normal[input_idx];
				normal = unit_normal_or_z(normal);
			}
			sample_normal[idx] = normal;
			return true;
		});
		return;
	}
	foreach_cell(samples, [&](PVertex v) {
		sample_normal[index_of(samples, v)] = Vec3(0, 0, 1);
		return true;
	});
}

// Step 3: removes the samples farther than alpha + filter_margin from the input (the input point cloud through its
// kd-tree, otherwise the input surface through its BVH). For the CPU fields the normals of the remaining samples are
// recomputed from the input; the neural samples keep the network gradient. Returns the number of removed samples.
template <typename FIELD, typename SURFACE, typename POINTS>
uint32 filter_samples(InputGeometry<SURFACE, POINTS>& input, Scalar alpha, POINTS& samples,
					  typename mesh_traits<POINTS>::template Attribute<Vec3>* sample_position,
					  typename mesh_traits<POINTS>::template Attribute<Vec3>* sample_normal)
{
	using PVertex = typename mesh_traits<POINTS>::Vertex;

	input.build_input_kdtree();
	if (!input.input_kdtree())
		input.build_surface_bvh();
	const uint32 removed =
		filter_samples_far_from_input(input, alpha + constants::filter_margin, samples, sample_position, sample_normal);
	if constexpr (FIELD::recompute_normals_from_input_after_filtering)
	{
		if (removed > 0 && nb_cells<PVertex>(samples) > 0)
			recompute_samples_normals_from_input(input, samples, *sample_position, *sample_normal);
	}
	return removed;
}

} // namespace offset_axis

} // namespace geometry

} // namespace cgogn

#endif // CGOGN_GEOMETRY_ALGOS_OFFSET_AXIS_SAMPLING_FILTER_H_
