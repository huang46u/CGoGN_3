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
#include <cgogn/geometry/types/vector_traits.h>

#include <utility>
#include <vector>

namespace cgogn
{

namespace geometry
{

namespace offset_axis
{

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

} // namespace offset_axis

} // namespace geometry

} // namespace cgogn

#endif // CGOGN_GEOMETRY_ALGOS_OFFSET_AXIS_SAMPLING_FILTER_H_
