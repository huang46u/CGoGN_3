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
#ifndef CGOGN_GEOMETRY_ALGOS_OFFSET_AXIS_SPATIAL_QUERY_H_
#define CGOGN_GEOMETRY_ALGOS_OFFSET_AXIS_SPATIAL_QUERY_H_

#include <cgogn/core/utils/numerics.h>
#include <cgogn/geometry/types/vector_traits.h>

#include <libacc/bvh_tree.h>
#include <libacc/kd_tree.h>

#include <utility>

namespace cgogn
{

namespace geometry
{

namespace offset_axis
{

inline bool query_surface_closest_point(const acc::BVHTree<uint32, Vec3>& bvh, const Vec3& query,
											uint32& primitive_index, Vec3& closest_point, Scalar& distance)
{
	std::pair<uint32, Vec3> result;
	if (!bvh.closest_point(query, &result))
		return false;

	primitive_index = result.first;
	closest_point = result.second;
	distance = (query - result.second).norm();
	return true;
}

inline bool query_point_cloud_nearest_point(const acc::KDTree<3, uint32>& kdtree, const Vec3& query,
											uint32& point_index, Scalar& distance)
{
	std::pair<uint32, Scalar> result;
	if (!kdtree.find_nn(query, &result))
		return false;

	point_index = result.first;
	distance = result.second;
	return true;
}

} // namespace offset_axis

} // namespace geometry

} // namespace cgogn

#endif // CGOGN_GEOMETRY_ALGOS_OFFSET_AXIS_SPATIAL_QUERY_H_
