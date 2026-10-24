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
#ifndef CGOGN_GEOMETRY_ALGOS_OFFSET_AXIS_RAY_LEVEL_SET_SAMPLER_TRAITS_H_
#define CGOGN_GEOMETRY_ALGOS_OFFSET_AXIS_RAY_LEVEL_SET_SAMPLER_TRAITS_H_

#include <cgogn/core/functions/mesh_info.h>
#include <cgogn/geometry/algos/offset_axis/spatial_query.h>
#include <cgogn/geometry/types/vector_traits.h>

#include <libacc/bvh_tree.h>
#include <libacc/kd_tree.h>

#include <cmath>
#include <limits>
#include <utility>
#include <vector>

namespace cgogn
{

namespace geometry
{

namespace offset_axis
{

struct RaySamplerNeural
{
};

struct RaySamplerSurface
{
};

struct RaySamplerPointCloud
{
};

template <typename RaySamplerTag, typename SURFACE, typename POINTS>
struct RaySamplerTraits
{
};

template <typename SURFACE>
struct SurfaceMeshRayTraits
{
	using Face = typename mesh_traits<SURFACE>::Face;
	template <typename T>
	using Attribute = typename mesh_traits<SURFACE>::template Attribute<T>;

	SURFACE* mesh = nullptr;
	Attribute<Vec3>* position = nullptr;
	acc::BVHTree<uint32, Vec3>* bvh = nullptr;
	std::vector<Face>* bvh_faces = nullptr;

	bool is_ready() const
	{
		return bvh != nullptr;
	}

	Scalar eval_distance(const Vec3& pos) const
	{
		if (!bvh)
			return std::numeric_limits<Scalar>::max();
		uint32 primitive_index;
		Vec3 closest_point;
		Scalar distance;
		if (!query_surface_closest_point(*bvh, pos, primitive_index, closest_point, distance))
			return std::numeric_limits<Scalar>::max();
		(void)primitive_index;
		(void)closest_point;
		return distance;
	}
};

template <typename SURFACE, typename POINTS>
struct RaySamplerTraits<RaySamplerSurface, SURFACE, POINTS>
{
	using SamplerTraits = SurfaceMeshRayTraits<SURFACE>;
	template <typename T>
	using Attribute = typename mesh_traits<SURFACE>::template Attribute<T>;
	using Face = typename mesh_traits<SURFACE>::Face;

	static SamplerTraits make(SURFACE& mesh, Attribute<Vec3>* position, acc::BVHTree<uint32, Vec3>* bvh,
							  std::vector<Face>* bvh_faces)
	{
		SamplerTraits traits;
		traits.mesh = &mesh;
		traits.position = position;
		traits.bvh = bvh;
		traits.bvh_faces = bvh_faces;
		return traits;
	}
};

template <typename POINTS>
struct PointCloudRayTraits
{
	using Vertex = typename mesh_traits<POINTS>::Vertex;
	template <typename T>
	using Attribute = typename mesh_traits<POINTS>::template Attribute<T>;

	POINTS* mesh = nullptr;
	Attribute<Vec3>* position = nullptr;
	Attribute<Vec3>* normal = nullptr;
	Attribute<std::vector<Vertex>>* knn = nullptr;
	acc::KDTree<3, uint32>* kdtree = nullptr;
	const std::vector<Vertex>* kdtree_vertices = nullptr;

	bool is_ready() const
	{
		return mesh && position && normal && knn && kdtree && kdtree_vertices;
	}

	Scalar eval_distance(const Vec3& pos) const
	{
		if (!kdtree)
			return std::numeric_limits<Scalar>::max();
		uint32 point_index;
		Scalar distance;
		if (!query_point_cloud_nearest_point(*kdtree, pos, point_index, distance))
			return std::numeric_limits<Scalar>::max();
		(void)point_index;
		return distance;

	}
};

template <typename SURFACE, typename POINTS>
struct RaySamplerTraits<RaySamplerPointCloud, SURFACE, POINTS>
{
	using SamplerTraits = PointCloudRayTraits<POINTS>;
	using Vertex = typename mesh_traits<POINTS>::Vertex;
	template <typename T>
	using Attribute = typename mesh_traits<POINTS>::template Attribute<T>;

	static SamplerTraits make(POINTS& mesh, Attribute<Vec3>* position, Attribute<Vec3>* normal,
							  Attribute<std::vector<Vertex>>* knn, acc::KDTree<3, uint32>* kdtree,
							  const std::vector<Vertex>* kdtree_vertices)
	{
		SamplerTraits traits;
		traits.mesh = &mesh;
		traits.position = position;
		traits.normal = normal;
		traits.knn = knn;
		traits.kdtree = kdtree;
		traits.kdtree_vertices = kdtree_vertices;
		return traits;
	}
};

} // namespace offset_axis

} // namespace geometry

} // namespace cgogn

#endif // CGOGN_GEOMETRY_ALGOS_OFFSET_AXIS_RAY_LEVEL_SET_SAMPLER_TRAITS_H_
