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
#ifndef CGOGN_GEOMETRY_ALGOS_OFFSET_AXIS_SKELETON_BUILD_H_
#define CGOGN_GEOMETRY_ALGOS_OFFSET_AXIS_SKELETON_BUILD_H_

#include <cgogn/core/functions/attributes.h>
#include <cgogn/core/functions/mesh_info.h>
#include <cgogn/core/functions/traversals/edge.h>
#include <cgogn/core/functions/traversals/face.h>
#include <cgogn/core/functions/traversals/vertex.h>
#include <cgogn/geometry/algos/offset_axis/skeleton_mesh_ops.h>
#include <cgogn/geometry/algos/offset_axis/skeleton_tets.h>
#include <cgogn/geometry/types/vector_traits.h>

#include <algorithm>
#include <array>
#include <cstddef>
#include <functional>
#include <map>
#include <set>
#include <unordered_map>
#include <utility>
#include <vector>

namespace cgogn
{

namespace geometry
{

namespace offset_axis
{

// Hash and equality of an unordered pair of vertex indices.
struct edge_hash
{
	std::size_t operator()(const std::pair<uint32, uint32>& edge) const
	{
		const uint32 a = std::min(edge.first, edge.second);
		const uint32 b = std::max(edge.first, edge.second);
		return std::hash<uint64>()((uint64(a) << 32) | uint64(b));
	}
};

struct edge_equal
{
	bool operator()(const std::pair<uint32, uint32>& edge1, const std::pair<uint32, uint32>& edge2) const
	{
		return ((edge1.first == edge2.first && edge1.second == edge2.second) ||
				(edge1.first == edge2.second && edge1.second == edge2.first));
	}
};

// Builds the skeleton of the spheres: one vertex per sphere, one edge per pair of neighbor spheres, one face per triple
// of pairwise neighbors, and one tet (4 faces, not a mesh cell) per 4 pairwise neighbors. The neighbor lists are sorted.
template <typename POINTS, typename NONMANIFOLD>
void build_skeleton_from_spheres(
	const POINTS& spheres, const typename mesh_traits<POINTS>::template Attribute<Vec3>* sphere_position,
	const typename mesh_traits<POINTS>::template Attribute<Scalar>* sphere_radius,
	const typename mesh_traits<POINTS>::template Attribute<std::vector<typename mesh_traits<POINTS>::Vertex>>*
		sphere_neighbors,
	NONMANIFOLD& skeleton, typename mesh_traits<NONMANIFOLD>::template Attribute<Vec3>* skeleton_position,
	typename mesh_traits<NONMANIFOLD>::template Attribute<Scalar>* skeleton_radius,
	typename mesh_traits<NONMANIFOLD>::template Attribute<typename mesh_traits<POINTS>::Vertex>* skeleton_source_sphere,
	typename mesh_traits<POINTS>::template Attribute<typename mesh_traits<NONMANIFOLD>::Vertex>* sphere_skeleton_vertex,
	typename mesh_traits<NONMANIFOLD>::template Attribute<uint32>* edge_degree,
	SkeletonTets<NONMANIFOLD>& tets)
{
	using SphereVertex = typename mesh_traits<POINTS>::Vertex;
	using SkeletonVertex = typename mesh_traits<NONMANIFOLD>::Vertex;
	using SkeletonEdge = typename mesh_traits<NONMANIFOLD>::Edge;
	using SkeletonFace = typename mesh_traits<NONMANIFOLD>::Face;
	using Tet = typename SkeletonTets<NONMANIFOLD>::Tet;
	using FaceKey = std::array<uint32, 3>;

	clear(skeleton);
	std::map<FaceKey, SkeletonFace> skeleton_faces_map;
	auto get_face_key = [](uint32 i1, uint32 i2, uint32 i3) -> FaceKey {
		std::array<uint32, 3> key = {i1, i2, i3};
		std::sort(key.begin(), key.end());
		return key;
	};
	std::vector<std::array<uint32, 4>> raw_tets;
	foreach_cell(spheres, [&](SphereVertex pv) -> bool {
		uint32 pv_index = index_of(spheres, pv);
		SkeletonVertex nmv = add_vertex(skeleton);
		const uint32 skeleton_index = index_of(skeleton, nmv);
		(*skeleton_position)[skeleton_index] = (*sphere_position)[pv_index];
		(*skeleton_radius)[skeleton_index] = (*sphere_radius)[pv_index];
		(*sphere_skeleton_vertex)[pv_index] = nmv;
		value<SphereVertex>(skeleton, skeleton_source_sphere, nmv) = pv;
		return true;
	});

	std::unordered_map<std::pair<uint32, uint32>, SkeletonEdge, edge_hash, edge_equal> edge_indices;

	foreach_cell(spheres, [&](SphereVertex pv) -> bool {
		uint32 pv_index = index_of(spheres, pv);
		SkeletonVertex nmv1 = (*sphere_skeleton_vertex)[pv_index];
		const std::vector<SphereVertex>& neighbors = (*sphere_neighbors)[pv_index];
		for (SphereVertex neighbor : neighbors)
		{
			uint32 n_index = index_of(spheres, neighbor);
			SkeletonVertex nmv2 = (*sphere_skeleton_vertex)[n_index];
			const std::pair<uint32, uint32> key{index_of(skeleton, nmv1), index_of(skeleton, nmv2)};
			if (edge_indices.find(key) == edge_indices.end())
				edge_indices[key] = add_edge(skeleton, nmv1, nmv2);
		}
		return true;
	});

	foreach_cell(spheres, [&](SphereVertex pv) -> bool {
		uint32 idx1 = index_of(spheres, pv);
		SkeletonVertex nmv1 = (*sphere_skeleton_vertex)[idx1];
		const std::vector<SphereVertex>& n_pv = (*sphere_neighbors)[idx1];
		for (const SphereVertex& ne1 : n_pv)
		{
			uint32 idx2 = index_of(spheres, ne1);
			if (idx1 >= idx2)
				continue;
			SkeletonVertex nmv2 = (*sphere_skeleton_vertex)[idx2];
			const std::vector<SphereVertex>& ne_ne1 = (*sphere_neighbors)[idx2];
			for (const SphereVertex& ne2 : ne_ne1)
			{
				if (!std::binary_search(n_pv.begin(), n_pv.end(), ne2))
					continue;
				uint32 idx3 = index_of(spheres, ne2);
				if (idx2 >= idx3)
					continue;

				SkeletonVertex nmv3 = (*sphere_skeleton_vertex)[idx3];

				// the three edges exist: the three spheres are pairwise neighbors
				std::vector<SkeletonEdge> edges;
				edges.reserve(3);
				edges.push_back(edge_indices.at({index_of(skeleton, nmv1), index_of(skeleton, nmv2)}));
				edges.push_back(edge_indices.at({index_of(skeleton, nmv2), index_of(skeleton, nmv3)}));
				edges.push_back(edge_indices.at({index_of(skeleton, nmv1), index_of(skeleton, nmv3)}));
				SkeletonFace new_face = add_face(skeleton, edges);

				skeleton_faces_map[get_face_key(idx1, idx2, idx3)] = new_face;

				const std::vector<SphereVertex>& ne_ne2 = (*sphere_neighbors)[idx3];

				for (const SphereVertex& ne3 : ne_ne2)
				{
					uint32 idx4 = index_of(spheres, ne3);
					if (idx3 >= idx4)
						continue;
					bool connected_v1 = std::binary_search(n_pv.begin(), n_pv.end(), ne3);
					bool connected_v2 = std::binary_search(ne_ne1.begin(), ne_ne1.end(), ne3);
					if (connected_v1 && connected_v2)
					{
						raw_tets.push_back({idx1, idx2, idx3, idx4});
					}
				}
			}
		}
		return true;
	});

	// Resolve Tets
	tets.clear();
	tets.reserve(raw_tets.size());
	foreach_cell(skeleton, [&](SkeletonFace f) -> bool {
		tets.reset_face(f);
		return true;
	});
	std::size_t tet_index = 0;
	for (const auto& rt : raw_tets)
	{
		Tet new_tet;
		uint32 v[4] = {rt[0], rt[1], rt[2], rt[3]};

		FaceKey keys[4] = {get_face_key(v[1], v[2], v[3]), get_face_key(v[0], v[2], v[3]),
						   get_face_key(v[0], v[1], v[3]), get_face_key(v[0], v[1], v[2])};

		// the four faces exist: every triple of pairwise neighbors with increasing indices got a face
		for (uint32 i = 0; i < 4; ++i)
		{
			const SkeletonFace f = skeleton_faces_map.at(keys[i]);
			new_tet.faces[i] = f;
		}
		tets.add(tet_index, new_tet);
		tet_index++;
	}
	compute_edge_degree(skeleton, edge_degree);
}

} // namespace offset_axis

} // namespace geometry

} // namespace cgogn

#endif // CGOGN_GEOMETRY_ALGOS_OFFSET_AXIS_SKELETON_BUILD_H_
