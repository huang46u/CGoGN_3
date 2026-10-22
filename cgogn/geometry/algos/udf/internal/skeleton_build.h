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
#ifndef CGOGN_GEOMETRY_ALGOS_UDF_INTERNAL_SKELETON_BUILD_H_
#define CGOGN_GEOMETRY_ALGOS_UDF_INTERNAL_SKELETON_BUILD_H_

#include <cgogn/core/functions/attributes.h>
#include <cgogn/core/functions/mesh_info.h>
#include <cgogn/core/functions/traversals/edge.h>
#include <cgogn/core/functions/traversals/face.h>
#include <cgogn/core/functions/traversals/vertex.h>
#include <cgogn/core/types/cell_marker.h>
#include <cgogn/geometry/algos/udf/internal/skeleton_mesh_ops.h>
#include <cgogn/geometry/algos/udf/internal/skeleton_topology_types.h>

#include <algorithm>
#include <array>
#include <cstddef>
#include <map>
#include <set>
#include <unordered_map>
#include <utility>
#include <vector>

namespace cgogn::geometry::skeleton_topology_detail
{

template <typename POINTS, typename NONMANIFOLD>
struct SkeletonBuild
{
	using Types = SkeletonTopologyTypes<POINTS, NONMANIFOLD>;
	using SphereVertex = typename Types::SphereVertex;
	using SkeletonVertex = typename Types::SkeletonVertex;
	using SkeletonEdge = typename Types::SkeletonEdge;
	using SkeletonFace = typename Types::SkeletonFace;
	using Context = typename Types::Context;
	using Tet = typename Types::Tet;
	using FaceKey = typename Types::FaceKey;
	using edge_hash = typename Types::edge_hash;
	using edge_equal = typename Types::edge_equal;
	using MeshOps = SkeletonMeshOps<POINTS, NONMANIFOLD>;

	static bool compute_skeleton(Context& p)
	{
		p.spheres_optimizer_->compute_sphere_neighbors();
		if (!p.samples_mesh_ || !p.samples_sphere_ || !p.samples_knn_)
			return false;
		clear(*p.skeleton_);
		std::map<FaceKey, SkeletonFace> skeleton_faces_map;
		auto get_face_key = [](uint32 i1, uint32 i2, uint32 i3) -> FaceKey {
			std::array<uint32, 3> key = {i1, i2, i3};
			std::sort(key.begin(), key.end());
			return key;
		};
		std::vector<std::array<uint32, 4>> raw_tets;
		if (p.spheres_skeleton_vertex_)
		{
			foreach_cell(*p.spheres_, [&](SphereVertex pv) -> bool {
				const uint32 sphere_index = index_of(*p.spheres_, pv);
				if (sphere_index != INVALID_INDEX)
					(*p.spheres_skeleton_vertex_)[sphere_index] = SkeletonVertex();
				return true;
			});
		}
		foreach_cell(*p.spheres_, [&](SphereVertex pv) -> bool {
			uint32 pv_index = index_of(*p.spheres_, pv);
			SkeletonVertex nmv = add_vertex(*p.skeleton_);
			const uint32 skeleton_index = index_of(*p.skeleton_, nmv);
			(*p.skeleton_position_)[skeleton_index] = (*p.spheres_position_)[pv_index];
			(*p.skeleton_radius_)[skeleton_index] = (*p.spheres_radius_)[pv_index];
			if (p.spheres_skeleton_vertex_)
				(*p.spheres_skeleton_vertex_)[pv_index] = nmv;
			if (p.skeleton_source_sphere_)
				value<SphereVertex>(*p.skeleton_, p.skeleton_source_sphere_, nmv) = pv;
			return true;
		});

		std::unordered_map<std::pair<uint32, uint32>, SkeletonEdge, edge_hash, edge_equal> edge_indices;

		foreach_cell(*p.spheres_, [&](SphereVertex pv) -> bool {
			uint32 pv_index = index_of(*p.spheres_, pv);
			SkeletonVertex nmv1 =
				p.spheres_skeleton_vertex_ ? (*p.spheres_skeleton_vertex_)[pv_index] : SkeletonVertex();
			const std::vector<SphereVertex>& neighbors = (*p.spheres_neighbor_clusters_)[pv_index];
			for (SphereVertex neighbor : neighbors)
			{
				uint32 n_index = index_of(*p.spheres_, neighbor);
				SkeletonVertex nmv2 =
					p.spheres_skeleton_vertex_ ? (*p.spheres_skeleton_vertex_)[n_index] : SkeletonVertex();
				std::vector<SkeletonVertex> av = adjacent_vertices_through_edge(*p.skeleton_, nmv1);
				if (std::find(av.begin(), av.end(), nmv2) == av.end())
				{
					SkeletonEdge e = add_edge(*p.skeleton_, nmv1, nmv2);
					edge_indices[{index_of(*p.skeleton_, nmv1), index_of(*p.skeleton_, nmv2)}] = e;
				}
			}
			return true;
		});

		foreach_cell(*p.spheres_, [&](SphereVertex pv) -> bool {
			uint32 idx1 = index_of(*p.spheres_, pv);
			SkeletonVertex nmv1 = p.spheres_skeleton_vertex_ ? (*p.spheres_skeleton_vertex_)[idx1] : SkeletonVertex();
			const std::vector<SphereVertex>& n_pv = (*p.spheres_neighbor_clusters_)[idx1];
			for (const SphereVertex& ne1 : n_pv)
			{
				uint32 idx2 = index_of(*p.spheres_, ne1);
				if (idx1 >= idx2)
					continue;
				SkeletonVertex nmv2 =
					p.spheres_skeleton_vertex_ ? (*p.spheres_skeleton_vertex_)[idx2] : SkeletonVertex();
				const std::vector<SphereVertex>& ne_ne1 = (*p.spheres_neighbor_clusters_)[idx2];
				for (const SphereVertex& ne2 : ne_ne1)
				{
					if (!std::binary_search(n_pv.begin(), n_pv.end(), ne2))
						continue;
					uint32 idx3 = index_of(*p.spheres_, ne2);
					if (idx2 >= idx3)
						continue;

					SkeletonVertex nmv3 =
						p.spheres_skeleton_vertex_ ? (*p.spheres_skeleton_vertex_)[idx3] : SkeletonVertex();

					const auto e12 = edge_indices.find({index_of(*p.skeleton_, nmv1), index_of(*p.skeleton_, nmv2)});
					const auto e23 = edge_indices.find({index_of(*p.skeleton_, nmv2), index_of(*p.skeleton_, nmv3)});
					const auto e13 = edge_indices.find({index_of(*p.skeleton_, nmv1), index_of(*p.skeleton_, nmv3)});
					if (e12 == edge_indices.end() || e23 == edge_indices.end() || e13 == edge_indices.end())
						continue;
					std::vector<SkeletonEdge> edges;
					edges.reserve(3);
					edges.push_back(e12->second);
					edges.push_back(e23->second);
					edges.push_back(e13->second);
					SkeletonFace new_face = add_face(*p.skeleton_, edges);

					skeleton_faces_map[get_face_key(idx1, idx2, idx3)] = new_face;

					const std::vector<SphereVertex>& ne_ne2 = (*p.spheres_neighbor_clusters_)[idx3];

					for (const SphereVertex& ne3 : ne_ne2)
					{
						uint32 idx4 = index_of(*p.spheres_, ne3);
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
		p.skeleton_tets_.clear();
		p.skeleton_tets_.reserve(raw_tets.size());
		foreach_cell(*p.skeleton_, [&](SkeletonFace f) -> bool {
			value<std::set<std::size_t>>(*p.skeleton_, p.incident_tets_, f).clear();
			return true;
		});
		std::size_t tet_index = 0;
		for (const auto& rt : raw_tets)
		{
			Tet new_tet;
			uint32 v[4] = {rt[0], rt[1], rt[2], rt[3]};

			FaceKey keys[4] = {get_face_key(v[1], v[2], v[3]), get_face_key(v[0], v[2], v[3]),
							   get_face_key(v[0], v[1], v[3]), get_face_key(v[0], v[1], v[2])};

			for (uint32 i = 0; i < 4; ++i)
			{
				auto it = skeleton_faces_map.find(keys[i]);
				if (it != skeleton_faces_map.end())
				{
					SkeletonFace f = it->second;
					new_tet.faces[i] = f;

					value<std::set<size_t>>(*p.skeleton_, p.incident_tets_, f).insert(tet_index);
				}
			}
			p.skeleton_tets_.insert({tet_index, new_tet});
			tet_index++;
		}
		MeshOps::compute_edge_degree(p);
		return true;
	}
};

} // namespace cgogn::geometry::skeleton_topology_detail

#endif // CGOGN_GEOMETRY_ALGOS_UDF_INTERNAL_SKELETON_BUILD_H_
