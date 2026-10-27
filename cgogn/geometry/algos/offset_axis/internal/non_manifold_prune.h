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
#ifndef CGOGN_GEOMETRY_ALGOS_OFFSET_AXIS_INTERNAL_NON_MANIFOLD_PRUNE_H_
#define CGOGN_GEOMETRY_ALGOS_OFFSET_AXIS_INTERNAL_NON_MANIFOLD_PRUNE_H_

#include <cgogn/core/functions/attributes.h>
#include <cgogn/core/functions/mesh_info.h>
#include <cgogn/core/functions/traversals/edge.h>
#include <cgogn/core/functions/traversals/face.h>
#include <cgogn/core/functions/traversals/vertex.h>
#include <cgogn/core/types/cell_marker.h>
#include <cgogn/geometry/algos/offset_axis/internal/skeleton_mesh_ops.h>
#include <cgogn/geometry/algos/offset_axis/internal/skeleton_topology_types.h>

#include <algorithm>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

namespace cgogn
{

namespace geometry
{

namespace offset_axis
{

namespace topology_detail
{

template <typename POINTS, typename NONMANIFOLD>
struct NonManifoldPrune
{
	using Types = SkeletonTopologyTypes<POINTS, NONMANIFOLD>;
	using SkeletonVertex = typename Types::SkeletonVertex;
	using SkeletonEdge = typename Types::SkeletonEdge;
	using SkeletonFace = typename Types::SkeletonFace;
	using Context = typename Types::Context;
	using RemovalStats = typename Types::RemovalStats;
	using ResidualPruneMode = typename Types::ResidualPruneMode;
	using MeshOps = SkeletonMeshOps<POINTS, NONMANIFOLD>;

	static bool prune_fully_non_manifold_triangles(Context& p)
	{
		bool removed_any_faces = false;

		for (;;)
		{
			std::vector<uint32> face_ids_to_remove;
			foreach_cell(p.data.skeleton, [&](SkeletonFace f) -> bool {
				if (!f.is_valid())
					return true;
				const uint32 face_id = index_of(p.data.skeleton, f);
				if (face_id == INVALID_INDEX)
					return true;

				const std::vector<SkeletonEdge> face_edges = incident_edges(p.data.skeleton, f);
				if (face_edges.size() != 3)
					return true;

				for (const SkeletonEdge& e : face_edges)
				{
					if (!e.is_valid() || incident_faces(p.data.skeleton, e).size() <= 2)
						return true;
				}

				face_ids_to_remove.push_back(face_id);
				return true;
			});

			if (face_ids_to_remove.empty())
				break;

			for (uint32 face_id : face_ids_to_remove)
			{
				const SkeletonFace f = of_index<SkeletonFace>(p.data.skeleton, face_id);
				if (!f.is_valid())
					continue;

				const std::vector<SkeletonEdge> affected_edges = incident_edges(p.data.skeleton, f);
				remove_face(p.data.skeleton, f);
				removed_any_faces = true;
				MeshOps::remove_orphan_edges_from_removed_face_edges(p, affected_edges);
			}
		}

		if (removed_any_faces)
			MeshOps::compute_edge_degree(p);

		return true;
	}

	static void run_residual_sheet_prune(Context& p, ResidualPruneMode mode)
	{
		MeshOps::compute_edge_degree(p);
		if (mode == ResidualPruneMode::current)
			prune_fully_non_manifold_triangles(p);

		std::vector<SkeletonFace> faces;
		faces.reserve(nb_cells<SkeletonFace>(p.data.skeleton));
		std::unordered_map<uint32, uint32> face_id_to_index;
		face_id_to_index.reserve(nb_cells<SkeletonFace>(p.data.skeleton));
		foreach_cell(p.data.skeleton, [&](SkeletonFace f) -> bool {
			if (!f.is_valid())
				return true;
			const uint32 face_id = index_of(p.data.skeleton, f);
			if (face_id == INVALID_INDEX)
				return true;
			face_id_to_index.emplace(face_id, static_cast<uint32>(faces.size()));
			faces.push_back(f);
			return true;
		});
		if (faces.empty())
			return;

		std::vector<uint32> parents(faces.size());
		std::vector<uint32> ranks(faces.size(), 0);
		for (uint32 i = 0; i < static_cast<uint32>(faces.size()); ++i)
			parents[i] = i;
		auto find_root = [&](uint32 index) {
			uint32 root = index;
			while (parents[root] != root)
				root = parents[root];
			while (parents[index] != index)
			{
				const uint32 next = parents[index];
				parents[index] = root;
				index = next;
			}
			return root;
		};
		auto unite = [&](uint32 first, uint32 second) {
			uint32 first_root = find_root(first);
			uint32 second_root = find_root(second);
			if (first_root == second_root)
				return;
			if (ranks[first_root] < ranks[second_root])
				std::swap(first_root, second_root);
			parents[second_root] = first_root;
			if (ranks[first_root] == ranks[second_root])
				++ranks[first_root];
		};
		foreach_cell(p.data.skeleton, [&](SkeletonEdge e) -> bool {
			if (!e.is_valid())
				return true;
			const auto incident = incident_faces(p.data.skeleton, e);
			if (incident.size() != 2)
				return true;
			const auto first = face_id_to_index.find(index_of(p.data.skeleton, incident[0]));
			const auto second = face_id_to_index.find(index_of(p.data.skeleton, incident[1]));
			if (first != face_id_to_index.end() && second != face_id_to_index.end())
				unite(first->second, second->second);
			return true;
		});

		std::unordered_map<uint32, uint32> root_to_sheet_label;
		root_to_sheet_label.reserve(faces.size());
		std::unordered_map<uint32, std::vector<uint32>> label_to_faces;
		label_to_faces.reserve(faces.size());
		for (uint32 i = 0; i < static_cast<uint32>(faces.size()); ++i)
		{
			const uint32 root = find_root(i);
			const auto [it, inserted] =
				root_to_sheet_label.emplace(root, static_cast<uint32>(root_to_sheet_label.size()));
			(void)inserted;
			const uint32 face_id = index_of(p.data.skeleton, faces[i]);
			if (face_id != INVALID_INDEX)
				label_to_faces[it->second].push_back(face_id);
		}

		std::vector<uint32> sorted_sheet_labels;
		sorted_sheet_labels.reserve(label_to_faces.size());
		for (const auto& kv : label_to_faces)
			sorted_sheet_labels.push_back(kv.first);
		std::sort(sorted_sheet_labels.begin(), sorted_sheet_labels.end());

		const size_t forced_delete_max_face_count = 2;
		auto sheet_has_vertex_connected_to_isolated_edge = [&](const std::vector<uint32>& face_ids) -> bool {
			std::unordered_set<uint32> visited_vertex_ids;
			visited_vertex_ids.reserve(face_ids.size() * 3);
			for (uint32 face_id : face_ids)
			{
				if (face_id == INVALID_INDEX)
					continue;
				const SkeletonFace f = of_index<SkeletonFace>(p.data.skeleton, face_id);
				if (!f.is_valid() || index_of(p.data.skeleton, f) != face_id)
					continue;
				const std::vector<SkeletonVertex> vertices = incident_vertices(p.data.skeleton, f);
				for (const SkeletonVertex& v : vertices)
				{
					if (!v.is_valid())
						continue;
					const uint32 vertex_id = index_of(p.data.skeleton, v);
					if (vertex_id == INVALID_INDEX || !visited_vertex_ids.insert(vertex_id).second)
						continue;
					for (const SkeletonEdge& e : incident_edges(p.data.skeleton, v))
					{
						if (!e.is_valid())
							continue;
						if (incident_faces(p.data.skeleton, e).empty())
							return true;
					}
				}
			}
			return false;
		};

		auto sheet_has_boundary_edge = [&](const std::vector<uint32>& face_ids) -> bool {
			for (uint32 face_id : face_ids)
			{
				if (face_id == INVALID_INDEX)
					continue;
				const SkeletonFace f = of_index<SkeletonFace>(p.data.skeleton, face_id);
				if (!f.is_valid() || index_of(p.data.skeleton, f) != face_id)
					continue;
				for (const SkeletonEdge& e : incident_edges(p.data.skeleton, f))
				{
					if (e.is_valid() && incident_faces(p.data.skeleton, e).size() == 1)
						return true;
				}
			}
			return false;
		};

		std::unordered_set<uint32> face_ids_to_delete;
		for (uint32 sheet_label : sorted_sheet_labels)
		{
			const auto it_faces = label_to_faces.find(sheet_label);
			if (it_faces == label_to_faces.end() || it_faces->second.empty())
				continue;

			const size_t face_count = it_faces->second.size();
			const bool delete_by_small_sheet = (face_count <= forced_delete_max_face_count);
			const bool keep_by_isolated_edge_guard =
				delete_by_small_sheet && sheet_has_vertex_connected_to_isolated_edge(it_faces->second);
			bool delete_sheet = delete_by_small_sheet && !keep_by_isolated_edge_guard;
			// a single-face sheet without a boundary edge is a fully non-manifold triangle
			if (delete_sheet && !sheet_has_boundary_edge(it_faces->second) &&
				(mode == ResidualPruneMode::collapsible_only ||
				 (mode == ResidualPruneMode::keep_fully_non_manifold && face_count == 1)))
				delete_sheet = false;
			if (!delete_sheet)
				continue;

			face_ids_to_delete.reserve(face_ids_to_delete.size() + it_faces->second.size());
			for (uint32 face_id : it_faces->second)
				face_ids_to_delete.insert(face_id);
		}

		if (face_ids_to_delete.empty())
		{
			MeshOps::compute_edge_degree(p);
			if (mode == ResidualPruneMode::current)
				prune_fully_non_manifold_triangles(p);
			return;
		}

		RemovalStats deletion_stats;
		if (!MeshOps::delete_skeleton_face_id_set(p, face_ids_to_delete, deletion_stats))
			return;

		MeshOps::compute_edge_degree(p);
		if (mode == ResidualPruneMode::current)
			prune_fully_non_manifold_triangles(p);
	}
};

} // namespace topology_detail

} // namespace offset_axis

} // namespace geometry

} // namespace cgogn

#endif // CGOGN_GEOMETRY_ALGOS_OFFSET_AXIS_INTERNAL_NON_MANIFOLD_PRUNE_H_
