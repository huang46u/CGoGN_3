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
#ifndef CGOGN_GEOMETRY_ALGOS_UDF_INTERNAL_NON_MANIFOLD_PRUNE_H_
#define CGOGN_GEOMETRY_ALGOS_UDF_INTERNAL_NON_MANIFOLD_PRUNE_H_

#include <cgogn/core/functions/attributes.h>
#include <cgogn/core/functions/mesh_info.h>
#include <cgogn/core/functions/traversals/edge.h>
#include <cgogn/core/functions/traversals/face.h>
#include <cgogn/core/functions/traversals/vertex.h>
#include <cgogn/core/types/cell_marker.h>
#include <cgogn/geometry/algos/udf/internal/skeleton_mesh_ops.h>
#include <cgogn/geometry/algos/udf/internal/skeleton_topology_types.h>

#include <algorithm>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

namespace cgogn::geometry::skeleton_topology_detail
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
		if (!p.skeleton_)
			return false;

		bool removed_any_faces = false;
		uint32 removed_face_count = 0;

		for (;;)
		{
			std::vector<uint32> face_ids_to_remove;
			foreach_cell(*p.skeleton_, [&](SkeletonFace f) -> bool {
				if (!f.is_valid())
					return true;
				const uint32 face_id = index_of(*p.skeleton_, f);
				if (face_id == INVALID_INDEX)
					return true;

				const std::vector<SkeletonEdge> face_edges = incident_edges(*p.skeleton_, f);
				if (face_edges.size() != 3)
					return true;

				for (const SkeletonEdge& e : face_edges)
				{
					if (!e.is_valid() || incident_faces(*p.skeleton_, e).size() <= 2)
						return true;
				}

				face_ids_to_remove.push_back(face_id);
				return true;
			});

			if (face_ids_to_remove.empty())
				break;

			for (uint32 face_id : face_ids_to_remove)
			{
				const SkeletonFace f = of_index<SkeletonFace>(*p.skeleton_, face_id);
				if (!f.is_valid())
					continue;

				const std::vector<SkeletonEdge> affected_edges = incident_edges(*p.skeleton_, f);
				remove_face(*p.skeleton_, f);
				removed_any_faces = true;
				++removed_face_count;

				std::vector<SkeletonVertex> candidate_vertices;
				candidate_vertices.reserve(affected_edges.size() * 2);
				for (const SkeletonEdge& e : affected_edges)
				{
					if (!e.is_valid())
						continue;
					for (const SkeletonVertex& v : incident_vertices(*p.skeleton_, e))
						candidate_vertices.push_back(v);
					if (incident_faces(*p.skeleton_, e).empty())
						remove_edge(*p.skeleton_, e);
				}
				for (const SkeletonVertex& v : candidate_vertices)
				{
					if (!v.is_valid())
						continue;
					const uint32 vertex_id = index_of(*p.skeleton_, v);
					if (vertex_id == INVALID_INDEX || !incident_edges(*p.skeleton_, v).empty())
						continue;
					MeshOps::remove_skeleton_vertex_and_linked_sphere(p, v);
				}
			}
		}

		p.metrics.residual_fully_nm_removed_faces += removed_face_count;
		if (p.edge_degree_ && removed_any_faces)
			MeshOps::compute_edge_degree(p);

		return true;
	}

	static int32 skeleton_euler_characteristic(const Context& p)
	{
		return static_cast<int32>(nb_cells<SkeletonVertex>(*p.skeleton_)) -
			   static_cast<int32>(nb_cells<SkeletonEdge>(*p.skeleton_)) +
			   static_cast<int32>(nb_cells<SkeletonFace>(*p.skeleton_));
	}

	static void run_residual_sheet_prune(Context& p)
	{
		if (!p.skeleton_ || !p.incident_tets_)
			return;

		p.metrics.residual_chi_before = skeleton_euler_characteristic(p);
		auto record_chi_after = [&]() { p.metrics.residual_chi_after = skeleton_euler_characteristic(p); };
		MeshOps::compute_edge_degree(p);
		if (p.residual_prune_mode == ResidualPruneMode::current)
			prune_fully_non_manifold_triangles(p);

		std::vector<SkeletonFace> faces;
		faces.reserve(nb_cells<SkeletonFace>(*p.skeleton_));
		std::unordered_map<uint32, uint32> face_id_to_index;
		face_id_to_index.reserve(nb_cells<SkeletonFace>(*p.skeleton_));
		foreach_cell(*p.skeleton_, [&](SkeletonFace f) -> bool {
			if (!f.is_valid())
				return true;
			const uint32 face_id = index_of(*p.skeleton_, f);
			if (face_id == INVALID_INDEX)
				return true;
			face_id_to_index.emplace(face_id, static_cast<uint32>(faces.size()));
			faces.push_back(f);
			return true;
		});
		if (faces.empty())
		{
			record_chi_after();
			return;
		}

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
		foreach_cell(*p.skeleton_, [&](SkeletonEdge e) -> bool {
			if (!e.is_valid())
				return true;
			const auto incident = incident_faces(*p.skeleton_, e);
			if (incident.size() != 2)
				return true;
			const auto first = face_id_to_index.find(index_of(*p.skeleton_, incident[0]));
			const auto second = face_id_to_index.find(index_of(*p.skeleton_, incident[1]));
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
			const uint32 face_id = index_of(*p.skeleton_, faces[i]);
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
				const SkeletonFace f = of_index<SkeletonFace>(*p.skeleton_, face_id);
				if (!f.is_valid() || index_of(*p.skeleton_, f) != face_id)
					continue;
				const std::vector<SkeletonVertex> vertices = incident_vertices(*p.skeleton_, f);
				for (const SkeletonVertex& v : vertices)
				{
					if (!v.is_valid())
						continue;
					const uint32 vertex_id = index_of(*p.skeleton_, v);
					if (vertex_id == INVALID_INDEX || !visited_vertex_ids.insert(vertex_id).second)
						continue;
					for (const SkeletonEdge& e : incident_edges(*p.skeleton_, v))
					{
						if (!e.is_valid())
							continue;
						if (incident_faces(*p.skeleton_, e).empty())
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
				const SkeletonFace f = of_index<SkeletonFace>(*p.skeleton_, face_id);
				if (!f.is_valid() || index_of(*p.skeleton_, f) != face_id)
					continue;
				for (const SkeletonEdge& e : incident_edges(*p.skeleton_, f))
				{
					if (e.is_valid() && incident_faces(*p.skeleton_, e).size() == 1)
						return true;
				}
			}
			return false;
		};

		uint32 deleted_sheet_count = 0;
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
			if (delete_by_small_sheet)
			{
				++p.metrics.residual_small_sheets;
				if (keep_by_isolated_edge_guard)
					++p.metrics.residual_small_sheets_guarded;
				else
				{
					const bool has_boundary_edge = sheet_has_boundary_edge(it_faces->second);
					if (has_boundary_edge)
						++p.metrics.residual_small_sheets_with_boundary_edge;
					else
						++p.metrics.residual_small_sheets_without_boundary_edge;
					// a single-face sheet without a boundary edge is a fully non-manifold triangle
					if (!has_boundary_edge &&
						(p.residual_prune_mode == ResidualPruneMode::collapsible_only ||
						 (p.residual_prune_mode == ResidualPruneMode::keep_fully_non_manifold && face_count == 1)))
						delete_sheet = false;
				}
			}
			if (!delete_sheet)
				continue;

			++deleted_sheet_count;
			face_ids_to_delete.reserve(face_ids_to_delete.size() + it_faces->second.size());
			for (uint32 face_id : it_faces->second)
				face_ids_to_delete.insert(face_id);
		}

		if (face_ids_to_delete.empty())
		{
			MeshOps::compute_edge_degree(p);
			if (p.residual_prune_mode == ResidualPruneMode::current)
				prune_fully_non_manifold_triangles(p);
			record_chi_after();
			return;
		}

		RemovalStats deletion_stats;
		if (!MeshOps::delete_skeleton_face_id_set(p, face_ids_to_delete, deletion_stats))
		{
			record_chi_after();
			return;
		}

		p.metrics.residual_removed_sheets += deleted_sheet_count;
		p.metrics.residual_removed_faces += deletion_stats.removed_faces;
		MeshOps::compute_edge_degree(p);
		if (p.residual_prune_mode == ResidualPruneMode::current)
			prune_fully_non_manifold_triangles(p);
		record_chi_after();
	}
};

} // namespace cgogn::geometry::skeleton_topology_detail

#endif // CGOGN_GEOMETRY_ALGOS_UDF_INTERNAL_NON_MANIFOLD_PRUNE_H_
