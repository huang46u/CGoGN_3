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
#include <cgogn/geometry/algos/offset_axis/internal/skeleton_mesh_ops.h>
#include <cgogn/geometry/algos/offset_axis/reconstruction_options.h>

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

// Removes the faces whose 3 edges have more than 2 faces, until there is none.
template <typename POINTS, typename NONMANIFOLD>
void remove_fully_non_manifold_triangles(
	NONMANIFOLD& skeleton, typename mesh_traits<NONMANIFOLD>::template Attribute<uint32>* edge_degree,
	SpheresOptimizer<POINTS>& spheres_optimizer,
	typename mesh_traits<NONMANIFOLD>::template Attribute<typename mesh_traits<POINTS>::Vertex>* skeleton_source_sphere,
	typename mesh_traits<POINTS>::template Attribute<typename mesh_traits<NONMANIFOLD>::Vertex>* sphere_skeleton_vertex)
{
	using SkeletonEdge = typename mesh_traits<NONMANIFOLD>::Edge;
	using SkeletonFace = typename mesh_traits<NONMANIFOLD>::Face;

	bool removed_any_faces = false;
	for (;;)
	{
		std::vector<SkeletonFace> faces_to_remove;
		foreach_cell(skeleton, [&](SkeletonFace f) -> bool {
			for (SkeletonEdge e : incident_edges(skeleton, f))
				if (incident_faces(skeleton, e).size() <= 2)
					return true;
			faces_to_remove.push_back(f);
			return true;
		});

		if (faces_to_remove.empty())
			break;

		for (SkeletonFace f : faces_to_remove)
			remove_faces(skeleton, {f}, spheres_optimizer, skeleton_source_sphere, sphere_skeleton_vertex);
		removed_any_faces = true;
	}

	if (removed_any_faces)
		compute_edge_degree(skeleton, edge_degree);
}

// Step 11: removes the small sheets (at most 2 faces connected by degree-2 edges) that have no vertex on an isolated
// edge, with their tets; the mode keeps some of them and removes the fully non-manifold triangles before and after.
template <typename POINTS, typename NONMANIFOLD>
void prune_residual_sheets(
	NONMANIFOLD& skeleton, SkeletonTets<NONMANIFOLD>& tets,
	typename mesh_traits<NONMANIFOLD>::template Attribute<uint32>* edge_degree,
	SpheresOptimizer<POINTS>& spheres_optimizer,
	typename mesh_traits<NONMANIFOLD>::template Attribute<typename mesh_traits<POINTS>::Vertex>* skeleton_source_sphere,
	typename mesh_traits<POINTS>::template Attribute<typename mesh_traits<NONMANIFOLD>::Vertex>* sphere_skeleton_vertex,
	ResidualPruneOption mode)
{
	using SkeletonVertex = typename mesh_traits<NONMANIFOLD>::Vertex;
	using SkeletonEdge = typename mesh_traits<NONMANIFOLD>::Edge;
	using SkeletonFace = typename mesh_traits<NONMANIFOLD>::Face;

	compute_edge_degree(skeleton, edge_degree);
	if (mode == ResidualPruneOption::current)
		remove_fully_non_manifold_triangles(skeleton, edge_degree, spheres_optimizer, skeleton_source_sphere,
											sphere_skeleton_vertex);

	std::vector<SkeletonFace> faces;
	faces.reserve(nb_cells<SkeletonFace>(skeleton));
	std::unordered_map<uint32, uint32> face_id_to_index;
	face_id_to_index.reserve(nb_cells<SkeletonFace>(skeleton));
	foreach_cell(skeleton, [&](SkeletonFace f) -> bool {
		face_id_to_index.emplace(index_of(skeleton, f), static_cast<uint32>(faces.size()));
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
	foreach_cell(skeleton, [&](SkeletonEdge e) -> bool {
		const std::vector<SkeletonFace> incident = incident_faces(skeleton, e);
		if (incident.size() != 2)
			return true;
		unite(face_id_to_index[index_of(skeleton, incident[0])], face_id_to_index[index_of(skeleton, incident[1])]);
		return true;
	});

	std::unordered_map<uint32, uint32> root_to_sheet_label;
	root_to_sheet_label.reserve(faces.size());
	std::unordered_map<uint32, std::vector<uint32>> label_to_faces;
	label_to_faces.reserve(faces.size());
	for (uint32 i = 0; i < static_cast<uint32>(faces.size()); ++i)
	{
		const auto it =
			root_to_sheet_label.emplace(find_root(i), static_cast<uint32>(root_to_sheet_label.size())).first;
		label_to_faces[it->second].push_back(index_of(skeleton, faces[i]));
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
			for (SkeletonVertex v : incident_vertices(skeleton, of_index<SkeletonFace>(skeleton, face_id)))
			{
				if (!visited_vertex_ids.insert(index_of(skeleton, v)).second)
					continue;
				for (SkeletonEdge e : incident_edges(skeleton, v))
					if (incident_faces(skeleton, e).empty())
						return true;
			}
		}
		return false;
	};

	auto sheet_has_boundary_edge = [&](const std::vector<uint32>& face_ids) -> bool {
		for (uint32 face_id : face_ids)
			for (SkeletonEdge e : incident_edges(skeleton, of_index<SkeletonFace>(skeleton, face_id)))
				if (incident_faces(skeleton, e).size() == 1)
					return true;
		return false;
	};

	std::unordered_set<uint32> face_ids_to_delete;
	for (uint32 sheet_label : sorted_sheet_labels)
	{
		const std::vector<uint32>& sheet_faces = label_to_faces[sheet_label];
		const size_t face_count = sheet_faces.size();
		const bool delete_by_small_sheet = (face_count <= forced_delete_max_face_count);
		const bool keep_by_isolated_edge_guard =
			delete_by_small_sheet && sheet_has_vertex_connected_to_isolated_edge(sheet_faces);
		bool delete_sheet = delete_by_small_sheet && !keep_by_isolated_edge_guard;
		// a single-face sheet without a boundary edge is a fully non-manifold triangle
		if (delete_sheet && !sheet_has_boundary_edge(sheet_faces) &&
			(mode == ResidualPruneOption::collapsible_only ||
			 (mode == ResidualPruneOption::keep_fully_non_manifold && face_count == 1)))
			delete_sheet = false;
		if (!delete_sheet)
			continue;

		face_ids_to_delete.reserve(face_ids_to_delete.size() + sheet_faces.size());
		for (uint32 face_id : sheet_faces)
			face_ids_to_delete.insert(face_id);
	}

	if (!face_ids_to_delete.empty())
		remove_faces_and_tets(skeleton, tets, face_ids_to_delete, spheres_optimizer, skeleton_source_sphere,
							  sphere_skeleton_vertex);

	compute_edge_degree(skeleton, edge_degree);
	if (mode == ResidualPruneOption::current)
		remove_fully_non_manifold_triangles(skeleton, edge_degree, spheres_optimizer, skeleton_source_sphere,
											sphere_skeleton_vertex);
}

} // namespace topology_detail

} // namespace offset_axis

} // namespace geometry

} // namespace cgogn

#endif // CGOGN_GEOMETRY_ALGOS_OFFSET_AXIS_INTERNAL_NON_MANIFOLD_PRUNE_H_
