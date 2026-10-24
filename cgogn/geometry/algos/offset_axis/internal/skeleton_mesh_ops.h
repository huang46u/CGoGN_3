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
#ifndef CGOGN_GEOMETRY_ALGOS_OFFSET_AXIS_INTERNAL_SKELETON_MESH_OPS_H_
#define CGOGN_GEOMETRY_ALGOS_OFFSET_AXIS_INTERNAL_SKELETON_MESH_OPS_H_

#include <cgogn/core/functions/attributes.h>
#include <cgogn/core/functions/mesh_info.h>
#include <cgogn/core/functions/traversals/edge.h>
#include <cgogn/core/functions/traversals/face.h>
#include <cgogn/core/functions/traversals/vertex.h>
#include <cgogn/core/types/cell_marker.h>
#include <cgogn/geometry/algos/offset_axis/internal/skeleton_topology_types.h>

#include <cstddef>
#include <set>
#include <unordered_set>
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
struct SkeletonMeshOps
{
	using Types = SkeletonTopologyTypes<POINTS, NONMANIFOLD>;
	using SphereVertex = typename Types::SphereVertex;
	using SkeletonVertex = typename Types::SkeletonVertex;
	using SkeletonEdge = typename Types::SkeletonEdge;
	using SkeletonFace = typename Types::SkeletonFace;
	using Context = typename Types::Context;
	using RemovalStats = typename Types::RemovalStats;
	using SkeletonFaceDeletionOptions = typename Types::SkeletonFaceDeletionOptions;
	using Tet = typename Types::Tet;

	static bool remove_skeleton_vertex_and_linked_sphere(Context& p, const SkeletonVertex& v)
	{
		if (!p.skeleton_ || !v.is_valid())
			return false;
		const uint32 skeleton_vertex_id = index_of(*p.skeleton_, v);
		if (skeleton_vertex_id == INVALID_INDEX)
			return false;

		SphereVertex linked_sphere;
		if (p.skeleton_source_sphere_)
			linked_sphere = value<SphereVertex>(*p.skeleton_, p.skeleton_source_sphere_, v);

		if (linked_sphere.is_valid() && p.spheres_)
		{
			const uint32 sphere_index = index_of(*p.spheres_, linked_sphere);
			if (sphere_index != INVALID_INDEX && p.spheres_skeleton_vertex_)
				(*p.spheres_skeleton_vertex_)[sphere_index] = SkeletonVertex();
			p.spheres_optimizer_->remove_sphere(linked_sphere);
		}

		if (p.skeleton_source_sphere_)
			value<SphereVertex>(*p.skeleton_, p.skeleton_source_sphere_, v) = SphereVertex();
		remove_vertex(*p.skeleton_, v);
		return true;
	}
	static void compute_edge_degree(Context& p)
	{
		parallel_foreach_cell(*p.skeleton_, [&](SkeletonEdge e) {
			auto in_face = incident_faces(*p.skeleton_, e);
			const uint32 ide = index_of(*p.skeleton_, e);
			const uint32 deg = static_cast<uint32>(in_face.size());

			if (ide != INVALID_INDEX && p.edge_degree_)
				(*p.edge_degree_)[ide] = deg;

			return true;
		});
	}
	static uint32 remove_orphan_edges_from_removed_face_edges(Context& p, const std::vector<SkeletonEdge>& face_edges)
	{
		uint32 removed_edges = 0;
		std::vector<SkeletonVertex> candidate_vertices;
		candidate_vertices.reserve(face_edges.size() * 2);
		for (SkeletonEdge e : face_edges)
		{
			if (!e.is_valid())
				continue;
			const std::vector<SkeletonVertex> edge_vertices = incident_vertices(*p.skeleton_, e);
			for (const SkeletonVertex& v : edge_vertices)
				candidate_vertices.push_back(v);
			if (incident_faces(*p.skeleton_, e).empty())
			{
				remove_edge(*p.skeleton_, e);
				++removed_edges;
			}
		}
		for (const SkeletonVertex& v : candidate_vertices)
		{
			if (!v.is_valid())
				continue;
			const uint32 idv = index_of(*p.skeleton_, v);
			if (idv == INVALID_INDEX)
				continue;
			if (!incident_edges(*p.skeleton_, v).empty())
				continue;
			remove_skeleton_vertex_and_linked_sphere(p, v);
		}
		return removed_edges;
	}
	static std::unordered_set<uint32> collect_current_orphan_edge_ids(Context& p)
	{
		std::unordered_set<uint32> orphan_edge_ids;
		if (!p.skeleton_)
			return orphan_edge_ids;
		orphan_edge_ids.reserve(nb_cells<SkeletonEdge>(*p.skeleton_));
		foreach_cell(*p.skeleton_, [&](SkeletonEdge e) -> bool {
			if (!e.is_valid())
				return true;
			const uint32 ide = index_of(*p.skeleton_, e);
			if (ide == INVALID_INDEX)
				return true;
			if (incident_faces(*p.skeleton_, e).empty())
				orphan_edge_ids.insert(ide);
			return true;
		});
		return orphan_edge_ids;
	}
	static void remove_global_orphan_skeleton_elements(
		Context& p, uint32* out_removed_edges = nullptr,
		const std::unordered_set<uint32>* preserved_orphan_edge_ids = nullptr)
	{
		uint32 removed_edges = 0;
		if (!p.skeleton_)
		{
			if (out_removed_edges)
				*out_removed_edges = 0;
			return;
		}

		std::vector<SkeletonEdge> orphan_edges;
		orphan_edges.reserve(nb_cells<SkeletonEdge>(*p.skeleton_));
		foreach_cell(*p.skeleton_, [&](SkeletonEdge e) -> bool {
			if (!e.is_valid())
				return true;
			if (!incident_faces(*p.skeleton_, e).empty())
				return true;
			const uint32 ide = index_of(*p.skeleton_, e);
			if (ide == INVALID_INDEX)
				return true;
			if (preserved_orphan_edge_ids && preserved_orphan_edge_ids->find(ide) != preserved_orphan_edge_ids->end())
				return true;
			orphan_edges.push_back(e);
			return true;
		});

		std::vector<SkeletonVertex> candidate_vertices;
		candidate_vertices.reserve(orphan_edges.size() * 2);
		for (const SkeletonEdge& e : orphan_edges)
		{
			if (!e.is_valid())
				continue;
			for (const SkeletonVertex& v : incident_vertices(*p.skeleton_, e))
				candidate_vertices.push_back(v);
			remove_edge(*p.skeleton_, e);
			++removed_edges;
		}

		foreach_cell(*p.skeleton_, [&](SkeletonVertex v) -> bool {
			if (!v.is_valid())
				return true;
			candidate_vertices.push_back(v);
			return true;
		});

		std::unordered_set<uint32> seen_vertex_ids;
		seen_vertex_ids.reserve(candidate_vertices.size());

		for (const SkeletonVertex& v : candidate_vertices)
		{
			if (!v.is_valid())
				continue;
			const uint32 idv = index_of(*p.skeleton_, v);
			if (idv == INVALID_INDEX)
				continue;
			if (!seen_vertex_ids.insert(idv).second)
				continue;
			if (!incident_edges(*p.skeleton_, v).empty())
				continue;
			remove_skeleton_vertex_and_linked_sphere(p, v);
		}

		if (out_removed_edges)
			*out_removed_edges = removed_edges;
	}
	static bool delete_skeleton_face_id_set(Context& p, const std::unordered_set<uint32>& face_ids_to_delete,
											RemovalStats& out_stats, const SkeletonFaceDeletionOptions& options = {})
	{
		out_stats = RemovalStats{};
		if (!p.skeleton_ || !p.incident_tets_ || face_ids_to_delete.empty())
			return false;
		const std::unordered_set<uint32> preexisting_orphan_edge_ids = options.remove_global_orphan_elements_after_batch
																		   ? collect_current_orphan_edge_ids(p)
																		   : std::unordered_set<uint32>{};

		for (uint32 face_id : face_ids_to_delete)
		{
			const SkeletonFace face = of_index<SkeletonFace>(*p.skeleton_, face_id);
			if (!face.is_valid())
				continue;
			const uint32 current_face_id = index_of(*p.skeleton_, face);
			if (current_face_id == INVALID_INDEX || current_face_id != face_id)
				continue;

			bool face_still_attached = false;
			for (SkeletonEdge e : incident_edges(*p.skeleton_, face))
			{
				if (!e.is_valid())
					continue;
				for (SkeletonFace ef : incident_faces(*p.skeleton_, e))
				{
					if (ef == face)
					{
						face_still_attached = true;
						break;
					}
				}
				if (face_still_attached)
					break;
			}
			if (!face_still_attached)
				continue;

			const std::set<std::size_t> in_tets = (*p.incident_tets_)[face_id];
			std::vector<SkeletonEdge> affected_edges = incident_edges(*p.skeleton_, face);
			remove_face(*p.skeleton_, face);
			++out_stats.removed_faces;

			if (options.remove_incident_orphan_edges_immediately)
				out_stats.removed_edges += remove_orphan_edges_from_removed_face_edges(p, affected_edges);

			for (std::size_t tet_id : in_tets)
			{
				auto it_tet = p.skeleton_tets_.find(tet_id);
				if (it_tet == p.skeleton_tets_.end())
					continue;
				const Tet old_tet = it_tet->second;
				for (uint32 i = 0; i < 4; ++i)
				{
					const SkeletonFace f = old_tet.faces[i];
					if (!f.is_valid() || f == face)
						continue;
					const uint32 idf = index_of(*p.skeleton_, f);
					if (idf != INVALID_INDEX)
						(*p.incident_tets_)[idf].erase(tet_id);
				}
				p.skeleton_tets_.erase(it_tet);
				++out_stats.removed_tets;
			}
		}

		if (out_stats.removed_faces == 0)
			return false;

		if (options.remove_global_orphan_elements_after_batch)
		{
			uint32 global_orphan_edges = 0;
			remove_global_orphan_skeleton_elements(p, &global_orphan_edges, &preexisting_orphan_edge_ids);
			out_stats.removed_edges += global_orphan_edges;
		}
		return true;
	}
	static bool find_skeleton_edge_between(Context& p, const SkeletonVertex& a, const SkeletonVertex& b,
										   SkeletonEdge& out_edge)
	{
		out_edge = SkeletonEdge();
		if (!p.skeleton_ || !a.is_valid() || !b.is_valid() || a == b)
			return false;
		for (SkeletonEdge edge : incident_edges(*p.skeleton_, a))
		{
			if (!edge.is_valid())
				continue;
			const std::vector<SkeletonVertex> edge_vertices = incident_vertices(*p.skeleton_, edge);
			if (edge_vertices.size() == 2 &&
				((edge_vertices[0] == a && edge_vertices[1] == b) || (edge_vertices[0] == b && edge_vertices[1] == a)))
			{
				out_edge = edge;
				return true;
			}
		}
		return false;
	}
};

} // namespace topology_detail

} // namespace offset_axis

} // namespace geometry

} // namespace cgogn

#endif // CGOGN_GEOMETRY_ALGOS_OFFSET_AXIS_INTERNAL_SKELETON_MESH_OPS_H_
