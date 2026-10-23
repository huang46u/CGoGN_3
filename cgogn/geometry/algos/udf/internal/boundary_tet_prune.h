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
#ifndef CGOGN_GEOMETRY_ALGOS_UDF_INTERNAL_BOUNDARY_TET_PRUNE_H_
#define CGOGN_GEOMETRY_ALGOS_UDF_INTERNAL_BOUNDARY_TET_PRUNE_H_

#include <cgogn/core/functions/attributes.h>
#include <cgogn/core/functions/mesh_info.h>
#include <cgogn/core/functions/traversals/edge.h>
#include <cgogn/core/functions/traversals/face.h>
#include <cgogn/core/functions/traversals/vertex.h>
#include <cgogn/core/types/cell_marker.h>
#include <cgogn/geometry/algos/udf/internal/quadrature_scores.h>
#include <cgogn/geometry/algos/udf/internal/skeleton_topology_types.h>

#include <array>
#include <cstddef>
#include <unordered_set>
#include <vector>

namespace cgogn::geometry::skeleton_topology_detail
{

template <typename POINTS, typename NONMANIFOLD>
struct BoundaryTetPrune
{
	using Types = SkeletonTopologyTypes<POINTS, NONMANIFOLD>;
	using SkeletonVertex = typename Types::SkeletonVertex;
	using SkeletonEdge = typename Types::SkeletonEdge;
	using SkeletonFace = typename Types::SkeletonFace;
	using Context = typename Types::Context;
	using Tet = typename Types::Tet;
	using RemovalStats = typename Types::RemovalStats;
	using TopologyFixScoreCache = typename Types::TopologyFixScoreCache;
	using QuadratureScores = skeleton_topology_detail::QuadratureScores<POINTS, NONMANIFOLD>;

	struct BoundaryEdgeChoice
	{
		SkeletonEdge edge;
		SkeletonFace f0;
		SkeletonFace f1;
		Scalar length = Scalar(-1);
	};

	// A tet is a boundary tet when one of its faces has at least two degree-2 edges.
	static bool is_boundary_tet(Context& p, const Tet& tet)
	{
		for (uint32 i = 0; i < 4; ++i)
		{
			const SkeletonFace f = tet.faces[i];
			uint32 deg2_edge_count = 0;
			for (SkeletonEdge e : incident_edges(*p.skeleton_, f))
			{
				if (incident_faces(*p.skeleton_, e).size() == 2)
					++deg2_edge_count;
			}
			if (deg2_edge_count >= 2)
				return true;
		}
		return false;
	}

	// Longest degree-2 edge among the faces of the tet. Returns false when there is none.
	static bool choose_longest_boundary_edge(Context& p, const Tet& tet, BoundaryEdgeChoice& best_choice)
	{
		std::unordered_set<uint32> seen_edges;
		seen_edges.reserve(8);

		for (uint32 i = 0; i < 4; ++i)
		{
			const SkeletonFace f = tet.faces[i];
			for (SkeletonEdge e : incident_edges(*p.skeleton_, f))
			{
				const uint32 ide = index_of(*p.skeleton_, e);
				if (!seen_edges.insert(ide).second)
					continue;

				const auto in_faces = incident_faces(*p.skeleton_, e);
				if (in_faces.size() != 2)
					continue;
				const SkeletonFace f0 = in_faces[0];
				const SkeletonFace f1 = in_faces[1];
				const uint32 idf0 = index_of(*p.skeleton_, f0);
				const uint32 idf1 = index_of(*p.skeleton_, f1);
				const std::vector<SkeletonVertex> edge_vertices = incident_vertices(*p.skeleton_, e);
				if (edge_vertices.size() != 2)
					continue;
				const Vec3 p0 = value<Vec3>(*p.skeleton_, p.skeleton_position_, edge_vertices[0]);
				const Vec3 p1 = value<Vec3>(*p.skeleton_, p.skeleton_position_, edge_vertices[1]);
				const Scalar edge_len = (p1 - p0).norm();
				if (edge_len > best_choice.length)
				{
					best_choice.edge = e;
					best_choice.f0 = f0;
					best_choice.f1 = f1;
					best_choice.length = edge_len;
				}
			}
		}

		return best_choice.length >= Scalar(0);
	}

	// Read-only: collects the faces of every current boundary tet and the longest degree-2 edge of each such tet.
	static void collect_boundary_tet_mask(Context& p, std::unordered_set<uint32>& face_ids,
										  std::unordered_set<uint32>& edge_ids)
	{
		for (const auto& kv : p.skeleton_tets_)
		{
			const Tet& tet = kv.second;
			if (!is_boundary_tet(p, tet))
				continue;
			for (uint32 i = 0; i < 4; ++i)
				face_ids.insert(index_of(*p.skeleton_, tet.faces[i]));
			BoundaryEdgeChoice best_choice;
			if (choose_longest_boundary_edge(p, tet, best_choice))
				edge_ids.insert(index_of(*p.skeleton_, best_choice.edge));
		}
	}


	static RemovalStats run_boundary_tet_face_deletion(Context& p, TopologyFixScoreCache* score_cache = nullptr)
	{
		RemovalStats stats;

		struct BoundaryTetCandidate
		{
			std::size_t tet_id = std::size_t(-1);
			SkeletonEdge edge;
			SkeletonFace f0;
			SkeletonFace f1;
		};

		std::vector<std::size_t> tet_ids;
		tet_ids.reserve(p.skeleton_tets_.size());
		for (const auto& kv : p.skeleton_tets_)
			tet_ids.push_back(kv.first);

		std::vector<BoundaryTetCandidate> candidates;
		candidates.reserve(tet_ids.size());

		// Phase 1: scan all current tets once and collect boundary-tet deletion candidates.
		for (std::size_t tet_id : tet_ids)
		{
			auto it_tet = p.skeleton_tets_.find(tet_id);
			if (it_tet == p.skeleton_tets_.end())
				continue;
			const Tet tet = it_tet->second;

			if (!is_boundary_tet(p, tet))
				continue;

			BoundaryEdgeChoice best_choice;
			if (choose_longest_boundary_edge(p, tet, best_choice))
				candidates.push_back({tet_id, best_choice.edge, best_choice.f0, best_choice.f1});
		}

		// Phase 2: delete using precomputed candidates only (no boundary re-check during deletion).
		for (const BoundaryTetCandidate& cand : candidates)
		{
			auto it_remove_tet = p.skeleton_tets_.find(cand.tet_id);
			if (it_remove_tet == p.skeleton_tets_.end())
				continue;
			const Tet old_tet = it_remove_tet->second;

			const std::array<SkeletonFace, 2> faces_to_remove = {cand.f0, cand.f1};
			std::vector<SkeletonEdge> affected_edges;
			affected_edges.reserve(6);
			for (const SkeletonFace& f_del : faces_to_remove)
			{
				if (!f_del.is_valid())
					continue;
				for (SkeletonEdge e : incident_edges(*p.skeleton_, f_del))
					affected_edges.push_back(e);
			}
			std::unordered_set<uint32> removed_face_ids;
			removed_face_ids.reserve(faces_to_remove.size() * 2 + 1);
			for (const SkeletonFace& f_del : faces_to_remove)
			{
				if (const uint32 idf_del = index_of(*p.skeleton_, f_del); idf_del != INVALID_INDEX)
					removed_face_ids.insert(idf_del);
				remove_face(*p.skeleton_, f_del);
				++stats.removed_faces;
			}
			QuadratureScores::erase_topology_fix_scores_for_removed_faces_and_orphan_edges(
				p, score_cache, removed_face_ids, affected_edges);
			if (cand.edge.is_valid())
			{
				const auto edge_in_faces = incident_faces(*p.skeleton_, cand.edge);
				if (edge_in_faces.empty())
				{
					remove_edge(*p.skeleton_, cand.edge);
					++stats.removed_edges;
				}
			}
			p.skeleton_tets_.erase(it_remove_tet);
			++stats.removed_tets;
			for (uint32 i = 0; i < 4; ++i)
			{
				const SkeletonFace f = old_tet.faces[i];
				const uint32 idf = index_of(*p.skeleton_, f);
				(*p.incident_tets_)[idf].erase(cand.tet_id);
			}
		}

		return stats;
	}
};

} // namespace cgogn::geometry::skeleton_topology_detail

#endif // CGOGN_GEOMETRY_ALGOS_UDF_INTERNAL_BOUNDARY_TET_PRUNE_H_
