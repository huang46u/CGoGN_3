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
#ifndef CGOGN_GEOMETRY_ALGOS_OFFSET_AXIS_INTERNAL_BOUNDARY_TET_PRUNE_H_
#define CGOGN_GEOMETRY_ALGOS_OFFSET_AXIS_INTERNAL_BOUNDARY_TET_PRUNE_H_

#include <cgogn/core/functions/attributes.h>
#include <cgogn/core/functions/mesh_info.h>
#include <cgogn/core/functions/traversals/edge.h>
#include <cgogn/core/functions/traversals/face.h>
#include <cgogn/core/functions/traversals/vertex.h>
#include <cgogn/core/types/cell_marker.h>
#include <cgogn/geometry/algos/offset_axis/internal/quadrature_scores.h>
#include <cgogn/geometry/algos/offset_axis/internal/skeleton_mesh_ops.h>
#include <cgogn/geometry/algos/offset_axis/internal/skeleton_topology_types.h>

#include <array>
#include <cstddef>
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
struct BoundaryTetPrune
{
	using Types = SkeletonTopologyTypes<POINTS, NONMANIFOLD>;
	using SkeletonVertex = typename Types::SkeletonVertex;
	using SkeletonEdge = typename Types::SkeletonEdge;
	using SkeletonFace = typename Types::SkeletonFace;
	using Context = typename Types::Context;
	using Tet = typename Types::Tet;
	using MeshOps = SkeletonMeshOps<POINTS, NONMANIFOLD>;

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
			for (SkeletonEdge e : incident_edges(p.data.skeleton, f))
			{
				if (incident_faces(p.data.skeleton, e).size() == 2)
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
			for (SkeletonEdge e : incident_edges(p.data.skeleton, f))
			{
				const uint32 ide = index_of(p.data.skeleton, e);
				if (!seen_edges.insert(ide).second)
					continue;

				const auto in_faces = incident_faces(p.data.skeleton, e);
				if (in_faces.size() != 2)
					continue;
				const SkeletonFace f0 = in_faces[0];
				const SkeletonFace f1 = in_faces[1];
				const uint32 idf0 = index_of(p.data.skeleton, f0);
				const uint32 idf1 = index_of(p.data.skeleton, f1);
				const std::vector<SkeletonVertex> edge_vertices = incident_vertices(p.data.skeleton, e);
				if (edge_vertices.size() != 2)
					continue;
				const Vec3 p0 = value<Vec3>(p.data.skeleton, p.data.skeleton_position, edge_vertices[0]);
				const Vec3 p1 = value<Vec3>(p.data.skeleton, p.data.skeleton_position, edge_vertices[1]);
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
		for (const auto& kv : p.skeleton_tets)
		{
			const Tet& tet = kv.second;
			if (!is_boundary_tet(p, tet))
				continue;
			for (uint32 i = 0; i < 4; ++i)
				face_ids.insert(index_of(p.data.skeleton, tet.faces[i]));
			BoundaryEdgeChoice best_choice;
			if (choose_longest_boundary_edge(p, tet, best_choice))
				edge_ids.insert(index_of(p.data.skeleton, best_choice.edge));
		}
	}


	static void run_boundary_tet_face_deletion(Context& p)
	{

		struct BoundaryTetCandidate
		{
			std::size_t tet_id = std::size_t(-1);
			SkeletonEdge edge;
			SkeletonFace f0;
			SkeletonFace f1;
		};

		std::vector<std::size_t> tet_ids;
		tet_ids.reserve(p.skeleton_tets.size());
		for (const auto& kv : p.skeleton_tets)
			tet_ids.push_back(kv.first);

		std::vector<BoundaryTetCandidate> candidates;
		candidates.reserve(tet_ids.size());

		// Phase 1: scan all current tets once and collect boundary-tet deletion candidates.
		for (std::size_t tet_id : tet_ids)
		{
			if (!p.skeleton_tets.contains(tet_id))
				continue;
			const Tet tet = p.skeleton_tets.tet(tet_id);

			if (!is_boundary_tet(p, tet))
				continue;

			BoundaryEdgeChoice best_choice;
			if (choose_longest_boundary_edge(p, tet, best_choice))
				candidates.push_back({tet_id, best_choice.edge, best_choice.f0, best_choice.f1});
		}

		// Phase 2: delete using precomputed candidates only (no boundary re-check during deletion).
		for (const BoundaryTetCandidate& cand : candidates)
		{
			if (!p.skeleton_tets.contains(cand.tet_id))
				continue;

			// the two faces of the boundary edge are faces of this tet only, and the other edges keep a face
			MeshOps::remove_faces(p, {cand.f0, cand.f1});
			p.skeleton_tets.erase(cand.tet_id);
		}
	}
};

} // namespace topology_detail

} // namespace offset_axis

} // namespace geometry

} // namespace cgogn

#endif // CGOGN_GEOMETRY_ALGOS_OFFSET_AXIS_INTERNAL_BOUNDARY_TET_PRUNE_H_
