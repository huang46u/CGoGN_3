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
#ifndef CGOGN_GEOMETRY_ALGOS_OFFSET_AXIS_BOUNDARY_TET_PRUNE_H_
#define CGOGN_GEOMETRY_ALGOS_OFFSET_AXIS_BOUNDARY_TET_PRUNE_H_

#include <cgogn/core/functions/attributes.h>
#include <cgogn/core/functions/mesh_info.h>
#include <cgogn/core/functions/traversals/edge.h>
#include <cgogn/core/functions/traversals/face.h>
#include <cgogn/core/functions/traversals/vertex.h>
#include <cgogn/core/types/cell_marker.h>
#include <cgogn/geometry/algos/offset_axis/skeleton_mesh_ops.h>

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

template <typename NONMANIFOLD>
struct BoundaryEdgeChoice
{
	typename mesh_traits<NONMANIFOLD>::Edge edge;
	typename mesh_traits<NONMANIFOLD>::Face f0;
	typename mesh_traits<NONMANIFOLD>::Face f1;
	Scalar length = Scalar(-1);
};

// A tet is a boundary tet when one of its faces has at least two degree-2 edges.
template <typename NONMANIFOLD>
bool is_boundary_tet(const NONMANIFOLD& skeleton, const typename SkeletonTets<NONMANIFOLD>::Tet& tet)
{
	using SkeletonEdge = typename mesh_traits<NONMANIFOLD>::Edge;
	using SkeletonFace = typename mesh_traits<NONMANIFOLD>::Face;

	for (SkeletonFace f : tet.faces)
	{
		uint32 deg2_edge_count = 0;
		for (SkeletonEdge e : incident_edges(skeleton, f))
		{
			if (incident_faces(skeleton, e).size() == 2)
				++deg2_edge_count;
		}
		if (deg2_edge_count >= 2)
			return true;
	}
	return false;
}

// Longest degree-2 edge among the faces of the tet. Returns false when there is none.
template <typename NONMANIFOLD>
bool choose_longest_boundary_edge(const NONMANIFOLD& skeleton,
								  const typename mesh_traits<NONMANIFOLD>::template Attribute<Vec3>* skeleton_position,
								  const typename SkeletonTets<NONMANIFOLD>::Tet& tet,
								  BoundaryEdgeChoice<NONMANIFOLD>& best_choice)
{
	using SkeletonVertex = typename mesh_traits<NONMANIFOLD>::Vertex;
	using SkeletonEdge = typename mesh_traits<NONMANIFOLD>::Edge;
	using SkeletonFace = typename mesh_traits<NONMANIFOLD>::Face;

	CellMarkerStore<NONMANIFOLD, SkeletonEdge> seen_edges(skeleton);

	for (SkeletonFace f : tet.faces)
	{
		for (SkeletonEdge e : incident_edges(skeleton, f))
		{
			if (seen_edges.is_marked(e))
				continue;
			seen_edges.mark(e);

			const std::vector<SkeletonFace> in_faces = incident_faces(skeleton, e);
			if (in_faces.size() != 2)
				continue;
			const std::vector<SkeletonVertex> edge_vertices = incident_vertices(skeleton, e);
			const Vec3 p0 = value<Vec3>(skeleton, skeleton_position, edge_vertices[0]);
			const Vec3 p1 = value<Vec3>(skeleton, skeleton_position, edge_vertices[1]);
			const Scalar edge_len = (p1 - p0).norm();
			if (edge_len > best_choice.length)
			{
				best_choice.edge = e;
				best_choice.f0 = in_faces[0];
				best_choice.f1 = in_faces[1];
				best_choice.length = edge_len;
			}
		}
	}

	return best_choice.length >= Scalar(0);
}

// Read-only: collects the faces of every current boundary tet and the longest degree-2 edge of each such tet.
template <typename NONMANIFOLD>
void boundary_tet_mask(const NONMANIFOLD& skeleton, const SkeletonTets<NONMANIFOLD>& tets,
					   const typename mesh_traits<NONMANIFOLD>::template Attribute<Vec3>* skeleton_position,
					   std::unordered_set<uint32>& face_ids, std::unordered_set<uint32>& edge_ids)
{
	for (const auto& [tet_id, tet] : tets)
	{
		if (!is_boundary_tet(skeleton, tet))
			continue;
		for (uint32 i = 0; i < 4; ++i)
			face_ids.insert(index_of(skeleton, tet.faces[i]));
		BoundaryEdgeChoice<NONMANIFOLD> best_choice;
		if (choose_longest_boundary_edge(skeleton, skeleton_position, tet, best_choice))
			edge_ids.insert(index_of(skeleton, best_choice.edge));
	}
}

// Removes the two faces of the longest degree-2 edge of each boundary tet, and the tet. The boundary tets are found
// before any removal.
template <typename POINTS, typename NONMANIFOLD>
void remove_boundary_tets(
	NONMANIFOLD& skeleton, SkeletonTets<NONMANIFOLD>& tets,
	const typename mesh_traits<NONMANIFOLD>::template Attribute<Vec3>* skeleton_position,
	SpheresOptimizer<POINTS>& spheres_optimizer,
	typename mesh_traits<NONMANIFOLD>::template Attribute<typename mesh_traits<POINTS>::Vertex>* skeleton_source_sphere,
	typename mesh_traits<POINTS>::template Attribute<typename mesh_traits<NONMANIFOLD>::Vertex>* sphere_skeleton_vertex)
{
	using SkeletonFace = typename mesh_traits<NONMANIFOLD>::Face;

	struct BoundaryTetCandidate
	{
		std::size_t tet_id;
		SkeletonFace f0;
		SkeletonFace f1;
	};

	std::vector<BoundaryTetCandidate> candidates;
	candidates.reserve(tets.size());
	for (const auto& [tet_id, tet] : tets)
	{
		if (!is_boundary_tet(skeleton, tet))
			continue;
		BoundaryEdgeChoice<NONMANIFOLD> best_choice;
		if (choose_longest_boundary_edge(skeleton, skeleton_position, tet, best_choice))
			candidates.push_back({tet_id, best_choice.f0, best_choice.f1});
	}

	for (const BoundaryTetCandidate& cand : candidates)
	{
		// the two faces of the boundary edge are faces of this tet only, and the other edges keep a face
		remove_faces(skeleton, {cand.f0, cand.f1}, spheres_optimizer, skeleton_source_sphere, sphere_skeleton_vertex);
		tets.erase(cand.tet_id);
	}
}

} // namespace offset_axis

} // namespace geometry

} // namespace cgogn

#endif // CGOGN_GEOMETRY_ALGOS_OFFSET_AXIS_BOUNDARY_TET_PRUNE_H_
