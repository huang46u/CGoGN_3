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
#ifndef CGOGN_GEOMETRY_ALGOS_OFFSET_AXIS_SKELETON_MESH_OPS_H_
#define CGOGN_GEOMETRY_ALGOS_OFFSET_AXIS_SKELETON_MESH_OPS_H_

#include <cgogn/core/functions/attributes.h>
#include <cgogn/core/functions/mesh_info.h>
#include <cgogn/core/functions/traversals/edge.h>
#include <cgogn/core/functions/traversals/face.h>
#include <cgogn/core/functions/traversals/vertex.h>
#include <cgogn/core/types/cell_marker.h>
#include <cgogn/geometry/algos/offset_axis/skeleton_tets.h>
#include <cgogn/geometry/algos/offset_axis/spheres_optimizer.h>

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

// Number of faces incident to each edge of the skeleton.
template <typename NONMANIFOLD>
void compute_edge_degree(const NONMANIFOLD& skeleton,
						 typename mesh_traits<NONMANIFOLD>::template Attribute<uint32>* edge_degree)
{
	using SkeletonEdge = typename mesh_traits<NONMANIFOLD>::Edge;
	parallel_foreach_cell(skeleton, [&](SkeletonEdge e) {
		(*edge_degree)[index_of(skeleton, e)] = static_cast<uint32>(incident_faces(skeleton, e).size());
		return true;
	});
}

// Removes the skeleton vertex and the sphere it comes from.
template <typename POINTS, typename NONMANIFOLD>
void remove_skeleton_vertex_and_sphere(
	NONMANIFOLD& skeleton, typename mesh_traits<NONMANIFOLD>::Vertex v, SpheresOptimizer<POINTS>& spheres_optimizer,
	typename mesh_traits<NONMANIFOLD>::template Attribute<typename mesh_traits<POINTS>::Vertex>* skeleton_source_sphere,
	typename mesh_traits<POINTS>::template Attribute<typename mesh_traits<NONMANIFOLD>::Vertex>* sphere_skeleton_vertex)
{
	using SphereVertex = typename mesh_traits<POINTS>::Vertex;
	using SkeletonVertex = typename mesh_traits<NONMANIFOLD>::Vertex;

	SphereVertex& linked_sphere = value<SphereVertex>(skeleton, skeleton_source_sphere, v);
	// invalid for a skeleton built without spheres (residual prune test)
	if (linked_sphere.is_valid())
	{
		value<SkeletonVertex>(spheres_optimizer.data().spheres, sphere_skeleton_vertex, linked_sphere) =
			SkeletonVertex();
		spheres_optimizer.remove_sphere(linked_sphere);
	}
	linked_sphere = SphereVertex();
	remove_vertex(skeleton, v);
}

// Removes the edges left without face and then their vertices left without edge (with their spheres).
template <typename POINTS, typename NONMANIFOLD>
void remove_orphan_edges_and_vertices(
	NONMANIFOLD& skeleton, const std::vector<typename mesh_traits<NONMANIFOLD>::Edge>& face_edges,
	SpheresOptimizer<POINTS>& spheres_optimizer,
	typename mesh_traits<NONMANIFOLD>::template Attribute<typename mesh_traits<POINTS>::Vertex>* skeleton_source_sphere,
	typename mesh_traits<POINTS>::template Attribute<typename mesh_traits<NONMANIFOLD>::Vertex>* sphere_skeleton_vertex)
{
	using SkeletonVertex = typename mesh_traits<NONMANIFOLD>::Vertex;
	using SkeletonEdge = typename mesh_traits<NONMANIFOLD>::Edge;

	// an edge shared by several removed faces and a vertex of several edges are visited once: removing a cell twice
	// releases its index twice
	CellMarkerStore<NONMANIFOLD, SkeletonEdge> visited_edges(skeleton);
	CellMarkerStore<NONMANIFOLD, SkeletonVertex> visited_vertices(skeleton);
	std::vector<SkeletonVertex> candidate_vertices;
	candidate_vertices.reserve(face_edges.size() * 2);
	for (SkeletonEdge e : face_edges)
	{
		if (visited_edges.is_marked(e))
			continue;
		visited_edges.mark(e);
		for (SkeletonVertex v : incident_vertices(skeleton, e))
		{
			if (visited_vertices.is_marked(v))
				continue;
			visited_vertices.mark(v);
			candidate_vertices.push_back(v);
		}
		if (incident_faces(skeleton, e).empty())
			remove_edge(skeleton, e);
	}
	for (SkeletonVertex v : candidate_vertices)
		if (incident_edges(skeleton, v).empty())
			remove_skeleton_vertex_and_sphere(skeleton, v, spheres_optimizer, skeleton_source_sphere,
											  sphere_skeleton_vertex);
}

// Removes the faces, then the edges and vertices left isolated (orphan edges and vertices). The tets are not
// changed.
template <typename POINTS, typename NONMANIFOLD>
void remove_faces(
	NONMANIFOLD& skeleton, const std::vector<typename mesh_traits<NONMANIFOLD>::Face>& faces,
	SpheresOptimizer<POINTS>& spheres_optimizer,
	typename mesh_traits<NONMANIFOLD>::template Attribute<typename mesh_traits<POINTS>::Vertex>* skeleton_source_sphere,
	typename mesh_traits<POINTS>::template Attribute<typename mesh_traits<NONMANIFOLD>::Vertex>* sphere_skeleton_vertex)
{
	using SkeletonEdge = typename mesh_traits<NONMANIFOLD>::Edge;
	using SkeletonFace = typename mesh_traits<NONMANIFOLD>::Face;

	std::vector<SkeletonEdge> face_edges;
	face_edges.reserve(faces.size() * 3);
	for (SkeletonFace f : faces)
		for (SkeletonEdge e : incident_edges(skeleton, f))
			face_edges.push_back(e);
	for (SkeletonFace f : faces)
		remove_face(skeleton, f);
	remove_orphan_edges_and_vertices(skeleton, face_edges, spheres_optimizer, skeleton_source_sphere,
									 sphere_skeleton_vertex);
}

// Removes every vertex without edge (with its sphere).
template <typename POINTS, typename NONMANIFOLD>
void remove_isolated_vertices(
	NONMANIFOLD& skeleton, SpheresOptimizer<POINTS>& spheres_optimizer,
	typename mesh_traits<NONMANIFOLD>::template Attribute<typename mesh_traits<POINTS>::Vertex>* skeleton_source_sphere,
	typename mesh_traits<POINTS>::template Attribute<typename mesh_traits<NONMANIFOLD>::Vertex>* sphere_skeleton_vertex)
{
	using SkeletonVertex = typename mesh_traits<NONMANIFOLD>::Vertex;

	std::vector<SkeletonVertex> isolated_vertices;
	foreach_cell(skeleton, [&](SkeletonVertex v) -> bool {
		if (incident_edges(skeleton, v).empty())
			isolated_vertices.push_back(v);
		return true;
	});
	for (SkeletonVertex v : isolated_vertices)
		remove_skeleton_vertex_and_sphere(skeleton, v, spheres_optimizer, skeleton_source_sphere,
										  sphere_skeleton_vertex);
}

// Removes the faces one by one with their tets, then every vertex without edge (with its sphere).
template <typename POINTS, typename NONMANIFOLD>
void remove_faces_and_tets(
	NONMANIFOLD& skeleton, SkeletonTets<NONMANIFOLD>& tets, const std::unordered_set<uint32>& face_ids,
	SpheresOptimizer<POINTS>& spheres_optimizer,
	typename mesh_traits<NONMANIFOLD>::template Attribute<typename mesh_traits<POINTS>::Vertex>* skeleton_source_sphere,
	typename mesh_traits<POINTS>::template Attribute<typename mesh_traits<NONMANIFOLD>::Vertex>* sphere_skeleton_vertex)
{
	using SkeletonFace = typename mesh_traits<NONMANIFOLD>::Face;

	for (uint32 face_id : face_ids)
	{
		const std::set<std::size_t> face_tets = tets.tets_of(face_id);
		remove_faces(skeleton, {of_index<SkeletonFace>(skeleton, face_id)}, spheres_optimizer, skeleton_source_sphere,
					 sphere_skeleton_vertex);
		for (std::size_t tet_id : face_tets)
			if (tets.contains(tet_id))
				tets.erase(tet_id);
	}
	remove_isolated_vertices(skeleton, spheres_optimizer, skeleton_source_sphere, sphere_skeleton_vertex);
}

// The edge between a and b, or an invalid edge if there is none.
template <typename NONMANIFOLD>
typename mesh_traits<NONMANIFOLD>::Edge find_edge(const NONMANIFOLD& skeleton,
												  typename mesh_traits<NONMANIFOLD>::Vertex a,
												  typename mesh_traits<NONMANIFOLD>::Vertex b)
{
	using SkeletonVertex = typename mesh_traits<NONMANIFOLD>::Vertex;
	using SkeletonEdge = typename mesh_traits<NONMANIFOLD>::Edge;

	for (SkeletonEdge e : incident_edges(skeleton, a))
	{
		const std::vector<SkeletonVertex> edge_vertices = incident_vertices(skeleton, e);
		if (edge_vertices[0] == b || edge_vertices[1] == b)
			return e;
	}
	return SkeletonEdge();
}

} // namespace offset_axis

} // namespace geometry

} // namespace cgogn

#endif // CGOGN_GEOMETRY_ALGOS_OFFSET_AXIS_SKELETON_MESH_OPS_H_
