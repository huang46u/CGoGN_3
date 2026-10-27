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
	using Tet = typename Types::Tet;

	static bool remove_skeleton_vertex_and_linked_sphere(Context& p, const SkeletonVertex& v)
	{
		if (!v.is_valid())
			return false;
		const uint32 skeleton_vertex_id = index_of(p.data.skeleton, v);
		if (skeleton_vertex_id == INVALID_INDEX)
			return false;

		SphereVertex linked_sphere;
		if (p.data.skeleton_source_sphere)
			linked_sphere = value<SphereVertex>(p.data.skeleton, p.data.skeleton_source_sphere, v);

		if (linked_sphere.is_valid())
		{
			const uint32 sphere_index = index_of(p.optimizer_data.spheres, linked_sphere);
			if (sphere_index != INVALID_INDEX && p.data.sphere_skeleton_vertex)
				(*p.data.sphere_skeleton_vertex)[sphere_index] = SkeletonVertex();
			p.data.spheres_optimizer.remove_sphere(linked_sphere);
		}

		if (p.data.skeleton_source_sphere)
			value<SphereVertex>(p.data.skeleton, p.data.skeleton_source_sphere, v) = SphereVertex();
		remove_vertex(p.data.skeleton, v);
		return true;
	}
	static void compute_edge_degree(Context& p)
	{
		offset_axis::compute_edge_degree(p.data.skeleton, p.data.edge_degree);
	}
	// Removes the edges left without face and then their vertices left without edge (with their spheres).
	static void remove_orphan_edges_and_vertices(Context& p, const std::vector<SkeletonEdge>& face_edges)
	{
		std::vector<SkeletonVertex> candidate_vertices;
		candidate_vertices.reserve(face_edges.size() * 2);
		for (SkeletonEdge e : face_edges)
		{
			if (!e.is_valid())
				continue;
			const std::vector<SkeletonVertex> edge_vertices = incident_vertices(p.data.skeleton, e);
			for (const SkeletonVertex& v : edge_vertices)
				candidate_vertices.push_back(v);
			if (incident_faces(p.data.skeleton, e).empty())
				remove_edge(p.data.skeleton, e);
		}
		for (const SkeletonVertex& v : candidate_vertices)
		{
			if (!v.is_valid())
				continue;
			const uint32 idv = index_of(p.data.skeleton, v);
			if (idv == INVALID_INDEX)
				continue;
			if (!incident_edges(p.data.skeleton, v).empty())
				continue;
			remove_skeleton_vertex_and_linked_sphere(p, v);
		}
	}

	// Removes the faces, then the edges and vertices left isolated (orphan edges and vertices). The tets are not
	// changed. Returns the edges of the removed faces.
	static std::vector<SkeletonEdge> remove_faces(Context& p, const std::vector<SkeletonFace>& faces)
	{
		std::vector<SkeletonEdge> face_edges;
		face_edges.reserve(faces.size() * 3);
		for (SkeletonFace f : faces)
			for (SkeletonEdge e : incident_edges(p.data.skeleton, f))
				face_edges.push_back(e);
		for (SkeletonFace f : faces)
			remove_face(p.data.skeleton, f);
		remove_orphan_edges_and_vertices(p, face_edges);
		return face_edges;
	}

	// Removes the faces one by one with their tets, then every vertex without edge (with its sphere).
	static void remove_faces_and_tets(Context& p, const std::unordered_set<uint32>& face_ids)
	{
		for (uint32 face_id : face_ids)
		{
			const std::set<std::size_t> face_tets = p.skeleton_tets.tets_of(face_id);
			remove_faces(p, {of_index<SkeletonFace>(p.data.skeleton, face_id)});
			for (std::size_t tet_id : face_tets)
				if (p.skeleton_tets.contains(tet_id))
					p.skeleton_tets.erase(tet_id);
		}
		remove_isolated_vertices(p);
	}

	static void remove_isolated_vertices(Context& p)
	{
		std::vector<SkeletonVertex> isolated_vertices;
		foreach_cell(p.data.skeleton, [&](SkeletonVertex v) -> bool {
			if (incident_edges(p.data.skeleton, v).empty())
				isolated_vertices.push_back(v);
			return true;
		});
		for (SkeletonVertex v : isolated_vertices)
			remove_skeleton_vertex_and_linked_sphere(p, v);
	}

	static bool find_skeleton_edge_between(Context& p, const SkeletonVertex& a, const SkeletonVertex& b,
										   SkeletonEdge& out_edge)
	{
		out_edge = SkeletonEdge();
		if (!a.is_valid() || !b.is_valid() || a == b)
			return false;
		for (SkeletonEdge edge : incident_edges(p.data.skeleton, a))
		{
			if (!edge.is_valid())
				continue;
			const std::vector<SkeletonVertex> edge_vertices = incident_vertices(p.data.skeleton, edge);
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
