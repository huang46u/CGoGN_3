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
#ifndef CGOGN_GEOMETRY_ALGOS_OFFSET_AXIS_SKELETON_FACE_COMPONENTS_H_
#define CGOGN_GEOMETRY_ALGOS_OFFSET_AXIS_SKELETON_FACE_COMPONENTS_H_

#include <cgogn/core/functions/attributes.h>
#include <cgogn/core/functions/mesh_info.h>
#include <cgogn/core/functions/traversals/edge.h>
#include <cgogn/core/functions/traversals/face.h>
#include <cgogn/geometry/algos/offset_axis/internal/skeleton_topology_types.h>

#include <cstddef>
#include <numeric>
#include <unordered_map>
#include <vector>

namespace cgogn
{

namespace geometry
{

namespace offset_axis
{

// Read-only: splits the faces of a skeleton into sheets. Two faces belong to the same sheet when they are
// connected through a manifold edge, i.e. an edge incident to exactly two faces.
// Returns a map from face index (index_of) to sheet id. Sheet ids are numbered from 0 in face traversal order.
// The skeleton is not modified.
template <typename NONMANIFOLD>
std::unordered_map<uint32, uint32> compute_skeleton_face_components(NONMANIFOLD& skeleton)
{
	using Face = typename mesh_traits<NONMANIFOLD>::Face;
	using Edge = typename mesh_traits<NONMANIFOLD>::Edge;

	std::vector<Face> faces;
	faces.reserve(nb_cells<Face>(skeleton));
	std::unordered_map<uint32, uint32> face_slot;
	face_slot.reserve(nb_cells<Face>(skeleton));
	foreach_cell(skeleton, [&](Face f) -> bool {
		const uint32 idf = index_of(skeleton, f);
		if (idf == INVALID_INDEX)
			return true;
		face_slot[idf] = static_cast<uint32>(faces.size());
		faces.push_back(f);
		return true;
	});

	std::vector<uint32> parent(faces.size());
	std::iota(parent.begin(), parent.end(), 0u);
	auto find_root = [&parent](uint32 x) {
		while (parent[x] != x)
		{
			parent[x] = parent[parent[x]];
			x = parent[x];
		}
		return x;
	};

	foreach_cell(skeleton, [&](Edge e) -> bool {
		const auto in_faces = incident_faces(skeleton, e);
		if (in_faces.size() != 2)
			return true;
		const auto it0 = face_slot.find(index_of(skeleton, in_faces[0]));
		const auto it1 = face_slot.find(index_of(skeleton, in_faces[1]));
		if (it0 == face_slot.end() || it1 == face_slot.end())
			return true;
		const uint32 r0 = find_root(it0->second);
		const uint32 r1 = find_root(it1->second);
		if (r0 != r1)
			parent[r1] = r0;
		return true;
	});

	std::unordered_map<uint32, uint32> sheet_of_face;
	sheet_of_face.reserve(faces.size());
	std::unordered_map<uint32, uint32> root_to_sheet;
	for (uint32 i = 0; i < static_cast<uint32>(faces.size()); ++i)
	{
		const uint32 root = find_root(i);
		const auto [it, inserted] = root_to_sheet.emplace(root, static_cast<uint32>(root_to_sheet.size()));
		(void)inserted;
		sheet_of_face[index_of(skeleton, faces[i])] = it->second;
	}
	return sheet_of_face;
}

} // namespace offset_axis

} // namespace geometry
} // namespace cgogn

#endif // CGOGN_GEOMETRY_ALGOS_OFFSET_AXIS_SKELETON_FACE_COMPONENTS_H_
