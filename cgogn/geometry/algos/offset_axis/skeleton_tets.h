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
#ifndef CGOGN_GEOMETRY_ALGOS_OFFSET_AXIS_SKELETON_TETS_H_
#define CGOGN_GEOMETRY_ALGOS_OFFSET_AXIS_SKELETON_TETS_H_

#include <cgogn/core/functions/attributes.h>
#include <cgogn/core/types/mesh_traits.h>
#include <cgogn/core/utils/numerics.h>

#include <cstddef>
#include <set>
#include <unordered_map>

namespace cgogn
{

namespace geometry
{

namespace offset_axis
{

// The tets of the skeleton: 4 faces each, not a cell of the mesh. The tets of each face are kept in a face attribute,
// written only here.
template <typename NONMANIFOLD>
class SkeletonTets
{
public:
	using Face = typename mesh_traits<NONMANIFOLD>::Face;
	template <typename T>
	using Attribute = typename mesh_traits<NONMANIFOLD>::template Attribute<T>;

	struct Tet
	{
		Face faces[4];
	};
	using Map = std::unordered_map<std::size_t, Tet>;

	SkeletonTets(const NONMANIFOLD& skeleton, Attribute<std::set<std::size_t>>* face_tets)
		: skeleton_(skeleton), face_tets_(face_tets)
	{
	}

	void clear()
	{
		tets_.clear();
	}

	void reserve(std::size_t n)
	{
		tets_.reserve(n);
	}

	// A new face may reuse the attribute slot of a removed face.
	void reset_face(Face f)
	{
		value<std::set<std::size_t>>(skeleton_, face_tets_, f).clear();
	}

	void add(std::size_t tet_id, const Tet& tet)
	{
		for (uint32 i = 0; i < 4; ++i)
			value<std::set<std::size_t>>(skeleton_, face_tets_, tet.faces[i]).insert(tet_id);
		tets_.insert({tet_id, tet});
	}

	void erase(std::size_t tet_id)
	{
		const auto it = tets_.find(tet_id);
		for (uint32 i = 0; i < 4; ++i)
			value<std::set<std::size_t>>(skeleton_, face_tets_, it->second.faces[i]).erase(tet_id);
		tets_.erase(it);
	}

	bool contains(std::size_t tet_id) const
	{
		return tets_.find(tet_id) != tets_.end();
	}

	const Tet& tet(std::size_t tet_id) const
	{
		return tets_.at(tet_id);
	}

	// Tets of the face of index face_index.
	const std::set<std::size_t>& tets_of(uint32 face_index) const
	{
		return (*face_tets_)[face_index];
	}

	std::size_t size() const
	{
		return tets_.size();
	}

	bool empty() const
	{
		return tets_.empty();
	}

	typename Map::const_iterator begin() const
	{
		return tets_.begin();
	}

	typename Map::const_iterator end() const
	{
		return tets_.end();
	}

private:
	const NONMANIFOLD& skeleton_;
	Attribute<std::set<std::size_t>>* face_tets_;
	Map tets_;
};

} // namespace offset_axis

} // namespace geometry

} // namespace cgogn

#endif // CGOGN_GEOMETRY_ALGOS_OFFSET_AXIS_SKELETON_TETS_H_
