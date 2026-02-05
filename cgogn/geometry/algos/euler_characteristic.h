/*******************************************************************************
 * CGoGN: Combinatorial and Geometric modeling with Generic N-dimensional Maps  *
 * Copyright (C), IGG Group, ICube, University of Strasbourg, France            *
 *                                                                              *
 * This library is free software; you can redistribute it and/or modify it      *
 * under the terms of the GNU Lesser General Public License as published by the *
 * Free Software Foundation; either version 2.1 of the License, or (at your     *
 * option) any later version.                                                   *
 *                                                                              *
 * This library is distributed in the hope that it will be useful, but WITHOUT  *
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or        *
 * FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License  *
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

#ifndef CGOGN_GEOMETRY_ALGOS_EULER_CHARACTERISTIC_H_
#define CGOGN_GEOMETRY_ALGOS_EULER_CHARACTERISTIC_H_

#include <cgogn/core/functions/mesh_info.h>
#include <cgogn/core/utils/tuples.h>

#include <cstdint>

namespace cgogn
{

namespace geometry
{

template <typename MESH>
inline int64 compute_euler_characteristic(const MESH& mesh)
{
	using Cells = typename mesh_traits<MESH>::Cells;
	int64 chi = 0;
	if constexpr (is_in_tuple_v<typename mesh_traits<MESH>::Vertex, Cells>)
		chi += static_cast<int64_t>(nb_cells<typename mesh_traits<MESH>::Vertex>(mesh));
	if constexpr (is_in_tuple_v<typename mesh_traits<MESH>::Edge, Cells>)
		chi -= static_cast<int64_t>(nb_cells<typename mesh_traits<MESH>::Edge>(mesh));
	if constexpr (is_in_tuple_v<typename mesh_traits<MESH>::Face, Cells>)
		chi += static_cast<int64_t>(nb_cells<typename mesh_traits<MESH>::Face>(mesh));
	if constexpr (is_in_tuple_v<typename mesh_traits<MESH>::Volume, Cells>)
		chi -= static_cast<int64_t>(nb_cells<typename mesh_traits<MESH>::Volume>(mesh));
	return chi;
}

} // namespace geometry

} // namespace cgogn

#endif // CGOGN_GEOMETRY_ALGOS_EULER_CHARACTERISTIC_H_
