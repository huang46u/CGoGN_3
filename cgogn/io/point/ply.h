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

#ifndef CGOGN_IO_POINT_PLY_H_
#define CGOGN_IO_POINT_PLY_H_

#include <cgogn/io/point/point_import.h>
#include <cgogn/io/utils.h>

#include <cgogn/core/functions/attributes.h>
#include <cgogn/core/functions/mesh_info.h>

#include <thirdparty/happly/happly.h>

namespace cgogn
{

namespace io
{

template <typename MESH>
typename std::enable_if<mesh_traits<MESH>::dimension == 0, bool>::type import_PLY(MESH& m, const std::string& filename)
{
	static_assert(mesh_traits<MESH>::dimension == 0, "MESH dimension should be 0");

	using Vertex = typename MESH::Vertex;

	Scoped_C_Locale loc;

	PointImportData point_data;

	happly::PLYData plyData(filename);
	std::vector<std::array<double, 3>> position = plyData.getVertexPositions();

	const uint32 nb_vertices = position.size();

	point_data.reserve(nb_vertices);

	for (uint32 i = 0u; i < nb_vertices; ++i)
	{
		const std::array<double, 3>& p = position[i];
		point_data.vertex_position_.push_back({p[0], p[1], p[2]});
	}

	import_point_data(m, point_data);

	return true;
}

template <typename MESH>
typename std::enable_if<mesh_traits<MESH>::dimension == 0, void>::type export_PLY(
	MESH& m, const typename mesh_traits<MESH>::template Attribute<geometry::Vec3>* vertex_position,
	const std::string& filename)
{
	static_assert(mesh_traits<MESH>::dimension == 0, "MESH dimension should be 0");

	using Vertex = typename mesh_traits<MESH>::Vertex;
	using Vec3 = geometry::Vec3;

	Scoped_C_Locale loc;

	std::vector<std::array<double, 3>> position;
	uint32 nb_vertices = nb_cells<Vertex>(m);
	position.reserve(nb_vertices);

	foreach_cell(m, [&](Vertex v) -> bool {
		const Vec3& p = value<geometry::Vec3>(m, vertex_position, v);
		position.push_back({p.x(), p.y(), p.z()});
		return true;
	});

	happly::PLYData plyOut;
	plyOut.addVertexPositions(position);
	plyOut.write(filename);
}

} // namespace io

} // namespace cgogn

#endif // CGOGN_IO_POINT_PLY_H_
