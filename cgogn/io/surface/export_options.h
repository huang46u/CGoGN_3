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

#ifndef CGOGN_IO_SURFACE_EXPORT_OPTIONS_H_
#define CGOGN_IO_SURFACE_EXPORT_OPTIONS_H_

#include <cgogn/core/types/mesh_traits.h>

#include <memory>
#include <vector>

namespace cgogn
{

namespace io
{

template <typename MESH>
struct SurfaceExportAttributeSelection
{
	using AttributeGen = typename mesh_traits<MESH>::AttributeGen;
	std::vector<std::shared_ptr<AttributeGen>> vertex_attributes;
	std::vector<std::shared_ptr<AttributeGen>> face_attributes;
};

} // namespace io

} // namespace cgogn

#endif // CGOGN_IO_SURFACE_EXPORT_OPTIONS_H_
