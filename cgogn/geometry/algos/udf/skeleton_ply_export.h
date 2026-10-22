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
#ifndef CGOGN_GEOMETRY_ALGOS_UDF_SKELETON_PLY_EXPORT_H_
#define CGOGN_GEOMETRY_ALGOS_UDF_SKELETON_PLY_EXPORT_H_

#include <cgogn/core/functions/attributes.h>
#include <cgogn/geometry/types/vector_traits.h>
#include <cgogn/io/surface/export_options.h>
#include <cgogn/io/surface/ply.h>

#include <filesystem>
#include <memory>
#include <system_error>

namespace cgogn::geometry
{

// Writes the skeleton with its per-vertex radius to a PLY file (directories are created, an existing file is replaced).
// Returns false if the file could not be written or is empty.
template <typename NONMANIFOLD>
bool export_skeleton_ply(NONMANIFOLD& skeleton,
						 const std::shared_ptr<typename mesh_traits<NONMANIFOLD>::template Attribute<Vec3>>& position,
						 const std::shared_ptr<typename mesh_traits<NONMANIFOLD>::template Attribute<Scalar>>& radius,
						 const std::filesystem::path& path)
{
	if (!position || !radius)
		return false;
	std::error_code error;
	if (!path.parent_path().empty())
	{
		std::filesystem::create_directories(path.parent_path(), error);
		if (error)
			return false;
	}
	if (std::filesystem::exists(path, error))
	{
		std::filesystem::remove(path, error);
		if (error)
			return false;
	}
	io::SurfaceExportAttributeSelection<NONMANIFOLD> attributes;
	attributes.vertex_attributes.push_back(radius);
	io::export_PLY(skeleton, position.get(), path.string(), &attributes);
	const bool exists = std::filesystem::exists(path, error);
	if (error || !exists)
		return false;
	const auto size = std::filesystem::file_size(path, error);
	return !error && size > 0;
}

} // namespace cgogn::geometry

#endif // CGOGN_GEOMETRY_ALGOS_UDF_SKELETON_PLY_EXPORT_H_
