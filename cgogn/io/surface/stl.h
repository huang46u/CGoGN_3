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

#ifndef CGOGN_IO_SURFACE_STL_H_
#define CGOGN_IO_SURFACE_STL_H_

#include <cgogn/io/surface/surface_import.h>
#include <cgogn/io/utils.h>

#include <cgogn/core/functions/attributes.h>
#include <cgogn/core/functions/mesh_info.h>

#include <fstream>

namespace cgogn
{

namespace io
{

template <typename MESH>
bool import_ascii_STL(MESH& m, const std::string& filename)
{
	static_assert(mesh_traits<MESH>::dimension == 2, "MESH dimension should be 2");

	using Vertex = typename mesh_traits<MESH>::Vertex;
	using Face = typename mesh_traits<MESH>::Face;

	Scoped_C_Locale loc;

	SurfaceImportData surface_data;

	std::ifstream fp(filename.c_str(), std::ios::in);

	std::string line;
	line.reserve(512u);

	// read STL header
	getline_safe(fp, line);
	if (line.find("solid") == std::string::npos)
	{
		std::cerr << "File \"" << filename << "\" is not a valid ASCII STL file." << std::endl;
		return false;
	}

	// STL files don't contain the number of vertices or faces explicitly
	// We will read the file until the end
	while (!fp.eof())
	{
		getline_safe(fp, line);
		if (line.find("facet normal") != std::string::npos)
		{
			// Read the facet normal
			std::istringstream iss(line.substr(15)); // Skip "facet normal"
			std::string normal_component;
			float64 nx, ny, nz;

			iss >> normal_component;
			nx =std::stod(normal_component);
			iss >> normal_component;
			ny = std::stod(normal_component);
			iss >> normal_component;
			nz = std::stod(normal_component);

			getline_safe(fp, line); // This should be the "outer loop" line
			std::vector<geometry::Vec3> vertices;
			for (int i = 0; i < 3; ++i) // Each face in STL is a triangle
			{
				getline_safe(fp, line);						   // This should be the "vertex" line
				std::istringstream vertex_iss(line.substr(13)); // Skip "vertex"
				std::string vertex_component;
				float64 x, y, z;

				vertex_iss >> vertex_component;
				x = std::stod(vertex_component);
				vertex_iss >> vertex_component;
				y = std::stod(vertex_component);
				vertex_iss >> vertex_component;
				z = std::stod(vertex_component);

				vertices.push_back({x, y, z});
			}

			surface_data.vertex_position_.insert(surface_data.vertex_position_.end(), vertices.begin(), vertices.end());
			// Skip the "endloop" and "endfacet" lines
			getline_safe(fp, line);
			getline_safe(fp, line);
		}
	}

	// Create faces from vertices (assuming each set of 3 vertices forms a face)
	for (uint32 i = 0; i < surface_data.vertex_position_.size(); i += 3)
	{
		std::vector<uint32> indices = {i, i + 1, i + 2};
		surface_data.faces_nb_vertices_.push_back(3);
		surface_data.faces_vertex_indices_.insert(surface_data.faces_vertex_indices_.end(), indices.begin(),
												  indices.end());
	}
	surface_data.reserve(surface_data.vertex_position_.size(), surface_data.faces_nb_vertices_.size());
	import_surface_data(m, surface_data);

	return true;
}
template <typename MESH>
bool import_binary_STL(MESH& m, const std::string& filename)
{
	static_assert(mesh_traits<MESH>::dimension == 2, "MESH dimension should be 2");
	using Vertex = typename mesh_traits<MESH>::Vertex;
	using Face = typename mesh_traits<MESH>::Face;
	Scoped_C_Locale loc;

	SurfaceImportData surface_data;
	std::ifstream file(filename, std::ios::binary);
	if (!file)
	{
		std::cerr << "Cannot open file: " << filename << std::endl;
		return false;
	}

	// Skip the header
	file.ignore(80);

	// Read the number of triangles
	uint32_t num_triangles;
	file.read(reinterpret_cast<char*>(&num_triangles), 4);
	std::vector<geometry::Vec3> vertices;
	// Read triangles
	for (uint32_t i = 0; i < num_triangles; ++i)
	{
		// Read and ignore normal vector
		file.ignore(12);

		// Read vertices
		for (int j = 0; j < 3; ++j)
		{
			float x, y, z;
			file.read(reinterpret_cast<char*>(&x), 4);
			file.read(reinterpret_cast<char*>(&y), 4);
			file.read(reinterpret_cast<char*>(&z), 4);
			
			surface_data.vertex_position_.push_back({x,y,z});
		}
		
		
		// Skip attribute byte count
		file.ignore(2);
	}

	if (!file)
	{
		std::cerr << "Error occurred while reading file: " << filename << std::endl;
		return false;
	}
	// Create faces from vertices (assuming each set of 3 vertices forms a face)
	for (uint32 i = 0; i < surface_data.vertex_position_.size(); i += 3)
	{
		std::vector<uint32> indices = {i, i + 1, i + 2};
		surface_data.faces_nb_vertices_.push_back(3);
		surface_data.faces_vertex_indices_.insert(surface_data.faces_vertex_indices_.end(), indices.begin(),
												  indices.end());
	}
	surface_data.reserve(surface_data.vertex_position_.size(), surface_data.faces_nb_vertices_.size());
	import_surface_data(m, surface_data);
	return true;
}

template <typename MESH>
bool import_STL(MESH& m, const std::string& filename)
{
	static_assert(mesh_traits<MESH>::dimension == 2, "MESH dimension should be 2");

	std::ifstream file(filename, std::ios::binary);
	if (!file)
	{
		std::cerr << "Cannot open file: " << filename << std::endl;
		return false;
	}

	// Read the first few bytes of the file to check the format
	char header[6];
	file.read(header, 5);
	header[5] = '\0';
	std::string headerStr(header);

	// Check if the file is ASCII or binary
	if (headerStr == "solid")
	{
		// Check next line to be sure
		std::string line;
		std::getline(file, line); // Finish the current line
		std::getline(file, line); // Read the next line
		if (line.find("facet normal") != std::string::npos)
		{
			// It's an ASCII STL file
			file.close();
			return import_ascii_STL(m, filename); // A separate function to handle ASCII STL
		}
	}

	// If not ASCII, it's assumed to be binary
	file.close();
	return import_binary_STL(m, filename); // A separate function to handle binary STL
}
template <typename MESH>
void export_STL(MESH& m, const typename mesh_traits<MESH>::template Attribute<geometry::Vec3>* vertex_position,
				const std::string& filename)
{
	static_assert(mesh_traits<MESH>::dimension == 2, "MESH dimension should be 2");

	using Vertex = typename mesh_traits<MESH>::Vertex;
	using Face = typename mesh_traits<MESH>::Face;

	std::ofstream out_file;
	out_file.open(filename);

	// STL header
	out_file << "solid mesh\n";

	foreach_cell(m, [&](Face f) -> bool {
		// Assuming each face is a triangle
		std::array<geometry::Vec3, 3> vertices;

		int i = 0;
		foreach_incident_vertex(m, f, [&](Vertex v) -> bool {
			if (i < 3)
			{
				vertices[i] = value<geometry::Vec3>(m, vertex_position, v);
				++i;
			}
			return true;
		});

		// Compute normal for the face
		geometry::Vec3 normal = compute_normal(vertices[0], vertices[1], vertices[2]);

		// Write face to STL
		out_file << "facet normal " << normal[0] << " " << normal[1] << " " << normal[2] << "\n";
		out_file << "    outer loop\n";
		for (const auto& vertex : vertices)
		{
			out_file << "        vertex " << vertex[0] << " " << vertex[1] << " " << vertex[2] << "\n";
		}
		out_file << "    endloop\n";
		out_file << "endfacet\n";

		return true;
	});

	out_file << "endsolid mesh\n";
	out_file.close();
}

} // namespace io

} // namespace cgogn

#endif // CGOGN_IO_SURFACE_OFF_H_
