/*******************************************************************************
 * CGoGN                                                                        *
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

#include <cgogn/core/types/mesh_traits.h>
#include <cgogn/core/functions/attributes.h>
#include <cgogn/core/functions/mesh_info.h>
#include <cgogn/geometry/types/vector_traits.h>
#include <cgogn/io/surface/surface_import.h>

#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

namespace cgogn
{
namespace io
{

template <typename MESH>
bool import_STL(MESH& m, const std::string& filename)
{
    static_assert(mesh_traits<MESH>::dimension == 2, "MESH dimension should be 2");
    
    using geometry::Vec3;
    
    SurfaceImportData surface_data;

    std::ifstream file(filename, std::ios::binary);
    if (!file.is_open())
    {
        std::cerr << "Cannot open file: " << filename << std::endl;
        return false;
    }    // Check if binary format
    char header[80];
    file.read(header, 80);
    
    uint32 num_triangles;
    file.read(reinterpret_cast<char*>(&num_triangles), sizeof(uint32));

    // Validate file size for binary STL format
    file.seekg(0, std::ios::end);
    size_t file_size = file.tellg();
    file.seekg(84, std::ios::beg); // Skip header

    bool is_binary = (file_size == 80 + 4 + num_triangles * 50);
    
    if (!is_binary)
    {
        // ASCII format
        file.close();
        file.open(filename, std::ios::in);
        if (!file.is_open())
            return false;

        std::string line;
        std::vector<Vec3> triangle_vertices;

        while (std::getline(file, line))
        {
            std::istringstream iss(line);
            std::string token;
            iss >> token;

            if (token == "vertex")
            {
                float x, y, z;
                iss >> x >> y >> z;
                triangle_vertices.push_back(Vec3(x, y, z));
            }
            else if (token == "endfacet")
            {
                if (triangle_vertices.size() >= 3)
                {
                    // Add the last 3 vertices as a triangle
                    size_t start = triangle_vertices.size() - 3;
                    for (size_t i = start; i < triangle_vertices.size(); ++i)
                    {
                        surface_data.vertex_position_.push_back(triangle_vertices[i]);
                    }
                    
                    // Add face with 3 vertices
                    surface_data.faces_nb_vertices_.push_back(3);
                    uint32 base_idx = static_cast<uint32>(surface_data.vertex_position_.size() - 3);
                    surface_data.faces_vertex_indices_.push_back(base_idx);
                    surface_data.faces_vertex_indices_.push_back(base_idx + 1);
                    surface_data.faces_vertex_indices_.push_back(base_idx + 2);
                }
            }
        }
    }
    else
    {        // Binary format
        for (uint32 i = 0; i < num_triangles; ++i)
        {
            // Skip normal vector (12 bytes)
            file.seekg(12, std::ios::cur);

            // Read 3 vertices
            for (int j = 0; j < 3; ++j)
            {
                float x, y, z;
                file.read(reinterpret_cast<char*>(&x), sizeof(float));
                file.read(reinterpret_cast<char*>(&y), sizeof(float));
                file.read(reinterpret_cast<char*>(&z), sizeof(float));

                surface_data.vertex_position_.push_back(Vec3(x, y, z));
            }

            // Skip attribute byte count (2 bytes)
            file.seekg(2, std::ios::cur);

            // Add face with 3 vertices
            surface_data.faces_nb_vertices_.push_back(3);
            uint32 base_idx = static_cast<uint32>(surface_data.vertex_position_.size() - 3);
            surface_data.faces_vertex_indices_.push_back(base_idx);
            surface_data.faces_vertex_indices_.push_back(base_idx + 1);
            surface_data.faces_vertex_indices_.push_back(base_idx + 2);
        }
    }

    surface_data.nb_vertices_ = static_cast<uint32>(surface_data.vertex_position_.size());
    surface_data.nb_faces_ = static_cast<uint32>(surface_data.faces_nb_vertices_.size());

    file.close();
    
    import_surface_data(m, surface_data);
    return true;
}

template <typename MESH>
bool export_STL(const MESH& m, const typename mesh_traits<MESH>::template Attribute<geometry::Vec3>* vertex_position,
                const std::string& filename);//TODO

} // namespace io
} // namespace cgogn

#endif // CGOGN_IO_SURFACE_STL_H_