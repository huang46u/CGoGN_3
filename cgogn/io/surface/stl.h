#ifndef CGOGN_IO_SURFACE_STL_H_
#define CGOGN_IO_SURFACE_STL_H_

#include <cgogn/io/surface/surface_import.h>
#include <cgogn/io/utils.h>

#include <algorithm>
#include <array>
#include <cctype>
#include <cstdint>
#include <fstream>
#include <string>

namespace cgogn
{

namespace io
{

template <typename MESH>
bool import_STL(MESH& m, const std::string& filename)
{
	static_assert(mesh_traits<MESH>::dimension == 2, "MESH dimension should be 2");

	Scoped_C_Locale loc;

	std::ifstream fp(filename, std::ios::binary);
	if (!fp.good())
	{
		std::cerr << "Error opening file " << filename << std::endl;
		return false;
	}

	fp.seekg(0, std::ios::end);
	uint64 file_size = static_cast<uint64>(fp.tellg());
	fp.seekg(0, std::ios::beg);

	if (file_size < 84u)
	{
		std::cerr << "File \"" << filename << "\" is too small to be valid STL." << std::endl;
		return false;
	}

	char header[80];
	fp.read(header, 80);
	uint32 tri_count = 0;
	fp.read(reinterpret_cast<char*>(&tri_count), sizeof(uint32));

	uint64 expected_binary_size = 84u + 50u * static_cast<uint64>(tri_count);
	bool is_binary = (expected_binary_size == file_size);

	SurfaceImportData surface_data;

	if (is_binary)
	{
		surface_data.reserve(tri_count * 3u, tri_count);
		surface_data.nb_vertices_ = 0u;
		surface_data.nb_faces_ = 0u;
		for (uint32 i = 0; i < tri_count; ++i)
		{
			float normal[3];
			float v[9];
			uint16 attr = 0;

			fp.read(reinterpret_cast<char*>(normal), sizeof(normal));
			fp.read(reinterpret_cast<char*>(v), sizeof(v));
			fp.read(reinterpret_cast<char*>(&attr), sizeof(attr));

			uint32 base = static_cast<uint32>(surface_data.vertex_position_.size());
			surface_data.vertex_position_.push_back(Vec3(v[0], v[1], v[2]));
			surface_data.vertex_position_.push_back(Vec3(v[3], v[4], v[5]));
			surface_data.vertex_position_.push_back(Vec3(v[6], v[7], v[8]));

			surface_data.faces_nb_vertices_.push_back(3);
			surface_data.faces_vertex_indices_.push_back(base);
			surface_data.faces_vertex_indices_.push_back(base + 1);
			surface_data.faces_vertex_indices_.push_back(base + 2);

			surface_data.nb_vertices_ += 3u;
			surface_data.nb_faces_ += 1u;
		}
	}
	else
	{
		fp.clear();
		fp.seekg(0, std::ios::beg);

		std::string token;
		std::array<Vec3, 3> tri;
		int tri_vertex_count = 0;

		while (fp >> token)
		{
			std::transform(token.begin(), token.end(), token.begin(),
						   [](unsigned char c) { return static_cast<char>(std::tolower(c)); });
			if (token == "vertex")
			{
				float64 x = 0.0;
				float64 y = 0.0;
				float64 z = 0.0;
				fp >> x >> y >> z;
				tri[tri_vertex_count] = Vec3(x, y, z);
				tri_vertex_count++;
				if (tri_vertex_count == 3)
				{
					uint32 base = static_cast<uint32>(surface_data.vertex_position_.size());
					surface_data.vertex_position_.push_back(tri[0]);
					surface_data.vertex_position_.push_back(tri[1]);
					surface_data.vertex_position_.push_back(tri[2]);

					surface_data.faces_nb_vertices_.push_back(3);
					surface_data.faces_vertex_indices_.push_back(base);
					surface_data.faces_vertex_indices_.push_back(base + 1);
					surface_data.faces_vertex_indices_.push_back(base + 2);

					surface_data.nb_vertices_ += 3u;
					surface_data.nb_faces_ += 1u;
					tri_vertex_count = 0;
				}
			}
		}
	}

	if (surface_data.nb_faces_ == 0u)
	{
		std::cerr << "File \"" << filename << " has no faces." << std::endl;
		return false;
	}

	import_surface_data(m, surface_data);
	return true;
}

} // namespace io

} // namespace cgogn

#endif // CGOGN_IO_SURFACE_STL_H_
