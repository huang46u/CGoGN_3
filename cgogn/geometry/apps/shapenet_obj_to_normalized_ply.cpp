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

#include <cgogn/core/functions/attributes.h>
#include <cgogn/core/functions/traversals/vertex.h>
#include <cgogn/core/types/cell_marker.h>
#include <cgogn/core/types/maps/cmap/cmap2.h>

#include <cgogn/geometry/functions/bounding_box.h>
#include <cgogn/geometry/types/vector_traits.h>

#include <cgogn/io/surface/obj.h>
#include <cgogn/io/surface/ply.h>

#include <algorithm>
#include <filesystem>
#include <iostream>
#include <string>
#include <vector>

using Surface = cgogn::CMap2;
using Vertex = cgogn::mesh_traits<Surface>::Vertex;
using cgogn::geometry::Vec3;

namespace
{

struct BatchConfig
{
	std::filesystem::path source_ply_dir =
		"D:/Code/CGoGN_3/data/skeleton/udf/shapenetcar_200/model";
	std::filesystem::path shapenet_category_dir =
		"D:/Code/Dataset/ShapeNet/02958343";
	std::filesystem::path output_ply_dir =
		"D:/Code/CGoGN_3/data/skeleton/udf/shapenetcar_200/model_ply";
};

std::vector<std::filesystem::path> collect_ply_files(const std::filesystem::path& directory)
{
	std::vector<std::filesystem::path> files;
	for (const auto& entry : std::filesystem::directory_iterator(directory))
	{
		if (!entry.is_regular_file())
			continue;
		if (entry.path().extension() == ".ply")
			files.push_back(entry.path());
	}
	std::sort(files.begin(), files.end());
	return files;
}

std::filesystem::path shapenet_obj_path_for_id(const std::filesystem::path& category_dir, const std::string& model_id)
{
	return category_dir / model_id / "models" / "model_normalized.obj";
}

bool normalize_obj_and_export_ply(
	const std::filesystem::path& input_obj,
	const std::filesystem::path& output_ply,
	Vec3& out_center,
	cgogn::geometry::Scalar& out_scale)
{
	Surface mesh;
	if (!cgogn::io::import_OBJ(mesh, input_obj.string()))
		return false;

	auto position = cgogn::get_attribute<Vec3, Vertex>(mesh, "position");
	if (!position)
		return false;

	std::tie(out_center, out_scale) = cgogn::geometry::normalize_centered(*position);
	cgogn::io::export_PLY(mesh, position.get(), output_ply.string());
	return true;
}

void print_usage(const char* exe_name)
{
	std::cout << "Usage: " << exe_name
			  << " [source_ply_dir] [shapenet_category_dir] [output_ply_dir]\n"
			  << "Defaults:\n"
			  << "  source_ply_dir        D:/Code/CGoGN_3/data/skeleton/udf/shapenetcar_200/model\n"
			  << "  shapenet_category_dir D:/Code/Dataset/ShapeNet/02958343\n"
			  << "  output_ply_dir        D:/Code/CGoGN_3/data/skeleton/udf/shapenetcar_200/model_ply\n"
			  << "Normalization matches geometry::normalize_centered used by udf_training.h:\n"
			  << "  v' = (v - bbox_center) / max_bbox_extent\n";
}

} // namespace

int main(int argc, char** argv)
{
	BatchConfig config;

	if (argc > 1)
	{
		const std::string first_arg = argv[1];
		if (first_arg == "-h" || first_arg == "--help")
		{
			print_usage(argv[0]);
			return 0;
		}
		config.source_ply_dir = argv[1];
	}
	if (argc > 2)
		config.shapenet_category_dir = argv[2];
	if (argc > 3)
		config.output_ply_dir = argv[3];
	if (argc > 4)
	{
		print_usage(argv[0]);
		return 1;
	}

	if (!std::filesystem::exists(config.source_ply_dir))
	{
		std::cerr << "Source PLY directory does not exist: " << config.source_ply_dir << std::endl;
		return 1;
	}
	if (!std::filesystem::exists(config.shapenet_category_dir))
	{
		std::cerr << "ShapeNet category directory does not exist: " << config.shapenet_category_dir << std::endl;
		return 1;
	}

	std::error_code ec;
	std::filesystem::create_directories(config.output_ply_dir, ec);
	if (ec)
	{
		std::cerr << "Failed to create output directory: " << config.output_ply_dir << "\n"
				  << "Reason: " << ec.message() << std::endl;
		return 1;
	}

	const std::vector<std::filesystem::path> ply_files = collect_ply_files(config.source_ply_dir);
	if (ply_files.empty())
	{
		std::cerr << "No PLY files found in: " << config.source_ply_dir << std::endl;
		return 1;
	}

	std::size_t converted_count = 0;
	std::size_t missing_obj_count = 0;
	std::size_t failed_count = 0;

	for (std::size_t i = 0; i < ply_files.size(); ++i)
	{
		const std::filesystem::path& ply_path = ply_files[i];
		const std::string model_id = ply_path.stem().string();
		const std::filesystem::path input_obj = shapenet_obj_path_for_id(config.shapenet_category_dir, model_id);
		const std::filesystem::path output_ply = config.output_ply_dir / ply_path.filename();

		std::cout << "[" << (i + 1) << "/" << ply_files.size() << "] " << model_id << std::endl;
		if (!std::filesystem::exists(input_obj))
		{
			++missing_obj_count;
			std::cerr << "  Missing OBJ: " << input_obj << std::endl;
			continue;
		}

		Vec3 center = Vec3::Zero();
		cgogn::geometry::Scalar scale = 0;
		if (!normalize_obj_and_export_ply(input_obj, output_ply, center, scale))
		{
			++failed_count;
			std::cerr << "  Failed to convert OBJ -> normalized PLY: " << input_obj << std::endl;
			continue;
		}

		++converted_count;
		std::cout << "  Wrote: " << output_ply << " | center=(" << center.transpose() << ") scale=" << scale
				  << std::endl;
	}

	std::cout << "\nSummary\n"
			  << "  Total source ids: " << ply_files.size() << "\n"
			  << "  Converted: " << converted_count << "\n"
			  << "  Missing OBJ: " << missing_obj_count << "\n"
			  << "  Failed: " << failed_count << std::endl;

	return (missing_obj_count == 0 && failed_count == 0) ? 0 : 1;
}
