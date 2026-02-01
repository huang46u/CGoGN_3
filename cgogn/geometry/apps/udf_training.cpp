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

#include <cgogn/core/types/incidence_graph/incidence_graph.h>
#include <cgogn/core/types/maps/cmap/cmap2.h>
#include <cgogn/core/types/maps/cmap/cmap0.h>

#include <cgogn/core/functions/attributes.h>

#include <cgogn/geometry/types/vector_traits.h>

#include <cgogn/ui/app.h>
#include <cgogn/ui/view.h>

#include <cgogn/core/ui_modules/mesh_provider.h>
#include <cgogn/geometry/ui_modules/shrinking_ball_debugger.h>
#include <cgogn/geometry/ui_modules/alpha_inside_sphere_debugger.h>
#include <cgogn/geometry/ui_modules/udf_training.h>
#include <cgogn/rendering/ui_modules/point_cloud_render.h>
#include <cgogn/rendering/ui_modules/surface_render.h>

#include <algorithm>
#include <cctype>

#define DEFAULT_MESH_PATH CGOGN_STR(CGOGN_DATA_PATH) "/meshes/"

using namespace cgogn::numerics;

using Surface = cgogn::CMap2;
using Points = cgogn::CMap0;
using NonManifold = cgogn::IncidenceGraph;

template <typename T>
using SAttribute = typename cgogn::mesh_traits<Surface>::Attribute<T>;
using PVertex = typename cgogn::mesh_traits<Points>::Vertex;

using cgogn::geometry::Vec3;

template <typename RaySamplerTag>
int run_udf_training_app(const std::string& filename, const std::string& model_path, bool use_neural_udf,
						 bool is_surface_file, bool neural_is_mf)
{
	cgogn::thread_start();

	cgogn::ui::App app;
	app.set_window_title("UDF Training");
	app.set_window_size(1000, 800);

	cgogn::ui::MeshProvider<Surface> mps(app);
	cgogn::ui::MeshProvider<Points> mpp(app);
	cgogn::ui::MeshProvider<NonManifold> mpnm(app);

	cgogn::ui::SurfaceRender<Surface> sr(app);
	cgogn::ui::PointCloudRender<Points> pcr(app);
	cgogn::ui::SurfaceRender<NonManifold> srnm(app);
	cgogn::ui::ShrinkingBallDebugger<Points> sbd(app);
	cgogn::ui::AlphaInsideSphereDebugger<Points> aisd(app);

	cgogn::ui::UDFTraining<Surface, Points, NonManifold, RaySamplerTag> udf(app);

	app.init_modules();

	cgogn::ui::View* v1 = app.current_view();
	v1->link_module(&mps);
	v1->link_module(&mpp);
	v1->link_module(&mpnm);
	v1->link_module(&pcr);
	v1->link_module(&sr);
	v1->link_module(&srnm);
	v1->link_module(&sbd);
	v1->link_module(&aisd);
	v1->link_module(&udf);

	Points* p = nullptr;
	if (is_surface_file)
	{
		std::cout << "Detected surface mesh file. Loading as Surface Mesh..." << std::endl;
		Surface* s = mps.load_surface_from_file(filename);
		if (!s)
		{
			std::cout << "Failed to load surface mesh." << std::endl;
			return 1;
		}

		using SVertex = typename cgogn::mesh_traits<Surface>::Vertex;
		auto s_pos = cgogn::get_attribute<Vec3, SVertex>(*s, "position");
		mps.set_mesh_bb_vertex_position(*s, s_pos);
		sr.set_vertex_position(*v1, *s, s_pos);
			udf.set_selected_surface(*s);

			p = mpp.add_mesh("input_points");
			udf.set_selected_points(*p);
		}
	else
	{
		std::cout << "Loading as Point Cloud..." << std::endl;
		p = mpp.load_points_from_file(filename);
		if (!p)
		{
			std::cout << "File could not be loaded" << std::endl;
			return 1;
		}
		auto p_vertex_position = cgogn::get_attribute<Vec3, PVertex>(*p, "position");
		mpp.set_mesh_bb_vertex_position(*p, p_vertex_position);

		udf.set_selected_points(*p);
		pcr.set_vertex_position(*v1, *p, p_vertex_position);
	}

	if (use_neural_udf)
	{
		std::cout << "=== Neural UDF Mode ===" << std::endl;
		auto p_vertex_position = cgogn::get_or_add_attribute<Vec3, PVertex>(*p, "position");
		mpp.set_mesh_bb_vertex_position(*p, p_vertex_position);
		udf.set_selected_points(*p);
		udf.load_neural_udf_model(*p, model_path,
								  neural_is_mf ? cgogn::ui::UDFTraining<Surface, Points, NonManifold, RaySamplerTag>::NEURAL_MODEL_MF
											   : cgogn::ui::UDFTraining<Surface, Points, NonManifold, RaySamplerTag>::NEURAL_MODEL_UDF);

		pcr.set_vertex_position(*v1, *p, p_vertex_position);
		auto [bb_min, bb_max] = mpp.meshes_bb();
		std::cout << "Loaded point cloud bounding box: min(" << bb_min.transpose() << "), max(" << bb_max.transpose()
				  << ")" << std::endl;
		std::cout << "Neural UDF model loaded. Use UI to sample alpha-level set." << std::endl;
	}

	app.background_color_ = cgogn::rendering::GLColor(0.5f, 0.5f, 0.5f, 1.0f);
	return app.launch();
}

int main(int argc, char** argv)
{
	auto print_usage = [&](const char* exe_name) {
		std::cout << "Usage: " << exe_name << " <mesh_path> [model.pt] [--neural-type udf|mf]\n"
				  << "  <mesh_path>           Input mesh or point cloud file.\n"
				  << "  model.pt              Optional Neural UDF model path.\n"
				  << "  --neural-type udf|mf   Required when model.pt is provided.\n"
				  << "  -h, --help             Show this help message.\n";
	};

	std::string model_path;
	bool use_neural_udf = false;
	bool neural_type_set = false;
	bool neural_is_mf = false;
	std::string filename;
	bool invalid_args = false;
	std::string invalid_reason;
	if (argc < 2)
		filename = std::string(DEFAULT_MESH_PATH) + std::string("ply/546.ply");
	else
	{
		std::string first_arg = std::string(argv[1]);
		if (first_arg == "-h" || first_arg == "--help")
		{
			print_usage(argv[0]);
			return 0;
		}
		if (!first_arg.empty() && first_arg[0] == '-')
		{
			invalid_args = true;
			invalid_reason = "Missing mesh path.";
		}
		else
		{
			filename = first_arg;
		}
		for (int i = 2; i < argc; ++i)
		{
			std::string arg = std::string(argv[i]);
			if (arg == "-h" || arg == "--help")
			{
				print_usage(argv[0]);
				return 0;
			}
			if (arg.rfind("--neural-type", 0) == 0)
			{
				std::string value;
				auto eq_pos = arg.find('=');
				if (eq_pos != std::string::npos)
					value = arg.substr(eq_pos + 1);
				else if (i + 1 < argc)
					value = std::string(argv[++i]);
				else
					value.clear();

				if (value == "udf")
				{
					neural_type_set = true;
					neural_is_mf = false;
				}
				else if (value == "mf")
				{
					neural_type_set = true;
					neural_is_mf = true;
				}
				else
				{
					invalid_args = true;
					invalid_reason = "Invalid value for --neural-type (use udf or mf).";
				}
				continue;
			}
			if (arg.find(".pt") != std::string::npos)
			{
				if (!model_path.empty())
				{
					invalid_args = true;
					invalid_reason = "Multiple model paths provided.";
				}
				else
				{
					model_path = arg;
					use_neural_udf = true;
					std::cout << "Neural UDF mode: will load model from " << model_path << std::endl;
				}
				continue;
			}
			if (!arg.empty() && arg[0] == '-')
			{
				invalid_args = true;
				invalid_reason = "Unknown option: " + arg;
				continue;
			}
			invalid_args = true;
			invalid_reason = "Unexpected argument: " + arg;
		}
	}
	if (use_neural_udf && !neural_type_set)
	{
		invalid_args = true;
		invalid_reason = "Missing --neural-type for neural model.";
	}
	if (!use_neural_udf && neural_type_set)
	{
		invalid_args = true;
		invalid_reason = "--neural-type provided without model.pt.";
	}
	if (invalid_args)
	{
		std::cerr << "Invalid arguments: " << invalid_reason << std::endl;
		print_usage(argv[0]);
		return 1;
	}
	std::string ext = filename.substr(filename.find_last_of(".") + 1);
	std::transform(ext.begin(), ext.end(), ext.begin(),
				   [](unsigned char c) { return static_cast<char>(std::tolower(c)); });
	const bool is_surface_file = (ext != "ply");

	if (use_neural_udf)
		return run_udf_training_app<cgogn::geometry::RaySamplerNeural>(filename, model_path, true, is_surface_file,
																	  neural_is_mf);
	if (is_surface_file)
		return run_udf_training_app<cgogn::geometry::RaySamplerSurface>(filename, model_path, false, true, false);
	return run_udf_training_app<cgogn::geometry::RaySamplerPointCloud>(filename, model_path, false, false, false);
}
