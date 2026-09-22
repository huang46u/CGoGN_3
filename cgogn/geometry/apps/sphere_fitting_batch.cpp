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
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS *
 * FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License for more *
 * details.                                                                      *
 *                                                                              *
 * You should have received a copy of the GNU Lesser General Public License     *
 * along with this library; if not, write to the Free Software Foundation,     *
 * Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA             *
 *                                                                              *
 * Web site: http://cgogn.unistra.fr/                                           *
 * Contact information: cgogn@unistra.fr                                       *
 *                                                                              *
 *******************************************************************************/

#include <cgogn/core/types/incidence_graph/incidence_graph.h>
#include <cgogn/core/types/maps/cmap/cmap2.h>
#include <cgogn/core/functions/attributes.h>
#include <cgogn/core/ui_modules/mesh_provider.h>
#include <cgogn/geometry/types/vector_traits.h>
#include <cgogn/geometry/ui_modules/sphere_fitting.h>
#include <cgogn/ui/app.h>

#include <filesystem>
#include <algorithm>
#include <chrono>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <string>
#include <vector>

#if defined(_WIN32)
#include <Windows.h>
#include <crtdbg.h>
#endif

namespace fs = std::filesystem;
using Surface = cgogn::CMap2;
using Points = cgogn::CMap0;
using NonManifold = cgogn::IncidenceGraph;
using cgogn::numerics::uint32;
using SVertex = cgogn::mesh_traits<Surface>::Vertex;
using Vec3 = cgogn::geometry::Vec3;
using cgogn::ui::SphereFittingBatchOptions;
using cgogn::ui::SphereFittingBatchResult;

namespace
{

struct BatchArguments
{
	fs::path input_dir = fs::path(CGOGN_STR(CGOGN_DATA_PATH)) / "meshes" / "off";
	fs::path input_file;
	fs::path output = "sphere_fitting_batch.csv";
	uint32 timeout_seconds = 10;
	uint32 max_iterations = 1000;
	uint32 surface_samples = 300000;
	uint32 skeleton_resolution = 50;
	uint32 random_seed = 1337;
};

bool parse_uint32(const char* text, uint32* value)
{
	try
	{
		const unsigned long parsed = std::stoul(text);
		if (parsed > std::numeric_limits<uint32>::max())
			return false;
		*value = static_cast<uint32>(parsed);
		return true;
	}
	catch (const std::exception&)
	{
		return false;
	}
}

bool parse_arguments(int argc, char** argv, BatchArguments* arguments)
{
	for (int i = 1; i < argc; ++i)
	{
		const std::string option = argv[i];
		if (option == "--input-dir" && i + 1 < argc)
			arguments->input_dir = argv[++i];
		else if (option == "--mesh" && i + 1 < argc)
			arguments->input_file = argv[++i];
		else if (option == "--output" && i + 1 < argc)
			arguments->output = argv[++i];
		else if (option == "--timeout" && i + 1 < argc && parse_uint32(argv[++i], &arguments->timeout_seconds))
			{}
		else if (option == "--max-iterations" && i + 1 < argc && parse_uint32(argv[++i], &arguments->max_iterations))
			{}
		else if (option == "--surface-samples" && i + 1 < argc && parse_uint32(argv[++i], &arguments->surface_samples))
			{}
		else if (option == "--skeleton-resolution" && i + 1 < argc &&
				 parse_uint32(argv[++i], &arguments->skeleton_resolution))
			{}
		else if (option == "--seed" && i + 1 < argc && parse_uint32(argv[++i], &arguments->random_seed))
			{}
		else
		{
			std::cerr << "Unknown or invalid option: " << option << std::endl;
			return false;
		}
	}
	return true;
}

void write_csv_value(std::ostream& stream, const std::string& value)
{
	stream << '"';
	for (const char character : value)
	{
		if (character == '"')
			stream << '"';
		stream << character;
	}
	stream << '"';
}

void write_result(std::ostream& stream, const std::string& mesh_name, uint32 target, const std::string& metric,
	bool use_line_quadric, const SphereFittingBatchOptions& options, const SphereFittingBatchResult& result,
	const std::string& status)
{
	write_csv_value(stream, mesh_name);
	stream << ',' << target << ',';
	write_csv_value(stream, metric);
	stream << ',' << (use_line_quadric ? 1 : 0) << ",1,1,1," << options.timeout_seconds << ','
		   << options.max_iterations << ',' << options.hausdorff_surface_samples << ','
		   << options.hausdorff_skeleton_resolution << ',' << options.random_seed << ',' << result.sphere_count << ','
		   << result.iterations << ',' << std::setprecision(17) << result.surface_to_skeleton << ','
		   << result.skeleton_to_surface << ',' << result.symmetric << ',' << result.surface_sample_count << ','
		   << result.skeleton_sample_count << ',' << result.preprocess_seconds << ','
		   << result.sphere_init_seconds << ',' << result.fit_seconds << ','
		   << result.eval_seconds << ',' << result.total_seconds << ',';
	write_csv_value(stream, status);
	stream << '\n';
}

std::string result_status(const SphereFittingBatchResult& result)
{
	if (result.converged)
		return "converged";
	if (result.timed_out)
		return "timeout";
	if (result.iteration_limited)
		return "max_iterations";
	return "no_result";
}

void write_failure_rows(std::ostream& stream, const std::string& mesh_name, const std::vector<uint32>& targets,
	const BatchArguments& arguments, const std::string& status)
{
	for (const uint32 target : targets)
		for (const bool use_line_quadric : {false, true})
		{
			SphereFittingBatchOptions options;
			options.target_spheres = target;
			options.timeout_seconds = arguments.timeout_seconds;
			options.max_iterations = arguments.max_iterations;
			options.hausdorff_surface_samples = arguments.surface_samples;
			options.hausdorff_skeleton_resolution = arguments.skeleton_resolution;
			options.random_seed = arguments.random_seed;
			options.use_line_quadric = use_line_quadric;
			write_result(stream, mesh_name, target,
				use_line_quadric ? "sqem_line_quadric" : "sqem_euclidean", use_line_quadric, options,
				SphereFittingBatchResult(), status);
		}
}

} // namespace

int main(int argc, char** argv)
{
#if defined(_WIN32)
	SetErrorMode(SEM_FAILCRITICALERRORS | SEM_NOGPFAULTERRORBOX | SEM_NOOPENFILEERRORBOX);
	_set_abort_behavior(0, _WRITE_ABORT_MSG | _CALL_REPORTFAULT);
#endif

	BatchArguments arguments;
	if (!parse_arguments(argc, argv, &arguments))
		return 2;
	if (!arguments.input_file.empty() && !fs::is_regular_file(arguments.input_file))
	{
		std::cerr << "Input mesh does not exist: " << arguments.input_file << std::endl;
		return 2;
	}
	if (arguments.input_file.empty() && !fs::is_directory(arguments.input_dir))
	{
		std::cerr << "Input directory does not exist: " << arguments.input_dir << std::endl;
		return 2;
	}
	if (arguments.output.has_parent_path())
		fs::create_directories(arguments.output.parent_path());

	std::ofstream output(arguments.output);
	if (!output)
	{
		std::cerr << "Could not open output CSV: " << arguments.output << std::endl;
		return 2;
	}
	output << "mesh,target_spheres,metric,use_line_quadric,auto_split,auto_stop,sphere_correction,timeout_seconds,"
		      "max_iterations,hausdorff_surface_samples,hausdorff_skeleton_resolution,random_seed,actual_spheres,"
		      "iterations,surface_to_skeleton,skeleton_to_surface,symmetric,surface_sample_count,"
		      "skeleton_sample_count,preprocess_seconds,sphere_init_seconds,fit_seconds,eval_seconds,total_seconds,status\n";

	cgogn::thread_start();
	cgogn::ui::App app;
	app.show_gui(false);
	cgogn::ui::MeshProvider<Surface> surface_provider(app);
	cgogn::ui::MeshProvider<Points> points_provider(app);
	cgogn::ui::MeshProvider<NonManifold> non_manifold_provider(app);
	cgogn::ui::SphereFitting<Surface, Points, NonManifold> sphere_fitting(app);
	app.init_modules();

	std::vector<fs::path> mesh_files;
	if (!arguments.input_file.empty())
		mesh_files.push_back(arguments.input_file);
	else
		for (const fs::directory_entry& entry : fs::directory_iterator(arguments.input_dir))
			if (entry.is_regular_file() && entry.path().extension() == ".off")
				mesh_files.push_back(entry.path());
	std::sort(mesh_files.begin(), mesh_files.end());

	const std::vector<uint32> targets = {50, 100, 200, 250};
	for (const fs::path& mesh_file : mesh_files)
	{
		std::cout << "Evaluating " << mesh_file.filename().string() << std::endl;
		Surface* surface = surface_provider.load_surface_from_file(mesh_file.string());
		if (!surface)
		{
			write_failure_rows(output, mesh_file.filename().string(), targets, arguments, "load_failed");
			continue;
		}

		auto position = cgogn::get_attribute<Vec3, SVertex>(*surface, "position");
		if (!position)
		{
			std::cerr << "Position attribute missing: " << mesh_file << std::endl;
			write_failure_rows(output, mesh_file.filename().string(), targets, arguments, "position_missing");
			continue;
		}
		sphere_fitting.set_selected_surface(*surface);
		sphere_fitting.set_surface_vertex_position(*surface, position);
		const auto preprocess_start = std::chrono::steady_clock::now();
		sphere_fitting.init_surface_data(*surface);
		const double preprocess_seconds =
			std::chrono::duration<double>(std::chrono::steady_clock::now() - preprocess_start).count();

		for (const uint32 target : targets)
			for (const bool use_line_quadric : {false, true})
			{
				SphereFittingBatchOptions options;
				options.target_spheres = target;
				options.timeout_seconds = arguments.timeout_seconds;
				options.max_iterations = arguments.max_iterations;
				options.hausdorff_surface_samples = arguments.surface_samples;
				options.hausdorff_skeleton_resolution = arguments.skeleton_resolution;
				options.random_seed = arguments.random_seed;
				options.use_line_quadric = use_line_quadric;
				SphereFittingBatchResult result;
				std::string status;
				try
				{
					result = sphere_fitting.run_batch(*surface, options);
					result.preprocess_seconds = preprocess_seconds;
					status = result_status(result);
				}
				catch (const std::exception& exception)
				{
					status = std::string("exception:") + exception.what();
				}
				write_result(output, mesh_file.filename().string(), target,
					use_line_quadric ? "sqem_line_quadric" : "sqem_euclidean", use_line_quadric, options, result,
					status);
				output.flush();
			}
	}
	return 0;
}
