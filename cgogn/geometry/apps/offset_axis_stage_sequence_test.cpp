#include <cgogn/core/types/incidence_graph/incidence_graph.h>
#include <cgogn/core/types/maps/cmap/cmap0.h>
#include <cgogn/core/types/maps/cmap/cmap2.h>
#include <cgogn/core/utils/thread.h>
#if !defined(CGOGN_UDF_WITH_TORCH) || CGOGN_UDF_WITH_TORCH
#include <cgogn/geometry/algos/offset_axis/neural/neural_distance_field.h>
#endif
#include <cgogn/geometry/algos/offset_axis/reconstruction.h>
#include <cgogn/geometry/algos/offset_axis/skeleton_ply_export.h>
#include <cgogn/geometry/algos/offset_axis/benchmark/benchmark_config.h>
#include <cgogn/geometry/algos/offset_axis/benchmark/offset_axis_benchmark_runner.h>

#include <exception>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <iterator>
#include <stdexcept>
#include <string>

namespace
{

using cgogn::geometry::ReconstructionCounts;
using cgogn::geometry::Scalar;
namespace benchmark = cgogn::geometry::benchmark;
using benchmark::BenchmarkConfig;
using benchmark::InputMode;

std::string ply_prefix(const BenchmarkConfig& config)
{
	if (config.output.skeleton_ply.empty())
		return "udf_stage_sequence";
	std::string prefix = config.output.skeleton_ply;
	const std::string suffix = ".ply";
	if (prefix.size() >= suffix.size() && prefix.compare(prefix.size() - suffix.size(), suffix.size(), suffix) == 0)
		prefix.erase(prefix.size() - suffix.size());
	return prefix;
}

std::string read_binary_file(const std::string& path)
{
	std::ifstream in(path, std::ios::binary);
	if (!in)
		throw std::runtime_error("could not open file: " + path);
	return std::string(std::istreambuf_iterator<char>(in), std::istreambuf_iterator<char>());
}

// Prints one line per differing count field and returns true when all fields are equal.
bool compare_counts(const ReconstructionCounts& staged, const ReconstructionCounts& run)
{
	bool same = true;
	auto check = [&](const char* name, auto a, auto b) {
		if (a != b)
		{
			std::cout << "count " << name << ": stages=" << a << " run=" << b << std::endl;
			same = false;
		}
	};
	check("input_vertices", staged.input_vertices, run.input_vertices);
	check("input_points", staged.input_points, run.input_points);
	check("samples_before_filtering", staged.samples_before_filtering, run.samples_before_filtering);
	check("samples", staged.samples, run.samples);
	check("spheres", staged.spheres, run.spheres);
	check("skeleton_vertices", staged.skeleton_vertices, run.skeleton_vertices);
	check("skeleton_edges", staged.skeleton_edges, run.skeleton_edges);
	check("skeleton_faces", staged.skeleton_faces, run.skeleton_faces);
	check("optimization_iterations", staged.optimization_iterations, run.optimization_iterations);
	return same;
}

template <typename RaySamplerTag>
bool run_test(const BenchmarkConfig& config)
{
	using Surface = cgogn::CMap2;
	using Points = cgogn::CMap0;
	using NonManifold = cgogn::IncidenceGraph;
	using Reconstruction = cgogn::geometry::UDFReconstruction<Surface, Points, NonManifold, RaySamplerTag>;
	using Status = typename Reconstruction::Status;
	using OptimizerStatus = typename Reconstruction::SpheresOptimizerType::Status;

	const bool neural = config.input.mode == InputMode::NeuralUDF;

	auto check = [](const char* step, Status status) {
		if (status == Status::success)
			return true;
		std::cerr << "stage failed: " << step << " (status " << static_cast<int>(status) << ")" << std::endl;
		return false;
	};

	auto run_once = [&](bool staged, const std::string& ply_path, ReconstructionCounts& counts) -> bool {
		Surface surface;
		Points input_points, samples, spheres;
		NonManifold skeleton;
		const bool use_surface = benchmark::load_benchmark_input(config, surface, input_points);

		typename Reconstruction::Data data;
		data.surface = use_surface ? &surface : nullptr;
		data.input_points = &input_points;
		data.samples = &samples;
		data.spheres = &spheres;
		data.skeleton = &skeleton;
		Reconstruction reconstruction(data, config.reconstruction);

		if (staged)
		{
			// Same stage order as the "Init spheres" / "Update spheres" / "Build Skeleton" / "Topology fix" buttons of
			// the UI.
			if (!check("initialize_data", reconstruction.initialize_data()))
				return false;
			if (neural &&
				!check("load_neural_model", reconstruction.load_neural_model(config.input.neural_udf_model_path,
																			 config.reconstruction.neural_model_type)))
				return false;
			if (!check("sample", reconstruction.sample()))
				return false;
			if (!check("compute_fitting_data", reconstruction.compute_fitting_data()))
				return false;
			if (!check("initialize_spheres", reconstruction.initialize_spheres()))
				return false;
			// The UI rebuilds the skeleton right after "Init spheres" and ignores its result; the same is done here.
			reconstruction.build_skeleton();
			reconstruction.spheres_optimizer().reset_optimization();
			OptimizerStatus status = OptimizerStatus::running;
			while (status == OptimizerStatus::running)
				status =
					reconstruction.update_spheres_once(Scalar(config.reconstruction.sqem_update_lambda_line_plane));
			if (status == OptimizerStatus::failed)
			{
				std::cerr << "stage failed: update_spheres_once (optimizer status failed)" << std::endl;
				return false;
			}
			if (!check("build_skeleton", reconstruction.build_skeleton()))
				return false;
			if (!check("fix_topology_and_prune", reconstruction.fix_topology_and_prune()))
				return false;
		}
		else
		{
			if (neural &&
				!check("load_neural_model", reconstruction.load_neural_model(config.input.neural_udf_model_path,
																			 config.reconstruction.neural_model_type)))
				return false;
			const auto result = reconstruction.run();
			if (!check("run", result.status))
				return false;
		}

		const std::filesystem::path output_path(ply_path);
		if (!output_path.parent_path().empty())
			std::filesystem::create_directories(output_path.parent_path());
		if (!(reconstruction.skeleton() && cgogn::geometry::export_skeleton_ply(
											   *reconstruction.skeleton(), reconstruction.skeleton_position_attribute(),
											   reconstruction.skeleton_radius_attribute(), ply_path)))
		{
			std::cerr << "Could not export skeleton: " << ply_path << std::endl;
			return false;
		}
		counts = reconstruction.counts();
		return true;
	};

	const std::string prefix = ply_prefix(config);
	const std::string staged_ply = prefix + "_stages.ply";
	const std::string run_ply = prefix + "_run.ply";

	ReconstructionCounts staged_counts;
	ReconstructionCounts run_counts;
	if (!run_once(true, staged_ply, staged_counts) || !run_once(false, run_ply, run_counts))
	{
		std::cout << "FAIL" << std::endl;
		return false;
	}

	const std::string staged_bytes = read_binary_file(staged_ply);
	const std::string run_bytes = read_binary_file(run_ply);
	const bool counts_same = compare_counts(staged_counts, run_counts);
	const bool bytes_same = staged_bytes == run_bytes;
	if (!bytes_same)
		std::cout << "PLY bytes differ" << std::endl;

	// a stage called before the stage producing its input reports invalid_data with the failed precondition
	bool precondition_ok = false;
	{
		Surface surface;
		Points input_points, samples, spheres;
		NonManifold skeleton;
		const bool use_surface = benchmark::load_benchmark_input(config, surface, input_points);
		typename Reconstruction::Data data;
		data.surface = use_surface ? &surface : nullptr;
		data.input_points = &input_points;
		data.samples = &samples;
		data.spheres = &spheres;
		data.skeleton = &skeleton;
		Reconstruction reconstruction(data, config.reconstruction);
		reconstruction.initialize_data();
		const auto status = reconstruction.compute_fitting_data();
		precondition_ok =
			status.code == Status::invalid_data && status.message.find("precondition failed") != std::string::npos;
		std::cout << "precondition check: " << (precondition_ok ? "ok" : "unexpected") << " (" << status.message << ")"
				  << std::endl;
	}

	const bool pass = counts_same && bytes_same && precondition_ok;
	std::cout << (pass ? "PASS" : "FAIL") << std::endl;
	return pass;
}

bool dispatch_test(const BenchmarkConfig& config)
{
	switch (config.input.mode)
	{
	case InputMode::PointCloud:
		return run_test<cgogn::geometry::RaySamplerPointCloud>(config);
	case InputMode::SurfaceMesh:
		return run_test<cgogn::geometry::RaySamplerSurface>(config);
	case InputMode::NeuralUDF:
#if !defined(CGOGN_UDF_WITH_TORCH) || CGOGN_UDF_WITH_TORCH
		return run_test<cgogn::geometry::RaySamplerNeural>(config);
#else
		std::cerr << "offset_axis_stage_sequence_test was built without LibTorch (CGOGN_USE_TORCH=OFF): neural_udf input is "
					 "not available"
				  << std::endl;
		return false;
#endif
	}
	return false;
}

} // namespace

int main(int argc, char** argv)
{
	if (argc != 2)
	{
		std::cout << "Usage: " << argv[0] << " <benchmark_config.json>\n"
				  << "  <benchmark_config.json>   Benchmark configuration file (single case).\n";
		return argc == 1 ? 1 : 0;
	}

	cgogn::thread_start(0);
	try
	{
		const BenchmarkConfig config = benchmark::load_benchmark_config(argv[1]);
		const bool pass = dispatch_test(config);
		cgogn::thread_stop();
		return pass ? 0 : 1;
	}
	catch (const std::exception& e)
	{
		std::cout << "error: " << e.what() << std::endl;
		std::cout << "FAIL" << std::endl;
		cgogn::thread_stop();
		return 2;
	}
}
