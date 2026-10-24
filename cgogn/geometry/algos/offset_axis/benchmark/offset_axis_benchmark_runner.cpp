#include <cgogn/geometry/algos/offset_axis/benchmark/offset_axis_benchmark_runner.h>

#include <cgogn/core/types/incidence_graph/incidence_graph.h>
#include <cgogn/core/types/maps/cmap/cmap0.h>
#include <cgogn/core/types/maps/cmap/cmap2.h>
#include <cgogn/core/utils/thread_pool.h>
#if !defined(CGOGN_UDF_WITH_TORCH) || CGOGN_UDF_WITH_TORCH
#include <cgogn/geometry/algos/offset_axis/neural/neural_distance_field.h>
#endif
#include <cgogn/geometry/algos/offset_axis/reconstruction.h>
#include <cgogn/geometry/algos/offset_axis/skeleton_ply_export.h>
#include <cgogn/geometry/functions/bounding_box.h>
#include <cgogn/io/point/ply.h>
#include <cgogn/io/surface/obj.h>
#include <cgogn/io/surface/off.h>
#include <cgogn/io/surface/ply.h>
#include <cgogn/io/surface/stl.h>

#include <boost/property_tree/json_parser.hpp>
#include <boost/property_tree/ptree.hpp>

#include <algorithm>
#include <cctype>
#include <chrono>
#include <exception>
#include <filesystem>
#include <iostream>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <vector>

namespace cgogn::geometry::benchmark
{
namespace
{

double seconds_from_milliseconds(double value_ms)
{
	return value_ms / 1000.0;
}

std::string geometry_input_path_for_run(const BenchmarkConfig& config)
{
	if (config.input.geometry_type == InputGeometryType::Mesh && !config.input.surface_path.empty())
		return config.input.surface_path;
	return config.input.input_path;
}

std::string lower_extension(const std::string& path)
{
	std::string extension = std::filesystem::path(path).extension().string();
	if (!extension.empty() && extension.front() == '.')
		extension.erase(extension.begin());
	std::transform(extension.begin(), extension.end(), extension.begin(),
				   [](unsigned char c) { return static_cast<char>(std::tolower(c)); });
	return extension;
}

bool import_points(CMap0& points, const std::string& filename)
{
	if (lower_extension(filename) == "ply")
		return io::import_PLY(points, filename);
	return false;
}

bool import_surface(CMap2& surface, const std::string& filename)
{
	const std::string extension = lower_extension(filename);
	if (extension == "ply")
		return io::import_PLY(surface, filename);
	if (extension == "obj")
		return io::import_OBJ(surface, filename);
	if (extension == "off")
		return io::import_OFF(surface, filename);
	if (extension == "stl")
		return io::import_STL(surface, filename);
	return false;
}

template <typename T>
void put_option_value(boost::property_tree::ptree& tree, const char* key, const T& value)
{
	if constexpr (std::is_enum_v<T>)
		tree.put(key, to_string(value));
	else
		tree.put(key, value);
}

template <typename RaySamplerTag>
BenchmarkResult run_impl(const BenchmarkConfig& config)
{
	using Surface = CMap2;
	using Points = CMap0;
	using NonManifold = IncidenceGraph;
	using Reconstruction = UDFReconstruction<Surface, Points, NonManifold, RaySamplerTag>;
	using Status = typename Reconstruction::Status;

	BenchmarkResult result;
	result.input_mode = to_string(config.input.mode);
	result.input_path = geometry_input_path_for_run(config);
	result.seed = config.reconstruction.seed;
	Surface surface;
	Points input_points, samples, spheres;
	NonManifold skeleton;
	const bool neural = config.input.mode == InputMode::NeuralUDF;
	const auto input_start = std::chrono::high_resolution_clock::now();
	const bool use_surface = load_benchmark_input(config, surface, input_points);
	const auto input_end = std::chrono::high_resolution_clock::now();
	result.extra_timing.input_loading_ms = std::chrono::duration<double, std::milli>(input_end - input_start).count();

	typename Reconstruction::Data data;
	data.surface = use_surface ? &surface : nullptr;
	data.input_points = &input_points;
	data.samples = &samples;
	data.spheres = &spheres;
	data.skeleton = &skeleton;
	typename Reconstruction::Options options = config.reconstruction;

	Reconstruction reconstruction(data, options);
	const auto preprocessing_start = std::chrono::high_resolution_clock::now();
	if (neural)
	{
		const auto model_status =
			reconstruction.load_neural_model(config.input.neural_udf_model_path, options.neural_model_type);
		if (model_status != Status::success)
			throw std::runtime_error("Could not load neural model: " + config.input.neural_udf_model_path);
	}
	const auto preprocessing_end = std::chrono::high_resolution_clock::now();
	result.extra_timing.preprocess_ms =
		std::chrono::duration<double, std::milli>(preprocessing_end - preprocessing_start).count();

	auto log_stage = [&](const char* stage) {
		if (config.benchmark.verbose)
			std::cerr << "[BenchmarkStage] " << stage << std::endl;
	};
	log_stage("reconstruction_begin");
	const auto reconstruction_result = reconstruction.run([&](const auto& metrics) {
		if (config.benchmark.verbose)
			std::cout << "[Optimize] iteration=" << metrics.iteration << " spheres=" << metrics.sphere_count
					  << " error=" << metrics.total_error << " diff=" << metrics.error_difference << std::endl;
	});
	if (reconstruction_result.status != Status::success)
		throw std::runtime_error(
			"UDF reconstruction failed with status " + std::to_string(static_cast<int>(reconstruction_result.status)) +
			(reconstruction_result.message.empty() ? std::string() : ": " + reconstruction_result.message));
	result.extra_timing.preprocess_ms += reconstruction_result.timing.preprocessing_ms;

	result.counts = reconstruction_result.counts;
	result.timing = reconstruction_result.timing;
	const auto& topology_metrics = reconstruction.skeleton_topology().metrics();
	result.repairs = reconstruction_result.repairs;
	result.residual_prune.residual_removed_sheets = topology_metrics.residual_removed_sheets;
	result.residual_prune.residual_removed_faces = topology_metrics.residual_removed_faces;
	result.residual_prune.residual_chi_before = topology_metrics.residual_chi_before;
	result.residual_prune.residual_chi_after = topology_metrics.residual_chi_after;
	result.residual_prune.residual_fully_nm_removed_faces = topology_metrics.residual_fully_nm_removed_faces;
	result.residual_prune.residual_small_sheets = topology_metrics.residual_small_sheets;
	result.residual_prune.residual_small_sheets_guarded = topology_metrics.residual_small_sheets_guarded;
	result.residual_prune.residual_small_sheets_with_boundary_edge =
		topology_metrics.residual_small_sheets_with_boundary_edge;
	result.residual_prune.residual_small_sheets_without_boundary_edge =
		topology_metrics.residual_small_sheets_without_boundary_edge;
	result.extra_timing.average_iteration_ms =
		reconstruction_result.counts.optimization_iterations > 0
			? reconstruction_result.timing.optimization_ms / reconstruction_result.counts.optimization_iterations
			: 0.0;

	const auto export_start = std::chrono::high_resolution_clock::now();
	if (!(reconstruction.skeleton() && geometry::export_skeleton_ply(*reconstruction.skeleton(),
			reconstruction.skeleton_position_attribute(), reconstruction.skeleton_radius_attribute(),
			config.output.skeleton_ply)))
		throw std::runtime_error("Could not export skeleton: " + config.output.skeleton_ply);
	const auto counts_after_export = reconstruction.counts();
	if (counts_after_export.skeleton_vertices != result.counts.skeleton_vertices ||
		counts_after_export.skeleton_edges != result.counts.skeleton_edges ||
		counts_after_export.skeleton_faces != result.counts.skeleton_faces)
		throw std::runtime_error("Skeleton export changed the reconstructed topology.");
	const auto export_end = std::chrono::high_resolution_clock::now();
	result.extra_timing.export_ms = std::chrono::duration<double, std::milli>(export_end - export_start).count();
	return result;
}
void write_timing_json_impl(const BenchmarkConfig& config, const BenchmarkResult& result)
{
	if (config.output.timing_json.empty())
		return;

	boost::property_tree::ptree root;
	root.put("input_mode", result.input_mode);
	root.put("input_path", result.input_path);
	root.put("seed", result.seed);

	boost::property_tree::ptree benchmark_config;
	benchmark_config.put("input.geometry_type", to_string(config.input.geometry_type));
	for_each_option(config.reconstruction, [&](const OptionInfo& info, const auto& value) {
		put_option_value(benchmark_config, info.json_key, value);
	});
	benchmark_config.put("profile", config.profile.empty() ? std::string("none") : config.profile);
	benchmark_config.put("benchmark.sequential", config.benchmark.sequential);
	root.add_child("benchmark_config", benchmark_config);

	boost::property_tree::ptree counts;
	counts.put("input_vertices", result.counts.input_vertices);
	counts.put("input_points", result.counts.input_points);
	counts.put("sample_points_before_filtering", result.counts.samples_before_filtering);
	counts.put("sample_points", result.counts.samples);
	counts.put("final_spheres", result.counts.spheres);
	counts.put("skeleton_vertices", result.counts.skeleton_vertices);
	counts.put("skeleton_edges", result.counts.skeleton_edges);
	counts.put("skeleton_faces", result.counts.skeleton_faces);
	counts.put("optimization_iterations", result.counts.optimization_iterations);
	root.add_child("counts", counts);

	boost::property_tree::ptree residual_prune;
	residual_prune.put("mode", to_string(config.reconstruction.residual_prune_mode));
	residual_prune.put("residual_removed_sheets", result.residual_prune.residual_removed_sheets);
	residual_prune.put("residual_removed_faces", result.residual_prune.residual_removed_faces);
	residual_prune.put("residual_chi_before", result.residual_prune.residual_chi_before);
	residual_prune.put("residual_chi_after", result.residual_prune.residual_chi_after);
	residual_prune.put("residual_fully_nm_removed_faces", result.residual_prune.residual_fully_nm_removed_faces);
	residual_prune.put("residual_small_sheets", result.residual_prune.residual_small_sheets);
	residual_prune.put("residual_small_sheets_guarded", result.residual_prune.residual_small_sheets_guarded);
	residual_prune.put("residual_small_sheets_with_boundary_edge",
					   result.residual_prune.residual_small_sheets_with_boundary_edge);
	residual_prune.put("residual_small_sheets_without_boundary_edge",
					   result.residual_prune.residual_small_sheets_without_boundary_edge);
	root.add_child("residual_prune", residual_prune);

	boost::property_tree::ptree repairs;
	repairs.put("ma_degenerate_normals", result.repairs.ma_degenerate_normals);
	repairs.put("ma_shrink_fallbacks", result.repairs.ma_shrink_fallbacks);
	repairs.put("ma_secondary_from_nn", result.repairs.ma_secondary_from_nn);
	repairs.put("ma_flip_flipped", result.repairs.ma_flip_flipped);
	repairs.put("ma_flip_deleted", result.repairs.ma_flip_deleted);
	repairs.put("sphere_update_nonfinite", result.repairs.sphere_update_nonfinite);
	root.add_child("repairs", repairs);

	boost::property_tree::ptree timing;
	for (const auto& field : kBenchmarkTimingFields)
		timing.put(field.json_key, seconds_from_milliseconds(field.value_ms(result)));
	root.add_child("timing_s", timing);

	const std::filesystem::path output_path(config.output.timing_json);
	if (!output_path.parent_path().empty())
		std::filesystem::create_directories(output_path.parent_path());
	boost::property_tree::write_json(config.output.timing_json, root);
	if (!std::filesystem::exists(config.output.timing_json))
		throw std::runtime_error("Timing export did not create file: " + config.output.timing_json);
}

} // namespace

bool load_benchmark_input(const BenchmarkConfig& config, CMap2& surface, CMap0& input_points)
{
	using Surface = CMap2;
	using Points = CMap0;
	const std::string input_path = geometry_input_path_for_run(config);
	const bool neural = config.input.mode == InputMode::NeuralUDF;
	const bool use_surface = config.input.mode == InputMode::SurfaceMesh ||
							 (neural && config.input.geometry_type == InputGeometryType::Mesh);
	if (use_surface)
	{
		if (!import_surface(surface, input_path))
			throw std::runtime_error("Could not import surface mesh: " + input_path);
		auto position = get_attribute<Vec3, mesh_traits<Surface>::Vertex>(surface, "position");
		if (!position)
			throw std::runtime_error("Imported surface has no position attribute.");
		geometry::rescale_centered(*position, 1);
		if (neural && config.reconstruction.neural_model_type == NeuralModelType::udf)
		{
			geometry::normalize_centered(*position);
		}
	}
	else
	{
		if (!import_points(input_points, input_path))
			throw std::runtime_error("Could not import point cloud (PLY required): " + input_path);
		auto position = get_attribute<Vec3, mesh_traits<Points>::Vertex>(input_points, "position");
		if (!position)
			throw std::runtime_error("Imported point cloud has no position attribute.");
		geometry::rescale_centered(*position, 1);
		if (neural && config.reconstruction.neural_model_type == NeuralModelType::udf)
		{
			geometry::normalize_centered(*position);
		}
	}
	return use_surface;
}

BenchmarkResult SkeletonBenchmarkRunner::run(const BenchmarkConfig& config) const
{
	const auto total_start = std::chrono::high_resolution_clock::now();
	// benchmark.sequential: no thread pool worker during the run, all workers again afterwards (also on exceptions)
	struct WorkerCountGuard
	{
		bool active;
		~WorkerCountGuard()
		{
			if (active)
				thread_pool()->set_nb_workers();
		}
	};
	const WorkerCountGuard worker_count_guard{config.benchmark.sequential};
	if (config.benchmark.sequential)
		thread_pool()->set_nb_workers(0);
	BenchmarkResult result;
	switch (config.input.mode)
	{
	case InputMode::PointCloud:
		result = run_impl<geometry::RaySamplerPointCloud>(config);
		break;
	case InputMode::SurfaceMesh:
		result = run_impl<geometry::RaySamplerSurface>(config);
		break;
	case InputMode::NeuralUDF:
#if !defined(CGOGN_UDF_WITH_TORCH) || CGOGN_UDF_WITH_TORCH
		result = run_impl<geometry::RaySamplerNeural>(config);
#else
		throw std::runtime_error("udf_benchmark was built without LibTorch (CGOGN_USE_TORCH=OFF): neural_udf input is not available");
#endif
		break;
	}
	const auto total_end = std::chrono::high_resolution_clock::now();
	result.extra_timing.total_ms += std::chrono::duration<double, std::milli>(total_end - total_start).count();
	return result;
}

void write_benchmark_timing_json(const BenchmarkConfig& config, const BenchmarkResult& result)
{
	write_timing_json_impl(config, result);
}

} // namespace cgogn::geometry::benchmark
