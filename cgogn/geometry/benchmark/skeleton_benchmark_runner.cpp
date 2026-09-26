#include <cgogn/geometry/benchmark/skeleton_benchmark_runner.h>

#include <cgogn/core/types/incidence_graph/incidence_graph.h>
#include <cgogn/core/types/maps/cmap/cmap0.h>
#include <cgogn/core/types/maps/cmap/cmap2.h>
#include <cgogn/geometry/algos/udf/reconstruction.h>
#include <cgogn/geometry/benchmark/benchmark_timer.h>
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
	result.seed = config.initialization.seed;
	Surface surface;
	Points input_points, samples, spheres;
	NonManifold skeleton;
	const bool neural = config.input.mode == InputMode::NeuralUDF;
	const bool use_surface = config.input.mode == InputMode::SurfaceMesh ||
							 (neural && config.input.geometry_type == InputGeometryType::Mesh);
	const auto input_start = std::chrono::high_resolution_clock::now();
	if (use_surface)
	{
		if (!import_surface(surface, result.input_path))
			throw std::runtime_error("Could not import surface mesh: " + result.input_path);
		auto position = get_attribute<Vec3, mesh_traits<Surface>::Vertex>(surface, "position");
		if (!position)
			throw std::runtime_error("Imported surface has no position attribute.");
		geometry::rescale_centered(*position, 1);
		if (neural && config.input.neural_model_type == NeuralModelType::UDF)
		{
			geometry::normalize_centered(*position);
		}
	}
	else
	{
		if (!import_points(input_points, result.input_path))
			throw std::runtime_error("Could not import point cloud (PLY required): " + result.input_path);
		auto position = get_attribute<Vec3, mesh_traits<Points>::Vertex>(input_points, "position");
		if (!position)
			throw std::runtime_error("Imported point cloud has no position attribute.");
		geometry::rescale_centered(*position, 1);
		if (neural && config.input.neural_model_type == NeuralModelType::UDF)
		{
			geometry::normalize_centered(*position);
		}
	}
	const auto input_end = std::chrono::high_resolution_clock::now();
	result.timing.input_loading_ms = std::chrono::duration<double, std::milli>(input_end - input_start).count();

	typename Reconstruction::Data data;
	data.surface = use_surface ? &surface : nullptr;
	data.input_points = &input_points;
	data.samples = &samples;
	data.spheres = &spheres;
	data.skeleton = &skeleton;
	typename Reconstruction::Options options;
	options.neural_model_type = config.input.neural_model_type == NeuralModelType::MF
									? Reconstruction::NeuralModelType::mf
									: Reconstruction::NeuralModelType::udf;
	options.ma_flip_prune = config.input.ma_flip_prune;
	options.ma_flip_prune_alpha_factor = config.input.ma_flip_prune_alpha_factor;
	options.alpha = config.sampling.alpha;
	options.knn_k = config.sampling.knn_k;
	options.apply_filtering = config.sampling.apply_filtering;
	options.recompute_normals_after_sampling = config.sampling.recompute_normals_after_sampling;
	options.ray_sampler_batch_size = config.sampling.ray_sampler_batch_size;
	options.batch_size = config.sampling.batch_size;
	options.tolerance = config.sampling.tol;
	options.udf_max_iterations = config.sampling.udf_max_iterations;
	options.sample_radius = config.sampling.bridson.sample_radius;
	options.sqem_update_lambda_line_plane = config.optimization.sqem_update_lambda_line_plane;
	options.seed = config.initialization.seed;
	options.init_dilation_constant = config.initialization.init_dilation_constant;
	options.residual_prune = config.postprocess.residual_prune;

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
	result.timing.preprocess_ms =
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
		throw std::runtime_error("UDF reconstruction failed with status " +
								 std::to_string(static_cast<int>(reconstruction_result.status)));
	result.timing.preprocess_ms += reconstruction_result.timing.preprocessing_ms;

	result.counts.input_vertices = reconstruction_result.counts.input_vertices;
	result.counts.input_points = reconstruction_result.counts.input_points;
	result.counts.sample_points_before_filtering = reconstruction_result.counts.samples_before_filtering;
	result.counts.sample_points = reconstruction_result.counts.samples;
	result.counts.final_spheres = reconstruction_result.counts.spheres;
	result.counts.optimization_iterations = reconstruction_result.counts.optimization_iterations;
	result.timing.average_iteration_ms =
		reconstruction_result.counts.optimization_iterations > 0
			? reconstruction_result.timing.optimization_ms / reconstruction_result.counts.optimization_iterations
			: 0.0;
	result.counts.skeleton_vertices = reconstruction_result.counts.skeleton_vertices;
	result.counts.skeleton_edges = reconstruction_result.counts.skeleton_edges;
	result.counts.skeleton_faces = reconstruction_result.counts.skeleton_faces;
	result.timing.sampling_ms = reconstruction_result.timing.sampling_ms;
	result.timing.sample_filtering_ms = reconstruction_result.timing.sample_filtering_ms;
	result.timing.kdtree_bvh_ms = reconstruction_result.timing.kdtree_and_normals_ms;
	result.timing.fitting_data_ms = reconstruction_result.timing.fitting_primitives_ms;
	result.timing.initial_medial_axis_ms = reconstruction_result.timing.initial_medial_axis_ms;
	result.timing.sphere_initialization_ms = reconstruction_result.timing.sphere_initialization_ms;
	result.timing.optimization_total_ms = reconstruction_result.timing.optimization_ms;
	result.timing.cluster_total_ms = reconstruction_result.timing.cluster_ms;
	result.timing.sphere_update_total_ms = reconstruction_result.timing.sphere_update_ms;
	result.timing.error_total_ms = reconstruction_result.timing.error_ms;
	result.timing.skeleton_construction_ms = reconstruction_result.timing.skeleton_construction_ms;
	result.timing.postprocess_ms = reconstruction_result.timing.topology_processing_ms;

	if (config.output.save_face_components)
	{
		const auto component_start = std::chrono::high_resolution_clock::now();
		if (reconstruction.prepare_face_components() != Status::success)
			throw std::runtime_error("Could not prepare skeleton face components.");
		result.timing.postprocess_ms +=
			std::chrono::duration<double, std::milli>(std::chrono::high_resolution_clock::now() - component_start)
				.count();
		const auto counts = reconstruction.counts();
		result.counts.skeleton_vertices = counts.skeleton_vertices;
		result.counts.skeleton_edges = counts.skeleton_edges;
		result.counts.skeleton_faces = counts.skeleton_faces;
	}
	const auto export_start = std::chrono::high_resolution_clock::now();
	if (!reconstruction.export_skeleton_ply(config.output.skeleton_ply, config.output.save_face_components))
		throw std::runtime_error("Could not export skeleton: " + config.output.skeleton_ply);
	const auto counts_after_export = reconstruction.counts();
	if (counts_after_export.skeleton_vertices != result.counts.skeleton_vertices ||
		counts_after_export.skeleton_edges != result.counts.skeleton_edges ||
		counts_after_export.skeleton_faces != result.counts.skeleton_faces)
		throw std::runtime_error("Skeleton export changed the reconstructed topology.");
	const auto export_end = std::chrono::high_resolution_clock::now();
	result.timing.export_ms = std::chrono::duration<double, std::milli>(export_end - export_start).count();
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
	benchmark_config.put("sampling.apply_filtering", config.sampling.apply_filtering);
	benchmark_config.put("sampling.recompute_normals_after_sampling", config.sampling.recompute_normals_after_sampling);
	benchmark_config.put("initialization.init_dilation_constant", config.initialization.init_dilation_constant);
	benchmark_config.put("input.geometry_type", to_string(config.input.geometry_type));
	benchmark_config.put("input.ma_flip_prune", config.input.ma_flip_prune);
	benchmark_config.put("input.ma_flip_prune_alpha_factor", config.input.ma_flip_prune_alpha_factor);
	benchmark_config.put("output.save_face_components", config.output.save_face_components);
	benchmark_config.put("postprocess.residual_prune", config.postprocess.residual_prune);
	root.add_child("benchmark_config", benchmark_config);

	boost::property_tree::ptree counts;
	counts.put("input_vertices", result.counts.input_vertices);
	counts.put("input_points", result.counts.input_points);
	counts.put("sample_points_before_filtering", result.counts.sample_points_before_filtering);
	counts.put("sample_points", result.counts.sample_points);
	counts.put("final_spheres", result.counts.final_spheres);
	counts.put("skeleton_vertices", result.counts.skeleton_vertices);
	counts.put("skeleton_edges", result.counts.skeleton_edges);
	counts.put("skeleton_faces", result.counts.skeleton_faces);
	counts.put("optimization_iterations", result.counts.optimization_iterations);
	root.add_child("counts", counts);

	boost::property_tree::ptree timing;
	timing.put("config_loading", seconds_from_milliseconds(result.timing.config_loading_ms));
	timing.put("input_loading", seconds_from_milliseconds(result.timing.input_loading_ms));
	timing.put("preprocess", seconds_from_milliseconds(result.timing.preprocess_ms));
	timing.put("sampling", seconds_from_milliseconds(result.timing.sampling_ms));
	timing.put("sample_filtering", seconds_from_milliseconds(result.timing.sample_filtering_ms));
	timing.put("kdtree_bvh", seconds_from_milliseconds(result.timing.kdtree_bvh_ms));
	timing.put("fitting_data", seconds_from_milliseconds(result.timing.fitting_data_ms));
	timing.put("initial_medial_axis", seconds_from_milliseconds(result.timing.initial_medial_axis_ms));
	timing.put("sphere_initialization", seconds_from_milliseconds(result.timing.sphere_initialization_ms));
	timing.put("optimization_total", seconds_from_milliseconds(result.timing.optimization_total_ms));
	timing.put("cluster_total", seconds_from_milliseconds(result.timing.cluster_total_ms));
	timing.put("sphere_update_total", seconds_from_milliseconds(result.timing.sphere_update_total_ms));
	timing.put("error_total", seconds_from_milliseconds(result.timing.error_total_ms));
	timing.put("skeleton_construction", seconds_from_milliseconds(result.timing.skeleton_construction_ms));
	timing.put("postprocess", seconds_from_milliseconds(result.timing.postprocess_ms));
	timing.put("export", seconds_from_milliseconds(result.timing.export_ms));
	timing.put("total", seconds_from_milliseconds(result.timing.total_ms));
	timing.put("average_iteration", seconds_from_milliseconds(result.timing.average_iteration_ms));
	root.add_child("timing_s", timing);

	const std::filesystem::path output_path(config.output.timing_json);
	if (!output_path.parent_path().empty())
		std::filesystem::create_directories(output_path.parent_path());
	boost::property_tree::write_json(config.output.timing_json, root);
	if (!std::filesystem::exists(config.output.timing_json))
		throw std::runtime_error("Timing export did not create file: " + config.output.timing_json);
}

} // namespace

BenchmarkResult SkeletonBenchmarkRunner::run(const BenchmarkConfig& config) const
{
	const auto total_start = std::chrono::high_resolution_clock::now();
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
		result = run_impl<geometry::RaySamplerNeural>(config);
		break;
	}
	const auto total_end = std::chrono::high_resolution_clock::now();
	result.timing.total_ms += std::chrono::duration<double, std::milli>(total_end - total_start).count();
	return result;
}

void write_benchmark_timing_json(const BenchmarkConfig& config, const BenchmarkResult& result)
{
	write_timing_json_impl(config, result);
}

} // namespace cgogn::geometry::benchmark
