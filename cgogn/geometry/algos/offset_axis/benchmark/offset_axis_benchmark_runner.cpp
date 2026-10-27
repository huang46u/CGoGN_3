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
#include <cgogn/geometry/algos/offset_axis/benchmark/offset_axis_benchmark_runner.h>

#include <cgogn/core/functions/attributes.h>
#include <cgogn/core/functions/mesh_info.h>
#include <cgogn/core/types/incidence_graph/incidence_graph.h>
#include <cgogn/core/types/maps/cmap/cmap0.h>
#include <cgogn/core/types/maps/cmap/cmap2.h>
#include <cgogn/core/utils/thread_pool.h>
#if !defined(CGOGN_OFFSET_AXIS_WITH_TORCH) || CGOGN_OFFSET_AXIS_WITH_TORCH
#include <cgogn/geometry/algos/offset_axis/neural/neural_distance_field.h>
#endif
#include <cgogn/geometry/algos/offset_axis/alpha_medial.h>
#include <cgogn/geometry/algos/offset_axis/alpha_sampling.h>
#include <cgogn/geometry/algos/offset_axis/input_geometry.h>
#include <cgogn/geometry/algos/offset_axis/topology_scores.h>
#include <cgogn/geometry/algos/offset_axis/point_cloud_distance_field.h>
#include <cgogn/geometry/algos/offset_axis/reconstruction_options.h>
#include <cgogn/geometry/algos/offset_axis/sample_processing.h>
#include <cgogn/geometry/algos/offset_axis/sampling_filter.h>
#include <cgogn/geometry/algos/offset_axis/skeleton_topology.h>
#include <cgogn/geometry/algos/offset_axis/spheres_optimizer.h>
#include <cgogn/geometry/algos/offset_axis/surface_distance_field.h>
#include <cgogn/geometry/types/line_quadric.h>
#include <cgogn/geometry/types/spherical_quadric.h>
#include <cgogn/io/surface/skeleton_ply.h>
#include <cgogn/geometry/functions/bounding_box.h>
#include <cgogn/io/point/ply.h>
#include <cgogn/io/surface/obj.h>
#include <cgogn/io/surface/off.h>
#include <cgogn/io/surface/ply.h>
#include <cgogn/io/surface/stl.h>

#include <boost/property_tree/json_parser.hpp>
#include <boost/property_tree/ptree.hpp>

#include <libacc/kd_tree.h>

#include <algorithm>
#include <cctype>
#include <chrono>
#include <exception>
#include <filesystem>
#include <iostream>
#include <memory>
#include <set>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <vector>

namespace cgogn
{

namespace geometry
{

namespace offset_axis
{

namespace benchmark
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

template <template <typename, typename> class FIELD>
BenchmarkResult run_impl(const BenchmarkConfig& config)
{
	using Surface = CMap2;
	using Points = CMap0;
	using NonManifold = IncidenceGraph;
	using PVertex = typename mesh_traits<Points>::Vertex;
	using SVertex = typename mesh_traits<Surface>::Vertex;
	using NMVertex = typename mesh_traits<NonManifold>::Vertex;
	using NMEdge = typename mesh_traits<NonManifold>::Edge;
	using NMFace = typename mesh_traits<NonManifold>::Face;
	using Field = FIELD<Surface, Points>;
	using Optimizer = SpheresOptimizer<Points>;
	using Topology = SkeletonTopology<Points, NonManifold>;
	using Clock = std::chrono::high_resolution_clock;

	auto elapsed_ms = [](Clock::time_point start) {
		return std::chrono::duration<double, std::milli>(Clock::now() - start).count();
	};

	BenchmarkResult result;
	result.input_mode = to_string(config.input.mode);
	result.input_path = geometry_input_path_for_run(config);
	BenchmarkTiming& timing = result.timing;
	Surface surface;
	Points input_points, samples, spheres;
	NonManifold skeleton;
	const bool neural = config.input.mode == InputMode::NeuralUDF;
	const auto loading_start = Clock::now();
	const bool use_surface = load_benchmark_input(config, surface, input_points);
	Surface* input_surface = use_surface ? &surface : nullptr;

	const ReconstructionOptions& options = config.reconstruction;
	InputGeometry<Surface, Points> input;
	Field field{input, options};
	std::unique_ptr<acc::KDTree<3, uint32>> sample_kdtree, sample_ma_kdtree;
	std::vector<PVertex> sample_vertices, sample_ma_vertices;

	if (neural && !field.load_model(config.input.neural_udf_model_path, options.neural_model_type))
		throw std::runtime_error("Could not load neural model: " + config.input.neural_udf_model_path);

	auto input_position = get_or_add_attribute<Vec3, PVertex>(input_points, "position");
	auto input_normal = get_or_add_attribute<Vec3, PVertex>(input_points, "normal");
	auto input_knn = get_or_add_attribute<std::vector<PVertex>, PVertex>(input_points, "knn");
	auto sample_position = get_or_add_attribute<Vec3, PVertex>(samples, "position");
	auto sample_normal = get_or_add_attribute<Vec3, PVertex>(samples, "normal");
	auto sample_area = get_or_add_attribute<Scalar, PVertex>(samples, "area");
	auto sample_knn = get_or_add_attribute<std::vector<PVertex>, PVertex>(samples, "knn");
	auto sample_quadric = get_or_add_attribute<Spherical_Quadric, PVertex>(samples, "quadric");
	auto sample_line_quadric = get_or_add_attribute<Line_Quadric, PVertex>(samples, "line_quadric");
	auto sample_ma_position = get_or_add_attribute<Vec3, PVertex>(samples, "ma_position");
	auto sample_ma_radius = get_or_add_attribute<Scalar, PVertex>(samples, "ma_radius");
	auto sample_ma_secondary = get_or_add_attribute<PVertex, PVertex>(samples, "ma_secondary_vertex");
	auto sample_sphere = get_or_add_attribute<PVertex, PVertex>(samples, "sphere");
	auto sample_error = get_or_add_attribute<Scalar, PVertex>(samples, "error");
	auto sphere_position = get_or_add_attribute<Vec3, PVertex>(spheres, "position");
	auto sphere_radius = get_or_add_attribute<Scalar, PVertex>(spheres, "radius");
	auto sphere_cluster = get_or_add_attribute<std::vector<PVertex>, PVertex>(spheres, "cluster");
	auto sphere_cluster_area = get_or_add_attribute<Scalar, PVertex>(spheres, "cluster_area");
	auto sphere_cluster_color = get_or_add_attribute<Vec4, PVertex>(spheres, "cluster_color");
	auto sphere_neighbors = get_or_add_attribute<std::vector<PVertex>, PVertex>(spheres, "neighbor_clusters");
	auto sphere_error = get_or_add_attribute<Scalar, PVertex>(spheres, "error");
	auto sphere_error_raw = get_or_add_attribute<Scalar, PVertex>(spheres, "error_not_normalized");
	auto sphere_skeleton_vertex = get_or_add_attribute<NMVertex, PVertex>(spheres, "skeleton_vertex");
	auto skeleton_position = get_or_add_attribute<Vec3, NMVertex>(skeleton, "position");
	auto skeleton_radius = get_or_add_attribute<Scalar, NMVertex>(skeleton, "radius");
	auto skeleton_source_sphere = get_or_add_attribute<PVertex, NMVertex>(skeleton, "source_sphere");
	auto skeleton_incident_tets = get_or_add_attribute<std::set<std::size_t>, NMFace>(skeleton, "incident_tets");
	auto skeleton_edge_degree = get_or_add_attribute<uint32, NMEdge>(skeleton, "degree");
	auto skeleton_face_score = get_or_add_attribute<Scalar, NMFace>(skeleton, "udf_score");
	auto skeleton_edge_score = get_or_add_attribute<Scalar, NMEdge>(skeleton, "udf_score");
	input.bind(input_surface, &input_points, input_position, input_normal);
	timing.loading_ms = elapsed_ms(loading_start);

	auto build_sample_kdtrees = [&]() {
		build_point_kdtree(samples, *sample_position, sample_kdtree, sample_vertices);
		build_medial_kdtree(samples, *sample_ma_position, *sample_ma_radius, sample_ma_kdtree, sample_ma_vertices);
	};

	auto stage_start = Clock::now();
	input.build_input_kdtree();
	const PrepareInputs<Surface, Points> inputs{input_surface, &input_points, input_position, input_normal, input_knn};
	if (field.prepare(inputs) != FieldPrepareResult::ok)
		throw std::runtime_error("Could not prepare the distance field.");
	if (sample_alpha_level_set(field, inputs, alpha_sampling_params(options), samples, sample_position.get(),
							   sample_normal.get()) != SamplingStatus::ok)
		throw std::runtime_error("Alpha level set sampling failed.");
	result.counts.samples_before_filtering = static_cast<uint32>(nb_cells<PVertex>(samples));
	build_sample_kdtrees();

	if (options.apply_filtering &&
		filter_samples<Field>(input, Scalar(options.alpha), samples, sample_position.get(), sample_normal.get()) > 0)
	{
		if (nb_cells<PVertex>(samples) == 0)
			throw std::runtime_error("Every sample is farther than alpha + filter_margin from the input (the field does "
									 "not match the input).");
		build_sample_kdtrees();
	}

	// PCA for every input: neural samples already carry the network gradient, and for MF models that gradient is not
	// orthogonal to the medial axis.
	if (options.recompute_normals_after_sampling)
		recompute_normals_pca(samples, *sample_position, *sample_normal, *sample_kdtree, sample_vertices,
							  options.knn_k);
	timing.sampling_ms = elapsed_ms(stage_start);

	// the medial flip prune changes normals and removes samples, so it runs before the primitives that use them
	stage_start = Clock::now();
	compute_initial_medial_axis(field, samples, *sample_position, *sample_normal, sample_kdtree, sample_vertices,
								*sample_ma_position, *sample_ma_radius, *sample_ma_secondary, Scalar(options.alpha),
								options.ma_flip_prune, Scalar(options.ma_flip_prune_alpha_factor), options.knn_k);
	build_sample_kdtrees();
	if (nb_cells<PVertex>(samples) == 0)
		throw std::runtime_error("The medial flip prune removed every sample.");

	compute_samples_area(samples, *sample_position, *sample_normal, *sample_area, *sample_knn, *sample_kdtree,
						 sample_vertices, options.knn_k, Scalar(options.alpha));
	compute_quadrics(samples, *sample_position, *sample_normal, *sample_area, *sample_knn, *sample_quadric,
					 *sample_line_quadric, options.knn_k);

	typename Optimizer::Data optimizer_data{samples, spheres, sample_kdtree, sample_vertices};
	optimizer_data.sample_position = sample_position.get();
	optimizer_data.sample_area = sample_area.get();
	optimizer_data.sample_knn = sample_knn.get();
	optimizer_data.sample_quadric = sample_quadric.get();
	optimizer_data.sample_line_quadric = sample_line_quadric.get();
	optimizer_data.sample_ma_position = sample_ma_position.get();
	optimizer_data.sample_ma_radius = sample_ma_radius.get();
	optimizer_data.sample_ma_secondary_vertex = sample_ma_secondary.get();
	optimizer_data.sample_sphere = sample_sphere.get();
	optimizer_data.sample_error = sample_error.get();
	optimizer_data.sphere_position = sphere_position.get();
	optimizer_data.sphere_radius = sphere_radius.get();
	optimizer_data.sphere_cluster = sphere_cluster.get();
	optimizer_data.sphere_cluster_area = sphere_cluster_area.get();
	optimizer_data.sphere_cluster_color = sphere_cluster_color.get();
	optimizer_data.sphere_neighbors = sphere_neighbors.get();
	optimizer_data.sphere_error = sphere_error.get();
	optimizer_data.sphere_error_not_normalized = sphere_error_raw.get();
	Optimizer optimizer(optimizer_data);
	optimizer.initialize_from_samples(Scalar(options.init_dilation_constant), options.init_owner_euclidean,
									  Scalar(options.sqem_update_lambda_line_plane));
	timing.initialization_ms = elapsed_ms(stage_start);
	if (nb_cells<PVertex>(spheres) == 0)
		throw std::runtime_error("The sphere initialization created no sphere.");

	stage_start = Clock::now();
	optimizer.reset_optimization();
	typename Optimizer::Status status = Optimizer::Status::running;
	while (status == Optimizer::Status::running)
	{
		status = optimizer.update_once(Scalar(options.sqem_update_lambda_line_plane));
		const auto& metrics = optimizer.metrics();
		if (config.benchmark.verbose)
			std::cout << "[Optimize] iteration=" << metrics.iteration << " spheres=" << metrics.sphere_count
					  << " error=" << metrics.total_error << " diff=" << metrics.error_difference << std::endl;
	}
	timing.optimization_ms = elapsed_ms(stage_start);
	if (status == Optimizer::Status::failed)
		throw std::runtime_error("The sphere optimization failed.");

	stage_start = Clock::now();
	optimizer.compute_sphere_neighbors();
	typename Topology::Data topology_data{skeleton, optimizer};
	topology_data.sphere_skeleton_vertex = sphere_skeleton_vertex.get();
	topology_data.skeleton_position = skeleton_position.get();
	topology_data.skeleton_radius = skeleton_radius.get();
	topology_data.skeleton_source_sphere = skeleton_source_sphere.get();
	topology_data.face_incident_tets = skeleton_incident_tets.get();
	topology_data.edge_degree = skeleton_edge_degree.get();
	topology_data.face_score = skeleton_face_score.get();
	topology_data.edge_score = skeleton_edge_score.get();
	Topology topology(topology_data);
	topology.build_from_spheres();
	if (nb_cells<NMVertex>(skeleton) == 0)
		throw std::runtime_error("The skeleton construction created no vertex.");

	if (!topology.fix_topology([&](const std::vector<Vec3>& points, std::vector<Scalar>& values) {
			return evaluate_topology_scores(field, samples, *sample_ma_position, *sample_ma_kdtree, sample_ma_vertices,
											points, values);
		}))
		throw std::runtime_error("The topology fix failed.");
	if (options.residual_prune)
		topology.prune_residual_sheets();
	timing.postprocess_ms = elapsed_ms(stage_start);

	ReconstructionCounts& counts = result.counts;
	auto count_skeleton = [&](ReconstructionCounts& c) {
		c.skeleton_vertices = nb_cells<NMVertex>(skeleton);
		c.skeleton_edges = nb_cells<NMEdge>(skeleton);
		c.skeleton_faces = nb_cells<NMFace>(skeleton);
	};
	counts.input_points = nb_cells<PVertex>(input_points);
	counts.input_vertices = use_surface ? nb_cells<SVertex>(surface) : counts.input_points;
	counts.samples = nb_cells<PVertex>(samples);
	counts.spheres = nb_cells<PVertex>(spheres);
	count_skeleton(counts);
	counts.optimization_iterations = optimizer.metrics().iteration;

	const auto export_start = Clock::now();
	if (!io::export_skeleton_ply(skeleton, skeleton_position, skeleton_radius, config.output.skeleton_ply))
		throw std::runtime_error("Could not export skeleton: " + config.output.skeleton_ply);
	ReconstructionCounts counts_after_export;
	count_skeleton(counts_after_export);
	if (counts_after_export.skeleton_vertices != counts.skeleton_vertices ||
		counts_after_export.skeleton_edges != counts.skeleton_edges ||
		counts_after_export.skeleton_faces != counts.skeleton_faces)
		throw std::runtime_error("Skeleton export changed the reconstructed topology.");
	timing.export_ms = elapsed_ms(export_start);
	return result;
}
void write_timing_json_impl(const BenchmarkConfig& config, const BenchmarkResult& result)
{
	if (config.output.timing_json.empty())
		return;

	boost::property_tree::ptree root;
	root.put("input_mode", result.input_mode);
	root.put("input_path", result.input_path);

	boost::property_tree::ptree benchmark_config;
	benchmark_config.put("input.geometry_type", to_string(config.input.geometry_type));
	for_each_option(config.reconstruction, [&](const OptionInfo& info, const auto& value) {
		put_option_value(benchmark_config, info.json_key, value);
	});
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

	boost::property_tree::ptree timing;
	for (const auto& field : benchmark_timing_fields)
		timing.put(field.json_key, seconds_from_milliseconds(result.timing.*field.value_ms));
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
		result = run_impl<offset_axis::PointCloudDistanceField>(config);
		break;
	case InputMode::SurfaceMesh:
		result = run_impl<offset_axis::SurfaceDistanceField>(config);
		break;
	case InputMode::NeuralUDF:
#if !defined(CGOGN_OFFSET_AXIS_WITH_TORCH) || CGOGN_OFFSET_AXIS_WITH_TORCH
		result = run_impl<offset_axis::NeuralDistanceField>(config);
#else
		throw std::runtime_error("offset_axis_benchmark was built without LibTorch (CGOGN_USE_TORCH=OFF): neural_udf input is not available");
#endif
		break;
	}
	const auto total_end = std::chrono::high_resolution_clock::now();
	result.timing.total_ms = std::chrono::duration<double, std::milli>(total_end - total_start).count();
	return result;
}

void write_benchmark_timing_json(const BenchmarkConfig& config, const BenchmarkResult& result)
{
	write_timing_json_impl(config, result);
}

} // namespace benchmark

} // namespace offset_axis

} // namespace geometry

} // namespace cgogn
