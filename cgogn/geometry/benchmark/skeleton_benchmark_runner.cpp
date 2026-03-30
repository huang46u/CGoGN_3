#include <cgogn/geometry/benchmark/skeleton_benchmark_runner.h>

#include <cgogn/geometry/algos/udf_pipeline_core.h>
#include <cgogn/geometry/benchmark/benchmark_timer.h>

#include <boost/property_tree/json_parser.hpp>
#include <boost/property_tree/ptree.hpp>

#include <chrono>
#include <filesystem>
#include <iostream>
#include <sstream>
#include <stdexcept>

namespace cgogn
{

namespace geometry
{

namespace benchmark
{

namespace
{

double seconds_from_milliseconds(double value_ms)
{
	return value_ms / 1000.0;
}

class ScopedStreamSilencer
{
public:
	explicit ScopedStreamSilencer(bool enabled) : enabled_(enabled)
	{
		if (!enabled_)
			return;
		old_cout_ = std::cout.rdbuf(null_stream_.rdbuf());
		old_cerr_ = std::cerr.rdbuf(null_stream_.rdbuf());
	}

	~ScopedStreamSilencer()
	{
		if (!enabled_)
			return;
		std::cout.rdbuf(old_cout_);
		std::cerr.rdbuf(old_cerr_);
	}

private:
	bool enabled_ = false;
	std::streambuf* old_cout_ = nullptr;
	std::streambuf* old_cerr_ = nullptr;
	std::ostringstream null_stream_;
};

template <typename Training>
typename Training::HeadlessBenchmarkOptions make_training_options(const BenchmarkConfig& config)
{
	typename Training::HeadlessBenchmarkOptions options;
	options.verbose_ = config.benchmark.verbose;
	options.initial_nb_spheres_ = config.initialization.initial_nb_spheres;
	switch (config.input.initial_ma_mode_override)
	{
	case InitialMAMode::Displacement:
		options.initial_ma_mode_override_ = Training::INITIAL_MA_DISPLACEMENT;
		break;
	case InitialMAMode::ShrinkingBall:
		options.initial_ma_mode_override_ = Training::INITIAL_MA_SHRINKING_BALL;
		break;
	case InitialMAMode::Auto:
	default:
		options.initial_ma_mode_override_ = Training::INITIAL_MA_AUTO;
		break;
	}
	options.sphere_correction_ = config.optimization.sphere_correction;
	options.sphere_correction_mode_ = (config.optimization.sphere_correction_mode == SphereCorrectionMode::OnSplit)
										  ? Training::CORRECT_ON_SPLIT
										  : Training::CORRECT_ALWAYS;
	options.lock_skeleton_connectivity_ = config.optimization.lock_skeleton_connectivity;
	options.distance_mode_ = (config.optimization.distance_mode == DistanceMode::LineQuadricDistanceFreeRadius)
								 ? Training::LINE_QUADRIC_DISTANCE_FREE_RADIUS
								 : Training::LINE_QUADRIC_DISTANCE;
	options.use_local_clusters_ = config.optimization.use_local_clusters;
	options.local_cluster_connectivity_refresh_interval_ =
		config.optimization.local_cluster_connectivity_refresh_interval;
	options.auto_stop_ = config.optimization.auto_stop;
	options.max_iterations_without_autosplit_ = config.optimization.max_iterations_without_autosplit;
	options.max_iterations_after_reaching_max_spheres_ = config.optimization.max_iterations_after_reaching_max_spheres;
	options.auto_split_ = config.auto_split.enabled;
	options.auto_split_mode_ =
		(config.auto_split.mode == AutoSplitMode::MaxNbSpheres) ? Training::MAX_NB_SPHERES : Training::ERROR_THRESHOLD;
	options.auto_split_error_threshold_ = config.auto_split.error_threshold;
	options.auto_split_max_nb_spheres_ = config.auto_split.max_nb_spheres;
	options.auto_split_ratio_ = config.auto_split.ratio;
	options.auto_split_max_per_iter_error_ = config.auto_split.max_per_iter_error;
	options.auto_split_max_per_iter_max_ = config.auto_split.max_per_iter_max;
	options.sqem_update_lambda_full_ = config.optimization.sqem_update_lambda_full;
	options.sqem_update_lambda_line_plane_ = config.optimization.sqem_update_lambda_line_plane;
	options.sqem_fix_radius_scale_ = config.optimization.sqem_fix_radius_scale;
	options.udf_center_enabled_ = config.optimization.udf_center_enabled;
	options.udf_center_lambda_ = config.optimization.udf_center_lambda;
	options.init_dilation_constant_ = config.initialization.init_dilation_constant;
	options.init_min_cover_points_ = config.initialization.init_min_cover_points;
	options.alpha_ = config.sampling.alpha;
	options.knn_k_ = config.sampling.knn_k;
	options.seed_ = config.initialization.seed;
	options.grid_cell_size_ = config.sampling.grid_cell_size;
	options.num_alpha_samples_ = config.sampling.num_alpha_samples;
	options.batch_size_ = config.sampling.batch_size;
	options.ray_sampler_batch_size_ = config.sampling.ray_sampler_batch_size;
	options.tol_ = config.sampling.tol;
	options.udf_max_iterations_ = config.sampling.udf_max_iterations;
	return options;
}

template <typename RaySamplerTag>
BenchmarkResult run_impl(const BenchmarkConfig& config)
{
	using Context = UDFPipelineContext<RaySamplerTag>;
	using Core = UDFPipelineCore<RaySamplerTag>;
	using Training = typename Context::Training;
	using Points = typename Context::Points;

	auto log_stage = [&](const char* stage) {
		if (config.benchmark.verbose)
			std::cerr << "[BenchmarkStage] " << stage << std::endl;
	};

	log_stage("construct_context_begin");
	Context context;
	log_stage("construct_context_end");

	Core core(context);
	BenchmarkResult result;
	result.input_mode = to_string(config.input.mode);
	result.input_path = !config.input.input_path.empty() ? config.input.input_path : config.input.surface_path;
	result.seed = config.initialization.seed;

	Points* points = nullptr;
	{
		ScopedStreamSilencer silence_internal_output(!config.benchmark.verbose);

		log_stage("input_loading_begin");
		{
			ScopedBenchmarkTimer timer(result.timing.input_loading_ms);
			const bool normalize_for_udf =
				(config.input.mode == InputMode::NeuralUDF && config.input.neural_model_type == NeuralModelType::UDF);
			switch (config.input.mode)
			{
			case InputMode::PointCloud:
				points = &core.load_point_cloud_input(config.input.input_path, normalize_for_udf);
				break;
			case InputMode::SurfaceMesh:
				core.load_surface_input(config.input.input_path, normalize_for_udf);
				points = &core.create_empty_points_input("input_points");
				break;
			case InputMode::NeuralUDF:
				if (!config.input.surface_path.empty())
				{
					core.load_surface_input(config.input.surface_path, normalize_for_udf);
					points = &core.create_empty_points_input("input_points");
				}
				else
				{
					points = &core.load_point_cloud_input(config.input.input_path, normalize_for_udf);
				}
				break;
			}
		}
		log_stage("input_loading_end");

		if (!points)
			throw std::runtime_error("Benchmark pipeline did not initialize a valid point container.");
		core.prepare_points(*points);

		const auto options = make_training_options<Training>(config);
		log_stage("preprocess_begin");
		{
			ScopedBenchmarkTimer timer(result.timing.preprocess_ms);
			core.apply_options_prepared(*points, options);
			if (config.input.mode == InputMode::NeuralUDF)
			{
				const auto model_type = (config.input.neural_model_type == NeuralModelType::MF)
											? Training::NEURAL_MODEL_MF
											: Training::NEURAL_MODEL_UDF;
				core.load_neural_model(*points, config.input.neural_udf_model_path, model_type);
			}
		}
		log_stage("preprocess_end");

		log_stage("sampling_begin");
		{
			ScopedBenchmarkTimer timer(result.timing.sampling_ms);
			core.sample_alpha_level_set_prepared(*points, options);
		}
		log_stage("sampling_end");

		{
			const auto sampled_counts = core.collect_counts(*points);
			result.counts.sample_points_before_filtering = sampled_counts.sample_points_;
		}

		if (config.sampling.apply_filtering)
		{
			log_stage("sample_filtering_begin");
			ScopedBenchmarkTimer timer(result.timing.sample_filtering_ms);
			core.apply_sampling_filtering_prepared(*points);
			log_stage("sample_filtering_end");
		}
		else
		{
			result.timing.sample_filtering_ms = 0.0;
		}

		log_stage("kdtree_begin");
		{
			ScopedBenchmarkTimer timer(result.timing.kdtree_bvh_ms);
			core.build_kdtree_and_normals_prepared(*points);
		}
		log_stage("kdtree_end");

		log_stage("fitting_begin");
		{
			ScopedBenchmarkTimer timer(result.timing.fitting_data_ms);
			core.compute_fitting_primitives_prepared(*points);
		}
		log_stage("fitting_end");

		log_stage("initial_ma_begin");
		{
			ScopedBenchmarkTimer timer(result.timing.initial_medial_axis_ms);
			core.compute_initial_medial_axis_prepared(*points);
		}
		log_stage("initial_ma_end");

		log_stage("sphere_init_begin");
		{
			ScopedBenchmarkTimer timer(result.timing.sphere_initialization_ms);
			core.init_spheres_prepared(*points, options.initial_nb_spheres_);
		}
		log_stage("sphere_init_end");

		log_stage("optimization_begin");
		{
			const auto stats = core.optimize_prepared(*points, config.benchmark.verbose);
			result.timing.optimization_total_ms = stats.optimization_total_ms_;
			result.timing.cluster_total_ms = stats.cluster_total_ms_;
			result.timing.sphere_update_total_ms = stats.sphere_update_total_ms_;
			result.timing.error_total_ms = stats.error_total_ms_;
			result.timing.split_total_ms = stats.split_total_ms_;
			result.timing.average_iteration_ms = stats.average_iteration_ms_;
		}
		log_stage("optimization_end");

		log_stage("skeleton_begin");
		{
			ScopedBenchmarkTimer timer(result.timing.skeleton_construction_ms);
			core.build_skeleton_prepared(*points);
		}
		log_stage("skeleton_end");

		if (config.postprocess.topology_fix || config.postprocess.deg_face_deletion ||
			config.postprocess.nm_two_layer_prune || config.postprocess.residual_prune ||
			config.postprocess.face_post_processing)
		{
			log_stage("postprocess_begin");
			ScopedBenchmarkTimer timer(result.timing.postprocess_ms);
			if (config.postprocess.topology_fix)
				core.run_topology_fix_prepared(*points, config.postprocess.deg_face_deletion);
			else if (config.postprocess.deg_face_deletion)
				core.run_deg_face_deletion_prepared(*points);
			if (config.postprocess.nm_two_layer_prune)
				core.run_nm_two_layer_prune_prepared(*points);
			if (config.postprocess.residual_prune)
				core.run_completion_residual_prune_prepared(*points, config.postprocess.residual_prune_threshold);
			if (config.postprocess.face_post_processing)
				core.run_face_post_processing_prepared(*points);
			log_stage("postprocess_end");
		}
		else
		{
			result.timing.postprocess_ms = 0.0;
		}

		log_stage("export_begin");
		{
			ScopedBenchmarkTimer timer(result.timing.export_ms);
			core.export_skeleton_ply_prepared(*points, config.output.skeleton_ply);
			if (!std::filesystem::exists(config.output.skeleton_ply))
				throw std::runtime_error("Skeleton export did not create file: " + config.output.skeleton_ply);
		}
		log_stage("export_end");

		const auto counts = core.collect_counts(*points);
		result.counts.input_vertices = counts.input_vertices_;
		result.counts.input_points = counts.input_points_;
		result.counts.sample_points = counts.sample_points_;
		if (result.counts.sample_points_before_filtering == 0)
			result.counts.sample_points_before_filtering = counts.sample_points_;
		result.counts.final_spheres = counts.final_spheres_;
		result.counts.skeleton_vertices = counts.skeleton_vertices_;
		result.counts.skeleton_edges = counts.skeleton_edges_;
		result.counts.skeleton_faces = counts.skeleton_faces_;
		result.counts.optimization_iterations = counts.optimization_iterations_;
	}

	std::cout << "Input mode: " << result.input_mode << std::endl;
	std::cout << "Input path: " << result.input_path << std::endl;
	std::cout << "Apply filtering: " << (config.sampling.apply_filtering ? "true" : "false") << std::endl;
	std::cout << "Initial spheres: " << config.initialization.initial_nb_spheres << std::endl;
	std::cout << "Topology fix: " << (config.postprocess.topology_fix ? "true" : "false") << std::endl;
	std::cout << "Deg face deletion: " << (config.postprocess.deg_face_deletion ? "true" : "false") << std::endl;
	std::cout << "NM two-layer prune: " << (config.postprocess.nm_two_layer_prune ? "true" : "false") << std::endl;
	std::cout << "Residual prune: " << (config.postprocess.residual_prune ? "true" : "false")
			  << " threshold=" << config.postprocess.residual_prune_threshold << std::endl;
	std::cout << "Face post processing: " << (config.postprocess.face_post_processing ? "true" : "false") << std::endl;
	std::cout << "Samples: " << result.counts.sample_points_before_filtering << " -> " << result.counts.sample_points
			  << std::endl;
	std::cout << "Final spheres: " << result.counts.final_spheres << std::endl;
	const long long euler_characteristic = static_cast<long long>(result.counts.skeleton_vertices) -
		static_cast<long long>(result.counts.skeleton_edges) +
		static_cast<long long>(result.counts.skeleton_faces);
	std::cout << "Skeleton V/E/F: " << result.counts.skeleton_vertices << "/" << result.counts.skeleton_edges << "/"
			  << result.counts.skeleton_faces << "  Euler X=" << euler_characteristic << std::endl;
	std::cout << "Optimization iterations: " << result.counts.optimization_iterations << std::endl;
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
	benchmark_config.put("initialization.initial_nb_spheres", config.initialization.initial_nb_spheres);
	benchmark_config.put("initialization.init_min_cover_points", config.initialization.init_min_cover_points);
	benchmark_config.put("input.initial_ma_mode_override", to_string(config.input.initial_ma_mode_override));
	benchmark_config.put("auto_split.enabled", config.auto_split.enabled);
	benchmark_config.put("optimization.max_iterations_without_autosplit",
						 config.optimization.max_iterations_without_autosplit);
	benchmark_config.put("optimization.max_iterations_after_reaching_max_spheres",
						 config.optimization.max_iterations_after_reaching_max_spheres);
	benchmark_config.put("postprocess.topology_fix", config.postprocess.topology_fix);
	benchmark_config.put("postprocess.deg_face_deletion", config.postprocess.deg_face_deletion);
	benchmark_config.put("postprocess.nm_two_layer_prune", config.postprocess.nm_two_layer_prune);
	benchmark_config.put("postprocess.residual_prune", config.postprocess.residual_prune);
	benchmark_config.put("postprocess.residual_prune_threshold", config.postprocess.residual_prune_threshold);
	benchmark_config.put("postprocess.face_post_processing", config.postprocess.face_post_processing);
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
	timing.put("split_total", seconds_from_milliseconds(result.timing.split_total_ms));
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

} // namespace benchmark

} // namespace geometry

} // namespace cgogn
