#include <cgogn/geometry/benchmark/benchmark_config.h>

#include <boost/property_tree/json_parser.hpp>
#include <boost/property_tree/ptree.hpp>

#include <filesystem>
#include <stdexcept>

namespace cgogn
{

namespace geometry
{

namespace benchmark
{

namespace
{

using boost::property_tree::ptree;

[[noreturn]] void fail_config(const std::string& message)
{
	throw std::runtime_error("Invalid benchmark config: " + message);
}

template <typename T>
T require_value(const ptree& tree, const std::string& key)
{
	auto value = tree.get_optional<T>(key);
	if (!value)
		fail_config("missing required field `" + key + "`");
	return *value;
}

template <typename T>
T get_value(const ptree& tree, const std::string& key, const T& default_value)
{
	return tree.get<T>(key, default_value);
}

std::string resolve_config_relative_path(const std::filesystem::path& config_directory, const std::string& raw_path)
{
	if (raw_path.empty())
		return raw_path;

	const std::filesystem::path path(raw_path);
	if (path.is_absolute())
		return path.lexically_normal().string();

	return std::filesystem::absolute(config_directory / path).lexically_normal().string();
}

InputMode parse_input_mode(const std::string& value)
{
	if (value == "point_cloud")
		return InputMode::PointCloud;
	if (value == "surface_mesh")
		return InputMode::SurfaceMesh;
	if (value == "neural_udf")
		return InputMode::NeuralUDF;
	fail_config("unsupported `input.mode`: " + value);
}

NeuralModelType parse_neural_model_type(const std::string& value)
{
	if (value == "udf")
		return NeuralModelType::UDF;
	if (value == "mf")
		return NeuralModelType::MF;
	fail_config("unsupported `input.neural_model_type`: " + value);
}

DistanceMode parse_distance_mode(const std::string& value)
{
	if (value == "line_quadric_distance")
		return DistanceMode::LineQuadricDistance;
	if (value == "line_quadric_distance_free_radius")
		return DistanceMode::LineQuadricDistanceFreeRadius;
	fail_config("unsupported `optimization.distance_mode`: " + value);
}

SphereCorrectionMode parse_correction_mode(const std::string& value)
{
	if (value == "always")
		return SphereCorrectionMode::Always;
	if (value == "on_split")
		return SphereCorrectionMode::OnSplit;
	fail_config("unsupported `optimization.sphere_correction_mode`: " + value);
}

AutoSplitMode parse_auto_split_mode(const std::string& value)
{
	if (value == "error_threshold")
		return AutoSplitMode::ErrorThreshold;
	if (value == "max_nb_spheres")
		return AutoSplitMode::MaxNbSpheres;
	fail_config("unsupported `auto_split.mode`: " + value);
}

InitialMAMode parse_initial_ma_mode(const std::string& value)
{
	if (value == "auto")
		return InitialMAMode::Auto;
	if (value == "displacement")
		return InitialMAMode::Displacement;
	if (value == "shrinking_ball")
		return InitialMAMode::ShrinkingBall;
	fail_config("unsupported `input.initial_ma_mode_override`: " + value);
}

} // namespace

std::string to_string(InputMode value)
{
	switch (value)
	{
	case InputMode::PointCloud:
		return "point_cloud";
	case InputMode::SurfaceMesh:
		return "surface_mesh";
	case InputMode::NeuralUDF:
	default:
		return "neural_udf";
	}
}

std::string to_string(NeuralModelType value)
{
	return value == NeuralModelType::MF ? "mf" : "udf";
}

std::string to_string(DistanceMode value)
{
	return value == DistanceMode::LineQuadricDistanceFreeRadius ? "line_quadric_distance_free_radius"
																 : "line_quadric_distance";
}

std::string to_string(SphereCorrectionMode value)
{
	return value == SphereCorrectionMode::OnSplit ? "on_split" : "always";
}

std::string to_string(AutoSplitMode value)
{
	return value == AutoSplitMode::MaxNbSpheres ? "max_nb_spheres" : "error_threshold";
}

std::string to_string(InitialMAMode value)
{
	switch (value)
	{
	case InitialMAMode::Displacement:
		return "displacement";
	case InitialMAMode::ShrinkingBall:
		return "shrinking_ball";
	case InitialMAMode::Auto:
	default:
		return "auto";
	}
}

BenchmarkConfig load_benchmark_config(const std::string& path)
{
	const std::filesystem::path config_path = std::filesystem::absolute(path).lexically_normal();
	if (!std::filesystem::exists(config_path))
		fail_config("config file does not exist: " + config_path.string());

	ptree root;
	boost::property_tree::read_json(config_path.string(), root);
	const std::filesystem::path config_directory = config_path.parent_path();

	BenchmarkConfig config;
	config.config_path = config_path.string();

	const ptree& input = root.get_child("input");
	config.input.mode = parse_input_mode(require_value<std::string>(input, "mode"));
	config.input.input_path = resolve_config_relative_path(config_directory, get_value<std::string>(input, "input_path", ""));
	config.input.surface_path = resolve_config_relative_path(config_directory, get_value<std::string>(input, "surface_path", ""));
	config.input.neural_udf_model_path =
		resolve_config_relative_path(config_directory, get_value<std::string>(input, "neural_udf_model_path", ""));
	config.input.initial_ma_mode_override =
		parse_initial_ma_mode(get_value<std::string>(input, "initial_ma_mode_override", "auto"));
	if (auto neural_type = input.get_optional<std::string>("neural_model_type"))
		config.input.neural_model_type = parse_neural_model_type(*neural_type);

	const ptree& sampling = root.get_child("sampling");
	config.sampling.alpha = require_value<float>(sampling, "alpha");
	config.sampling.grid_cell_size = require_value<float>(sampling, "grid_cell_size");
	config.sampling.knn_k = require_value<int>(sampling, "knn_k");
	config.sampling.num_alpha_samples = require_value<int>(sampling, "num_alpha_samples");
	config.sampling.apply_filtering = get_value<bool>(sampling, "apply_filtering", false);
	config.sampling.ray_sampler_batch_size = require_value<int>(sampling, "ray_sampler_batch_size");
	config.sampling.batch_size = require_value<int>(sampling, "batch_size");
	config.sampling.tol = require_value<float>(sampling, "tol");
	config.sampling.udf_max_iterations = require_value<int>(sampling, "udf_max_iterations");

	const ptree& optimization = root.get_child("optimization");
	config.optimization.distance_mode =
		parse_distance_mode(require_value<std::string>(optimization, "distance_mode"));
	config.optimization.use_local_clusters = get_value<bool>(optimization, "use_local_clusters", false);
	config.optimization.local_cluster_connectivity_refresh_interval =
		get_value<unsigned int>(optimization, "local_cluster_connectivity_refresh_interval", 10u);
	config.optimization.max_iterations_without_autosplit =
		get_value<unsigned int>(optimization, "max_iterations_without_autosplit", 300u);
	config.optimization.max_iterations_after_reaching_max_spheres =
		get_value<unsigned int>(optimization, "max_iterations_after_reaching_max_spheres", 100u);
	config.optimization.sqem_update_lambda_full = require_value<float>(optimization, "sqem_update_lambda_full");
	config.optimization.sqem_update_lambda_line_plane =
		require_value<float>(optimization, "sqem_update_lambda_line_plane");
	config.optimization.sqem_fix_radius_scale = require_value<float>(optimization, "sqem_fix_radius_scale");
	config.optimization.udf_center_enabled = get_value<bool>(optimization, "udf_center_enabled", false);
	config.optimization.udf_center_lambda = get_value<float>(optimization, "udf_center_lambda", 0.10f);
	config.optimization.sphere_correction = get_value<bool>(optimization, "sphere_correction", false);
	config.optimization.sphere_correction_mode = parse_correction_mode(
		get_value<std::string>(optimization, "sphere_correction_mode", "always"));
	config.optimization.lock_skeleton_connectivity =
		get_value<bool>(optimization, "lock_skeleton_connectivity", false);
	config.optimization.auto_stop = get_value<bool>(optimization, "auto_stop", false);

	const ptree& auto_split = root.get_child("auto_split");
	config.auto_split.enabled = get_value<bool>(auto_split, "enabled", false);
	if (auto mode = auto_split.get_optional<std::string>("mode"))
		config.auto_split.mode = parse_auto_split_mode(*mode);
	config.auto_split.error_threshold = get_value<float>(auto_split, "error_threshold", 0.00025f);
	config.auto_split.max_nb_spheres = get_value<unsigned int>(auto_split, "max_nb_spheres", 500u);
	config.auto_split.ratio = get_value<float>(auto_split, "ratio", 0.2f);
	config.auto_split.max_per_iter_error = get_value<unsigned int>(auto_split, "max_per_iter_error", 10u);
	config.auto_split.max_per_iter_max = get_value<unsigned int>(auto_split, "max_per_iter_max", 100u);

	const ptree& initialization = root.get_child("initialization");
	config.initialization.seed = require_value<int>(initialization, "seed");
	config.initialization.init_dilation_constant = require_value<float>(initialization, "init_dilation_constant");
	config.initialization.init_min_cover_points = require_value<unsigned int>(initialization, "init_min_cover_points");
	config.initialization.initial_nb_spheres =
		get_value<unsigned int>(initialization, "initial_nb_spheres", 1u);

	if (auto postprocess = root.get_child_optional("postprocess"))
	{
		config.postprocess.topology_fix = get_value<bool>(*postprocess, "topology_fix", false);
		config.postprocess.deg_face_deletion = get_value<bool>(*postprocess, "deg_face_deletion", false);
		config.postprocess.face_post_processing = get_value<bool>(*postprocess, "face_post_processing", false);
	}

	const ptree& benchmark_tree = root.get_child("benchmark");
	config.benchmark.num_measure_runs = get_value<int>(benchmark_tree, "num_measure_runs", 1);
	config.benchmark.verbose = get_value<bool>(benchmark_tree, "verbose", true);

	const ptree& output = root.get_child("output");
	config.output.skeleton_ply = resolve_config_relative_path(
		config_directory, require_value<std::string>(output, "skeleton_ply"));
	config.output.timing_json =
		resolve_config_relative_path(config_directory, get_value<std::string>(output, "timing_json", ""));

	if (config.benchmark.num_measure_runs != 1)
		fail_config("`benchmark.num_measure_runs` is reserved in v1 and must be `1`");
	if (config.input.mode == InputMode::NeuralUDF)
	{
		if (config.input.neural_udf_model_path.empty())
			fail_config("`input.neural_udf_model_path` is required when `input.mode` is `neural_udf`");
		if (!std::filesystem::exists(config.input.neural_udf_model_path))
			fail_config("neural model does not exist: " + config.input.neural_udf_model_path);
	}
	if (config.input.input_path.empty() && config.input.surface_path.empty())
		fail_config("one of `input.input_path` or `input.surface_path` must be provided");
	if (config.initialization.initial_nb_spheres == 0)
		fail_config("`initialization.initial_nb_spheres` must be > 0");
	if (config.optimization.max_iterations_without_autosplit == 0)
		fail_config("`optimization.max_iterations_without_autosplit` must be > 0");
	if (config.optimization.max_iterations_after_reaching_max_spheres == 0)
		fail_config("`optimization.max_iterations_after_reaching_max_spheres` must be > 0");
	if (config.auto_split.enabled)
	{
		if (config.auto_split.mode == AutoSplitMode::ErrorThreshold)
		{
			if (config.auto_split.error_threshold <= 0.0f)
				fail_config("`auto_split.error_threshold` must be > 0 when `auto_split.mode=error_threshold`");
		}
		else if (config.auto_split.max_nb_spheres == 0)
		{
			fail_config("`auto_split.max_nb_spheres` must be > 0 when `auto_split.mode=max_nb_spheres`");
		}
	}

	return config;
}

} // namespace benchmark

} // namespace geometry

} // namespace cgogn
