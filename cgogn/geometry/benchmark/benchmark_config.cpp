#include <cgogn/geometry/benchmark/benchmark_config.h>

#include <boost/property_tree/json_parser.hpp>
#include <boost/property_tree/ptree.hpp>

#include <algorithm>
#include <filesystem>
#include <stdexcept>
#include <vector>

namespace cgogn
{

namespace geometry
{

namespace benchmark
{

namespace
{

using boost::property_tree::ptree;
namespace fs = std::filesystem;

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

bool has_matching_extension(const fs::path& path, const std::string& extension)
{
	if (extension.empty())
		return true;
	return path.extension().string() == extension;
}

std::vector<fs::path> collect_input_files(const fs::path& directory, const std::string& extension, bool recursive)
{
	std::vector<fs::path> files;
	if (!fs::exists(directory))
		fail_config("batch input directory does not exist: " + directory.string());
	if (!fs::is_directory(directory))
		fail_config("batch input directory is not a directory: " + directory.string());

	if (recursive)
	{
		for (const auto& entry : fs::recursive_directory_iterator(directory))
		{
			if (entry.is_regular_file() && has_matching_extension(entry.path(), extension))
				files.push_back(entry.path().lexically_normal());
		}
	}
	else
	{
		for (const auto& entry : fs::directory_iterator(directory))
		{
			if (entry.is_regular_file() && has_matching_extension(entry.path(), extension))
				files.push_back(entry.path().lexically_normal());
		}
	}

	std::sort(files.begin(), files.end(), [](const fs::path& lhs, const fs::path& rhs) {
		return lhs.filename().string() < rhs.filename().string();
	});
	return files;
}

fs::path find_surface_file_for_case(const BenchmarkBatchConfig& batch, const std::string& stem)
{
	if (batch.surface_directory.empty())
		return {};

	const fs::path surface_dir(batch.surface_directory);
	if (!fs::exists(surface_dir))
		fail_config("batch surface directory does not exist: " + surface_dir.string());

	if (!batch.surface_extension.empty())
	{
		const fs::path candidate = surface_dir / (stem + batch.surface_extension);
		if (!fs::exists(candidate))
			fail_config("missing batch surface file for case `" + stem + "`: " + candidate.string());
		return candidate.lexically_normal();
	}

	static const char* kDefaultSurfaceExtensions[] = {".obj", ".ply", ".off", ".stl"};
	for (const char* ext : kDefaultSurfaceExtensions)
	{
		const fs::path candidate = surface_dir / (stem + ext);
		if (fs::exists(candidate))
			return candidate.lexically_normal();
	}

	fail_config("missing batch surface file for case `" + stem + "` in directory: " + surface_dir.string());
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
	if (value == "line_quadric_distance_fix_r" || value == "line_quadric_distance")
		return DistanceMode::LineQuadricDistance;
	if (value == "line_quadric_distance_free_r" || value == "line_quadric_distance_free_radius")
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
	return value == DistanceMode::LineQuadricDistanceFreeRadius ? "line_quadric_distance_free_r"
																 : "line_quadric_distance_fix_r";
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
	const fs::path config_path = fs::absolute(path).lexically_normal();
	if (!fs::exists(config_path))
		fail_config("config file does not exist: " + config_path.string());

	ptree root;
	boost::property_tree::read_json(config_path.string(), root);
	const fs::path config_directory = config_path.parent_path();

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

	if (auto batch = root.get_child_optional("batch"))
	{
		config.batch.enabled = get_value<bool>(*batch, "enabled", false);
		config.batch.input_directory =
			resolve_config_relative_path(config_directory, get_value<std::string>(*batch, "input_directory", ""));
		config.batch.surface_directory =
			resolve_config_relative_path(config_directory, get_value<std::string>(*batch, "surface_directory", ""));
		config.batch.neural_udf_model_directory = resolve_config_relative_path(
			config_directory, get_value<std::string>(*batch, "neural_udf_model_directory", ""));
		config.batch.input_extension = get_value<std::string>(*batch, "input_extension", ".ply");
		config.batch.surface_extension = get_value<std::string>(*batch, "surface_extension", "");
		config.batch.neural_udf_model_extension =
			get_value<std::string>(*batch, "neural_udf_model_extension", ".pt");
		config.batch.recursive = get_value<bool>(*batch, "recursive", false);
		config.batch.output_directory =
			resolve_config_relative_path(config_directory, get_value<std::string>(*batch, "output_directory", ""));
		config.batch.timing_directory =
			resolve_config_relative_path(config_directory, get_value<std::string>(*batch, "timing_directory", ""));
	}

	const ptree& output = root.get_child("output");
	config.output.skeleton_ply =
		resolve_config_relative_path(config_directory, get_value<std::string>(output, "skeleton_ply", ""));
	config.output.timing_json =
		resolve_config_relative_path(config_directory, get_value<std::string>(output, "timing_json", ""));

	if (config.benchmark.num_measure_runs != 1)
		fail_config("`benchmark.num_measure_runs` is reserved in v1 and must be `1`");
	if (config.input.mode == InputMode::NeuralUDF && !config.batch.enabled)
	{
		if (config.input.neural_udf_model_path.empty())
			fail_config("`input.neural_udf_model_path` is required when `input.mode` is `neural_udf`");
		if (!fs::exists(config.input.neural_udf_model_path))
			fail_config("neural model does not exist: " + config.input.neural_udf_model_path);
	}
	if (config.batch.enabled)
	{
		if (config.batch.input_directory.empty())
			fail_config("`batch.input_directory` is required when `batch.enabled` is true");
		if (config.batch.output_directory.empty())
			fail_config("`batch.output_directory` is required when `batch.enabled` is true");
		if (config.input.mode == InputMode::NeuralUDF && config.batch.neural_udf_model_directory.empty())
			fail_config("`batch.neural_udf_model_directory` is required when `input.mode` is `neural_udf`");
	}
	else if (config.input.input_path.empty() && config.input.surface_path.empty())
	{
		fail_config("one of `input.input_path` or `input.surface_path` must be provided");
	}
	if (!config.batch.enabled && config.output.skeleton_ply.empty())
		fail_config("`output.skeleton_ply` is required in single-run mode");
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

bool is_batch_benchmark_config(const BenchmarkConfig& config)
{
	return config.batch.enabled;
}

std::vector<BenchmarkConfig> expand_batch_benchmark_configs(const BenchmarkConfig& base_config)
{
	if (!base_config.batch.enabled)
		return {base_config};

	const fs::path input_directory(base_config.batch.input_directory);
	const auto input_files =
		collect_input_files(input_directory, base_config.batch.input_extension, base_config.batch.recursive);
	if (input_files.empty())
		fail_config("no batch input files found in: " + input_directory.string());

	std::vector<BenchmarkConfig> expanded_configs;
	expanded_configs.reserve(input_files.size());

	for (const fs::path& input_file : input_files)
	{
		BenchmarkConfig config = base_config;
		const std::string stem = input_file.stem().string();
		config.input.input_path = input_file.string();
		config.input.surface_path.clear();
		config.input.neural_udf_model_path.clear();

		if (!config.batch.surface_directory.empty())
			config.input.surface_path = find_surface_file_for_case(config.batch, stem).string();

		if (config.input.mode == InputMode::NeuralUDF)
		{
			const fs::path model_path =
				fs::path(config.batch.neural_udf_model_directory) / (stem + config.batch.neural_udf_model_extension);
			if (!fs::exists(model_path))
				fail_config("missing batch neural model for case `" + stem + "`: " + model_path.string());
			config.input.neural_udf_model_path = model_path.lexically_normal().string();
		}

		config.output.skeleton_ply =
			(fs::path(config.batch.output_directory) / (stem + ".ply")).lexically_normal().string();
		if (!config.batch.timing_directory.empty())
		{
			config.output.timing_json =
				(fs::path(config.batch.timing_directory) / (stem + "_timing.json")).lexically_normal().string();
		}
		else
		{
			config.output.timing_json.clear();
		}

		expanded_configs.push_back(std::move(config));
	}

	return expanded_configs;
}

} // namespace benchmark

} // namespace geometry

} // namespace cgogn
