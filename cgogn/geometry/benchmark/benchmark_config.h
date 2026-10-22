#ifndef CGOGN_GEOMETRY_BENCHMARK_CONFIG_H_
#define CGOGN_GEOMETRY_BENCHMARK_CONFIG_H_

#include <cgogn/geometry/algos/udf/neural_model_type.h>

#include <optional>
#include <string>
#include <vector>

namespace cgogn
{

namespace geometry
{

namespace benchmark
{

using cgogn::geometry::NeuralModelType;

enum class InputMode
{
	PointCloud,
	SurfaceMesh,
	NeuralUDF
};

enum class InputGeometryType
{
	PointCloud,
	Mesh
};

struct BenchmarkInputConfig
{
	InputMode mode = InputMode::PointCloud;
	InputGeometryType geometry_type = InputGeometryType::PointCloud;
	std::string input_path;
	std::string surface_path;
	std::string neural_udf_model_path;
	NeuralModelType neural_model_type = NeuralModelType::udf;
	bool ma_flip_prune = true;
	float ma_flip_prune_alpha_factor = 1.0f;
};

struct BenchmarkSamplingConfig
{
	float alpha = 0.005f;
	int knn_k = 10;
	float sample_radius = 0.0025f;
	int candidates_per_parent = 12;
	int parent_batch_size = 1024;
	// unset: default of the neural model type (default_apply_filtering / default_recompute_normals_after_sampling)
	std::optional<bool> apply_filtering;
	std::optional<bool> recompute_normals_after_sampling;
};

struct BenchmarkOptimizationConfig
{
	float sqem_update_lambda_line_plane = 0.20f;
};

struct BenchmarkInitializationConfig
{
	int seed = 42;
	float init_dilation_constant = 0.001f;
};

struct BenchmarkPostprocessConfig
{
	bool residual_prune = false;
	std::string residual_prune_mode = "collapsible_only";
};

struct BenchmarkRuntimeConfig
{
	int num_measure_runs = 1;
	bool verbose = true;
};

struct BenchmarkOutputConfig
{
	std::string skeleton_ply;
	std::string timing_json;
};

struct BenchmarkBatchConfig
{
	bool enabled = false;
	std::string input_directory;
	std::string surface_directory;
	std::string neural_udf_model_directory;
	std::string input_extension = ".ply";
	std::string surface_extension;
	std::string neural_udf_model_extension = ".pt";
	bool recursive = false;
	bool resume_from_existing = false;
	std::string output_directory;
	std::string timing_directory;
};

struct BenchmarkConfig
{
	std::string config_path;
	BenchmarkInputConfig input;
	BenchmarkSamplingConfig sampling;
	BenchmarkOptimizationConfig optimization;
	BenchmarkInitializationConfig initialization;
	BenchmarkPostprocessConfig postprocess;
	BenchmarkRuntimeConfig benchmark;
	BenchmarkOutputConfig output;
	BenchmarkBatchConfig batch;
};

BenchmarkConfig load_benchmark_config(const std::string& path);
std::vector<BenchmarkConfig> expand_batch_benchmark_configs(const BenchmarkConfig& base_config);
bool is_batch_benchmark_config(const BenchmarkConfig& config);
bool is_batch_case_completed(const BenchmarkConfig& config);

std::string to_string(InputMode value);
std::string to_string(InputGeometryType value);
std::string to_string(NeuralModelType value);

} // namespace benchmark

} // namespace geometry

} // namespace cgogn

#endif // CGOGN_GEOMETRY_BENCHMARK_CONFIG_H_
