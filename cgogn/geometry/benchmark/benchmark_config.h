#ifndef CGOGN_GEOMETRY_BENCHMARK_CONFIG_H_
#define CGOGN_GEOMETRY_BENCHMARK_CONFIG_H_

#include <string>
#include <vector>

namespace cgogn
{

namespace geometry
{

namespace benchmark
{

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

enum class NeuralModelType
{
	UDF,
	MF
};

struct BenchmarkInputConfig
{
	InputMode mode = InputMode::PointCloud;
	InputGeometryType geometry_type = InputGeometryType::PointCloud;
	std::string input_path;
	std::string surface_path;
	std::string neural_udf_model_path;
	NeuralModelType neural_model_type = NeuralModelType::UDF;
	bool ma_flip_prune = true;
	float ma_flip_prune_alpha_factor = 1.0f;
};

struct BenchmarkSamplingConfig
{
	struct BridsonConfig
	{
		float sample_radius = 0.0025f;
	};

	float alpha = 0.005f;
	int knn_k = 10;
	bool apply_filtering = false;
	bool recompute_normals_after_sampling = false;
	int ray_sampler_batch_size = 4096;
	int batch_size = 1310640;
	float tol = 1e-5f;
	int udf_max_iterations = 3000;
	BridsonConfig bridson;
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
	bool topology_fix = false;
	bool deg_face_deletion = false;
	bool nm_two_layer_prune = false;
	bool residual_prune = false;
	bool face_post_processing = false;
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
	bool save_face_components = false;
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
