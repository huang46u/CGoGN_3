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

enum class NeuralModelType
{
	UDF,
	MF
};

enum class DistanceMode
{
	LineQuadricDistance,
	LineQuadricDistanceFreeRadius
};

enum class SphereCorrectionMode
{
	Always,
	OnSplit
};

enum class AutoSplitMode
{
	ErrorThreshold,
	MaxNbSpheres
};

enum class InitialMAMode
{
	Auto,
	Displacement,
	ShrinkingBall
};

enum class BridsonCandidateMode
{
	Shell3D,
	Plane2D
};

struct BenchmarkInputConfig
{
	InputMode mode = InputMode::PointCloud;
	std::string input_path;
	std::string surface_path;
	std::string neural_udf_model_path;
	NeuralModelType neural_model_type = NeuralModelType::UDF;
	InitialMAMode initial_ma_mode_override = InitialMAMode::Auto;
};

struct BenchmarkSamplingConfig
{
	struct BridsonConfig
	{
		float sample_radius = 0.0025f;
		BridsonCandidateMode candidate_mode = BridsonCandidateMode::Shell3D;
		float outer_radius_scale = 2.0f;
		int warmup_iterations = 2;
		int sample_iterations = 30;
		int parent_batch_size = 8192;
		int max_samples = 4000000;
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
	DistanceMode distance_mode = DistanceMode::LineQuadricDistance;
	bool use_local_clusters = false;
	unsigned int local_cluster_connectivity_refresh_interval = 10;
	unsigned int max_iterations_without_autosplit = 300;
	unsigned int max_iterations_after_reaching_max_spheres = 100;
	float sqem_update_lambda_full = 0.05f;
	float sqem_update_lambda_line_plane = 0.20f;
	float sqem_fix_radius_scale = 1.0f;
	bool udf_center_enabled = false;
	float udf_center_lambda = 0.10f;
	bool sphere_correction = false;
	SphereCorrectionMode sphere_correction_mode = SphereCorrectionMode::Always;
	bool lock_skeleton_connectivity = false;
	bool auto_stop = false;
};

struct BenchmarkAutoSplitConfig
{
	bool enabled = false;
	AutoSplitMode mode = AutoSplitMode::ErrorThreshold;
	float error_threshold = 0.00025f;
	unsigned int max_nb_spheres = 500;
	float ratio = 0.2f;
	unsigned int max_per_iter_error = 10;
	unsigned int max_per_iter_max = 100;
};

struct BenchmarkInitializationConfig
{
	int seed = 42;
	float init_dilation_constant = 0.001f;
	unsigned int init_min_cover_points = 10;
	unsigned int initial_nb_spheres = 1;
};

struct BenchmarkPostprocessConfig
{
	bool topology_fix = false;
	bool deg_face_deletion = false;
	bool nm_two_layer_prune = false;
	bool residual_prune = false;
	float residual_prune_threshold = -1.0f;
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
	BenchmarkAutoSplitConfig auto_split;
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
std::string to_string(NeuralModelType value);
std::string to_string(DistanceMode value);
std::string to_string(SphereCorrectionMode value);
std::string to_string(AutoSplitMode value);
std::string to_string(InitialMAMode value);
std::string to_string(BridsonCandidateMode value);

} // namespace benchmark

} // namespace geometry

} // namespace cgogn

#endif // CGOGN_GEOMETRY_BENCHMARK_CONFIG_H_
