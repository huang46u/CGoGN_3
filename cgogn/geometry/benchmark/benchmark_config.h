#ifndef CGOGN_GEOMETRY_BENCHMARK_CONFIG_H_
#define CGOGN_GEOMETRY_BENCHMARK_CONFIG_H_

#include <cgogn/geometry/algos/udf/neural_model_type.h>
#include <cgogn/geometry/algos/udf/reconstruction_options.h>

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
};

struct BenchmarkRuntimeConfig
{
	int num_measure_runs = 1;
	bool verbose = true;
	// runs the reconstruction without thread pool workers (every parallel loop on the calling thread); used to check
	// that the sequential and the parallel runs give the same skeleton
	bool sequential = false;
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
	cgogn::geometry::ReconstructionOptions reconstruction;
	std::string profile;
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

} // namespace benchmark

} // namespace geometry

} // namespace cgogn

#endif // CGOGN_GEOMETRY_BENCHMARK_CONFIG_H_
