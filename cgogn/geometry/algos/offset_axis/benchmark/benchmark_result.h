#ifndef CGOGN_GEOMETRY_ALGOS_OFFSET_AXIS_BENCHMARK_RESULT_H_
#define CGOGN_GEOMETRY_ALGOS_OFFSET_AXIS_BENCHMARK_RESULT_H_

#include <cgogn/geometry/algos/offset_axis/reconstruction_stats.h>

#include <string>

namespace cgogn
{

namespace geometry
{

namespace benchmark
{

// timings measured by the benchmark around UDFReconstruction
struct BenchmarkExtraTiming
{
	double config_loading_ms = 0.0;
	double input_loading_ms = 0.0;
	double preprocess_ms = 0.0;
	double export_ms = 0.0;
	double total_ms = 0.0;
	double average_iteration_ms = 0.0;
};

// residual sheet prune statistics copied from SkeletonTopology::metrics() after reconstruction
struct BenchmarkResidualPruneStats
{
	uint32 residual_removed_sheets = 0;
	uint32 residual_removed_faces = 0;
	int32 residual_chi_before = 0;
	int32 residual_chi_after = 0;
	uint32 residual_fully_nm_removed_faces = 0;
	uint32 residual_small_sheets = 0;
	uint32 residual_small_sheets_guarded = 0;
	uint32 residual_small_sheets_with_boundary_edge = 0;
	uint32 residual_small_sheets_without_boundary_edge = 0;
};

struct BenchmarkResult
{
	std::string input_mode;
	std::string input_path;
	int seed = 0;
	ReconstructionCounts counts;
	ReconstructionTiming timing;
	BenchmarkExtraTiming extra_timing;
	BenchmarkResidualPruneStats residual_prune;
	ReconstructionRepairs repairs;
};

// single table shared by the console summary and the timing JSON export; order here is the JSON key order
struct BenchmarkTimingField
{
	const char* json_key;
	const char* print_label; // nullptr: not printed on the console
	double (*value_ms)(const BenchmarkResult& result);
	bool config_loading_only = false; // printed only when config loading is part of the run
};

inline constexpr BenchmarkTimingField kBenchmarkTimingFields[] = {
	{"config_loading", "Config loading", [](const BenchmarkResult& r) { return r.extra_timing.config_loading_ms; }, true},
	{"input_loading", "Input loading", [](const BenchmarkResult& r) { return r.extra_timing.input_loading_ms; }},
	{"preprocess", "Preprocess", [](const BenchmarkResult& r) { return r.extra_timing.preprocess_ms; }},
	{"sampling", "Sampling", [](const BenchmarkResult& r) { return r.timing.sampling_ms; }},
	{"sample_filtering", nullptr, [](const BenchmarkResult& r) { return r.timing.sample_filtering_ms; }},
	{"kdtree_bvh", nullptr, [](const BenchmarkResult& r) { return r.timing.kdtree_and_normals_ms; }},
	{"fitting_data", nullptr, [](const BenchmarkResult& r) { return r.timing.fitting_primitives_ms; }},
	{"initial_medial_axis", nullptr, [](const BenchmarkResult& r) { return r.timing.initial_medial_axis_ms; }},
	{"sphere_initialization", nullptr, [](const BenchmarkResult& r) { return r.timing.sphere_initialization_ms; }},
	{"optimization_total", "Optimization total", [](const BenchmarkResult& r) { return r.timing.optimization_ms; }},
	{"cluster_total", nullptr, [](const BenchmarkResult& r) { return r.timing.cluster_ms; }},
	{"sphere_update_total", nullptr, [](const BenchmarkResult& r) { return r.timing.sphere_update_ms; }},
	{"error_total", nullptr, [](const BenchmarkResult& r) { return r.timing.error_ms; }},
	{"skeleton_construction", "Skeleton construction",
	 [](const BenchmarkResult& r) { return r.timing.skeleton_construction_ms; }},
	{"postprocess", "Postprocess", [](const BenchmarkResult& r) { return r.timing.topology_processing_ms; }},
	{"export", "Export", [](const BenchmarkResult& r) { return r.extra_timing.export_ms; }},
	{"total", "Total", [](const BenchmarkResult& r) { return r.extra_timing.total_ms; }},
	{"average_iteration", nullptr, [](const BenchmarkResult& r) { return r.extra_timing.average_iteration_ms; }},
};

} // namespace benchmark

} // namespace geometry

} // namespace cgogn

#endif // CGOGN_GEOMETRY_ALGOS_OFFSET_AXIS_BENCHMARK_RESULT_H_
