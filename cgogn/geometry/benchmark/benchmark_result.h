#ifndef CGOGN_GEOMETRY_BENCHMARK_RESULT_H_
#define CGOGN_GEOMETRY_BENCHMARK_RESULT_H_

#include <string>

namespace cgogn
{

namespace geometry
{

namespace benchmark
{

struct BenchmarkCounts
{
	unsigned int input_vertices = 0;
	unsigned int input_points = 0;
	unsigned int sample_points_before_filtering = 0;
	unsigned int sample_points = 0;
	unsigned int final_spheres = 0;
	unsigned int skeleton_vertices = 0;
	unsigned int skeleton_edges = 0;
	unsigned int skeleton_faces = 0;
	unsigned int optimization_iterations = 0;
};

struct BenchmarkTiming
{
	double config_loading_ms = 0.0;
	double input_loading_ms = 0.0;
	double preprocess_ms = 0.0;
	double sampling_ms = 0.0;
	double sample_filtering_ms = 0.0;
	double kdtree_bvh_ms = 0.0;
	double fitting_data_ms = 0.0;
	double initial_medial_axis_ms = 0.0;
	double sphere_initialization_ms = 0.0;
	double optimization_total_ms = 0.0;
	double cluster_total_ms = 0.0;
	double sphere_update_total_ms = 0.0;
	double error_total_ms = 0.0;
	double split_total_ms = 0.0;
	double skeleton_construction_ms = 0.0;
	double postprocess_ms = 0.0;
	double export_ms = 0.0;
	double total_ms = 0.0;
	double average_iteration_ms = 0.0;
};

struct BenchmarkResult
{
	std::string input_mode;
	std::string input_path;
	int seed = 0;
	BenchmarkCounts counts;
	BenchmarkTiming timing;
};

} // namespace benchmark

} // namespace geometry

} // namespace cgogn

#endif // CGOGN_GEOMETRY_BENCHMARK_RESULT_H_
