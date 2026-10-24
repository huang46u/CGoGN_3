#ifndef CGOGN_GEOMETRY_ALGOS_OFFSET_AXIS_BENCHMARK_RUNNER_H_
#define CGOGN_GEOMETRY_ALGOS_OFFSET_AXIS_BENCHMARK_RUNNER_H_

#include <cgogn/geometry/algos/offset_axis/benchmark/benchmark_config.h>
#include <cgogn/geometry/algos/offset_axis/benchmark/benchmark_result.h>
#include <cgogn/core/types/maps/cmap/cmap0.h>
#include <cgogn/core/types/maps/cmap/cmap2.h>

namespace cgogn
{

namespace geometry
{

namespace offset_axis
{

namespace benchmark
{

class SkeletonBenchmarkRunner
{
public:
	BenchmarkResult run(const BenchmarkConfig& config) const;
};

// Imports the input geometry of the config (surface or point cloud), rescaled exactly as for a benchmark run.
// Returns true when the surface is used (data.surface must point to it), false for the point cloud. Throws on error.
bool load_benchmark_input(const BenchmarkConfig& config, CMap2& surface, CMap0& input_points);

void write_benchmark_timing_json(const BenchmarkConfig& config, const BenchmarkResult& result);

} // namespace benchmark

} // namespace offset_axis

} // namespace geometry

} // namespace cgogn

#endif // CGOGN_GEOMETRY_ALGOS_OFFSET_AXIS_BENCHMARK_RUNNER_H_

