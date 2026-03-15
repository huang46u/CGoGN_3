#ifndef CGOGN_GEOMETRY_BENCHMARK_SKELETON_BENCHMARK_RUNNER_H_
#define CGOGN_GEOMETRY_BENCHMARK_SKELETON_BENCHMARK_RUNNER_H_

#include <cgogn/geometry/benchmark/benchmark_config.h>
#include <cgogn/geometry/benchmark/benchmark_result.h>

namespace cgogn
{

namespace geometry
{

namespace benchmark
{

class SkeletonBenchmarkRunner
{
public:
	BenchmarkResult run(const BenchmarkConfig& config) const;
};

void write_benchmark_timing_json(const BenchmarkConfig& config, const BenchmarkResult& result);

} // namespace benchmark

} // namespace geometry

} // namespace cgogn

#endif // CGOGN_GEOMETRY_BENCHMARK_SKELETON_BENCHMARK_RUNNER_H_

