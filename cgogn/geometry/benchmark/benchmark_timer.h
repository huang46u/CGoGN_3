#ifndef CGOGN_GEOMETRY_BENCHMARK_TIMER_H_
#define CGOGN_GEOMETRY_BENCHMARK_TIMER_H_

#include <chrono>

namespace cgogn
{

namespace geometry
{

namespace benchmark
{

class ScopedBenchmarkTimer
{
public:
	explicit ScopedBenchmarkTimer(double& output_ms)
		: output_ms_(output_ms), start_(std::chrono::high_resolution_clock::now())
	{
	}

	~ScopedBenchmarkTimer()
	{
		const auto end = std::chrono::high_resolution_clock::now();
		output_ms_ += std::chrono::duration<double, std::milli>(end - start_).count();
	}

private:
	double& output_ms_;
	std::chrono::time_point<std::chrono::high_resolution_clock> start_;
};

} // namespace benchmark

} // namespace geometry

} // namespace cgogn

#endif // CGOGN_GEOMETRY_BENCHMARK_TIMER_H_
