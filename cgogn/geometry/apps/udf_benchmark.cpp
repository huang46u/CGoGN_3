#include <cgogn/core/utils/thread.h>
#include <cgogn/geometry/benchmark/benchmark_config.h>
#include <cgogn/geometry/benchmark/skeleton_benchmark_runner.h>

#include <chrono>
#include <exception>
#include <iostream>

namespace
{

double seconds_from_milliseconds(double value_ms)
{
	return value_ms / 1000.0;
}

} // namespace

int main(int argc, char** argv)
{
	using cgogn::geometry::benchmark::BenchmarkResult;
	using cgogn::geometry::benchmark::SkeletonBenchmarkRunner;
	using cgogn::geometry::benchmark::load_benchmark_config;
	using cgogn::geometry::benchmark::write_benchmark_timing_json;

	auto print_usage = [&](const char* exe_name) {
		std::cout << "Usage: " << exe_name << " <config.json>\n"
				  << "  <config.json>   Benchmark configuration file.\n";
	};

	if (argc != 2)
	{
		print_usage(argv[0]);
		return argc == 1 ? 1 : 0;
	}

	const std::string arg = argv[1];
	if (arg == "-h" || arg == "--help")
	{
		print_usage(argv[0]);
		return 0;
	}

	try
	{
		cgogn::thread_start(0);
		try
		{
			const auto config_start = std::chrono::high_resolution_clock::now();
			auto config = load_benchmark_config(arg);
			const auto config_end = std::chrono::high_resolution_clock::now();

			SkeletonBenchmarkRunner runner;
			BenchmarkResult result = runner.run(config);
			result.timing.config_loading_ms =
				std::chrono::duration<double, std::milli>(config_end - config_start).count();
			result.timing.total_ms += result.timing.config_loading_ms;

			write_benchmark_timing_json(config, result);

			std::cout << "Config loading (s): " << seconds_from_milliseconds(result.timing.config_loading_ms) << std::endl;
			std::cout << "Input loading (s): " << seconds_from_milliseconds(result.timing.input_loading_ms) << std::endl;
			std::cout << "Preprocess (s): " << seconds_from_milliseconds(result.timing.preprocess_ms) << std::endl;
			std::cout << "Sampling (s): " << seconds_from_milliseconds(result.timing.sampling_ms) << std::endl;
			std::cout << "Optimization total (s): " << seconds_from_milliseconds(result.timing.optimization_total_ms) << std::endl;
			std::cout << "Skeleton construction (s): " << seconds_from_milliseconds(result.timing.skeleton_construction_ms) << std::endl;
			std::cout << "Postprocess (s): " << seconds_from_milliseconds(result.timing.postprocess_ms) << std::endl;
			std::cout << "Export (s): " << seconds_from_milliseconds(result.timing.export_ms) << std::endl;
			std::cout << "Total (s): " << seconds_from_milliseconds(result.timing.total_ms) << std::endl;
			std::cout << "Skeleton saved to: " << config.output.skeleton_ply << std::endl;
			if (!config.output.timing_json.empty())
				std::cout << "Timing saved to: " << config.output.timing_json << std::endl;
			cgogn::thread_stop();
			return 0;
		}
		catch (...)
		{
			cgogn::thread_stop();
			throw;
		}
	}
	catch (const std::exception& e)
	{
		std::cerr << "udf_benchmark failed: " << e.what() << std::endl;
		return 1;
	}
}
