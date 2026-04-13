#include <cgogn/core/utils/thread.h>
#include <cgogn/geometry/benchmark/benchmark_config.h>
#include <cgogn/geometry/benchmark/skeleton_benchmark_runner.h>

#include <chrono>
#include <exception>
#include <filesystem>
#include <iostream>
#include <string>
#include <vector>

namespace
{

double seconds_from_milliseconds(double value_ms)
{
	return value_ms / 1000.0;
}

} // namespace

int main(int argc, char** argv)
{
	using cgogn::geometry::benchmark::BenchmarkConfig;
	using cgogn::geometry::benchmark::BenchmarkResult;
	using cgogn::geometry::benchmark::SkeletonBenchmarkRunner;
	using cgogn::geometry::benchmark::expand_batch_benchmark_configs;
	using cgogn::geometry::benchmark::is_batch_case_completed;
	using cgogn::geometry::benchmark::is_batch_benchmark_config;
	using cgogn::geometry::benchmark::load_benchmark_config;
	using cgogn::geometry::benchmark::write_benchmark_timing_json;

	auto print_usage = [&](const char* exe_name) {
		std::cout << "Usage: " << exe_name << " <config.json>\n"
				  << "  <config.json>   Benchmark configuration file.\n";
	};

	auto ensure_output_directory = [](const std::string& output_path) {
		if (output_path.empty())
			return;
		const std::filesystem::path path(output_path);
		if (!path.parent_path().empty())
			std::filesystem::create_directories(path.parent_path());
	};

	auto print_result = [&](const BenchmarkConfig& config, const BenchmarkResult& result, bool include_config_loading) {
		if (include_config_loading)
			std::cout << "Config loading (s): " << seconds_from_milliseconds(result.timing.config_loading_ms) << std::endl;
		std::cout << "Input loading (s): " << seconds_from_milliseconds(result.timing.input_loading_ms) << std::endl;
		std::cout << "Preprocess (s): " << seconds_from_milliseconds(result.timing.preprocess_ms) << std::endl;
		std::cout << "Sampling (s): " << seconds_from_milliseconds(result.timing.sampling_ms) << std::endl;
		std::cout << "Optimization total (s): " << seconds_from_milliseconds(result.timing.optimization_total_ms) << std::endl;
		std::cout << "Skeleton construction (s): "
				  << seconds_from_milliseconds(result.timing.skeleton_construction_ms) << std::endl;
		std::cout << "Postprocess (s): " << seconds_from_milliseconds(result.timing.postprocess_ms) << std::endl;
		std::cout << "Export (s): " << seconds_from_milliseconds(result.timing.export_ms) << std::endl;
		std::cout << "Total (s): " << seconds_from_milliseconds(result.timing.total_ms) << std::endl;
		std::cout << "Skeleton saved to: " << config.output.skeleton_ply << std::endl;
		if (!config.output.timing_json.empty())
			std::cout << "Timing saved to: " << config.output.timing_json << std::endl;
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
			const double config_loading_ms =
				std::chrono::duration<double, std::milli>(config_end - config_start).count();

			if (!is_batch_benchmark_config(config))
			{
				ensure_output_directory(config.output.skeleton_ply);
				ensure_output_directory(config.output.timing_json);
				BenchmarkResult result = runner.run(config);
				result.timing.config_loading_ms = config_loading_ms;
				result.timing.total_ms += result.timing.config_loading_ms;
				write_benchmark_timing_json(config, result);
				print_result(config, result, true);
			}
			else
			{
				const auto expanded_case_configs = expand_batch_benchmark_configs(config);
				std::vector<BenchmarkConfig> case_configs;
				std::vector<std::string> skipped_completed_cases;
				case_configs.reserve(expanded_case_configs.size());
				skipped_completed_cases.reserve(expanded_case_configs.size());
				for (const auto& case_config : expanded_case_configs)
				{
					const std::string case_name = std::filesystem::path(case_config.input.input_path).stem().string();
					if (config.batch.resume_from_existing && is_batch_case_completed(case_config))
					{
						skipped_completed_cases.push_back(case_name);
						continue;
					}
					case_configs.push_back(case_config);
				}
				std::vector<std::string> failed_cases;
				std::cout << "Batch cases: total=" << expanded_case_configs.size() << " pending=" << case_configs.size();
				if (config.batch.resume_from_existing)
					std::cout << " skipped_completed=" << skipped_completed_cases.size();
				std::cout << std::endl;
				if (config.batch.resume_from_existing && config.benchmark.verbose)
				{
					for (const auto& skipped_case : skipped_completed_cases)
						std::cout << "Skip completed case: " << skipped_case << std::endl;
				}
				for (std::size_t i = 0; i < case_configs.size(); ++i)
				{
					const auto& case_config = case_configs[i];
					const std::string case_name = std::filesystem::path(case_config.input.input_path).stem().string();
					std::cout << "===== Case [" << (i + 1) << "/" << case_configs.size() << "] " << case_name
							  << " =====" << std::endl;
					try
					{
						ensure_output_directory(case_config.output.skeleton_ply);
						ensure_output_directory(case_config.output.timing_json);
						BenchmarkResult result = runner.run(case_config);
						write_benchmark_timing_json(case_config, result);
						print_result(case_config, result, false);
					}
					catch (const std::exception& e)
					{
						failed_cases.push_back(case_name);
						std::cerr << "Case failed: " << case_name << " -> " << e.what() << std::endl;
					}
				}

				std::cout << "Batch finished: success=" << (case_configs.size() - failed_cases.size())
						  << " failed=" << failed_cases.size();
				if (config.batch.resume_from_existing)
					std::cout << " skipped_completed=" << skipped_completed_cases.size();
				std::cout << std::endl;
				if (!failed_cases.empty())
				{
					std::string failed_summary;
					for (std::size_t i = 0; i < failed_cases.size(); ++i)
					{
						if (i > 0)
							failed_summary += ", ";
						failed_summary += failed_cases[i];
					}
					throw std::runtime_error("batch benchmark failed for cases: " + failed_summary);
				}
			}
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
