/*******************************************************************************
 * CGoGN: Combinatorial and Geometric modeling with Generic N-dimensional Maps  *
 * Copyright (C), IGG Group, ICube, University of Strasbourg, France            *
 *                                                                              *
 * This library is free software; you can redistribute it and/or modify it      *
 * under the terms of the GNU Lesser General Public License as published by the *
 * Free Software Foundation; either version 2.1 of the License, or (at your     *
 * option) any later version.                                                   *
 *                                                                              *
 * This library is distributed in the hope that it will be useful, but WITHOUT *
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or        *
 * FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License *
 * for more details.                                                            *
 *                                                                              *
 * You should have received a copy of the GNU Lesser General Public License     *
 * along with this library; if not, write to the Free Software Foundation,      *
 * Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301 USA.           *
 *                                                                              *
 * Web site: http://cgogn.unistra.fr/                                           *
 * Contact information: cgogn@unistra.fr                                        *
 *                                                                              *
 *******************************************************************************/
#ifndef CGOGN_GEOMETRY_ALGOS_OFFSET_AXIS_BENCHMARK_RESULT_H_
#define CGOGN_GEOMETRY_ALGOS_OFFSET_AXIS_BENCHMARK_RESULT_H_

#include <cgogn/core/utils/numerics.h>

#include <string>

namespace cgogn
{

namespace geometry
{

namespace offset_axis
{

namespace benchmark
{

// time of each reconstruction step
struct ReconstructionTiming
{
	float64 preprocessing_ms = 0.0;
	float64 sampling_ms = 0.0;
	float64 sample_filtering_ms = 0.0;
	float64 kdtree_and_normals_ms = 0.0;
	float64 initial_medial_axis_ms = 0.0;
	float64 fitting_primitives_ms = 0.0;
	float64 sphere_initialization_ms = 0.0;
	float64 optimization_ms = 0.0;
	float64 skeleton_construction_ms = 0.0;
	float64 topology_processing_ms = 0.0;
};

// sizes of the meshes after the reconstruction
struct ReconstructionCounts
{
	uint32 input_vertices = 0;
	uint32 input_points = 0;
	uint32 samples_before_filtering = 0;
	uint32 samples = 0;
	uint32 spheres = 0;
	uint32 skeleton_vertices = 0;
	uint32 skeleton_edges = 0;
	uint32 skeleton_faces = 0;
	uint32 optimization_iterations = 0;
};

// timings of the benchmark outside the reconstruction steps
struct BenchmarkExtraTiming
{
	double config_loading_ms = 0.0;
	double input_loading_ms = 0.0;
	double preprocess_ms = 0.0;
	double export_ms = 0.0;
	double total_ms = 0.0;
	double average_iteration_ms = 0.0;
};

struct BenchmarkResult
{
	std::string input_mode;
	std::string input_path;
	int seed = 0;
	ReconstructionCounts counts;
	ReconstructionTiming timing;
	BenchmarkExtraTiming extra_timing;
};

// single table shared by the console summary and the timing JSON export; order here is the JSON key order
struct BenchmarkTimingField
{
	const char* json_key;
	const char* print_label; // nullptr: not printed on the console
	double (*value_ms)(const BenchmarkResult& result);
	bool config_loading_only = false; // printed only when config loading is part of the run
};

inline constexpr BenchmarkTimingField benchmark_timing_fields[] = {
	{"config_loading", "Config loading", [](const BenchmarkResult& r) { return r.extra_timing.config_loading_ms; }, true},
	{"input_loading", "Input loading", [](const BenchmarkResult& r) { return r.extra_timing.input_loading_ms; }},
	{"preprocess", "Preprocess", [](const BenchmarkResult& r) { return r.extra_timing.preprocess_ms; }},
	{"sampling", "Sampling", [](const BenchmarkResult& r) { return r.timing.sampling_ms; }},
	{"sample_filtering", nullptr, [](const BenchmarkResult& r) { return r.timing.sample_filtering_ms; }},
	{"kdtree_bvh", nullptr, [](const BenchmarkResult& r) { return r.timing.kdtree_and_normals_ms; }},
	{"initial_medial_axis", nullptr, [](const BenchmarkResult& r) { return r.timing.initial_medial_axis_ms; }},
	{"fitting_data", nullptr, [](const BenchmarkResult& r) { return r.timing.fitting_primitives_ms; }},
	{"sphere_initialization", nullptr, [](const BenchmarkResult& r) { return r.timing.sphere_initialization_ms; }},
	{"optimization_total", "Optimization total", [](const BenchmarkResult& r) { return r.timing.optimization_ms; }},
	{"skeleton_construction", "Skeleton construction",
	 [](const BenchmarkResult& r) { return r.timing.skeleton_construction_ms; }},
	{"postprocess", "Postprocess", [](const BenchmarkResult& r) { return r.timing.topology_processing_ms; }},
	{"export", "Export", [](const BenchmarkResult& r) { return r.extra_timing.export_ms; }},
	{"total", "Total", [](const BenchmarkResult& r) { return r.extra_timing.total_ms; }},
	{"average_iteration", nullptr, [](const BenchmarkResult& r) { return r.extra_timing.average_iteration_ms; }},
};

} // namespace benchmark

} // namespace offset_axis

} // namespace geometry

} // namespace cgogn

#endif // CGOGN_GEOMETRY_ALGOS_OFFSET_AXIS_BENCHMARK_RESULT_H_
