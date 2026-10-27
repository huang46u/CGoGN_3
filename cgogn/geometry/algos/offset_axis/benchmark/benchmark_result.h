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

// wall-clock times of the benchmark steps
struct BenchmarkTiming
{
	double loading_ms = 0.0;        // config (single runs only), input, neural model and attributes
	double sampling_ms = 0.0;       // level set sampling, sample filtering and PCA normals
	double initialization_ms = 0.0; // initial medial axis with its flip prune, fitting data and initial spheres
	double optimization_ms = 0.0;
	double postprocess_ms = 0.0; // skeleton construction, topology fix and residual prune
	double export_ms = 0.0;
	double total_ms = 0.0;
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

struct BenchmarkResult
{
	std::string input_mode;
	std::string input_path;
	ReconstructionCounts counts;
	BenchmarkTiming timing;
};

// single table shared by the console summary and the timing JSON export; order here is the JSON key order
struct BenchmarkTimingField
{
	const char* json_key;
	const char* print_label;
	double BenchmarkTiming::*value_ms;
};

inline constexpr BenchmarkTimingField benchmark_timing_fields[] = {
	{"loading", "Loading", &BenchmarkTiming::loading_ms},
	{"sampling", "Sampling", &BenchmarkTiming::sampling_ms},
	{"initialization", "Initialization", &BenchmarkTiming::initialization_ms},
	{"optimization", "Optimization", &BenchmarkTiming::optimization_ms},
	{"postprocess", "Postprocess", &BenchmarkTiming::postprocess_ms},
	{"export", "Export", &BenchmarkTiming::export_ms},
	{"total", "Total", &BenchmarkTiming::total_ms},
};

} // namespace benchmark

} // namespace offset_axis

} // namespace geometry

} // namespace cgogn

#endif // CGOGN_GEOMETRY_ALGOS_OFFSET_AXIS_BENCHMARK_RESULT_H_
