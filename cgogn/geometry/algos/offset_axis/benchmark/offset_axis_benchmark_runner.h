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

