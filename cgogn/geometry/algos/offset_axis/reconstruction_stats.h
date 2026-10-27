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
#ifndef CGOGN_GEOMETRY_ALGOS_OFFSET_AXIS_RECONSTRUCTION_STATS_H_
#define CGOGN_GEOMETRY_ALGOS_OFFSET_AXIS_RECONSTRUCTION_STATS_H_

#include <cgogn/core/utils/numerics.h>

namespace cgogn
{

namespace geometry
{

namespace offset_axis
{

struct ReconstructionTiming
{
	float64 preprocessing_ms = 0.0;
	float64 sampling_ms = 0.0;
	float64 sample_filtering_ms = 0.0;
	float64 kdtree_and_normals_ms = 0.0;
	float64 fitting_primitives_ms = 0.0;
	float64 initial_medial_axis_ms = 0.0;
	float64 sphere_initialization_ms = 0.0;
	float64 sphere_init_cover_ms = 0.0;
	float64 sphere_init_owner_ms = 0.0;
	float64 sphere_init_cluster_ms = 0.0;
	float64 optimization_ms = 0.0;
	float64 cluster_ms = 0.0;
	float64 sphere_update_ms = 0.0;
	float64 error_ms = 0.0;
	float64 skeleton_construction_ms = 0.0;
	float64 topology_processing_ms = 0.0;
};

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

// Silent repairs: places where the pipeline replaces an invalid value by a fallback instead of failing.
struct ReconstructionRepairs
{
	// sphere updates skipped because the solve gave non-finite values (sum over all optimizer iterations)
	uint32 sphere_update_nonfinite = 0;
};

} // namespace offset_axis

} // namespace geometry

} // namespace cgogn

#endif // CGOGN_GEOMETRY_ALGOS_OFFSET_AXIS_RECONSTRUCTION_STATS_H_
