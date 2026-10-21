#ifndef CGOGN_GEOMETRY_ALGOS_UDF_RECONSTRUCTION_STATS_H_
#define CGOGN_GEOMETRY_ALGOS_UDF_RECONSTRUCTION_STATS_H_

#include <cgogn/core/utils/numerics.h>

namespace cgogn::geometry
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

} // namespace cgogn::geometry

#endif // CGOGN_GEOMETRY_ALGOS_UDF_RECONSTRUCTION_STATS_H_
