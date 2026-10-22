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
#ifndef CGOGN_GEOMETRY_ALGOS_UDF_RECONSTRUCTION_STATE_H_
#define CGOGN_GEOMETRY_ALGOS_UDF_RECONSTRUCTION_STATE_H_

#include <cgogn/core/functions/attributes.h>
#include <cgogn/core/functions/mesh_info.h>
#include <cgogn/core/functions/traversals/edge.h>
#include <cgogn/core/types/cell_marker.h>
#include <cgogn/core/utils/thread_pool.h>
#include <cgogn/geometry/algos/area.h>
#include <cgogn/geometry/algos/fitting.h>
#include <cgogn/geometry/algos/medial_axis.h>
#include <cgogn/geometry/algos/normal.h>
#include <cgogn/geometry/algos/udf/alpha_projection.h>
#include <cgogn/geometry/algos/udf/alpha_sampling.h>
#include <cgogn/geometry/algos/udf/distance_field.h>
#include <cgogn/geometry/algos/udf/input_geometry.h>
#include <cgogn/geometry/algos/udf/neural_model_type.h>
#include <cgogn/geometry/algos/udf/normal_utils.h>
#include <cgogn/geometry/algos/udf/parallel_for.h>
#include <cgogn/geometry/algos/udf/point_cloud_distance_field.h>
#include <cgogn/geometry/algos/udf/reconstruction_options.h>
#include <cgogn/geometry/algos/udf/reconstruction_stats.h>
#include <cgogn/geometry/algos/udf/sample_processing.h>
#include <cgogn/geometry/algos/udf/skeleton_topology.h>
#include <cgogn/geometry/algos/udf/spatial_query.h>
#include <cgogn/geometry/algos/udf/spheres_optimizer.h>
#include <cgogn/geometry/algos/udf/surface_distance_field.h>
#include <cgogn/geometry/functions/bounding_box.h>
#include <cgogn/geometry/functions/normal.h>
#include <cgogn/geometry/types/line_quadric.h>
#include <cgogn/geometry/types/ray_level_set_sampler.h>
#include <cgogn/geometry/types/ray_level_set_sampler_traits.h>
#include <cgogn/geometry/types/spatial_grid.h>
#include <cgogn/geometry/types/spherical_quadric.h>
#include <cgogn/geometry/types/vector_traits.h>

#include <libacc/bvh_tree.h>
#include <libacc/kd_tree.h>

#include <Eigen/Dense>

#include <algorithm>
#include <atomic>
#include <chrono>
#include <cmath>
#include <filesystem>
#include <functional>
#include <future>
#include <memory>
#include <set>
#include <stdexcept>
#include <string>
#include <tuple>
#include <type_traits>
#include <utility>
#include <vector>

namespace cgogn::geometry
{

// All the state of the UDF reconstruction: the types and the data members of the former UDFReconstruction.
// The members are public and are accessed by the UDFReconstruction facade as s_.name.
// Their declaration order is the initialization order: field refers to input and options, so it is declared
// after them.
template <typename SURFACE, typename POINTS, typename NONMANIFOLD, typename RaySamplerTag>
struct ReconstructionState
{
	using PVertex = typename mesh_traits<POINTS>::Vertex;
	using SVertex = typename mesh_traits<SURFACE>::Vertex;
	using SFace = typename mesh_traits<SURFACE>::Face;
	using NMVertex = typename mesh_traits<NONMANIFOLD>::Vertex;
	using NMEdge = typename mesh_traits<NONMANIFOLD>::Edge;
	using NMFace = typename mesh_traits<NONMANIFOLD>::Face;
	using SpheresOptimizerType = SpheresOptimizer<POINTS>;
	using SpheresOptimizerMetrics = typename SpheresOptimizerType::Metrics;
	using TopologyType = SkeletonTopology<POINTS, NONMANIFOLD>;
	using Field = typename DistanceFieldFor<RaySamplerTag, SURFACE, POINTS>::type;
	using RayConfig = typename Field::RayConfig;
	using RayTraits = typename RayConfig::SamplerTraits;
	using RaySampler = typename Field::RaySampler;
	using RayParams = typename RaySampler::Params;
	using PAttribute = typename mesh_traits<POINTS>::template Attribute<Vec3>;
	template <typename T>
	using PointAttribute = typename mesh_traits<POINTS>::template Attribute<T>;
	template <typename T>
	using SurfaceAttribute = typename mesh_traits<SURFACE>::template Attribute<T>;
	template <typename T>
	using SkeletonAttribute = typename mesh_traits<NONMANIFOLD>::template Attribute<T>;

	using NeuralModelType = cgogn::geometry::NeuralModelType;
	enum class Status
	{
		success,
		invalid_data,
		model_load_failed,
		sampling_failed,
		fitting_failed,
		medial_axis_failed,
		sphere_initialization_failed,
		optimization_failed,
		skeleton_build_failed,
		topology_fix_failed
	};

	struct Data
	{
		SURFACE* surface = nullptr;
		POINTS* input_points = nullptr;
		POINTS* samples = nullptr;
		POINTS* spheres = nullptr;
		NONMANIFOLD* skeleton = nullptr;
	};

	using Options = ReconstructionOptions;

	using Timing = ReconstructionTiming;
	using Counts = ReconstructionCounts;

	struct Result
	{
		Status status = Status::invalid_data;
		Timing timing;
		Counts counts;
	};

	using IterationCallback = std::function<void(const SpheresOptimizerMetrics&)>;

	Data data;
	Options options;
	std::unique_ptr<RaySampler> ray_sampler;
	InputGeometry<SURFACE, POINTS> input;
	// distance field of the input kind (CPU fields or the neural model); it refers to input and options, so it is
	// declared after them
	Field field{input, options};
	std::unique_ptr<acc::KDTree<3, uint32>> sample_kdtree;
	std::unique_ptr<acc::KDTree<3, uint32>> sample_ma_kdtree;
	std::vector<PVertex> sample_vertices;
	std::vector<PVertex> sample_ma_vertices;
	std::unique_ptr<SpheresOptimizerType> optimizer;
	std::unique_ptr<TopologyType> topology;
	bool fitting_data_computed = false;
	bool initialized = false;
	uint32 samples_before_filtering = 0;

	std::shared_ptr<PointAttribute<Vec3>> input_position, input_normal;
	std::shared_ptr<PointAttribute<std::vector<PVertex>>> input_knn;
	std::shared_ptr<PointAttribute<Vec3>> sample_position, sample_normal;
	std::shared_ptr<PointAttribute<Scalar>> sample_area;
	std::shared_ptr<PointAttribute<std::vector<PVertex>>> sample_knn;
	std::shared_ptr<PointAttribute<Spherical_Quadric>> sample_quadric;
	std::shared_ptr<PointAttribute<Line_Quadric>> sample_line_quadric;
	std::shared_ptr<PointAttribute<Vec3>> sample_ma_position;
	std::shared_ptr<PointAttribute<Scalar>> sample_ma_radius;
	std::shared_ptr<PointAttribute<PVertex>> sample_ma_secondary, sample_sphere;
	std::shared_ptr<PointAttribute<Scalar>> sample_error;
	std::shared_ptr<PointAttribute<Vec3>> sphere_position;
	std::shared_ptr<PointAttribute<Scalar>> sphere_radius;
	std::shared_ptr<PointAttribute<std::vector<PVertex>>> sphere_cluster;
	std::shared_ptr<PointAttribute<Scalar>> sphere_cluster_area;
	std::shared_ptr<PointAttribute<Vec4>> sphere_cluster_color;
	std::shared_ptr<PointAttribute<std::vector<PVertex>>> sphere_neighbors;
	std::shared_ptr<PointAttribute<Scalar>> sphere_error, sphere_error_raw;
	std::shared_ptr<PointAttribute<NMVertex>> sphere_skeleton_vertex;
	std::shared_ptr<SkeletonAttribute<Vec3>> skeleton_position;
	std::shared_ptr<SkeletonAttribute<Scalar>> skeleton_radius;
	std::shared_ptr<SkeletonAttribute<PVertex>> skeleton_source_sphere;
	std::shared_ptr<SkeletonAttribute<std::set<std::size_t>>> skeleton_incident_tets;
	std::shared_ptr<SkeletonAttribute<uint32>> skeleton_edge_degree;
};

} // namespace cgogn::geometry

#endif // CGOGN_GEOMETRY_ALGOS_UDF_RECONSTRUCTION_STATE_H_
