/*******************************************************************************
 * CGoGN: Combinatorial and Geometric modeling with Generic N-dimensional Maps  *
 * Copyright (C), IGG Group, ICube, University of Strasbourg, France            *
 *******************************************************************************/
#ifndef CGOGN_GEOMETRY_ALGOS_UDF_RECONSTRUCTION_H_
#define CGOGN_GEOMETRY_ALGOS_UDF_RECONSTRUCTION_H_

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
#include <cgogn/geometry/algos/udf/reconstruction_stages.h>
#include <cgogn/geometry/algos/udf/reconstruction_state.h>
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

template <typename SURFACE, typename POINTS, typename NONMANIFOLD, typename RaySamplerTag>
class UDFReconstruction
{
public:
	using State = ReconstructionState<SURFACE, POINTS, NONMANIFOLD, RaySamplerTag>;

	using PVertex = typename State::PVertex;
	using SVertex = typename State::SVertex;
	using SFace = typename State::SFace;
	using NMVertex = typename State::NMVertex;
	using NMEdge = typename State::NMEdge;
	using NMFace = typename State::NMFace;
	using SpheresOptimizerType = typename State::SpheresOptimizerType;
	using SpheresOptimizerMetrics = typename State::SpheresOptimizerMetrics;
	using TopologyType = typename State::TopologyType;
	using Field = typename State::Field;
	using RayConfig = typename State::RayConfig;
	using RayTraits = typename State::RayTraits;
	using RaySampler = typename State::RaySampler;
	using RayParams = typename State::RayParams;
	using PAttribute = typename State::PAttribute;
	template <typename T>
	using PointAttribute = typename State::template PointAttribute<T>;
	template <typename T>
	using SurfaceAttribute = typename State::template SurfaceAttribute<T>;
	template <typename T>
	using SkeletonAttribute = typename State::template SkeletonAttribute<T>;

	using NeuralModelType = typename State::NeuralModelType;
	using Status = typename State::Status;
	using Data = typename State::Data;
	using Options = typename State::Options;
	using Timing = typename State::Timing;
	using Counts = typename State::Counts;
	using Result = typename State::Result;
	using IterationCallback = typename State::IterationCallback;

	explicit UDFReconstruction(Data data, Options options = {})
	{
		s_.data = data;
		s_.options = options;
	}
	UDFReconstruction(const UDFReconstruction&) = delete;
	UDFReconstruction& operator=(const UDFReconstruction&) = delete;

	State& state()
	{
		return s_;
	}
	const State& state() const
	{
		return s_;
	}

	Data& data()
	{
		s_.input.invalidate();
		return s_.data;
	}
	const Data& data() const
	{
		return s_.data;
	}
	Options& options()
	{
		return s_.options;
	}
	const Options& options() const
	{
		return s_.options;
	}

	void reset()
	{
		udf_stages::reset(s_);
	}
	Status initialize_data()
	{
		return udf_stages::initialize_data(s_);
	}
	Status load_neural_model(const std::filesystem::path& path, NeuralModelType type)
	{
		return udf_stages::load_neural_model(s_, path, type);
	}
	Status sample_alpha_level_set()
	{
		return udf_stages::sample_alpha_level_set(s_);
	}
	Status apply_sampling_filtering()
	{
		return udf_stages::apply_sampling_filtering(s_);
	}
	Status build_kdtree_and_normals()
	{
		return udf_stages::build_kdtree_and_normals(s_);
	}
	Status compute_fitting_primitives()
	{
		return udf_stages::compute_fitting_primitives(s_);
	}
	Status compute_initial_medial_axis()
	{
		return udf_stages::compute_initial_medial_axis(s_);
	}
	Status initialize_spheres()
	{
		return udf_stages::initialize_spheres(s_);
	}
	Status optimize_spheres(const IterationCallback& iteration_callback = {})
	{
		return udf_stages::optimize_spheres(s_, iteration_callback);
	}
	Status build_skeleton()
	{
		return udf_stages::build_skeleton(s_);
	}
	Status fix_topology()
	{
		return udf_stages::fix_topology(s_);
	}
	Status prune_residual_sheets()
	{
		return udf_stages::prune_residual_sheets(s_);
	}
	Status sample()
	{
		return udf_stages::sample(s_);
	}
	Status compute_fitting_data()
	{
		return udf_stages::compute_fitting_data(s_);
	}
	typename SpheresOptimizerType::Status update_spheres_once(Scalar sqem_update_lambda_line_plane)
	{
		return udf_stages::update_spheres_once(s_, sqem_update_lambda_line_plane);
	}
	Status fix_topology_and_prune()
	{
		return udf_stages::fix_topology_and_prune(s_);
	}
	Result run(const IterationCallback& iteration_callback = {})
	{
		return udf_stages::reconstruct(s_, iteration_callback);
	}

	Counts counts() const
	{
		return udf_stages::counts(s_);
	}
	ReconstructionRepairs repairs() const
	{
		return udf_stages::repairs(s_);
	}

	NONMANIFOLD* skeleton() const
	{
		return s_.data.skeleton;
	}
	const std::shared_ptr<SkeletonAttribute<Vec3>>& skeleton_position_attribute() const
	{
		return s_.skeleton_position;
	}
	const std::shared_ptr<SkeletonAttribute<Scalar>>& skeleton_radius_attribute() const
	{
		return s_.skeleton_radius;
	}
	SpheresOptimizerType& spheres_optimizer()
	{
		return udf_stages::optimizer_for(s_);
	}
	TopologyType& skeleton_topology()
	{
		return udf_stages::topology_for(s_);
	}

private:
	State s_;
};

} // namespace cgogn::geometry

#endif // CGOGN_GEOMETRY_ALGOS_UDF_RECONSTRUCTION_H_
