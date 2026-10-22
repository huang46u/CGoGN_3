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
#include <cgogn/geometry/algos/udf/reconstruction_options.h>
#include <cgogn/geometry/algos/udf/reconstruction_state.h>
#include <cgogn/geometry/algos/udf/reconstruction_stats.h>
#include <cgogn/geometry/algos/udf/sample_processing.h>
#include <cgogn/geometry/algos/udf/skeleton_topology.h>
#include <cgogn/geometry/algos/udf/spatial_query.h>
#include <cgogn/geometry/algos/udf/surface_distance_field.h>
#include <cgogn/geometry/algos/udf/neural_model_type.h>
#include <cgogn/geometry/algos/udf/normal_utils.h>
#include <cgogn/geometry/algos/udf/parallel_for.h>
#include <cgogn/geometry/algos/udf/point_cloud_distance_field.h>
#include <cgogn/geometry/algos/udf/spheres_optimizer.h>
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

// Fixed constants of the reconstruction (not user options). Lengths are absolute, in model units (the inputs are
// rescaled to [0,1] before the reconstruction).
namespace udf_constants
{
// Padding added on each side of the input bounding box. The padded box is the sampling domain and the box the
// projected points are clamped to.
inline constexpr Scalar sampling_padding = Scalar(0.1);
// Relative enlargement of the padded sampling box from which the seed rays are cast (its half size is multiplied by
// 1 + ray_box_relative_expand). Overrides the ray sampler default (0.1).
inline constexpr Scalar ray_box_relative_expand = Scalar(0.05);
// Samples farther than alpha + filter_margin from the input are removed by the sampling filter.
inline constexpr Scalar filter_margin = Scalar(1e-2);
// MAFlipPrune, neural UDF: a medial centre is retried when its UDF value exceeds this factor times alpha.
inline constexpr Scalar udf_flip_threshold_factor = Scalar(1.2);
// MAFlipPrune, MF: range the ma_flip_prune_alpha_factor option is clamped to.
inline constexpr Scalar mf_flip_alpha_factor_min = Scalar(1);
inline constexpr Scalar mf_flip_alpha_factor_max = Scalar(5);
} // namespace udf_constants

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
		s_.input.invalidate();
		s_.sample_kdtree.reset();
		s_.sample_ma_kdtree.reset();
		s_.sample_vertices.clear();
		s_.sample_ma_vertices.clear();
		s_.ray_sampler.reset();
		s_.optimizer.reset();
		s_.topology.reset();
		s_.field.reset();
		s_.fitting_data_computed = false;
		s_.samples_before_filtering = 0;
		s_.initialized = false;
	}

	Status initialize_data()
	{
		if (s_.initialized)
			return Status::success;
		if (!s_.data.input_points || !s_.data.samples || !s_.data.spheres || !s_.data.skeleton)
			return Status::invalid_data;
		s_.input_position = get_or_add_attribute<Vec3, PVertex>(*s_.data.input_points, "position");
		s_.input_normal = get_or_add_attribute<Vec3, PVertex>(*s_.data.input_points, "normal");
		s_.input_knn = get_or_add_attribute<std::vector<PVertex>, PVertex>(*s_.data.input_points, "knn");
		s_.sample_position = get_or_add_attribute<Vec3, PVertex>(*s_.data.samples, "position");
		s_.sample_normal = get_or_add_attribute<Vec3, PVertex>(*s_.data.samples, "normal");
		s_.sample_area = get_or_add_attribute<Scalar, PVertex>(*s_.data.samples, "area");
		s_.sample_knn = get_or_add_attribute<std::vector<PVertex>, PVertex>(*s_.data.samples, "knn");
		s_.sample_quadric = get_or_add_attribute<Spherical_Quadric, PVertex>(*s_.data.samples, "quadric");
		s_.sample_line_quadric = get_or_add_attribute<Line_Quadric, PVertex>(*s_.data.samples, "line_quadric");
		s_.sample_ma_position = get_or_add_attribute<Vec3, PVertex>(*s_.data.samples, "ma_position");
		s_.sample_ma_radius = get_or_add_attribute<Scalar, PVertex>(*s_.data.samples, "ma_radius");
		s_.sample_ma_secondary = get_or_add_attribute<PVertex, PVertex>(*s_.data.samples, "ma_secondary_vertex");
		s_.sample_sphere = get_or_add_attribute<PVertex, PVertex>(*s_.data.samples, "sphere");
		s_.sample_error = get_or_add_attribute<Scalar, PVertex>(*s_.data.samples, "error");
		s_.sphere_position = get_or_add_attribute<Vec3, PVertex>(*s_.data.spheres, "position");
		s_.sphere_radius = get_or_add_attribute<Scalar, PVertex>(*s_.data.spheres, "radius");
		s_.sphere_cluster = get_or_add_attribute<std::vector<PVertex>, PVertex>(*s_.data.spheres, "cluster");
		s_.sphere_cluster_area = get_or_add_attribute<Scalar, PVertex>(*s_.data.spheres, "cluster_area");
		s_.sphere_cluster_color = get_or_add_attribute<Vec4, PVertex>(*s_.data.spheres, "cluster_color");
		s_.sphere_neighbors = get_or_add_attribute<std::vector<PVertex>, PVertex>(*s_.data.spheres, "neighbor_clusters");
		s_.sphere_error = get_or_add_attribute<Scalar, PVertex>(*s_.data.spheres, "error");
		s_.sphere_error_raw = get_or_add_attribute<Scalar, PVertex>(*s_.data.spheres, "error_not_normalized");
		s_.sphere_skeleton_vertex = get_or_add_attribute<NMVertex, PVertex>(*s_.data.spheres, "skeleton_vertex");
		s_.skeleton_position = get_or_add_attribute<Vec3, NMVertex>(*s_.data.skeleton, "position");
		s_.skeleton_radius = get_or_add_attribute<Scalar, NMVertex>(*s_.data.skeleton, "radius");
		s_.skeleton_source_sphere = get_or_add_attribute<PVertex, NMVertex>(*s_.data.skeleton, "source_sphere");
		s_.skeleton_incident_tets = get_or_add_attribute<std::set<std::size_t>, NMFace>(*s_.data.skeleton, "incident_tets");
		s_.skeleton_edge_degree = get_or_add_attribute<uint32, NMEdge>(*s_.data.skeleton, "degree");
		s_.initialized = s_.input_position && s_.sample_position && s_.sample_normal && s_.sphere_position && s_.skeleton_position;
		s_.input.bind(s_.data.surface, s_.data.input_points, s_.input_position, s_.input_normal);
		return s_.initialized ? Status::success : Status::invalid_data;
	}

	Status load_neural_model(const std::filesystem::path& path, NeuralModelType type)
	{
		if (!s_.field.load_model(path, type))
			return Status::model_load_failed;
		s_.options.neural_model_type = type;
		return Status::success;
	}

	Status sample_alpha_level_set();
	Status apply_sampling_filtering();
	Status build_kdtree_and_normals();
	Status compute_fitting_primitives();
	Status compute_initial_medial_axis();
	Status initialize_spheres();
	Status optimize_spheres(const IterationCallback& iteration_callback = {});
	Status build_skeleton();
	Status fix_topology();
	Status prune_residual_sheets();
	Result run(const IterationCallback& iteration_callback = {});

	Counts counts() const
	{
		Counts c;
		if (s_.data.input_points)
			c.input_points = nb_cells<PVertex>(*s_.data.input_points);
		if (s_.data.surface)
			c.input_vertices = nb_cells<SVertex>(*s_.data.surface);
		else
			c.input_vertices = c.input_points;
		if (s_.data.samples)
			c.samples = nb_cells<PVertex>(*s_.data.samples);
		c.samples_before_filtering = s_.samples_before_filtering;
		if (s_.data.spheres)
			c.spheres = nb_cells<PVertex>(*s_.data.spheres);
		if (s_.data.skeleton)
		{
			c.skeleton_vertices = nb_cells<NMVertex>(*s_.data.skeleton);
			c.skeleton_edges = nb_cells<NMEdge>(*s_.data.skeleton);
			c.skeleton_faces = nb_cells<NMFace>(*s_.data.skeleton);
		}
		if (s_.optimizer)
			c.optimization_iterations = s_.optimizer->metrics().iteration;
		return c;
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
		return optimizer_for();
	}
	TopologyType& skeleton_topology()
	{
		return topology_for();
	}

private:
	State s_;

	void build_sample_kdtree();
	void recompute_samples_normals_from_input();
	bool recompute_sample_normal_pca(PVertex vertex);
	SpheresOptimizerType& optimizer_for();
	TopologyType& topology_for();
	bool evaluate_topology_scores(const std::vector<Vec3>& points, std::vector<Scalar>& values);
	static Status status_of(FieldPrepareResult result)
	{
		switch (result)
		{
		case FieldPrepareResult::ok:
			return Status::success;
		case FieldPrepareResult::invalid_data:
			return Status::invalid_data;
		case FieldPrepareResult::model_not_loaded:
			return Status::model_load_failed;
		}
		return Status::invalid_data;
	}
};

template <typename SURFACE, typename POINTS, typename NONMANIFOLD, typename RaySamplerTag>
void UDFReconstruction<SURFACE, POINTS, NONMANIFOLD, RaySamplerTag>::build_sample_kdtree()
{
	s_.sample_kdtree.reset();
	s_.sample_ma_kdtree.reset();
	s_.sample_vertices.clear();
	s_.sample_ma_vertices.clear();
	if (!s_.data.samples || !s_.sample_position)
		return;
	std::vector<Vec3> positions;
	std::vector<Vec3> ma_positions;
	foreach_cell(*s_.data.samples, [&](PVertex v) {
		const uint32 idx = index_of(*s_.data.samples, v);
		positions.push_back((*s_.sample_position)[idx]);
		s_.sample_vertices.push_back(v);
		const Vec3 ma = (*s_.sample_ma_position)[idx];
		if (ma.allFinite() && (*s_.sample_ma_radius)[idx] > Scalar(0))
		{
			ma_positions.push_back(ma);
			s_.sample_ma_vertices.push_back(v);
		}
		return true;
	});
	if (!positions.empty())
		s_.sample_kdtree = std::make_unique<acc::KDTree<3, uint32>>(positions);
	if (!ma_positions.empty())
		s_.sample_ma_kdtree = std::make_unique<acc::KDTree<3, uint32>>(ma_positions);
}

template <typename SURFACE, typename POINTS, typename NONMANIFOLD, typename RaySamplerTag>
void UDFReconstruction<SURFACE, POINTS, NONMANIFOLD, RaySamplerTag>::recompute_samples_normals_from_input()
{
	s_.input.bind(s_.data.surface, s_.data.input_points, s_.input_position, s_.input_normal);
	if (!s_.input.input_kdtree() && s_.data.surface)
	{
		s_.input.build_surface_bvh();
		if (s_.input.surface_bvh())
		{
			foreach_cell(*s_.data.samples, [&](PVertex v) {
				const uint32 idx = index_of(*s_.data.samples, v);
				const Vec3& position = (*s_.sample_position)[idx];
				std::pair<uint32, Vec3> closest;
				Vec3 normal(0, 0, 1);
				if (s_.input.surface_bvh()->closest_point(position, &closest))
				{
					normal = unit_normal_or_z(position - closest.second);
				}
				(*s_.sample_normal)[idx] = normal;
				return true;
			});
			return;
		}
	}
	if (s_.input.input_kdtree() && s_.data.input_points)
	{
		foreach_cell(*s_.data.samples, [&](PVertex v) {
			const uint32 idx = index_of(*s_.data.samples, v);
			const Vec3& position = (*s_.sample_position)[idx];
			std::pair<uint32, Scalar> nearest;
			Vec3 normal(0, 0, 1);
			if (s_.input.input_kdtree()->find_nn(position, &nearest) && nearest.first < s_.input.input_vertices().size())
			{
				const uint32 input_idx = index_of(*s_.data.input_points, s_.input.input_vertices()[nearest.first]);
				if (input_idx != INVALID_INDEX)
				{
					normal = position - (*s_.input_position)[input_idx];
					if (normal.squaredNorm() < Scalar(1e-12))
						normal = (*s_.input_normal)[input_idx];
					normal = unit_normal_or_z(normal);
				}
			}
			(*s_.sample_normal)[idx] = normal;
			return true;
		});
		return;
	}
	foreach_cell(*s_.data.samples, [&](PVertex v) {
		(*s_.sample_normal)[index_of(*s_.data.samples, v)] = Vec3(0, 0, 1);
		return true;
	});
}

template <typename SURFACE, typename POINTS, typename NONMANIFOLD, typename RaySamplerTag>
bool UDFReconstruction<SURFACE, POINTS, NONMANIFOLD, RaySamplerTag>::recompute_sample_normal_pca(PVertex vertex)
{
	if (!s_.data.samples || !s_.sample_kdtree || !vertex.is_valid())
		return false;
	const int k = std::max(3, s_.options.knn_k);
	const uint32 vertex_index = index_of(*s_.data.samples, vertex);
	if (vertex_index == INVALID_INDEX)
		return false;
	const Vec3& center = (*s_.sample_position)[vertex_index];
	std::vector<std::pair<uint32, Scalar>> neighbors_index;
	s_.sample_kdtree->find_nns(center, k + 1, &neighbors_index);
	std::vector<Vec3> neighbors;
	neighbors.reserve(k);
	for (const auto& item : neighbors_index)
	{
		const PVertex neighbor = s_.sample_vertices[item.first];
		if (neighbor == vertex)
			continue;
		const uint32 neighbor_index = index_of(*s_.data.samples, neighbor);
		neighbors.push_back((*s_.sample_position)[neighbor_index]);
		if (static_cast<int>(neighbors.size()) >= k)
			break;
	}
	if (neighbors.size() < 3)
		return false;
	Vec3 mean(0, 0, 0);
	for (const Vec3& p : neighbors)
		mean += p;
	mean /= Scalar(neighbors.size());
	Eigen::Matrix<Scalar, 3, 3> covariance = Eigen::Matrix<Scalar, 3, 3>::Zero();
	for (const Vec3& p : neighbors)
	{
		const Vec3 d = p - mean;
		covariance(0, 0) += d.x() * d.x();
		covariance(0, 1) += d.x() * d.y();
		covariance(0, 2) += d.x() * d.z();
		covariance(1, 1) += d.y() * d.y();
		covariance(1, 2) += d.y() * d.z();
		covariance(2, 2) += d.z() * d.z();
	}
	covariance(1, 0) = covariance(0, 1);
	covariance(2, 0) = covariance(0, 2);
	covariance(2, 1) = covariance(1, 2);
	Eigen::SelfAdjointEigenSolver<Eigen::Matrix<Scalar, 3, 3>> solver(covariance);
	if (solver.info() != Eigen::Success)
		return false;
	const auto eigenvector = solver.eigenvectors().col(0);
	Vec3 normal(eigenvector(0), eigenvector(1), eigenvector(2));
	if (normal.squaredNorm() < Scalar(1e-12))
		return false;
	const Vec3 current = (*s_.sample_normal)[vertex_index];
	if (current.squaredNorm() > Scalar(1e-12) && normal.dot(current) < Scalar(0))
		normal = -normal;
	normal.normalize();
	(*s_.sample_normal)[vertex_index] = normal;
	return true;
}

template <typename SURFACE, typename POINTS, typename NONMANIFOLD, typename RaySamplerTag>
typename UDFReconstruction<SURFACE, POINTS, NONMANIFOLD, RaySamplerTag>::SpheresOptimizerType& UDFReconstruction<
	SURFACE, POINTS, NONMANIFOLD, RaySamplerTag>::optimizer_for()
{
	typename SpheresOptimizerType::Data d;
	d.samples_mesh = s_.data.samples;
	d.spheres = s_.data.spheres;
	d.sample_position = s_.sample_position;
	d.sample_area = s_.sample_area;
	d.sample_knn = s_.sample_knn;
	d.sample_quadric = s_.sample_quadric;
	d.sample_line_quadric = s_.sample_line_quadric;
	d.sample_ma_position = s_.sample_ma_position;
	d.sample_ma_radius = s_.sample_ma_radius;
	d.sample_ma_secondary_vertex = s_.sample_ma_secondary;
	d.sample_sphere = s_.sample_sphere;
	d.sample_error = s_.sample_error;
	d.sphere_position = s_.sphere_position;
	d.sphere_radius = s_.sphere_radius;
	d.sphere_cluster = s_.sphere_cluster;
	d.sphere_cluster_area = s_.sphere_cluster_area;
	d.sphere_cluster_color = s_.sphere_cluster_color;
	d.sphere_neighbors = s_.sphere_neighbors;
	d.sphere_error = s_.sphere_error;
	d.sphere_error_not_normalized = s_.sphere_error_raw;
	d.sample_kdtree = s_.sample_kdtree.get();
	d.sample_kdtree_vertices = &s_.sample_vertices;
	d.sqem_update_lambda_line_plane = Scalar(s_.options.sqem_update_lambda_line_plane);
	if (!s_.optimizer)
		s_.optimizer = std::make_unique<SpheresOptimizerType>(std::move(d));
	else
		s_.optimizer->data() = std::move(d);
	return *s_.optimizer;
}

template <typename SURFACE, typename POINTS, typename NONMANIFOLD, typename RaySamplerTag>
typename UDFReconstruction<SURFACE, POINTS, NONMANIFOLD, RaySamplerTag>::TopologyType& UDFReconstruction<
	SURFACE, POINTS, NONMANIFOLD, RaySamplerTag>::topology_for()
{
	typename TopologyType::Data d;
	d.skeleton = s_.data.skeleton;
	d.spheres_optimizer = &optimizer_for();
	d.sphere_skeleton_vertex = s_.sphere_skeleton_vertex;
	d.skeleton_position = s_.skeleton_position;
	d.skeleton_radius = s_.skeleton_radius;
	d.skeleton_source_sphere = s_.skeleton_source_sphere;
	d.face_incident_tets = s_.skeleton_incident_tets;
	d.edge_degree = s_.skeleton_edge_degree;
	if (!s_.topology)
		s_.topology = std::make_unique<TopologyType>(std::move(d));
	else
		s_.topology->data() = std::move(d);
	return *s_.topology;
}

template <typename SURFACE, typename POINTS, typename NONMANIFOLD, typename RaySamplerTag>
typename UDFReconstruction<SURFACE, POINTS, NONMANIFOLD, RaySamplerTag>::Status UDFReconstruction<
	SURFACE, POINTS, NONMANIFOLD, RaySamplerTag>::sample_alpha_level_set()
{
	if (initialize_data() != Status::success)
		return Status::invalid_data;
	s_.input.bind(s_.data.surface, s_.data.input_points, s_.input_position, s_.input_normal);
	s_.input.invalidate();
	s_.input.build_input_kdtree();
	const PrepareInputs<SURFACE, POINTS> inputs{s_.data.surface, s_.data.input_points, s_.input_position, s_.input_normal,
												s_.input_knn};
	if (const FieldPrepareResult prepared = s_.field.prepare(inputs); prepared != FieldPrepareResult::ok)
		return status_of(prepared);

	Vec3 bbox_min(0, 0, 0), bbox_max(1, 1, 1);
	if (const FieldPrepareResult boxed = s_.field.sampling_box(inputs, bbox_min, bbox_max);
		boxed != FieldPrepareResult::ok)
		return status_of(boxed);
	const Scalar padding = udf_constants::sampling_padding;
	bbox_min -= Vec3(padding, padding, padding);
	bbox_max += Vec3(padding, padding, padding);

	RayTraits traits = s_.field.ray_traits(inputs);
	if (!traits.is_ready())
		return Status::invalid_data;
	RayParams ray_params;
	ray_params.bbox_expand = udf_constants::ray_box_relative_expand;
	ray_params.alpha = Scalar(s_.options.alpha);
	ray_params.tol = Scalar(s_.options.tolerance);
	ray_params.step_bound = Scalar(2.0);
	ray_params.batch_size = std::max(1, s_.options.ray_sampler_batch_size);
	ray_params.max_iterations = s_.options.udf_max_iterations;
	ray_params.min_active_rays = s_.options.udf_min_active_rays;
	ray_params.max_outer_iterations = 10;
	ray_params.seed = s_.options.seed;
	if (!s_.ray_sampler)
		s_.ray_sampler = s_.field.make_ray_sampler(ray_params);
	else
		s_.field.configure_ray_sampler(*s_.ray_sampler, ray_params);
	AlphaSamplingParameters params;
	params.sample_radius = Scalar(s_.options.sample_radius);
	params.seed = static_cast<uint32>(s_.options.seed);
	params.ray_sampler_batch_size = static_cast<size_t>(std::max(1, s_.options.ray_sampler_batch_size));
	params.candidates_per_parent = static_cast<size_t>(std::max(1, s_.options.alpha_sampling_candidates_per_parent));
	params.parent_batch_size = static_cast<size_t>(std::max(1, s_.options.alpha_sampling_parent_batch_size));
	params.bbox_min = bbox_min;
	params.bbox_max = bbox_max;
	auto compute_normals = [&](const std::vector<Vec3>& points, std::vector<Vec3>& normals,
							   std::vector<uint8_t>& keep) {
		keep.assign(points.size(), uint8_t(0));
		normals.assign(points.size(), Vec3(0, 0, 1));
		return s_.field.normals(points, normals, keep);
	};
	const AlphaProjectionParameters pp{Scalar(s_.options.alpha), Scalar(s_.options.tolerance), bbox_min, bbox_max, 1};
	auto project_points = [&](const std::vector<Vec3>& points, AlphaProjectionResult& projection) {
		return s_.field.project(points, pp, projection);
	};
	AlphaSamplingResult result =
		geometry::sample_alpha_level_set(traits, *s_.ray_sampler, params, compute_normals, project_points);
	if (!result.success || result.positions.empty())
		return Status::sampling_failed;
	clear(*s_.data.samples);
	s_.samples_before_filtering = static_cast<uint32>(result.positions.size());
	for (size_t i = 0; i < result.positions.size(); ++i)
	{
		PVertex v = add_vertex(*s_.data.samples);
		const uint32 idx = index_of(*s_.data.samples, v);
		(*s_.sample_position)[idx] = result.positions[i];
		(*s_.sample_normal)[idx] = result.normals[i];
	}
	build_sample_kdtree();
	return Status::success;
}

template <typename SURFACE, typename POINTS, typename NONMANIFOLD, typename RaySamplerTag>
typename UDFReconstruction<SURFACE, POINTS, NONMANIFOLD, RaySamplerTag>::Status UDFReconstruction<
	SURFACE, POINTS, NONMANIFOLD, RaySamplerTag>::apply_sampling_filtering()
{
	if (!s_.options.apply_filtering || !s_.data.samples || nb_cells<PVertex>(*s_.data.samples) == 0)
		return Status::success;
	s_.input.bind(s_.data.surface, s_.data.input_points, s_.input_position, s_.input_normal);
	s_.input.build_input_kdtree();
	if (!s_.input.input_kdtree())
		s_.input.build_surface_bvh();
	std::vector<Vec3> kept, kept_normals;
	kept.reserve(nb_cells<PVertex>(*s_.data.samples));
	foreach_cell(*s_.data.samples, [&](PVertex v) {
		const Vec3 p = (*s_.sample_position)[index_of(*s_.data.samples, v)];
		bool within_distance = false;
		const Scalar max_distance = Scalar(s_.options.alpha) + udf_constants::filter_margin;
		if (s_.input.input_kdtree())
		{
			std::pair<uint32, Scalar> nearest;
			within_distance = s_.input.input_kdtree()->find_nn(p, &nearest, max_distance);
		}
		else if (s_.input.surface_bvh())
		{
			std::pair<uint32, Vec3> closest;
			within_distance = s_.input.surface_bvh()->closest_point(p, &closest, max_distance);
		}
		if (within_distance)
		{
			kept.push_back(p);
			kept_normals.push_back((*s_.sample_normal)[index_of(*s_.data.samples, v)]);
		}
		return true;
	});
	const size_t before = nb_cells<PVertex>(*s_.data.samples);
	if (kept.size() == before)
		return Status::success;
	clear(*s_.data.samples);
	for (size_t i = 0; i < kept.size(); ++i)
	{
		PVertex v = add_vertex(*s_.data.samples);
		(*s_.sample_position)[index_of(*s_.data.samples, v)] = kept[i];
		(*s_.sample_normal)[index_of(*s_.data.samples, v)] = kept_normals[i];
	}
	s_.fitting_data_computed = false;
	if (kept.empty())
	{
		s_.sample_kdtree.reset();
		s_.sample_vertices.clear();
		return Status::success;
	}
	// neural samples keep their projection normals (the network gradient); for point cloud / surface input the
	// normal is recomputed at the projected position from the input
	if constexpr (Field::recompute_normals_from_input_after_filtering)
		recompute_samples_normals_from_input();
	build_sample_kdtree();
	return Status::success;
}

template <typename SURFACE, typename POINTS, typename NONMANIFOLD, typename RaySamplerTag>
typename UDFReconstruction<SURFACE, POINTS, NONMANIFOLD, RaySamplerTag>::Status UDFReconstruction<
	SURFACE, POINTS, NONMANIFOLD, RaySamplerTag>::build_kdtree_and_normals()
{
	if (!s_.data.samples || nb_cells<PVertex>(*s_.data.samples) == 0)
		return Status::invalid_data;
	if (!s_.sample_kdtree)
		return Status::invalid_data;
	// PCA for every input: neural samples already carry the network gradient, and for MF models that gradient
	// is not orthogonal to the medial axis.
	if (s_.options.recompute_normals_after_sampling)
		parallel_foreach_cell(*s_.data.samples, [&](PVertex v) {
			recompute_sample_normal_pca(v);
			return true;
		});
	return Status::success;
}

template <typename SURFACE, typename POINTS, typename NONMANIFOLD, typename RaySamplerTag>
typename UDFReconstruction<SURFACE, POINTS, NONMANIFOLD, RaySamplerTag>::Status UDFReconstruction<
	SURFACE, POINTS, NONMANIFOLD, RaySamplerTag>::compute_fitting_primitives()
{
	if (!s_.data.samples || !s_.sample_kdtree)
		return Status::invalid_data;
	compute_samples_area(*s_.data.samples, *s_.sample_position, *s_.sample_normal, *s_.sample_area, *s_.sample_knn,
						 *s_.sample_kdtree, s_.sample_vertices, s_.options.knn_k, Scalar(s_.options.alpha));
	compute_quadrics(*s_.data.samples, *s_.sample_position, *s_.sample_normal, *s_.sample_area, *s_.sample_knn, *s_.sample_quadric,
					 *s_.sample_line_quadric, s_.options.knn_k);
	s_.fitting_data_computed = true;
	return Status::success;
}

template <typename SURFACE, typename POINTS, typename NONMANIFOLD, typename RaySamplerTag>
typename UDFReconstruction<SURFACE, POINTS, NONMANIFOLD, RaySamplerTag>::Status UDFReconstruction<
	SURFACE, POINTS, NONMANIFOLD, RaySamplerTag>::compute_initial_medial_axis()
{
	if (!s_.data.samples || !s_.sample_kdtree)
		return Status::invalid_data;
	const Scalar fallback_radius = Scalar(s_.options.alpha);
	const Scalar min_norm = Scalar(1e-12);
	const bool neural_input = Field::is_neural_type && s_.field.loaded();
	const bool mf_model = neural_input && s_.options.neural_model_type == NeuralModelType::mf;
	const bool udf_model = neural_input && s_.options.neural_model_type == NeuralModelType::udf;
	const Scalar initial_radius = std::max<Scalar>(fallback_radius * Scalar(10), Scalar(0));
	auto run_shrinking_ball_for_vertex = [&](PVertex v) -> bool {
		if (!v.is_valid())
			return false;
		const uint32 idx = index_of(*s_.data.samples, v);
		if (idx == INVALID_INDEX)
			return false;
		const Vec3& pt = (*s_.sample_position)[idx];
		Vec3 n = (*s_.sample_normal)[idx];
		const bool normal_finite = n.allFinite();
		const bool can_shrink = normal_finite && n.squaredNorm() > min_norm;
		if (can_shrink)
			n.normalize();
		else
			n = Vec3(0, 0, 1);
		Vec3 center = pt - n * fallback_radius;
		Scalar radius = fallback_radius;
		PVertex secondary;
		if (can_shrink)
		{
			auto [c1, r1, q1] =
				geometry::shrinking_ball_center<PVertex>(pt, n, s_.sample_kdtree.get(), s_.sample_vertices, initial_radius);
			if (c1.allFinite() && std::isfinite(static_cast<double>(r1)) && r1 > Scalar(0))
			{
				center = c1;
				radius = r1;
				secondary = q1;
			}
		}
		if (!center.allFinite())
			center = pt - n * fallback_radius;
		if (!std::isfinite(static_cast<double>(radius)) || radius <= Scalar(0))
			radius = fallback_radius;
		if (!secondary.is_valid() && s_.sample_kdtree && !s_.sample_vertices.empty())
		{
			std::pair<uint32, Scalar> nn;
			if (s_.sample_kdtree->find_nn(center - n * radius, &nn) && nn.first < s_.sample_vertices.size())
				secondary = s_.sample_vertices[nn.first];
		}
		(*s_.sample_ma_position)[idx] = center;
		(*s_.sample_ma_radius)[idx] = std::max(radius, fallback_radius);
		(*s_.sample_ma_secondary)[idx] = secondary;
		return true;
	};
	auto run_all = [&] {
		parallel_foreach_cell(*s_.data.samples, [&](PVertex v) {
			run_shrinking_ball_for_vertex(v);
			return true;
		});
	};
	run_all();
	if (s_.options.ma_flip_prune)
	{
		uint32 flipped = 0, deleted = 0;
		std::vector<uint32> flipped_indices;
		const Scalar udf_threshold =
			udf_model ? udf_constants::udf_flip_threshold_factor * fallback_radius : fallback_radius;
		const Scalar mf_threshold =
			fallback_radius * std::max(udf_constants::mf_flip_alpha_factor_min,
									   std::min(udf_constants::mf_flip_alpha_factor_max,
												Scalar(s_.options.ma_flip_prune_alpha_factor)));
		const bool supported = mf_model || !Field::is_neural_type || neural_input;
		auto scores = [&](const std::vector<Vec3>& points, std::vector<Scalar>& values,
						  std::vector<Scalar>& sdf) -> bool {
			sdf.clear();
			if (!mf_model)
				return evaluate_topology_scores(points, values);
			return s_.field.values_sdf(points, values, sdf);
		};
		if (supported)
		{
			std::vector<PVertex> vertices;
			std::vector<Vec3> centers;
			vertices.reserve(nb_cells<PVertex>(*s_.data.samples));
			centers.reserve(vertices.capacity());
			foreach_cell(*s_.data.samples, [&](PVertex v) {
				const uint32 i = index_of(*s_.data.samples, v);
				vertices.push_back(v);
				centers.push_back((*s_.sample_ma_position)[i]);
				return true;
			});
			std::vector<Scalar> values, sdf;
			if (scores(centers, values, sdf) && values.size() == centers.size())
			{
				std::vector<PVertex> retry;
				for (size_t i = 0; i < vertices.size(); ++i)
					if ((mf_model && i < sdf.size() && sdf[i] > Scalar(0)) ||
						values[i] > (mf_model ? mf_threshold : udf_threshold))
						retry.push_back(vertices[i]);
				flipped_indices.reserve(retry.size());
				for (PVertex v : retry)
				{
					const uint32 i = index_of(*s_.data.samples, v);
					if (i == INVALID_INDEX)
						continue;
					Vec3 n = (*s_.sample_normal)[i];
					if (n.squaredNorm() > min_norm)
					{
						n.normalize();
						(*s_.sample_normal)[i] = -n;
						flipped_indices.push_back(i);
						++flipped;
					}
					run_shrinking_ball_for_vertex(v);
				}
				auto prune = [&](const std::vector<PVertex>& candidates) {
					std::vector<Vec3> p;
					std::vector<PVertex> active;
					for (PVertex v : candidates)
					{
						uint32 i = index_of(*s_.data.samples, v);
						if (i != INVALID_INDEX)
						{
							active.push_back(v);
							p.push_back((*s_.sample_ma_position)[i]);
						}
					}
					std::vector<Scalar> val, sd;
					if (!scores(p, val, sd) || val.size() != p.size())
						return uint32(0);
					uint32 removed = 0;
					for (size_t i = 0; i < active.size(); ++i)
						if ((mf_model && i < sd.size() && sd[i] > Scalar(0)) ||
							val[i] > (mf_model ? mf_threshold : udf_threshold))
						{
							remove_vertex(*s_.data.samples, active[i]);
							++removed;
						}
					return removed;
				};
				deleted += prune(retry);
				if (deleted > 0 && !flipped_indices.empty() && nb_cells<PVertex>(*s_.data.samples) > 0)
				{
					build_sample_kdtree();
					std::vector<PVertex> survivors;
					survivors.reserve(flipped_indices.size());
					for (uint32 i : flipped_indices)
					{
						PVertex v = of_index<PVertex>(*s_.data.samples, i);
						if (v.is_valid())
							survivors.push_back(v);
					}
					for (PVertex v : survivors)
						recompute_sample_normal_pca(v);
					for (PVertex v : survivors)
						run_shrinking_ball_for_vertex(v);
					deleted += prune(survivors);
				}
			}
		}
	}
	build_sample_kdtree();
	return Status::success;
}

template <typename SURFACE, typename POINTS, typename NONMANIFOLD, typename RaySamplerTag>
typename UDFReconstruction<SURFACE, POINTS, NONMANIFOLD, RaySamplerTag>::Status UDFReconstruction<
	SURFACE, POINTS, NONMANIFOLD, RaySamplerTag>::initialize_spheres()
{
	if (!s_.fitting_data_computed || !s_.data.samples || !s_.data.spheres)
		return Status::invalid_data;
	if (!optimizer_for().initialize_from_samples(Scalar(s_.options.init_dilation_constant)))
		return Status::sphere_initialization_failed;
	return Status::success;
}

template <typename SURFACE, typename POINTS, typename NONMANIFOLD, typename RaySamplerTag>
typename UDFReconstruction<SURFACE, POINTS, NONMANIFOLD, RaySamplerTag>::Status UDFReconstruction<
	SURFACE, POINTS, NONMANIFOLD, RaySamplerTag>::optimize_spheres(const IterationCallback& callback)
{
	auto& optimizer = optimizer_for();
	optimizer.reset_optimization();
	typename SpheresOptimizerType::Status status = SpheresOptimizerType::Status::running;
	while (status == SpheresOptimizerType::Status::running)
	{
		status = optimizer.update_once(Scalar(s_.options.sqem_update_lambda_line_plane));
		if (callback)
			callback(optimizer.metrics());
	}
	return status == SpheresOptimizerType::Status::failed ? Status::optimization_failed : Status::success;
}

template <typename SURFACE, typename POINTS, typename NONMANIFOLD, typename RaySamplerTag>
typename UDFReconstruction<SURFACE, POINTS, NONMANIFOLD, RaySamplerTag>::Status UDFReconstruction<
	SURFACE, POINTS, NONMANIFOLD, RaySamplerTag>::build_skeleton()
{
	return topology_for().build_from_spheres() == TopologyType::Status::success ? Status::success
																				: Status::skeleton_build_failed;
}

template <typename SURFACE, typename POINTS, typename NONMANIFOLD, typename RaySamplerTag>
typename UDFReconstruction<SURFACE, POINTS, NONMANIFOLD, RaySamplerTag>::Status UDFReconstruction<
	SURFACE, POINTS, NONMANIFOLD, RaySamplerTag>::fix_topology()
{
	return topology_for().fix_topology([this](const std::vector<Vec3>& p, std::vector<Scalar>& v) {
		return evaluate_topology_scores(p, v);
	}) == TopologyType::Status::success
			   ? Status::success
			   : Status::topology_fix_failed;
}

template <typename SURFACE, typename POINTS, typename NONMANIFOLD, typename RaySamplerTag>
typename UDFReconstruction<SURFACE, POINTS, NONMANIFOLD, RaySamplerTag>::Status UDFReconstruction<
	SURFACE, POINTS, NONMANIFOLD, RaySamplerTag>::prune_residual_sheets()
{
	topology_for().set_residual_prune_mode(
		static_cast<typename TopologyType::ResidualPruneMode>(static_cast<int>(s_.options.residual_prune_mode)));
	return topology_for().prune_residual_sheets() == TopologyType::Status::success ? Status::success
																	   : Status::topology_fix_failed;
}

template <typename SURFACE, typename POINTS, typename NONMANIFOLD, typename RaySamplerTag>
bool UDFReconstruction<SURFACE, POINTS, NONMANIFOLD, RaySamplerTag>::evaluate_topology_scores(
	const std::vector<Vec3>& points, std::vector<Scalar>& values)
{
	s_.input.bind(s_.data.surface, s_.data.input_points, s_.input_position, s_.input_normal);
	if (Field::is_neural_type && s_.options.neural_model_type == NeuralModelType::mf)
	{
		// MF: mf - |sdf| is not a stable medial axis (its minimum is above zero and varies a lot over a shape),
		// so the score is the distance to the medial point of the nearest sample instead of the network value
		if (!s_.sample_ma_kdtree)
			return false;
		values.resize(points.size(), Scalar(0));
		detail::parallel_for_index(points.size(), [&](std::size_t i) {
			std::pair<uint32, Scalar> nearest;
			if (!s_.sample_ma_kdtree->find_nn(points[i], &nearest) || nearest.first >= s_.sample_ma_vertices.size())
				return;
			const uint32 sid = index_of(*s_.data.samples, s_.sample_ma_vertices[nearest.first]);
			if (sid != INVALID_INDEX)
				values[i] = (points[i] - (*s_.sample_ma_position)[sid]).norm();
		});
		return true;
	}
	return s_.field.distances(points, values);
}

template <typename SURFACE, typename POINTS, typename NONMANIFOLD, typename RaySamplerTag>
typename UDFReconstruction<SURFACE, POINTS, NONMANIFOLD, RaySamplerTag>::Result UDFReconstruction<
	SURFACE, POINTS, NONMANIFOLD, RaySamplerTag>::run(const IterationCallback& callback)
{
	Result r;
	s_.input.invalidate();
	if (s_.data.samples)
		clear(*s_.data.samples);
	if (s_.data.spheres)
		clear(*s_.data.spheres);
	if (s_.data.skeleton)
		clear(*s_.data.skeleton);
	s_.sample_kdtree.reset();
	s_.sample_ma_kdtree.reset();
	s_.sample_vertices.clear();
	s_.sample_ma_vertices.clear();
	s_.optimizer.reset();
	s_.topology.reset();
	s_.fitting_data_computed = false;
	s_.samples_before_filtering = 0;
	auto run_stage = [](float64& time, auto&& stage) {
		const auto start = std::chrono::high_resolution_clock::now();
		const Status status = stage();
		time = std::chrono::duration<float64, std::milli>(std::chrono::high_resolution_clock::now() - start).count();
		return status;
	};
	if ((r.status = run_stage(r.timing.preprocessing_ms, [&] { return initialize_data(); })) != Status::success)
		return r;
	if ((r.status = run_stage(r.timing.sampling_ms, [&] { return sample_alpha_level_set(); })) != Status::success)
		return r;
	r.counts.samples_before_filtering = s_.samples_before_filtering;
	if ((r.status = run_stage(r.timing.sample_filtering_ms, [&] { return apply_sampling_filtering(); })) !=
		Status::success)
		return r;
	if ((r.status = run_stage(r.timing.kdtree_and_normals_ms, [&] { return build_kdtree_and_normals(); })) !=
		Status::success)
		return r;
	if ((r.status = run_stage(r.timing.fitting_primitives_ms, [&] { return compute_fitting_primitives(); })) !=
		Status::success)
		return r;
	if ((r.status = run_stage(r.timing.initial_medial_axis_ms, [&] { return compute_initial_medial_axis(); })) !=
		Status::success)
		return r;
	if ((r.status = run_stage(r.timing.sphere_initialization_ms, [&] { return initialize_spheres(); })) !=
		Status::success)
		return r;
	if ((r.status = run_stage(r.timing.optimization_ms, [&] { return optimize_spheres(callback); })) != Status::success)
		return r;
	if (s_.optimizer)
	{
		r.timing.cluster_ms = s_.optimizer->metrics().cluster_total_ms;
		r.timing.sphere_update_ms = s_.optimizer->metrics().sphere_update_total_ms;
		r.timing.error_ms = s_.optimizer->metrics().error_total_ms;
	}
	if ((r.status = run_stage(r.timing.skeleton_construction_ms, [&] { return build_skeleton(); })) != Status::success)
		return r;
	if ((r.status = run_stage(r.timing.topology_processing_ms, [&] { return fix_topology(); })) != Status::success)
		return r;
	if (s_.options.residual_prune)
	{
		float64 residual_ms = 0.0;
		if ((r.status = run_stage(residual_ms, [&] { return prune_residual_sheets(); })) != Status::success)
			return r;
		r.timing.topology_processing_ms += residual_ms;
	}
	r.status = Status::success;
	r.counts = counts();
	return r;
}

} // namespace cgogn::geometry

#endif // CGOGN_GEOMETRY_ALGOS_UDF_RECONSTRUCTION_H_
