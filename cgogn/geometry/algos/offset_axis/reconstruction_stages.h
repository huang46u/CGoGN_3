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
#ifndef CGOGN_GEOMETRY_ALGOS_OFFSET_AXIS_RECONSTRUCTION_STAGES_H_
#define CGOGN_GEOMETRY_ALGOS_OFFSET_AXIS_RECONSTRUCTION_STAGES_H_

#include <cgogn/core/functions/attributes.h>
#include <cgogn/core/functions/mesh_info.h>
#include <cgogn/core/functions/traversals/edge.h>
#include <cgogn/core/types/cell_marker.h>
#include <cgogn/core/utils/thread_pool.h>
#include <cgogn/geometry/algos/area.h>
#include <cgogn/geometry/algos/fitting.h>
#include <cgogn/geometry/algos/medial_axis.h>
#include <cgogn/geometry/algos/normal.h>
#include <cgogn/geometry/algos/offset_axis/alpha_projection.h>
#include <cgogn/geometry/algos/offset_axis/alpha_sampling.h>
#include <cgogn/geometry/algos/offset_axis/distance_field.h>
#include <cgogn/geometry/algos/offset_axis/input_geometry.h>
#include <cgogn/geometry/algos/offset_axis/neural_model_type.h>
#include <cgogn/geometry/algos/offset_axis/normal_utils.h>
#include <cgogn/geometry/algos/offset_axis/parallel_for.h>
#include <cgogn/geometry/algos/offset_axis/point_cloud_distance_field.h>
#include <cgogn/geometry/algos/offset_axis/reconstruction_options.h>
#include <cgogn/geometry/algos/offset_axis/reconstruction_state.h>
#include <cgogn/geometry/algos/offset_axis/reconstruction_stats.h>
#include <cgogn/geometry/algos/offset_axis/sample_processing.h>
#include <cgogn/geometry/algos/offset_axis/skeleton_topology.h>
#include <cgogn/geometry/algos/offset_axis/spatial_query.h>
#include <cgogn/geometry/algos/offset_axis/spheres_optimizer.h>
#include <cgogn/geometry/algos/offset_axis/surface_distance_field.h>
#include <cgogn/geometry/functions/bounding_box.h>
#include <cgogn/geometry/functions/normal.h>
#include <cgogn/geometry/types/line_quadric.h>
#include <cgogn/geometry/algos/offset_axis/ray_level_set_sampler.h>
#include <cgogn/geometry/algos/offset_axis/ray_level_set_sampler_traits.h>
#include <cgogn/geometry/algos/offset_axis/spatial_grid.h>
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

namespace cgogn
{

namespace geometry
{

namespace offset_axis
{

// Fixed constants of the reconstruction (not user options). Lengths are absolute, in model units (the inputs are
// rescaled to [0,1] before the reconstruction).
namespace constants
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
} // namespace constants

} // namespace offset_axis

} // namespace geometry

} // namespace cgogn

// Stage precondition: when cond does not hold, the enclosing stage function (template parameter State) returns code
// with a message naming the stage and the condition.
#define OFFSET_AXIS_PRECONDITION(cond, code)                                                                                   \
	do                                                                                                                 \
	{                                                                                                                  \
		if (!(cond))                                                                                                   \
			return typename State::StatusInfo{(code), std::string(__func__) + ": precondition failed: " #cond};        \
	} while (false)

namespace cgogn
{

namespace geometry
{

namespace offset_axis
{

namespace stages
{

// Stage free functions of the UDF reconstruction. They operate on a State (ReconstructionState) and are
// called by the OffsetAxisReconstruction facade. reconstruct() is the orchestration of the whole pipeline.

template <typename State>
void reset(State& s);

template <typename State>
typename State::StatusInfo initialize_data(State& s);

template <typename State>
typename State::StatusInfo load_neural_model(State& s, const std::filesystem::path& path,
											 typename State::NeuralModelType type);

template <typename State>
typename State::Counts counts(State& s);

template <typename State>
typename State::Repairs repairs(State& s);

template <typename State>
void build_sample_kdtree(State& s);

template <typename State>
void recompute_samples_normals_from_input(State& s);

template <typename State>
bool recompute_sample_normal_pca(State& s, typename State::PVertex vertex);

template <typename State>
typename State::SpheresOptimizerType& optimizer_for(State& s);

template <typename State>
typename State::TopologyType& topology_for(State& s);

// true if the field is neural and holds a model of the given type (CPU fields: false)
template <typename State>
bool neural_model_is(const State& s, NeuralModelType type)
{
	if constexpr (State::Field::is_neural_type)
		return s.field.model_type() == type;
	else
		return false;
}

template <typename State>
bool evaluate_topology_scores(State& s, const std::vector<Vec3>& points, std::vector<Scalar>& values);

template <typename State>
typename State::Status status_of(State&, FieldPrepareResult result);

template <typename State>
typename State::StatusInfo sample_alpha_level_set(State& s);

template <typename State>
typename State::StatusInfo apply_sampling_filtering(State& s);

template <typename State>
typename State::StatusInfo build_kdtree_and_normals(State& s);

template <typename State>
typename State::StatusInfo compute_fitting_primitives(State& s);

template <typename State>
typename State::StatusInfo compute_initial_medial_axis(State& s);

template <typename State>
typename State::StatusInfo initialize_spheres(State& s);

template <typename State>
typename State::StatusInfo optimize_spheres(State& s, const typename State::IterationCallback& callback = {});

template <typename State>
typename State::StatusInfo build_skeleton(State& s);

template <typename State>
typename State::StatusInfo fix_topology(State& s);

template <typename State>
typename State::StatusInfo prune_residual_sheets(State& s);

// Composite stages shared by the UI buttons and the stage sequence test.
template <typename State>
typename State::StatusInfo sample(State& s);

template <typename State>
typename State::StatusInfo compute_fitting_data(State& s);

template <typename State>
typename State::SpheresOptimizerType::Status update_spheres_once(State& s, Scalar sqem_update_lambda_line_plane);

template <typename State>
typename State::StatusInfo fix_topology_and_prune(State& s);

template <typename State>
typename State::Result reconstruct(State& s, const typename State::IterationCallback& callback = {});

template <typename State>
void reset(State& s)
{
	s.input.invalidate();
	s.sample_kdtree.reset();
	s.sample_ma_kdtree.reset();
	s.sample_vertices.clear();
	s.sample_ma_vertices.clear();
	s.ray_sampler.reset();
	s.optimizer.reset();
	s.topology.reset();
	s.field.reset();
	s.fitting_data_computed = false;
	s.samples_before_filtering = 0;
	s.ma_repairs = {};
	s.initialized = false;
}

template <typename State>
typename State::StatusInfo initialize_data(State& s)
{
	using PVertex = typename State::PVertex;
	using NMVertex = typename State::NMVertex;
	using NMEdge = typename State::NMEdge;
	using NMFace = typename State::NMFace;
	using Status = typename State::Status;

	if (s.initialized)
		return Status::success;
	OFFSET_AXIS_PRECONDITION(s.data.input_points && s.data.samples && s.data.spheres && s.data.skeleton, Status::invalid_data);
	s.input_position = get_or_add_attribute<Vec3, PVertex>(*s.data.input_points, "position");
	s.input_normal = get_or_add_attribute<Vec3, PVertex>(*s.data.input_points, "normal");
	s.input_knn = get_or_add_attribute<std::vector<PVertex>, PVertex>(*s.data.input_points, "knn");
	s.sample_position = get_or_add_attribute<Vec3, PVertex>(*s.data.samples, "position");
	s.sample_normal = get_or_add_attribute<Vec3, PVertex>(*s.data.samples, "normal");
	s.sample_area = get_or_add_attribute<Scalar, PVertex>(*s.data.samples, "area");
	s.sample_knn = get_or_add_attribute<std::vector<PVertex>, PVertex>(*s.data.samples, "knn");
	s.sample_quadric = get_or_add_attribute<Spherical_Quadric, PVertex>(*s.data.samples, "quadric");
	s.sample_line_quadric = get_or_add_attribute<Line_Quadric, PVertex>(*s.data.samples, "line_quadric");
	s.sample_ma_position = get_or_add_attribute<Vec3, PVertex>(*s.data.samples, "ma_position");
	s.sample_ma_radius = get_or_add_attribute<Scalar, PVertex>(*s.data.samples, "ma_radius");
	s.sample_ma_secondary = get_or_add_attribute<PVertex, PVertex>(*s.data.samples, "ma_secondary_vertex");
	s.sample_sphere = get_or_add_attribute<PVertex, PVertex>(*s.data.samples, "sphere");
	s.sample_error = get_or_add_attribute<Scalar, PVertex>(*s.data.samples, "error");
	s.sphere_position = get_or_add_attribute<Vec3, PVertex>(*s.data.spheres, "position");
	s.sphere_radius = get_or_add_attribute<Scalar, PVertex>(*s.data.spheres, "radius");
	s.sphere_cluster = get_or_add_attribute<std::vector<PVertex>, PVertex>(*s.data.spheres, "cluster");
	s.sphere_cluster_area = get_or_add_attribute<Scalar, PVertex>(*s.data.spheres, "cluster_area");
	s.sphere_cluster_color = get_or_add_attribute<Vec4, PVertex>(*s.data.spheres, "cluster_color");
	s.sphere_neighbors = get_or_add_attribute<std::vector<PVertex>, PVertex>(*s.data.spheres, "neighbor_clusters");
	s.sphere_error = get_or_add_attribute<Scalar, PVertex>(*s.data.spheres, "error");
	s.sphere_error_raw = get_or_add_attribute<Scalar, PVertex>(*s.data.spheres, "error_not_normalized");
	s.sphere_skeleton_vertex = get_or_add_attribute<NMVertex, PVertex>(*s.data.spheres, "skeleton_vertex");
	s.skeleton_position = get_or_add_attribute<Vec3, NMVertex>(*s.data.skeleton, "position");
	s.skeleton_radius = get_or_add_attribute<Scalar, NMVertex>(*s.data.skeleton, "radius");
	s.skeleton_source_sphere = get_or_add_attribute<PVertex, NMVertex>(*s.data.skeleton, "source_sphere");
	s.skeleton_incident_tets = get_or_add_attribute<std::set<std::size_t>, NMFace>(*s.data.skeleton, "incident_tets");
	s.skeleton_edge_degree = get_or_add_attribute<uint32, NMEdge>(*s.data.skeleton, "degree");
	s.initialized =
		s.input_position && s.sample_position && s.sample_normal && s.sphere_position && s.skeleton_position;
	s.input.bind(s.data.surface, s.data.input_points, s.input_position, s.input_normal);
	return s.initialized ? Status::success : Status::invalid_data;
}

template <typename State>
typename State::StatusInfo load_neural_model(State& s, const std::filesystem::path& path,
											 typename State::NeuralModelType type)
{
	using Status = typename State::Status;

	if (!s.field.load_model(path, type))
		return typename State::StatusInfo{Status::model_load_failed, "load_neural_model: " + s.field.load_error()};
	return Status::success;
}

template <typename State>
typename State::Counts counts(State& s)
{
	using PVertex = typename State::PVertex;
	using SVertex = typename State::SVertex;
	using NMVertex = typename State::NMVertex;
	using NMEdge = typename State::NMEdge;
	using NMFace = typename State::NMFace;
	using Counts = typename State::Counts;

	Counts c;
	if (s.data.input_points)
		c.input_points = nb_cells<PVertex>(*s.data.input_points);
	if (s.data.surface)
		c.input_vertices = nb_cells<SVertex>(*s.data.surface);
	else
		c.input_vertices = c.input_points;
	if (s.data.samples)
		c.samples = nb_cells<PVertex>(*s.data.samples);
	c.samples_before_filtering = s.samples_before_filtering;
	if (s.data.spheres)
		c.spheres = nb_cells<PVertex>(*s.data.spheres);
	if (s.data.skeleton)
	{
		c.skeleton_vertices = nb_cells<NMVertex>(*s.data.skeleton);
		c.skeleton_edges = nb_cells<NMEdge>(*s.data.skeleton);
		c.skeleton_faces = nb_cells<NMFace>(*s.data.skeleton);
	}
	if (s.optimizer)
		c.optimization_iterations = s.optimizer->metrics().iteration;
	return c;
}

template <typename State>
typename State::Repairs repairs(State& s)
{
	typename State::Repairs r = s.ma_repairs;
	if (s.optimizer)
		r.sphere_update_nonfinite = s.optimizer->metrics().sphere_update_nonfinite;
	return r;
}

template <typename State>
void build_sample_kdtree(State& s)
{
	using PVertex = typename State::PVertex;

	s.sample_kdtree.reset();
	s.sample_ma_kdtree.reset();
	s.sample_vertices.clear();
	s.sample_ma_vertices.clear();
	if (!s.data.samples || !s.sample_position)
		return;
	std::vector<Vec3> positions;
	std::vector<Vec3> ma_positions;
	foreach_cell(*s.data.samples, [&](PVertex v) {
		const uint32 idx = index_of(*s.data.samples, v);
		positions.push_back((*s.sample_position)[idx]);
		s.sample_vertices.push_back(v);
		const Vec3 ma = (*s.sample_ma_position)[idx];
		if (ma.allFinite() && (*s.sample_ma_radius)[idx] > Scalar(0))
		{
			ma_positions.push_back(ma);
			s.sample_ma_vertices.push_back(v);
		}
		return true;
	});
	if (!positions.empty())
		s.sample_kdtree = std::make_unique<acc::KDTree<3, uint32>>(positions);
	if (!ma_positions.empty())
		s.sample_ma_kdtree = std::make_unique<acc::KDTree<3, uint32>>(ma_positions);
}

template <typename State>
void recompute_samples_normals_from_input(State& s)
{
	using PVertex = typename State::PVertex;

	s.input.bind(s.data.surface, s.data.input_points, s.input_position, s.input_normal);
	if (!s.input.input_kdtree() && s.data.surface)
	{
		s.input.build_surface_bvh();
		if (s.input.surface_bvh())
		{
			foreach_cell(*s.data.samples, [&](PVertex v) {
				const uint32 idx = index_of(*s.data.samples, v);
				const Vec3& position = (*s.sample_position)[idx];
				std::pair<uint32, Vec3> closest;
				Vec3 normal(0, 0, 1);
				if (s.input.surface_bvh()->closest_point(position, &closest))
				{
					normal = unit_normal_or_z(position - closest.second);
				}
				(*s.sample_normal)[idx] = normal;
				return true;
			});
			return;
		}
	}
	if (s.input.input_kdtree() && s.data.input_points)
	{
		foreach_cell(*s.data.samples, [&](PVertex v) {
			const uint32 idx = index_of(*s.data.samples, v);
			const Vec3& position = (*s.sample_position)[idx];
			std::pair<uint32, Scalar> nearest;
			Vec3 normal(0, 0, 1);
			if (s.input.input_kdtree()->find_nn(position, &nearest) && nearest.first < s.input.input_vertices().size())
			{
				const uint32 input_idx = index_of(*s.data.input_points, s.input.input_vertices()[nearest.first]);
				if (input_idx != INVALID_INDEX)
				{
					normal = position - (*s.input_position)[input_idx];
					if (normal.squaredNorm() < Scalar(1e-12))
						normal = (*s.input_normal)[input_idx];
					normal = unit_normal_or_z(normal);
				}
			}
			(*s.sample_normal)[idx] = normal;
			return true;
		});
		return;
	}
	foreach_cell(*s.data.samples, [&](PVertex v) {
		(*s.sample_normal)[index_of(*s.data.samples, v)] = Vec3(0, 0, 1);
		return true;
	});
}

template <typename State>
bool recompute_sample_normal_pca(State& s, typename State::PVertex vertex)
{
	using PVertex = typename State::PVertex;

	if (!s.data.samples || !s.sample_kdtree || !vertex.is_valid())
		return false;
	const int k = std::max(3, s.options.knn_k);
	const uint32 vertex_index = index_of(*s.data.samples, vertex);
	if (vertex_index == INVALID_INDEX)
		return false;
	const Vec3& center = (*s.sample_position)[vertex_index];
	std::vector<std::pair<uint32, Scalar>> neighbors_index;
	s.sample_kdtree->find_nns(center, k + 1, &neighbors_index);
	std::vector<Vec3> neighbors;
	neighbors.reserve(k);
	for (const auto& item : neighbors_index)
	{
		const PVertex neighbor = s.sample_vertices[item.first];
		if (neighbor == vertex)
			continue;
		const uint32 neighbor_index = index_of(*s.data.samples, neighbor);
		neighbors.push_back((*s.sample_position)[neighbor_index]);
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
	const Vec3 current = (*s.sample_normal)[vertex_index];
	if (current.squaredNorm() > Scalar(1e-12) && normal.dot(current) < Scalar(0))
		normal = -normal;
	normal.normalize();
	(*s.sample_normal)[vertex_index] = normal;
	return true;
}

template <typename State>
typename State::SpheresOptimizerType& optimizer_for(State& s)
{
	using SpheresOptimizerType = typename State::SpheresOptimizerType;

	typename SpheresOptimizerType::Data d;
	d.samples_mesh = s.data.samples;
	d.spheres = s.data.spheres;
	d.sample_position = s.sample_position;
	d.sample_area = s.sample_area;
	d.sample_knn = s.sample_knn;
	d.sample_quadric = s.sample_quadric;
	d.sample_line_quadric = s.sample_line_quadric;
	d.sample_ma_position = s.sample_ma_position;
	d.sample_ma_radius = s.sample_ma_radius;
	d.sample_ma_secondary_vertex = s.sample_ma_secondary;
	d.sample_sphere = s.sample_sphere;
	d.sample_error = s.sample_error;
	d.sphere_position = s.sphere_position;
	d.sphere_radius = s.sphere_radius;
	d.sphere_cluster = s.sphere_cluster;
	d.sphere_cluster_area = s.sphere_cluster_area;
	d.sphere_cluster_color = s.sphere_cluster_color;
	d.sphere_neighbors = s.sphere_neighbors;
	d.sphere_error = s.sphere_error;
	d.sphere_error_not_normalized = s.sphere_error_raw;
	d.sample_kdtree = s.sample_kdtree.get();
	d.sample_kdtree_vertices = &s.sample_vertices;
	d.sqem_update_lambda_line_plane = Scalar(s.options.sqem_update_lambda_line_plane);
	if (!s.optimizer)
		s.optimizer = std::make_unique<SpheresOptimizerType>(std::move(d));
	else
		s.optimizer->data() = std::move(d);
	return *s.optimizer;
}

template <typename State>
typename State::TopologyType& topology_for(State& s)
{
	using TopologyType = typename State::TopologyType;

	typename TopologyType::Data d;
	d.skeleton = s.data.skeleton;
	d.spheres_optimizer = &optimizer_for(s);
	d.sphere_skeleton_vertex = s.sphere_skeleton_vertex;
	d.skeleton_position = s.skeleton_position;
	d.skeleton_radius = s.skeleton_radius;
	d.skeleton_source_sphere = s.skeleton_source_sphere;
	d.face_incident_tets = s.skeleton_incident_tets;
	d.edge_degree = s.skeleton_edge_degree;
	if (!s.topology)
		s.topology = std::make_unique<TopologyType>(std::move(d));
	else
		s.topology->data() = std::move(d);
	return *s.topology;
}

template <typename State>
bool evaluate_topology_scores(State& s, const std::vector<Vec3>& points, std::vector<Scalar>& values)
{
	using NeuralModelType = typename State::NeuralModelType;

	s.input.bind(s.data.surface, s.data.input_points, s.input_position, s.input_normal);
	if (neural_model_is(s, NeuralModelType::mf))
	{
		// MF: mf - |sdf| is not a stable medial axis (its minimum is above zero and varies a lot over a shape),
		// so the score is the distance to the medial point of the nearest sample instead of the network value
		if (!s.sample_ma_kdtree)
			return false;
		values.resize(points.size(), Scalar(0));
		detail::parallel_for_index(points.size(), [&](std::size_t i) {
			std::pair<uint32, Scalar> nearest;
			if (!s.sample_ma_kdtree->find_nn(points[i], &nearest) || nearest.first >= s.sample_ma_vertices.size())
				return;
			const uint32 sid = index_of(*s.data.samples, s.sample_ma_vertices[nearest.first]);
			if (sid != INVALID_INDEX)
				values[i] = (points[i] - (*s.sample_ma_position)[sid]).norm();
		});
		return true;
	}
	return s.field.distances(points, values);
}

template <typename State>
typename State::Status status_of(State&, FieldPrepareResult result)
{
	using Status = typename State::Status;

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

template <typename State>
typename State::StatusInfo sample_alpha_level_set(State& s)
{
	using PVertex = typename State::PVertex;
	using RayTraits = typename State::RayTraits;
	using RayParams = typename State::RayParams;
	using Status = typename State::Status;

	if (const typename State::StatusInfo initialized = initialize_data(s); initialized != Status::success)
		return initialized;
	s.input.bind(s.data.surface, s.data.input_points, s.input_position, s.input_normal);
	s.input.invalidate();
	s.input.build_input_kdtree();
	const PrepareInputs<typename State::Surface, typename State::Points> inputs{
		s.data.surface, s.data.input_points, s.input_position, s.input_normal, s.input_knn};
	if (const FieldPrepareResult prepared = s.field.prepare(inputs); prepared != FieldPrepareResult::ok)
		return status_of(s, prepared);

	Vec3 bbox_min(0, 0, 0), bbox_max(1, 1, 1);
	if (const FieldPrepareResult boxed = s.field.sampling_box(inputs, bbox_min, bbox_max);
		boxed != FieldPrepareResult::ok)
		return status_of(s, boxed);
	const Scalar padding = constants::sampling_padding;
	bbox_min -= Vec3(padding, padding, padding);
	bbox_max += Vec3(padding, padding, padding);

	RayTraits traits = s.field.ray_traits(inputs);
	OFFSET_AXIS_PRECONDITION(traits.is_ready(), Status::invalid_data);
	RayParams ray_params;
	ray_params.bbox_expand = constants::ray_box_relative_expand;
	ray_params.alpha = Scalar(s.options.alpha);
	ray_params.tol = Scalar(s.options.tolerance);
	ray_params.step_bound = Scalar(2.0);
	ray_params.batch_size = std::max(1, s.options.ray_sampler_batch_size);
	ray_params.max_iterations = s.options.udf_max_iterations;
	ray_params.min_active_rays = s.options.udf_min_active_rays;
	ray_params.max_outer_iterations = 10;
	ray_params.seed = s.options.seed;
	if (!s.ray_sampler)
		s.ray_sampler = s.field.make_ray_sampler(ray_params);
	else
		s.field.configure_ray_sampler(*s.ray_sampler, ray_params);
	AlphaSamplingParameters params;
	params.sample_radius = Scalar(s.options.sample_radius);
	params.seed = static_cast<uint32>(s.options.seed);
	params.ray_sampler_batch_size = static_cast<size_t>(std::max(1, s.options.ray_sampler_batch_size));
	params.candidates_per_parent = static_cast<size_t>(std::max(1, s.options.alpha_sampling_candidates_per_parent));
	params.parent_batch_size = static_cast<size_t>(std::max(1, s.options.alpha_sampling_parent_batch_size));
	params.bbox_min = bbox_min;
	params.bbox_max = bbox_max;
	auto compute_normals = [&](const std::vector<Vec3>& points, std::vector<Vec3>& normals,
							   std::vector<uint8_t>& keep) {
		keep.assign(points.size(), uint8_t(0));
		normals.assign(points.size(), Vec3(0, 0, 1));
		return s.field.normals(points, normals, keep);
	};
	const AlphaProjectionParameters pp{Scalar(s.options.alpha), Scalar(s.options.tolerance), bbox_min, bbox_max, 1};
	auto project_points = [&](const std::vector<Vec3>& points, AlphaProjectionResult& projection) {
		return s.field.project(points, pp, projection);
	};
	AlphaSamplingResult result =
		offset_axis::sample_alpha_level_set(traits, *s.ray_sampler, params, compute_normals, project_points);
	if (!result.success || result.positions.empty())
		return Status::sampling_failed;
	clear(*s.data.samples);
	s.samples_before_filtering = static_cast<uint32>(result.positions.size());
	for (size_t i = 0; i < result.positions.size(); ++i)
	{
		PVertex v = add_vertex(*s.data.samples);
		const uint32 idx = index_of(*s.data.samples, v);
		(*s.sample_position)[idx] = result.positions[i];
		(*s.sample_normal)[idx] = result.normals[i];
	}
	build_sample_kdtree(s);
	return Status::success;
}

template <typename State>
typename State::StatusInfo apply_sampling_filtering(State& s)
{
	using PVertex = typename State::PVertex;
	using Field = typename State::Field;
	using Status = typename State::Status;

	// attributes come from initialize_data
	OFFSET_AXIS_PRECONDITION(s.initialized, Status::invalid_data);
	if (!s.options.apply_filtering || !s.data.samples || nb_cells<PVertex>(*s.data.samples) == 0)
		return Status::success;
	s.input.bind(s.data.surface, s.data.input_points, s.input_position, s.input_normal);
	s.input.build_input_kdtree();
	if (!s.input.input_kdtree())
		s.input.build_surface_bvh();
	std::vector<Vec3> kept, kept_normals;
	kept.reserve(nb_cells<PVertex>(*s.data.samples));
	foreach_cell(*s.data.samples, [&](PVertex v) {
		const Vec3 p = (*s.sample_position)[index_of(*s.data.samples, v)];
		bool within_distance = false;
		const Scalar max_distance = Scalar(s.options.alpha) + constants::filter_margin;
		if (s.input.input_kdtree())
		{
			std::pair<uint32, Scalar> nearest;
			within_distance = s.input.input_kdtree()->find_nn(p, &nearest, max_distance);
		}
		else if (s.input.surface_bvh())
		{
			std::pair<uint32, Vec3> closest;
			within_distance = s.input.surface_bvh()->closest_point(p, &closest, max_distance);
		}
		if (within_distance)
		{
			kept.push_back(p);
			kept_normals.push_back((*s.sample_normal)[index_of(*s.data.samples, v)]);
		}
		return true;
	});
	const size_t before = nb_cells<PVertex>(*s.data.samples);
	if (kept.size() == before)
		return Status::success;
	clear(*s.data.samples);
	for (size_t i = 0; i < kept.size(); ++i)
	{
		PVertex v = add_vertex(*s.data.samples);
		(*s.sample_position)[index_of(*s.data.samples, v)] = kept[i];
		(*s.sample_normal)[index_of(*s.data.samples, v)] = kept_normals[i];
	}
	s.fitting_data_computed = false;
	if (kept.empty())
	{
		s.sample_kdtree.reset();
		s.sample_vertices.clear();
		return Status::success;
	}
	// neural samples keep their projection normals (the network gradient); for point cloud / surface input the
	// normal is recomputed at the projected position from the input
	if constexpr (Field::recompute_normals_from_input_after_filtering)
		recompute_samples_normals_from_input(s);
	build_sample_kdtree(s);
	return Status::success;
}

template <typename State>
typename State::StatusInfo build_kdtree_and_normals(State& s)
{
	using PVertex = typename State::PVertex;
	using Status = typename State::Status;

	// samples and their kd-tree come from sample_alpha_level_set / apply_sampling_filtering
	OFFSET_AXIS_PRECONDITION(s.data.samples && nb_cells<PVertex>(*s.data.samples) > 0, Status::invalid_data);
	OFFSET_AXIS_PRECONDITION(s.sample_kdtree, Status::invalid_data);
	// PCA for every input: neural samples already carry the network gradient, and for MF models that gradient
	// is not orthogonal to the medial axis.
	if (s.options.recompute_normals_after_sampling)
		parallel_foreach_cell(*s.data.samples, [&](PVertex v) {
			recompute_sample_normal_pca(s, v);
			return true;
		});
	return Status::success;
}

template <typename State>
typename State::StatusInfo compute_fitting_primitives(State& s)
{
	using PVertex = typename State::PVertex;
	using Status = typename State::Status;

	// samples and their kd-tree come from sample_alpha_level_set / apply_sampling_filtering
	OFFSET_AXIS_PRECONDITION(s.data.samples && nb_cells<PVertex>(*s.data.samples) > 0 && s.sample_kdtree, Status::invalid_data);
	compute_samples_area(*s.data.samples, *s.sample_position, *s.sample_normal, *s.sample_area, *s.sample_knn,
						 *s.sample_kdtree, s.sample_vertices, s.options.knn_k, Scalar(s.options.alpha));
	compute_quadrics(*s.data.samples, *s.sample_position, *s.sample_normal, *s.sample_area, *s.sample_knn,
					 *s.sample_quadric, *s.sample_line_quadric, s.options.knn_k);
	s.fitting_data_computed = true;
	return Status::success;
}

template <typename State>
typename State::StatusInfo compute_initial_medial_axis(State& s)
{
	using PVertex = typename State::PVertex;
	using Field = typename State::Field;
	using NeuralModelType = typename State::NeuralModelType;
	using Status = typename State::Status;

	// samples and their kd-tree come from sample_alpha_level_set / apply_sampling_filtering
	OFFSET_AXIS_PRECONDITION(s.data.samples && nb_cells<PVertex>(*s.data.samples) > 0 && s.sample_kdtree, Status::invalid_data);
	const Scalar fallback_radius = Scalar(s.options.alpha);
	const Scalar min_norm = Scalar(1e-12);
	const bool neural_input = Field::is_neural_type && s.field.loaded();
	const bool mf_model = neural_input && neural_model_is(s, NeuralModelType::mf);
	const bool udf_model = neural_input && neural_model_is(s, NeuralModelType::udf);
	const Scalar initial_radius = std::max<Scalar>(fallback_radius * Scalar(10), Scalar(0));
	// repair counters (atomic: run_shrinking_ball_for_vertex runs in parallel); they do not affect the result
	s.ma_repairs = {};
	std::atomic<uint32> degenerate_normals{0}, shrink_fallbacks{0}, secondary_from_nn{0};
	auto run_shrinking_ball_for_vertex = [&](PVertex v) -> bool {
		if (!v.is_valid())
			return false;
		const uint32 idx = index_of(*s.data.samples, v);
		if (idx == INVALID_INDEX)
			return false;
		const Vec3& pt = (*s.sample_position)[idx];
		Vec3 n = (*s.sample_normal)[idx];
		const bool normal_finite = n.allFinite();
		const bool can_shrink = normal_finite && n.squaredNorm() > min_norm;
		if (can_shrink)
			n.normalize();
		else
		{
			n = Vec3(0, 0, 1);
			degenerate_normals.fetch_add(1, std::memory_order_relaxed);
		}
		Vec3 center = pt - n * fallback_radius;
		Scalar radius = fallback_radius;
		PVertex secondary;
		if (can_shrink)
		{
			auto [c1, r1, q1] = geometry::shrinking_ball_center<PVertex>(pt, n, s.sample_kdtree.get(),
																		 s.sample_vertices, initial_radius);
			if (c1.allFinite() && std::isfinite(static_cast<double>(r1)) && r1 > Scalar(0))
			{
				center = c1;
				radius = r1;
				secondary = q1;
			}
			else
				shrink_fallbacks.fetch_add(1, std::memory_order_relaxed);
		}
		if (!center.allFinite())
			center = pt - n * fallback_radius;
		if (!std::isfinite(static_cast<double>(radius)) || radius <= Scalar(0))
			radius = fallback_radius;
		if (!secondary.is_valid() && s.sample_kdtree && !s.sample_vertices.empty())
		{
			std::pair<uint32, Scalar> nn;
			if (s.sample_kdtree->find_nn(center - n * radius, &nn) && nn.first < s.sample_vertices.size())
			{
				secondary = s.sample_vertices[nn.first];
				secondary_from_nn.fetch_add(1, std::memory_order_relaxed);
			}
		}
		(*s.sample_ma_position)[idx] = center;
		(*s.sample_ma_radius)[idx] = std::max(radius, fallback_radius);
		(*s.sample_ma_secondary)[idx] = secondary;
		return true;
	};
	auto run_all = [&] {
		parallel_foreach_cell(*s.data.samples, [&](PVertex v) {
			run_shrinking_ball_for_vertex(v);
			return true;
		});
	};
	run_all();
	if (s.options.ma_flip_prune)
	{
		uint32 flipped = 0, deleted = 0;
		std::vector<uint32> flipped_indices;
		const Scalar udf_threshold =
			udf_model ? constants::udf_flip_threshold_factor * fallback_radius : fallback_radius;
		const Scalar mf_threshold = fallback_radius * std::max(constants::mf_flip_alpha_factor_min,
															   std::min(constants::mf_flip_alpha_factor_max,
																		Scalar(s.options.ma_flip_prune_alpha_factor)));
		const bool supported = mf_model || !Field::is_neural_type || neural_input;
		auto scores = [&](const std::vector<Vec3>& points, std::vector<Scalar>& values,
						  std::vector<Scalar>& sdf) -> bool {
			sdf.clear();
			if (!mf_model)
				return evaluate_topology_scores(s, points, values);
			return s.field.values_sdf(points, values, sdf);
		};
		if (supported)
		{
			std::vector<PVertex> vertices;
			std::vector<Vec3> centers;
			vertices.reserve(nb_cells<PVertex>(*s.data.samples));
			centers.reserve(vertices.capacity());
			foreach_cell(*s.data.samples, [&](PVertex v) {
				const uint32 i = index_of(*s.data.samples, v);
				vertices.push_back(v);
				centers.push_back((*s.sample_ma_position)[i]);
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
					const uint32 i = index_of(*s.data.samples, v);
					if (i == INVALID_INDEX)
						continue;
					Vec3 n = (*s.sample_normal)[i];
					if (n.squaredNorm() > min_norm)
					{
						n.normalize();
						(*s.sample_normal)[i] = -n;
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
						uint32 i = index_of(*s.data.samples, v);
						if (i != INVALID_INDEX)
						{
							active.push_back(v);
							p.push_back((*s.sample_ma_position)[i]);
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
							remove_vertex(*s.data.samples, active[i]);
							++removed;
						}
					return removed;
				};
				deleted += prune(retry);
				if (deleted > 0 && !flipped_indices.empty() && nb_cells<PVertex>(*s.data.samples) > 0)
				{
					build_sample_kdtree(s);
					std::vector<PVertex> survivors;
					survivors.reserve(flipped_indices.size());
					for (uint32 i : flipped_indices)
					{
						PVertex v = of_index<PVertex>(*s.data.samples, i);
						if (v.is_valid())
							survivors.push_back(v);
					}
					for (PVertex v : survivors)
						recompute_sample_normal_pca(s, v);
					for (PVertex v : survivors)
						run_shrinking_ball_for_vertex(v);
					deleted += prune(survivors);
				}
			}
		}
		s.ma_repairs.ma_flip_flipped = flipped;
		s.ma_repairs.ma_flip_deleted = deleted;
	}
	s.ma_repairs.ma_degenerate_normals = degenerate_normals.load();
	s.ma_repairs.ma_shrink_fallbacks = shrink_fallbacks.load();
	s.ma_repairs.ma_secondary_from_nn = secondary_from_nn.load();
	build_sample_kdtree(s);
	return Status::success;
}

template <typename State>
typename State::StatusInfo initialize_spheres(State& s)
{
	using Status = typename State::Status;

	OFFSET_AXIS_PRECONDITION(s.fitting_data_computed && s.data.samples && s.data.spheres, Status::invalid_data);
	if (!optimizer_for(s).initialize_from_samples(Scalar(s.options.init_dilation_constant)))
		return Status::sphere_initialization_failed;
	return Status::success;
}

template <typename State>
typename State::StatusInfo optimize_spheres(State& s, const typename State::IterationCallback& callback)
{
	using PVertex = typename State::PVertex;
	using SpheresOptimizerType = typename State::SpheresOptimizerType;
	using Status = typename State::Status;

	// spheres come from initialize_spheres
	OFFSET_AXIS_PRECONDITION(s.data.spheres && nb_cells<PVertex>(*s.data.spheres) > 0, Status::invalid_data);
	auto& optimizer = optimizer_for(s);
	optimizer.reset_optimization();
	typename SpheresOptimizerType::Status status = SpheresOptimizerType::Status::running;
	while (status == SpheresOptimizerType::Status::running)
	{
		status = update_spheres_once(s, Scalar(s.options.sqem_update_lambda_line_plane));
		if (callback)
			callback(optimizer.metrics());
	}
	return status == SpheresOptimizerType::Status::failed ? Status::optimization_failed : Status::success;
}

template <typename State>
typename State::StatusInfo build_skeleton(State& s)
{
	using PVertex = typename State::PVertex;
	using TopologyType = typename State::TopologyType;
	using Status = typename State::Status;

	// spheres come from initialize_spheres / optimize_spheres
	OFFSET_AXIS_PRECONDITION(s.data.spheres && nb_cells<PVertex>(*s.data.spheres) > 0, Status::invalid_data);
	return topology_for(s).build_from_spheres() == TopologyType::Status::success ? Status::success
																				 : Status::skeleton_build_failed;
}

template <typename State>
typename State::StatusInfo fix_topology(State& s)
{
	using NMVertex = typename State::NMVertex;
	using TopologyType = typename State::TopologyType;
	using Status = typename State::Status;

	// skeleton comes from build_skeleton
	OFFSET_AXIS_PRECONDITION(s.data.skeleton && nb_cells<NMVertex>(*s.data.skeleton) > 0, Status::invalid_data);
	return topology_for(s).fix_topology([&](const std::vector<Vec3>& p, std::vector<Scalar>& v) {
		return evaluate_topology_scores(s, p, v);
	}) == TopologyType::Status::success
			   ? Status::success
			   : Status::topology_fix_failed;
}

template <typename State>
typename State::StatusInfo prune_residual_sheets(State& s)
{
	using NMVertex = typename State::NMVertex;
	using TopologyType = typename State::TopologyType;
	using Status = typename State::Status;

	// skeleton comes from build_skeleton
	OFFSET_AXIS_PRECONDITION(s.data.skeleton && nb_cells<NMVertex>(*s.data.skeleton) > 0, Status::invalid_data);
	topology_for(s).set_residual_prune_mode(
		static_cast<typename TopologyType::ResidualPruneMode>(static_cast<int>(s.options.residual_prune_mode)));
	return topology_for(s).prune_residual_sheets() == TopologyType::Status::success ? Status::success
																					: Status::topology_fix_failed;
}

template <typename State>
typename State::StatusInfo sample(State& s)
{
	using Status = typename State::Status;

	const typename State::StatusInfo status = sample_alpha_level_set(s);
	if (status != Status::success)
		return status;
	return apply_sampling_filtering(s);
}

template <typename State>
typename State::StatusInfo compute_fitting_data(State& s)
{
	using Status = typename State::Status;

	typename State::StatusInfo status = build_kdtree_and_normals(s);
	if (status != Status::success)
		return status;
	status = compute_fitting_primitives(s);
	if (status != Status::success)
		return status;
	return compute_initial_medial_axis(s);
}

// One optimizer iteration; the caller passes lambda explicitly so that a UI thread can feed a live value.
template <typename State>
typename State::SpheresOptimizerType::Status update_spheres_once(State& s, Scalar sqem_update_lambda_line_plane)
{
	using PVertex = typename State::PVertex;
	using SpheresOptimizerType = typename State::SpheresOptimizerType;

	// spheres come from initialize_spheres
	if (!s.data.spheres || nb_cells<PVertex>(*s.data.spheres) == 0)
		return SpheresOptimizerType::Status::failed;
	return optimizer_for(s).update_once(sqem_update_lambda_line_plane);
}

template <typename State>
typename State::StatusInfo fix_topology_and_prune(State& s)
{
	using Status = typename State::Status;

	const typename State::StatusInfo status = fix_topology(s);
	if (status != Status::success || !s.options.residual_prune)
		return status;
	return prune_residual_sheets(s);
}

template <typename State>
typename State::Result reconstruct(State& s, const typename State::IterationCallback& callback)
{
	using Status = typename State::Status;
	using Result = typename State::Result;

	Result r;
	s.input.invalidate();
	if (s.data.samples)
		clear(*s.data.samples);
	if (s.data.spheres)
		clear(*s.data.spheres);
	if (s.data.skeleton)
		clear(*s.data.skeleton);
	s.sample_kdtree.reset();
	s.sample_ma_kdtree.reset();
	s.sample_vertices.clear();
	s.sample_ma_vertices.clear();
	s.optimizer.reset();
	s.topology.reset();
	s.fitting_data_computed = false;
	s.samples_before_filtering = 0;
	s.ma_repairs = {};
	auto run_stage = [&r](float64& time, auto&& stage) -> Status {
		const auto start = std::chrono::high_resolution_clock::now();
		const typename State::StatusInfo status = stage();
		time = std::chrono::duration<float64, std::milli>(std::chrono::high_resolution_clock::now() - start).count();
		if (status.code != Status::success)
			r.message = status.message;
		return status.code;
	};
	if ((r.status = run_stage(r.timing.preprocessing_ms, [&] { return initialize_data(s); })) != Status::success)
		return r;
	if ((r.status = run_stage(r.timing.sampling_ms, [&] { return sample_alpha_level_set(s); })) != Status::success)
		return r;
	r.counts.samples_before_filtering = s.samples_before_filtering;
	if ((r.status = run_stage(r.timing.sample_filtering_ms, [&] { return apply_sampling_filtering(s); })) !=
		Status::success)
		return r;
	if ((r.status = run_stage(r.timing.kdtree_and_normals_ms, [&] { return build_kdtree_and_normals(s); })) !=
		Status::success)
		return r;
	if ((r.status = run_stage(r.timing.fitting_primitives_ms, [&] { return compute_fitting_primitives(s); })) !=
		Status::success)
		return r;
	if ((r.status = run_stage(r.timing.initial_medial_axis_ms, [&] { return compute_initial_medial_axis(s); })) !=
		Status::success)
		return r;
	if ((r.status = run_stage(r.timing.sphere_initialization_ms, [&] { return initialize_spheres(s); })) !=
		Status::success)
		return r;
	if ((r.status = run_stage(r.timing.optimization_ms, [&] { return optimize_spheres(s, callback); })) !=
		Status::success)
		return r;
	if (s.optimizer)
	{
		r.timing.cluster_ms = s.optimizer->metrics().cluster_total_ms;
		r.timing.sphere_update_ms = s.optimizer->metrics().sphere_update_total_ms;
		r.timing.error_ms = s.optimizer->metrics().error_total_ms;
	}
	if ((r.status = run_stage(r.timing.skeleton_construction_ms, [&] { return build_skeleton(s); })) != Status::success)
		return r;
	if ((r.status = run_stage(r.timing.topology_processing_ms, [&] { return fix_topology(s); })) != Status::success)
		return r;
	if (s.options.residual_prune)
	{
		float64 residual_ms = 0.0;
		if ((r.status = run_stage(residual_ms, [&] { return prune_residual_sheets(s); })) != Status::success)
			return r;
		r.timing.topology_processing_ms += residual_ms;
	}
	r.status = Status::success;
	r.counts = counts(s);
	r.repairs = repairs(s);
	return r;
}

} // namespace stages

} // namespace offset_axis

} // namespace geometry

} // namespace cgogn

#undef OFFSET_AXIS_PRECONDITION

#endif // CGOGN_GEOMETRY_ALGOS_OFFSET_AXIS_RECONSTRUCTION_STAGES_H_
