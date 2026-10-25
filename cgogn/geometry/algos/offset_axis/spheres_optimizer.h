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
 * along with this library; if not, write to the Free Software Foundation,     *
 * Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301 USA.           *
 *                                                                              *
 * Web site: http://cgogn.unistra.fr/                                           *
 * Contact information: cgogn@unistra.fr                                        *
 *                                                                              *
 *******************************************************************************/
#ifndef CGOGN_GEOMETRY_ALGOS_OFFSET_AXIS_SPHERES_OPTIMIZER_H_
#define CGOGN_GEOMETRY_ALGOS_OFFSET_AXIS_SPHERES_OPTIMIZER_H_

#include <cgogn/core/functions/attributes.h>
#include <cgogn/core/functions/mesh_info.h>
#include <cgogn/core/types/cell_marker.h>
#include <cgogn/geometry/types/line_quadric.h>
#include <cgogn/geometry/types/spherical_quadric.h>
#include <cgogn/geometry/types/vector_traits.h>
#include <libacc/kd_tree.h>

#include <algorithm>
#include <atomic>
#include <chrono>
#include <cmath>
#include <cstdlib>
#include <limits>
#include <memory>
#include <set>
#include <unordered_set>
#include <utility>
#include <vector>

namespace cgogn
{

namespace geometry
{

namespace offset_axis
{

// sorted, duplicate-free vector used as a set: same element order and membership as std::set
template <typename T>
void sorted_vector_insert(std::vector<T>& v, const T& item)
{
	auto it = std::lower_bound(v.begin(), v.end(), item);
	if (it == v.end() || item < *it)
		v.insert(it, item);
}

template <typename T>
void sorted_vector_erase(std::vector<T>& v, const T& item)
{
	auto it = std::lower_bound(v.begin(), v.end(), item);
	if (it != v.end() && !(item < *it))
		v.erase(it);
}

template <typename POINTS>
class SpheresOptimizer
{
public:
	using Vertex = typename mesh_traits<POINTS>::Vertex;
	template <typename T>
	using Attribute = typename mesh_traits<POINTS>::template Attribute<T>;

	struct Data
	{
		POINTS* samples_mesh = nullptr;
		POINTS* spheres = nullptr;
		std::shared_ptr<Attribute<Vec3>> sample_position;
		std::shared_ptr<Attribute<Scalar>> sample_area;
		std::shared_ptr<Attribute<std::vector<Vertex>>> sample_knn;
		std::shared_ptr<Attribute<Spherical_Quadric>> sample_quadric;
		std::shared_ptr<Attribute<Line_Quadric>> sample_line_quadric;
		std::shared_ptr<Attribute<Vec3>> sample_ma_position;
		std::shared_ptr<Attribute<Scalar>> sample_ma_radius;
		std::shared_ptr<Attribute<Vertex>> sample_ma_secondary_vertex;
		std::shared_ptr<Attribute<Vertex>> sample_sphere;
		std::shared_ptr<Attribute<Scalar>> sample_error;
		std::shared_ptr<Attribute<Vec3>> sphere_position;
		std::shared_ptr<Attribute<Scalar>> sphere_radius;
		std::shared_ptr<Attribute<std::vector<Vertex>>> sphere_cluster;
		std::shared_ptr<Attribute<Scalar>> sphere_cluster_area;
		std::shared_ptr<Attribute<Vec4>> sphere_cluster_color;
		std::shared_ptr<Attribute<std::vector<Vertex>>> sphere_neighbors;
		std::shared_ptr<Attribute<Scalar>> sphere_error;
		std::shared_ptr<Attribute<Scalar>> sphere_error_not_normalized;
		const acc::KDTree<3, uint32>* sample_kdtree = nullptr;
		const std::vector<Vertex>* sample_kdtree_vertices = nullptr;
		Scalar sqem_update_lambda_line_plane = Scalar(0.20);
	};

	struct Metrics
	{
		uint32 sphere_count = 0;
		uint32 iteration = 0;
		Scalar total_error = Scalar(0);
		Scalar total_error_not_normalized = Scalar(0);
		Scalar last_total_error = std::numeric_limits<Scalar>::max();
		Scalar error_difference = Scalar(0);
		Scalar minimum_error = Scalar(0);
		Scalar maximum_error = Scalar(0);
		float64 cluster_total_ms = 0.0;
		float64 sphere_update_total_ms = 0.0;
		float64 error_total_ms = 0.0;
		// sphere initialization split: cover loop, Euclidean owner assignment, first cluster pass
		float64 init_cover_ms = 0.0;
		float64 init_owner_ms = 0.0;
		float64 init_cluster_ms = 0.0;
		bool sphere_topology_changed = false;
		// sphere updates skipped because the solve gave non-finite values (sum over the iterations)
		uint32 sphere_update_nonfinite = 0;
	};

	enum class Status
	{
		running,
		converged,
		max_iterations,
		failed
	};

	explicit SpheresOptimizer(Data data) : data_(std::move(data))
	{
		sync_sphere_count();
	}

	SpheresOptimizer(const SpheresOptimizer&) = delete;
	SpheresOptimizer& operator=(const SpheresOptimizer&) = delete;
	SpheresOptimizer(SpheresOptimizer&&) = delete;
	SpheresOptimizer& operator=(SpheresOptimizer&&) = delete;

	Data& data() { return data_; }
	const Data& data() const { return data_; }
	const Metrics& metrics() const { return metrics_; }

	bool initialize_from_samples(Scalar init_dilation_constant, bool owner_euclidean = false, uint32 owner_knn = 16)
	{
		const auto init_start = std::chrono::steady_clock::now();
		if (!ready() || !data_.sample_position || !data_.sample_knn ||
			!data_.sample_ma_position || !data_.sample_ma_radius || !data_.sample_ma_secondary_vertex ||
			!data_.sphere_position || !data_.sphere_radius || !data_.sphere_cluster_color)
			return false;

		clear(*data_.spheres);
		constexpr uint32 min_cover_points = 10;
		constexpr uint32 max_nb_spheres = 100000;
		std::vector<Vertex> sorted_vertices;
		uint32 max_sample_index = 0;
		foreach_cell(*data_.samples_mesh, [&](Vertex v) {
			sorted_vertices.push_back(v);
			const uint32 idx = index_of(*data_.samples_mesh, v);
			if (idx != INVALID_INDEX)
				max_sample_index = std::max(max_sample_index, idx);
			return true;
		});
		// same comparison outcomes as sorting the vertices directly, so ties end up in the same order
		std::vector<std::pair<Scalar, Vertex>> keyed_vertices;
		keyed_vertices.reserve(sorted_vertices.size());
		for (Vertex v : sorted_vertices)
		{
			const uint32 idx = index_of(*data_.samples_mesh, v);
			const Scalar r = idx != INVALID_INDEX ? (*data_.sample_ma_radius)[idx] : Scalar(-1);
			keyed_vertices.emplace_back(std::isfinite(static_cast<double>(r)) ? r : Scalar(-1), v);
		}
		std::sort(keyed_vertices.begin(), keyed_vertices.end(),
				  [](const std::pair<Scalar, Vertex>& a, const std::pair<Scalar, Vertex>& b) { return a.first > b.first; });
		for (std::size_t i = 0; i < keyed_vertices.size(); ++i)
			sorted_vertices[i] = keyed_vertices[i].second;
		auto covered = get_or_add_attribute<bool, Vertex>(*data_.samples_mesh, "__covered");
		covered->fill(false);
		const std::size_t marks_size = sorted_vertices.empty() ? 0 : static_cast<std::size_t>(max_sample_index) + 1;
		std::vector<uint32> candidate_marks(marks_size, 0);
		uint32 candidate_mark_token = 1;
		uint32 sphere_count = 0;
		std::vector<Vertex> candidate_cover;
		std::vector<Vertex> stack;
		for (Vertex v : sorted_vertices)
		{
			const uint32 v_index = index_of(*data_.samples_mesh, v);
			if (v_index == INVALID_INDEX || sphere_count >= max_nb_spheres || (*covered)[v_index])
				continue;
			const Vec3& vp = (*data_.sample_ma_position)[v_index];
			const Scalar vr = (*data_.sample_ma_radius)[v_index];
			if (!vp.allFinite() || !std::isfinite(static_cast<double>(vr)) || vr <= Scalar(0))
				continue;
			const Scalar dilation_radius = std::max<Scalar>(vr + init_dilation_constant, Scalar(0));
			const Scalar dilation_radius_sq = dilation_radius * dilation_radius;
			if (candidate_mark_token == std::numeric_limits<uint32>::max())
			{
				std::fill(candidate_marks.begin(), candidate_marks.end(), 0u);
				candidate_mark_token = 1;
			}
			const uint32 current_mark = candidate_mark_token++;
			candidate_cover.clear();
			auto flood_cover = [&](Vertex seed) {
				if (!seed.is_valid())
					return;
				const uint32 seed_idx = index_of(*data_.samples_mesh, seed);
				if (seed_idx == INVALID_INDEX || seed_idx >= candidate_marks.size() || (*covered)[seed_idx] ||
					candidate_marks[seed_idx] == current_mark)
					return;
				stack.clear();
				candidate_marks[seed_idx] = current_mark;
				candidate_cover.push_back(seed);
				stack.push_back(seed);
				while (!stack.empty())
				{
					const Vertex w = stack.back();
					stack.pop_back();
					const uint32 w_idx = index_of(*data_.samples_mesh, w);
					if (w_idx == INVALID_INDEX)
						continue;
					for (Vertex u : (*data_.sample_knn)[w_idx])
					{
						const uint32 u_idx = index_of(*data_.samples_mesh, u);
						if (u_idx == INVALID_INDEX || u_idx >= candidate_marks.size())
							continue;
						if (!(*covered)[u_idx] && candidate_marks[u_idx] != current_mark &&
							((*data_.sample_position)[u_idx] - vp).squaredNorm() < dilation_radius_sq)
						{
							candidate_marks[u_idx] = current_mark;
							candidate_cover.push_back(u);
							stack.push_back(u);
						}
					}
				}
			};
			flood_cover(v);
			flood_cover((*data_.sample_ma_secondary_vertex)[v_index]);
			for (Vertex covered_vertex : candidate_cover)
			{
				const uint32 covered_index = index_of(*data_.samples_mesh, covered_vertex);
				if (covered_index != INVALID_INDEX)
					(*covered)[covered_index] = true;
			}
			if (candidate_cover.size() < min_cover_points && (sphere_count != 0 || candidate_cover.empty()))
				continue;
			const Vertex sphere = add_vertex(*data_.spheres);
			const uint32 sphere_index = index_of(*data_.spheres, sphere);
			++sphere_count;
			(*data_.sphere_position)[sphere_index] = vp;
			(*data_.sphere_radius)[sphere_index] = vr;
			(*data_.sphere_cluster_color)[sphere_index] =
				Vec4(0.5 + 0.5 * (rand() % 256) / 256.0, 0.5 + 0.5 * (rand() % 256) / 256.0,
					 0.5 + 0.5 * (rand() % 256) / 256.0, 1.0);
		}
		remove_attribute<Vertex>(*data_.samples_mesh, covered);
		if (data_.sample_sphere)
			data_.sample_sphere->fill(Vertex());
		metrics_ = Metrics{};
		metrics_.sphere_count = sphere_count;
		const auto cover_end = std::chrono::steady_clock::now();
		metrics_.init_cover_ms = std::chrono::duration<float64, std::milli>(cover_end - init_start).count();
		auto owner_end = cover_end;
		if (owner_euclidean && sphere_count > 0)
		{
			assign_owners_euclidean(owner_knn);
			owner_end = std::chrono::steady_clock::now();
		}
		metrics_.init_owner_ms = std::chrono::duration<float64, std::milli>(owner_end - cover_end).count();
		compute_clusters_local();
		const auto cluster_end = std::chrono::steady_clock::now();
		metrics_.init_cluster_ms = std::chrono::duration<float64, std::milli>(cluster_end - owner_end).count();
		return true;
	}

	void reset_optimization()
	{
		const Metrics previous = metrics_;
		metrics_ = Metrics{};
		metrics_.sphere_count = previous.sphere_count;
		metrics_.init_cover_ms = previous.init_cover_ms;
		metrics_.init_owner_ms = previous.init_owner_ms;
		metrics_.init_cluster_ms = previous.init_cluster_ms;
		convergence_reached_ = false;
		post_convergence_iterations_ = 0;
	}

	Status update_once(Scalar sqem_update_lambda_line_plane)
	{
		if (!ready())
			return Status::failed;
		data_.sqem_update_lambda_line_plane = sqem_update_lambda_line_plane;
		metrics_.sphere_topology_changed = false;
		sync_sphere_count();
		const auto cluster_start = std::chrono::high_resolution_clock::now();
		if (metrics_.iteration % 10 == 0)
			compute_sphere_neighbors();
		compute_clusters_local();
		prune_empty_clusters();
		const auto cluster_end = std::chrono::high_resolution_clock::now();
		metrics_.cluster_total_ms += std::chrono::duration<float64, std::milli>(cluster_end - cluster_start).count();

		const auto update_start = std::chrono::high_resolution_clock::now();
		parallel_foreach_cell(*data_.spheres, [&](Vertex v) {
			update_sphere_line_quadric_distance_free_radius(v);
			return true;
		});
		metrics_.sphere_update_nonfinite += nonfinite_updates_.exchange(0);
		const auto update_end = std::chrono::high_resolution_clock::now();
		metrics_.sphere_update_total_ms += std::chrono::duration<float64, std::milli>(update_end - update_start).count();

		const auto error_start = std::chrono::high_resolution_clock::now();
		compute_spheres_error();
		const auto error_end = std::chrono::high_resolution_clock::now();
		metrics_.error_total_ms += std::chrono::duration<float64, std::milli>(error_end - error_start).count();
		++metrics_.iteration;
		return update_stop_status();
	}

	void compute_sphere_neighbors()
	{
		if (!data_.samples_mesh || !data_.spheres || !data_.sample_knn || !data_.sample_sphere || !data_.sphere_neighbors)
			return;
		parallel_foreach_cell(*data_.spheres, [&](Vertex v) {
			(*data_.sphere_neighbors)[index_of(*data_.spheres, v)].clear();
			return true;
		});
		foreach_cell(*data_.samples_mesh, [&](Vertex v) {
			const uint32 v_index = index_of(*data_.samples_mesh, v);
			const Vertex v_sphere = (*data_.sample_sphere)[v_index];
			for (Vertex w : (*data_.sample_knn)[v_index])
			{
				const Vertex w_sphere = (*data_.sample_sphere)[index_of(*data_.samples_mesh, w)];
				if (v_sphere.is_valid() && w_sphere.is_valid() && v_sphere != w_sphere)
				{
					(*data_.sphere_neighbors)[index_of(*data_.spheres, v_sphere)].push_back(w_sphere);
					(*data_.sphere_neighbors)[index_of(*data_.spheres, w_sphere)].push_back(v_sphere);
				}
			}
			return true;
		});
		parallel_foreach_cell(*data_.spheres, [&](Vertex v) {
			std::vector<Vertex>& neighbors = (*data_.sphere_neighbors)[index_of(*data_.spheres, v)];
			std::sort(neighbors.begin(), neighbors.end());
			neighbors.erase(std::unique(neighbors.begin(), neighbors.end()), neighbors.end());
			return true;
		});
	}

	void compute_clusters_local()
	{
		if (!ready())
			return;
		if (metrics_.sphere_count == 0)
		{
			data_.sample_sphere->fill(Vertex());
			return;
		}
		parallel_foreach_cell(*data_.spheres, [&](Vertex v) {
			const uint32 idx = index_of(*data_.spheres, v);
			(*data_.sphere_cluster)[idx].clear();
			(*data_.sphere_cluster_area)[idx] = Scalar(0);
			return true;
		});
		parallel_foreach_cell(*data_.samples_mesh, [&](Vertex v) {
			const uint32 v_index = index_of(*data_.samples_mesh, v);
			const Vertex owner = (*data_.sample_sphere)[v_index];
			const uint32 owner_index = owner.is_valid() ? index_of(*data_.spheres, owner) : INVALID_INDEX;
			Scalar min_distance = std::numeric_limits<Scalar>::max();
			Vertex closest_sphere;
			uint32 closest_index = INVALID_INDEX;
			auto evaluate_candidate = [&](Vertex candidate) {
				const uint32 idx = index_of(*data_.spheres, candidate);
				if (!is_valid_sphere_for_clustering(idx))
				{
					return;
				}
				const Vec3& center = (*data_.sphere_position)[idx];
				const Scalar radius = (*data_.sphere_radius)[idx];
				const Scalar distance = (*data_.sample_quadric)[v_index].eval(Vec4(center.x(), center.y(), center.z(), radius)) +
					data_.sqem_update_lambda_line_plane * (*data_.sample_line_quadric)[v_index].eval(center);
				if (!std::isfinite(static_cast<double>(distance)))
				{
					return;
				}
				if (distance < min_distance)
				{
					min_distance = distance;
					closest_sphere = candidate;
					closest_index = idx;
				}
			};
			if (is_valid_sphere_for_clustering(owner_index))
			{
				const std::vector<Vertex>& neighbors = (*data_.sphere_neighbors)[owner_index];
				bool owner_evaluated = false;
				for (Vertex candidate : neighbors)
				{
					if (!owner_evaluated && owner < candidate)
					{
						evaluate_candidate(owner);
						owner_evaluated = true;
					}
					if (candidate == owner)
						owner_evaluated = true;
					evaluate_candidate(candidate);
				}
				if (!owner_evaluated)
					evaluate_candidate(owner);
			}
			else
				foreach_cell(*data_.spheres, [&](Vertex candidate) { evaluate_candidate(candidate); return true; });
			(*data_.sample_sphere)[v_index] = closest_index == INVALID_INDEX ? Vertex() : closest_sphere;
			return true;
		});
		// clusters are filled in sample order, so the sums over a cluster do not depend on thread scheduling
		foreach_cell(*data_.samples_mesh, [&](Vertex v) {
			const uint32 v_index = index_of(*data_.samples_mesh, v);
			const Vertex sphere = (*data_.sample_sphere)[v_index];
			if (!sphere.is_valid())
				return true;
			const uint32 sphere_index = index_of(*data_.spheres, sphere);
			(*data_.sphere_cluster)[sphere_index].push_back(v);
			(*data_.sphere_cluster_area)[sphere_index] += (*data_.sample_area)[v_index];
			return true;
		});
	}

	void compute_spheres_error()
	{
		if (!ready())
			return;
		parallel_foreach_cell(*data_.spheres, [&](Vertex v) {
			const uint32 idx = index_of(*data_.spheres, v);
			const Vec3& center = (*data_.sphere_position)[idx];
			const Scalar radius = (*data_.sphere_radius)[idx];
			const std::vector<Vertex>& cluster = (*data_.sphere_cluster)[idx];
			Scalar cluster_error = Scalar(0);
			for (Vertex sample : cluster)
			{
				const uint32 sample_idx = index_of(*data_.samples_mesh, sample);
				const Scalar distance = (*data_.sample_quadric)[sample_idx].eval(Vec4(center.x(), center.y(), center.z(), radius)) +
					data_.sqem_update_lambda_line_plane * (*data_.sample_line_quadric)[sample_idx].eval(center);
				if (data_.sample_error)
					(*data_.sample_error)[sample_idx] = distance;
				cluster_error += distance;
			}
			(*data_.sphere_error)[idx] = (*data_.sphere_cluster_area)[idx] > Scalar(0)
				? cluster_error / (*data_.sphere_cluster_area)[idx] : Scalar(0);
			(*data_.sphere_error_not_normalized)[idx] = cluster_error;
			return true;
		});
		metrics_.minimum_error = std::numeric_limits<Scalar>::max();
		metrics_.maximum_error = Scalar(0);
		metrics_.total_error = Scalar(0);
		metrics_.total_error_not_normalized = Scalar(0);
		foreach_cell(*data_.spheres, [&](Vertex v) {
			const uint32 idx = index_of(*data_.spheres, v);
			const Scalar error = (*data_.sphere_error)[idx];
			const Scalar unnormalized = (*data_.sphere_error_not_normalized)[idx];
			if (error < metrics_.minimum_error)
				metrics_.minimum_error = error;
			if (error > metrics_.maximum_error)
				metrics_.maximum_error = error;
			metrics_.total_error += error;
			metrics_.total_error_not_normalized += unnormalized;
			return true;
		});
		metrics_.error_difference = std::abs(metrics_.total_error - metrics_.last_total_error);
		metrics_.last_total_error = metrics_.total_error;
		sync_sphere_count();
	}

	void remove_sphere(Vertex sphere)
	{
		if (!data_.samples_mesh || !data_.spheres || !data_.sample_sphere || !data_.sphere_cluster ||
			!sphere.is_valid() || sphere.dart_.index_ >= data_.spheres->darts_.maximum_index())
			return;
		const uint32 sphere_index = index_of(*data_.spheres, sphere);
		if (sphere_index == INVALID_INDEX)
			return;
		const std::vector<Vertex> removed_cluster = (*data_.sphere_cluster)[sphere_index];
		const std::vector<Vertex> removed_neighbors =
			data_.sphere_neighbors ? (*data_.sphere_neighbors)[sphere_index] : std::vector<Vertex>();
		for (Vertex sample : removed_cluster)
		{
			const uint32 idx = index_of(*data_.samples_mesh, sample);
			if (idx != INVALID_INDEX)
				(*data_.sample_sphere)[idx] = Vertex();
		}
		if (data_.sphere_neighbors)
		{
			for (Vertex neighbor : (*data_.sphere_neighbors)[sphere_index])
			{
				const uint32 idx = neighbor.is_valid() ? index_of(*data_.spheres, neighbor) : INVALID_INDEX;
				if (idx != INVALID_INDEX)
					sorted_vector_erase((*data_.sphere_neighbors)[idx], sphere);
			}
			(*data_.sphere_neighbors)[sphere_index].clear();
		}
		remove_vertex(*data_.spheres, sphere);
		redistribute_removed_sphere_cluster(removed_cluster, removed_neighbors);
		metrics_.sphere_topology_changed = true;
		sync_sphere_count();
	}

	Vertex add_sphere(const Vec3& position, Scalar radius)
	{
		if (!data_.spheres || !data_.sphere_position || !data_.sphere_radius)
			return Vertex();

		const Vertex sphere = add_vertex(*data_.spheres);
		const uint32 sphere_id = index_of(*data_.spheres, sphere);
		if (sphere_id == INVALID_INDEX)
			return Vertex();

		(*data_.sphere_position)[sphere_id] = position;
		(*data_.sphere_radius)[sphere_id] = radius;
		if (data_.sphere_cluster)
			(*data_.sphere_cluster)[sphere_id].clear();
		if (data_.sphere_cluster_area)
			(*data_.sphere_cluster_area)[sphere_id] = Scalar(0);
		if (data_.sphere_neighbors)
			(*data_.sphere_neighbors)[sphere_id].clear();
		if (data_.sphere_error)
			(*data_.sphere_error)[sphere_id] = Scalar(0);
		if (data_.sphere_error_not_normalized)
			(*data_.sphere_error_not_normalized)[sphere_id] = Scalar(0);

		metrics_.sphere_topology_changed = true;
		sync_sphere_count();
		return sphere;
	}

private:
	// Initial owner of every sample: the sphere minimizing |p - c| - r among the knn nearest centers.
	void assign_owners_euclidean(uint32 owner_knn)
	{
		std::vector<Vec3> centers;
		std::vector<Vertex> center_spheres;
		foreach_cell(*data_.spheres, [&](Vertex s) {
			const uint32 idx = index_of(*data_.spheres, s);
			if (is_valid_sphere_for_clustering(idx))
			{
				centers.push_back((*data_.sphere_position)[idx]);
				center_spheres.push_back(s);
			}
			return true;
		});
		if (centers.empty())
			return;
		const acc::KDTree<3, uint32> tree(centers);
		const std::size_t k = std::min<std::size_t>(std::max<uint32>(owner_knn, 1u), centers.size());
		parallel_foreach_cell(*data_.samples_mesh, [&](Vertex v) {
			const uint32 v_index = index_of(*data_.samples_mesh, v);
			const Vec3& p = (*data_.sample_position)[v_index];
			std::vector<std::pair<uint32, double>> nearest;
			tree.find_nns(p, k, &nearest);
			Scalar best = std::numeric_limits<Scalar>::max();
			Vertex best_sphere;
			for (const auto& item : nearest)
			{
				const Vertex s = center_spheres[item.first];
				const uint32 idx = index_of(*data_.spheres, s);
				const Scalar gap = ((*data_.sphere_position)[idx] - p).norm() - (*data_.sphere_radius)[idx];
				if (std::isfinite(static_cast<double>(gap)) && gap < best)
				{
					best = gap;
					best_sphere = s;
				}
			}
			(*data_.sample_sphere)[v_index] = best_sphere;
			return true;
		});
		compute_sphere_neighbors();
	}

	bool is_valid_sphere_for_clustering(uint32 sphere_index) const
	{
		if (!data_.sphere_position || !data_.sphere_radius || sphere_index == INVALID_INDEX)
			return false;
		const Vec3& center = (*data_.sphere_position)[sphere_index];
		const Scalar radius = (*data_.sphere_radius)[sphere_index];
		return center.allFinite() && std::isfinite(static_cast<double>(radius)) && radius > Scalar(0);
	}

	bool ready() const
	{
		return data_.samples_mesh && data_.spheres && data_.sample_position && data_.sample_area && data_.sample_quadric &&
			data_.sample_line_quadric && data_.sample_sphere && data_.sphere_position && data_.sphere_radius &&
			data_.sphere_cluster && data_.sphere_cluster_area && data_.sphere_neighbors && data_.sphere_error &&
			data_.sphere_error_not_normalized;
	}

	void sync_sphere_count()
	{
		metrics_.sphere_count = data_.spheres ? nb_cells<Vertex>(*data_.spheres) : 0;
	}

	void prune_empty_clusters()
	{
		std::vector<Vertex> to_remove;
		foreach_cell(*data_.spheres, [&](Vertex sphere) {
			if ((*data_.sphere_cluster)[index_of(*data_.spheres, sphere)].empty())
				to_remove.push_back(sphere);
			return true;
		});
		for (Vertex sphere : to_remove)
			remove_sphere(sphere);
	}

	bool try_get_nearest_sample_ma_radius(const Vec3& query, Scalar& out_radius) const
	{
		if (!data_.samples_mesh || !data_.sample_ma_radius || !data_.sample_kdtree || !data_.sample_kdtree_vertices)
			return false;
		std::pair<uint32, Scalar> result;
		if (!data_.sample_kdtree->find_nn(query, &result) || result.first >= data_.sample_kdtree_vertices->size())
			return false;
		const uint32 sample_index = index_of(*data_.samples_mesh, (*data_.sample_kdtree_vertices)[result.first]);
		if (sample_index == INVALID_INDEX)
			return false;
		const Scalar radius = (*data_.sample_ma_radius)[sample_index];
		if (!std::isfinite(radius) || radius <= Scalar(0))
			return false;
		out_radius = radius;
		return true;
	}

	void update_sphere_line_quadric_distance_fix_current_radius(Vertex sphere, Scalar fixed_radius)
	{
		if (!std::isfinite(fixed_radius) || fixed_radius <= Scalar(0))
			return;
		const uint32 sphere_index = index_of(*data_.spheres, sphere);
		const std::vector<Vertex>& cluster = (*data_.sphere_cluster)[sphere_index];
		if (cluster.empty())
			return;
		Vec3 center = (*data_.sphere_position)[sphere_index];
		Spherical_Quadric q;
		Line_Quadric lq;
		for (Vertex sample : cluster)
		{
			const uint32 idx = index_of(*data_.samples_mesh, sample);
			const Scalar weight = value<Scalar>(*data_.samples_mesh, data_.sample_area, sample);
			if (weight <= Scalar(0))
				continue;
			q += (*data_.sample_quadric)[idx] * weight;
			lq += (*data_.sample_line_quadric)[idx] * weight;
		}
		const Mat4 Ql = lq.get_quadric().matrix();
		const Mat3 Al = Ql.block<3, 3>(0, 0);
		const Vec3 bl = -Ql.block<3, 1>(0, 3);
		const Mat3 As = q._A.block<3, 3>(0, 0);
		const Vec3 bs = q._b.head<3>();
		const Vec3 Asr = q._A.block<3, 1>(0, 3);
		const Mat3 A = As + data_.sqem_update_lambda_line_plane * Al;
		const Vec3 b = (bs + data_.sqem_update_lambda_line_plane * bl) - Asr * fixed_radius;
		center = A.ldlt().solve(b);
		if (!center.allFinite())
		{
			nonfinite_updates_.fetch_add(1, std::memory_order_relaxed);
			return;
		}
		(*data_.sphere_position)[sphere_index] = center;
		(*data_.sphere_radius)[sphere_index] = fixed_radius;
	}

	void update_sphere_line_quadric_distance_free_radius(Vertex sphere)
	{
		const uint32 sphere_index = index_of(*data_.spheres, sphere);
		const std::vector<Vertex>& cluster = (*data_.sphere_cluster)[sphere_index];
		const Scalar radius = (*data_.sphere_radius)[sphere_index];
		if (cluster.empty())
			return;
		Spherical_Quadric q;
		Line_Quadric lq;
		Scalar weight_sum = Scalar(0);
		for (Vertex sample : cluster)
		{
			const uint32 idx = index_of(*data_.samples_mesh, sample);
			const Scalar weight = value<Scalar>(*data_.samples_mesh, data_.sample_area, sample);
			if (weight <= Scalar(0))
				continue;
			q += (*data_.sample_quadric)[idx] * weight;
			lq += (*data_.sample_line_quadric)[idx] * weight;
			weight_sum += weight;
		}
		if (weight_sum <= Scalar(0))
			return;
		const Mat4 Ql = lq.get_quadric().matrix();
		const Mat3 Al = Ql.block<3, 3>(0, 0);
		const Vec3 bl = -Ql.block<3, 1>(0, 3);
		Mat4 Al_ext = Mat4::Zero();
		Al_ext.block<3, 3>(0, 0) = Al;
		Vec4 bl_ext = Vec4::Zero();
		bl_ext.head<3>() = bl;
		const Mat4 A = q._A + data_.sqem_update_lambda_line_plane * Al_ext;
		const Vec4 b = q._b + data_.sqem_update_lambda_line_plane * bl_ext;
		const Vec4 s = A.completeOrthogonalDecomposition().solve(b);
		if (!s.allFinite())
		{
			nonfinite_updates_.fetch_add(1, std::memory_order_relaxed);
			return;
		}
		Scalar nearest_ma_radius = Scalar(0);
		const bool too_large = try_get_nearest_sample_ma_radius(s.head<3>(), nearest_ma_radius) &&
			(s[3] > nearest_ma_radius * Scalar(1.5));
		if (too_large || s[3] <= Scalar(0))
		{
			update_sphere_line_quadric_distance_fix_current_radius(sphere, radius);
			return;
		}
		(*data_.sphere_position)[sphere_index] = s.head<3>();
		(*data_.sphere_radius)[sphere_index] = s[3];
	}

	void redistribute_removed_sphere_cluster(const std::vector<Vertex>& removed_cluster,
										 const std::vector<Vertex>& removed_neighbors)
	{
		if (removed_cluster.empty() || metrics_.sphere_count == 0)
			return;
		std::vector<Vertex> candidates;
		std::unordered_set<uint32> candidate_indices;
		auto add_candidate = [&](Vertex sphere) {
			if (!sphere.is_valid())
				return;
			const uint32 idx = index_of(*data_.spheres, sphere);
			if (idx != INVALID_INDEX && candidate_indices.insert(idx).second)
				candidates.push_back(sphere);
		};
		for (Vertex neighbor : removed_neighbors)
			add_candidate(neighbor);
		if (candidates.empty())
			foreach_cell(*data_.spheres, [&](Vertex sphere) { add_candidate(sphere); return true; });
		for (Vertex sample : removed_cluster)
		{
			const uint32 sample_index = index_of(*data_.samples_mesh, sample);
			if (sample_index == INVALID_INDEX)
				continue;
			Scalar min_distance = std::numeric_limits<Scalar>::max();
			Vertex closest_sphere;
			uint32 closest_index = INVALID_INDEX;
			for (Vertex candidate : candidates)
			{
				const uint32 idx = index_of(*data_.spheres, candidate);
				const Vec3& center = (*data_.sphere_position)[idx];
				const Scalar radius = (*data_.sphere_radius)[idx];
				const Scalar distance =
					(*data_.sample_quadric)[sample_index].eval(Vec4(center.x(), center.y(), center.z(), radius)) +
					data_.sqem_update_lambda_line_plane * (*data_.sample_line_quadric)[sample_index].eval(center);
				if (distance < min_distance)
				{
					min_distance = distance;
					closest_sphere = candidate;
					closest_index = idx;
				}
			}
			if (closest_index == INVALID_INDEX)
				continue;
			(*data_.sample_sphere)[sample_index] = closest_sphere;
			(*data_.sphere_cluster)[closest_index].push_back(sample);
			(*data_.sphere_cluster_area)[closest_index] += (*data_.sample_area)[sample_index];
		}
	}

	Status update_stop_status()
	{
		constexpr uint32 max_iterations = 150;
		constexpr uint32 max_post_convergence_iterations = 10;
		constexpr Scalar convergence_epsilon = Scalar(1e-10);
		if (metrics_.error_difference < convergence_epsilon)
		{
			if (!convergence_reached_)
			{
				convergence_reached_ = true;
				post_convergence_iterations_ = 0;
			}
			else
				++post_convergence_iterations_;
			if (post_convergence_iterations_ >= max_post_convergence_iterations)
				return Status::converged;
		}
		else if (convergence_reached_)
		{
			convergence_reached_ = false;
			post_convergence_iterations_ = 0;
		}
		return metrics_.iteration >= max_iterations ? Status::max_iterations : Status::running;
	}

	Data data_;
	Metrics metrics_;
	bool convergence_reached_ = false;
	// non-finite solves of the current update pass (sphere updates run in parallel)
	std::atomic<uint32> nonfinite_updates_{0};
	uint32 post_convergence_iterations_ = 0;
};

} // namespace offset_axis

} // namespace geometry
} // namespace cgogn

#endif // CGOGN_GEOMETRY_ALGOS_OFFSET_AXIS_SPHERES_OPTIMIZER_H_
