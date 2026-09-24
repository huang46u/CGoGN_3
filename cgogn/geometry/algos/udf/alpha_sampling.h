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
#ifndef CGOGN_GEOMETRY_ALGOS_UDF_ALPHA_SAMPLING_H_
#define CGOGN_GEOMETRY_ALGOS_UDF_ALPHA_SAMPLING_H_

#include <cgogn/geometry/types/spatial_grid.h>
#include <cgogn/geometry/types/vector_traits.h>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <memory>
#include <random>
#include <unordered_map>
#include <vector>

namespace cgogn
{

namespace geometry
{

struct AlphaSamplingParameters
{
	Scalar sample_radius = Scalar(0);
	uint32 seed = 0;
	size_t ray_sampler_batch_size = 1;
	Vec3 bbox_min = Vec3(0, 0, 0);
	Vec3 bbox_max = Vec3(1, 1, 1);
};

struct AlphaSamplingResult
{
	std::vector<Vec3> positions;
	std::vector<Vec3> normals;
	std::unique_ptr<SpatialGrid> spatial_grid;
	bool success = true;
};

namespace detail
{

inline constexpr Scalar alpha_sampling_outer_radius_scale = Scalar(1.25);
inline constexpr size_t alpha_sampling_warmup_iterations = 2;
inline constexpr size_t alpha_sampling_candidates_per_parent = 12;
inline constexpr size_t alpha_sampling_parent_batch_size = 1024;
inline constexpr size_t alpha_sampling_max_samples = 10000000;

inline void build_stable_tangent_basis(const Vec3& normal_in, Vec3& tangent_u, Vec3& tangent_v)
{
	Vec3 n = normal_in;
	if (n.squaredNorm() < Scalar(1e-12))
		n = Vec3(0, 0, 1);
	else
		n.normalize();
	Vec3 ref = (std::abs(n.z()) < Scalar(0.9)) ? Vec3(0, 0, 1) : Vec3(0, 1, 0);
	tangent_u = n.cross(ref);
	if (tangent_u.squaredNorm() < Scalar(1e-12))
		tangent_u = n.cross(Vec3(1, 0, 0));
	tangent_u.normalize();
	tangent_v = n.cross(tangent_u);
	if (tangent_v.squaredNorm() < Scalar(1e-12))
		tangent_v = Vec3(0, 1, 0);
	else
		tangent_v.normalize();
}

inline Vec3 sample_tangent_annulus_candidate(const Vec3& center, const Vec3& normal, Scalar radius,
												std::uniform_real_distribution<Scalar>& uniform, std::mt19937& generator)
{
	Vec3 tangent_u, tangent_v;
	build_stable_tangent_basis(normal, tangent_u, tangent_v);
	const Scalar outer_radius_sq = alpha_sampling_outer_radius_scale * alpha_sampling_outer_radius_scale;
	const Scalar annulus_radius = radius *
		std::sqrt(Scalar(1.0) + (outer_radius_sq - Scalar(1.0)) * uniform(generator));
	const Scalar theta = Scalar(2.0 * M_PI) * uniform(generator);
	return center + annulus_radius * (std::cos(theta) * tangent_u + std::sin(theta) * tangent_v);
}

inline void sample_active_parent_indices(const std::vector<uint32>& active_indices, size_t batch_size,
											 std::mt19937& generator, std::vector<uint32>& selected_samples)
{
	selected_samples.clear();
	const size_t active_count = active_indices.size();
	if (active_count == 0)
		return;

	const size_t pick_count = std::min(batch_size, active_count);
	selected_samples.reserve(pick_count);
	std::unordered_map<size_t, size_t> remap;
	remap.reserve(pick_count * 2);
	for (size_t i = 0; i < pick_count; ++i)
	{
		std::uniform_int_distribution<size_t> distribution(i, active_count - 1);
		const size_t j = distribution(generator);
		const size_t mapped_j = remap.count(j) ? remap[j] : j;
		const size_t mapped_i = remap.count(i) ? remap[i] : i;
		remap[j] = mapped_i;
		remap[i] = mapped_j;
		selected_samples.push_back(active_indices[mapped_j]);
	}
}

} // namespace detail

template <typename Traits, typename RaySampler, typename NormalFn, typename ProjectionFn>
AlphaSamplingResult sample_alpha_level_set(Traits& traits, RaySampler& ray_sampler,
											   const AlphaSamplingParameters& parameters, NormalFn&& compute_normals,
											   ProjectionFn&& project_points)
{
	AlphaSamplingResult result;
	if (parameters.sample_radius <= Scalar(0) || !traits.is_ready())
	{
		result.success = false;
		return result;
	}

	const Scalar grid_cell_size = std::max(parameters.sample_radius / std::sqrt(Scalar(3.0)), Scalar(1e-8));
	const size_t parent_batch_size = detail::alpha_sampling_parent_batch_size;
	const size_t max_samples = detail::alpha_sampling_max_samples;
	const size_t seed_rays_per_round =
		std::max<size_t>(1, parameters.ray_sampler_batch_size) * detail::alpha_sampling_warmup_iterations;
	const size_t initial_reserve = std::max(seed_rays_per_round, parent_batch_size) * size_t(4);

	std::mt19937 generator(parameters.seed);
	result.spatial_grid = std::make_unique<SpatialGrid>(grid_cell_size);
	std::vector<uint32> active_indices;
	std::vector<uint32> active_slot_of_sample;
	std::vector<uint32> selected_parents;
	std::vector<Vec3> raw_candidates;
	std::vector<Vec3> projected_candidates;
	std::vector<Vec3> projected_normals;
	std::vector<uint8_t> projected_keep_mask;
	result.positions.reserve(initial_reserve);
	result.normals.reserve(initial_reserve);
	active_indices.reserve(initial_reserve);
	active_slot_of_sample.reserve(initial_reserve);
	raw_candidates.reserve(parent_batch_size * detail::alpha_sampling_candidates_per_parent);
	projected_candidates.reserve(parent_batch_size * detail::alpha_sampling_candidates_per_parent);
	projected_normals.reserve(parent_batch_size * detail::alpha_sampling_candidates_per_parent);

	auto add_sample = [&](const Vec3& position, const Vec3& normal) {
		const uint32 sample_index = static_cast<uint32>(result.positions.size());
		Vec3 normalized_normal = normal;
		if (normalized_normal.squaredNorm() < Scalar(1e-12))
			normalized_normal = Vec3(0, 0, 1);
		else
			normalized_normal.normalize();
		result.positions.push_back(position);
		result.normals.push_back(normalized_normal);
		active_slot_of_sample.push_back(static_cast<uint32>(active_indices.size()));
		active_indices.push_back(sample_index);
		result.spatial_grid->insert(position, sample_index);
	};

	auto remove_active_sample = [&](uint32 sample_index) {
		if (sample_index >= active_slot_of_sample.size())
			return;
		const uint32 slot = active_slot_of_sample[sample_index];
		if (slot == std::numeric_limits<uint32>::max() || slot >= active_indices.size())
			return;
		const uint32 last_sample = active_indices.back();
		active_indices[slot] = last_sample;
		active_slot_of_sample[last_sample] = slot;
		active_indices.pop_back();
		active_slot_of_sample[sample_index] = std::numeric_limits<uint32>::max();
	};

	if (result.positions.size() < max_samples)
	{
		SpatialGrid working_grid(grid_cell_size);
		for (uint32 sample_index = 0; sample_index < result.positions.size(); ++sample_index)
			working_grid.insert(result.positions[sample_index], sample_index);

		bool used_sdf_filter = false;
		std::vector<Vec3> seed_points = ray_sampler.sample_alpha_level_set_rays_fixed_iterations(
			traits, static_cast<int>(detail::alpha_sampling_warmup_iterations), &working_grid, parameters.sample_radius,
			parameters.bbox_min, parameters.bbox_max, &result.positions, &used_sdf_filter);
		if (!seed_points.empty())
		{
			const size_t remaining_slots = max_samples - result.positions.size();
			if (seed_points.size() > remaining_slots)
				seed_points.resize(remaining_slots);

			std::vector<Vec3> seed_normals;
			std::vector<uint8_t> seed_keep_mask;
			if (!compute_normals(seed_points, seed_normals, seed_keep_mask))
			{
				result.success = false;
				return result;
			}
			if (seed_normals.size() != seed_points.size() || seed_keep_mask.size() != seed_points.size())
			{
				result.success = false;
				return result;
			}
			std::vector<Vec3> accepted_seed_points;
			std::vector<Vec3> accepted_seed_normals;
			accepted_seed_points.reserve(seed_points.size());
			accepted_seed_normals.reserve(seed_points.size());
			for (size_t i = 0; i < seed_points.size(); ++i)
			{
				if (!seed_points[i].allFinite() || seed_keep_mask[i] == 0 || !seed_normals[i].allFinite())
					continue;
				if (result.positions.size() >= max_samples)
					break;
				accepted_seed_points.push_back(seed_points[i]);
				accepted_seed_normals.push_back(seed_normals[i]);
			}
			for (size_t i = 0; i < accepted_seed_points.size(); ++i)
				add_sample(accepted_seed_points[i], accepted_seed_normals[i]);
		}
	}

	while (result.positions.size() < max_samples && !active_indices.empty())
	{
		detail::sample_active_parent_indices(active_indices, parent_batch_size, generator, selected_parents);
		raw_candidates.clear();
		std::uniform_real_distribution<Scalar> uniform(Scalar(0), Scalar(1));
		for (uint32 parent_index : selected_parents)
			for (size_t candidate_index = 0; candidate_index < detail::alpha_sampling_candidates_per_parent;
				 candidate_index++)
				raw_candidates.push_back(detail::sample_tangent_annulus_candidate(
					result.positions[parent_index], result.normals[parent_index], parameters.sample_radius, uniform, generator));

		if (!project_points(raw_candidates, projected_candidates, projected_normals, projected_keep_mask))
		{
			result.success = false;
			break;
		}
		if (projected_candidates.size() != raw_candidates.size() || projected_normals.size() != raw_candidates.size() ||
			projected_keep_mask.size() != raw_candidates.size())
		{
			result.success = false;
			break;
		}

		for (size_t parent_batch_index = 0; parent_batch_index < selected_parents.size(); ++parent_batch_index)
		{
			if (result.positions.size() >= max_samples)
				break;
			const uint32 parent_index = selected_parents[parent_batch_index];
			bool parent_found_new_sample = false;
			const size_t candidate_begin = parent_batch_index * detail::alpha_sampling_candidates_per_parent;
			for (size_t candidate_offset = 0; candidate_offset < detail::alpha_sampling_candidates_per_parent;
				 candidate_offset++)
			{
				if (result.positions.size() >= max_samples)
					break;
				const size_t candidate_index = candidate_begin + candidate_offset;
				const Vec3& child_position = projected_candidates[candidate_index];
				const Vec3& child_normal = projected_normals[candidate_index];
				const Vec3& raw_candidate = raw_candidates[candidate_index];
				bool accept = projected_keep_mask[candidate_index] != 0;
				if (accept && (child_position - raw_candidate).norm() > parameters.sample_radius)
					accept = false;
				if (accept && !result.spatial_grid->is_valid_sample(child_position, parameters.sample_radius,
															 result.positions))
					accept = false;
				if (!accept)
					continue;
				add_sample(child_position, child_normal);
				parent_found_new_sample = true;
				break;
			}
			if (!parent_found_new_sample)
				remove_active_sample(parent_index);
		}
	}

	return result;
}

} // namespace geometry

} // namespace cgogn

#endif // CGOGN_GEOMETRY_ALGOS_UDF_ALPHA_SAMPLING_H_
