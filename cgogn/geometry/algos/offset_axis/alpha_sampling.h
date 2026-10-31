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
#ifndef CGOGN_GEOMETRY_ALGOS_OFFSET_AXIS_ALPHA_SAMPLING_H_
#define CGOGN_GEOMETRY_ALGOS_OFFSET_AXIS_ALPHA_SAMPLING_H_

#include <cgogn/core/types/mesh_traits.h>
#include <cgogn/core/utils/assert.h>
#include <cgogn/geometry/algos/offset_axis/alpha_projection.h>
#include <cgogn/geometry/algos/offset_axis/distance_field.h>
#include <cgogn/geometry/algos/offset_axis/normal_utils.h>
#include <cgogn/geometry/algos/offset_axis/ray_level_set_sampler.h>
#include <cgogn/geometry/algos/offset_axis/reconstruction_options.h>
#include <cgogn/geometry/algos/offset_axis/spatial_grid.h>
#include <cgogn/geometry/types/vector_traits.h>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <optional>
#include <random>
#include <unordered_map>
#include <vector>

namespace cgogn
{

namespace geometry
{

namespace offset_axis
{

namespace constants
{

// Padding added on each side of the input bounding box. The padded box is the sampling domain and the box the
// projected points are clamped to.
inline constexpr Scalar sampling_padding = Scalar(0.1);
// Relative enlargement of the padded sampling box from which the seed rays are cast (its half size is multiplied by
// 1 + ray_box_relative_expand). Overrides the ray sampler default (0.1).
inline constexpr Scalar ray_box_relative_expand = Scalar(0.05);
// Final reprojection of the samples: lengths relative to the sample radius. The Newton direction is the central
// difference gradient of the field values, since the gradient predicted by the MF networks is not the gradient of
// their values.
inline constexpr Scalar reprojection_difference_step = Scalar(0.05);
inline constexpr Scalar reprojection_max_step = Scalar(0.5);
// from the sampled position
inline constexpr Scalar reprojection_max_displacement = Scalar(3);
// on |f - alpha| / |grad f|
inline constexpr Scalar reprojection_tolerance = Scalar(0.01);

} // namespace constants

enum class ReprojectionState : uint8_t
{
	active,
	converged,
	rejected // farther than the maximum displacement from its sampled position, or no gradient
};

struct AlphaSamplingParams
{
	Scalar sample_radius = Scalar(0);
	size_t candidates_per_parent = 1;
	size_t parent_batch_size = 1;
	int reprojection_iterations = 0;
	// field values queried per call during the reprojection
	size_t query_batch_size = 1;
	RayLevelSetParams rays;
};

inline AlphaSamplingParams alpha_sampling_params(const ReconstructionOptions& options)
{
	AlphaSamplingParams params;
	params.sample_radius = Scalar(options.sample_radius);
	params.candidates_per_parent = static_cast<size_t>(options.alpha_sampling_candidates_per_parent);
	params.parent_batch_size = static_cast<size_t>(options.alpha_sampling_parent_batch_size);
	params.reprojection_iterations = options.reprojection_iterations;
	params.query_batch_size = static_cast<size_t>(options.query_batch_size);
	params.rays.bbox_expand = constants::ray_box_relative_expand;
	params.rays.alpha = Scalar(options.alpha);
	params.rays.tol = Scalar(options.tolerance);
	params.rays.step_bound = Scalar(2.0);
	params.rays.batch_size = options.ray_sampler_batch_size;
	params.rays.max_iterations = options.udf_max_iterations;
	params.rays.min_active_rays = options.udf_min_active_rays;
	params.rays.seed = options.seed;
	return params;
}

struct AlphaSamplingResult
{
	std::vector<Vec3> positions;
	std::vector<Vec3> normals;
	bool success = true;
};

enum class SamplingStatus
{
	ok,
	invalid_data,
	failed
};

namespace detail
{

inline constexpr Scalar alpha_sampling_outer_radius_scale = Scalar(1.25);
inline constexpr size_t alpha_sampling_warmup_iterations = 2;
inline constexpr size_t alpha_sampling_max_samples = 10000000;

inline void build_stable_tangent_basis(const Vec3& normal_in, Vec3& tangent_u, Vec3& tangent_v)
{
	Vec3 n = unit_normal_or_z(normal_in);
	Vec3 ref = (std::abs(n.z()) < Scalar(0.9)) ? Vec3(0, 0, 1) : Vec3(0, 1, 0);
	// |n x ref| >= sqrt(1 - 0.9^2) since ref is not the axis n is closest to
	tangent_u = n.cross(ref).normalized();
	tangent_v = n.cross(tangent_u).normalized();
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

template <typename Traits, typename NormalFn, typename ProjectionFn>
AlphaSamplingResult sample_alpha_level_set(Traits& traits, const AlphaSamplingParams& parameters, const Vec3& bbox_min,
										   const Vec3& bbox_max, NormalFn&& compute_normals,
										   ProjectionFn&& project_points)
{
	AlphaSamplingResult result;

	// cell = radius, so every neighbor closer than radius lies in the 27 surrounding cells
	const Scalar grid_cell_size = parameters.sample_radius;
	const size_t parent_batch_size = parameters.parent_batch_size;
	const size_t max_samples = detail::alpha_sampling_max_samples;
	const size_t seed_rays_per_round =
		static_cast<size_t>(parameters.rays.batch_size) * detail::alpha_sampling_warmup_iterations;
	const size_t initial_reserve = std::max(seed_rays_per_round, parent_batch_size) * size_t(4);

	std::mt19937 generator(static_cast<uint32>(parameters.rays.seed));
	SpatialGrid spatial_grid(grid_cell_size);
	std::vector<uint32> active_indices;
	std::vector<uint32> active_slot_of_sample;
	std::vector<uint32> selected_parents;
	std::vector<Vec3> raw_candidates;
	AlphaProjectionResult projection;
	result.positions.reserve(initial_reserve);
	result.normals.reserve(initial_reserve);
	active_indices.reserve(initial_reserve);
	active_slot_of_sample.reserve(initial_reserve);
	raw_candidates.reserve(parent_batch_size * parameters.candidates_per_parent);

	auto add_sample = [&](const Vec3& position, const Vec3& normal) {
		const uint32 sample_index = static_cast<uint32>(result.positions.size());
		result.positions.push_back(position);
		result.normals.push_back(unit_normal_or_z(normal));
		active_slot_of_sample.push_back(static_cast<uint32>(active_indices.size()));
		active_indices.push_back(sample_index);
		spatial_grid.insert(position, sample_index);
	};

	// the parents of a batch are distinct active samples
	auto remove_active_sample = [&](uint32 sample_index) {
		const uint32 slot = active_slot_of_sample[sample_index];
		cgogn_message_assert(slot < active_indices.size(), "removed sample is not active");
		const uint32 last_sample = active_indices.back();
		active_indices[slot] = last_sample;
		active_slot_of_sample[last_sample] = slot;
		active_indices.pop_back();
		active_slot_of_sample[sample_index] = std::numeric_limits<uint32>::max();
	};

	// empty grid: no samples exist yet at this point, so nothing is inserted
	SpatialGrid working_grid(grid_cell_size);

	std::optional<std::vector<Vec3>> sampled_seeds = sample_level_set_rays(
		traits, parameters.rays, static_cast<int>(detail::alpha_sampling_warmup_iterations), &working_grid,
		parameters.sample_radius, bbox_min, bbox_max, &result.positions);
	if (!sampled_seeds)
	{
		result.success = false;
		return result;
	}
	std::vector<Vec3>& seed_points = *sampled_seeds;
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
		cgogn_message_assert(seed_normals.size() == seed_points.size() && seed_keep_mask.size() == seed_points.size(),
							 "compute_normals output size mismatch");
		// a neural gradient can be non finite
		for (size_t i = 0; i < seed_points.size(); ++i)
			if (seed_keep_mask[i] != 0 && seed_normals[i].allFinite())
				add_sample(seed_points[i], seed_normals[i]);
	}

	while (result.positions.size() < max_samples && !active_indices.empty())
	{
		detail::sample_active_parent_indices(active_indices, parent_batch_size, generator, selected_parents);
		raw_candidates.clear();
		std::uniform_real_distribution<Scalar> uniform(Scalar(0), Scalar(1));
		for (uint32 parent_index : selected_parents)
			for (size_t candidate_index = 0; candidate_index < parameters.candidates_per_parent;
				 candidate_index++)
				raw_candidates.push_back(detail::sample_tangent_annulus_candidate(
					result.positions[parent_index], result.normals[parent_index], parameters.sample_radius, uniform, generator));

		if (!project_points(raw_candidates, projection))
		{
			result.success = false;
			break;
		}
		cgogn_message_assert(projection.projected_points.size() == raw_candidates.size() &&
								 projection.normals.size() == raw_candidates.size() &&
								 projection.keep_mask.size() == raw_candidates.size(),
							 "project_points output size mismatch");

		for (size_t parent_batch_index = 0; parent_batch_index < selected_parents.size(); ++parent_batch_index)
		{
			if (result.positions.size() >= max_samples)
				break;
			const uint32 parent_index = selected_parents[parent_batch_index];
			bool parent_found_new_sample = false;
			const size_t candidate_begin = parent_batch_index * parameters.candidates_per_parent;
			for (size_t candidate_offset = 0; candidate_offset < parameters.candidates_per_parent;
				 candidate_offset++)
			{
				if (result.positions.size() >= max_samples)
					break;
				const size_t candidate_index = candidate_begin + candidate_offset;
				const Vec3& child_position = projection.projected_points[candidate_index];
				const Vec3& child_normal = projection.normals[candidate_index];
				const Vec3& raw_candidate = raw_candidates[candidate_index];
				bool accept = projection.keep_mask[candidate_index] != 0;
				if (accept && (child_position - raw_candidate).norm() > parameters.sample_radius)
					accept = false;
				if (accept && !spatial_grid.is_valid_sample(child_position, parameters.sample_radius,
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

// Evaluates the MF field and its central difference gradient at the active samples, stores the normalized gradient as
// their normal and their distance to the level set (f - alpha) / |grad f| as their residual, marks the converged ones
// (rejected instead on the positive sdf side, the alpha level set outside the shape) and, when move is set, moves the
// others by one Newton step. Each sample takes 7 queries, so the field is queried at most query_batch_size / 7 samples
// at a time.
template <typename FIELD, typename POINTS>
bool reprojection_iteration(FIELD& field, const POINTS& samples,
							typename mesh_traits<POINTS>::template Attribute<Vec3>* sample_position,
							typename mesh_traits<POINTS>::template Attribute<Vec3>* sample_normal,
							const typename mesh_traits<POINTS>::template Attribute<Vec3>* sample_origin,
							typename mesh_traits<POINTS>::template Attribute<ReprojectionState>* sample_state,
							typename mesh_traits<POINTS>::template Attribute<Scalar>* sample_residual, Scalar alpha,
							Scalar sample_radius, size_t query_batch_size, bool move)
{
	using PVertex = typename mesh_traits<POINTS>::Vertex;

	std::vector<uint32> active;
	foreach_cell(samples, [&](PVertex v) {
		const uint32 i = index_of(samples, v);
		if ((*sample_state)[i] == ReprojectionState::active)
			active.push_back(i);
		return true;
	});
	if (active.empty())
		return true;

	const Scalar h = constants::reprojection_difference_step * sample_radius;
	const Scalar tolerance = constants::reprojection_tolerance * sample_radius;
	const Scalar max_step = constants::reprojection_max_step * sample_radius;
	const Scalar max_displacement = constants::reprojection_max_displacement * sample_radius;
	const size_t chunk = std::max<size_t>(1, query_batch_size / 7);
	std::vector<Vec3> points;
	std::vector<Scalar> values, sdf;
	for (size_t begin = 0; begin < active.size(); begin += chunk)
	{
		const size_t count = std::min(chunk, active.size() - begin);
		points.resize(count * 7);
		for (size_t a = 0; a < count; ++a)
		{
			const Vec3& x = (*sample_position)[active[begin + a]];
			points[7 * a] = x;
			for (int k = 0; k < 3; ++k)
			{
				Vec3 e(0, 0, 0);
				e[k] = h;
				points[7 * a + 1 + 2 * k] = x + e;
				points[7 * a + 2 + 2 * k] = x - e;
			}
		}
		if (!field.values_sdf(points, values, sdf))
			return false;
		for (size_t a = 0; a < count; ++a)
		{
			const uint32 i = active[begin + a];
			const Scalar* f = &values[7 * a];
			const Vec3 g((f[1] - f[2]) / (2 * h), (f[3] - f[4]) / (2 * h), (f[5] - f[6]) / (2 * h));
			const Scalar g_norm = g.norm();
			// the network can return non finite values far from the input
			if (!(g_norm > Scalar(1e-12)) || !std::isfinite(f[0]))
			{
				(*sample_state)[i] = ReprojectionState::rejected;
				continue;
			}
			const Vec3 n = g / g_norm;
			const Scalar residual = (f[0] - alpha) / g_norm;
			(*sample_normal)[i] = n;
			(*sample_residual)[i] = residual;
			if (std::abs(residual) <= tolerance)
			{
				(*sample_state)[i] = sdf[7 * a] > 0 ? ReprojectionState::rejected : ReprojectionState::converged;
				continue;
			}
			if (!move)
				continue;
			Vec3& x = (*sample_position)[i];
			x -= std::clamp(residual, -max_step, max_step) * n;
			if ((x - (*sample_origin)[i]).norm() > max_displacement)
				(*sample_state)[i] = ReprojectionState::rejected;
		}
	}
	return true;
}

} // namespace detail

// Starts the reprojection of the samples onto the alpha level set: every sample is active and its position is kept in
// sample_origin.
template <typename POINTS>
void start_samples_reprojection(const POINTS& samples,
								const typename mesh_traits<POINTS>::template Attribute<Vec3>* sample_position,
								typename mesh_traits<POINTS>::template Attribute<Vec3>* sample_origin,
								typename mesh_traits<POINTS>::template Attribute<ReprojectionState>* sample_state,
								typename mesh_traits<POINTS>::template Attribute<Scalar>* sample_residual)
{
	using PVertex = typename mesh_traits<POINTS>::Vertex;
	foreach_cell(samples, [&](PVertex v) {
		const uint32 i = index_of(samples, v);
		(*sample_origin)[i] = (*sample_position)[i];
		(*sample_state)[i] = ReprojectionState::active;
		(*sample_residual)[i] = Scalar(0);
		return true;
	});
}

// One Newton iteration of the active samples (see detail::reprojection_iteration). Each step is at most
// reprojection_max_step; a sample that ends farther than reprojection_max_displacement from its sampled position is
// rejected.
template <typename FIELD, typename POINTS>
bool reproject_samples_step(FIELD& field, const POINTS& samples,
							typename mesh_traits<POINTS>::template Attribute<Vec3>* sample_position,
							typename mesh_traits<POINTS>::template Attribute<Vec3>* sample_normal,
							const typename mesh_traits<POINTS>::template Attribute<Vec3>* sample_origin,
							typename mesh_traits<POINTS>::template Attribute<ReprojectionState>* sample_state,
							typename mesh_traits<POINTS>::template Attribute<Scalar>* sample_residual, Scalar alpha,
							Scalar sample_radius, size_t query_batch_size)
{
	return detail::reprojection_iteration(field, samples, sample_position, sample_normal, sample_origin, sample_state,
										  sample_residual, alpha, sample_radius, query_batch_size, true);
}

// Ends the reprojection: the samples still active are evaluated once more, then only the converged samples are kept.
template <typename FIELD, typename POINTS>
bool finish_samples_reprojection(FIELD& field, POINTS& samples,
								 typename mesh_traits<POINTS>::template Attribute<Vec3>* sample_position,
								 typename mesh_traits<POINTS>::template Attribute<Vec3>* sample_normal,
								 const typename mesh_traits<POINTS>::template Attribute<Vec3>* sample_origin,
								 typename mesh_traits<POINTS>::template Attribute<ReprojectionState>* sample_state,
								 typename mesh_traits<POINTS>::template Attribute<Scalar>* sample_residual, Scalar alpha,
								 Scalar sample_radius, size_t query_batch_size)
{
	using PVertex = typename mesh_traits<POINTS>::Vertex;

	if (!detail::reprojection_iteration(field, samples, sample_position, sample_normal, sample_origin, sample_state,
										sample_residual, alpha, sample_radius, query_batch_size, false))
		return false;
	std::vector<Vec3> positions, normals;
	foreach_cell(samples, [&](PVertex v) {
		const uint32 i = index_of(samples, v);
		if ((*sample_state)[i] != ReprojectionState::converged)
			return true;
		positions.push_back((*sample_position)[i]);
		normals.push_back((*sample_normal)[i]);
		return true;
	});
	clear(samples);
	for (size_t i = 0; i < positions.size(); ++i)
	{
		PVertex v = add_vertex(samples);
		(*sample_position)[index_of(samples, v)] = positions[i];
		(*sample_normal)[index_of(samples, v)] = normals[i];
	}
	return true;
}

// Samples the alpha level set of the field and stores the samples in the samples point cloud (cleared first), with
// their positions and normals in sample_position and sample_normal. field.prepare(inputs) must have succeeded.
template <typename FIELD, typename INPUTS, typename POINTS>
SamplingStatus sample_alpha_level_set(FIELD& field, const INPUTS& inputs, const AlphaSamplingParams& params,
									  POINTS& samples,
									  typename mesh_traits<POINTS>::template Attribute<Vec3>* sample_position,
									  typename mesh_traits<POINTS>::template Attribute<Vec3>* sample_normal)
{
	using PVertex = typename mesh_traits<POINTS>::Vertex;

	Vec3 bbox_min(0, 0, 0), bbox_max(1, 1, 1);
	if (field.sampling_box(inputs, bbox_min, bbox_max) != FieldPrepareResult::ok)
		return SamplingStatus::invalid_data;
	const Scalar padding = constants::sampling_padding;
	bbox_min -= Vec3(padding, padding, padding);
	bbox_max += Vec3(padding, padding, padding);

	auto traits = field.ray_traits(inputs);
	cgogn_message_assert(traits.is_ready(), "field.prepare must succeed before sampling");
	auto compute_normals = [&](const std::vector<Vec3>& points, std::vector<Vec3>& normals,
							   std::vector<uint8_t>& keep) {
		keep.assign(points.size(), uint8_t(0));
		normals.assign(points.size(), Vec3(0, 0, 1));
		return field.normals(points, normals, keep);
	};
	const AlphaProjectionParameters pp{params.rays.alpha, params.rays.tol, bbox_min, bbox_max, 1};
	auto project_points = [&](const std::vector<Vec3>& points, AlphaProjectionResult& projection) {
		return field.project(points, pp, projection);
	};
	AlphaSamplingResult result = detail::sample_alpha_level_set(traits, params, bbox_min, bbox_max, compute_normals, project_points);
	if (!result.success || result.positions.empty())
		return SamplingStatus::failed;
	clear(samples);
	for (size_t i = 0; i < result.positions.size(); ++i)
	{
		PVertex v = add_vertex(samples);
		const uint32 idx = index_of(samples, v);
		(*sample_position)[idx] = result.positions[i];
		(*sample_normal)[idx] = result.normals[i];
	}
	// the Bridson children are projected with one Newton iteration; the samples of an MF field are projected again
	if (params.reprojection_iterations > 0 && is_neural_model(field, NeuralModelType::mf))
	{
		auto origin = add_attribute<Vec3, PVertex>(samples, "__reprojection_origin");
		auto state = add_attribute<ReprojectionState, PVertex>(samples, "__reprojection_state");
		auto residual = add_attribute<Scalar, PVertex>(samples, "__reprojection_residual");
		start_samples_reprojection(samples, sample_position, origin.get(), state.get(), residual.get());
		bool ok = true;
		for (int it = 0; ok && it < params.reprojection_iterations; ++it)
			ok = reproject_samples_step(field, samples, sample_position, sample_normal, origin.get(), state.get(),
										residual.get(), params.rays.alpha, params.sample_radius,
										params.query_batch_size);
		ok = ok && finish_samples_reprojection(field, samples, sample_position, sample_normal, origin.get(),
											   state.get(), residual.get(), params.rays.alpha, params.sample_radius,
											   params.query_batch_size);
		remove_attribute<PVertex>(samples, origin);
		remove_attribute<PVertex>(samples, state);
		remove_attribute<PVertex>(samples, residual);
		if (!ok || nb_cells<PVertex>(samples) == 0)
			return SamplingStatus::failed;
	}
	return SamplingStatus::ok;
}

} // namespace offset_axis

} // namespace geometry

} // namespace cgogn

#endif // CGOGN_GEOMETRY_ALGOS_OFFSET_AXIS_ALPHA_SAMPLING_H_
