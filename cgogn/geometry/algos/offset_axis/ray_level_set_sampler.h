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
#ifndef CGOGN_GEOMETRY_ALGOS_OFFSET_AXIS_RAY_LEVEL_SET_SAMPLER_H_
#define CGOGN_GEOMETRY_ALGOS_OFFSET_AXIS_RAY_LEVEL_SET_SAMPLER_H_

#include <cgogn/core/utils/numerics.h>
#include <cgogn/geometry/algos/offset_axis/parallel_for.h>
#include <cgogn/geometry/algos/offset_axis/spatial_grid.h>
#include <cgogn/geometry/types/vector_traits.h>

#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>
#include <random>
#include <utility>
#include <vector>

namespace cgogn
{

namespace geometry
{

namespace offset_axis
{

// A ray of the level set sampler: origin + t * direction for t in [t_min, t_max].
struct LevelSetRay
{
	Vec3 origin;
	Vec3 direction;
	Scalar t_min;
	Scalar t_max;
};

// Parameters of the level set sampler.
struct RayLevelSetParams
{
	Scalar bbox_expand = Scalar(0.1);
	Scalar alpha = Scalar(0.0);
	Scalar tol = Scalar(1e-5);
	Scalar step_bound = Scalar(2.0);
	int batch_size = 0;
	int max_iterations = 0;
	// GPU marching: end a batch once fewer rays than this are still marching (0 = march until max_iterations)
	int min_active_rays = 0;
	int seed = 42;
};

namespace detail
{

// CPU and GPU samplers share these helpers

inline Scalar convergence_tolerance(const RayLevelSetParams& params)
{
	return std::max(Scalar(0), params.tol);
}

inline Vec3 clamp_point(const Vec3& p, const Vec3& bbox_min, const Vec3& bbox_max)
{
	return p.cwiseMax(bbox_min).cwiseMin(bbox_max);
}

inline std::pair<Scalar, Scalar> intersect_bbox(const LevelSetRay& ray, const Vec3& bbox_min, const Vec3& bbox_max)
{
	Scalar t_min = 0.0;
	Scalar t_max = std::numeric_limits<Scalar>::max();
	const Scalar eps = Scalar(1e-8);

	for (int i = 0; i < 3; ++i)
	{
		Scalar o = ray.origin[i];
		Scalar d = ray.direction[i];

		if (std::abs(d) < eps)
		{
			// Ray parallel to slab
			if (o < bbox_min[i] || o > bbox_max[i])
				return {Scalar(1), Scalar(0)}; // no intersection
			continue;
		}

		Scalar t1 = (bbox_min[i] - o) / d;
		Scalar t2 = (bbox_max[i] - o) / d;
		if (t1 > t2)
			std::swap(t1, t2);

		t_min = std::max(t_min, t1);
		t_max = std::min(t_max, t2);

		if (t_max < t_min)
			return {Scalar(1), Scalar(0)};
	}

	return {t_min, t_max};
}

inline std::vector<LevelSetRay> generate_rays(const RayLevelSetParams& params, const Vec3& bbox_min,
											  const Vec3& bbox_max, int num_rays, std::mt19937& gen)
{
	std::vector<LevelSetRay> rays;
	rays.reserve(num_rays);

	Vec3 center = (bbox_min + bbox_max) * 0.5;
	Vec3 half_size = (bbox_max - bbox_min) * 0.5 * (1.0 + params.bbox_expand);
	Vec3 expanded_min = center - half_size;
	Vec3 expanded_max = center + half_size;
	const Scalar diag = (expanded_max - expanded_min).norm();

	std::uniform_real_distribution<Scalar> uniform(0.0, 1.0);
	auto uniform_signed = [&](Scalar a) {
		return (uniform(gen) * 2.0 - 1.0) * a; // [-a, a]
	};
	for (int i = 0; i < num_rays; ++i)
	{
		LevelSetRay ray;

		// Uniform direction on sphere (Marsaglia method)
		Scalar u1, u2, s;
		do
		{
			u1 = uniform(gen) * 2.0 - 1.0;
			u2 = uniform(gen) * 2.0 - 1.0;
			s = u1 * u1 + u2 * u2;
		} while (s >= 1.0);

		Scalar factor = 2.0 * std::sqrt(1.0 - s);
		ray.direction = Vec3(u1 * factor, u2 * factor, 1.0 - 2.0 * s);

		// Generate orthogonal basis
		Vec3 l = ray.direction;

		// Chose a vector not aligned with l
		Vec3 a, n, b;
		Vec3 ad = l.cwiseAbs();
		if (ad.x() <= ad.y() && ad.x() <= ad.z())
			a = Vec3(1, 0, 0);
		else if (ad.y() <= ad.x() && ad.y() <= ad.z())
			a = Vec3(0, 1, 0);
		else
			a = Vec3(0, 0, 1);

		n = l.cross(a).normalized();
		b = l.cross(n);

		// Sample on plane perpendicular to major axis
		while (true)
		{
			Scalar u = uniform_signed(diag * 0.5);
			Scalar v = uniform_signed(diag * 0.5);

			Vec3 q = center + n * u + b * v; // candidate point on plane

			LevelSetRay tmp;
			tmp.origin = q;
			tmp.direction = l;

			auto [t_enter, t_exit] = intersect_bbox(tmp, expanded_min, expanded_max);
			if (!(t_exit > t_enter))
				continue;

			// Unique ray: start at entry point
			ray.origin = q + t_enter * l;
			ray.t_min = 0.0;
			ray.t_max = t_exit - t_enter;
			break;
		}
		rays.push_back(ray);
	}

	return rays;
}

} // namespace detail

// CPU level set sampler: the rays of each outer iteration are marched in lockstep, and the points where the udf of the
// traits reaches params.alpha are returned. With a spatial grid, a point is kept only if it is not closer than
// sample_spacing to the other samples (the existing samples are the neighbors that count).
template <typename TRAITS>
std::vector<Vec3> sample_level_set_rays(const TRAITS& traits, const RayLevelSetParams& params, int outer_iterations,
										SpatialGrid* spatial_grid, Scalar sample_spacing, const Vec3& bbox_min,
										const Vec3& bbox_max, const std::vector<Vec3>* existing_samples = nullptr)
{
	const int max_outer_iterations = std::max(1, outer_iterations);
	const size_t reserve_hint = std::max<size_t>(static_cast<size_t>(std::max(1, params.batch_size)) *
													 static_cast<size_t>(std::max(1, max_outer_iterations)),
												 size_t(1024));
	std::vector<Vec3> neighbor_samples;
	if (existing_samples)
		neighbor_samples = *existing_samples;
	neighbor_samples.reserve(neighbor_samples.size() + reserve_hint);
	std::vector<Vec3> all_samples;
	all_samples.reserve(reserve_hint);

	std::mt19937 gen(params.seed);
	const Scalar eps = detail::convergence_tolerance(params);
	const int max_steps = std::max(1, params.max_iterations);
	const int check_interval = 10;

	for (int iter = 0; iter < max_outer_iterations; ++iter)
	{
		auto rays = detail::generate_rays(params, bbox_min, bbox_max, params.batch_size, gen);
		const int R = static_cast<int>(rays.size());
		if (R == 0)
			break;

		std::cout << "Iteration " << iter << ": Processing " << R << " rays on CPU..." << std::endl;

		std::vector<Scalar> t(static_cast<size_t>(R));
		std::vector<Scalar> t_max(static_cast<size_t>(R));
		std::vector<Scalar> step_udf(static_cast<size_t>(R));
		for (int i = 0; i < R; ++i)
		{
			t[i] = rays[i].t_min;
			t_max[i] = rays[i].t_max;
		}
		int active_count = R;
		std::vector<Vec3> iter_samples;
		for (int step = 0; step < max_steps; ++step)
		{
			if (step % check_interval == 0)
			{
				if (active_count == 0)
					break;
				std::cout << "    Step " << step << ": " << active_count << " active\r" << std::flush;
			}

			// stage 1 (parallel): udf of all active rays; reads t/rays only, writes only step_udf[i]
			// (serial once few rays are left, the long tail would otherwise pay a pool dispatch per step)
			auto eval_ray = [&](std::size_t i) {
				if (t[i] < t_max[i])
				{
					const LevelSetRay& ray = rays[i];
					step_udf[i] = traits.eval_distance(ray.origin + t[i] * ray.direction);
				}
			};
			if (active_count >= 1024)
				detail::parallel_for_index(static_cast<std::size_t>(R), eval_ray);
			else
				for (std::size_t i = 0; i < static_cast<std::size_t>(R); ++i)
					eval_ray(i);

			// stage 2 (serial): same ordered logic as before, using the precomputed udf
			for (int i = 0; i < R; ++i)
			{
				if (t[i] >= t_max[i])
					continue;

				const LevelSetRay& ray = rays[i];
				const Vec3 X = ray.origin + t[i] * ray.direction;
				const Scalar udf = step_udf[i];
				if (!std::isfinite(static_cast<double>(udf)))
				{
					t[i] = t_max[i];
					continue;
				}

				const Scalar delta = std::abs(udf - params.alpha);
				if (delta < eps)
				{
					Vec3 pt = detail::clamp_point(X, bbox_min, bbox_max);
					if (!spatial_grid || spatial_grid->is_valid_sample(pt, sample_spacing, neighbor_samples))
					{
						iter_samples.push_back(pt);
						all_samples.push_back(pt);
						neighbor_samples.push_back(pt);
						if (spatial_grid)
							spatial_grid->insert(pt, static_cast<uint32>(neighbor_samples.size() - 1));
					}

					const Scalar step_size = std::max(delta, eps * Scalar(0.5));
					t[i] += step_size;
				}
				else
				{
					const Scalar step_size = std::max(delta / params.step_bound, eps * Scalar(0.5));
					t[i] += step_size;
				}
				if (t[i] >= t_max[i])
					--active_count;
			}
		}

		std::cout << std::endl;
		std::cout << "  Collected " << iter_samples.size() << " unique samples this iteration" << std::endl;
		std::cout << "  Total: " << all_samples.size() << std::endl;
	}

	return all_samples;
}

} // namespace offset_axis

} // namespace geometry

} // namespace cgogn

#endif // CGOGN_GEOMETRY_ALGOS_OFFSET_AXIS_RAY_LEVEL_SET_SAMPLER_H_
