#ifndef CGOGN_GEOMETRY_TYPES_RAY_LEVEL_SET_SAMPLER_H_
#define CGOGN_GEOMETRY_TYPES_RAY_LEVEL_SET_SAMPLER_H_

#include <cgogn/core/utils/numerics.h>
#include <cgogn/core/utils/thread_pool.h>
#include <cgogn/geometry/types/spatial_grid.h>
#include <cgogn/geometry/types/vector_traits.h>

#include <algorithm>
#include <cmath>
#include <future>
#include <iostream>
#include <limits>
#include <random>
#include <utility>
#include <vector>

namespace cgogn
{

namespace geometry
{

template <typename Traits>
class RayLevelSetSampler
{
public:
	struct Ray
	{
		Vec3 origin;
		Vec3 direction;
		Scalar t_min;
		Scalar t_max;
	};

	struct Params
	{
		Scalar bbox_expand = Scalar(0.1);
		Scalar alpha = Scalar(0.0);
		Scalar tol = Scalar(1e-5);
		Scalar step_bound = Scalar(2.0);
		int batch_size = 0;
		int max_iterations = 0;
		// GPU marching: end a batch once fewer rays than this are still marching (0 = march until max_iterations)
		int min_active_rays = 0;
		int max_outer_iterations = 500;
		int seed = 42;
	};

	struct SampleResult
	{
		std::vector<Vec3> surface_samples;
		std::vector<Vec3> inside_samples;
	};

	explicit RayLevelSetSampler(const Params& params) : params_(params) {}

	void set_params(const Params& params)
	{
		params_ = params;
	}

	std::vector<Vec3> sample_alpha_level_set_rays(Traits& traits, size_t target_num_points,
												  SpatialGrid* spatial_grid, Scalar sample_spacing,
												  const Vec3& bbox_min, const Vec3& bbox_max,
												  bool* used_sdf_filter = nullptr)
	{
		if (used_sdf_filter)
			*used_sdf_filter = false;
		return sample_alpha_level_set_rays_cpu(traits, target_num_points, spatial_grid, sample_spacing, bbox_min,
											   bbox_max, params_.max_outer_iterations, nullptr);
	}

	std::vector<Vec3> sample_alpha_level_set_rays_fixed_iterations(Traits& traits, int outer_iterations,
																   SpatialGrid* spatial_grid, Scalar sample_spacing,
																   const Vec3& bbox_min, const Vec3& bbox_max,
																   const std::vector<Vec3>* existing_samples = nullptr,
																   bool* used_sdf_filter = nullptr)
	{
		const int clamped_outer_iterations = std::max(1, outer_iterations);
		if (used_sdf_filter)
			*used_sdf_filter = false;
		return sample_alpha_level_set_rays_cpu(traits, unlimited_sample_target_, spatial_grid, sample_spacing,
											   bbox_min, bbox_max, clamped_outer_iterations, existing_samples);
	}

	SampleResult sample_alpha_level_set_rays_with_inside(Traits& traits, size_t target_surface_points,
														 size_t target_inside_points, SpatialGrid* surface_grid,
														 SpatialGrid* inside_grid, Scalar grid_cell_size,
														 const Vec3& bbox_min, const Vec3& bbox_max,
														 bool* used_sdf_filter = nullptr)
	{
		if (target_inside_points == 0)
		{
			SampleResult result;
			result.surface_samples = sample_alpha_level_set_rays(traits, target_surface_points, surface_grid,
																 grid_cell_size, bbox_min, bbox_max,
																 used_sdf_filter);
			return result;
		}

		if (used_sdf_filter)
			*used_sdf_filter = false;
		return sample_alpha_level_set_rays_with_inside_cpu(traits, target_surface_points, target_inside_points,
														   surface_grid, inside_grid, grid_cell_size, bbox_min,
														   bbox_max);
	}

protected:
	static constexpr size_t unlimited_sample_target_ = std::numeric_limits<size_t>::max();

	static Vec3 clamp_point(const Vec3& p, const Vec3& bbox_min, const Vec3& bbox_max)
	{
		return p.cwiseMax(bbox_min).cwiseMin(bbox_max);
	}

	Scalar convergence_tolerance() const
	{
		return std::max(Scalar(0), params_.tol);
	}

	// runs f(i) for i in [0, n), in chunks on the thread pool; must not be called from a pool worker (nested wait)
	template <typename FUNC>
	static void parallel_for_rays(std::size_t n, const FUNC& f)
	{
		constexpr std::size_t min_parallel_size = 1024;
		ThreadPool* pool = thread_pool();
		const std::size_t nb_workers = pool->nb_workers();
		if (n < min_parallel_size || nb_workers == 0)
		{
			for (std::size_t i = 0; i < n; ++i)
				f(i);
			return;
		}

		const std::size_t nb_chunks = nb_workers * 4;
		const std::size_t chunk_size = (n + nb_chunks - 1) / nb_chunks;
		std::vector<std::future<void>> futures;
		futures.reserve(nb_chunks);
		for (std::size_t begin = 0; begin < n; begin += chunk_size)
		{
			const std::size_t end = std::min(n, begin + chunk_size);
			futures.push_back(pool->enqueue([&f, begin, end]() {
				for (std::size_t i = begin; i < end; ++i)
					f(i);
			}));
		}
		for (auto& future : futures)
			future.wait();
	}

	struct InsideSegment
	{
		Vec3 origin;
		Vec3 direction;
		Scalar t0;
		Scalar t1;
		Scalar length;
	};

	void append_segments_from_hits(const std::vector<Ray>& rays,
								   const std::vector<std::vector<Scalar>>& ray_hits,
								   const std::vector<uint8_t>& ray_inside0,
								   std::vector<InsideSegment>& segments, Scalar& total_length)
	{
		const Scalar eps = convergence_tolerance();

		for (size_t i = 0; i < rays.size(); ++i)
		{
			const Ray& ray = rays[i];
			if (ray.t_max <= ray.t_min + eps)
				continue;

			std::vector<Scalar> hits = ray_hits[i];
			if (!hits.empty())
				std::sort(hits.begin(), hits.end());

			bool inside = (i < ray_inside0.size() && ray_inside0[i] != 0);
			Scalar prev_t = ray.t_min;
			Scalar last_t = ray.t_min - eps * Scalar(2);
			for (Scalar t_hit : hits)
			{
				Scalar t_clamped = std::max(ray.t_min, std::min(ray.t_max, t_hit));
				if (t_clamped <= last_t + eps)
					continue;
				last_t = t_clamped;
				if (t_clamped <= prev_t + eps)
				{
					inside = !inside;
					prev_t = t_clamped;
					continue;
				}
				if (inside)
				{
					InsideSegment seg{ray.origin, ray.direction, prev_t, t_clamped, t_clamped - prev_t};
					if (seg.length > eps)
					{
						segments.push_back(seg);
						total_length += seg.length;
					}
				}
				inside = !inside;
				prev_t = t_clamped;
			}
			if (inside && ray.t_max > prev_t + eps)
			{
				InsideSegment seg{ray.origin, ray.direction, prev_t, ray.t_max, ray.t_max - prev_t};
				if (seg.length > eps)
				{
					segments.push_back(seg);
					total_length += seg.length;
				}
			}
		}
	}

	bool sample_inside_from_segments(const std::vector<InsideSegment>& segments, Scalar total_length,
									 size_t target_inside_points, SpatialGrid* inside_grid, Scalar grid_cell_size,
									 const Vec3& bbox_min, const Vec3& bbox_max,
									 std::vector<Vec3>& inside_samples, std::mt19937& gen,
									 size_t max_reject_attempts)
	{
		(void)total_length;
		if (inside_samples.size() >= target_inside_points)
			return false;
		if (segments.empty())
			return false;
		const Scalar eps = convergence_tolerance();
		if (total_length <= eps)
			return false;

		std::vector<Scalar> prefix;
		prefix.reserve(segments.size());
		Scalar acc_length = Scalar(0);
		for (const InsideSegment& seg : segments)
		{
			acc_length += seg.length;
			prefix.push_back(acc_length);
		}
		if (acc_length <= eps)
			return false;

		std::uniform_real_distribution<Scalar> uniform(0.0, 1.0);
		size_t reject_count = 0;
		while (inside_samples.size() < target_inside_points && reject_count < max_reject_attempts)
		{
			Scalar r = uniform(gen) * acc_length;
			auto it = std::lower_bound(prefix.begin(), prefix.end(), r);
			size_t idx = static_cast<size_t>(it - prefix.begin());
			if (idx >= segments.size())
				idx = segments.size() - 1;
			const InsideSegment& seg = segments[idx];
			Scalar prev = (idx == 0) ? Scalar(0) : prefix[idx - 1];
			Scalar local = r - prev;
			Scalar t = seg.t0 + local;
			Vec3 pt = seg.origin + seg.direction * t;
			pt = clamp_point(pt, bbox_min, bbox_max);
			if (!inside_grid || inside_grid->is_valid_sample(pt, grid_cell_size, inside_samples))
			{
				inside_samples.push_back(pt);
				if (inside_grid)
					inside_grid->insert(pt, static_cast<uint32>(inside_samples.size() - 1));
				reject_count = 0;
			}
			else
			{
				reject_count++;
			}
		}
		return reject_count >= max_reject_attempts;
	}

	std::vector<Vec3> sample_alpha_level_set_rays_cpu(Traits& traits, size_t target_num_points,
													  SpatialGrid* spatial_grid, Scalar sample_spacing,
													  const Vec3& bbox_min, const Vec3& bbox_max,
													  int max_outer_iterations,
													  const std::vector<Vec3>* existing_samples)
	{
		const bool bounded_target = (target_num_points != unlimited_sample_target_);
		const size_t reserve_hint =
			bounded_target
				? std::max<size_t>(target_num_points + target_num_points / 2, size_t(1024))
				: std::max<size_t>(static_cast<size_t>(std::max(1, params_.batch_size)) *
									   static_cast<size_t>(std::max(1, max_outer_iterations)),
								   size_t(1024));
		std::vector<Vec3> neighbor_samples;
		if (existing_samples)
			neighbor_samples = *existing_samples;
		neighbor_samples.reserve(neighbor_samples.size() + reserve_hint);
		std::vector<Vec3> all_samples;
		all_samples.reserve(reserve_hint);

		std::mt19937 gen(params_.seed);
		const Scalar eps = convergence_tolerance();
		const int max_steps = std::max(1, params_.max_iterations);
		const int check_interval = 10;

		for (int iter = 0; iter < max_outer_iterations && (!bounded_target || all_samples.size() < target_num_points);
			 ++iter)
		{
			auto rays = generate_rays(bbox_min, bbox_max, params_.batch_size, gen);
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
						const Ray& ray = rays[i];
						step_udf[i] = traits.eval_distance(ray.origin + t[i] * ray.direction);
					}
				};
				if (active_count >= 1024)
					parallel_for_rays(static_cast<std::size_t>(R), eval_ray);
				else
					for (std::size_t i = 0; i < static_cast<std::size_t>(R); ++i)
						eval_ray(i);

				// stage 2 (serial): same ordered logic as before, using the precomputed udf
				for (int i = 0; i < R; ++i)
				{
					if (t[i] >= t_max[i])
						continue;

					const Ray& ray = rays[i];
					const Vec3 X = ray.origin + t[i] * ray.direction;
					const Scalar udf = step_udf[i];
					if (!std::isfinite(static_cast<double>(udf)))
					{
						t[i] = t_max[i];
						continue;
					}

					const Scalar delta = std::abs(udf - params_.alpha);
					if (delta < eps)
					{
						Vec3 pt = clamp_point(X, bbox_min, bbox_max);
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
						const Scalar step_size = std::max(delta / params_.step_bound, eps * Scalar(0.5));
						t[i] += step_size;
					}
					if (t[i] >= t_max[i])
						--active_count;
					if (bounded_target && all_samples.size() >= target_num_points)
						break;
				}
				if (bounded_target && all_samples.size() >= target_num_points)
					break;
			}

			std::cout << std::endl;
			std::cout << "  Collected " << iter_samples.size() << " unique samples this iteration" << std::endl;
			if (bounded_target)
				std::cout << "  Total: " << all_samples.size() << " / " << target_num_points << std::endl;
			else
				std::cout << "  Total: " << all_samples.size() << std::endl;
		}

		if (bounded_target && all_samples.size() > target_num_points)
		{
			std::shuffle(all_samples.begin(), all_samples.end(), gen);
			all_samples.resize(target_num_points);
		}

		return all_samples;
	}

	SampleResult sample_alpha_level_set_rays_with_inside_cpu(Traits& traits, size_t target_surface_points,
															 size_t target_inside_points, SpatialGrid* surface_grid,
															 SpatialGrid* inside_grid, Scalar grid_cell_size,
															 const Vec3& bbox_min, const Vec3& bbox_max)
	{
		SampleResult result;
		result.surface_samples.reserve(target_surface_points * 1.5);
		result.inside_samples.reserve(target_inside_points * 1.5);

		std::mt19937 gen(params_.seed);
		const Scalar eps = convergence_tolerance();
		const int max_steps = std::max(1, params_.max_iterations);
		const int check_interval = 10;
		const size_t max_inside_reject = 1000;

		std::vector<InsideSegment> all_segments;
		Scalar segments_total_length = Scalar(0);

		for (int iter = 0; iter < params_.max_outer_iterations; ++iter)
		{
			if (result.surface_samples.size() >= target_surface_points)
				break;
			auto rays = generate_rays(bbox_min, bbox_max, params_.batch_size, gen);
			const int R = static_cast<int>(rays.size());
			if (R == 0)
				break;

			std::cout << "Iteration " << iter << ": Processing " << R << " rays on CPU..." << std::endl;

			std::vector<Scalar> t(static_cast<size_t>(R));
			std::vector<Scalar> t_max(static_cast<size_t>(R));
			std::vector<std::vector<Scalar>> ray_hits(static_cast<size_t>(R));
			std::vector<uint8_t> ray_inside0(static_cast<size_t>(R), 0);
			std::vector<uint8_t> inside_prev(static_cast<size_t>(R), 0);
			std::vector<Scalar> last_hit(static_cast<size_t>(R), -std::numeric_limits<Scalar>::max());
			for (int i = 0; i < R; ++i)
			{
				t[i] = rays[i].t_min;
				t_max[i] = rays[i].t_max;
				const Vec3 X = rays[i].origin + t[i] * rays[i].direction;
				const Scalar udf = traits.eval_distance(X);
				const bool inside = (udf <= params_.alpha);
				ray_inside0[static_cast<size_t>(i)] = inside ? uint8_t(1) : uint8_t(0);
				inside_prev[static_cast<size_t>(i)] = ray_inside0[static_cast<size_t>(i)];
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

				for (int i = 0; i < R; ++i)
				{
					if (t[i] >= t_max[i])
						continue;

					const Ray& ray = rays[i];
					const Vec3 X = ray.origin + t[i] * ray.direction;
					const Scalar udf = traits.eval_distance(X);
					if (!std::isfinite(static_cast<double>(udf)))
					{
						t[i] = t_max[i];
						continue;
					}

					const bool inside = (udf <= params_.alpha);
					if (inside != (inside_prev[static_cast<size_t>(i)] != 0))
					{
						Scalar& last = last_hit[static_cast<size_t>(i)];
						if (t[i] - last > eps)
						{
							ray_hits[static_cast<size_t>(i)].push_back(t[i]);
							last = t[i];
						}
						inside_prev[static_cast<size_t>(i)] = inside ? uint8_t(1) : uint8_t(0);
					}

					const Scalar delta = std::abs(udf - params_.alpha);
					if (delta < eps)
					{
						Vec3 pt = clamp_point(X, bbox_min, bbox_max);
						if (result.surface_samples.size() < target_surface_points)
						{
							if (!surface_grid || surface_grid->is_valid_sample(pt, grid_cell_size, result.surface_samples))
							{
								iter_samples.push_back(pt);
								result.surface_samples.push_back(pt);
								if (surface_grid)
									surface_grid->insert(pt, static_cast<uint32>(result.surface_samples.size() - 1));
							}
						}

						const Scalar step_size = std::max(delta, eps * Scalar(0.5));
						t[i] += step_size;
					}
					else
					{
						const Scalar step_size = std::max(delta / params_.step_bound, eps * Scalar(0.5));
						t[i] += step_size;
					}
					if (t[i] >= t_max[i])
						--active_count;
					if (result.surface_samples.size() >= target_surface_points &&
						result.inside_samples.size() >= target_inside_points)
						break;
				}
				if (result.surface_samples.size() >= target_surface_points &&
					result.inside_samples.size() >= target_inside_points)
					break;
			}

			std::cout << std::endl;
			std::cout << "  Collected " << iter_samples.size() << " surface samples this iteration" << std::endl;
			std::cout << "  Total surface: " << result.surface_samples.size() << " / " << target_surface_points
					  << std::endl;

			if (target_inside_points > 0)
				append_segments_from_hits(rays, ray_hits, ray_inside0, all_segments, segments_total_length);
		}

		if (target_inside_points > 0 && result.inside_samples.size() < target_inside_points)
		{
			sample_inside_from_segments(
				all_segments, segments_total_length, target_inside_points, inside_grid, grid_cell_size, bbox_min,
				bbox_max, result.inside_samples, gen, max_inside_reject);
			std::cout << "  Total inside: " << result.inside_samples.size() << " / " << target_inside_points
					  << std::endl;
		}

		if (result.surface_samples.size() > target_surface_points)
		{
			std::shuffle(result.surface_samples.begin(), result.surface_samples.end(), gen);
			result.surface_samples.resize(target_surface_points);
		}
		if (result.inside_samples.size() > target_inside_points)
		{
			std::shuffle(result.inside_samples.begin(), result.inside_samples.end(), gen);
			result.inside_samples.resize(target_inside_points);
		}

		return result;
	}

	std::vector<Ray> generate_rays(const Vec3& bbox_min, const Vec3& bbox_max, int num_rays, std::mt19937& gen) const
	{
		std::vector<Ray> rays;
		rays.reserve(num_rays);

		Vec3 center = (bbox_min + bbox_max) * 0.5;
		Vec3 half_size = (bbox_max - bbox_min) * 0.5 * (1.0 + params_.bbox_expand);
		Vec3 expanded_min = center - half_size;
		Vec3 expanded_max = center + half_size;
		const Scalar diag = (expanded_max - expanded_min).norm();

		std::uniform_real_distribution<Scalar> uniform(0.0, 1.0);
		auto uniform_signed = [&](Scalar a) {
			return (uniform(gen) * 2.0 - 1.0) * a; // [-a, a]
		};
		for (int i = 0; i < num_rays; ++i)
		{
			Ray ray;

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

				Ray tmp;
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

	static std::pair<Scalar, Scalar> intersect_bbox(const Ray& ray, const Vec3& bbox_min, const Vec3& bbox_max)
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

	Params params_;
};

} // namespace geometry

} // namespace cgogn

#endif // CGOGN_GEOMETRY_TYPES_RAY_LEVEL_SET_SAMPLER_H_
