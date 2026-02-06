#ifndef CGOGN_GEOMETRY_TYPES_RAY_LEVEL_SET_SAMPLER_H_
#define CGOGN_GEOMETRY_TYPES_RAY_LEVEL_SET_SAMPLER_H_

#include <cgogn/core/utils/numerics.h>
#include <cgogn/geometry/types/spatial_grid.h>
#include <cgogn/geometry/types/vector_traits.h>

#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>
#include <random>
#include <utility>
#include <vector>
#include <torch/torch.h>

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
		int max_outer_iterations = 500;
		int seed = 42;
	};

	struct SampleResult
	{
		std::vector<Vec3> surface_samples;
		std::vector<Vec3> inside_samples;
	};

	explicit RayLevelSetSampler(const Params& params, const torch::Device& device)
		: params_(params), device_(device)
	{
	}

	void set_params(const Params& params)
	{
		params_ = params;
	}

	void set_device(const torch::Device& device)
	{
		device_ = device;
	}

	std::vector<Vec3> sample_alpha_level_set_rays(Traits& traits, size_t target_num_points,
												  SpatialGrid* spatial_grid, Scalar grid_cell_size,
												  const Vec3& bbox_min, const Vec3& bbox_max,
												  bool* used_sdf_filter = nullptr)
	{
	
		if constexpr (Traits::kUsesTorch)
		{
			set_device(traits.device());
			return sample_alpha_level_set_rays_gpu(traits, target_num_points, spatial_grid, grid_cell_size, bbox_min,
												   bbox_max, used_sdf_filter);
		}
		else
		{
			if (used_sdf_filter)
				*used_sdf_filter = false;
			return sample_alpha_level_set_rays_cpu(traits, target_num_points, spatial_grid, grid_cell_size, bbox_min,
												   bbox_max);
		}
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

		if constexpr (Traits::kUsesTorch)
		{
			set_device(traits.device());
			return sample_alpha_level_set_rays_with_inside_gpu(traits, target_surface_points, target_inside_points,
															   surface_grid, inside_grid, grid_cell_size, bbox_min,
															   bbox_max, used_sdf_filter);
		}
		else
		{
			if (used_sdf_filter)
				*used_sdf_filter = false;
			return sample_alpha_level_set_rays_with_inside_cpu(traits, target_surface_points, target_inside_points,
															   surface_grid, inside_grid, grid_cell_size, bbox_min,
															   bbox_max);
		}
	}

private:
	static Vec3 clamp_point(const Vec3& p, const Vec3& bbox_min, const Vec3& bbox_max)
	{
		return p.cwiseMax(bbox_min).cwiseMin(bbox_max);
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
		const Scalar eps = std::max(Scalar(1e-6), params_.tol);

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
		const Scalar eps = std::max(Scalar(1e-6), params_.tol);
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

	std::vector<Vec3> sample_alpha_level_set_rays_gpu(Traits& traits, size_t target_num_points,
													  SpatialGrid* spatial_grid, Scalar grid_cell_size,
													  const Vec3& bbox_min, const Vec3& bbox_max,
													  bool* used_sdf_filter)
	{
		std::vector<Vec3> all_samples;
		all_samples.reserve(target_num_points * 1.5);

		std::mt19937 gen(params_.seed);
		const Scalar eps = std::max(Scalar(1e-5), params_.tol);
		const int max_steps = std::max(1, params_.max_iterations);
		const int check_interval = 500;
		int total_rays = 0;

		bool has_sdf_any = false;
		for (int iter = 0; iter < params_.max_outer_iterations && all_samples.size() < target_num_points; ++iter)
		{
			auto rays = generate_rays(bbox_min, bbox_max, params_.batch_size, gen);
			total_rays += static_cast<int>(rays.size());
			const int R = static_cast<int>(rays.size());
			if (R == 0)
				break;

			std::cout << "Iteration " << iter << ": Processing " << R << " rays on GPU..." << std::endl;

			auto buffers = prepare_ray_batch_buffers(rays);

			auto& O = buffers.O;
			auto& D = buffers.D;
			auto& t = buffers.t;
			auto& t_max = buffers.t_max;

			torch::Tensor X = O + t.unsqueeze(1) * D;
			auto d_sdf = traits.eval_values_sdf(X);
			torch::Tensor d = d_sdf.first;
			torch::Tensor sdf = d_sdf.second;
			const bool has_sdf = sdf.defined();
			if (has_sdf)
				has_sdf_any = true;
			if (!d.defined())
				break;

			d = d - params_.alpha;
			torch::Tensor delta = torch::abs(d);

			torch::Tensor not_converged = t < t_max;

			std::vector<Vec3> iter_samples;
			for (int step = 0; step < max_steps; ++step)
			{
				if (step % check_interval == 0)
				{
					int active_count = torch::sum(not_converged).item<int>();
					if (active_count == 0)
						break;
					std::cout << "    Step " << step << ": " << active_count << " active\r" << std::flush;
				}

				torch::Tensor inside_mask;
				if (has_sdf)
					inside_mask = sdf < 0;
				torch::Tensor near_mask = not_converged & (delta < eps);
				if (has_sdf)
					near_mask = near_mask & inside_mask;
				torch::Tensor hit_idx = torch::nonzero(near_mask).squeeze(1);
				if (hit_idx.numel() > 0)
				{
					torch::Tensor hit_X = X.index_select(0, hit_idx).to(torch::kCPU);
					auto hit_acc = hit_X.accessor<float, 2>();

					for (int64_t i = 0; i < hit_X.size(0); ++i)
					{
						Vec3 pt(hit_acc[i][0], hit_acc[i][1], hit_acc[i][2]);
						pt = clamp_point(pt, bbox_min, bbox_max);

						if (!spatial_grid || spatial_grid->is_valid_sample(pt, grid_cell_size, all_samples))
						{
							iter_samples.push_back(pt);
							all_samples.push_back(pt);
							if (spatial_grid)
								spatial_grid->insert(pt, static_cast<uint32>(all_samples.size() - 1));
						}
					}

					// One-step cross to avoid repeated near-surface iterations.
					torch::Tensor p_near = X.index_select(0, hit_idx);
					torch::Tensor delta_near = delta.index_select(0, hit_idx);
					torch::Tensor t_near = t.index_select(0, hit_idx);
					torch::Tensor D_near = D.index_select(0, hit_idx);

					torch::Tensor min_step = torch::full_like(delta_near, eps * 0.5f);
					torch::Tensor step_size = torch::maximum(delta_near, min_step);

					p_near = p_near + step_size.unsqueeze(1) * D_near;
					t_near = t_near + step_size;

					auto d_sdf_near = traits.eval_values_sdf(p_near);
					torch::Tensor d_near = d_sdf_near.first;
					torch::Tensor sdf_near = d_sdf_near.second;
					if (!d_near.defined() || (has_sdf && !sdf_near.defined()))
						break;

					d_near = d_near - params_.alpha;
					delta_near = torch::abs(d_near);

					X.index_copy_(0, hit_idx, p_near);
					t.index_copy_(0, hit_idx, t_near);
					delta.index_copy_(0, hit_idx, delta_near);
					d.index_copy_(0, hit_idx, d_near);
					if (has_sdf)
						sdf.index_copy_(0, hit_idx, sdf_near);
				}

				torch::Tensor march_mask = not_converged & (~near_mask);
				torch::Tensor nc_idx = torch::nonzero(march_mask).squeeze(1);
				if (nc_idx.numel() == 0)
				{
					not_converged = t < t_max;
					if (all_samples.size() >= target_num_points)
						break;
					continue;
				}

				torch::Tensor p_far = X.index_select(0, nc_idx);
				torch::Tensor delta_far = delta.index_select(0, nc_idx);
				torch::Tensor t_far = t.index_select(0, nc_idx);
				torch::Tensor D_far = D.index_select(0, nc_idx);

				torch::Tensor step_size;
				if (has_sdf)
				{
					torch::Tensor sdf_far_current = sdf.index_select(0, nc_idx);
					torch::Tensor inside_far = sdf_far_current < 0;
					torch::Tensor sdf_abs = torch::abs(sdf_far_current);
					torch::Tensor step_raw = torch::where(inside_far, delta_far, sdf_abs);
					torch::Tensor min_step = torch::full_like(step_raw, eps * 0.5f);
					step_size = torch::maximum(step_raw / params_.step_bound, min_step);
					torch::Tensor cross_mask = (~inside_far) & (sdf_abs < eps);
					step_size = torch::where(cross_mask, torch::full_like(step_size, eps), step_size);
				}
				else
				{
					step_size = delta_far / params_.step_bound;
				}
				p_far = p_far + step_size.unsqueeze(1) * D_far;
				t_far = t_far + step_size;

				auto d_sdf_far = traits.eval_values_sdf(p_far);
				torch::Tensor d_far = d_sdf_far.first;
				torch::Tensor sdf_far = d_sdf_far.second;
				if (!d_far.defined() || (has_sdf && !sdf_far.defined()))
					break;

				d_far = d_far - params_.alpha;
				torch::Tensor delta_far_new = torch::abs(d_far);

				X.index_copy_(0, nc_idx, p_far);
				t.index_copy_(0, nc_idx, t_far);
				delta.index_copy_(0, nc_idx, delta_far_new);
				d.index_copy_(0, nc_idx, d_far);
				if (has_sdf)
					sdf.index_copy_(0, nc_idx, sdf_far);

				not_converged = t < t_max;

				if (all_samples.size() >= target_num_points)
					break;
			}

			std::cout << std::endl;
			std::cout << "  Collected " << iter_samples.size() << " unique samples this iteration" << std::endl;
			std::cout << "  Total: " << all_samples.size() << " / " << target_num_points << std::endl;
		}

		if (all_samples.size() > target_num_points)
		{
			std::shuffle(all_samples.begin(), all_samples.end(), gen);
			all_samples.resize(target_num_points);
		}
		if (used_sdf_filter)
			*used_sdf_filter = has_sdf_any;
		return all_samples;
	}

	std::vector<Vec3> sample_alpha_level_set_rays_cpu(Traits& traits, size_t target_num_points,
													  SpatialGrid* spatial_grid, Scalar grid_cell_size,
													  const Vec3& bbox_min, const Vec3& bbox_max)
	{
		std::vector<Vec3> all_samples;
		all_samples.reserve(target_num_points * 1.5);

		std::mt19937 gen(params_.seed);
		const Scalar eps = std::max(Scalar(1e-5), params_.tol);
		const int max_steps = std::max(1, params_.max_iterations);
		const int check_interval = 10;

		for (int iter = 0; iter < params_.max_outer_iterations && all_samples.size() < target_num_points; ++iter)
		{
			auto rays = generate_rays(bbox_min, bbox_max, params_.batch_size, gen);
			const int R = static_cast<int>(rays.size());
			if (R == 0)
				break;

			std::cout << "Iteration " << iter << ": Processing " << R << " rays on CPU..." << std::endl;

			std::vector<Scalar> t(static_cast<size_t>(R));
			std::vector<Scalar> t_max(static_cast<size_t>(R));
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

					const Scalar delta = std::abs(udf - params_.alpha);
					if (delta < eps)
					{
						Vec3 pt = clamp_point(X, bbox_min, bbox_max);
						if (!spatial_grid || spatial_grid->is_valid_sample(pt, grid_cell_size, all_samples))
						{
							iter_samples.push_back(pt);
							all_samples.push_back(pt);
							if (spatial_grid)
								spatial_grid->insert(pt, static_cast<uint32>(all_samples.size() - 1));
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
					if (all_samples.size() >= target_num_points)
						break;
				}
				if (all_samples.size() >= target_num_points)
					break;
			}

			std::cout << std::endl;
			std::cout << "  Collected " << iter_samples.size() << " unique samples this iteration" << std::endl;
			std::cout << "  Total: " << all_samples.size() << " / " << target_num_points << std::endl;
		}

		if (all_samples.size() > target_num_points)
		{
			std::shuffle(all_samples.begin(), all_samples.end(), gen);
			all_samples.resize(target_num_points);
		}

		return all_samples;
	}

	SampleResult sample_alpha_level_set_rays_with_inside_gpu(Traits& traits, size_t target_surface_points,
															 size_t target_inside_points, SpatialGrid* surface_grid,
															 SpatialGrid* inside_grid, Scalar grid_cell_size,
															 const Vec3& bbox_min, const Vec3& bbox_max,
															 bool* used_sdf_filter)
	{
		SampleResult result;
		result.surface_samples.reserve(target_surface_points * 1.5);
		result.inside_samples.reserve(target_inside_points * 1.5);

		std::mt19937 gen(params_.seed);
		const Scalar eps = std::max(Scalar(1e-5), params_.tol);
		const int max_steps = std::max(1, params_.max_iterations);
		const int check_interval = 500;
		const size_t max_inside_reject = 1000;

		std::vector<InsideSegment> all_segments;
		Scalar segments_total_length = Scalar(0);

		bool has_sdf_any = false;
		for (int iter = 0; iter < params_.max_outer_iterations; ++iter)
		{
			if (result.surface_samples.size() >= target_surface_points)
				break;

			auto rays = generate_rays(bbox_min, bbox_max, params_.batch_size, gen);
			const int R = static_cast<int>(rays.size());
			if (R == 0)
				break;

			std::cout << "Iteration " << iter << ": Processing " << R << " rays on GPU..." << std::endl;

			auto buffers = prepare_ray_batch_buffers(rays);

			auto& O = buffers.O;
			auto& D = buffers.D;
			auto& t = buffers.t;
			auto& t_max = buffers.t_max;

			torch::Tensor X = O + t.unsqueeze(1) * D;
			auto d_sdf = traits.eval_values_sdf(X);
			torch::Tensor f = d_sdf.first;
			torch::Tensor sdf = d_sdf.second;
			const bool has_sdf = sdf.defined();
			if (has_sdf)
				has_sdf_any = true;
			if (!f.defined())
				break;

			if (f.dim() == 2 && f.size(1) == 1)
				f = f.squeeze(1);
			if (has_sdf && sdf.dim() == 2 && sdf.size(1) == 1)
				sdf = sdf.squeeze(1);

			torch::Tensor d = f - params_.alpha;
			torch::Tensor delta = torch::abs(d);

			torch::Tensor not_converged = t < t_max;

			std::vector<std::vector<Scalar>> ray_hits;
			std::vector<uint8_t> ray_inside0;
			std::vector<Scalar> last_hit;
			torch::Tensor inside_prev;
			if (target_inside_points > 0)
			{
				ray_hits.resize(static_cast<size_t>(R));
				ray_inside0.resize(static_cast<size_t>(R));
				last_hit.assign(static_cast<size_t>(R), -std::numeric_limits<Scalar>::max());

				torch::Tensor inside_init = d <= 0;
				if (has_sdf)
					inside_init = inside_init & (sdf < 0);
				inside_prev = inside_init.clone();

				torch::Tensor inside_cpu = inside_init.to(torch::kCPU);
				auto inside_acc = inside_cpu.accessor<bool, 1>();
				for (int i = 0; i < R; ++i)
					ray_inside0[static_cast<size_t>(i)] = inside_acc[i] ? uint8_t(1) : uint8_t(0);
			}

			std::vector<Vec3> iter_samples;
			for (int step = 0; step < max_steps; ++step)
			{
				if (step % check_interval == 0)
				{
					int active_count = torch::sum(not_converged).item<int>();
					if (active_count == 0)
						break;
					std::cout << "    Step " << step << ": " << active_count << " active\r" << std::flush;
				}

				torch::Tensor near_mask = not_converged & (delta < eps);

				torch::Tensor near_mask_surface = near_mask;
				if (has_sdf)
				{
					torch::Tensor inside_mask = sdf < 0;
					near_mask_surface = near_mask_surface & inside_mask;
				}

				torch::Tensor hit_idx_surface = torch::nonzero(near_mask_surface).squeeze(1);
				if (hit_idx_surface.numel() > 0 && result.surface_samples.size() < target_surface_points)
				{
					torch::Tensor hit_X = X.index_select(0, hit_idx_surface).to(torch::kCPU);
					auto hit_acc = hit_X.accessor<float, 2>();

					for (int64_t i = 0; i < hit_X.size(0); ++i)
					{
						Vec3 pt(hit_acc[i][0], hit_acc[i][1], hit_acc[i][2]);
						pt = clamp_point(pt, bbox_min, bbox_max);

						if (!surface_grid || surface_grid->is_valid_sample(pt, grid_cell_size, result.surface_samples))
						{
							iter_samples.push_back(pt);
							result.surface_samples.push_back(pt);
							if (surface_grid)
								surface_grid->insert(pt, static_cast<uint32>(result.surface_samples.size() - 1));
						}
					}
				}

				torch::Tensor hit_idx = torch::nonzero(near_mask).squeeze(1);
				if (hit_idx.numel() > 0)
				{
					torch::Tensor p_near = X.index_select(0, hit_idx);
					torch::Tensor delta_near = delta.index_select(0, hit_idx);
					torch::Tensor t_near = t.index_select(0, hit_idx);
					torch::Tensor D_near = D.index_select(0, hit_idx);

					torch::Tensor min_step = torch::full_like(delta_near, eps * 0.5f);
					torch::Tensor step_size = torch::maximum(delta_near, min_step);

					p_near = p_near + step_size.unsqueeze(1) * D_near;
					t_near = t_near + step_size;

					auto d_sdf_near = traits.eval_values_sdf(p_near);
					torch::Tensor f_near = d_sdf_near.first;
					torch::Tensor sdf_near = d_sdf_near.second;
					if (!f_near.defined() || (has_sdf && !sdf_near.defined()))
						break;

					if (f_near.dim() == 2 && f_near.size(1) == 1)
						f_near = f_near.squeeze(1);
					if (has_sdf && sdf_near.dim() == 2 && sdf_near.size(1) == 1)
						sdf_near = sdf_near.squeeze(1);
					torch::Tensor d_near = f_near - params_.alpha;
					delta_near = torch::abs(d_near);

					X.index_copy_(0, hit_idx, p_near);
					t.index_copy_(0, hit_idx, t_near);
					delta.index_copy_(0, hit_idx, delta_near);
					d.index_copy_(0, hit_idx, d_near);
					if (has_sdf)
						sdf.index_copy_(0, hit_idx, sdf_near);
				}

				torch::Tensor march_mask = not_converged & (~near_mask);
				torch::Tensor nc_idx = torch::nonzero(march_mask).squeeze(1);
				if (nc_idx.numel() > 0)
				{
					torch::Tensor p_far = X.index_select(0, nc_idx);
					torch::Tensor delta_far = delta.index_select(0, nc_idx);
					torch::Tensor t_far = t.index_select(0, nc_idx);
					torch::Tensor D_far = D.index_select(0, nc_idx);

					torch::Tensor step_size;
					if (has_sdf)
					{
						torch::Tensor sdf_far_current = sdf.index_select(0, nc_idx);
						torch::Tensor inside_far = sdf_far_current < 0;
						torch::Tensor sdf_abs = torch::abs(sdf_far_current);
						torch::Tensor step_raw = torch::where(inside_far, delta_far, sdf_abs);
						torch::Tensor min_step = torch::full_like(step_raw, eps * 0.5f);
						step_size = torch::maximum(step_raw / params_.step_bound, min_step);
						torch::Tensor cross_mask = (~inside_far) & (sdf_abs < eps);
						step_size = torch::where(cross_mask, torch::full_like(step_size, eps), step_size);
					}
					else
					{
						step_size = delta_far / params_.step_bound;
					}
					p_far = p_far + step_size.unsqueeze(1) * D_far;
					t_far = t_far + step_size;

					auto d_sdf_far = traits.eval_values_sdf(p_far);
					torch::Tensor f_far = d_sdf_far.first;
					torch::Tensor sdf_far = d_sdf_far.second;
					if (!f_far.defined() || (has_sdf && !sdf_far.defined()))
						break;

					if (f_far.dim() == 2 && f_far.size(1) == 1)
						f_far = f_far.squeeze(1);
					if (has_sdf && sdf_far.dim() == 2 && sdf_far.size(1) == 1)
						sdf_far = sdf_far.squeeze(1);
					torch::Tensor d_far = f_far - params_.alpha;
					torch::Tensor delta_far_new = torch::abs(d_far);

					X.index_copy_(0, nc_idx, p_far);
					t.index_copy_(0, nc_idx, t_far);
					delta.index_copy_(0, nc_idx, delta_far_new);
					d.index_copy_(0, nc_idx, d_far);
					if (has_sdf)
						sdf.index_copy_(0, nc_idx, sdf_far);
				}

				not_converged = t < t_max;

				if (target_inside_points > 0)
				{
					torch::Tensor inside_current = d <= 0;
					if (has_sdf)
						inside_current = inside_current & (sdf < 0);
					torch::Tensor change_mask = inside_current != inside_prev;
					torch::Tensor change_idx = torch::nonzero(change_mask & not_converged).squeeze(1);
					if (change_idx.numel() > 0)
					{
						torch::Tensor hit_t = t.index_select(0, change_idx).to(torch::kCPU);
						torch::Tensor hit_idx_cpu = change_idx.to(torch::kCPU);
						auto hit_t_acc = hit_t.accessor<float, 1>();
						auto hit_idx_acc = hit_idx_cpu.accessor<int64_t, 1>();
						for (int64_t i = 0; i < hit_t.size(0); ++i)
						{
							int64_t rid = hit_idx_acc[i];
							if (rid < 0 || rid >= R)
								continue;
							Scalar t_val = static_cast<Scalar>(hit_t_acc[i]);
							Scalar& last = last_hit[static_cast<size_t>(rid)];
							if (t_val - last > eps)
							{
								ray_hits[static_cast<size_t>(rid)].push_back(t_val);
								last = t_val;
							}
						}
					}
					inside_prev = inside_current;
				}

				if (nc_idx.numel() == 0)
				{
					if (result.surface_samples.size() >= target_surface_points &&
						result.inside_samples.size() >= target_inside_points)
						break;
					continue;
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
		if (used_sdf_filter)
			*used_sdf_filter = has_sdf_any;
		return result;
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
		const Scalar eps = std::max(Scalar(1e-5), params_.tol);
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

	struct RayBatchBuffers
	{
		torch::Tensor O;
		torch::Tensor D;
		torch::Tensor t;
		torch::Tensor t_max;
	};

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

	void init_ray_buffers(int R)
	{
		if (R <= 0)
			return;

		const bool need_realloc = (ray_capacity_ < R) || (!ray_origins_.defined()) || (ray_device_ != device_);
		if (!need_realloc)
			return;

		ray_capacity_ = R;
		ray_device_ = device_;

		auto cput_opts = torch::TensorOptions().dtype(torch::kFloat32).device(torch::kCPU);
		if (device_.is_cuda())
			cput_opts = cput_opts.pinned_memory(true);

		ray_origins_cpu_ = torch::empty({ray_capacity_, 3}, cput_opts);
		ray_directions_cpu_ = torch::empty({ray_capacity_, 3}, cput_opts);
		ray_t_cpu_ = torch::empty({ray_capacity_}, cput_opts);
		ray_tmax_cpu_ = torch::empty({ray_capacity_}, cput_opts);

		auto gput_opts = torch::TensorOptions().dtype(torch::kFloat32).device(device_);

		ray_origins_ = torch::empty({ray_capacity_, 3}, gput_opts);
		ray_directions_ = torch::empty({ray_capacity_, 3}, gput_opts);
		ray_t_ = torch::empty({ray_capacity_}, gput_opts);
		ray_tmax_ = torch::empty({ray_capacity_}, gput_opts);
	}

	RayBatchBuffers prepare_ray_batch_buffers(const std::vector<Ray>& rays)
	{
		const int R = static_cast<int>(rays.size());
		init_ray_buffers(R);

		auto O_cpu = ray_origins_cpu_.narrow(0, 0, R);
		auto D_cpu = ray_directions_cpu_.narrow(0, 0, R);
		auto t_cpu = ray_t_cpu_.narrow(0, 0, R);
		auto tmax_cpu = ray_tmax_cpu_.narrow(0, 0, R);

		auto O_acc = O_cpu.accessor<float, 2>();
		auto D_acc = D_cpu.accessor<float, 2>();
		auto t_acc = t_cpu.accessor<float, 1>();
		auto tmax_acc = tmax_cpu.accessor<float, 1>();

		for (int i = 0; i < R; ++i)
		{
			O_acc[i][0] = static_cast<float>(rays[i].origin.x());
			O_acc[i][1] = static_cast<float>(rays[i].origin.y());
			O_acc[i][2] = static_cast<float>(rays[i].origin.z());
			D_acc[i][0] = static_cast<float>(rays[i].direction.x());
			D_acc[i][1] = static_cast<float>(rays[i].direction.y());
			D_acc[i][2] = static_cast<float>(rays[i].direction.z());
			t_acc[i] = static_cast<float>(rays[i].t_min);
			tmax_acc[i] = static_cast<float>(rays[i].t_max);
		}

		RayBatchBuffers buffers;
		buffers.O = ray_origins_.narrow(0, 0, R);
		buffers.D = ray_directions_.narrow(0, 0, R);
		buffers.t = ray_t_.narrow(0, 0, R);
		buffers.t_max = ray_tmax_.narrow(0, 0, R);

		const bool nonblocking = device_.is_cuda();
		buffers.O.copy_(O_cpu, nonblocking);
		buffers.D.copy_(D_cpu, nonblocking);
		buffers.t.copy_(t_cpu, nonblocking);
		buffers.t_max.copy_(tmax_cpu, nonblocking);

		return buffers;
	}

	Params params_;
	torch::Device device_;
	int ray_capacity_ = 0;
	torch::Device ray_device_ = torch::kCPU;

	torch::Tensor ray_origins_cpu_;
	torch::Tensor ray_directions_cpu_;
	torch::Tensor ray_t_cpu_;
	torch::Tensor ray_tmax_cpu_;

	torch::Tensor ray_origins_;
	torch::Tensor ray_directions_;
	torch::Tensor ray_t_;
	torch::Tensor ray_tmax_;
};

} // namespace geometry

} // namespace cgogn

#endif // CGOGN_GEOMETRY_TYPES_RAY_LEVEL_SET_SAMPLER_H_
