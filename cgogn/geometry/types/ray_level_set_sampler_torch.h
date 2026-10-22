#ifndef CGOGN_GEOMETRY_TYPES_RAY_LEVEL_SET_SAMPLER_TORCH_H_
#define CGOGN_GEOMETRY_TYPES_RAY_LEVEL_SET_SAMPLER_TORCH_H_

#include <cgogn/geometry/types/ray_level_set_sampler.h>

#include <torch/torch.h>

#include <algorithm>
#include <cstdint>
#include <iostream>
#include <limits>
#include <random>
#include <utility>
#include <vector>

namespace cgogn
{

namespace geometry
{

// GPU (LibTorch) version of RayLevelSetSampler: the rays are marched in batches on the device of the traits. The
// public entry points below hide the CPU ones of the base class; the base members used by the GPU code are imported
// with using-declarations.
template <typename Traits>
class TorchRayLevelSetSampler : public RayLevelSetSampler<Traits>
{
	using Base = RayLevelSetSampler<Traits>;
	using Base::append_segments_from_hits;
	using Base::clamp_point;
	using Base::convergence_tolerance;
	using Base::generate_rays;
	using Base::params_;
	using Base::sample_inside_from_segments;
	using Base::unlimited_sample_target_;
	using typename Base::InsideSegment;
	using typename Base::Ray;

public:
	using typename Base::Params;
	using typename Base::SampleResult;

	TorchRayLevelSetSampler(const Params& params, const torch::Device& device) : Base(params), device_(device)
	{
	}

	void set_device(const torch::Device& device)
	{
		device_ = device;
	}

	std::vector<Vec3> sample_alpha_level_set_rays(Traits& traits, size_t target_num_points, SpatialGrid* spatial_grid,
												  Scalar sample_spacing, const Vec3& bbox_min, const Vec3& bbox_max,
												  bool* used_sdf_filter = nullptr)
	{
		set_device(traits.device());
		return sample_alpha_level_set_rays_gpu(traits, target_num_points, spatial_grid, sample_spacing, bbox_min,
											   bbox_max, used_sdf_filter, params_.max_outer_iterations, nullptr);
	}

	std::vector<Vec3> sample_alpha_level_set_rays_fixed_iterations(Traits& traits, int outer_iterations,
																   SpatialGrid* spatial_grid, Scalar sample_spacing,
																   const Vec3& bbox_min, const Vec3& bbox_max,
																   const std::vector<Vec3>* existing_samples = nullptr,
																   bool* used_sdf_filter = nullptr)
	{
		const int clamped_outer_iterations = std::max(1, outer_iterations);
		set_device(traits.device());
		return sample_alpha_level_set_rays_gpu(traits, unlimited_sample_target_, spatial_grid, sample_spacing, bbox_min,
											   bbox_max, used_sdf_filter, clamped_outer_iterations, existing_samples);
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
																 grid_cell_size, bbox_min, bbox_max, used_sdf_filter);
			return result;
		}

		set_device(traits.device());
		return sample_alpha_level_set_rays_with_inside_gpu(traits, target_surface_points, target_inside_points,
														   surface_grid, inside_grid, grid_cell_size, bbox_min,
														   bbox_max, used_sdf_filter);
	}

private:
	struct RayBatchBuffers
	{
		torch::Tensor O;
		torch::Tensor D;
		torch::Tensor t;
		torch::Tensor t_max;
	};

	std::vector<Vec3> sample_alpha_level_set_rays_gpu(Traits& traits, size_t target_num_points,
													  SpatialGrid* spatial_grid, Scalar sample_spacing,
													  const Vec3& bbox_min, const Vec3& bbox_max, bool* used_sdf_filter,
													  int max_outer_iterations,
													  const std::vector<Vec3>* existing_samples)
	{
		const bool bounded_target = (target_num_points != unlimited_sample_target_);
		const size_t reserve_hint = bounded_target
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
		const int check_interval = 500;
		int total_rays = 0;

		bool has_sdf_any = false;
		for (int iter = 0; iter < max_outer_iterations && (!bounded_target || all_samples.size() < target_num_points);
			 ++iter)
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
				// the last rays run nearly tangent to the level set and creep forward by tiny steps: they take most
				// of the marching time but find few seeds, which the Poisson expansion reaches anyway
				if (params_.min_active_rays > 0 && torch::sum(not_converged).item<int>() < params_.min_active_rays)
					break;

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

						if (!spatial_grid || spatial_grid->is_valid_sample(pt, sample_spacing, neighbor_samples))
						{
							iter_samples.push_back(pt);
							all_samples.push_back(pt);
							neighbor_samples.push_back(pt);
							if (spatial_grid)
								spatial_grid->insert(pt, static_cast<uint32>(neighbor_samples.size() - 1));
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
					if (bounded_target && all_samples.size() >= target_num_points)
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
		if (used_sdf_filter)
			*used_sdf_filter = has_sdf_any;
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
		const Scalar eps = convergence_tolerance();
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
			sample_inside_from_segments(all_segments, segments_total_length, target_inside_points, inside_grid,
										grid_cell_size, bbox_min, bbox_max, result.inside_samples, gen,
										max_inside_reject);
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
		// the page-locked buffers below are uploaded with non_blocking copies and reused by every batch: wait until
		// the previous batch's copies are done before overwriting them
		if (device_.is_cuda())
			torch::cuda::synchronize();

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

#endif // CGOGN_GEOMETRY_TYPES_RAY_LEVEL_SET_SAMPLER_TORCH_H_
