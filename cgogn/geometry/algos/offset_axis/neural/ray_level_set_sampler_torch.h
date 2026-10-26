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
#ifndef CGOGN_GEOMETRY_ALGOS_OFFSET_AXIS_NEURAL_RAY_LEVEL_SET_SAMPLER_TORCH_H_
#define CGOGN_GEOMETRY_ALGOS_OFFSET_AXIS_NEURAL_RAY_LEVEL_SET_SAMPLER_TORCH_H_

#include <cgogn/geometry/algos/offset_axis/ray_level_set_sampler.h>

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

namespace offset_axis
{

// GPU (LibTorch) version of RayLevelSetSampler: the rays are marched in batches on the device of the traits. The
// public entry points below hide the CPU ones of the base class; the base members used by the GPU code are imported
// with using-declarations.
template <typename Traits>
class TorchRayLevelSetSampler : public RayLevelSetSampler<Traits>
{
	using Base = RayLevelSetSampler<Traits>;
	using Base::clamp_point;
	using Base::convergence_tolerance;
	using Base::generate_rays;
	using Base::params_;
	using typename Base::Ray;

public:
	using typename Base::Params;

	TorchRayLevelSetSampler(const Params& params, const torch::Device& device) : Base(params), device_(device)
	{
	}

	void set_device(const torch::Device& device)
	{
		device_ = device;
	}

	std::vector<Vec3> sample_alpha_level_set_rays_fixed_iterations(Traits& traits, int outer_iterations,
																   SpatialGrid* spatial_grid, Scalar sample_spacing,
																   const Vec3& bbox_min, const Vec3& bbox_max,
																   const std::vector<Vec3>* existing_samples = nullptr)
	{
		const int clamped_outer_iterations = std::max(1, outer_iterations);
		set_device(traits.device());
		return sample_alpha_level_set_rays_gpu(traits, spatial_grid, sample_spacing, bbox_min, bbox_max,
											   clamped_outer_iterations, existing_samples);
	}

private:
	struct RayBatchBuffers
	{
		torch::Tensor O;
		torch::Tensor D;
		torch::Tensor t;
		torch::Tensor t_max;
	};

	std::vector<Vec3> sample_alpha_level_set_rays_gpu(Traits& traits, SpatialGrid* spatial_grid,
													  Scalar sample_spacing, const Vec3& bbox_min,
													  const Vec3& bbox_max, int max_outer_iterations,
													  const std::vector<Vec3>* existing_samples)
	{
		const size_t reserve_hint = std::max<size_t>(static_cast<size_t>(std::max(1, params_.batch_size)) *
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

		for (int iter = 0; iter < max_outer_iterations; ++iter)
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
			}

			std::cout << std::endl;
			std::cout << "  Collected " << iter_samples.size() << " unique samples this iteration" << std::endl;
			std::cout << "  Total: " << all_samples.size() << std::endl;
		}

		return all_samples;
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

} // namespace offset_axis

} // namespace geometry

} // namespace cgogn

#endif // CGOGN_GEOMETRY_ALGOS_OFFSET_AXIS_NEURAL_RAY_LEVEL_SET_SAMPLER_TORCH_H_
