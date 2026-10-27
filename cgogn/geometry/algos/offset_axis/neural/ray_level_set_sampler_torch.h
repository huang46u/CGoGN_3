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

#include <cgogn/geometry/algos/offset_axis/neural/neural_ray_traits.h>
#include <cgogn/geometry/algos/offset_axis/ray_level_set_sampler.h>

#include <torch/torch.h>

#include <algorithm>
#include <cstdint>
#include <iostream>
#include <limits>
#include <random>
#include <utility>
#include <optional>
#include <vector>

namespace cgogn
{

namespace geometry
{

namespace offset_axis
{

// Device and reusable buffers of the GPU level set sampler. The neural field owns one of them and points the traits to
// it (NeuralFieldRayTraits::ray_buffers), so the buffers are kept between the sampling calls.
struct TorchRayBuffers
{
	torch::Device device = torch::kCPU;
	int ray_capacity = 0;
	torch::Device ray_device = torch::kCPU;

	torch::Tensor ray_origins_cpu;
	torch::Tensor ray_directions_cpu;
	torch::Tensor ray_t_cpu;
	torch::Tensor ray_tmax_cpu;

	torch::Tensor ray_origins;
	torch::Tensor ray_directions;
	torch::Tensor ray_t;
	torch::Tensor ray_tmax;
};

// The rays of one batch, as device tensors (views into the buffers)
struct TorchRayBatch
{
	torch::Tensor O;
	torch::Tensor D;
	torch::Tensor t;
	torch::Tensor t_max;
};

inline void init_ray_buffers(TorchRayBuffers& buffers, int R)
{
	const bool need_realloc = buffers.ray_capacity < R || buffers.ray_device != buffers.device;
	if (!need_realloc)
		return;

	buffers.ray_capacity = R;
	buffers.ray_device = buffers.device;

	auto cput_opts = torch::TensorOptions().dtype(torch::kFloat32).device(torch::kCPU);
	if (buffers.device.is_cuda())
		cput_opts = cput_opts.pinned_memory(true);

	buffers.ray_origins_cpu = torch::empty({buffers.ray_capacity, 3}, cput_opts);
	buffers.ray_directions_cpu = torch::empty({buffers.ray_capacity, 3}, cput_opts);
	buffers.ray_t_cpu = torch::empty({buffers.ray_capacity}, cput_opts);
	buffers.ray_tmax_cpu = torch::empty({buffers.ray_capacity}, cput_opts);

	auto gput_opts = torch::TensorOptions().dtype(torch::kFloat32).device(buffers.device);

	buffers.ray_origins = torch::empty({buffers.ray_capacity, 3}, gput_opts);
	buffers.ray_directions = torch::empty({buffers.ray_capacity, 3}, gput_opts);
	buffers.ray_t = torch::empty({buffers.ray_capacity}, gput_opts);
	buffers.ray_tmax = torch::empty({buffers.ray_capacity}, gput_opts);
}

inline TorchRayBatch prepare_ray_batch(TorchRayBuffers& buffers, const std::vector<LevelSetRay>& rays)
{
	const int R = static_cast<int>(rays.size());
	init_ray_buffers(buffers, R);
	// the page-locked buffers below are uploaded with non_blocking copies and reused by every batch: wait until
	// the previous batch's copies are done before overwriting them
	if (buffers.device.is_cuda())
		torch::cuda::synchronize();

	auto O_cpu = buffers.ray_origins_cpu.narrow(0, 0, R);
	auto D_cpu = buffers.ray_directions_cpu.narrow(0, 0, R);
	auto t_cpu = buffers.ray_t_cpu.narrow(0, 0, R);
	auto tmax_cpu = buffers.ray_tmax_cpu.narrow(0, 0, R);

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

	TorchRayBatch batch;
	batch.O = buffers.ray_origins.narrow(0, 0, R);
	batch.D = buffers.ray_directions.narrow(0, 0, R);
	batch.t = buffers.ray_t.narrow(0, 0, R);
	batch.t_max = buffers.ray_tmax.narrow(0, 0, R);

	const bool nonblocking = buffers.device.is_cuda();
	batch.O.copy_(O_cpu, nonblocking);
	batch.D.copy_(D_cpu, nonblocking);
	batch.t.copy_(t_cpu, nonblocking);
	batch.t_max.copy_(tmax_cpu, nonblocking);

	return batch;
}

// GPU level set sampler: the rays are marched in batches on the device of the traits (the buffers are the ones of the
// field, traits.ray_buffers). Same contract as the CPU sample_level_set_rays; this non-template overload is preferred
// for NeuralFieldRayTraits.
inline std::optional<std::vector<Vec3>> sample_level_set_rays(const NeuralFieldRayTraits& traits, const RayLevelSetParams& params,
											   int outer_iterations, SpatialGrid* spatial_grid, Scalar sample_spacing,
											   const Vec3& bbox_min, const Vec3& bbox_max,
											   const std::vector<Vec3>* existing_samples = nullptr)
{
	TorchRayBuffers& ray_buffers = *traits.ray_buffers;

	const size_t reserve_hint =
		std::max(static_cast<size_t>(params.batch_size) * static_cast<size_t>(outer_iterations), size_t(1024));
	std::vector<Vec3> neighbor_samples;
	if (existing_samples)
		neighbor_samples = *existing_samples;
	neighbor_samples.reserve(neighbor_samples.size() + reserve_hint);
	std::vector<Vec3> all_samples;
	all_samples.reserve(reserve_hint);

	std::mt19937 gen(params.seed);
	const Scalar eps = params.tol;
	const int check_interval = 500;

	for (int iter = 0; iter < outer_iterations; ++iter)
	{
		auto rays = detail::generate_rays(params, bbox_min, bbox_max, params.batch_size, gen);
		const int R = static_cast<int>(rays.size());

		std::cout << "Iteration " << iter << ": Processing " << R << " rays on GPU..." << std::endl;

		auto batch = prepare_ray_batch(ray_buffers, rays);

		auto& O = batch.O;
		auto& D = batch.D;
		auto& t = batch.t;
		auto& t_max = batch.t_max;

		torch::Tensor X = O + t.unsqueeze(1) * D;
		auto d_sdf = traits.eval_values_sdf(X);
		torch::Tensor d = d_sdf.first;
		torch::Tensor sdf = d_sdf.second;
		const bool has_sdf = sdf.defined();
		if (!d.defined())
			return std::nullopt;

		d = d - params.alpha;
		torch::Tensor delta = torch::abs(d);

		torch::Tensor not_converged = t < t_max;

		std::vector<Vec3> iter_samples;
		for (int step = 0; step < params.max_iterations; ++step)
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
			if (params.min_active_rays > 0 && torch::sum(not_converged).item<int>() < params.min_active_rays)
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
					pt = detail::clamp_point(pt, bbox_min, bbox_max);

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
					return std::nullopt;

				d_near = d_near - params.alpha;
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
				step_size = torch::maximum(step_raw / params.step_bound, min_step);
				torch::Tensor cross_mask = (~inside_far) & (sdf_abs < eps);
				step_size = torch::where(cross_mask, torch::full_like(step_size, eps), step_size);
			}
			else
			{
				step_size = delta_far / params.step_bound;
			}
			p_far = p_far + step_size.unsqueeze(1) * D_far;
			t_far = t_far + step_size;

			auto d_sdf_far = traits.eval_values_sdf(p_far);
			torch::Tensor d_far = d_sdf_far.first;
			torch::Tensor sdf_far = d_sdf_far.second;
			if (!d_far.defined() || (has_sdf && !sdf_far.defined()))
				return std::nullopt;

			d_far = d_far - params.alpha;
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

} // namespace offset_axis

} // namespace geometry

} // namespace cgogn

#endif // CGOGN_GEOMETRY_ALGOS_OFFSET_AXIS_NEURAL_RAY_LEVEL_SET_SAMPLER_TORCH_H_
