/*******************************************************************************
 * CGoGN: Combinatorial and Geometric modeling with Generic N-dimensional Maps  *
 * Copyright (C), IGG Group, ICube, University of Strasbourg, France            *
 *                                                                              *
 * This library is free software; you can redistribute it and/or modify it      *
 * under the terms of the GNU Lesser General Public License as published by the *
 * Free Software Foundation; either version 2.1 of the License, or (at your     *
 * option) any later version.                                                   *
 *                                                                              *
 * This library is distributed in the hope that it will be useful, but WITHOUT  *
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or        *
 * FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License  *
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
#ifndef CGOGN_GEOMETRY_ALGOS_UDF_NEURAL_NEURAL_ALPHA_PROJECTION_H_
#define CGOGN_GEOMETRY_ALGOS_UDF_NEURAL_NEURAL_ALPHA_PROJECTION_H_

#include <cgogn/geometry/algos/udf/alpha_projection.h>
#include <cgogn/geometry/algos/udf/normal_utils.h>
#include <cgogn/geometry/algos/udf/neural/neural_field_forward.h>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <iostream>
#include <vector>

namespace cgogn
{

namespace geometry
{

struct NeuralProjectionWorkspace
{
	size_t capacity_ = 0;
	torch::Tensor points_cpu_;
	torch::Tensor points_gpu_;
	torch::Tensor active_gpu_;
	torch::Tensor grad_cpu_;
	// projected xyz, gradient, value and (MF only) sdf per point, so a projection needs one device-to-host copy
	torch::Tensor packed_cpu_;
	torch::Tensor clamp_min_gpu_;
	torch::Tensor clamp_max_gpu_;
	Vec3 clamp_min_cache_ = Vec3(0, 0, 0);
	Vec3 clamp_max_cache_ = Vec3(0, 0, 0);
	bool clamp_bounds_valid_ = false;
};

namespace detail
{

inline void ensure_neural_projection_workspace(NeuralProjectionWorkspace& workspace, size_t capacity,
													  const torch::Device& device)
{
	const bool device_changed = workspace.points_gpu_.defined() && workspace.points_gpu_.device() != device;
	if (workspace.capacity_ >= capacity && workspace.points_gpu_.defined() && !device_changed)
		return;

	size_t new_capacity = workspace.capacity_;
	if (new_capacity < capacity)
		new_capacity = std::max(capacity, new_capacity * 2);
	auto cpu_opts = torch::TensorOptions().dtype(torch::kFloat32).device(torch::kCPU);
	if (device.is_cuda())
		cpu_opts = cpu_opts.pinned_memory(true);
	auto gpu_float_opts = torch::TensorOptions().dtype(torch::kFloat32).device(device);
	auto gpu_bool_opts = torch::TensorOptions().dtype(torch::kBool).device(device);

	workspace.points_cpu_ = torch::empty({static_cast<int64_t>(new_capacity), 3}, cpu_opts);
	workspace.points_gpu_ = torch::empty({static_cast<int64_t>(new_capacity), 3}, gpu_float_opts);
	workspace.active_gpu_ = torch::empty({static_cast<int64_t>(new_capacity)}, gpu_bool_opts);
	workspace.grad_cpu_ = torch::empty({static_cast<int64_t>(new_capacity), 3}, cpu_opts);
	workspace.packed_cpu_ = torch::empty({static_cast<int64_t>(new_capacity) * 8}, cpu_opts);
	workspace.clamp_bounds_valid_ = false;
	workspace.capacity_ = new_capacity;
}

inline void ensure_neural_projection_clamp_tensors(NeuralProjectionWorkspace& workspace, const Vec3& bbox_min,
													 const Vec3& bbox_max, const torch::Device& device)
{
	if (workspace.clamp_bounds_valid_ && workspace.clamp_min_cache_.isApprox(bbox_min) &&
		workspace.clamp_max_cache_.isApprox(bbox_max))
		return;

	auto gpu_opts = torch::TensorOptions().dtype(torch::kFloat32).device(device);
	workspace.clamp_min_gpu_ = torch::tensor(
		{static_cast<float>(bbox_min.x()), static_cast<float>(bbox_min.y()), static_cast<float>(bbox_min.z())}, gpu_opts);
	workspace.clamp_max_gpu_ = torch::tensor(
		{static_cast<float>(bbox_max.x()), static_cast<float>(bbox_max.y()), static_cast<float>(bbox_max.z())}, gpu_opts);
	workspace.clamp_min_cache_ = bbox_min;
	workspace.clamp_max_cache_ = bbox_max;
	workspace.clamp_bounds_valid_ = true;
}

// writes points[begin, begin + count) as float xyz rows into dst (a CPU tensor with at least count rows)
inline void fill_points_tensor(const std::vector<Vec3>& points, size_t begin, size_t count, torch::Tensor& dst)
{
	auto acc = dst.accessor<float, 2>();
	for (size_t i = 0; i < count; ++i)
	{
		const Vec3& p = points[begin + i];
		acc[static_cast<long>(i)][0] = static_cast<float>(p.x());
		acc[static_cast<long>(i)][1] = static_cast<float>(p.y());
		acc[static_cast<long>(i)][2] = static_cast<float>(p.z());
	}
}

// points[begin, begin + count) in a new CPU tensor, page-locked when the field runs on CUDA (faster upload)
inline torch::Tensor make_points_tensor(const std::vector<Vec3>& points, size_t begin, size_t count,
										const torch::Device& device)
{
	auto opts = torch::TensorOptions().dtype(torch::kFloat32).device(torch::kCPU);
	if (device.is_cuda())
		opts = opts.pinned_memory(true);
	torch::Tensor t = torch::empty({static_cast<int64_t>(count), 3}, opts);
	fill_points_tensor(points, begin, count, t);
	return t;
}

} // namespace detail

inline bool project_points_to_alpha(const NeuralFieldForward& udf, const std::vector<Vec3>& points,
										const AlphaProjectionParameters& parameters, bool filter_positive,
										NeuralProjectionWorkspace& workspace, AlphaProjectionResult& result)
{
	const size_t nb_points = points.size();
	result.projected_points.assign(nb_points, Vec3(0, 0, 0));
	result.normals.assign(nb_points, Vec3(0, 0, 1));
	result.keep_mask.assign(nb_points, uint8_t(1));
	if (nb_points == 0)
		return true;

	detail::ensure_neural_projection_workspace(workspace, nb_points, udf.device());
	detail::ensure_neural_projection_clamp_tensors(workspace, parameters.bbox_min, parameters.bbox_max, udf.device());
	try
	{
		const int64_t count = static_cast<int64_t>(nb_points);
		torch::Tensor cpu_points = workspace.points_cpu_.narrow(0, 0, count);
		detail::fill_points_tensor(points, 0, nb_points, cpu_points);

		torch::Tensor X = workspace.points_gpu_.narrow(0, 0, count);
		X.copy_(cpu_points);
		torch::Tensor active = workspace.active_gpu_.narrow(0, 0, count);
		active.fill_(true);
		const torch::Tensor& bbox_min = workspace.clamp_min_gpu_;
		const torch::Tensor& bbox_max = workspace.clamp_max_gpu_;

		for (int iter = 0; iter < parameters.max_iterations; ++iter)
		{
			auto [values, grad] = udf.forward_values_grad_gpu(X);
			if (!values.defined() || !grad.defined())
				return false;
			torch::Tensor residual = values - parameters.alpha;
			torch::Tensor next_active = active & (torch::abs(residual) > parameters.tolerance);
			active.copy_(next_active);
			// the early exit costs a host sync; on the last iteration the same skip is done on the device below
			const bool last_iteration = iter + 1 == parameters.max_iterations;
			if (!last_iteration && torch::sum(active).item<int64_t>() == 0)
				break;

			torch::Tensor grad_norm2 = (grad * grad).sum(1, true);
			torch::Tensor step = residual.unsqueeze(1) / (grad_norm2 + 1e-12f);
			step = torch::where(active.unsqueeze(1), step, torch::zeros_like(step));
			torch::Tensor moved = torch::maximum(torch::minimum(X - step * grad, bbox_max), bbox_min);
			X.copy_(last_iteration ? torch::where(active.any(), moved, X) : moved);
		}

		torch::Tensor sdf;
		auto [values, grad] = udf.forward_values_grad_gpu(X, &sdf);
		if (!values.defined() || !grad.defined() || values.numel() != count)
			return false;
		const bool has_sdf = sdf.defined() && sdf.numel() == count;
		std::vector<torch::Tensor> columns = {X, grad.to(torch::kFloat32), values.to(torch::kFloat32).reshape({count, 1})};
		if (has_sdf)
			columns.push_back(sdf.to(torch::kFloat32).reshape({count, 1}));
		torch::Tensor packed_gpu = torch::cat(columns, 1);
		torch::Tensor packed_cpu =
			workspace.packed_cpu_.narrow(0, 0, count * packed_gpu.size(1)).view({count, packed_gpu.size(1)});
		packed_cpu.copy_(packed_gpu);

		// columns: projected xyz, gradient xyz, value, sdf; residual = |udf(x) - alpha| at the projected position
		auto packed_acc = packed_cpu.accessor<float, 2>();
		result.residuals.resize(nb_points);
		for (size_t i = 0; i < nb_points; ++i)
		{
			const long k = static_cast<long>(i);
			if (has_sdf)
				result.keep_mask[i] = ((packed_acc[k][7] <= 0.0f) || !filter_positive) ? uint8_t(1) : uint8_t(0);
			result.projected_points[i] = Vec3(packed_acc[k][0], packed_acc[k][1], static_cast<Scalar>(packed_acc[k][2]));
			result.residuals[i] = std::abs(static_cast<Scalar>(packed_acc[k][6]) - parameters.alpha);
			Vec3 normal(static_cast<Scalar>(packed_acc[k][3]), static_cast<Scalar>(packed_acc[k][4]),
						static_cast<Scalar>(packed_acc[k][5]));
			result.normals[i] = unit_normal_or_z(normal);
		}
		return true;
	}
	catch (const c10::Error& e)
	{
		std::cerr << "Neural alpha projection failed: " << e.what() << std::endl;
		return false;
	}
}

} // namespace geometry

} // namespace cgogn

#endif // CGOGN_GEOMETRY_ALGOS_UDF_NEURAL_NEURAL_ALPHA_PROJECTION_H_
