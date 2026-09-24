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
#ifndef CGOGN_GEOMETRY_ALGOS_UDF_NEURAL_ALPHA_PROJECTION_H_
#define CGOGN_GEOMETRY_ALGOS_UDF_NEURAL_ALPHA_PROJECTION_H_

#include <cgogn/geometry/algos/udf/alpha_projection.h>
#include <cgogn/geometry/types/neural_field_forward.h>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
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
	torch::Tensor projected_cpu_;
	torch::Tensor grad_cpu_;
	torch::Tensor sdf_cpu_;
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
	workspace.projected_cpu_ = torch::empty({static_cast<int64_t>(new_capacity), 3}, cpu_opts);
	workspace.grad_cpu_ = torch::empty({static_cast<int64_t>(new_capacity), 3}, cpu_opts);
	workspace.sdf_cpu_ = torch::empty({static_cast<int64_t>(new_capacity)}, cpu_opts);
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
		auto x_acc = cpu_points.accessor<float, 2>();
		for (size_t i = 0; i < nb_points; ++i)
		{
			x_acc[static_cast<long>(i)][0] = static_cast<float>(points[i].x());
			x_acc[static_cast<long>(i)][1] = static_cast<float>(points[i].y());
			x_acc[static_cast<long>(i)][2] = static_cast<float>(points[i].z());
		}

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
			torch::Tensor f = values;
			if (f.dim() == 2 && f.size(1) == 1)
				f = f.squeeze(1);
			torch::Tensor residual = f - parameters.alpha;
			torch::Tensor next_active = active & (torch::abs(residual) > parameters.tolerance);
			active.copy_(next_active);
			if (torch::sum(active).item<int64_t>() == 0)
				break;

			torch::Tensor grad_norm2 = (grad * grad).sum(1, true);
			torch::Tensor step = residual.unsqueeze(1) / (grad_norm2 + 1e-12f);
			step = torch::where(active.unsqueeze(1), step, torch::zeros_like(step));
			X.sub_(step * grad);
			X.copy_(torch::maximum(torch::minimum(X, bbox_max), bbox_min));
		}

		auto [values, grad] = udf.forward_values_grad_gpu(X);
		if (!values.defined() || !grad.defined())
			return false;
		torch::Tensor projected_cpu = workspace.projected_cpu_.narrow(0, 0, count);
		torch::Tensor grad_cpu = workspace.grad_cpu_.narrow(0, 0, count);
		projected_cpu.copy_(X);
		grad_cpu.copy_(grad);
		auto [values_sdf, sdf] = udf.forward_values_sdf_gpu(X);
		(void)values_sdf;
		torch::Tensor sdf_cpu;
		if (sdf.defined())
		{
			if (sdf.dim() == 2 && sdf.size(1) == 1)
				sdf = sdf.squeeze(1);
			sdf_cpu = workspace.sdf_cpu_.narrow(0, 0, count);
			sdf_cpu.copy_(sdf);
		}

		auto projected_acc = projected_cpu.accessor<float, 2>();
		auto grad_acc = grad_cpu.accessor<float, 2>();
		if (sdf_cpu.defined() && sdf_cpu.dim() == 1 && sdf_cpu.size(0) == count)
		{
			auto sdf_acc = sdf_cpu.accessor<float, 1>();
			for (size_t i = 0; i < nb_points; ++i)
			{
				result.keep_mask[i] = ((sdf_acc[static_cast<long>(i)] <= 0.0f) || !filter_positive) ? uint8_t(1)
																 : uint8_t(0);
				result.projected_points[i] =
					Vec3(projected_acc[static_cast<long>(i)][0], projected_acc[static_cast<long>(i)][1],
						 static_cast<Scalar>(projected_acc[static_cast<long>(i)][2]));
				Vec3 normal(static_cast<Scalar>(grad_acc[static_cast<long>(i)][0]),
						 static_cast<Scalar>(grad_acc[static_cast<long>(i)][1]),
						 static_cast<Scalar>(grad_acc[static_cast<long>(i)][2]));
				if (normal.squaredNorm() < Scalar(1e-12))
					normal = Vec3(0, 0, 1);
				else
					normal.normalize();
				result.normals[i] = normal;
			}
		}
		else
		{
			for (size_t i = 0; i < nb_points; ++i)
			{
				result.projected_points[i] =
					Vec3(projected_acc[static_cast<long>(i)][0], projected_acc[static_cast<long>(i)][1],
						 static_cast<Scalar>(projected_acc[static_cast<long>(i)][2]));
				Vec3 normal(static_cast<Scalar>(grad_acc[static_cast<long>(i)][0]),
						 static_cast<Scalar>(grad_acc[static_cast<long>(i)][1]),
						 static_cast<Scalar>(grad_acc[static_cast<long>(i)][2]));
				if (normal.squaredNorm() < Scalar(1e-12))
					normal = Vec3(0, 0, 1);
				else
					normal.normalize();
				result.normals[i] = normal;
			}
		}
		return true;
	}
	catch (const c10::Error&)
	{
		return false;
	}
}

} // namespace geometry

} // namespace cgogn

#endif // CGOGN_GEOMETRY_ALGOS_UDF_NEURAL_ALPHA_PROJECTION_H_
