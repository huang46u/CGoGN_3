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
#ifndef CGOGN_GEOMETRY_ALGOS_OFFSET_AXIS_NEURAL_NEURAL_FIELD_QUERY_H_
#define CGOGN_GEOMETRY_ALGOS_OFFSET_AXIS_NEURAL_NEURAL_FIELD_QUERY_H_

#include <cgogn/geometry/algos/offset_axis/neural/neural_field_forward.h>
#include <cgogn/geometry/algos/offset_axis/neural/neural_alpha_projection.h>
#include <cgogn/geometry/algos/offset_axis/normal_utils.h>

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <vector>

namespace cgogn
{

namespace geometry
{

namespace offset_axis
{

inline bool evaluate_udf_values(const NeuralFieldForward& field, const std::vector<Vec3>& points, size_t batch_size,
								std::vector<Scalar>& out_values)
{
	out_values.clear();
	out_values.resize(points.size(), Scalar(0));
	if (points.empty())
		return true;
	if (!field.is_loaded())
		return false;

	const size_t batch = std::max<size_t>(1, batch_size);
	for (size_t offset = 0; offset < points.size(); offset += batch)
	{
		const size_t count = std::min(batch, points.size() - offset);
		torch::Tensor points_gpu = detail::make_points_tensor(points, offset, count, field.device()).to(field.device());
		torch::Tensor values = field.forward_values_gpu(points_gpu);
		if (!values.defined())
			return false;
		torch::Tensor values_cpu = values.to(torch::kCPU).contiguous();
		auto values_acc = values_cpu.accessor<float, 1>();
		for (size_t i = 0; i < count; ++i)
			out_values[offset + i] = static_cast<Scalar>(values_acc[static_cast<long>(i)]);
	}
	return true;
}

inline bool evaluate_udf_normals(const NeuralFieldForward& field, const std::vector<Vec3>& points, size_t batch_size,
								NeuralProjectionWorkspace& workspace, std::vector<Vec3>& normals)
{
	normals.clear();
	if (points.empty())
		return true;
	if (!field.is_loaded())
		return false;

	const size_t batch = std::max<size_t>(1, batch_size);
	detail::ensure_neural_projection_workspace(workspace, std::min(points.size(), batch), field.device());
	normals.assign(points.size(), Vec3(0, 0, 1));
	for (size_t offset = 0; offset < points.size(); offset += batch)
	{
		const size_t count = std::min(batch, points.size() - offset);
		const int64_t tensor_count = static_cast<int64_t>(count);
		torch::Tensor points_cpu = workspace.points_cpu_.narrow(0, 0, tensor_count);
		detail::fill_points_tensor(points, offset, count, points_cpu);

		torch::Tensor points_gpu = workspace.points_gpu_.narrow(0, 0, tensor_count);
		points_gpu.copy_(points_cpu);
		auto [values, gradients] = field.forward_values_grad_gpu(points_gpu);
		if (!values.defined() || !gradients.defined())
			return false;
		if (gradients.dim() != 2 || gradients.size(0) != tensor_count || gradients.size(1) != 3)
			return false;

		torch::Tensor gradients_cpu = workspace.grad_cpu_.narrow(0, 0, tensor_count);
		gradients_cpu.copy_(gradients);
		auto gradients_acc = gradients_cpu.accessor<float, 2>();
		for (size_t i = 0; i < count; ++i)
		{
			Vec3 normal(static_cast<Scalar>(gradients_acc[static_cast<long>(i)][0]),
						static_cast<Scalar>(gradients_acc[static_cast<long>(i)][1]),
						static_cast<Scalar>(gradients_acc[static_cast<long>(i)][2]));
			normals[offset + i] = unit_normal_or_z(normal);
		}
	}
	return true;
}

// udf value (mf - |sdf| for MF models) and sdf of every point; false if the model has no sdf output
inline bool evaluate_udf_values_sdf(const NeuralFieldForward& field, const std::vector<Vec3>& points,
									size_t batch_size, std::vector<Scalar>& out_values, std::vector<Scalar>& out_sdf)
{
	out_values.assign(points.size(), Scalar(0));
	out_sdf.assign(points.size(), Scalar(0));
	if (points.empty())
		return true;
	if (!field.is_loaded())
		return false;

	const size_t batch = std::max<size_t>(1, batch_size);
	for (size_t offset = 0; offset < points.size(); offset += batch)
	{
		const size_t count = std::min(batch, points.size() - offset);
		auto [values, sdf] = field.forward_values_sdf_gpu(detail::make_points_tensor(points, offset, count, field.device()));
		if (!values.defined() || !sdf.defined() || values.numel() != static_cast<int64_t>(count) ||
			sdf.numel() != static_cast<int64_t>(count))
			return false;
		torch::Tensor both = torch::stack({values.to(torch::kFloat32), sdf.to(torch::kFloat32)}, 1).to(torch::kCPU);
		auto acc = both.accessor<float, 2>();
		for (size_t i = 0; i < count; ++i)
		{
			out_values[offset + i] = static_cast<Scalar>(acc[static_cast<long>(i)][0]);
			out_sdf[offset + i] = static_cast<Scalar>(acc[static_cast<long>(i)][1]);
		}
	}
	return true;
}

} // namespace offset_axis

} // namespace geometry

} // namespace cgogn

#endif // CGOGN_GEOMETRY_ALGOS_OFFSET_AXIS_NEURAL_NEURAL_FIELD_QUERY_H_
