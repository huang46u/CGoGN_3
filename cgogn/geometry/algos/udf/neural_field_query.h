#ifndef CGOGN_GEOMETRY_ALGOS_UDF_NEURAL_FIELD_QUERY_H_
#define CGOGN_GEOMETRY_ALGOS_UDF_NEURAL_FIELD_QUERY_H_

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
		auto cpu_opts = torch::TensorOptions().dtype(torch::kFloat32).device(torch::kCPU);
		if (field.device().is_cuda())
			cpu_opts = cpu_opts.pinned_memory(true);

		torch::Tensor points_cpu = torch::empty({static_cast<int64_t>(count), 3}, cpu_opts);
		auto points_acc = points_cpu.accessor<float, 2>();
		for (size_t i = 0; i < count; ++i)
		{
			const Vec3& point = points[offset + i];
			points_acc[static_cast<long>(i)][0] = static_cast<float>(point.x());
			points_acc[static_cast<long>(i)][1] = static_cast<float>(point.y());
			points_acc[static_cast<long>(i)][2] = static_cast<float>(point.z());
		}

		torch::Tensor points_gpu = points_cpu.to(field.device());
		torch::Tensor values = field.forward_values_gpu(points_gpu);
		if (!values.defined())
			return false;
		if (values.dim() == 2 && values.size(1) == 1)
			values = values.squeeze(1);
		torch::Tensor values_cpu = values.to(torch::kCPU).contiguous();
		auto values_acc = values_cpu.accessor<float, 1>();
		for (size_t i = 0; i < count; ++i)
			out_values[offset + i] = static_cast<Scalar>(values_acc[static_cast<long>(i)]);
	}
	return true;
}

inline bool evaluate_udf_value_and_gradient(const NeuralFieldForward& field, const Vec3& query_point, Scalar& value,
												Vec3& gradient)
{
	if (!field.is_loaded())
		return false;
	const auto result = field.forward_point_with_grad(query_point);
	value = result.first;
	gradient = result.second;
	return std::isfinite(value) && gradient.allFinite();
}

} // namespace geometry

} // namespace cgogn

#endif // CGOGN_GEOMETRY_ALGOS_UDF_NEURAL_FIELD_QUERY_H_
