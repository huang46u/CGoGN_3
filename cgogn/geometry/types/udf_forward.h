#ifndef CGOGN_GEOMETRY_TYPES_UDF_FORWARD_H_
#define CGOGN_GEOMETRY_TYPES_UDF_FORWARD_H_

#include <cgogn/geometry/types/vector_traits.h>

#include <torch/script.h>
#include <torch/torch.h>

#include <iostream>
#include <limits>
#include <utility>
#include <vector>

namespace cgogn
{

namespace geometry
{

struct BatchUDFResult
{
	std::vector<Scalar> values;	 // size N
	std::vector<Vec3> gradients; // size N
	bool ok = false;
};

class UDFForward
{
public:
	UDFForward(torch::jit::Module* module, bool loaded, const torch::Device& device)
		: module_(module), loaded_(loaded), device_(device)
	{
	}

	bool is_loaded() const
	{
		return loaded_ && module_;
	}

	const torch::Device& device() const
	{
		return device_;
	}

	Scalar forward_point(const Vec3& query_point) const
	{
		return forward_point_impl(query_point, device_);
	}

	Scalar forward_point_cpu(const Vec3& query_point) const
	{
		return forward_point_impl(query_point, torch::kCPU);
	}

	Scalar forward_point_gpu(const Vec3& query_point) const
	{
		return forward_point_impl(query_point, device_);
	}

	std::pair<Scalar, Vec3> forward_point_with_grad(const Vec3& query_point) const
	{
		return forward_point_with_grad_impl(query_point, device_);
	}

	std::pair<Scalar, Vec3> forward_point_with_grad_cpu(const Vec3& query_point) const
	{
		return forward_point_with_grad_impl(query_point, torch::kCPU);
	}

	std::pair<Scalar, Vec3> forward_point_with_grad_gpu(const Vec3& query_point) const
	{
		return forward_point_with_grad_impl(query_point, device_);
	}

	BatchUDFResult forward_batch(const std::vector<Vec3>& query_points) const
	{
		return forward_batch_impl(query_points, device_, false);
	}

	BatchUDFResult forward_batch_cpu(const std::vector<Vec3>& query_points) const
	{
		return forward_batch_impl(query_points, torch::kCPU, false);
	}

	BatchUDFResult forward_batch_gpu(const std::vector<Vec3>& query_points) const
	{
		return forward_batch_impl(query_points, device_, false);
	}

	BatchUDFResult forward_batch_with_grad(const std::vector<Vec3>& query_points) const
	{
		return forward_batch_impl(query_points, device_, true);
	}

	BatchUDFResult forward_batch_with_grad_cpu(const std::vector<Vec3>& query_points) const
	{
		return forward_batch_impl(query_points, torch::kCPU, true);
	}

	BatchUDFResult forward_batch_with_grad_gpu(const std::vector<Vec3>& query_points) const
	{
		return forward_batch_impl(query_points, device_, true);
	}

	torch::Tensor forward_values_gpu(const torch::Tensor& points_gpu) const
	{
		if (!is_loaded())
			return torch::Tensor();
		try
		{
			ModuleDeviceScope module_scope(module_, device_, device_);
			torch::NoGradGuard no_grad;
			torch::Tensor x = points_gpu;
			if (x.device() != device_)
				x = x.to(device_);

			torch::Tensor output = module_->forward({x}).toTensor();
			if (output.dim() == 2 && output.size(1) == 1)
				output = output.squeeze(1);
			return output.detach();
		}
		catch (const c10::Error& e)
		{
			std::cerr << "UDF forward (values) failed: " << e.what() << std::endl;
			return torch::Tensor();
		}
	}

	std::pair<torch::Tensor, torch::Tensor> forward_values_grad_gpu(const torch::Tensor& points_gpu) const
	{
		if (!is_loaded())
			return {torch::Tensor(), torch::Tensor()};
		try
		{
			ModuleDeviceScope module_scope(module_, device_, device_);

			torch::Tensor x = points_gpu;
			if (x.device() != device_)
				x = x.to(device_);
			x = x.detach().requires_grad_(true);

			torch::Tensor y = module_->forward({x}).toTensor();
			if (y.dim() == 2 && y.size(1) == 1)
				y = y.squeeze(1);

			torch::Tensor grad_outputs = torch::ones_like(y);
			std::vector<torch::Tensor> grads = torch::autograd::grad(
				/*outputs=*/{y},
				/*inputs=*/{x},
				/*grad_outputs=*/{grad_outputs},
				/*retain_graph=*/false,
				/*create_graph=*/false,
				/*allow_unused=*/false);

			torch::Tensor dx = grads[0];
			return {y.detach(), dx.detach()};
		}
		catch (const c10::Error& e)
		{
			std::cerr << "UDF forward (values+grad) failed: " << e.what() << std::endl;
			return {torch::Tensor(), torch::Tensor()};
		}
	}

private:
	class ModuleDeviceScope
	{
	public:
		ModuleDeviceScope(torch::jit::Module* module, const torch::Device& current_device,
						  const torch::Device& target_device)
			: module_(module), current_device_(current_device), moved_(false)
		{
			if (module_ && current_device_ != target_device)
			{
				module_->to(target_device);
				moved_ = true;
			}
		}

		~ModuleDeviceScope()
		{
			if (module_ && moved_)
				module_->to(current_device_);
		}

	private:
		torch::jit::Module* module_;
		torch::Device current_device_;
		bool moved_;
	};

	Scalar forward_point_impl(const Vec3& query_point, const torch::Device& device) const
	{
		BatchUDFResult r = forward_batch_impl({query_point}, device, false);
		if (!r.ok || r.values.empty())
			return std::numeric_limits<Scalar>::max();
		return r.values[0];
	}

	std::pair<Scalar, Vec3> forward_point_with_grad_impl(const Vec3& query_point, const torch::Device& device) const
	{
		BatchUDFResult r = forward_batch_impl({query_point}, device, true);
		if (!r.ok || r.values.empty() || r.gradients.empty())
			return {std::numeric_limits<Scalar>::max(), Vec3(0, 0, 0)};
		return {r.values[0], r.gradients[0]};
	}

	BatchUDFResult forward_batch_impl(const std::vector<Vec3>& query_points, const torch::Device& device,
									  bool with_grad) const
	{
		BatchUDFResult r;
		r.ok = false;
		const size_t N = query_points.size();
		if (!is_loaded() || N == 0)
			return r;
		try
		{
			ModuleDeviceScope module_scope(module_, device_, device);
			if (!with_grad)
			{
				torch::NoGradGuard no_grad;
				return forward_batch_impl_no_grad(query_points, device);
			}
			return forward_batch_impl_with_grad(query_points, device);
		}
		catch (const c10::Error& e)
		{
			std::cerr << "UDF batch forward failed: " << e.what() << std::endl;
			return r;
		}
	}

	BatchUDFResult forward_batch_impl_no_grad(const std::vector<Vec3>& query_points, const torch::Device& device) const
	{
		BatchUDFResult r;
		r.ok = false;
		const size_t N = query_points.size();

		auto cpu_opts = torch::TensorOptions().dtype(torch::kFloat32).device(torch::kCPU);
		if (device.is_cuda())
			cpu_opts = cpu_opts.pinned_memory(true);

		torch::Tensor points_cpu = torch::empty({static_cast<long>(N), 3}, cpu_opts);
		auto acc = points_cpu.accessor<float, 2>();
		for (size_t i = 0; i < N; ++i)
		{
			acc[i][0] = query_points[i].x();
			acc[i][1] = query_points[i].y();
			acc[i][2] = query_points[i].z();
		}

		torch::Tensor points = points_cpu.to(device);
		torch::Tensor output = module_->forward({points}).toTensor();
		if (output.dim() == 2 && output.size(1) == 1)
			output = output.squeeze(1);
		if (output.dim() != 1 || output.size(0) != static_cast<long>(N))
			return r;

		torch::Tensor output_cpu = output.detach().to(torch::kCPU).contiguous();
		auto out_acc = output_cpu.accessor<float, 1>();

		r.values.resize(N);
		r.gradients.clear();
		for (size_t i = 0; i < N; ++i)
			r.values[i] = static_cast<Scalar>(out_acc[(long)i]);

		r.ok = true;
		return r;
	}

	BatchUDFResult forward_batch_impl_with_grad(const std::vector<Vec3>& query_points, const torch::Device& device) const
	{
		BatchUDFResult r;
		r.ok = false;
		const size_t N = query_points.size();

		auto cpu_opts = torch::TensorOptions().dtype(torch::kFloat32).device(torch::kCPU);
		if (device.is_cuda())
			cpu_opts = cpu_opts.pinned_memory(true);

		torch::Tensor points_cpu = torch::empty({static_cast<long>(N), 3}, cpu_opts);
		auto acc = points_cpu.accessor<float, 2>();
		for (size_t i = 0; i < N; ++i)
		{
			acc[i][0] = query_points[i].x();
			acc[i][1] = query_points[i].y();
			acc[i][2] = query_points[i].z();
		}

		torch::Tensor points = points_cpu.to(device);
		points.set_requires_grad(true);

		torch::Tensor output = module_->forward({points}).toTensor();
		if (output.dim() == 2 && output.size(1) == 1)
			output = output.squeeze(1);
		if (output.dim() != 1 || output.size(0) != static_cast<long>(N))
			return r;

		torch::Tensor grad_outputs = torch::ones_like(output);
		std::vector<torch::Tensor> grads = torch::autograd::grad(
			/*outputs=*/{output},
			/*inputs=*/{points},
			/*grad_outputs=*/{grad_outputs},
			/*retain_graph=*/false,
			/*create_graph=*/false,
			/*allow_unused=*/false);

		torch::Tensor grad = grads[0];
		if (!grad.defined() || grad.dim() != 2 || grad.size(0) != static_cast<long>(N) || grad.size(1) != 3)
			return r;

		torch::Tensor output_cpu = output.detach().to(torch::kCPU).contiguous();
		torch::Tensor grad_cpu = grad.detach().to(torch::kCPU).contiguous();

		auto out_acc = output_cpu.accessor<float, 1>();
		auto grad_acc = grad_cpu.accessor<float, 2>();

		r.values.resize(N);
		r.gradients.resize(N);
		for (size_t i = 0; i < N; ++i)
		{
			r.values[i] = static_cast<Scalar>(out_acc[(long)i]);
			r.gradients[i] = Vec3(static_cast<Scalar>(grad_acc[(long)i][0]), static_cast<Scalar>(grad_acc[(long)i][1]),
								  static_cast<Scalar>(grad_acc[(long)i][2]));
		}

		r.ok = true;
		return r;
	}

	torch::jit::Module* module_;
	bool loaded_;
	torch::Device device_;
};

} // namespace geometry

} // namespace cgogn

#endif // CGOGN_GEOMETRY_TYPES_UDF_FORWARD_H_
