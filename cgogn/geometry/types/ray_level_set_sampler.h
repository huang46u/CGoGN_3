#ifndef CGOGN_GEOMETRY_TYPES_RAY_LEVEL_SET_SAMPLER_H_
#define CGOGN_GEOMETRY_TYPES_RAY_LEVEL_SET_SAMPLER_H_

#include <cgogn/core/utils/numerics.h>
#include <cgogn/geometry/types/spatial_grid.h>
#include <cgogn/geometry/types/neural_field_forward.h>
#include <cgogn/geometry/types/vector_traits.h>

#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>
#include <random>
#include <utility>
#include <vector>

namespace cgogn
{

namespace geometry
{

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

	std::vector<Vec3> sample_alpha_level_set_rays(NeuralFieldForward& udf, size_t target_num_points,
												  SpatialGrid* spatial_grid, Scalar grid_cell_size,
												  const Vec3& bbox_min, const Vec3& bbox_max,
												  bool* used_sdf_filter = nullptr)
	{
		if (!udf.is_loaded())
		{
			std::cerr << "Neural UDF model not loaded." << std::endl;
			return {};
		}

		set_device(udf.device());

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
			auto d_sdf = udf.forward_values_sdf_gpu(X);
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
						pt = pt.cwiseMax(0.0).cwiseMin(1.0);

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

					auto d_sdf_near = udf.forward_values_sdf_gpu(p_near);
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

				auto d_sdf_far = udf.forward_values_sdf_gpu(p_far);
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

private:
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
