#ifndef CGOGN_GEOMETRY_TORCH_VMAS_TROCH_H_
#define CGOGN_GEOMETRY_TORCH_VMAS_TROCH_H_

#include <torch/torch.h>

#include <algorithm>
#include <limits>
#include <tuple>
#include <utility>

namespace cgogn
{

namespace geometry
{

namespace vmas_troch
{

inline torch::Tensor eval_sqem_batch(const torch::Tensor& A, // [N,4,4]
									 const torch::Tensor& b, // [N,4]
									 const torch::Tensor& c, // [N]
									 const torch::Tensor& v  // [B,4]
)
{
	torch::Tensor v_exp = v.unsqueeze(0).expand({A.size(0), v.size(0), 4}); // [N,B,4]
	torch::Tensor A_exp = A.unsqueeze(1);									// [N,1,4,4]

	torch::Tensor Av = torch::matmul(A_exp, v_exp.unsqueeze(-1)).squeeze(-1); // [N,B,4]

	torch::Tensor vTAv = (v_exp * Av).sum(-1);									 // [N,B]
	torch::Tensor b_exp = b.unsqueeze(1);										 // [N,1,4]
	torch::Tensor dist = 0.5f * vTAv - (b_exp * v_exp).sum(-1) + c.unsqueeze(1); // [N,B]
	return dist;
}

inline torch::Tensor eval_line_quadric_batch(const torch::Tensor& Q, // [N,4,4]
											 const torch::Tensor& v) // [B,4], v = [cx,cy,cz,1]
{
	torch::Tensor v_exp = v.unsqueeze(0).expand({Q.size(0), v.size(0), 4});
	torch::Tensor Q_exp = Q.unsqueeze(1);
	torch::Tensor Qv = torch::matmul(Q_exp, v_exp.unsqueeze(-1)).squeeze(-1);
	torch::Tensor vTQv = (v_exp * Qv).sum(-1);
	return vTQv;
}

inline torch::Tensor eval_sqem_single(const torch::Tensor& A, // [N,4,4]
									  const torch::Tensor& b, // [N,4]
									  const torch::Tensor& c, // [N]
									  const torch::Tensor& v  // [N,4]
)
{
	torch::Tensor Av = torch::matmul(A, v.unsqueeze(-1)).squeeze(-1); // [N,4]
	torch::Tensor vTAv = (v * Av).sum(-1);								// [N]
	torch::Tensor dist = 0.5f * vTAv - (b * v).sum(-1) + c;			// [N]
	return dist;
}

inline torch::Tensor eval_line_quadric_single(const torch::Tensor& Q, // [N,4,4]
											  const torch::Tensor& v) // [N,4]
{
	torch::Tensor Qv = torch::matmul(Q, v.unsqueeze(-1)).squeeze(-1); // [N,4]
	return (v * Qv).sum(-1);										 // [N]
}

inline std::pair<torch::Tensor, torch::Tensor> compute_clusters_gpu(
	const torch::Tensor& samples_pos,   // [N,3]
	const torch::Tensor& samples_area,  // [N]
	const torch::Tensor& centers,	   // [M,3]
	const torch::Tensor& radius,		   // [M]
	const torch::Tensor& samples_sqem_A, // [N,4,4]
	const torch::Tensor& samples_sqem_b, // [N,4]
	const torch::Tensor& samples_sqem_c, // [N]
	const torch::Tensor& samples_line_Q, // [N,4,4]
	int64_t distance_mode,
	float sqem_clustering_lambda,
	int64_t chunk = 128)
{
	const auto device = samples_pos.device();
	const auto fopts = torch::TensorOptions().dtype(torch::kFloat32).device(device);
	const auto iopts = torch::TensorOptions().dtype(torch::kInt64).device(device);

	const int64_t N = samples_pos.size(0);
	const int64_t M = centers.size(0);

	torch::Tensor min_dist = torch::full({N}, std::numeric_limits<float>::infinity(), fopts);
	torch::Tensor argmin = torch::full({N}, -1, iopts);

	for (int64_t k = 0; k < M; k += chunk)
	{
		const int64_t B = std::min(chunk, M - k);
		auto centers_c = centers.narrow(0, k, B); // [B,3]
		auto radius_c = radius.narrow(0, k, B);	 // [B]

		torch::Tensor dist;

		if (distance_mode == 2) // PURE_EUCLIDEAN_DISTANCE
		{
			torch::Tensor diff = samples_pos.unsqueeze(1) - centers_c.unsqueeze(0); // [N,B,3]
			torch::Tensor d = torch::sqrt(torch::sum(diff * diff, 2));				// [N,B]

			dist = d - radius_c.unsqueeze(0);
			dist = dist * dist;
			dist = dist * samples_area.unsqueeze(1);
		}
		else //SPHERE_EUCLIDEAN_DISTANCE
		{
			torch::Tensor v_sphere = torch::cat({centers_c, radius_c.unsqueeze(1)}, 1); // [B,4]
			torch::Tensor dist_sqem =
				eval_sqem_batch(samples_sqem_A, samples_sqem_b, samples_sqem_c, v_sphere); // [N,B]

			torch::Tensor dist_other;
			if (distance_mode == 0) // SPHERE_EUCLIDEAN_DISTANCE
			{
				torch::Tensor diff = samples_pos.unsqueeze(1) - centers_c.unsqueeze(0); // [N,B,3]
				torch::Tensor d = torch::sqrt(torch::sum(diff * diff, 2));				// [N,B]

				dist_other = d - radius_c.unsqueeze(0);
				dist_other = dist_other * dist_other * samples_area.unsqueeze(1);
			}
			else // LINE_QUADRIC_DISTANCE
			{
				torch::Tensor v_sphere_line = torch::cat({centers_c, torch::ones({B, 1}, fopts)}, 1); // [B,4]
				dist_other = eval_line_quadric_batch(samples_line_Q, v_sphere_line);				   // [N,B]
			}

			dist = dist_sqem + sqem_clustering_lambda * dist_other;
		}

		auto min_pair = dist.min(1);
		torch::Tensor dist_min = std::get<0>(min_pair); // [N]
		torch::Tensor dist_idx = std::get<1>(min_pair); // [N] in [0,B)

		torch::Tensor global_idx = dist_idx + k;
		torch::Tensor mask = dist_min < min_dist;

		min_dist = torch::where(mask, dist_min, min_dist);
		argmin = torch::where(mask, global_idx, argmin);
	}

	return {argmin, min_dist};
}

inline std::pair<torch::Tensor, torch::Tensor> update_spheres_euclidean_gpu(
	const torch::Tensor& samples_pos,  // [N,3]
	const torch::Tensor& samples_area, // [N]
	const torch::Tensor& cluster_idx,  // [N]
	const torch::Tensor& centers,	  // [B,3]
	float alpha,
	int64_t max_iters = 10,
	float eps = 1e-8f)
{
	const auto device = samples_pos.device();
	const auto fopts = torch::TensorOptions().dtype(torch::kFloat32).device(device);

	const int64_t N = samples_pos.size(0);
	const int64_t B = centers.size(0);

	torch::Tensor centers_cur = centers.clone();

	torch::Tensor valid = cluster_idx >= 0;
	torch::Tensor valid_f = valid.to(fopts);
	torch::Tensor cluster_idx_safe = torch::where(valid, cluster_idx, torch::zeros_like(cluster_idx));

	torch::Tensor counts = torch::zeros({B}, fopts);
	counts.index_add_(0, cluster_idx_safe, valid_f);
	torch::Tensor has_samples = counts > 0;

	torch::Tensor eye = torch::eye(3, fopts).unsqueeze(0);

	for (int64_t iter = 0; iter < max_iters; ++iter)
	{
		torch::Tensor c_i = centers_cur.index_select(0, cluster_idx_safe);
		torch::Tensor d = samples_pos - c_i;
		torch::Tensor l = torch::sqrt(torch::sum(d * d, 1)).clamp_min(eps);

		torch::Tensor w = torch::sqrt(samples_area).clamp_min(0) * valid_f;
		torch::Tensor J = (-d / l.unsqueeze(1)) * w.unsqueeze(1); // [N,3]
		torch::Tensor b = -(l - alpha) * w;					   // [N]

		torch::Tensor outer = J.unsqueeze(2) * J.unsqueeze(1); // [N,3,3]
		torch::Tensor H = torch::zeros({B, 3, 3}, fopts);
		H.view({B, 9}).index_add_(0, cluster_idx_safe, outer.view({N, 9}));
		H = H + eps * eye;

		torch::Tensor g = torch::zeros({B, 3}, fopts);
		g.index_add_(0, cluster_idx_safe, J * b.unsqueeze(1));

		torch::Tensor delta = torch::linalg_solve(H, g.unsqueeze(-1)).squeeze(-1);
		delta = torch::where(has_samples.unsqueeze(1), delta, torch::zeros_like(delta));

		centers_cur = centers_cur + delta;
	}

	torch::Tensor radii = torch::full({B}, alpha, fopts);
	return {centers_cur, radii};
}

inline std::pair<torch::Tensor, torch::Tensor> update_spheres_sqem_gpu(
	const torch::Tensor& samples_pos,	 // [N,3]
	const torch::Tensor& samples_normal, // [N,3]
	const torch::Tensor& samples_area,	 // [N]
	const torch::Tensor& knn_pos,		 // [N,K,3]
	const torch::Tensor& knn_normal,	 // [N,K,3]
	const torch::Tensor& knn_area,		 // [N,K]
	const torch::Tensor& cluster_idx,	 // [N]
	const torch::Tensor& centers,		 // [B,3]
	const torch::Tensor& radii,			 // [B]
	int64_t knn_k,
	float sqem_update_lambda,
	int64_t max_iters = 10,
	float eps = 1e-8f)
{
	const auto device = samples_pos.device();
	const auto fopts = torch::TensorOptions().dtype(torch::kFloat32).device(device);

	const int64_t N = samples_pos.size(0);
	const int64_t B = centers.size(0);
	const int64_t K = knn_pos.defined() ? knn_pos.size(1) : 0;

	torch::Tensor centers_cur = centers.clone();
	torch::Tensor radii_cur = radii.clone();

	torch::Tensor valid = cluster_idx >= 0;
	torch::Tensor valid_f = valid.to(fopts);
	torch::Tensor cluster_idx_safe = torch::where(valid, cluster_idx, torch::zeros_like(cluster_idx));

	torch::Tensor counts = torch::zeros({B}, fopts);
	counts.index_add_(0, cluster_idx_safe, valid_f);
	torch::Tensor has_samples = counts > 0;

	torch::Tensor eye = torch::eye(4, fopts).unsqueeze(0);
	const float knn_den = static_cast<float>(knn_k) + 1.0f;

	for (int64_t iter = 0; iter < max_iters; ++iter)
	{
		torch::Tensor c_i = centers_cur.index_select(0, cluster_idx_safe); // [N,3]
		torch::Tensor r_i = radii_cur.index_select(0, cluster_idx_safe);	// [N]

		torch::Tensor a_self = torch::sqrt(samples_area / knn_den).clamp_min(0) * valid_f;
		torch::Tensor ones = torch::ones({N, 1}, fopts);
		torch::Tensor n4_self = torch::cat({samples_normal, ones}, 1); // [N,4]
		torch::Tensor lhs = -n4_self * a_self.unsqueeze(1);			 // [N,4]

		torch::Tensor dot_self = ((samples_pos - c_i) * samples_normal).sum(1); // [N]
		torch::Tensor rhs = -(dot_self - r_i) * a_self;						 // [N]

		if (K > 0 && knn_pos.defined() && knn_normal.defined() && knn_area.defined())
		{
			torch::Tensor a_knn = torch::sqrt(knn_area / knn_den).clamp_min(0); // [N,K]
			torch::Tensor ones_knn = torch::ones({N, K, 1}, fopts);
			torch::Tensor n4_knn = torch::cat({knn_normal, ones_knn}, 2); // [N,K,4]

			lhs += -(n4_knn * a_knn.unsqueeze(-1)).sum(1);

			torch::Tensor dot_knn = ((knn_pos - c_i.unsqueeze(1)) * knn_normal).sum(-1); // [N,K]
			rhs += -((dot_knn - r_i.unsqueeze(1)) * a_knn).sum(1);
		}

		torch::Tensor d = samples_pos - c_i;
		torch::Tensor l = torch::sqrt(torch::sum(d * d, 1)).clamp_min(eps);
		torch::Tensor a_dist = torch::sqrt(samples_area).clamp_min(0) * valid_f;
		torch::Tensor scale = (a_dist * sqem_update_lambda).unsqueeze(1);

		torch::Tensor dist_row = torch::cat({-(d / l.unsqueeze(1)), -torch::ones({N, 1}, fopts)}, 1) * scale; // [N,4]
		torch::Tensor b_dist = -(l - r_i) * a_dist * sqem_update_lambda;										 // [N]

		torch::Tensor H = torch::zeros({B, 4, 4}, fopts);
		torch::Tensor outer_lhs = lhs.unsqueeze(2) * lhs.unsqueeze(1);
		torch::Tensor outer_dist = dist_row.unsqueeze(2) * dist_row.unsqueeze(1);
		H.view({B, 16}).index_add_(0, cluster_idx_safe, (outer_lhs + outer_dist).view({N, 16}));
		H = H + eps * eye;

		torch::Tensor g = torch::zeros({B, 4}, fopts);
		torch::Tensor g_rows = lhs * rhs.unsqueeze(1) + dist_row * b_dist.unsqueeze(1);
		g.index_add_(0, cluster_idx_safe, g_rows);

		torch::Tensor delta = torch::linalg_solve(H, g.unsqueeze(-1)).squeeze(-1); // [B,4]
		delta = torch::where(has_samples.unsqueeze(1), delta, torch::zeros_like(delta));

		centers_cur = centers_cur + delta.narrow(1, 0, 3);
		radii_cur = radii_cur + delta.select(1, 3);
	}

	return {centers_cur, radii_cur};
}

inline std::pair<torch::Tensor, torch::Tensor> update_spheres_line_quadric_gpu(
	const torch::Tensor& samples_area,   // [N]
	const torch::Tensor& samples_sqem_A, // [N,4,4]
	const torch::Tensor& samples_sqem_b, // [N,4]
	const torch::Tensor& samples_line_Q, // [N,4,4]
	const torch::Tensor& cluster_idx,	// [N]
	const torch::Tensor& centers,		// [B,3]
	float alpha,
	float sqem_update_lambda,
	float eps = 1e-8f)
{
	const auto device = samples_area.device();
	const auto fopts = torch::TensorOptions().dtype(torch::kFloat32).device(device);

	const int64_t N = samples_area.size(0);
	const int64_t B = centers.size(0);

	torch::Tensor valid = cluster_idx >= 0;
	torch::Tensor valid_f = valid.to(fopts);
	torch::Tensor cluster_idx_safe = torch::where(valid, cluster_idx, torch::zeros_like(cluster_idx));

	torch::Tensor counts = torch::zeros({B}, fopts);
	counts.index_add_(0, cluster_idx_safe, valid_f);
	torch::Tensor has_samples = counts > 0;

	torch::Tensor w = samples_area * valid_f;
	torch::Tensor A_w = samples_sqem_A * w.view({N, 1, 1});
	torch::Tensor b_w = samples_sqem_b * w.view({N, 1});
	torch::Tensor Q_w = samples_line_Q * w.view({N, 1, 1});

	torch::Tensor A_sum = torch::zeros({B, 4, 4}, fopts);
	torch::Tensor b_sum = torch::zeros({B, 4}, fopts);
	torch::Tensor Q_sum = torch::zeros({B, 4, 4}, fopts);

	A_sum.view({B, 16}).index_add_(0, cluster_idx_safe, A_w.view({N, 16}));
	b_sum.index_add_(0, cluster_idx_safe, b_w);
	Q_sum.view({B, 16}).index_add_(0, cluster_idx_safe, Q_w.view({N, 16}));

	torch::Tensor As = A_sum.narrow(1, 0, 3).narrow(2, 0, 3); // [B,3,3]
	torch::Tensor Asr = A_sum.narrow(1, 0, 3).select(2, 3);	// [B,3]
	torch::Tensor bs = b_sum.narrow(1, 0, 3);				// [B,3]

	torch::Tensor Al = Q_sum.narrow(1, 0, 3).narrow(2, 0, 3); // [B,3,3]
	torch::Tensor bl = -Q_sum.narrow(1, 0, 3).select(2, 3);	  // [B,3]

	torch::Tensor A = As + sqem_update_lambda * Al;
	torch::Tensor b = (bs + sqem_update_lambda * bl) - Asr * alpha;

	torch::Tensor eye = torch::eye(3, fopts).unsqueeze(0);
	A = A + eps * eye;

	torch::Tensor x = torch::linalg_solve(A, b.unsqueeze(-1)).squeeze(-1); // [B,3]
	x = torch::where(has_samples.unsqueeze(1), x, centers);

	torch::Tensor r = torch::full({B}, alpha, fopts);
	return {x, r};
}

inline std::tuple<torch::Tensor, torch::Tensor, torch::Tensor, torch::Tensor, torch::Tensor, torch::Tensor, torch::Tensor,
				  torch::Tensor>
compute_spheres_error_gpu(const torch::Tensor& samples_pos,	 // [N,3]
						  const torch::Tensor& samples_area,	 // [N]
						  const torch::Tensor& samples_sqem_A, // [N,4,4]
						  const torch::Tensor& samples_sqem_b, // [N,4]
						  const torch::Tensor& samples_sqem_c, // [N]
						  const torch::Tensor& samples_line_Q, // [N,4,4]
						  const torch::Tensor& centers,		 // [B,3]
						  const torch::Tensor& radius,		 // [B]
						  const torch::Tensor& cluster_idx,	 // [N]
						  int64_t distance_mode,
						  float sqem_clustering_lambda,
						  float eps = 1e-8f)
{
	const auto device = samples_pos.device();
	const auto fopts = torch::TensorOptions().dtype(torch::kFloat32).device(device);

	const int64_t N = samples_pos.size(0);
	const int64_t B = centers.size(0);

	if (N == 0 || B == 0)
	{
		torch::Tensor empty = torch::zeros({0}, fopts);
		torch::Tensor zero = torch::zeros({}, fopts);
		torch::Tensor zero_i = torch::zeros({}, torch::TensorOptions().dtype(torch::kInt64).device(device));
		return {empty, empty, empty, zero, zero, zero, zero, zero_i};
	}

	torch::Tensor valid = cluster_idx >= 0;
	torch::Tensor valid_f = valid.to(fopts);
	torch::Tensor cluster_idx_safe = torch::where(valid, cluster_idx, torch::zeros_like(cluster_idx));

	torch::Tensor c_i = centers.index_select(0, cluster_idx_safe); // [N,3]
	torch::Tensor r_i = radius.index_select(0, cluster_idx_safe);	 // [N]

	torch::Tensor dist;
	if (distance_mode == 2) // PURE_EUCLIDEAN_DISTANCE
	{
		torch::Tensor d = torch::sqrt(torch::sum((samples_pos - c_i) * (samples_pos - c_i), 1));
		torch::Tensor diff = d - r_i;
		dist = diff * diff * samples_area;
	}
	else
	{
		torch::Tensor v_sphere = torch::cat({c_i, r_i.unsqueeze(1)}, 1); // [N,4]
		torch::Tensor dist_sqem = eval_sqem_single(samples_sqem_A, samples_sqem_b, samples_sqem_c, v_sphere);

		torch::Tensor dist_other;
		if (distance_mode == 0) // SPHERE_EUCLIDEAN_DISTANCE
		{
			torch::Tensor d = torch::sqrt(torch::sum((samples_pos - c_i) * (samples_pos - c_i), 1));
			torch::Tensor diff = d - r_i;
			dist_other = diff * diff * samples_area;
		}
		else // LINE_QUADRIC_DISTANCE
		{
			torch::Tensor v_sphere_line = torch::cat({c_i, torch::ones({N, 1}, fopts)}, 1);
			dist_other = eval_line_quadric_single(samples_line_Q, v_sphere_line);
		}

		dist = dist_sqem + sqem_clustering_lambda * dist_other;
	}

	dist = dist * valid_f;
	torch::Tensor area = samples_area * valid_f;

	torch::Tensor cluster_error = torch::zeros({B}, fopts);
	cluster_error.index_add_(0, cluster_idx_safe, dist);

	torch::Tensor cluster_area = torch::zeros({B}, fopts);
	cluster_area.index_add_(0, cluster_idx_safe, area);

	torch::Tensor sphere_error = torch::where(cluster_area > 0, cluster_error / cluster_area, torch::zeros_like(cluster_area));
	torch::Tensor sphere_error_nn = cluster_error;

	torch::Tensor total_error = sphere_error.sum();
	torch::Tensor total_error_not_norm = cluster_error.sum();
	torch::Tensor min_error = sphere_error.min();
	auto max_pair = sphere_error.max(0);
	torch::Tensor max_error = std::get<0>(max_pair);
	torch::Tensor max_idx = std::get<1>(max_pair);

	return {sphere_error, sphere_error_nn, cluster_area, total_error, total_error_not_norm, min_error, max_error, max_idx};
}

} // namespace vmas_troch

} // namespace geometry

} // namespace cgogn

#endif // CGOGN_GEOMETRY_TORCH_VMAS_TROCH_H_
