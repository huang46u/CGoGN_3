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
 * along with this library; if not, write to the Free Software Foundation,     *
 * Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301 USA.           *
 *                                                                              *
 * Web site: http://cgogn.unistra.fr/                                           *
 * Contact information: cgogn@unistra.fr                                        *
 *                                                                              *
 *******************************************************************************/
#ifndef CGOGN_GEOMETRY_ALGOS_UDF_SAMPLE_PROCESSING_H_
#define CGOGN_GEOMETRY_ALGOS_UDF_SAMPLE_PROCESSING_H_

#include <cgogn/core/functions/mesh_info.h>
#include <cgogn/core/functions/attributes.h>

#include <cgogn/geometry/types/line_quadric.h>
#include <cgogn/geometry/types/spherical_quadric.h>
#include <cgogn/geometry/types/vector_traits.h>

#include <Eigen/Eigenvalues>

#include <libacc/kd_tree.h>

#include <cmath>
#include <cstdint>
#include <utility>
#include <vector>

namespace cgogn
{

namespace geometry
{

namespace detail
{

template <typename POINTS>
Vec3 compute_pca_normal(
	const POINTS& points,
	const typename mesh_traits<POINTS>::template Attribute<Vec3>& position,
	const std::vector<uint32>& neighbor_indices,
	const std::vector<typename mesh_traits<POINTS>::Vertex>& kdtree_vertices)
{
	using Vertex = typename mesh_traits<POINTS>::Vertex;
	if (neighbor_indices.size() < 3)
		return Vec3(0, 0, 1);

	Vec3 centroid(0, 0, 0);
	for (uint32 idx : neighbor_indices)
		centroid += position[index_of(points, kdtree_vertices[idx])];
	centroid /= Scalar(neighbor_indices.size());

	Eigen::Matrix3d covariance = Eigen::Matrix3d::Zero();
	for (uint32 idx : neighbor_indices)
	{
		const Vertex v = kdtree_vertices[idx];
		const Vec3 diff = position[index_of(points, v)] - centroid;
		const Eigen::Vector3d point(diff[0], diff[1], diff[2]);
		covariance += point * point.transpose();
	}
	covariance /= Scalar(neighbor_indices.size());
	const Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> solver(covariance);
	const Eigen::Vector3d normal = solver.eigenvectors().col(0);
	return Vec3(normal[0], normal[1], normal[2]).normalized();
}

} // namespace detail

template <typename POINTS>
void compute_input_normals(
	POINTS& points,
	const typename mesh_traits<POINTS>::template Attribute<Vec3>& position,
	typename mesh_traits<POINTS>::template Attribute<Vec3>& normal,
	typename mesh_traits<POINTS>::template Attribute<std::vector<typename mesh_traits<POINTS>::Vertex>>& knn,
	const acc::KDTree<3, uint32>& kdtree,
	const std::vector<typename mesh_traits<POINTS>::Vertex>& kdtree_vertices,
	int knn_k)
{
	using Vertex = typename mesh_traits<POINTS>::Vertex;
	parallel_foreach_cell(points, [&](Vertex v) -> bool {
		const uint32 v_idx = index_of(points, v);
		const Vec3& pt = position[v_idx];
		std::vector<std::pair<uint32, Scalar>> knn_res;
		kdtree.find_nns(pt, knn_k + 1, &knn_res);
		std::vector<uint32> indices;
		knn[v_idx].clear();
		for (const auto& res : knn_res)
		{
			indices.push_back(res.first);
			knn[v_idx].push_back(kdtree_vertices[res.first]);
		}
		normal[v_idx] = detail::compute_pca_normal(points, position, indices, kdtree_vertices);
		return true;
	});
}

template <typename POINTS>
void compute_samples_area(
	POINTS& samples,
	const typename mesh_traits<POINTS>::template Attribute<Vec3>& position,
	const typename mesh_traits<POINTS>::template Attribute<Vec3>& normal,
	typename mesh_traits<POINTS>::template Attribute<Scalar>& area,
	typename mesh_traits<POINTS>::template Attribute<std::vector<typename mesh_traits<POINTS>::Vertex>>& knn,
	const acc::KDTree<3, uint32>& kdtree,
	const std::vector<typename mesh_traits<POINTS>::Vertex>& kdtree_vertices,
	int knn_k,
	Scalar alpha)
{
	using Vertex = typename mesh_traits<POINTS>::Vertex;
	parallel_foreach_cell(samples, [&](Vertex v) -> bool {
		const uint32 v_idx = index_of(samples, v);
		const Vec3& pt = position[v_idx];
		std::vector<std::pair<uint32, Scalar>> knn_res;
		kdtree.find_nns(pt, knn_k + 10, &knn_res);

		knn[v_idx].clear();
		Scalar sum_dist = 0.0;
		const Scalar eps = Scalar(1e-12);
		const Scalar band = alpha;
		Vec3 n = normal[v_idx];
		const Scalar n2 = n.squaredNorm();
		if (n2 > eps)
			n /= std::sqrt(n2);
		else
			n = Vec3(0, 0, 1);
		int kept = 0;
		for (auto& res : knn_res)
		{
			if (kdtree_vertices[res.first] != v)
			{
				const Vertex nb = kdtree_vertices[res.first];
				const uint32 nb_idx = index_of(samples, nb);
				const Vec3& q = position[nb_idx];
				const Scalar dn = (q - pt).dot(n);
				if (std::abs(dn) > band)
					continue;
				knn[v_idx].push_back(nb);
				sum_dist += res.second;
				++kept;
				if (kept >= knn_k + 1)
					break;
			}
		}
		if (kept == 0)
		{
			for (auto& res : knn_res)
			{
				if (kdtree_vertices[res.first] != v)
				{
					const Vertex nb = kdtree_vertices[res.first];
					knn[v_idx].push_back(nb);
					sum_dist += res.second;
				}
			}
		}
		area[v_idx] = (sum_dist * sum_dist) / (2.0 * knn_k);
		return true;
	});
}

template <typename POINTS>
void compute_quadrics(
	POINTS& samples,
	const typename mesh_traits<POINTS>::template Attribute<Vec3>& position,
	const typename mesh_traits<POINTS>::template Attribute<Vec3>& normal,
	const typename mesh_traits<POINTS>::template Attribute<Scalar>& area,
	const typename mesh_traits<POINTS>::template Attribute<std::vector<typename mesh_traits<POINTS>::Vertex>>& knn,
	typename mesh_traits<POINTS>::template Attribute<Spherical_Quadric>& spherical_quadric,
	typename mesh_traits<POINTS>::template Attribute<Line_Quadric>& line_quadric,
	int knn_k)
{
	using Vertex = typename mesh_traits<POINTS>::Vertex;
	parallel_foreach_cell(samples, [&](Vertex v) -> bool {
		const uint32 v_idx = index_of(samples, v);
		Spherical_Quadric& q = spherical_quadric[v_idx];
		Line_Quadric& lq = line_quadric[v_idx];
		q.clear();
		lq.clear();
		const Vec3& pos = position[v_idx];
		const Vec3& n = normal[v_idx];
		const Scalar a = area[v_idx] / (knn_k + 1.0);
		q += Spherical_Quadric(Vec4(pos.x(), pos.y(), pos.z(), 0), Vec4(n.x(), n.y(), n.z(), 1)) * a;
		lq += Line_Quadric(pos, n) * a;
		for (Vertex vn : knn[v_idx])
		{
			const uint32 vn_idx = index_of(samples, vn);
			const Vec3& pn = position[vn_idx];
			const Vec3& nn = normal[vn_idx];
			const Scalar an = area[vn_idx] / (knn_k + 1.0);
			q += Spherical_Quadric(Vec4(pn.x(), pn.y(), pn.z(), 0), Vec4(nn.x(), nn.y(), nn.z(), 1)) * an;
			lq += Line_Quadric(pn, nn) * an;
		}
		return true;
	});
}

} // namespace geometry

} // namespace cgogn

#endif // CGOGN_GEOMETRY_ALGOS_UDF_SAMPLE_PROCESSING_H_
