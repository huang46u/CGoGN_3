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
#ifndef CGOGN_GEOMETRY_ALGOS_OFFSET_AXIS_SAMPLE_PROCESSING_H_
#define CGOGN_GEOMETRY_ALGOS_OFFSET_AXIS_SAMPLE_PROCESSING_H_

#include <cgogn/core/functions/mesh_info.h>
#include <cgogn/core/functions/attributes.h>

#include <cgogn/geometry/algos/medial_axis.h>
#include <cgogn/geometry/types/line_quadric.h>
#include <cgogn/geometry/types/spherical_quadric.h>
#include <cgogn/geometry/types/vector_traits.h>

#include <Eigen/Eigenvalues>

#include <libacc/kd_tree.h>

#include <algorithm>
#include <cmath>
#include <memory>
#include <cstdint>
#include <utility>
#include <vector>

namespace cgogn
{

namespace geometry
{

namespace offset_axis
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
			int fallback_count = 0;
			for (auto& res : knn_res)
			{
				if (kdtree_vertices[res.first] != v)
				{
					const Vertex nb = kdtree_vertices[res.first];
					knn[v_idx].push_back(nb);
					sum_dist += res.second;
					++fallback_count;
					if (fallback_count >= knn_k + 1)
						break;
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

// kd-tree of the positions of the points; kdtree_vertices maps the kd-tree indices to the vertices (no tree when the
// point cloud is empty)
template <typename POINTS>
void build_point_kdtree(POINTS& points, const typename mesh_traits<POINTS>::template Attribute<Vec3>& position,
						std::unique_ptr<acc::KDTree<3, uint32>>& kdtree,
						std::vector<typename mesh_traits<POINTS>::Vertex>& kdtree_vertices)
{
	using Vertex = typename mesh_traits<POINTS>::Vertex;
	kdtree.reset();
	kdtree_vertices.clear();
	std::vector<Vec3> positions;
	foreach_cell(points, [&](Vertex v) {
		positions.push_back(position[index_of(points, v)]);
		kdtree_vertices.push_back(v);
		return true;
	});
	if (!positions.empty())
		kdtree = std::make_unique<acc::KDTree<3, uint32>>(positions);
}

// kd-tree of the medial points of the samples that have one (finite medial point, positive radius)
template <typename POINTS>
void build_medial_kdtree(POINTS& samples, const typename mesh_traits<POINTS>::template Attribute<Vec3>& ma_position,
						 const typename mesh_traits<POINTS>::template Attribute<Scalar>& ma_radius,
						 std::unique_ptr<acc::KDTree<3, uint32>>& kdtree,
						 std::vector<typename mesh_traits<POINTS>::Vertex>& kdtree_vertices)
{
	using Vertex = typename mesh_traits<POINTS>::Vertex;
	kdtree.reset();
	kdtree_vertices.clear();
	std::vector<Vec3> positions;
	foreach_cell(samples, [&](Vertex v) {
		const uint32 idx = index_of(samples, v);
		const Vec3 ma = ma_position[idx];
		if (ma.allFinite() && ma_radius[idx] > Scalar(0))
		{
			positions.push_back(ma);
			kdtree_vertices.push_back(v);
		}
		return true;
	});
	if (!positions.empty())
		kdtree = std::make_unique<acc::KDTree<3, uint32>>(positions);
}

namespace detail
{

template <typename POINTS>
void shrinking_ball_center(POINTS& samples, const typename mesh_traits<POINTS>::template Attribute<Vec3>& position,
						   const typename mesh_traits<POINTS>::template Attribute<Vec3>& normal,
						   const acc::KDTree<3, uint32>& kdtree,
						   const std::vector<typename mesh_traits<POINTS>::Vertex>& kdtree_vertices,
						   Scalar initial_radius, Scalar min_radius, typename mesh_traits<POINTS>::Vertex vertex,
						   typename mesh_traits<POINTS>::template Attribute<Vec3>& ma_position,
						   typename mesh_traits<POINTS>::template Attribute<Scalar>& ma_radius,
						   typename mesh_traits<POINTS>::template Attribute<typename mesh_traits<POINTS>::Vertex>&
							   ma_secondary)
{
	using Vertex = typename mesh_traits<POINTS>::Vertex;

	const Scalar min_norm = Scalar(1e-12);
	const uint32 idx = index_of(samples, vertex);
	const Vec3& pt = position[idx];
	Vec3 n = normal[idx];
	const bool can_shrink = n.allFinite() && n.squaredNorm() > min_norm;
	if (can_shrink)
		n.normalize();
	else
		n = Vec3(0, 0, 1);
	Vec3 center = pt - n * min_radius;
	Scalar radius = min_radius;
	Vertex secondary;
	if (can_shrink)
	{
		auto [c1, r1, q1] = geometry::shrinking_ball_center<Vertex>(pt, n, &kdtree, kdtree_vertices, initial_radius);
		if (c1.allFinite() && std::isfinite(static_cast<double>(r1)) && r1 > Scalar(0))
		{
			center = c1;
			radius = r1;
			secondary = q1;
		}
	}
	if (!secondary.is_valid())
	{
		std::pair<uint32, Scalar> nn;
		kdtree.find_nn(center - n * radius, &nn);
		secondary = kdtree_vertices[nn.first];
	}
	ma_position[idx] = center;
	ma_radius[idx] = std::max(radius, min_radius);
	ma_secondary[idx] = secondary;
}

} // namespace detail

// Shrinking ball centers of the given samples, whose positions are in the kd-tree (kdtree_vertices maps the kd-tree
// indices to the vertices). A ball that cannot be computed (degenerate normal, non finite or non positive result)
// falls back to the ball of radius min_radius tangent at the point, and the radius is at least min_radius. Without a
// second contact vertex, the vertex of the kd-tree nearest to the bottom of the ball is used.
template <typename POINTS>
void shrinking_ball_centers(POINTS& samples, const typename mesh_traits<POINTS>::template Attribute<Vec3>& position,
							const typename mesh_traits<POINTS>::template Attribute<Vec3>& normal,
							const acc::KDTree<3, uint32>& kdtree,
							const std::vector<typename mesh_traits<POINTS>::Vertex>& kdtree_vertices,
							Scalar initial_radius, Scalar min_radius,
							const std::vector<typename mesh_traits<POINTS>::Vertex>& vertices,
							typename mesh_traits<POINTS>::template Attribute<Vec3>& ma_position,
							typename mesh_traits<POINTS>::template Attribute<Scalar>& ma_radius,
							typename mesh_traits<POINTS>::template Attribute<typename mesh_traits<POINTS>::Vertex>&
								ma_secondary)
{
	for (auto v : vertices)
		detail::shrinking_ball_center(samples, position, normal, kdtree, kdtree_vertices, initial_radius, min_radius, v,
									  ma_position, ma_radius, ma_secondary);
}

// Shrinking ball centers of all the samples (see above), computed in parallel
template <typename POINTS>
void shrinking_ball_centers(POINTS& samples, const typename mesh_traits<POINTS>::template Attribute<Vec3>& position,
							const typename mesh_traits<POINTS>::template Attribute<Vec3>& normal,
							const acc::KDTree<3, uint32>& kdtree,
							const std::vector<typename mesh_traits<POINTS>::Vertex>& kdtree_vertices,
							Scalar initial_radius, Scalar min_radius,
							typename mesh_traits<POINTS>::template Attribute<Vec3>& ma_position,
							typename mesh_traits<POINTS>::template Attribute<Scalar>& ma_radius,
							typename mesh_traits<POINTS>::template Attribute<typename mesh_traits<POINTS>::Vertex>&
								ma_secondary)
{
	using Vertex = typename mesh_traits<POINTS>::Vertex;
	parallel_foreach_cell(samples, [&](Vertex v) {
		detail::shrinking_ball_center(samples, position, normal, kdtree, kdtree_vertices, initial_radius, min_radius, v,
									  ma_position, ma_radius, ma_secondary);
		return true;
	});
}

namespace detail
{

// PCA normal of the sample from its k nearest samples, oriented like its current normal (unchanged with fewer than 3
// neighbors)
template <typename POINTS>
void recompute_normal_pca(POINTS& samples, const typename mesh_traits<POINTS>::template Attribute<Vec3>& position,
						  typename mesh_traits<POINTS>::template Attribute<Vec3>& normal,
						  const acc::KDTree<3, uint32>& kdtree,
						  const std::vector<typename mesh_traits<POINTS>::Vertex>& kdtree_vertices, int knn_k,
						  typename mesh_traits<POINTS>::Vertex vertex)
{
	using Vertex = typename mesh_traits<POINTS>::Vertex;

	const int k = knn_k;
	const uint32 vertex_index = index_of(samples, vertex);
	const Vec3& center = position[vertex_index];
	std::vector<std::pair<uint32, Scalar>> neighbors_index;
	kdtree.find_nns(center, k + 1, &neighbors_index);
	std::vector<Vec3> neighbors;
	neighbors.reserve(k);
	for (const auto& item : neighbors_index)
	{
		const Vertex neighbor = kdtree_vertices[item.first];
		if (neighbor == vertex)
			continue;
		const uint32 neighbor_index = index_of(samples, neighbor);
		neighbors.push_back(position[neighbor_index]);
		if (static_cast<int>(neighbors.size()) >= k)
			break;
	}
	if (neighbors.size() < 3)
		return;
	Vec3 mean(0, 0, 0);
	for (const Vec3& p : neighbors)
		mean += p;
	mean /= Scalar(neighbors.size());
	Eigen::Matrix<Scalar, 3, 3> covariance = Eigen::Matrix<Scalar, 3, 3>::Zero();
	for (const Vec3& p : neighbors)
	{
		const Vec3 d = p - mean;
		covariance(0, 0) += d.x() * d.x();
		covariance(0, 1) += d.x() * d.y();
		covariance(0, 2) += d.x() * d.z();
		covariance(1, 1) += d.y() * d.y();
		covariance(1, 2) += d.y() * d.z();
		covariance(2, 2) += d.z() * d.z();
	}
	covariance(1, 0) = covariance(0, 1);
	covariance(2, 0) = covariance(0, 2);
	covariance(2, 1) = covariance(1, 2);
	Eigen::SelfAdjointEigenSolver<Eigen::Matrix<Scalar, 3, 3>> solver(covariance);
	const auto eigenvector = solver.eigenvectors().col(0);
	Vec3 n(eigenvector(0), eigenvector(1), eigenvector(2));
	if (n.dot(normal[vertex_index]) < Scalar(0))
		n = -n;
	n.normalize();
	normal[vertex_index] = n;
}

} // namespace detail

// PCA normals of all the samples from their k nearest samples, oriented like their current normals
template <typename POINTS>
void recompute_normals_pca(POINTS& samples, const typename mesh_traits<POINTS>::template Attribute<Vec3>& position,
						   typename mesh_traits<POINTS>::template Attribute<Vec3>& normal,
						   const acc::KDTree<3, uint32>& kdtree,
						   const std::vector<typename mesh_traits<POINTS>::Vertex>& kdtree_vertices, int knn_k)
{
	using Vertex = typename mesh_traits<POINTS>::Vertex;
	parallel_foreach_cell(samples, [&](Vertex v) {
		detail::recompute_normal_pca(samples, position, normal, kdtree, kdtree_vertices, knn_k, v);
		return true;
	});
}

} // namespace offset_axis

} // namespace geometry

} // namespace cgogn

#endif // CGOGN_GEOMETRY_ALGOS_OFFSET_AXIS_SAMPLE_PROCESSING_H_
