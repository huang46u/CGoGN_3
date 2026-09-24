/*******************************************************************************
 * CGoGN: Combinatorial and Geometric modeling with Generic N-dimensional Maps  *
 * Copyright (C), IGG Group, ICube, University of Strasbourg, France            *
 *                                                                              *
 * This library is free software; you can redistribute it and/or modify it      *
 * under the terms of the GNU Lesser General Public License as published by    *
 * the Free Software Foundation; either version 2.1 of the License, or (at     *
 * your option) any later version.                                              *
 *******************************************************************************/
#ifndef CGOGN_GEOMETRY_ALGOS_UDF_SAMPLE_GEOMETRY_H_
#define CGOGN_GEOMETRY_ALGOS_UDF_SAMPLE_GEOMETRY_H_

#include <cgogn/geometry/types/line_quadric.h>
#include <cgogn/geometry/types/quadric.h>
#include <cgogn/geometry/types/spherical_quadric.h>
#include <cgogn/geometry/types/vector_traits.h>

#include <Eigen/Eigenvalues>

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <utility>
#include <vector>

namespace cgogn
{
namespace geometry
{

struct SampleNeighborhoodParameters
{
	int knn_k = 10;
	Scalar alpha = Scalar(0);
};

template <typename MESH, typename POSITION, typename VERTICES>
Vec3 estimate_pca_normal(const MESH& mesh, const POSITION& position, const std::vector<uint32>& indices,
						 const VERTICES& kdtree_vertices)
{
	using Vertex = typename mesh_traits<MESH>::Vertex;
	if (indices.size() < 3)
		return Vec3(0, 0, 1);

	Vec3 centroid(0, 0, 0);
	for (uint32 idx : indices)
		centroid += position[index_of(mesh, kdtree_vertices[idx])];
	centroid /= Scalar(indices.size());

	Eigen::Matrix3d covariance = Eigen::Matrix3d::Zero();
	for (uint32 idx : indices)
	{
		const Vertex v = kdtree_vertices[idx];
		const Vec3 diff = position[index_of(mesh, v)] - centroid;
		const Eigen::Vector3d point(diff[0], diff[1], diff[2]);
		covariance += point * point.transpose();
	}
	covariance /= Scalar(indices.size());
	const Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> solver(covariance);
	const Eigen::Vector3d normal = solver.eigenvectors().col(0);
	return Vec3(normal[0], normal[1], normal[2]).normalized();
}

template <typename MESH, typename POSITION, typename NORMAL, typename KD_TREE, typename VERTICES>
void recompute_pca_normals(const MESH& mesh, const POSITION& position, NORMAL& normals, KD_TREE& kdtree,
						  const VERTICES& kdtree_vertices, int knn_k,
						  const std::vector<typename mesh_traits<MESH>::Vertex>* subset = nullptr)
{
	using Vertex = typename mesh_traits<MESH>::Vertex;
	const auto recompute = [&](Vertex v) {
		if (!v.is_valid())
			return;
		const uint32 v_idx = index_of(mesh, v);
		const Vec3 original = normals[v_idx];
		const Vec3& center = position[v_idx];
		std::vector<std::pair<uint32, Scalar>> neighbors;
		const int k = std::max(3, knn_k);
		kdtree.find_nns(center, k + 1, &neighbors);
		std::vector<Vec3> neighbor_positions;
		neighbor_positions.reserve(k);
		for (const auto& neighbor : neighbors)
		{
			const Vertex adjacent = kdtree_vertices[neighbor.first];
			if (adjacent == v)
				continue;
			neighbor_positions.push_back(position[index_of(mesh, adjacent)]);
			if (static_cast<int>(neighbor_positions.size()) >= k)
				break;
		}
		if (neighbor_positions.size() < 3)
			return;

		Vec3 mean(0, 0, 0);
		for (const Vec3& point : neighbor_positions)
			mean += point;
		mean /= Scalar(neighbor_positions.size());
		Eigen::Matrix<Scalar, 3, 3> covariance = Eigen::Matrix<Scalar, 3, 3>::Zero();
		for (const Vec3& point : neighbor_positions)
		{
			const Vec3 diff = point - mean;
			covariance(0, 0) += diff.x() * diff.x();
			covariance(0, 1) += diff.x() * diff.y();
			covariance(0, 2) += diff.x() * diff.z();
			covariance(1, 1) += diff.y() * diff.y();
			covariance(1, 2) += diff.y() * diff.z();
			covariance(2, 2) += diff.z() * diff.z();
		}
		covariance(1, 0) = covariance(0, 1);
		covariance(2, 0) = covariance(0, 2);
		covariance(2, 1) = covariance(1, 2);
		Eigen::SelfAdjointEigenSolver<Eigen::Matrix<Scalar, 3, 3>> solver(covariance);
		if (solver.info() != Eigen::Success)
			return;
		const Eigen::Matrix<Scalar, 3, 1> eigenvector = solver.eigenvectors().col(0);
		Vec3 normal(eigenvector(0), eigenvector(1), eigenvector(2));
		if (normal.squaredNorm() < Scalar(1e-12))
			return;
		if (original.squaredNorm() > Scalar(1e-12) && normal.dot(original) < Scalar(0))
			normal = -normal;
		normals[v_idx] = normal.normalized();
	};
	if (subset)
	{
		for (Vertex v : *subset)
			if (v.is_valid() && index_of(mesh, v) != INVALID_INDEX)
				recompute(v);
	}
	else
		foreach_cell(mesh, [&](Vertex v) { recompute(v); return true; });
}

template <typename MESH, typename POSITION, typename NORMAL, typename AREA, typename KNN, typename KD_TREE,
		  typename VERTICES>
void compute_sample_neighborhoods_and_areas(const MESH& mesh, const POSITION& positions, const NORMAL& normals,
												AREA& areas, KNN& knn, KD_TREE& kdtree,
												const VERTICES& kdtree_vertices,
												const SampleNeighborhoodParameters& parameters)
{
	using Vertex = typename mesh_traits<MESH>::Vertex;
	parallel_foreach_cell(mesh, [&](Vertex v) {
		const uint32 v_idx = index_of(mesh, v);
		const Vec3& point = positions[v_idx];
		std::vector<std::pair<uint32, Scalar>> neighbors;
		kdtree.find_nns(point, parameters.knn_k + 10, &neighbors);
		knn[v_idx].clear();
		Scalar sum_dist = Scalar(0);
		const Scalar eps = Scalar(1e-12);
		Vec3 normal = normals[v_idx];
		const Scalar norm2 = normal.squaredNorm();
		if (norm2 > eps)
			normal /= std::sqrt(norm2);
		else
			normal = Vec3(0, 0, 1);
		int kept = 0;
		for (const auto& neighbor : neighbors)
		{
			if (kdtree_vertices[neighbor.first] == v)
				continue;
			const Vertex adjacent = kdtree_vertices[neighbor.first];
			const Vec3& q = positions[index_of(mesh, adjacent)];
			if (std::abs((q - point).dot(normal)) > parameters.alpha)
				continue;
			knn[v_idx].push_back(adjacent);
			sum_dist += neighbor.second;
			if (++kept >= parameters.knn_k + 1)
				break;
		}
		if (kept == 0)
			for (const auto& neighbor : neighbors)
				if (kdtree_vertices[neighbor.first] != v)
				{
					knn[v_idx].push_back(kdtree_vertices[neighbor.first]);
					sum_dist += neighbor.second;
				}
		areas[v_idx] = (sum_dist * sum_dist) / (Scalar(2) * parameters.knn_k);
		return true;
	});
}

template <typename MESH, typename POSITION, typename NORMAL, typename AREA, typename KNN, typename QUADRIC,
		  typename LINE_QUADRIC>
void compute_sample_quadrics(const MESH& mesh, const POSITION& positions, const NORMAL& normals, const AREA& areas,
							 const KNN& knn, QUADRIC& quadrics, LINE_QUADRIC& line_quadrics, int knn_k)
{
	using Vertex = typename mesh_traits<MESH>::Vertex;
	parallel_foreach_cell(mesh, [&](Vertex v) {
		const uint32 v_idx = index_of(mesh, v);
		Spherical_Quadric& quadric = quadrics[v_idx];
		Line_Quadric& line_quadric = line_quadrics[v_idx];
		quadric.clear();
		line_quadric.clear();
		const Vec3& position = positions[v_idx];
		const Vec3& normal = normals[v_idx];
		const Scalar weight = areas[v_idx] / (knn_k + 1.0);
		quadric += Spherical_Quadric(Vec4(position.x(), position.y(), position.z(), 0),
									Vec4(normal.x(), normal.y(), normal.z(), 1)) * weight;
		line_quadric += Line_Quadric(position, normal) * weight;
		for (Vertex neighbor : knn[v_idx])
		{
			const uint32 neighbor_idx = index_of(mesh, neighbor);
			const Vec3& neighbor_position = positions[neighbor_idx];
			const Vec3& neighbor_normal = normals[neighbor_idx];
			const Scalar neighbor_weight = areas[neighbor_idx] / (knn_k + 1.0);
			quadric += Spherical_Quadric(Vec4(neighbor_position.x(), neighbor_position.y(), neighbor_position.z(), 0),
										 Vec4(neighbor_normal.x(), neighbor_normal.y(), neighbor_normal.z(), 1)) * neighbor_weight;
			line_quadric += Line_Quadric(neighbor_position, neighbor_normal) * neighbor_weight;
		}
		return true;
	});
}

} // namespace geometry
} // namespace cgogn

#endif // CGOGN_GEOMETRY_ALGOS_UDF_SAMPLE_GEOMETRY_H_
