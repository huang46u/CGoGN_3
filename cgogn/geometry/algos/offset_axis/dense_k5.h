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
#ifndef CGOGN_GEOMETRY_ALGOS_OFFSET_AXIS_DENSE_K5_H_
#define CGOGN_GEOMETRY_ALGOS_OFFSET_AXIS_DENSE_K5_H_

#include <Eigen/Dense>
#include <Eigen/Eigenvalues>
#include <cgogn/core/functions/attributes.h>
#include <cgogn/core/functions/mesh_info.h>
#include <cgogn/core/functions/traversals/edge.h>
#include <cgogn/core/functions/traversals/face.h>
#include <cgogn/core/functions/traversals/vertex.h>
#include <cgogn/core/types/cell_marker.h>
#include <cgogn/geometry/algos/offset_axis/skeleton_mesh_ops.h>
#include <cgogn/geometry/algos/offset_axis/skeleton_tets.h>
#include <cgogn/geometry/algos/offset_axis/spheres_optimizer.h>
#include <cgogn/geometry/types/vector_traits.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <limits>
#include <set>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

namespace cgogn
{

namespace geometry
{

namespace offset_axis
{

// The 5 spheres of a dense region of the skeleton (5 vertices pairwise adjacent, or with one missing edge).
template <typename POINTS>
struct K5Region
{
	std::array<typename mesh_traits<POINTS>::Vertex, 5> spheres;
};

struct DenseFiveVertexDetectionResult
{
	std::vector<std::array<uint32, 5>> independent_exact_k5_regions;
	std::vector<std::array<uint32, 5>> independent_k5_minus_1_regions;
};

template <typename NONMANIFOLD>
DenseFiveVertexDetectionResult detect_k5_and_k5_minus_1_cells(const NONMANIFOLD& skeleton,
															  const SkeletonTets<NONMANIFOLD>& tets)
{
	using SkeletonVertex = typename mesh_traits<NONMANIFOLD>::Vertex;
	using SkeletonEdge = typename mesh_traits<NONMANIFOLD>::Edge;
	using SkeletonFace = typename mesh_traits<NONMANIFOLD>::Face;

	DenseFiveVertexDetectionResult result;
	auto make_edge_key = [](uint32 a, uint32 b) -> uint64 {
		if (a > b)
			std::swap(a, b);
		return (uint64(a) << 32) | uint64(b);
	};
	std::vector<std::array<uint32, 5>> exact_k5_regions;
	std::vector<std::array<uint32, 5>> k5_minus_1_regions;

	std::unordered_map<uint32, std::unordered_set<uint32>> adjacency;
	adjacency.reserve(nb_cells<SkeletonVertex>(skeleton));
	std::unordered_set<uint64> edge_pair_to_id;
	edge_pair_to_id.reserve(nb_cells<SkeletonEdge>(skeleton) * 2 + 1);

	foreach_cell(skeleton, [&](SkeletonEdge e) -> bool {
		const std::vector<SkeletonVertex> vv = incident_vertices(skeleton, e);
		const uint32 a = index_of(skeleton, vv[0]);
		const uint32 b = index_of(skeleton, vv[1]);
		adjacency[a].insert(b);
		adjacency[b].insert(a);
		edge_pair_to_id.insert(make_edge_key(a, b));
		return true;
	});

	std::set<std::array<uint32, 5>> seen_region_keys;
	for (const auto& [tet_id, tet] : tets)
	{
		std::array<uint32, 4> base_vertices = {INVALID_INDEX, INVALID_INDEX, INVALID_INDEX, INVALID_INDEX};
		uint32 base_vertex_count = 0;
		for (SkeletonFace face : tet.faces)
		{
			for (SkeletonVertex vertex : incident_vertices(skeleton, face))
			{
				const uint32 vertex_id = index_of(skeleton, vertex);
				if (std::find(base_vertices.begin(), base_vertices.begin() + base_vertex_count, vertex_id) !=
					base_vertices.begin() + base_vertex_count)
					continue;
				if (base_vertex_count < 4)
					base_vertices[base_vertex_count++] = vertex_id;
			}
		}
		if (base_vertex_count != 4)
			continue;
		std::sort(base_vertices.begin(), base_vertices.end());

		std::unordered_map<uint32, uint32> fifth_vertex_connection_count;
		for (uint32 base_vertex : base_vertices)
		{
			const auto adjacency_it = adjacency.find(base_vertex);
			if (adjacency_it == adjacency.end())
				continue;
			for (uint32 neighbor : adjacency_it->second)
			{
				if (!std::binary_search(base_vertices.begin(), base_vertices.end(), neighbor))
					++fifth_vertex_connection_count[neighbor];
			}
		}

		for (const auto& candidate : fifth_vertex_connection_count)
		{
			if (candidate.second < 3)
				continue;

			std::array<uint32, 5> region = {base_vertices[0], base_vertices[1], base_vertices[2], base_vertices[3],
											candidate.first};
			std::sort(region.begin(), region.end());
			if (!seen_region_keys.insert(region).second)
				continue;

			const bool is_exact_k5 = candidate.second == 4;
			if (is_exact_k5)
				exact_k5_regions.push_back(region);
			else
				k5_minus_1_regions.push_back(region);
		}
	}

	std::unordered_map<uint64, uint32> edge_region_count;
	edge_region_count.reserve((exact_k5_regions.size() + k5_minus_1_regions.size()) * 10);
	auto count_region_edges = [&](const std::array<uint32, 5>& region) {
		for (uint32 i = 0; i < 5; ++i)
		{
			for (uint32 j = i + 1; j < 5; ++j)
			{
				const uint64 edge_key = make_edge_key(region[i], region[j]);
				if (edge_pair_to_id.find(edge_key) != edge_pair_to_id.end())
					++edge_region_count[edge_key];
			}
		}
	};
	for (const auto& region : exact_k5_regions)
		count_region_edges(region);
	for (const auto& region : k5_minus_1_regions)
		count_region_edges(region);

	auto append_independent_region = [&](const std::array<uint32, 5>& region, bool is_exact_k5) {
		for (uint32 i = 0; i < 5; ++i)
		{
			for (uint32 j = i + 1; j < 5; ++j)
			{
				const uint64 edge_key = make_edge_key(region[i], region[j]);
				const auto count_it = edge_region_count.find(edge_key);
				if (count_it != edge_region_count.end() && count_it->second > 1)
					return;
			}
		}

		if (is_exact_k5)
			result.independent_exact_k5_regions.push_back(region);
		else
			result.independent_k5_minus_1_regions.push_back(region);
	};
	for (const auto& region : exact_k5_regions)
		append_independent_region(region, true);
	for (const auto& region : k5_minus_1_regions)
		append_independent_region(region, false);

	return result;
}

template <typename POINTS>
bool is_live_sphere_vertex(const POINTS& spheres, typename mesh_traits<POINTS>::Vertex sphere, uint32& sphere_index)
{
	sphere_index = INVALID_INDEX;
	if (!sphere.is_valid() || sphere.dart_.index_ >= spheres.darts_.maximum_index())
		return false;
	sphere_index = index_of(spheres, sphere);
	return sphere_index != INVALID_INDEX;
}

template <typename NONMANIFOLD>
bool is_live_skeleton_vertex(const NONMANIFOLD& skeleton, typename mesh_traits<NONMANIFOLD>::Vertex vertex,
							 uint32& vertex_index)
{
	vertex_index = INVALID_INDEX;
	if (!vertex.is_valid())
		return false;
	vertex_index = index_of(skeleton, vertex);
	return vertex_index != INVALID_INDEX;
}

// The dense regions that share no edge with another dense region, found at the start of the topology fix and
// repaired at its end by repair_k5_regions.
template <typename POINTS, typename NONMANIFOLD>
std::vector<K5Region<POINTS>> find_k5_regions(
	const NONMANIFOLD& skeleton, const SkeletonTets<NONMANIFOLD>& tets, const POINTS& spheres,
	const typename mesh_traits<NONMANIFOLD>::template Attribute<typename mesh_traits<POINTS>::Vertex>*
		skeleton_source_sphere)
{
	using SphereVertex = typename mesh_traits<POINTS>::Vertex;
	using SkeletonVertex = typename mesh_traits<NONMANIFOLD>::Vertex;

	const DenseFiveVertexDetectionResult detection = detect_k5_and_k5_minus_1_cells(skeleton, tets);
	std::vector<K5Region<POINTS>> seeds;
	seeds.reserve(detection.independent_exact_k5_regions.size() + detection.independent_k5_minus_1_regions.size());
	auto append_regions = [&](const std::vector<std::array<uint32, 5>>& regions) {
		for (const auto& region : regions)
		{
			K5Region<POINTS> seed;
			bool valid = true;
			for (uint32 i = 0; i < 5; ++i)
			{
				const SkeletonVertex vertex = of_index<SkeletonVertex>(skeleton, region[i]);
				uint32 current_vertex_index = INVALID_INDEX;
				if (!is_live_skeleton_vertex(skeleton, vertex, current_vertex_index) ||
					current_vertex_index != region[i])
				{
					valid = false;
					break;
				}
				const SphereVertex sphere = value<SphereVertex>(skeleton, skeleton_source_sphere, vertex);
				uint32 sphere_index = INVALID_INDEX;
				if (!is_live_sphere_vertex(spheres, sphere, sphere_index))
				{
					valid = false;
					break;
				}
				seed.spheres[i] = sphere;
			}
			if (valid)
				seeds.push_back(seed);
		}
	};
	append_regions(detection.independent_exact_k5_regions);
	append_regions(detection.independent_k5_minus_1_regions);
	return seeds;
}

template <typename NONMANIFOLD>
bool face_is_inside_vertex_set(const NONMANIFOLD& skeleton, typename mesh_traits<NONMANIFOLD>::Face face,
							   const CellMarkerStore<NONMANIFOLD, typename mesh_traits<NONMANIFOLD>::Vertex>& vertices)
{
	using SkeletonVertex = typename mesh_traits<NONMANIFOLD>::Vertex;

	for (SkeletonVertex vertex : incident_vertices(skeleton, face))
		if (!vertices.is_marked(vertex))
			return false;
	return true;
}

// Replaces the faces inside the region by a fan around a new sphere at the centroid of the 5 spheres, when the region
// still exists and is flat enough. Returns true if the region is repaired.
template <typename POINTS, typename NONMANIFOLD>
bool repair_k5_region(
	NONMANIFOLD& skeleton, SkeletonTets<NONMANIFOLD>& tets, SpheresOptimizer<POINTS>& spheres_optimizer,
	typename mesh_traits<NONMANIFOLD>::template Attribute<Vec3>* skeleton_position,
	typename mesh_traits<NONMANIFOLD>::template Attribute<Scalar>* skeleton_radius,
	typename mesh_traits<NONMANIFOLD>::template Attribute<typename mesh_traits<POINTS>::Vertex>* skeleton_source_sphere,
	typename mesh_traits<POINTS>::template Attribute<typename mesh_traits<NONMANIFOLD>::Vertex>* sphere_skeleton_vertex,
	const K5Region<POINTS>& seed, std::vector<typename mesh_traits<POINTS>::Vertex>& added_spheres)
{
	using SphereVertex = typename mesh_traits<POINTS>::Vertex;
	using SkeletonVertex = typename mesh_traits<NONMANIFOLD>::Vertex;
	using SkeletonEdge = typename mesh_traits<NONMANIFOLD>::Edge;
	using SkeletonFace = typename mesh_traits<NONMANIFOLD>::Face;

	auto& optimizer_data = spheres_optimizer.data();
	std::array<uint32, 5> sphere_ids;
	std::array<SkeletonVertex, 5> vertices;
	std::array<uint32, 5> vertex_ids;
	std::array<Vec3, 5> positions;
	Scalar average_radius = Scalar(0);
	Vec3 centroid = Vec3::Zero();
	CellMarkerStore<NONMANIFOLD, SkeletonVertex> region_vertices(skeleton);
	for (uint32 i = 0; i < 5; ++i)
	{
		if (!is_live_sphere_vertex(optimizer_data.spheres, seed.spheres[i], sphere_ids[i]))
			return false;
		vertices[i] = (*sphere_skeleton_vertex)[sphere_ids[i]];
		if (!is_live_skeleton_vertex(skeleton, vertices[i], vertex_ids[i]))
			return false;
		if (value<SphereVertex>(skeleton, skeleton_source_sphere, vertices[i]) != seed.spheres[i])
			return false;
		positions[i] = (*optimizer_data.sphere_position)[sphere_ids[i]];
		centroid += positions[i];
		average_radius += (*optimizer_data.sphere_radius)[sphere_ids[i]];
		// the 5 vertices are distinct
		if (region_vertices.is_marked(vertices[i]))
			return false;
		region_vertices.mark(vertices[i]);
	}
	centroid /= Scalar(5);
	average_radius /= Scalar(5);
	if (!centroid.allFinite() || !std::isfinite(static_cast<double>(average_radius)) || average_radius <= Scalar(0))
		return false;

	auto make_edge_key = [](uint32 a, uint32 b) -> uint64 {
		if (a > b)
			std::swap(a, b);
		return (uint64(a) << 32) | uint64(b);
	};
	std::unordered_map<uint64, SkeletonEdge> region_edges;
	region_edges.reserve(16);
	std::unordered_set<uint32> region_nm_edge_ids;
	region_nm_edge_ids.reserve(16);
	for (uint32 i = 0; i < 5; ++i)
	{
		for (uint32 j = i + 1; j < 5; ++j)
		{
			const SkeletonEdge edge = find_edge(skeleton, vertices[i], vertices[j]);
			if (!edge.is_valid())
				continue;
			region_edges[make_edge_key(vertex_ids[i], vertex_ids[j])] = edge;
			if (incident_faces(skeleton, edge).size() > 2)
				region_nm_edge_ids.insert(index_of(skeleton, edge));
		}
	}
	if (region_nm_edge_ids.empty())
		return false;

	for (uint32 nm_edge_id : region_nm_edge_ids)
	{
		const SkeletonEdge nm_edge = of_index<SkeletonEdge>(skeleton, nm_edge_id);
		for (SkeletonVertex endpoint : incident_vertices(skeleton, nm_edge))
		{
			for (SkeletonEdge adjacent_edge : incident_edges(skeleton, endpoint))
			{
				if (incident_faces(skeleton, adjacent_edge).size() <= 2)
					continue;
				if (region_nm_edge_ids.find(index_of(skeleton, adjacent_edge)) == region_nm_edge_ids.end())
					return false;
			}
		}
	}

	Eigen::Matrix<Scalar, 3, 3> covariance = Eigen::Matrix<Scalar, 3, 3>::Zero();
	for (const Vec3& position : positions)
	{
		const Vec3 delta = position - centroid;
		covariance += delta * delta.transpose();
	}
	Eigen::SelfAdjointEigenSolver<Eigen::Matrix<Scalar, 3, 3>> eigen_solver(covariance);
	if (eigen_solver.info() != Eigen::Success ||
		eigen_solver.eigenvalues()[1] <= std::numeric_limits<Scalar>::epsilon())
		return false;
	const Vec3 axis_u = eigen_solver.eigenvectors().col(2);
	const Vec3 axis_v = eigen_solver.eigenvectors().col(1);
	std::array<uint32, 5> cyclic_order = {0, 1, 2, 3, 4};
	std::sort(cyclic_order.begin(), cyclic_order.end(), [&](uint32 lhs, uint32 rhs) {
		const Vec3 lhs_delta = positions[lhs] - centroid;
		const Vec3 rhs_delta = positions[rhs] - centroid;
		const Scalar lhs_angle = std::atan2(lhs_delta.dot(axis_v), lhs_delta.dot(axis_u));
		const Scalar rhs_angle = std::atan2(rhs_delta.dot(axis_v), rhs_delta.dot(axis_u));
		return lhs_angle < rhs_angle;
	});
	Scalar max_radial_squared = Scalar(0);
	for (const Vec3& position : positions)
		max_radial_squared = std::max(max_radial_squared, (position - centroid).squaredNorm());
	const Scalar fan_area_epsilon =
		std::max(std::numeric_limits<Scalar>::epsilon(), max_radial_squared * Scalar(1e-10));
	for (uint32 i = 0; i < 5; ++i)
	{
		const Vec3 radial_a = positions[cyclic_order[i]] - centroid;
		const Vec3 radial_b = positions[cyclic_order[(i + 1) % 5]] - centroid;
		if (radial_a.cross(radial_b).norm() <= fan_area_epsilon)
			return false;
	}

	std::unordered_set<uint64> cycle_edge_keys;
	cycle_edge_keys.reserve(8);
	for (uint32 i = 0; i < 5; ++i)
	{
		const uint32 a = cyclic_order[i];
		const uint32 b = cyclic_order[(i + 1) % 5];
		cycle_edge_keys.insert(make_edge_key(vertex_ids[a], vertex_ids[b]));
	}

	std::unordered_set<uint32> internal_face_ids;
	internal_face_ids.reserve(16);
	foreach_cell(skeleton, [&](SkeletonFace face) -> bool {
		if (face_is_inside_vertex_set(skeleton, face, region_vertices))
			internal_face_ids.insert(index_of(skeleton, face));
		return true;
	});
	if (internal_face_ids.empty())
		return false;

	for (const auto& edge_entry : region_edges)
	{
		uint32 external_face_count = 0;
		for (SkeletonFace face : incident_faces(skeleton, edge_entry.second))
		{
			if (internal_face_ids.find(index_of(skeleton, face)) == internal_face_ids.end())
				++external_face_count;
		}
		if (cycle_edge_keys.find(edge_entry.first) == cycle_edge_keys.end())
		{
			if (external_face_count != 0)
				return false;
		}
		else if (external_face_count > 1)
		{
			return false;
		}
	}

	// the faces and their tets only: the region edges are handled below
	for (uint32 face_id : internal_face_ids)
	{
		const std::set<std::size_t> face_tets = tets.tets_of(face_id);
		remove_face(skeleton, of_index<SkeletonFace>(skeleton, face_id));
		for (std::size_t tet_id : face_tets)
			tets.erase(tet_id);
	}

	for (const auto& edge_entry : region_edges)
	{
		if (cycle_edge_keys.find(edge_entry.first) != cycle_edge_keys.end())
			continue;
		if (incident_faces(skeleton, edge_entry.second).empty())
			remove_edge(skeleton, edge_entry.second);
	}

	for (uint32 i = 0; i < 5; ++i)
	{
		for (uint32 j = i + 1; j < 5; ++j)
		{
			const uint64 edge_key = make_edge_key(vertex_ids[i], vertex_ids[j]);
			const bool is_cycle_edge = cycle_edge_keys.find(edge_key) != cycle_edge_keys.end();
			if (is_cycle_edge)
			{
				sorted_vector_insert((*optimizer_data.sphere_neighbors)[sphere_ids[i]], seed.spheres[j]);
				sorted_vector_insert((*optimizer_data.sphere_neighbors)[sphere_ids[j]], seed.spheres[i]);
			}
			else
			{
				sorted_vector_erase((*optimizer_data.sphere_neighbors)[sphere_ids[i]], seed.spheres[j]);
				sorted_vector_erase((*optimizer_data.sphere_neighbors)[sphere_ids[j]], seed.spheres[i]);
			}
		}
	}

	const SphereVertex center_sphere = spheres_optimizer.add_sphere(centroid, average_radius);
	if (!center_sphere.is_valid())
		return false;
	const uint32 center_sphere_id = index_of(optimizer_data.spheres, center_sphere);
	added_spheres.push_back(center_sphere);

	const SkeletonVertex center_vertex = add_vertex(skeleton);
	value<Vec3>(skeleton, skeleton_position, center_vertex) = centroid;
	value<Scalar>(skeleton, skeleton_radius, center_vertex) = average_radius;
	value<SphereVertex>(skeleton, skeleton_source_sphere, center_vertex) = center_sphere;
	(*sphere_skeleton_vertex)[center_sphere_id] = center_vertex;

	std::array<SkeletonEdge, 5> cycle_edges;
	std::array<SkeletonEdge, 5> spoke_edges;
	for (uint32 i = 0; i < 5; ++i)
	{
		const uint32 a = cyclic_order[i];
		const uint32 b = cyclic_order[(i + 1) % 5];
		cycle_edges[i] = find_edge(skeleton, vertices[a], vertices[b]);
		if (!cycle_edges[i].is_valid())
			cycle_edges[i] = add_edge(skeleton, vertices[a], vertices[b]);
		spoke_edges[i] = add_edge(skeleton, center_vertex, vertices[a]);
		sorted_vector_insert((*optimizer_data.sphere_neighbors)[center_sphere_id], seed.spheres[a]);
		sorted_vector_insert((*optimizer_data.sphere_neighbors)[sphere_ids[a]], center_sphere);
	}

	for (uint32 i = 0; i < 5; ++i)
	{
		std::vector<SkeletonEdge> face_edges = {cycle_edges[i], spoke_edges[i], spoke_edges[(i + 1) % 5]};
		const SkeletonFace new_face = add_face(skeleton, face_edges);
		tets.reset_face(new_face);
	}

	return true;
}

// Repairs the regions of find_k5_regions that still exist; the new spheres are appended to added_spheres.
template <typename POINTS, typename NONMANIFOLD>
void repair_k5_regions(
	NONMANIFOLD& skeleton, SkeletonTets<NONMANIFOLD>& tets, SpheresOptimizer<POINTS>& spheres_optimizer,
	typename mesh_traits<NONMANIFOLD>::template Attribute<Vec3>* skeleton_position,
	typename mesh_traits<NONMANIFOLD>::template Attribute<Scalar>* skeleton_radius,
	typename mesh_traits<NONMANIFOLD>::template Attribute<typename mesh_traits<POINTS>::Vertex>* skeleton_source_sphere,
	typename mesh_traits<POINTS>::template Attribute<typename mesh_traits<NONMANIFOLD>::Vertex>* sphere_skeleton_vertex,
	typename mesh_traits<NONMANIFOLD>::template Attribute<uint32>* edge_degree,
	const std::vector<K5Region<POINTS>>& seeds, std::vector<typename mesh_traits<POINTS>::Vertex>& added_spheres)
{
	using SphereVertex = typename mesh_traits<POINTS>::Vertex;
	using SkeletonVertex = typename mesh_traits<NONMANIFOLD>::Vertex;

	auto& optimizer_data = spheres_optimizer.data();
	std::vector<K5Region<POINTS>> live_seeds;
	live_seeds.reserve(seeds.size());
	for (const auto& seed : seeds)
	{
		bool all_vertices_live = true;
		CellMarkerStore<NONMANIFOLD, SkeletonVertex> seed_vertices(skeleton);
		for (const SphereVertex& sphere : seed.spheres)
		{
			uint32 sphere_id = INVALID_INDEX;
			if (!is_live_sphere_vertex(optimizer_data.spheres, sphere, sphere_id))
			{
				all_vertices_live = false;
				break;
			}
			const SkeletonVertex vertex = (*sphere_skeleton_vertex)[sphere_id];
			uint32 vertex_id = INVALID_INDEX;
			if (!is_live_skeleton_vertex(skeleton, vertex, vertex_id) ||
				value<SphereVertex>(skeleton, skeleton_source_sphere, vertex) != sphere ||
				seed_vertices.is_marked(vertex))
			{
				all_vertices_live = false;
				break;
			}
			seed_vertices.mark(vertex);
		}
		if (all_vertices_live)
			live_seeds.push_back(seed);
	}

	bool repaired_any = false;
	for (const auto& seed : live_seeds)
		if (repair_k5_region(skeleton, tets, spheres_optimizer, skeleton_position, skeleton_radius,
							 skeleton_source_sphere, sphere_skeleton_vertex, seed, added_spheres))
			repaired_any = true;

	if (repaired_any)
		compute_edge_degree(skeleton, edge_degree);
}

} // namespace offset_axis

} // namespace geometry

} // namespace cgogn

#endif // CGOGN_GEOMETRY_ALGOS_OFFSET_AXIS_DENSE_K5_H_
