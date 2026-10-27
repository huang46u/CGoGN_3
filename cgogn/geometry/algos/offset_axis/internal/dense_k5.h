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
#ifndef CGOGN_GEOMETRY_ALGOS_OFFSET_AXIS_INTERNAL_DENSE_K5_H_
#define CGOGN_GEOMETRY_ALGOS_OFFSET_AXIS_INTERNAL_DENSE_K5_H_

#include <Eigen/Dense>
#include <Eigen/Eigenvalues>
#include <cgogn/core/functions/attributes.h>
#include <cgogn/core/functions/mesh_info.h>
#include <cgogn/core/functions/traversals/edge.h>
#include <cgogn/core/functions/traversals/face.h>
#include <cgogn/core/functions/traversals/vertex.h>
#include <cgogn/core/types/cell_marker.h>
#include <cgogn/geometry/algos/offset_axis/internal/skeleton_mesh_ops.h>
#include <cgogn/geometry/algos/offset_axis/internal/skeleton_topology_types.h>

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

namespace topology_detail
{

template <typename POINTS, typename NONMANIFOLD>
struct DenseK5Repair
{
	using Types = SkeletonTopologyTypes<POINTS, NONMANIFOLD>;
	using SphereVertex = typename Types::SphereVertex;
	using SkeletonVertex = typename Types::SkeletonVertex;
	using SkeletonEdge = typename Types::SkeletonEdge;
	using SkeletonFace = typename Types::SkeletonFace;
	using Context = typename Types::Context;
	using Tet = typename Types::Tet;
	using RemovalStats = typename Types::RemovalStats;
	using SkeletonFaceDeletionOptions = typename Types::SkeletonFaceDeletionOptions;
	using MeshOps = topology_detail::SkeletonMeshOps<POINTS, NONMANIFOLD>;

	struct DenseFiveVertexDetectionResult
	{
		std::vector<std::array<uint32, 5>> independent_exact_k5_regions;
		std::vector<std::array<uint32, 5>> independent_k5_minus_1_regions;
	};

	static DenseFiveVertexDetectionResult detect_k5_and_k5_minus_1_cells(Context& p)
	{
		DenseFiveVertexDetectionResult result;
		auto make_edge_key = [](uint32 a, uint32 b) -> uint64 {
			if (a > b)
				std::swap(a, b);
			return (uint64(a) << 32) | uint64(b);
		};
		std::vector<std::array<uint32, 5>> exact_k5_regions;
		std::vector<std::array<uint32, 5>> k5_minus_1_regions;

		std::unordered_map<uint32, std::unordered_set<uint32>> adjacency;
		adjacency.reserve(nb_cells<SkeletonVertex>(p.data.skeleton));
		std::unordered_set<uint64> edge_pair_to_id;
		edge_pair_to_id.reserve(nb_cells<SkeletonEdge>(p.data.skeleton) * 2 + 1);

		foreach_cell(p.data.skeleton, [&](SkeletonEdge e) -> bool {
			const uint32 ide = index_of(p.data.skeleton, e);
			if (ide == INVALID_INDEX)
				return true;
			const std::vector<SkeletonVertex> vv = incident_vertices(p.data.skeleton, e);
			if (vv.size() != 2)
				return true;
			const uint32 a = index_of(p.data.skeleton, vv[0]);
			const uint32 b = index_of(p.data.skeleton, vv[1]);
			if (a == INVALID_INDEX || b == INVALID_INDEX || a == b)
				return true;
			adjacency[a].insert(b);
			adjacency[b].insert(a);
			edge_pair_to_id.insert(make_edge_key(a, b));
			return true;
		});

		std::set<std::array<uint32, 5>> seen_region_keys;
		for (const auto& tet_entry : p.skeleton_tets)
		{
			const Tet& tet = tet_entry.second;
			std::array<uint32, 4> base_vertices = {INVALID_INDEX, INVALID_INDEX, INVALID_INDEX, INVALID_INDEX};
			uint32 base_vertex_count = 0;
			for (uint32 face_index = 0; face_index < 4; ++face_index)
			{
				const SkeletonFace face = tet.faces[face_index];
				if (!face.is_valid())
					continue;
				for (SkeletonVertex vertex : incident_vertices(p.data.skeleton, face))
				{
					const uint32 vertex_id = index_of(p.data.skeleton, vertex);
					if (vertex_id == INVALID_INDEX ||
						std::find(base_vertices.begin(), base_vertices.begin() + base_vertex_count, vertex_id) !=
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

	struct IndependentDenseFiveVertexRepairSeed
	{
		std::array<SphereVertex, 5> spheres;
	};


	static bool is_live_sphere_vertex(Context& p, const SphereVertex& sphere, uint32& sphere_index)
	{
		sphere_index = INVALID_INDEX;
		if (!sphere.is_valid() || sphere.dart_.index_ >= p.optimizer_data.spheres.darts_.maximum_index())
			return false;
		sphere_index = index_of(p.optimizer_data.spheres, sphere);
		return sphere_index != INVALID_INDEX;
	}

	static bool is_live_skeleton_vertex(Context& p, const SkeletonVertex& vertex, uint32& vertex_index)
	{
		vertex_index = INVALID_INDEX;
		if (!vertex.is_valid())
			return false;
		vertex_index = index_of(p.data.skeleton, vertex);
		return vertex_index != INVALID_INDEX;
	}

	static std::vector<IndependentDenseFiveVertexRepairSeed> collect_independent_dense_five_vertex_repair_seeds(
		Context& p, const DenseFiveVertexDetectionResult& detection)
	{
		std::vector<IndependentDenseFiveVertexRepairSeed> seeds;
		seeds.reserve(detection.independent_exact_k5_regions.size() + detection.independent_k5_minus_1_regions.size());
		auto append_regions = [&](const std::vector<std::array<uint32, 5>>& regions) {
			for (const auto& region : regions)
			{
				IndependentDenseFiveVertexRepairSeed seed;
				bool valid = true;
				for (uint32 i = 0; i < 5; ++i)
				{
					const SkeletonVertex vertex = of_index<SkeletonVertex>(p.data.skeleton, region[i]);
					uint32 current_vertex_index = INVALID_INDEX;
					if (!is_live_skeleton_vertex(p, vertex, current_vertex_index) || current_vertex_index != region[i])
					{
						valid = false;
						break;
					}
					const SphereVertex sphere = value<SphereVertex>(p.data.skeleton, p.data.skeleton_source_sphere, vertex);
					uint32 sphere_index = INVALID_INDEX;
					if (!is_live_sphere_vertex(p, sphere, sphere_index))
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

	static bool face_is_inside_vertex_set(Context& p, const SkeletonFace& face,
										  const std::unordered_set<uint32>& vertex_ids)
	{
		if (!face.is_valid())
			return false;
		const std::vector<SkeletonVertex> vertices = incident_vertices(p.data.skeleton, face);
		if (vertices.size() != 3)
			return false;
		for (SkeletonVertex vertex : vertices)
		{
			const uint32 vertex_id = index_of(p.data.skeleton, vertex);
			if (vertex_id == INVALID_INDEX || vertex_ids.find(vertex_id) == vertex_ids.end())
				return false;
		}
		return true;
	}

	// Returns true if the region is repaired.
	static bool repair_independent_dense_five_vertex_region(Context& p, const IndependentDenseFiveVertexRepairSeed& seed)
	{
		std::array<uint32, 5> sphere_ids;
		std::array<SkeletonVertex, 5> vertices;
		std::array<uint32, 5> vertex_ids;
		std::array<Vec3, 5> positions;
		Scalar average_radius = Scalar(0);
		Vec3 centroid = Vec3::Zero();
		std::unordered_set<uint32> region_vertex_ids;
		region_vertex_ids.reserve(8);
		for (uint32 i = 0; i < 5; ++i)
		{
			if (!is_live_sphere_vertex(p, seed.spheres[i], sphere_ids[i]))
				return false;
			vertices[i] = (*p.data.sphere_skeleton_vertex)[sphere_ids[i]];
			if (!is_live_skeleton_vertex(p, vertices[i], vertex_ids[i]))
				return false;
			if (value<SphereVertex>(p.data.skeleton, p.data.skeleton_source_sphere, vertices[i]) != seed.spheres[i])
				return false;
			positions[i] = (*p.optimizer_data.sphere_position)[sphere_ids[i]];
			centroid += positions[i];
			average_radius += (*p.optimizer_data.sphere_radius)[sphere_ids[i]];
			region_vertex_ids.insert(vertex_ids[i]);
		}
		if (region_vertex_ids.size() != 5)
			return false;
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
				SkeletonEdge edge;
				if (!MeshOps::find_skeleton_edge_between(p, vertices[i], vertices[j], edge))
					continue;
				region_edges[make_edge_key(vertex_ids[i], vertex_ids[j])] = edge;
				if (incident_faces(p.data.skeleton, edge).size() > 2)
				{
					const uint32 edge_id = index_of(p.data.skeleton, edge);
					if (edge_id != INVALID_INDEX)
						region_nm_edge_ids.insert(edge_id);
				}
			}
		}
		if (region_nm_edge_ids.empty())
			return false;

		for (uint32 nm_edge_id : region_nm_edge_ids)
		{
			const SkeletonEdge nm_edge = of_index<SkeletonEdge>(p.data.skeleton, nm_edge_id);
			if (!nm_edge.is_valid())
				continue;
			for (SkeletonVertex endpoint : incident_vertices(p.data.skeleton, nm_edge))
			{
				for (SkeletonEdge adjacent_edge : incident_edges(p.data.skeleton, endpoint))
				{
					if (!adjacent_edge.is_valid() || incident_faces(p.data.skeleton, adjacent_edge).size() <= 2)
						continue;
					const uint32 adjacent_edge_id = index_of(p.data.skeleton, adjacent_edge);
					if (adjacent_edge_id != INVALID_INDEX &&
						region_nm_edge_ids.find(adjacent_edge_id) == region_nm_edge_ids.end())
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
		foreach_cell(p.data.skeleton, [&](SkeletonFace face) -> bool {
			if (face_is_inside_vertex_set(p, face, region_vertex_ids))
			{
				const uint32 face_id = index_of(p.data.skeleton, face);
				if (face_id != INVALID_INDEX)
					internal_face_ids.insert(face_id);
			}
			return true;
		});
		if (internal_face_ids.empty())
			return false;

		for (const auto& edge_entry : region_edges)
		{
			uint32 external_face_count = 0;
			for (SkeletonFace face : incident_faces(p.data.skeleton, edge_entry.second))
			{
				const uint32 face_id = index_of(p.data.skeleton, face);
				if (face_id == INVALID_INDEX || internal_face_ids.find(face_id) == internal_face_ids.end())
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

		RemovalStats deletion_stats;
		SkeletonFaceDeletionOptions deletion_options;
		deletion_options.remove_incident_orphan_edges_immediately = false;
		deletion_options.remove_global_orphan_elements_after_batch = false;
		if (!MeshOps::delete_skeleton_face_id_set(p, internal_face_ids, deletion_stats, deletion_options))
			return false;

		for (const auto& edge_entry : region_edges)
		{
			if (cycle_edge_keys.find(edge_entry.first) != cycle_edge_keys.end())
				continue;
			if (edge_entry.second.is_valid() && incident_faces(p.data.skeleton, edge_entry.second).empty())
				remove_edge(p.data.skeleton, edge_entry.second);
		}

		for (uint32 i = 0; i < 5; ++i)
		{
			for (uint32 j = i + 1; j < 5; ++j)
			{
				const uint64 edge_key = make_edge_key(vertex_ids[i], vertex_ids[j]);
				const bool is_cycle_edge = cycle_edge_keys.find(edge_key) != cycle_edge_keys.end();
				if (is_cycle_edge)
				{
					sorted_vector_insert((*p.optimizer_data.sphere_neighbors)[sphere_ids[i]], seed.spheres[j]);
					sorted_vector_insert((*p.optimizer_data.sphere_neighbors)[sphere_ids[j]], seed.spheres[i]);
				}
				else
				{
					sorted_vector_erase((*p.optimizer_data.sphere_neighbors)[sphere_ids[i]], seed.spheres[j]);
					sorted_vector_erase((*p.optimizer_data.sphere_neighbors)[sphere_ids[j]], seed.spheres[i]);
				}
			}
		}

		const SphereVertex center_sphere = p.data.spheres_optimizer.add_sphere(centroid, average_radius);
		if (!center_sphere.is_valid())
			return false;
		const uint32 center_sphere_id = index_of(p.optimizer_data.spheres, center_sphere);
		p.added_spheres.push_back(center_sphere);

		const SkeletonVertex center_vertex = add_vertex(p.data.skeleton);
		value<Vec3>(p.data.skeleton, p.data.skeleton_position, center_vertex) = centroid;
		value<Scalar>(p.data.skeleton, p.data.skeleton_radius, center_vertex) = average_radius;
		value<SphereVertex>(p.data.skeleton, p.data.skeleton_source_sphere, center_vertex) = center_sphere;
		(*p.data.sphere_skeleton_vertex)[center_sphere_id] = center_vertex;

		std::array<SkeletonEdge, 5> cycle_edges;
		std::array<SkeletonEdge, 5> spoke_edges;
		for (uint32 i = 0; i < 5; ++i)
		{
			const uint32 a = cyclic_order[i];
			const uint32 b = cyclic_order[(i + 1) % 5];
			if (!MeshOps::find_skeleton_edge_between(p, vertices[a], vertices[b], cycle_edges[i]))
				cycle_edges[i] = add_edge(p.data.skeleton, vertices[a], vertices[b]);
			spoke_edges[i] = add_edge(p.data.skeleton, center_vertex, vertices[a]);
			sorted_vector_insert((*p.optimizer_data.sphere_neighbors)[center_sphere_id], seed.spheres[a]);
			sorted_vector_insert((*p.optimizer_data.sphere_neighbors)[sphere_ids[a]], center_sphere);
		}

		for (uint32 i = 0; i < 5; ++i)
		{
			std::vector<SkeletonEdge> face_edges = {cycle_edges[i], spoke_edges[i], spoke_edges[(i + 1) % 5]};
			const SkeletonFace new_face = add_face(p.data.skeleton, face_edges);
			p.skeleton_tets.reset_face(new_face);
		}

		return true;
	}

	static void repair_isolated_independent_dense_five_vertex_regions(
		Context& p, const std::vector<IndependentDenseFiveVertexRepairSeed>& seeds)
	{
		std::vector<IndependentDenseFiveVertexRepairSeed> live_seeds;
		live_seeds.reserve(seeds.size());
		for (const auto& seed : seeds)
		{
			bool all_vertices_live = true;
			std::unordered_set<uint32> unique_skeleton_vertex_ids;
			unique_skeleton_vertex_ids.reserve(8);
			for (const SphereVertex& sphere : seed.spheres)
			{
				uint32 sphere_id = INVALID_INDEX;
				if (!is_live_sphere_vertex(p, sphere, sphere_id))
				{
					all_vertices_live = false;
					break;
				}
				const SkeletonVertex vertex = (*p.data.sphere_skeleton_vertex)[sphere_id];
				uint32 vertex_id = INVALID_INDEX;
				if (!is_live_skeleton_vertex(p, vertex, vertex_id) ||
					value<SphereVertex>(p.data.skeleton, p.data.skeleton_source_sphere, vertex) != sphere)
				{
					all_vertices_live = false;
					break;
				}
				unique_skeleton_vertex_ids.insert(vertex_id);
			}
			if (all_vertices_live && unique_skeleton_vertex_ids.size() == 5)
				live_seeds.push_back(seed);
		}

		bool repaired_any = false;
		for (const auto& seed : live_seeds)
			if (repair_independent_dense_five_vertex_region(p, seed))
				repaired_any = true;

		if (repaired_any)
			MeshOps::compute_edge_degree(p);
	}
};

} // namespace topology_detail

} // namespace offset_axis

} // namespace geometry

} // namespace cgogn

#endif // CGOGN_GEOMETRY_ALGOS_OFFSET_AXIS_INTERNAL_DENSE_K5_H_
