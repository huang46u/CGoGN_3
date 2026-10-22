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
#ifndef CGOGN_GEOMETRY_ALGOS_UDF_SKELETON_TOPOLOGY_H_
#define CGOGN_GEOMETRY_ALGOS_UDF_SKELETON_TOPOLOGY_H_

#include <cgogn/core/functions/attributes.h>
#include <cgogn/core/functions/mesh_info.h>
#include <cgogn/core/functions/traversals/edge.h>
#include <cgogn/core/functions/traversals/face.h>
#include <cgogn/core/functions/traversals/vertex.h>
#include <cgogn/core/types/cell_marker.h>
#include <cgogn/core/utils/numerics.h>
#include <cgogn/geometry/algos/area.h>
#include <cgogn/geometry/algos/udf/internal/skeleton_topology_types.h>
#include <cgogn/geometry/algos/udf/internal/greedy_tet_removal.h>
#include <cgogn/geometry/algos/udf/internal/non_manifold_prune.h>
#include <cgogn/geometry/algos/udf/internal/quadrature_scores.h>
#include <cgogn/geometry/algos/udf/internal/skeleton_build.h>
#include <cgogn/geometry/algos/udf/internal/skeleton_mesh_ops.h>
#include <cgogn/geometry/algos/udf/spheres_optimizer.h>
#include <cgogn/geometry/types/vector_traits.h>

#include <Eigen/Dense>
#include <Eigen/Eigenvalues>
#include <libacc/kd_tree.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <limits>
#include <map>
#include <memory>
#include <queue>
#include <set>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

namespace cgogn
{
namespace geometry
{

template <typename POINTS, typename NONMANIFOLD>
class SkeletonTopology
{
public:
	using Types = skeleton_topology_detail::SkeletonTopologyTypes<POINTS, NONMANIFOLD>;
	using SphereVertex = typename Types::SphereVertex;
	using SkeletonVertex = typename Types::SkeletonVertex;
	using SkeletonEdge = typename Types::SkeletonEdge;
	using SkeletonFace = typename Types::SkeletonFace;
	template <typename T> using SphereAttribute = typename Types::template SphereAttribute<T>;
	template <typename T> using SkeletonAttribute = typename Types::template SkeletonAttribute<T>;
	using Data = typename Types::Data;
	using Metrics = typename Types::Metrics;
	using Status = typename Types::Status;
	using ResidualPruneMode = typename Types::ResidualPruneMode;

	explicit SkeletonTopology(Data data) : data_(std::move(data)) {}
	SkeletonTopology(const SkeletonTopology&) = delete;
	SkeletonTopology& operator=(const SkeletonTopology&) = delete;
	Data& data() { return data_; }
	const Data& data() const { return data_; }
	const Metrics& metrics() const { return metrics_; }
	void set_residual_prune_mode(ResidualPruneMode mode) { residual_prune_mode_ = mode; }

	Status build_from_spheres()
	{
		if (!valid_data())
			return Status::invalid_data;
		metrics_ = Metrics{};
		Context p(data_, skeleton_tets_, metrics_, residual_prune_mode_);
		if (!SkeletonBuild::compute_skeleton(p))
			return Status::invalid_data;
		metrics_.tet_count = static_cast<uint32>(skeleton_tets_.size());
		return Status::success;
	}

	template <typename ScoreEvaluator>
	Status fix_topology(ScoreEvaluator&& evaluate_scores)
	{
		if (!valid_data())
			return Status::invalid_data;
		metrics_ = Metrics{};
		if (skeleton_tets_.empty())
			return Status::success;
		Context p(data_, skeleton_tets_, metrics_, residual_prune_mode_);
		return run_topology_fix_pipeline(p, evaluate_scores);
	}

	Status prune_fully_non_manifold_triangles()
	{
		if (!valid_data())
			return Status::invalid_data;
		Context p(data_, skeleton_tets_, metrics_, residual_prune_mode_);
		return NonManifoldPrune::prune_fully_non_manifold_triangles(p) ? Status::success : Status::invalid_data;
	}

	void compute_edge_degree()
	{
		if (!valid_data())
			return;
		Context p(data_, skeleton_tets_, metrics_, residual_prune_mode_);
		MeshOps::compute_edge_degree(p);
	}

	Status prune_residual_sheets()
	{
		if (!valid_data())
			return Status::invalid_data;
		metrics_ = Metrics{};
		Context p(data_, skeleton_tets_, metrics_, residual_prune_mode_);
		NonManifoldPrune::run_residual_sheet_prune(p);
		metrics_.tet_count = static_cast<uint32>(skeleton_tets_.size());
		return Status::success;
	}

private:
	using SpheresOptimizerType = typename Types::SpheresOptimizerType;
	using FaceKey = typename Types::FaceKey;
	using Tet = typename Types::Tet;
	using TetMap = typename Types::TetMap;
	using Context = typename Types::Context;
	using edge_hash = typename Types::edge_hash;
	using edge_equal = typename Types::edge_equal;
	using RemovalStats = typename Types::RemovalStats;
	using SkeletonFaceDeletionOptions = typename Types::SkeletonFaceDeletionOptions;
	using TopologyFixScoreCache = typename Types::TopologyFixScoreCache;
	using EdgeTetDeleteMode = typename Types::EdgeTetDeleteMode;
	using MeshOps = skeleton_topology_detail::SkeletonMeshOps<POINTS, NONMANIFOLD>;
	using SkeletonBuild = skeleton_topology_detail::SkeletonBuild<POINTS, NONMANIFOLD>;
	using QuadratureScores = skeleton_topology_detail::QuadratureScores<POINTS, NONMANIFOLD>;
	using NonManifoldPrune = skeleton_topology_detail::NonManifoldPrune<POINTS, NONMANIFOLD>;
	using GreedyTetRemoval = skeleton_topology_detail::GreedyTetRemoval<POINTS, NONMANIFOLD>;

	bool valid_data() const
	{
		return data_.skeleton && data_.spheres_optimizer && data_.skeleton_position && data_.skeleton_radius &&
			data_.skeleton_source_sphere && data_.face_incident_tets && data_.edge_degree &&
			data_.sphere_skeleton_vertex;
	}



	//------------------------------//
	//-----Topology correction------//
	//------------------------------//



	struct DenseFiveVertexDetectionResult
	{
		std::vector<std::array<uint32, 5>> independent_exact_k5_regions;
		std::vector<std::array<uint32, 5>> independent_k5_minus_1_regions;
	};

	DenseFiveVertexDetectionResult detect_k5_and_k5_minus_1_cells(Context& p)
	{
		DenseFiveVertexDetectionResult result;
		if (!p.skeleton_)
			return result;

		auto make_edge_key = [](uint32 a, uint32 b) -> uint64 {
			if (a > b)
				std::swap(a, b);
			return (uint64(a) << 32) | uint64(b);
		};
		std::vector<std::array<uint32, 5>> exact_k5_regions;
		std::vector<std::array<uint32, 5>> k5_minus_1_regions;

		std::unordered_map<uint32, std::unordered_set<uint32>> adjacency;
		adjacency.reserve(nb_cells<SkeletonVertex>(*p.skeleton_));
		std::unordered_set<uint64> edge_pair_to_id;
		edge_pair_to_id.reserve(nb_cells<SkeletonEdge>(*p.skeleton_) * 2 + 1);

		foreach_cell(*p.skeleton_, [&](SkeletonEdge e) -> bool {
			const uint32 ide = index_of(*p.skeleton_, e);
			if (ide == INVALID_INDEX)
				return true;
			const std::vector<SkeletonVertex> vv = incident_vertices(*p.skeleton_, e);
			if (vv.size() != 2)
				return true;
			const uint32 a = index_of(*p.skeleton_, vv[0]);
			const uint32 b = index_of(*p.skeleton_, vv[1]);
			if (a == INVALID_INDEX || b == INVALID_INDEX || a == b)
				return true;
			adjacency[a].insert(b);
			adjacency[b].insert(a);
			edge_pair_to_id.insert(make_edge_key(a, b));
			return true;
		});

		std::set<std::array<uint32, 5>> seen_region_keys;
		for (const auto& tet_entry : p.skeleton_tets_)
		{
			const Tet& tet = tet_entry.second;
			std::array<uint32, 4> base_vertices = {INVALID_INDEX, INVALID_INDEX, INVALID_INDEX, INVALID_INDEX};
			uint32 base_vertex_count = 0;
			for (uint32 face_index = 0; face_index < 4; ++face_index)
			{
				const SkeletonFace face = tet.faces[face_index];
				if (!face.is_valid())
					continue;
				for (SkeletonVertex vertex : incident_vertices(*p.skeleton_, face))
				{
					const uint32 vertex_id = index_of(*p.skeleton_, vertex);
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

				std::array<uint32, 5> region = {
					base_vertices[0], base_vertices[1], base_vertices[2], base_vertices[3], candidate.first};
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

	RemovalStats run_boundary_tet_face_deletion(
		Context& p, TopologyFixScoreCache* score_cache = nullptr)
	{
		RemovalStats stats;

		struct BoundaryTetCandidate
		{
			std::size_t tet_id = std::size_t(-1);
			SkeletonEdge edge;
			SkeletonFace f0;
			SkeletonFace f1;
		};

		std::vector<std::size_t> tet_ids;
		tet_ids.reserve(p.skeleton_tets_.size());
		for (const auto& kv : p.skeleton_tets_)
			tet_ids.push_back(kv.first);

		std::vector<BoundaryTetCandidate> candidates;
		candidates.reserve(tet_ids.size());

		// Phase 1: scan all current tets once and collect boundary-tet deletion candidates.
		for (std::size_t tet_id : tet_ids)
		{
			auto it_tet = p.skeleton_tets_.find(tet_id);
			if (it_tet == p.skeleton_tets_.end())
				continue;
			const Tet tet = it_tet->second;

			bool is_boundary_tet = false;
			for (uint32 i = 0; i < 4; ++i)
			{
				const SkeletonFace f = tet.faces[i];
				uint32 deg2_edge_count = 0;
				for (SkeletonEdge e : incident_edges(*p.skeleton_, f))
				{
					if (incident_faces(*p.skeleton_, e).size() == 2)
						++deg2_edge_count;
				}
				if (deg2_edge_count >= 2)
				{
					is_boundary_tet = true;
					break;
				}
			}
			if (!is_boundary_tet)
				continue;

			struct BoundaryEdgeChoice
			{
				SkeletonEdge edge;
				SkeletonFace f0;
				SkeletonFace f1;
				Scalar length = Scalar(-1);
			};
			BoundaryEdgeChoice best_choice;
			std::unordered_set<uint32> seen_edges;
			seen_edges.reserve(8);

			for (uint32 i = 0; i < 4; ++i)
			{
				const SkeletonFace f = tet.faces[i];
				for (SkeletonEdge e : incident_edges(*p.skeleton_, f))
				{
					const uint32 ide = index_of(*p.skeleton_, e);
					if (!seen_edges.insert(ide).second)
						continue;

					const auto in_faces = incident_faces(*p.skeleton_, e);
					if (in_faces.size() != 2)
						continue;
					const SkeletonFace f0 = in_faces[0];
					const SkeletonFace f1 = in_faces[1];
					const uint32 idf0 = index_of(*p.skeleton_, f0);
					const uint32 idf1 = index_of(*p.skeleton_, f1);
					const std::vector<SkeletonVertex> edge_vertices = incident_vertices(*p.skeleton_, e);
					if (edge_vertices.size() != 2)
						continue;
					const Vec3 p0 = value<Vec3>(*p.skeleton_, p.skeleton_position_, edge_vertices[0]);
					const Vec3 p1 = value<Vec3>(*p.skeleton_, p.skeleton_position_, edge_vertices[1]);
					const Scalar edge_len = (p1 - p0).norm();
					if (edge_len > best_choice.length)
					{
						best_choice.edge = e;
						best_choice.f0 = f0;
						best_choice.f1 = f1;
						best_choice.length = edge_len;
					}
				}
			}

			if (best_choice.length >= Scalar(0))
				candidates.push_back({tet_id, best_choice.edge, best_choice.f0, best_choice.f1});
		}

		// Phase 2: delete using precomputed candidates only (no boundary re-check during deletion).
		for (const BoundaryTetCandidate& cand : candidates)
		{
			auto it_remove_tet = p.skeleton_tets_.find(cand.tet_id);
			if (it_remove_tet == p.skeleton_tets_.end())
				continue;
			const Tet old_tet = it_remove_tet->second;

			const std::array<SkeletonFace, 2> faces_to_remove = {cand.f0, cand.f1};
			std::vector<SkeletonEdge> affected_edges;
			affected_edges.reserve(6);
			for (const SkeletonFace& f_del : faces_to_remove)
			{
				if (!f_del.is_valid())
					continue;
				for (SkeletonEdge e : incident_edges(*p.skeleton_, f_del))
					affected_edges.push_back(e);
			}
			std::unordered_set<uint32> removed_face_ids;
			removed_face_ids.reserve(faces_to_remove.size() * 2 + 1);
			for (const SkeletonFace& f_del : faces_to_remove)
			{
				if (const uint32 idf_del = index_of(*p.skeleton_, f_del); idf_del != INVALID_INDEX)
					removed_face_ids.insert(idf_del);
				remove_face(*p.skeleton_, f_del);
				++stats.removed_faces;
			}
			QuadratureScores::erase_topology_fix_scores_for_removed_faces_and_orphan_edges(p, score_cache, removed_face_ids, affected_edges);
			if (cand.edge.is_valid())
			{
				const auto edge_in_faces = incident_faces(*p.skeleton_, cand.edge);
				if (edge_in_faces.empty())
				{
					remove_edge(*p.skeleton_, cand.edge);
					++stats.removed_edges;
				}
			}
			p.skeleton_tets_.erase(it_remove_tet);
			++stats.removed_tets;
			for (uint32 i = 0; i < 4; ++i)
			{
				const SkeletonFace f = old_tet.faces[i];
				const uint32 idf = index_of(*p.skeleton_, f);
				(*p.incident_tets_)[idf].erase(cand.tet_id);
			}
		}

		return stats;
	}



	struct IndependentDenseFiveVertexRepairSeed
	{
		std::array<SphereVertex, 5> spheres;
	};

	struct IndependentDenseFiveVertexRepairStats
	{
		uint32 repaired_regions = 0;
		uint32 added_faces = 0;
	};

	bool is_live_sphere_vertex(Context& p, const SphereVertex& sphere, uint32& sphere_index)
	{
		sphere_index = INVALID_INDEX;
		if (!p.spheres_ || !sphere.is_valid() || sphere.dart_.index_ >= p.spheres_->darts_.maximum_index())
			return false;
		sphere_index = index_of(*p.spheres_, sphere);
		return sphere_index != INVALID_INDEX;
	}

	bool is_live_skeleton_vertex(Context& p, const SkeletonVertex& vertex, uint32& vertex_index)
	{
		vertex_index = INVALID_INDEX;
		if (!p.skeleton_ || !vertex.is_valid())
			return false;
		vertex_index = index_of(*p.skeleton_, vertex);
		return vertex_index != INVALID_INDEX;
	}

	std::vector<IndependentDenseFiveVertexRepairSeed> collect_independent_dense_five_vertex_repair_seeds(
		Context& p, const DenseFiveVertexDetectionResult& detection)
	{
		std::vector<IndependentDenseFiveVertexRepairSeed> seeds;
		seeds.reserve(detection.independent_exact_k5_regions.size() +
					  detection.independent_k5_minus_1_regions.size());
		if (!p.skeleton_ || !p.spheres_ || !p.skeleton_source_sphere_)
			return seeds;

		auto append_regions = [&](const std::vector<std::array<uint32, 5>>& regions) {
			for (const auto& region : regions)
			{
				IndependentDenseFiveVertexRepairSeed seed;
				bool valid = true;
				for (uint32 i = 0; i < 5; ++i)
				{
					const SkeletonVertex vertex = of_index<SkeletonVertex>(*p.skeleton_, region[i]);
					uint32 current_vertex_index = INVALID_INDEX;
					if (!is_live_skeleton_vertex(p, vertex, current_vertex_index) || current_vertex_index != region[i])
					{
						valid = false;
						break;
					}
					const SphereVertex sphere = value<SphereVertex>(*p.skeleton_, p.skeleton_source_sphere_, vertex);
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

	bool face_is_inside_vertex_set(Context& p, const SkeletonFace& face,
								   const std::unordered_set<uint32>& vertex_ids)
	{
		if (!p.skeleton_ || !face.is_valid())
			return false;
		const std::vector<SkeletonVertex> vertices = incident_vertices(*p.skeleton_, face);
		if (vertices.size() != 3)
			return false;
		for (SkeletonVertex vertex : vertices)
		{
			const uint32 vertex_id = index_of(*p.skeleton_, vertex);
			if (vertex_id == INVALID_INDEX || vertex_ids.find(vertex_id) == vertex_ids.end())
				return false;
		}
		return true;
	}

	enum class IndependentDenseFiveVertexRepairResult
	{
		Repaired,
		MissingVertices,
		NoNonManifoldEdges,
		AdjacentNonManifoldEdges,
		UnsafeExternalFaces,
		DegenerateGeometry
	};

	IndependentDenseFiveVertexRepairResult repair_independent_dense_five_vertex_region(
		Context& p, const IndependentDenseFiveVertexRepairSeed& seed, uint32& out_added_faces)
	{
		out_added_faces = 0;
		if (!p.spheres_ || !p.spheres_position_ || !p.spheres_radius_ || !p.spheres_neighbor_clusters_ ||
			!p.spheres_skeleton_vertex_ || !p.skeleton_ || !p.skeleton_position_ || !p.skeleton_radius_ ||
			!p.skeleton_source_sphere_ || !p.incident_tets_)
			return IndependentDenseFiveVertexRepairResult::MissingVertices;

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
				return IndependentDenseFiveVertexRepairResult::MissingVertices;
			vertices[i] = (*p.spheres_skeleton_vertex_)[sphere_ids[i]];
			if (!is_live_skeleton_vertex(p, vertices[i], vertex_ids[i]))
				return IndependentDenseFiveVertexRepairResult::MissingVertices;
			if (value<SphereVertex>(*p.skeleton_, p.skeleton_source_sphere_, vertices[i]) != seed.spheres[i])
				return IndependentDenseFiveVertexRepairResult::MissingVertices;
			positions[i] = (*p.spheres_position_)[sphere_ids[i]];
			centroid += positions[i];
			average_radius += (*p.spheres_radius_)[sphere_ids[i]];
			region_vertex_ids.insert(vertex_ids[i]);
		}
		if (region_vertex_ids.size() != 5)
			return IndependentDenseFiveVertexRepairResult::MissingVertices;
		centroid /= Scalar(5);
		average_radius /= Scalar(5);
		if (!centroid.allFinite() || !std::isfinite(static_cast<double>(average_radius)) || average_radius <= Scalar(0))
			return IndependentDenseFiveVertexRepairResult::DegenerateGeometry;

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
				if (incident_faces(*p.skeleton_, edge).size() > 2)
				{
					const uint32 edge_id = index_of(*p.skeleton_, edge);
					if (edge_id != INVALID_INDEX)
						region_nm_edge_ids.insert(edge_id);
				}
			}
		}
		if (region_nm_edge_ids.empty())
			return IndependentDenseFiveVertexRepairResult::NoNonManifoldEdges;

		for (uint32 nm_edge_id : region_nm_edge_ids)
		{
			const SkeletonEdge nm_edge = of_index<SkeletonEdge>(*p.skeleton_, nm_edge_id);
			if (!nm_edge.is_valid())
				continue;
			for (SkeletonVertex endpoint : incident_vertices(*p.skeleton_, nm_edge))
			{
				for (SkeletonEdge adjacent_edge : incident_edges(*p.skeleton_, endpoint))
				{
					if (!adjacent_edge.is_valid() || incident_faces(*p.skeleton_, adjacent_edge).size() <= 2)
						continue;
					const uint32 adjacent_edge_id = index_of(*p.skeleton_, adjacent_edge);
					if (adjacent_edge_id != INVALID_INDEX &&
						region_nm_edge_ids.find(adjacent_edge_id) == region_nm_edge_ids.end())
						return IndependentDenseFiveVertexRepairResult::AdjacentNonManifoldEdges;
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
			return IndependentDenseFiveVertexRepairResult::DegenerateGeometry;
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
				return IndependentDenseFiveVertexRepairResult::DegenerateGeometry;
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
		foreach_cell(*p.skeleton_, [&](SkeletonFace face) -> bool {
			if (face_is_inside_vertex_set(p, face, region_vertex_ids))
			{
				const uint32 face_id = index_of(*p.skeleton_, face);
				if (face_id != INVALID_INDEX)
					internal_face_ids.insert(face_id);
			}
			return true;
		});
		if (internal_face_ids.empty())
			return IndependentDenseFiveVertexRepairResult::UnsafeExternalFaces;

		for (const auto& edge_entry : region_edges)
		{
			uint32 external_face_count = 0;
			for (SkeletonFace face : incident_faces(*p.skeleton_, edge_entry.second))
			{
				const uint32 face_id = index_of(*p.skeleton_, face);
				if (face_id == INVALID_INDEX || internal_face_ids.find(face_id) == internal_face_ids.end())
					++external_face_count;
			}
			if (cycle_edge_keys.find(edge_entry.first) == cycle_edge_keys.end())
			{
				if (external_face_count != 0)
					return IndependentDenseFiveVertexRepairResult::UnsafeExternalFaces;
			}
			else if (external_face_count > 1)
			{
				return IndependentDenseFiveVertexRepairResult::UnsafeExternalFaces;
			}
		}

		RemovalStats deletion_stats;
		SkeletonFaceDeletionOptions deletion_options;
		deletion_options.remove_incident_orphan_edges_immediately = false;
		deletion_options.remove_global_orphan_elements_after_batch = false;
		if (!MeshOps::delete_skeleton_face_id_set(p, internal_face_ids, deletion_stats, deletion_options))
			return IndependentDenseFiveVertexRepairResult::UnsafeExternalFaces;

		for (const auto& edge_entry : region_edges)
		{
			if (cycle_edge_keys.find(edge_entry.first) != cycle_edge_keys.end())
				continue;
			if (edge_entry.second.is_valid() && incident_faces(*p.skeleton_, edge_entry.second).empty())
				remove_edge(*p.skeleton_, edge_entry.second);
		}

		for (uint32 i = 0; i < 5; ++i)
		{
			for (uint32 j = i + 1; j < 5; ++j)
			{
				const uint64 edge_key = make_edge_key(vertex_ids[i], vertex_ids[j]);
				const bool is_cycle_edge = cycle_edge_keys.find(edge_key) != cycle_edge_keys.end();
				if (is_cycle_edge)
				{
					sorted_vector_insert((*p.spheres_neighbor_clusters_)[sphere_ids[i]], seed.spheres[j]);
					sorted_vector_insert((*p.spheres_neighbor_clusters_)[sphere_ids[j]], seed.spheres[i]);
				}
				else
				{
					sorted_vector_erase((*p.spheres_neighbor_clusters_)[sphere_ids[i]], seed.spheres[j]);
					sorted_vector_erase((*p.spheres_neighbor_clusters_)[sphere_ids[j]], seed.spheres[i]);
				}
			}
		}

		const SphereVertex center_sphere = p.spheres_optimizer_->add_sphere(centroid, average_radius);
		if (!center_sphere.is_valid())
			return IndependentDenseFiveVertexRepairResult::MissingVertices;
		const uint32 center_sphere_id = index_of(*p.spheres_, center_sphere);
		p.metrics.added_spheres.push_back(center_sphere);

		const SkeletonVertex center_vertex = add_vertex(*p.skeleton_);
		value<Vec3>(*p.skeleton_, p.skeleton_position_, center_vertex) = centroid;
		value<Scalar>(*p.skeleton_, p.skeleton_radius_, center_vertex) = average_radius;
		value<SphereVertex>(*p.skeleton_, p.skeleton_source_sphere_, center_vertex) = center_sphere;
		(*p.spheres_skeleton_vertex_)[center_sphere_id] = center_vertex;

		std::array<SkeletonEdge, 5> cycle_edges;
		std::array<SkeletonEdge, 5> spoke_edges;
		for (uint32 i = 0; i < 5; ++i)
		{
			const uint32 a = cyclic_order[i];
			const uint32 b = cyclic_order[(i + 1) % 5];
			if (!MeshOps::find_skeleton_edge_between(p, vertices[a], vertices[b], cycle_edges[i]))
				cycle_edges[i] = add_edge(*p.skeleton_, vertices[a], vertices[b]);
			spoke_edges[i] = add_edge(*p.skeleton_, center_vertex, vertices[a]);
			sorted_vector_insert((*p.spheres_neighbor_clusters_)[center_sphere_id], seed.spheres[a]);
			sorted_vector_insert((*p.spheres_neighbor_clusters_)[sphere_ids[a]], center_sphere);
		}

		for (uint32 i = 0; i < 5; ++i)
		{
			std::vector<SkeletonEdge> face_edges = {cycle_edges[i], spoke_edges[i], spoke_edges[(i + 1) % 5]};
			const SkeletonFace new_face = add_face(*p.skeleton_, face_edges);
			value<std::set<std::size_t>>(*p.skeleton_, p.incident_tets_, new_face).clear();
			++out_added_faces;
		}

		return IndependentDenseFiveVertexRepairResult::Repaired;
	}

	IndependentDenseFiveVertexRepairStats repair_isolated_independent_dense_five_vertex_regions(
		Context& p, const std::vector<IndependentDenseFiveVertexRepairSeed>& seeds)
	{
		IndependentDenseFiveVertexRepairStats stats;
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
				const SkeletonVertex vertex = (*p.spheres_skeleton_vertex_)[sphere_id];
				uint32 vertex_id = INVALID_INDEX;
				if (!is_live_skeleton_vertex(p, vertex, vertex_id) ||
					value<SphereVertex>(*p.skeleton_, p.skeleton_source_sphere_, vertex) != sphere)
				{
					all_vertices_live = false;
					break;
				}
				unique_skeleton_vertex_ids.insert(vertex_id);
			}
			if (all_vertices_live && unique_skeleton_vertex_ids.size() == 5)
				live_seeds.push_back(seed);
		}

		for (const auto& seed : live_seeds)
		{
			uint32 added_faces = 0;
			const IndependentDenseFiveVertexRepairResult result =
				repair_independent_dense_five_vertex_region(p, seed, added_faces);
			switch (result)
			{
			case IndependentDenseFiveVertexRepairResult::Repaired:
				++stats.repaired_regions;
				stats.added_faces += added_faces;
				break;
			case IndependentDenseFiveVertexRepairResult::NoNonManifoldEdges:
			case IndependentDenseFiveVertexRepairResult::MissingVertices:
			case IndependentDenseFiveVertexRepairResult::AdjacentNonManifoldEdges:
			case IndependentDenseFiveVertexRepairResult::UnsafeExternalFaces:
			case IndependentDenseFiveVertexRepairResult::DegenerateGeometry:
				break;
			}
		}

		if (stats.repaired_regions != 0)
			MeshOps::compute_edge_degree(p);
		return stats;
	}
	template <typename ScoreEvaluator>
	Status run_topology_fix_pipeline(Context& p, ScoreEvaluator& evaluator)
	{
		if (!p.skeleton_ || !p.incident_tets_)
			return Status::invalid_data;
		const DenseFiveVertexDetectionResult dense_five_vertex_detection = detect_k5_and_k5_minus_1_cells(p);
		const std::vector<IndependentDenseFiveVertexRepairSeed> independent_dense_repair_seeds =
			collect_independent_dense_five_vertex_repair_seeds(p, dense_five_vertex_detection);
		TopologyFixScoreCache score_cache;
		if (!QuadratureScores::initialize_topology_fix_score_cache(p, evaluator, score_cache))
			return Status::score_evaluation_failed;
		while (!p.skeleton_tets_.empty())
		{
			const uint32 round_start_tets = static_cast<uint32>(p.skeleton_tets_.size());
			const RemovalStats boundary_stats = run_boundary_tet_face_deletion(p, &score_cache);
			RemovalStats simple_stats;
			if (!p.skeleton_tets_.empty())
				simple_stats = GreedyTetRemoval::run_edge_score_tet_mode_topology_fix(
					p, EdgeTetDeleteMode::SimpleTet, score_cache, evaluator);
			RemovalStats nonsimple_stats;
			if (!p.skeleton_tets_.empty())
				nonsimple_stats = GreedyTetRemoval::run_edge_score_tet_mode_topology_fix(
					p, EdgeTetDeleteMode::NonSimpleTet, score_cache, evaluator);
			const uint32 round_end_tets = static_cast<uint32>(p.skeleton_tets_.size());
			const uint32 round_removed_tets = round_start_tets - round_end_tets;
			++p.metrics.topology_rounds;
			p.metrics.boundary_removed_tets += boundary_stats.removed_tets;
			p.metrics.simple_removed_tets += simple_stats.removed_tets;
			p.metrics.nonsimple_removed_tets += nonsimple_stats.removed_tets;
			RemovalStats round_stats = boundary_stats;
			round_stats += simple_stats;
			round_stats += nonsimple_stats;
			p.metrics.removed_faces += round_stats.removed_faces;
			p.metrics.removed_edges += round_stats.removed_edges;
			p.metrics.removed_tets += round_removed_tets;
			if (round_removed_tets == 0)
				break;
		}
		const IndependentDenseFiveVertexRepairStats repair_stats = repair_isolated_independent_dense_five_vertex_regions(p, independent_dense_repair_seeds);
		p.metrics.repaired_dense_regions += repair_stats.repaired_regions;
		p.metrics.added_faces += repair_stats.added_faces;
		p.metrics.tet_count = static_cast<uint32>(p.skeleton_tets_.size());
		return Status::success;
	}

	Data data_;
	TetMap skeleton_tets_;
	Metrics metrics_;
	ResidualPruneMode residual_prune_mode_ = ResidualPruneMode::collapsible_only;
};

} // namespace geometry
} // namespace cgogn

#endif // CGOGN_GEOMETRY_ALGOS_UDF_SKELETON_TOPOLOGY_H_
