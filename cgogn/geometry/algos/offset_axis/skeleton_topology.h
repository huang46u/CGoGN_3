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
#ifndef CGOGN_GEOMETRY_ALGOS_OFFSET_AXIS_SKELETON_TOPOLOGY_H_
#define CGOGN_GEOMETRY_ALGOS_OFFSET_AXIS_SKELETON_TOPOLOGY_H_

#include <cgogn/core/utils/numerics.h>
#include <cgogn/geometry/algos/offset_axis/internal/boundary_tet_prune.h>
#include <cgogn/geometry/algos/offset_axis/internal/dense_k5.h>
#include <cgogn/geometry/algos/offset_axis/internal/greedy_tet_removal.h>
#include <cgogn/geometry/algos/offset_axis/internal/non_manifold_prune.h>
#include <cgogn/geometry/algos/offset_axis/internal/quadrature_scores.h>
#include <cgogn/geometry/algos/offset_axis/internal/skeleton_build.h>
#include <cgogn/geometry/algos/offset_axis/internal/skeleton_mesh_ops.h>
#include <cgogn/geometry/algos/offset_axis/internal/skeleton_topology_types.h>

#include <algorithm>
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

template <typename POINTS, typename NONMANIFOLD>
class SkeletonTopology
{
public:
	using Types = topology_detail::SkeletonTopologyTypes<POINTS, NONMANIFOLD>;
	using SphereVertex = typename Types::SphereVertex;
	using SkeletonVertex = typename Types::SkeletonVertex;
	using SkeletonEdge = typename Types::SkeletonEdge;
	using SkeletonFace = typename Types::SkeletonFace;
	template <typename T>
	using SphereAttribute = typename Types::template SphereAttribute<T>;
	template <typename T>
	using SkeletonAttribute = typename Types::template SkeletonAttribute<T>;
	using Data = typename Types::Data;
	using Metrics = typename Types::Metrics;
	using Status = typename Types::Status;
	using ResidualPruneMode = typename Types::ResidualPruneMode;

	explicit SkeletonTopology(Data data) : data_(std::move(data))
	{
	}
	SkeletonTopology(const SkeletonTopology&) = delete;
	SkeletonTopology& operator=(const SkeletonTopology&) = delete;
	Data& data()
	{
		return data_;
	}
	const Data& data() const
	{
		return data_;
	}
	const Metrics& metrics() const
	{
		return metrics_;
	}
	bool has_tets() const
	{
		return !skeleton_tets_.empty();
	}
	void set_residual_prune_mode(ResidualPruneMode mode)
	{
		residual_prune_mode_ = mode;
	}

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

	// Read-only: ids (index_of) of the faces of the current boundary tets, and of the longest degree-2 edge of each
	// such tet. Intended for display only; the skeleton and the tet set are not modified.
	struct BoundaryTetMask
	{
		std::vector<uint32> face_ids;
		std::vector<uint32> edge_ids;
	};

	BoundaryTetMask boundary_tet_mask()
	{
		BoundaryTetMask mask;
		if (!valid_data())
			return mask;
		Context p(data_, skeleton_tets_, metrics_, residual_prune_mode_);
		std::unordered_set<uint32> face_ids;
		std::unordered_set<uint32> edge_ids;
		BoundaryTetPrune::collect_boundary_tet_mask(p, face_ids, edge_ids);
		mask.face_ids.assign(face_ids.begin(), face_ids.end());
		mask.edge_ids.assign(edge_ids.begin(), edge_ids.end());
		std::sort(mask.face_ids.begin(), mask.face_ids.end());
		std::sort(mask.edge_ids.begin(), mask.edge_ids.end());
		return mask;
	}

	// Read-only: UDF score of each face adjacent to a current tet, keyed by index_of(face). Same integral as the
	// topology-fix cache (QuadratureScores::compute_skeleton_face_scores), optionally divided by the face area.
	// The skeleton and the tet set are not modified. Returns false if the scores cannot be evaluated.
	template <typename ScoreEvaluator>
	bool tet_face_udf_scores(ScoreEvaluator& evaluate_scores, bool normalize_by_area,
							 std::unordered_map<uint32, Scalar>& face_scores)
	{
		face_scores.clear();
		if (!valid_data())
			return false;
		Context p(data_, skeleton_tets_, metrics_, residual_prune_mode_);
		const std::unordered_set<uint32> face_ids = QuadratureScores::collect_current_tet_face_id_whitelist(p);
		return QuadratureScores::compute_skeleton_face_scores(p, evaluate_scores, face_ids, face_scores,
															  normalize_by_area);
	}

	// Read-only: 3-point Gauss UDF average along each edge adjacent to a current tet, keyed by index_of(edge). Same
	// values as the topology-fix cache (QuadratureScores::compute_skeleton_edge_scores_gauss3_normalized).
	template <typename ScoreEvaluator>
	bool tet_edge_gauss3_scores(ScoreEvaluator& evaluate_scores, std::unordered_map<uint32, Scalar>& edge_scores)
	{
		edge_scores.clear();
		if (!valid_data())
			return false;
		Context p(data_, skeleton_tets_, metrics_, residual_prune_mode_);
		const std::unordered_set<uint32> edge_ids = QuadratureScores::collect_current_tet_edge_ids(p);
		return QuadratureScores::compute_skeleton_edge_scores_gauss3_normalized(p, evaluate_scores, edge_ids,
																				edge_scores);
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
	using TetMap = typename Types::TetMap;
	using Context = typename Types::Context;
	using RemovalStats = typename Types::RemovalStats;
	using TopologyFixScoreCache = typename Types::TopologyFixScoreCache;
	using MeshOps = topology_detail::SkeletonMeshOps<POINTS, NONMANIFOLD>;
	using SkeletonBuild = topology_detail::SkeletonBuild<POINTS, NONMANIFOLD>;
	using QuadratureScores = topology_detail::QuadratureScores<POINTS, NONMANIFOLD>;
	using NonManifoldPrune = topology_detail::NonManifoldPrune<POINTS, NONMANIFOLD>;
	using GreedyTetRemoval = topology_detail::GreedyTetRemoval<POINTS, NONMANIFOLD>;
	using BoundaryTetPrune = topology_detail::BoundaryTetPrune<POINTS, NONMANIFOLD>;
	using DenseK5Repair = topology_detail::DenseK5Repair<POINTS, NONMANIFOLD>;
	using DenseFiveVertexDetectionResult = typename DenseK5Repair::DenseFiveVertexDetectionResult;
	using IndependentDenseFiveVertexRepairSeed = typename DenseK5Repair::IndependentDenseFiveVertexRepairSeed;
	using IndependentDenseFiveVertexRepairStats = typename DenseK5Repair::IndependentDenseFiveVertexRepairStats;

	bool valid_data() const
	{
		return data_.skeleton && data_.spheres_optimizer && data_.skeleton_position && data_.skeleton_radius &&
			   data_.skeleton_source_sphere && data_.face_incident_tets && data_.edge_degree &&
			   data_.sphere_skeleton_vertex;
	}

	//------------------------------//
	//-----Topology correction------//
	//------------------------------//

	template <typename ScoreEvaluator>
	Status run_topology_fix_pipeline(Context& p, ScoreEvaluator& evaluator)
	{
		if (!p.skeleton_ || !p.incident_tets_)
			return Status::invalid_data;
		const DenseFiveVertexDetectionResult dense_five_vertex_detection =
			DenseK5Repair::detect_k5_and_k5_minus_1_cells(p);
		const std::vector<IndependentDenseFiveVertexRepairSeed> independent_dense_repair_seeds =
			DenseK5Repair::collect_independent_dense_five_vertex_repair_seeds(p, dense_five_vertex_detection);
		TopologyFixScoreCache score_cache;
		if (!QuadratureScores::initialize_topology_fix_score_cache(p, evaluator, score_cache))
			return Status::score_evaluation_failed;
		while (!p.skeleton_tets_.empty())
		{
			const uint32 round_start_tets = static_cast<uint32>(p.skeleton_tets_.size());
			const RemovalStats boundary_stats = BoundaryTetPrune::run_boundary_tet_face_deletion(p, &score_cache);
			RemovalStats simple_stats;
			if (!p.skeleton_tets_.empty())
				simple_stats = GreedyTetRemoval::template run_edge_score_tet_mode_topology_fix<
					topology_detail::SimpleTetPolicy>(p, score_cache, evaluator);
			RemovalStats nonsimple_stats;
			if (!p.skeleton_tets_.empty())
				nonsimple_stats = GreedyTetRemoval::template run_edge_score_tet_mode_topology_fix<
					topology_detail::NonSimpleTetPolicy>(p, score_cache, evaluator);
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
		const IndependentDenseFiveVertexRepairStats repair_stats =
			DenseK5Repair::repair_isolated_independent_dense_five_vertex_regions(p, independent_dense_repair_seeds);
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

} // namespace offset_axis

} // namespace geometry
} // namespace cgogn

#endif // CGOGN_GEOMETRY_ALGOS_OFFSET_AXIS_SKELETON_TOPOLOGY_H_
