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
#include <cgogn/geometry/algos/udf/internal/boundary_tet_prune.h>
#include <cgogn/geometry/algos/udf/internal/dense_k5.h>
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
	using BoundaryTetPrune = skeleton_topology_detail::BoundaryTetPrune<POINTS, NONMANIFOLD>;
	using DenseK5Repair = skeleton_topology_detail::DenseK5Repair<POINTS, NONMANIFOLD>;
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
		const DenseFiveVertexDetectionResult dense_five_vertex_detection = DenseK5Repair::detect_k5_and_k5_minus_1_cells(p);
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
		const IndependentDenseFiveVertexRepairStats repair_stats = DenseK5Repair::repair_isolated_independent_dense_five_vertex_regions(p, independent_dense_repair_seeds);
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
