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
#include <cgogn/geometry/algos/offset_axis/internal/skeleton_build.h>
#include <cgogn/geometry/algos/offset_axis/internal/skeleton_mesh_ops.h>
#include <cgogn/geometry/algos/offset_axis/internal/skeleton_topology_types.h>
#include <cgogn/geometry/algos/offset_axis/internal/topology_scores.h>

#include <algorithm>
#include <limits>
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
	using ResidualPruneMode = typename Types::ResidualPruneMode;

	explicit SkeletonTopology(Data data)
		: data_(std::move(data)), skeleton_tets_(data_.skeleton, data_.face_incident_tets)
	{
	}
	SkeletonTopology(const SkeletonTopology&) = delete;
	SkeletonTopology& operator=(const SkeletonTopology&) = delete;
	// Spheres added by the dense K5 repair of the last fix_topology.
	const std::vector<SphereVertex>& added_spheres() const
	{
		return added_spheres_;
	}
	bool has_tets() const
	{
		return !skeleton_tets_.empty();
	}
	void build_from_spheres()
	{
		added_spheres_.clear();
		const auto& spheres = data_.spheres_optimizer.data();
		build_skeleton_from_spheres(spheres.spheres, spheres.sphere_position, spheres.sphere_radius,
									spheres.sphere_neighbors, data_.skeleton, data_.skeleton_position,
									data_.skeleton_radius, data_.skeleton_source_sphere, data_.sphere_skeleton_vertex,
									data_.edge_degree, skeleton_tets_);
	}

	// Returns false if the scores cannot be evaluated.
	template <typename ScoreEvaluator>
	bool fix_topology(ScoreEvaluator&& evaluate_scores)
	{
		added_spheres_.clear();
		if (skeleton_tets_.empty())
			return true;
		Context p(data_, skeleton_tets_, added_spheres_);
		return run_topology_fix_pipeline(p, evaluate_scores);
	}

	void compute_edge_degree()
	{
		offset_axis::compute_edge_degree(data_.skeleton, data_.edge_degree);
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
		Context p(data_, skeleton_tets_, added_spheres_);
		std::unordered_set<uint32> face_ids;
		std::unordered_set<uint32> edge_ids;
		BoundaryTetPrune::collect_boundary_tet_mask(p, face_ids, edge_ids);
		mask.face_ids.assign(face_ids.begin(), face_ids.end());
		mask.edge_ids.assign(edge_ids.begin(), edge_ids.end());
		std::sort(mask.face_ids.begin(), mask.face_ids.end());
		std::sort(mask.edge_ids.begin(), mask.edge_ids.end());
		return mask;
	}

	// For display: sets face_score to the UDF score of the topology fix (optionally divided by the face area) on
	// the faces of the current tets and to NaN on the other faces. Returns false if the scores cannot be evaluated.
	template <typename ScoreEvaluator>
	bool compute_face_scores(ScoreEvaluator& evaluate_scores, bool normalize_by_area)
	{
		data_.face_score->fill(std::numeric_limits<Scalar>::quiet_NaN());
		Context p(data_, skeleton_tets_, added_spheres_);
		return TopologyScores::compute_face_scores(p, evaluate_scores, normalize_by_area);
	}

	// For display: sets edge_score to the UDF score of the topology fix on the edges of the current tets and to NaN
	// on the other edges. Returns false if the scores cannot be evaluated.
	template <typename ScoreEvaluator>
	bool compute_edge_scores(ScoreEvaluator& evaluate_scores)
	{
		data_.edge_score->fill(std::numeric_limits<Scalar>::quiet_NaN());
		Context p(data_, skeleton_tets_, added_spheres_);
		return TopologyScores::compute_edge_scores(p, evaluate_scores);
	}

	void prune_residual_sheets(ResidualPruneMode mode)
	{
		// the prune can remove spheres
		added_spheres_.clear();
		Context p(data_, skeleton_tets_, added_spheres_);
		NonManifoldPrune::run_residual_sheet_prune(p, mode);
	}

private:
	using Tets = typename Types::Tets;
	using Context = typename Types::Context;
	using MeshOps = topology_detail::SkeletonMeshOps<POINTS, NONMANIFOLD>;
	using TopologyScores = topology_detail::TopologyScores<POINTS, NONMANIFOLD>;
	using NonManifoldPrune = topology_detail::NonManifoldPrune<POINTS, NONMANIFOLD>;
	using GreedyTetRemoval = topology_detail::GreedyTetRemoval<POINTS, NONMANIFOLD>;
	using BoundaryTetPrune = topology_detail::BoundaryTetPrune<POINTS, NONMANIFOLD>;
	using DenseK5Repair = topology_detail::DenseK5Repair<POINTS, NONMANIFOLD>;
	using DenseFiveVertexDetectionResult = typename DenseK5Repair::DenseFiveVertexDetectionResult;
	using IndependentDenseFiveVertexRepairSeed = typename DenseK5Repair::IndependentDenseFiveVertexRepairSeed;

	//------------------------------//
	//-----Topology correction------//
	//------------------------------//

	template <typename ScoreEvaluator>
	bool run_topology_fix_pipeline(Context& p, ScoreEvaluator& evaluator)
	{
		const DenseFiveVertexDetectionResult dense_five_vertex_detection =
			DenseK5Repair::detect_k5_and_k5_minus_1_cells(p);
		const std::vector<IndependentDenseFiveVertexRepairSeed> independent_dense_repair_seeds =
			DenseK5Repair::collect_independent_dense_five_vertex_repair_seeds(p, dense_five_vertex_detection);
		// the cells that are not in a tet keep a zero score
		data_.edge_score->fill(Scalar(0));
		data_.face_score->fill(Scalar(0));
		if (!TopologyScores::compute_edge_scores(p, evaluator) || !TopologyScores::compute_face_scores(p, evaluator, false))
			return false;
		while (!p.skeleton_tets.empty())
		{
			const uint32 round_start_tets = static_cast<uint32>(p.skeleton_tets.size());
			BoundaryTetPrune::run_boundary_tet_face_deletion(p);
			if (!p.skeleton_tets.empty())
				GreedyTetRemoval::template run_edge_score_tet_mode_topology_fix<topology_detail::SimpleTetPolicy>(p);
			if (!p.skeleton_tets.empty())
				GreedyTetRemoval::template run_edge_score_tet_mode_topology_fix<topology_detail::NonSimpleTetPolicy>(p);
			const uint32 round_end_tets = static_cast<uint32>(p.skeleton_tets.size());
			const uint32 round_removed_tets = round_start_tets - round_end_tets;
			if (round_removed_tets == 0)
				break;
		}
		DenseK5Repair::repair_isolated_independent_dense_five_vertex_regions(p, independent_dense_repair_seeds);
		return true;
	}

	Data data_;
	Tets skeleton_tets_;
	std::vector<SphereVertex> added_spheres_;
};

} // namespace offset_axis

} // namespace geometry
} // namespace cgogn

#endif // CGOGN_GEOMETRY_ALGOS_OFFSET_AXIS_SKELETON_TOPOLOGY_H_
