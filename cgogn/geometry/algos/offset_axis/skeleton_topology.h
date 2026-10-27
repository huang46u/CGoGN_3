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
#include <cgogn/geometry/algos/offset_axis/internal/topology_scores.h>
#include <cgogn/geometry/algos/offset_axis/reconstruction_options.h>
#include <cgogn/geometry/algos/offset_axis/spheres_optimizer.h>

#include <algorithm>
#include <limits>
#include <set>
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
	using SphereVertex = typename mesh_traits<POINTS>::Vertex;
	using SkeletonVertex = typename mesh_traits<NONMANIFOLD>::Vertex;
	using SkeletonEdge = typename mesh_traits<NONMANIFOLD>::Edge;
	using SkeletonFace = typename mesh_traits<NONMANIFOLD>::Face;
	template <typename T>
	using SphereAttribute = typename mesh_traits<POINTS>::template Attribute<T>;
	template <typename T>
	using SkeletonAttribute = typename mesh_traits<NONMANIFOLD>::template Attribute<T>;
	using ResidualPruneMode = ResidualPruneOption;

	// Borrowed: the meshes and the attributes belong to the caller.
	struct Data
	{
		NONMANIFOLD& skeleton;
		SpheresOptimizer<POINTS>& spheres_optimizer;
		SphereAttribute<SkeletonVertex>* sphere_skeleton_vertex = nullptr;
		SkeletonAttribute<Vec3>* skeleton_position = nullptr;
		SkeletonAttribute<Scalar>* skeleton_radius = nullptr;
		SkeletonAttribute<SphereVertex>* skeleton_source_sphere = nullptr;
		SkeletonAttribute<std::set<std::size_t>>* face_incident_tets = nullptr;
		SkeletonAttribute<uint32>* edge_degree = nullptr;
		SkeletonAttribute<Scalar>* face_score = nullptr;
		SkeletonAttribute<Scalar>* edge_score = nullptr;
	};

	explicit SkeletonTopology(Data data) : data_(std::move(data)), tets_(data_.skeleton, data_.face_incident_tets)
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
		return !tets_.empty();
	}
	void build_from_spheres()
	{
		added_spheres_.clear();
		const auto& spheres = data_.spheres_optimizer.data();
		build_skeleton_from_spheres(spheres.spheres, spheres.sphere_position, spheres.sphere_radius,
									spheres.sphere_neighbors, data_.skeleton, data_.skeleton_position,
									data_.skeleton_radius, data_.skeleton_source_sphere, data_.sphere_skeleton_vertex,
									data_.edge_degree, tets_);
	}

	// Step 10: removes the tets of the skeleton (boundary tets, then greedy by score, until no tet is removed) and
	// repairs the dense K5 regions found at the start. Returns false if the scores cannot be evaluated.
	template <typename ScoreEvaluator>
	bool fix_topology(ScoreEvaluator&& evaluate_scores)
	{
		using namespace topology_detail;

		added_spheres_.clear();
		if (tets_.empty())
			return true;
		NONMANIFOLD& skeleton = data_.skeleton;
		const std::vector<K5Region<POINTS>> k5_regions =
			find_k5_regions(skeleton, tets_, data_.spheres_optimizer.data().spheres, data_.skeleton_source_sphere);
		// the cells that are not in a tet keep a zero score
		data_.edge_score->fill(Scalar(0));
		data_.face_score->fill(Scalar(0));
		if (!topology_detail::compute_edge_scores(skeleton, data_.skeleton_position, tet_edge_ids(skeleton, tets_),
												  evaluate_scores, data_.edge_score) ||
			!topology_detail::compute_face_scores(skeleton, data_.skeleton_position, tet_face_ids(skeleton, tets_),
												  evaluate_scores, false, data_.face_score))
			return false;
		while (!tets_.empty())
		{
			const std::size_t nb_tets = tets_.size();
			remove_boundary_tets(skeleton, tets_, data_.skeleton_position, data_.spheres_optimizer,
								 data_.skeleton_source_sphere, data_.sphere_skeleton_vertex);
			if (!tets_.empty())
				remove_tets_greedy<SimpleTetPolicy>(skeleton, tets_, data_.face_score, data_.edge_score,
													data_.spheres_optimizer, data_.skeleton_source_sphere,
													data_.sphere_skeleton_vertex);
			if (!tets_.empty())
				remove_tets_greedy<NonSimpleTetPolicy>(skeleton, tets_, data_.face_score, data_.edge_score,
													   data_.spheres_optimizer, data_.skeleton_source_sphere,
													   data_.sphere_skeleton_vertex);
			if (tets_.size() == nb_tets)
				break;
		}
		repair_k5_regions(skeleton, tets_, data_.spheres_optimizer, data_.skeleton_position, data_.skeleton_radius,
						  data_.skeleton_source_sphere, data_.sphere_skeleton_vertex, data_.edge_degree, k5_regions,
						  added_spheres_);
		return true;
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
		std::unordered_set<uint32> face_ids;
		std::unordered_set<uint32> edge_ids;
		topology_detail::boundary_tet_mask(data_.skeleton, tets_, data_.skeleton_position, face_ids, edge_ids);
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
		return topology_detail::compute_face_scores(data_.skeleton, data_.skeleton_position,
													topology_detail::tet_face_ids(data_.skeleton, tets_),
													evaluate_scores, normalize_by_area, data_.face_score);
	}

	// For display: sets edge_score to the UDF score of the topology fix on the edges of the current tets and to NaN
	// on the other edges. Returns false if the scores cannot be evaluated.
	template <typename ScoreEvaluator>
	bool compute_edge_scores(ScoreEvaluator& evaluate_scores)
	{
		data_.edge_score->fill(std::numeric_limits<Scalar>::quiet_NaN());
		return topology_detail::compute_edge_scores(data_.skeleton, data_.skeleton_position,
													topology_detail::tet_edge_ids(data_.skeleton, tets_),
													evaluate_scores, data_.edge_score);
	}

	// Step 11.
	void prune_residual_sheets(ResidualPruneMode mode)
	{
		topology_detail::prune_residual_sheets(data_.skeleton, tets_, data_.edge_degree, data_.spheres_optimizer,
											   data_.skeleton_source_sphere, data_.sphere_skeleton_vertex, mode);
		// the prune can remove added spheres: a removed sphere is unlinked from its skeleton vertex
		const POINTS& spheres = data_.spheres_optimizer.data().spheres;
		added_spheres_.erase(
			std::remove_if(added_spheres_.begin(), added_spheres_.end(),
						   [&](SphereVertex sphere) {
							   return !value<SkeletonVertex>(spheres, data_.sphere_skeleton_vertex, sphere).is_valid();
						   }),
			added_spheres_.end());
	}

private:
	Data data_;
	SkeletonTets<NONMANIFOLD> tets_;
	std::vector<SphereVertex> added_spheres_;
};

} // namespace offset_axis

} // namespace geometry
} // namespace cgogn

#endif // CGOGN_GEOMETRY_ALGOS_OFFSET_AXIS_SKELETON_TOPOLOGY_H_
