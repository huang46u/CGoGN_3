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
#ifndef CGOGN_GEOMETRY_ALGOS_OFFSET_AXIS_INTERNAL_SKELETON_TOPOLOGY_TYPES_H_
#define CGOGN_GEOMETRY_ALGOS_OFFSET_AXIS_INTERNAL_SKELETON_TOPOLOGY_TYPES_H_

#include <cgogn/core/functions/attributes.h>
#include <cgogn/core/utils/numerics.h>
#include <cgogn/geometry/algos/offset_axis/reconstruction_options.h>
#include <cgogn/geometry/algos/offset_axis/spheres_optimizer.h>
#include <cgogn/geometry/types/vector_traits.h>

#include <libacc/kd_tree.h>

#include <algorithm>
#include <array>
#include <cstddef>
#include <functional>
#include <memory>
#include <set>
#include <unordered_map>
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
struct SkeletonTopologyTypes
{
	using SphereVertex = typename mesh_traits<POINTS>::Vertex;
	using SkeletonVertex = typename mesh_traits<NONMANIFOLD>::Vertex;
	using SkeletonEdge = typename mesh_traits<NONMANIFOLD>::Edge;
	using SkeletonFace = typename mesh_traits<NONMANIFOLD>::Face;
	template <typename T>
	using SphereAttribute = typename mesh_traits<POINTS>::template Attribute<T>;
	template <typename T>
	using SkeletonAttribute = typename mesh_traits<NONMANIFOLD>::template Attribute<T>;

	struct Data
	{
		NONMANIFOLD* skeleton = nullptr;
		SpheresOptimizer<POINTS>* spheres_optimizer = nullptr;
		std::shared_ptr<SphereAttribute<SkeletonVertex>> sphere_skeleton_vertex;
		std::shared_ptr<SkeletonAttribute<Vec3>> skeleton_position;
		std::shared_ptr<SkeletonAttribute<Scalar>> skeleton_radius;
		std::shared_ptr<SkeletonAttribute<SphereVertex>> skeleton_source_sphere;
		std::shared_ptr<SkeletonAttribute<std::set<std::size_t>>> face_incident_tets;
		std::shared_ptr<SkeletonAttribute<uint32>> edge_degree;
	};

	struct Metrics
	{
		uint32 residual_removed_sheets = 0;
		uint32 residual_removed_faces = 0;
		int32 residual_chi_before = 0;
		int32 residual_chi_after = 0;
		uint32 residual_fully_nm_removed_faces = 0;
		uint32 residual_small_sheets = 0;
		uint32 residual_small_sheets_guarded = 0;
		uint32 residual_small_sheets_with_boundary_edge = 0;
		uint32 residual_small_sheets_without_boundary_edge = 0;
		std::vector<SphereVertex> added_spheres;
	};

	enum class Status
	{
		success,
		invalid_data,
		score_evaluation_failed
	};
	using ResidualPruneMode = offset_axis::ResidualPruneOption;

	using SpheresOptimizerType = SpheresOptimizer<POINTS>;
	using FaceKey = std::array<uint32, 3>;

	struct Tet
	{
		SkeletonFace faces[4];
	};
	using TetMap = std::unordered_map<std::size_t, Tet>;

	// What the topology functions work on: the topology data, the data of the spheres optimizer (samples and spheres),
	// the current tets, the metrics of the running step and the residual prune mode.
	struct Context
	{
		Data& data;
		typename SpheresOptimizerType::Data& optimizer_data;
		TetMap& skeleton_tets;
		Metrics& metrics;
		ResidualPruneMode residual_prune_mode;

		Context(Data& d, TetMap& t, Metrics& m, ResidualPruneMode mode)
			: data(d), optimizer_data(d.spheres_optimizer->data()), skeleton_tets(t), metrics(m), residual_prune_mode(mode)
		{
		}
	};

	struct edge_hash
	{
		std::size_t operator()(const std::pair<uint32, uint32>& edge) const
		{
			const uint32 a = std::min(edge.first, edge.second);
			const uint32 b = std::max(edge.first, edge.second);
			return std::hash<uint64>()((uint64(a) << 32) | uint64(b));
		}
	};

	struct edge_equal
	{
		bool operator()(const std::pair<uint32, uint32>& edge1, const std::pair<uint32, uint32>& edge2) const
		{
			return ((edge1.first == edge2.first && edge1.second == edge2.second) ||
					(edge1.first == edge2.second && edge1.second == edge2.first));
		}
	};

	struct RemovalStats
	{
		uint32 removed_faces = 0;
		uint32 removed_edges = 0;
		uint32 removed_tets = 0;

		RemovalStats& operator+=(const RemovalStats& o)
		{
			removed_faces += o.removed_faces;
			removed_edges += o.removed_edges;
			removed_tets += o.removed_tets;
			return *this;
		}
	};

	struct SkeletonFaceDeletionOptions
	{
		bool remove_incident_orphan_edges_immediately = true;
		bool remove_global_orphan_elements_after_batch = true;
	};

	struct TopologyFixScoreCache
	{
		std::unordered_map<uint32, Scalar> edge_scores;
		std::unordered_map<uint32, Scalar> face_scores;
		bool initialized = false;
	};
};

} // namespace topology_detail

} // namespace offset_axis

} // namespace geometry

} // namespace cgogn

#endif // CGOGN_GEOMETRY_ALGOS_OFFSET_AXIS_INTERNAL_SKELETON_TOPOLOGY_TYPES_H_
