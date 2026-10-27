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
#include <cgogn/geometry/algos/offset_axis/internal/skeleton_tets.h>
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
	};

	using ResidualPruneMode = offset_axis::ResidualPruneOption;

	using SpheresOptimizerType = SpheresOptimizer<POINTS>;
	using FaceKey = std::array<uint32, 3>;

	using Tets = SkeletonTets<NONMANIFOLD>;
	using Tet = typename Tets::Tet;

	// What the topology functions work on: the topology data, the data of the spheres optimizer (samples and spheres),
	// the current tets and the spheres added by the dense K5 repair.
	struct Context
	{
		Data& data;
		typename SpheresOptimizerType::Data& optimizer_data;
		Tets& skeleton_tets;
		std::vector<SphereVertex>& added_spheres;

		Context(Data& d, Tets& t, std::vector<SphereVertex>& added)
			: data(d), optimizer_data(d.spheres_optimizer.data()), skeleton_tets(t), added_spheres(added)
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
