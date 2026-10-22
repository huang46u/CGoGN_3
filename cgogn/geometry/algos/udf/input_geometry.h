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
#ifndef CGOGN_GEOMETRY_ALGOS_UDF_INPUT_GEOMETRY_H_
#define CGOGN_GEOMETRY_ALGOS_UDF_INPUT_GEOMETRY_H_

#include <cgogn/core/functions/attributes.h>
#include <cgogn/core/functions/mesh_info.h>
#include <cgogn/core/functions/traversals/face.h>
#include <cgogn/core/functions/traversals/vertex.h>
#include <cgogn/core/utils/numerics.h>
#include <cgogn/geometry/algos/normal.h>
#include <cgogn/geometry/algos/udf/alpha_projection.h>
#include <cgogn/geometry/algos/udf/normal_utils.h>
#include <cgogn/geometry/algos/udf/spatial_query.h>
#include <cgogn/geometry/types/vector_traits.h>

#include <libacc/bvh_tree.h>
#include <libacc/kd_tree.h>

#include <cmath>
#include <memory>
#include <vector>

namespace cgogn::geometry
{

// Caches of the input geometry used by the sampling, the projections and the topology scores of the UDF
// reconstruction: the k-d tree and the vertex list of the input point cloud, and the BVH of the input surface with its
// face and vertex lists. The caches are built on demand and dropped by invalidate(). The input surface, the input
// point cloud and their attributes are owned by the caller: bind() only stores pointers to them.
template <typename SURFACE, typename POINTS>
class InputGeometry
{
public:
	using PVertex = typename mesh_traits<POINTS>::Vertex;
	using SVertex = typename mesh_traits<SURFACE>::Vertex;
	using SFace = typename mesh_traits<SURFACE>::Face;
	using PAttribute = typename mesh_traits<POINTS>::template Attribute<Vec3>;
	template <typename T>
	using PointAttribute = typename mesh_traits<POINTS>::template Attribute<T>;
	template <typename T>
	using SurfaceAttribute = typename mesh_traits<SURFACE>::template Attribute<T>;

	// non-owning pointers to the input (the owner is ReconstructionState::data)
	struct
	{
		SURFACE* surface = nullptr;
		POINTS* input_points = nullptr;
	} data_;

	void bind(SURFACE* surface, POINTS* input_points, std::shared_ptr<PAttribute> input_position,
			  std::shared_ptr<PAttribute> input_normal)
	{
		data_.surface = surface;
		data_.input_points = input_points;
		input_position_ = input_position;
		input_normal_ = input_normal;
	}

	void invalidate();
	void build_input_kdtree();
	void build_surface_bvh();
	bool query_surface_projection_info(const Vec3& query, AlphaProjectionInfo& info);
	bool query_point_cloud_projection_info(const Vec3& query, AlphaProjectionInfo& info);

	const std::unique_ptr<acc::KDTree<3, uint32>>& input_kdtree() const
	{
		return input_kdtree_;
	}
	std::vector<PVertex>& input_vertices()
	{
		return input_vertices_;
	}
	const std::unique_ptr<acc::BVHTree<uint32, Vec3>>& surface_bvh() const
	{
		return surface_bvh_;
	}
	std::vector<SFace>& surface_faces()
	{
		return surface_faces_;
	}
	const std::shared_ptr<SurfaceAttribute<Vec3>>& surface_vertex_normal() const
	{
		return surface_vertex_normal_;
	}

private:
	std::shared_ptr<PointAttribute<Vec3>> input_position_, input_normal_;
	std::unique_ptr<acc::KDTree<3, uint32>> input_kdtree_;
	std::vector<PVertex> input_vertices_;
	std::unique_ptr<acc::BVHTree<uint32, Vec3>> surface_bvh_;
	std::vector<SFace> surface_faces_;
	std::vector<SVertex> surface_vertices_;
	std::vector<Vec3> surface_positions_;
	std::shared_ptr<SurfaceAttribute<Vec3>> surface_vertex_normal_;
};

template <typename SURFACE, typename POINTS>
void InputGeometry<SURFACE, POINTS>::invalidate()
{
	input_kdtree_.reset();
	input_vertices_.clear();
	surface_bvh_.reset();
	surface_faces_.clear();
	surface_vertices_.clear();
	surface_positions_.clear();
	surface_vertex_normal_.reset();
}

template <typename SURFACE, typename POINTS>
void InputGeometry<SURFACE, POINTS>::build_input_kdtree()
{
	if (input_kdtree_)
		return;
	input_kdtree_.reset();
	input_vertices_.clear();
	if (!data_.input_points || !input_position_)
		return;
	std::vector<Vec3> positions;
	foreach_cell(*data_.input_points, [&](PVertex v) {
		const uint32 idx = index_of(*data_.input_points, v);
		positions.push_back((*input_position_)[idx]);
		input_vertices_.push_back(v);
		return true;
	});
	if (!positions.empty())
		input_kdtree_ = std::make_unique<acc::KDTree<3, uint32>>(positions);
}

template <typename SURFACE, typename POINTS>
void InputGeometry<SURFACE, POINTS>::build_surface_bvh()
{
	if (surface_bvh_)
		return;
	surface_bvh_.reset();
	surface_faces_.clear();
	surface_vertices_.clear();
	surface_positions_.clear();
	if (!data_.surface)
		return;
	auto position = get_attribute<Vec3, SVertex>(*data_.surface, "position");
	if (!position)
		return;
	surface_vertex_normal_ = get_or_add_attribute<Vec3, SVertex>(*data_.surface, "normal");
	geometry::compute_normal<SVertex>(*data_.surface, position.get(), surface_vertex_normal_.get());
	auto vertex_index = get_or_add_attribute<uint32, SVertex>(*data_.surface, "__udf_bvh_vertex_index");
	uint32 next = 0;
	foreach_cell(*data_.surface, [&](SVertex v) {
		surface_vertices_.push_back(v);
		surface_positions_.push_back((*position)[index_of(*data_.surface, v)]);
		value<uint32>(*data_.surface, vertex_index, v) = next++;
		return true;
	});
	std::vector<uint32> face_indices;
	foreach_cell(*data_.surface, [&](typename mesh_traits<SURFACE>::Face f) {
		surface_faces_.push_back(f);
		foreach_incident_vertex(*data_.surface, f, [&](SVertex v) {
			face_indices.push_back(value<uint32>(*data_.surface, vertex_index, v));
			return true;
		});
		return true;
	});
	remove_attribute<SVertex>(*data_.surface, vertex_index);
	if (!face_indices.empty() && !surface_positions_.empty())
		surface_bvh_ = std::make_unique<acc::BVHTree<uint32, Vec3>>(face_indices, surface_positions_);
}

template <typename SURFACE, typename POINTS>
bool InputGeometry<SURFACE, POINTS>::query_surface_projection_info(const Vec3& query, AlphaProjectionInfo& info)
{
	info = AlphaProjectionInfo{};
	if (!data_.surface || !surface_bvh_ || !surface_vertex_normal_)
		return false;
	uint32 primitive;
	if (!query_surface_closest_point(*surface_bvh_, query, primitive, info.closest_point, info.distance))
		return false;
	info.normal = query - info.closest_point;
	if (info.normal.squaredNorm() < Scalar(1e-12) && primitive < surface_faces_.size())
		info.normal = geometry::normal(*data_.surface, surface_faces_[primitive], surface_vertex_normal_.get());
	info.normal = unit_normal_or_z(info.normal);
	info.valid =
		info.closest_point.allFinite() && info.normal.allFinite() && std::isfinite(static_cast<double>(info.distance));
	return info.valid;
}

template <typename SURFACE, typename POINTS>
bool InputGeometry<SURFACE, POINTS>::query_point_cloud_projection_info(const Vec3& query, AlphaProjectionInfo& info)
{
	info = AlphaProjectionInfo{};
	if (!input_kdtree_ || !data_.input_points || input_vertices_.empty())
		return false;
	uint32 point_index;
	if (!query_point_cloud_nearest_point(*input_kdtree_, query, point_index, info.distance) ||
		point_index >= input_vertices_.size())
		return false;
	const uint32 vertex_index = index_of(*data_.input_points, input_vertices_[point_index]);
	if (vertex_index == INVALID_INDEX)
		return false;
	info.closest_point = (*input_position_)[vertex_index];
	info.normal = query - info.closest_point;
	if (info.normal.squaredNorm() < Scalar(1e-12))
		info.normal = (*input_normal_)[vertex_index];
	info.normal = unit_normal_or_z(info.normal);
	info.valid =
		info.closest_point.allFinite() && info.normal.allFinite() && std::isfinite(static_cast<double>(info.distance));
	return info.valid;
}

} // namespace cgogn::geometry

#endif // CGOGN_GEOMETRY_ALGOS_UDF_INPUT_GEOMETRY_H_
