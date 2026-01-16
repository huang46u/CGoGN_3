/*******************************************************************************
 * CGoGN: Combinatorial and Geometric modeling with Generic N-dimensional Maps  *
 * Copyright (C), IGG Group, ICube, University of Strasbourg, France            *
 *                                                                              *
 * This library is free software; you can redistribute it and/or modify it      *
 * under the terms of the GNU Lesser General Public License as published by the *
 * Free Software Foundation; either version 2.1 of the License, or (at your     *
 * option) any later version.                                                   *
 *                                                                              *
 * This library is distributed in the hope that it will be useful, but WITHOUT  *
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or        *
 * FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License  *
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

#ifndef CGOGN_GEOMETRY_TYPES_FAST_WINDING_NUMBER_TRAITS_H_
#define CGOGN_GEOMETRY_TYPES_FAST_WINDING_NUMBER_TRAITS_H_

#include <algorithm>
#include <cmath>
#include <numeric>
#include <vector>

#include <cgogn/core/functions/attributes.h>
#include <cgogn/core/utils/numerics.h>
#include <cgogn/core/functions/traversals/vertex.h>
#include <cgogn/geometry/types/vector_traits.h>
#include <libacc/bvh_tree.h>
#include <libacc/bvh_tree_spheres.h>

namespace cgogn
{

namespace geometry
{

template <typename MESH>
class FWN_Triangle_Traits
{
	template <typename T>
	using Attribute = typename mesh_traits<MESH>::template Attribute<T>;

	using Face = typename mesh_traits<MESH>::Face;
	using Vertex = typename mesh_traits<MESH>::Vertex;

public:
	using BVH = acc::BVHTree<uint32, Vec3>;
	using PrimID = std::size_t;

	FWN_Triangle_Traits(MESH& mesh, const BVH* bvh_tree, const std::vector<Face>& faces, const Attribute<Vec3>* vpos,
						const Attribute<Vec3>* fnorm, const Attribute<Scalar>* farea,
						const Attribute<Vec3>* fcentroid)
		: m_(mesh), bvh_tree_(bvh_tree), bvh_faces_(faces), vertex_position_(vpos), face_normal_(fnorm),
		  face_area_(farea), face_centroid_(fcentroid)
	{
		// Fix: traits now have consistent names and explicit BVH accessors.
	}

	inline void collect_primitives_id(std::size_t node_id, std::vector<PrimID>& primitives) const
	{
		bvh_tree_->collect_primitives_id(node_id, primitives);
	}

	inline void collect_all_primitives(std::vector<PrimID>& primitives) const
	{
		primitives.resize(bvh_faces_.size());
		std::iota(primitives.begin(), primitives.end(), 0);
	}

	inline void primitive_data(PrimID pid, Scalar& a, Vec3& c, Vec3& n) const
	{
		const Face f = bvh_faces_[pid];
		a = value<Scalar>(m_, face_area_, f);
		c = value<Vec3>(m_, face_centroid_, f);
		n = value<Vec3>(m_, face_normal_, f);
	}

	inline Scalar primitive_extent(PrimID pid, const Vec3& p_tilde) const
	{
		const Face f = bvh_faces_[pid];

		auto iv = incident_vertices(m_, f);
		Scalar max_norm = 0;
		for (Vertex v : iv)
		{
			Vec3 pos = value<Vec3>(m_, vertex_position_, v);
			Scalar norm = (pos - p_tilde).norm();
			max_norm = std::max(max_norm, norm);
		}
		return max_norm;
	}

	inline Scalar primitive_solid_angle(const Vec3& q, PrimID pid) const
	{
		const Face f = bvh_faces_[pid];
		return solid_angle(q, f) / (4 * M_PI);
	}

	template <int ORDER, class Coeff>
	inline void accumulate(Coeff& coeff, PrimID pid, const Vec3& p_tilde) const
	{
		Scalar a;
		Vec3 c, n;
		primitive_data(pid, a, c, n);

		const Vec3 r = c - p_tilde;

		if constexpr (ORDER < 3)
			coeff.add(r, n, a);
		else
		{
			const Face f = bvh_faces_[pid];
			auto iv = incident_vertices(m_, f);

			const Vec3 p1 = value<Vec3>(m_, vertex_position_, iv[0]);
			const Vec3 p2 = value<Vec3>(m_, vertex_position_, iv[1]);
			const Vec3 p3 = value<Vec3>(m_, vertex_position_, iv[2]);

			const Vec3 v1 = 0.5 * (p1 + p2) - p_tilde;
			const Vec3 v2 = 0.5 * (p2 + p3) - p_tilde;
			const Vec3 v3 = 0.5 * (p3 + p1) - p_tilde;

			coeff.add(v1, v2, v3, r, n, a);
		}
	}

private:
	Scalar solid_angle(const Vec3& q, const Face& f) const
	{
		auto iv = incident_vertices(m_, f);

		Vec3 v1 = value<Vec3>(m_, vertex_position_, iv[0]) - q;
		Vec3 v2 = value<Vec3>(m_, vertex_position_, iv[1]) - q;
		Vec3 v3 = value<Vec3>(m_, vertex_position_, iv[2]) - q;

		Scalar l1 = v1.norm();
		Scalar l2 = v2.norm();
		Scalar l3 = v3.norm();

		if (l1 == 0 || l2 == 0 || l3 == 0)
			return 0;

		v1 /= l1;
		v2 /= l2;
		v3 /= l3;

		const Scalar numerator = v1.dot((v2 - v1).cross(v3 - v1));
		if (numerator == 0)
			return 0;

		const Scalar denominator = 1 + v1.dot(v2) + v2.dot(v3) + v3.dot(v1);
		return 2 * std::atan2(numerator, denominator);
	}

private:
	MESH& m_;
	const BVH* bvh_tree_;
	const std::vector<Face>& bvh_faces_;

	const Attribute<Vec3>* vertex_position_;
	const Attribute<Vec3>* face_normal_;
	const Attribute<Scalar>* face_area_;
	const Attribute<Vec3>* face_centroid_;
};

template <typename MESH>
class FWN_Point_Traits
{
	template <typename T>
	using Attribute = typename mesh_traits<MESH>::template Attribute<T>;
	using Vertex = typename mesh_traits<MESH>::Vertex;

public:
	using BVH = acc::BVHTreeSpheres<uint32, Vec3>;
	using PrimID = std::size_t;

	FWN_Point_Traits(MESH& mesh, const BVH* bvh_tree, const std::vector<Vertex>& vertices, const Attribute<Vec3>* vpos,
					const Attribute<Vec3>* vnorm, const Attribute<Scalar>* varea)
		: m_(mesh), bvh_tree_(bvh_tree), bvh_vertices_(vertices), vertex_position_(vpos), vertex_normal_(vnorm),
		  vertex_area_(varea)
	{
	}

	Scalar area_to_radius(Scalar area) const
	{
		return std::sqrt(area / M_PI);
	}

	inline void collect_primitives_id(std::size_t node_id, std::vector<PrimID>& primitives) const
	{
		bvh_tree_->collect_primitives_id(node_id, primitives);
	}

	inline void collect_all_primitives(std::vector<PrimID>& primitives) const
	{
		primitives.resize(bvh_vertices_.size());
		std::iota(primitives.begin(), primitives.end(), 0);
	}

	inline void primitive_data(PrimID pid, Scalar& a, Vec3& c, Vec3& n) const
	{
		const Vertex v = bvh_vertices_[pid];
		a = value<Scalar>(m_, vertex_area_, v);
		c = value<Vec3>(m_, vertex_position_, v);
		n = value<Vec3>(m_, vertex_normal_, v);
	}

	inline Scalar primitive_extent(PrimID pid, const Vec3& p_tilde) const
	{
		const Vertex v = bvh_vertices_[pid];
		Vec3 pos = value<Vec3>(m_, vertex_position_, v);
		Scalar radius = area_to_radius(value<Scalar>(m_, vertex_area_, v));
		return (pos - p_tilde).norm() + radius;
	}

	inline Scalar primitive_solid_angle(const Vec3& q, PrimID pid) const
	{
		Scalar a;
		Vec3 c, n;
		primitive_data(pid, a, c, n);

		const Vec3 dir = c - q;
		const Scalar R = std::max(dir.norm(), Scalar(1e-20));
		const Scalar inv_r3 = Scalar(1.0) / (R * R * R);
		return (a * n.dot(dir)) * inv_r3 / (4 * M_PI);
	}

	template <int ORDER, class Coeff>
	inline void accumulate(Coeff& coeff, PrimID pid, const Vec3& p_tilde) const
	{
		Scalar a;
		Vec3 c, n;
		primitive_data(pid, a, c, n);

		const Vec3 r = c - p_tilde;

		if constexpr (ORDER < 3)
			coeff.add(r, n, a);
		else
			coeff.add(r, r, r, r, n, a);
	}

private:
	MESH& m_;
	const BVH* bvh_tree_;
	const std::vector<Vertex>& bvh_vertices_;

	const Attribute<Vec3>* vertex_position_;
	const Attribute<Vec3>* vertex_normal_;
	const Attribute<Scalar>* vertex_area_;
};

} // namespace geometry

} // namespace cgogn

#endif // CGOGN_GEOMETRY_TYPES_FAST_WINDING_NUMBER_TRAITS_H_
