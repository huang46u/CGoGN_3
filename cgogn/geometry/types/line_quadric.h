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

#ifndef CGOGN_GEOMETRY_TYPES_LINE_Quadric_H_
#define CGOGN_GEOMETRY_TYPES_LINE_Quadric_H_

#include <cgogn/core/utils/numerics.h>
#include <cgogn/geometry/types/quadric.h>
#include <cgogn/geometry/types/vector_traits.h>

namespace cgogn
{

namespace geometry
{

class Line_Quadric
{
public:
	inline Line_Quadric()
	{
		q_.zero();
	}

	inline Line_Quadric(const Vec3& p, const Vec3& normal)
	{
		auto [e2, e3] = gram_schimdt(normal);
		Quadric q1 = Quadric(p, e2);
		Quadric q2 = Quadric(p, e3);

		q_ += q1 ;
		q_ += q2 ;
	}

	Line_Quadric(const Line_Quadric& q)
	{
		q_ = q.q_;
	}

	inline void zero()
	{
		q_.zero();
	}

	Line_Quadric& operator=(const Line_Quadric& q)
	{
		q_ = q.q_;
		return *this;
	}

	Line_Quadric& operator+=(const Line_Quadric& q)
	{
		q_ += q.q_;
		return *this;
	}
	Line_Quadric operator+(const Line_Quadric& q)
	{
		Line_Quadric res(*this);
		res += q;
		return res;
	}
	Line_Quadric& operator*=(Scalar s)
	{
		q_ *= s;
		return *this;
	}
	Line_Quadric operator*(Scalar s)
	{
		Line_Quadric res(*this);
		res *= s;
		return res;
	}
	Scalar eval(const Vec3& v)
	{
		return eval(Vec4{v[0], v[1], v[2], 1.});
	}

	inline Scalar eval(const Vec4& v)
	{
		return q_.eval(v);
	}

	bool optimized(Vec3& v)
	{
		return q_.optimized(v);
	}

	Quadric get_quadric() const
	{
		return q_;
	}

private:
	std::pair<Vec3, Vec3> gram_schimdt(const Vec3& n)
	{
		Vec3 e1 = n.normalized();
		Vec3 v2 = (std::abs(e1.x()) < std::abs(e1.y())) ? Vec3(1, 0, 0) : Vec3(0, 1, 0);
		Vec3 e2 = (v2 - v2.dot(e1) * e1).normalized();
		Vec3 e3 = e1.cross(e2);
		return {e2, e3};
	}

private:
	Quadric q_;
	
};

} // namespace geometry

} // namespace cgogn

#endif // CGOGN_GEOMETRY_TYPES_Line_Quadric_H_