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
 *******************************************************************************/

#ifndef CGOGN_GEOMETRY_TYPES_LINE_QUADRIC_H_
#define CGOGN_GEOMETRY_TYPES_LINE_QUADRIC_H_

#include <cmath>

#include <cgogn/geometry/types/vector_traits.h>

namespace cgogn
{

namespace geometry
{

struct Line_Quadric
{
	Line_Quadric()
	{
		clear();
	}

	Line_Quadric(const Vec3& p, const Vec3& n)
	{
		clear();

		Vec3 direction = n;
		const Scalar squared_norm = direction.squaredNorm();
		if (!(squared_norm > 0.0) || !std::isfinite(squared_norm))
			return;
		direction /= std::sqrt(squared_norm);

		const Mat3 projection = Mat3::Identity() - direction * direction.transpose();
		_A.template block<3, 3>(0, 0) = 2.0 * projection;
		const Vec4 position(p.x(), p.y(), p.z(), 0.0);
		_b = _A * position;
		_c = Scalar(0.5 * position.transpose() * _A * position);
	}

	void clear()
	{
		_A.setZero();
		_b.setZero();
		_c = 0.0;
	}

	Line_Quadric& operator=(const Line_Quadric& q)
	{
		_A = q._A;
		_b = q._b;
		_c = q._c;
		return *this;
	}

	Line_Quadric& operator+=(const Line_Quadric& q)
	{
		_A += q._A;
		_b += q._b;
		_c += q._c;
		return *this;
	}

	Line_Quadric& operator*=(Scalar s)
	{
		_A *= s;
		_b *= s;
		_c *= s;
		return *this;
	}

	friend Line_Quadric operator*(const Line_Quadric& q, Scalar s)
	{
		Line_Quadric result(q);
		result *= s;
		return result;
	}

	friend Line_Quadric operator+(const Line_Quadric& lhs, const Line_Quadric& rhs)
	{
		Line_Quadric result(lhs);
		result += rhs;
		return result;
	}

	Scalar eval(const Vec3& p) const
	{
		return eval(Vec4(p.x(), p.y(), p.z(), 0.0));
	}

	Scalar eval(const Vec4& p) const
	{
		return Scalar(0.5 * p.transpose() * _A * p) - Scalar(_b.transpose() * p) + _c;
	}

	Mat4 _A;
	Vec4 _b;
	Scalar _c = 0.0;
};

} // namespace geometry

} // namespace cgogn

#endif // CGOGN_GEOMETRY_TYPES_LINE_QUADRIC_H_
