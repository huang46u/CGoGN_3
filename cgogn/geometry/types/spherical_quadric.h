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

#ifndef CGOGN_GEOMETRY_TYPES_SPHEREICAL_QUADRIC_H_
#define CGOGN_GEOMETRY_TYPES_SPHEREICAL_QUADRIC_H_

#include <cgogn/geometry/types/vector_traits.h>

namespace cgogn
{

namespace geometry
{

enum class SQEM_CASE
{
	Case1_Full,
	Case2_Line,
	Case3_Plane,
	Case4_Degenerate
};

struct Spherical_Quadric
{
	

	Spherical_Quadric()
	{
		this->clear();
	}

	void clear()
	{
		_A.setZero();
		_b.setZero();
		_c = 0;
	}

	inline Spherical_Quadric(const Vec4& p, const Vec4& n)
	{
		_A = 2 * (n * n.transpose());
		_b = 2 * (n.transpose() * p * n.transpose());
		_c = n.transpose() * p * p.transpose() * n;
	}

	Spherical_Quadric& operator=(const Spherical_Quadric& q)
	{
		_A = q._A;
		_b = q._b;
		_c = q._c;
		return *this;
	}

	Spherical_Quadric& operator+=(const Spherical_Quadric& q)
	{
		_A += q._A;
		_b += q._b;
		_c += q._c;
		return *this;
	}

	Spherical_Quadric& operator*=(Scalar s)
	{
		_A *= s;
		_b *= s;
		_c *= s;
		return *this;
	}

	friend Spherical_Quadric operator*(const Spherical_Quadric& q, Scalar s)
	{
		Spherical_Quadric result(q);
		result *= s;
		return result;
	}

	friend Spherical_Quadric operator+(const Spherical_Quadric& lhs, const Spherical_Quadric& rhs)
	{
		Spherical_Quadric result(lhs);
		result += rhs;
		return result;
	}

	Scalar eval(const Vec4& p) const
	{
		return (0.5 * p.transpose() * _A * p) - (Scalar)(_b.transpose() * p) + _c;
	}

	Vec4 gradient(const Vec4& p) const
	{
		return _A * p - _b;
	}

	SQEM_CASE well_conditioned(Scalar& r_out) const
	{
		Eigen::JacobiSVD<Mat4> svd(_A, Eigen::ComputeFullU | Eigen::ComputeFullV);
		Mat4 U = svd.matrixU();
		Vec4 S = svd.singularValues();
		Mat4 V = svd.matrixV();
		std::vector<double> sorted_sv{S(0), S(1), S(2), S(3)};

		std::sort(sorted_sv.begin(), sorted_sv.end(), std::greater<>());
		int rank = 0;
		for(int i = 0; i < 4; ++i)
		{
			if(sorted_sv[i] > 1e-6)
				rank++;
		}
		if(rank>1){
			Mat4 Sp = Mat4::Zero();
			for(int i = 0; i < 4; ++i)
			{
				if(S(i) > 1e-6)
					Sp(i,i) = 1.0 / S(i);
				else 
					Sp(i,i) = 0;	
			}
			Vec4 m = V * Sp * U.transpose() * _b;
			r_out = m(3);
		}
		switch (rank)
		{
		case 4:
			return SQEM_CASE::Case1_Full;
		case 3:
			return SQEM_CASE::Case2_Line;
		case 2:
			return SQEM_CASE::Case3_Plane;
		default:
			return SQEM_CASE::Case4_Degenerate;
		}
	}

	bool optimized(Vec4& sphere)
	{
		Mat4 inverse;
		inverse.setZero();
		Scalar determinant;
		bool invertible;
		_A.computeInverseAndDetWithCheck(inverse, determinant, invertible, 0.01);
		if (invertible)
			sphere = inverse * _b;
		return invertible;
	}

	friend std::ostream& operator<<(std::ostream& os, const Spherical_Quadric& quad)
	{
		os << quad._A << ", " << quad._b << ", " << quad._c;
		return os;
	}

	Mat4 _A;
	Vec4 _b;
	Scalar _c = 0;
};

} // namespace geometry

} // namespace cgogn

#endif // CGOGN_GEOMETRY_TYPES_SPHERICAL_QUADRIC_H_
