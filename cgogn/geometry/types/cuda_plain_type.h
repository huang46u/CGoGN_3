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

#ifndef CGOGN_GEOMETRY_TYPES_CUDA_PLAIN_TYPE_H_
#define CGOGN_GEOMETRY_TYPES_CUDA_PLAIN_TYPE_H_

#include <cgogn/core/types/cuda_plain_traits.h>
#include <cgogn/geometry/types/spherical_quadric.h>
#include <cgogn/geometry/types/vector_traits.h>
#include <vector_types.h>
namespace cgogn
{
namespace cuda
{
	using Scalar = cgogn::geometry::Scalar;
	using Vec3 = cgogn::geometry::Vec3;
	using Spherical_Quadric = cgogn::geometry::Spherical_Quadric; 

template <>
struct CudaPlainValueTraits<Scalar>
{
	using PlainType = float;

	static PlainType to_plain(const Scalar& value)
	{
		return static_cast<PlainType>(value);
	}

	static Scalar from_plain(const PlainType& plain)
	{
		return static_cast<Scalar>(plain);
	}

	static constexpr bool can_write_back = true;
};

template <>
struct CudaPlainValueTraits<cgogn::geometry::Vec3>
{
	using PlainType = float4;

	static PlainType to_plain(const Vec3& value)
	{
		return make_float4(static_cast<float>(value[0]), static_cast<float>(value[1]), static_cast<float>(value[2]), 0.f);
	}

	static Vec3 from_plain(const PlainType& plain)
	{
		return Vec3(static_cast<Scalar>(plain.x), static_cast<Scalar>(plain.y), static_cast<Scalar>(plain.z));
	}

	static constexpr bool can_write_back = true;
};
struct PlainSphericalQuadric
{
	float A[16];
	float b[4];
	float c;
};

template <>
struct CudaPlainValueTraits<Spherical_Quadric>
{
	using PlainType = PlainSphericalQuadric;
	static PlainType to_plain(const Spherical_Quadric& value)
	{
		PlainType plain;
		auto A = value.A();
		auto b = value.b();
		for (int r = 0; r < 4; ++r)
			for (int c = 0; c < 4; ++c)
				plain.A[r * 4 + c] = static_cast<float>(A(r, c));
		for (int i = 0; i < 4; ++i)
			plain.b[i] = static_cast<float>(b(i));
		plain.c = static_cast<float>(value.c());
		return plain;
	}

	static Spherical_Quadric from_plain(const PlainType& plain)
	{
		Spherical_Quadric q;
		int idx = 0;
		for (int r = 0; r < 4; ++r)
			for (int c = 0; c < 4; ++c)
				q._A(r, c) = plain.A[idx++];
		for (int i = 0; i < 4; ++i)
			q._b(i) = plain.b[i];
		q._c = plain.c;
		return q;
	}

	static constexpr bool can_write_back = true;
};
} // namespace cuda
} // namespace cgogn
#endif /* CGOGN_GEOMETRY_TYPES_CUDA_PLAIN_TYPE_H_ */