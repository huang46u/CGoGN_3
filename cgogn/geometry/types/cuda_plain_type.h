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
#include <cuda_runtime.h>
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
    float4 A[4];   
    float4 b;      
    float c;      
};

template <>
struct CudaPlainValueTraits<Spherical_Quadric>
{
    using PlainType = PlainSphericalQuadric;

    static PlainType to_plain(const Spherical_Quadric& value)
    {
        PlainType plain;
        const auto A = value.A();
        const auto b = value.b();

        for (int r = 0; r < 4; ++r)
        {
            plain.A[r] = make_float4(
                static_cast<float>(A(r, 0)),
                static_cast<float>(A(r, 1)),
                static_cast<float>(A(r, 2)),
                static_cast<float>(A(r, 3)));
        }
        plain.b = make_float4(
            static_cast<float>(b(0)),
            static_cast<float>(b(1)),
            static_cast<float>(b(2)),
            static_cast<float>(b(3)));
        plain.c = static_cast<float>(value.c());
        return plain;
    }

    static Spherical_Quadric from_plain(const PlainType& plain)
    {
        Spherical_Quadric q;
        for (int r = 0; r < 4; ++r)
        {
            q._A(r,0) = plain.A[r].x;
            q._A(r,1) = plain.A[r].y;
            q._A(r,2) = plain.A[r].z;
            q._A(r,3) = plain.A[r].w;
        }
        q._b(0) = plain.b.x;
        q._b(1) = plain.b.y;
        q._b(2) = plain.b.z;
        q._b(3) = plain.b.w;
        q._c    = plain.c;
        return q;
    }

    static constexpr bool can_write_back = true;
};
} // namespace cuda
} // namespace cgogn
#endif /* CGOGN_GEOMETRY_TYPES_CUDA_PLAIN_TYPE_H_ */