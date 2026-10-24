#ifndef CGOGN_GEOMETRY_ALGOS_OFFSET_AXIS_NORMAL_UTILS_H_
#define CGOGN_GEOMETRY_ALGOS_OFFSET_AXIS_NORMAL_UTILS_H_

#include <cgogn/geometry/types/vector_traits.h>

namespace cgogn::geometry
{

// unit vector along v, or +z when v is (numerically) zero
inline Vec3 unit_normal_or_z(const Vec3& v)
{
	if (v.squaredNorm() < Scalar(1e-12))
		return Vec3(0, 0, 1);
	return v.normalized();
}

} // namespace cgogn::geometry

#endif // CGOGN_GEOMETRY_ALGOS_OFFSET_AXIS_NORMAL_UTILS_H_
