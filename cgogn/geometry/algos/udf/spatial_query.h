#ifndef CGOGN_GEOMETRY_ALGOS_UDF_SPATIAL_QUERY_H_
#define CGOGN_GEOMETRY_ALGOS_UDF_SPATIAL_QUERY_H_

#include <cgogn/geometry/types/vector_traits.h>

#include <libacc/bvh_tree.h>
#include <libacc/kd_tree.h>

#include <utility>

namespace cgogn
{

namespace geometry
{

inline bool query_surface_closest_point(const acc::BVHTree<uint32, Vec3>& bvh, const Vec3& query,
											uint32& primitive_index, Vec3& closest_point, Scalar& distance)
{
	std::pair<uint32, Vec3> result;
	if (!bvh.closest_point(query, &result))
		return false;

	primitive_index = result.first;
	closest_point = result.second;
	distance = (query - result.second).norm();
	return true;
}

inline bool query_point_cloud_nearest_point(const acc::KDTree<3, uint32>& kdtree, const Vec3& query,
											uint32& point_index, Scalar& distance)
{
	std::pair<uint32, Scalar> result;
	if (!kdtree.find_nn(query, &result))
		return false;

	point_index = result.first;
	distance = result.second;
	return true;
}

} // namespace geometry

} // namespace cgogn

#endif // CGOGN_GEOMETRY_ALGOS_UDF_SPATIAL_QUERY_H_
