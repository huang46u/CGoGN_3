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
#ifndef CGOGN_GEOMETRY_ALGOS_UDF_ALPHA_PROJECTION_H_
#define CGOGN_GEOMETRY_ALGOS_UDF_ALPHA_PROJECTION_H_

#include <cgogn/geometry/types/vector_traits.h>

#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <vector>

namespace cgogn
{

namespace geometry
{

struct AlphaProjectionParameters
{
	Scalar alpha = Scalar(0);
	Scalar tolerance = Scalar(1e-6);
	Vec3 bbox_min = Vec3(0, 0, 0);
	Vec3 bbox_max = Vec3(1, 1, 1);
	int max_iterations = 1;
};

struct AlphaProjectionResult
{
	std::vector<Vec3> projected_points;
	std::vector<Vec3> normals;
	std::vector<uint8_t> keep_mask;
};

struct AlphaProjectionInfo
{
	Vec3 closest_point = Vec3(0, 0, 0);
	Vec3 normal = Vec3(0, 0, 1);
	Scalar distance = std::numeric_limits<Scalar>::max();
	bool valid = false;
};

template <typename QueryFn>
bool project_points_to_alpha(const std::vector<Vec3>& points, const AlphaProjectionParameters& parameters,
							 QueryFn&& query, AlphaProjectionResult& result)
{
	const size_t nb_points = points.size();
	result.projected_points.assign(nb_points, Vec3(0, 0, 0));
	result.normals.assign(nb_points, Vec3(0, 0, 1));
	result.keep_mask.assign(nb_points, uint8_t(0));
	if (nb_points == 0)
		return true;

	for (size_t i = 0; i < nb_points; ++i)
	{
		Vec3 x = points[i].cwiseMax(parameters.bbox_min).cwiseMin(parameters.bbox_max);
		AlphaProjectionInfo info;
		bool valid = false;
		for (int iter = 0; iter < parameters.max_iterations; ++iter)
		{
			if (!query(x, info) || !info.valid)
			{
				valid = false;
				break;
			}
			valid = true;
			const Scalar residual = info.distance - parameters.alpha;
			if (std::abs(residual) <= parameters.tolerance)
				break;
			if (info.normal.squaredNorm() < Scalar(1e-12) || !info.normal.allFinite())
			{
				valid = false;
				break;
			}
			x -= residual * info.normal;
			x = x.cwiseMax(parameters.bbox_min).cwiseMin(parameters.bbox_max);
		}

		if (!valid || !query(x, info) || !info.valid)
			continue;

		const Scalar residual = std::abs(info.distance - parameters.alpha);
		if (residual > parameters.tolerance)
			continue;

		result.projected_points[i] = info.closest_point + info.normal * parameters.alpha;
		result.normals[i] = info.normal;
		result.keep_mask[i] = uint8_t(1);
	}
	return true;
}

} // namespace geometry

} // namespace cgogn

#endif // CGOGN_GEOMETRY_ALGOS_UDF_ALPHA_PROJECTION_H_
