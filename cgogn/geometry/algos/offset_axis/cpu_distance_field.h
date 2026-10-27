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
#ifndef CGOGN_GEOMETRY_ALGOS_OFFSET_AXIS_CPU_DISTANCE_FIELD_H_
#define CGOGN_GEOMETRY_ALGOS_OFFSET_AXIS_CPU_DISTANCE_FIELD_H_

#include <cgogn/geometry/algos/offset_axis/alpha_projection.h>
#include <cgogn/geometry/algos/offset_axis/distance_field.h>
#include <cgogn/geometry/algos/offset_axis/neural_model_type.h>
#include <cgogn/geometry/algos/offset_axis/normal_utils.h>
#include <cgogn/geometry/algos/offset_axis/parallel_for.h>
#include <cgogn/geometry/algos/offset_axis/ray_level_set_sampler.h>
#include <cgogn/geometry/types/vector_traits.h>

#include <atomic>
#include <cstdint>
#include <filesystem>
#include <memory>
#include <string>
#include <vector>

namespace cgogn
{

namespace geometry
{

namespace offset_axis
{

// Common part of the CPU distance fields (SurfaceDistanceField, PointCloudDistanceField): the candidate normals, the
// alpha projection and the topology scores all go through one closest-point query of the input, given by the derived
// field; there is no neural model and no MF value.
template <typename RAY_TRAITS>
class CpuDistanceField
{
public:
	using RayTraits = RAY_TRAITS;
	static constexpr bool is_neural_type = false;
	// the normals of the samples are recomputed from the input after the sampling filter
	static constexpr bool recompute_normals_from_input_after_filtering = true;

	// MF values are only defined for neural fields; this method is not called for the CPU fields
	bool values_sdf(const std::vector<Vec3>& /*points*/, std::vector<Scalar>& /*values*/, std::vector<Scalar>& /*sdf*/)
	{
		return false;
	}

	// the CPU fields never hold a neural model (see NeuralDistanceField)
	bool load_model(const std::filesystem::path& /*path*/, NeuralModelType /*type*/)
	{
		return false;
	}
	const std::string& load_error() const
	{
		static const std::string no_model = "this distance field holds no neural model";
		return no_model;
	}
	void reset()
	{
	}
	bool loaded() const
	{
		return false;
	}

protected:
	// normals of the candidate points from query(point, AlphaProjectionInfo&) -> bool
	template <typename QUERY>
	bool normals_from(const std::vector<Vec3>& points, std::vector<Vec3>& normals, std::vector<uint8_t>& keep,
					  const QUERY& query)
	{
		detail::parallel_for_index(points.size(), [&](std::size_t i) {
			AlphaProjectionInfo info;
			if (!query(points[i], info))
				return;
			normals[i] = unit_normal_or_z(info.normal);
			keep[i] = uint8_t(1);
		});
		return true;
	}

	// projection onto the alpha level set with the same query
	template <typename QUERY>
	bool project_with(const std::vector<Vec3>& points, const AlphaProjectionParameters& pp,
					  AlphaProjectionResult& projection, const QUERY& query)
	{
		return offset_axis::project_points_to_alpha(points, pp, query, projection);
	}

	// topology scores: closest(point, Scalar& distance) -> bool for every point; false if one point has no answer
	template <typename CLOSEST>
	bool distances_from(const std::vector<Vec3>& points, std::vector<Scalar>& values, const CLOSEST& closest)
	{
		values.resize(points.size());
		std::atomic<bool> all_found{true};
		detail::parallel_for_index(points.size(), [&](std::size_t i) {
			Scalar distance;
			if (!closest(points[i], distance))
			{
				all_found = false;
				return;
			}
			values[i] = distance;
		});
		if (!all_found)
			return false;
		return true;
	}
};

} // namespace offset_axis

} // namespace geometry

} // namespace cgogn

#endif // CGOGN_GEOMETRY_ALGOS_OFFSET_AXIS_CPU_DISTANCE_FIELD_H_
