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

#ifndef CGOGN_GEOMETRY_ALGOS_SURFACE_SAMPLING_H_
#define CGOGN_GEOMETRY_ALGOS_SURFACE_SAMPLING_H_

#include <cgogn/core/functions/attributes.h>
#include <cgogn/core/functions/mesh_info.h>

#include <cgogn/geometry/types/vector_traits.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <random>
#include <vector>

namespace cgogn
{

namespace geometry
{

template <typename MESH>
std::vector<Vec3> sample_surface_area_weighted(
	const MESH& mesh, const typename mesh_traits<MESH>::template Attribute<Vec3>* vertex_position,
	const typename mesh_traits<MESH>::template Attribute<Scalar>* face_area, uint32 sample_count, uint32 random_seed)
{
	using Vertex = typename mesh_traits<MESH>::Vertex;
	using Face = typename mesh_traits<MESH>::Face;

	std::vector<std::array<Vec3, 3>> triangles;
	std::vector<Scalar> cumulative_area;
	Scalar total_area = 0.0;

	foreach_cell(mesh, [&](Face f) -> bool {
		const std::vector<Vertex> vertices = incident_vertices(mesh, f);
		const Scalar stored_face_area = value<Scalar>(mesh, face_area, f);
		for (uint32 i = 1; i + 1 < vertices.size(); ++i)
		{
			const Vec3& a = value<Vec3>(mesh, vertex_position, vertices[0]);
			const Vec3& b = value<Vec3>(mesh, vertex_position, vertices[i]);
			const Vec3& c = value<Vec3>(mesh, vertex_position, vertices[i + 1]);
			const Scalar triangle_area = Scalar(0.5) * (b - a).cross(c - a).norm();
			const Scalar area = vertices.size() == 3 ? stored_face_area : triangle_area;
			if (area > Scalar(0))
			{
				triangles.push_back({a, b, c});
				total_area += area;
				cumulative_area.push_back(total_area);
			}
		}
		return true;
	});

	std::vector<Vec3> samples;
	if (triangles.empty() || sample_count == 0)
		return samples;

	samples.reserve(sample_count);
	std::mt19937 generator(random_seed);
	std::uniform_real_distribution<Scalar> area_distribution(Scalar(0), total_area);
	std::uniform_real_distribution<Scalar> unit_distribution(Scalar(0), Scalar(1));

	for (uint32 i = 0; i < sample_count; ++i)
	{
		const Scalar target_area = area_distribution(generator);
		const auto triangle_it = std::lower_bound(cumulative_area.begin(), cumulative_area.end(), target_area);
		const uint32 triangle_index =
			std::min(uint32(triangle_it - cumulative_area.begin()), uint32(triangles.size() - 1));
		const Scalar sqrt_u = std::sqrt(unit_distribution(generator));
		const Scalar v = unit_distribution(generator);
		const Scalar w0 = Scalar(1) - sqrt_u;
		const Scalar w1 = sqrt_u * (Scalar(1) - v);
		const Scalar w2 = sqrt_u * v;
		const auto& triangle = triangles[triangle_index];
		samples.push_back(w0 * triangle[0] + w1 * triangle[1] + w2 * triangle[2]);
	}

	return samples;
}

} // namespace geometry

} // namespace cgogn

#endif // CGOGN_GEOMETRY_ALGOS_SURFACE_SAMPLING_H_
