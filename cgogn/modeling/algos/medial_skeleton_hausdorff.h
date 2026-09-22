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

#ifndef CGOGN_MODELING_ALGOS_MEDIAL_SKELETON_HAUSDORFF_H_
#define CGOGN_MODELING_ALGOS_MEDIAL_SKELETON_HAUSDORFF_H_

#include <cgogn/geometry/algos/surface_sampling.h>
#include <cgogn/modeling/skeleton_sampling.h>

#include <libacc/bvh_tree.h>

#include <algorithm>
#include <cmath>
#include <vector>

namespace cgogn
{

namespace modeling
{

using geometry::Scalar;
using geometry::Vec3;
using geometry::Vec4;

struct MedialSkeletonHausdorffOptions
{
	uint32 surface_sample_count = 10000;
	uint32 skeleton_sample_resolution = 50;
	uint32 random_seed = 1337;
};

struct MedialSkeletonHausdorffResult
{
	Scalar surface_to_skeleton = 0.0;
	Scalar skeleton_to_surface = 0.0;
	Scalar symmetric = 0.0;
	uint32 surface_sample_count = 0;
	uint32 skeleton_sample_count = 0;
};

namespace internal
{

inline Vec4 medial_sphere_vector(const Vec3& center, Scalar radius)
{
	return Vec4(center[0], center[1], center[2], radius);
}

} // namespace internal

template <typename SURFACE, typename SKELETON>
MedialSkeletonHausdorffResult compute_medial_skeleton_hausdorff(
	const SURFACE& surface,
	const typename mesh_traits<SURFACE>::template Attribute<Vec3>* surface_vertex_position,
	const typename mesh_traits<SURFACE>::template Attribute<Scalar>* surface_face_area, const SKELETON& skeleton,
	const typename mesh_traits<SKELETON>::template Attribute<Vec3>* skeleton_vertex_position,
	const typename mesh_traits<SKELETON>::template Attribute<Scalar>* skeleton_vertex_radius,
	acc::BVHTree<uint32, Vec3>* surface_bvh, const MedialSkeletonHausdorffOptions& options)
{
	using SkeletonVertex = typename mesh_traits<SKELETON>::Vertex;
	using SkeletonEdge = typename mesh_traits<SKELETON>::Edge;
	using SkeletonFace = typename mesh_traits<SKELETON>::Face;
	using Sampler = SkeletonSampler<Vec4, Vec3, Scalar>;

	Sampler sampler;
	foreach_cell(skeleton, [&](SkeletonVertex v) -> bool {
		sampler.add_vertex(value<Vec3>(skeleton, skeleton_vertex_position, v),
			value<Scalar>(skeleton, skeleton_vertex_radius, v));
		return true;
	});
	foreach_cell(skeleton, [&](SkeletonEdge e) -> bool {
		const std::vector<SkeletonVertex> vertices = incident_vertices(skeleton, e);
		if (vertices.size() == 2)
			sampler.add_edge(value<Vec3>(skeleton, skeleton_vertex_position, vertices[0]),
				value<Scalar>(skeleton, skeleton_vertex_radius, vertices[0]),
				value<Vec3>(skeleton, skeleton_vertex_position, vertices[1]),
				value<Scalar>(skeleton, skeleton_vertex_radius, vertices[1]));
		return true;
	});
	foreach_cell(skeleton, [&](SkeletonFace f) -> bool {
		const std::vector<SkeletonVertex> vertices = incident_vertices(skeleton, f);
		if (vertices.size() == 3)
		{
			const Vec3& c0 = value<Vec3>(skeleton, skeleton_vertex_position, vertices[0]);
			const Vec3& c1 = value<Vec3>(skeleton, skeleton_vertex_position, vertices[1]);
			const Vec3& c2 = value<Vec3>(skeleton, skeleton_vertex_position, vertices[2]);
			sampler.add_triangle(
				internal::medial_sphere_vector(c0, value<Scalar>(skeleton, skeleton_vertex_radius, vertices[0])),
				internal::medial_sphere_vector(c1, value<Scalar>(skeleton, skeleton_vertex_radius, vertices[1])),
				internal::medial_sphere_vector(c2, value<Scalar>(skeleton, skeleton_vertex_radius, vertices[2])));
		}
		return true;
	});

	const std::vector<Vec3> surface_samples = geometry::sample_surface_area_weighted(
		surface, surface_vertex_position, surface_face_area, options.surface_sample_count, options.random_seed);

	MedialSkeletonHausdorffResult result;
	result.surface_sample_count = uint32(surface_samples.size());
	for (const Vec3& point : surface_samples)
		result.surface_to_skeleton = std::max(result.surface_to_skeleton, std::abs(sampler.eval_skeleton(point)));

	const uint32 skeleton_sample_resolution = std::max(uint32(1), options.skeleton_sample_resolution);
	const Scalar skeleton_sample_step = sampler.BBwidth().norm() / Scalar(skeleton_sample_resolution);
	sampler.sample(skeleton_sample_step, Scalar(0), options.random_seed);

	const std::vector<Vec3> skeleton_samples = sampler.samples();
	result.skeleton_sample_count = uint32(skeleton_samples.size());
	for (const Vec3& point : skeleton_samples)
		result.skeleton_to_surface =
			std::max(result.skeleton_to_surface, (surface_bvh->closest_point(point) - point).norm());

	result.symmetric = std::max(result.surface_to_skeleton, result.skeleton_to_surface);
	return result;
}

} // namespace modeling

} // namespace cgogn

#endif // CGOGN_MODELING_ALGOS_MEDIAL_SKELETON_HAUSDORFF_H_
