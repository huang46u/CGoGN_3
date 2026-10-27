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
#ifndef CGOGN_GEOMETRY_ALGOS_OFFSET_AXIS_ALPHA_MEDIAL_H_
#define CGOGN_GEOMETRY_ALGOS_OFFSET_AXIS_ALPHA_MEDIAL_H_

#include <cgogn/core/functions/attributes.h>
#include <cgogn/core/functions/mesh_info.h>
#include <cgogn/core/types/mesh_traits.h>
#include <cgogn/core/utils/assert.h>
#include <cgogn/geometry/algos/offset_axis/distance_field.h>
#include <cgogn/geometry/algos/offset_axis/neural_model_type.h>
#include <cgogn/geometry/algos/offset_axis/sample_processing.h>
#include <cgogn/geometry/types/vector_traits.h>

#include <libacc/kd_tree.h>

#include <algorithm>
#include <memory>
#include <vector>

namespace cgogn
{

namespace geometry
{

namespace offset_axis
{

namespace constants
{

// MAFlipPrune, neural UDF: a medial centre is retried when its UDF value exceeds this factor times alpha.
inline constexpr Scalar udf_flip_threshold_factor = Scalar(1.2);

} // namespace constants

// Medial flip prune. scores(points, values, sdf) -> bool scores the medial points (sdf is filled for MF models
// only): a medial point whose score is above threshold (or, for MF models, whose sdf is positive) is on the wrong
// side, so the normal of its sample is flipped and its shrinking ball recomputed. The samples whose new medial
// point still scores above threshold are removed; if some were, the kd-tree is rebuilt and the flipped samples
// that remain get a PCA normal and a new shrinking ball, and are pruned again.
template <typename POINTS, typename SCORES>
void flip_prune_medial_axis(POINTS& samples, const typename mesh_traits<POINTS>::template Attribute<Vec3>& position,
							typename mesh_traits<POINTS>::template Attribute<Vec3>& normal,
							std::unique_ptr<acc::KDTree<3, uint32>>& kdtree,
							std::vector<typename mesh_traits<POINTS>::Vertex>& kdtree_vertices,
							typename mesh_traits<POINTS>::template Attribute<Vec3>& ma_position,
							typename mesh_traits<POINTS>::template Attribute<Scalar>& ma_radius,
							typename mesh_traits<POINTS>::template Attribute<typename mesh_traits<POINTS>::Vertex>&
								ma_secondary,
							bool mf_model, Scalar threshold, Scalar alpha, int knn_k, SCORES&& scores)
{
	using Vertex = typename mesh_traits<POINTS>::Vertex;

	const Scalar initial_radius = alpha * Scalar(10);
	auto shrink = [&](const std::vector<Vertex>& vertices) {
		shrinking_ball_centers(samples, position, normal, *kdtree, kdtree_vertices, initial_radius, alpha, vertices,
							   ma_position, ma_radius, ma_secondary);
	};

	std::vector<Vertex> vertices;
	std::vector<Vec3> centers;
	vertices.reserve(nb_cells<Vertex>(samples));
	centers.reserve(vertices.capacity());
	foreach_cell(samples, [&](Vertex v) {
		const uint32 i = index_of(samples, v);
		vertices.push_back(v);
		centers.push_back(ma_position[i]);
		return true;
	});
	std::vector<Scalar> values, sdf;
	if (!scores(centers, values, sdf))
		return;
	cgogn_message_assert(values.size() == centers.size() && (!mf_model || sdf.size() == centers.size()),
						 "scores output size mismatch");

	std::vector<Vertex> retry;
	for (size_t i = 0; i < vertices.size(); ++i)
		if ((mf_model && sdf[i] > Scalar(0)) || values[i] > threshold)
			retry.push_back(vertices[i]);
	std::vector<uint32> flipped_indices;
	flipped_indices.reserve(retry.size());
	for (Vertex v : retry)
	{
		const uint32 i = index_of(samples, v);
		normal[i] = -normal[i].normalized();
		flipped_indices.push_back(i);
	}
	shrink(retry);

	auto prune = [&](const std::vector<Vertex>& candidates) {
		std::vector<Vec3> p;
		p.reserve(candidates.size());
		for (Vertex v : candidates)
			p.push_back(ma_position[index_of(samples, v)]);
		std::vector<Scalar> val, sd;
		if (!scores(p, val, sd))
			return uint32(0);
		cgogn_message_assert(val.size() == p.size() && (!mf_model || sd.size() == p.size()),
							 "scores output size mismatch");
		uint32 removed = 0;
		for (size_t i = 0; i < candidates.size(); ++i)
			if ((mf_model && sd[i] > Scalar(0)) || val[i] > threshold)
			{
				remove_vertex(samples, candidates[i]);
				++removed;
			}
		return removed;
	};
	const uint32 deleted = prune(retry);
	if (deleted > 0 && !flipped_indices.empty() && nb_cells<Vertex>(samples) > 0)
	{
		build_point_kdtree(samples, position, kdtree, kdtree_vertices);
		std::vector<Vertex> survivors;
		survivors.reserve(flipped_indices.size());
		for (uint32 i : flipped_indices)
		{
			Vertex v = of_index<Vertex>(samples, i);
			if (v.is_valid())
				survivors.push_back(v);
		}
		for (Vertex v : survivors)
			detail::recompute_normal_pca(samples, position, normal, *kdtree, kdtree_vertices, knn_k, v);
		shrink(survivors);
		prune(survivors);
	}
}

// The shrinking ball medial point of every sample, then (ma_flip_prune) the medial flip prune. Its threshold is
// ma_flip_prune_alpha_factor * alpha on the network value for MF models, udf_flip_threshold_factor * alpha on the UDF
// for UDF models and alpha on the distance for the CPU fields. The caller rebuilds the kd-trees of the samples.
template <typename FIELD, typename POINTS>
void compute_initial_medial_axis(
	FIELD& field, POINTS& samples, const typename mesh_traits<POINTS>::template Attribute<Vec3>& position,
	typename mesh_traits<POINTS>::template Attribute<Vec3>& normal, std::unique_ptr<acc::KDTree<3, uint32>>& kdtree,
	std::vector<typename mesh_traits<POINTS>::Vertex>& kdtree_vertices,
	typename mesh_traits<POINTS>::template Attribute<Vec3>& ma_position,
	typename mesh_traits<POINTS>::template Attribute<Scalar>& ma_radius,
	typename mesh_traits<POINTS>::template Attribute<typename mesh_traits<POINTS>::Vertex>& ma_secondary, Scalar alpha,
	bool ma_flip_prune, Scalar ma_flip_prune_alpha_factor, int knn_k)
{
	const bool mf_model = is_neural_model(field, NeuralModelType::mf);
	const bool udf_model = is_neural_model(field, NeuralModelType::udf);
	shrinking_ball_centers(samples, position, normal, *kdtree, kdtree_vertices, alpha * Scalar(10), alpha, ma_position,
						   ma_radius, ma_secondary);
	if (!ma_flip_prune)
		return;
	const Scalar threshold = mf_model	 ? alpha * ma_flip_prune_alpha_factor
							 : udf_model ? constants::udf_flip_threshold_factor * alpha
										 : alpha;
	auto scores = [&](const std::vector<Vec3>& points, std::vector<Scalar>& values, std::vector<Scalar>& sdf) -> bool {
		sdf.clear();
		if (!mf_model)
			return field.distances(points, values);
		return field.values_sdf(points, values, sdf);
	};
	flip_prune_medial_axis(samples, position, normal, kdtree, kdtree_vertices, ma_position, ma_radius, ma_secondary,
						   mf_model, threshold, alpha, knn_k, scores);
}

} // namespace offset_axis

} // namespace geometry

} // namespace cgogn

#endif // CGOGN_GEOMETRY_ALGOS_OFFSET_AXIS_ALPHA_MEDIAL_H_
