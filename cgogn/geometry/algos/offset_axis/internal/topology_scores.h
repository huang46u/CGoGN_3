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
#ifndef CGOGN_GEOMETRY_ALGOS_OFFSET_AXIS_INTERNAL_TOPOLOGY_SCORES_H_
#define CGOGN_GEOMETRY_ALGOS_OFFSET_AXIS_INTERNAL_TOPOLOGY_SCORES_H_

#include <cgogn/core/functions/attributes.h>
#include <cgogn/core/functions/mesh_info.h>
#include <cgogn/core/functions/traversals/edge.h>
#include <cgogn/core/functions/traversals/face.h>
#include <cgogn/core/functions/traversals/vertex.h>
#include <cgogn/core/types/cell_marker.h>
#include <cgogn/geometry/algos/area.h>
#include <cgogn/geometry/algos/offset_axis/internal/skeleton_topology_types.h>

#include <algorithm>
#include <array>
#include <limits>
#include <unordered_set>
#include <vector>

namespace cgogn
{

namespace geometry
{

namespace offset_axis
{

namespace topology_detail
{

template <typename POINTS, typename NONMANIFOLD>
struct TopologyScores
{
	using Types = SkeletonTopologyTypes<POINTS, NONMANIFOLD>;
	using SkeletonVertex = typename Types::SkeletonVertex;
	using SkeletonEdge = typename Types::SkeletonEdge;
	using SkeletonFace = typename Types::SkeletonFace;
	using Context = typename Types::Context;
	using Tet = typename Types::Tet;

	// Sets the face_score of each face of a current tet to the UDF integral over the face (adaptive Dunavant-7
	// quadrature), divided by the face area if normalize_by_area. The other faces keep their score.
	// Returns false if the scores cannot be evaluated.
	template <typename ScoreEvaluator>
	static bool compute_face_scores(Context& p, ScoreEvaluator& evaluator, bool normalize_by_area)
	{
		const std::unordered_set<uint32> face_ids = collect_current_tet_face_ids(p);
		if (face_ids.empty())
			return true;
		auto& face_scores = *p.data.face_score;

		struct AdaptiveTriangleTask
		{
			uint32 face_id = INVALID_INDEX;
			Vec3 a;
			Vec3 b;
			Vec3 c;
			uint32 depth = 0;
		};

		static const std::array<std::array<Scalar, 3>, 7> dunavant7_bary = {{
			{{Scalar(1.0 / 3.0), Scalar(1.0 / 3.0), Scalar(1.0 / 3.0)}},
			{{Scalar(0.470142064105115), Scalar(0.470142064105115), Scalar(0.059715871789770)}},
			{{Scalar(0.470142064105115), Scalar(0.059715871789770), Scalar(0.470142064105115)}},
			{{Scalar(0.059715871789770), Scalar(0.470142064105115), Scalar(0.470142064105115)}},
			{{Scalar(0.101286507323456), Scalar(0.101286507323456), Scalar(0.797426985353087)}},
			{{Scalar(0.101286507323456), Scalar(0.797426985353087), Scalar(0.101286507323456)}},
			{{Scalar(0.797426985353087), Scalar(0.101286507323456), Scalar(0.101286507323456)}},
		}};
		static const std::array<Scalar, 7> dunavant7_w = {
			Scalar(0.225000000000000), Scalar(0.132394152788506), Scalar(0.132394152788506), Scalar(0.132394152788506),
			Scalar(0.125939180544827), Scalar(0.125939180544827), Scalar(0.125939180544827)};

		auto subdivide_triangle = [&](const AdaptiveTriangleTask& tri, std::vector<AdaptiveTriangleTask>& out) {
			const Vec3 ab = (tri.a + tri.b) * Scalar(0.5);
			const Vec3 bc = (tri.b + tri.c) * Scalar(0.5);
			const Vec3 ca = (tri.c + tri.a) * Scalar(0.5);
			const uint32 next_depth = tri.depth + 1;
			out.push_back({tri.face_id, tri.a, ab, ca, next_depth});
			out.push_back({tri.face_id, ab, tri.b, bc, next_depth});
			out.push_back({tri.face_id, ca, bc, tri.c, next_depth});
			out.push_back({tri.face_id, ab, bc, ca, next_depth});
		};

		std::vector<AdaptiveTriangleTask> pending_tris;
		pending_tris.reserve(face_ids.size() * 2);
		foreach_cell(p.data.skeleton, [&](SkeletonFace f) -> bool {
			const uint32 idf = index_of(p.data.skeleton, f);
			if (idf == INVALID_INDEX)
				return true;
			if (face_ids.find(idf) == face_ids.end())
				return true;
			face_scores[idf] = Scalar(0);
			if (!f.is_valid())
				return true;
			const std::vector<SkeletonVertex> vertices = incident_vertices(p.data.skeleton, f);
			if (vertices.size() < 3)
				return true;
			const Vec3 p0 = value<Vec3>(p.data.skeleton, p.data.skeleton_position, vertices[0]);
			for (uint32 i = 1; i + 1 < vertices.size(); ++i)
			{
				const Vec3 p1 = value<Vec3>(p.data.skeleton, p.data.skeleton_position, vertices[i]);
				const Vec3 p2 = value<Vec3>(p.data.skeleton, p.data.skeleton_position, vertices[i + 1]);
				if (!(geometry::area(p0, p1, p2) > Scalar(0)))
					continue;
				pending_tris.push_back({idf, p0, p1, p2, 0});
			}
			return true;
		});

		while (!pending_tris.empty())
		{
			std::vector<AdaptiveTriangleTask> next_tris;
			next_tris.reserve(pending_tris.size() * 2);

			struct TriBatchMeta
			{
				size_t tri_idx = 0;
				size_t sample_offset = 0;
				Scalar area = Scalar(0);
			};
			std::vector<TriBatchMeta> tri_batch;
			tri_batch.reserve(pending_tris.size());

			std::vector<Vec3> sample_points;
			sample_points.reserve(pending_tris.size() * 7);

			for (size_t tri_idx = 0; tri_idx < pending_tris.size(); ++tri_idx)
			{
				const AdaptiveTriangleTask& tri = pending_tris[tri_idx];
				const Scalar tri_area = geometry::area(tri.a, tri.b, tri.c);
				if (!(tri_area > Scalar(0)))
					continue;
				const size_t sample_offset = sample_points.size();
				for (uint32 k = 0; k < 7; ++k)
				{
					const Scalar l0 = dunavant7_bary[k][0];
					const Scalar l1 = dunavant7_bary[k][1];
					const Scalar l2 = dunavant7_bary[k][2];
					sample_points.push_back(tri.a * l0 + tri.b * l1 + tri.c * l2);
				}
				tri_batch.push_back({tri_idx, sample_offset, tri_area});
			}

			if (sample_points.empty())
			{
				pending_tris.clear();
				break;
			}

			std::vector<Scalar> score_values;
			if (!evaluator(sample_points, score_values))
				return false;

			for (const TriBatchMeta& tri_eval : tri_batch)
			{
				const AdaptiveTriangleTask& tri = pending_tris[tri_eval.tri_idx];
				Scalar weighted_mean = Scalar(0);
				Scalar min_udf = std::numeric_limits<Scalar>::max();
				Scalar max_udf = std::numeric_limits<Scalar>::lowest();
				for (uint32 k = 0; k < 7; ++k)
				{
					const Scalar u = score_values[tri_eval.sample_offset + k];
					weighted_mean += dunavant7_w[k] * u;
					min_udf = std::min(min_udf, u);
					max_udf = std::max(max_udf, u);
				}

				const Scalar tri_integral = tri_eval.area * weighted_mean;
				const Scalar udf_range = max_udf - min_udf;
				const Scalar refine_threshold =
					std::max(Scalar(1e-4), Scalar(0.35) * std::max(weighted_mean, Scalar(0)));
				const bool should_refine = (tri.depth < 2) && (udf_range > refine_threshold);

				if (should_refine)
				{
					subdivide_triangle(tri, next_tris);
				}
				else
				{
					face_scores[tri.face_id] += tri_integral;
				}
			}
			pending_tris.swap(next_tris);
		}

		if (normalize_by_area)
		{
			foreach_cell(p.data.skeleton, [&](SkeletonFace f) -> bool {
				const uint32 face_id = index_of(p.data.skeleton, f);
				if (face_id == INVALID_INDEX || face_ids.find(face_id) == face_ids.end())
					return true;
				const std::vector<SkeletonVertex> vertices = incident_vertices(p.data.skeleton, f);
				if (vertices.size() < 3)
				{
					face_scores[face_id] = Scalar(0);
					return true;
				}
				const Vec3 p0 = value<Vec3>(p.data.skeleton, p.data.skeleton_position, vertices[0]);
				Scalar area = Scalar(0);
				for (uint32 i = 1; i + 1 < vertices.size(); ++i)
				{
					const Vec3 p1 = value<Vec3>(p.data.skeleton, p.data.skeleton_position, vertices[i]);
					const Vec3 p2 = value<Vec3>(p.data.skeleton, p.data.skeleton_position, vertices[i + 1]);
					area += geometry::area(p0, p1, p2);
				}
				face_scores[face_id] = (area > Scalar(0)) ? (face_scores[face_id] / area) : Scalar(0);
				return true;
			});
		}
		return true;
	}

	static std::unordered_set<uint32> collect_current_tet_edge_ids(Context& p)
	{
		std::unordered_set<uint32> tet_edge_ids;
		tet_edge_ids.reserve(nb_cells<SkeletonEdge>(p.data.skeleton));
		for (const auto& kv : p.skeleton_tets)
		{
			const Tet& tet = kv.second;
			for (uint32 i = 0; i < 4; ++i)
			{
				const SkeletonFace f = tet.faces[i];
				if (!f.is_valid())
					continue;
				for (SkeletonEdge e : incident_edges(p.data.skeleton, f))
				{
					if (!e.is_valid())
						continue;
					const uint32 ide = index_of(p.data.skeleton, e);
					if (ide != INVALID_INDEX)
						tet_edge_ids.insert(ide);
				}
			}
		}
		return tet_edge_ids;
	}

	// Sets the edge_score of each edge of a current tet to the mean UDF along the edge (3-point Gauss quadrature).
	// The other edges keep their score. Returns false if the scores cannot be evaluated.
	template <typename ScoreEvaluator>
	static bool compute_edge_scores(Context& p, ScoreEvaluator& evaluator)
	{
		const std::unordered_set<uint32> edge_ids = collect_current_tet_edge_ids(p);
		if (edge_ids.empty())
			return true;
		auto& edge_scores = *p.data.edge_score;

		static const Scalar gauss3_xi[3] = {Scalar(-0.7745966692414834), Scalar(0.0), Scalar(0.7745966692414834)};
		static const Scalar gauss3_w_ref[3] = {Scalar(0.5555555555555556), Scalar(0.8888888888888888),
											   Scalar(0.5555555555555556)};
		static const Scalar gauss3_w_01[3] = {Scalar(0.5) * gauss3_w_ref[0], Scalar(0.5) * gauss3_w_ref[1],
											  Scalar(0.5) * gauss3_w_ref[2]};

		struct EdgeBatchMeta
		{
			uint32 edge_id = INVALID_INDEX;
			size_t sample_offset = 0;
		};

		std::vector<EdgeBatchMeta> edge_batch;
		edge_batch.reserve(edge_ids.size());
		std::vector<Vec3> sample_points;
		sample_points.reserve(edge_ids.size() * 3);

		foreach_cell(p.data.skeleton, [&](SkeletonEdge e) -> bool {
			if (!e.is_valid())
				return true;
			const uint32 ide = index_of(p.data.skeleton, e);
			if (ide == INVALID_INDEX)
				return true;
			if (edge_ids.find(ide) == edge_ids.end())
				return true;
			edge_scores[ide] = Scalar(0);

			const std::vector<SkeletonVertex> verts = incident_vertices(p.data.skeleton, e);
			if (verts.size() != 2)
				return true;

			const Vec3 p0 = value<Vec3>(p.data.skeleton, p.data.skeleton_position, verts[0]);
			const Vec3 p1 = value<Vec3>(p.data.skeleton, p.data.skeleton_position, verts[1]);
			const Scalar edge_len = (p1 - p0).norm();
			if (!(edge_len > Scalar(0)))
				return true;

			const size_t sample_offset = sample_points.size();
			for (uint32 k = 0; k < 3; ++k)
			{
				const Scalar t = Scalar(0.5) * (gauss3_xi[k] + Scalar(1.0));
				sample_points.push_back(p0 * (Scalar(1.0) - t) + p1 * t);
			}
			edge_batch.push_back({ide, sample_offset});
			return true;
		});

		if (sample_points.empty())
			return true;

		std::vector<Scalar> score_values;
		if (!evaluator(sample_points, score_values))
			return false;

		for (const EdgeBatchMeta& edge_eval : edge_batch)
		{
			Scalar avg_udf = Scalar(0);
			for (uint32 k = 0; k < 3; ++k)
				avg_udf += gauss3_w_01[k] * score_values[edge_eval.sample_offset + k];
			edge_scores[edge_eval.edge_id] = avg_udf;
		}
		return true;
	}

	static std::unordered_set<uint32> collect_current_tet_face_ids(Context& p)
	{
		std::unordered_set<uint32> tet_face_ids;
		tet_face_ids.reserve(nb_cells<SkeletonFace>(p.data.skeleton));
		for (const auto& kv : p.skeleton_tets)
		{
			const Tet& tet = kv.second;
			for (uint32 i = 0; i < 4; ++i)
			{
				const SkeletonFace f = tet.faces[i];
				if (!f.is_valid())
					continue;
				const uint32 idf = index_of(p.data.skeleton, f);
				if (idf != INVALID_INDEX)
					tet_face_ids.insert(idf);
			}
		}
		return tet_face_ids;
	}
};

} // namespace topology_detail

} // namespace offset_axis

} // namespace geometry

} // namespace cgogn

#endif // CGOGN_GEOMETRY_ALGOS_OFFSET_AXIS_INTERNAL_TOPOLOGY_SCORES_H_
