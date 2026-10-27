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
#ifndef CGOGN_GEOMETRY_ALGOS_OFFSET_AXIS_INTERNAL_GREEDY_TET_REMOVAL_H_
#define CGOGN_GEOMETRY_ALGOS_OFFSET_AXIS_INTERNAL_GREEDY_TET_REMOVAL_H_

#include <cgogn/core/functions/attributes.h>
#include <cgogn/core/functions/mesh_info.h>
#include <cgogn/core/functions/traversals/edge.h>
#include <cgogn/core/functions/traversals/face.h>
#include <cgogn/core/functions/traversals/vertex.h>
#include <cgogn/core/types/cell_marker.h>
#include <cgogn/geometry/algos/offset_axis/internal/skeleton_mesh_ops.h>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <queue>
#include <set>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

namespace cgogn
{

namespace geometry
{

namespace offset_axis
{

namespace topology_detail
{

// Greedy removal of tets that have a simple face (a face owned by this tet only): the first face removed must be
// simple, and at most one face that becomes dangling is removed in the second step.
struct SimpleTetPolicy
{
	static bool accepts_tet(bool has_simple_face)
	{
		return has_simple_face;
	}
	static bool accepts_first_face(std::size_t face_tet_count)
	{
		return face_tet_count == 1;
	}
	static constexpr std::size_t max_second_step_faces = 1;
};

// Greedy removal of the remaining (non-simple) tets: any face of the tet may go first, and every face that becomes
// dangling is removed in the second step.
struct NonSimpleTetPolicy
{
	static bool accepts_tet(bool has_simple_face)
	{
		return !has_simple_face;
	}
	static bool accepts_first_face(std::size_t /*face_tet_count*/)
	{
		return true;
	}
	static constexpr std::size_t max_second_step_faces = std::numeric_limits<std::size_t>::max();
};

// Lazy greedy removal of skeleton tets for the topology fix. The tets are kept in a lazy priority queue ordered by
// edge score (descending), then face score (descending), then tet id (ascending). Each tet has a version number;
// bumping it invalidates the older queue entries, which are re-validated when popped. Each tet may lose at most two
// faces (face budget). A step removes one first face, then the faces that became dangling (how many is set by
// Policy::max_second_step_faces). The tets touched by the removed edges are re-queued.
template <typename POINTS, typename NONMANIFOLD, typename Policy>
class LazyGreedyTetRemover
{
public:
	using SphereVertex = typename mesh_traits<POINTS>::Vertex;
	using SkeletonVertex = typename mesh_traits<NONMANIFOLD>::Vertex;
	using SkeletonEdge = typename mesh_traits<NONMANIFOLD>::Edge;
	using SkeletonFace = typename mesh_traits<NONMANIFOLD>::Face;
	template <typename T>
	using SkeletonAttribute = typename mesh_traits<NONMANIFOLD>::template Attribute<T>;
	template <typename T>
	using SphereAttribute = typename mesh_traits<POINTS>::template Attribute<T>;

	LazyGreedyTetRemover(NONMANIFOLD& skeleton, SkeletonTets<NONMANIFOLD>& tets,
						 const SkeletonAttribute<Scalar>* face_score, const SkeletonAttribute<Scalar>* edge_score,
						 SpheresOptimizer<POINTS>& spheres_optimizer,
						 SkeletonAttribute<SphereVertex>* skeleton_source_sphere,
						 SphereAttribute<SkeletonVertex>* sphere_skeleton_vertex)
		: skeleton(skeleton), tets(tets), face_score(face_score), edge_score(edge_score),
		  spheres_optimizer(spheres_optimizer), skeleton_source_sphere(skeleton_source_sphere),
		  sphere_skeleton_vertex(sphere_skeleton_vertex)
	{
	}

	void run()
	{
		face_owner_tets.reserve(nb_cells<SkeletonFace>(skeleton));
		deleted_face_count_per_tet.reserve(tets.size());
		for (const auto& [tet_id, tet] : tets)
		{
			for (SkeletonFace f : tet.faces)
				face_owner_tets[index_of(skeleton, f)].push_back(tet_id);
			deleted_face_count_per_tet[tet_id] = 0;
		}

		tet_versions.reserve(tets.size() * 2 + 1);
		for (const auto& [tet_id, tet] : tets)
		{
			tet_versions[tet_id] = 0;
			queue_tet_candidate(tet_id, false);
		}

		while (true)
		{
			StepCandidate best;
			bool has_best = false;
			const Scalar edge_eps = Scalar(1e-12);
			const Scalar face_eps = Scalar(1e-12);
			while (!candidate_queue.empty())
			{
				const TetCandidateQueueEntry entry = candidate_queue.top();
				candidate_queue.pop();

				// the version of a removed tet is erased
				auto it_version = tet_versions.find(entry.tet_id);
				if (it_version == tet_versions.end() || it_version->second != entry.version)
					continue;

				StepCandidate refreshed;
				if (!build_step_candidate(entry.tet_id, refreshed))
				{
					++(it_version->second);
					continue;
				}

				const bool matches_cached = (refreshed.edge_id == entry.edge_id) &&
											(refreshed.first_face_id == entry.face_id) &&
											(std::abs(refreshed.edge_score - entry.edge_score) <= edge_eps) &&
											(std::abs(refreshed.first_face_score - entry.face_score) <= face_eps);
				if (!matches_cached)
				{
					++(it_version->second);
					candidate_queue.push({refreshed.tet_id, refreshed.edge_id, refreshed.first_face_id,
										  refreshed.edge_score, refreshed.first_face_score, it_version->second});
					continue;
				}

				best = refreshed;
				has_best = true;
				break;
			}

			if (!has_best)
				break;

			// its order is the order of the second step faces before their sort, which is not stable
			std::unordered_set<std::size_t> first_removed_tets;
			for (std::size_t tet_id : tets.tets_of(best.first_face_id))
				first_removed_tets.insert(tet_id);

			std::unordered_set<std::size_t> step_tets_to_erase;
			std::vector<SkeletonEdge> step_affected_edges;
			remove_face_collect_tets(best.first_face, step_tets_to_erase, step_affected_edges);

			std::vector<FaceChoice> created_deg1_faces;
			created_deg1_faces.reserve(16);
			CellMarkerStore<NONMANIFOLD, SkeletonFace> seen_created_faces(skeleton);
			for (std::size_t tet_id : first_removed_tets)
			{
				for (SkeletonFace f : tets.tet(tet_id).faces)
				{
					const uint32 idf = index_of(skeleton, f);
					if (idf == best.first_face_id)
						continue;
					if (seen_created_faces.is_marked(f))
						continue;
					seen_created_faces.mark(f);
					if (!has_face_budget(idf))
						continue;
					if (!face_has_degree1_edge(f))
						continue;
					if (face_has_edge_with_tet_face_count_gt2(f, step_tets_to_erase))
						continue;
					created_deg1_faces.push_back({f, (*face_score)[idf]});
				}
			}
			std::sort(created_deg1_faces.begin(), created_deg1_faces.end(),
					  [](const FaceChoice& a, const FaceChoice& b) { return a.score > b.score; });
			if (created_deg1_faces.size() > Policy::max_second_step_faces)
				created_deg1_faces.resize(Policy::max_second_step_faces);

			for (const FaceChoice& fc : created_deg1_faces)
			{
				if (!has_face_budget(index_of(skeleton, fc.face)))
					continue;
				if (!face_has_degree1_edge(fc.face))
					continue;
				if (face_has_edge_with_tet_face_count_gt2(fc.face, step_tets_to_erase))
					continue;
				remove_face_collect_tets(fc.face, step_tets_to_erase, step_affected_edges);
			}
			for (std::size_t tet_id : step_tets_to_erase)
			{
				tets.erase(tet_id);
				tet_versions.erase(tet_id);
			}

			for (std::size_t tet_id : dirty_tets(step_affected_edges))
				queue_tet_candidate(tet_id, true);
		}
	}

private:
	struct StepCandidate
	{
		std::size_t tet_id = 0;
		uint32 edge_id = INVALID_INDEX;
		Scalar edge_score = Scalar(0);
		SkeletonFace first_face;
		uint32 first_face_id = INVALID_INDEX;
		Scalar first_face_score = Scalar(0);
	};
	struct FaceChoice
	{
		SkeletonFace face;
		Scalar score = Scalar(0);
	};
	struct TetCandidateQueueEntry
	{
		std::size_t tet_id = 0;
		uint32 edge_id = INVALID_INDEX;
		uint32 face_id = INVALID_INDEX;
		Scalar edge_score = Scalar(0);
		Scalar face_score = Scalar(0);
		uint32 version = 0;
	};
	struct TetCandidateQueueCompare
	{
		bool operator()(const TetCandidateQueueEntry& a, const TetCandidateQueueEntry& b) const
		{
			if (a.edge_score != b.edge_score)
				return a.edge_score < b.edge_score;
			if (a.face_score != b.face_score)
				return a.face_score < b.face_score;
			return a.tet_id > b.tet_id;
		}
	};
	using TetCandidateQueue =
		std::priority_queue<TetCandidateQueueEntry, std::vector<TetCandidateQueueEntry>, TetCandidateQueueCompare>;

	NONMANIFOLD& skeleton;
	SkeletonTets<NONMANIFOLD>& tets;
	const SkeletonAttribute<Scalar>* face_score;
	const SkeletonAttribute<Scalar>* edge_score;
	SpheresOptimizer<POINTS>& spheres_optimizer;
	SkeletonAttribute<SphereVertex>* skeleton_source_sphere;
	SphereAttribute<SkeletonVertex>* sphere_skeleton_vertex;
	std::unordered_map<uint32, std::vector<std::size_t>> face_owner_tets;
	std::unordered_map<std::size_t, uint32> deleted_face_count_per_tet;
	std::unordered_map<std::size_t, uint32> tet_versions;
	TetCandidateQueue candidate_queue;

	bool has_face_budget(uint32 face_id)
	{
		auto it_owner = face_owner_tets.find(face_id);
		if (it_owner == face_owner_tets.end())
			return false;
		for (std::size_t tet_id : it_owner->second)
		{
			if (!tets.contains(tet_id))
				continue;
			auto it_count = deleted_face_count_per_tet.find(tet_id);
			if (it_count != deleted_face_count_per_tet.end() && it_count->second >= 2)
				return false;
		}
		return true;
	}

	void account_face_deletion_budget(uint32 face_id)
	{
		auto it_owner = face_owner_tets.find(face_id);
		if (it_owner == face_owner_tets.end())
			return;
		for (std::size_t tet_id : it_owner->second)
		{
			if (!tets.contains(tet_id))
				continue;
			auto it_count = deleted_face_count_per_tet.find(tet_id);
			if (it_count != deleted_face_count_per_tet.end())
				++(it_count->second);
		}
	}

	bool face_belongs_to_tet(SkeletonFace f, std::size_t tet_id)
	{
		const std::set<std::size_t>& in_tets = tets.tets_of(index_of(skeleton, f));
		return in_tets.find(tet_id) != in_tets.end();
	}

	bool face_has_degree1_edge(SkeletonFace f)
	{
		for (SkeletonEdge e : incident_edges(skeleton, f))
			if (incident_faces(skeleton, e).size() == 1)
				return true;
		return false;
	}

	// True if an edge of the face has more than 2 faces that bound a current tet not in ignored_tets.
	bool face_has_edge_with_tet_face_count_gt2(SkeletonFace f, const std::unordered_set<std::size_t>& ignored_tets)
	{
		for (SkeletonEdge e : incident_edges(skeleton, f))
		{
			size_t tet_face_count = 0;
			for (SkeletonFace ef : incident_faces(skeleton, e))
			{
				for (std::size_t tet_id : tets.tets_of(index_of(skeleton, ef)))
				{
					if (ignored_tets.find(tet_id) == ignored_tets.end())
					{
						++tet_face_count;
						break;
					}
				}
			}
			if (tet_face_count > 2)
				return true;
		}
		return false;
	}

	bool tet_has_simple_face(std::size_t tet_id)
	{
		for (SkeletonFace f : tets.tet(tet_id).faces)
			if (tets.tets_of(index_of(skeleton, f)).size() == 1)
				return true;
		return false;
	}

	bool pick_tet_best_edge_and_face(std::size_t tet_id, uint32& out_edge_id, Scalar& out_edge_score,
									 SkeletonFace& out_face, uint32& out_face_id, Scalar& out_face_score)
	{
		CellMarkerStore<NONMANIFOLD, SkeletonEdge> seen_edges(skeleton);
		bool found = false;
		const Scalar edge_eps = Scalar(1e-12);
		const Scalar face_eps = Scalar(1e-12);
		for (SkeletonFace f : tets.tet(tet_id).faces)
		{
			for (SkeletonEdge e : incident_edges(skeleton, f))
			{
				if (seen_edges.is_marked(e))
					continue;
				seen_edges.mark(e);
				const uint32 ide = index_of(skeleton, e);
				SkeletonFace best_face_on_edge;
				uint32 best_face_on_edge_id = INVALID_INDEX;
				Scalar best_face_on_edge_score = Scalar(0);
				bool edge_has_candidate = false;
				for (SkeletonFace ef : incident_faces(skeleton, e))
				{
					if (!face_belongs_to_tet(ef, tet_id))
						continue;
					const uint32 idf = index_of(skeleton, ef);
					if (!has_face_budget(idf))
						continue;
					if (!Policy::accepts_first_face(tets.tets_of(idf).size()))
						continue;
					const Scalar score = (*face_score)[idf];
					if (!edge_has_candidate || score > best_face_on_edge_score + face_eps)
					{
						edge_has_candidate = true;
						best_face_on_edge = ef;
						best_face_on_edge_id = idf;
						best_face_on_edge_score = score;
					}
				}
				if (!edge_has_candidate)
					continue;
				const Scalar score = (*edge_score)[ide];
				const bool better_edge = (!found) || (score > out_edge_score + edge_eps);
				const bool tie_edge = found && (std::abs(score - out_edge_score) <= edge_eps);
				const bool better_face_on_tie = tie_edge && (best_face_on_edge_score > out_face_score + face_eps);
				if (better_edge || better_face_on_tie)
				{
					found = true;
					out_edge_id = ide;
					out_edge_score = score;
					out_face = best_face_on_edge;
					out_face_id = best_face_on_edge_id;
					out_face_score = best_face_on_edge_score;
				}
			}
		}
		return found;
	}

	// Removes the face, a face of a current tet, and appends its edges to affected_edges. Its tets are added to
	// tets_to_erase: they are erased at the end of the step.
	void remove_face_collect_tets(SkeletonFace f, std::unordered_set<std::size_t>& tets_to_erase,
								  std::vector<SkeletonEdge>& affected_edges)
	{
		const uint32 idf = index_of(skeleton, f);
		const std::set<std::size_t>& face_tets = tets.tets_of(idf);
		tets_to_erase.insert(face_tets.begin(), face_tets.end());
		account_face_deletion_budget(idf);
		for (SkeletonEdge e : incident_edges(skeleton, f))
			affected_edges.push_back(e);
		remove_faces(skeleton, {f}, spheres_optimizer, skeleton_source_sphere, sphere_skeleton_vertex);
	}

	bool build_step_candidate(std::size_t tet_id, StepCandidate& out_candidate)
	{
		if (!Policy::accepts_tet(tet_has_simple_face(tet_id)))
			return false;
		out_candidate.tet_id = tet_id;
		return pick_tet_best_edge_and_face(tet_id, out_candidate.edge_id, out_candidate.edge_score,
										   out_candidate.first_face, out_candidate.first_face_id,
										   out_candidate.first_face_score);
	}

	void queue_tet_candidate(std::size_t tet_id, bool bump_version)
	{
		uint32& version = tet_versions[tet_id];
		if (bump_version)
			++version;
		StepCandidate candidate;
		if (!build_step_candidate(tet_id, candidate))
			return;
		candidate_queue.push({candidate.tet_id, candidate.edge_id, candidate.first_face_id, candidate.edge_score,
							  candidate.first_face_score, version});
	}

	// The current tets of the faces of the edges.
	std::unordered_set<std::size_t> dirty_tets(const std::vector<SkeletonEdge>& affected_edges)
	{
		std::unordered_set<std::size_t> result;
		for (SkeletonEdge e : affected_edges)
			for (SkeletonFace f : incident_faces(skeleton, e))
				for (std::size_t tet_id : tets.tets_of(index_of(skeleton, f)))
					result.insert(tet_id);
		return result;
	}
};

// Steps b and c of the topology fix loop: greedy removal of the tets accepted by Policy (SimpleTetPolicy or
// NonSimpleTetPolicy), the best edge and face first by score.
template <typename Policy, typename POINTS, typename NONMANIFOLD>
void remove_tets_greedy(
	NONMANIFOLD& skeleton, SkeletonTets<NONMANIFOLD>& tets,
	const typename mesh_traits<NONMANIFOLD>::template Attribute<Scalar>* face_score,
	const typename mesh_traits<NONMANIFOLD>::template Attribute<Scalar>* edge_score,
	SpheresOptimizer<POINTS>& spheres_optimizer,
	typename mesh_traits<NONMANIFOLD>::template Attribute<typename mesh_traits<POINTS>::Vertex>* skeleton_source_sphere,
	typename mesh_traits<POINTS>::template Attribute<typename mesh_traits<NONMANIFOLD>::Vertex>* sphere_skeleton_vertex)
{
	LazyGreedyTetRemover<POINTS, NONMANIFOLD, Policy>(skeleton, tets, face_score, edge_score, spheres_optimizer,
													  skeleton_source_sphere, sphere_skeleton_vertex)
		.run();
}

} // namespace topology_detail

} // namespace offset_axis

} // namespace geometry

} // namespace cgogn

#endif // CGOGN_GEOMETRY_ALGOS_OFFSET_AXIS_INTERNAL_GREEDY_TET_REMOVAL_H_
