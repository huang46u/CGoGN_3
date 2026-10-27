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
#include <cgogn/geometry/algos/offset_axis/internal/quadrature_scores.h>
#include <cgogn/geometry/algos/offset_axis/internal/skeleton_mesh_ops.h>
#include <cgogn/geometry/algos/offset_axis/internal/skeleton_topology_types.h>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <queue>
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
template <typename POINTS, typename NONMANIFOLD, typename Policy, typename ScoreEvaluator>
class LazyGreedyTetRemover
{
public:
	using Types = SkeletonTopologyTypes<POINTS, NONMANIFOLD>;
	using SkeletonEdge = typename Types::SkeletonEdge;
	using SkeletonFace = typename Types::SkeletonFace;
	using Context = typename Types::Context;
	using Tet = typename Types::Tet;
	using RemovalStats = typename Types::RemovalStats;
	using TopologyFixScoreCache = typename Types::TopologyFixScoreCache;
	using MeshOps = topology_detail::SkeletonMeshOps<POINTS, NONMANIFOLD>;
	using QuadratureScores = topology_detail::QuadratureScores<POINTS, NONMANIFOLD>;

	LazyGreedyTetRemover(Context& p, TopologyFixScoreCache& score_cache, ScoreEvaluator& evaluator)
		: p(p), score_cache(score_cache), evaluator(evaluator), edge_score_cache(score_cache.edge_scores),
		  face_score_cache(score_cache.face_scores)
	{
	}

	RemovalStats run()
	{
		RemovalStats stats;
		if (!QuadratureScores::initialize_topology_fix_score_cache(p, evaluator, score_cache))
			return stats;

		face_owner_tets.reserve(nb_cells<SkeletonFace>(p.data.skeleton));
		deleted_face_count_per_tet.reserve(p.skeleton_tets.size());
		for (const auto& kv : p.skeleton_tets)
		{
			const std::size_t tet_id = kv.first;
			const Tet& tet = kv.second;
			uint32 initial_deleted_faces = 0;
			for (uint32 i = 0; i < 4; ++i)
			{
				const SkeletonFace f = tet.faces[i];
				if (!f.is_valid())
				{
					++initial_deleted_faces;
					continue;
				}
				const uint32 idf = index_of(p.data.skeleton, f);
				if (idf == INVALID_INDEX)
				{
					++initial_deleted_faces;
					continue;
				}
				face_owner_tets[idf].push_back(tet_id);
			}
			deleted_face_count_per_tet[tet_id] = initial_deleted_faces;
		}

		uint32 removed_faces = 0;
		uint32 removed_edges = 0;
		uint32 removed_tets = 0;
		tet_versions.reserve(p.skeleton_tets.size() * 2 + 1);
		for (const auto& kv : p.skeleton_tets)
		{
			tet_versions[kv.first] = 0;
			queue_tet_candidate(kv.first, false);
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

				auto it_version = tet_versions.find(entry.tet_id);
				if (it_version == tet_versions.end() || it_version->second != entry.version)
					continue;
				if (!p.skeleton_tets.contains(entry.tet_id))
				{
					tet_versions.erase(entry.tet_id);
					continue;
				}

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
			{
				break;
			}

			std::unordered_set<std::size_t> first_removed_tets;
			for (std::size_t tet_id : p.skeleton_tets.tets_of(best.first_face_id))
				if (p.skeleton_tets.contains(tet_id))
					first_removed_tets.insert(tet_id);

			uint32 step_removed_faces = 0;
			uint32 step_removed_edges = 0;
			uint32 step_removed_tets = 0;
			std::unordered_set<std::size_t> step_tets_to_erase;
			std::unordered_set<std::size_t> first_tets_to_erase;
			std::vector<SkeletonEdge> step_affected_edges;
			if (!remove_faces_collect_tets({best.first_face}, first_tets_to_erase, step_removed_faces,
										   step_removed_edges, &step_affected_edges))
			{
				queue_tet_candidate(best.tet_id, true);
				continue;
			}
			step_tets_to_erase.insert(first_tets_to_erase.begin(), first_tets_to_erase.end());

			std::vector<FaceChoice> created_deg1_faces;
			created_deg1_faces.reserve(16);
			std::unordered_set<uint32> seen_created_ids;
			seen_created_ids.reserve(32);
			for (std::size_t tet_id : first_removed_tets)
			{
				if (!p.skeleton_tets.contains(tet_id))
					continue;
				const Tet& tet = p.skeleton_tets.tet(tet_id);
				for (uint32 i = 0; i < 4; ++i)
				{
					const SkeletonFace f = tet.faces[i];
					const uint32 idf = index_of(p.data.skeleton, f);
					if (idf == INVALID_INDEX)
					{
						continue;
					}
					if (idf == best.first_face_id)
						continue;
					if (!seen_created_ids.insert(idf).second)
						continue;
					if (!has_face_budget(idf))
						continue;
					if (!face_has_degree1_edge(f))
						continue;
					if (face_has_edge_with_tet_face_count_gt2(f, &step_tets_to_erase))
						continue;
					created_deg1_faces.push_back({f, get_face_score(idf)});
				}
			}
			std::sort(created_deg1_faces.begin(), created_deg1_faces.end(),
					  [](const FaceChoice& a, const FaceChoice& b) { return a.score > b.score; });
			if (created_deg1_faces.size() > Policy::max_second_step_faces)
				created_deg1_faces.resize(Policy::max_second_step_faces);

			for (const FaceChoice& fc : created_deg1_faces)
			{
				if (!fc.face.is_valid())
					continue;
				const uint32 idf_now = index_of(p.data.skeleton, fc.face);
				if (idf_now == INVALID_INDEX)
					continue;
				if (!has_face_budget(idf_now))
					continue;
				if (!face_has_degree1_edge(fc.face))
					continue;
				if (face_has_edge_with_tet_face_count_gt2(fc.face, &step_tets_to_erase))
					continue;
				std::unordered_set<std::size_t> second_tets_to_erase;
				uint32 rf = 0, re = 0;
				if (remove_faces_collect_tets({fc.face}, second_tets_to_erase, rf, re, &step_affected_edges))
				{
					step_removed_faces += rf;
					step_removed_edges += re;
					step_tets_to_erase.insert(second_tets_to_erase.begin(), second_tets_to_erase.end());
				}
			}
			finalize_erase_tets(step_tets_to_erase, step_removed_tets);
			for (std::size_t tet_id : step_tets_to_erase)
				tet_versions.erase(tet_id);

			std::unordered_set<std::size_t> dirty_tets;
			collect_dirty_tets_from_edges(step_affected_edges, dirty_tets);
			for (std::size_t tet_id : dirty_tets)
				queue_tet_candidate(tet_id, true);

			removed_faces += step_removed_faces;
			removed_edges += step_removed_edges;
			removed_tets += step_removed_tets;
		}

		stats.removed_faces = removed_faces;
		stats.removed_edges = removed_edges;
		stats.removed_tets = removed_tets;
		return stats;
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

	Context& p;
	TopologyFixScoreCache& score_cache;
	ScoreEvaluator& evaluator;
	std::unordered_map<uint32, Scalar>& edge_score_cache;
	std::unordered_map<uint32, Scalar>& face_score_cache;
	std::unordered_map<uint32, std::vector<std::size_t>> face_owner_tets;
	std::unordered_map<std::size_t, uint32> deleted_face_count_per_tet;
	std::unordered_map<std::size_t, uint32> tet_versions;
	TetCandidateQueue candidate_queue;

	Scalar get_face_score(uint32 idf)
	{
		auto it = face_score_cache.find(idf);
		return (it != face_score_cache.end()) ? it->second : Scalar(0);
	}

	bool has_face_budget(uint32 face_id)
	{
		auto it_owner = face_owner_tets.find(face_id);
		if (it_owner == face_owner_tets.end())
			return false;
		for (std::size_t tet_id : it_owner->second)
		{
			if (!p.skeleton_tets.contains(tet_id))
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
			if (!p.skeleton_tets.contains(tet_id))
				continue;
			auto it_count = deleted_face_count_per_tet.find(tet_id);
			if (it_count != deleted_face_count_per_tet.end())
				++(it_count->second);
		}
	}

	bool face_belongs_to_tet(const SkeletonFace& f, std::size_t tet_id)
	{
		if (!f.is_valid())
			return false;
		const uint32 idf = index_of(p.data.skeleton, f);
		if (idf == INVALID_INDEX)
			return false;
		const auto& in_tets = p.skeleton_tets.tets_of(idf);
		return in_tets.find(tet_id) != in_tets.end();
	}

	bool face_has_degree1_edge(const SkeletonFace& f)
	{
		if (!f.is_valid())
			return false;
		for (SkeletonEdge e : incident_edges(p.data.skeleton, f))
		{
			if (!e.is_valid())
				continue;
			if (incident_faces(p.data.skeleton, e).size() == 1)
				return true;
		}
		return false;
	}

	bool face_has_edge_with_tet_face_count_gt2(const SkeletonFace& f,
											   const std::unordered_set<std::size_t>* ignored_tets)
	{
		if (!f.is_valid())
			return false;
		for (SkeletonEdge e : incident_edges(p.data.skeleton, f))
		{
			if (!e.is_valid())
				continue;
			size_t tet_face_count = 0;
			for (SkeletonFace ef : incident_faces(p.data.skeleton, e))
			{
				if (!ef.is_valid())
					continue;
				const uint32 idef = index_of(p.data.skeleton, ef);
				if (idef == INVALID_INDEX)
					continue;
				bool has_alive_owner = false;
				for (std::size_t tet_id : p.skeleton_tets.tets_of(idef))
				{
					if (!p.skeleton_tets.contains(tet_id))
						continue;
					if (ignored_tets && ignored_tets->find(tet_id) != ignored_tets->end())
						continue;
					has_alive_owner = true;
					break;
				}
				if (has_alive_owner)
					++tet_face_count;
			}
			if (tet_face_count > 2)
				return true;
		}
		return false;
	}

	bool tet_has_simple_face(std::size_t tet_id)
	{
		if (!p.skeleton_tets.contains(tet_id))
			return false;
		const Tet& tet = p.skeleton_tets.tet(tet_id);
		for (uint32 i = 0; i < 4; ++i)
		{
			const SkeletonFace f = tet.faces[i];
			if (!f.is_valid())
				continue;
			const uint32 idf = index_of(p.data.skeleton, f);
			if (idf == INVALID_INDEX)
				continue;
			const auto& in_tets = p.skeleton_tets.tets_of(idf);
			if (in_tets.find(tet_id) == in_tets.end())
				continue;
			if (in_tets.size() == 1)
				return true;
		}
		return false;
	}

	bool pick_tet_best_edge_and_face(std::size_t tet_id, SkeletonEdge& out_edge, uint32& out_edge_id,
									 Scalar& out_edge_score, SkeletonFace& out_face, uint32& out_face_id,
									 Scalar& out_face_score)
	{
		if (!p.skeleton_tets.contains(tet_id))
			return false;
		const Tet& tet = p.skeleton_tets.tet(tet_id);
		std::unordered_set<uint32> seen_edges;
		seen_edges.reserve(8);
		bool found = false;
		const Scalar edge_eps = Scalar(1e-12);
		const Scalar face_eps = Scalar(1e-12);
		for (uint32 i = 0; i < 4; ++i)
		{
			const SkeletonFace f = tet.faces[i];
			if (!face_belongs_to_tet(f, tet_id))
				continue;
			for (SkeletonEdge e : incident_edges(p.data.skeleton, f))
			{
				if (!e.is_valid())
					continue;
				const uint32 ide = index_of(p.data.skeleton, e);
				if (ide == INVALID_INDEX)
					continue;
				if (!seen_edges.insert(ide).second)
					continue;
				SkeletonFace best_face_on_edge;
				uint32 best_face_on_edge_id = INVALID_INDEX;
				Scalar best_face_on_edge_score = Scalar(0);
				bool edge_has_candidate = false;
				for (const SkeletonFace& ef : incident_faces(p.data.skeleton, e))
				{
					if (!face_belongs_to_tet(ef, tet_id))
						continue;
					const uint32 idf = index_of(p.data.skeleton, ef);
					if (idf == INVALID_INDEX || !has_face_budget(idf))
						continue;
					if (!Policy::accepts_first_face(p.skeleton_tets.tets_of(idf).size()))
						continue;
					const Scalar face_score = get_face_score(idf);
					if (!edge_has_candidate || face_score > best_face_on_edge_score + face_eps)
					{
						edge_has_candidate = true;
						best_face_on_edge = ef;
						best_face_on_edge_id = idf;
						best_face_on_edge_score = face_score;
					}
				}
				if (!edge_has_candidate)
					continue;
				const auto it_score = edge_score_cache.find(ide);
				const Scalar score = (it_score != edge_score_cache.end()) ? it_score->second : Scalar(0);
				const bool better_edge = (!found) || (score > out_edge_score + edge_eps);
				const bool tie_edge = found && (std::abs(score - out_edge_score) <= edge_eps);
				const bool better_face_on_tie = tie_edge && (best_face_on_edge_score > out_face_score + face_eps);
				if (better_edge || better_face_on_tie)
				{
					found = true;
					out_edge = e;
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

	bool remove_faces_collect_tets(const std::vector<SkeletonFace>& faces_to_remove,
								   std::unordered_set<std::size_t>& out_tets_to_erase, uint32& out_removed_faces,
								   uint32& out_removed_edges, std::vector<SkeletonEdge>* out_affected_edges = nullptr)
	{
		out_tets_to_erase.clear();
		std::vector<std::pair<uint32, SkeletonFace>> unique_faces;
		unique_faces.reserve(faces_to_remove.size());
		std::unordered_set<uint32> seen_face_ids;
		seen_face_ids.reserve(faces_to_remove.size() * 2 + 1);
		for (const SkeletonFace& f : faces_to_remove)
		{
			if (!f.is_valid())
				continue;
			const uint32 idf = index_of(p.data.skeleton, f);
			if (idf == INVALID_INDEX)
				continue;
			if (!seen_face_ids.insert(idf).second)
				continue;
			unique_faces.push_back({idf, f});
		}
		if (unique_faces.empty())
			return false;

		std::unordered_set<std::size_t> tets_to_remove_local;
		tets_to_remove_local.reserve(16);
		for (const auto& fpair : unique_faces)
		{
			const uint32 idf = fpair.first;
			for (std::size_t tet_id : p.skeleton_tets.tets_of(idf))
				if (p.skeleton_tets.contains(tet_id))
					tets_to_remove_local.insert(tet_id);
		}
		if (tets_to_remove_local.empty())
			return false;

		std::vector<SkeletonEdge> affected_edges;
		affected_edges.reserve(unique_faces.size() * 3);
		for (const auto& fpair : unique_faces)
		{
			const SkeletonFace f = fpair.second;
			for (SkeletonEdge e : incident_edges(p.data.skeleton, f))
				affected_edges.push_back(e);
		}
		if (out_affected_edges)
			out_affected_edges->insert(out_affected_edges->end(), affected_edges.begin(), affected_edges.end());

		std::unordered_set<uint32> removed_face_ids;
		removed_face_ids.reserve(unique_faces.size() * 2 + 1);
		for (const auto& fpair : unique_faces)
		{
			const SkeletonFace f = fpair.second;
			if (!f.is_valid())
				continue;
			const uint32 idf_now = index_of(p.data.skeleton, f);
			if (idf_now == INVALID_INDEX)
				continue;
			account_face_deletion_budget(idf_now);
			removed_face_ids.insert(idf_now);
			remove_face(p.data.skeleton, f);
			++out_removed_faces;
		}

		QuadratureScores::erase_topology_fix_scores_for_removed_faces_and_orphan_edges(
			p, &score_cache, removed_face_ids, affected_edges);
		out_removed_edges += MeshOps::remove_orphan_edges_from_removed_face_edges(p, affected_edges);
		out_tets_to_erase = std::move(tets_to_remove_local);
		return true;
	}

	void finalize_erase_tets(const std::unordered_set<std::size_t>& tets_to_erase, uint32& out_removed_tets)
	{
		for (std::size_t tet_id : tets_to_erase)
		{
			if (!p.skeleton_tets.contains(tet_id))
				continue;
			p.skeleton_tets.erase(tet_id);
			++out_removed_tets;
		}
	}

	bool build_step_candidate(std::size_t tet_id, StepCandidate& out_candidate)
	{
		const bool is_simple_tet = tet_has_simple_face(tet_id);
		if (!Policy::accepts_tet(is_simple_tet))
			return false;

		SkeletonEdge best_edge_local;
		uint32 best_edge_local_id = INVALID_INDEX;
		Scalar best_edge_local_score = Scalar(0);
		SkeletonFace first_face_local;
		uint32 first_face_local_id = INVALID_INDEX;
		Scalar first_face_local_score = Scalar(0);
		if (!pick_tet_best_edge_and_face(tet_id, best_edge_local, best_edge_local_id, best_edge_local_score,
										 first_face_local, first_face_local_id, first_face_local_score))
			return false;

		out_candidate.tet_id = tet_id;
		out_candidate.edge_id = best_edge_local_id;
		out_candidate.edge_score = best_edge_local_score;
		out_candidate.first_face = first_face_local;
		out_candidate.first_face_id = first_face_local_id;
		out_candidate.first_face_score = first_face_local_score;
		return true;
	}

	void queue_tet_candidate(std::size_t tet_id, bool bump_version)
	{
		if (!p.skeleton_tets.contains(tet_id))
		{
			tet_versions.erase(tet_id);
			return;
		}
		uint32& version = tet_versions[tet_id];
		if (bump_version)
			++version;
		StepCandidate candidate;
		if (!build_step_candidate(tet_id, candidate))
			return;
		candidate_queue.push({candidate.tet_id, candidate.edge_id, candidate.first_face_id, candidate.edge_score,
							  candidate.first_face_score, version});
	}

	void collect_dirty_tets_from_edges(const std::vector<SkeletonEdge>& affected_edges,
									   std::unordered_set<std::size_t>& dirty_tets)
	{
		dirty_tets.clear();
		for (SkeletonEdge e : affected_edges)
		{
			if (!e.is_valid())
				continue;
			for (SkeletonFace f : incident_faces(p.data.skeleton, e))
			{
				if (!f.is_valid())
					continue;
				const uint32 idf = index_of(p.data.skeleton, f);
				if (idf == INVALID_INDEX)
					continue;
				for (std::size_t tet_id : p.skeleton_tets.tets_of(idf))
					if (p.skeleton_tets.contains(tet_id))
						dirty_tets.insert(tet_id);
			}
		}
	}
};

template <typename POINTS, typename NONMANIFOLD>
struct GreedyTetRemoval
{
	using Types = SkeletonTopologyTypes<POINTS, NONMANIFOLD>;
	using SkeletonEdge = typename Types::SkeletonEdge;
	using SkeletonFace = typename Types::SkeletonFace;
	using Context = typename Types::Context;
	using Tet = typename Types::Tet;
	using RemovalStats = typename Types::RemovalStats;
	using TopologyFixScoreCache = typename Types::TopologyFixScoreCache;
	using MeshOps = topology_detail::SkeletonMeshOps<POINTS, NONMANIFOLD>;
	using QuadratureScores = topology_detail::QuadratureScores<POINTS, NONMANIFOLD>;

	template <typename Policy, typename ScoreEvaluator>
	static RemovalStats run_edge_score_tet_mode_topology_fix(Context& p, TopologyFixScoreCache& score_cache,
															 ScoreEvaluator& evaluator)
	{
		return LazyGreedyTetRemover<POINTS, NONMANIFOLD, Policy, ScoreEvaluator>(p, score_cache, evaluator).run();
	}
};

} // namespace topology_detail

} // namespace offset_axis

} // namespace geometry

} // namespace cgogn

#endif // CGOGN_GEOMETRY_ALGOS_OFFSET_AXIS_INTERNAL_GREEDY_TET_REMOVAL_H_
