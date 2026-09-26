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
#ifndef CGOGN_GEOMETRY_ALGOS_UDF_SKELETON_TOPOLOGY_H_
#define CGOGN_GEOMETRY_ALGOS_UDF_SKELETON_TOPOLOGY_H_

#include <cgogn/core/functions/attributes.h>
#include <cgogn/core/functions/mesh_info.h>
#include <cgogn/core/functions/traversals/edge.h>
#include <cgogn/geometry/algos/area.h>
#include <cgogn/geometry/algos/udf/spheres_optimizer.h>
#include <Eigen/Dense>
#include <Eigen/Eigenvalues>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <limits>
#include <map>
#include <memory>
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

template <typename POINTS, typename NONMANIFOLD>
class SkeletonTopology
{
public:
	using SphereVertex = typename mesh_traits<POINTS>::Vertex;
	using SkeletonVertex = typename mesh_traits<NONMANIFOLD>::Vertex;
	using SkeletonEdge = typename mesh_traits<NONMANIFOLD>::Edge;
	using SkeletonFace = typename mesh_traits<NONMANIFOLD>::Face;
	template <typename T> using SphereAttribute = typename mesh_traits<POINTS>::template Attribute<T>;
	template <typename T> using SkeletonAttribute = typename mesh_traits<NONMANIFOLD>::template Attribute<T>;

	struct Data
	{
		NONMANIFOLD* skeleton = nullptr;
		SpheresOptimizer<POINTS>* spheres_optimizer = nullptr;
		std::shared_ptr<SphereAttribute<SkeletonVertex>> sphere_skeleton_vertex;
		std::shared_ptr<SkeletonAttribute<Vec3>> skeleton_position;
		std::shared_ptr<SkeletonAttribute<Scalar>> skeleton_radius;
		std::shared_ptr<SkeletonAttribute<SphereVertex>> skeleton_source_sphere;
		std::shared_ptr<SkeletonAttribute<std::set<std::size_t>>> face_incident_tets;
		std::shared_ptr<SkeletonAttribute<uint32>> face_component_id;
		std::shared_ptr<SkeletonAttribute<uint32>> edge_degree;
	};

	struct Metrics
	{
		uint32 tet_count = 0;
		uint32 topology_rounds = 0;
		uint32 boundary_removed_tets = 0;
		uint32 simple_removed_tets = 0;
		uint32 nonsimple_removed_tets = 0;
		uint32 removed_faces = 0;
		uint32 removed_edges = 0;
		uint32 removed_tets = 0;
		uint32 repaired_dense_regions = 0;
		uint32 added_faces = 0;
		uint32 residual_removed_sheets = 0;
		uint32 residual_removed_faces = 0;
		std::vector<SphereVertex> added_spheres;
	};

	enum class Status { success, invalid_data, score_evaluation_failed };

	explicit SkeletonTopology(Data data) : data_(std::move(data)) {}
	SkeletonTopology(const SkeletonTopology&) = delete;
	SkeletonTopology& operator=(const SkeletonTopology&) = delete;
	Data& data() { return data_; }
	const Data& data() const { return data_; }
	const Metrics& metrics() const { return metrics_; }

	Status build_from_spheres()
	{
		if (!valid_data())
			return Status::invalid_data;
		metrics_ = Metrics{};
		TopologyParameters p(data_, skeleton_tets_);
		if (!compute_skeleton(p))
			return Status::invalid_data;
		metrics_.tet_count = static_cast<uint32>(skeleton_tets_.size());
		return Status::success;
	}

	template <typename ScoreEvaluator>
	Status fix_topology(ScoreEvaluator&& evaluate_scores)
	{
		if (!valid_data())
			return Status::invalid_data;
		metrics_ = Metrics{};
		if (skeleton_tets_.empty())
			return Status::success;
		TopologyParameters p(data_, skeleton_tets_);
		return fix_topology_impl(p, evaluate_scores);
	}

	Status prune_fully_non_manifold_triangles()
	{
		if (!valid_data())
			return Status::invalid_data;
		TopologyParameters p(data_, skeleton_tets_);
		return prune_fully_non_manifold_triangles(p) ? Status::success : Status::invalid_data;
	}

	Status compute_skeleton_face_components_union_find()
	{
		if (!valid_data())
			return Status::invalid_data;
		TopologyParameters p(data_, skeleton_tets_);
		return compute_skeleton_face_components_union_find(p) ? Status::success : Status::invalid_data;
	}

	void compute_edge_degree()
	{
		if (!valid_data())
			return;
		TopologyParameters p(data_, skeleton_tets_);
		compute_edge_degree(p);
	}

	Status prune_residual_sheets()
	{
		if (!valid_data())
			return Status::invalid_data;
		metrics_ = Metrics{};
		TopologyParameters p(data_, skeleton_tets_);
		run_residual_sheet_prune(p);
		metrics_.tet_count = static_cast<uint32>(skeleton_tets_.size());
		return Status::success;
	}

private:
	using SpheresOptimizerType = SpheresOptimizer<POINTS>;
	using SpheresOptimizerData = typename SpheresOptimizerType::Data;
	using FaceKey = std::array<uint32, 3>;

	struct Tet { SkeletonFace faces[4]; };
	using TetMap = std::unordered_map<std::size_t, Tet>;
	struct TopologyParameters
	{
		Data& data;
		SpheresOptimizerData& optimizer_data;
		TetMap& skeleton_tets_;
		NONMANIFOLD*& skeleton_;
		POINTS*& samples_mesh_;
		std::shared_ptr<SphereAttribute<Vec3>>& samples_position_;
		std::shared_ptr<SphereAttribute<Scalar>>& samples_area_;
		std::shared_ptr<SphereAttribute<std::vector<SphereVertex>>>& samples_knn_;
		std::shared_ptr<SphereAttribute<Spherical_Quadric>>& samples_quadric_;
		std::shared_ptr<SphereAttribute<Line_Quadric>>& samples_line_quadric_;
		std::shared_ptr<SphereAttribute<Vec3>>& samples_ma_position_;
		std::shared_ptr<SphereAttribute<Scalar>>& samples_ma_radius_;
		std::shared_ptr<SphereAttribute<SphereVertex>>& samples_ma_secondary_vertex_;
		std::shared_ptr<SphereAttribute<SphereVertex>>& samples_sphere_;
		std::shared_ptr<SphereAttribute<Scalar>>& samples_error_;
		const acc::KDTree<3, uint32>*& samples_kdtree_;
		const std::vector<SphereVertex>*& samples_kdtree_vertices_;
		POINTS*& spheres_;
		std::shared_ptr<SphereAttribute<Vec3>>& spheres_position_;
		std::shared_ptr<SphereAttribute<Scalar>>& spheres_radius_;
		std::shared_ptr<SphereAttribute<std::vector<SphereVertex>>>& spheres_cluster_;
		std::shared_ptr<SphereAttribute<Scalar>>& spheres_cluster_area_;
		std::shared_ptr<SphereAttribute<std::set<SphereVertex>>>& spheres_neighbor_clusters_;
		std::shared_ptr<SphereAttribute<Scalar>>& spheres_error_;
		std::shared_ptr<SphereAttribute<Scalar>>& spheres_error_not_normalized_;
		std::shared_ptr<SphereAttribute<SkeletonVertex>>& spheres_skeleton_vertex_;
		std::shared_ptr<SkeletonAttribute<Vec3>>& skeleton_position_;
		std::shared_ptr<SkeletonAttribute<Scalar>>& skeleton_radius_;
		std::shared_ptr<SkeletonAttribute<std::set<std::size_t>>>& incident_tets_;
		std::shared_ptr<SkeletonAttribute<uint32>>& skeleton_face_component_id_;
		std::shared_ptr<SkeletonAttribute<SphereVertex>>& skeleton_source_sphere_;
		std::shared_ptr<SkeletonAttribute<uint32>>& edge_degree_;
		SpheresOptimizerType*& spheres_optimizer_;
		Scalar sqem_update_lambda_line_plane_ = Scalar(0.20);

		TopologyParameters(Data& d, TetMap& t) : data(d), optimizer_data(d.spheres_optimizer->data()), skeleton_tets_(t),
			skeleton_(d.skeleton), samples_mesh_(optimizer_data.samples_mesh), samples_position_(optimizer_data.sample_position),
			samples_area_(optimizer_data.sample_area), samples_knn_(optimizer_data.sample_knn),
			samples_quadric_(optimizer_data.sample_quadric), samples_line_quadric_(optimizer_data.sample_line_quadric),
			samples_ma_position_(optimizer_data.sample_ma_position), samples_ma_radius_(optimizer_data.sample_ma_radius),
			samples_ma_secondary_vertex_(optimizer_data.sample_ma_secondary_vertex), samples_sphere_(optimizer_data.sample_sphere),
			samples_error_(optimizer_data.sample_error), samples_kdtree_(optimizer_data.sample_kdtree),
			samples_kdtree_vertices_(optimizer_data.sample_kdtree_vertices),
			spheres_(optimizer_data.spheres), spheres_position_(optimizer_data.sphere_position), spheres_radius_(optimizer_data.sphere_radius),
			spheres_cluster_(optimizer_data.sphere_cluster), spheres_cluster_area_(optimizer_data.sphere_cluster_area),
			spheres_neighbor_clusters_(optimizer_data.sphere_neighbors),
			spheres_error_(optimizer_data.sphere_error), spheres_error_not_normalized_(optimizer_data.sphere_error_not_normalized),
			spheres_skeleton_vertex_(d.sphere_skeleton_vertex), skeleton_position_(d.skeleton_position),
			skeleton_radius_(d.skeleton_radius), incident_tets_(d.face_incident_tets),
			skeleton_face_component_id_(d.face_component_id), skeleton_source_sphere_(d.skeleton_source_sphere),
			edge_degree_(d.edge_degree), spheres_optimizer_(d.spheres_optimizer) {}
	};

	bool valid_data() const
	{
		return data_.skeleton && data_.spheres_optimizer && data_.skeleton_position && data_.skeleton_radius &&
			data_.skeleton_source_sphere && data_.face_incident_tets && data_.face_component_id && data_.edge_degree &&
			data_.sphere_skeleton_vertex;
	}

	template <typename ScoreEvaluator>
	Status fix_topology_impl(TopologyParameters& p, ScoreEvaluator& evaluate_scores)
	{
		return run_topology_fix_pipeline(p, evaluate_scores);
	}


	bool remove_skeleton_vertex_and_linked_sphere(TopologyParameters& p, const SkeletonVertex& v)
	{
		if (!p.skeleton_ || !v.is_valid())
			return false;
		const uint32 skeleton_vertex_id = index_of(*p.skeleton_, v);
		if (skeleton_vertex_id == INVALID_INDEX)
			return false;

		SphereVertex linked_sphere;
		if (p.skeleton_source_sphere_)
			linked_sphere = value<SphereVertex>(*p.skeleton_, p.skeleton_source_sphere_, v);

		if (linked_sphere.is_valid() && p.spheres_)
		{
			const uint32 sphere_index = index_of(*p.spheres_, linked_sphere);
			if (sphere_index != INVALID_INDEX && p.spheres_skeleton_vertex_)
				(*p.spheres_skeleton_vertex_)[sphere_index] = SkeletonVertex();
			p.spheres_optimizer_->remove_sphere(linked_sphere);
		}

		if (p.skeleton_source_sphere_)
			value<SphereVertex>(*p.skeleton_, p.skeleton_source_sphere_, v) = SphereVertex();
		remove_vertex(*p.skeleton_, v);
		return true;
	}

	struct edge_hash
	{
		std::size_t operator()(const std::pair<uint32, uint32>& edge) const
		{
			const uint32 a = std::min(edge.first, edge.second);
			const uint32 b = std::max(edge.first, edge.second);
			return std::hash<uint64>()((uint64(a) << 32) | uint64(b));
		}
	};

	struct edge_equal
	{
		bool operator()(const std::pair<uint32, uint32>& edge1, const std::pair<uint32, uint32>& edge2) const
		{
			return ((edge1.first == edge2.first && edge1.second == edge2.second) ||
					(edge1.first == edge2.second && edge1.second == edge2.first));
		}
	};

	bool compute_skeleton(TopologyParameters& p)
	{
		p.spheres_optimizer_->compute_sphere_neighbors();
		if (!p.samples_mesh_ || !p.samples_sphere_ || !p.samples_knn_)
			return false;
		clear(*p.skeleton_);
		std::map<FaceKey, SkeletonFace> skeleton_faces_map;
		auto get_face_key = [](uint32 i1, uint32 i2, uint32 i3) -> FaceKey {
			std::array<uint32, 3> key = {i1, i2, i3};
			std::sort(key.begin(), key.end());
			return key;
		};
		std::vector<std::array<uint32, 4>> raw_tets;
		if (p.spheres_skeleton_vertex_)
		{
			foreach_cell(*p.spheres_, [&](SphereVertex pv) -> bool {
				const uint32 sphere_index = index_of(*p.spheres_, pv);
				if (sphere_index != INVALID_INDEX)
					(*p.spheres_skeleton_vertex_)[sphere_index] = SkeletonVertex();
				return true;
			});
		}
		foreach_cell(*p.spheres_, [&](SphereVertex pv) -> bool {
			uint32 pv_index = index_of(*p.spheres_, pv);
			SkeletonVertex nmv = add_vertex(*p.skeleton_);
			const uint32 skeleton_index = index_of(*p.skeleton_, nmv);
			(*p.skeleton_position_)[skeleton_index] = (*p.spheres_position_)[pv_index];
			(*p.skeleton_radius_)[skeleton_index] = (*p.spheres_radius_)[pv_index];
			if (p.spheres_skeleton_vertex_)
				(*p.spheres_skeleton_vertex_)[pv_index] = nmv;
			if (p.skeleton_source_sphere_)
				value<SphereVertex>(*p.skeleton_, p.skeleton_source_sphere_, nmv) = pv;
			return true;
		});

		std::unordered_map<std::pair<uint32, uint32>, SkeletonEdge, edge_hash, edge_equal> edge_indices;

		foreach_cell(*p.spheres_, [&](SphereVertex pv) -> bool {
			uint32 pv_index = index_of(*p.spheres_, pv);
			SkeletonVertex nmv1 = p.spheres_skeleton_vertex_ ? (*p.spheres_skeleton_vertex_)[pv_index] : SkeletonVertex();
			const std::set<SphereVertex>& neighbors = (*p.spheres_neighbor_clusters_)[pv_index];
			for (SphereVertex neighbor : neighbors)
			{
				uint32 n_index = index_of(*p.spheres_, neighbor);
				SkeletonVertex nmv2 = p.spheres_skeleton_vertex_ ? (*p.spheres_skeleton_vertex_)[n_index] : SkeletonVertex();
				std::vector<SkeletonVertex> av = adjacent_vertices_through_edge(*p.skeleton_, nmv1);
				if (std::find(av.begin(), av.end(), nmv2) == av.end())
				{
					SkeletonEdge e = add_edge(*p.skeleton_, nmv1, nmv2);
					edge_indices[{index_of(*p.skeleton_, nmv1), index_of(*p.skeleton_, nmv2)}] = e;
				}
			}
			return true;
		});

		foreach_cell(*p.spheres_, [&](SphereVertex pv) -> bool {
			uint32 idx1 = index_of(*p.spheres_, pv);
			SkeletonVertex nmv1 = p.spheres_skeleton_vertex_ ? (*p.spheres_skeleton_vertex_)[idx1] : SkeletonVertex();
			const std::set<SphereVertex>& n_pv = (*p.spheres_neighbor_clusters_)[idx1];
			for (const SphereVertex& ne1 : n_pv)
			{
				uint32 idx2 = index_of(*p.spheres_, ne1);
				if (idx1 >= idx2)
					continue;
				SkeletonVertex nmv2 = p.spheres_skeleton_vertex_ ? (*p.spheres_skeleton_vertex_)[idx2] : SkeletonVertex();
				const std::set<SphereVertex>& ne_ne1 = (*p.spheres_neighbor_clusters_)[idx2];
				for (const SphereVertex& ne2 : ne_ne1)
				{
					if (n_pv.find(ne2) == n_pv.end())
						continue;
					uint32 idx3 = index_of(*p.spheres_, ne2);
					if (idx2 >= idx3)
						continue;

					SkeletonVertex nmv3 = p.spheres_skeleton_vertex_ ? (*p.spheres_skeleton_vertex_)[idx3] : SkeletonVertex();

					std::vector<SkeletonEdge> edges;
					edges.reserve(3);
					edges.push_back(edge_indices[{index_of(*p.skeleton_, nmv1), index_of(*p.skeleton_, nmv2)}]);
					edges.push_back(edge_indices[{index_of(*p.skeleton_, nmv2), index_of(*p.skeleton_, nmv3)}]);
					edges.push_back(edge_indices[{index_of(*p.skeleton_, nmv1), index_of(*p.skeleton_, nmv3)}]);
					SkeletonFace new_face = add_face(*p.skeleton_, edges);

					skeleton_faces_map[get_face_key(idx1, idx2, idx3)] = new_face;

					const std::set<SphereVertex>& ne_ne2 = (*p.spheres_neighbor_clusters_)[idx3];

					for (const SphereVertex& ne3 : ne_ne2)
					{
						uint32 idx4 = index_of(*p.spheres_, ne3);
						if (idx3 >= idx4)
							continue;
						bool connected_v1 = (n_pv.find(ne3) != n_pv.end());
						bool connected_v2 = (ne_ne1.find(ne3) != ne_ne1.end());
						if (connected_v1 && connected_v2)
						{
							raw_tets.push_back({idx1, idx2, idx3, idx4});
						}
					}
				}
			}
			return true;
		});

		// Resolve Tets
		p.skeleton_tets_.clear();
		p.skeleton_tets_.reserve(raw_tets.size());
		foreach_cell(*p.skeleton_, [&](SkeletonFace f) -> bool {
			value<std::set<std::size_t>>(*p.skeleton_, p.incident_tets_, f).clear();
			return true;
		});
		std::size_t tet_index = 0;
		for (const auto& rt : raw_tets)
		{
			Tet new_tet;
			uint32 v[4] = {rt[0], rt[1], rt[2], rt[3]};

			FaceKey keys[4] = {get_face_key(v[1], v[2], v[3]), get_face_key(v[0], v[2], v[3]),
								 get_face_key(v[0], v[1], v[3]), get_face_key(v[0], v[1], v[2])};

			for (uint32 i = 0; i < 4; ++i)
			{
				auto it = skeleton_faces_map.find(keys[i]);
				if (it != skeleton_faces_map.end())
				{
					SkeletonFace f = it->second;
					new_tet.faces[i] = f;

					value<std::set<size_t>>(*p.skeleton_, p.incident_tets_, f).insert(tet_index);

				}
			}
			p.skeleton_tets_.insert({tet_index, new_tet});
			tet_index++;
		}
		compute_edge_degree(p);
		return true;
	}


	//------------------------------//
	//-----Topology correction------//
	//------------------------------//

	void compute_edge_degree(TopologyParameters& p)
	{
		parallel_foreach_cell(*p.skeleton_, [&](SkeletonEdge e) {
			auto in_face = incident_faces(*p.skeleton_, e);
			const uint32 ide = index_of(*p.skeleton_, e);
			const uint32 deg = static_cast<uint32>(in_face.size());

			if (ide != INVALID_INDEX && p.edge_degree_)
				(*p.edge_degree_)[ide] = deg;

			return true;
		});
	}

	void update_edge_degree_for_edge_ids(TopologyParameters& p, const std::unordered_set<uint32>& edge_ids)
	{
		if (!p.skeleton_ || !p.edge_degree_)
			return;
		for (uint32 edge_id : edge_ids)
		{
			if (edge_id == INVALID_INDEX)
				continue;
			const SkeletonEdge e = of_index<SkeletonEdge>(*p.skeleton_, edge_id);
			if (!e.is_valid())
				continue;
			if (index_of(*p.skeleton_, e) != edge_id)
				continue;

			const uint32 deg = static_cast<uint32>(incident_faces(*p.skeleton_, e).size());
			(*p.edge_degree_)[edge_id] = deg;
		}
	}

	struct FaceComponentUnionFind
	{
		std::vector<uint32> parent_;
		std::vector<uint32> rank_;

		FaceComponentUnionFind() = default;

		explicit FaceComponentUnionFind(size_t n)
		{
			reset(n);
		}

		void reset(size_t n)
		{
			parent_.resize(n);
			rank_.assign(n, 0);
			for (uint32 i = 0; i < static_cast<uint32>(n); ++i)
				parent_[i] = i;
		}

		uint32 find(uint32 x)
		{
			uint32 root = x;
			while (parent_[root] != root)
				root = parent_[root];
			while (parent_[x] != x)
			{
				const uint32 next = parent_[x];
				parent_[x] = root;
				x = next;
			}
			return root;
		}

		void unite(uint32 a, uint32 b)
		{
			uint32 ra = find(a);
			uint32 rb = find(b);
			if (ra == rb)
				return;
			if (rank_[ra] < rank_[rb])
				std::swap(ra, rb);
			parent_[rb] = ra;
			if (rank_[ra] == rank_[rb])
				++rank_[ra];
		}
	};



	bool prune_fully_non_manifold_triangles(TopologyParameters& p)
	{
		if (!p.skeleton_)
			return false;

		bool removed_any_faces = false;

		for (;;)
		{
			std::vector<uint32> face_ids_to_remove;
			foreach_cell(*p.skeleton_, [&](SkeletonFace f) -> bool {
				if (!f.is_valid())
					return true;
				const uint32 face_id = index_of(*p.skeleton_, f);
				if (face_id == INVALID_INDEX)
					return true;

				const std::vector<SkeletonEdge> face_edges = incident_edges(*p.skeleton_, f);
				if (face_edges.size() != 3)
					return true;

				for (const SkeletonEdge& e : face_edges)
				{
					if (!e.is_valid() || incident_faces(*p.skeleton_, e).size() <= 2)
						return true;
				}

				face_ids_to_remove.push_back(face_id);
				return true;
			});

			if (face_ids_to_remove.empty())
				break;

			for (uint32 face_id : face_ids_to_remove)
			{
				const SkeletonFace f = of_index<SkeletonFace>(*p.skeleton_, face_id);
				if (!f.is_valid())
					continue;

				const std::vector<SkeletonEdge> affected_edges = incident_edges(*p.skeleton_, f);
				remove_face(*p.skeleton_, f);
				removed_any_faces = true;

				std::vector<SkeletonVertex> candidate_vertices;
				candidate_vertices.reserve(affected_edges.size() * 2);
				for (const SkeletonEdge& e : affected_edges)
				{
					if (!e.is_valid())
						continue;
					for (const SkeletonVertex& v : incident_vertices(*p.skeleton_, e))
						candidate_vertices.push_back(v);
					if (incident_faces(*p.skeleton_, e).empty())
						remove_edge(*p.skeleton_, e);
				}
				for (const SkeletonVertex& v : candidate_vertices)
				{
					if (!v.is_valid())
						continue;
					const uint32 vertex_id = index_of(*p.skeleton_, v);
					if (vertex_id == INVALID_INDEX || !incident_edges(*p.skeleton_, v).empty())
						continue;
					remove_skeleton_vertex_and_linked_sphere(p, v);
				}
			}
		}

		if (p.edge_degree_ && removed_any_faces)
			compute_edge_degree(p);

		return true;
	}

	bool compute_skeleton_face_components_union_find(TopologyParameters& p)
	{
		if (!p.skeleton_ || !p.skeleton_face_component_id_)
			return false;

		std::vector<SkeletonFace> faces;
		faces.reserve(nb_cells<SkeletonFace>(*p.skeleton_));
		std::unordered_map<uint32, uint32> face_id_to_uf_index;
		face_id_to_uf_index.reserve(nb_cells<SkeletonFace>(*p.skeleton_));
		foreach_cell(*p.skeleton_, [&](SkeletonFace f) -> bool {
			if (!f.is_valid())
				return true;
			const uint32 face_id = index_of(*p.skeleton_, f);
			if (face_id == INVALID_INDEX)
				return true;
			(*p.skeleton_face_component_id_)[face_id] = INVALID_INDEX;
			face_id_to_uf_index.emplace(face_id, static_cast<uint32>(faces.size()));
			faces.push_back(f);
			return true;
		});

		FaceComponentUnionFind uf(faces.size());
		foreach_cell(*p.skeleton_, [&](SkeletonEdge e) -> bool {
			if (!e.is_valid())
				return true;
			const auto incident = incident_faces(*p.skeleton_, e);
			if (incident.size() != 2)
				return true;
			const uint32 f0 = index_of(*p.skeleton_, incident[0]);
			const uint32 f1 = index_of(*p.skeleton_, incident[1]);
			const auto i0 = face_id_to_uf_index.find(f0);
			const auto i1 = face_id_to_uf_index.find(f1);
			if (i0 != face_id_to_uf_index.end() && i1 != face_id_to_uf_index.end())
				uf.unite(i0->second, i1->second);
			return true;
		});

		std::unordered_map<uint32, uint32> root_to_component_id;
		root_to_component_id.reserve(faces.size());
		for (uint32 i = 0; i < static_cast<uint32>(faces.size()); ++i)
		{
			const uint32 root = uf.find(i);
			const auto [it, inserted] =
				root_to_component_id.emplace(root, static_cast<uint32>(root_to_component_id.size()));
			(void)inserted;
			const uint32 face_id = index_of(*p.skeleton_, faces[i]);
			if (face_id != INVALID_INDEX)
				(*p.skeleton_face_component_id_)[face_id] = it->second;
		}
		return true;
	}


	void run_residual_sheet_prune(TopologyParameters& p)
	{
		if (!p.skeleton_ || !p.skeleton_face_component_id_ || !p.incident_tets_)
			return;

		compute_edge_degree(p);
		prune_fully_non_manifold_triangles(p);
		if (!compute_skeleton_face_components_union_find(p))
			return;

		std::unordered_map<uint32, std::vector<uint32>> label_to_faces;
		label_to_faces.reserve(nb_cells<SkeletonFace>(*p.skeleton_));
		std::unordered_set<uint32> all_face_ids;
		all_face_ids.reserve(nb_cells<SkeletonFace>(*p.skeleton_));
		foreach_cell(*p.skeleton_, [&](SkeletonFace f) -> bool {
			if (!f.is_valid())
				return true;
			const uint32 face_id = index_of(*p.skeleton_, f);
			if (face_id == INVALID_INDEX)
				return true;
			const uint32 label = (*p.skeleton_face_component_id_)[face_id];
			if (label == INVALID_INDEX)
				return true;
			label_to_faces[label].push_back(face_id);
			all_face_ids.insert(face_id);
			return true;
		});
		if (all_face_ids.empty())
			return;

		std::vector<uint32> sorted_sheet_labels;
		sorted_sheet_labels.reserve(label_to_faces.size());
		for (const auto& kv : label_to_faces)
			sorted_sheet_labels.push_back(kv.first);
		std::sort(sorted_sheet_labels.begin(), sorted_sheet_labels.end());

		const size_t forced_delete_max_face_count = 2;
		auto sheet_has_vertex_connected_to_isolated_edge = [&](const std::vector<uint32>& face_ids) -> bool {
			std::unordered_set<uint32> visited_vertex_ids;
			visited_vertex_ids.reserve(face_ids.size() * 3);
			for (uint32 face_id : face_ids)
			{
				if (face_id == INVALID_INDEX)
					continue;
				const SkeletonFace f = of_index<SkeletonFace>(*p.skeleton_, face_id);
				if (!f.is_valid() || index_of(*p.skeleton_, f) != face_id)
					continue;
				const std::vector<SkeletonVertex> vertices = incident_vertices(*p.skeleton_, f);
				for (const SkeletonVertex& v : vertices)
				{
					if (!v.is_valid())
						continue;
					const uint32 vertex_id = index_of(*p.skeleton_, v);
					if (vertex_id == INVALID_INDEX || !visited_vertex_ids.insert(vertex_id).second)
						continue;
					for (const SkeletonEdge& e : incident_edges(*p.skeleton_, v))
					{
						if (!e.is_valid())
							continue;
						if (incident_faces(*p.skeleton_, e).empty())
							return true;
					}
				}
			}
			return false;
		};

		uint32 deleted_sheet_count = 0;
		std::unordered_set<uint32> face_ids_to_delete;
		for (uint32 sheet_label : sorted_sheet_labels)
		{
			const auto it_faces = label_to_faces.find(sheet_label);
			if (it_faces == label_to_faces.end() || it_faces->second.empty())
				continue;

			const size_t face_count = it_faces->second.size();
			const bool delete_by_small_sheet = (face_count <= forced_delete_max_face_count);
			const bool keep_by_isolated_edge_guard =
				delete_by_small_sheet && sheet_has_vertex_connected_to_isolated_edge(it_faces->second);
			const bool delete_sheet = delete_by_small_sheet && !keep_by_isolated_edge_guard;
			if (!delete_sheet)
				continue;

			++deleted_sheet_count;
			face_ids_to_delete.reserve(face_ids_to_delete.size() + it_faces->second.size());
			for (uint32 face_id : it_faces->second)
				face_ids_to_delete.insert(face_id);
		}

		if (face_ids_to_delete.empty())
		{
			compute_edge_degree(p);
			prune_fully_non_manifold_triangles(p);
			compute_skeleton_face_components_union_find(p);
			return;
		}

		SkeletonFaceDeletionStats deletion_stats;
		if (!delete_skeleton_face_id_set(p, face_ids_to_delete, deletion_stats))
			return;

		metrics_.residual_removed_sheets += deleted_sheet_count;
		metrics_.residual_removed_faces += deletion_stats.removed_faces;
		compute_edge_degree(p);
		prune_fully_non_manifold_triangles(p);
		compute_skeleton_face_components_union_find(p);
	}
	uint32 remove_orphan_edges_from_removed_face_edges(TopologyParameters& p, const std::vector<SkeletonEdge>& face_edges,
													   uint32* out_removed_vertices = nullptr)
	{
		uint32 removed_edges = 0;
		std::vector<SkeletonVertex> candidate_vertices;
		candidate_vertices.reserve(face_edges.size() * 2);
		for (SkeletonEdge e : face_edges)
		{
			if (!e.is_valid())
				continue;
			const std::vector<SkeletonVertex> edge_vertices = incident_vertices(*p.skeleton_, e);
			for (const SkeletonVertex& v : edge_vertices)
				candidate_vertices.push_back(v);
			if (incident_faces(*p.skeleton_, e).empty())
			{
				remove_edge(*p.skeleton_, e);
				++removed_edges;
			}
		}
		uint32 removed_vertices = 0;
		for (const SkeletonVertex& v : candidate_vertices)
		{
			if (!v.is_valid())
				continue;
			const uint32 idv = index_of(*p.skeleton_, v);
			if (idv == INVALID_INDEX)
				continue;
			if (!incident_edges(*p.skeleton_, v).empty())
				continue;
			if (remove_skeleton_vertex_and_linked_sphere(p, v))
				++removed_vertices;
		}
		if (out_removed_vertices)
			*out_removed_vertices = removed_vertices;
		return removed_edges;
	}

	std::unordered_set<uint32> collect_current_orphan_edge_ids(TopologyParameters& p)
	{
		std::unordered_set<uint32> orphan_edge_ids;
		if (!p.skeleton_)
			return orphan_edge_ids;
		orphan_edge_ids.reserve(nb_cells<SkeletonEdge>(*p.skeleton_));
		foreach_cell(*p.skeleton_, [&](SkeletonEdge e) -> bool {
			if (!e.is_valid())
				return true;
			const uint32 ide = index_of(*p.skeleton_, e);
			if (ide == INVALID_INDEX)
				return true;
			if (incident_faces(*p.skeleton_, e).empty())
				orphan_edge_ids.insert(ide);
			return true;
		});
		return orphan_edge_ids;
	}

	void remove_global_orphan_skeleton_elements(TopologyParameters& p, uint32* out_removed_edges = nullptr,
												uint32* out_removed_vertices = nullptr,
												const std::unordered_set<uint32>* preserved_orphan_edge_ids = nullptr)
	{
		uint32 removed_edges = 0;
		uint32 removed_vertices = 0;
		if (!p.skeleton_)
		{
			if (out_removed_edges)
				*out_removed_edges = 0;
			if (out_removed_vertices)
				*out_removed_vertices = 0;
			return;
		}

		std::vector<SkeletonEdge> orphan_edges;
		orphan_edges.reserve(nb_cells<SkeletonEdge>(*p.skeleton_));
		foreach_cell(*p.skeleton_, [&](SkeletonEdge e) -> bool {
			if (!e.is_valid())
				return true;
			if (!incident_faces(*p.skeleton_, e).empty())
				return true;
			const uint32 ide = index_of(*p.skeleton_, e);
			if (ide == INVALID_INDEX)
				return true;
			if (preserved_orphan_edge_ids && preserved_orphan_edge_ids->find(ide) != preserved_orphan_edge_ids->end())
				return true;
			orphan_edges.push_back(e);
			return true;
		});

		std::vector<SkeletonVertex> candidate_vertices;
		candidate_vertices.reserve(orphan_edges.size() * 2);
		for (const SkeletonEdge& e : orphan_edges)
		{
			if (!e.is_valid())
				continue;
			for (const SkeletonVertex& v : incident_vertices(*p.skeleton_, e))
				candidate_vertices.push_back(v);
			remove_edge(*p.skeleton_, e);
			++removed_edges;
		}

		foreach_cell(*p.skeleton_, [&](SkeletonVertex v) -> bool {
			if (!v.is_valid())
				return true;
			candidate_vertices.push_back(v);
			return true;
		});

		std::unordered_set<uint32> seen_vertex_ids;
		seen_vertex_ids.reserve(candidate_vertices.size());

		for (const SkeletonVertex& v : candidate_vertices)
		{
			if (!v.is_valid())
				continue;
			const uint32 idv = index_of(*p.skeleton_, v);
			if (idv == INVALID_INDEX)
				continue;
			if (!seen_vertex_ids.insert(idv).second)
				continue;
			if (!incident_edges(*p.skeleton_, v).empty())
				continue;
			if (remove_skeleton_vertex_and_linked_sphere(p, v))
				++removed_vertices;
		}

		if (out_removed_edges)
			*out_removed_edges = removed_edges;
		if (out_removed_vertices)
			*out_removed_vertices = removed_vertices;
	}

	struct SkeletonFaceDeletionStats
	{
		uint32 removed_faces = 0;
		uint32 removed_edges = 0;
		uint32 removed_vertices = 0;
		uint32 removed_tets = 0;
	};

	struct SkeletonFaceDeletionOptions
	{
		bool remove_incident_orphan_edges_immediately = true;
		bool remove_global_orphan_elements_after_batch = true;
	};

	bool delete_skeleton_face_id_set(TopologyParameters& p, const std::unordered_set<uint32>& face_ids_to_delete,
									 SkeletonFaceDeletionStats& out_stats,
									 const SkeletonFaceDeletionOptions& options = {})
	{
		out_stats = SkeletonFaceDeletionStats{};
		if (!p.skeleton_ || !p.incident_tets_ || face_ids_to_delete.empty())
			return false;
		const std::unordered_set<uint32> preexisting_orphan_edge_ids =
			options.remove_global_orphan_elements_after_batch ? collect_current_orphan_edge_ids(p) : std::unordered_set<uint32>{};

		for (uint32 face_id : face_ids_to_delete)
		{
			const SkeletonFace face = of_index<SkeletonFace>(*p.skeleton_, face_id);
			if (!face.is_valid())
				continue;
			const uint32 current_face_id = index_of(*p.skeleton_, face);
			if (current_face_id == INVALID_INDEX || current_face_id != face_id)
				continue;

			bool face_still_attached = false;
			for (SkeletonEdge e : incident_edges(*p.skeleton_, face))
			{
				if (!e.is_valid())
					continue;
				for (SkeletonFace ef : incident_faces(*p.skeleton_, e))
				{
					if (ef == face)
					{
						face_still_attached = true;
						break;
					}
				}
				if (face_still_attached)
					break;
			}
			if (!face_still_attached)
				continue;

			const std::set<std::size_t> in_tets = (*p.incident_tets_)[face_id];
			std::vector<SkeletonEdge> affected_edges = incident_edges(*p.skeleton_, face);
			remove_face(*p.skeleton_, face);
			++out_stats.removed_faces;

			if (options.remove_incident_orphan_edges_immediately)
			{
				uint32 removed_vertices_this_face = 0;
				out_stats.removed_edges +=
					remove_orphan_edges_from_removed_face_edges(p, affected_edges, &removed_vertices_this_face);
				out_stats.removed_vertices += removed_vertices_this_face;
			}

			for (std::size_t tet_id : in_tets)
			{
				auto it_tet = p.skeleton_tets_.find(tet_id);
				if (it_tet == p.skeleton_tets_.end())
					continue;
				const Tet old_tet = it_tet->second;
				for (uint32 i = 0; i < 4; ++i)
				{
					const SkeletonFace f = old_tet.faces[i];
					if (!f.is_valid() || f == face)
						continue;
					const uint32 idf = index_of(*p.skeleton_, f);
					if (idf != INVALID_INDEX)
						(*p.incident_tets_)[idf].erase(tet_id);
				}
				p.skeleton_tets_.erase(it_tet);
				++out_stats.removed_tets;
			}
		}

		if (out_stats.removed_faces == 0)
			return false;

		if (options.remove_global_orphan_elements_after_batch)
		{
			uint32 global_orphan_edges = 0;
			uint32 global_orphan_vertices = 0;
			remove_global_orphan_skeleton_elements(
				p, &global_orphan_edges, &global_orphan_vertices, &preexisting_orphan_edge_ids);
			out_stats.removed_edges += global_orphan_edges;
			out_stats.removed_vertices += global_orphan_vertices;
		}
		return true;
	}

	template <typename ScoreEvaluator>
	bool compute_skeleton_face_scores(TopologyParameters& p, ScoreEvaluator& evaluator,
							  const std::unordered_set<uint32>& face_ids,
								  std::unordered_map<uint32, Scalar>& face_scores, bool normalize_by_area)
	{
		if (!p.skeleton_ || !p.skeleton_position_)
			return false;
		if (face_ids.empty())
		{
			face_scores.clear();
			return true;
		}

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
			Scalar(0.225000000000000), Scalar(0.132394152788506), Scalar(0.132394152788506),
			Scalar(0.132394152788506), Scalar(0.125939180544827), Scalar(0.125939180544827),
			Scalar(0.125939180544827)};

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

		face_scores.clear();
		face_scores.reserve(face_ids.size());
		std::vector<AdaptiveTriangleTask> pending_tris;
		pending_tris.reserve(face_ids.size() * 2);
		foreach_cell(*p.skeleton_, [&](SkeletonFace f) -> bool {
			const uint32 idf = index_of(*p.skeleton_, f);
			if (idf == INVALID_INDEX)
				return true;
			if (face_ids.find(idf) == face_ids.end())
				return true;
			face_scores[idf] = Scalar(0);
			if (!f.is_valid())
				return true;
			const std::vector<SkeletonVertex> vertices = incident_vertices(*p.skeleton_, f);
			if (vertices.size() < 3)
				return true;
			const Vec3 p0 = value<Vec3>(*p.skeleton_, p.skeleton_position_, vertices[0]);
			for (uint32 i = 1; i + 1 < vertices.size(); ++i)
			{
				const Vec3 p1 = value<Vec3>(*p.skeleton_, p.skeleton_position_, vertices[i]);
				const Vec3 p2 = value<Vec3>(*p.skeleton_, p.skeleton_position_, vertices[i + 1]);
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
			foreach_cell(*p.skeleton_, [&](SkeletonFace f) -> bool {
				const uint32 face_id = index_of(*p.skeleton_, f);
				auto it_score = face_scores.find(face_id);
				if (face_id == INVALID_INDEX || it_score == face_scores.end())
					return true;
				const std::vector<SkeletonVertex> vertices = incident_vertices(*p.skeleton_, f);
				if (vertices.size() < 3)
				{
					it_score->second = Scalar(0);
					return true;
				}
				const Vec3 p0 = value<Vec3>(*p.skeleton_, p.skeleton_position_, vertices[0]);
				Scalar area = Scalar(0);
				for (uint32 i = 1; i + 1 < vertices.size(); ++i)
				{
					const Vec3 p1 = value<Vec3>(*p.skeleton_, p.skeleton_position_, vertices[i]);
					const Vec3 p2 = value<Vec3>(*p.skeleton_, p.skeleton_position_, vertices[i + 1]);
					area += geometry::area(p0, p1, p2);
				}
				it_score->second = (area > Scalar(0)) ? (it_score->second / area) : Scalar(0);
				return true;
			});
		}
		return true;
	}

	std::unordered_set<uint32> collect_current_tet_edge_ids(TopologyParameters& p)
	{
		std::unordered_set<uint32> tet_edge_ids;
		if (!p.skeleton_)
			return tet_edge_ids;
		tet_edge_ids.reserve(nb_cells<SkeletonEdge>(*p.skeleton_));
		for (const auto& kv : p.skeleton_tets_)
		{
			const Tet& tet = kv.second;
			for (uint32 i = 0; i < 4; ++i)
			{
				const SkeletonFace f = tet.faces[i];
				if (!f.is_valid())
					continue;
				for (SkeletonEdge e : incident_edges(*p.skeleton_, f))
				{
					if (!e.is_valid())
						continue;
					const uint32 ide = index_of(*p.skeleton_, e);
					if (ide != INVALID_INDEX)
						tet_edge_ids.insert(ide);
				}
			}
		}
		return tet_edge_ids;
	}

	template <typename ScoreEvaluator>
	bool compute_skeleton_edge_scores_gauss3_normalized(TopologyParameters& p, ScoreEvaluator& evaluator,
										 const std::unordered_set<uint32>& edge_ids,
										std::unordered_map<uint32, Scalar>& edge_scores)
	{
		if (!p.skeleton_ || !p.skeleton_position_)
			return false;
		edge_scores.clear();
		if (edge_ids.empty())
			return true;

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

		edge_scores.reserve(edge_ids.size());
		std::vector<EdgeBatchMeta> edge_batch;
		edge_batch.reserve(edge_ids.size());
		std::vector<Vec3> sample_points;
		sample_points.reserve(edge_ids.size() * 3);

		foreach_cell(*p.skeleton_, [&](SkeletonEdge e) -> bool {
			if (!e.is_valid())
				return true;
			const uint32 ide = index_of(*p.skeleton_, e);
			if (ide == INVALID_INDEX)
				return true;
			if (edge_ids.find(ide) == edge_ids.end())
				return true;
			edge_scores[ide] = Scalar(0);

			const std::vector<SkeletonVertex> verts = incident_vertices(*p.skeleton_, e);
			if (verts.size() != 2)
				return true;

			const Vec3 p0 = value<Vec3>(*p.skeleton_, p.skeleton_position_, verts[0]);
			const Vec3 p1 = value<Vec3>(*p.skeleton_, p.skeleton_position_, verts[1]);
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

	struct TopologyFixScoreCache
	{
		std::unordered_map<uint32, Scalar> edge_scores;
		std::unordered_map<uint32, Scalar> face_scores;
		bool initialized = false;
	};

	template <typename ScoreEvaluator>
	bool initialize_topology_fix_score_cache(
		TopologyParameters& p, ScoreEvaluator& evaluator, TopologyFixScoreCache& cache)
	{
		if (cache.initialized)
			return true;
		const std::unordered_set<uint32> tet_face_ids = collect_current_tet_face_id_whitelist(p);
		const std::unordered_set<uint32> tet_edge_ids = collect_current_tet_edge_ids(p);
		if (!compute_skeleton_edge_scores_gauss3_normalized(p, evaluator, tet_edge_ids, cache.edge_scores))
			return false;
		if (!compute_skeleton_face_scores(p, evaluator, tet_face_ids, cache.face_scores, false))
			return false;
		cache.initialized = true;
		return true;
	}




	void erase_topology_fix_scores_for_removed_faces_and_orphan_edges(
		TopologyParameters& p, TopologyFixScoreCache* cache, const std::unordered_set<uint32>& removed_face_ids,
		const std::vector<SkeletonEdge>& affected_edges)
	{
		if (!cache)
			return;
		for (uint32 face_id : removed_face_ids)
			cache->face_scores.erase(face_id);

		std::unordered_set<uint32> seen_edge_ids;
		seen_edge_ids.reserve(affected_edges.size() * 2 + 1);
		for (SkeletonEdge e : affected_edges)
		{
			if (!e.is_valid())
				continue;
			const uint32 ide = index_of(*p.skeleton_, e);
			if (ide == INVALID_INDEX || !seen_edge_ids.insert(ide).second)
				continue;
			if (incident_faces(*p.skeleton_, e).empty())
				cache->edge_scores.erase(ide);
		}
	}

	enum class EdgeTetDeleteMode
	{
		SimpleTet,
		NonSimpleTet
	};

	struct EdgeTetModeRunStats
	{
		uint32 removed_faces = 0;
		uint32 removed_edges = 0;
		uint32 removed_tets = 0;
	};

	template <typename ScoreEvaluator>
	EdgeTetModeRunStats run_edge_score_tet_mode_topology_fix(
		TopologyParameters& p, EdgeTetDeleteMode mode, TopologyFixScoreCache& score_cache, ScoreEvaluator& evaluator)
	{
		EdgeTetModeRunStats stats;
		if (!p.skeleton_ || !p.incident_tets_)
			return stats;

		if (!initialize_topology_fix_score_cache(p, evaluator, score_cache))
			return stats;
		auto& edge_score_cache = score_cache.edge_scores;
		auto& face_score_cache = score_cache.face_scores;

		std::unordered_map<uint32, std::vector<std::size_t>> face_owner_tets;
		face_owner_tets.reserve(nb_cells<SkeletonFace>(*p.skeleton_));
		std::unordered_map<std::size_t, uint32> deleted_face_count_per_tet;
		deleted_face_count_per_tet.reserve(p.skeleton_tets_.size());
		for (const auto& kv : p.skeleton_tets_)
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
				const uint32 idf = index_of(*p.skeleton_, f);
				if (idf == INVALID_INDEX)
				{
					++initial_deleted_faces;
					continue;
				}
				face_owner_tets[idf].push_back(tet_id);
			}
			deleted_face_count_per_tet[tet_id] = initial_deleted_faces;
		}

		auto get_face_score = [&](uint32 idf) -> Scalar {
			auto it = face_score_cache.find(idf);
			return (it != face_score_cache.end()) ? it->second : Scalar(0);
		};
		auto has_face_budget = [&](uint32 face_id) -> bool {
			auto it_owner = face_owner_tets.find(face_id);
			if (it_owner == face_owner_tets.end())
				return false;
			for (std::size_t tet_id : it_owner->second)
			{
				if (p.skeleton_tets_.find(tet_id) == p.skeleton_tets_.end())
					continue;
				auto it_count = deleted_face_count_per_tet.find(tet_id);
				if (it_count != deleted_face_count_per_tet.end() && it_count->second >= 2)
					return false;
			}
			return true;
		};
		auto account_face_deletion_budget = [&](uint32 face_id) {
			auto it_owner = face_owner_tets.find(face_id);
			if (it_owner == face_owner_tets.end())
				return;
			for (std::size_t tet_id : it_owner->second)
			{
				if (p.skeleton_tets_.find(tet_id) == p.skeleton_tets_.end())
					continue;
				auto it_count = deleted_face_count_per_tet.find(tet_id);
				if (it_count != deleted_face_count_per_tet.end())
					++(it_count->second);
			}
		};
		auto face_belongs_to_tet = [&](const SkeletonFace& f, std::size_t tet_id) -> bool {
			if (!f.is_valid())
				return false;
			const uint32 idf = index_of(*p.skeleton_, f);
			if (idf == INVALID_INDEX)
				return false;
			const auto& in_tets = (*p.incident_tets_)[idf];
			return in_tets.find(tet_id) != in_tets.end();
		};
		auto face_has_degree1_edge = [&](const SkeletonFace& f) -> bool {
			if (!f.is_valid())
				return false;
			for (SkeletonEdge e : incident_edges(*p.skeleton_, f))
			{
				if (!e.is_valid())
					continue;
				if (incident_faces(*p.skeleton_, e).size() == 1)
					return true;
			}
			return false;
		};
		auto face_has_edge_with_tet_face_count_gt2 =
			[&](const SkeletonFace& f, const std::unordered_set<std::size_t>* ignored_tets) -> bool {
			if (!f.is_valid())
				return false;
			for (SkeletonEdge e : incident_edges(*p.skeleton_, f))
			{
				if (!e.is_valid())
					continue;
				size_t tet_face_count = 0;
				for (SkeletonFace ef : incident_faces(*p.skeleton_, e))
				{
					if (!ef.is_valid())
						continue;
					const uint32 idef = index_of(*p.skeleton_, ef);
					if (idef == INVALID_INDEX)
						continue;
					bool has_alive_owner = false;
					for (std::size_t tet_id : (*p.incident_tets_)[idef])
					{
						if (p.skeleton_tets_.find(tet_id) == p.skeleton_tets_.end())
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
		};
		auto tet_has_simple_face = [&](std::size_t tet_id) -> bool {
			auto it_tet = p.skeleton_tets_.find(tet_id);
			if (it_tet == p.skeleton_tets_.end())
				return false;
			const Tet& tet = it_tet->second;
			for (uint32 i = 0; i < 4; ++i)
			{
				const SkeletonFace f = tet.faces[i];
				if (!f.is_valid())
					continue;
				const uint32 idf = index_of(*p.skeleton_, f);
				if (idf == INVALID_INDEX)
					continue;
				const auto& in_tets = (*p.incident_tets_)[idf];
				if (in_tets.find(tet_id) == in_tets.end())
					continue;
				if (in_tets.size() == 1)
					return true;
			}
			return false;
		};
		auto pick_tet_best_edge_and_face = [&](std::size_t tet_id, SkeletonEdge& out_edge, uint32& out_edge_id,
											  Scalar& out_edge_score, SkeletonFace& out_face, uint32& out_face_id,
											  Scalar& out_face_score) -> bool {
			auto it_tet = p.skeleton_tets_.find(tet_id);
			if (it_tet == p.skeleton_tets_.end())
				return false;
			const Tet& tet = it_tet->second;
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
				for (SkeletonEdge e : incident_edges(*p.skeleton_, f))
				{
					if (!e.is_valid())
						continue;
					const uint32 ide = index_of(*p.skeleton_, e);
					if (ide == INVALID_INDEX)
						continue;
					if (!seen_edges.insert(ide).second)
						continue;
					SkeletonFace best_face_on_edge;
					uint32 best_face_on_edge_id = INVALID_INDEX;
					Scalar best_face_on_edge_score = Scalar(0);
					bool edge_has_candidate = false;
					for (const SkeletonFace& ef : incident_faces(*p.skeleton_, e))
					{
						if (!face_belongs_to_tet(ef, tet_id))
							continue;
						const uint32 idf = index_of(*p.skeleton_, ef);
						if (idf == INVALID_INDEX || !has_face_budget(idf))
							continue;
						if (mode == EdgeTetDeleteMode::SimpleTet && (*p.incident_tets_)[idf].size() != 1)
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
		};
		auto remove_faces_collect_tets =
			[&](const std::vector<SkeletonFace>& faces_to_remove, std::unordered_set<std::size_t>& out_tets_to_erase,
				uint32& out_removed_faces, uint32& out_removed_edges,
				std::vector<SkeletonEdge>* out_affected_edges = nullptr) -> bool {
			out_tets_to_erase.clear();
			std::vector<std::pair<uint32, SkeletonFace>> unique_faces;
			unique_faces.reserve(faces_to_remove.size());
			std::unordered_set<uint32> seen_face_ids;
			seen_face_ids.reserve(faces_to_remove.size() * 2 + 1);
			for (const SkeletonFace& f : faces_to_remove)
			{
				if (!f.is_valid())
					continue;
				const uint32 idf = index_of(*p.skeleton_, f);
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
				for (std::size_t tet_id : (*p.incident_tets_)[idf])
					if (p.skeleton_tets_.find(tet_id) != p.skeleton_tets_.end())
						tets_to_remove_local.insert(tet_id);
			}
			if (tets_to_remove_local.empty())
				return false;

			std::vector<SkeletonEdge> affected_edges;
			affected_edges.reserve(unique_faces.size() * 3);
			for (const auto& fpair : unique_faces)
			{
				const SkeletonFace f = fpair.second;
				for (SkeletonEdge e : incident_edges(*p.skeleton_, f))
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
				const uint32 idf_now = index_of(*p.skeleton_, f);
				if (idf_now == INVALID_INDEX)
					continue;
				account_face_deletion_budget(idf_now);
				removed_face_ids.insert(idf_now);
				remove_face(*p.skeleton_, f);
				++out_removed_faces;
			}

			erase_topology_fix_scores_for_removed_faces_and_orphan_edges(p, &score_cache, removed_face_ids, affected_edges);
			out_removed_edges += remove_orphan_edges_from_removed_face_edges(p, affected_edges);
			out_tets_to_erase = std::move(tets_to_remove_local);
			return true;
		};
		auto finalize_erase_tets = [&](const std::unordered_set<std::size_t>& tets_to_erase, uint32& out_removed_tets) {
			for (std::size_t tet_id : tets_to_erase)
			{
				auto it_tet = p.skeleton_tets_.find(tet_id);
				if (it_tet == p.skeleton_tets_.end())
					continue;
				const Tet old_tet = it_tet->second;
				for (uint32 i = 0; i < 4; ++i)
				{
					const SkeletonFace f = old_tet.faces[i];
					if (!f.is_valid())
						continue;
					const uint32 idf = index_of(*p.skeleton_, f);
					if (idf == INVALID_INDEX)
						continue;
					(*p.incident_tets_)[idf].erase(tet_id);
				}
				p.skeleton_tets_.erase(it_tet);
				++out_removed_tets;
			}
		};
		struct StepCandidate
		{
			std::size_t tet_id = 0;
			SkeletonEdge edge;
			uint32 edge_id = INVALID_INDEX;
			Scalar edge_score = Scalar(0);
			SkeletonFace first_face;
			uint32 first_face_id = INVALID_INDEX;
			Scalar first_face_score = Scalar(0);
		};
		struct FaceChoice
		{
			SkeletonFace face;
			uint32 face_id = INVALID_INDEX;
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

		uint32 removed_faces = 0;
		uint32 removed_edges = 0;
		uint32 removed_tets = 0;
		std::unordered_map<std::size_t, uint32> tet_versions;
		tet_versions.reserve(p.skeleton_tets_.size() * 2 + 1);
		TetCandidateQueue candidate_queue;

		auto build_step_candidate = [&](std::size_t tet_id, StepCandidate& out_candidate) -> bool {
			const bool is_simple_tet = tet_has_simple_face(tet_id);
			if (mode == EdgeTetDeleteMode::SimpleTet && !is_simple_tet)
				return false;
			if (mode == EdgeTetDeleteMode::NonSimpleTet && is_simple_tet)
				return false;

			SkeletonEdge best_edge_local;
			uint32 best_edge_local_id = INVALID_INDEX;
			Scalar best_edge_local_score = Scalar(0);
			SkeletonFace first_face_local;
			uint32 first_face_local_id = INVALID_INDEX;
			Scalar first_face_local_score = Scalar(0);
			if (!pick_tet_best_edge_and_face(
					tet_id, best_edge_local, best_edge_local_id, best_edge_local_score, first_face_local,
					first_face_local_id, first_face_local_score))
				return false;

			out_candidate.tet_id = tet_id;
			out_candidate.edge = best_edge_local;
			out_candidate.edge_id = best_edge_local_id;
			out_candidate.edge_score = best_edge_local_score;
			out_candidate.first_face = first_face_local;
			out_candidate.first_face_id = first_face_local_id;
			out_candidate.first_face_score = first_face_local_score;
			return true;
		};
		auto queue_tet_candidate = [&](std::size_t tet_id, bool bump_version) {
			if (p.skeleton_tets_.find(tet_id) == p.skeleton_tets_.end())
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
			candidate_queue.push(
				{candidate.tet_id, candidate.edge_id, candidate.first_face_id, candidate.edge_score,
				 candidate.first_face_score, version});
		};
		auto collect_dirty_tets_from_edges =
			[&](const std::vector<SkeletonEdge>& affected_edges, std::unordered_set<std::size_t>& dirty_tets) {
				dirty_tets.clear();
				for (SkeletonEdge e : affected_edges)
				{
					if (!e.is_valid())
						continue;
					for (SkeletonFace f : incident_faces(*p.skeleton_, e))
					{
						if (!f.is_valid())
							continue;
						const uint32 idf = index_of(*p.skeleton_, f);
						if (idf == INVALID_INDEX)
							continue;
						for (std::size_t tet_id : (*p.incident_tets_)[idf])
							if (p.skeleton_tets_.find(tet_id) != p.skeleton_tets_.end())
								dirty_tets.insert(tet_id);
					}
				}
			};

		for (const auto& kv : p.skeleton_tets_)
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
				if (p.skeleton_tets_.find(entry.tet_id) == p.skeleton_tets_.end())
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

				const bool matches_cached =
					(refreshed.edge_id == entry.edge_id) && (refreshed.first_face_id == entry.face_id) &&
					(std::abs(refreshed.edge_score - entry.edge_score) <= edge_eps) &&
					(std::abs(refreshed.first_face_score - entry.face_score) <= face_eps);
				if (!matches_cached)
				{
					++(it_version->second);
					candidate_queue.push(
						{refreshed.tet_id, refreshed.edge_id, refreshed.first_face_id, refreshed.edge_score,
						 refreshed.first_face_score, it_version->second});
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
			for (std::size_t tet_id : (*p.incident_tets_)[best.first_face_id])
				if (p.skeleton_tets_.find(tet_id) != p.skeleton_tets_.end())
					first_removed_tets.insert(tet_id);

			uint32 step_removed_faces = 0;
			uint32 step_removed_edges = 0;
			uint32 step_removed_tets = 0;
			std::unordered_set<std::size_t> step_tets_to_erase;
			std::unordered_set<std::size_t> first_tets_to_erase;
			std::vector<SkeletonEdge> step_affected_edges;
			if (!remove_faces_collect_tets(
					{best.first_face}, first_tets_to_erase, step_removed_faces, step_removed_edges, &step_affected_edges))
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
				auto it_tet = p.skeleton_tets_.find(tet_id);
				if (it_tet == p.skeleton_tets_.end())
					continue;
				const Tet& tet = it_tet->second;
				for (uint32 i = 0; i < 4; ++i)
				{
					const SkeletonFace f = tet.faces[i];
					const uint32 idf = index_of(*p.skeleton_, f);
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
					created_deg1_faces.push_back({f, idf, get_face_score(idf)});
				}
			}
			std::sort(created_deg1_faces.begin(), created_deg1_faces.end(),
					  [](const FaceChoice& a, const FaceChoice& b) { return a.score > b.score; });
			if (mode == EdgeTetDeleteMode::SimpleTet && created_deg1_faces.size() > 1)
				created_deg1_faces.resize(1);

			for (const FaceChoice& fc : created_deg1_faces)
			{
				if (!fc.face.is_valid())
					continue;
				const uint32 idf_now = index_of(*p.skeleton_, fc.face);
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
				if (remove_faces_collect_tets(
						{fc.face}, second_tets_to_erase, rf, re, &step_affected_edges))
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

	struct DenseFiveVertexDetectionResult
	{
		std::vector<std::array<uint32, 5>> independent_exact_k5_regions;
		std::vector<std::array<uint32, 5>> independent_k5_minus_1_regions;
	};

	DenseFiveVertexDetectionResult detect_k5_and_k5_minus_1_cells(TopologyParameters& p)
	{
		DenseFiveVertexDetectionResult result;
		if (!p.skeleton_)
			return result;

		auto make_edge_key = [](uint32 a, uint32 b) -> uint64 {
			if (a > b)
				std::swap(a, b);
			return (uint64(a) << 32) | uint64(b);
		};
		std::vector<std::array<uint32, 5>> exact_k5_regions;
		std::vector<std::array<uint32, 5>> k5_minus_1_regions;

		std::unordered_map<uint32, std::unordered_set<uint32>> adjacency;
		adjacency.reserve(nb_cells<SkeletonVertex>(*p.skeleton_));
		std::unordered_map<uint64, uint32> edge_pair_to_id;
		edge_pair_to_id.reserve(nb_cells<SkeletonEdge>(*p.skeleton_) * 2 + 1);

		foreach_cell(*p.skeleton_, [&](SkeletonEdge e) -> bool {
			const uint32 ide = index_of(*p.skeleton_, e);
			if (ide == INVALID_INDEX)
				return true;
			const std::vector<SkeletonVertex> vv = incident_vertices(*p.skeleton_, e);
			if (vv.size() != 2)
				return true;
			const uint32 a = index_of(*p.skeleton_, vv[0]);
			const uint32 b = index_of(*p.skeleton_, vv[1]);
			if (a == INVALID_INDEX || b == INVALID_INDEX || a == b)
				return true;
			adjacency[a].insert(b);
			adjacency[b].insert(a);
			edge_pair_to_id[make_edge_key(a, b)] = ide;
			return true;
		});

		std::set<std::array<uint32, 5>> seen_region_keys;
		for (const auto& tet_entry : p.skeleton_tets_)
		{
			const Tet& tet = tet_entry.second;
			std::array<uint32, 4> base_vertices = {INVALID_INDEX, INVALID_INDEX, INVALID_INDEX, INVALID_INDEX};
			uint32 base_vertex_count = 0;
			for (uint32 face_index = 0; face_index < 4; ++face_index)
			{
				const SkeletonFace face = tet.faces[face_index];
				if (!face.is_valid())
					continue;
				for (SkeletonVertex vertex : incident_vertices(*p.skeleton_, face))
				{
					const uint32 vertex_id = index_of(*p.skeleton_, vertex);
					if (vertex_id == INVALID_INDEX ||
						std::find(base_vertices.begin(), base_vertices.begin() + base_vertex_count, vertex_id) !=
							base_vertices.begin() + base_vertex_count)
						continue;
					if (base_vertex_count < 4)
						base_vertices[base_vertex_count++] = vertex_id;
				}
			}
			if (base_vertex_count != 4)
				continue;
			std::sort(base_vertices.begin(), base_vertices.end());

			std::unordered_map<uint32, uint32> fifth_vertex_connection_count;
			for (uint32 base_vertex : base_vertices)
			{
				const auto adjacency_it = adjacency.find(base_vertex);
				if (adjacency_it == adjacency.end())
					continue;
				for (uint32 neighbor : adjacency_it->second)
				{
					if (!std::binary_search(base_vertices.begin(), base_vertices.end(), neighbor))
						++fifth_vertex_connection_count[neighbor];
				}
			}

			for (const auto& candidate : fifth_vertex_connection_count)
			{
				if (candidate.second < 3)
					continue;

				std::array<uint32, 5> region = {
					base_vertices[0], base_vertices[1], base_vertices[2], base_vertices[3], candidate.first};
				std::sort(region.begin(), region.end());
				if (!seen_region_keys.insert(region).second)
					continue;

				const bool is_exact_k5 = candidate.second == 4;
				if (is_exact_k5)
					exact_k5_regions.push_back(region);
				else
					k5_minus_1_regions.push_back(region);
			}
		}

		std::unordered_map<uint64, uint32> edge_region_count;
		edge_region_count.reserve((exact_k5_regions.size() + k5_minus_1_regions.size()) * 10);
		auto count_region_edges = [&](const std::array<uint32, 5>& region) {
			for (uint32 i = 0; i < 5; ++i)
			{
				for (uint32 j = i + 1; j < 5; ++j)
				{
					const uint64 edge_key = make_edge_key(region[i], region[j]);
					if (edge_pair_to_id.find(edge_key) != edge_pair_to_id.end())
						++edge_region_count[edge_key];
				}
			}
		};
		for (const auto& region : exact_k5_regions)
			count_region_edges(region);
		for (const auto& region : k5_minus_1_regions)
			count_region_edges(region);

		auto append_independent_region = [&](const std::array<uint32, 5>& region, bool is_exact_k5) {
			for (uint32 i = 0; i < 5; ++i)
			{
				for (uint32 j = i + 1; j < 5; ++j)
				{
					const uint64 edge_key = make_edge_key(region[i], region[j]);
					const auto count_it = edge_region_count.find(edge_key);
					if (count_it != edge_region_count.end() && count_it->second > 1)
						return;
				}
			}

			if (is_exact_k5)
				result.independent_exact_k5_regions.push_back(region);
			else
				result.independent_k5_minus_1_regions.push_back(region);
		};
		for (const auto& region : exact_k5_regions)
			append_independent_region(region, true);
		for (const auto& region : k5_minus_1_regions)
			append_independent_region(region, false);

		return result;
	}

	struct BoundaryTetPrepassStats
	{
		uint32 removed_faces = 0;
		uint32 removed_tets = 0;
		uint32 removed_edges = 0;
	};

	BoundaryTetPrepassStats run_boundary_tet_face_deletion(
		TopologyParameters& p, TopologyFixScoreCache* score_cache = nullptr)
	{
		BoundaryTetPrepassStats stats;

		struct BoundaryTetCandidate
		{
			std::size_t tet_id = std::size_t(-1);
			SkeletonEdge edge;
			SkeletonFace f0;
			SkeletonFace f1;
		};

		std::vector<std::size_t> tet_ids;
		tet_ids.reserve(p.skeleton_tets_.size());
		for (const auto& kv : p.skeleton_tets_)
			tet_ids.push_back(kv.first);

		std::vector<BoundaryTetCandidate> candidates;
		candidates.reserve(tet_ids.size());

		// Phase 1: scan all current tets once and collect boundary-tet deletion candidates.
		for (std::size_t tet_id : tet_ids)
		{
			auto it_tet = p.skeleton_tets_.find(tet_id);
			if (it_tet == p.skeleton_tets_.end())
				continue;
			const Tet tet = it_tet->second;

			bool is_boundary_tet = false;
			std::unordered_set<uint32> tet_face_ids;
			tet_face_ids.reserve(4);
			for (uint32 i = 0; i < 4; ++i)
			{
				const SkeletonFace f = tet.faces[i];
				const uint32 idf = index_of(*p.skeleton_, f);
				tet_face_ids.insert(idf);

				uint32 deg2_edge_count = 0;
				for (SkeletonEdge e : incident_edges(*p.skeleton_, f))
				{
					if (incident_faces(*p.skeleton_, e).size() == 2)
						++deg2_edge_count;
				}
				if (deg2_edge_count >= 2)
				{
					is_boundary_tet = true;
					break;
				}
			}
			if (!is_boundary_tet || tet_face_ids.empty())
				continue;

			struct BoundaryEdgeChoice
			{
				SkeletonEdge edge;
				SkeletonFace f0;
				SkeletonFace f1;
				Scalar length = Scalar(-1);
			};
			BoundaryEdgeChoice best_choice;
			std::unordered_set<uint32> seen_edges;
			seen_edges.reserve(8);

			for (uint32 i = 0; i < 4; ++i)
			{
				const SkeletonFace f = tet.faces[i];
				for (SkeletonEdge e : incident_edges(*p.skeleton_, f))
				{
					const uint32 ide = index_of(*p.skeleton_, e);
					if (!seen_edges.insert(ide).second)
						continue;

					const auto in_faces = incident_faces(*p.skeleton_, e);
					if (in_faces.size() != 2)
						continue;
					const SkeletonFace f0 = in_faces[0];
					const SkeletonFace f1 = in_faces[1];
					const uint32 idf0 = index_of(*p.skeleton_, f0);
					const uint32 idf1 = index_of(*p.skeleton_, f1);
					const std::vector<SkeletonVertex> edge_vertices = incident_vertices(*p.skeleton_, e);
					if (edge_vertices.size() != 2)
						continue;
					const Vec3 p0 = value<Vec3>(*p.skeleton_, p.skeleton_position_, edge_vertices[0]);
					const Vec3 p1 = value<Vec3>(*p.skeleton_, p.skeleton_position_, edge_vertices[1]);
					const Scalar edge_len = (p1 - p0).norm();
					if (edge_len > best_choice.length)
					{
						best_choice.edge = e;
						best_choice.f0 = f0;
						best_choice.f1 = f1;
						best_choice.length = edge_len;
					}
				}
			}

			if (best_choice.length >= Scalar(0))
				candidates.push_back({tet_id, best_choice.edge, best_choice.f0, best_choice.f1});
		}

		// Phase 2: delete using precomputed candidates only (no boundary re-check during deletion).
		for (const BoundaryTetCandidate& cand : candidates)
		{
			auto it_remove_tet = p.skeleton_tets_.find(cand.tet_id);
			if (it_remove_tet == p.skeleton_tets_.end())
				continue;
			const Tet old_tet = it_remove_tet->second;

			const std::array<SkeletonFace, 2> faces_to_remove = {cand.f0, cand.f1};
			std::vector<SkeletonEdge> affected_edges;
			affected_edges.reserve(6);
			for (const SkeletonFace& f_del : faces_to_remove)
			{
				if (!f_del.is_valid())
					continue;
				for (SkeletonEdge e : incident_edges(*p.skeleton_, f_del))
					affected_edges.push_back(e);
			}
			std::unordered_set<uint32> removed_face_ids;
			removed_face_ids.reserve(faces_to_remove.size() * 2 + 1);
			for (const SkeletonFace& f_del : faces_to_remove)
			{
				if (const uint32 idf_del = index_of(*p.skeleton_, f_del); idf_del != INVALID_INDEX)
					removed_face_ids.insert(idf_del);
				remove_face(*p.skeleton_, f_del);
				++stats.removed_faces;
			}
			erase_topology_fix_scores_for_removed_faces_and_orphan_edges(p, score_cache, removed_face_ids, affected_edges);
			if (cand.edge.is_valid())
			{
				const auto edge_in_faces = incident_faces(*p.skeleton_, cand.edge);
				if (edge_in_faces.empty())
				{
					remove_edge(*p.skeleton_, cand.edge);
					++stats.removed_edges;
				}
			}
			p.skeleton_tets_.erase(it_remove_tet);
			++stats.removed_tets;
			for (uint32 i = 0; i < 4; ++i)
			{
				const SkeletonFace f = old_tet.faces[i];
				const uint32 idf = index_of(*p.skeleton_, f);
				(*p.incident_tets_)[idf].erase(cand.tet_id);
			}
		}

		return stats;
	}



	std::unordered_set<uint32> collect_current_tet_face_id_whitelist(TopologyParameters& p)
	{
		std::unordered_set<uint32> tet_face_ids;
		tet_face_ids.reserve(nb_cells<SkeletonFace>(*p.skeleton_));
		for (const auto& kv : p.skeleton_tets_)
		{
			const Tet& tet = kv.second;
			for (uint32 i = 0; i < 4; ++i)
			{
				const SkeletonFace f = tet.faces[i];
				if (!f.is_valid())
					continue;
				const uint32 idf = index_of(*p.skeleton_, f);
				if (idf != INVALID_INDEX)
					tet_face_ids.insert(idf);
			}
		}
		return tet_face_ids;
	}

	struct IndependentDenseFiveVertexRepairSeed
	{
		std::array<SphereVertex, 5> spheres;
	};

	struct IndependentDenseFiveVertexRepairStats
	{
		uint32 repaired_regions = 0;
		uint32 added_faces = 0;
	};

	bool is_live_sphere_vertex(TopologyParameters& p, const SphereVertex& sphere, uint32& sphere_index)
	{
		sphere_index = INVALID_INDEX;
		if (!p.spheres_ || !sphere.is_valid() || sphere.dart_.index_ >= p.spheres_->darts_.maximum_index())
			return false;
		sphere_index = index_of(*p.spheres_, sphere);
		return sphere_index != INVALID_INDEX;
	}

	bool is_live_skeleton_vertex(TopologyParameters& p, const SkeletonVertex& vertex, uint32& vertex_index)
	{
		vertex_index = INVALID_INDEX;
		if (!p.skeleton_ || !vertex.is_valid())
			return false;
		vertex_index = index_of(*p.skeleton_, vertex);
		return vertex_index != INVALID_INDEX;
	}

	std::vector<IndependentDenseFiveVertexRepairSeed> collect_independent_dense_five_vertex_repair_seeds(
		TopologyParameters& p, const DenseFiveVertexDetectionResult& detection)
	{
		std::vector<IndependentDenseFiveVertexRepairSeed> seeds;
		seeds.reserve(detection.independent_exact_k5_regions.size() +
					  detection.independent_k5_minus_1_regions.size());
		if (!p.skeleton_ || !p.spheres_ || !p.skeleton_source_sphere_)
			return seeds;

		auto append_regions = [&](const std::vector<std::array<uint32, 5>>& regions) {
			for (const auto& region : regions)
			{
				IndependentDenseFiveVertexRepairSeed seed;
				bool valid = true;
				for (uint32 i = 0; i < 5; ++i)
				{
					const SkeletonVertex vertex = of_index<SkeletonVertex>(*p.skeleton_, region[i]);
					uint32 current_vertex_index = INVALID_INDEX;
					if (!is_live_skeleton_vertex(p, vertex, current_vertex_index) || current_vertex_index != region[i])
					{
						valid = false;
						break;
					}
					const SphereVertex sphere = value<SphereVertex>(*p.skeleton_, p.skeleton_source_sphere_, vertex);
					uint32 sphere_index = INVALID_INDEX;
					if (!is_live_sphere_vertex(p, sphere, sphere_index))
					{
						valid = false;
						break;
					}
					seed.spheres[i] = sphere;
				}
				if (valid)
					seeds.push_back(seed);
			}
		};
		append_regions(detection.independent_exact_k5_regions);
		append_regions(detection.independent_k5_minus_1_regions);
		return seeds;
	}

	bool find_skeleton_edge_between(TopologyParameters& p, const SkeletonVertex& a, const SkeletonVertex& b, SkeletonEdge& out_edge)
	{
		out_edge = SkeletonEdge();
		if (!p.skeleton_ || !a.is_valid() || !b.is_valid() || a == b)
			return false;
		for (SkeletonEdge edge : incident_edges(*p.skeleton_, a))
		{
			if (!edge.is_valid())
				continue;
			const std::vector<SkeletonVertex> edge_vertices = incident_vertices(*p.skeleton_, edge);
			if (edge_vertices.size() == 2 &&
				((edge_vertices[0] == a && edge_vertices[1] == b) ||
				 (edge_vertices[0] == b && edge_vertices[1] == a)))
			{
				out_edge = edge;
				return true;
			}
		}
		return false;
	}

	bool face_is_inside_vertex_set(TopologyParameters& p, const SkeletonFace& face,
								   const std::unordered_set<uint32>& vertex_ids)
	{
		if (!p.skeleton_ || !face.is_valid())
			return false;
		const std::vector<SkeletonVertex> vertices = incident_vertices(*p.skeleton_, face);
		if (vertices.size() != 3)
			return false;
		for (SkeletonVertex vertex : vertices)
		{
			const uint32 vertex_id = index_of(*p.skeleton_, vertex);
			if (vertex_id == INVALID_INDEX || vertex_ids.find(vertex_id) == vertex_ids.end())
				return false;
		}
		return true;
	}

	enum class IndependentDenseFiveVertexRepairResult
	{
		Repaired,
		MissingVertices,
		NoNonManifoldEdges,
		AdjacentNonManifoldEdges,
		UnsafeExternalFaces,
		DegenerateGeometry
	};

	IndependentDenseFiveVertexRepairResult repair_independent_dense_five_vertex_region(
		TopologyParameters& p, const IndependentDenseFiveVertexRepairSeed& seed, uint32& out_added_faces)
	{
		out_added_faces = 0;
		if (!p.spheres_ || !p.spheres_position_ || !p.spheres_radius_ || !p.spheres_neighbor_clusters_ ||
			!p.spheres_skeleton_vertex_ || !p.skeleton_ || !p.skeleton_position_ || !p.skeleton_radius_ ||
			!p.skeleton_source_sphere_ || !p.incident_tets_)
			return IndependentDenseFiveVertexRepairResult::MissingVertices;

		std::array<uint32, 5> sphere_ids;
		std::array<SkeletonVertex, 5> vertices;
		std::array<uint32, 5> vertex_ids;
		std::array<Vec3, 5> positions;
		Scalar average_radius = Scalar(0);
		Vec3 centroid = Vec3::Zero();
		std::unordered_set<uint32> region_vertex_ids;
		region_vertex_ids.reserve(8);
		for (uint32 i = 0; i < 5; ++i)
		{
			if (!is_live_sphere_vertex(p, seed.spheres[i], sphere_ids[i]))
				return IndependentDenseFiveVertexRepairResult::MissingVertices;
			vertices[i] = (*p.spheres_skeleton_vertex_)[sphere_ids[i]];
			if (!is_live_skeleton_vertex(p, vertices[i], vertex_ids[i]))
				return IndependentDenseFiveVertexRepairResult::MissingVertices;
			if (value<SphereVertex>(*p.skeleton_, p.skeleton_source_sphere_, vertices[i]) != seed.spheres[i])
				return IndependentDenseFiveVertexRepairResult::MissingVertices;
			positions[i] = (*p.spheres_position_)[sphere_ids[i]];
			centroid += positions[i];
			average_radius += (*p.spheres_radius_)[sphere_ids[i]];
			region_vertex_ids.insert(vertex_ids[i]);
		}
		if (region_vertex_ids.size() != 5)
			return IndependentDenseFiveVertexRepairResult::MissingVertices;
		centroid /= Scalar(5);
		average_radius /= Scalar(5);
		if (!centroid.allFinite() || !std::isfinite(static_cast<double>(average_radius)) || average_radius <= Scalar(0))
			return IndependentDenseFiveVertexRepairResult::DegenerateGeometry;

		auto make_edge_key = [](uint32 a, uint32 b) -> uint64 {
			if (a > b)
				std::swap(a, b);
			return (uint64(a) << 32) | uint64(b);
		};
		std::unordered_map<uint64, SkeletonEdge> region_edges;
		region_edges.reserve(16);
		std::unordered_set<uint32> region_nm_edge_ids;
		region_nm_edge_ids.reserve(16);
		for (uint32 i = 0; i < 5; ++i)
		{
			for (uint32 j = i + 1; j < 5; ++j)
			{
				SkeletonEdge edge;
				if (!find_skeleton_edge_between(p, vertices[i], vertices[j], edge))
					continue;
				region_edges[make_edge_key(vertex_ids[i], vertex_ids[j])] = edge;
				if (incident_faces(*p.skeleton_, edge).size() > 2)
				{
					const uint32 edge_id = index_of(*p.skeleton_, edge);
					if (edge_id != INVALID_INDEX)
						region_nm_edge_ids.insert(edge_id);
				}
			}
		}
		if (region_nm_edge_ids.empty())
			return IndependentDenseFiveVertexRepairResult::NoNonManifoldEdges;

		for (uint32 nm_edge_id : region_nm_edge_ids)
		{
			const SkeletonEdge nm_edge = of_index<SkeletonEdge>(*p.skeleton_, nm_edge_id);
			if (!nm_edge.is_valid())
				continue;
			for (SkeletonVertex endpoint : incident_vertices(*p.skeleton_, nm_edge))
			{
				for (SkeletonEdge adjacent_edge : incident_edges(*p.skeleton_, endpoint))
				{
					if (!adjacent_edge.is_valid() || incident_faces(*p.skeleton_, adjacent_edge).size() <= 2)
						continue;
					const uint32 adjacent_edge_id = index_of(*p.skeleton_, adjacent_edge);
					if (adjacent_edge_id != INVALID_INDEX &&
						region_nm_edge_ids.find(adjacent_edge_id) == region_nm_edge_ids.end())
						return IndependentDenseFiveVertexRepairResult::AdjacentNonManifoldEdges;
				}
			}
		}

		Eigen::Matrix<Scalar, 3, 3> covariance = Eigen::Matrix<Scalar, 3, 3>::Zero();
		for (const Vec3& position : positions)
		{
			const Vec3 delta = position - centroid;
			covariance += delta * delta.transpose();
		}
		Eigen::SelfAdjointEigenSolver<Eigen::Matrix<Scalar, 3, 3>> eigen_solver(covariance);
		if (eigen_solver.info() != Eigen::Success ||
			eigen_solver.eigenvalues()[1] <= std::numeric_limits<Scalar>::epsilon())
			return IndependentDenseFiveVertexRepairResult::DegenerateGeometry;
		const Vec3 axis_u = eigen_solver.eigenvectors().col(2);
		const Vec3 axis_v = eigen_solver.eigenvectors().col(1);
		std::array<uint32, 5> cyclic_order = {0, 1, 2, 3, 4};
		std::sort(cyclic_order.begin(), cyclic_order.end(), [&](uint32 lhs, uint32 rhs) {
			const Vec3 lhs_delta = positions[lhs] - centroid;
			const Vec3 rhs_delta = positions[rhs] - centroid;
			const Scalar lhs_angle = std::atan2(lhs_delta.dot(axis_v), lhs_delta.dot(axis_u));
			const Scalar rhs_angle = std::atan2(rhs_delta.dot(axis_v), rhs_delta.dot(axis_u));
			return lhs_angle < rhs_angle;
		});
		Scalar max_radial_squared = Scalar(0);
		for (const Vec3& position : positions)
			max_radial_squared = std::max(max_radial_squared, (position - centroid).squaredNorm());
		const Scalar fan_area_epsilon =
			std::max(std::numeric_limits<Scalar>::epsilon(), max_radial_squared * Scalar(1e-10));
		for (uint32 i = 0; i < 5; ++i)
		{
			const Vec3 radial_a = positions[cyclic_order[i]] - centroid;
			const Vec3 radial_b = positions[cyclic_order[(i + 1) % 5]] - centroid;
			if (radial_a.cross(radial_b).norm() <= fan_area_epsilon)
				return IndependentDenseFiveVertexRepairResult::DegenerateGeometry;
		}

		std::unordered_set<uint64> cycle_edge_keys;
		cycle_edge_keys.reserve(8);
		for (uint32 i = 0; i < 5; ++i)
		{
			const uint32 a = cyclic_order[i];
			const uint32 b = cyclic_order[(i + 1) % 5];
			cycle_edge_keys.insert(make_edge_key(vertex_ids[a], vertex_ids[b]));
		}

		std::unordered_set<uint32> internal_face_ids;
		internal_face_ids.reserve(16);
		foreach_cell(*p.skeleton_, [&](SkeletonFace face) -> bool {
			if (face_is_inside_vertex_set(p, face, region_vertex_ids))
			{
				const uint32 face_id = index_of(*p.skeleton_, face);
				if (face_id != INVALID_INDEX)
					internal_face_ids.insert(face_id);
			}
			return true;
		});
		if (internal_face_ids.empty())
			return IndependentDenseFiveVertexRepairResult::UnsafeExternalFaces;

		for (const auto& edge_entry : region_edges)
		{
			uint32 external_face_count = 0;
			for (SkeletonFace face : incident_faces(*p.skeleton_, edge_entry.second))
			{
				const uint32 face_id = index_of(*p.skeleton_, face);
				if (face_id == INVALID_INDEX || internal_face_ids.find(face_id) == internal_face_ids.end())
					++external_face_count;
			}
			if (cycle_edge_keys.find(edge_entry.first) == cycle_edge_keys.end())
			{
				if (external_face_count != 0)
					return IndependentDenseFiveVertexRepairResult::UnsafeExternalFaces;
			}
			else if (external_face_count > 1)
			{
				return IndependentDenseFiveVertexRepairResult::UnsafeExternalFaces;
			}
		}

		SkeletonFaceDeletionStats deletion_stats;
		SkeletonFaceDeletionOptions deletion_options;
		deletion_options.remove_incident_orphan_edges_immediately = false;
		deletion_options.remove_global_orphan_elements_after_batch = false;
		if (!delete_skeleton_face_id_set(p, internal_face_ids, deletion_stats, deletion_options))
			return IndependentDenseFiveVertexRepairResult::UnsafeExternalFaces;

		for (const auto& edge_entry : region_edges)
		{
			if (cycle_edge_keys.find(edge_entry.first) != cycle_edge_keys.end())
				continue;
			if (edge_entry.second.is_valid() && incident_faces(*p.skeleton_, edge_entry.second).empty())
				remove_edge(*p.skeleton_, edge_entry.second);
		}

		for (uint32 i = 0; i < 5; ++i)
		{
			for (uint32 j = i + 1; j < 5; ++j)
			{
				const uint64 edge_key = make_edge_key(vertex_ids[i], vertex_ids[j]);
				const bool is_cycle_edge = cycle_edge_keys.find(edge_key) != cycle_edge_keys.end();
				if (is_cycle_edge)
				{
					(*p.spheres_neighbor_clusters_)[sphere_ids[i]].insert(seed.spheres[j]);
					(*p.spheres_neighbor_clusters_)[sphere_ids[j]].insert(seed.spheres[i]);
				}
				else
				{
					(*p.spheres_neighbor_clusters_)[sphere_ids[i]].erase(seed.spheres[j]);
					(*p.spheres_neighbor_clusters_)[sphere_ids[j]].erase(seed.spheres[i]);
				}
			}
		}

		const SphereVertex center_sphere = p.spheres_optimizer_->add_sphere(centroid, average_radius);
		if (!center_sphere.is_valid())
			return IndependentDenseFiveVertexRepairResult::MissingVertices;
		const uint32 center_sphere_id = index_of(*p.spheres_, center_sphere);
		metrics_.added_spheres.push_back(center_sphere);

		const SkeletonVertex center_vertex = add_vertex(*p.skeleton_);
		value<Vec3>(*p.skeleton_, p.skeleton_position_, center_vertex) = centroid;
		value<Scalar>(*p.skeleton_, p.skeleton_radius_, center_vertex) = average_radius;
		value<SphereVertex>(*p.skeleton_, p.skeleton_source_sphere_, center_vertex) = center_sphere;
		(*p.spheres_skeleton_vertex_)[center_sphere_id] = center_vertex;

		std::array<SkeletonEdge, 5> cycle_edges;
		std::array<SkeletonEdge, 5> spoke_edges;
		for (uint32 i = 0; i < 5; ++i)
		{
			const uint32 a = cyclic_order[i];
			const uint32 b = cyclic_order[(i + 1) % 5];
			if (!find_skeleton_edge_between(p, vertices[a], vertices[b], cycle_edges[i]))
				cycle_edges[i] = add_edge(*p.skeleton_, vertices[a], vertices[b]);
			spoke_edges[i] = add_edge(*p.skeleton_, center_vertex, vertices[a]);
			(*p.spheres_neighbor_clusters_)[center_sphere_id].insert(seed.spheres[a]);
			(*p.spheres_neighbor_clusters_)[sphere_ids[a]].insert(center_sphere);
		}

		for (uint32 i = 0; i < 5; ++i)
		{
			std::vector<SkeletonEdge> face_edges = {cycle_edges[i], spoke_edges[i], spoke_edges[(i + 1) % 5]};
			const SkeletonFace new_face = add_face(*p.skeleton_, face_edges);
			value<std::set<std::size_t>>(*p.skeleton_, p.incident_tets_, new_face).clear();
			++out_added_faces;
		}

		return IndependentDenseFiveVertexRepairResult::Repaired;
	}

	IndependentDenseFiveVertexRepairStats repair_isolated_independent_dense_five_vertex_regions(
		TopologyParameters& p, const std::vector<IndependentDenseFiveVertexRepairSeed>& seeds)
	{
		IndependentDenseFiveVertexRepairStats stats;
		std::vector<IndependentDenseFiveVertexRepairSeed> live_seeds;
		live_seeds.reserve(seeds.size());
		for (const auto& seed : seeds)
		{
			bool all_vertices_live = true;
			std::unordered_set<uint32> unique_skeleton_vertex_ids;
			unique_skeleton_vertex_ids.reserve(8);
			for (const SphereVertex& sphere : seed.spheres)
			{
				uint32 sphere_id = INVALID_INDEX;
				if (!is_live_sphere_vertex(p, sphere, sphere_id))
				{
					all_vertices_live = false;
					break;
				}
				const SkeletonVertex vertex = (*p.spheres_skeleton_vertex_)[sphere_id];
				uint32 vertex_id = INVALID_INDEX;
				if (!is_live_skeleton_vertex(p, vertex, vertex_id) ||
					value<SphereVertex>(*p.skeleton_, p.skeleton_source_sphere_, vertex) != sphere)
				{
					all_vertices_live = false;
					break;
				}
				unique_skeleton_vertex_ids.insert(vertex_id);
			}
			if (all_vertices_live && unique_skeleton_vertex_ids.size() == 5)
				live_seeds.push_back(seed);
		}

		for (const auto& seed : live_seeds)
		{
			uint32 added_faces = 0;
			const IndependentDenseFiveVertexRepairResult result =
				repair_independent_dense_five_vertex_region(p, seed, added_faces);
			switch (result)
			{
			case IndependentDenseFiveVertexRepairResult::Repaired:
				++stats.repaired_regions;
				stats.added_faces += added_faces;
				break;
			case IndependentDenseFiveVertexRepairResult::NoNonManifoldEdges:
			case IndependentDenseFiveVertexRepairResult::MissingVertices:
			case IndependentDenseFiveVertexRepairResult::AdjacentNonManifoldEdges:
			case IndependentDenseFiveVertexRepairResult::UnsafeExternalFaces:
			case IndependentDenseFiveVertexRepairResult::DegenerateGeometry:
				break;
			}
		}

		if (stats.repaired_regions != 0)
			compute_edge_degree(p);
		return stats;
	}
	template <typename ScoreEvaluator>
	Status run_topology_fix_pipeline(TopologyParameters& p, ScoreEvaluator& evaluator)
	{
		if (!p.skeleton_ || !p.incident_tets_)
			return Status::invalid_data;
		const DenseFiveVertexDetectionResult dense_five_vertex_detection = detect_k5_and_k5_minus_1_cells(p);
		const std::vector<IndependentDenseFiveVertexRepairSeed> independent_dense_repair_seeds =
			collect_independent_dense_five_vertex_repair_seeds(p, dense_five_vertex_detection);
		TopologyFixScoreCache score_cache;
		if (!initialize_topology_fix_score_cache(p, evaluator, score_cache))
			return Status::score_evaluation_failed;
		while (!p.skeleton_tets_.empty())
		{
			const uint32 round_start_tets = static_cast<uint32>(p.skeleton_tets_.size());
			const BoundaryTetPrepassStats boundary_stats = run_boundary_tet_face_deletion(p, &score_cache);
			EdgeTetModeRunStats simple_stats;
			if (!p.skeleton_tets_.empty())
				simple_stats = run_edge_score_tet_mode_topology_fix(
					p, EdgeTetDeleteMode::SimpleTet, score_cache, evaluator);
			EdgeTetModeRunStats nonsimple_stats;
			if (!p.skeleton_tets_.empty())
				nonsimple_stats = run_edge_score_tet_mode_topology_fix(
					p, EdgeTetDeleteMode::NonSimpleTet, score_cache, evaluator);
			const uint32 round_end_tets = static_cast<uint32>(p.skeleton_tets_.size());
			const uint32 round_removed_tets = round_start_tets - round_end_tets;
			++metrics_.topology_rounds;
			metrics_.boundary_removed_tets += boundary_stats.removed_tets;
			metrics_.simple_removed_tets += simple_stats.removed_tets;
			metrics_.nonsimple_removed_tets += nonsimple_stats.removed_tets;
			metrics_.removed_faces += boundary_stats.removed_faces + simple_stats.removed_faces + nonsimple_stats.removed_faces;
			metrics_.removed_edges += boundary_stats.removed_edges + simple_stats.removed_edges + nonsimple_stats.removed_edges;
			metrics_.removed_tets += round_removed_tets;
			if (round_removed_tets == 0)
				break;
		}
		const IndependentDenseFiveVertexRepairStats repair_stats = repair_isolated_independent_dense_five_vertex_regions(p, independent_dense_repair_seeds);
		metrics_.repaired_dense_regions += repair_stats.repaired_regions;
		metrics_.added_faces += repair_stats.added_faces;
		metrics_.tet_count = static_cast<uint32>(p.skeleton_tets_.size());
		return Status::success;
	}

	Data data_;
	TetMap skeleton_tets_;
	Metrics metrics_;
};

} // namespace geometry
} // namespace cgogn

#endif // CGOGN_GEOMETRY_ALGOS_UDF_SKELETON_TOPOLOGY_H_
