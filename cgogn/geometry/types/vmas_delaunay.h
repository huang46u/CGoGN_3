#ifndef CGOGN_DELAUNAY_MANAGER_H_
#define CGOGN_DELAUNAY_MANAGER_H_

#include <CGAL/Delaunay_triangulation_3.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Triangulation_vertex_base_with_info_3.h>
#include <queue>
#include <unordered_map>
#include <vector>
#include <map>
#include <unordered_set>
namespace cgogn
{
namespace geometry
{
using Scalar = cgogn::geometry::Scalar;
using Vec3 = cgogn::geometry::Vec3;
using K = CGAL::Exact_predicates_inexact_constructions_kernel;
using Point_3 = K::Point_3;

enum class Point_Type
{
	VOLUME_SAMPLE,
	SURFACE_VERTEX,
	SURFACE_PROJECTION
};
struct VertexInfo
{

	uint32 cgogn_index;
	Point_Type type;
	uint32 cluster_id;

	VertexInfo() : cgogn_index(0), type(Point_Type::VOLUME_SAMPLE), cluster_id((std::numeric_limits<uint32>::max)())
	{
	}
	VertexInfo(uint32 idx, Point_Type t, uint32 cluster_id) : cgogn_index(idx), type(t), cluster_id(cluster_id)
	{
	}
	VertexInfo(uint32 idx, Point_Type t) : VertexInfo(idx, t, std::numeric_limits<uint32>::max())
	{
	}
};

using Vb = CGAL::Triangulation_vertex_base_with_info_3<VertexInfo, K>;
using Cb = CGAL::Triangulation_cell_base_3<K>;
using Tds = CGAL::Triangulation_data_structure_3<Vb, Cb>;
using Delaunay = CGAL::Delaunay_triangulation_3<K, Tds>;

using Vertex_handle = Delaunay::Vertex_handle;
using Cell_handle = Delaunay::Cell_handle;
using Cell_circulatror = Delaunay::Cell_circulator;
class Vmas_Delaunay_Manager
{

public:
	Vertex_handle insert_sample(const Vec3& position, uint32 original_index, Point_Type type)
	{
		if (type == Point_Type::VOLUME_SAMPLE)
		{
			if (auto it = volume_sample_map_.find(original_index); it != volume_sample_map_.end())
				return it->second;
		}
		Vertex_handle vh = delaunay_.insert(to_cgal_point(position));
		vh->info() = VertexInfo(original_index, type);
		if (type == Point_Type::VOLUME_SAMPLE)
			volume_sample_map_[original_index] = vh;

		return vh;
	}
	void clear()
	{
		delaunay_.clear();
		volume_sample_map_.clear();
	}

	size_t number_of_vertices() const
	{
		return delaunay_.number_of_vertices();
	}
	Delaunay& get_delaunay()
	{
		return delaunay_;
	}
	const Delaunay& get_delaunay() const
	{
		return delaunay_;
	}
	double compute_cell_volume(const Cell_handle& cell) const
	{
		if (delaunay_.is_infinite(cell))
		{
			return 0.0;
		}

		Point_3 p0 = cell->vertex(0)->point();
		Point_3 p1 = cell->vertex(1)->point();
		Point_3 p2 = cell->vertex(2)->point();
		Point_3 p3 = cell->vertex(3)->point();

		CGAL::Tetrahedron_3<K> tet(p0, p1, p2, p3);
		return CGAL::to_double(tet.volume());
	}
	std::unordered_map<uint32, double> compute_samples_volume() const
	{
		std::unordered_map<uint32, double> samples_volumes;
		for (auto cit = delaunay_.finite_cells_begin(); cit != delaunay_.finite_cells_end(); ++cit)
		{
			double cell_volume = compute_cell_volume(cit);
			int volume_count = 0;
			for (int i = 0; i < 4; ++i)
				if (cit->vertex(i)->info().type == Point_Type::VOLUME_SAMPLE)
					++volume_count;

			if (volume_count == 0)
				continue;

			double share = cell_volume / volume_count;
			for (int i = 0; i < 4; ++i)
			{
				Vertex_handle vh = cit->vertex(i);
				if (vh->info().type != Point_Type::VOLUME_SAMPLE)
					continue;
				samples_volumes[vh->info().cgogn_index] += share;
			}
		}
		return samples_volumes;
	}


	bool set_cluster_id(uint32 cgogn_index, uint32 cluster_id)
	{
		auto it = volume_sample_map_.find(cgogn_index);
		if (it != volume_sample_map_.end())
		{
			Vertex_handle vh = it->second;
			vh->info().cluster_id = cluster_id;
			return true;
		}
		return false;
	}

	uint32 get_sample_cluster_id(uint32 cgogn_index) const
	{
		auto it = volume_sample_map_.find(cgogn_index);
		if (it == volume_sample_map_.end())
			return (std::numeric_limits<uint32>::max)();
		return it->second->info().cluster_id;
	}

	std::unordered_map<uint32, std::set<uint32>> compute_adjacent_clusters() const
	{
		std::unordered_map<uint32, std::set<uint32>> adjacent_clusters;

		// Build volume-sample list and cluster labels
		std::vector<uint32> local_to_index;
		local_to_index.reserve(volume_sample_map_.size());
		std::unordered_map<uint32, uint32> index_to_local;
		for (const auto& [idx, vh] : volume_sample_map_)
		{
			uint32 cid = vh->info().cluster_id;
			if (cid == (std::numeric_limits<uint32>::max)())
				continue;
			uint32 lid = uint32(local_to_index.size());
			index_to_local[idx] = lid;
			local_to_index.push_back(idx);
		}
		const uint32 n = uint32(local_to_index.size());
		if (n == 0)
			return adjacent_clusters;

		std::vector<uint32> cluster_of(n);
		for (uint32 lid = 0; lid < n; ++lid)
		{
			uint32 idx = local_to_index[lid];
			cluster_of[lid] = volume_sample_map_.at(idx)->info().cluster_id;
		}

		// adjacency list on local ids (only volume samples with valid cluster)
		std::vector<std::vector<uint32>> adj(n);
		for (auto eit = delaunay_.finite_edges_begin(); eit != delaunay_.finite_edges_end(); ++eit)
		{
			Vertex_handle v0 = eit->first->vertex(eit->second);
			Vertex_handle v1 = eit->first->vertex(eit->third);
			if (v0->info().type != Point_Type::VOLUME_SAMPLE || v1->info().type != Point_Type::VOLUME_SAMPLE)
				continue;
			uint32 c0 = v0->info().cluster_id;
			uint32 c1 = v1->info().cluster_id;
			if (c0 == (std::numeric_limits<uint32>::max)() || c1 == (std::numeric_limits<uint32>::max)())
				continue;
			auto it0 = index_to_local.find(v0->info().cgogn_index);
			auto it1 = index_to_local.find(v1->info().cgogn_index);
			if (it0 == index_to_local.end() || it1 == index_to_local.end())
				continue;
			uint32 l0 = it0->second;
			uint32 l1 = it1->second;
			adj[l0].push_back(l1);
			adj[l1].push_back(l0);
		}

		// find largest CC per cluster
		std::vector<uint8> visited(n, 0);
		std::vector<uint8> in_main(n, 0);
		std::unordered_map<uint32, uint32> cluster_best_size;
		for (uint32 lid = 0; lid < n; ++lid)
		{
			if (visited[lid])
				continue;
			uint32 cid = cluster_of[lid];
			uint32 comp_size = 0;
			std::queue<uint32> q;
			q.push(lid);
			visited[lid] = 1;
			std::vector<uint32> nodes;
			while (!q.empty())
			{
				uint32 u = q.front();
				q.pop();
				++comp_size;
				nodes.push_back(u);
				for (uint32 nb : adj[u])
				{
					if (!visited[nb] && cluster_of[nb] == cid)
					{
						visited[nb] = 1;
						q.push(nb);
					}
				}
			}
			if (comp_size > cluster_best_size[cid])
			{
				cluster_best_size[cid] = comp_size;
				for (uint32 k = 0; k < n; ++k)
					if (cluster_of[k] == cid)
						in_main[k] = 0;
				for (uint32 u : nodes)
					in_main[u] = 1;
			}
		}

		// collect adjacent clusters using only vertices in largest CC of each cluster
		for (auto cit = delaunay_.finite_cells_begin(); cit != delaunay_.finite_cells_end(); ++cit)
		{
			std::set<uint32> cell_clusters;
			for (int i = 0; i < 4; ++i)
			{
				auto vh = cit->vertex(i);
				if (vh->info().type != Point_Type::VOLUME_SAMPLE)
					continue;
				uint32 cid = vh->info().cluster_id;
				if (cid == (std::numeric_limits<uint32>::max)())
					continue;
				auto itl = index_to_local.find(vh->info().cgogn_index);
				if (itl == index_to_local.end())
					continue;
				uint32 lid = itl->second;
				if (!in_main[lid])
					continue;
				cell_clusters.insert(cid);
			}
			if (cell_clusters.size() < 2)
				continue;
			for (auto it1 = cell_clusters.begin(); it1 != cell_clusters.end(); ++it1)
				for (auto it2 = std::next(it1); it2 != cell_clusters.end(); ++it2)
				{
					uint32 c1 = *it1;
					uint32 c2 = *it2;
					if (c1 > c2)
						std::swap(c1, c2);
					adjacent_clusters[c1].insert(c2);
					adjacent_clusters[c2].insert(c1);
				}
		}
		return adjacent_clusters;
	}

	// adjacency using only the largest connected component of each cluster (ignores tiny CCs)
	std::unordered_map<uint32, std::set<uint32>> compute_adjacent_clusters_main_cc() const
	{
		std::unordered_map<uint32, std::set<uint32>> adjacent_clusters;

		// Build volume-sample list and cluster labels
		std::vector<uint32> local_to_index;
		local_to_index.reserve(volume_sample_map_.size());
		std::unordered_map<uint32, uint32> index_to_local;
		for (const auto& [idx, vh] : volume_sample_map_)
		{
			uint32 cid = vh->info().cluster_id;
			if (cid == (std::numeric_limits<uint32>::max)())
				continue;
			uint32 lid = uint32(local_to_index.size());
			index_to_local[idx] = lid;
			local_to_index.push_back(idx);
		}
		const uint32 n = uint32(local_to_index.size());
		if (n == 0)
			return adjacent_clusters;

		std::vector<uint32> cluster_of(n);
		for (uint32 lid = 0; lid < n; ++lid)
		{
			uint32 idx = local_to_index[lid];
			cluster_of[lid] = volume_sample_map_.at(idx)->info().cluster_id;
		}

		// adjacency list on local ids (only volume samples with valid cluster)
		std::vector<std::vector<uint32>> adj(n);
		for (auto eit = delaunay_.finite_edges_begin(); eit != delaunay_.finite_edges_end(); ++eit)
		{
			Vertex_handle v0 = eit->first->vertex(eit->second);
			Vertex_handle v1 = eit->first->vertex(eit->third);
			if (v0->info().type != Point_Type::VOLUME_SAMPLE || v1->info().type != Point_Type::VOLUME_SAMPLE)
				continue;
			uint32 c0 = v0->info().cluster_id;
			uint32 c1 = v1->info().cluster_id;
			if (c0 == (std::numeric_limits<uint32>::max)() || c1 == (std::numeric_limits<uint32>::max)())
				continue;
			auto it0 = index_to_local.find(v0->info().cgogn_index);
			auto it1 = index_to_local.find(v1->info().cgogn_index);
			if (it0 == index_to_local.end() || it1 == index_to_local.end())
				continue;
			uint32 l0 = it0->second;
			uint32 l1 = it1->second;
			adj[l0].push_back(l1);
			adj[l1].push_back(l0);
		}

		// find largest CC per cluster
		std::vector<uint8> visited(n, 0);
		std::vector<uint8> in_main(n, 0);
		std::unordered_map<uint32, uint32> cluster_best_size;
		for (uint32 lid = 0; lid < n; ++lid)
		{
			if (visited[lid])
				continue;
			uint32 cid = cluster_of[lid];
			uint32 comp_size = 0;
			std::queue<uint32> q;
			q.push(lid);
			visited[lid] = 1;
			std::vector<uint32> nodes;
			while (!q.empty())
			{
				uint32 u = q.front();
				q.pop();
				++comp_size;
				nodes.push_back(u);
				for (uint32 nb : adj[u])
				{
					if (!visited[nb] && cluster_of[nb] == cid)
					{
						visited[nb] = 1;
						q.push(nb);
					}
				}
			}
			if (comp_size > cluster_best_size[cid])
			{
				cluster_best_size[cid] = comp_size;
				for (uint32 u : nodes)
					in_main[u] = 1;
			}
		}

		// collect adjacent clusters using only main CC vertices via edges
		for (uint32 l0 = 0; l0 < n; ++l0)
		{
			if (!in_main[l0])
				continue;
			for (uint32 l1 : adj[l0])
			{
				if (!in_main[l1])
					continue;
				uint32 c0 = cluster_of[l0];
				uint32 c1 = cluster_of[l1];
				if (c0 == c1)
					continue;
				uint32 a = std::min(c0, c1);
				uint32 b = std::max(c0, c1);
				adjacent_clusters[a].insert(b);
				adjacent_clusters[b].insert(a);
			}
		}

		return adjacent_clusters;
	}

	// number of connected components within each cluster (using Delaunay edge adjacency between volume samples)
	// optionally returns sizes of each connected component per cluster in comps_sizes (vector of counts)
	std::unordered_map<uint32, uint32>
	compute_clusters_components(std::unordered_map<uint32, std::vector<uint32>>* comps_sizes = nullptr) const
	{
		std::unordered_map<uint32, uint32> components_count;
		if (volume_sample_map_.empty())
			return components_count;

		// map cgogn index -> local id
		std::vector<uint32> local_to_index;
		local_to_index.reserve(volume_sample_map_.size());
		std::unordered_map<uint32, uint32> index_to_local;
		for (const auto& [idx, vh] : volume_sample_map_)
		{
			uint32 cid = vh->info().cluster_id;
			if (cid == (std::numeric_limits<uint32>::max)())
				continue;
			uint32 lid = uint32(local_to_index.size());
			index_to_local[idx] = lid;
			local_to_index.push_back(idx);
		}
		const uint32 n = uint32(local_to_index.size());
		if (n == 0)
			return components_count;

		std::vector<uint32> cluster_of(n);
		for (uint32 lid = 0; lid < n; ++lid)
		{
			uint32 idx = local_to_index[lid];
			cluster_of[lid] = volume_sample_map_.at(idx)->info().cluster_id;
		}

		// adjacency list on local ids (only volume samples with valid cluster)
		std::vector<std::vector<uint32>> adj(n);
		for (auto eit = delaunay_.finite_edges_begin(); eit != delaunay_.finite_edges_end(); ++eit)
		{
			Vertex_handle v0 = eit->first->vertex(eit->second);
			Vertex_handle v1 = eit->first->vertex(eit->third);
			if (v0->info().type != Point_Type::VOLUME_SAMPLE || v1->info().type != Point_Type::VOLUME_SAMPLE)
				continue;
			uint32 c0 = v0->info().cluster_id;
			uint32 c1 = v1->info().cluster_id;
			if (c0 == (std::numeric_limits<uint32>::max)() || c1 == (std::numeric_limits<uint32>::max)())
				continue;
			auto it0 = index_to_local.find(v0->info().cgogn_index);
			auto it1 = index_to_local.find(v1->info().cgogn_index);
			if (it0 == index_to_local.end() || it1 == index_to_local.end())
				continue;
			uint32 l0 = it0->second;
			uint32 l1 = it1->second;
			adj[l0].push_back(l1);
			adj[l1].push_back(l0);
		}

		// group nodes by cluster
		std::unordered_map<uint32, std::vector<uint32>> cluster_nodes;
		for (uint32 lid = 0; lid < n; ++lid)
			cluster_nodes[cluster_of[lid]].push_back(lid);

		std::vector<uint8> visited(n, 0);
		for (const auto& [cid, nodes] : cluster_nodes)
		{
			uint32 comps = 0;
			if (comps_sizes)
				(*comps_sizes)[cid].clear();
			for (uint32 lid : nodes)
			{
				if (visited[lid])
					continue;
				++comps;
				// BFS restricted to same cluster
				std::queue<uint32> q;
				q.push(lid);
				visited[lid] = 1;
				uint32 comp_size = 0;
				while (!q.empty())
				{
					uint32 u = q.front();
					q.pop();
					++comp_size;
					for (uint32 nb : adj[u])
					{
						if (!visited[nb] && cluster_of[nb] == cid)
						{
							visited[nb] = 1;
							q.push(nb);
						}
					}
				}
				if (comps_sizes)
					(*comps_sizes)[cid].push_back(comp_size);
			}
			components_count[cid] = comps;
		}

		return components_count;
	}

	// number of connected components when considering only vertices from two adjacent clusters together
	using ClusterPair = std::pair<uint32, uint32>;
	struct PairHash
	{
		std::size_t operator()(const ClusterPair& p) const noexcept
		{
			return (std::size_t(p.first) << 32) ^ std::size_t(p.second);
		}
	};

	std::unordered_map<ClusterPair, uint32, PairHash> compute_adjacent_clusters_components() const
	{
		std::unordered_map<ClusterPair, uint32, PairHash> result;
		if (volume_sample_map_.empty())
			return result;

		// reuse cluster grouping and adjacency
		std::vector<uint32> local_to_index;
		local_to_index.reserve(volume_sample_map_.size());
		std::unordered_map<uint32, uint32> index_to_local;
		for (const auto& [idx, vh] : volume_sample_map_)
		{
			uint32 cid = vh->info().cluster_id;
			if (cid == (std::numeric_limits<uint32>::max)())
				continue;
			uint32 lid = uint32(local_to_index.size());
			index_to_local[idx] = lid;
			local_to_index.push_back(idx);
		}
		const uint32 n = uint32(local_to_index.size());
		if (n == 0)
			return result;

		std::vector<uint32> cluster_of(n);
		for (uint32 lid = 0; lid < n; ++lid)
		{
			uint32 idx = local_to_index[lid];
			cluster_of[lid] = volume_sample_map_.at(idx)->info().cluster_id;
		}

		std::vector<std::vector<uint32>> adj(n);
		for (auto eit = delaunay_.finite_edges_begin(); eit != delaunay_.finite_edges_end(); ++eit)
		{
			Vertex_handle v0 = eit->first->vertex(eit->second);
			Vertex_handle v1 = eit->first->vertex(eit->third);
			if (v0->info().type != Point_Type::VOLUME_SAMPLE || v1->info().type != Point_Type::VOLUME_SAMPLE)
				continue;
			uint32 c0 = v0->info().cluster_id;
			uint32 c1 = v1->info().cluster_id;
			if (c0 == (std::numeric_limits<uint32>::max)() || c1 == (std::numeric_limits<uint32>::max)())
				continue;
			auto it0 = index_to_local.find(v0->info().cgogn_index);
			auto it1 = index_to_local.find(v1->info().cgogn_index);
			if (it0 == index_to_local.end() || it1 == index_to_local.end())
				continue;
			uint32 l0 = it0->second;
			uint32 l1 = it1->second;
			adj[l0].push_back(l1);
			adj[l1].push_back(l0);
		}

		// collect adjacency pairs
		auto adj_clusters = compute_adjacent_clusters();
		std::vector<ClusterPair> pairs;
		for (const auto& [c1, set] : adj_clusters)
			for (uint32 c2 : set)
			{
				ClusterPair p = c1 < c2 ? ClusterPair{c1, c2} : ClusterPair{c2, c1};
				pairs.push_back(p);
			}
		std::sort(pairs.begin(), pairs.end(), [](const ClusterPair& a, const ClusterPair& b) {
			if (a.first != b.first)
				return a.first < b.first;
			return a.second < b.second;
		});
		pairs.erase(std::unique(pairs.begin(), pairs.end()), pairs.end());

		std::vector<uint8> visited(n, 0);
		for (const auto& p : pairs)
		{
			uint32 c1 = p.first;
			uint32 c2 = p.second;
			// reset visited
			std::fill(visited.begin(), visited.end(), 0);
			uint32 comps = 0;
			for (uint32 lid = 0; lid < n; ++lid)
			{
				if (visited[lid])
					continue;
				uint32 cid = cluster_of[lid];
				if (cid != c1 && cid != c2)
					continue;
				++comps;
				std::queue<uint32> q;
				q.push(lid);
				visited[lid] = 1;
				while (!q.empty())
				{
					uint32 u = q.front();
					q.pop();
					for (uint32 nb : adj[u])
					{
						if (!visited[nb] && (cluster_of[nb] == c1 || cluster_of[nb] == c2))
						{
							visited[nb] = 1;
							q.push(nb);
						}
					}
				}
			}
			result[p] = comps;
		}

		return result;
	}

	

private:
	Point_3 to_cgal_point(const Vec3& position)
	{
		return Point_3(position[0], position[1], position[2]);
	}

private:
	Delaunay delaunay_;
	std::unordered_map<uint32, Vertex_handle> volume_sample_map_;

public:
	using ClusterAdjacency = std::unordered_map<uint32, std::set<uint32>>;
	using ClusterPair = std::pair<uint32, uint32>;
	using ClusterPairComponents = std::unordered_map<ClusterPair, uint32, PairHash>;
	using ClusterComponentsSizes = std::unordered_map<uint32, std::vector<uint32>>;
};

#endif // CGOGN_DELAUNAY_MANAGER_H_
} // namespace geometry
} // namespace cgogn
