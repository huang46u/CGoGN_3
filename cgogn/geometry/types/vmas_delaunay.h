#ifndef CGOGN_DELAUNAY_MANAGER_H_
#define CGOGN_DELAUNAY_MANAGER_H_

#include <CGAL/Delaunay_triangulation_3.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Triangulation_vertex_base_with_info_3.h>
#include <unordered_map>
#include <vector>
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

	std::unordered_map<uint32, std::set<uint32>> compute_adjacent_clusters() const
	{
		std::unordered_map<uint32, std::set<uint32>> adjacent_clusters;
		for (auto cit = delaunay_.finite_cells_begin(); cit != delaunay_.finite_cells_end(); ++cit)
		{
			std::set<uint32> cell_clusters;
			for (int i = 0; i < 4; ++i)
			{
				auto vh = cit->vertex(i);
				if (vh->info().type == Point_Type::VOLUME_SAMPLE)
				{
					uint32 cid = vh->info().cluster_id;
					if (cid != (std::numeric_limits<uint32>::max)())
					{
						cell_clusters.insert(cid);
					}
				}
			}
			const std::size_t n = cell_clusters.size();
			if (n < 2 ||n > 3)
				continue;
			for (auto it1 = cell_clusters.begin(); it1 != cell_clusters.end(); ++it1)
			{
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
			
		}
		return adjacent_clusters;
	}

	

private:
	Point_3 to_cgal_point(const Vec3& position)
	{
		return Point_3(position[0], position[1], position[2]);
	}

private:
	Delaunay delaunay_;
	std::unordered_map<uint32, Vertex_handle> volume_sample_map_;
};

#endif // CGOGN_DELAUNAY_MANAGER_H_
} // namespace geometry
} // namespace cgogn