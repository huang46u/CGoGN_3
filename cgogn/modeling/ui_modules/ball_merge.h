/*******************************************************************************
 * CGoGN                                                                        *
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

#ifndef CGOGN_MODULE_SKELETON_EXTRACTOR_H_
#define CGOGN_MODULE_SKELETON_EXTRACTOR_H_

#include <cgogn/core/ui_modules/mesh_provider.h>
#include <cgogn/ui/app.h>
#include <cgogn/ui/module.h>

#include <cgogn/core/functions/convert.h>

#include <cgogn/core/functions/traversals/face.h>
#include <cgogn/core/functions/traversals/vertex.h>

#include <cgogn/geometry/types/vector_traits.h>

#include <CGAL/Delaunay_triangulation_3.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/IO/read_off_points.h>
#include <CGAL/Object.h>
#include <CGAL/Polygon_mesh_processing/distance.h>
#include <CGAL/Side_of_triangle_mesh.h>
#include <CGAL/Surface_mesh.h>
#include <libacc/bvh_tree.h>
#include <libacc/kd_tree.h>
#include <cgogn/geometry/algos/medial_axis.h>
#include <CGAL/Triangulation_cell_base_with_info_3.h>
#include <CGAL/Triangulation_vertex_base_with_info_3.h>
#include <CGAL/double.h>
#include <CGAL/optimal_bounding_box.h>

#include <boost/synapse/connect.hpp>

#include <unordered_map>

namespace cgogn
{

namespace ui
{

using geometry::Mat3;
using geometry::Scalar;
using geometry::Vec3;
using geometry::Vec4;

template <typename SURFACE, typename NONMANIFOLD> 
class BallMerge : public Module
{
	using K = CGAL::Exact_predicates_inexact_constructions_kernel;
	using Point = K::Point_3;
	using Weight_Point = K::Weighted_point_3;
	struct DelaunayCellInfo
	{
		int id = -1;
		bool inside = false;
		Point centroid;
		double radius2; // radius square
		int32 group_id = -1;
	};

	struct VertexInfo
	{
		int32 id = -1;
		bool is_medial = false;
	};
	using Vb = CGAL::Triangulation_vertex_base_with_info_3<VertexInfo, K>;
	// Delaunay
	using Cb = CGAL::Triangulation_cell_base_with_info_3<DelaunayCellInfo, K>;
	using Tds = CGAL::Triangulation_data_structure_3<Vb, Cb>;
	using Delaunay = CGAL::Delaunay_triangulation_3<K, Tds, CGAL::Fast_location>;
	using Delaunay_Cell_handle = typename Delaunay::Cell_handle;
	using Delaunay_Cell_circulator = typename Delaunay::Cell_circulator;
	using Delaunay_Vertex_handle = typename Delaunay::Vertex_handle;
	using Facet = typename Delaunay::Facet;

	using CGAL_Surface_mesh = CGAL::Surface_mesh<Point>;
	using Point_inside = CGAL::Side_of_triangle_mesh<CGAL_Surface_mesh, K>;
	using Primitive = CGAL::AABB_face_graph_triangle_primitive<CGAL_Surface_mesh>;
	using Traits = CGAL::AABB_traits_3<K, Primitive>;
	using Tree = CGAL::AABB_tree<Traits>;

	template <typename T>
	using SAttribute = typename mesh_traits<SURFACE>::template Attribute<T>;
	template <typename T>
	using NAttribute = typename mesh_traits<NONMANIFOLD>::template Attribute<T>;

	using SVertex = typename mesh_traits<SURFACE>::Vertex;
	using SEdge = typename mesh_traits<SURFACE>::Edge;
	using SFace = typename mesh_traits<SURFACE>::Face;

	using NVertex = typename mesh_traits<NONMANIFOLD>::Vertex;

	struct SurfaceParameters
	{
		SURFACE* mesh_;
		std::shared_ptr<SAttribute<Vec3>> surface_positions_ = nullptr;
		std::shared_ptr<SAttribute<Vec3>> surface_normals_ = nullptr;
		std::shared_ptr<SAttribute<Vec3>> medial_axis_position_ = nullptr;
		std::shared_ptr<SAttribute<Scalar>> medial_axis_radius_ = nullptr;
		std::shared_ptr<SAttribute<SVertex>> medial_axis_secondary_vertex_ = nullptr;

		acc::BVHTree<uint32, Vec3>* surface_bvh_ = nullptr;
		std::vector<SFace> surface_bvh_faces_;
		acc::KDTree<3, uint32>* surface_kdt_ = nullptr;
		std::vector<SVertex> surface_kdt_vertices_;

		NONMANIFOLD* result_mesh_ = nullptr;
		std::shared_ptr<NAttribute<Vec3>> result_mesh_positions_ = nullptr;

		NONMANIFOLD* result_mesh_2 = nullptr;
		std::shared_ptr<NAttribute<Vec3>> result_mesh_positions_2 = nullptr;

		CGAL_Surface_mesh csm;
		Delaunay delaunay;
		Tree tree;
		std::unique_ptr<Point_inside> inside_tester;
		uint32 cell_count = 0;
		float ball_merge_threshold_ = 2.0f;
	};

public:
	BallMerge(const App& app) : Module(app, "BallMerge"), selected_surface_(nullptr)
	{
	}

	~BallMerge()
	{
	}

	void init_surface_mesh(SURFACE* s)
	{
		SurfaceParameters& p = surface_parameters_[s];
		p.mesh_ = s;
		p.surface_positions_ = get_or_add_attribute<Vec3, SVertex>(*s, "position");
		p.surface_normals_ = get_or_add_attribute<Vec3, SVertex>(*s, "normal");
		geometry::compute_normal<SVertex>(*s, p.surface_positions_.get(), p.surface_normals_.get());

		MeshData<SURFACE>& md = surface_provider_->mesh_data(*s);
		uint32 nb_vertices = md.template nb_cells<SVertex>();
		uint32 nb_faces = md.template nb_cells<SFace>();

		auto bvh_vertex_index = get_or_add_attribute<uint32, SVertex>(*s, "__bvh_vertex_index");

		p.surface_kdt_vertices_.clear();
		p.surface_kdt_vertices_.reserve(nb_vertices);
		std::vector<Vec3> vertex_position_vector;
		vertex_position_vector.reserve(nb_vertices);
		uint32 idx = 0;
		foreach_cell(*s, [&](SVertex v) -> bool {
			p.surface_kdt_vertices_.push_back(v);
			value<uint32>(*s, bvh_vertex_index, v) = idx++;
			vertex_position_vector.push_back(value<Vec3>(*s, p.surface_positions_, v));
			return true;
		});

		p.surface_bvh_faces_.clear();
		p.surface_bvh_faces_.reserve(nb_faces);
		std::vector<uint32> face_vertex_indices;
		face_vertex_indices.reserve(nb_faces * 3);
		foreach_cell(*s, [&](SFace f) -> bool {
			p.surface_bvh_faces_.push_back(f);
			foreach_incident_vertex(*s, f, [&](SVertex v) -> bool {
				face_vertex_indices.push_back(value<uint32>(*s, bvh_vertex_index, v));
				return true;
			});
			return true;
		});

		if (p.surface_bvh_)
			delete p.surface_bvh_;
		p.surface_bvh_ = new acc::BVHTree<uint32, Vec3>(face_vertex_indices, vertex_position_vector);

		if (p.surface_kdt_)
			delete p.surface_kdt_;
		p.surface_kdt_ = new acc::KDTree<3, uint32>(vertex_position_vector);

		remove_attribute<SVertex>(*s, bvh_vertex_index);
		p.medial_axis_position_ = get_or_add_attribute<Vec3, SVertex>(*s, "medial_axis_position");
		p.medial_axis_radius_ = get_or_add_attribute<Scalar, SVertex>(*s, "medial_axis_radius");
		p.medial_axis_secondary_vertex_ = get_or_add_attribute<SVertex, SVertex>(*s, "medial_axis_secondary_vertex_");

		parallel_foreach_cell(*s, [&](SVertex v) -> bool {
			uint32 v_index = index_of(*s, v);
			auto [c, r, q] = cgogn::geometry::shrinking_ball_center(
				*s, (*p.surface_positions_)[v_index], (*p.surface_normals_)[v_index], p.surface_positions_.get(),
				p.surface_bvh_, p.surface_bvh_faces_, p.surface_kdt_,
				p.surface_kdt_vertices_);
			(*p.medial_axis_position_)[v_index] = c;
			(*p.medial_axis_radius_)[v_index] = r;
			(*p.medial_axis_secondary_vertex_)[v_index] = q;
			return true;
		});

		if (!p.result_mesh_)
			p.result_mesh_ = nonmanifold_provider_->add_mesh(surface_provider_->mesh_name(*s) + "_ball_merge");
		p.result_mesh_positions_ = get_or_add_attribute<Vec3, NVertex>(*p.result_mesh_, "position");

		if (!p.result_mesh_2)
			p.result_mesh_2 = nonmanifold_provider_->add_mesh(surface_provider_->mesh_name(*s) + "_ball_merge_2");
		p.result_mesh_positions_2 = get_or_add_attribute<Vec3, NVertex>(*p.result_mesh_2, "position");

		std::string filename = surface_provider_->mesh_filename(*s);
		if (!filename.empty())
		{
			std::ifstream input(filename);
			if (!input || !(input >> p.csm))
			{
				std::cerr << "Error: input file could not be read" << std::endl;
				return;
			}
		}
		normalize_surface_mesh(p.csm);
		// Convert to CGAL Surface_mesh
		p.tree = Tree(faces(p.csm).first, faces(p.csm).second, p.csm);
		p.tree.accelerate_distance_queries();
		p.inside_tester = std::make_unique<Point_inside>(p.tree);
		compute_delaunay(p);
	}

	void global_ball_merge(SurfaceParameters& p)
	{

		for (auto cit = p.delaunay.finite_cells_begin(); cit != p.delaunay.finite_cells_end(); ++cit)
		{
			cit->info().group_id = -1;
		}

		int32 group_id = -1;
		uint32 max_group_id = 0;
		uint32 second_group_id = 0;
		uint32 max_first = 0;
		uint32 max_second = 0;

		for (auto cit = p.delaunay.finite_cells_begin(); cit != p.delaunay.finite_cells_end(); ++cit)
		{
			if (cit->info().group_id == -1)
			{
				group_id++;
				uint32 gcount = Group(cit, group_id, p.ball_merge_threshold_);
				if (max_first < gcount)
				{
					second_group_id = max_group_id;
					max_second = max_first;
					max_group_id = group_id;
					max_first = gcount;
				}
				else if (max_second < gcount && gcount != max_first)
				{
					second_group_id = group_id;
					max_second = gcount;
				}
			}
		}
		std::cout << "Largest group size: " << max_first << std::endl;
		std::cout << "Second largest group size: " << max_second << std::endl;
		
		construct_result_mesh(p, p.result_mesh_, max_group_id);
		construct_result_mesh(p, p.result_mesh_2, second_group_id);
	}

private:
	void normalize_surface_mesh(CGAL_Surface_mesh& mesh)
	{
		if (mesh.is_empty())
			return;

		CGAL::Bbox_3 bbox;
		for (auto v : mesh.vertices())
			bbox = bbox + mesh.point(v).bbox();

		double cx = (bbox.xmin() + bbox.xmax()) / 2.0;
		double cy = (bbox.ymin() + bbox.ymax()) / 2.0;
		double cz = (bbox.zmin() + bbox.zmax()) / 2.0;
		double max_dim = std::max({bbox.xmax() - bbox.xmin(), bbox.ymax() - bbox.ymin(), bbox.zmax() - bbox.zmin()});

		for (auto v : mesh.vertices())
		{
			Point p = mesh.point(v);
			double nx = (p.x() - bbox.xmin()) / max_dim;
			double ny = (p.y() - bbox.ymin()) / max_dim;
			double nz = (p.z() - bbox.zmin()) / max_dim;
			mesh.point(v) = Point(nx, ny, nz);
		}
	}

	void compute_delaunay(SurfaceParameters& p)
	{
		p.delaunay.clear();
		p.cell_count = 0;
		foreach_cell(*p.mesh_, [&](SVertex v) {
			auto index = index_of(*p.mesh_, v);
			auto pos = (*p.surface_positions_)[index];
			auto medial_pos = (*p.medial_axis_position_)[index];

			Delaunay_Vertex_handle vh = p.delaunay.insert(Point(pos[0], pos[1], pos[2]));
			vh->info().is_medial = false;

			Delaunay_Vertex_handle mh = p.delaunay.insert(Point(medial_pos[0], medial_pos[1], medial_pos[2]));
			mh->info().is_medial = true;
			
			return true;
		});
		
		

		for (auto cit = p.delaunay.finite_cells_begin(); cit != p.delaunay.finite_cells_end(); ++cit)
		{
			cit->info().centroid = CGAL::circumcenter(p.delaunay.tetrahedron(cit));
			cit->info().radius2 = CGAL::squared_distance(cit->info().centroid, cit->vertex(0)->point());
			cit->info().inside = ((*p.inside_tester)(cit->info().centroid) == CGAL::ON_BOUNDED_SIDE);
			cit->info().id = p.cell_count;
			p.cell_count++;
		}
		std::cout << "Cell count: " << p.cell_count << std::endl;
	}

	bool can_merge(Delaunay_Cell_handle c1, Delaunay_Cell_handle c2, float threshold)
	{
		double r0 = CGAL::sqrt(c1->info().radius2);
		double r1 = CGAL::sqrt(c2->info().radius2);
		double d = CGAL::sqrt(CGAL::squared_distance(c1->info().centroid, c2->info().centroid));
		if (std::isnan(d))
			return false;
		double ir = std::max((r0 + r1 - d) / r0, (r0 + r1 - d) / r1);
		return ir >= threshold;
	}

	void construct_result_mesh(SurfaceParameters& p, NONMANIFOLD* s, uint32 group_id)
	{
		clear(*s);
		cgogn::io::IncidenceGraphImportData ball_merge_non_manifold_data;
		std::unordered_map<std::pair<uint32, uint32>, uint32, edge_hash, edge_equal> edge_indices;
		auto poisition = get_or_add_attribute<Vec3, NVertex>(*s, "position");

		uint32 global_count = 0;
		for (auto vit = p.delaunay.finite_vertices_begin(); vit != p.delaunay.finite_vertices_end(); ++vit)
		{
			if (!vit->info().is_medial)
				continue;

			Point point = vit->point();
			vit->info().id = global_count++;
			ball_merge_non_manifold_data.vertex_position_.push_back(Vec3(point.x(), point.y(), point.z()));
		}
	
		for (auto cit = p.delaunay.finite_cells_begin(); cit != p.delaunay.finite_cells_end(); ++cit)
		{
			for (int i = 0; i < 4; i++)
			{
				if (cit->info().group_id == group_id &&
					(cit->neighbor(i)->info().group_id != group_id || p.delaunay.is_infinite(cit->neighbor(i))))
				{
					std::vector<int> indices(3);
					bool is_medial = true;
					for (size_t idx = 0; idx < 3; ++idx)
					{
						auto v_info = cit->vertex((i + 1 + idx) % 4)->info();
						is_medial &= v_info.is_medial;
						indices[idx] = v_info.id;
					}
					if (!is_medial)
						continue;
					ball_merge_non_manifold_data.faces_nb_edges_.push_back(3);
					for (size_t i = 0; i < 3; ++i)
					{
						if (edge_indices.find({indices[i], indices[(i + 1) % 3]}) == edge_indices.end())
						{
							edge_indices[{indices[i], indices[(i + 1) % 3]}] = uint32(edge_indices.size());
							ball_merge_non_manifold_data.edges_vertex_indices_.push_back(indices[i]);
							ball_merge_non_manifold_data.edges_vertex_indices_.push_back(indices[(i + 1) % 3]);
						}

						ball_merge_non_manifold_data.faces_edge_indices_.push_back(
							edge_indices[{indices[i], indices[(i + 1) % 3]}]);
					}
				}
					
				//}
			}
		}
		uint32 nb_vertices = ball_merge_non_manifold_data.vertex_position_.size();
		uint32 nb_edges = ball_merge_non_manifold_data.edges_vertex_indices_.size() / 2;
		uint32 nb_faces = ball_merge_non_manifold_data.faces_nb_edges_.size();
		ball_merge_non_manifold_data.reserve(nb_vertices,nb_edges,nb_faces);
		cgogn::io::import_incidence_graph_data(*s, ball_merge_non_manifold_data);
		nonmanifold_provider_->emit_connectivity_changed(*s);
		nonmanifold_provider_->emit_attribute_changed(*s, poisition.get());
	}
	uint32 Group(Delaunay_Cell_handle ch, uint32 group, float threshold)
	{
		std::queue<Delaunay_Cell_handle> cell_queue;
		cell_queue.push(ch);
		uint32 gcount = 0;
		while (!cell_queue.empty())
		{
			Delaunay_Cell_handle current = cell_queue.front();
			cell_queue.pop();
			current->info().group_id = group;
			gcount++;
			for (int i = 0; i < 4; i++)
			{
				Delaunay_Cell_handle neighbor = current->neighbor(i);
				if (neighbor->info().group_id == -1 && can_merge(current, neighbor, threshold))
				{
					neighbor->info().group_id = group;
					cell_queue.push(neighbor);
				}
			}
		}
		return gcount;
	}
	
	struct edge_hash
	{
		std::size_t operator()(const std::pair<uint32, uint32>& edge) const
		{
			return std::hash<uint32>()(edge.first) + std::hash<uint32>()(edge.second);
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

protected:
	void init() override
	{

		surface_provider_ = static_cast<ui::MeshProvider<SURFACE>*>(
			app_.module("MeshProvider (" + std::string{mesh_traits<SURFACE>::name} + ")"));
		nonmanifold_provider_ = static_cast<ui::MeshProvider<NONMANIFOLD>*>(
			app_.module("MeshProvider (" + std::string{mesh_traits<NONMANIFOLD>::name} + ")"));
		surface_provider_->foreach_mesh([this](SURFACE& m, const std::string&) { init_surface_mesh(&m); });
	}

	void left_panel() override
	{
		imgui_mesh_selector(surface_provider_, selected_surface_, "Surface", [&](SURFACE& m) {
			selected_surface_ = &m;
			surface_provider_->mesh_data(m).outlined_until_ = App::frame_time_ + 1.0;
		});

		if (selected_surface_)
		{
			SurfaceParameters& p = surface_parameters_[selected_surface_];

			ImGui::Separator();
			ImGui::SliderFloat("Ball Merge threshold", &p.ball_merge_threshold_, 0.0f, 2.0f);

			if (ImGui::IsItemDeactivatedAfterEdit())
			{
				global_ball_merge(p);
			}
		}
	}

private:
	MeshProvider<SURFACE>* surface_provider_ = nullptr;
	MeshProvider<NONMANIFOLD>* nonmanifold_provider_ = nullptr;

	SURFACE* selected_surface_ = nullptr;
	std::unordered_map<const SURFACE*, SurfaceParameters> surface_parameters_;

	std::vector<std::shared_ptr<boost::synapse::connection>> connections_;
};

} // namespace ui

} // namespace cgogn

#endif // CGOGN_MODULE_SKELETON_EXTRACTOR_H_
