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

template <typename SURFACE>
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
	using SurfaceAttribute = typename mesh_traits<SURFACE>::template Attribute<T>;

	using SurfaceVertex = typename mesh_traits<SURFACE>::Vertex;
	using SurfaceEdge = typename mesh_traits<SURFACE>::Edge;
	using SurfaceFace = typename mesh_traits<SURFACE>::Face;

	struct SurfaceParameters
	{
		SURFACE* mesh_;
		std::shared_ptr<SurfaceAttribute<Vec3>> surface_positions_ = nullptr;

		SURFACE* result_mesh_ = nullptr;
		std::shared_ptr<SurfaceAttribute<Vec3>> result_mesh_positions_ = nullptr;

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
		p.surface_positions_ = get_or_add_attribute<Vec3, SurfaceVertex>(*s, "position");
		if (!p.result_mesh_)
			p.result_mesh_ = surface_provider_->add_mesh(surface_provider_->mesh_name(*s) + "_ball_merge");
		p.result_mesh_positions_ = get_or_add_attribute<Vec3, SurfaceVertex>(*p.result_mesh_, "position");

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
		foreach_cell(*p.mesh_, [&](SurfaceVertex v) {
			auto index = index_of(*p.mesh_, v);
			auto pos = (*p.surface_positions_)[index];
			p.delaunay.insert(Point(pos[0], pos[1], pos[2]));
			return true;
		});

		uint32 inside_counter_ = 0;
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
			return 0;
		double ir = std::max((r0 + r1 - d) / r0, (r0 + r1 - d) / r1);
		return ir >= threshold;
	}

	void construct_result_mesh(SurfaceParameters& p, std::vector<Delaunay_Cell_handle> largest_group)
	{
		clear(*p.result_mesh_);
		for (auto vit = p.delaunay.finite_vertices_begin(); vit != p.delaunay.finite_vertices_end(); ++vit)
		{
			vit->info().id = -1;
		}
		uint32 global_count = 0;
		cgogn::io::SurfaceImportData ball_merge_surface_data;
		for (const auto& cell : largest_group)
		{
			for (int i = 0; i < 4; i++)
			{
				Facet f(cell, i);
				auto g_id = cell->info().group_id;
				auto neighbour = p.delaunay.mirror_facet(f).first;
				int neighbour_g_id = -1;
				if (!p.delaunay.is_infinite(neighbour))
					neighbour_g_id = neighbour->info().group_id;

				if (g_id != neighbour_g_id)
				{
					Delaunay_Vertex_handle v[3];
					std::vector<uint32> indices;
					uint32 count = 0;
					for (size_t idx = 0; idx < 4; ++idx)
					{
						if (idx != f.second)
						{
							v[count] = f.first->vertex(idx);
							if (v[count]->info().id == -1)
							{
								ball_merge_surface_data.nb_vertices_++;
								ball_merge_surface_data.vertex_position_.push_back(
									Vec3(v[count]->point().x(), v[count]->point().y(), v[count]->point().z()));
								v[count]->info().id = global_count;
								global_count++;
							}
							indices.push_back(v[count]->info().id);
							count++;
						}
					}
					ball_merge_surface_data.nb_faces_++;
					ball_merge_surface_data.faces_nb_vertices_.push_back(3);
					ball_merge_surface_data.faces_vertex_indices_.insert(
						ball_merge_surface_data.faces_vertex_indices_.end(), indices.begin(), indices.end());
				}
			}

			// for (int i = 0; i < 4; i++)
			// {
			// 	Facet f(cell, i);
			// 	Delaunay_Vertex_handle v[3];
			// 	std::vector<uint32> indices;
			// 	uint32 count = 0;
			// 	for (size_t idx = 0; idx < 4; ++idx)
			// 	{
			// 		if (idx != f.second)
			// 		{
			// 			v[count] = f.first->vertex(idx);
			// 			if (v[count]->info().id == -1)
			// 			{
			// 				ball_merge_surface_data.nb_vertices_++;
			// 				ball_merge_surface_data.vertex_position_.push_back(
			// 					Vec3(v[count]->point().x(), v[count]->point().y(), v[count]->point().z()));
			// 				v[count]->info().id = global_count;
			// 				global_count++;
			// 			}
			// 			indices.push_back(v[count]->info().id);
			// 			count++;
			// 		}
			// 	}
			// 	ball_merge_surface_data.nb_faces_++;
			// 	ball_merge_surface_data.faces_nb_vertices_.push_back(3);
			// 	ball_merge_surface_data.faces_vertex_indices_.insert(
			// 		ball_merge_surface_data.faces_vertex_indices_.end(), indices.begin(), indices.end());

			// }
		}
		cgogn::io::import_surface_data(*p.result_mesh_, ball_merge_surface_data);
		surface_provider_->emit_connectivity_changed(*p.result_mesh_);
		surface_provider_->emit_attribute_changed(*p.result_mesh_, p.result_mesh_positions_.get());
	}

	void global_ball_merge(SurfaceParameters& p)
	{

		for (auto cit = p.delaunay.finite_cells_begin(); cit != p.delaunay.finite_cells_end(); ++cit)
		{
			cit->info().group_id = -1;
		}

		int32 group_id = -1;
		std::queue<Delaunay_Cell_handle> cell_queue;
		std::vector<Delaunay_Cell_handle> largest_group;
		for (auto cit = p.delaunay.finite_cells_begin(); cit != p.delaunay.finite_cells_end(); ++cit)
		{

			if (cit->info().id != -1 && cit->info().group_id == -1)
			{
				std::vector<Delaunay_Cell_handle> new_group;
				group_id++;
				cell_queue.push(cit);
				while (!cell_queue.empty())
				{
					Delaunay_Cell_handle current = cell_queue.front();
					cell_queue.pop();
					new_group.push_back(current);
					current->info().group_id = group_id;
					for (int i = 0; i < 4; i++)
					{
						Delaunay_Cell_handle neighbor = current->neighbor(i);
						if (p.delaunay.is_infinite(neighbor))
							continue;
						if (neighbor->info().group_id == -1 && can_merge(current, neighbor, p.ball_merge_threshold_))
						{
							neighbor->info().group_id = group_id;
							cell_queue.push(neighbor);
						}
					}
					// std::cout << "Cell queue size: " << cell_queue.size() << std::endl;
				}
				if (new_group.size() > largest_group.size())
					largest_group = new_group;
			}
		}
		std::cout << "Largest group size: " << largest_group.size() << std::endl;
		construct_result_mesh(p, largest_group);
	}

protected:
	void init() override
	{

		surface_provider_ = static_cast<ui::MeshProvider<SURFACE>*>(
			app_.module("MeshProvider (" + std::string{mesh_traits<SURFACE>::name} + ")"));

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

	SURFACE* selected_surface_ = nullptr;
	std::unordered_map<const SURFACE*, SurfaceParameters> surface_parameters_;

	std::vector<std::shared_ptr<boost::synapse::connection>> connections_;
};

} // namespace ui

} // namespace cgogn

#endif // CGOGN_MODULE_SKELETON_EXTRACTOR_H_
