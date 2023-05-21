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

#ifndef CGOGN_MODULE_POWER_SHAPE_H_
#define CGOGN_MODULE_POWER_SHAPE_H_

#include <cgogn/core/ui_modules/mesh_provider.h>
#include <cgogn/ui/app.h>
#include <cgogn/ui/module.h>

#include <cgogn/geometry/types/vector_traits.h>

// import CGAL
#include <CGAL/Delaunay_triangulation_3.h>
#include <CGAL/Delaunay_triangulation_cell_base_3.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/IO/read_off_points.h>
#include <CGAL/Object.h>
#include <CGAL/Polygon_mesh_processing/distance.h>
#include <CGAL/Regular_triangulation_3.h>
#include <CGAL/Regular_triangulation_cell_base_3.h>
#include <CGAL/Regular_triangulation_vertex_base_3.h>
#include <CGAL/Side_of_triangle_mesh.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/Triangulation_vertex_base_with_info_3.h>
#include <CGAL/optimal_bounding_box.h>

namespace cgogn
{

namespace ui
{

template <typename SURFACE, typename NONMANIFOLD>
class PowerShape : public Module
{
	static_assert(mesh_traits<SURFACE>::dimension >= 2, "PowerShape can only be used with meshes of dimension >= 2");

	// Kernel for construct Delaunay
	using K = CGAL::Exact_predicates_inexact_constructions_kernel;
	using Vb = CGAL::Triangulation_vertex_base_with_info_3<uint32_t, K>;
	// Delaunay
	using Cb = CGAL::Delaunay_triangulation_cell_base_3<K>;
	using Tds = CGAL::Triangulation_data_structure_3<Vb, Cb>;
	using Delaunay = CGAL::Delaunay_triangulation_3<K, Tds, CGAL::Fast_location>;
	// Regular
	using Vb0 = CGAL::Regular_triangulation_vertex_base_3<K>;
	using RVb = CGAL::Triangulation_vertex_base_with_info_3<std::pair<uint32_t, bool>, K, Vb0>;
	using RCb = CGAL::Regular_triangulation_cell_base_3<K>;
	using RTds = CGAL::Triangulation_data_structure_3<RVb, RCb>;
	using Regular = CGAL::Regular_triangulation_3<K, RTds, CGAL::Fast_location>;

	using Point = K::Point_3;
	using Weight = K::FT;
	using Weight_Point = K::Weighted_point_3;

	using Cgal_Surface_mesh = CGAL::Surface_mesh<Point>;
	using Point_inside = CGAL::Side_of_triangle_mesh<Cgal_Surface_mesh, K>;
	using Primitive = CGAL::AABB_face_graph_triangle_primitive<Cgal_Surface_mesh>;
	using Traits = CGAL::AABB_traits<K, Primitive>;
	using Tree = CGAL::AABB_tree<Traits>;

	template <typename T>
	using NonManifoldAttribute = typename mesh_traits<NONMANIFOLD>::template Attribute<T>;
	template <typename T>
	using SurfaceAttribute = typename mesh_traits<SURFACE>::template Attribute<T>;

	using NonManifoldVertex = typename mesh_traits<NONMANIFOLD>::Vertex;
	using NonManifoldEdge = typename mesh_traits<NONMANIFOLD>::Edge;
	using NonManifoldFace = typename mesh_traits<NONMANIFOLD>::Face;

	using SurfaceVertex = typename mesh_traits<SURFACE>::Vertex;
	using SurfaceEdge = typename mesh_traits<SURFACE>::Edge;
	using SurfaceFace = typename mesh_traits<SURFACE>::Face;

	using Vec3 = geometry::Vec3;
	using Scalar = geometry::Scalar;

public:
	PowerShape(const App& app)
		: Module(app, "PowerShape"), selected_surface_mesh_(nullptr)
	{
	}
	~PowerShape()
	{
	}

private:
	
	struct point_hash
	{
		std::size_t operator()(const Point& p) const
		{
			return ((std::hash<double>()(p.x()) ^ (std::hash<double>()(p.y()) << 1)) >> 1) ^
				   (std::hash<double>()(p.z()) << 1);
		}
	};

	struct edge_hash
	{
		std::size_t operator()(const std::pair<uint32, uint32>& edge) const
		{
			return std::hash<uint32>()(edge.first) + std::hash<uint32>()(edge.second);
		}
	};

	bool pointInside(Tree& tree, Point& query)
	{
		// Initialize the point-in-polyhedron tester
		Point_inside inside_tester(tree);

		// Determine the side and return true if inside!
		return inside_tester(query) == CGAL::ON_BOUNDED_SIDE;
	}

public:
	Delaunay compute_delaunay_tredrahedron(Cgal_Surface_mesh& csm, Tree& tree)
	{
		std::vector<Point> Delaunay_tri_point;
		std::array<Point, 8> obb_points;
		Point acc(0, 0, 0);
		CGAL::oriented_bounding_box(csm, obb_points, CGAL::parameters::use_convex_hull(true));
		for (size_t i = 0; i < obb_points.size(); i++)
		{
			acc += K::Vector_3(obb_points[i].x(), obb_points[i].y(), obb_points[i].z());
		}
		std::array<double, 3> center{acc.x() / 8, acc.y() / 8, acc.z() / 8}; // TODO
		// Create a large box surrounding object so that the Voronoi vertices are bounded
		double offset = 2.0; // TODO
		std::array<double, 3> offset_array = {offset, offset, offset};
		std::array<std::array<double, 3>, 8> cube_corners = {
			{{center[0] - offset_array[0], center[1] - offset_array[1], center[2] - offset_array[2]},
			 {center[0] - offset_array[0], center[1] - offset_array[1], center[2] + offset_array[2]},
			 {center[0] - offset_array[0], center[1] + offset_array[1], center[2] - offset_array[2]},
			 {center[0] - offset_array[0], center[1] + offset_array[1], center[2] + offset_array[2]},
			 {center[0] + offset_array[0], center[1] - offset_array[1], center[2] - offset_array[2]},
			 {center[0] + offset_array[0], center[1] - offset_array[1], center[2] + offset_array[2]},
			 {center[0] + offset_array[0], center[1] + offset_array[1], center[2] - offset_array[2]},
			 {center[0] + offset_array[0], center[1] + offset_array[1], center[2] + offset_array[2]}}};

		// Sampling the mesh surface
		std::vector<Point> mesh_samples;
		CGAL::Polygon_mesh_processing::sample_triangle_mesh(
			csm, std::back_inserter(mesh_samples),
			// CGAL::parameters::use_monte_carlo_sampling(true).number_of_points_per_area_unit(0.01));
			CGAL::parameters::use_grid_sampling(true).grid_spacing(0.1));

		// 	Add bounding box vertices in the sample points set
		for (auto& p : cube_corners)
		{
			Delaunay_tri_point.emplace_back(p[0], p[1], p[2]);
		}

		// Add sampled vertices into the volume data to construct the delauney tredrahedron
		for (auto& s : mesh_samples)
		{
			Delaunay_tri_point.emplace_back(s[0], s[1], s[2]);
		}

		uint32 nb_vertices = Delaunay_tri_point.size();

		//		auto start_timer = std::chrono::high_resolution_clock::now();

		// Indices info for constructing volume data in Cgogn
		std::vector<unsigned> indices;
		indices.reserve(nb_vertices);
		for (unsigned int i = 0; i < nb_vertices; ++i)
			indices.push_back(i);

		// Construct delauney tredrahedron using CGAL
		Delaunay tri(boost::make_zip_iterator(boost::make_tuple(Delaunay_tri_point.begin(), indices.begin())),
					 boost::make_zip_iterator(boost::make_tuple(Delaunay_tri_point.end(), indices.end())));
		return tri;
	}
	
	void compute_inner_voronoi(Delaunay& tri, Tree& tree, std::vector<Weight_Point>& power_point,
							   std::vector<std::pair<uint32, bool>>& point_info,
							   std::unordered_map<uint32, uint32>& inside_indices)
	{
		uint32 count = 0, inside_vertices_count = 0;
		double dis;
		std::unordered_set<Point, point_hash> voronoi_vertices;
		std::vector<Delaunay::Cell_handle> cells;
		for (auto vit = tri.finite_vertices_begin(); vit != tri.finite_vertices_end(); ++vit)
		{
			cells.clear();
			tri.finite_incident_cells(vit, std::back_inserter(cells));
			if (cells.size())
			{
				for (auto c = cells.begin(); c != cells.end(); ++c)
				{
					Point centroid = tri.dual(*c);
					auto dis = CGAL::squared_distance(centroid, vit->point());
					if (voronoi_vertices.find(centroid) == voronoi_vertices.end())
					{
						voronoi_vertices.insert(centroid);
						power_point.push_back(Weight_Point(centroid, dis));
						if (pointInside(tree, centroid))
						{
							point_info.push_back({count, true});
							inside_indices.insert({count++, inside_vertices_count++});
						}
						else if (!pointInside(tree, centroid))
						{
							point_info.push_back({count++, false});
						
						}
					}
				}
			}
		}
	}

	void compute_inner_poles(Delaunay& tri, Tree& tree, std::vector<Weight_Point>& power_point,
							 std::vector<std::pair<uint32, bool>>& point_info,
							 std::unordered_map<uint32, uint32>& inside_indices)
	{
		std::vector<std::pair<uint32, bool>> power_indices;
		// Find inside and outside poles
		uint32 count = 0, inside_poles_count = 0;
		double dis;
		Point farthest_inside_point, farthest_outside_point;
		double farthest_inside_distance, farthest_outside_distance;
		std::unordered_set<Point, point_hash> poles;
		 // Use for the construction of medial axis
		std::vector<Delaunay::Cell_handle> cells;
		for (auto vit = tri.finite_vertices_begin(); vit != tri.finite_vertices_end(); ++vit)
		{
			farthest_inside_distance = 0;
			farthest_outside_distance = 0;
			cells.clear();
			tri.finite_incident_cells(vit, std::back_inserter(cells));
			if (cells.size())
			{
				for (auto c = cells.begin(); c != cells.end(); ++c)
				{
					Point centroid = tri.dual(*c);
					dis = CGAL::squared_distance(centroid, vit->point());
					if (pointInside(tree, centroid) && dis > farthest_inside_distance)
					{
						farthest_inside_point = centroid;
						farthest_inside_distance = dis;
					}
					else if (!pointInside(tree, centroid) && dis > farthest_outside_distance)
					{
						farthest_outside_point = centroid;
						farthest_outside_distance = dis;
					}
				}
			}
			if (farthest_inside_distance != 0)
			{
				if (poles.find(farthest_inside_point) == poles.end())
				{
					poles.insert(farthest_inside_point);
					power_point.push_back(Weight_Point(farthest_inside_point, farthest_inside_distance));
					point_info.push_back({count, true});
					inside_indices.insert({count++, inside_poles_count++});
				}
			}
			if (farthest_outside_distance != 0)
			{
				if (poles.find(farthest_outside_point) == poles.end())
				{
					poles.insert(farthest_outside_point);
					power_point.push_back(Weight_Point(farthest_outside_point, farthest_outside_distance));
					point_info.push_back({count++, false});
				}
			}
		}
	}
	
	void constrcut_power_diagram(NONMANIFOLD* mv, SURFACE& surface, std::vector<Weight_Point>& power_point,
								 std::vector<std::pair<uint32, bool>>& point_info,
								 std::unordered_map<uint32, uint32>& inside_indices)
	{
		
		cgogn::io::IncidenceGraphImportData Power_shape_data;
		std::unordered_map<std::pair<uint32, uint32>, uint32, edge_hash> edge_indices;
		uint32 edge_count = 0;

		medial_axis = Regular(boost::make_zip_iterator(boost::make_tuple(power_point.begin(), point_info.begin())),
					boost::make_zip_iterator(boost::make_tuple(power_point.end(), point_info.end())));
		// Add vertices
		for (size_t idx = 0; idx < power_point.size(); ++idx)
		{
			if (point_info[idx].second)
			{
				Power_shape_data.vertex_position_.push_back(
					{power_point[idx].x(), power_point[idx].y(), power_point[idx].z()});
				
			}	
		}
		bool inside;
		uint32 v, v1, v2;
		for (auto fit = medial_axis.finite_facets_begin(); fit != medial_axis.finite_facets_end(); ++fit)
		{
			inside = true;
			v = fit->second;
			// If face is inside
			for (size_t idx = 0; idx < 4; ++idx)
			{
				if (idx != v)
				{
					inside &= fit->first->vertex(idx)->info().second;
				}
			}
			if (inside)
			{
				Power_shape_data.faces_nb_edges_.push_back(3);

				for (size_t i = 0; i < 4; i++)
				{
					if (i != v)
					{
						for (size_t j = i + 1; j < 4; j++)
						{
							if (j != v)
							{

								// Add edge
								v1 = inside_indices[fit->first->vertex(i)->info().first];
								v2 = inside_indices[fit->first->vertex(j)->info().first];
								if (edge_indices.find({v1, v2}) == edge_indices.end())
								{
									edge_indices.insert({{v1, v2}, edge_count});
									Power_shape_data.edges_vertex_indices_.push_back(v1);
									Power_shape_data.edges_vertex_indices_.push_back(v2);
									edge_count++;
								}
								// Add face
								Power_shape_data.faces_edge_indices_.push_back(edge_indices[{v1, v2}]);
							}
						}
					}
				}
			}
		}
		uint32 power_nb_vertices = Power_shape_data.vertex_position_.size();
		uint32 power_nb_edges = Power_shape_data.edges_vertex_indices_.size() / 2;
		uint32 power_nb_faces = Power_shape_data.faces_nb_edges_.size();

		Power_shape_data.set_parameter(power_nb_vertices, power_nb_edges, power_nb_faces);

		import_incidence_graph_data(*mv, Power_shape_data);

		auto sphere_raidus = add_attribute<double, NonManifoldVertex>(*mv, "radius");
		for (uint32 i = 0u; i < power_point.size(); ++i)
		{
			uint32 vertex_id = inside_indices[point_info[i].first];
			(*sphere_raidus)[vertex_id] = power_point[i].weight();
		}

		std::shared_ptr<NonManifoldAttribute<Vec3>> mv_vertex_position =
			get_attribute<Vec3, NonManifoldVertex>(*mv, "position");
		if (mv_vertex_position)
			nonmanifold_provider_->set_mesh_bb_vertex_position(*mv, mv_vertex_position);

		nonmanifold_provider_->emit_connectivity_changed(*mv);
	}

	void compute_power_shape(SURFACE& surface)
	{
		NONMANIFOLD* mv = nonmanifold_provider_->add_mesh(surface_provider_->mesh_name(surface) + "_power_shape");
		Cgal_Surface_mesh csm;

		std::string filename = surface_provider_->mesh_filename(surface);
		if (!filename.empty())
		{
			std::ifstream input(filename);
			if (!input || !(input >> csm))
			{
				std::cerr << "Error: input file could not be read" << std::endl;
				return;
			}
		}
		Tree tree(faces(csm).first, faces(csm).second, csm);
		tree.accelerate_distance_queries();
		Delaunay tri = compute_delaunay_tredrahedron(csm,tree);
		std::vector<Weight_Point> Power_point;
		std::vector<std::pair<uint32, bool>> Point_info;
		std::unordered_map<uint32, uint32> Inside_indices;
		compute_inner_poles(tri, tree, Power_point, Point_info, Inside_indices);
		constrcut_power_diagram(mv, surface, Power_point, Point_info, Inside_indices);
	}

	void compute_original_power_diagram(SURFACE& surface)
	{
		NONMANIFOLD* mv = nonmanifold_provider_->add_mesh(surface_provider_->mesh_name(surface) + "_medial_axis");
		Cgal_Surface_mesh csm;

		std::string filename = surface_provider_->mesh_filename(surface);
		if (!filename.empty())
		{
			std::ifstream input(filename);
			if (!input || !(input >> csm))
			{
				std::cerr << "Error: input file could not be read" << std::endl;
				return;
			}
		}
		Tree tree(faces(csm).first, faces(csm).second, csm);
		tree.accelerate_distance_queries();
		Delaunay tri = compute_delaunay_tredrahedron(csm, tree);
		std::vector<Weight_Point> Power_point;
		std::vector<std::pair<uint32, bool>> Point_info;
		std::unordered_map<uint32, uint32> Inside_indices;
		compute_inner_voronoi(tri, tree, Power_point, Point_info, Inside_indices);
		constrcut_power_diagram(mv, surface, Power_point, Point_info, Inside_indices);
	}

	double compute_stability_ratio(IncidenceGraph& ig, NonManifoldEdge e)
	{
		//To do
		auto sphere_radius = get_attribute<double, NonManifoldVertex>(ig, "radius");
		auto position = get_attribute<double, NonManifoldVertex>(ig, "position");
		NonManifoldVertex [v1, v2] = (*ig.edge_incident_vertices_)[e.index_];
		Vec3 v1_p = value<Vec3>(ig, position, v1);
		Vec3 v2_p = value<Vec3>(ig, position, v2);
		double r1 = value<double>(ig, sphere_radius, v1);
		double r2 = value<double>(ig, sphere_radius, v2);
		double dis = std::max(0.0, ((v1_p - v2_p).squaredNorm() - std::abs(r1 - r2)));
		double stability = dis / (v1_p - v2_p).squaredNorm();
		return stability;
	}

 	void compute_stability_ratio(NONMANIFOLD& mv)
	{
		auto stability_ratio = add_attribute<double, NonManifoldEdge>(mv, "stability_ratio");
		auto stability_color = add_attribute<Vec3, NonManifoldEdge>(mv, "stability_color");
		foreach_cell(mv, [&](NonManifoldEdge e) -> bool {
			double stability = compute_stability_ratio(mv, e); 
			(*stability_ratio)[e.index_] = stability;
			(*stability_color)[e.index_] = Vec3(stability, (1 - stability), 0);
			return true;
		});
	}
	void collapse_non_manifold_using_QMat()
	{
		using EdgeQueue = std::multimap<Scalar, NonManifoldEdge>;
		using EdgeQueueIt = typename EdgeQueue::const_iterator;
		using EdgeInfo = std::pair<bool, EdgeQueueIt>; // {valid, iterator}

		EdgeQueue queue;
		auto edge_queue_it = add_attribute<EdgeInfo, NonManifoldEdge>(*non_manifold_, "__non_manifold_edge_queue_it");
		foreach_cell(*non_manifold_, [&](NonManifoldEdge e) -> bool {
			value<EdgeInfo>(*non_manifold_, edge_queue_it, e) = {
				true, queue.emplace(geometry::length(*non_manifold_, e, non_manifold_vertex_position_.get()), e)};
			return true;
		});

		using PositionAccu = std::vector<Vec3>;
		auto vertex_position_accu =
			add_attribute<PositionAccu, NonManifoldVertex>(*non_manifold_, "__non_manifold_vertex_position_accu");
		foreach_cell(*non_manifold_, [&](NonManifoldVertex v) -> bool {
			value<PositionAccu>(*non_manifold_, vertex_position_accu,
								v) = {value<Vec3>(*non_manifold_, non_manifold_vertex_position_, v)};
			return true;
		});

		while (!queue.empty())
		{
			auto it = queue.begin();
			NonManifoldEdge e = (*it).second;
			queue.erase(it);
			value<EdgeInfo>(*non_manifold_, edge_queue_it, e).first = false;

			std::vector<NonManifoldFace> ifaces = incident_faces(*non_manifold_, e);
			if (ifaces.size() == 0)
				continue;

			// iv[0] will be removed and iv[1] will survive
			std::vector<NonManifoldVertex> iv = incident_vertices(*non_manifold_, e);
			PositionAccu& accu0 = value<PositionAccu>(*non_manifold_, vertex_position_accu, iv[0]);
			PositionAccu& accu1 = value<PositionAccu>(*non_manifold_, vertex_position_accu, iv[1]);
			accu1.insert(accu1.end(), accu0.begin(), accu0.end());

			auto [v, removed_edges] = collapse_edge(*non_manifold_, e);
			for (NonManifoldEdge re : removed_edges)
			{
				EdgeInfo einfo = value<EdgeInfo>(*non_manifold_, edge_queue_it, re);
				if (einfo.first)
					queue.erase(einfo.second);
			}

			foreach_incident_edge(*non_manifold_, v, [&](NonManifoldEdge ie) -> bool {
				EdgeInfo einfo = value<EdgeInfo>(*non_manifold_, edge_queue_it, ie);
				if (einfo.first)
					queue.erase(einfo.second);
				value<EdgeInfo>(*non_manifold_, edge_queue_it, ie) = {
					true, queue.emplace(geometry::length(*non_manifold_, ie, non_manifold_vertex_position_.get()), ie)};
				return true;
			});
		}

		foreach_cell(*non_manifold_, [&](NonManifoldVertex v) -> bool {
			Vec3 mean{0, 0, 0};
			uint32 count = 0;
			for (Vec3& p : value<PositionAccu>(*non_manifold_, vertex_position_accu, v))
			{
				mean += p;
				++count;
			}
			mean /= count;
			value<Vec3>(*non_manifold_, non_manifold_vertex_position_, v) = mean;
			return true;
		});

		remove_attribute<NonManifoldEdge>(*non_manifold_, edge_queue_it);
		remove_attribute<NonManifoldVertex>(*non_manifold_, vertex_position_accu);

		non_manifold_provider_->emit_connectivity_changed(*non_manifold_);
		non_manifold_provider_->emit_attribute_changed(*non_manifold_, non_manifold_vertex_position_.get());
	}


protected:
	void init() override
	{
		surface_provider_ = static_cast<ui::MeshProvider<SURFACE>*>(
			app_.module("MeshProvider (" + std::string{mesh_traits<SURFACE>::name} + ")"));
		nonmanifold_provider_ = static_cast<ui::MeshProvider<NONMANIFOLD>*>(
			app_.module("MeshProvider (" + std::string{mesh_traits<NONMANIFOLD>::name} + ")"));
	}

	void left_panel() override
	{
		imgui_mesh_selector(surface_provider_, selected_surface_mesh_, "Surface", [&](SURFACE& m) {
			selected_surface_mesh_ = &m;
			surface_provider_->mesh_data(m).outlined_until_ = App::frame_time_ + 1.0;
		});
		imgui_mesh_selector(nonmanifold_provider_, selected_medial_axis, "Medial_axis", [&](NONMANIFOLD& nm) {
			selected_medial_axis = &nm;
			nonmanifold_provider_->mesh_data(nm).outlined_until_ = App::frame_time_ + 1.0;
		});


		if (selected_surface_mesh_)
		{
			if (ImGui::Button("Power shape"))
			{
				compute_power_shape(*selected_surface_mesh_);
				if (selected_medial_axis)
				{
					if (ImGui::Button("Compute stablility ratio"))
					{
						compute_stability_ratio(*selected_medial_axis);
					}
				}
			}
			if (ImGui::Button("Original Medial Axis"))
				compute_original_power_diagram(*selected_surface_mesh_);
		}
	}

private:
	SURFACE* selected_surface_mesh_;
	NONMANIFOLD* selected_medial_axis;
	MeshProvider<SURFACE>* surface_provider_;
	MeshProvider<NONMANIFOLD>* nonmanifold_provider_;
	Regular medial_axis;
};

} // namespace ui

} // namespace cgogn

#endif // CGOGN_MODULE_POWER_SHAPE_H_
