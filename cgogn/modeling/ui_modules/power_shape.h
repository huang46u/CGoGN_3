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
#include <filesystem>
#include <cmath>
#include <nlopt/src/api/nlopt.h>
#include <cgogn/core/ui_modules/mesh_provider.h>
#include <cgogn/ui/app.h>
#include <cgogn/ui/module.h>
#include <cgogn/core/types/maps/gmap/gmap_base.h>
#include <cgogn/geometry/types/slab_quadric.h>

#include <cgogn/io/point/point_import.h>
#include <cgogn/geometry/algos/medial_axis.h>
#include <cgogn/modeling/algos/decimation/SQEM_helper.h>
#include <cgogn/modeling/algos/decimation/QEM_helper.h>
#include <libacc/kd_tree.h>
#include <libacc/bvh_tree.h>
#include <cgogn/core/types/maps/cmap/cmap0.h>
// import CGAL
#include <CGAL/Bbox_3.h>
#include <CGAL/Delaunay_triangulation_3.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/IO/read_off_points.h>
#include <CGAL/Object.h>
#include <CGAL/Polygon_mesh_processing/distance.h>
#include <CGAL/Regular_triangulation_3.h>
#include <CGAL/Regular_triangulation_cell_base_3.h>
#include <CGAL/Regular_triangulation_vertex_base_3.h>
#include <CGAL/Side_of_triangle_mesh.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/Triangulation_cell_base_with_info_3.h>
#include <CGAL/Triangulation_vertex_base_with_info_3.h>
#include <CGAL/Min_sphere_of_points_d_traits_3.h>
#include <CGAL/Min_sphere_of_spheres_d.h>
#include <CGAL/double.h>
#include <CGAL/optimal_bounding_box.h>

#include <Highs.h>

#include <iomanip>
#include <limits>

namespace cgogn
{

namespace ui
{

template <typename POINT, typename SURFACE, typename NONMANIFOLD>
class PowerShape : public Module
{
	static_assert(mesh_traits<SURFACE>::dimension >= 2, "PowerShape can only be used with meshes of dimension >= 2");
	// Kernel for construct Delaunay
	using K = CGAL::Exact_predicates_inexact_constructions_kernel;
	using Point = K::Point_3;
	using Weight_Point = K::Weighted_point_3;

	class VertexInfo
	{
	public:
		bool inside = false;
		uint32 id = -1;
		Point inside_pole;
		Point outside_pole;
		double inside_pole_distance = 0.0;
		double outside_pole_distance = 0.0;
	};

	class DelaunayCellInfo
	{
	public:
		int32 id = -1;
		bool is_pole = false;
		bool inside = false;
		bool angle_flag = false;
		bool radius_flag = false;
		bool distance_flag = true;
		bool selected = false;
		Point centroid;
		double angle;
		double radius2; // radius square
	};
	class RegularVertexInfo
	{
	public:
		bool inside = false;
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
	// Regular
	using Vb0 = CGAL::Regular_triangulation_vertex_base_3<K>;
	using RVb = CGAL::Triangulation_vertex_base_with_info_3<RegularVertexInfo, K, Vb0>;
	using RCb = CGAL::Regular_triangulation_cell_base_3<K>;
	using RTds = CGAL::Triangulation_data_structure_3<RVb, RCb>;
	using Regular = CGAL::Regular_triangulation_3<K, RTds, CGAL::Fast_location>;
	using Regular_Cell_handle = typename Regular::Cell_handle;

	using Cgal_Surface_mesh = CGAL::Surface_mesh<Point>;
	using Point_inside = CGAL::Side_of_triangle_mesh<Cgal_Surface_mesh, K>;
	using Primitive = CGAL::AABB_face_graph_triangle_primitive<Cgal_Surface_mesh>;
	using Tree_Traits = CGAL::AABB_traits<K, Primitive>;
	using Tree = CGAL::AABB_tree<Tree_Traits>;

	using Min_Sphere_Tratis =  CGAL::Min_sphere_of_points_d_traits_3<K, K::FT>;
	using Min_sphere =  CGAL::Min_sphere_of_spheres_d<Min_Sphere_Tratis>;

	template <typename T>
	using NonManifoldAttribute = typename mesh_traits<NONMANIFOLD>::template Attribute<T>;
	template <typename T>
	using SurfaceAttribute = typename mesh_traits<SURFACE>::template Attribute<T>;

	using PointVertex = typename mesh_traits<POINT>::Vertex;

	using NonManifoldVertex = typename mesh_traits<NONMANIFOLD>::Vertex;
	using NonManifoldEdge = typename mesh_traits<NONMANIFOLD>::Edge;
	using NonManifoldFace = typename mesh_traits<NONMANIFOLD>::Face;

	using SurfaceVertex = typename mesh_traits<SURFACE>::Vertex;
	using SurfaceEdge = typename mesh_traits<SURFACE>::Edge;
	using SurfaceFace = typename mesh_traits<SURFACE>::Face;


	using Vec4 = geometry::Vec4;
	using Vec3 = geometry::Vec3;
	using Scalar = geometry::Scalar;

public:
	PowerShape(const App& app) : Module(app, "PowerShape"), tree(nullptr), inside_tester(nullptr)
	{
		clustering_mode = 0;
	}
	~PowerShape()
	{
		if (!tree)
			delete tree;
		if (!inside_tester)
			delete inside_tester;
		
	}

private:
	
	template <typename T>
	struct T3_hash
	{
		std::size_t operator()(const T& p) const
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
	struct edge_equal
	{
		bool operator()(const std::pair<uint32, uint32>& edge1, const std::pair<uint32, uint32>& edge2) const
		{
			return ((edge1.first == edge2.first && edge1.second == edge2.second) ||
					(edge1.first == edge2.second && edge1.second == edge2.first));
		}
	};
	bool inside_sphere(const Vec3& point, const Vec3& center, double radius)
	{
		return (point - center).norm() <= radius;
	}
	bool pointInside(Tree& tree, Point& query)
	{
		// Initialize the point-in-polyhedron tester
		Point_inside inside_tester(tree);

		// Determine the side and return true if inside!
		return inside_tester(query) == CGAL::ON_BOUNDED_SIDE;
	}
	bool inside(Point_inside& inside_tester, Point& query)
	{
		// Determine the side and return true if inside!
		return inside_tester(query) == CGAL::ON_BOUNDED_SIDE;
	}

	void load_model_in_cgal(SURFACE& surface, Cgal_Surface_mesh& csm)
	{
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
	}

	struct Cell_handle_hash
	{
		std::size_t operator()(const Delaunay_Cell_handle& cell) const
		{
			return std::hash<ptrdiff_t>()(cell->info().id);
		}
	};

	struct Cell_handle_equal
	{
		bool operator()(const Delaunay_Cell_handle& cell1, const Delaunay_Cell_handle& cell2) const
		{
			return cell1->info().id == cell2->info().id;
		}
	};

	void mark_poles(Delaunay& tri)
	{
		std::vector<Delaunay_Cell_handle> vertex_incident_cells;
		for (auto vit = tri.finite_vertices_begin(); vit != tri.finite_vertices_end(); ++vit)
		{
			vertex_incident_cells.clear();
			Delaunay_Cell_handle ch_inside, ch_outside;
			bool inside = false, outside = false;
			tri.finite_incident_cells(vit, std::back_inserter(vertex_incident_cells));
			for (auto cell : vertex_incident_cells)
			{
				if (cell->info().inside)
				{
					if (cell->info().radius2 > vit->info().inside_pole_distance)
					{
						vit->info().inside_pole_distance = cell->info().radius2;
						inside = true;
						ch_inside = cell;
					}
				}
				else
				{
					if (cell->info().radius2 > vit->info().outside_pole_distance)
					{
						vit->info().outside_pole_distance = cell->info().radius2;
						outside = true;
						ch_outside = cell;
					}
				}
			}
			
			if (inside)
			{
				ch_inside->info().is_pole = true;
				vit->info().inside_pole = ch_inside->info().centroid;
				vit->info().inside_pole_distance = std::sqrt(ch_inside->info().radius2);
			}
			if (outside)
			{
				ch_outside->info().is_pole = true;
				vit->info().outside_pole = ch_outside->info().centroid;
				vit->info().outside_pole_distance = std::sqrt(ch_outside->info().radius2);
			}
		}
	}

	void filter_by_angle(Delaunay& tri, double angle_threshold)
	{
		uint32 count = 0;
		for (auto cit = tri.finite_cells_begin(); cit != tri.finite_cells_end(); ++cit)
		{
			double max_angle = 0;
			for (size_t i = 0; i < 4; i++)
			{
				for (size_t j = i + 1; j < 4; j++)
				{
					Point p1 = cit->vertex(i)->point();
					Point p2 = cit->vertex(j)->point();
					Point p3 = cit->info().centroid;
					Vec3 v1 = Vec3(p1.x(), p1.y(), p1.z()) - Vec3(p3.x(), p3.y(), p3.z());
					Vec3 v2 = Vec3(p2.x(), p2.y(), p2.z()) - Vec3(p3.x(), p3.y(), p3.z());
					double angle = std::acos(v1.dot(v2) / (std::sqrt(v1.dot(v1)) * std::sqrt(v2.dot(v2))));
					max_angle = std::max(angle, max_angle);
				}
			}
			cit->info().angle = max_angle;
			cit->info().angle_flag = max_angle > angle_threshold ? true : false;
			if (!cit->info().angle_flag)
				count++;
			
		}
		std::cout << "angle deletion " << count << std::endl;
	}

	void filter_by_circumradius(Delaunay& tri, double radius_threshold)
	{
		uint32 count = 0;
		for (auto cit = tri.finite_cells_begin(); cit != tri.finite_cells_end(); ++cit)
		{
			cit->info().radius_flag = std::sqrt(cit->info().radius2) > radius_threshold
										  ? (std::sqrt(cit->info().radius2) < max_radius_ ? true : false)
										  : false;
			if (!cit->info().radius_flag)
				count++;
		}
		std::cout << "circumradius deletion: " << count << std::endl;
	}

	void filter_by_distance(Delaunay& tri, double distance_threshold)
	{

		double min_dist = std::numeric_limits<double>::max();
		double max_dist = std::numeric_limits<double>::min();
		int count = 0;
		for (auto fit = tri.finite_facets_begin(); fit != tri.finite_facets_end(); ++fit)
		{
			auto opposite = tri.mirror_facet(*fit);
			if (fit->first->info().distance_flag && opposite.first->info().distance_flag && fit->first->info().inside &&
				opposite.first->info().inside)
			{
				double distance =
					std::sqrt(CGAL::squared_distance(fit->first->info().centroid, opposite.first->info().centroid));
				min_dist = std::min(distance, min_dist);
				max_dist = std::max(max_dist, distance);
				if (distance <= distance_threshold)
				{
					if (fit->first->info().radius2 >= opposite.first->info().radius2)
					{
						opposite.first->info().distance_flag = false;
					}
					else
					{
						fit->first->info().distance_flag = false;
					}
					count++;
				}
			}
		}
		std::cout << "distance deletion: " << count << std::endl;
	}
	
	void normalise_scalar(SURFACE& surface, std::shared_ptr<SurfaceAttribute<Scalar>> attribute)
	{
		Scalar min = 1e30;
		Scalar max = -1e30;
		foreach_cell(surface, [&](SurfaceVertex sv) {
			Scalar s = value<Scalar>(surface, attribute, sv);
			if (s > 1e5)
			{
				return true;
			}
			if (s < min)
				min = s;
			if (s > max)
				max = s;
			return true;
		});
		std::cout << "min: " << min << ", max: " << max << std::endl;
		parallel_foreach_cell(surface, [&](SurfaceVertex sv) {
			Scalar s = value<Scalar>(surface, attribute, sv);
			
			if (s > 1e5)
			{
				value<Scalar>(surface, attribute, sv) = 1;
				return true;
			}
			value<Scalar>(surface, attribute, sv) = (s - min) / (max - min);
			return true;
		});
		surface_provider_->emit_attribute_changed(surface, attribute.get());
	}

	std::pair<Vec3, Scalar> sphere_fitting_algo(SURFACE& surface, std::vector<SurfaceVertex> cluster_vertices)
	{
		auto position = get_attribute<Vec3, SurfaceVertex>(surface, "position");
		Eigen::Matrix3d A;
		Eigen::Vector3d B;
		Scalar x_avg = 0;
		Scalar y_avg = 0;
		Scalar z_avg = 0;
		Scalar xy_avg = 0;
		Scalar xz_avg = 0;
		Scalar yz_avg = 0;
		Scalar x2_avg = 0;
		Scalar y2_avg = 0;
		Scalar z2_avg = 0;
		Scalar x2y_avg = 0;
		Scalar x2z_avg = 0;
		Scalar xy2_avg = 0;
		Scalar y2z_avg = 0;
		Scalar xz2_avg = 0;
		Scalar yz2_avg = 0;
		Scalar zy2_avg = 0;
		Scalar x3_avg = 0;
		Scalar y3_avg = 0;
		Scalar z3_avg = 0;
		for (auto sv : cluster_vertices)
		{
			Vec3 pos = value<Vec3>(surface, position, sv);
			x_avg += pos.x();
			y_avg += pos.y();
			z_avg += pos.z();
			xy_avg += pos.x() * pos.y();
			xz_avg += pos.x() * pos.z();
			yz_avg += pos.y() * pos.z();
			x2_avg += pos.x() * pos.x();
			y2_avg += pos.y() * pos.y();
			z2_avg += pos.z() * pos.z();
			x2y_avg += pos.x() * pos.x() * pos.y();
			x2z_avg += pos.x() * pos.x() * pos.z();
			xy2_avg += pos.x() * pos.y() * pos.y();
			y2z_avg += pos.y() * pos.y() * pos.z();
			xz2_avg += pos.x() * pos.z() * pos.z();
			yz2_avg += pos.y() * pos.z() * pos.z();
			zy2_avg += pos.z() * pos.y() * pos.y();
			x3_avg += pos.x() * pos.x() * pos.x();
			y3_avg += pos.y() * pos.y() * pos.y();
			z3_avg += pos.z() * pos.z() * pos.z();
		}
		x_avg /= cluster_vertices.size();
		y_avg /= cluster_vertices.size();
		z_avg /= cluster_vertices.size();
		xy_avg /= cluster_vertices.size();
		xz_avg /= cluster_vertices.size();
		yz_avg /= cluster_vertices.size();
		x2_avg /= cluster_vertices.size();
		y2_avg /= cluster_vertices.size();
		z2_avg /= cluster_vertices.size();
		x2y_avg /= cluster_vertices.size();
		x2z_avg /= cluster_vertices.size();
		xy2_avg /= cluster_vertices.size();
		y2z_avg /= cluster_vertices.size();
		xz2_avg /= cluster_vertices.size();
		yz2_avg /= cluster_vertices.size();
		zy2_avg /= cluster_vertices.size();
		x3_avg /= cluster_vertices.size();
		y3_avg /= cluster_vertices.size();
		z3_avg /= cluster_vertices.size();
		A << x2_avg - x_avg * x_avg, xy_avg - x_avg * y_avg, xz_avg - x_avg * z_avg,
			xy_avg - x_avg * y_avg, y2_avg - y_avg * y_avg, yz_avg - y_avg * z_avg,
			xz_avg - x_avg * z_avg, yz_avg - y_avg * z_avg, z2_avg - z_avg * z_avg;
		B << (x3_avg - x_avg * x2_avg) +
			(x_avg * y2_avg - x_avg * y2_avg) + 
			(xz2_avg - x_avg*z2_avg),
			(x2y_avg - x2_avg * y_avg)
			+ (y3_avg - y_avg * y2_avg)
			+ yz2_avg-y_avg*z2_avg,
			(x2z_avg- x2_avg * z_avg)
			 +(zy2_avg - z_avg*y2_avg)
			+(z3_avg - z_avg *z2_avg);
		Eigen::Vector3d center = A.inverse() *0.5 * B;
		Scalar radius = std::sqrt((x2_avg - 2 * x_avg * center.x() + center.x() * center.x() +
								  y2_avg - 2 * y_avg * center.y() + center.y() * center.y() +
								  z2_avg - 2 * z_avg * center.z() + center.z() * center.z()));
		return {Vec3(center[0], center[1], center[2]), radius};
	}

 public:

	std::array<std::array<double, 3>, 8> compute_big_box(SURFACE& surface, Cgal_Surface_mesh& csm)
	{
		std::array<Point, 8> obb_points;
		Point acc(0, 0, 0);
		CGAL::oriented_bounding_box(csm, obb_points, CGAL::parameters::use_convex_hull(true));
		for (size_t i = 0; i < obb_points.size(); i++)
		{
			acc += K::Vector_3(obb_points[i].x(), obb_points[i].y(), obb_points[i].z());
		}
		std::array<double, 3> center{acc.x() / 8, acc.y() / 8, acc.z() / 8};
		// Create a large box surrounding object so that the Voronoi vertices are bounded
		MeshData<SURFACE>& md = surface_provider_->mesh_data(surface);
		double offset = (md.bb_max_ - md.bb_min_).norm();
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
		return cube_corners;
	}

	void plot_surface_samples(SURFACE& surface, std::vector<Point>& mesh_samples)
	{
		cgogn::io::PointImportData samples;
		surface_sample_ = point_provider_->add_mesh(point_provider_->mesh_name(surface) + "surface_samples");
		uint32 nb_vertices = mesh_samples.size();
		samples.reserve(nb_vertices);
		for (auto& s : mesh_samples)
		{
			samples.vertex_position_.emplace_back(s[0], s[1], s[2]);
		}
		cgogn::io::import_point_data(*surface_sample_, samples);
		auto position = get_attribute<Vec3, PointVertex>(*surface_sample_, "position");
		if (position)
			point_provider_->set_mesh_bb_vertex_position(*surface_sample_, position);
	}

	

	Delaunay compute_delaunay_tredrahedron(SURFACE& surface, Cgal_Surface_mesh& csm, Tree& tree)
	{

		std::vector<Point> Delaunay_tri_point;
		// Sampling the mesh surface
		std::vector<Point> mesh_samples;
		CGAL::Polygon_mesh_processing::sample_triangle_mesh(
			csm, std::back_inserter(mesh_samples),
			// CGAL::parameters::use_monte_carlo_sampling(true).number_of_points_per_area_unit(50));
			CGAL::parameters::use_grid_sampling(true).grid_spacing(1));


 		std::array<std::array<double, 3>, 8> cube_corners = compute_big_box(surface, csm);
 		// 	Add bounding box vertices in the sample points set
 		for (auto& p : cube_corners)
 		{
 			Delaunay_tri_point.emplace_back(p[0], p[1], p[2]);
			mesh_samples.emplace_back(p[0], p[1], p[2]);
 		}

		// Add sampled vertices into the volume data to construct the delauney tredrahedron
		for (auto& s : mesh_samples)
		{
			Delaunay_tri_point.emplace_back(s[0], s[1], s[2]);
		}

		plot_surface_samples(surface, mesh_samples);

		// auto start_timer = std::chrono::high_resolution_clock::now();

		// Construct delauney tredrahedron using CGAL
		int count = 0;
		Delaunay tri;
		for (Point p : Delaunay_tri_point)
		{
			Delaunay_Vertex_handle vh = tri.insert(p);
			vh->info().id = count;
			count++;
		}
		int cell_count = 0;
		for (auto cit = tri.finite_cells_begin(); cit != tri.finite_cells_end(); ++cit)
		{
			cit->info().id = -1;
			cit->info().centroid = CGAL::circumcenter(tri.tetrahedron(cit));
			cit->info().radius2 = CGAL::squared_distance(cit->info().centroid, cit->vertex(0)->point());

			if (pointInside(tree, cit->info().centroid))
			{
				cit->info().inside = true;
				cit->info().id = cell_count;
				min_radius_ = std::min(min_radius_, std::sqrt(cit->info().radius2));
				max_radius_ = std::max(max_radius_, std::sqrt(cit->info().radius2));
			}
			else
			{
				cit->info().inside = false;
				cit->info().id = cell_count;
			}
			cell_count++;
		}
		mark_poles(tri);
		return tri;
	}


	void construct_voronoi_cell(SURFACE& surface, Delaunay& tri)
	{
		NONMANIFOLD* voronoi_cell =
			nonmanifold_provider_->add_mesh(surface_provider_->mesh_name(surface) + "voronoi cells");
		cgogn::io::IncidenceGraphImportData Cell_non_manifold;
		std::unordered_map<std::pair<uint32, uint32>, size_t, edge_hash, edge_equal> edge_indices;
		std::vector<double> sphere_radius;
		std::vector<Vec3> sphere_color;
		bool valid;
		bool inside;
		bool outside;
		int vertex_count = 0;
		int face_number = 0;
		std::vector<Delaunay_Cell_handle> incells;
		
		for (auto eit = tri.finite_edges_begin(); eit != tri.finite_edges_end(); ++eit)
		{
			valid = false;
			incells.clear();
			inside = false;
			outside = false;
			Delaunay_Cell_circulator cc = tri.incident_cells(*eit);
			do
			{
				
				if (cc->info().radius_flag)
				{
					if (cc->info().inside && (distance_filtering_ ? cc->info().distance_flag : true) &&
						(angle_filtering_ ? cc->info().angle_flag : true) &&
						(pole_filtering_ ? cc->info().is_pole : true))
					{
						cc->info().selected = true;
						valid = true;
					}
					incells.push_back(cc);
					if (cc->info().inside)
						inside = true;
					else
						outside = true;
				}
				
				cc++;
			} while (cc != tri.incident_cells(*eit));
			if (!valid || !inside || !outside || incells.size()<3)
				continue;
			for (size_t k = 0; k < incells.size(); k++)
			{
				if (incells[k]->info().selected)
				{
					sphere_color.push_back(Vec3(1.0, 0.0, 0.0));
				}
				else
				{
					sphere_color.push_back(Vec3(0.0, 1.0, 0.0));
				}
				Point centroid = incells[k]->info().centroid;
				double radius = std::sqrt(incells[k]->info().radius2);
				Cell_non_manifold.vertex_position_.emplace_back(centroid[0], centroid[1], centroid[2]);
				sphere_radius.push_back(radius);
				incells[k]->info().id = vertex_count;
				vertex_count++;
				
			}
			Cell_non_manifold.faces_nb_edges_.push_back(incells.size());
			for (size_t k = 0; k < incells.size(); k++)
			{
				uint32 ev1 = incells[k]->info().id;
				uint32 ev2 = incells[(k + 1) % incells.size()]->info().id;
				if (edge_indices.find({ev1, ev2}) == edge_indices.end())
				{
					Cell_non_manifold.edges_vertex_indices_.push_back(ev1);
					Cell_non_manifold.edges_vertex_indices_.push_back(ev2);
					edge_indices.insert({{ev1, ev2}, edge_indices.size()});
				}
				
				Cell_non_manifold.faces_edge_indices_.push_back(edge_indices[{ev1,ev2}]);
			}
			face_number++;
		}

		uint32 Cell_non_manifold_nb_vertices = Cell_non_manifold.vertex_position_.size();
		uint32 Cell_non_manifold_nb_edges = Cell_non_manifold.edges_vertex_indices_.size() / 2;
		uint32 Cell_non_manifold_nb_faces = face_number;
		Cell_non_manifold.reserve(Cell_non_manifold_nb_vertices, Cell_non_manifold_nb_edges,
										   Cell_non_manifold_nb_faces);

		import_incidence_graph_data(*voronoi_cell, Cell_non_manifold);

		auto sphere_raidus_att = add_attribute<double, NonManifoldVertex>(*voronoi_cell, "sphere_radius");
		auto sphere_color_att = add_attribute<Vec3, NonManifoldVertex>(*voronoi_cell, "sphere_color");
		for (size_t idx = 0; idx < sphere_radius.size(); idx++)
		{
			(*sphere_raidus_att)[idx] = sphere_radius[idx];
			(*sphere_color_att)[idx] = sphere_color[idx];
		}
		std::shared_ptr<NonManifoldAttribute<Vec3>> mv_vertex_position =
			get_attribute<Vec3, NonManifoldVertex>(*voronoi_cell, "position");
		if (mv_vertex_position)
			nonmanifold_provider_->set_mesh_bb_vertex_position(*voronoi_cell, mv_vertex_position);

		nonmanifold_provider_->emit_connectivity_changed(*voronoi_cell);

	}

	void construct_candidates_points(SURFACE& surface, Delaunay& tri)
	{
		cgogn::io::PointImportData candidates;
		cgogn::io::PointImportData outside_poles_data;
		POINT* candidates_point = point_provider_->add_mesh(point_provider_->mesh_name(surface) + "candidates");
		POINT* outside_poles_point = point_provider_->add_mesh(point_provider_->mesh_name(surface) + "outside_poles");
		

		std::vector<double> candidates_radius;
		std::vector<double> angle;

		std::vector<Vec3> oustside_poles_center;
		std::vector<double> oustside_poles_radius;

		if (distance_filtering_)
			filter_by_distance(tri, distance_threshold_);
		if (circumradius_filtering_)
			filter_by_circumradius(tri, radius_threshold_);
		if (angle_filtering_)
			filter_by_angle(tri, angle_threshold_);

		construct_voronoi_cell(surface, tri);

		for (auto cit = tri.finite_cells_begin(); cit != tri.finite_cells_end(); ++cit)
		{
			if (cit->info().inside && (distance_filtering_ ? cit->info().distance_flag : true) &&
				(circumradius_filtering_ ? cit->info().radius_flag : true) &&
				(angle_filtering_ ? cit->info().angle_flag : true) && (pole_filtering_ ? cit->info().is_pole : true))
			{
				candidates.vertex_position_.emplace_back(cit->info().centroid.x(), cit->info().centroid.y(),
														 cit->info().centroid.z());
				candidates_radius.push_back(std::sqrt(cit->info().radius2));
				angle.push_back(cit->info().angle);

				double min_dist = std::numeric_limits<double>::max();
				Vec3 nearest_pole(0,0,0);
				for (size_t i = 0; i < 4; ++i)
				{
					if (min_dist > cit->vertex(i)->info().outside_pole_distance)
					{
						min_dist = cit->vertex(i)->info().outside_pole_distance;
						nearest_pole =
							Vec3(cit->vertex(i)->info().outside_pole.x(), 
												cit->vertex(i)->info().outside_pole.y(), 
												cit->vertex(i)->info().outside_pole.z());
					}
				}	
				oustside_poles_radius.emplace_back(min_dist);
				oustside_poles_center.push_back(nearest_pole);
				outside_poles_data.vertex_position_.emplace_back(nearest_pole.x(), nearest_pole.y(), nearest_pole.z());
			}
		}
		candidates.reserve(candidates_radius.size());
		outside_poles_data.reserve(oustside_poles_radius.size());
		cgogn::io::import_point_data(*candidates_point, candidates);
		cgogn::io::import_point_data(*outside_poles_point, outside_poles_data);
		auto candidates_position = get_attribute<Vec3, PointVertex>(*candidates_point, "position");
		auto outside_poles_position = get_attribute<Vec3, PointVertex>(*outside_poles_point, "position");
		if (candidates_position)
		{
			point_provider_->set_mesh_bb_vertex_position(*candidates_point, candidates_position);
		}
		if (outside_poles_position)
		{
			point_provider_->set_mesh_bb_vertex_position(*outside_poles_point, outside_poles_position);
		}
		auto sphere_radius_att = add_attribute<double, PointVertex>(*candidates_point, "sphere_radius");
		auto point_angle_att = add_attribute<Vec3, PointVertex>(*candidates_point, "angle");
		auto outside_pole_radius_att = add_attribute<double, PointVertex>(*candidates_point, "outside_pole_radius");
		auto outside_pole_center_att = add_attribute<Vec3, PointVertex>(*candidates_point, "outside_pole_center");
		for (size_t idx = 0; idx < candidates_radius.size(); idx++)
		{
			(*sphere_radius_att)[idx] = candidates_radius[idx];
			(*point_angle_att)[idx] = Vec3(((angle[idx] - angle_threshold_) / (M_PI - angle_threshold_)), 0, 0);
			(*outside_pole_radius_att)[idx] = oustside_poles_radius[idx];
			(*outside_pole_center_att)[idx] = oustside_poles_center[idx];
		}
		construct_candidates_delaunay(surface, *candidates_point);
		Delaunay RVD = compute_restricted_voronoi_diagram(surface, *candidates_point);
		construct_inner_restricted_voronoi_diagram(RVD, "constrained_voronoi_diagram");
		std::cout << "point size: " << candidates_radius.size() << std::endl;
	}

	void construct_candidates_delaunay(SURFACE& surface, POINT& candidates)
	{
		Delaunay tri;
		auto position = get_attribute<Vec3, PointVertex>(candidates, "position");
		auto [bb_min, bb_max] = geometry::bounding_box(*position);
		std::vector<Vec3> bb_box;
		bb_box.emplace_back(bb_min[0], bb_max[1], bb_max[2]);
		bb_box.emplace_back(bb_max[0], bb_max[1], bb_max[2]);
		bb_box.emplace_back(bb_min[0], bb_min[1], bb_max[2]);
		bb_box.emplace_back(bb_max[0], bb_min[1], bb_max[2]);
		bb_box.emplace_back(bb_min[0], bb_max[1], bb_min[2]);
		bb_box.emplace_back(bb_max[0], bb_max[1], bb_min[2]);
		bb_box.emplace_back(bb_min[0], bb_min[1], bb_min[2]);
		bb_box.emplace_back(bb_max[0], bb_min[1], bb_min[2]);
		for (Vec3 vec : bb_box)
			tri.insert(Point(vec.x(), vec.y(), vec.z()));
		foreach_cell(candidates, [&](PointVertex v) {
			tri.insert(Point(value<Vec3>(candidates, position, v).x(), value<Vec3>(candidates, position, v).y(),
							 value<Vec3>(candidates, position, v).z()));
			return true;
		});

		int count = 0;
		for (auto cit = tri.finite_cells_begin(); cit != tri.finite_cells_end(); ++cit)
		{
			cit->info().centroid = CGAL::circumcenter(tri.tetrahedron(cit));
			cit->info().radius2 = CGAL::squared_distance(cit->info().centroid, cit->vertex(0)->point());
			cit->info().id = count;
			count++;
		}
		mark_poles(tri);
		NONMANIFOLD* candidates_voronoi =
			nonmanifold_provider_->add_mesh(surface_provider_->mesh_name(surface) + "candidates_voronoi");
		cgogn::io::IncidenceGraphImportData Candidates_voronoi_non_manifold;
		std::unordered_map<std::pair<uint32, uint32>, size_t, edge_hash, edge_equal> edge_indices;
		std::vector<Delaunay_Cell_handle> incells;
		bool all_finite_inside = true;
		int face_number = 0, vertex_count = 0;
		for (auto vit = tri.finite_vertices_begin(); vit != tri.finite_vertices_end(); ++vit)
		{
			if (vit->info().outside_pole_distance < 0.2)
			{

				Candidates_voronoi_non_manifold.vertex_position_.emplace_back(vit->point().x(), vit->point().y(),
																			  vit->point().z());
				Candidates_voronoi_non_manifold.vertex_position_.emplace_back(
					vit->info().outside_pole.x(), vit->info().outside_pole.y(), vit->info().outside_pole.z());
				Candidates_voronoi_non_manifold.edges_vertex_indices_.emplace_back(vertex_count);
				Candidates_voronoi_non_manifold.edges_vertex_indices_.emplace_back(vertex_count + 1);
				vertex_count += 2;
			}
		}

		uint32 Cell_non_manifold_nb_vertices = Candidates_voronoi_non_manifold.vertex_position_.size();
		uint32 Cell_non_manifold_nb_edges = Candidates_voronoi_non_manifold.edges_vertex_indices_.size() / 2;
		uint32 Cell_non_manifold_nb_faces = face_number;
		Candidates_voronoi_non_manifold.reserve(Cell_non_manifold_nb_vertices, Cell_non_manifold_nb_edges,
										Cell_non_manifold_nb_faces);

		import_incidence_graph_data(*candidates_voronoi, Candidates_voronoi_non_manifold);
		std::shared_ptr<NonManifoldAttribute<Vec3>> mv_vertex_position =
			get_attribute<Vec3, NonManifoldVertex>(*candidates_voronoi, "position");
		if (mv_vertex_position)
			nonmanifold_provider_->set_mesh_bb_vertex_position(*candidates_voronoi, mv_vertex_position);

		nonmanifold_provider_->emit_connectivity_changed(*candidates_voronoi);

	}
	NONMANIFOLD& compute_initial_non_manifold(Delaunay& tri, string name)
	{
		std::vector<double> sphere_radius;
		std::vector<Point> sphere_center;
		cgogn::io::IncidenceGraphImportData Initial_non_manifold;
		std::unordered_map<std::pair<uint32, uint32>, size_t, edge_hash, edge_equal> edge_indices;
		NONMANIFOLD* mv =
			nonmanifold_provider_->add_mesh(name + std::to_string(nonmanifold_provider_->number_of_meshes()));
		// Add vertices
		int count = 0;
		for (auto cit = tri.finite_cells_begin(); cit != tri.finite_cells_end(); ++cit)
		{
			if (cit->info().inside)
			{
				cit->info().id = count;
				Point centroid = cit->info().centroid;
				double radius = std::sqrt(cit->info().radius2);
				Initial_non_manifold.vertex_position_.emplace_back(centroid[0], centroid[1], centroid[2]);
				sphere_radius.push_back(radius);
				sphere_center.push_back(centroid);
				count++;
			}
		}
		// Add edges
		for (auto fit = tri.finite_facets_begin(); fit != tri.finite_facets_end(); ++fit)
		{
			if (fit->first->info().inside && tri.mirror_facet(*fit).first->info().inside)
			{
				uint32 v1_ind = fit->first->info().id;
				uint32 v2_ind = tri.mirror_facet(*fit).first->info().id;
				Initial_non_manifold.edges_vertex_indices_.push_back(v1_ind);
				Initial_non_manifold.edges_vertex_indices_.push_back(v2_ind);
				edge_indices.insert({{v1_ind, v2_ind}, edge_indices.size()});
			}
		}
		bool all_finite_inside;

		std::vector<Delaunay_Cell_handle> incells;
		for (auto eit = tri.finite_edges_begin(); eit != tri.finite_edges_end(); ++eit)
		{
			all_finite_inside = true;
			incells.clear();
			Delaunay_Cell_circulator cc = tri.incident_cells(*eit);
			do
			{
				if (tri.is_infinite(cc))
				{
					all_finite_inside = false;
					break;
				}
				else if (cc->info().inside == false)
				{
					all_finite_inside = false;
					break;
				}
				incells.push_back(cc);
				cc++;
			} while (cc != tri.incident_cells(*eit));
			if (!all_finite_inside)
				continue;
			for (size_t k = 2; k < incells.size() - 1; k++)
			{
				uint32 ev1 = incells[0]->info().id;
				uint32 ev2 = incells[k]->info().id;
				// Check if the edge is already added
				if (edge_indices.find({ev1, ev2}) == edge_indices.end())
				{
					Initial_non_manifold.edges_vertex_indices_.push_back(ev1);
					Initial_non_manifold.edges_vertex_indices_.push_back(ev2);
					edge_indices.insert({{ev1, ev2}, edge_indices.size()});
				}
			}
			for (size_t k = 1; k < incells.size() - 1; k++)
			{
				uint32 v1 = incells[0]->info().id;
				uint32 v2 = incells[k]->info().id;
				uint32 v3 = incells[k + 1]->info().id;
				uint32 e1, e2, e3;
				e1 = edge_indices[{v1, v2}];
				e2 = edge_indices[{v2, v3}];
				e3 = edge_indices[{v3, v1}];

				Initial_non_manifold.faces_nb_edges_.push_back(3);
				Initial_non_manifold.faces_edge_indices_.push_back(e1);
				Initial_non_manifold.faces_edge_indices_.push_back(e2);
				Initial_non_manifold.faces_edge_indices_.push_back(e3);
			}
		}

		uint32 Initial_non_manifold_nb_vertices = Initial_non_manifold.vertex_position_.size();
		uint32 Initial_non_manifold_nb_edges = Initial_non_manifold.edges_vertex_indices_.size() / 2;
		uint32 Initial_non_manifold_nb_faces = Initial_non_manifold.faces_nb_edges_.size();
		Initial_non_manifold.reserve(Initial_non_manifold_nb_vertices, Initial_non_manifold_nb_edges,
									 Initial_non_manifold_nb_faces);

		import_incidence_graph_data(*mv, Initial_non_manifold);

		auto sphere_raidus = add_attribute<double, NonManifoldVertex>(*mv, "sphere_radius");
		auto sphere_info = add_attribute<Vec4, NonManifoldVertex>(*mv, "sphere_info");
		for (size_t idx = 0; idx < sphere_center.size(); idx++)
		{
			(*sphere_raidus)[idx] = sphere_radius[idx];
			(*sphere_info)[idx] =
				Vec4(sphere_center[idx].x(), sphere_center[idx].y(), sphere_center[idx].z(), sphere_radius[idx]);
		}
		std::shared_ptr<NonManifoldAttribute<Vec3>> mv_vertex_position =
			get_attribute<Vec3, NonManifoldVertex>(*mv, "position");
		if (mv_vertex_position)
			nonmanifold_provider_->set_mesh_bb_vertex_position(*mv, mv_vertex_position);

		nonmanifold_provider_->emit_connectivity_changed(*mv);
		return *mv;
	}

	Regular compute_regular_tredrahedron(Tree& tree, Delaunay& tri)
	{
		Regular power_shape;
		for (auto cit = tri.finite_cells_begin(); cit != tri.finite_cells_end(); ++cit)
		{
			if (cit->info().inside && (distance_filtering_ ? cit->info().distance_flag : true) &&
				(circumradius_filtering_ ? cit->info().radius_flag : true) &&
				(angle_filtering_ ? cit->info().angle_flag : true) && cit->info().is_pole)
			{
				Weight_Point wp(cit->info().centroid, cit->info().radius2);
				power_shape.insert(wp);
				
			}
		}
		int count = 0;
		for (auto vit = power_shape.finite_vertices_begin(); vit != power_shape.finite_vertices_end(); ++vit)
		{
			// if the point is inside
			Point p = vit->point().point();
			if (pointInside(tree, p))
			{
				vit->info().id = count;
				vit->info().inside = true;
				count++;
			}
		}
		return power_shape;
	}

	Delaunay compute_restricted_voronoi_diagram(SURFACE& surface, POINT& medial_points)
	{
		Cgal_Surface_mesh csm;
		load_model_in_cgal(surface, csm);
		auto inner_position = get_attribute<Vec3, PointVertex>(medial_points, "position");
		auto sample_position = get_attribute<Vec3, Vertex>(surface, "position");
		auto sphere_radius = get_attribute<double, PointVertex>(medial_points, "sphere_radius");
		Delaunay constrained_voronoi_diagram;
		uint32 count = 0;
		foreach_cell(medial_points, [&](PointVertex nv){
			Vec3 pos = value<Vec3>(medial_points, inner_position, nv);
			Delaunay_Vertex_handle v = constrained_voronoi_diagram.insert(Point(pos[0], pos[1], pos[2]));
			v->info().inside = true;
			v->info().id = count;
			count++;
			return true;
		});

		foreach_cell(surface, [&](Vertex v) {
			Vec3 pos = value<Vec3>(surface, sample_position, v);
			Delaunay_Vertex_handle vhd = constrained_voronoi_diagram.insert(Point(pos[0], pos[1], pos[2]));
			vhd->info().inside = false;
			vhd->info().id = count;
			count++;
			return true;
		});

		return constrained_voronoi_diagram;
	}
	
	void construct_complete_constrained_voronoi_diagram(Delaunay& delaunay, std::string& name)
	{
		NONMANIFOLD* mv =
			nonmanifold_provider_->add_mesh(std::to_string(nonmanifold_provider_->number_of_meshes()) + "_" + name);
		cgogn::io::IncidenceGraphImportData Inner_RVD_data;
		std::unordered_map<std::pair<uint32, uint32>, uint32, edge_hash, edge_equal> edge_indices;
		uint32 edge_count = 0;
		std::vector<Point> RVD_point;
		for (auto vit = delaunay.finite_vertices_begin(); vit != delaunay.finite_vertices_end(); ++vit)
		{
			RVD_point.push_back(vit->point());
			Inner_RVD_data.vertex_position_.push_back({vit->point().x(), vit->point().y(), vit->point().z()});
		}

		for (auto eit = delaunay.finite_edges_begin(); eit != delaunay.finite_edges_end(); ++eit)
		{
			Delaunay_Vertex_handle v1 = eit->first->vertex(eit->second);
			Delaunay_Vertex_handle v2 = eit->first->vertex(eit->third);
			uint32 v1_ind = v1->info().id;
			uint32 v2_ind = v2->info().id;
			edge_indices.insert({{v1_ind, v2_ind}, edge_count});
			Inner_RVD_data.edges_vertex_indices_.push_back(v1_ind);
			Inner_RVD_data.edges_vertex_indices_.push_back(v2_ind);
			edge_count++;
		}

		for (auto fit = delaunay.finite_facets_begin(); fit != delaunay.finite_facets_end(); ++fit)
		{
			Delaunay_Vertex_handle v[3];
			int count = 0;
			// If face is inside
			for (size_t idx = 0; idx < 4; ++idx)
			{
				if (idx != fit->second)
				{
					v[count] = fit->first->vertex(idx);
					count++;
				}
			}
			Inner_RVD_data.faces_nb_edges_.push_back(3);
			Inner_RVD_data.faces_edge_indices_.push_back(edge_indices[{v[0]->info().id, v[1]->info().id}]);
			Inner_RVD_data.faces_edge_indices_.push_back(edge_indices[{v[1]->info().id, v[2]->info().id}]);
			Inner_RVD_data.faces_edge_indices_.push_back(edge_indices[{v[2]->info().id, v[0]->info().id}]);
		}
		uint32 inner_power_nb_vertices = Inner_RVD_data.vertex_position_.size();
		uint32 inner_power_nb_edges = Inner_RVD_data.edges_vertex_indices_.size() / 2;
		uint32 inner_power_nb_faces = Inner_RVD_data.faces_nb_edges_.size();

		Inner_RVD_data.reserve(inner_power_nb_vertices, inner_power_nb_edges, inner_power_nb_faces);

		import_incidence_graph_data(*mv, Inner_RVD_data);

		std::shared_ptr<NonManifoldAttribute<Vec3>> mv_vertex_position =
			get_attribute<Vec3, NonManifoldVertex>(*mv, "position");
		if (mv_vertex_position)
			nonmanifold_provider_->set_mesh_bb_vertex_position(*mv, mv_vertex_position);

		nonmanifold_provider_->emit_connectivity_changed(*mv);
	}
	
	NONMANIFOLD& construct_inner_restricted_voronoi_diagram(Delaunay& RVD, std::string name)
	{
		
		NONMANIFOLD* mv =
			nonmanifold_provider_->add_mesh(std::to_string(nonmanifold_provider_->number_of_meshes()) + "_" + name);
		cgogn::io::IncidenceGraphImportData Inner_RVD_data;
		std::unordered_map<std::pair<uint32, uint32>, uint32, edge_hash, edge_equal> edge_indices;
		uint32 edge_count = 0;
		std::vector<Point> RVD_point;
		for (auto vit = RVD.finite_vertices_begin(); vit != RVD.finite_vertices_end(); ++vit)
		{
			// if the point is inside
			if (vit->info().inside)
			{
				RVD_point.push_back(vit->point());
				Inner_RVD_data.vertex_position_.push_back({vit->point().x(), vit->point().y(), vit->point().z()});
			}
		}

		for (auto eit = RVD.finite_edges_begin(); eit != RVD.finite_edges_end(); ++eit)
		{
			Delaunay_Vertex_handle v1 = eit->first->vertex(eit->second);
			Delaunay_Vertex_handle v2 = eit->first->vertex(eit->third);
			if (v1->info().inside && v2->info().inside)
			{
				// Add edge
				uint32 v1_ind = v1->info().id;
				uint32 v2_ind = v2->info().id;
				edge_indices.insert({{v1_ind, v2_ind}, edge_count});
				Inner_RVD_data.edges_vertex_indices_.push_back(v1_ind);
				Inner_RVD_data.edges_vertex_indices_.push_back(v2_ind);
				edge_count++;
			}
		}

		for (auto fit = RVD.finite_facets_begin(); fit != RVD.finite_facets_end(); ++fit)
		{
			Delaunay_Vertex_handle v[3];
			int count = 0;
			bool inside = true;
			// If face is inside
			for (size_t idx = 0; idx < 4; ++idx)
			{
				if (idx != fit->second)
				{
					v[count] = fit->first->vertex(idx);
					inside &= v[count]->info().inside;
					count++;
				}
			}
			if (inside)
			{
				Inner_RVD_data.faces_nb_edges_.push_back(3);
				Inner_RVD_data.faces_edge_indices_.push_back(edge_indices[{v[0]->info().id, v[1]->info().id}]);
				Inner_RVD_data.faces_edge_indices_.push_back(edge_indices[{v[1]->info().id, v[2]->info().id}]);
				Inner_RVD_data.faces_edge_indices_.push_back(edge_indices[{v[2]->info().id, v[0]->info().id}]);
			}
		}
		uint32 inner_power_nb_vertices = Inner_RVD_data.vertex_position_.size();
		uint32 inner_power_nb_edges = Inner_RVD_data.edges_vertex_indices_.size() / 2;
		uint32 inner_power_nb_faces = Inner_RVD_data.faces_nb_edges_.size();

		Inner_RVD_data.reserve(inner_power_nb_vertices, inner_power_nb_edges, inner_power_nb_faces);

		import_incidence_graph_data(*mv, Inner_RVD_data);

		std::shared_ptr<NonManifoldAttribute<Vec3>> mv_vertex_position =
			get_attribute<Vec3, NonManifoldVertex>(*mv, "position");
		if (mv_vertex_position)
			nonmanifold_provider_->set_mesh_bb_vertex_position(*mv, mv_vertex_position);

		nonmanifold_provider_->emit_connectivity_changed(*mv);
		return *mv;
	}
	
	NONMANIFOLD& constrcut_power_shape_non_manifold(Regular& power_shape, string name)
	{
		NONMANIFOLD* mv =
			nonmanifold_provider_->add_mesh(std::to_string(nonmanifold_provider_->number_of_meshes()) + "_" +name);
		cgogn::io::IncidenceGraphImportData Inner_Power_shape_data;
		std::unordered_map<std::pair<uint32, uint32>, uint32, edge_hash, edge_equal> edge_indices;
		uint32 edge_count = 0;
		std::vector<Weight_Point> power_point;
		for (auto vit = power_shape.finite_vertices_begin(); vit != power_shape.finite_vertices_end(); ++vit)
		{
			// if the point is inside
			
			if (vit->info().inside)
			{
				power_point.push_back(vit->point());
				Inner_Power_shape_data.vertex_position_.push_back(
					{vit->point().x(), vit->point().y(), vit->point().z()});
			}
		}
	
		for (auto eit = power_shape.finite_edges_begin(); eit != power_shape.finite_edges_end(); ++eit)
		{
			Regular::Vertex_handle v1 = eit->first->vertex(eit->second);
			Regular::Vertex_handle v2 = eit->first->vertex(eit->third);
			if (v1->info().inside && v2->info().inside)
			{
				// Add edge
				uint32 v1_ind = v1->info().id;
				uint32 v2_ind = v2->info().id;
				edge_indices.insert({{v1_ind, v2_ind}, edge_count});
				Inner_Power_shape_data.edges_vertex_indices_.push_back(v1_ind);
				Inner_Power_shape_data.edges_vertex_indices_.push_back(v2_ind);
				edge_count++;
			}
		}

		for (auto fit = power_shape.finite_facets_begin(); fit != power_shape.finite_facets_end(); ++fit)
		{
			Regular::Vertex_handle v[3];
			int count = 0;
			bool inside = true;
			// If face is inside
			for (size_t idx = 0; idx < 4; ++idx)
			{
				if (idx != fit->second)
				{
					v[count]= fit->first->vertex(idx);
					inside &= v[count]->info().inside;
					count++;
				}
			}
			if (inside)
			{
				Inner_Power_shape_data.faces_nb_edges_.push_back(3);
				Inner_Power_shape_data.faces_edge_indices_.push_back(edge_indices[{v[0]->info().id, v[1]->info().id}]);
				Inner_Power_shape_data.faces_edge_indices_.push_back(edge_indices[{v[1]->info().id, v[2]->info().id}]);
				Inner_Power_shape_data.faces_edge_indices_.push_back(edge_indices[{v[2]->info().id, v[0]->info().id}]);
				
			}
		}
		uint32 inner_power_nb_vertices = Inner_Power_shape_data.vertex_position_.size();
		uint32 inner_power_nb_edges = Inner_Power_shape_data.edges_vertex_indices_.size() / 2;
		uint32 inner_power_nb_faces = Inner_Power_shape_data.faces_nb_edges_.size();

		Inner_Power_shape_data.reserve(inner_power_nb_vertices, inner_power_nb_edges, inner_power_nb_faces);

		import_incidence_graph_data(*mv, Inner_Power_shape_data);
		auto sphere_raidus = add_attribute<double, NonManifoldVertex>(*mv, "sphere_radius");
		auto sphere_info = add_attribute<Vec4, NonManifoldVertex>(*mv, "sphere_info");
		for (uint32 idx = 0u; idx < power_point.size(); ++idx)
		{
			(*sphere_raidus)[idx] = std::sqrt(power_point[idx].weight());
			(*sphere_info)[idx] =
				Vec4(power_point[idx].x(), power_point[idx].y(), power_point[idx].z(), power_point[idx].weight());
			
		}

		std::shared_ptr<NonManifoldAttribute<Vec3>> mv_vertex_position =
			get_attribute<Vec3, NonManifoldVertex>(*mv, "position");
		if (mv_vertex_position)
			nonmanifold_provider_->set_mesh_bb_vertex_position(*mv, mv_vertex_position);

		nonmanifold_provider_->emit_connectivity_changed(*mv);
		return *mv;
	}
	
	void compute_stability_ratio_edge(NONMANIFOLD& nm, NonManifoldEdge e,
									  std::shared_ptr<NonManifoldAttribute<double>>& stability_ratio,
									  std::shared_ptr<NonManifoldAttribute<Vec3>>& stability_color,
									  std::shared_ptr<NonManifoldAttribute<Vec4>>& sphere_info)
	{
		auto iv = incident_vertices(nm, e);
		NonManifoldVertex v1 = iv[0];
		NonManifoldVertex v2 = iv[1];
		const Vec3& v1_p = value<Vec4>(nm, sphere_info, v1).head<3>();
		const Vec3& v2_p = value<Vec4>(nm, sphere_info, v2).head<3>();
		const double& r1 = value<Vec4>(nm, sphere_info, v1).w();
		const double& r2 = value<Vec4>(nm, sphere_info, v2).w();
		const double center_dist = (v1_p - v2_p).norm();
		double dis = std::max(0.0, (center_dist - std::abs(r1 - r2)));
		if (center_dist == 0.0)
		{
			(*stability_ratio)[e.index_] = 0.0;
			(*stability_color)[e.index_] = Vec3(0, 0, 0.5);
			return;
		}
		double stability = dis / center_dist;
		value<double>(nm, stability_ratio,e) = stability;
		value<Vec3>(nm, stability_color, e) =
			(stability <= 0.5) ? Vec3(0, stability, (0.5 - stability)) : Vec3(stability - 0.5, (1 - stability), 0);
	}

	void compute_stability_ratio(NONMANIFOLD& nm)
	{
		auto stability_ratio = add_attribute<double, NonManifoldEdge>(nm, "stability_ratio");
		auto stability_color = add_attribute<Vec3, NonManifoldEdge>(nm, "stability_color");
		auto sphere_info =
			get_attribute<Vec4, NonManifoldVertex>(nm, "sphere_info"); // {center, radius} = {x, y, z, r}
		parallel_foreach_cell(nm, [&](NonManifoldEdge e) -> bool {
			compute_stability_ratio_edge(nm, e, stability_ratio, stability_color, sphere_info);
			return true;
		});
	}

	void collapse_non_manifold_using_QMat(NONMANIFOLD& nm, uint32 number_vertices_remain, float k)
	{
		using QMatHelper = modeling::DecimationSQEM_Helper<NONMANIFOLD>;
		using Slab_Quadric = geometry::Slab_Quadric;

		auto sphere_radius = get_attribute<double, NonManifoldVertex>(nm, "sphere_radius");
		auto stability_ratio = get_attribute<double, NonManifoldEdge>(nm, "stability_ratio");
		auto stability_color = get_attribute<Vec3, NonManifoldEdge>(nm, "stability_color");
		auto position = get_attribute<Vec3, NonManifoldVertex>(nm, "position");
		auto sphere_info = get_attribute<Vec4, NonManifoldVertex>(nm, "sphere_info");
		auto fixed_vertex = add_attribute<bool, NonManifoldVertex>(nm, "fixed_vertex");
		foreach_cell(nm, [&](NonManifoldVertex v) -> bool {
			value<bool>(nm, fixed_vertex, v) = false;
			return true;
		});

		QMatHelper helper(k, nm, position, sphere_info, stability_color, stability_ratio, sphere_radius, fixed_vertex);

		helper.initial_slab_mesh();
		helper.initial_boundary_mesh();
		helper.initial_collapse_queue();
		helper.simplify(number_vertices_remain, true);

		remove_attribute<NonManifoldVertex>(nm, fixed_vertex);
		nonmanifold_provider_->emit_connectivity_changed(nm);
		nonmanifold_provider_->emit_attribute_changed(nm, position.get());
		nonmanifold_provider_->emit_attribute_changed(nm, sphere_radius.get());
		nonmanifold_provider_->emit_attribute_changed(nm, stability_ratio.get());
		nonmanifold_provider_->emit_attribute_changed(nm, stability_color.get());
	}

	
	 void coverage_axis_RVD(SURFACE& surface, POINT& mv, HighsSolution& solution)
	{
		Cgal_Surface_mesh csm;
		load_model_in_cgal(surface, csm);
		auto inner_position = get_attribute<Vec3, PointVertex>(mv, "position");
		auto sample_position = get_attribute<Vec3, SurfaceVertex>(surface, "position");
		auto sphere_radius = get_attribute<double, PointVertex>(mv, "sphere_radius");
		Delaunay constrained_voronoi_diagram;
		uint32 count = 0;
		foreach_cell(mv, [&](PointVertex nv) {
			if (solution.col_value[index_of(mv, nv)] > 1e-5)
			{
				Vec3 pos = value<Vec3>(mv, inner_position, nv);
				Delaunay_Vertex_handle v = constrained_voronoi_diagram.insert(Point(pos[0], pos[1], pos[2]));
				v->info().inside = true;
				v->info().id = count;
				count++;
			}
			return true;
		});

		foreach_cell(surface, [&](SurfaceVertex v) {
			Vec3 pos = value<Vec3>(surface, sample_position, v);
			Delaunay_Vertex_handle vhd = constrained_voronoi_diagram.insert(Point(pos[0], pos[1], pos[2]));
			vhd->info().inside = false;
			vhd->info().id = count;
			count++;
			return true;
		});
		construct_inner_restricted_voronoi_diagram(constrained_voronoi_diagram,
													"inner_RVD_" + surface_provider_->mesh_name(surface));
		construct_complete_constrained_voronoi_diagram(constrained_voronoi_diagram,
													   "complete_RVD_" + surface_provider_->mesh_name(surface));
	}

	void coverage_axis_PD(SURFACE& surface, POINT& selected_points, HighsSolution& solution, double dilation_factor)
	{
		auto inner_position = get_attribute<Vec3, PointVertex>(selected_points, "position");
		auto sphere_radius = get_attribute<double, PointVertex>(selected_points, "sphere_radius");
		auto dilated_radius = get_attribute<double, PointVertex>(selected_points, "dilated_radius");
		auto oustide_pole_center_att = get_attribute<Vec3, PointVertex>(selected_points, "outside_pole_center");
		auto oustide_pole_radius_att = get_attribute<double, PointVertex>(selected_points, "outside_pole_radius");

		auto coverage_infos = get_attribute<std::vector<uint32>, Vertex>(surface, "coverage_infos");
		auto sample_position = get_attribute<Vec3, Vertex>(surface, "position");
	
		POINT* outside_poles_point = point_provider_->add_mesh(point_provider_->mesh_name(surface) + "filtered_outside_poles");

		cgogn::io::PointImportData outside_poles;

		std::vector<double> outside_poles_radius;
		std::vector<double> selected_points_radius(nb_cells<PointVertex>(selected_points), 0.0);
		Regular power_shape;
		foreach_cell(selected_points, [&](PointVertex nv) {
			if (solution.col_value[index_of(selected_points, nv)] > 1e-5)
			{
				Vec3 pos = value<Vec3>(selected_points, inner_position, nv);

				selected_points_radius[index_of(selected_points, nv)] =
					value<double>(selected_points, dilated_radius, nv) -
					value<double>(selected_points, sphere_radius, nv);
				power_shape.insert(
					Weight_Point(Point(pos[0], pos[1], pos[2]), value<double>(selected_points, dilated_radius, nv) *
													value<double>(selected_points, dilated_radius, nv)));

				Vec3 outside_pole_pos = value<Vec3>(selected_points, oustide_pole_center_att, nv);
				outside_poles.vertex_position_.emplace_back(outside_pole_pos[0], outside_pole_pos[1],
															outside_pole_pos[2]);
				outside_poles_radius.push_back(value<double>(selected_points, oustide_pole_radius_att, nv));
			}
			return true;
		});
		outside_poles.reserve(outside_poles_radius.size());

		cgogn::io::import_point_data(*outside_poles_point, outside_poles);
		auto outside_poles_position = get_attribute<Vec3, PointVertex>(*outside_poles_point, "position");
		if (outside_poles_position)
		{
			point_provider_->set_mesh_bb_vertex_position(*outside_poles_point, outside_poles_position);
		}
		auto sphere_radius_att = add_attribute<double, PointVertex>(*outside_poles_point, "sphere_radius");
		for (size_t idx = 0; idx < outside_poles_radius.size(); idx++)
		{
			(*sphere_radius_att)[idx] = outside_poles_radius[idx];
		}

		foreach_cell(surface, [&](SurfaceVertex v) {
			Vec3 pos = value<Vec3>(surface, sample_position, v);
			std::vector<uint32> covered_sphere = value<std::vector<uint32>>(surface, coverage_infos, v);
			double sample_sphere_radius = std::numeric_limits<double>::min();
			for (uint32 idx : covered_sphere)
			{
				if (solution.col_value[idx] > 1e-5)
				{
					sample_sphere_radius = std::max(sample_sphere_radius, selected_points_radius[idx]);
				}
			}
			power_shape.insert(
				Weight_Point(Point(pos[0], pos[1], pos[2]), sample_sphere_radius * sample_sphere_radius));
			return true;
		});
		int count = 0;
	
		for (auto vit = power_shape.finite_vertices_begin(); vit != power_shape.finite_vertices_end(); ++vit)
		{
			// if the point is inside
			Point p = vit->point().point();
			if (pointInside(*tree, p))
			{
				vit->info().id = count;
				vit->info().inside = true;
				count++;
			}
		}
		
		constrcut_power_shape_non_manifold(power_shape, "_coverage_axis" + surface_provider_->mesh_name(surface));
	}

	void coverage_axis_collapse(SURFACE& surface, POINT& medial_points, HighsSolution& solution)
	{
		Delaunay RVD = compute_restricted_voronoi_diagram(surface, medial_points);
		NONMANIFOLD& nm = construct_inner_restricted_voronoi_diagram(RVD, "constrained_voronoi_diagram");

		using QMatHelper = modeling::DecimationSQEM_Helper<NONMANIFOLD>;
		using Slab_Quadric = geometry::Slab_Quadric;

		compute_stability_ratio(nm);
		auto sphere_radius = get_attribute<double, NonManifoldVertex>(nm, "sphere_radius");
		auto stability_ratio = get_attribute<double, NonManifoldEdge>(nm, "stability_ratio");
		auto stability_color = get_attribute<Vec3, NonManifoldEdge>(nm, "stability_color");
		auto position = get_attribute<Vec3, NonManifoldVertex>(nm, "position");
		auto sphere_info = get_attribute<Vec4, NonManifoldVertex>(nm, "sphere_info");
		auto fixed_vertex = add_attribute<bool, NonManifoldVertex>(nm, "fixed_vertex");
		int number_vertex_remain = 0;
		foreach_cell(nm, [&](NonManifoldVertex nv) -> bool {
			Vec3 pos = value<Vec3>(nm, position, nv);
			if (solution.col_value[index_of(nm, nv)] > 1e-5)
			{
				value<bool>(nm, fixed_vertex, nv) = true;
				number_vertex_remain++;
			}
			else
				value<bool>(nm, fixed_vertex, nv) = false;

			return true;
		});
		QMatHelper helper(0.00001, nm, position, sphere_info, stability_color, stability_ratio, sphere_radius,
						  fixed_vertex);

		helper.initial_slab_mesh();
		helper.initial_boundary_mesh();
		helper.initial_collapse_queue();
		helper.simplify(number_vertex_remain, false);

		foreach_cell(nm, [&](NonManifoldVertex nv) {
			Vec3 pos = value<Vec3>(nm, position, nv);
			return true;
		});

		remove_attribute<NonManifoldVertex>(nm, fixed_vertex);
		nonmanifold_provider_->emit_connectivity_changed(nm);
		nonmanifold_provider_->emit_attribute_changed(nm, position.get());
		nonmanifold_provider_->emit_attribute_changed(nm, sphere_radius.get());
		nonmanifold_provider_->emit_attribute_changed(nm, stability_ratio.get());
		nonmanifold_provider_->emit_attribute_changed(nm, stability_color.get());
	}

	HighsSolution point_selection_by_coverage_axis(SURFACE& surface, POINT& candidates, double dilation_factor)
	{
		typedef Eigen::SparseMatrix<double> SpMat;
		typedef Eigen::Triplet<double> T;
		std::vector<T> triplets;

		auto dilated_radius = add_attribute<double, PointVertex>(candidates, "dilated_radius");
		auto coverage_infos = add_attribute<std::vector<uint32>, SurfaceVertex>(surface, "coverage_infos");

		auto inner_position = get_attribute<Vec3, PointVertex>(candidates, "position");
		auto sphere_radius = get_attribute<double, PointVertex>(candidates, "sphere_radius");
		auto outside_pole_radius = get_attribute<double, PointVertex>(candidates, "outside_pole_radius");
		auto sample_position = get_attribute<Vec3, SurfaceVertex>(surface, "position");
		auto inner_point_nb = nb_cells<PointVertex>(candidates);
		auto sample_point_nb = nb_cells<SurfaceVertex>(surface);

		std::vector<Vec3> sample_points(sample_point_nb, Vec3(0, 0, 0));
		std::vector<Vec3> inner_points(inner_point_nb, Vec3(0, 0, 0));
		std::vector<double> weights(inner_point_nb, 0);
		parallel_foreach_cell(candidates, [&](PointVertex nv) {
			uint32 index = index_of(candidates, nv);
			double radius = value<double>(candidates, sphere_radius, nv);
			double outside_radius = value<double>(candidates, outside_pole_radius, nv);
			if (radius * dilation_factor > outside_radius)
			{
				weights[index] = radius + 0.8 * outside_radius;
			}
			else
			{
				weights[index] = radius * (1 + dilation_factor);
			}
			value<double>(candidates, dilated_radius, nv) = weights[index];
			inner_points[index] = value<Vec3>(candidates, inner_position, nv);

			return true;
		});
		foreach_cell(surface, [&](SurfaceVertex v) {
			for (size_t idx = 0; idx < inner_point_nb; idx++)
			{
				if (inside_sphere(value<Vec3>(surface, sample_position, v), inner_points[idx], weights[idx]))
				{
					triplets.push_back(T(index_of(surface, v), idx, 1.0));
					value<std::vector<uint32>>(surface, coverage_infos, v).push_back(idx);
				}
			}
			return true;
		});
		std::cout << triplets.size() << std::endl;

		SpMat A(sample_point_nb, inner_point_nb);
		A.setFromTriplets(triplets.begin(), triplets.end());
		A.makeCompressed();
		HighsModel model;
		model.lp_.a_matrix_.format_ = MatrixFormat::kColwise;
		model.lp_.num_col_ = A.cols();
		model.lp_.num_row_ = A.rows();
		model.lp_.sense_ = ObjSense::kMinimize;
		// Adding decision variable bounds
		HighsVarType type = HighsVarType::kInteger;
		model.lp_.col_cost_ = std::vector<double>(model.lp_.num_col_, 1.0);
		model.lp_.col_lower_ = std::vector<double>(model.lp_.num_col_, 0.0);
		model.lp_.col_upper_ = std::vector<double>(model.lp_.num_col_, 1.0);
		model.lp_.row_lower_ = std::vector<double>(model.lp_.num_row_, 1.0);
		model.lp_.row_upper_ = std::vector<double>(model.lp_.num_row_, 1e30);
		model.lp_.integrality_ = std::vector<HighsVarType>(model.lp_.num_col_, type);

		model.lp_.a_matrix_.num_col_ = model.lp_.num_col_;
		model.lp_.a_matrix_.num_row_ = model.lp_.num_row_;
		model.lp_.a_matrix_.start_.resize(A.cols() + 1);
		model.lp_.a_matrix_.index_.resize(A.nonZeros());
		model.lp_.a_matrix_.value_.resize(A.nonZeros());

		// Copy the data from A to the vectors
		std::copy(A.outerIndexPtr(), A.outerIndexPtr() + A.cols() + 1, model.lp_.a_matrix_.start_.begin());
		std::copy(A.innerIndexPtr(), A.innerIndexPtr() + A.nonZeros(), model.lp_.a_matrix_.index_.begin());
		std::copy(A.valuePtr(), A.valuePtr() + A.nonZeros(), model.lp_.a_matrix_.value_.begin());

		Highs highs;
		HighsStatus status = highs.passModel(model);
		HighsSolution solution;
		highs.setOptionValue("time_limit", 200);
		if (status == HighsStatus::kOk)
		{
			highs.run();

			assert(status == HighsStatus::kOk);
			solution = highs.getSolution();
		}
		return solution;
	}
	void compute_power_shape(SURFACE& surface)
	{
		Cgal_Surface_mesh csm;
		load_model_in_cgal(*selected_surface_mesh_, csm);
		Tree tree(faces(csm).first, faces(csm).second, csm);
		tree.accelerate_distance_queries();
		Regular reg = compute_regular_tredrahedron(tree, tri_);
		constrcut_power_shape_non_manifold(reg, surface_provider_->mesh_name(surface) + "_inner_power_shape");
	}

	void compute_original_power_diagram(SURFACE& surface)
	{
		compute_initial_non_manifold(tri_, surface_provider_->mesh_name(surface) + "_inner_voronoi_diagram");
	}

	void construct_surface_sample(SURFACE& surface)
	{
		auto position = get_or_add_attribute<Vec3,SurfaceVertex>(surface, "position");
		auto bvh_vertex_index = get_or_add_attribute<uint32, SurfaceVertex>(surface, "bvh_vertex_index");
		uint32 id = 0;
		load_model_in_cgal(surface, csm);
		std::vector<Point> mesh_samples;
		CGAL::Polygon_mesh_processing::sample_triangle_mesh(
			csm, std::back_inserter(mesh_samples),
			CGAL::parameters::use_grid_sampling(true).grid_spacing(0.01));
		foreach_cell(surface, [&](SurfaceVertex v) {
			value<uint32>(surface, bvh_vertex_index, v) = id++;
			surface_kdt_vertices.push_back(v);
			surface_vertex_position_vector.push_back(value<Vec3>(surface, position, v));
			return true;
		});
		uint32 nb_faces = nb_cells<SurfaceFace>(surface);
		std::vector<SurfaceFace> bvh_faces;
		bvh_faces.reserve(nb_faces);
		std::vector<uint32> face_vertex_indices;
		face_vertex_indices.reserve(nb_faces * 3);

		foreach_cell(surface, [&](SurfaceFace f) -> bool {
			bvh_faces.push_back(f);
			foreach_incident_vertex(surface, f, [&](SurfaceVertex v) -> bool {
				face_vertex_indices.push_back(value<uint32>(surface, bvh_vertex_index, v));
				return true;
			});
			return true;
		});

		surface_bvh = std::make_unique<acc::BVHTree<uint32, Vec3>>(face_vertex_indices, surface_vertex_position_vector);
		POINT* surface_sample = point_provider_->add_mesh(point_provider_->mesh_name(surface) + "surface_sample");
		auto sample_position = get_or_add_attribute<Vec3, PointVertex>(*surface_sample, "position");
		auto sample_normal = get_or_add_attribute<Vec3, PointVertex>(*surface_sample, "normal");
		for (Point& p : mesh_samples)
		{
			Vec3 pos(p.x(), p.y(), p.z());
			std::pair<uint32, Vec3> bvh_res;
			bool found = surface_bvh->closest_point(pos, &bvh_res);
			if (!found)
				std::cout << "not found"<<std::endl;
			SurfaceFace& f = bvh_faces[bvh_res.first];
			PointVertex new_v = add_vertex(*surface_sample);
			value<Vec3>(*surface_sample, sample_position, new_v) = pos;
			value<Vec3>(*surface_sample, sample_normal, new_v) = geometry::normal(surface, f, position.get());
		}
		if (sample_position)
			point_provider_->set_mesh_bb_vertex_position(*surface_sample, sample_position);
		point_provider_->emit_connectivity_changed(*surface_sample);
		point_provider_->emit_attribute_changed(*surface_sample, sample_position.get());
		point_provider_->emit_attribute_changed(*surface_sample, sample_normal.get());

	}
	/*void transfer_file()
	{
		using fs = std::filesystem;
		std::string in_file_dir = "C:/Users/huang/Devel/CGoGN_3/data/meshes/Thingi10K/Thingi10K/rawmesh";
		std::string out_file_dir = "C:/Users/huang/Devel/CGoGN_3/data/meshes/off/Thingi10K";

		for (const auto& entry : fs::directory_iterator(in_file_dir))
		{
			std::string in_file = entry.path().string();
			std::string off_file = out_file_dir + "/" + entry.path().filename().string() + ;
		}

	}*/
	void sample_medial_axis(SURFACE& s)
	{
		auto vertex_position = get_or_add_attribute<Vec3, SurfaceVertex>(s, "position");
		auto vertex_normal = get_or_add_attribute<Vec3, SurfaceVertex>(s, "normal");
		auto medial_axis_samples_position_ =
			get_or_add_attribute<Vec3, SurfaceVertex>(s, "medial_axis_samples_position");
		auto medial_axis_samples_radius_ = get_or_add_attribute<Scalar, SurfaceVertex>(s, "medial_axis_samples_radius");
		auto medial_axis_samples_closest_points_ =
			get_or_add_attribute<std::pair<SurfaceVertex, SurfaceVertex>, SurfaceVertex>(
				s, "medial_axis_samples_closest_points");
		auto medial_axis_samples_angle_ = get_or_add_attribute<Scalar, SurfaceVertex>(s, "medial_axis_samples_angle");
		auto medial_axis_samples_feature_value_ =
			get_or_add_attribute<Scalar, SurfaceVertex>(s, "medial_axis_samples_feature_value");
		auto medial_axis_samples_weight_ = get_or_add_attribute<Scalar, SurfaceVertex>(s, "medial_axis_samples_weight");
		auto kmax = get_or_add_attribute<Scalar, SurfaceVertex>(s, "kmax");
		geometry::shrinking_ball_centers<SURFACE, false>(
			s, vertex_position.get(), vertex_normal.get(),
										 medial_axis_samples_position_.get(),
										 medial_axis_samples_radius_.get(),
										 medial_axis_samples_closest_points_.get());
		MeshData<SURFACE>& md = surface_provider_->mesh_data(s);
		auto filtered_medial_axis_samples_set_ = &md.template get_or_add_cells_set<SurfaceVertex>("filtered_medial_axis_samples");
		filtered_medial_axis_samples_set_->select_if([&](SurfaceVertex v) { return true; });
		foreach_cell(s, [&](SurfaceVertex v) -> bool {
			const Vec3& c = value<Vec3>(s, medial_axis_samples_position_, v);
			const auto& [v1, v2] =
				value<std::pair<SurfaceVertex, SurfaceVertex>>(s, medial_axis_samples_closest_points_, v);
			const auto c1 = value<Vec3>(s, vertex_position, v1);
			const auto c2 = value<Vec3>(s, vertex_position, v2);
			const Scalar r = value<Scalar>(s, medial_axis_samples_radius_, v);
			max_radius_ = std::max(max_radius_, r);
			min_radius_ = std::min(min_radius_, r);
			const Scalar angle = geometry::angle(c1 - c, c2 - c);
			max_angle_ = std::max(max_angle_, angle);
			min_angle_ = std::min(min_angle_, angle);
			value<Scalar>(s, medial_axis_samples_angle_, v) = angle;
			value<Scalar>(s, medial_axis_samples_feature_value_, v) =
				(1000/(1 + std::exp(1.9*angle / M_PI))) * (1/(1 + std::exp(30*r)));
			
			/*Scalar sum = 0;
			uint32_t count = 0;
			for_n_ring<SURFACE,SurfaceVertex>(s, v, 3, [&](SurfaceVertex iv)
			{
				count++;
				sum += value<Scalar>(s, medial_axis_samples_feature_value_, iv);
				return true;
			});*/
			value<Scalar>(s, medial_axis_samples_weight_, v) = value<Scalar>(s, kmax, v);
			/*std::cout << "feature value: "<< value<Scalar>(s, medial_axis_samples_feature_value_, v) 
				<< ", weight: " << value<Scalar>(s, medial_axis_samples_weight_, v) << std::endl;*/

			return true;
		});
		normalise_scalar(s, medial_axis_samples_weight_);
		normalise_scalar(s, medial_axis_samples_feature_value_);
		surface_provider_->emit_attribute_changed(s, medial_axis_samples_position_.get());
		surface_provider_->emit_attribute_changed(s, medial_axis_samples_radius_.get());
		surface_provider_->emit_attribute_changed(s, medial_axis_samples_angle_.get());
		surface_provider_->emit_cells_set_changed(s, filtered_medial_axis_samples_set_);
	}

	void filter_medial_axis_samples(SURFACE& s)
	{
		auto vertex_position = get_attribute<Vec3, SurfaceVertex>(s, "position");
		MeshData<SURFACE>& md = surface_provider_->mesh_data(s);
		auto filtered_medial_axis_samples_set_ = &md.template get_or_add_cells_set<SurfaceVertex>("filtered_medial_axis_samples");
		filtered_medial_axis_samples_set_->clear();
		auto medial_axis_samples_position_ = get_attribute<Vec3, SurfaceVertex>(s, "medial_axis_samples_position");
		auto medial_axis_samples_radius_ = get_attribute<Scalar, SurfaceVertex>(s, "medial_axis_samples_radius");
		auto medial_axis_samples_closest_points_ =
			get_attribute<std::pair<SurfaceVertex, SurfaceVertex>, SurfaceVertex>(s,
																				  "medial_axis_samples_closest_points");
		auto medial_axis_samples_angle_ = get_or_add_attribute<Scalar, SurfaceVertex>(s, "medial_axis_samples_angle");
		foreach_cell(s, [&](SurfaceVertex v) -> bool {
			const Scalar r = value<Scalar>(s, medial_axis_samples_radius_, v);
			const Scalar angle = value<Scalar>(s, medial_axis_samples_angle_, v);
			if (r > radius_threshold_ && angle > angle_threshold_)
				filtered_medial_axis_samples_set_->select(v);
			return true;
		});

		surface_provider_->emit_cells_set_changed(s, filtered_medial_axis_samples_set_);
	}

	void Twin_medial_axis_samples(SURFACE& s)
	{
		auto vertex_position = get_attribute<Vec3, SurfaceVertex>(s, "position");
		MeshData<SURFACE>& md = surface_provider_->mesh_data(s);
		auto twin_medial_axis_samples_set_ =
			&md.template get_or_add_cells_set<SurfaceVertex>("twin_medial_axis_samples_set_");
		twin_medial_axis_samples_set_->clear();
		auto medial_axis_samples_position_ = get_attribute<Vec3, SurfaceVertex>(s, "medial_axis_samples_position");
		auto medial_axis_samples_radius_ = get_attribute<Scalar, SurfaceVertex>(s, "medial_axis_samples_radius");
		auto medial_axis_samples_closest_points_ =
			get_attribute<std::pair<SurfaceVertex, SurfaceVertex>, SurfaceVertex>(s,
																				  "medial_axis_samples_closest_points");
		auto medial_axis_samples_angle_ = get_or_add_attribute<Scalar, SurfaceVertex>(s, "medial_axis_samples_angle");
		foreach_cell(s, [&](SurfaceVertex v) -> bool {
			auto [v1, v2] = value<std::pair<SurfaceVertex, SurfaceVertex>>(s, medial_axis_samples_closest_points_, v);
			auto [v3, v4] = value<std::pair<SurfaceVertex, SurfaceVertex>>(s, medial_axis_samples_closest_points_, v2);
			if (index_of(s, v1) == index_of(s,v4))
				twin_medial_axis_samples_set_->select(v);
			return true;
		});

		surface_provider_->emit_cells_set_changed(s, twin_medial_axis_samples_set_);
	}


	HighsSolution shrinking_balls_coverage_axis(SURFACE& surface, Scalar dilation_factor)
	{
		typedef Eigen::SparseMatrix<Scalar> SpMat;
		typedef Eigen::Triplet<Scalar> T;
		std::unordered_map<uint32, uint32> local_index;
		std::vector<T> triplets;

		auto medial_axis_samples_dilated_radius =
			add_attribute<Scalar, SurfaceVertex>(surface, "medial_axis_samples_dilated_radius");
		auto medial_axis_samples_position = get_attribute<Vec3, SurfaceVertex>(surface, "medial_axis_samples_position");
		auto medial_axis_samples_radius = get_attribute<Scalar, SurfaceVertex>(surface, "medial_axis_samples_radius");
		auto sample_position = get_attribute<Vec3, SurfaceVertex>(surface, "position");
		auto filtered_medial_axis_samples_set =
			&surface_provider_->mesh_data(surface).template get_or_add_cells_set<SurfaceVertex>(
				"filtered_medial_axis_samples");
		auto selected_medial_axis_samples_set =
			&surface_provider_->mesh_data(surface).template get_or_add_cells_set<SurfaceVertex>(
				"selected_medial_axis_samples");
		uint32 nb_candidates = 0;
		uint32 nb_surface_samples = 0;
		std::vector<Vec3> sample_points;
		std::vector<Vec3> inner_points;
		std::vector<Scalar> weights;
		filtered_medial_axis_samples_set->foreach_cell([&](SurfaceVertex sv) 
			{
			value<Scalar>(surface, medial_axis_samples_dilated_radius, sv) = value<Scalar>(surface, medial_axis_samples_radius, sv) /** (1 */+ 0.02;
			sample_points.push_back(value<Vec3>(surface, sample_position, sv));
			inner_points.push_back(value<Vec3>(surface, medial_axis_samples_position, sv));
			weights.push_back(value<Scalar>(surface, medial_axis_samples_dilated_radius, sv));
			local_index.insert({index_of(surface, sv), nb_candidates});
			nb_candidates++;
			nb_surface_samples++;
			return true;
		});
		for (size_t idx1 = 0; idx1 < nb_surface_samples; idx1++)
		{
			for (size_t idx2 = 0; idx2 < nb_candidates; idx2++)
			{
				if (inside_sphere(sample_points[idx1], inner_points[idx2], weights[idx2]))
				{
					triplets.push_back(T(idx1, idx2, 1.0));
				}
			}
		}
		std::cout << triplets.size() << std::endl;

		SpMat A(nb_surface_samples, nb_candidates);
		A.setFromTriplets(triplets.begin(), triplets.end());
		A.makeCompressed();
		HighsModel model;
		model.lp_.a_matrix_.format_ = MatrixFormat::kColwise;
		model.lp_.num_col_ = A.cols();
		model.lp_.num_row_ = A.rows();
		model.lp_.sense_ = ObjSense::kMinimize;
		// Adding decision variable bounds
		HighsVarType type = HighsVarType::kInteger;
		model.lp_.col_cost_ = std::vector<double>(model.lp_.num_col_, 1.0);
		model.lp_.col_lower_ = std::vector<double>(model.lp_.num_col_, 0.0);
		model.lp_.col_upper_ = std::vector<double>(model.lp_.num_col_, 1.0);
		model.lp_.row_lower_ = std::vector<double>(model.lp_.num_row_, 1.0);
		model.lp_.row_upper_ = std::vector<double>(model.lp_.num_row_, 1e30);
		model.lp_.integrality_ = std::vector<HighsVarType>(model.lp_.num_col_, type);

		model.lp_.a_matrix_.num_col_ = model.lp_.num_col_;
		model.lp_.a_matrix_.num_row_ = model.lp_.num_row_;
		model.lp_.a_matrix_.start_.resize(A.cols() + 1);
		model.lp_.a_matrix_.index_.resize(A.nonZeros());
		model.lp_.a_matrix_.value_.resize(A.nonZeros());

		// Copy the data from A to the vectors
		std::copy(A.outerIndexPtr(), A.outerIndexPtr() + A.cols() + 1, model.lp_.a_matrix_.start_.begin());
		std::copy(A.innerIndexPtr(), A.innerIndexPtr() + A.nonZeros(), model.lp_.a_matrix_.index_.begin());
		std::copy(A.valuePtr(), A.valuePtr() + A.nonZeros(), model.lp_.a_matrix_.value_.begin());

		Highs highs;
		HighsStatus status = highs.passModel(model);
		HighsSolution solution;
		highs.setOptionValue("time_limit", 1000);
		highs.setOptionValue("parallel", "on");
		/*highs.setHighsOptionValue("solver", "simplex");*/
		if (status == HighsStatus::kOk)
		{
			highs.run();

			assert(status == HighsStatus::kOk);
			solution = highs.getSolution();
		}
		filtered_medial_axis_samples_set->foreach_cell([&](SurfaceVertex sv) -> bool {
			if (solution.col_value[local_index[index_of(surface, sv)]] > 1e-5)
				selected_medial_axis_samples_set->select(sv);
			return true;
		});
		shrinking_balls_coverage_color(surface);
		return solution;
	}

	void shrinking_balls_coverage_color(SURFACE& surface)
	{
		auto coverage_color = add_attribute<Vec3, SurfaceVertex>(surface, "coverage_color");
		foreach_cell(surface, [&](SurfaceVertex v) {
			value<Vec3>(surface, coverage_color, v) = Vec3(0, 0, 0);
			return true;
		});
		auto medial_axis_samples_dilated_radius =
			get_attribute<Scalar, SurfaceVertex>(surface, "medial_axis_samples_dilated_radius");
		auto medial_axis_samples_position = get_attribute<Vec3, SurfaceVertex>(surface, "medial_axis_samples_position");
		auto medial_axis_samples_radius = get_attribute<Scalar, SurfaceVertex>(surface, "medial_axis_samples_radius");
		auto sample_position = get_attribute<Vec3, SurfaceVertex>(surface, "position");
		auto selected_medial_axis_samples_set =
			&surface_provider_->mesh_data(surface).template get_or_add_cells_set<SurfaceVertex>(
				"selected_medial_axis_samples");
		auto filter_medial_axis_samples =
			&surface_provider_->mesh_data(surface).template get_or_add_cells_set<SurfaceVertex>(
				"filtered_medial_axis_samples");
		CellMarker<Surface, SurfaceVertex> marker(surface);
		selected_medial_axis_samples_set->foreach_cell([&](SurfaceVertex svc) {
			Vec3 color = Vec3(rand()/double(RAND_MAX), rand()/double(RAND_MAX), rand()/double(RAND_MAX));
			foreach_cell(surface, [&](SurfaceVertex svs) {
				if (marker.is_marked(svs))
					return true;
				if (inside_sphere(value<Vec3>(surface, sample_position, svs),
													  value<Vec3>(surface, medial_axis_samples_position, svc),
								  value<Scalar>(surface, medial_axis_samples_dilated_radius, svc)))
				{
					value<Vec3>(surface, coverage_color, svs) = color;
					marker.mark(svs);
				}
				return true;
				});
			return true;
		});

	}

	void shrinking_balls_coverage_axis_PD(SURFACE& surface)
	{
		load_model_in_cgal(surface, csm);
		Tree tree(faces(csm).first, faces(csm).second, csm);
		tree.accelerate_distance_queries();
		auto medial_axis_samples_dilated_radius =
			get_attribute<Scalar, SurfaceVertex>(surface, "medial_axis_samples_dilated_radius");
		auto medial_axis_samples_position = get_attribute<Vec3, SurfaceVertex>(surface, "medial_axis_samples_position");
		auto medial_axis_samples_radius = get_attribute<Scalar, SurfaceVertex>(surface, "medial_axis_samples_radius");
		auto sample_position = get_attribute<Vec3, SurfaceVertex>(surface, "position");
		auto selected_medial_axis_samples_set =
			&surface_provider_->mesh_data(surface).template get_or_add_cells_set<SurfaceVertex>(
				"selected_medial_axis_samples");
		Regular power_shape;
		std::cout << "collect selected points info" << std::endl;
		selected_medial_axis_samples_set->foreach_cell([&](SurfaceVertex sv) {
			Vec3 pos = value<Vec3>(surface, medial_axis_samples_position, sv);
			power_shape.insert(
				Weight_Point(Point(pos[0], pos[1], pos[2]),
											(value<double>(surface, medial_axis_samples_radius, sv)+0.02) *
												(value<double>(surface, medial_axis_samples_radius, sv)+0.02)));

			return true;
		});
		std::cout << "collect surface samples info " << std::endl;
		foreach_cell(surface, [&](SurfaceVertex v) {
			Vec3 pos = value<Vec3>(surface, sample_position, v);
			power_shape.insert(
				Weight_Point(Point(pos[0], pos[1], pos[2]), 0.02*0.02));
			return true;
		});
		std::cout << "determine if the point is outside" << std::endl;
		int count = 0;
		for (auto vit = power_shape.finite_vertices_begin(); vit != power_shape.finite_vertices_end(); ++vit)
		{
			// if the point is inside
			Point p = vit->point().point();
			if (pointInside(tree, p))
			{
				vit->info().id = count;
				vit->info().inside = true;
				count++;
			}
		}
		std::cout << "construct non manifold PD" << std::endl;
		constrcut_power_shape_non_manifold(power_shape, "_coverage_axis" + surface_provider_->mesh_name(surface));
	}

	void compute_samples_inside_distance(SURFACE& surface)
	{
		std::cout << "compute distance" << std::endl;
		std::vector<Vec3> sample_point_position(nb_cells<SurfaceVertex>(surface), Vec3(0, 0, 0));
		std::vector<Vec3> medial_point_position(nb_cells<SurfaceVertex>(surface), Vec3(0, 0, 0));
		std::vector<Scalar> medial_point_radius(nb_cells<SurfaceVertex>(surface), 0);
		std::cout << "get necessary info" << std::endl;
		MeshData<SURFACE>& md = surface_provider_->mesh_data(surface);
		auto filtered_medial_axis_samples =
			&md.template get_or_add_cells_set<SurfaceVertex>("filtered_medial_axis_samples");
		auto sample_position = get_attribute<Vec3, SurfaceVertex>(surface, "position");
		auto medial_axis_samples_position_ =
			get_attribute<Vec3, SurfaceVertex>(surface, "medial_axis_samples_position");
		auto medial_axis_samples_radius_ = get_attribute<Scalar, SurfaceVertex>(surface, "medial_axis_samples_radius");
		std::cout << "generate necessary vector" << std::endl;
		dis_matrix.resize(nb_cells<SurfaceVertex>(surface), nb_cells<SurfaceVertex>(surface));
		dis_matrix.setZero();
		parallel_foreach_cell(surface, [&](SurfaceVertex sv) {
			sample_point_position[index_of(surface, sv)] = value<Vec3>(surface, sample_position, sv);
			return true;
			});
		filtered_medial_axis_samples->foreach_cell( [&](SurfaceVertex sv) {
			medial_point_position[index_of(surface, sv)] = value<Vec3>(surface, medial_axis_samples_position_, sv);
			medial_point_radius[index_of(surface, sv)] = value<Scalar>(surface, medial_axis_samples_radius_, sv);
			return true;
		});
		std::cout << "loop begin" << std::endl;
		for (size_t idx1 = 0; idx1 < sample_point_position.size(); idx1++)
		{
			for (size_t idx2 = 0; idx2 < medial_point_position.size(); idx2++)
			{
				Vec3 sample_point_minus_medial_point = sample_point_position[idx1] - medial_point_position[idx2];
				dis_matrix(idx1, idx2) = std::sqrt(sample_point_minus_medial_point.dot(sample_point_minus_medial_point))-medial_point_radius[idx2];
			}
		}
		std::cout << "loop end" << std::endl;
	}

	struct Cluster_Info
	{
		std::vector<SurfaceVertex> cluster_vertices;
		Scalar cluster_variance;
		
	};
	
	using SphereQueue = std::multimap<Scalar, SurfaceVertex,std::greater<Scalar>>;
	using SphereQueueIt = typename SphereQueue::const_iterator;
	using SphereInfo = std::pair<bool, SphereQueueIt>;

	void iniialise_cluster(SURFACE& surface)
	{
		auto sample_position = get_or_add_attribute<Vec3, SurfaceVertex>(surface, "position");
		auto sample_normal = get_or_add_attribute<Vec3, SurfaceVertex>(surface, "normal");
		auto vertex_cluster_info =
			get_or_add_attribute<std::pair<uint32, PointVertex>, SurfaceVertex>(surface, "vertex_cluster_info");
		auto vertex_cluster_color = get_or_add_attribute<Vec3, SurfaceVertex>(surface, "cluster_color");
		auto medial_axis_samples_position =
			get_or_add_attribute<Vec3, SurfaceVertex>(surface, "medial_axis_samples_position");
		auto medial_axis_samples_radius =
			get_or_add_attribute<Scalar, SurfaceVertex>(surface, "medial_axis_samples_radius");
		auto medial_axis_samples_closest_points_ =
			get_attribute<std::pair<SurfaceVertex, SurfaceVertex>, SurfaceVertex>(surface, "medial_axis_samples_closest_points");
		auto distance_to_cluster = get_or_add_attribute<Scalar, SurfaceVertex>(surface, "distance_to_cluster");
		auto bvh_vertex_index = get_or_add_attribute<uint32, SurfaceVertex>(surface, "bvh_vertex_index");
		auto cloest_point_color = get_or_add_attribute<Vec3, SurfaceVertex>(surface, "cloest_point_color");
		MeshData<SURFACE>& md = surface_provider_->mesh_data(surface);
		auto filtered_medial_axis_samples =
			&md.template get_or_add_cells_set<SurfaceVertex>("filtered_medial_axis_samples");
	
		POINT* clusters = point_provider_->add_mesh(point_provider_->mesh_name(surface) + "clusters");
		auto cluster_color = get_or_add_attribute<Vec3, PointVertex>(*clusters, "color");
		auto cluster_positions = get_or_add_attribute<Vec3, PointVertex>(*clusters, "position");
		auto cluster_radius = get_or_add_attribute<Scalar, PointVertex>(*clusters, "clusters_radius");
		auto clusters_infos = get_or_add_attribute<Cluster_Info, PointVertex>(*clusters, "clusters_infos");
		auto cluster_cloest_sample = get_or_add_attribute<std::pair<SurfaceVertex,SurfaceVertex>,PointVertex>(*clusters, "cluster_cloest_sample");
		
		parallel_foreach_cell(surface, [&](SurfaceVertex sv) {
			value<Vec3>(surface, cloest_point_color, sv) = Vec3(0, 0, 0);
			return true;
		});
		
		load_model_in_cgal(surface, csm);
		tree = new Tree(faces(csm).first, faces(csm).second, csm);
		tree->accelerate_distance_queries();
		inside_tester = new Point_inside(*tree);

		SphereQueue sphere_queue;
		auto sphere_info = get_or_add_attribute<SphereInfo, SurfaceVertex>(surface, "sphere_info");
		//Push all the medial axis sphere into the maximum queue 
		filtered_medial_axis_samples->foreach_cell([&](SurfaceVertex v) {
			value<SphereInfo>(surface, sphere_info, v) = {
				true, sphere_queue.emplace(value<Scalar>(surface, medial_axis_samples_radius, v), v)};
			medial_vertex_position_vector.push_back(value<Vec3>(surface, medial_axis_samples_position, v));
			return true;
		});
		uint32 id = 0;
		foreach_cell(surface, [&](SurfaceVertex v) {
			value<uint32>(surface, bvh_vertex_index, v) = id++;
			surface_kdt_vertices.push_back(v);
			surface_vertex_position_vector.push_back(value<Vec3>(surface, sample_position, v));
			
			return true;
			});
		uint32 nb_faces = nb_cells<SurfaceFace>(surface);
		std::vector<SurfaceFace> bvh_faces;
		bvh_faces.reserve(nb_faces);
		std::vector<uint32> face_vertex_indices;
		face_vertex_indices.reserve(nb_faces * 3);
		
		foreach_cell(surface, [&](SurfaceFace f) -> bool {
			bvh_faces.push_back(f);
			foreach_incident_vertex(surface, f, [&](SurfaceVertex v) -> bool {
				face_vertex_indices.push_back(value<uint32>(surface, bvh_vertex_index, v));
				return true;
			});
			return true;
		});

		sqem_helper =
			std::make_unique<modeling::ClusteringSQEM_Helper<Surface>>(surface, sample_position.get(), sample_normal.get());
		surface_bvh = std::make_unique<acc::BVHTree<uint32, Vec3>>(face_vertex_indices, surface_vertex_position_vector);
		surface_kdt = std::make_unique<acc::KDTree<3, uint32>>(surface_vertex_position_vector);
		medial_kdt = std::make_unique<acc::KDTree<3, uint32>>(medial_vertex_position_vector);
		CellMarker<SURFACE, SurfaceVertex> bfs_marker(surface);
		bfs_marker.unmark_all();
		std::queue<SurfaceVertex> vertex_queue;
		//Start from the biggest sphere, affect the surface points to the sphere if the distance 
		// between surface point and sphere is less than threshold
		while (sphere_queue.size() > 0)
		{
			SphereQueueIt it = sphere_queue.begin();
			auto [radius, v] = *it;
			if (!value<SphereInfo>(surface, sphere_info, v).first)
			{
				sphere_queue.erase(it);
				continue;
			}
			auto [v1, v2] =
				value<std::pair<SurfaceVertex, SurfaceVertex>>(surface, medial_axis_samples_closest_points_, v);
			Vec3& sphere_center = value<Vec3>(surface, medial_axis_samples_position, v1);
			value<SphereInfo>(surface, sphere_info, v1).first = false;
			if (!value<SphereInfo>(surface, sphere_info, v2).first)
			{
				continue;
			}
			
			//Set new clusters info
			PointVertex new_cluster = add_vertex(*clusters);
			value<Vec3>(*clusters, cluster_positions, new_cluster) = sphere_center;
			value<std::pair<SurfaceVertex, SurfaceVertex>>(*clusters, cluster_cloest_sample, new_cluster) = {v1, v2};
			Cluster_Info& cf = value<Cluster_Info>(*clusters, clusters_infos, new_cluster);
			value<Vec3>(*clusters, cluster_color, new_cluster) =
				Vec3(rand() / (double)RAND_MAX, rand() / (double)RAND_MAX, rand() / (double)RAND_MAX);
			cf.cluster_variance = 0;
			cf.cluster_vertices.clear();
			value<Scalar>(*clusters, cluster_radius, new_cluster) = radius;
			value<Vec3>(surface, cloest_point_color, v1) = value<Vec3>(*clusters, cluster_color, new_cluster);
			value<Vec3>(surface, cloest_point_color, v2) = value<Vec3>(*clusters, cluster_color, new_cluster);

			// Add v1 into cluster and remove it from the queue
			vertex_queue.push(v1);
			value<std::pair<uint32, PointVertex>>(surface, vertex_cluster_info, v1) =
				std::make_pair(index_of(*clusters, new_cluster), new_cluster);
			value<Vec3>(surface, vertex_cluster_color, v1) = value<Vec3>(*clusters, cluster_color, new_cluster);
			value<Scalar>(surface, distance_to_cluster, v1) = 0;
			/*value<SphereInfo>(surface, sphere_info, v1).first = false;*/
			
			// Add v2 into cluster and remove it from the queue
			vertex_queue.push(v2);
			value<std::pair<uint32, PointVertex>>(surface, vertex_cluster_info, v2) =
				std::make_pair(index_of(*clusters, new_cluster), new_cluster);
			value<Vec3>(surface, vertex_cluster_color, v2) = value<Vec3>(*clusters, cluster_color, new_cluster);
			value<Scalar>(surface, distance_to_cluster, v2) = 0;
			value<SphereInfo>(surface, sphere_info, v2).first = false;
			
			//Start BFS
			bfs_marker.unmark_all();
			while (vertex_queue.size() > 0)
			{
				SurfaceVertex v = vertex_queue.front();
				vertex_queue.pop();
				// Affect the surface points to the cluster
				foreach_adjacent_vertex_through_edge(surface, v, [&](SurfaceVertex sv) {
					if (bfs_marker.is_marked(sv))
						return true;
					Vec3& pos = value<Vec3>(surface, sample_position, sv);
					Vec4 sphere = Vec4(sphere_center.x(), sphere_center.y(), sphere_center.z(), radius);
					// Scalar distance = (pos - sphere_center).norm() - radius;
					// Scalar error = sqem_helper->vertex_cost(sv, sphere);
					// error = distance  + error;
					//Scalar cosine = (pos-sphere_center).normalized().dot(value<Vec3>(surface, sample_normal, sv));
					// if the distance is less than threshold
					Scalar error = (pos - sphere_center).norm() - radius;
					if (error < dilation_factor /*&& cosine>=0*/)
					{
						vertex_queue.push(sv);
						if (!value<SphereInfo>(surface,sphere_info,sv).first)
						{
							if (error < value<Scalar>(surface, distance_to_cluster, sv))
							{
								value<std::pair<uint32, PointVertex>>(surface, vertex_cluster_info, sv) =
									std::make_pair(index_of(*clusters, new_cluster), new_cluster);
								value<Vec3>(surface, vertex_cluster_color, sv) =
									value<Vec3>(*clusters, cluster_color, new_cluster);
								value<Scalar>(surface, distance_to_cluster, sv) = error;
							}
						}
						else
						{
							value<std::pair<uint32, PointVertex>>(surface, vertex_cluster_info, sv) =
								std::make_pair(index_of(*clusters, new_cluster), new_cluster);
							value<Vec3>(surface, vertex_cluster_color, sv) =
								value<Vec3>(*clusters, cluster_color, new_cluster);
							value<Scalar>(surface, distance_to_cluster, sv) = error;
							value<SphereInfo>(surface, sphere_info, sv).first = false;
							
						}
					}
					bfs_marker.mark(sv);
					return true;
				});
			}
		}
		assign_cluster(surface, *clusters);
		if (cluster_positions)
			point_provider_->set_mesh_bb_vertex_position(*clusters, cluster_positions);
		point_provider_->emit_connectivity_changed(*clusters);
		surface_provider_->emit_attribute_changed(surface, vertex_cluster_color.get());
	}
	
	void iniialise_cluster2(SURFACE& surface)
	{
		auto sample_position = get_or_add_attribute<Vec3, SurfaceVertex>(surface, "position");
		auto sample_normal = get_or_add_attribute<Vec3, SurfaceVertex>(surface, "normal");
		auto vertex_cluster_info =
			get_or_add_attribute<std::pair<uint32, PointVertex>, SurfaceVertex>(surface, "vertex_cluster_info");
		auto vertex_cluster_color = get_or_add_attribute<Vec3, SurfaceVertex>(surface, "cluster_color");
		auto medial_axis_samples_position =
			get_or_add_attribute<Vec3, SurfaceVertex>(surface, "medial_axis_samples_position");
		auto medial_axis_samples_radius =
			get_or_add_attribute<Scalar, SurfaceVertex>(surface, "medial_axis_samples_radius");
		auto medial_axis_samples_closest_points_ =
			get_attribute<std::pair<SurfaceVertex, SurfaceVertex>, SurfaceVertex>(surface,
																				  "medial_axis_samples_closest_points");
		auto distance_to_cluster = get_or_add_attribute<Scalar, SurfaceVertex>(surface, "distance_to_cluster");
		auto bvh_vertex_index = get_or_add_attribute<uint32, SurfaceVertex>(surface, "bvh_vertex_index");
		auto cloest_point_color = get_or_add_attribute<Vec3, SurfaceVertex>(surface, "cloest_point_color");
		MeshData<SURFACE>& md = surface_provider_->mesh_data(surface);
		auto filtered_medial_axis_samples =
			&md.template get_or_add_cells_set<SurfaceVertex>("filtered_medial_axis_samples");

		POINT* clusters = point_provider_->add_mesh(point_provider_->mesh_name(surface) + "clusters");
		auto cluster_color = get_or_add_attribute<Vec3, PointVertex>(*clusters, "color");
		auto cluster_positions = get_or_add_attribute<Vec3, PointVertex>(*clusters, "position");
		auto cluster_radius = get_or_add_attribute<Scalar, PointVertex>(*clusters, "clusters_radius");
		auto clusters_infos = get_or_add_attribute<Cluster_Info, PointVertex>(*clusters, "clusters_infos");
		auto cluster_cloest_sample = get_or_add_attribute<std::pair<SurfaceVertex, SurfaceVertex>, PointVertex>(
			*clusters, "cluster_cloest_sample");

		parallel_foreach_cell(surface, [&](SurfaceVertex sv) {
			value<Vec3>(surface, cloest_point_color, sv) = Vec3(0, 0, 0);
			return true;
		});

		load_model_in_cgal(surface, csm);
		tree = new Tree(faces(csm).first, faces(csm).second, csm);
		tree->accelerate_distance_queries();
		inside_tester = new Point_inside(*tree);

		SphereQueue sphere_queue;
		auto sphere_info = get_or_add_attribute<SphereInfo, SurfaceVertex>(surface, "sphere_info");
		// Push all the medial axis sphere into the maximum queue
		filtered_medial_axis_samples->foreach_cell([&](SurfaceVertex v) {
			value<SphereInfo>(surface, sphere_info, v) = {
				true, sphere_queue.emplace(value<Scalar>(surface, medial_axis_samples_radius, v), v)};
			medial_vertex_position_vector.push_back(value<Vec3>(surface, medial_axis_samples_position, v));
			return true;
		});
		uint32 id = 0;
		foreach_cell(surface, [&](SurfaceVertex v) {
			value<uint32>(surface, bvh_vertex_index, v) = id++;
			surface_kdt_vertices.push_back(v);
			surface_vertex_position_vector.push_back(value<Vec3>(surface, sample_position, v));

			return true;
		});
		uint32 nb_faces = nb_cells<SurfaceFace>(surface);
		std::vector<SurfaceFace> bvh_faces;
		bvh_faces.reserve(nb_faces);
		std::vector<uint32> face_vertex_indices;
		face_vertex_indices.reserve(nb_faces * 3);

		foreach_cell(surface, [&](SurfaceFace f) -> bool {
			bvh_faces.push_back(f);
			foreach_incident_vertex(surface, f, [&](SurfaceVertex v) -> bool {
				face_vertex_indices.push_back(value<uint32>(surface, bvh_vertex_index, v));
				return true;
			});
			return true;
		});

		sqem_helper = std::make_unique<modeling::ClusteringSQEM_Helper<Surface>>(surface, sample_position.get(),
																				 sample_normal.get());
		surface_bvh = std::make_unique<acc::BVHTree<uint32, Vec3>>(face_vertex_indices, surface_vertex_position_vector);
		surface_kdt = std::make_unique<acc::KDTree<3, uint32>>(surface_vertex_position_vector);
		medial_kdt = std::make_unique<acc::KDTree<3, uint32>>(medial_vertex_position_vector);
		CellMarker<SURFACE, SurfaceVertex> bfs_marker(surface);
		bfs_marker.unmark_all();
		std::queue<SurfaceVertex> vertex_queue;
		// Start from the biggest sphere, affect the surface points to the sphere if the distance
		//  between surface point and sphere is less than threshold
		while (sphere_queue.size() > 0)
		{
			SphereQueueIt it = sphere_queue.begin();
			auto [radius, v] = *it;
			if (!value<SphereInfo>(surface, sphere_info, v).first)
			{
				sphere_queue.erase(it);
				continue;
			}
			auto [v1, v2] =
				value<std::pair<SurfaceVertex, SurfaceVertex>>(surface, medial_axis_samples_closest_points_, v);
			Vec3& sphere_center = value<Vec3>(surface, medial_axis_samples_position, v1);
			value<SphereInfo>(surface, sphere_info, v1).first = false;

			// Set new clusters info
			PointVertex new_cluster = add_vertex(*clusters);
			value<Vec3>(*clusters, cluster_positions, new_cluster) = sphere_center;
			value<std::pair<SurfaceVertex, SurfaceVertex>>(*clusters, cluster_cloest_sample, new_cluster) = {v1, v2};
			Cluster_Info& cf = value<Cluster_Info>(*clusters, clusters_infos, new_cluster);
			value<Vec3>(*clusters, cluster_color, new_cluster) =
				Vec3(rand() / (double)RAND_MAX, rand() / (double)RAND_MAX, rand() / (double)RAND_MAX);
			cf.cluster_variance = 0;
			cf.cluster_vertices.clear();
			value<Scalar>(*clusters, cluster_radius, new_cluster) = radius;
			value<Vec3>(surface, cloest_point_color, v1) = value<Vec3>(*clusters, cluster_color, new_cluster);
			value<Vec3>(surface, cloest_point_color, v2) = value<Vec3>(*clusters, cluster_color, new_cluster);

			foreach_cell(surface, [&](SurfaceVertex sv) {
				Vec3& center = value<Vec3>(surface, medial_axis_samples_position, sv);
				if (((center - sphere_center).norm() - radius) < 0.06)
				{
					value<SphereInfo>(surface, sphere_info, sv).first = false;
					value<std::pair<uint32, PointVertex>>(surface, vertex_cluster_info, sv) =
						std::make_pair(index_of(*clusters, new_cluster), new_cluster);
					value<Vec3>(surface, vertex_cluster_color, sv) = value<Vec3>(*clusters, cluster_color, new_cluster);
					
				}
				return true;
				});
		}
		assign_cluster(surface, *clusters);
		if (cluster_positions)
			point_provider_->set_mesh_bb_vertex_position(*clusters, cluster_positions);
		point_provider_->emit_connectivity_changed(*clusters);
		surface_provider_->emit_attribute_changed(surface, vertex_cluster_color.get());
	}
	void assign_cluster(SURFACE& surface, POINT& clusters)
	{
		auto sample_position = get_or_add_attribute<Vec3, SurfaceVertex>(surface, "position");
		auto sample_normal = get_or_add_attribute<Vec3, SurfaceVertex>(surface, "normal");
		auto vertex_cluster_info =
			get_or_add_attribute<std::pair<uint32, PointVertex>, SurfaceVertex>(surface, "vertex_cluster_info");
		auto vertex_cluster_color = get_or_add_attribute<Vec3, SurfaceVertex>(surface, "cluster_color");
		auto distance_to_cluster = get_or_add_attribute<Scalar, SurfaceVertex>(surface, "distance_to_cluster");
		auto euclidean_distance = get_or_add_attribute<Scalar, SurfaceVertex>(surface, "euclidean_distance");
		auto medial_axis_samples_closest_points_ =
			get_attribute<std::pair<SurfaceVertex, SurfaceVertex>, SurfaceVertex>(surface,
																				  "medial_axis_samples_closest_points");
		auto neighbours_set =
			get_or_add_attribute<std::unordered_map<uint32, PointVertex>, PointVertex>(clusters, "neighbours_set");

		MeshData<SURFACE>& md = surface_provider_->mesh_data(surface);
		auto filtered_medial_axis_samples =
			&md.template get_or_add_cells_set<SurfaceVertex>("filtered_medial_axis_samples");
		auto cluster_color = get_or_add_attribute<Vec3, PointVertex>(clusters, "color");
		auto cluster_positions = get_or_add_attribute<Vec3, PointVertex>(clusters, "position");
		auto cluster_radius = get_or_add_attribute<Scalar, PointVertex>(clusters, "clusters_radius");
		auto clusters_infos = get_or_add_attribute<Cluster_Info, PointVertex>(clusters, "clusters_infos");
		auto cluster_cloest_sample =
			get_or_add_attribute<std::pair<SurfaceVertex,SurfaceVertex>, PointVertex>(clusters, "cluster_cloest_sample");

		CellMarker<SURFACE, SurfaceVertex> bfs_marker(surface);
		CellMarker<SURFACE, SurfaceVertex> cluster_marker(surface);
		bfs_marker.unmark_all();
		cluster_marker.unmark_all();
		std::queue<std::pair<SurfaceVertex, PointVertex>> vertex_queue;
		foreach_cell(clusters, [&](PointVertex pv) {
			value<Cluster_Info>(clusters, clusters_infos, pv).cluster_vertices.clear();
			return true;
		});
		/*foreach_cell(clusters, [&](PointVertex pv) {
			Cluster_Info& cf = value<Cluster_Info>(clusters, clusters_infos, pv);
			foreach_cell(surface, [&](SurfaceVertex sv) {
				Vec3& sphere_center = value<Vec3>(clusters, cluster_positions, pv);
				Vec3& pos = value<Vec3>(surface, sample_position, sv);
				Scalar distance = (pos - sphere_center).norm() - value<Scalar>(clusters, cluster_radius, pv);
				if (cluster_marker.is_marked(sv)) {
					if (distance < value<Scalar>(surface, distance_to_cluster, sv))
					{
						value<std::pair<uint32, PointVertex>>(surface, vertex_cluster_info, sv) =
							std::make_pair(index_of(clusters, pv), pv);
						value<Vec3>(surface, vertex_cluster_color, sv) = value<Vec3>(clusters, cluster_color, pv);
						value<Scalar>(surface, distance_to_cluster, sv) = distance;
					}
					else
					{
						value<std::unordered_map<uint32, PointVertex>>(clusters, neighbours_set, pv)
							.insert(value<std::pair<uint32, PointVertex>>(surface, vertex_cluster_info, sv));
					}
				}
				else {
					value<std::pair<uint32, PointVertex>>(surface, vertex_cluster_info, sv) =
						std::make_pair(index_of(clusters, pv), pv);
					value<Vec3>(surface, vertex_cluster_color, sv) = value<Vec3>(clusters, cluster_color, pv);
					value<Scalar>(surface, distance_to_cluster, sv) = distance;
					cluster_marker.mark(sv);
				}
				return true;
				});


			return true;
			});*/
		 foreach_cell(clusters, [&](PointVertex pv) {
			Cluster_Info& cf = value<Cluster_Info>(clusters, clusters_infos, pv);

			auto [v1, v2] = value<std::pair<SurfaceVertex, SurfaceVertex>>(clusters, cluster_cloest_sample, pv);
			Vec3& sphere_center = value<Vec3>(clusters, cluster_positions, pv);
			Vec3& pos_v1 = value<Vec3>(surface, sample_position, v1);
			Scalar distance_v1 = (pos_v1 - sphere_center).norm() - value<Scalar>(clusters, cluster_radius, pv);
			if (cluster_marker.is_marked(v1))
			{
				if (distance_v1 < value<Scalar>(surface, distance_to_cluster, v1))
				{
					value<std::pair<uint32, PointVertex>>(surface, vertex_cluster_info, v1) =
						std::make_pair(index_of(clusters, pv), pv);
					value<Vec3>(surface, vertex_cluster_color, v1) = value<Vec3>(clusters, cluster_color, pv);
					value<Scalar>(surface, distance_to_cluster, v1) = distance_v1;
				}
			}
			else
			{
				value<std::pair<uint32, PointVertex>>(surface, vertex_cluster_info, v1) =
					std::make_pair(index_of(clusters, pv), pv);
				value<Vec3>(surface, vertex_cluster_color, v1) = value<Vec3>(clusters, cluster_color, pv);
				value<Scalar>(surface, distance_to_cluster, v1) = distance_v1;
				cluster_marker.mark(v1);
			}

			vertex_queue.push({v1, pv});
			Vec3& pos_v2 = value<Vec3>(surface, sample_position, v2);
			Scalar distance_v2 = (pos_v2 - sphere_center).norm() - value<Scalar>(clusters, cluster_radius, pv);

			if (cluster_marker.is_marked(v2))
			{
				if (distance_v2 < value<Scalar>(surface, distance_to_cluster, v2))
				{
					value<std::pair<uint32, PointVertex>>(surface, vertex_cluster_info, v2) =
						std::make_pair(index_of(clusters, pv), pv);
					value<Vec3>(surface, vertex_cluster_color, v2) = value<Vec3>(clusters, cluster_color, pv);
					value<Scalar>(surface, distance_to_cluster, v2) = distance_v2;
				}
			}
			else
			{
				value<std::pair<uint32, PointVertex>>(surface, vertex_cluster_info, v2) =
					std::make_pair(index_of(clusters, pv), pv);
				value<Vec3>(surface, vertex_cluster_color, v2) = value<Vec3>(clusters, cluster_color, pv);
				value<Scalar>(surface, distance_to_cluster, v2) = distance_v2;
				cluster_marker.mark(v2);
			}

			vertex_queue.push({v2, pv});
			return true;
		});
		// Start BFS
		//bfs_marker.unmark_all();
		while (vertex_queue.size() > 0)
		{
			auto [v, pv] = vertex_queue.front();
			vertex_queue.pop();
			// Affect the surface points to the cluster
			foreach_adjacent_vertex_through_edge(surface, v, [&](SurfaceVertex sv) {
				
				Vec3& sphere_center = value<Vec3>(clusters, cluster_positions, pv);
				Vec3& pos = value<Vec3>(surface, sample_position, sv);
				Scalar distance = (pos - sphere_center).norm() - value<Scalar>(clusters, cluster_radius, pv);
				
				//if encounter a vertex who has already been affected
				if (cluster_marker.is_marked(sv))
				{
					PointVertex pv2 = value<std::pair<uint32, PointVertex>>(surface, vertex_cluster_info, sv).second;
					if (index_of(clusters, pv) != index_of(clusters, pv2))
					{

						// if the point-cluster distance is smaller, change the label
						if (distance < value<Scalar>(surface, distance_to_cluster, sv))
						{
							vertex_queue.push({sv, pv});
							value<std::pair<uint32, PointVertex>>(surface, vertex_cluster_info, sv) =
								std::make_pair(index_of(clusters, pv), pv);
							value<Vec3>(surface, vertex_cluster_color, sv) = value<Vec3>(clusters, cluster_color, pv);
							value<Scalar>(surface, distance_to_cluster, sv) = distance;
						}
						// else add adjacency info for the current cluster
						else
						{
							value<std::unordered_map<uint32, PointVertex>>(clusters, neighbours_set, pv)
								.insert(value<std::pair<uint32, PointVertex>>(surface, vertex_cluster_info, sv));
						}
					}
				}
				else
				{
					vertex_queue.push({sv, pv});
					value<std::pair<uint32, PointVertex>>(surface, vertex_cluster_info, sv) =
						std::make_pair(index_of(clusters, pv), pv);
					value<Vec3>(surface, vertex_cluster_color, sv) = value<Vec3>(clusters, cluster_color, pv);
					value<Scalar>(surface, distance_to_cluster, sv) = distance;
					cluster_marker.mark(sv);
				}
				return true;
			});
		}

		foreach_cell(surface, [&](SurfaceVertex sv) { 
			std::pair<uint32, PointVertex>&cluster_info = value<std::pair<uint32, PointVertex>>(surface, vertex_cluster_info, sv);
			PointVertex cluster_point = cluster_info.second;
			Cluster_Info& cf = value<Cluster_Info>(clusters, clusters_infos, cluster_point);
			cf.cluster_vertices.push_back(sv);
			return true;
		});
		std::vector<PointVertex> clusters_to_remove;
		foreach_cell(clusters, [&](PointVertex pv) {
			Cluster_Info& cf = value<Cluster_Info>(clusters, clusters_infos, pv);
			if (cf.cluster_vertices.size() == 0)
			{
				clusters_to_remove.push_back(pv);
			}
			return true;
		});
		for (PointVertex& sv : clusters_to_remove)
		{
			std::cout << "delete cluster" << index_of(clusters, sv) << std::endl;
			remove_vertex(clusters, sv);
		}
		Scalar error = 0;
		foreach_cell(surface, [&](SurfaceVertex sv) { 
			error += value<Scalar>(surface, distance_to_cluster, sv);
			//std::cout << "distance: " << value<Scalar>(surface, distance_to_cluster, sv) << std::endl;
			return true;
		});
		std::cout << "Global distance: " << error << std::endl; 
		
		surface_provider_->emit_attribute_changed(surface, vertex_cluster_color.get());
		surface_provider_->emit_attribute_changed(surface, distance_to_cluster.get());
	}
	
	
	
	void update_filtered_cluster(SURFACE& surface, POINT& clusters)
	{
		
		auto sample_position = get_or_add_attribute<Vec3, SurfaceVertex>(surface, "position");
		auto sample_normal = get_or_add_attribute<Vec3, SurfaceVertex>(surface, "normal");
		auto cluster_color = get_or_add_attribute<Vec3, PointVertex>(clusters, "color");
		auto cluster_position = get_attribute<Vec3, PointVertex>(clusters, "position");
		auto clusters_infos = get_or_add_attribute<Cluster_Info, PointVertex>(clusters, "clusters_infos");
		auto clusters_radius = get_or_add_attribute<Scalar, PointVertex>(clusters, "clusters_radius");
		auto vertex_cluster_info =
			get_or_add_attribute<std::pair<uint32, PointVertex>, SurfaceVertex>(surface, "vertex_cluster_info");
		auto distance_to_cluster = get_or_add_attribute<Scalar, SurfaceVertex>(surface, "distance_to_cluster");
		auto cluster_cloest_sample = get_or_add_attribute<std::pair<SurfaceVertex, SurfaceVertex>, PointVertex>(
			clusters, "cluster_cloest_sample");
		auto medial_axis_samples_position =
			get_or_add_attribute<Vec3, SurfaceVertex>(surface, "medial_axis_samples_position");
		auto medial_axis_sample_radius_ =
			get_or_add_attribute<Scalar, SurfaceVertex>(surface, "medial_axis_samples_radius");
		auto medial_axis_samples_closest_points =
			get_attribute<std::pair<SurfaceVertex, SurfaceVertex>, SurfaceVertex>(surface,
																				  "medial_axis_samples_closest_points");
		auto medial_axis_samples_feature_value = get_or_add_attribute<Scalar, SurfaceVertex>(surface, "medial_axis_samples_feature_value");
		auto medial_axis_samples_weight = get_or_add_attribute<Scalar, SurfaceVertex>(surface, "medial_axis_samples_weight");
		auto cloest_point_color = get_or_add_attribute<Vec3, SurfaceVertex>(surface, "cloest_point_color");
		parallel_foreach_cell(surface, [&](SurfaceVertex sv) { 
			value<Vec3>(surface, cloest_point_color, sv) = Vec3(0, 0, 0);
			return true;
		});
		auto kmax = get_or_add_attribute<Scalar, SurfaceVertex>(surface, "kmax");
		//For each cluster, find the nearest medial point
		foreach_cell(clusters, [&](PointVertex pv) {
			Cluster_Info& cf = value<Cluster_Info>(clusters, clusters_infos, pv);
			Vec3 opti_coord;
			switch (clustering_mode)
			{
			case 0: {
				Vec3 sum_coord = Vec3(0, 0, 0);
				Scalar weight = 0;
				for (SurfaceVertex sv1 : cf.cluster_vertices)
				{
					Scalar w = 1;
					sum_coord += value<Vec3>(surface, medial_axis_samples_position, sv1) * w;
					weight += w;
				}

				opti_coord = sum_coord / weight;
				Vec3 original_coord = opti_coord;
				auto [radius, v1, v2] = geometry::non_linear_solver(surface, sample_position.get(), sample_normal.get(),
																	cf.cluster_vertices, opti_coord, surface_kdt.get(),
																	medial_kdt.get(), surface_bvh.get());
 				/*auto [radius, v1, v2] = geometry::move_point_to_medial_axis(
					surface, sample_position.get(), sample_normal.get(), surface_kdt_vertices, opti_coord,
					surface_kdt.get(), medial_kdt.get(), surface_bvh.get());*/
				//PointVertex pv = value<std::pair<uint32, PointVertex>>(surface, vertex_cluster_info, v1).second;*/
				/*SurfaceVertex v2 =
					value<std::pair<SurfaceVertex, SurfaceVertex>>(surface, medial_axis_samples_closest_points, v1)
						.second;*/
				//std::cout << "cloest medial vertex index: " << index_of(surface, v) << std::endl;
				/*auto [v1, v2] =
					value<std::pair<SurfaceVertex, SurfaceVertex>>(surface, medial_axis_samples_closest_points, v);*/
				/*value<Vec3>(surface, cloest_point_color, v1) = value<Vec3>(clusters, cluster_color, pv);
				value<Vec3>(surface, cloest_point_color, v2) = value<Vec3>(clusters, cluster_color, pv);*/
				value<Vec3>(clusters, cluster_position, pv) = opti_coord;
				value<Scalar>(clusters, clusters_radius, pv) = radius;
				//value<std::pair<SurfaceVertex, SurfaceVertex>>(clusters, cluster_cloest_sample, pv) = {v1, v2};
				
			}
			break;
			case 1: {
				Vec3 sum_coord = Vec3(0, 0, 0);
				Scalar weight = 0;
				for (SurfaceVertex sv1 : cf.cluster_vertices)
				{
					Scalar w = 1;
					sum_coord += value<Vec3>(surface, sample_position, sv1) * w;
					weight += w;
				}

				opti_coord = sum_coord / weight;
				Vec3 original_coord = opti_coord;
				auto [radius, v1, v2] = geometry::move_point_to_medial_axis(
					surface, sample_position.get(), sample_normal.get(), surface_kdt_vertices, opti_coord,
					surface_kdt.get(), medial_kdt.get(), surface_bvh.get());
				// PointVertex pv = value<std::pair<uint32, PointVertex>>(surface, vertex_cluster_info, v1).second;*/
				/*SurfaceVertex v2 =
					value<std::pair<SurfaceVertex, SurfaceVertex>>(surface, medial_axis_samples_closest_points, v1)
						.second;*/
				// std::cout << "cloest medial vertex index: " << index_of(surface, v) << std::endl;
				/* auto [v1, v2] =
					value<std::pair<SurfaceVertex, SurfaceVertex>>(surface, medial_axis_samples_closest_points, v);*/
				value<Vec3>(surface, cloest_point_color, v1) = value<Vec3>(clusters, cluster_color, pv);
				value<Vec3>(surface, cloest_point_color, v2) = value<Vec3>(clusters, cluster_color, pv);
				value<Vec3>(clusters, cluster_position, pv) = opti_coord;
				value<Scalar>(clusters, clusters_radius, pv) = radius;
				value<std::pair<SurfaceVertex, SurfaceVertex>>(clusters, cluster_cloest_sample, pv) = {v1, v2};
			}
			break;
			case 2: {
				auto [center, r] = sphere_fitting_algo(surface, cf.cluster_vertices);
				
				auto [radius, v1, v2] = geometry::move_point_to_medial_axis(
					surface, sample_position.get(), sample_normal.get(), surface_kdt_vertices, center,
					surface_kdt.get(), medial_kdt.get(), surface_bvh.get());
				value<Vec3>(surface, cloest_point_color, v1) = value<Vec3>(clusters, cluster_color, pv);
				value<Vec3>(surface, cloest_point_color, v2) = value<Vec3>(clusters, cluster_color, pv);
				value<Vec3>(clusters, cluster_position, pv) = center;
				value<Scalar>(clusters, clusters_radius, pv) = radius;
				value<std::pair<SurfaceVertex, SurfaceVertex>>(clusters, cluster_cloest_sample, pv) = {v1, v2};
			}
			break;
			case 3: {
				Vec3 sum_coord = Vec3(0, 0, 0);
				Scalar weight = 0;
				for (SurfaceVertex sv1 : cf.cluster_vertices)
				{
					Scalar w = 100 * value<Scalar>(surface, medial_axis_sample_radius_, sv1) + 
						value<Scalar>(surface, medial_axis_samples_weight, sv1);
					sum_coord += value<Vec3>(surface, medial_axis_samples_position, sv1) * w;
					weight += w;
				}

				opti_coord = sum_coord / weight;
				Vec3 original_coord = opti_coord;
				auto [radius, v1,v2] = geometry::move_point_to_medial_axis(
					surface, sample_position.get(), sample_normal.get(), surface_kdt_vertices, opti_coord,
					surface_kdt.get(), medial_kdt.get(), surface_bvh.get());
				//PointVertex pv = value<std::pair<uint32, PointVertex>>(surface, vertex_cluster_info, v1).second;
				/*SurfaceVertex v2 =
					value<std::pair<SurfaceVertex, SurfaceVertex>>(surface, medial_axis_samples_closest_points, v1)
						.second;*/
				//std::cout << "cloest medial vertex index: " << index_of(surface, v) << std::endl;
				/*auto [v1, v2] =
					value<std::pair<SurfaceVertex, SurfaceVertex>>(surface, medial_axis_samples_closest_points, v);*/
				value<Vec3>(surface, cloest_point_color, v1) = value<Vec3>(clusters, cluster_color, pv);
				value<Vec3>(surface, cloest_point_color, v2) = value<Vec3>(clusters, cluster_color, pv);
				value<Vec3>(clusters, cluster_position, pv) = opti_coord;
				value<Scalar>(clusters, clusters_radius, pv) = radius;
				value<std::pair<SurfaceVertex, SurfaceVertex>>(clusters, cluster_cloest_sample, pv) = {v1, v2};
			}
			break;
			case 4: {
				Vec3 sum_coord = Vec3(0, 0, 0);
				Scalar weight = 0;
				for (SurfaceVertex sv1 : cf.cluster_vertices)
				{
					/*Scalar w = 100 * value<Scalar>(surface, medial_axis_sample_radius_, sv1) +
							   100 * value<Scalar>(surface, medial_axis_samples_feature_value, sv1);*/
					Scalar w = value<Scalar>(surface, medial_axis_samples_weight, sv1) +
							   value<Scalar>(surface, medial_axis_sample_radius_, sv1);
					
					sum_coord += value<Vec3>(surface, sample_position, sv1) * w;
					weight += w;
				}

				opti_coord = sum_coord / weight;
				Vec3 original_coord = opti_coord;
				auto [radius, v1, v2] = geometry::move_point_to_medial_axis(
					surface, sample_position.get(), sample_normal.get(), surface_kdt_vertices, opti_coord,
					surface_kdt.get(), medial_kdt.get(), surface_bvh.get());
				PointVertex pv = value<std::pair<uint32, PointVertex>>(surface, vertex_cluster_info, v1).second;
				/*SurfaceVertex v2 =
					value<std::pair<SurfaceVertex, SurfaceVertex>>(surface, medial_axis_samples_closest_points, v1)
						.second;*/
				// std::cout << "cloest medial vertex index: " << index_of(surface, v) << std::endl;
				/*auto [v1, v2] =
					value<std::pair<SurfaceVertex, SurfaceVertex>>(surface, medial_axis_samples_closest_points, v);*/
				
				value<Vec3>(surface, cloest_point_color, v1) = value<Vec3>(clusters, cluster_color, pv);
				value<Vec3>(surface, cloest_point_color, v2) = value<Vec3>(clusters, cluster_color, pv);
				value<Vec3>(clusters, cluster_position, pv) = opti_coord;
				value<Scalar>(clusters, clusters_radius, pv) = radius;
				value<std::pair<SurfaceVertex, SurfaceVertex>>(clusters, cluster_cloest_sample, pv) = {v1, v2};
			}
			break;
			case 5: {
				std::vector<Point> surface_points;
				for (SurfaceVertex sv1 : cf.cluster_vertices)
				{
					Vec3 pos = value<Vec3>(surface, sample_position, sv1);
					surface_points.push_back(Point(pos.x(), pos.y(), pos.z()));
				}
				Min_sphere min_sphere(surface_points.begin(), surface_points.end());
				assert(min_sphere.is_valid());
				auto it = min_sphere.center_cartesian_begin();
				K::FT coord[3];
				for (size_t idx = 0; idx < 3 && it != min_sphere.center_cartesian_end(); idx++)
				{
					coord[idx] = *it;
					it++;
				}
				opti_coord = Vec3(coord[0], coord[1], coord[2]);
				K::FT radius = min_sphere.radius();
				auto [r, v1, v2] = geometry::move_point_to_medial_axis(
					surface, sample_position.get(), sample_normal.get(), surface_kdt_vertices, opti_coord,
					surface_kdt.get(), medial_kdt.get(), surface_bvh.get());
				/*auto [v1, v2] =
					value<std::pair<SurfaceVertex, SurfaceVertex>>(surface, medial_axis_samples_closest_points, v);*/
				value<Vec3>(clusters, cluster_position, pv) = opti_coord;
				value<Scalar>(clusters, clusters_radius, pv) = radius;
				value<Vec3>(surface, cloest_point_color, v1) = value<Vec3>(clusters, cluster_color, pv);
				value<Vec3>(surface, cloest_point_color, v2) = value<Vec3>(clusters, cluster_color, pv);
				value<std::pair<SurfaceVertex, SurfaceVertex>>(clusters, cluster_cloest_sample, pv) = {v1, v2};

				break;
			}
			case 6: {
				std::vector<Point> surface_points;
				for (SurfaceVertex sv1 : cf.cluster_vertices)
				{
					Vec3 pos = value<Vec3>(surface, medial_axis_samples_position, sv1);
					surface_points.push_back(Point(pos.x(), pos.y(), pos.z()));
				}
				Min_sphere min_sphere(surface_points.begin(), surface_points.end());
				assert(min_sphere.is_valid());
				auto it = min_sphere.center_cartesian_begin();
				K::FT coord[3];
				for (size_t idx = 0; idx < 3 && it != min_sphere.center_cartesian_end(); idx++)
				{
					coord[idx] = *it;
					it++;
				}
				opti_coord = Vec3(coord[0], coord[1], coord[2]);
				K::FT radius = min_sphere.radius();
				auto [r, v1, v2] = geometry::move_point_to_medial_axis(
					surface, sample_position.get(), sample_normal.get(), surface_kdt_vertices, opti_coord,
					surface_kdt.get(), medial_kdt.get(), surface_bvh.get());
				/*auto [v1, v2] =
					value<std::pair<SurfaceVertex, SurfaceVertex>>(surface, medial_axis_samples_closest_points, v);*/
				value<Vec3>(surface, cloest_point_color, v1) = value<Vec3>(clusters, cluster_color, pv);
				value<Vec3>(surface, cloest_point_color, v2) = value<Vec3>(clusters, cluster_color, pv);
				value<Vec3>(clusters, cluster_position, pv) = opti_coord;
				value<Scalar>(clusters, clusters_radius, pv) = r;
				value<std::pair<SurfaceVertex, SurfaceVertex>>(clusters, cluster_cloest_sample, pv) = {v1, v2};

				break;
			}
			
			default:
				break;
			}
			
			
			return true;
		});
		//assign_cluster(surface, clusters);
		point_provider_->emit_connectivity_changed(clusters);
		surface_provider_->emit_attribute_changed(surface, cloest_point_color.get());
		point_provider_->emit_attribute_changed(clusters, cluster_position.get());
		point_provider_->emit_attribute_changed(clusters, clusters_radius.get());
	}
	
	void sphere_fitting(SURFACE& surface)
	{
		using ClusteringSQEM_Helper = modeling::ClusteringSQEM_Helper<SURFACE>;
		auto sample_position = get_or_add_attribute<Vec3, SurfaceVertex>(surface, "position");
		auto sample_normal = get_or_add_attribute<Vec3, SurfaceVertex>(surface, "normal");
		auto bvh_vertex_index = get_or_add_attribute<uint32, SurfaceVertex>(surface, "bvh_vertex_index");
		ClusteringSQEM_Helper helper(surface, sample_position.get(), sample_normal.get());
		std::vector<SurfaceVertex> vertex_vec;
		uint32 id = 0;
		foreach_cell(surface, [&](SurfaceVertex sv) {
			value<uint32>(surface, bvh_vertex_index, sv) = id++;
			vertex_vec.push_back(sv);
			surface_vertex_position_vector.push_back(value<Vec3>(surface, sample_position,sv));
			return true;
		});
		uint32 nb_faces = nb_cells<SurfaceFace>(surface);
		std::vector<SurfaceFace> bvh_faces;
		bvh_faces.reserve(nb_faces);
		std::vector<uint32> face_vertex_indices;
		face_vertex_indices.reserve(nb_faces * 3);

		foreach_cell(surface, [&](SurfaceFace f) -> bool {
			bvh_faces.push_back(f);
			foreach_incident_vertex(surface, f, [&](SurfaceVertex v) -> bool {
				face_vertex_indices.push_back(value<uint32>(surface, bvh_vertex_index, v));
				return true;
			});
			return true;
		});
		surface_bvh = std::make_unique<acc::BVHTree<uint32, Vec3>>(face_vertex_indices, surface_vertex_position_vector);
		surface_kdt = std::make_unique<acc::KDTree<3, uint32>>(surface_vertex_position_vector);
		Vec4 opt_sphere = Vec4(0,0,0,0);
		helper.optimal_sphere(vertex_vec, opt_sphere);
		POINT* sphere = point_provider_->add_mesh(point_provider_->mesh_name(surface) + "fitting sphere");
		auto sphere_position = get_or_add_attribute<Vec3, PointVertex>(*sphere, "position");
		auto sphere_radius = get_or_add_attribute<Scalar, PointVertex>(*sphere, "radius");
		auto sphere_color = get_or_add_attribute<Vec3, PointVertex>(*sphere, "color");
		/*PointVertex v = add_vertex(*sphere);
		value<Vec3>(*sphere, sphere_position, v) = Vec3(opt_sphere[0], opt_sphere[1], opt_sphere[2]);
		value<Scalar>(*sphere, sphere_radius, v) = opt_sphere[3];
		std::cout << opt_sphere << std::endl;
		value<Vec3>(*sphere, sphere_color, v) =
			Vec3(rand() / (double)RAND_MAX, rand() / (double)RAND_MAX, rand() / (double)RAND_MAX);*/

		PointVertex new_v = add_vertex(*sphere);
		Vec3 pos = Vec3(opt_sphere[0], opt_sphere[1], opt_sphere[2]);
		auto [new_radius, v1,v2] = geometry::move_point_to_medial_axis(
			surface, sample_position.get(),sample_normal.get(), vertex_vec, pos, surface_kdt.get(), medial_kdt.get(), surface_bvh.get());
		
		std::cout << new_radius << std::endl;
		value<Vec3>(*sphere, sphere_position, new_v) = pos;
		value<Scalar>(*sphere, sphere_radius, new_v) = new_radius;
		value<Vec3>(*sphere, sphere_color, new_v) =
			Vec3(rand() / (double)RAND_MAX, rand() / (double)RAND_MAX, rand() / (double)RAND_MAX);
		
		if (sphere_position)
			point_provider_->set_mesh_bb_vertex_position(*sphere, sphere_position);
		point_provider_->emit_connectivity_changed(*sphere);
	}
	
	void split_filtered_cluster(SURFACE& surface, POINT& clusters)
	{
			auto sample_position = get_or_add_attribute<Vec3, SurfaceVertex>(surface, "position");
			auto sample_normal = get_or_add_attribute<Vec3, SurfaceVertex>(surface, "normal");
			auto cluster_position = get_attribute<Vec3, PointVertex>(clusters, "position");
			auto clusters_infos = get_or_add_attribute<Cluster_Info, PointVertex>(clusters, "clusters_infos");
			auto cluster_color = get_or_add_attribute<Vec3, PointVertex>(clusters, "color");
			auto vertex_cluster_color = get_or_add_attribute<Vec3, SurfaceVertex>(surface, "cluster_color");
			auto clusters_radius = get_or_add_attribute<Scalar, PointVertex>(clusters, "clusters_radius");
			auto distance_to_cluster = get_or_add_attribute<Scalar, SurfaceVertex>(surface, "distance_to_cluster");
			auto medial_axis_samples_position_ =
				get_attribute<Vec3, SurfaceVertex>(surface, "medial_axis_samples_position");
			auto medial_axis_sample_radius_ =
				get_attribute<Scalar, SurfaceVertex>(surface, "medial_axis_samples_radius");
			auto medial_axis_samples_closest_points_ =
				get_attribute<std::pair<SurfaceVertex, SurfaceVertex>, SurfaceVertex>(
					surface, "medial_axis_samples_closest_points");
			auto meidal_axis_samples_weight =
				get_or_add_attribute<Scalar, SurfaceVertex>(surface, "medial_axis_samples_weight");
			auto medial_axis_samples_feature_value_ =
				get_or_add_attribute<Scalar, SurfaceVertex>(surface, "medial_axis_samples_feature_value");
			auto cluster_cloest_sample = get_or_add_attribute<std::pair<SurfaceVertex, SurfaceVertex>, PointVertex>(
				clusters, "cluster_cloest_sample");
			auto neighbours_set =
				get_or_add_attribute<std::unordered_map<uint32, PointVertex>, PointVertex>(clusters, "neighbours_set");

			MeshData<POINT>& md_cluster = point_provider_->mesh_data(clusters);
			bool split = false;
			CellMarkerStore<POINT, PointVertex> marker(clusters);
			std::cout << "split cluster ------------------------------------" << std::endl;
			//Test if the cluster center is outside
			foreach_cell(clusters, [&](PointVertex pv) { 
				auto [v1, v2] = value<std::pair<SurfaceVertex, SurfaceVertex>>(clusters, cluster_cloest_sample, pv);
				Vec3 n1 = value<Vec3>(surface, sample_normal, v1);
				Vec3 vec = value<Vec3>(surface, sample_position, v1) - value<Vec3>(clusters, cluster_position, pv);
				Scalar cosine = vec.dot(n1);
				if (cosine < 0)
				{
					Vec3 pos1 = value<Vec3>(surface, medial_axis_samples_position_, v1);
					Vec3 pos2 = value<Vec3>(surface, medial_axis_samples_position_, v2);
					auto [r1, v3, v4 ] = geometry::move_point_to_medial_axis(
						surface, sample_position.get(), sample_normal.get(), surface_kdt_vertices, pos1,
						surface_kdt.get(), medial_kdt.get(), surface_bvh.get());
					auto [r2, v5, v6] = geometry::move_point_to_medial_axis(
						surface, sample_position.get(), sample_normal.get(), surface_kdt_vertices, pos2,
						surface_kdt.get(), medial_kdt.get(), surface_bvh.get());
					value<Vec3>(clusters, cluster_position, pv) = pos1;
					value<Scalar>(clusters, clusters_radius, pv) = r1;
					value<std::pair<SurfaceVertex, SurfaceVertex>>(clusters, cluster_cloest_sample, pv) = {v3, v4};
					
					PointVertex new_cluster = add_vertex(clusters);
					value<Vec3>(clusters, cluster_position, new_cluster) = pos2;
					value<Scalar>(clusters, clusters_radius, new_cluster) = r2;
					value<std::pair<SurfaceVertex, SurfaceVertex>>(clusters, cluster_cloest_sample, new_cluster) = {v5,
																													v6};
					value<Vec3>(clusters, cluster_color, new_cluster) =
						Vec3(rand() / (double)RAND_MAX, rand() / (double)RAND_MAX, rand() / (double)RAND_MAX);
					
					split = true;
				
				}
				
				return true;
			});
			struct CompareScalar
			{
			bool operator()(const std::pair<PointVertex, Scalar>& a, const std::pair<PointVertex, Scalar>& b)
			{
				// Assuming higher variance is given higher priority
				return a.second < b.second;
			}
			};
			std::priority_queue<std::pair<PointVertex, Scalar>, 
				std::vector<std::pair<PointVertex, Scalar>>,
				CompareScalar> pq;
			if (!split)
			{
				// Test if the cluster can be split
				switch (clustering_mode)
				{
				case 0: 
				{ // variance
					Scalar max_error = -1e30;
					PointVertex candidate_cluster;
					foreach_cell(clusters, [&](PointVertex pv) {
						Cluster_Info& cf = value<Cluster_Info>(clusters, clusters_infos, pv);
						Scalar variance = geometry::surface_medial_distance_variance<SURFACE, POINT>(
							surface, clusters, pv, sample_position.get(), sample_normal.get(), cluster_position.get(),
							meidal_axis_samples_weight.get(), cf.cluster_vertices,
							value<Scalar>(clusters, clusters_radius, pv));
						Scalar error = variance ;
						std::cout << std::setiosflags(std::ios::fixed);
						std::cout << "variance: " << std::setprecision(9) 
							<< variance
								  << ", radius: " << 0.01 * value<Scalar>(clusters, clusters_radius, pv)
								  << ", error: " << error
								  << std::endl;
						if (error > split_variance_threshold_)
						{
							pq.push({pv, error});
						}
						cf.cluster_variance = error;
						return true;
					});
					
					
					break;
				}
				case 1: { // distance
					Scalar max_error = -1e30;
					PointVertex candidate_cluster;
					SurfaceVertex max_dist_vertex, global_max_dist_vertex;
					foreach_cell(clusters, [&](PointVertex pv) {
						Cluster_Info& cf = value<Cluster_Info>(clusters, clusters_infos, pv);
						Scalar error = 0;
						Scalar cluster_max_dist = -1e30;
						for (SurfaceVertex sv: cf.cluster_vertices)
						{
							Scalar dist =
								value<Scalar>(surface, distance_to_cluster, sv) * value<Scalar>(surface, meidal_axis_samples_weight, sv);
							error += dist;
						}
						error /= cf.cluster_vertices.size();
						std::cout << std::setiosflags(std::ios::fixed);
						std::cout << "distance: " << std::setprecision(9) << error
								  << ", radius: " << 0.01 * value<Scalar>(clusters, clusters_radius, pv)
								  << ", error: " << error << std::endl;
						if (error > split_distance_threshold_)
						{
							pq.push({pv, error});
						}
						return true;
					});
					break;
				}
				default:
					break;
				}
				while (pq.size() > 0)
				{
					auto& [candidate_cluster, error] = pq.top();
					pq.pop();
					if (marker.is_marked(candidate_cluster))
					continue;
					for (auto& [id, neighbour] :
						 value<std::unordered_map<uint32, PointVertex>>(clusters, neighbours_set, candidate_cluster))
					{
						marker.mark(neighbour);
					}
					Cluster_Info& cf = value<Cluster_Info>(clusters, clusters_infos, candidate_cluster);
					Scalar cluster_max_dist = -1e30;
					SurfaceVertex max_dist_vertex;
					for (SurfaceVertex sv : cf.cluster_vertices)
					{
						Scalar dist = value<Scalar>(surface, distance_to_cluster, sv) * 1000*
									  value<Scalar>(surface, meidal_axis_samples_weight, sv);
						if (dist > cluster_max_dist)
						{
							max_dist_vertex = sv;
						}
					}
					Vec3 pos = value<Vec3>(surface, medial_axis_samples_position_, max_dist_vertex);
					auto [r, v1, v2] = geometry::move_point_to_medial_axis(
						surface, sample_position.get(), sample_normal.get(), surface_kdt_vertices, pos,
						surface_kdt.get(), medial_kdt.get(), surface_bvh.get());

					PointVertex new_cluster = add_vertex(clusters);
					value<Vec3>(clusters, cluster_position, new_cluster) = pos;
					value<Scalar>(clusters, clusters_radius, new_cluster) = r;
					value<std::pair<SurfaceVertex, SurfaceVertex>>(clusters, cluster_cloest_sample, new_cluster) = {v1,
																													v2};
					value<Vec3>(clusters, cluster_color, new_cluster) =
						Vec3(rand() / (double)RAND_MAX, rand() / (double)RAND_MAX, rand() / (double)RAND_MAX);
				}
			}
			point_provider_->emit_connectivity_changed(clusters);
			point_provider_->emit_attribute_changed(clusters, cluster_position.get());
			point_provider_->emit_attribute_changed(clusters, clusters_radius.get());
			point_provider_->emit_attribute_changed(clusters, cluster_color.get());
			assign_cluster(surface, clusters);
			
	}

	void shrinking_balls_clustering_connectivity(SURFACE& surface, POINT& clusters)
	{
		NONMANIFOLD* mv = nonmanifold_provider_->add_mesh(
				std::to_string(nonmanifold_provider_->number_of_meshes()) + "_" + "skeleton");
		auto sample_position = get_or_add_attribute<Vec3, SurfaceVertex>(surface, "position");
		auto sample_normal = get_or_add_attribute<Vec3, SurfaceVertex>(surface, "normal");
		auto cluster_position = get_attribute<Vec3, PointVertex>(clusters, "position");
		auto clusters_infos = get_or_add_attribute<Cluster_Info, PointVertex>(clusters, "clusters_infos");
		auto cluster_color = get_or_add_attribute<Vec3, SurfaceVertex>(surface, "cluster_color");
		auto clusters_radius = get_or_add_attribute<Scalar, PointVertex>(clusters, "clusters_radius");
		auto cluster_cloest_sample = get_or_add_attribute<std::pair<SurfaceVertex, SurfaceVertex>, PointVertex>(
			clusters, "cluster_cloest_sample");

		auto medial_axis_samples_position_ =
			get_attribute<Vec3, SurfaceVertex>(surface, "medial_axis_samples_position");
		auto medial_axis_sample_radius_ =
			get_attribute<Scalar, SurfaceVertex>(surface, "medial_axis_samples_radius");
		
		auto vertex_cluster_info =
			get_or_add_attribute<std::pair<uint32, PointVertex>, SurfaceVertex>(surface, "vertex_cluster_info");
		auto neighbours_set = get_or_add_attribute<std::unordered_map<uint32, PointVertex>, PointVertex>(clusters, "neighbours_set");

		//Reindex the cluster
		uint32 index = 0;
		cluster_map.clear();
		foreach_cell(clusters, [&](PointVertex pv)
		{
			value<std::unordered_map<uint32, PointVertex>>(clusters, neighbours_set, pv).clear();
			cluster_map[index_of(clusters, pv)] = index;
			index++;
			return true;
		});

		foreach_cell(surface, [&](SurfaceVertex sv)
		{
				PointVertex pv = value<std::pair<uint32, PointVertex>>(surface, vertex_cluster_info, sv).second;
				foreach_adjacent_vertex_through_edge(surface, sv, [&](SurfaceVertex v) {
					if (value<std::pair<uint32, PointVertex>>(surface, vertex_cluster_info, v).first !=
						value<std::pair<uint32, PointVertex>>(surface, vertex_cluster_info, sv).first)
					{
						value<std::unordered_map<uint32, PointVertex>>(clusters, neighbours_set, pv)
							.insert(value<std::pair<uint32, PointVertex>>(surface, vertex_cluster_info, v));
					}
					return true;
				});
				return true;
		});

		/* foreach_cell(clusters, [&](PointVertex pv) {
			CellMarker<SURFACE, SurfaceVertex> visited(surface);
			std::queue<SurfaceVertex> queue;
			auto [v1, v2] = value<std::pair<SurfaceVertex, SurfaceVertex>>(clusters, cluster_cloest_sample, pv);
			queue.push(v1);
			queue.push(v2);
			visited.mark(v1);
			visited.mark(v2);
			while (queue.size() > 0)
			{
				SurfaceVertex current = queue.front();
				queue.pop();
				foreach_adjacent_vertex_through_edge(surface, current, [&](SurfaceVertex v){
					
					if (visited.is_marked(v))
					{
						return true;
					}
					visited.mark(v);
					if (value<std::pair<uint32, PointVertex>>(surface, vertex_cluster_info, v).first ==
						value<std::pair<uint32, PointVertex>>(surface, vertex_cluster_info, current).first)
					{
						queue.push(v);
					}
					else
					{
						value<std::unordered_map<uint32, PointVertex>>(clusters, neighbours_set, pv)
							.insert(value<std::pair<uint32, PointVertex>>(surface, vertex_cluster_info, v));
					}
					return true;
				});
				
			}
			return true;
		});*/
		
		
		cgogn::io::IncidenceGraphImportData skeleton_data;
		std::unordered_map<std::pair<uint32, uint32>, uint32, edge_hash, edge_equal> edge_indices;
		skeleton_data.vertex_position_ = std::vector<Vec3>(nb_cells<PointVertex>(clusters), Vec3(0, 0, 0));
		uint32 edge_count = 0;
		//Add vertex and edge
		foreach_cell(clusters, [&](PointVertex pv) {
			uint32 cluster_index_1 = cluster_map[index_of(clusters, pv)];
			Vec3 pos = value<Vec3>(clusters, cluster_position, pv);
			skeleton_data.vertex_position_[cluster_index_1] = pos;
			std::unordered_map<uint32, PointVertex>& neighbours_pv =
				value<std::unordered_map<uint32, PointVertex>>(clusters, neighbours_set, pv);
			for (auto& vertex_neighbours_1 : neighbours_pv)
			{
				uint32 cluster_index_2 = cluster_map[index_of(clusters, vertex_neighbours_1.second)];
				if (edge_indices.find({cluster_index_1, cluster_index_2}) == edge_indices.end())
				{
					edge_indices[{cluster_index_1, cluster_index_2}] = edge_count;
					skeleton_data.edges_vertex_indices_.push_back(cluster_index_1);
					skeleton_data.edges_vertex_indices_.push_back(cluster_index_2);
					edge_count++;
				}
			}
			return true;
		});
		//Add face
		foreach_cell(clusters, [&](PointVertex pv) {
			uint32 cluster_index_1 = cluster_map[index_of(clusters, pv)];
			std::unordered_map<uint32, PointVertex>& neighbours_pv =
				value<std::unordered_map<uint32, PointVertex>>(clusters, neighbours_set, pv);
			for (auto& vertex_neighbours_1 : neighbours_pv)
			{
				uint32 cluster_index_2 = cluster_map[index_of(clusters, vertex_neighbours_1.second)]; 
				std::unordered_map<uint32, PointVertex>& neighbours_1 =
					value<std::unordered_map<uint32, PointVertex>>(clusters, neighbours_set, vertex_neighbours_1.second);
				for (auto& vertex_neighbours_2 : neighbours_1)
				{
					if (neighbours_pv.find(vertex_neighbours_2.first) != neighbours_pv.end())
					{
						uint32 cluster_index_3 = cluster_map[index_of(clusters, vertex_neighbours_2.second)]; 
						skeleton_data.faces_nb_edges_.push_back(3);
						skeleton_data.faces_edge_indices_.push_back(
							edge_indices[{cluster_index_1, cluster_index_2}]);
						skeleton_data.faces_edge_indices_.push_back(
							edge_indices[{cluster_index_1, cluster_index_3}]);
						skeleton_data.faces_edge_indices_.push_back(
							edge_indices[{cluster_index_2, cluster_index_3}]);
					}
				}
			}
			return true;
		});
		uint32 skeleton_nb_vertices = skeleton_data.vertex_position_.size();
		uint32 skeleton_nb_edges = skeleton_data.edges_vertex_indices_.size() / 2;
		uint32 skeleton_nb_faces = skeleton_data.faces_nb_edges_.size();

		skeleton_data.reserve(skeleton_nb_vertices, skeleton_nb_edges, skeleton_nb_faces);

		import_incidence_graph_data(*mv, skeleton_data);
		std::shared_ptr<NonManifoldAttribute<Vec3>> mv_vertex_position =
			get_attribute<Vec3, NonManifoldVertex>(*mv, "position");
		if (mv_vertex_position)
			nonmanifold_provider_->set_mesh_bb_vertex_position(*mv, mv_vertex_position);
		
		nonmanifold_provider_->emit_connectivity_changed(*mv);
		std::cout << "finish creating skeleton"<< std::endl;
		remove_attribute<PointVertex>(clusters, neighbours_set.get());
	}


 protected:
	void init() override
	{
		point_provider_ = static_cast<ui::MeshProvider<POINT>*>(
			app_.module("MeshProvider (" + std::string{mesh_traits<POINT>::name} + ")"));
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
		
		imgui_mesh_selector(point_provider_, selected_clusters, "Clusters", [&](POINT& m) {
			selected_clusters = &m;
			point_provider_->mesh_data(m).outlined_until_ = App::frame_time_ + 1.0;
		});
	

		if (selected_surface_mesh_){
		
			/*if (ImGui::Button("fitting sphere"))
				sphere_fitting(*selected_surface_mesh_);*/
			/*if (ImGui::Button("Sample surface"))
				construct_surface_sample(*selected_surface_mesh_);*/
			if (ImGui::Button("Sample medial axis"))
				sample_medial_axis(*selected_surface_mesh_);
			if (ImGui::SliderFloat("Min radius (log)", &radius_threshold_, min_radius_, max_radius_, "%.4f"))
				filter_medial_axis_samples(*selected_surface_mesh_);
			if (ImGui::SliderFloat("Min angle (log)", &angle_threshold_, min_angle_, max_angle_, "%.4f"))
				filter_medial_axis_samples(*selected_surface_mesh_);
			if (ImGui::Button("Filter medial axis samples"))
				filter_medial_axis_samples(*selected_surface_mesh_);
			if (ImGui::Button("Twin point"))
				Twin_medial_axis_samples(*selected_surface_mesh_);
			if (ImGui::Button("Initialise Cluster"))
				iniialise_cluster2(*selected_surface_mesh_);
			if (selected_clusters)
			{
				ImGui::SliderInt("Update time", &(int)update_times, 1, 100);
				if (ImGui::Button("Average of medial points"))
				{
					clustering_mode = 0;
					for (uint32 i = 0; i < update_times; i++)
						update_filtered_cluster(*selected_surface_mesh_, *selected_clusters);
				}
				if (ImGui::Button("Average of sample points"))
				{
					clustering_mode = 1;
					for (uint32 i = 0; i < update_times; i++)
						update_filtered_cluster(*selected_surface_mesh_, *selected_clusters);
				}
				if (ImGui::Button("Sphere fitting"))
				{
					clustering_mode = 2;
					for (uint32 i = 0; i < update_times; i++)
						update_filtered_cluster(*selected_surface_mesh_, *selected_clusters);
				
				}
				if (ImGui::Button("Average of medial points weighted with curvature"))
				{
					clustering_mode = 3;
					for (uint32 i = 0; i < update_times; i++)
						update_filtered_cluster(*selected_surface_mesh_, *selected_clusters);
				}
				if (ImGui::Button("Average of medial points weighted with feature value"))
				{
					clustering_mode = 4;
					for (uint32 i = 0; i < update_times; i++)
						update_filtered_cluster(*selected_surface_mesh_, *selected_clusters);
				}
				
				if (ImGui::Button("Update Cluster by minimum enclosing ball of surface points"))
				{
					clustering_mode = 5;
					for (uint32 i = 0; i < update_times; i++)
						update_filtered_cluster(*selected_surface_mesh_, *selected_clusters);
				}
				if (ImGui::Button("Update Cluster by minimum enclosing ball of medial points"))
				{
					clustering_mode = 6;
					for (uint32 i = 0; i < update_times; i++)
						update_filtered_cluster(*selected_surface_mesh_, *selected_clusters);
				}
				
				
				ImGui::DragFloat("Split vairance threshold", &split_variance_threshold_, 0.00000001f, 0.0f, 0.000001f,
								 "%.6f");
				ImGui::DragFloat("Split distance threshold", &split_distance_threshold_, 0.00000001f, 0.0f, 0.000001f,
								 "%.6f");
				if (ImGui::Button("Split Cluster by variance"))
				{
					clustering_mode = 0;
					split_filtered_cluster(*selected_surface_mesh_,*selected_clusters);
				}
				if (ImGui::Button("Split Cluster by distance"))
				{
					clustering_mode = 1;
					split_filtered_cluster(*selected_surface_mesh_, *selected_clusters);
				}
				if (ImGui::Button("Construct skeleton"))
					shrinking_balls_clustering_connectivity(*selected_surface_mesh_, *selected_clusters);
			}
			ImGui::DragFloat("Dilation factor", &dilation_factor, 0.0001f, 0.0f, 1.0f, "%.6f");
			if (ImGui::Button("Shrinking Ball Coverage Axis"))
			{
				solution = shrinking_balls_coverage_axis(*selected_surface_mesh_, dilation_factor);
			}
			if (solution.col_value.size() > 0)
			{
				if (ImGui::Button("Shrinking balls PD"))
					shrinking_balls_coverage_axis_PD(*selected_surface_mesh_);
			}
			if (ImGui::Button("Compute delaunay"))
			{
				load_model_in_cgal(*selected_surface_mesh_, csm);
				Tree tree(faces(csm).first, faces(csm).second, csm);
				tree.accelerate_distance_queries();
				tri_ = compute_delaunay_tredrahedron(*selected_surface_mesh_, csm, tree);
			}

			imgui_mesh_selector(point_provider_, selected_candidates_, "Candidates", [&](POINT& m) {
				selected_candidates_ = &m;
				point_provider_->mesh_data(m).outlined_until_ = App::frame_time_ + 1.0;
			});
			if (selected_candidates_)
			{
				ImGui::Checkbox("Preseve only poles", &pole_filtering_);
				ImGui::Checkbox("Filter by Angle", &angle_filtering_);
				if (angle_filtering_)
				{
					ImGui::DragFloat("angle", &(float)angle_threshold_, 0.01, 0.0, M_PI, "%.2f");
				}

				ImGui::Checkbox("Filter by circumradius", &circumradius_filtering_);
				if (circumradius_filtering_)
				{
					ImGui::DragFloat("circumradius", &(float)radius_threshold_, (max_radius_ - min_radius_) / 100,
									 min_radius_,
									 max_radius_, "%.3f");
				}
				ImGui::Checkbox("Filter by distance", &distance_filtering_);
				if (distance_filtering_)
				{
					ImGui::DragFloat("distance", &distance_threshold_, 1e-2, 0, 1, "%.4f");
				}

				if (ImGui::Button("Generate candidates"))
				{
					construct_candidates_points(*selected_surface_mesh_, tri_);
				}
				if (ImGui::Button("Power shape"))
				{
					compute_power_shape(*selected_surface_mesh_);
				}
				if (ImGui::Button("Original Medial Axis"))
					compute_original_power_diagram(*selected_surface_mesh_);
				imgui_mesh_selector(nonmanifold_provider_, selected_medial_axis_, "Medial_axis", [&](NONMANIFOLD& nm) {
					selected_medial_axis_ = &nm;
					nonmanifold_provider_->mesh_data(nm).outlined_until_ = App::frame_time_ + 1.0;
				});
				if (selected_medial_axis_)
				{
					if (ImGui::Button("Compute stability ratio"))
						compute_stability_ratio(*selected_medial_axis_);
					static int32 number_vertex_remain = 1;
					static float k = 1e-5;
					ImGui::DragInt("Vertices to delete", &number_vertex_remain, 1, 0,
								   nb_cells<NonManifoldVertex>(*selected_medial_axis_));
					ImGui::DragFloat("K", &k, 1e-5, 0.0f, 1.0f, "%.5f");
					if (ImGui::Button("QMAT"))
					{
						collapse_non_manifold_using_QMat(*selected_medial_axis_, number_vertex_remain, k);
					}
				}
				if (selected_candidates_)
				{
					ImGui::DragFloat("Dilation factor", &dilation_factor, 0.001f, 0.0f, 1.0f, "%.4f");
					if (ImGui::Button("Coverage Axis"))
					{
						solution = point_selection_by_coverage_axis(*selected_surface_mesh_, *selected_candidates_,
																	dilation_factor);
					}
					if (solution.col_value.size() > 0)
					{
						if (ImGui::Button("Collpase"))
							coverage_axis_collapse(*selected_surface_mesh_, *selected_candidates_, solution);
						if (ImGui::Button("PD"))
							coverage_axis_PD(*selected_surface_mesh_, *selected_candidates_, solution, dilation_factor);
						if (ImGui::Button("RVD"))
							coverage_axis_RVD(*selected_surface_mesh_, *selected_candidates_, solution);
					}
				}
			}
		}
	}

private:
	uint32 clustering_mode = 0;//0: least square, 1: QEM
	
	POINT* selected_candidates_ = nullptr;
	POINT* surface_sample_ = nullptr;
	POINT* selected_clusters = nullptr;
	SURFACE* selected_surface_mesh_ = nullptr;
	NONMANIFOLD* selected_medial_axis_ = nullptr;
	MeshProvider<POINT>* point_provider_;
	MeshProvider<SURFACE>* surface_provider_;
	MeshProvider<NONMANIFOLD>* nonmanifold_provider_;
	HighsSolution solution;
	Delaunay tri_;
	Tree* tree;
	Point_inside* inside_tester;
	Cgal_Surface_mesh csm;
	std::unordered_map<uint32, uint32> cluster_map;
	std::vector<Vec3> colors;
	float dilation_factor = 0.1f;
	bool angle_filtering_ = true;
	bool circumradius_filtering_ = true;
	bool distance_filtering_ = true;
	bool pole_filtering_ = true;
	float distance_threshold_ = 0.001;
	float angle_threshold_ = 1.9;
	float radius_threshold_ = 0.030;
	float split_variance_threshold_ = 0.000001;
	float split_distance_threshold_ = 0.003;
	uint32 update_times = 5;
	double min_radius_ = std::numeric_limits<double>::max();
	double max_radius_ = std::numeric_limits<double>::min();
	double min_angle_ = std::numeric_limits<double>::max();
	double max_angle_ = std::numeric_limits<double>::min();
	Eigen::MatrixXd dis_matrix;
	std::vector<SurfaceVertex> surface_kdt_vertices;
	std::vector<Vec3> surface_vertex_position_vector;
	std::vector<Vec3> medial_vertex_position_vector;
	std::unique_ptr<acc::KDTree<3, uint32>> surface_kdt;
	std::unique_ptr<acc::BVHTree<uint32, Vec3>> surface_bvh;
	std::unique_ptr<acc::KDTree<3, uint32>> medial_kdt;
	std::unique_ptr<CellMarker<SURFACE, SurfaceVertex>> unsplittable;
	std::unique_ptr<modeling::ClusteringSQEM_Helper<SURFACE>> sqem_helper;
	};

} // namespace ui

} // namespace cgogn

#endif // CGOGN_MODULE_POWER_SHAPE_H_
