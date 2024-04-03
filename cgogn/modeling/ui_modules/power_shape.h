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
#include <libacc/bvh_trees_spheres.h>

#include <cgogn/io/point/point_import.h>
#include <cgogn/geometry/algos/medial_axis.h>
#include <cgogn/modeling/algos/decimation/SQEM_helper.h>
#include <cgogn/modeling/algos/decimation/QEM_helper.h>
#include <cgogn/rendering/skelshape.h>
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
#include <GLFW/glfw3.h>
#include <Highs.h>

#include <iomanip>
#include <limits>

namespace cgogn
{

namespace ui
{
enum UpdateMethod
{
	MEAN_MEDIAL_SAMPLE,
	MEAN_SURFACE_SAMPLE,
	MEAN_MEDIAL_SAMPLE_WEIGHTED_RAYON,
	MEAN_MEDIAL_SAMPLE_WEIGHTED_CURVATURE,
	SPHERE_FITTING,
	ENCLOSING_SPHERE_OF_MEDIAL_SAMPLE,
	ENCLOSING_SPHERE_OF_SURFACE_SAMPLE,
	MASS_SPRING,
	SQEM,
	OPTIMAL_SPHERE_UNDER_CONSTRAINT
};

enum SplitMethod
{
	DEVIATION,
	DISTANCE_POINT_SPHERE,
	DISTANCE_TO_MEDIAL_AXIS,
	HAUSDORFF_DISTANCE,
};

enum InitMethod
{
	CONSTANT,
	FACTOR
};
enum CandidateGenerationMethod
{
	RANDOM,
	SHRINKING_BALL,
	DELAUNAY

};


template <typename POINT, typename SURFACE, typename NONMANIFOLD>
class PowerShape : public ViewModule
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
	template <typename T>
	using PointAttribute = typename mesh_traits<POINT>::template Attribute<T>;

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
	PowerShape(const App& app) : ViewModule(app, "PowerShape"), tree(nullptr), inside_tester(nullptr)
	{
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

	std::pair<Vec3, Scalar> sphere_fitting_algo(SURFACE& surface, std::vector<SurfaceVertex> clusters_surface_vertices_)
	{
		ClusterAxisParameter& p = cluster_axis_parameters_[&surface];

		Eigen::MatrixXd A(clusters_surface_vertices_.size(), 4);
		Eigen::VectorXd b(clusters_surface_vertices_.size());
		uint32 idx = 0;
		for (SurfaceVertex v : clusters_surface_vertices_)
		{
			const Vec3& pos = value<Vec3>(surface, p.surface_vertex_position_, v);
			A.row(idx) = Eigen::Vector4d(-2.0 * pos[0], -2.0 * pos[1], -2.0 * pos[2], 1.0);
			b(idx) = -(pos[0] * pos[0]) - (pos[1] * pos[1]) - (pos[2] * pos[2]);
			++idx;
		};
		Eigen::LDLT<Eigen::MatrixXd> solver(A.transpose() * A);
		Eigen::MatrixXd s1 = solver.solve(A.transpose() * b);
		s1(3) = std::sqrt(s1(0) * s1(0) + s1(1) * s1(1) + s1(2) * s1(2) - s1(3));

		Vec3 s1c = Vec3(s1(0), s1(1), s1(2));
		Scalar s1r = s1(3);

		Eigen::MatrixXd J(clusters_surface_vertices_.size(), 4);
		Eigen::VectorXd r(clusters_surface_vertices_.size());
		Eigen::VectorXd s2(4);
		s2 << s1(0), s1(1), s1(2), s1(3);
		for (uint32 i = 0; i < 5; ++i) // TODO: check number of iterations
		{
			idx = 0;
			for (SurfaceVertex v : clusters_surface_vertices_)
			{
				const Vec3& pos = value<Vec3>(surface, p.surface_vertex_position_, v);
				Vec3 d = pos - Vec3(s2(0), s2(1), s2(2));
				Scalar l = d.norm();
				J.row(idx) = Eigen::Vector4d(-(d[0] / l), -(d[1] / l), -(d[2] / l), -1.0);
				r(idx) = -(l - s2(3));
				++idx;
			}
			Eigen::LDLT<Eigen::MatrixXd> solver(J.transpose() * J);
			s2 += solver.solve(J.transpose() * r);
		}

		Vec3 s2c = Vec3(s2(0), s2(1), s2(2));
		Scalar s2r = s2(3);
		return {s2c, s2r};
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
		auto surface_vertex_position_ = get_attribute<Vec3, Vertex>(surface, "position");
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
			Vec3 pos = value<Vec3>(surface, surface_vertex_position_, v);
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
		load_model_in_cgal(surface,csm);
		auto inner_position = get_attribute<Vec3, PointVertex>(mv, "position");
		auto surface_vertex_position_ = get_attribute<Vec3, SurfaceVertex>(surface, "position");
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
			Vec3 pos = value<Vec3>(surface, surface_vertex_position_, v);
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
		auto surface_vertex_position_ = get_attribute<Vec3, Vertex>(surface, "position");
	
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
			Vec3 pos = value<Vec3>(surface, surface_vertex_position_, v);
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
		auto surface_vertex_position_ = get_attribute<Vec3, SurfaceVertex>(surface, "position");
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
				if (inside_sphere(value<Vec3>(surface, surface_vertex_position_, v), inner_points[idx], weights[idx]))
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

	void sample_medial_axis(SURFACE& s)
	{
		auto vertex_position = get_or_add_attribute<Vec3, SurfaceVertex>(s, "position");
		auto vertex_normal = get_or_add_attribute<Vec3, SurfaceVertex>(s, "normal");
		auto medial_axis_position_ =
			get_or_add_attribute<Vec3, SurfaceVertex>(s, "medial_axis_position_");
		auto medial_axis_samples_radius_ = get_or_add_attribute<Scalar, SurfaceVertex>(s, "medial_axis_samples_radius");
		auto medial_axis_closest_points_ =
			get_or_add_attribute<std::pair<SurfaceVertex, SurfaceVertex>, SurfaceVertex>(
				s, "medial_axis_samples_closest_points");
		auto medial_axis_samples_angle_ = get_or_add_attribute<Scalar, SurfaceVertex>(s, "medial_axis_samples_angle");
		auto medial_axis_samples_weight_ = get_or_add_attribute<Scalar, SurfaceVertex>(s, "medial_axis_samples_weight");
		auto kmax = get_or_add_attribute<Scalar, SurfaceVertex>(s, "kmax");
		geometry::shrinking_ball_centers<SURFACE>(
			s, vertex_position.get(), vertex_normal.get(),
										 medial_axis_position_.get(),
										 medial_axis_samples_radius_.get(),
										 medial_axis_closest_points_.get());
		MeshData<SURFACE>& md = surface_provider_->mesh_data(s);
		auto filtered_medial_axis_samples_set_ = &md.template get_or_add_cells_set<SurfaceVertex>("filtered_medial_axis_samples");
		filtered_medial_axis_samples_set_->select_if([&](SurfaceVertex v) { return true; });
		foreach_cell(s, [&](SurfaceVertex v) -> bool {
			const Vec3& c = value<Vec3>(s, medial_axis_position_, v);
			const auto& [v1, v2] =
				value<std::pair<SurfaceVertex, SurfaceVertex>>(s, medial_axis_closest_points_, v);
			const auto c1 = value<Vec3>(s, vertex_position, v1);
			const auto c2 = value<Vec3>(s, vertex_position, v2);
			const Scalar r = value<Scalar>(s, medial_axis_samples_radius_, v);
			max_radius_ = std::max(max_radius_, r);
			min_radius_ = std::min(min_radius_, r);
			const Scalar angle = geometry::angle(c1 - c, c2 - c);
			max_angle_ = std::max(max_angle_, angle);
			min_angle_ = std::min(min_angle_, angle);
			value<Scalar>(s, medial_axis_samples_angle_, v) = angle;
			value<Scalar>(s, medial_axis_samples_weight_, v) = value<Scalar>(s, kmax, v);
		
			return true;
		});
		normalise_scalar(s, medial_axis_samples_weight_);
		surface_provider_->emit_attribute_changed(s, medial_axis_position_.get());
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
		auto medial_axis_position_ = get_attribute<Vec3, SurfaceVertex>(s, "medial_axis_position_");
		auto medial_axis_samples_radius_ = get_attribute<Scalar, SurfaceVertex>(s, "medial_axis_samples_radius");
		auto medial_axis_closest_points_ =
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
		auto medial_axis_position_ = get_attribute<Vec3, SurfaceVertex>(s, "medial_axis_position_");
		auto medial_axis_samples_radius_ = get_attribute<Scalar, SurfaceVertex>(s, "medial_axis_samples_radius");
		auto medial_axis_closest_points_ =
			get_attribute<std::pair<SurfaceVertex, SurfaceVertex>, SurfaceVertex>(s,
																				  "medial_axis_samples_closest_points");
		auto medial_axis_samples_angle_ = get_or_add_attribute<Scalar, SurfaceVertex>(s, "medial_axis_samples_angle");
		foreach_cell(s, [&](SurfaceVertex v) -> bool {
			auto [v1, v2] = value<std::pair<SurfaceVertex, SurfaceVertex>>(s, medial_axis_closest_points_, v);
			auto [v3, v4] = value<std::pair<SurfaceVertex, SurfaceVertex>>(s, medial_axis_closest_points_, v2);
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
		auto medial_axis_position_ = get_attribute<Vec3, SurfaceVertex>(surface, "medial_axis_position_");
		auto medial_axis_samples_radius = get_attribute<Scalar, SurfaceVertex>(surface, "medial_axis_samples_radius");
		auto surface_vertex_position_ = get_attribute<Vec3, SurfaceVertex>(surface, "position");
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
			sample_points.push_back(value<Vec3>(surface, surface_vertex_position_, sv));
			inner_points.push_back(value<Vec3>(surface, medial_axis_position_, sv));
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
		auto medial_axis_position_ = get_attribute<Vec3, SurfaceVertex>(surface, "medial_axis_position_");
		auto medial_axis_samples_radius = get_attribute<Scalar, SurfaceVertex>(surface, "medial_axis_samples_radius");
		auto surface_vertex_position_ = get_attribute<Vec3, SurfaceVertex>(surface, "position");
		auto selected_medial_axis_samples_set =
			&surface_provider_->mesh_data(surface).template get_or_add_cells_set<SurfaceVertex>(
				"selected_medial_axis_samples");
		auto filter_medial_axis_samples =
			&surface_provider_->mesh_data(surface).template get_or_add_cells_set<SurfaceVertex>(
				"filtered_medial_axis_samples");
		CellMarker<SURFACE, SurfaceVertex> marker(surface);
		selected_medial_axis_samples_set->foreach_cell([&](SurfaceVertex svc) {
			Vec3 color = Vec3(rand()/double(RAND_MAX), rand()/double(RAND_MAX), rand()/double(RAND_MAX));
			foreach_cell(surface, [&](SurfaceVertex svs) {
				if (marker.is_marked(svs))
					return true;
				if (inside_sphere(value<Vec3>(surface, surface_vertex_position_, svs),
													  value<Vec3>(surface, medial_axis_position_, svc),
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
		auto medial_axis_position_ = get_attribute<Vec3, SurfaceVertex>(surface, "medial_axis_position_");
		auto medial_axis_samples_radius = get_attribute<Scalar, SurfaceVertex>(surface, "medial_axis_samples_radius");
		auto surface_vertex_position_ = get_attribute<Vec3, SurfaceVertex>(surface, "position");
		auto selected_medial_axis_samples_set =
			&surface_provider_->mesh_data(surface).template get_or_add_cells_set<SurfaceVertex>(
				"selected_medial_axis_samples");
		Regular power_shape;
		std::cout << "collect selected points info" << std::endl;
		selected_medial_axis_samples_set->foreach_cell([&](SurfaceVertex sv) {
			Vec3 pos = value<Vec3>(surface, medial_axis_position_, sv);
			power_shape.insert(
				Weight_Point(Point(pos[0], pos[1], pos[2]),
											(value<double>(surface, medial_axis_samples_radius, sv)+0.02) *
												(value<double>(surface, medial_axis_samples_radius, sv)+0.02)));

			return true;
		});
		std::cout << "collect surface samples info " << std::endl;
		foreach_cell(surface, [&](SurfaceVertex v) {
			Vec3 pos = value<Vec3>(surface, surface_vertex_position_, v);
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

	
	struct CoverageAxisParameter
	{
		bool initialized_ = false;
		bool candidates_valid = false;
		Eigen::MatrixXi coverage_matrix;
		SURFACE* surface_;
		std::shared_ptr<SurfaceAttribute<Vec3>> surface_vertex_position_ = nullptr;
		std::shared_ptr<SurfaceAttribute<Vec3>> surface_vertex_normal_ = nullptr;
		std::shared_ptr<SurfaceAttribute<Vec3>> medial_axis_position_ = nullptr;
		std::shared_ptr<SurfaceAttribute<Scalar>> medial_axis_radius_ = nullptr;
		std::shared_ptr<SurfaceAttribute<std::pair<SurfaceVertex, SurfaceVertex>>> medial_axis_closest_points_ =
			nullptr;
		std::shared_ptr<SurfaceAttribute<Scalar>> medial_axis_angle_ = nullptr;


		POINT* candidate_points_ = nullptr;
		
		std::shared_ptr<PointAttribute<Vec3>> candidate_points_position_ = nullptr;
		std::shared_ptr<PointAttribute<Scalar>> candidate_points_radius_ = nullptr;
		std::shared_ptr<PointAttribute<Scalar>> candidate_points_score_ = nullptr;
		std::shared_ptr<PointAttribute<Scalar>> candidate_points_coverage_score_ = nullptr;
		std::shared_ptr<PointAttribute<Scalar>> candidate_points_uniformity_score_ = nullptr;
		std::shared_ptr<PointAttribute<NonManifoldVertex>> candidate_points_associated_vertex_ = nullptr;

		POINT* selected_points_ = nullptr;
		std::shared_ptr<PointAttribute<Vec3>> selected_points_position_ = nullptr;
		std::shared_ptr<PointAttribute<Scalar>> selected_points_radius_ = nullptr;
		std::shared_ptr<PointAttribute<NonManifoldVertex>> selected_points_associated_vertex_ = nullptr;

		NONMANIFOLD* voronoi_ = nullptr;
		std::shared_ptr<NonManifoldAttribute<Vec3>> voronoi_position_ = nullptr;
		std::shared_ptr<NonManifoldAttribute<Scalar>> voronoi_radius_ = nullptr;
		std::shared_ptr<NonManifoldAttribute<Scalar>> voronoi_stability_ratio_ = nullptr;
		std::shared_ptr<NonManifoldAttribute<Vec3>> voronoi_stability_color_ = nullptr;
		std::shared_ptr<NonManifoldAttribute<Vec4>> voronoi_sphere_info_ = nullptr;
		std::shared_ptr<NonManifoldAttribute<bool>> voronoi_fixed_vertices_ = nullptr;

		std::vector<NonManifoldVertex> voronoi_kdt_vertices_;
		std::shared_ptr<acc::KDTree<3, uint32>> voronoi_kdt_ = nullptr;

		float dilation_factor = 0.02;
		CandidateGenerationMethod candidate_generation_method = SHRINKING_BALL;

		uint32 max_selected_number = 200;
		std::vector<PointVertex> selected_inner_points;
		std::vector<PointVertex> candidate_inner_points;

		Cgal_Surface_mesh csm_;
		std::shared_ptr<Tree> tree_ = nullptr;
		std::shared_ptr<Point_inside> inside_tester_ = nullptr;
		Delaunay tri;
		Scalar min_radius_ = std::numeric_limits<Scalar>::max();
		Scalar max_radius_ = std::numeric_limits<Scalar>::min();
		Scalar min_angle_ = std::numeric_limits<Scalar>::max();
		Scalar max_angle_ = std::numeric_limits<Scalar>::min();
	};

	

	void compute_initial_non_manifold(CoverageAxisParameter& p)
	{
		foreach_cell(*p.surface_, [&](SurfaceVertex sv) {
			Delaunay_Vertex_handle vh =
				p.tri.insert(Point(value<Vec3>(*p.surface_, p.surface_vertex_position_, sv)[0],
								   value<Vec3>(*p.surface_, p.surface_vertex_position_, sv)[1],
								   value<Vec3>(*p.surface_, p.surface_vertex_position_, sv)[2]));
			return true;
		});
		// auto start_timer = std::chrono::high_resolution_clock::now();

		// Construct delauney
		for (auto cit = p.tri.finite_cells_begin(); cit != p.tri.finite_cells_end(); ++cit)
		{
			cit->info().id = -1;
			cit->info().centroid = CGAL::circumcenter(p.tri.tetrahedron(cit));
			cit->info().radius2 = CGAL::squared_distance(cit->info().centroid, cit->vertex(0)->point());

			if (inside(*(p.inside_tester_.get()), cit->info().centroid))
			{
				cit->info().inside = true;
			}
			else
			{
				cit->info().inside = false;
				
			}
		}
		cgogn::io::IncidenceGraphImportData Initial_non_manifold;
		std::unordered_map<std::pair<uint32, uint32>, size_t, edge_hash, edge_equal> edge_indices;
		std::vector<Scalar> sphere_radius;
		std::vector<Vec3> sphere_center;
		int count = 0;

		// Add Vertex
		for (auto cit = p.tri.finite_cells_begin(); cit != p.tri.finite_cells_end(); ++cit)
		{
			if (cit->info().inside)
			{
				cit->info().id = count;
				Point centroid = cit->info().centroid;
				Scalar radius = std::sqrt(cit->info().radius2);
				Initial_non_manifold.vertex_position_.emplace_back(centroid[0], centroid[1], centroid[2]);
				sphere_radius.push_back(radius);
				sphere_center.push_back(Vec3(centroid[0], centroid[1], centroid[2]));
				count++;
			}
		}
		// Add edges
		for (auto fit = p.tri.finite_facets_begin(); fit != p.tri.finite_facets_end(); ++fit)
		{
			if (fit->first->info().inside && p.tri.mirror_facet(*fit).first->info().inside)
			{
				uint32 v1_ind = fit->first->info().id;
				uint32 v2_ind = p.tri.mirror_facet(*fit).first->info().id;
				Initial_non_manifold.edges_vertex_indices_.push_back(v1_ind);
				Initial_non_manifold.edges_vertex_indices_.push_back(v2_ind);
				edge_indices.insert({{v1_ind, v2_ind}, edge_indices.size()});
			}
		}
		bool all_finite_inside;

		std::vector<Delaunay_Cell_handle> incells;
		for (auto eit = p.tri.finite_edges_begin(); eit != p.tri.finite_edges_end(); ++eit)
		{
			all_finite_inside = true;
			incells.clear();
			Delaunay_Cell_circulator cc = p.tri.incident_cells(*eit);
			do
			{
				if (p.tri.is_infinite(cc))
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
			} while (cc != p.tri.incident_cells(*eit));
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

		import_incidence_graph_data(*p.voronoi_, Initial_non_manifold);
		
		for (size_t idx = 0; idx < sphere_center.size(); idx++)
		{
			(*p.voronoi_radius_)[idx] = sphere_radius[idx];
			(*p.voronoi_sphere_info_)[idx] =
				Vec4(sphere_center[idx].x(), sphere_center[idx].y(), sphere_center[idx].z(), sphere_radius[idx]);
		}
		if (p.voronoi_position_)
			nonmanifold_provider_->set_mesh_bb_vertex_position(*p.voronoi_, p.voronoi_position_);
		nonmanifold_provider_->emit_connectivity_changed(*p.voronoi_);
		nonmanifold_provider_->emit_attribute_changed(*p.voronoi_, p.voronoi_position_.get());
		nonmanifold_provider_->emit_attribute_changed(*p.voronoi_, p.voronoi_radius_.get());
		nonmanifold_provider_->emit_attribute_changed(*p.voronoi_, p.voronoi_sphere_info_.get());
		
	}

	void init_coverage_axis_plus_plus(SURFACE& surface)
	{
		CoverageAxisParameter& p = coverage_axis_parameters_[&surface];
		p.surface_ = &surface;
		p.initialized_ = true;
		//Surface attributes
		p.medial_axis_position_ = get_or_add_attribute<Vec3, SurfaceVertex>(surface, "medial_axis_position");
		p.medial_axis_radius_ = get_or_add_attribute<Scalar, SurfaceVertex>(surface, "medial_axis_radius");
		p.medial_axis_closest_points_ = get_or_add_attribute<std::pair<SurfaceVertex, SurfaceVertex>, SurfaceVertex>(surface, "medial_axis_closest_points");
		p.medial_axis_angle_ = get_or_add_attribute<Scalar, SurfaceVertex>(surface, "medial_axis_angle");
		p.surface_vertex_normal_ = get_or_add_attribute<Vec3, SurfaceVertex>(surface, "normal");
		geometry::compute_normal<SurfaceVertex>(surface, p.surface_vertex_position_.get(),
												p.surface_vertex_normal_.get());


		//Point attributes
		p.candidate_points_ = point_provider_->add_mesh(point_provider_->mesh_name(surface) + "candidates");
		p.candidate_points_position_ = get_or_add_attribute<Vec3, PointVertex>(*p.candidate_points_, "position");
		p.candidate_points_radius_ = get_or_add_attribute<Scalar, PointVertex>(*p.candidate_points_, "radius");
		p.candidate_points_score_ = get_or_add_attribute<Scalar, PointVertex>(*p.candidate_points_, "score");
		p.candidate_points_coverage_score_ = get_or_add_attribute<Scalar, PointVertex>(*p.candidate_points_, "coverage_score");
		p.candidate_points_uniformity_score_ = get_or_add_attribute<Scalar, PointVertex>(*p.candidate_points_, "uniformity_score");
		p.candidate_points_associated_vertex_ = get_or_add_attribute<NonManifoldVertex, PointVertex>(*p.candidate_points_, "associated_vertex");


		p.selected_points_ = point_provider_->add_mesh(point_provider_->mesh_name(surface) + "selected");
		p.selected_points_position_ = get_or_add_attribute<Vec3, PointVertex>(*p.selected_points_, "position");
		p.selected_points_radius_ = get_or_add_attribute<Scalar, PointVertex>(*p.selected_points_, "radius");
		p.selected_points_associated_vertex_ = get_or_add_attribute<NonManifoldVertex, PointVertex>(*p.selected_points_, "associated_vertex");

		//Non-manifold attributes
		p.voronoi_ = nonmanifold_provider_->add_mesh(surface_provider_->mesh_name(surface) + "_coverage_axis");
		p.voronoi_position_ = get_or_add_attribute<Vec3, NonManifoldVertex>(*p.voronoi_, "position");
		p.voronoi_radius_ = get_or_add_attribute<Scalar, NonManifoldVertex>(*p.voronoi_, "radius");
		p.voronoi_stability_ratio_ = get_or_add_attribute<Scalar, NonManifoldEdge>(*p.voronoi_, "stability_ratio");
		p.voronoi_stability_color_ = get_or_add_attribute<Vec3, NonManifoldEdge>(*p.voronoi_, "stability_color");
		p.voronoi_sphere_info_ = get_or_add_attribute<Vec4, NonManifoldVertex>(*p.voronoi_, "sphere_info");
		p.voronoi_fixed_vertices_ = get_or_add_attribute<bool, NonManifoldVertex>(*p.voronoi_, "fixed_vertices");
		

		//Compute initalize non-manifold for Qmat

		load_model_in_cgal(surface, p.csm_);
		p.tree_ = std::make_shared<Tree>(faces(p.csm_).first, faces(p.csm_).second, p.csm_);
		p.tree_->accelerate_distance_queries();
		p.inside_tester_ = std::make_shared<Point_inside>(*(p.tree_.get()));

		//Construct non-manifold
		compute_initial_non_manifold(p);
		uint32 nb_vertices = nb_cells<NonManifoldVertex>(*p.voronoi_);
		//build Kd-tree for the voronoi diagram
		p.voronoi_kdt_vertices_.clear();
		p.voronoi_kdt_vertices_.reserve(nb_vertices);
		std::vector<Vec3> vertex_position_vector;
		vertex_position_vector.reserve(nb_vertices);
		uint32 idx = 0;
		foreach_cell(*p.voronoi_, [&](NonManifoldVertex v) -> bool {
			p.voronoi_kdt_vertices_.push_back(v);
			vertex_position_vector.push_back(value<Vec3>(*p.voronoi_, p.voronoi_position_, v));
			return true;
		});
		p.voronoi_kdt_ = std::make_shared<acc::KDTree<3, uint32>>(vertex_position_vector);
	
	}
	void random_sampling(CoverageAxisParameter& p)
	{

		//TODO
	}

	void delaunay_sampling(CoverageAxisParameter& p)
	{
		foreach_cell(*p.voronoi_, [&](NonManifoldVertex nv) {
			PointVertex new_candidate = add_vertex(*p.candidate_points_);
			p.candidate_inner_points.push_back(new_candidate);
			value<Vec3>(*p.candidate_points_, p.candidate_points_position_, new_candidate) =
				value<Vec3>(*p.voronoi_, p.voronoi_position_, nv);
			value<Scalar>(*p.candidate_points_, p.candidate_points_radius_, new_candidate) =
				value<Scalar>(*p.voronoi_, p.voronoi_radius_, nv);
			//bind non-manifold vertex to candidate point
			value<NonManifoldVertex>(*p.candidate_points_, p.candidate_points_associated_vertex_, new_candidate) = nv;

			return true;
		});
	}

	void shrinking_ball_sampling(CoverageAxisParameter& p)
	{
		geometry::shrinking_ball_centers<SURFACE>(*p.surface_, p.surface_vertex_position_.get(),
												  p.surface_vertex_normal_.get(), p.medial_axis_position_.get(),
												  p.medial_axis_radius_.get(), p.medial_axis_closest_points_.get());
		foreach_cell(*p.surface_, [&](SurfaceVertex v) -> bool {
			const Vec3& c = value<Vec3>(*p.surface_, p.medial_axis_position_, v);
			const auto& [v1, v2] =
				value<std::pair<SurfaceVertex, SurfaceVertex>>(*p.surface_, p.medial_axis_closest_points_, v);
			if (v2.is_valid())
			{
				const auto& c1 = value<Vec3>(*p.surface_, p.surface_vertex_position_, v1);
				const auto& c2 = value<Vec3>(*p.surface_, p.surface_vertex_position_, v2);
				const Scalar r = value<Scalar>(*p.surface_, p.medial_axis_radius_, v);
				p.max_radius_ = std::max(p.max_radius_, r);
				p.min_radius_ = std::min(p.min_radius_, r);
				const Scalar angle = geometry::angle(c1 - c, c2 - c);
				p.max_angle_ = std::max(p.max_angle_, angle);
				p.min_angle_ = std::min(p.min_angle_, angle);
				value<Scalar>(*p.surface_, p.medial_axis_angle_, v) = angle;

				PointVertex new_candidate = add_vertex(*p.candidate_points_);
				p.candidate_inner_points.push_back(new_candidate);
				value<Vec3>(*p.candidate_points_, p.candidate_points_position_, new_candidate) = c;
				value<Scalar>(*p.candidate_points_, p.candidate_points_radius_, new_candidate) = r;
				// bind non-manifold vertex to candidate point

				std::pair<uint32, Scalar> k_res;
				if (!p.voronoi_kdt_->find_nn(c, &k_res))
				{
					std::cout << "closest point not found !!!";
					return true;
				}
				value<NonManifoldVertex>(*p.candidate_points_, p.candidate_points_associated_vertex_, new_candidate) =
					p.voronoi_kdt_vertices_[k_res.first];

			}
			return true;
		});

		
	}
	void generate_candidates(CoverageAxisParameter& p)
	{
		switch (p.candidate_generation_method)
		{
		case SHRINKING_BALL: {
			shrinking_ball_sampling(p);
		}
		break;
		case RANDOM: {
		}
		break;
		case DELAUNAY: {
			delaunay_sampling(p);
		}
		break;
		default:
			break;
		}
		p.candidates_valid = true;
		point_provider_->emit_connectivity_changed(*p.candidate_points_);
		point_provider_->emit_attribute_changed(*p.candidate_points_, p.candidate_points_position_.get());
		point_provider_->emit_attribute_changed(*p.candidate_points_, p.candidate_points_radius_.get());
	}

	void compute_coverage_matrix(CoverageAxisParameter& p)
	{
		p.coverage_matrix.resize(nb_cells<SurfaceVertex>(*p.surface_), nb_cells<PointVertex>(*p.candidate_points_));
		p.coverage_matrix.setZero();
		foreach_cell(*p.surface_, [&](SurfaceVertex sample) {
			Vec3 pos = value<Vec3>(*p.surface_, p.surface_vertex_position_, sample);
			foreach_cell(*p.candidate_points_, [&](PointVertex candidate) {
				Vec3 center = value<Vec3>(*p.candidate_points_, p.candidate_points_position_, candidate);
				Scalar radius = value<Scalar>(*p.candidate_points_, p.candidate_points_radius_, candidate);
				if (inside_sphere(pos, center, radius + p.dilation_factor))
				{
					p.coverage_matrix(index_of(*p.surface_, sample), index_of(*p.candidate_points_, candidate)) = 1;
				}
				return true;
			});
			return true;
		});
	}

	Scalar compute_coverage_score(CoverageAxisParameter& p, PointVertex pv)
	{
		return p.coverage_matrix.col(index_of(*p.candidate_points_, pv)).sum();
	}

	Scalar compute_uniform_score(CoverageAxisParameter& p, PointVertex pv)
	{
		Scalar min_dist = std::numeric_limits<Scalar>::max();
		Vec3 pos = value<Vec3>(*p.candidate_points_, p.candidate_points_position_, pv);
		for (PointVertex& pi : p.selected_inner_points)
		{
			Scalar dis = (pos - value<Vec3>(*p.candidate_points_, p.candidate_points_position_, pi)).norm();
			min_dist = std::min(dis, min_dist);
		}
		return min_dist;
	}

	void compute_score(CoverageAxisParameter& p, std::vector<PointVertex>& candidate_points)
	{
		Scalar coverage_sum = 0;
		Scalar uniform_sum = 0;
		for (PointVertex& p_minus : candidate_points)
		{
			value<Scalar>(*p.candidate_points_, p.candidate_points_coverage_score_, p_minus) =
				compute_coverage_score(p, p_minus);
			coverage_sum += value<Scalar>(*p.candidate_points_, p.candidate_points_coverage_score_, p_minus);
			value<Scalar>(*p.candidate_points_, p.candidate_points_uniformity_score_, p_minus) =
				compute_uniform_score(p, p_minus);
			uniform_sum += value<Scalar>(*p.candidate_points_, p.candidate_points_uniformity_score_, p_minus);
		}
		Scalar mean_coverage = coverage_sum / candidate_points.size();
		Scalar mean_uniform = uniform_sum / candidate_points.size();
		Scalar deviation_coverage = 0;
		Scalar deviation_uniform = 0;
		
		for (PointVertex& p_minus : candidate_points)
		{
			Scalar& coverage_score = value<Scalar>(*p.candidate_points_, p.candidate_points_coverage_score_, p_minus);
			Scalar& uniform_score = value<Scalar>(*p.candidate_points_, p.candidate_points_uniformity_score_, p_minus);
			deviation_coverage += (coverage_score - mean_coverage) * (coverage_score - mean_coverage);
			deviation_uniform += (uniform_score - mean_uniform) * (uniform_score - mean_uniform);
		}
		deviation_coverage = sqrt(deviation_coverage);
		deviation_uniform = sqrt(deviation_uniform);

		for (PointVertex& p_minus : candidate_points)
		{
			Scalar coverage_normalised=  
				(value<Scalar>(*p.candidate_points_, p.candidate_points_coverage_score_, p_minus) - mean_coverage) /
				deviation_coverage;
			Scalar uniform_normalised =
				(value<Scalar>(*p.candidate_points_, p.candidate_points_uniformity_score_, p_minus) - mean_uniform) /
				deviation_uniform;
			value<Scalar>(*p.candidate_points_, p.candidate_points_score_, p_minus) =
				coverage_normalised + uniform_normalised;
		}
	}

	void coverage_axis_plus_plus(CoverageAxisParameter& p)
	{
		clear(*p.selected_points_);
		compute_coverage_matrix(p);
		p.selected_inner_points.clear();
		auto copie_candidates = std::vector<PointVertex>(p.candidate_inner_points);
		uint32 covered_number = 0;
		while (covered_number < nb_cells<SurfaceVertex>(*p.surface_) && 
			p.selected_inner_points.size()<p.max_selected_number)
		{
			compute_score(p, copie_candidates);
			auto& max_element_it =
				std::max_element(copie_candidates.begin(), copie_candidates.end(),
								 [&p](PointVertex& v1, PointVertex& v2) {
									 return value<Scalar>(*p.candidate_points_, p.candidate_points_score_, v1) <
											value<Scalar>(*p.candidate_points_, p.candidate_points_score_, v2);
								 });
			if (max_element_it != copie_candidates.end())
			{
				PointVertex max_score_vertex = *max_element_it;
				copie_candidates.erase(max_element_it);
				p.selected_inner_points.push_back(max_score_vertex);
				
				//Update coverage matrix;
				uint32 index = index_of(*p.candidate_points_, max_score_vertex);
				
				for (size_t idx = 0; idx < p.coverage_matrix.rows(); ++idx)
				{
					if (p.coverage_matrix(idx, index) != 0)
					{
						covered_number++;
						p.coverage_matrix.row(idx).setZero();
					}
				}
				PointVertex new_selected = add_vertex(*p.selected_points_);

				value<Vec3>(*p.selected_points_, p.selected_points_position_, new_selected) =
					value<Vec3>(*p.candidate_points_, p.candidate_points_position_, max_score_vertex);
				value<Scalar>(*p.selected_points_, p.selected_points_radius_, new_selected) =
					value<Scalar>(*p.candidate_points_, p.candidate_points_radius_, max_score_vertex);
				value<NonManifoldVertex>(*p.selected_points_, p.selected_points_associated_vertex_, new_selected) =
					value<NonManifoldVertex>(*p.candidate_points_, p.candidate_points_associated_vertex_, max_score_vertex);
			}
			
		}
		point_provider_->emit_connectivity_changed(*p.selected_points_);
		point_provider_->emit_attribute_changed(*p.selected_points_, p.selected_points_position_.get());
		point_provider_->emit_attribute_changed(*p.selected_points_, p.selected_points_radius_.get());
	}

	void coverage_axis_plus_plus_connectivity(CoverageAxisParameter& p)
	{
		using QMatHelper = modeling::DecimationSQEM_Helper<NONMANIFOLD>;
		using Slab_Quadric = geometry::Slab_Quadric;

		foreach_cell(*p.voronoi_, [&](NonManifoldEdge e) -> bool {
			auto iv = incident_vertices(*p.voronoi_, e);
			NonManifoldVertex v1 = iv[0];
			NonManifoldVertex v2 = iv[1];
			const Vec3& v1_p = value<Vec4>(*p.voronoi_, p.voronoi_sphere_info_, v1).head<3>();
			const Vec3& v2_p = value<Vec4>(*p.voronoi_, p.voronoi_sphere_info_, v2).head<3>();
			const Scalar & r1 = value<Vec4>(*p.voronoi_, p.voronoi_sphere_info_, v1).w();
			const Scalar& r2 = value<Vec4>(*p.voronoi_, p.voronoi_sphere_info_, v2).w();
			const Scalar center_dist = (v1_p - v2_p).norm();
			Scalar dis = std::max(0.0, (center_dist - std::abs(r1 - r2)));
			if (center_dist == 0.0)
			{
				(*p.voronoi_stability_ratio_)[e.index_] = 0.0;
				(*p.voronoi_stability_color_)[e.index_] = Vec3(0, 0, 0.5);
				return true;
			}
			Scalar stability = dis / center_dist;
			value<Scalar>(*p.voronoi_, p.voronoi_stability_ratio_, e) = stability;
			value<Vec3>(*p.voronoi_, p.voronoi_stability_color_, e) =
				(stability <= 0.5) ? Vec3(0, stability, (0.5 - stability)) : Vec3(stability - 0.5, (1 - stability), 0);
			return true;
		});
		parallel_foreach_cell(*p.voronoi_, [&](NonManifoldVertex nv) -> bool {
			value<bool>(*p.voronoi_, p.voronoi_fixed_vertices_, nv) = false;
			return true;
			});
		
		int number_vertex_remain = 0;
		foreach_cell(*p.selected_points_, [&](PointVertex pv) -> bool {
			NonManifoldVertex nv = value<NonManifoldVertex>(*p.selected_points_, p.selected_points_associated_vertex_, pv);
			value<bool>(*p.voronoi_, p.voronoi_fixed_vertices_, nv) = true;
			number_vertex_remain++;
			return true;
		});
		QMatHelper helper(0.00001, *p.voronoi_, p.voronoi_position_, p.voronoi_sphere_info_, p.voronoi_stability_color_,
						  p.voronoi_stability_ratio_, p.voronoi_radius_, p.voronoi_fixed_vertices_);

		helper.initial_slab_mesh();
		helper.initial_boundary_mesh();
		helper.initial_collapse_queue();
		helper.simplify(number_vertex_remain, true);

		nonmanifold_provider_->emit_connectivity_changed(*p.voronoi_);
		nonmanifold_provider_->emit_attribute_changed(*p.voronoi_, p.voronoi_position_.get());
		nonmanifold_provider_->emit_attribute_changed(*p.voronoi_, p.voronoi_stability_color_.get());
		nonmanifold_provider_->emit_attribute_changed(*p.voronoi_, p.voronoi_stability_ratio_.get());
		nonmanifold_provider_->emit_attribute_changed(*p.voronoi_, p.voronoi_radius_.get());
		
	}

	using SphereQueue = std::multimap<Scalar, SurfaceVertex,std::greater<Scalar>>;
	using SphereQueueIt = typename SphereQueue::const_iterator;
	using SphereInfo = std::pair<bool, SphereQueueIt>;

	
	struct ClusterAxisParameter
	{
		bool initialized_ = false;

		SURFACE* surface_;
		std::shared_ptr<SurfaceAttribute<Vec3>> surface_vertex_position_ = nullptr;
		std::shared_ptr<SurfaceAttribute<Vec3>> surface_vertex_normal_ = nullptr;
		std::shared_ptr<SurfaceAttribute<Vec3>> surface_face_normal_ = nullptr;
		std::shared_ptr<SurfaceAttribute<Vec3>> surface_vertex_color_ = nullptr;
		std::shared_ptr<SurfaceAttribute<Scalar>> surface_vertex_area_ = nullptr; 
		std::shared_ptr<SurfaceAttribute<Vec3>> medial_axis_position_ = nullptr;
		std::shared_ptr<SurfaceAttribute<Scalar>> medial_axis_radius_ = nullptr;
		std::shared_ptr<SurfaceAttribute<Scalar>> surface_distance_to_cluster_ = nullptr;
		std::shared_ptr<SurfaceAttribute<Scalar>> surface_distance_to_enveloppe= nullptr;
		std::shared_ptr<SurfaceAttribute<PointVertex>> surface_cluster_info_ = nullptr;
		std::shared_ptr<SurfaceAttribute<uint32>> surface_clusters_cc_idx_ = nullptr;
		std::shared_ptr<SurfaceAttribute<std::vector<std::pair<Scalar, Vec3>>>> clusters_fuzzy_cluster_color_ = nullptr;
		std::shared_ptr<SurfaceAttribute<std::pair<SurfaceVertex,SurfaceVertex>>> medial_axis_closest_points_ = nullptr;
		std::shared_ptr<SurfaceAttribute<Vec3>> medial_axis_cloest_point_color_ = nullptr;
		

		std::vector<SurfaceFace> surface_bvh_faces_;
		std::vector<SurfaceVertex> surface_kdt_vertices_;
		std::shared_ptr<acc::BVHTree<uint32, Vec3>> surface_bvh_ = nullptr;
		std::shared_ptr < acc::KDTree<3, uint32>> surface_kdt_ = nullptr;
		Cgal_Surface_mesh csm_;
		std::shared_ptr<Tree> tree_ = nullptr;
		std::shared_ptr<Point_inside> inside_tester_ = nullptr;
		std::shared_ptr<modeling::ClusteringSQEM_Helper<SURFACE>> sqem_helper_ = nullptr;
		
		POINT* clusters_ = nullptr;
		
		std::shared_ptr<PointAttribute<Vec3>> clusters_position_ = nullptr;
		std::shared_ptr<PointAttribute<Scalar>> clusters_radius_ = nullptr;
		std::shared_ptr<PointAttribute<Vec4>> clusters_color_ = nullptr;
		std::shared_ptr<PointAttribute<Vec4>> clusters_surface_color_ = nullptr;
		std::shared_ptr<PointAttribute<std::vector<SurfaceVertex>>> clusters_surface_vertices_ = nullptr;
		std::shared_ptr<PointAttribute<std::vector<SurfaceVertex>>> clusters_fuzzy_surface_vertices_ = nullptr;
		std::shared_ptr<PointAttribute<Scalar>> clusters_error_ = nullptr;
		std::shared_ptr<PointAttribute<Scalar>> clusters_deviation_ = nullptr;
		std::shared_ptr<PointAttribute<std::set<PointVertex>>> clusters_neighbours_ = nullptr;
		std::shared_ptr<PointAttribute<std::unordered_map<PointVertex*, std::set<uint32>>>> clusters_adjacency_info_ =
			nullptr;
		std::shared_ptr<PointAttribute<std::pair<SurfaceVertex, SurfaceVertex>>> clusters_cloest_sample_ = nullptr;
		std::shared_ptr<PointAttribute<Scalar>> clusters_max_distance_ = nullptr;
		std::shared_ptr<PointAttribute<SurfaceVertex>> clusters_max_vertex_ = nullptr;
		std::shared_ptr<PointAttribute<Scalar>> correction_error_ = nullptr;
		std::shared_ptr<PointAttribute<Vec3>> clusters_without_correction_position_ = nullptr;
		std::shared_ptr<PointAttribute<Scalar>> clusters_without_correction_radius_ = nullptr;
		std::shared_ptr<PointAttribute<uint32>> clusters_cc_number_ = nullptr;

		NONMANIFOLD* non_manifold_ = nullptr;
		std::shared_ptr<NonManifoldAttribute<Vec3>> non_manifold_vertex_position_ = nullptr;
		std::shared_ptr<NonManifoldAttribute<Scalar>> non_manifold_vertex_radius_ = nullptr;
		std::shared_ptr<NonManifoldAttribute<Vec4>> non_manifold_sphere_info_ = nullptr;
		std::shared_ptr<NonManifoldAttribute<std::vector<SurfaceVertex>>> non_manifold_cluster_vertices_ = nullptr;
																						 
		PointAttribute<Scalar>* selected_clusters_error_ = nullptr;
		std::shared_ptr<CellMarkerStore<SURFACE, SurfaceVertex>> no_ball;
		float init_min_radius_ = 0.01;
		float init_cover_dist_ = 0.06;

		UpdateMethod update_method_ = SPHERE_FITTING;
		SplitMethod split_method_  = DEVIATION;
		InitMethod init_method_ = CONSTANT;
		float mean_update_curvature_weight_ = 0.2;
		bool auto_split_outside_spheres_ = false;
		float split_variance_threshold_ = 0.0001f;
		float split_distance_threshold_ = 0.01f;
		float split_distance_to_medial_axis_threshold_ = 0.0001f;
		float split_hausdorff_distance_threshold_ = 0.0001f;
		float fuzzy_distance_ = 0.001f;
		Scalar total_error_ = 0.0;
		Scalar min_error_ = 0.0;
		Scalar max_error_ = 0.0;
		float init_distance_ = 0.1f;
		float init_factor_ = 0.8f;
		Scalar hausdorff_distance_enveloppe_ = 0.0;
		Scalar huasdorff_distance_cluster_ = 0.0;

		std::mutex mutex_;
		bool running_ = false;
		bool stopping_ = false;
		bool slow_down_ = true;
		bool fuzzy_clustering_ = false;
		bool connectivity_surgery = false;
		bool compute_enveloppe_distance_ = false;
		Scalar dilation_factor = 0.02;

		uint32 update_rate_ = 20;
		uint32 max_split_number = 1;

		cgogn::rendering::SkelShapeDrawer skel_drawer_;
		bool draw_enveloppe = false;
	};

	void cluster_axis_init(SURFACE& surface)
	{
		ClusterAxisParameter& p = cluster_axis_parameters_[&surface];
		p.surface_ = &surface;
		if (p.initialized_)
		{
			std::cout << "Surface data already initialized" << std::endl;
			return;
		}

		// create BVH and KDTree for the surface
		p.no_ball = std::make_shared<CellMarkerStore<SURFACE, SurfaceVertex>>(surface);
		MeshData<SURFACE>& md = surface_provider_->mesh_data(surface);
		uint32 nb_vertices = md.template nb_cells<SurfaceVertex>();
		uint32 nb_faces = md.template nb_cells<SurfaceFace>();

		auto bvh_vertex_index = add_attribute<uint32, SurfaceVertex>(surface, "__bvh_vertex_index");

		p.surface_kdt_vertices_.clear();
		p.surface_kdt_vertices_.reserve(nb_vertices);
		std::vector<Vec3> vertex_position_vector;
		vertex_position_vector.reserve(nb_vertices);
		uint32 idx = 0;
		foreach_cell(surface, [&](SurfaceVertex v) -> bool {
			p.surface_kdt_vertices_.push_back(v);
			value<uint32>(surface, bvh_vertex_index, v) = idx++;
			vertex_position_vector.push_back(value<Vec3>(surface, p.surface_vertex_position_, v));
			return true;
		});

		p.surface_bvh_faces_.clear();
		p.surface_bvh_faces_.reserve(nb_faces);
		std::vector<uint32> face_vertex_indices;
		face_vertex_indices.reserve(nb_faces * 3);
		foreach_cell(surface, [&](SurfaceFace f) -> bool {
			p.surface_bvh_faces_.push_back(f);
			foreach_incident_vertex(surface, f, [&](SurfaceVertex v) -> bool {
				face_vertex_indices.push_back(value<uint32>(surface, bvh_vertex_index, v));
				return true;
			});
			return true;
		});

		p.surface_bvh_ = std::make_shared<acc::BVHTree<uint32, Vec3>>(face_vertex_indices, vertex_position_vector);
		p.surface_kdt_ = std::make_shared < acc::KDTree<3, uint32>>(vertex_position_vector);
		
		//Create inside tester
		
		load_model_in_cgal(surface, p.csm_);
		p.tree_ = std::make_shared<Tree>(faces(p.csm_).first, faces(p.csm_).second, p.csm_);
		p.tree_->accelerate_distance_queries();
		p.inside_tester_ = std::make_shared<Point_inside>(*(p.tree_.get()));

		remove_attribute<SurfaceVertex>(surface, bvh_vertex_index);

		// compute normals

		p.surface_face_normal_ = get_or_add_attribute<Vec3, SurfaceFace>(surface, "normal");
		geometry::compute_normal<SurfaceFace>(surface, p.surface_vertex_position_.get(), p.surface_face_normal_.get());

		p.surface_vertex_normal_ = get_or_add_attribute<Vec3, SurfaceVertex>(surface, "normal");
		geometry::compute_normal<SurfaceVertex>(surface, p.surface_vertex_position_.get(),
												p.surface_vertex_normal_.get());
		p.surface_vertex_area_ = get_or_add_attribute<Scalar, SurfaceVertex>(surface, "area");
		geometry::compute_area<SurfaceVertex>(surface, p.surface_vertex_position_.get(), p.surface_vertex_area_.get());

		p.sqem_helper_ = std::make_shared<modeling::ClusteringSQEM_Helper<SURFACE>>(
			surface, p.surface_vertex_position_.get(), p.surface_vertex_normal_.get());
		// compute shrinking balls for the surface vertices

		p.medial_axis_position_ = get_or_add_attribute<Vec3, SurfaceVertex>(surface, "medial_axis_position");
		p.medial_axis_radius_ = get_or_add_attribute<Scalar, SurfaceVertex>(surface, "medial_axis_radius");
		p.surface_vertex_color_ = get_or_add_attribute<Vec3, SurfaceVertex>(surface, "surface_vertex_color");
		p.medial_axis_cloest_point_color_ = get_or_add_attribute<Vec3, SurfaceVertex>(surface, "medial_axis_cloest_point_color_");
		p.surface_distance_to_cluster_ = get_or_add_attribute<Scalar, SurfaceVertex>(surface, "surface_distance_to_cluster_");
		p.surface_distance_to_enveloppe = get_or_add_attribute<Scalar, SurfaceVertex>(surface, "surface_distance_to_enveloppe");
		p.surface_cluster_info_ = get_or_add_attribute<PointVertex, SurfaceVertex>(surface, "surface_cluster_info_");
		p.surface_clusters_cc_idx_ = get_or_add_attribute<uint32, SurfaceVertex>(surface, "surface_clusters_cc_idx_");
		p.medial_axis_closest_points_ = get_or_add_attribute<std::pair<SurfaceVertex, SurfaceVertex>, SurfaceVertex>(
			surface, "medial_axis_closest_points_");
		p.clusters_fuzzy_cluster_color_ = get_or_add_attribute<std::vector<std::pair<Scalar, Vec3>>, SurfaceVertex>(surface, "clusters_fuzzy_cluster_color_");
		
		parallel_foreach_cell(surface, [&](SurfaceVertex v) -> bool {
			uint32 index = index_of(surface, v);
			auto [c, r, q] = geometry::shrinking_ball_center(
				surface, value<Vec3>(surface, p.surface_vertex_position_, v),
				value<Vec3>(surface, p.surface_vertex_normal_, v), p.surface_vertex_position_.get(),
				p.surface_bvh_.get(), p.surface_bvh_faces_, p.surface_kdt_.get(),
				p.surface_kdt_vertices_);
			if (!q.is_valid())
			{
				std::cout << "surface vertex : " << index_of(surface, v)
						  << " has no intersection in the mesh"<<std::endl;
				p.no_ball->mark(v);
				return true;
			}value<Vec3>(surface, p.medial_axis_position_, v) = c;
			value<Scalar>(surface, p.medial_axis_radius_, v) = r;
		
			value<std::pair<SurfaceVertex, SurfaceVertex>>(surface, p.medial_axis_closest_points_, v) = {v, q};
			return true;
		});

		p.clusters_ = point_provider_->add_mesh(point_provider_->mesh_name(*p.surface_) + "clusters");
		
		p.clusters_position_ = get_or_add_attribute<Vec3, PointVertex>(*p.clusters_, "clusters_position");
		p.clusters_radius_ = get_or_add_attribute<Scalar, PointVertex>(*p.clusters_, "clusters_radius_");
		p.clusters_color_ = get_or_add_attribute<Vec4, PointVertex>(*p.clusters_, "clusters_color");
		p.clusters_surface_color_ = get_or_add_attribute<Vec4, PointVertex>(*p.clusters_, "cluster_surface_color");
		p.clusters_surface_vertices_ = get_or_add_attribute<std::vector<SurfaceVertex>, PointVertex>(
			*p.clusters_, "clusters_surface_vertices");
		p.clusters_fuzzy_surface_vertices_ = get_or_add_attribute<std::vector<SurfaceVertex>, PointVertex>(
			*p.clusters_, "clusters_fuzzy_surface_vertices_");

		p.clusters_error_ = get_or_add_attribute<Scalar, PointVertex>(*p.clusters_, "clusters_error");
		p.clusters_deviation_ = get_or_add_attribute<Scalar, PointVertex>(*p.clusters_, "clusters_deviation");
		p.clusters_neighbours_ = get_or_add_attribute<std::set<PointVertex>, PointVertex>(
			*p.clusters_, "clusters_neighbours");
		p.clusters_cloest_sample_ = get_or_add_attribute<std::pair<SurfaceVertex, SurfaceVertex>, PointVertex>(
			*p.clusters_, "clusters_cloest_sample_");
		p.clusters_max_distance_ = get_or_add_attribute<Scalar, PointVertex>(
			*p.clusters_, "clusters_max_distance");
		p.clusters_max_vertex_ = get_or_add_attribute<SurfaceVertex, PointVertex>(*p.clusters_, "cluster_max_vertex");
		p.correction_error_ = get_or_add_attribute<Scalar, PointVertex>(*p.clusters_, "correction_error");
		p.clusters_cc_number_ = get_or_add_attribute<uint32, PointVertex>(*p.clusters_, "clusters_cc_number");
		p.clusters_adjacency_info_ = get_or_add_attribute<std::unordered_map<PointVertex*, std::set<uint32>>, PointVertex>(
			*p.clusters_, "clusters_adjacency_info");
		p.selected_clusters_error_ = p.clusters_error_.get();

		
		p.clusters_without_correction_position_ =
			get_or_add_attribute<Vec3, PointVertex>(*p.clusters_, "clusters_without_correction_position");
		p.clusters_without_correction_radius_ =
			get_or_add_attribute<Scalar, PointVertex>(*p.clusters_, "clusters_without_correction_radius");


		p.non_manifold_ =
			nonmanifold_provider_->add_mesh("skeleton" + std::to_string(nonmanifold_provider_->number_of_meshes()));
		p.non_manifold_vertex_position_ = get_or_add_attribute<Vec3, NonManifoldVertex>(*p.non_manifold_, "position");
		p.non_manifold_sphere_info_ = get_or_add_attribute<Vec4, NonManifoldVertex>(*p.non_manifold_, "sphere_info");
		p.non_manifold_vertex_radius_ = get_or_add_attribute<Scalar, NonManifoldVertex>(*p.non_manifold_, "radius");
		p.non_manifold_cluster_vertices_ =
			get_or_add_attribute<std::vector<SurfaceVertex>, NonManifoldVertex>(*p.non_manifold_, "cluster_vertices");

		p.skel_drawer_.set_color({1, 0, 0, 0.5});
		p.skel_drawer_.set_subdiv(40);

		p.initialized_ = true;

		
	}

	void set_selected_clusters_error(ClusterAxisParameter& p, PointAttribute<Scalar>* attribute)
	{
		p.selected_clusters_error_ = attribute;
		update_clusters_color(p);
		if (!p.running_)
			update_render_data(p);
	}
	void initialise_cluster(ClusterAxisParameter& p)
	{
		clear(*p.clusters_);
		
		foreach_cell(*p.surface_, [&](SurfaceVertex sv) {
			value<Vec3>(*p.surface_, p.surface_vertex_color_, sv) = Vec3(0, 0, 0);
			return true;
		});
		
		SphereQueue sphere_queue;
		auto sphere_info = get_or_add_attribute<SphereInfo, SurfaceVertex>(*p.surface_, "sphere_info");
		//Push all the medial axis sphere into the maximum queue 
		foreach_cell(*p.surface_, [&](SurfaceVertex v) {
			if (p.no_ball->is_marked(v))
				return true;
			value<SphereInfo>(*p.surface_, sphere_info, v) = {
				true, sphere_queue.emplace(value<Scalar>(*p.surface_, p.medial_axis_radius_, v), v)};
			return true;
		});
		std::cout << "finish pushing sph��res" << std::endl;

		CellMarker<SURFACE, SurfaceVertex> bfs_marker(*p.surface_);
		bfs_marker.unmark_all();
		std::queue<SurfaceVertex> vertex_queue;
		//Start from the biggest sphere, affect the surface points to the sphere if the distance 
		// between surface point and sphere is less than threshold
		std::cout << "start BFS " << std::endl;
		while (sphere_queue.size() > 0)
		{
			uint32 size = sphere_queue.size();
			SphereQueueIt it = sphere_queue.begin();
			auto [radius, v] = *it;
			if (!value<SphereInfo>(*p.surface_, sphere_info, v).first)
			{
				sphere_queue.erase(it);
				continue;
			}
			auto [v1, v2] = value<std::pair<SurfaceVertex, SurfaceVertex>>(*p.surface_, p.medial_axis_closest_points_, v);
			uint32 index = index_of(*p.surface_, v1);
			Vec3& sphere_center = value<Vec3>(*p.surface_, p.medial_axis_position_, v1);
			Scalar& sphere_radius = value<Scalar>(*p.surface_, p.medial_axis_radius_, v1);
			
			value<SphereInfo>(*p.surface_, sphere_info, v1).first = false;
	
			if (!value<SphereInfo>(*p.surface_, sphere_info, v2).first)
			{
				continue;
			}
			
			//Set new clusters info
			PointVertex new_cluster = add_vertex(*p.clusters_);
			value<Vec3>(*p.clusters_, p.clusters_position_, new_cluster) = sphere_center;
			value<std::pair<SurfaceVertex, SurfaceVertex>>(*p.clusters_, p.clusters_cloest_sample_,
														   new_cluster) = {v1, v2};
			
			value<Vec4>(*p.clusters_, p.clusters_surface_color_, new_cluster) =
				Vec4(rand() / (double)RAND_MAX, rand() / (double)RAND_MAX, rand() / (double)RAND_MAX,1.0);
			value<std::vector<SurfaceVertex>>(*p.clusters_, p.clusters_surface_vertices_, new_cluster).clear();
			value<Scalar>(*p.clusters_, p.clusters_radius_, new_cluster) = radius;
			value<Vec3>(*p.surface_, p.surface_vertex_color_, v1) = value<Vec4>(*p.clusters_, p.clusters_surface_color_, new_cluster).head<3>();
			value<Vec3>(*p.surface_, p.surface_vertex_color_, v2) =
				value<Vec4>(*p.clusters_, p.clusters_surface_color_, new_cluster).head<3>();

			// Add v1 into cluster and remove it from the queue
			vertex_queue.push(v1);
			Vec3& v1_pos = value<Vec3>(*p.surface_, p.surface_vertex_position_, v1);
			Vec3& v1_center = value<Vec3>(*p.surface_, p.medial_axis_position_, v1);
			Scalar v1_radius = value<Scalar>(*p.surface_, p.medial_axis_radius_, v1);
			value<PointVertex>(*p.surface_, p.surface_cluster_info_, v1) = new_cluster;
			value<Vec3>(*p.surface_, p.surface_vertex_color_, v1) =
				value<Vec4>(*p.clusters_, p.clusters_surface_color_, new_cluster).head<3>();
			value<Scalar>(*p.surface_, p.surface_distance_to_cluster_, v1) = (v1_pos - v1_center).norm() - v1_radius;
			/*value<SphereInfo>(*p.surface_, sphere_info, v1).first = false;*/
			
			// Add v2 into cluster and remove it from the queue
			vertex_queue.push(v2);
			Vec3& v2_pos = value<Vec3>(*p.surface_, p.surface_vertex_position_, v2);
			Vec3& v2_center = value<Vec3>(*p.surface_, p.medial_axis_position_, v2);
			Scalar v2_radius = value<Scalar>(*p.surface_, p.medial_axis_radius_, v2);
			value<PointVertex>(*p.surface_, p.surface_cluster_info_, v2) = new_cluster;
			value<Vec3>(*p.surface_, p.surface_vertex_color_, v2) =
				value<Vec4>(*p.clusters_, p.clusters_surface_color_, new_cluster).head<3>();
			value<Scalar>(*p.surface_, p.surface_distance_to_cluster_, v2) = (v2_pos - v2_center).norm() - v2_radius;
			value<SphereInfo>(*p.surface_, sphere_info, v2).first = false;
			
			//Start BFS
			bfs_marker.unmark_all();
			while (vertex_queue.size() > 0)
			{
				SurfaceVertex vf = vertex_queue.front();
				vertex_queue.pop();
				// Affect the surface points to the cluster
				foreach_adjacent_vertex_through_edge(*p.surface_, vf, [&](SurfaceVertex sv) {
					if (bfs_marker.is_marked(sv))
						return true;
					Vec3& pos = value<Vec3>(*p.surface_, p.surface_vertex_position_, sv);
					Vec4 sphere = Vec4(sphere_center.x(), sphere_center.y(), sphere_center.z(), radius);
					// if the distance is less than threshold
					Scalar error = (pos - sphere_center).norm() - radius;
					switch (p.init_method_)
					{
					case FACTOR: {
						error = error / radius;
						if (error < p.init_factor_ )
						{
							vertex_queue.push(sv);
							if (!value<SphereInfo>(*p.surface_, sphere_info, sv).first)
							{
								if (error < value<Scalar>(*p.surface_, p.surface_distance_to_cluster_, sv))
								{
									value<PointVertex>(*p.surface_, p.surface_cluster_info_, sv) = new_cluster;
									value<Vec3>(*p.surface_, p.surface_vertex_color_, sv) =
										value<Vec4>(*p.clusters_, p.clusters_surface_color_, new_cluster).head<3>();
									value<Scalar>(*p.surface_, p.surface_distance_to_cluster_, sv) = error;
								}
							}
							else
							{
								value<PointVertex>(*p.surface_, p.surface_cluster_info_, sv) = new_cluster;
								value<Vec3>(*p.surface_, p.surface_vertex_color_, sv) =
									value<Vec4>(*p.clusters_, p.clusters_surface_color_, new_cluster).head<3>();
								value<Scalar>(*p.surface_, p.surface_distance_to_cluster_, sv) = error;
								value<SphereInfo>(*p.surface_, sphere_info, sv).first = false;
							}
						}
						bfs_marker.mark(sv);
					}
					break;
					case CONSTANT: {
						if (error < p.init_distance_)
						{
							vertex_queue.push(sv);
							if (!value<SphereInfo>(*p.surface_, sphere_info, sv).first)
							{
								if (error < value<Scalar>(*p.surface_, p.surface_distance_to_cluster_, sv))
								{
									value<PointVertex>(*p.surface_, p.surface_cluster_info_, sv) = new_cluster;
									value<Vec3>(*p.surface_, p.surface_vertex_color_, sv) =
										value<Vec4>(*p.clusters_, p.clusters_surface_color_, new_cluster).head<3>();
									value<Scalar>(*p.surface_, p.surface_distance_to_cluster_, sv) = error;
								}
							}
							else
							{
								value<PointVertex>(*p.surface_, p.surface_cluster_info_, sv) = new_cluster;
								value<Vec3>(*p.surface_, p.surface_vertex_color_, sv) =
									value<Vec4>(*p.clusters_, p.clusters_surface_color_, new_cluster).head<3>();
								value<Scalar>(*p.surface_, p.surface_distance_to_cluster_, sv) = error;
								value<SphereInfo>(*p.surface_, sphere_info, sv).first = false;
							}
						}
						bfs_marker.mark(sv);
					}
					break;
					}
					return true;
					
				});
			}
		}
		remove_attribute<SurfaceVertex>(*p.surface_, sphere_info.get());
		                          
		assign_cluster(p);
		if (p.connectivity_surgery)
		{
			compute_clusters_cc_numbers(p);
			compute_correct_clusters_neighbor(p);
		}
		else
		{

			update_clusters_neighbor(p);
		}
		if (p.fuzzy_clustering_)
		{
			assign_fuzzy_clusters(p);
		}
		update_clusters_data(p);
		update_clusters_color(p);
		if (!p.running_)
			update_render_data(p);
		
	}
	
	void assign_cluster(ClusterAxisParameter& p)
	{
		CellMarker<SURFACE, SurfaceVertex> cluster_marker(*p.surface_);
		cluster_marker.unmark_all();
		std::queue<std::pair<SurfaceVertex, PointVertex>> vertex_queue;
		parallel_foreach_cell(*p.clusters_, [&](PointVertex pv) {
			value<std::vector<SurfaceVertex>>(*p.clusters_, p.clusters_surface_vertices_, pv).clear();
			value<std::set<PointVertex>>(*p.clusters_, p.clusters_neighbours_, pv).clear();
			value<uint32>(*p.clusters_, p.clusters_cc_number_, pv) = 0;
			value<std::unordered_map<PointVertex*, std::set<uint32>>>(*p.clusters_, p.clusters_adjacency_info_,
																	  pv).clear();
			return true;
		});
		MeshData<POINT>& md = point_provider_->mesh_data(*p.clusters_);
		uint32 nb_spheres = md.template nb_cells<PointVertex>();
		std::vector<Vec3> sphere_centers;
		sphere_centers.reserve(nb_spheres);
		std::vector<Scalar> sphere_radii;
		sphere_radii.reserve(nb_spheres);
		std::vector<PointVertex> spheres_bvh_vertices;
		spheres_bvh_vertices.reserve(nb_spheres);
		foreach_cell(*p.clusters_, [&](PointVertex v) -> bool {
			spheres_bvh_vertices.push_back(v);
			sphere_centers.push_back(value<Vec3>(*p.clusters_, p.clusters_position_, v));
			sphere_radii.push_back(value<Scalar>(*p.clusters_, p.clusters_radius_, v));
			return true;
		});
		acc::BVHTreeSpheres<uint32, Vec3> spheres_bvh(sphere_centers, sphere_radii);
		parallel_foreach_cell(*p.surface_, [&](SurfaceVertex v) -> bool {
	
			const Vec3& vp = value<Vec3>(*p.surface_, p.surface_vertex_position_, v);
			PointVertex closest_sphere;

			std::pair<uint32, Vec3> bvh_res;
			spheres_bvh.closest_point(vp, &bvh_res);
			closest_sphere = spheres_bvh_vertices[bvh_res.first];
			auto& s_center = sphere_centers[bvh_res.first];
			auto& s_radius = sphere_radii[bvh_res.first];
			value<PointVertex>(*p.surface_, p.surface_cluster_info_, v) = closest_sphere;
			value<Vec3>(*p.surface_, p.surface_vertex_color_, v) =
				value<Vec4>(*p.clusters_, p.clusters_surface_color_, closest_sphere).head<3>();
			value<Scalar>(*p.surface_, p.surface_distance_to_cluster_, v) = (vp-s_center).norm() - s_radius;
			std::lock_guard<std::mutex> lock(spheres_mutex_[bvh_res.first % spheres_mutex_.size()]);
			value<std::vector<SurfaceVertex>>(*p.clusters_, p.clusters_surface_vertices_, closest_sphere).push_back(v);

			return true;
		});
		foreach_cell(*p.clusters_, [&](PointVertex pv) {
			auto& clusters_surface_vertices =
				value<std::vector<SurfaceVertex>>(*p.clusters_, p.clusters_surface_vertices_, pv);
			if (clusters_surface_vertices.size() < 4)
			{
				const std::vector<SurfaceVertex>& cluster = value<std::vector<SurfaceVertex>>(*p.clusters_, p.clusters_surface_vertices_, pv);
				for (SurfaceVertex sv : cluster)
					value<PointVertex>(*p.surface_, p.surface_cluster_info_, sv) = PointVertex();
				std::cout << "delete cluster:" << index_of(*p.clusters_, pv) << std::endl;
				remove_vertex(*p.clusters_, pv);
			}
			return true;
		});
		
	}
	
	void compute_clusters_cc_numbers(ClusterAxisParameter& p)
	{
		
		CellMarker<SURFACE, SurfaceVertex> visited(*p.surface_);
		foreach_cell(*p.surface_, [&](SurfaceVertex sv) {
			uint32 nb_vertex = 0;
			if (visited.is_marked(sv))
				return true;
			std::queue<SurfaceVertex> queue;
			PointVertex cluster = value<PointVertex>(*p.surface_, p.surface_cluster_info_, sv);
			if (!cluster.is_valid())
				return true;
			uint32 idx = value<uint32>(*p.clusters_, p.clusters_cc_number_, cluster);
			queue.push(sv);
			visited.mark(sv);
			while (!queue.empty())
			{
				SurfaceVertex v = queue.front();
				value<uint32>(*p.surface_, p.surface_clusters_cc_idx_, v) = idx;
				queue.pop();
				nb_vertex++;
				foreach_adjacent_vertex_through_edge(*p.surface_, v, [&](SurfaceVertex n) {
					PointVertex adjacent_cluster = value<PointVertex>(*p.surface_, p.surface_cluster_info_, n);
					if (!visited.is_marked(n) && index_of(*p.clusters_, cluster) == index_of(*p.clusters_, adjacent_cluster))
					{
						queue.push(n);
						visited.mark(n);
					}
					return true;
				});
			}
			if (nb_vertex > 2)
				value<uint32>(*p.clusters_, p.clusters_cc_number_, cluster) +=1;
			return true;
		});
	}

	void compute_correct_clusters_neighbor(ClusterAxisParameter& p)
	{
		foreach_cell(*p.surface_, [&](SurfaceEdge e) -> bool {
			std::vector<SurfaceVertex> vertices = incident_vertices(*p.surface_, e);

			PointVertex& v1_cluster = value<PointVertex>(*p.surface_, p.surface_cluster_info_, vertices[0]);
			PointVertex& v2_cluster = value<PointVertex>(*p.surface_, p.surface_cluster_info_, vertices[1]);
			if (v1_cluster.is_valid() && v2_cluster.is_valid() && v1_cluster != v2_cluster)
			{
				uint32 v1_cc_id = value<uint32>(*p.surface_, p.surface_clusters_cc_idx_, vertices[0]);
				uint32 v2_cc_id = value<uint32>(*p.surface_, p.surface_clusters_cc_idx_, vertices[1]);
				value<std::unordered_map<PointVertex*, std::set<uint32>>>(*p.clusters_, p.clusters_adjacency_info_,
																		  v1_cluster)[&v2_cluster]
					.insert(v2_cc_id);
				value<std::unordered_map<PointVertex*, std::set<uint32>>>(*p.clusters_, p.clusters_adjacency_info_,
																		  v2_cluster)[&v1_cluster]
					.insert(v1_cc_id);
			}
			return true;
		});
		foreach_cell(*p.clusters_, [&](PointVertex pv) {
			auto& adjacency_info = value<std::unordered_map<PointVertex*, std::set<uint32>>>(*p.clusters_, p.clusters_adjacency_info_, pv);
			auto& cc_numbers = value<uint32>(*p.clusters_, p.clusters_cc_number_, pv);
			
			for (auto& it : adjacency_info)
			{
				PointVertex* cluster = it.first;
				auto& cc_number_other = value<uint32>(*p.clusters_, p.clusters_cc_number_, *cluster);
				if (cc_numbers<=cc_number_other)
				{
					if (it.second.size() == cc_number_other)
					{
						value<std::set<PointVertex>>(*p.clusters_, p.clusters_neighbours_, pv).insert(*cluster);
					}
				}
				
			}
			return true;
		});


	}

	// Update the neighborhood of the clusters
	void update_clusters_neighbor(ClusterAxisParameter& p)
	{
		foreach_cell(*p.surface_, [&](SurfaceEdge e) -> bool {
			std::vector<SurfaceVertex> vertices = incident_vertices(*p.surface_, e);

			PointVertex v1_cluster = value<PointVertex>(*p.surface_, p.surface_cluster_info_, vertices[0]);
			PointVertex v2_cluster = value<PointVertex>(*p.surface_, p.surface_cluster_info_, vertices[1]);
			
			if (v1_cluster.is_valid() && v2_cluster.is_valid() && v1_cluster != v2_cluster)
			{
				value<std::set<PointVertex>>(*p.clusters_, p.clusters_neighbours_, v1_cluster).insert(v2_cluster);
				value<std::set<PointVertex>>(*p.clusters_, p.clusters_neighbours_, v2_cluster).insert(v1_cluster);
			}
			return true;
		});
		
		
	}

	void assign_fuzzy_clusters(ClusterAxisParameter& p)
	{
		parallel_foreach_cell(*p.clusters_, [&](PointVertex pv)
			{ 
			// copy crisp vertices
			value<std::vector<SurfaceVertex>>(*p.clusters_, p.clusters_fuzzy_surface_vertices_, pv) =
				value<std::vector<SurfaceVertex>>(*p.clusters_, p.clusters_surface_vertices_, pv);
				return true;
			}
		);
		foreach_cell(*p.surface_, [&](SurfaceVertex sv) {
			if (value<Scalar>(*p.surface_, p.surface_distance_to_cluster_, sv) == 0)
				return true;
			auto& clusters_color =
				value<std::vector<std::pair<Scalar, Vec3>>>(*p.surface_, p.clusters_fuzzy_cluster_color_, sv);
			clusters_color.clear();
			PointVertex pv = value<PointVertex>(*p.surface_, p.surface_cluster_info_, sv);
			
			Scalar dis = value<Scalar>(*p.surface_, p.surface_distance_to_cluster_, sv);
			clusters_color.push_back({dis, value<Vec4>(*p.clusters_, p.clusters_surface_color_, pv).head<3>()});
			for (auto& cluster :
				 value<std::set<PointVertex>>(*p.clusters_, p.clusters_neighbours_, pv))
			{
				Vec3& pos = value<Vec3>(*p.surface_, p.surface_vertex_position_, sv);
				Vec3& sphere_center = value<Vec3>(*p.clusters_, p.clusters_position_, cluster);
				Scalar distance =
					(pos - sphere_center).norm() - value<Scalar>(*p.clusters_, p.clusters_radius_, cluster);
				
				if (std::fabs(dis - distance) <= p.fuzzy_distance_)
				{
					value<std::vector<std::pair<Scalar, Vec3>>>(*p.surface_, p.clusters_fuzzy_cluster_color_, sv)
						.push_back({distance, value<Vec4>(*p.clusters_, p.clusters_surface_color_, cluster).head<3>()});
					value<std::vector<SurfaceVertex>>(*p.clusters_, p.clusters_fuzzy_surface_vertices_, cluster).push_back(sv);
				}
			}
			return true;
		});

		foreach_cell(*p.surface_, [&](SurfaceVertex sv) {
			Vec3 color = Vec3(0, 0, 0);
			Scalar w = 0;
			for (auto& [weight, c] :
				 value<std::vector<std::pair<Scalar, Vec3>>>(*p.surface_, p.clusters_fuzzy_cluster_color_, sv))
			{
			
					w += 1 / weight;
					color += 1 / weight * c;
				
			}
			color /= w;
			value<Vec3>(*p.surface_, p.surface_vertex_color_, sv) = color;
			return true;
		});
	}

	void update_clusters_data(ClusterAxisParameter& p)
	{
		p.huasdorff_distance_cluster_ = 0.0;
		parallel_foreach_cell(*p.clusters_, [&](PointVertex pv) {
			Scalar error = 0.0;
			Scalar cluster_max_error = 0.0;
			SurfaceVertex max_error_vertex;
			const Vec3& center = value<Vec3>(*p.clusters_, p.clusters_position_, pv);
			Scalar radius = value<Scalar>(*p.clusters_, p.clusters_radius_, pv);
			Scalar weight = 0;
			const std::vector<SurfaceVertex>& cluster =
				value<std::vector<SurfaceVertex>>(*p.clusters_, p.clusters_surface_vertices_, pv);
			for (SurfaceVertex sv : cluster)
			{
				Scalar d = value<Scalar>(*p.surface_, p.surface_distance_to_cluster_, sv);
				Scalar area = value<Scalar>(*p.surface_, p.surface_vertex_area_, sv);
				if (cluster_max_error < d)
				{
					cluster_max_error = d;
					max_error_vertex = sv;
				}
				error += d *area;
				weight += area;
			}
			value<Scalar>(*p.clusters_, p.clusters_error_, pv) = error;
			value<Scalar>(*p.clusters_, p.clusters_max_distance_,
													pv) = cluster_max_error;
			value<SurfaceVertex>(*p.clusters_, p.clusters_max_vertex_, pv) = max_error_vertex;
			Scalar mean_error = error / weight;
			
			Scalar variance = 0.0;
			for (SurfaceVertex sv : cluster)
			{
				Scalar d = value<Scalar>(*p.surface_, p.surface_distance_to_cluster_, sv);
				Scalar area = value<Scalar>(*p.surface_, p.surface_vertex_area_, sv);
				variance += area * (d - mean_error) * (d - mean_error);
			}
			Scalar denominator = ((Scalar)(cluster.size()-1) / (Scalar)cluster.size()) * weight;
			variance /= denominator;
			Scalar deviation = std::sqrt(variance);
			value<Scalar>(*p.clusters_, p.clusters_deviation_, pv) = deviation;
			
			
			return true;
		});
		
		p.min_error_ = std::numeric_limits<Scalar>::max();
		p.max_error_ = std::numeric_limits<Scalar>::min();
		p.total_error_ = 0.0;
		for (Scalar e : *p.clusters_error_)
		{
			p.min_error_ = std::min(p.min_error_, e);
			p.max_error_ = std::max(p.max_error_, e);
			p.total_error_ += e;
		}
		for (Scalar max_distance : *p.clusters_max_distance_)
		{
			p.huasdorff_distance_cluster_ = std::max(p.huasdorff_distance_cluster_, max_distance);
		}
		
	} 
		
	void update_clusters_color(ClusterAxisParameter& p)
	{

		Scalar min = std::numeric_limits<Scalar>::max();
		Scalar max = std::numeric_limits<Scalar>::min();
		for (Scalar e : *p.selected_clusters_error_)
		{
			min = std::min(min, e);
			max = std::max(max, e);
		}
		parallel_foreach_cell(*p.clusters_, [&](PointVertex v) {
			value<Vec4>(*p.clusters_, p.clusters_color_, v) =
				color_map(value<Scalar>(*p.clusters_, p.selected_clusters_error_, v), min, max);
			return true;
		});
	}

	Vec4 color_map(Scalar x, Scalar min, Scalar max)
	{
		x = (x - min) / (max - min);
		x = std::clamp(x, 0.0, 1.0);

		Scalar x2 = 2.0 * x;
		switch (int(std::floor(std::max(0.0, x2 + 1.0))))
		{
		case 0:
			return Vec4(0.0, 0.0, 1.0, 0.5);
		case 1:
			return Vec4(x2, x2, 1.0, 0.5);
		case 2:
			return Vec4(1.0, 2.0 - x2, 2.0 - x2, 0.5);
		}
		return Vec4(1.0, 0.0, 0.0, 0.5);
	}

	void update_clusters(ClusterAxisParameter& p)
	{
		std::lock_guard<std::mutex> lock(p.mutex_);
		auto medial_axis_samples_weight =
			get_or_add_attribute<Scalar, SurfaceVertex>(*p.surface_, "medial_axis_samples_weight");
		auto kmax = get_or_add_attribute<Scalar, SurfaceVertex>(*p.surface_, "kmax");
		uint32 split_count = 0;
		parallel_foreach_cell(*p.surface_, [&](SurfaceVertex sv) {
			value<Vec3>(*p.surface_, p.medial_axis_cloest_point_color_, sv) = Vec3(0, 0, 0);
			return true;
		});
		
		// For each cluster, update it's position and radius based on the surface points
		parallel_foreach_cell(*p.clusters_, [&](PointVertex pv) {
			Vec3 opti_coord;
			Scalar rad;
			switch (p.update_method_)
			{
			case MEAN_MEDIAL_SAMPLE: {
				Vec3 sum_coord = Vec3(0, 0, 0);
				Scalar weight = 0;
				for (SurfaceVertex sv1 :
					 value<std::vector<SurfaceVertex>>(*p.clusters_, p.clusters_surface_vertices_, pv))
				{
					Scalar w = 1;
					sum_coord += value<Vec3>(*p.surface_, p.medial_axis_position_, sv1) * w;
					weight += w;
				}
				opti_coord = sum_coord / weight;
			}
			break;
			case MEAN_SURFACE_SAMPLE: {
				Vec3 sum_coord = Vec3(0, 0, 0);
				Scalar weight = 0;
				for (SurfaceVertex sv1 :
					 value<std::vector<SurfaceVertex>>(*p.clusters_, p.clusters_surface_vertices_, pv))
				{
					Scalar w = 1;
					sum_coord += value<Vec3>(*p.surface_, p.surface_vertex_position_, sv1) * w;
					weight += w;
				}

				opti_coord = sum_coord / weight;
			}
			break;
			case SPHERE_FITTING: {
				if (p.fuzzy_clustering_)
				{
					auto [center, r] = sphere_fitting_algo(
						*p.surface_, value<std::vector<SurfaceVertex>>(*p.clusters_, p.clusters_fuzzy_surface_vertices_, pv));
					
					opti_coord = center;
					rad = r;
				}
				else
				{
					auto [center, r] = sphere_fitting_algo(
						*p.surface_, value<std::vector<SurfaceVertex>>(*p.clusters_, p.clusters_surface_vertices_, pv));
					
					opti_coord = center;
					rad = r;
				}
			}
			break;
			case MEAN_MEDIAL_SAMPLE_WEIGHTED_RAYON: {
				Vec3 sum_coord = Vec3(0, 0, 0);
				Scalar weight = 0;
				for (SurfaceVertex sv1 :
					 value<std::vector<SurfaceVertex>>(*p.clusters_, p.clusters_surface_vertices_, pv))
				{
					Scalar w = 2 * value<Scalar>(*p.surface_, p.medial_axis_radius_, sv1) /*+
						value<Scalar>(*p.surface_, medial_axis_samples_weight, sv1)*/
						;
					sum_coord += value<Vec3>(*p.surface_, p.medial_axis_position_, sv1) * w;
					weight += w;
				}
				opti_coord = sum_coord / weight;
			}
			break;
			case MEAN_MEDIAL_SAMPLE_WEIGHTED_CURVATURE: {
				Vec3 sum_coord = Vec3(0, 0, 0);
				Scalar weight = 0;
				for (SurfaceVertex sv1 :
					 value<std::vector<SurfaceVertex>>(*p.clusters_, p.clusters_surface_vertices_, pv))
				{
					Scalar w = value<Scalar>(*p.surface_, medial_axis_samples_weight, sv1);
					sum_coord += value<Vec3>(*p.surface_, p.surface_vertex_position_, sv1) * w;
					weight += w;
				}
				opti_coord = sum_coord / weight;
			}
			break;
			case ENCLOSING_SPHERE_OF_SURFACE_SAMPLE: {
				std::vector<Point> surface_points;
				if (p.fuzzy_clustering_)
				{
					for (SurfaceVertex sv1 :
						 value<std::vector<SurfaceVertex>>(*p.clusters_, p.clusters_fuzzy_surface_vertices_, pv))
					{
						Vec3 pos = value<Vec3>(*p.surface_, p.surface_vertex_position_, sv1);
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
					rad = min_sphere.radius();
				}
				else
				{
					for (SurfaceVertex sv1 :
						 value<std::vector<SurfaceVertex>>(*p.clusters_, p.clusters_surface_vertices_, pv))
					{
						Vec3 pos = value<Vec3>(*p.surface_, p.surface_vertex_position_, sv1);
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
					rad = min_sphere.radius();
				}
				break;
			}
			case ENCLOSING_SPHERE_OF_MEDIAL_SAMPLE: {
				std::vector<Point> surface_points;
				if (p.fuzzy_clustering_)
				{
					for (SurfaceVertex sv1 :
						 value<std::vector<SurfaceVertex>>(*p.clusters_, p.clusters_fuzzy_surface_vertices_, pv))
					{
						Vec3 pos = value<Vec3>(*p.surface_, p.medial_axis_position_, sv1);
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
					rad = min_sphere.radius();
				}
				else
				{
					for (SurfaceVertex sv1 :
						 value<std::vector<SurfaceVertex>>(*p.clusters_, p.clusters_surface_vertices_, pv))
					{
						Vec3 pos = value<Vec3>(*p.surface_, p.medial_axis_position_, sv1);
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
					rad = min_sphere.radius();
				}
				
				break;
			}
			case MASS_SPRING: {

				Vec3 current_center = value<Vec3>(*p.clusters_, p.clusters_position_, pv);
				Vec3 opt_coord = current_center;
				Vec3 new_coord;
				uint32 iteration = 0;
				Scalar radius = value<Scalar>(*p.clusters_, p.clusters_radius_, pv);
				Eigen::Vector3d f;
				Eigen::Matrix3d H;
				do
				{
					iteration++;
					current_center = opt_coord;
					f.setZero();
					H.setZero();
					Eigen::Matrix3d I = Eigen::Matrix3d::Identity(3, 3);
					for (SurfaceVertex sv1 :
						 value<std::vector<SurfaceVertex>>(*p.clusters_, p.clusters_surface_vertices_, pv))
					{
						Vec3 x01 = value<Vec3>(*p.surface_, p.surface_vertex_position_, sv1) - opt_coord;
						Scalar norm = x01.norm();
						Vec3 x01_normalised = x01.normalized();

						f -= 2 * (norm - radius) * x01_normalised;
						Eigen::Matrix3d x01_m = x01 * x01.transpose() / (norm * norm);
						H += 2 * x01_m + 2 * (1 - radius / norm) * (I - x01_m);
					}

					// Newtons method
					new_coord = opt_coord - H.inverse() * f;
					opt_coord = new_coord;
					std::cout << "old_coord: " << current_center.x() << ", " << current_center.y() << ", "
							  << current_center.z() << ", new_coord: " << new_coord.x() << ", " << new_coord.y() << ", "
							  << new_coord.z() << std::endl;
				} while (iteration < 20 && (new_coord - current_center).norm() > 1e-4);
				std::cout << "------------------------------------------------" << std::endl;
				rad = radius;
				break;
			}
			case SQEM: {
				Vec4 opti_sphere;
				/*	std::cout << "cluster " << index_of(*p.clusters_, svc) << std::endl;*/
				bool find_optimal = p.sqem_helper_->optimal_sphere(
					value<std::vector<SurfaceVertex>>(*p.clusters_, p.clusters_surface_vertices_, pv), opti_sphere);
				if (find_optimal)
				{
					std::cout << "opti_sphere radius: " << opti_sphere[3] << std::endl;
					opti_coord = Vec3(opti_sphere[0], opti_sphere[1], opti_sphere[2]);
					rad = opti_sphere[3];
					break;
				}
				else
				{
					std::cout << "can't find optimal sphere" << std::endl;
				}
				
			}
			break;
			case OPTIMAL_SPHERE_UNDER_CONSTRAINT: {
				opti_coord = value<Vec3>(*p.clusters_, p.clusters_position_, pv);
				std::cout << "before optimization" << opti_coord.x() << " " << opti_coord.y() << " " << opti_coord.z()
						  << std::endl;
				auto [radius, v1, v2] = geometry::non_linear_solver(
					*p.surface_, p.surface_vertex_position_.get(), p.surface_vertex_normal_.get(),
					value<std::vector<SurfaceVertex>>(*p.clusters_, p.clusters_surface_vertices_, pv), opti_coord,
					value<Scalar>(*p.clusters_, p.clusters_radius_, pv));
				value<Vec3>(*p.clusters_, p.clusters_position_, pv) = opti_coord;
				value<Scalar>(*p.clusters_, p.clusters_radius_, pv) = radius;

				break;
			}

			default:
				break;
			}
			Vec3 last_coord = opti_coord;
			std::pair<uint32, Vec3> bvh_res;
			p.surface_bvh_->closest_point(opti_coord, &bvh_res);
			Vec3 closest_point_position = bvh_res.second;
			Vec3 closest_point_dir = (closest_point_position - opti_coord).normalized();
			const Vec3& closest_face_normal =
				value<Vec3>(*p.surface_, p.surface_face_normal_, p.surface_bvh_faces_[bvh_res.first]);
			SurfaceVertex v1;
			std::vector<std::pair<uint32, Scalar>> kdt_res;
			if (!p.surface_kdt_->find_nns(opti_coord,10,&kdt_res))
			{
				std::cout << "Didn't find cloest vertex of update position" << std::endl;
			}
			else
			{
				for (auto& [idx, dist] : kdt_res)
				{
					if (p.no_ball->is_marked(p.surface_kdt_vertices_[idx]))
						continue;
					v1 = p.surface_kdt_vertices_[idx];
					break;
				}
			}
			if (!inside(*(p.inside_tester_.get()), K::Point_3(opti_coord.x(), opti_coord.y(), opti_coord.z())))
			{
			
				if (p.auto_split_outside_spheres_ && split_count<2)
				{
					split_outside_cluster(p, opti_coord, pv);
					split_count++;
					return true;
				}
				else {
					closest_point_dir = -closest_point_dir;
				}
			}

			auto [c, r, v2] = geometry::shrinking_ball_center(
				*p.surface_, closest_point_position, closest_point_dir, p.surface_vertex_position_.get(),
				p.surface_bvh_.get(), p.surface_bvh_faces_, p.surface_kdt_.get(), p.surface_kdt_vertices_);

			if (v2.is_valid())
			{
				value<Scalar>(*p.clusters_, p.correction_error_, pv) = (c - last_coord).norm();

				value<Vec3>(*p.surface_, p.medial_axis_cloest_point_color_, v1) =
					value<Vec4>(*p.clusters_, p.clusters_surface_color_, pv).head<3>();
				value<Vec3>(*p.surface_, p.medial_axis_cloest_point_color_, v2) =
					value<Vec4>(*p.clusters_, p.clusters_surface_color_, pv).head<3>();
				value<Vec3>(*p.clusters_, p.clusters_without_correction_position_, pv) = opti_coord;
				value<Scalar>(*p.clusters_, p.clusters_without_correction_radius_, pv) = rad;

				value<Vec3>(*p.clusters_, p.clusters_position_, pv) = c;
				value<Scalar>(*p.clusters_, p.clusters_radius_, pv) = r;
				value<std::pair<SurfaceVertex, SurfaceVertex>>(*p.clusters_, p.clusters_cloest_sample_, pv) = {v1, v2};
			}
			else
			{
				std::cout << "v2 is invalid" << std::endl;
			}
			
			return true;
		});

		assign_cluster(p);
		if (p.connectivity_surgery)
		{
			compute_clusters_cc_numbers(p);
			compute_correct_clusters_neighbor(p);
		}
		else
		{

			update_clusters_neighbor(p);
		}
		if (p.fuzzy_clustering_)
		{
			assign_fuzzy_clusters(p);
		}
		update_clusters_data(p);
		update_clusters_color(p);
		if (!p.running_)
			update_render_data(p);
	}
	
	void split_outside_cluster(ClusterAxisParameter& p, const Vec3& center, PointVertex cluster)
	{
		std::vector<std::pair<uint32, Scalar>> k_res;
		p.surface_kdt_->find_nns(center, 10, &k_res);

		// look for the closest vertex with maximal angle with the first closest vertex
		Scalar max_angle = 0.0;
		SurfaceVertex v, max_angle_vertex;
		uint32 idx = 0;
		//ensure v has medial ball
		for (idx = 0; idx < k_res.size(); ++idx)
		{
			if (p.no_ball->is_marked(p.surface_kdt_vertices_[k_res[idx].first]))
				continue;
			else
				v = p.surface_kdt_vertices_[k_res[idx].first];
				break;
			
		}
		for(uint32 i = idx+1; i<k_res.size(); ++i)
		{
			if (p.no_ball->is_marked(p.surface_kdt_vertices_[k_res[i].first]))
				continue;
			Scalar a = geometry::angle(p.surface_kdt_->vertex(k_res[idx].first) - center,
									   p.surface_kdt_->vertex(k_res[i].first) - center);
			if (a > max_angle)
			{
				max_angle = a;
				max_angle_vertex = p.surface_kdt_vertices_[k_res[i].first];
			}
		}
		if (!v.is_valid() || !max_angle_vertex.is_valid())
		{
			std::cout << "failed split outside cluster" << std::endl;
			return;
		}
		auto& [v1, v2] = value<std::pair<SurfaceVertex, SurfaceVertex>>(*p.surface_, p.medial_axis_closest_points_, v);
		auto& [v3, v4] = value<std::pair<SurfaceVertex, SurfaceVertex>>(*p.surface_, p.medial_axis_closest_points_, max_angle_vertex);

		value<Vec3>(*p.clusters_, p.clusters_position_, cluster) = value<Vec3>(*p.surface_, p.medial_axis_position_, v1);
		value<Scalar>(*p.clusters_, p.clusters_radius_, cluster) = value<Scalar>(*p.surface_, p.medial_axis_radius_, v1);
		value<std::pair<SurfaceVertex, SurfaceVertex>>(*p.clusters_, p.clusters_cloest_sample_, cluster) = {v1, v2};
		value<Vec3>(*p.surface_, p.medial_axis_cloest_point_color_, v1) =
			value<Vec4>(*p.clusters_, p.clusters_surface_color_, cluster).head<3>();
		value<Vec3>(*p.surface_, p.medial_axis_cloest_point_color_, v2) =
			value<Vec4>(*p.clusters_, p.clusters_surface_color_, cluster).head<3>();

		PointVertex new_cluster = add_vertex(*p.clusters_);
		value<Vec3>(*p.clusters_, p.clusters_position_, new_cluster) = value<Vec3>(*p.surface_, p.medial_axis_position_, v3);
		value<Scalar>(*p.clusters_, p.clusters_radius_, new_cluster) = value<Scalar>(*p.surface_, p.medial_axis_radius_, v3);

		value<std::pair<SurfaceVertex, SurfaceVertex>>(*p.clusters_, p.clusters_cloest_sample_,
													   new_cluster) = {v3, v4};

		value<Vec4>(*p.clusters_, p.clusters_surface_color_, new_cluster) =
			Vec4(rand() / (double)RAND_MAX, rand() / (double)RAND_MAX, rand() / (double)RAND_MAX,1);
		value<std::vector<SurfaceVertex>>(*p.clusters_, p.clusters_surface_vertices_, new_cluster).clear();
		
		value<Vec3>(*p.surface_, p.medial_axis_cloest_point_color_,
					v3) = value<Vec4>(*p.clusters_, p.clusters_surface_color_, new_cluster).head<3>();
		value<Vec3>(*p.surface_, p.medial_axis_cloest_point_color_,
					v4) = value<Vec4>(*p.clusters_, p.clusters_surface_color_, new_cluster).head<3>();
	}

	void split(ClusterAxisParameter& p, PointVertex& candidate_cluster)
	{
	
		
		SurfaceVertex max_dist_vertex = value<SurfaceVertex>(*p.clusters_, p.clusters_max_vertex_, candidate_cluster);
		Vec3 pos = value<Vec3>(*p.surface_, p.medial_axis_position_, max_dist_vertex);
		PointVertex new_cluster = add_vertex(*p.clusters_);
		value<Vec3>(*p.clusters_, p.clusters_position_, new_cluster) = pos;
		value<Scalar>(*p.clusters_, p.clusters_radius_, new_cluster) =
			value<Scalar>(*p.surface_, p.medial_axis_radius_, max_dist_vertex);
		value<std::pair<SurfaceVertex, SurfaceVertex>>(*p.clusters_, p.clusters_cloest_sample_, new_cluster) =
			value<std::pair<SurfaceVertex, SurfaceVertex>>(*p.surface_, p.medial_axis_closest_points_, max_dist_vertex);
		value<Vec4>(*p.clusters_, p.clusters_surface_color_, new_cluster) =
			Vec4(rand() / (double)RAND_MAX, rand() / (double)RAND_MAX, rand() / (double)RAND_MAX, 1);

	}

	void split_equal(ClusterAxisParameter& p, PointVertex& candidate_cluster)
	{
		auto& [v1, v2] = value<std::pair<SurfaceVertex, SurfaceVertex>>(*p.clusters_, p.clusters_cloest_sample_, candidate_cluster);
		Vec3 pos1 = value<Vec3>(*p.surface_, p.medial_axis_position_, v1);
		Vec3 pos2 = value<Vec3>(*p.surface_, p.medial_axis_position_, v2);
		Scalar radius1 = value<Scalar>(*p.surface_, p.medial_axis_radius_, v1);
		Scalar radius2 = value<Scalar>(*p.surface_, p.medial_axis_radius_, v2);
		PointVertex new_cluster = add_vertex(*p.clusters_);
		value<Vec3>(*p.clusters_, p.clusters_position_, new_cluster) = pos2;
		value<Scalar>(*p.clusters_, p.clusters_radius_, new_cluster) = radius2;
		value<Vec3>(*p.clusters_, p.clusters_position_, candidate_cluster) = pos1;
		value<Scalar>(*p.clusters_, p.clusters_radius_, candidate_cluster) = radius1;

		value<std::pair<SurfaceVertex, SurfaceVertex>>(*p.clusters_, p.clusters_cloest_sample_, candidate_cluster) =
			value<std::pair<SurfaceVertex, SurfaceVertex>>(*p.surface_, p.medial_axis_closest_points_, v1);
		value<std::pair<SurfaceVertex, SurfaceVertex>>(*p.clusters_, p.clusters_cloest_sample_, new_cluster) =
			value<std::pair<SurfaceVertex, SurfaceVertex>>(*p.surface_, p.medial_axis_closest_points_, v2);
		value<Vec4>(*p.clusters_, p.clusters_surface_color_, new_cluster) =
			Vec4(rand() / (double)RAND_MAX, rand() / (double)RAND_MAX, rand() / (double)RAND_MAX, 1);

	}
	void remove_one_cluster(ClusterAxisParameter& p, PointVertex v)
	{

		std::lock_guard<std::mutex> lock(p.mutex_);

		const std::vector<SurfaceVertex>& clusters_vertices = value<std::vector<SurfaceVertex>>(*p.clusters_, p.clusters_surface_vertices_, v);
		for (SurfaceVertex sv : clusters_vertices)
			value<PointVertex>(*p.surface_, p.surface_cluster_info_, sv) = PointVertex();
		remove_vertex(*p.clusters_, v);

		assign_cluster(p);
		if (p.connectivity_surgery)
		{
			compute_clusters_cc_numbers(p);
			compute_correct_clusters_neighbor(p);
		}
		else
		{

			update_clusters_neighbor(p);
		}
		if (p.fuzzy_clustering_)
		{
			assign_fuzzy_clusters(p);
		}
		update_clusters_data(p);
		update_clusters_color(p);
		if (!p.running_)
			update_render_data(p);
	}


	void split_one_cluster(ClusterAxisParameter& p, PointVertex& candidate_cluster)
	{
		std::lock_guard<std::mutex> lock(p.mutex_);
		split(p, candidate_cluster);
		assign_cluster(p);
		if (p.connectivity_surgery)
		{
			compute_clusters_cc_numbers(p);
			compute_correct_clusters_neighbor(p);
		}
		else
		{

			update_clusters_neighbor(p);
		}
		if (p.fuzzy_clustering_)
		{
			assign_fuzzy_clusters(p);
		}
		update_clusters_data(p);
		update_clusters_color(p);
		

		if (!p.running_)
			update_render_data(p);
	
	}

	struct CompareScalar
	{
		bool operator()(const std::pair<PointVertex, Scalar>& a, const std::pair<PointVertex, Scalar>& b)
		{
			return a.second < b.second;
		}
	};
	void split_cluster(ClusterAxisParameter& p)
	{
		std::lock_guard<std::mutex> lock(p.mutex_);
		
		MeshData<POINT>& md_cluster = point_provider_->mesh_data(*p.clusters_);
		CellMarkerStore<POINT, PointVertex> marker(*p.clusters_);
		Scalar threshold;
		uint32 split_count = 0;
		
		std::priority_queue<std::pair<PointVertex, Scalar>, std::vector<std::pair<PointVertex, Scalar>>, CompareScalar>
			pq;
		// Test if the cluster can be split
		switch (p.split_method_)
		{
		case DEVIATION:
			threshold = p.split_variance_threshold_;
			break;
		case DISTANCE_POINT_SPHERE:
			threshold = p.split_distance_threshold_;
			break;
		case DISTANCE_TO_MEDIAL_AXIS:
			threshold = p.split_distance_to_medial_axis_threshold_;
			break;
		case HAUSDORFF_DISTANCE:
			threshold = p.split_hausdorff_distance_threshold_;
			break;
		default:
			break;
		}
		/*foreach_cell(*p.clusters_, [&](PointVertex pv) {
			std::cout << std::setiosflags(std::ios::fixed);
			std::cout << "cluster: " << index_of(*p.clusters_, pv) << ", error : " << std::setprecision(9)
					  << value<Scalar>(*p.clusters_, p.selected_clusters_error_, pv)
			<< std::endl;
			return true;
		});*/
		foreach_cell(*p.clusters_, [&](PointVertex pv) {
			if (marker.is_marked(pv))
			{
				return true;
			}
			Scalar error = value<Scalar>(*p.clusters_, p.selected_clusters_error_, pv);
			if (error > threshold)
			{
				pq.push({pv, error});
			}
			return true;
		});
			
		while (pq.size() > 0 && split_count < p.max_split_number)
		{
			auto [candidate_cluster, error] = pq.top();
			pq.pop();
			if (marker.is_marked(candidate_cluster))
			{
				continue;
			}
			const std::set<PointVertex>& neighbours =
				value<std::set<PointVertex>>(*p.clusters_, p.clusters_neighbours_, candidate_cluster);
			for (const PointVertex& neighbour : neighbours)
			{
				marker.mark(neighbour);
			}
			split(p, candidate_cluster);
			split_count++;
		}

		assign_cluster(p);
		if (p.connectivity_surgery)
		{
			compute_clusters_cc_numbers(p);
			compute_correct_clusters_neighbor(p);
		}
		else
		{

			update_clusters_neighbor(p);
		}
		if (p.fuzzy_clustering_)
		{
			assign_fuzzy_clusters(p);
		}
		update_clusters_data(p);
		update_clusters_color(p);

		if (!p.running_)
			update_render_data(p);
	}
	void start_clusters_update(ClusterAxisParameter& p)
	{
		p.running_ = true;
		launch_thread([&]() {
			while (true)
			{
				update_clusters(p);
				if (p.slow_down_)
					std::this_thread::sleep_for(std::chrono::microseconds(1000000 / p.update_rate_));
				else
					std::this_thread::yield();
				if (p.stopping_)
				{
					p.running_ = false;
					p.stopping_ = false;
					break;
				}
			}
			
		});

		
		app_.start_timer(100, [&]() -> bool { return !p.running_; });
	}

	void stop_clusters_update(ClusterAxisParameter& p)
	{
		p.stopping_ = true;
	}

	void update_render_data(ClusterAxisParameter& p)
	{
		if (p.running_)
		{
			std::lock_guard<std::mutex> lock(p.mutex_);
			point_provider_->emit_connectivity_changed(*p.clusters_);
			point_provider_->emit_attribute_changed(*p.clusters_, p.clusters_position_.get());
			point_provider_->emit_attribute_changed(*p.clusters_, p.clusters_radius_.get());
			point_provider_->emit_attribute_changed(*p.clusters_, p.clusters_surface_color_.get());
			point_provider_->emit_attribute_changed(*p.clusters_, p.clusters_color_.get());
			point_provider_->emit_attribute_changed(*p.clusters_, p.clusters_without_correction_position_.get());
			point_provider_->emit_attribute_changed(*p.clusters_, p.clusters_without_correction_radius_.get());
			compute_skeleton(p);

			surface_provider_->emit_attribute_changed(*p.surface_, p.surface_vertex_color_.get());
			surface_provider_->emit_attribute_changed(*p.surface_, p.surface_distance_to_cluster_.get());
			surface_provider_->emit_attribute_changed(*p.surface_, p.surface_distance_to_enveloppe.get());
			
			nonmanifold_provider_->emit_connectivity_changed(*p.non_manifold_);
			nonmanifold_provider_->emit_attribute_changed(*p.non_manifold_, p.non_manifold_vertex_position_.get());
			
		}
		else
		{
			point_provider_->emit_connectivity_changed(*p.clusters_);
			point_provider_->emit_attribute_changed(*p.clusters_, p.clusters_position_.get());
			point_provider_->emit_attribute_changed(*p.clusters_, p.clusters_radius_.get());
			point_provider_->emit_attribute_changed(*p.clusters_, p.clusters_surface_color_.get());
			point_provider_->emit_attribute_changed(*p.clusters_, p.clusters_color_.get());
			point_provider_->emit_attribute_changed(*p.clusters_, p.clusters_without_correction_position_.get());
			point_provider_->emit_attribute_changed(*p.clusters_, p.clusters_without_correction_radius_.get());
			compute_skeleton(p);
			
			surface_provider_->emit_attribute_changed(*p.surface_, p.surface_vertex_color_.get());
			surface_provider_->emit_attribute_changed(*p.surface_, p.surface_distance_to_cluster_.get());
			surface_provider_->emit_attribute_changed(*p.surface_, p.surface_distance_to_enveloppe.get());
			
			nonmanifold_provider_->emit_connectivity_changed(*p.non_manifold_);
			nonmanifold_provider_->emit_attribute_changed(*p.non_manifold_, p.non_manifold_vertex_position_.get());
		}
	}
	void compute_skeleton(ClusterAxisParameter& p)
	{
		//Reindex the cluster
		uint32 index = 0;
		clear(*p.non_manifold_);
		if (p.draw_enveloppe)
			p.skel_drawer_.clear();
		auto spheres_skeleton_vertex_map =
			add_attribute<NonManifoldVertex, PointVertex>(*p.clusters_, "__spheres_skeleton_vertex_map");
		//Add vertex
		foreach_cell(*p.clusters_, [&](PointVertex pv) -> bool {
			NonManifoldVertex nmv = add_vertex(*p.non_manifold_);
			Vec3& center = value<Vec3>(*p.clusters_, p.clusters_position_, pv);
			Scalar radius = value<Scalar>(*p.clusters_, p.clusters_radius_, pv);
			value<Vec3>(*p.non_manifold_, p.non_manifold_vertex_position_, nmv) = center;
			value<Scalar>(*p.non_manifold_, p.non_manifold_vertex_radius_, nmv) = radius;
			value<Vec4>(*p.non_manifold_, p.non_manifold_sphere_info_, nmv) = Vec4(center.x(), center.y(), center.z(), radius);
			value<std::vector<SurfaceVertex>>(*p.non_manifold_, p.non_manifold_cluster_vertices_, nmv) =
				value<std::vector<SurfaceVertex>>(*p.clusters_, p.clusters_surface_vertices_, pv);
			value<NonManifoldVertex>(*p.clusters_, spheres_skeleton_vertex_map, pv) = nmv;
			if (p.draw_enveloppe)
			{
				p.skel_drawer_.add_vertex(value<Vec3>(*p.clusters_, p.clusters_position_, pv),
										  value<Scalar>(*p.clusters_, p.clusters_radius_, pv));
			}
			return true;
		});

		//Add edge
		std::unordered_map<std::pair<uint32, uint32>, NonManifoldEdge, edge_hash, edge_equal> edge_indices;
		foreach_cell(*p.clusters_, [&](PointVertex pv) -> bool {
			NonManifoldVertex nmv1 = value<NonManifoldVertex>(*p.clusters_, spheres_skeleton_vertex_map, pv);
			const std::set<PointVertex>& neighbors =
				value<std::set<PointVertex>>(*p.clusters_, p.clusters_neighbours_, pv);
			for (const PointVertex& neighbor : neighbors)
			{
				NonManifoldVertex nmv2 = value<NonManifoldVertex>(*p.clusters_, spheres_skeleton_vertex_map, neighbor);
				std::vector<NonManifoldVertex> av = adjacent_vertices_through_edge(*p.non_manifold_, nmv1);
				if (std::find(av.begin(), av.end(), nmv2) == av.end())
				{
					NonManifoldEdge e = add_edge(*p.non_manifold_, nmv1, nmv2);
					edge_indices[{index_of(*p.non_manifold_, nmv1), index_of(*p.non_manifold_, nmv2)}] = e;
					if (p.draw_enveloppe)
					{
					
					p.skel_drawer_.add_edge(value<Vec3>(*p.clusters_, p.clusters_position_, pv),
												value<Scalar>(*p.clusters_, p.clusters_radius_, pv),
												value<Vec3>(*p.clusters_, p.clusters_position_, neighbor),
												value<Scalar>(*p.clusters_, p.clusters_radius_, neighbor));
					}
				}
			}
			return true;
		});

		//Add face
		foreach_cell(*p.clusters_, [&](PointVertex pv) -> bool {
			NonManifoldVertex nmv1 = value<NonManifoldVertex>(*p.clusters_, spheres_skeleton_vertex_map, pv);
			const std::set<PointVertex>& n_pv =
				value<std::set<PointVertex>>(*p.clusters_, p.clusters_neighbours_, pv);
			for (const PointVertex& ne1 : n_pv)
			{
				NonManifoldVertex nmv2 = value<NonManifoldVertex>(*p.clusters_, spheres_skeleton_vertex_map, ne1);
				const std::set<PointVertex>& ne_ne1 =
					value<std::set<PointVertex>>(*p.clusters_, p.clusters_neighbours_, ne1);
				for (const PointVertex& ne2 : ne_ne1)
				{
					if (n_pv.find(ne2) != n_pv.end())
					{
						NonManifoldVertex nmv3 =
							value<NonManifoldVertex>(*p.clusters_, spheres_skeleton_vertex_map, ne2);
						if (index_of(*p.non_manifold_, nmv1) < index_of(*p.non_manifold_, nmv2) &&
							index_of(*p.non_manifold_, nmv2) < index_of(*p.non_manifold_, nmv3))
						{
							std::vector<NonManifoldEdge> edges;
							edges.push_back(
								edge_indices[{index_of(*p.non_manifold_, nmv1), index_of(*p.non_manifold_, nmv2)}]);
							edges.push_back(
								edge_indices[{index_of(*p.non_manifold_, nmv2), index_of(*p.non_manifold_, nmv3)}]);
							edges.push_back(
								edge_indices[{index_of(*p.non_manifold_, nmv1), index_of(*p.non_manifold_, nmv3)}]);
							add_face(*p.non_manifold_, edges);
							if (p.draw_enveloppe)
							{
								p.skel_drawer_.add_triangle(value<Vec3>(*p.clusters_, p.clusters_position_, pv),
														value<Scalar>(*p.clusters_, p.clusters_radius_, pv),
														value<Vec3>(*p.clusters_, p.clusters_position_, ne1),
														value<Scalar>(*p.clusters_, p.clusters_radius_, ne1),
														value<Vec3>(*p.clusters_, p.clusters_position_, ne2),
														value<Scalar>(*p.clusters_, p.clusters_radius_, ne2));
							}
								
						}
					}
				}
			}
			
			return true;
		});
		if (p.compute_enveloppe_distance_)
		{
			modeling::SphereMeshConstructor<SURFACE, NONMANIFOLD> sphere_mesh_constructor(
				*p.surface_, *p.non_manifold_, p.surface_vertex_position_, p.non_manifold_sphere_info_,
				p.non_manifold_cluster_vertices_);
			Scalar max_dist = 0.0;
			foreach_cell(*p.non_manifold_, [&](NonManifoldVertex nv){
				for (SurfaceVertex sv :
					 value<std::vector<SurfaceVertex>>(*p.non_manifold_, p.non_manifold_cluster_vertices_, nv))
				{
					Scalar min_dist = sphere_mesh_constructor.min_distance_to_enveloppe(nv, sv);
					value<Scalar>(*p.surface_, p.surface_distance_to_enveloppe, sv) = min_dist;
					max_dist = std::max(max_dist, min_dist);
				}
				return true;
			});
			p.hausdorff_distance_enveloppe_ = max_dist;

		}
		
		if (p.draw_enveloppe)
		{
			p.skel_drawer_.update();
		}
		remove_attribute<PointVertex>(*p.clusters_, spheres_skeleton_vertex_map);
		
	}

	void init_delaunay(SURFACE& surface)
	{
		
		m.initialized_ = true;
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
		timer_connection_ = boost::synapse::connect<App::timer_tick>(&app_, [this]() {
			ClusterAxisParameter& p = cluster_axis_parameters_[selected_surface_mesh_];
			update_render_data(p);
		});
	}

	void draw(View* view) override
	{
		ClusterAxisParameter& p = cluster_axis_parameters_[selected_surface_mesh_];
		auto& proj_matrix = view->projection_matrix();
		auto& view_matrix = view->modelview_matrix();
		if(p.draw_enveloppe)
			p.skel_drawer_.draw(proj_matrix, view_matrix);
	}

	void key_press_event(View* view, int32 keycode) override{
		ClusterAxisParameter& p= cluster_axis_parameters_[selected_surface_mesh_];
		if (keycode == GLFW_KEY_I || keycode == GLFW_KEY_S || keycode == GLFW_KEY_D)
		{
			int32 x = view->mouse_x();
			int32 y = view->mouse_y();

			rendering::GLVec3d near_d = view->unproject(x, y, 0.0);
			rendering::GLVec3d far_d = view->unproject(x, y, 1.0);
			Vec3 A{near_d.x(), near_d.y(), near_d.z()};
			Vec3 B{far_d.x(), far_d.y(), far_d.z()};

			Vec3 picked_sphere_center;
			foreach_cell(*p.clusters_, [&](PointVertex v) -> bool {
				if (!picked_sphere_.is_valid())
				{
					picked_sphere_ = v;
					picked_sphere_center = value<Vec3>(*p.clusters_, p.clusters_position_, picked_sphere_);
					return true;
				}
				const Vec3& sp = value<Vec3>(*p.clusters_, p.clusters_position_, v);
				// Scalar sr = value<Scalar>(*clusters, clusters_radius_, v);
				if (geometry::squared_distance_line_point(A, B, sp) <
					geometry::squared_distance_line_point(A, B, picked_sphere_center))
				{
					picked_sphere_ = v;
					picked_sphere_center = sp;
				}
				return true;
			});

		}
		if (keycode == GLFW_KEY_S && picked_sphere_.is_valid())
			split_one_cluster(p, picked_sphere_);
		if (keycode == GLFW_KEY_D && picked_sphere_.is_valid())
			remove_one_cluster(p, picked_sphere_);
	}

	

	void popups() override
	{
		ClusterAxisParameter& p = cluster_axis_parameters_[selected_surface_mesh_];
		if (sphere_info_popup_)
			ImGui::OpenPopup("Sphere Info");
		
		if (ImGui::BeginPopup("Sphere Info"))
		{
			if (picked_sphere_.is_valid())
			{
					ImGui::Text("Picked sphere:");
					const Vec3& sp = value<Vec3>(*p.clusters_, p.clusters_position_, picked_sphere_);
					ImGui::Text("Index: %d", index_of(*p.clusters_, picked_sphere_));
					ImGui::Text("Radius: %f", value<Scalar>(*p.clusters_, p.clusters_radius_, picked_sphere_));
			}
			else
			{
					ImGui::Text("No sphere picked");
			}
			ImGui::EndPopup();
		}
	}

	void left_panel() override
	{
		imgui_mesh_selector(surface_provider_, selected_surface_mesh_, "Surface", [&](SURFACE& m) {
			selected_surface_mesh_ = &m;
			surface_provider_->mesh_data(m).outlined_until_ = App::frame_time_ + 1.0;
		});
	

		if (selected_surface_mesh_){	
		
			ClusterAxisParameter& p = cluster_axis_parameters_[selected_surface_mesh_];
			CoverageAxisParameter& c = coverage_axis_parameters_[selected_surface_mesh_];

			imgui_combo_attribute<SurfaceVertex, Vec3>(*selected_surface_mesh_, p.surface_vertex_position_, "Position",
													   [&](const std::shared_ptr<SurfaceAttribute<Vec3>>& attribute) {
														   p.surface_vertex_position_ = attribute;
														   c.surface_vertex_position_ = attribute;
													   });
			if (ImGui::Button("Cluster Axis"))
			{
				cluster_axis_init(*selected_surface_mesh_);
			}
			if (p.initialized_)
			{
				
				ImGui::RadioButton("Init method: Dilation Constant", (int*)&p.init_method_, CONSTANT);
				ImGui::SameLine();
				ImGui::DragFloat("Init Distance", &p.init_distance_, 0.01f,0.0f, 5.0f, "%.2f" );
				ImGui::RadioButton("Init method: Dilation Factor", (int*)&p.init_method_, FACTOR);
				ImGui::SameLine();
				ImGui::DragFloat("Init factor", &p.init_factor_, 0.01f, 0.0f, 5.0f, "%.2f");
				if (ImGui::Button("Initialise Cluster"))
					initialise_cluster(p);

				if (picked_sphere_.is_valid())
				{
					ImGui::Text("Picked sphere:");
					const Vec3& sp = value<Vec3>(*p.clusters_, p.clusters_position_, picked_sphere_);
					ImGui::Text("Index: %d", index_of(*p.clusters_, picked_sphere_));
					ImGui::Text("Radius: %f", value<Scalar>(*p.clusters_, p.clusters_radius_, picked_sphere_));
				}
				else
				{
					ImGui::Text("No sphere picked");
				}
				ImGui::SliderInt("Update time", &(int)update_times, 1, 100);
				ImGui::RadioButton("Average of medial points", (int*)&p.update_method_, MEAN_MEDIAL_SAMPLE);
				ImGui::RadioButton("Average of sample points", (int*)&p.update_method_, MEAN_SURFACE_SAMPLE);
				ImGui::RadioButton("Sphere fitting", (int*)&p.update_method_, SPHERE_FITTING);
				ImGui::RadioButton("Average of medial points weighted with radius", (int*)&p.update_method_,
								   MEAN_MEDIAL_SAMPLE_WEIGHTED_RAYON);
				ImGui::RadioButton("Average of medial points weighted with curvature", (int*)&p.update_method_, MEAN_MEDIAL_SAMPLE_WEIGHTED_CURVATURE);
				ImGui::RadioButton("Update Cluster by minimum enclosing ball of surface points", (int*)&p.update_method_, ENCLOSING_SPHERE_OF_SURFACE_SAMPLE);
				ImGui::RadioButton("Update Cluster by minimum enclosing ball of medial points", (int*)&p.update_method_, ENCLOSING_SPHERE_OF_MEDIAL_SAMPLE);
				ImGui::RadioButton("Mass spring optimization", (int*)&p.update_method_, MASS_SPRING);
				ImGui::RadioButton("SQEM", (int*)&p.update_method_, SQEM);
				ImGui::RadioButton("Non linear optimization", (int*)&p.update_method_, OPTIMAL_SPHERE_UNDER_CONSTRAINT);

				if (ImGui::Button("Update clusters"))
				{
					for (uint32 i = 0; i < update_times; i++)
						update_clusters(p);
				}
				ImGui::Checkbox("slow down update", &p.slow_down_);
				ImGui::SliderInt("Update rate", (int*)&p.update_rate_, 1, 1000);
				if (ImGui::Checkbox("Fuzzy clustering", &p.fuzzy_clustering_))
				{
					if(p.fuzzy_clustering_)
						assign_fuzzy_clusters(p);
				}
				ImGui::DragFloat("Fuzzy distance", &p.fuzzy_distance_, 0.00001f, 0.0f, 1.0f, "%.5f");
				
				if (!p.running_)
				{
					if (ImGui::Button("Start clusters update"))
						start_clusters_update(p);
				}
				else
				{
					if (ImGui::Button("Stop clusters update"))
						stop_clusters_update(p);
				}
				ImGui::Separator();

				if (ImGui::RadioButton("Distance               ", (int*)&p.split_method_, DISTANCE_POINT_SPHERE))
				{
					set_selected_clusters_error(p, p.clusters_error_.get());

				}
				ImGui::SameLine();
				ImGui::DragFloat("distance threshold", &p.split_distance_threshold_, 0.0001f, 0.0f, 0.1f, "%.6f");
								
				if(ImGui::RadioButton("Deviation               ", (int*)&p.split_method_, DEVIATION))
				{
					set_selected_clusters_error(p, p.clusters_deviation_.get());
				}
				ImGui::SameLine();
				ImGui::DragFloat("Deviation threshold", &p.split_variance_threshold_, 0.000001f, 0.0f, 0.1f,
								 "%.6f");
				
				if(ImGui::RadioButton("Correction distance     ", (int*)&p.split_method_, DISTANCE_TO_MEDIAL_AXIS))
				{
					set_selected_clusters_error(p, p.correction_error_.get());
				}
				ImGui::SameLine();
				ImGui::DragFloat("Correction distance threshold", &p.split_distance_to_medial_axis_threshold_,
								 0.000001f, 0.0f, 0.1f, "%.6f");

				if(ImGui::RadioButton("Hausdorff distance ", (int*)&p.split_method_, HAUSDORFF_DISTANCE))
				{
					set_selected_clusters_error(p, p.clusters_max_distance_.get());
				}
				ImGui::SameLine();
				ImGui::DragFloat("Hausdorff Distance threshold", &p.split_hausdorff_distance_threshold_, 0.000001f, 0.0f, 0.1f, "%.6f");

				
				ImGui::Checkbox("Auto split outside spheres", &p.auto_split_outside_spheres_);
				ImGui::DragInt("Max split number", &(int)p.max_split_number, 1, 0, 100);
				ImGui::Checkbox("Connectivity surgery", &p.connectivity_surgery);
				if (ImGui::Button("Split clusters"))
					split_cluster(p);
				ImGui::Separator();
				ImGui::Checkbox("Compute one sided hausdorff distance", &p.compute_enveloppe_distance_);
				if (ImGui::Checkbox("Draw reconstruction", &p.draw_enveloppe))
				{
					update_render_data(p);
				}
			

				ImGui::Text("Total error: %f", p.total_error_);
				ImGui::Text("Min error: %f", p.min_error_);
				ImGui::Text("Max error: %f", p.max_error_);
				MeshData<SURFACE>& md = surface_provider_->mesh_data(*selected_surface_mesh_);
				ImGui::Text("One sided hausdorff distance to enveloppe: %f \% ",
							p.hausdorff_distance_enveloppe_ / md.diangonal_length() * 100);
				ImGui::Text("One sided hausdorff distance to sphere: %f \% ",
							p.huasdorff_distance_cluster_ / md.diangonal_length() * 100);
				ImGui::Separator();

				ImGui::Text("Pick the sphere under the mouse with I, split it with S and delete it with D");
				if (picked_sphere_.is_valid())
				{
					ImGui::Text("Picked sphere:");
					const Vec3& sp = value<Vec3>(*p.clusters_, p.clusters_position_, picked_sphere_);
					ImGui::Text("Center: (%f, %f, %f)", sp[0], sp[1], sp[2]);
					ImGui::Text("Radius: %f", value<Scalar>(*p.clusters_, p.clusters_radius_, picked_sphere_));
				}
				
			}

			ImGui::Separator();
			if (ImGui::Button("Coverage Axis++"))
			{
				init_coverage_axis_plus_plus(*selected_surface_mesh_);
			}
			if (c.initialized_)
			{
				ImGui::RadioButton("Generate Random candidates", (int*)&c.candidate_generation_method, RANDOM);
				ImGui::RadioButton("Generate Shrinking ball candidates", (int*)&c.candidate_generation_method, SHRINKING_BALL);
				ImGui::RadioButton("Generate Delaunay candidates", (int*)&c.candidate_generation_method, DELAUNAY);
				if (ImGui::Button("Generate candidates"))
				{
					generate_candidates(c);
				}
				if (c.candidates_valid)
				{
					ImGui::DragFloat("Dilation factor", &c.dilation_factor, 0.001f, 0.0f, 1.0f, "%.3f");
					ImGui::DragInt("Max Vertices", (int*)&c.max_selected_number,1 ,0, 500);
					if (ImGui::Button("Compute Coverage axis++"))
					{
						coverage_axis_plus_plus(c);
					}
					if (ImGui::Button("Generate connectivity by Qmat"))
					{
						coverage_axis_plus_plus_connectivity(c);
					}
				}
			}
			ImGui::Separator();
			
			if (ImGui::Button("Coverage_axis"))
			{
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
			}
			if (ImGui::Button("Delaunay based method"))
			{

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
				if (ImGui::Button("Coverage Axis"))
				{
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
				}
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
										 min_radius_, max_radius_, "%.3f");
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
					imgui_mesh_selector(
						nonmanifold_provider_, selected_medial_axis_, "Medial_axis", [&](NONMANIFOLD& nm) {
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
								coverage_axis_PD(*selected_surface_mesh_, *selected_candidates_, solution,
												 dilation_factor);
							if (ImGui::Button("RVD"))
								coverage_axis_RVD(*selected_surface_mesh_, *selected_candidates_, solution);
						}
					}
				}
			}
		}
	}
	

private :
	
	SURFACE* selected_surface_mesh_ = nullptr;
	NONMANIFOLD* selected_medial_axis_ = nullptr;
	POINT* selected_candidates_ = nullptr;
	POINT* surface_sample_ = nullptr;
	std::unordered_map<const SURFACE*, ClusterAxisParameter> cluster_axis_parameters_;
	std::unordered_map<const SURFACE*, CoverageAxisParameter> coverage_axis_parameters_;
	PointVertex picked_sphere_;
	bool sphere_info_popup_ = false;

	MeshProvider<POINT>* point_provider_;
	MeshProvider<SURFACE>* surface_provider_;
	MeshProvider<NONMANIFOLD>* nonmanifold_provider_;
	HighsSolution solution;
	Delaunay tri_;

	Cgal_Surface_mesh csm;
	Tree* tree;
	Point_inside* inside_tester;
	std::array<std::mutex, 43> spheres_mutex_;
	
	std::vector<Vec3> colors;
	
	bool angle_filtering_ = true;
	bool circumradius_filtering_ = true;
	bool distance_filtering_ = true;
	bool pole_filtering_ = true;
	float distance_threshold_ = 0.001;
	float angle_threshold_ = 1.9;
	float radius_threshold_ = 0.030;
	float dilation_factor = 0.1f;
	uint32 update_times = 5;
	double min_radius_ = std::numeric_limits<double>::max();
	double max_radius_ = std::numeric_limits<double>::min();
	double min_angle_ = std::numeric_limits<double>::max();
	double max_angle_ = std::numeric_limits<double>::min();

	std::shared_ptr<boost::synapse::connection> timer_connection_;
	};

} // namespace ui

} // namespace cgogn

#endif // CGOGN_MODULE_POWER_SHAPE_H_
