
#ifndef CGOGN_MODULE_COVERAGE_AXIS_H_
#define CGOGN_MODULE_COVERAGE_AXIS_H_
#include <cgogn/core/types/maps/gmap/gmap_base.h>
#include <cgogn/core/ui_modules/mesh_provider.h>
#include <cgogn/geometry/types/slab_quadric.h>
#include <cgogn/modeling/skeleton_sampling.h>
#include <cgogn/ui/app.h>
#include <cgogn/ui/module.h>
#include <cmath>
#include <filesystem>
#include <libacc/bvh_trees_spheres.h>

#include <cgogn/core/types/maps/cmap/cmap0.h>
#include <cgogn/geometry/algos/medial_axis.h>
#include <cgogn/io/point/point_import.h>
#include <cgogn/modeling/algos/decimation/QEM_helper.h>
#include <cgogn/modeling/algos/decimation/SQEM_helper.h>
#include <cgogn/rendering/skelshape.h>
#include <libacc/bvh_tree.h>
#include <libacc/kd_tree.h>
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
#include <CGAL/double.h>
#include <CGAL/optimal_bounding_box.h>
#include <CGAL/point_generators_3.h>

#include <GLFW/glfw3.h>
#include <Highs.h>

#include <iomanip>
#include <limits>

namespace cgogn
{

namespace ui
{
template <typename POINT, typename SURFACE, typename NONMANIFOLD>
class CoverageAxis : public ViewModule
{
	static_assert(mesh_traits<SURFACE>::dimension >= 2, "Coverage axis can only be used with meshes of dimension >= 2");
	// Kernel for construct Delaunay
	using K = CGAL::Exact_predicates_inexact_constructions_kernel;
	using Point = K::Point_3;
	using Weight_Point = K::Weighted_point_3;
	enum CandidateGenerationMethod
	{
		RANDOM,
		SHRINKING_BALL,
		DELAUNAY
	};
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
	CoverageAxis(const App& app) : ViewModule(app, "CoverageAxis")
	{
	}
	~CoverageAxis()
	{
		
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
public:
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
		std::shared_ptr<PointAttribute<Scalar>> candidate_points_centrality_score_ = nullptr;
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

		std::vector<SurfaceFace> surface_bvh_faces_;
		std::vector<SurfaceVertex> surface_kdt_vertices_;
		std::shared_ptr<acc::BVHTree<uint32, Vec3>> surface_bvh_ = nullptr;
		std::shared_ptr<acc::KDTree<3, uint32>> surface_kdt_ = nullptr;

		float dilation_factor = 0.02;
		CandidateGenerationMethod candidate_generation_method = SHRINKING_BALL;

		uint32 surface_samples_number = 1500;
		uint32 candidates_number = 20000;
		uint32 max_selected_number = 200;
		std::vector<PointVertex> selected_inner_points;
		std::vector<PointVertex> candidate_inner_points;
		std::vector<Vec3> sample_points;

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

	NONMANIFOLD& constrcut_power_shape_non_manifold(Regular& power_shape, string name)
	{
		NONMANIFOLD* mv =
			nonmanifold_provider_->add_mesh(std::to_string(nonmanifold_provider_->number_of_meshes()) + "_" + name);
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
				Inner_Power_shape_data.vertex_radius_.push_back(std::sqrt(vit->point().weight()));
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
					v[count] = fit->first->vertex(idx);
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
		auto sphere_info = add_attribute<Vec4, NonManifoldVertex>(*mv, "sphere_info");
		for (uint32 idx = 0u; idx < power_point.size(); ++idx)
		{
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
	void init_coverage_axis_plus_plus(SURFACE& surface)
	{
		CoverageAxisParameter& p = coverage_axis_parameters_[&surface];
		p.surface_ = &surface;
		p.initialized_ = true;
		// Surface attributes
		p.medial_axis_position_ = get_or_add_attribute<Vec3, SurfaceVertex>(surface, "medial_axis_position");
		p.medial_axis_radius_ = get_or_add_attribute<Scalar, SurfaceVertex>(surface, "medial_axis_radius");
		p.medial_axis_closest_points_ = get_or_add_attribute<std::pair<SurfaceVertex, SurfaceVertex>, SurfaceVertex>(
			surface, "medial_axis_closest_points");
		p.medial_axis_angle_ = get_or_add_attribute<Scalar, SurfaceVertex>(surface, "medial_axis_angle");
		p.surface_vertex_normal_ = get_or_add_attribute<Vec3, SurfaceVertex>(surface, "normal");
		geometry::compute_normal<SurfaceVertex>(surface, p.surface_vertex_position_.get(),
												p.surface_vertex_normal_.get());

		// Point attributes
		p.candidate_points_ = point_provider_->add_mesh(point_provider_->mesh_name(surface) + "candidates");
		p.candidate_points_position_ = get_or_add_attribute<Vec3, PointVertex>(*p.candidate_points_, "position");
		p.candidate_points_radius_ = get_or_add_attribute<Scalar, PointVertex>(*p.candidate_points_, "radius");
		p.candidate_points_score_ = get_or_add_attribute<Scalar, PointVertex>(*p.candidate_points_, "score");
		p.candidate_points_coverage_score_ =
			get_or_add_attribute<Scalar, PointVertex>(*p.candidate_points_, "coverage_score");
		p.candidate_points_uniformity_score_ =
			get_or_add_attribute<Scalar, PointVertex>(*p.candidate_points_, "uniformity_score");
		p.candidate_points_centrality_score_ =
			get_or_add_attribute<Scalar, PointVertex>(*p.candidate_points_, "centrality_score");
		p.candidate_points_associated_vertex_ =
			get_or_add_attribute<NonManifoldVertex, PointVertex>(*p.candidate_points_, "associated_vertex");

		p.selected_points_ = point_provider_->add_mesh(point_provider_->mesh_name(surface) + "selected");
		p.selected_points_position_ = get_or_add_attribute<Vec3, PointVertex>(*p.selected_points_, "position");
		p.selected_points_radius_ = get_or_add_attribute<Scalar, PointVertex>(*p.selected_points_, "radius");
		p.selected_points_associated_vertex_ =
			get_or_add_attribute<NonManifoldVertex, PointVertex>(*p.selected_points_, "associated_vertex");

		// Non-manifold attributes
		p.voronoi_ = nonmanifold_provider_->add_mesh(surface_provider_->mesh_name(surface) + "_coverage_axis");
		p.voronoi_position_ = get_or_add_attribute<Vec3, NonManifoldVertex>(*p.voronoi_, "position");
		p.voronoi_radius_ = get_or_add_attribute<Scalar, NonManifoldVertex>(*p.voronoi_, "radius");
		p.voronoi_stability_ratio_ = get_or_add_attribute<Scalar, NonManifoldEdge>(*p.voronoi_, "stability_ratio");
		p.voronoi_stability_color_ = get_or_add_attribute<Vec3, NonManifoldEdge>(*p.voronoi_, "stability_color");
		p.voronoi_sphere_info_ = get_or_add_attribute<Vec4, NonManifoldVertex>(*p.voronoi_, "sphere_info");
		p.voronoi_fixed_vertices_ = get_or_add_attribute<bool, NonManifoldVertex>(*p.voronoi_, "fixed_vertices");
		MeshData<SURFACE>& md = surface_provider_->mesh_data(surface);
		// Build bvh and kdt for surface
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
		p.surface_kdt_ = std::make_shared<acc::KDTree<3, uint32>>(vertex_position_vector);
		// Compute initalize non-manifold for Qmat

		load_model_in_cgal(surface, p.csm_);
		p.tree_ = std::make_shared<Tree>(faces(p.csm_).first, faces(p.csm_).second, p.csm_);
		p.tree_->accelerate_distance_queries();
		p.inside_tester_ = std::make_shared<Point_inside>(*(p.tree_.get()));

		// Construct non-manifold
		compute_initial_non_manifold(p);
		nb_vertices = nb_cells<NonManifoldVertex>(*p.voronoi_);
		// build Kd-tree for the voronoi diagram
		p.voronoi_kdt_vertices_.clear();
		p.voronoi_kdt_vertices_.reserve(nb_vertices);
		std::vector<Vec3> voronoi_vertex_position_vector;
		voronoi_vertex_position_vector.reserve(nb_vertices);
		idx = 0;
		foreach_cell(*p.voronoi_, [&](NonManifoldVertex v) -> bool {
			p.voronoi_kdt_vertices_.push_back(v);
			voronoi_vertex_position_vector.push_back(value<Vec3>(*p.voronoi_, p.voronoi_position_, v));
			return true;
		});
		p.voronoi_kdt_ = std::make_shared<acc::KDTree<3, uint32>>(voronoi_vertex_position_vector);
	}
	void random_sampling(CoverageAxisParameter& p)
	{
		Vec3 bias = Vec3(0.5, 0.5, 0.5);
		std::vector<Vec3> samples;
		CGAL::Random_points_in_cube_3<Point> generator(0.5);
		while (samples.size() < p.candidates_number)
		{
			Point s = *generator++;
			Vec3 pos = Vec3(s.x(), s.y(), s.z()) + bias;
			if (inside(*(p.inside_tester_.get()), Point(pos.x(), pos.y(), pos.z())))
			{
				samples.push_back(pos);
			}
		}
		for (auto& pos : samples)
		{
			PointVertex new_candidate = add_vertex(*p.candidate_points_);
			p.candidate_inner_points.push_back(new_candidate);
			std::pair<uint32, Vec3> bvh_res;
			p.surface_bvh_->closest_point(pos, &bvh_res);
			Vec3 closest_surface_position = bvh_res.second;
			value<Vec3>(*p.candidate_points_, p.candidate_points_position_, new_candidate) = pos;
			value<Scalar>(*p.candidate_points_, p.candidate_points_radius_, new_candidate) =
				(pos - closest_surface_position).norm();
			// bind non-manifold vertex to candidate point
			std::pair<uint32, Scalar> k_res;
			if (!p.voronoi_kdt_->find_nn(pos, &k_res))
			{
				std::cout << "closest point not found !!!";
				continue;
			}
			value<NonManifoldVertex>(*p.candidate_points_, p.candidate_points_associated_vertex_, new_candidate) =
				p.voronoi_kdt_vertices_[k_res.first];
		}
		// TODO
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
			// bind non-manifold vertex to candidate point
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
		clear(*p.candidate_points_);
		switch (p.candidate_generation_method)
		{
		case SHRINKING_BALL: {
			shrinking_ball_sampling(p);
		}
		break;
		case RANDOM: {
			random_sampling(p);
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
		p.coverage_matrix.resize(p.surface_samples_number, nb_cells<PointVertex>(*p.candidate_points_));
		p.sample_points.clear();
		std::vector<Point> points;
		CGAL::Random_points_in_triangle_mesh_3<Cgal_Surface_mesh> generator(p.csm_);
		std::copy_n(generator, p.surface_samples_number, std::back_inserter(points));
		for (auto& pos : points)
		{
			p.sample_points.push_back(Vec3(pos.x(), pos.y(), pos.z()));
		}
		p.coverage_matrix.setZero();
		uint32 idx = 0;
		for (auto& pos : p.sample_points)
		{
			foreach_cell(*p.candidate_points_, [&](PointVertex candidate) {
				Vec3 center = value<Vec3>(*p.candidate_points_, p.candidate_points_position_, candidate);
				Scalar radius = value<Scalar>(*p.candidate_points_, p.candidate_points_radius_, candidate);
				if (inside_sphere(pos, center, radius + p.dilation_factor))
				{
					p.coverage_matrix(idx, index_of(*p.candidate_points_, candidate)) = 1;
				}
				return true;
			});
			idx++;
		}
	}

	Scalar compute_centrality_score(CoverageAxisParameter& p, PointVertex pv)
	{
		return -1 / value<Scalar>(*p.candidate_points_, p.candidate_points_radius_, pv);
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
		Scalar centrality_sum = 0;
		for (PointVertex& p_minus : candidate_points)
		{
			uint32 v_index = index_of(*p.candidate_points_, p_minus);
			(*p.candidate_points_coverage_score_)[v_index] = compute_coverage_score(p, p_minus);
			coverage_sum += (*p.candidate_points_coverage_score_)[v_index];
			(*p.candidate_points_uniformity_score_)[v_index] = compute_uniform_score(p, p_minus);
			uniform_sum += (*p.candidate_points_uniformity_score_)[v_index];
			(*p.candidate_points_centrality_score_)[v_index] = compute_centrality_score(p, p_minus);
			centrality_sum += (*p.candidate_points_centrality_score_)[v_index];
		}
		Scalar mean_coverage = coverage_sum / candidate_points.size();
		Scalar mean_uniform = uniform_sum / candidate_points.size();
		Scalar mean_centrality = centrality_sum / candidate_points.size();
		Scalar deviation_coverage = 0;
		Scalar deviation_uniform = 0;
		Scalar deviation_centrality = 0;

		for (PointVertex& p_minus : candidate_points)
		{
			uint32 v_index = index_of(*p.candidate_points_, p_minus);
			Scalar& coverage_score = (*p.candidate_points_coverage_score_)[v_index];
			Scalar& uniform_score = (*p.candidate_points_uniformity_score_)[v_index];
			Scalar& centrality_score = (*p.candidate_points_centrality_score_)[v_index];
			deviation_coverage += (coverage_score - mean_coverage) * (coverage_score - mean_coverage);
			deviation_uniform += (uniform_score - mean_uniform) * (uniform_score - mean_uniform);
			deviation_centrality += (centrality_score - mean_centrality) * (centrality_score - mean_centrality);
		}
		deviation_coverage = sqrt(deviation_coverage);
		deviation_uniform = sqrt(deviation_uniform);
		deviation_centrality = sqrt(deviation_centrality);

		for (PointVertex& p_minus : candidate_points)
		{
			uint32 v_index = index_of(*p.candidate_points_, p_minus);
			Scalar coverage_normalised =
				((*p.candidate_points_coverage_score_)[v_index] - mean_coverage) / deviation_coverage;
			Scalar uniform_normalised =
				((*p.candidate_points_uniformity_score_)[v_index] - mean_uniform) / deviation_uniform;
			Scalar centrality_normalised =
				((*p.candidate_points_centrality_score_)[v_index] - mean_centrality) / deviation_centrality;
			value<Scalar>(*p.candidate_points_, p.candidate_points_score_, p_minus) =
				coverage_normalised + uniform_normalised + centrality_normalised;
		}
	}

	void coverage_axis(CoverageAxisParameter& p)
	{
		typedef Eigen::SparseMatrix<Scalar> SpMat;
		typedef Eigen::Triplet<Scalar> T;
		std::vector<T> triplets;
		clear(*p.selected_points_);
		p.selected_inner_points.clear();
		std::vector<Vec3> sample_points;
		std::vector<Vec3> inner_points;
		std::vector<Scalar> weights;
		std::vector<PointVertex> candidates;

		std::vector<Point> points;
		CGAL::Random_points_in_triangle_mesh_3<Cgal_Surface_mesh> generator(p.csm_);
		std::copy_n(generator, p.surface_samples_number, std::back_inserter(points));
		for (auto& pos : points)
		{
			sample_points.push_back(Vec3(pos.x(), pos.y(), pos.z()));
		}
		foreach_cell(*p.candidate_points_, [&](PointVertex pv) {
			candidates.push_back(pv);
			inner_points.push_back(value<Vec3>(*p.candidate_points_, p.candidate_points_position_, pv));
			weights.push_back(value<Scalar>(*p.candidate_points_, p.candidate_points_radius_, pv));
			return true;
		});
		uint32 nb_candidates = inner_points.size();
		uint32 nb_surface_samples = sample_points.size();
		for (size_t idx1 = 0; idx1 < nb_surface_samples; idx1++)
		{
			for (size_t idx2 = 0; idx2 < nb_candidates; idx2++)
			{
				if (inside_sphere(sample_points[idx1], inner_points[idx2], weights[idx2] + p.dilation_factor))
				{
					triplets.push_back(T(idx1, idx2, 1.0));
				}
			}
		}

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
		highs.setOptionValue("time_limit", 1200);
		highs.setOptionValue("parallel", "on");
		/*highs.setHighsOptionValue("solver", "simplex");*/
		if (status == HighsStatus::kOk)
		{
			highs.run();

			assert(status == HighsStatus::kOk);
			solution = highs.getSolution();
		}
		for (size_t idx = 0; idx < nb_candidates; idx++)
		{
			if (solution.col_value[idx] > 0.5)
			{
				PointVertex new_selected = add_vertex(*p.selected_points_);
				value<Vec3>(*p.selected_points_, p.selected_points_position_, new_selected) = inner_points[idx];
				value<Scalar>(*p.selected_points_, p.selected_points_radius_, new_selected) = weights[idx];
				value<NonManifoldVertex>(*p.selected_points_, p.selected_points_associated_vertex_, new_selected) =
					value<NonManifoldVertex>(*p.candidate_points_, p.candidate_points_associated_vertex_,
											 candidates[idx]);
			}
		}
		point_provider_->emit_connectivity_changed(*p.selected_points_);
		point_provider_->emit_attribute_changed(*p.selected_points_, p.selected_points_position_.get());
		point_provider_->emit_attribute_changed(*p.selected_points_, p.selected_points_radius_.get());
	}

	void coverage_axis_plus_plus(CoverageAxisParameter& p)
	{
		clear(*p.selected_points_);
		compute_coverage_matrix(p);
		p.selected_inner_points.clear();
		auto copie_candidates = std::vector<PointVertex>(p.candidate_inner_points);
		uint32 covered_number = 0;
		while (covered_number < p.sample_points.size() && p.selected_inner_points.size() < p.max_selected_number)
		{
			compute_score(p, copie_candidates);
			auto& max_element_it = std::max_element(
				copie_candidates.begin(), copie_candidates.end(), [&p](PointVertex& v1, PointVertex& v2) {
					return value<Scalar>(*p.candidate_points_, p.candidate_points_score_, v1) <
						   value<Scalar>(*p.candidate_points_, p.candidate_points_score_, v2);
				});
			if (max_element_it != copie_candidates.end())
			{
				PointVertex max_score_vertex = *max_element_it;
				copie_candidates.erase(max_element_it);
				p.selected_inner_points.push_back(max_score_vertex);

				// Update coverage matrix;
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
					value<NonManifoldVertex>(*p.candidate_points_, p.candidate_points_associated_vertex_,
											 max_score_vertex);
			}
		}
		point_provider_->emit_connectivity_changed(*p.selected_points_);
		point_provider_->emit_attribute_changed(*p.selected_points_, p.selected_points_position_.get());
		point_provider_->emit_attribute_changed(*p.selected_points_, p.selected_points_radius_.get());
	}

	void coverage_axis_connect_by_power_diagram(CoverageAxisParameter& p)
	{
		Regular power_shape;
		foreach_cell(*p.selected_points_, [&](PointVertex pv) {
			Vec3 pos = value<Vec3>(*p.selected_points_, p.selected_points_position_, pv);
			Scalar radius = value<Scalar>(*p.selected_points_, p.selected_points_radius_, pv);
			power_shape.insert(Weight_Point(Point(pos[0], pos[1], pos[2]),
											(radius + p.dilation_factor) * (radius + p.dilation_factor)));

			return true;
		});
		foreach_cell(*p.surface_, [&](SurfaceVertex v) {
			Vec3 pos = value<Vec3>(*p.surface_, p.surface_vertex_position_, v);
			power_shape.insert(Weight_Point(Point(pos[0], pos[1], pos[2]), p.dilation_factor * p.dilation_factor));
			return true;
		});
		std::cout << "determine if the point is outside" << std::endl;
		int count = 0;
		for (auto vit = power_shape.finite_vertices_begin(); vit != power_shape.finite_vertices_end(); ++vit)
		{
			// if the point is inside
			Point po = vit->point().point();
			if (inside(*(p.inside_tester_.get()), po))
			{
				vit->info().id = count;
				vit->info().inside = true;
				count++;
			}
		}
		std::cout << "construct non manifold PD" << std::endl;
		constrcut_power_shape_non_manifold(power_shape, "_coverage_axis" + surface_provider_->mesh_name(*p.surface_));
	}
	void coverage_axis_connect_by_Qmat(CoverageAxisParameter& p)
	{
		using QMatHelper = modeling::DecimationSQEM_Helper<NONMANIFOLD>;
		using Slab_Quadric = geometry::Slab_Quadric;

		foreach_cell(*p.voronoi_, [&](NonManifoldEdge e) -> bool {
			auto iv = incident_vertices(*p.voronoi_, e);
			NonManifoldVertex v1 = iv[0];
			NonManifoldVertex v2 = iv[1];
			const Vec3& v1_p = value<Vec4>(*p.voronoi_, p.voronoi_sphere_info_, v1).head<3>();
			const Vec3& v2_p = value<Vec4>(*p.voronoi_, p.voronoi_sphere_info_, v2).head<3>();
			const Scalar& r1 = value<Vec4>(*p.voronoi_, p.voronoi_sphere_info_, v1).w();
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
			NonManifoldVertex nv =
				value<NonManifoldVertex>(*p.selected_points_, p.selected_points_associated_vertex_, pv);
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
	void interpolate_skeleton(NONMANIFOLD& non_manifold)
	{
		auto& non_manifold_vertex_position = get_or_add_attribute<Vec3, NonManifoldVertex>(non_manifold, "position");
		auto& non_manifold_vertex_radius = get_or_add_attribute<Scalar, NonManifoldVertex>(non_manifold, "radius");
		skeleton_drawer.clear();
		skeleton_drawer.set_color({1.0, 1.0, 1.0, 0.5});
		skeleton_drawer.set_subdiv(40);
		parallel_foreach_cell(non_manifold, [&](NonManifoldVertex nv) {
			skeleton_drawer.add_vertex(value<Vec3>(non_manifold, non_manifold_vertex_position, nv),
									   value<Scalar>(non_manifold, non_manifold_vertex_radius, nv));
			return true;
		});
		parallel_foreach_cell(non_manifold, [&](NonManifoldEdge ne) {
			auto& v_vec = incident_vertices(non_manifold, ne);
			skeleton_drawer.add_edge(value<Vec3>(non_manifold, non_manifold_vertex_position, v_vec[0]),
									 value<Scalar>(non_manifold, non_manifold_vertex_radius, v_vec[0]),
									 value<Vec3>(non_manifold, non_manifold_vertex_position, v_vec[1]),
									 value<Scalar>(non_manifold, non_manifold_vertex_radius, v_vec[1]));
			return true;
		});
		parallel_foreach_cell(non_manifold, [&](NonManifoldFace nf) {
			auto& v_vec = incident_vertices(non_manifold, nf);
			skeleton_drawer.add_triangle(value<Vec3>(non_manifold, non_manifold_vertex_position, v_vec[0]),
										 value<Scalar>(non_manifold, non_manifold_vertex_radius, v_vec[0]),
										 value<Vec3>(non_manifold, non_manifold_vertex_position, v_vec[1]),
										 value<Scalar>(non_manifold, non_manifold_vertex_radius, v_vec[1]),
										 value<Vec3>(non_manifold, non_manifold_vertex_position, v_vec[2]),
										 value<Scalar>(non_manifold, non_manifold_vertex_radius, v_vec[2]));
			return true;
		});
		skeleton_drawer.update();
	}
	void export_spheres_OBJ(CoverageAxisParameter& p)
	{
		auto& mesh_name = surface_provider_->mesh_name(*p.surface_);
		std::ofstream file(mesh_name + "_spheres.obj");
		if (!file.is_open())
		{
			std::cerr << "Error opening file" << std::endl;
			return;
		}
		foreach_cell(*p.selected_inner_points, [&](PointVertex pv) {
			Vec3 center = value<Vec3>(*p.selected_inner_points, p.selected_points_position_, pv);
			Scalar radius = value<Scalar>(*p.selected_inner_points, p.selected_points_radius_, pv);
			file << "v " << center.x() << " " << center.y() << " " << center.z() << " " << radius << std::endl;
			return true;
		});
		file.close();
	}

	acc::BVHTree<uint32, Vec3>* build_bvh(SURFACE& surface)
	{
		auto surface_vertex_position = get_or_add_attribute<Vec3, SurfaceVertex>(surface, "position");
		auto bvh_vertex_index = add_attribute<uint32, SurfaceVertex>(surface, "__bvh_vertex_index");

		std::vector<Vec3> vertex_position_vector;
		vertex_position_vector.reserve(nb_cells<SurfaceVertex>(surface));
		uint32 idx = 0;
		foreach_cell(surface, [&](SurfaceVertex v) -> bool {
			value<uint32>(surface, bvh_vertex_index, v) = idx++;
			vertex_position_vector.push_back(value<Vec3>(surface, surface_vertex_position, v));
			return true;
		});

		vector<SurfaceFace> surface_bvh_faces_;
		uint32 nb_faces = nb_cells<SurfaceFace>(surface);
		surface_bvh_faces_.reserve(nb_faces);
		std::vector<uint32> face_vertex_indices;
		face_vertex_indices.reserve(nb_faces * 3);
		foreach_cell(surface, [&](SurfaceFace f) -> bool {
			surface_bvh_faces_.push_back(f);
			foreach_incident_vertex(surface, f, [&](SurfaceVertex v) -> bool {
				face_vertex_indices.push_back(value<uint32>(surface, bvh_vertex_index, v));
				return true;
			});
			return true;
		});

		acc::BVHTree<uint32, Vec3>* surface_bvh_ =
			new acc::BVHTree<uint32, Vec3>(face_vertex_indices, vertex_position_vector);
		remove_attribute<SurfaceVertex>(surface, bvh_vertex_index);
		return surface_bvh_;
	}

	void compute_hausdorff_distance(NONMANIFOLD& non_manifold, SURFACE& surface)
	{
		auto& surface_vertex_position = get_or_add_attribute<Vec3, SurfaceVertex>(surface, "position");

		auto& non_manifold_vertex_position = get_or_add_attribute<Vec3, NonManifoldVertex>(non_manifold, "position");
		auto& non_manifold_vertex_radius = get_or_add_attribute<Scalar, NonManifoldVertex>(non_manifold, "radius");
		auto& non_manifold_sphere_info =
			get_or_add_attribute<Vec4, NonManifoldVertex>(non_manifold, "non_manifold_sphere_info");
		auto& nonmanifold_cluster_vertices = get_or_add_attribute<std::vector<SurfaceVertex>, NonManifoldVertex>(
			non_manifold, "nonmanifold_cluster_vertices");
		foreach_cell(surface, [&](SurfaceVertex sv) {
			Scalar min_dist = std::numeric_limits<Scalar>::max();
			Vec3& pos = value<Vec3>(surface, surface_vertex_position, sv);
			NonManifoldVertex non_manifold_vertex;
			foreach_cell(non_manifold, [&](NonManifoldVertex nv) {
				Vec3 center = value<Vec3>(non_manifold, non_manifold_vertex_position, nv);
				Scalar radius = value<Scalar>(non_manifold, non_manifold_vertex_radius, nv);
				Scalar dist = (center - pos).norm() - radius;
				if (dist < min_dist)
				{
					min_dist = dist;
					non_manifold_vertex = nv;
				}
				return true;
			});
			value<std::vector<SurfaceVertex>>(non_manifold, nonmanifold_cluster_vertices, non_manifold_vertex)
				.push_back(sv);
			return true;
		});
		foreach_cell(non_manifold, [&](NonManifoldVertex nv) {
			Vec3 center = value<Vec3>(non_manifold, non_manifold_vertex_position, nv);
			Scalar radius = value<Scalar>(non_manifold, non_manifold_vertex_radius, nv);
			value<Vec4>(non_manifold, non_manifold_sphere_info, nv) = Vec4(center.x(), center.y(), center.z(), radius);
			return true;
		});

		// Compute the distance from the enveloppe to the shape
		modeling::SphereMeshConstructor<SURFACE, NONMANIFOLD> sphere_mesh_constructor(
			surface, non_manifold, surface_vertex_position, non_manifold_sphere_info, nonmanifold_cluster_vertices);
		Scalar max_dist_shpae_to_enveloppe = 0.0;
		foreach_cell(non_manifold, [&](NonManifoldVertex nv) {
			for (SurfaceVertex sv : value<std::vector<SurfaceVertex>>(non_manifold, nonmanifold_cluster_vertices, nv))
			{
				Scalar min_dist = sphere_mesh_constructor.min_distance_to_enveloppe(nv, sv);
				max_dist_shpae_to_enveloppe = std::max(max_dist_shpae_to_enveloppe, min_dist);
			}
			return true;
		});

		Scalar max_dist_enveloppe_to_shape = 0.0;
		cgogn::modeling::SkeletonSampler<Vec4, Vec3, Scalar> skeleton_sampler;

		parallel_foreach_cell(non_manifold, [&](NonManifoldVertex nv) {
			skeleton_sampler.add_vertex(value<Vec3>(non_manifold, non_manifold_vertex_position, nv),
										value<Scalar>(non_manifold, non_manifold_vertex_radius, nv));
			return true;
		});
		parallel_foreach_cell(non_manifold, [&](NonManifoldEdge ne) {
			auto& v_vec = incident_vertices(non_manifold, ne);
			skeleton_sampler.add_edge(value<Vec3>(non_manifold, non_manifold_vertex_position, v_vec[0]),
									  value<Scalar>(non_manifold, non_manifold_vertex_radius, v_vec[0]),
									  value<Vec3>(non_manifold, non_manifold_vertex_position, v_vec[1]),
									  value<Scalar>(non_manifold, non_manifold_vertex_radius, v_vec[1]));
			return true;
		});
		parallel_foreach_cell(non_manifold, [&](NonManifoldFace nf) {
			auto& v_vec = incident_vertices(non_manifold, nf);
			skeleton_sampler.add_triangle(value<Vec3>(non_manifold, non_manifold_vertex_position, v_vec[0]),
										  value<Scalar>(non_manifold, non_manifold_vertex_radius, v_vec[0]),
										  value<Vec3>(non_manifold, non_manifold_vertex_position, v_vec[1]),
										  value<Scalar>(non_manifold, non_manifold_vertex_radius, v_vec[1]),
										  value<Vec3>(non_manifold, non_manifold_vertex_position, v_vec[2]),
										  value<Scalar>(non_manifold, non_manifold_vertex_radius, v_vec[2]));
			return true;
		});

		Vec3 bbw = skeleton_sampler.BBwidth();
		Scalar dia_length = bbw.norm();
		float step = std::min(std::min(bbw.x(), bbw.y()), bbw.z()) / 200;
		skeleton_sampler.sample(step);
		std::vector<Vec3> enveloppe_points = skeleton_sampler.samples();

		// Build bvh
		auto surface_bvh = build_bvh(surface);
		std::pair<uint32, Vec3> bvh_res;
		for (Vec3 v : enveloppe_points)
		{
			surface_bvh->closest_point(v, &bvh_res);
			Vec3 closest_point_position = bvh_res.second;
			Scalar dist = (closest_point_position - v).norm();
			max_dist_enveloppe_to_shape = std::max(max_dist_enveloppe_to_shape, dist);
		}
		std::cout << "hausdorff distance shape to enveloppe:" << max_dist_shpae_to_enveloppe / dia_length * 100 << "%"
				  << std::endl;
		std::cout << "hausdorff distance enveloppe to shape:" << max_dist_enveloppe_to_shape / dia_length * 100 << "%"
				  << std::endl;
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
	void draw(View* view) override
	{
		auto& proj_matrix = view->projection_matrix();
		auto& view_matrix = view->modelview_matrix();
		if (draw_enveloppe)
			skeleton_drawer.draw(proj_matrix, view_matrix);
	}
	void left_panel() override
	{
		imgui_mesh_selector(surface_provider_, selected_surface_mesh_, "Surface", [&](SURFACE& m) {
			selected_surface_mesh_ = &m;
			surface_provider_->mesh_data(m).outlined_until_ = App::frame_time_ + 1.0;
		});

		ImGui::Separator();
		if (selected_surface_mesh_)
		{
			imgui_mesh_selector(nonmanifold_provider_, selected_medial_axis_, "Skeleton", [&](NONMANIFOLD& m) {
				selected_medial_axis_ = &m;
				nonmanifold_provider_->mesh_data(m).outlined_until_ = App::frame_time_ + 1.0;
			});
			if (selected_medial_axis_)
			{
				if (ImGui::Button("Compute hausdorff distance"))
				{
					compute_hausdorff_distance(*selected_medial_axis_, *selected_surface_mesh_);
				}
				if (ImGui::Button("Draw skeleton"))
				{
					interpolate_skeleton(*selected_medial_axis_);
				}
				ImGui::Checkbox("Show skeleton", &draw_enveloppe);
			}
		
			CoverageAxisParameter& c = coverage_axis_parameters_[selected_surface_mesh_];
			imgui_combo_attribute<SurfaceVertex, Vec3>(*selected_surface_mesh_, c.surface_vertex_position_, "Position",
													   [&](const std::shared_ptr<SurfaceAttribute<Vec3>>& attribute) {
														   c.surface_vertex_position_ = attribute;
													   });
			if (ImGui::Button("Coverage Axis"))
			{
				init_coverage_axis_plus_plus(*selected_surface_mesh_);
			}
			if (c.initialized_)
			{
				ImGui::RadioButton("Generate Random candidates", (int*)&c.candidate_generation_method, RANDOM);
				ImGui::RadioButton("Generate Shrinking ball candidates", (int*)&c.candidate_generation_method,
								   SHRINKING_BALL);
				ImGui::RadioButton("Generate Delaunay candidates", (int*)&c.candidate_generation_method, DELAUNAY);
				ImGui::DragInt("Candidated Number", (int*)&c.candidates_number, 100, 0, 50000);
				if (ImGui::Button("Generate candidates"))
				{
					generate_candidates(c);
				}
				if (c.candidates_valid)
				{
					ImGui::DragFloat("Dilation factor", &c.dilation_factor, 0.001f, 0.0f, 1.0f, "%.3f");
					ImGui::DragInt("Surface sample number", (int*)&c.surface_samples_number, 100, 0, 10000);
					if (ImGui::Button("Compute Coverage axis"))
					{
						coverage_axis(c);
					}
					ImGui::DragInt("Max Vertices", (int*)&c.max_selected_number, 1, 0, 500);
					if (ImGui::Button("Compute Coverage axis++"))
					{
						coverage_axis_plus_plus(c);
					}
					if (ImGui::Button("Generate connectivity by Power diagram"))
					{
						coverage_axis_connect_by_power_diagram(c);
					}
					if (ImGui::Button("Generate connectivity by Qmat"))
					{
						coverage_axis_connect_by_Qmat(c);
					}
				}
			}
		}
	}

private:
	cgogn::rendering::SkelShapeDrawer skeleton_drawer;
	SURFACE* selected_surface_mesh_ = nullptr;
	NONMANIFOLD* selected_medial_axis_ = nullptr;

	std::unordered_map<const SURFACE*, CoverageAxisParameter> coverage_axis_parameters_;
	PointVertex picked_sphere_;
	bool sphere_info_popup_ = false;
	bool draw_enveloppe = false;

	MeshProvider<POINT>* point_provider_;
	MeshProvider<SURFACE>* surface_provider_;
	MeshProvider<NONMANIFOLD>* nonmanifold_provider_;
	HighsSolution solution;
	Delaunay tri_;


};

} // namespace ui

} // namespace cgogn

#endif // CGOGN_MODULE_COVERAGE_AXIS_H_