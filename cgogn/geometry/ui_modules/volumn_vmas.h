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

#ifndef CGOGN_MODULE_VOLUMN_VMAS_H_
#define CGOGN_MODULE_VOLUMN_VMAS_H_

#include <cgogn/ui/app.h>
#include <cgogn/ui/module.h>

#include <cgogn/core/ui_modules/mesh_provider.h>

#include <cgogn/geometry/algos/area.h>
#include <cgogn/geometry/algos/fitting.h>
#include <cgogn/geometry/algos/length.h>
#include <cgogn/geometry/algos/medial_axis.h>
#include <cgogn/geometry/functions/angle.h>
#include <cgogn/geometry/functions/distance.h>
#include <cgogn/geometry/types/spherical_quadric.h>

#include <Eigen/Sparse>
#include <libacc/bvh_tree.h>
#include <libacc/bvh_tree_spheres.h>
#include <libacc/kd_tree.h>

#include <GLFW/glfw3.h>
// import CGAL
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Side_of_triangle_mesh.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/point_generators_3.h>

#include <boost/synapse/connect.hpp>
#include <random>
#include <set>

namespace cgogn
{

namespace ui
{

using geometry::Mat3;
using geometry::Mat4;
using geometry::Scalar;
using geometry::Spherical_Quadric;
using geometry::Vec3;
using geometry::Vec4;

template <typename SURFACE, typename POINTS, typename NONMANIFOLD>
class Volumn_VMAS : public ViewModule
{
	using K = CGAL::Exact_predicates_inexact_constructions_kernel;
	using Point_3 = K::Point_3;
	using CGAL_SurfaceMesh = CGAL::Surface_mesh<Point_3>;
	using Side_tester = CGAL::Side_of_triangle_mesh<CGAL_SurfaceMesh, K>;

	template <typename T>
	using SAttribute = typename mesh_traits<SURFACE>::template Attribute<T>;
	using SVertex = typename mesh_traits<SURFACE>::Vertex;
	using SEdge = typename mesh_traits<SURFACE>::Edge;
	using SFace = typename mesh_traits<SURFACE>::Face;

	template <typename T>
	using PAttribute = typename mesh_traits<POINTS>::template Attribute<T>;
	using PVertex = typename mesh_traits<POINTS>::Vertex;

	template <typename T>
	using NMAttribute = typename mesh_traits<NONMANIFOLD>::template Attribute<T>;
	using NMVertex = typename mesh_traits<NONMANIFOLD>::Vertex;
	using NMEdge = typename mesh_traits<NONMANIFOLD>::Edge;

	// enum UpdateMethod : uint32
	// {
	// 	FIT,
	// 	SQEM
	// };

	enum AutoSplitMode : uint32
	{
		MAX_NB_SPHERES,
		ERROR_THRESHOLD
	};
	enum DistanceMode : uint32
	{
		SPHERE_EUCLIDEAN_DISTANCE,
		SPHERE_CENTER_DISTANCE,
		SPHERE_POWER_DISTANCE,
		SPHERE_POWER_DISTANCE_SQUARED
	};

	enum CorrectionMode : uint32
	{
		CORRECT_ALWAYS,
		CORRECT_ON_SPLIT
	};

	const uint32 k = 10;

	struct SurfaceParameters
	{
		bool initialized_ = false;

		SURFACE* surface_;
		std::shared_ptr<SAttribute<Vec3>> surface_vertex_position_ = nullptr;
		std::shared_ptr<SAttribute<Vec3>> surface_vertex_position_original_ = nullptr;
		std::shared_ptr<SAttribute<Vec3>> surface_vertex_normal_ = nullptr;
		std::shared_ptr<SAttribute<Vec3>> surface_vertex_color_ = nullptr;
		std::shared_ptr<SAttribute<Vec3>> surface_face_normal_ = nullptr;
		std::shared_ptr<SAttribute<Scalar>> surface_vertex_area_pc_ = nullptr;
		std::shared_ptr<SAttribute<Scalar>> surface_vertex_area_surf_ = nullptr;
		std::shared_ptr<SAttribute<Scalar>> surface_face_area_ = nullptr;

		// std::shared_ptr<SAttribute<bool>> medial_axis_selected_ = nullptr;

		acc::BVHTree<uint32, Vec3>* surface_bvh_ = nullptr;
		std::vector<SFace> surface_bvh_faces_;
		acc::KDTree<3, uint32>* projected_surface_kdt_ = nullptr;
		std::vector<PVertex> projected_kdt_vertices_;
		acc::KDTree<3, uint32>* samples_kdt_ = nullptr;
		std::vector<PVertex> samples_kdt_vertices_;

		float32 noise_factor_ = 0.01f;

		bool point_cloud_mode_ = false;

		float32 sqem_update_lambda_ = 0.2f;
		float32 sqem_clustering_lambda_ = 0.2f; // initialized with mean edge length

		POINTS* spheres_;
		uint32 nb_spheres_ = 0;
		std::shared_ptr<PAttribute<Vec3>> spheres_position_ = nullptr;
		std::shared_ptr<PAttribute<Scalar>> spheres_radius_ = nullptr;
		std::shared_ptr<PAttribute<Vec4>> spheres_color_ = nullptr;
		std::shared_ptr<PAttribute<std::vector<PVertex>>> spheres_cluster_ = nullptr;
		std::shared_ptr<PAttribute<Scalar>> spheres_cluster_area_ = nullptr;
		std::shared_ptr<PAttribute<Vec3>> spheres_cluster_color_ = nullptr;
		std::shared_ptr<PAttribute<std::set<PVertex>>> spheres_neighbor_clusters_ = nullptr;
		std::shared_ptr<PAttribute<bool>> spheres_do_not_split_ = nullptr;
		std::shared_ptr<PAttribute<Scalar>> spheres_error_ = nullptr;
		std::shared_ptr<PAttribute<Scalar>> spheres_error_not_normalized_ = nullptr;
		// PAttribute<Scalar>* selected_spheres_error_ = nullptr;

		//Volumn Samples
		POINTS* samples_;
		std::shared_ptr<PAttribute<Vec3>> samples_position_ = nullptr;
		std::shared_ptr<PAttribute<Vec4>> samples_vertex_color_ = nullptr;
		std::shared_ptr<PAttribute<Vec3>> projected_samples_position_ = nullptr;
		std::shared_ptr<PAttribute<Vec3>> projected_samples_normal_ = nullptr;
		std::shared_ptr<PAttribute<Spherical_Quadric>> samples_quadric_ = nullptr;
		std::shared_ptr<PAttribute<Vec3>> medial_axis_position_ = nullptr;
		std::shared_ptr<PAttribute<Scalar>> medial_axis_radius_ = nullptr;
		std::shared_ptr<PAttribute<PVertex>> medial_axis_secondary_vertex_ = nullptr;
		std::shared_ptr<PAttribute<PVertex>> samples_vertex_sphere_ = nullptr;
		std::shared_ptr<PAttribute<Scalar>> samples_vertex_error_ = nullptr;
		std::shared_ptr<PAttribute<std::vector<PVertex>>> samples_vertex_knn_ = nullptr;
		uint32 sampels_numbers_ = 100000;
		std::unique_ptr<Side_tester> inside_tester_;
		CGAL_SurfaceMesh cgal_surface_mesh_;

		NONMANIFOLD* skeleton_;
		std::shared_ptr<NMAttribute<Vec3>> skeleton_position_ = nullptr;

		float32 filter_radius_threshold_ = 0.0f;
		float32 filter_angle_threshold_ = 0.0f;

		float32 init_dilation_factor_ = 3.0f;

		// UpdateMethod update_method_ = SQEM;

		// bool auto_split_outside_spheres_ = false;

		bool sphere_correction_ = true;
		CorrectionMode sphere_correction_mode_ = CORRECT_ALWAYS;

		bool auto_stop_ = false;
		bool auto_split_ = true;
		AutoSplitMode auto_split_mode_ = ERROR_THRESHOLD;
		DistanceMode distance_mode_ = SPHERE_EUCLIDEAN_DISTANCE;

		// bool auto_simplify_ = false;
		float32 auto_split_error_threshold_ = 0.00025f;
		uint32 auto_split_max_nb_spheres_ = 50;

		Scalar total_error_ = 0.0;
		Scalar total_error_not_normalized_ = 0.0;
		Scalar last_total_error_ = 0.0;
		Scalar total_error_diff_ = 0.0;
		Scalar min_error_ = 0.0;
		Scalar max_error_ = 0.0;
		PVertex max_error_sphere_ = PVertex();

		bool error_as_spheres_color_ = false;
		float32 spheres_transparency_ = 0.5f;

		//fuzzy membership
		std::shared_ptr<PAttribute<std::unordered_map<uint32, Scalar>>> samples_membership_ = nullptr;
		float tau_ = 0.01f;
		uint32 top_k_spheres_ = 8;
		Scalar eps_ = 1e-12;

		float theta_gap_log_ = 0.1f;

		uint32 iteration_count_ = 0;
		std::mutex mutex_;
		bool running_ = false;
		bool stopping_ = false;
		bool slow_down_ = true;
		uint32 update_rate_ = 20;
		// poisson disk sampling
		Scalar r_ = 0.015;
		uint32 K_ = 30;
	};

public:
	Volumn_VMAS(const App& app)
		: ViewModule(app, "Volumn_VMAS (" + std::string{mesh_traits<SURFACE>::name} + "," +
							  std::string{mesh_traits<POINTS>::name} + ")")
	{
	}
	~Volumn_VMAS()
	{
	}

	void normalize_surface_mesh(CGAL_SurfaceMesh& mesh)
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
			Point_3 p = mesh.point(v);
			double nx = (p.x() - bbox.xmin()) / max_dim;
			double ny = (p.y() - bbox.ymin()) / max_dim;
			double nz = (p.z() - bbox.zmin()) / max_dim;
			mesh.point(v) = Point_3(nx, ny, nz);
		}
	}

	bool is_inside(const SurfaceParameters& p, const Vec3& pos)
	{
		Point_3 query(pos.x(), pos.y(), pos.z());
		CGAL::Bounded_side res = (*p.inside_tester_)(query);
		return res == CGAL::ON_BOUNDED_SIDE;
	}
	Vec3 random_sample_around(const Vec3& p, Scalar R)
	{
		static thread_local std::mt19937 rng{std::random_device{}()};
		static thread_local std::uniform_real_distribution<Scalar> U(0.0, 1.0);

		Scalar u = U(rng);
		Scalar v = U(rng);
		Scalar w = U(rng);

		Scalar R3 = R * R * R;
		Scalar r = std::cbrt(R3 + u * (8 * R3 - R3)); // r = pow((R^3 + u(8R^3 - R^3)), 1/3)

		Scalar phi = v * 2.0 * M_PI;
		Scalar theta = std::acos(1.0 - 2.0 * w);

		Scalar x = r * std::sin(theta) * std::cos(phi);
		Scalar y = r * std::sin(theta) * std::sin(phi);
		Scalar z = r * std::cos(theta);

		return p + Vec3(x, y, z);
	}

	std::tuple<uint32, uint32, uint32> pos_to_grid_cell(const Vec3& pos, Scalar cell_size, uint32 grid_size)
	{
		uint32 x = std::min(uint32(pos.x() / cell_size), grid_size - 1);
		uint32 y = std::min(uint32(pos.y() / cell_size), grid_size - 1);
		uint32 z = std::min(uint32(pos.z() / cell_size), grid_size - 1);
		return {x, y, z};
	}
	void poisson_disk_sampling(SurfaceParameters& p)
	{
		using SampleIndex = std::size_t;
		const SampleIndex INVALID = (std::numeric_limits<SampleIndex>::max)();
		Scalar cell_size = p.r_ / std::sqrt(3.0);
		uint32 grid_size = uint32(std::ceil(1.0 / cell_size));
		std::vector<std::vector<std::vector<SampleIndex>>> grid(
			grid_size, std::vector<std::vector<SampleIndex>>(grid_size, std::vector<SampleIndex>(grid_size, INVALID)));
		std::vector<SampleIndex> active_list;
		std::vector<Vec3> samples;
		std::mt19937 rng{std::random_device{}()};
		auto pick_active_index = [&](void) -> size_t {
			std::uniform_int_distribution<size_t> dist(0, active_list.size() - 1);
			return dist(rng);
		};
		Vec3 bias = Vec3(0.5, 0.5, 0.5);
		CGAL::Random_points_in_cube_3<Point_3> generator(0.5); // sample first point
		while (active_list.empty())
		{
			Point_3 s = *generator++;
			Vec3 pos = Vec3(s.x(), s.y(), s.z()) + bias;
			if (!is_inside(p, pos))
				continue;
			auto [x, y, z] = pos_to_grid_cell(pos, cell_size, grid_size);
			SampleIndex new_index = (SampleIndex)samples.size();
			samples.push_back(pos);
			grid[x][y][z] = new_index;
			active_list.push_back(new_index);
		}
		while (active_list.size() != 0)
		{ // Randomly pick a point from the active list
			size_t idx = pick_active_index();
			SampleIndex current_index = active_list[idx];
			Vec3 current_sample = samples[current_index];
			bool permenantly_remove = true;
			for (uint32 i = 0; i < p.K_; i++)
			{
				Vec3 next_pos;
				next_pos = random_sample_around(current_sample, p.r_);
				if (!is_inside(p, next_pos))
					continue;
				auto [x, y, z] = pos_to_grid_cell(next_pos, cell_size,
												  grid_size); // Test if next_pos is at least R far from other samples
				bool valid = true;
				if (grid[x][y][z] != INVALID)
				{
					valid = false;
				}
				for (int32 dx = -2; dx <= 2 && valid; dx++)
				{
					for (int32 dy = -2; dy <= 2 && valid; dy++)
					{
						for (int32 dz = -2; dz <= 2 && valid; dz++)
						{
							int32 nx = int32(x) + dx;
							int32 ny = int32(y) + dy;
							int32 nz = int32(z) + dz;
							if (nx >= 0 && nx < int32(grid_size) && ny >= 0 && ny < int32(grid_size) && nz >= 0 &&
								nz < int32(grid_size))
							{
								SampleIndex si = grid[nx][ny][nz];
								if (si != INVALID)
								{
									const Vec3& neighbor_pos = samples[si];
									Vec3 cp = next_pos - neighbor_pos;
									if (cp.dot(cp) < p.r_ * p.r_)
									{
										valid = false;
									}
								}
							}
						}
					}
				}
				if (valid)
				{
					SampleIndex new_index = (SampleIndex)samples.size();
					samples.push_back(next_pos);
					active_list.push_back(new_index);
					grid[x][y][z] = new_index;
					permenantly_remove = false;
					break;
				}
			}
			if (permenantly_remove)
			{
				active_list[idx] = active_list.back();
				active_list.pop_back();
			}
		}
		for (const Vec3& pos : samples)
		{
			PVertex new_sample = add_vertex(*p.samples_);
			std::pair<uint32, Vec3> bvh_res;
			p.surface_bvh_->closest_point(pos, &bvh_res);
			Vec3 closest_surface_position = bvh_res.second;
			SFace closest_face = p.surface_bvh_faces_[bvh_res.first];
			Vec3 closest_face_normal = value<Vec3>(*p.surface_, p.surface_face_normal_, closest_face);
			value<Vec3>(*p.samples_, p.samples_position_, new_sample) = pos;
			value<Vec3>(*p.samples_, p.projected_samples_position_, new_sample) = closest_surface_position;
			value<Vec3>(*p.samples_, p.projected_samples_normal_, new_sample) = closest_face_normal;
		}
		points_provider_->emit_connectivity_changed(*p.samples_);
	}

	void set_selected_surface(SURFACE& s)
	{
		selected_surface_ = &s;
		picked_sphere_ = PVertex();
	}

	void set_surface_vertex_position(SURFACE& s, const std::shared_ptr<SAttribute<Vec3>>& surface_vertex_position)
	{
		SurfaceParameters& p = surface_parameters_[&s];
		p.surface_vertex_position_ = surface_vertex_position;
	}

	void add_surface_noise(SurfaceParameters& p)
	{
		if (!p.surface_vertex_position_ || !p.surface_vertex_normal_)
		{
			std::cout << "No surface vertex position or normal attribute set" << std::endl;
			return;
		}

		const MeshData<SURFACE>& md = surface_provider_->mesh_data(*p.surface_);
		Scalar bb_diag = (md.bb_max_ - md.bb_min_).norm() / 2.0;

		parallel_foreach_cell(*p.surface_, [&](SVertex v) -> bool {
			uint32 v_index = index_of(*p.surface_, v);
			float32 r = float32(rand()) / float32(RAND_MAX);
			int i = rand() % 2;
			if (i == 0)
				r *= -1.0f;
			(*p.surface_vertex_position_)[v_index] +=
				r * p.noise_factor_ * bb_diag * (*p.surface_vertex_normal_)[v_index];
			return true;
		});

		surface_provider_->emit_attribute_changed(*p.surface_, p.surface_vertex_position_.get());
	}

	void restore_surface_position(SurfaceParameters& p)
	{
		if (!p.surface_vertex_position_ || !p.surface_vertex_position_original_)
		{
			std::cout << "No surface vertex position or original position attribute set" << std::endl;
			return;
		}

		p.surface_vertex_position_->copy(p.surface_vertex_position_original_.get());

		surface_provider_->emit_attribute_changed(*p.surface_, p.surface_vertex_position_.get());
	}

	void compute_quadrics(SurfaceParameters& p)
	{
		parallel_foreach_cell(*p.samples_, [&](PVertex v) -> bool {
			uint32 v_index = index_of(*p.samples_, v);
			Spherical_Quadric& q = (*p.samples_quadric_)[v_index];
			q.clear();
			const Vec3& pos = (*p.projected_samples_position_)[v_index];
			const Vec3& n = (*p.projected_samples_normal_)[v_index];
			q = Spherical_Quadric(Vec4(pos.x(), pos.y(), pos.z(), 0), Vec4(n.x(), n.y(), n.z(), 1));
			return true;
		});
	}

	void init_surface_data(SURFACE& s)
	{
		SurfaceParameters& p = surface_parameters_[&s];
		p.surface_ = &s;

		if (!p.surface_vertex_position_)
		{
			std::cout << "No surface vertex position attribute set" << std::endl;
			return;
		}

		// set signal connections to update the data when the surface connectivity or position changes

		if (surface_connections_.find(&s) == surface_connections_.end())
		{
			surface_connections_[&s].push_back(
				boost::synapse::connect<typename MeshProvider<SURFACE>::connectivity_changed>(&s, [this, &s = s]() {
					SurfaceParameters& p = surface_parameters_[&s];
					init_surface_data(s);
				}));
			surface_connections_[&s].push_back(
				boost::synapse::connect<typename MeshProvider<SURFACE>::template attribute_changed_t<Vec3>>(
					&s, [this, &s = s](SAttribute<Vec3>* attribute) {
						SurfaceParameters& p = surface_parameters_[&s];
						if (attribute == p.surface_vertex_position_.get())
							init_surface_data(s);
					}));
		}

		// save original positions

		p.surface_vertex_position_original_ = get_attribute<Vec3, SVertex>(s, "position_original");
		if (!p.surface_vertex_position_original_)
		{
			p.surface_vertex_position_original_ = add_attribute<Vec3, SVertex>(s, "position_original");
			p.surface_vertex_position_original_->copy(p.surface_vertex_position_.get());
		}

		// compute normals
		p.surface_face_normal_ = get_or_add_attribute<Vec3, SFace>(s, "normal");
		geometry::compute_normal<SFace>(s, p.surface_vertex_position_.get(), p.surface_face_normal_.get());

		p.surface_vertex_normal_ = get_or_add_attribute<Vec3, SVertex>(s, "normal");
		geometry::compute_normal<SVertex>(s, p.surface_vertex_position_.get(), p.surface_vertex_normal_.get());

		// vertex colors

		p.surface_vertex_color_ = get_or_add_attribute<Vec3, SVertex>(s, "color");

		// compute areas

		p.surface_face_area_ = get_or_add_attribute<Scalar, SFace>(s, "area");
		geometry::compute_area<SFace>(s, p.surface_vertex_position_.get(), p.surface_face_area_.get());

		/* p.surface_vertex_area_pc_ = get_or_add_attribute<Scalar, SVertex>(s, "area");
		foreach_cell(s, [&](SVertex v) -> bool {
			uint32 v_index = index_of(s, v);
			Scalar sum = 0.0;
			const Vec3& pos = (*p.surface_vertex_position_)[v_index];
			for (SVertex u : (*p.surface_vertex_knn_)[v_index])
				sum += (value<Vec3>(s, p.surface_vertex_position_, u) - pos).norm();
			(*p.surface_vertex_area_pc_)[v_index] = (sum * sum) / (2.0 * k);
			return true;
		});*/

		p.surface_vertex_area_surf_ = get_or_add_attribute<Scalar, SVertex>(s, "area");
		geometry::compute_area<SVertex>(s, p.surface_vertex_position_.get(), p.surface_vertex_area_surf_.get(),
										geometry::VertexAreaPolicy::THIRD);

		MeshData<SURFACE>& md = surface_provider_->mesh_data(s);
		uint32 nb_vertices = md.template nb_cells<SVertex>();
		uint32 nb_faces = md.template nb_cells<SFace>();

		// create BVH for the surface vertices
		auto bvh_vertex_index = get_or_add_attribute<uint32, SVertex>(s, "__bvh_vertex_index");

		std::vector<Vec3> vertex_position_vector;
		vertex_position_vector.reserve(nb_vertices);
		uint32 idx = 0;
		foreach_cell(s, [&](SVertex v) -> bool {
			value<uint32>(s, bvh_vertex_index, v) = idx++;
			vertex_position_vector.push_back(value<Vec3>(s, p.surface_vertex_position_, v));
			return true;
		});

		p.surface_bvh_faces_.clear();
		p.surface_bvh_faces_.reserve(nb_faces);
		std::vector<uint32> face_vertex_indices;
		face_vertex_indices.reserve(nb_faces * 3);
		foreach_cell(s, [&](SFace f) -> bool {
			p.surface_bvh_faces_.push_back(f);
			foreach_incident_vertex(s, f, [&](SVertex v) -> bool {
				face_vertex_indices.push_back(value<uint32>(s, bvh_vertex_index, v));
				return true;
			});
			return true;
		});
		if (p.surface_bvh_)
			delete p.surface_bvh_;
		p.surface_bvh_ = new acc::BVHTree<uint32, Vec3>(face_vertex_indices, vertex_position_vector);
		remove_attribute<SVertex>(s, bvh_vertex_index);

		// load CGAL surface mesh
		std::string filename = surface_provider_->mesh_filename(s);
		if (!filename.empty())
		{
			if (!CGAL::IO::read_polygon_mesh(filename, p.cgal_surface_mesh_) || p.cgal_surface_mesh_.is_empty())
			{
				std::cout << "Error loading CGAL surface mesh from file: " << filename << std::endl;
				return;
			}
		}
		normalize_surface_mesh(p.cgal_surface_mesh_);
		// CGAL side tester
		p.inside_tester_ = std::make_unique<Side_tester>(p.cgal_surface_mesh_);

		if (!p.samples_)
			p.samples_ = points_provider_->add_mesh(surface_provider_->mesh_name(s) + "_samples");
		// Volumn Samples
		p.samples_position_ = get_or_add_attribute<Vec3, PVertex>(*p.samples_, "position");
		p.samples_vertex_color_ = get_or_add_attribute<Vec4, PVertex>(*p.samples_, "color");
		p.projected_samples_position_ = get_or_add_attribute<Vec3, PVertex>(*p.samples_, "projected_position");
		p.projected_samples_normal_ = get_or_add_attribute<Vec3, PVertex>(*p.samples_, "projected_normal");
		p.samples_membership_ =
			get_or_add_attribute<std::unordered_map<uint32, Scalar>, PVertex>(*p.samples_, "membership");
		// initialize volumn samples
		poisson_disk_sampling(p);

		// create KDTree for the projected surface vertices
		p.projected_kdt_vertices_.clear();
		p.samples_kdt_vertices_.clear();
		p.projected_kdt_vertices_.reserve(p.sampels_numbers_);
		p.samples_kdt_vertices_.reserve(p.sampels_numbers_);
		std::vector<Vec3> projected_surface_position_vector;
		std::vector<Vec3> samples_position_vector;
		foreach_cell(*p.samples_, [&](PVertex v) {
			p.projected_kdt_vertices_.push_back(v);
			p.samples_kdt_vertices_.push_back(v);
			samples_position_vector.push_back(value<Vec3>(*p.samples_, p.samples_position_, v));
			projected_surface_position_vector.push_back(value<Vec3>(*p.samples_, p.projected_samples_position_, v));
			return true;
		});
		if (p.projected_surface_kdt_)
			delete p.projected_surface_kdt_;
		if (p.samples_kdt_)
			delete p.samples_kdt_;
		p.samples_kdt_ = new acc::KDTree<3, uint32>(samples_position_vector);
		p.projected_surface_kdt_ = new acc::KDTree<3, uint32>(projected_surface_position_vector);

		// compute knn graph
		p.samples_vertex_knn_ = get_or_add_attribute<std::vector<PVertex>, PVertex>(*p.samples_, "knn");
		foreach_cell(*p.samples_, [&](PVertex v) -> bool {
			uint32 v_index = index_of(*p.samples_, v);
			const Vec3& pos = (*p.samples_position_)[v_index];
			std::vector<std::pair<uint32, Scalar>> k_res;
			p.samples_kdt_->find_nns(pos, k, &k_res);
			(*p.samples_vertex_knn_)[v_index].clear();
			(*p.samples_vertex_knn_)[v_index].reserve(k_res.size());
			for (const auto& [idx, dist] : k_res)
			{
				PVertex nv = p.samples_kdt_vertices_[idx];
				if (nv != v)
					(*p.samples_vertex_knn_)[v_index].push_back(nv);
			}
			return true;
		});

		// initialize SQEM quadrics
		p.samples_quadric_ = get_or_add_attribute<Spherical_Quadric, PVertex>(*p.samples_, "quadric");
		compute_quadrics(p);

		// compute shrinking balls for the surface vertices

		p.medial_axis_position_ = get_or_add_attribute<Vec3, PVertex>(*p.samples_, "medial_axis_position");
		p.medial_axis_radius_ = get_or_add_attribute<Scalar, PVertex>(*p.samples_, "medial_axis_radius");
		p.medial_axis_secondary_vertex_ =
			get_or_add_attribute<PVertex, PVertex>(*p.samples_, "medial_axis_secondary_vertex_");

		parallel_foreach_cell(*p.samples_, [&](PVertex v) -> bool {
			uint32 v_index = index_of(*p.samples_, v);
			auto [c, r, q] = geometry::shrinking_ball_center(
				*p.samples_, (*p.projected_samples_position_)[v_index], (*p.projected_samples_normal_)[v_index],
				p.projected_samples_position_.get(), p.projected_surface_kdt_, p.projected_kdt_vertices_);
			(*p.medial_axis_position_)[v_index] = c;
			(*p.medial_axis_radius_)[v_index] = r;
			(*p.medial_axis_secondary_vertex_)[v_index] = q;
			return true;
		});

		/* parallel_foreach_cell(*p.samples_, [&](PVertex v) -> bool {
			uint32 v_index = index_of(*p.samples_, v);
			auto [c, r] = geometry::shrinking_ball_center(
				*p.surface_,
				(*p.projected_samples_position_)[v_index],
				(*p.projected_samples_normal_)[v_index],
				p.surface_bvh_,
				p.surface_bvh_faces_);
			(*p.medial_axis_position_)[v_index] = c;
			(*p.medial_axis_radius_)[v_index] = r;
			return true;
		});*/

		// create the spheres mesh

		if (!p.spheres_)
			p.spheres_ = points_provider_->add_mesh(surface_provider_->mesh_name(s) + "_spheres");

		p.spheres_position_ = get_or_add_attribute<Vec3, PVertex>(*p.spheres_, "position");
		p.spheres_radius_ = get_or_add_attribute<Scalar, PVertex>(*p.spheres_, "radius");
		p.spheres_color_ = get_or_add_attribute<Vec4, PVertex>(*p.spheres_, "color");

		p.spheres_cluster_ = get_or_add_attribute<std::vector<PVertex>, PVertex>(
			*p.spheres_, "cluster"); // surface vertices in the cluster
		p.spheres_cluster_color_ = get_or_add_attribute<Vec3, PVertex>(*p.spheres_, "cluster_color");
		p.spheres_cluster_area_ = get_or_add_attribute<Scalar, PVertex>(*p.spheres_, "cluster_area");

		p.spheres_error_ = get_or_add_attribute<Scalar, PVertex>(*p.spheres_, "error");
		p.spheres_error_not_normalized_ = get_or_add_attribute<Scalar, PVertex>(*p.spheres_, "error_not_normalized");
		p.spheres_error_->fill(0.0);
		p.spheres_error_not_normalized_->fill(0.0);

		p.spheres_do_not_split_ = get_or_add_attribute<bool, PVertex>(*p.spheres_, "do_not_split");

		p.samples_vertex_error_ =
			get_or_add_attribute<Scalar, PVertex>(*p.samples_, "error"); // error of a vertex w.r.t. its sphere

		p.samples_vertex_sphere_ =
			get_or_add_attribute<PVertex, PVertex>(*p.samples_, "sphere"); // cluster of the surface vertex
		p.spheres_neighbor_clusters_ =
			get_or_add_attribute<std::set<PVertex>, PVertex>(*p.spheres_, "neighbor_clusters"); // neighbor clusters

		// create the skeleton mesh

		if (!p.skeleton_)
			p.skeleton_ = non_manifold_provider_->add_mesh(surface_provider_->mesh_name(s) + "_skeleton");
		p.skeleton_position_ = get_or_add_attribute<Vec3, NMVertex>(*p.skeleton_, "position");

		// if we already have spheres, we need to recompute the clusters and errors
		// (cleans out surface vertex sphere data)

		//compute_clusters(p);
		compute_membership_soft(p);
		materialize_top1_labels(p);

		// update the render data (spheres and skeleton)
		// (the skeleton is reconstructed)

		if (!p.running_)
			update_render_data(p);

		p.initialized_ = true;
	}

	void init_spheres(SURFACE& s, uint32 max_nb_spheres)
	{
		SurfaceParameters& p = surface_parameters_[&s];

		clear(*p.spheres_);

		std::vector<PVertex> sorted_vertices;
		foreach_cell(*p.samples_, [&](PVertex v) -> bool {
			sorted_vertices.push_back(v);
			return true;
		});
		std::sort(sorted_vertices.begin(), sorted_vertices.end(), [&](PVertex a, PVertex b) {
			// sort candidate spheres by decreasing radius
			return value<Scalar>(*p.samples_, p.medial_axis_radius_, a) >
				   value<Scalar>(*p.samples_, p.medial_axis_radius_, b);
			// sort candidate spheres by decreasing number of covered vertices
			// return value<uint32>(*p.surface_, nb_covered, a) > value<uint32>(*p.surface_, nb_covered, b);
		});

		p.nb_spheres_ = 0;

		PVertex v = sorted_vertices[0];

		uint32 v_index = index_of(*p.samples_, v);

		const Vec3& vp = (*p.medial_axis_position_)[v_index];
		Scalar vr = (*p.medial_axis_radius_)[v_index];

		PVertex sphere = add_vertex(*p.spheres_);
		p.nb_spheres_++;
		uint32 sphere_index = index_of(*p.spheres_, sphere);

		(*p.spheres_position_)[sphere_index] = vp;
		(*p.spheres_radius_)[sphere_index] = vr;
		(*p.spheres_cluster_color_)[sphere_index] = Vec3(
			0.5 + 0.5 * (rand() % 256) / 256.0, 0.5 + 0.5 * (rand() % 256) / 256.0, 0.5 + 0.5 * (rand() % 256) / 256.0);
		// Vec3(rand() % 256 / 255.0f, rand() % 256 / 255.0f, rand() % 256 / 255.0f);

		//compute_clusters(p);
		compute_membership_soft(p);
		materialize_top1_labels(p);

		if (!p.running_)
			update_render_data(p);
	}

	// void compute_clusters(SurfaceParameters& p)
	// {
	// 	// clean cluster affectation
	// 	parallel_foreach_cell(*p.spheres_, [&](PVertex v) -> bool {
	// 		uint32 v_index = index_of(*p.spheres_, v);
	// 		(*p.spheres_cluster_)[v_index].clear();
	// 		(*p.spheres_cluster_area_)[v_index] = 0.0;
	// 		return true;
	// 	});
	// 	p.samples_vertex_sphere_->fill(PVertex());

	// 	MeshData<POINTS>& md = points_provider_->mesh_data(*p.spheres_);
	// 	if (p.nb_spheres_ == 0)
	// 		return;

	// 	// auto start = std::chrono::high_resolution_clock::now();

	// 	parallel_foreach_cell(*p.samples_, [&](PVertex v) -> bool {
	// 		uint32 v_index = index_of(*p.samples_, v);
	// 		// if (!(*p.medial_axis_selected_)[v_index])
	// 		// 	return true;

	// 		// Each volumn sample has weight 1
	// 		Scalar a = 1.0;

	// 		const Vec3& vp = (*p.projected_samples_position_)[v_index];
	// 		Scalar min_distance = std::numeric_limits<Scalar>::max();
	// 		PVertex closest_sphere;
	// 		uint32 closest_sphere_index;

	// 		foreach_cell(*p.spheres_, [&](PVertex pv) {
	// 			uint32 pv_index = index_of(*p.spheres_, pv);
	// 			const Vec3& center = (*p.spheres_position_)[pv_index];
	// 			Scalar radius = (*p.spheres_radius_)[pv_index];
	// 			Scalar dist_sqem =
	// 				(*p.samples_quadric_)[v_index].eval(Vec4(center.x(), center.y(), center.z(), radius));
	// 			Scalar dist_other = 0.0;
	// 			switch (p.distance_mode_)
	// 			{
	// 			case SPHERE_EUCLIDEAN_DISTANCE: {
	// 				dist_other = ((*p.projected_samples_position_)[v_index] - center).norm() - radius;
	// 				dist_other *= dist_other;
	// 			}
	// 			break;
	// 			case SPHERE_CENTER_DISTANCE: {
	// 				dist_other = ((*p.samples_position_)[v_index] - center).norm();
	// 				dist_other *= dist_other;
	// 			}
	// 			break;
	// 			case SPHERE_POWER_DISTANCE: {
	// 				dist_other =
	// 					((*p.samples_position_)[v_index] - center).dot((*p.samples_position_)[v_index] - center) -
	// 					radius * radius;
	// 			}
	// 			break;
	// 			case SPHERE_POWER_DISTANCE_SQUARED: {
	// 				dist_other =
	// 					((*p.samples_position_)[v_index] - center).dot((*p.samples_position_)[v_index] - center) -
	// 					radius * radius;
	// 				dist_other *= dist_other;
	// 			}
	// 			break;
	// 			default:
	// 				break;
	// 			}
	// 			dist_other *= a;
	// 			Scalar dist = dist_sqem + p.sqem_clustering_lambda_ * dist_other;
	// 			if (dist < min_distance)
	// 			{
	// 				min_distance = dist;
	// 				closest_sphere = pv;
	// 				closest_sphere_index = pv_index;
	// 			}
	// 			return true;
	// 		});

	// 		value<PVertex>(*p.samples_, p.samples_vertex_sphere_, v) = closest_sphere;

	// 		std::lock_guard<std::mutex> lock(spheres_mutex_[closest_sphere_index % spheres_mutex_.size()]);
	// 		value<std::vector<PVertex>>(*p.spheres_, p.spheres_cluster_, closest_sphere).push_back(v);
	// 		value<Scalar>(*p.spheres_, p.spheres_cluster_area_, closest_sphere) += a;

	// 		return true;
	// 	});
	// 	// 	break;
	// 	// }
	// 	// }

	// 	// auto end = std::chrono::high_resolution_clock::now();

	// 	// std::cout << "Cluster computation time: " << std::chrono::duration<Scalar>(end - start).count() << "s"
	// 	// 		  << std::endl;

	// 	//remove small clusters
	// 	foreach_cell(*p.spheres_, [&](PVertex v) -> bool {
	// 		std::vector<PVertex>& cluster = value<std::vector<PVertex>>(*p.spheres_, p.spheres_cluster_, v);
	// 		if (cluster.size() < 4)
	// 		{
	// 			for (PVertex sv : cluster)
	// 				value<PVertex>(*p.samples_, p.samples_vertex_sphere_, sv) = PVertex();
	// 			remove_vertex(*p.spheres_, v);
	// 			p.nb_spheres_--;
	// 		}
	// 		return true;
	// 	});
	// }

	void compute_membership_soft(SurfaceParameters& p)
	{
		parallel_foreach_cell(*p.spheres_, [&](PVertex v) -> bool {
			uint32 v_index = index_of(*p.spheres_, v);
			(*p.spheres_cluster_)[v_index].clear();
			(*p.spheres_cluster_area_)[v_index] = 0.0;
			return true;
		});
		p.samples_vertex_sphere_->fill(PVertex());

		MeshData<POINTS>& md = points_provider_->mesh_data(*p.spheres_);
		if (p.nb_spheres_ == 0)
			return;
		struct Cand
		{
			uint32 sidx;
			Scalar E;
		};
		foreach_cell(*p.samples_, [&](PVertex v) -> bool {
			uint32 v_index = index_of(*p.samples_, v);
			
			std::vector<Cand> candidates;
			candidates.reserve(p.nb_spheres_);
			// Each volumn sample has weight 1
			Scalar a = 1.0;

			const Vec3& vp = (*p.projected_samples_position_)[v_index];
			Scalar total_weight = 0.0;
			std::vector<std::pair<uint32, Scalar>> weights;
			foreach_cell(*p.spheres_, [&](PVertex pv) {
				uint32 pv_index = index_of(*p.spheres_, pv);
				const Vec3& center = (*p.spheres_position_)[pv_index];
				Scalar radius = (*p.spheres_radius_)[pv_index];
				Scalar dist_sqem =
					(*p.samples_quadric_)[v_index].eval(Vec4(center.x(), center.y(), center.z(), radius));
				Scalar dist_other = 0.0;
				switch (p.distance_mode_)
				{
				case SPHERE_EUCLIDEAN_DISTANCE: {
					dist_other = ((*p.projected_samples_position_)[v_index] - center).norm() - radius;
					dist_other *= dist_other;
				}
				break;
				case SPHERE_CENTER_DISTANCE: {
					dist_other = ((*p.samples_position_)[v_index] - center).norm();
					dist_other *= dist_other;
				}
				break;
				case SPHERE_POWER_DISTANCE: {
					dist_other =
						((*p.samples_position_)[v_index] - center).dot((*p.samples_position_)[v_index] - center) -
						radius * radius;
				}
				break;
				case SPHERE_POWER_DISTANCE_SQUARED: {
					dist_other =
						((*p.samples_position_)[v_index] - center).dot((*p.samples_position_)[v_index] - center) -
						radius * radius;
					dist_other *= dist_other;
				}
				break;
				default:
					break;
				}
				dist_other *= a;
				Scalar dist = dist_sqem + p.sqem_clustering_lambda_ * dist_other;
				candidates.push_back({ pv_index, dist });
				return true;
			});
			const uint32 K = std::min(p.top_k_spheres_, uint32(candidates.size()));
			std::nth_element(candidates.begin(), candidates.begin() + K, candidates.end(), [](const Cand& a, const Cand& b) {
				return a.E < b.E;
			});
			candidates.resize(K);
			std::sort(candidates.begin(), candidates.end(), [](const Cand& a, const Cand& b) { return a.E < b.E; });
			
			
			
			Scalar Emin = candidates[0].E;
			Scalar tau = p.tau_;
			if (candidates.size() == 0)
				return true;
			std::vector<Scalar> log_e;
			log_e.reserve(candidates.size());
			for (const auto& c : candidates)
				log_e.push_back(std::log(c.E + 1e-16));
			
		
			
			const Scalar eps = std::max<Scalar>(p.eps_, 1e-16);
			std::vector<Scalar> gap;
			gap.reserve((K > 1) ? (K - 1): 0);

			for (uint32 t = 0; t+1<K; ++t)
			{ 
				Scalar x_t = std::log(candidates[t].E + eps);
				Scalar x_t1 = std::log(candidates[t + 1].E + eps);
				gap.push_back(x_t1 - x_t);
			
			}
			const Scalar theta_gap = p.theta_gap_log_;
			uint32 pike =K - 1;
			for (uint32 t = 0; t < gap.size(); ++t)
			{
				if (gap[t] >= p.theta_gap_log_)
				{
					pike = t;
					break;
				}
			}
			Scalar denom = 0.0;
			std::vector<Cand> keep;
			keep.reserve(pike+1);
			for (uint32 t = 0; t <= pike; ++t)
				keep.push_back(candidates[t]);
			if (keep.empty())
				keep.push_back(candidates.front());
			/*if (keep.size() > 1)
			{
				std::cout << "candidates: ";
				for (const auto& c : candidates)
					std::cout << c.E << " ";
				std::cout << std::endl;
				std::cout << "log_e: ";
				for (const auto& le : log_e)
					std::cout << le << " ";
				std::cout << std::endl;
				std::cout << " keep: ";
				for (auto& c : keep)
					std::cout << c.E << " ";
				std::cout << std::endl;
				std::cout << " ----------------------------------------" << std::endl;
			}*/
			const auto nK = keep.size();
			std::vector<Scalar> ww(nK);
			for (uint32 t = 0; t < nK; ++t)
			{
				Scalar w = std::exp(-keep[t].E/ tau);
				denom += w;
				ww[t] = w;
			}
			auto& mmap = (*p.samples_membership_)[v_index];
			mmap.clear();
			mmap.reserve(nK);
			const Scalar inv = (denom > 0) ? (1.0/denom) : 1.0;
			for (uint32 t = 0; t < nK; ++t)
			{
				Scalar w = ww[t] * inv;
				if (w > 0) mmap[keep[t].sidx] = w;
			}
			for (const auto& [sidx, w] : mmap)
			{
				std::lock_guard<std::mutex> lock(spheres_mutex_[sidx % spheres_mutex_.size()]);
				(*p.spheres_cluster_)[sidx].push_back(v);
				(*p.spheres_cluster_area_)[sidx] += a * w;
			}
			return true;
		});
	}

	void materialize_top1_labels(SurfaceParameters& p)
	{
		if (p.nb_spheres_ == 0)
			return;
		foreach_cell(*p.samples_, [&](PVertex vi)->bool {
			const uint32 vi_idx = index_of(*p.samples_, vi);
			const auto& mmap = (*p.samples_membership_)[vi_idx];
			assert(!mmap.empty());
			uint32 best_idx = 0; 
			Scalar bw = -1.0;
			for (const auto& kv : mmap) {
				if (kv.second > bw)
				{
					bw = kv.second;
					best_idx = kv.first;
				}
			}
			value<PVertex>(*p.samples_, p.samples_vertex_sphere_, vi) = of_index<PVertex>(*p.spheres_, best_idx);
			(*p.spheres_cluster_)[best_idx].push_back(vi);
			return true;
		});
	}

	void compute_spheres_error(SurfaceParameters& p)
	{
		parallel_foreach_cell(*p.spheres_, [&](PVertex v) {
			uint32 v_index = index_of(*p.spheres_, v);

			const Vec3& center = (*p.spheres_position_)[v_index];
			Scalar radius = (*p.spheres_radius_)[v_index];
			const std::vector<PVertex>& cluster = (*p.spheres_cluster_)[v_index];

			Scalar cluster_error = 0.0;
			foreach_cell(*p.samples_, [&](PVertex sv) -> bool {
				uint32 sv_index = index_of(*p.samples_, sv);
				auto& mmap = (*p.samples_membership_)[sv_index];
				auto it = mmap.find(v_index);
				if (it == mmap.end())
					return true;
				Scalar membership = it->second;

				Scalar a = 1.0;
				Scalar dist_sqem =
					(*p.samples_quadric_)[sv_index].eval(Vec4(center.x(), center.y(), center.z(), radius));
				Scalar dist_other = 0.0;
				switch (p.distance_mode_)
				{
				case SPHERE_EUCLIDEAN_DISTANCE: {
					dist_other = ((*p.projected_samples_position_)[sv_index] - center).norm() - radius;
					dist_other *= dist_other;
				}
				break;
				case SPHERE_CENTER_DISTANCE: {
					dist_other = ((*p.samples_position_)[sv_index] - center).norm();
					dist_other *= dist_other;
				}
				break;
				case SPHERE_POWER_DISTANCE: {
					dist_other =
						((*p.samples_position_)[sv_index] - center).dot((*p.samples_position_)[sv_index] - center) -
						radius * radius;
				}
				break;
				case SPHERE_POWER_DISTANCE_SQUARED: {
					dist_other =
						((*p.samples_position_)[sv_index] - center).dot((*p.samples_position_)[sv_index] - center) -
						radius * radius;
					dist_other *= dist_other;
				}
				break;
				default:
					break;
				}
				dist_other *= a;
				Scalar dist = membership * (dist_sqem + p.sqem_clustering_lambda_ * dist_other);
				(*p.samples_vertex_error_)[sv_index] = dist;
				cluster_error += dist;
				return true;
			});
			(*p.spheres_error_)[v_index] = cluster_error / (*p.spheres_cluster_area_)[v_index];
			(*p.spheres_error_not_normalized_)[v_index] = cluster_error;

			return true;
		});

		p.min_error_ = std::numeric_limits<Scalar>::max();
		p.max_error_ = std::numeric_limits<Scalar>::min();
		p.max_error_sphere_ = PVertex();
		p.total_error_ = 0.0;
		p.total_error_not_normalized_ = 0.0;

		foreach_cell(*p.spheres_, [&](PVertex v) -> bool {
			Scalar error = value<Scalar>(*p.spheres_, p.spheres_error_, v);
			Scalar error_not_normalized = value<Scalar>(*p.spheres_, p.spheres_error_not_normalized_, v);
			p.min_error_ = std::min(p.min_error_, error);
			if (error > p.max_error_)
			{
				p.max_error_ = error;
				p.max_error_sphere_ = v;
			}
			p.total_error_ += error;
			p.total_error_not_normalized_ += error_not_normalized;
			return true;
		});

		p.total_error_diff_ = fabs(p.total_error_ - p.last_total_error_);
		p.last_total_error_ = p.total_error_;
	}

	void update_spheres_color(SurfaceParameters& p)
	{
		parallel_foreach_cell(*p.spheres_, [&](PVertex v) -> bool {
			uint32 v_index = index_of(*p.spheres_, v);
			if (p.error_as_spheres_color_)
				(*p.spheres_color_)[v_index] =
					color_map((*p.spheres_error_)[v_index], p.min_error_, p.max_error_, p.spheres_transparency_);
			else
			{
				const Vec3& c = (*p.spheres_cluster_color_)[v_index];
				(*p.spheres_color_)[v_index] = Vec4(c.x(), c.y(), c.z(), p.spheres_transparency_);
			}
			return true;
		});
	}

	Vec4 color_map(Scalar x, Scalar min, Scalar max, float32 transparency = 1.0)
	{
		x = (x - min) / (max - min);
		x = std::clamp(x, 0.0, 1.0);

		Scalar x2 = 2.0 * x;
		switch (int(std::floor(std::max(0.0, x2 + 1.0))))
		{
		case 0:
			return Vec4(0.0, 0.0, 1.0, transparency);
		case 1:
			return Vec4(x2, x2, 1.0, transparency);
		case 2:
			return Vec4(1.0, 2.0 - x2, 2.0 - x2, transparency);
		}
		return Vec4(1.0, 0.0, 0.0, transparency);
	}

	void update_sphere_euclidean_distance(SurfaceParameters& p, PVertex sphere)
	{
		uint32 sphere_index = index_of(*p.spheres_, sphere);

		const std::vector<PVertex>& cluster = (*p.spheres_cluster_)[sphere_index];
		Vec3 c = (*p.spheres_position_)[sphere_index];
		Scalar r = (*p.spheres_radius_)[sphere_index];
		// Verify if the SQEM is well conditioned
		Spherical_Quadric q;
		for (PVertex v : cluster)
		{
			uint32 v_index = index_of(*p.samples_, v);
			q += (*p.samples_quadric_)[v_index];
		}
		if (q.well_conditioned())
		{
			Eigen::MatrixXd J(2 * cluster.size(), 4);
			J.setZero();
			Eigen::VectorXd b(2 * cluster.size());
			b.setZero();
			uint32 idx = 0;
			Eigen::VectorXd s(4);
			s << c[0], c[1], c[2], r;
			for (uint32 i = 0; i < 10; ++i)
			{
				idx = 0;
				for (PVertex v : cluster)
				{
					uint32 v_index = index_of(*p.samples_, v);
					const Vec3& pos = (*p.projected_samples_position_)[v_index];

					// SQEM energy
					Eigen::Vector4d lhs = Eigen::Vector4d::Zero();
					Scalar rhs = 0.0;
					const Vec3& n = (*p.projected_samples_normal_)[v_index];
					Vec4 n4 = Vec4(n.x(), n.y(), n.z(), 1.0);
					Scalar a = 1;
					lhs += -n4 * a;
					rhs += -1.0 * ((pos - Vec3(s(0), s(1), s(2))).dot(n) - s(3)) * a;
					J.row(idx) = lhs;
					b(idx) = rhs;
					//}
					++idx;

					// distance energy
					Vec3 d = pos - Vec3(s(0), s(1), s(2));
					Scalar l = d.norm();
					if (p.point_cloud_mode_)
					{
						Scalar a = sqrt((*p.surface_vertex_area_pc_)[v_index]);
						J.row(idx) =
							Eigen::Vector4d(-(d[0] / l), -(d[1] / l), -(d[2] / l), -1.0) * a * p.sqem_update_lambda_;
						b(idx) = -(l - s(3)) * a * p.sqem_update_lambda_; // scale the row by the update lambda
					}
					else
					{
						Scalar a = 1;
						J.row(idx) =
							Eigen::Vector4d(-(d[0] / l), -(d[1] / l), -(d[2] / l), -1.0) * a * p.sqem_update_lambda_;
						b(idx) = -(l - s(3)) * a * p.sqem_update_lambda_; // scale the row by the update lambda
					}
					++idx;
				};

				Eigen::LDLT<Eigen::MatrixXd> solver(J.transpose() * J);
				Eigen::VectorXd delta_s = solver.solve(J.transpose() * b);
				s += delta_s;
				if (delta_s.norm() < 1e-6) // stop early if converged
					break;
			}

			c = s.head<3>();
			r = s[3];
		}
		else
		{
			std::cout << "Sphere " << sphere_index << " is not well conditioned, using shrinking ball" << std::endl;
			// apply shrinking ball
			std::pair<uint32, Vec3> bvh_res;
			p.surface_bvh_->closest_point(c, &bvh_res);
			Vec3 closest_point_position = bvh_res.second;
			Vec3 closest_point_dir = (closest_point_position - c).normalized();

			const Vec3& closest_face_normal =
				value<Vec3>(*p.surface_, p.surface_face_normal_, p.surface_bvh_faces_[bvh_res.first]);
			// TODO: exterior detection is not reliable
			if ((*p.inside_tester_)(Point_3(c.x(), c.y(), c.z())) == CGAL::ON_UNBOUNDED_SIDE)
				closest_point_dir = -closest_point_dir;

			auto [center, radius, secondary] = geometry::shrinking_ball_center(
				*p.samples_, closest_point_position, closest_point_dir, p.projected_samples_position_.get(),
				p.projected_surface_kdt_, p.projected_kdt_vertices_, p.point_cloud_mode_, 0.25f);
			c = center;
			r = radius;
		}
		(*p.spheres_position_)[sphere_index] = c;
		(*p.spheres_radius_)[sphere_index] = r;
	}

	
	void update_fuzzy_sphere_euclidean_distance(SurfaceParameters& p, PVertex sphere)
	{
		uint32 sphere_index = index_of(*p.spheres_, sphere);

		const std::vector<PVertex>& cluster = (*p.spheres_cluster_)[sphere_index];
		Vec3 c = (*p.spheres_position_)[sphere_index];
		Scalar r = (*p.spheres_radius_)[sphere_index];
		Vec4 s;
		s << c(0), c(1), c(2), r;

		for (uint32 i = 0; i < 10; ++i)
		{
			std::vector<Eigen::Vector4d> rows;
			std::vector<Scalar> rhs;
			rows.reserve(2048);
			rhs.reserve(2048);
			foreach_cell(*p.samples_, [&](PVertex v) -> bool {
				uint32 v_index = index_of(*p.samples_, v);
				auto& mmap = (*p.samples_membership_)[v_index];
				auto it = mmap.find(sphere_index);
				if (it == mmap.end())
					return true;
				const Scalar weight = std::max(0.0, it->second);
				const Scalar sw = std::sqrt(weight);
				const Vec3& pos = (*p.projected_samples_position_)[v_index];

				// SQEM energy
				const Vec3& n = (*p.projected_samples_normal_)[v_index];
				Vec4 n4 = Vec4(n.x(), n.y(), n.z(), 1.0);
				Scalar a = 1;
				Vec4 lhs=  -n4* a;
				Scalar r1 = - 1.0 * ((pos - Vec3(s(0), s(1), s(2))).dot(n) - s(3)) * a;
				rows.push_back(sw* lhs);
				rhs.push_back(sw * r1);

				// distance energy
				Vec3 d = pos - Vec3(s(0), s(1), s(2));
				Scalar l = d.norm();

				lhs = Vec4(-(d[0] / l), -(d[1] / l), -(d[2] / l), -1.0) * a * p.sqem_update_lambda_;
				Scalar r2 =  -(l - s(3)) * a * p.sqem_update_lambda_; // scale the row by the update lambda
				rows.push_back(lhs * sw);
				rhs.push_back(sw * r2);
				return true;
			});
			Eigen::MatrixXd J(rows.size(), 4);
			Eigen::VectorXd b(rhs.size());
			for (uint32 j = 0; j < rows.size(); ++j)
			{
				J.row(j) = rows[j];
				b(j) = rhs[j];
			}
			Eigen::LDLT<Eigen::MatrixXd> solver(J.transpose() * J);
			Eigen::VectorXd delta_s = solver.solve(J.transpose() * b);
			s += delta_s;
			if (delta_s.norm() < 1e-6) // stop early if converged
				break;
		}

		(*p.spheres_position_)[sphere_index] = s.head<3>();
		(*p.spheres_radius_)[sphere_index] = s[3];
	}

	void update_sphere_center_distance(SurfaceParameters& p, PVertex sphere)
	{
		uint32 sphere_index = index_of(*p.spheres_, sphere);

		const std::vector<PVertex>& cluster = (*p.spheres_cluster_)[sphere_index];
		Vec3 c = (*p.spheres_position_)[sphere_index];
		Scalar r = (*p.spheres_radius_)[sphere_index];
		// Verify if the SQEM is well conditioned
		Spherical_Quadric q;
		Scalar area = 0.0;
		Vec3 h;
		h.setZero();
		for (PVertex v : cluster)
		{
			uint32 v_index = index_of(*p.samples_, v);
			q += (*p.samples_quadric_)[v_index];
			h += 1.0 * (*p.samples_position_)[v_index];
			area += 1.0;
		}
		if (q.well_conditioned())
		{
			Mat4 As = Mat4::Identity();
			As.block<3, 3>(0, 0) = 2 * Mat3::Identity() * area;
			As(3, 3) = 0;
			Vec4 bs;
			bs.head<3>() = 2 * h;
			bs(3) = 0;
			Mat4 A = q._A + p.sqem_clustering_lambda_ * As;
			Vec4 b = q._b + p.sqem_clustering_lambda_ * bs;
			Vec4 s = A.ldlt().solve(b);
			c = s.head<3>();
			r = s(3);
		}
		else
		{
			std::cout << "Sphere " << sphere_index << " is not well conditioned, using shrinking ball" << std::endl;
			// apply shrinking ball
			std::pair<uint32, Vec3> bvh_res;
			p.surface_bvh_->closest_point(c, &bvh_res);
			Vec3 closest_point_position = bvh_res.second;
			Vec3 closest_point_dir = (closest_point_position - c).normalized();

			const Vec3& closest_face_normal =
				value<Vec3>(*p.surface_, p.surface_face_normal_, p.surface_bvh_faces_[bvh_res.first]);
			// TODO: exterior detection is not reliable
			if ((*p.inside_tester_)(Point_3(c.x(), c.y(), c.z())) == CGAL::ON_UNBOUNDED_SIDE)
				closest_point_dir = -closest_point_dir;

			auto [center, radius, secondary] = geometry::shrinking_ball_center(
				*p.samples_, closest_point_position, closest_point_dir, p.projected_samples_position_.get(),
				p.projected_surface_kdt_, p.projected_kdt_vertices_, p.point_cloud_mode_, 0.25f);
			c = center;
			r = radius;
		}
		(*p.spheres_position_)[sphere_index] = c;
		(*p.spheres_radius_)[sphere_index] = r;
	}

	void update_power_distance(SurfaceParameters& p, PVertex sphere)
	{
		uint32 sphere_index = index_of(*p.spheres_, sphere);

		const std::vector<PVertex>& cluster = (*p.spheres_cluster_)[sphere_index];
		Vec3 c = (*p.spheres_position_)[sphere_index];
		Scalar r = (*p.spheres_radius_)[sphere_index];
		// Verify if the SQEM is well conditioned
		Spherical_Quadric q;
		Scalar area = 0.0;
		Vec3 h;
		h.setZero();
		for (PVertex v : cluster)
		{
			uint32 v_index = index_of(*p.samples_, v);
			q += (*p.samples_quadric_)[v_index];
			h += 1.0 * (*p.samples_position_)[v_index];
			area += 1.0;
		}
		if (q.well_conditioned())
		{
			Mat4 As = Mat4::Zero();
			As.block<3, 3>(0, 0) = 2 * Mat3::Identity() * area;
			As(3, 3) = -2 * area;
			Vec4 bs;
			bs.head<3>() = 2 * h;
			bs(3) = 0;
			Mat4 A = q._A + p.sqem_clustering_lambda_ * As;
			Vec4 b = q._b + p.sqem_clustering_lambda_ * bs;
			Vec4 s = A.ldlt().solve(b);
			c = s.head<3>();
			r = s(3);
		}
		else
		{
			std::cout << "Sphere " << sphere_index << " is not well conditioned, using shrinking ball" << std::endl;
			// apply shrinking ball
			std::pair<uint32, Vec3> bvh_res;
			p.surface_bvh_->closest_point(c, &bvh_res);
			Vec3 closest_point_position = bvh_res.second;
			Vec3 closest_point_dir = (closest_point_position - c).normalized();

			const Vec3& closest_face_normal =
				value<Vec3>(*p.surface_, p.surface_face_normal_, p.surface_bvh_faces_[bvh_res.first]);
			// TODO: exterior detection is not reliable
			if ((*p.inside_tester_)(Point_3(c.x(), c.y(), c.z())) == CGAL::ON_UNBOUNDED_SIDE)
				closest_point_dir = -closest_point_dir;

			auto [center, radius, secondary] = geometry::shrinking_ball_center(
				*p.samples_, closest_point_position, closest_point_dir, p.projected_samples_position_.get(),
				p.projected_surface_kdt_, p.projected_kdt_vertices_, p.point_cloud_mode_, 0.25f);
			c = center;
			r = radius;
		}
		(*p.spheres_position_)[sphere_index] = c;
		(*p.spheres_radius_)[sphere_index] = r;
	}

	void update_power_distance_squared(SurfaceParameters& p, PVertex sphere)
	{
		uint32 sphere_index = index_of(*p.spheres_, sphere);

		const std::vector<PVertex>& cluster = (*p.spheres_cluster_)[sphere_index];
		Vec3 c = (*p.spheres_position_)[sphere_index];
		Scalar r = (*p.spheres_radius_)[sphere_index];
		// Verify if the SQEM is well conditioned
		Spherical_Quadric q;
		for (PVertex v : cluster)
		{
			uint32 v_index = index_of(*p.samples_, v);
			q += (*p.samples_quadric_)[v_index];
		}
		if (q.well_conditioned())
		{
			Eigen::MatrixXd J(2 * cluster.size(), 4);
			J.setZero();
			Eigen::VectorXd b(2 * cluster.size());
			b.setZero();
			uint32 idx = 0;
			Eigen::VectorXd s(4);
			s << c[0], c[1], c[2], r;
			for (uint32 i = 0; i < 10; ++i)
			{
				idx = 0;
				for (PVertex v : cluster)
				{
					uint32 v_index = index_of(*p.samples_, v);
					const Vec3& pos = (*p.projected_samples_position_)[v_index];

					// SQEM energy
					Eigen::Vector4d lhs = Eigen::Vector4d::Zero();
					Scalar rhs = 0.0;
					const Vec3& n = (*p.projected_samples_normal_)[v_index];
					Vec4 n4 = Vec4(n.x(), n.y(), n.z(), 1.0);
					Scalar a = 1;
					lhs += -n4 * a;
					rhs += -1.0 * ((pos - Vec3(s(0), s(1), s(2))).dot(n) - s(3)) * a;
					J.row(idx) = lhs;
					b(idx) = rhs;
					//}
					++idx;

					// distance energy
					const Vec3 vol_pos = (*p.samples_position_)[v_index];
					Vec3 d = Vec3(s(0), s(1), s(2)) - vol_pos;
					Scalar l = d.norm();
					J.row(idx) = Eigen::Vector4d(2 * d(0) / l, 2 * d(1) / l, 2 * d(2) / l, -2.0 * s(3)) * a *
								 p.sqem_update_lambda_;
					b(idx) = -(l * l - s(3) * s(3)) * a * p.sqem_update_lambda_; // scale the row by the update lambda

					++idx;
				}

				Eigen::LDLT<Eigen::MatrixXd> solver(J.transpose() * J);
				Eigen::VectorXd delta_s = solver.solve(J.transpose() * b);
				s += delta_s;
				if (delta_s.norm() < 1e-6) // stop early if converged
					break;
			}

			c = s.head<3>();
			r = s[3];
		}
		else
		{
			std::cout << "Sphere " << sphere_index << " is not well conditioned, using shrinking ball" << std::endl;
			// apply shrinking ball
			std::pair<uint32, Vec3> bvh_res;
			p.surface_bvh_->closest_point(c, &bvh_res);
			Vec3 closest_point_position = bvh_res.second;
			Vec3 closest_point_dir = (closest_point_position - c).normalized();

			const Vec3& closest_face_normal =
				value<Vec3>(*p.surface_, p.surface_face_normal_, p.surface_bvh_faces_[bvh_res.first]);
			// TODO: exterior detection is not reliable
			if ((*p.inside_tester_)(Point_3(c.x(), c.y(), c.z())) == CGAL::ON_UNBOUNDED_SIDE)
				closest_point_dir = -closest_point_dir;

			auto [center, radius, secondary] = geometry::shrinking_ball_center(
				*p.samples_, closest_point_position, closest_point_dir, p.projected_samples_position_.get(),
				p.projected_surface_kdt_, p.projected_kdt_vertices_, p.point_cloud_mode_, 0.25f);
			c = center;
			r = radius;
		}
		(*p.spheres_position_)[sphere_index] = c;
		(*p.spheres_radius_)[sphere_index] = r;
	}
	void correct_sphere(SurfaceParameters& p, PVertex v)
	{
		uint32 v_index = index_of(*p.spheres_, v);

		Vec3& c = (*p.spheres_position_)[v_index];
		Scalar& r = (*p.spheres_radius_)[v_index];

		Vec3 closest_point_position;
		Vec3 closest_point_dir;

		if (p.point_cloud_mode_)
		{
			std::pair<uint32, Scalar> k_res;
			p.projected_surface_kdt_->find_nn(c, &k_res);
			PVertex closest_vertex = p.projected_kdt_vertices_[k_res.first];
			closest_point_position = p.projected_surface_kdt_->vertex(k_res.first);
			closest_point_dir = (closest_point_position - c).normalized();

			const Vec3& closest_vertex_normal = value<Vec3>(*p.samples_, p.projected_samples_normal_, closest_vertex);
			// TODO: exterior detection is not reliable
			if (closest_point_dir.dot(closest_vertex_normal) <= 0.0)
				closest_point_dir = -closest_point_dir;
		}
		else
		{
			std::pair<uint32, Vec3> bvh_res;
			p.surface_bvh_->closest_point(c, &bvh_res);
			closest_point_position = bvh_res.second;
			closest_point_dir = (closest_point_position - c).normalized();

			const Vec3& closest_face_normal =
				value<Vec3>(*p.surface_, p.surface_face_normal_, p.surface_bvh_faces_[bvh_res.first]);
			// TODO: exterior detection is not reliable
			if (closest_point_dir.dot(closest_face_normal) <= 0.0)
				closest_point_dir = -closest_point_dir;
		}

		auto [center, radius, secondary] = geometry::shrinking_ball_center(
			*p.samples_, closest_point_position, closest_point_dir, p.projected_samples_position_.get(),
			p.projected_surface_kdt_, p.projected_kdt_vertices_, p.point_cloud_mode_, 0.25f);

		c = center;
		r = radius;
	}

	void update_spheres(SurfaceParameters& p)
	{
		// auto start = std::chrono::high_resolution_clock::now();

		//compute_clusters(p);
		compute_membership_soft(p);
		materialize_top1_labels(p);
		// switch (p.update_method_)
		// {
		// case FIT: {
		// 	parallel_foreach_cell(*p.spheres_, [&](PVertex v) -> bool {
		// 		update_sphere_fit(p, v);
		// 		return true;
		// 	});
		// 	break;
		// }
		// case SQEM: {
		parallel_foreach_cell(*p.spheres_, [&](PVertex v) -> bool {
			switch (p.distance_mode_)
			{
			case SPHERE_EUCLIDEAN_DISTANCE:
				update_fuzzy_sphere_euclidean_distance(p, v);
				break;
			case SPHERE_CENTER_DISTANCE:
				update_sphere_center_distance(p, v);
				break;
			case SPHERE_POWER_DISTANCE:
				update_power_distance(p, v);
				break;
			case SPHERE_POWER_DISTANCE_SQUARED:
				update_power_distance_squared(p, v);
				break;
			}

			if (p.sphere_correction_ && p.sphere_correction_mode_ == CORRECT_ALWAYS)
				correct_sphere(p, v);
			value<bool>(*p.spheres_, p.spheres_do_not_split_, v) = false;
			return true;
		});
		// 	break;
		// }
		// }

		compute_spheres_error(p); // compute spheres error
		// std::cout << p.total_error_not_normalized_ << std::endl;

		if (p.auto_split_ &&
			(p.total_error_diff_ < 1e-5 || p.iteration_count_ % 10 == 0)) // wait for convergence or max 10 iterations
		{
			switch (p.auto_split_mode_)
			{
			case ERROR_THRESHOLD: {
				if (p.max_error_ > p.auto_split_error_threshold_)
				{
					compute_cluster_neighbour(p);

					if (p.sphere_correction_ && p.sphere_correction_mode_ == CORRECT_ON_SPLIT)
					{
						parallel_foreach_cell(*p.spheres_, [&](PVertex v) -> bool {
							correct_sphere(p, v);
							return true;
						});
					}

					// std::cout << 0.0 << std::endl;

					std::vector<PVertex> sorted_spheres;
					sorted_spheres.reserve(p.nb_spheres_);
					foreach_cell(*p.spheres_, [&](PVertex v) -> bool {
						sorted_spheres.push_back(v);
						return true;
					});
					std::sort(sorted_spheres.begin(), sorted_spheres.end(), [&](PVertex a, PVertex b) {
						return value<Scalar>(*p.spheres_, p.spheres_error_, a) >
							   value<Scalar>(*p.spheres_, p.spheres_error_, b);
					});
					// do not split more than 20% of the spheres at once (and 10 at most)
					uint32 to_split_max = std::min(uint32(std::ceil(p.nb_spheres_ * 0.2)), 10u);
					for (PVertex sphere : sorted_spheres)
					{
						uint32 s_index = index_of(*p.spheres_, sphere);
						Scalar error = (*p.spheres_error_)[s_index];
						if (error < p.auto_split_error_threshold_)
							break;
						if (to_split_max > 0 && !(*p.spheres_do_not_split_)[s_index]) // do not split neighbor spheres
						{
							for (PVertex neighbor : (*p.spheres_neighbor_clusters_)[s_index])
								value<bool>(*p.spheres_, p.spheres_do_not_split_, neighbor) = true;
							split_sphere(p, sphere);
							--to_split_max;
						}
					}
				}
			}
			break;
			case MAX_NB_SPHERES: {
				if (p.nb_spheres_ < p.auto_split_max_nb_spheres_)
				{
					compute_cluster_neighbour(p);

					if (p.sphere_correction_ && p.sphere_correction_mode_ == CORRECT_ON_SPLIT)
					{
						parallel_foreach_cell(*p.spheres_, [&](PVertex v) -> bool {
							correct_sphere(p, v);
							return true;
						});
					}

					// std::cout << 0.0 << std::endl;

					std::vector<PVertex> sorted_spheres;
					sorted_spheres.reserve(p.nb_spheres_);
					foreach_cell(*p.spheres_, [&](PVertex v) -> bool {
						sorted_spheres.push_back(v);
						return true;
					});
					std::sort(sorted_spheres.begin(), sorted_spheres.end(), [&](PVertex a, PVertex b) {
						return value<Scalar>(*p.spheres_, p.spheres_error_, a) >
							   value<Scalar>(*p.spheres_, p.spheres_error_, b);
					});
					// do not split more than 20% of the spheres at once (and 10 at most)
					uint32 to_split_max = std::min(uint32(std::ceil(p.nb_spheres_ * 0.2)), 10u);
					for (PVertex sphere : sorted_spheres)
					{
						uint32 s_index = index_of(*p.spheres_, sphere);
						if (p.auto_split_max_nb_spheres_ - p.nb_spheres_ <= 0)
							break;
						if (to_split_max > 0 && !(*p.spheres_do_not_split_)[s_index]) // do not split neighbor spheres
						{
							for (PVertex neighbor : (*p.spheres_neighbor_clusters_)[s_index])
								value<bool>(*p.spheres_, p.spheres_do_not_split_, neighbor) = true;
							split_sphere(p, sphere);
							--to_split_max;
						}
					}
				}
			}
			break;
			}
		}

		// auto end = std::chrono::high_resolution_clock::now();
		// std::cout << "Update spheres: " << std::chrono::duration<Scalar>(end - start).count() << "s" << std::endl;

		if (!p.running_)
			update_render_data(p);
	}

	void split_sphere(SurfaceParameters& p, PVertex v)
	{
		uint32 v_index = index_of(*p.spheres_, v);

		// find the vertex of the cluster with maximal error w.r.t. the sphere

		PVertex max_error_vertex;
		Scalar max_error_vertex_error = 0.0;
		const std::vector<PVertex>& cluster = (*p.spheres_cluster_)[v_index];
		for (PVertex sv : cluster)
		{
			Scalar error = value<Scalar>(*p.samples_, p.samples_vertex_error_, sv);
			if (error > max_error_vertex_error)
			{
				max_error_vertex_error = error;
				max_error_vertex = sv;
			}
		}
		uint32 max_error_vertex_index = index_of(*p.samples_, max_error_vertex);

		// insert a new sphere for this vertex

		PVertex new_sphere = add_vertex(*p.spheres_);
		p.nb_spheres_++;
		uint32 new_sphere_index = index_of(*p.spheres_, new_sphere);
		(*p.spheres_position_)[new_sphere_index] = (*p.medial_axis_position_)[max_error_vertex_index];
		(*p.spheres_radius_)[new_sphere_index] = (*p.medial_axis_radius_)[max_error_vertex_index];

		// sample a random bright pastel color for the new sphere
		(*p.spheres_cluster_color_)[new_sphere_index] = Vec3(
			0.5 + 0.5 * (rand() % 256) / 256.0, 0.5 + 0.5 * (rand() % 256) / 256.0, 0.5 + 0.5 * (rand() % 256) / 256.0);
	}

	void remove_sphere(SurfaceParameters& p, PVertex v)
	{
		const std::vector<PVertex>& cluster = value<std::vector<PVertex>>(*p.spheres_, p.spheres_cluster_, v);
		for (PVertex sv : cluster)
			value<PVertex>(*p.samples_, p.samples_vertex_sphere_, sv) = PVertex();
		remove_vertex(*p.spheres_, v);
		p.nb_spheres_--;
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

	
	void compute_cluster_neighbour(SurfaceParameters& p)
	{
		// clean neighbor clusters sets
		parallel_foreach_cell(*p.spheres_, [&](PVertex v) -> bool {
			uint32 v_index = index_of(*p.spheres_, v);
			(*p.spheres_neighbor_clusters_)[v_index].clear();
			return true;
		});

		// compute neighbor clusters
		foreach_cell(*p.samples_, [&](PVertex v) -> bool {
			uint32 v_index = index_of(*p.samples_, v);
			PVertex v_sphere = (*p.samples_vertex_sphere_)[v_index];
			for (PVertex w : (*p.samples_vertex_knn_)[v_index])
			{
				PVertex w_sphere = value<PVertex>(*p.samples_, p.samples_vertex_sphere_, w);
				if (v_sphere.is_valid() && w_sphere.is_valid() && v_sphere != w_sphere)
				{
					value<std::set<PVertex>>(*p.spheres_, p.spheres_neighbor_clusters_, v_sphere).insert(w_sphere);
					value<std::set<PVertex>>(*p.spheres_, p.spheres_neighbor_clusters_, w_sphere).insert(v_sphere);
				}
			}
			return true;
		});

	}

	void compute_skeleton(SurfaceParameters& p, bool compute_neighbor_clusters_only = false)
	{
		// clean neighbor clusters sets
		parallel_foreach_cell(*p.spheres_, [&](PVertex v) -> bool {
			uint32 v_index = index_of(*p.spheres_, v);
			(*p.spheres_neighbor_clusters_)[v_index].clear();
			(*p.spheres_cluster_)[v_index].clear();
			(*p.spheres_cluster_area_)[v_index] = 0.0;
			return true;
		});

		p.samples_vertex_sphere_->fill(PVertex());

		MeshData<POINTS>& md = points_provider_->mesh_data(*p.spheres_);
		if (p.nb_spheres_ == 0)
			return;

		// auto start = std::chrono::high_resolution_clock::now();

		parallel_foreach_cell(*p.samples_, [&](PVertex v) -> bool {
			uint32 v_index = index_of(*p.samples_, v);

			// Each volumn sample has weight 1
			Scalar a = 1.0;

			Scalar min_distance = std::numeric_limits<Scalar>::max();
			PVertex closest_sphere;
			uint32 closest_sphere_index;

			foreach_cell(*p.spheres_, [&](PVertex pv) {
				uint32 pv_index = index_of(*p.spheres_, pv);
				const Vec3& center = (*p.spheres_position_)[pv_index];
				Scalar radius = (*p.spheres_radius_)[pv_index];

				Scalar power_distance =
					((*p.samples_position_)[v_index] - center).dot((*p.samples_position_)[v_index] - center)-radius * radius;
				
				power_distance *= a;
				
				if (power_distance < min_distance)
				{
					min_distance = power_distance;
					closest_sphere = pv;
					closest_sphere_index = pv_index;
				}
				return true;
			});

			value<PVertex>(*p.samples_, p.samples_vertex_sphere_, v) = closest_sphere;

			std::lock_guard<std::mutex> lock(spheres_mutex_[closest_sphere_index % spheres_mutex_.size()]);
			value<std::vector<PVertex>>(*p.spheres_, p.spheres_cluster_, closest_sphere).push_back(v);
			value<Scalar>(*p.spheres_, p.spheres_cluster_area_, closest_sphere) += a;

			return true;
		});
	
		// compute neighbor clusters
		foreach_cell(*p.samples_, [&](PVertex v) -> bool {
			uint32 v_index = index_of(*p.samples_, v);
			PVertex v_sphere = (*p.samples_vertex_sphere_)[v_index];
			for (PVertex w : (*p.samples_vertex_knn_)[v_index])
			{
				PVertex w_sphere = value<PVertex>(*p.samples_, p.samples_vertex_sphere_, w);
				if (v_sphere.is_valid() && w_sphere.is_valid() && v_sphere != w_sphere)
				{
					value<std::set<PVertex>>(*p.spheres_, p.spheres_neighbor_clusters_, v_sphere).insert(w_sphere);
					value<std::set<PVertex>>(*p.spheres_, p.spheres_neighbor_clusters_, w_sphere).insert(v_sphere);
				}
			}
			return true;
		});


		clear(*p.skeleton_);

		auto spheres_skeleton_vertex_map =
			add_attribute<NMVertex, PVertex>(*p.spheres_, "__spheres_skeleton_vertex_map");

		foreach_cell(*p.spheres_, [&](PVertex pv) -> bool {
			uint32 pv_index = index_of(*p.spheres_, pv);
			NMVertex nmv = add_vertex(*p.skeleton_);
			value<Vec3>(*p.skeleton_, p.skeleton_position_, nmv) = (*p.spheres_position_)[pv_index];
			(*spheres_skeleton_vertex_map)[pv_index] = nmv;
			return true;
		});

		std::unordered_map<std::pair<uint32, uint32>, NMEdge, edge_hash, edge_equal> edge_indices;

		foreach_cell(*p.spheres_, [&](PVertex pv) -> bool {
			uint32 pv_index = index_of(*p.spheres_, pv);
			NMVertex nmv1 = (*spheres_skeleton_vertex_map)[pv_index];
			const std::set<PVertex>& neighbors = (*p.spheres_neighbor_clusters_)[pv_index];
			for (PVertex neighbor : neighbors)
			{
				NMVertex nmv2 = value<NMVertex>(*p.spheres_, spheres_skeleton_vertex_map, neighbor);
				std::vector<NMVertex> av = adjacent_vertices_through_edge(*p.skeleton_, nmv1);
				if (std::find(av.begin(), av.end(), nmv2) == av.end())
				{
					NMEdge e = add_edge(*p.skeleton_, nmv1, nmv2);
					edge_indices[{index_of(*p.skeleton_, nmv1), index_of(*p.skeleton_, nmv2)}] = e;
				}
			}
			return true;
		});

		foreach_cell(*p.spheres_, [&](PVertex pv) -> bool {
			NMVertex nmv1 = value<NMVertex>(*p.spheres_, spheres_skeleton_vertex_map, pv);
			const std::set<PVertex>& n_pv = value<std::set<PVertex>>(*p.spheres_, p.spheres_neighbor_clusters_, pv);
			for (const PVertex& ne1 : n_pv)
			{
				NMVertex nmv2 = value<NMVertex>(*p.spheres_, spheres_skeleton_vertex_map, ne1);
				const std::set<PVertex>& ne_ne1 =
					value<std::set<PVertex>>(*p.spheres_, p.spheres_neighbor_clusters_, ne1);
				for (const PVertex& ne2 : ne_ne1)
				{
					if (n_pv.find(ne2) != n_pv.end())
					{
						NMVertex nmv3 = value<NMVertex>(*p.spheres_, spheres_skeleton_vertex_map, ne2);
						if (index_of(*p.skeleton_, nmv1) < index_of(*p.skeleton_, nmv2) &&
							index_of(*p.skeleton_, nmv2) < index_of(*p.skeleton_, nmv3))
						{
							std::vector<NMEdge> edges;
							edges.reserve(3);
							edges.push_back(edge_indices[{index_of(*p.skeleton_, nmv1), index_of(*p.skeleton_, nmv2)}]);
							edges.push_back(edge_indices[{index_of(*p.skeleton_, nmv2), index_of(*p.skeleton_, nmv3)}]);
							edges.push_back(edge_indices[{index_of(*p.skeleton_, nmv1), index_of(*p.skeleton_, nmv3)}]);
							add_face(*p.skeleton_, edges);
						}
					}
				}
			}
			return true;
		});

		remove_attribute<PVertex>(*p.spheres_, spheres_skeleton_vertex_map);
	}
	void update_samples_color(SurfaceParameters& p)
	{
		parallel_foreach_cell(*p.samples_, [&](PVertex v) -> bool {
			uint32 v_index = index_of(*p.samples_, v);
			auto& mmap = (*p.samples_membership_)[v_index];
			Vec4 color = Vec4(0, 0, 0, 1);
			for (auto it = mmap.begin(); it != mmap.end(); ++it)
			{
				auto c = (*p.spheres_color_)[it->first];
				color += it->second * c;
			}
			color(3) = 0.8;
			(*p.samples_vertex_color_)[v_index] = color;
			return true;
		});
		
	}

protected:
	void init() override
	{
		surface_provider_ = static_cast<ui::MeshProvider<SURFACE>*>(
			app_.module("MeshProvider (" + std::string{mesh_traits<SURFACE>::name} + ")"));

		points_provider_ = static_cast<ui::MeshProvider<POINTS>*>(
			app_.module("MeshProvider (" + std::string{mesh_traits<POINTS>::name} + ")"));

		non_manifold_provider_ = static_cast<ui::MeshProvider<NONMANIFOLD>*>(
			app_.module("MeshProvider (" + std::string{mesh_traits<NONMANIFOLD>::name} + ")"));

		timer_connection_ = boost::synapse::connect<App::timer_tick>(&app_, [this]() {
			SurfaceParameters& p = surface_parameters_[selected_surface_];
			update_render_data(p);
		});
	}

	void update_render_data(SurfaceParameters& p)
	{
		if (p.running_)
		{
			std::lock_guard<std::mutex> lock(p.mutex_);
			points_provider_->emit_connectivity_changed(*p.spheres_);
			points_provider_->emit_attribute_changed(*p.spheres_, p.spheres_position_.get());
			points_provider_->emit_attribute_changed(*p.spheres_, p.spheres_radius_.get());

			update_spheres_color(p);
			points_provider_->emit_attribute_changed(*p.spheres_, p.spheres_color_.get());
			update_samples_color(p);
			points_provider_->emit_attribute_changed(*p.samples_, p.samples_vertex_color_.get());

			compute_skeleton(p);
		}
		else
		{
			points_provider_->emit_connectivity_changed(*p.spheres_);
			points_provider_->emit_attribute_changed(*p.spheres_, p.spheres_position_.get());
			points_provider_->emit_attribute_changed(*p.spheres_, p.spheres_radius_.get());

			update_spheres_color(p);
			points_provider_->emit_attribute_changed(*p.spheres_, p.spheres_color_.get());
			update_samples_color(p);
			points_provider_->emit_attribute_changed(*p.samples_, p.samples_vertex_color_.get());

			compute_skeleton(p);
		}

		non_manifold_provider_->emit_connectivity_changed(*p.skeleton_);
		non_manifold_provider_->emit_attribute_changed(*p.skeleton_, p.skeleton_position_.get());
	}

	void start_spheres_update(SurfaceParameters& p)
	{
		p.running_ = true;
		p.iteration_count_ = 0;
		p.total_error_diff_ = 0.0;
		p.last_total_error_ = std::numeric_limits<Scalar>::max();

		launch_thread([&]() {
			auto start = std::chrono::high_resolution_clock::now();
			while (true)
			{
				{
					std::lock_guard<std::mutex> lock(p.mutex_);
					update_spheres(p);
					p.iteration_count_++;
				}
				if (p.slow_down_)
					std::this_thread::sleep_for(std::chrono::microseconds(1000000 / p.update_rate_));
				else
					std::this_thread::yield();

				if (p.auto_stop_)
				{
					if (p.total_error_diff_ < 1e-5)
					{
						switch (p.auto_split_mode_)
						{
						case ERROR_THRESHOLD:
							if (p.max_error_ < p.auto_split_error_threshold_)
								p.stopping_ = true;
							break;
						case MAX_NB_SPHERES:
							if (p.nb_spheres_ >= p.auto_split_max_nb_spheres_)
								p.stopping_ = true;
							break;
						}
					}
				}
				if (p.stopping_)
				{
					p.stopping_ = false;
					p.running_ = false;
					break;
				}
			}
			auto end = std::chrono::high_resolution_clock::now();
			std::cout << "Sphere optimizations time: " << std::chrono::duration<Scalar>(end - start).count() << "s"
					  << std::endl;
			std::cout << "Nb iterations: " << p.iteration_count_ << std::endl;
		});

		app_.start_timer(100, [&]() -> bool { return !p.running_; });
	}

	void stop_spheres_update(SurfaceParameters& p)
	{
		p.stopping_ = true;
	}

	void key_press_event(View* view, int32 key_code) override
	{
		SurfaceParameters& p = surface_parameters_[selected_surface_];

		if (key_code == GLFW_KEY_U)
		{
			if (!p.running_)
			{
				std::lock_guard<std::mutex> lock(p.mutex_);
				update_render_data(p);
			}
		}
		else if (key_code == GLFW_KEY_I || key_code == GLFW_KEY_S || key_code == GLFW_KEY_D)
		{
			int32 x = view->mouse_x();
			int32 y = view->mouse_y();

			// minwindef.h (included by MSVC) defines near and far macros, which is why the underscore is needed
			rendering::GLVec3d near_ = view->unproject(x, y, 0.0);
			rendering::GLVec3d far_d = view->unproject(x, y, 1.0);
			Vec3 A{near_.x(), near_.y(), near_.z()};
			Vec3 B{far_d.x(), far_d.y(), far_d.z()};

			Vec3 picked_sphere_center;
			foreach_cell(*p.spheres_, [&](PVertex v) -> bool {
				if (!picked_sphere_.is_valid())
				{
					picked_sphere_ = v;
					picked_sphere_center = value<Vec3>(*p.spheres_, p.spheres_position_, picked_sphere_);
					return true;
				}
				const Vec3& sp = value<Vec3>(*p.spheres_, p.spheres_position_, v);
				if (geometry::squared_distance_line_point(A, B, sp) <
					geometry::squared_distance_line_point(A, B, picked_sphere_center))
				{
					picked_sphere_ = v;
					picked_sphere_center = sp;
				}
				return true;
			});

			if (key_code == GLFW_KEY_S && picked_sphere_.is_valid())
			{
				std::lock_guard<std::mutex> lock(p.mutex_);
				split_sphere(p, picked_sphere_);
				//compute_clusters(p);
				compute_membership_soft(p);
				materialize_top1_labels(p);
				compute_spheres_error(p);
				if (!p.running_)
					update_render_data(p);
			}
			else if (key_code == GLFW_KEY_D && picked_sphere_.is_valid())
			{
				std::lock_guard<std::mutex> lock(p.mutex_);
				remove_sphere(p, picked_sphere_);
				//compute_clusters(p);
				compute_membership_soft(p);
				materialize_top1_labels(p);
				compute_spheres_error(p);
				if (!p.running_)
					update_render_data(p);
			}
		}
	}

	void left_panel() override
	{
		imgui_mesh_selector(surface_provider_, selected_surface_, "Surface",
							[&](SURFACE& s) { set_selected_surface(s); });

		if (selected_surface_)
		{
			SurfaceParameters& p = surface_parameters_[selected_surface_];

			imgui_combo_attribute<SVertex, Vec3>(
				*selected_surface_, p.surface_vertex_position_, "Position",
				[&](const std::shared_ptr<SAttribute<Vec3>>& attribute) { p.surface_vertex_position_ = attribute; });

			if (p.surface_vertex_position_ && !p.initialized_)
			{
				if (ImGui::Button("Init surface data"))
					init_surface_data(*selected_surface_);
			}

			if (p.initialized_)
			{
				ImGui::SliderFloat("Noise factor", &p.noise_factor_, 0.0f, 0.1f, "%.6f");
				if (ImGui::Button("Add noise"))
				{
					std::lock_guard<std::mutex> lock(p.mutex_);
					add_surface_noise(p);
				}
				ImGui::SameLine();
				if (ImGui::Button("Restore position"))
				{
					std::lock_guard<std::mutex> lock(p.mutex_);
					restore_surface_position(p);
				}

				ImGui::SliderFloat("Init dilation factor", &p.init_dilation_factor_, 1.0, 4.0);
				static uint32 init_max_nb_spheres = 1;
				ImGui::InputScalar("Init nb spheres", ImGuiDataType_U32, &init_max_nb_spheres);
				if (ImGui::Button("Init spheres"))
				{
					std::lock_guard<std::mutex> lock(p.mutex_);
					init_spheres(*selected_surface_, init_max_nb_spheres);
				}

				// if (ImGui::SliderFloat("Samples min radius", &p.filter_radius_threshold_, 0.0f, 1.0f, "%.6f",
				// 					   ImGuiSliderFlags_Logarithmic))
				// {
				// 	std::lock_guard<std::mutex> lock(p.mutex_);
				// 	filter_medial_samples(p);
				// }
				// if (ImGui::SliderFloat("Samples min angle", &p.filter_angle_threshold_, 0.0f, M_PI, "%.6f"))

				// {
				// 	std::lock_guard<std::mutex> lock(p.mutex_);
				// 	filter_medial_samples(p);
				// }

				if (ImGui::Checkbox("Point cloud mode", &p.point_cloud_mode_))
				{
					std::lock_guard<std::mutex> lock(p.mutex_);
					compute_quadrics(p);
				}

				// ImGui::RadioButton("Fit", (int*)&p.update_method_, FIT);
				// ImGui::SameLine();
				// ImGui::RadioButton("SQEM", (int*)&p.update_method_, SQEM);

				// if (p.update_method_ == SQEM)
				// {
				static bool sync_lambda = true;
				ImGui::Checkbox("Sync lambda", &sync_lambda);
				if (ImGui::SliderFloat("update lambda", &p.sqem_update_lambda_, 0.0f, 2.0f, "%.6f"))
				{
					if (sync_lambda)
						p.sqem_clustering_lambda_ = p.sqem_update_lambda_;
				}
				if (ImGui::SliderFloat("clustering lambda", &p.sqem_clustering_lambda_, 0.0f, 2.0f, "%.6f"))
				{
					if (sync_lambda)
						p.sqem_update_lambda_ = p.sqem_clustering_lambda_;
				}
				ImGui::SliderFloat("Tau", &p.tau_, 0.0001f, 1.0f, "%.6f");
				ImGui::SliderFloat("Threshold", &p.theta_gap_log_, 0.01f,5.0f, "%.2f");
				if (ImGui::Button("Print weight"))
				{
					uint32 count = 0;
					std::unordered_map<uint32, Scalar> weight_sum;
					foreach_cell(*p.samples_, [&](PVertex v) -> bool {
						uint32 v_index = index_of(*p.samples_, v);
						auto& mmap = (*p.samples_membership_)[v_index];
						Scalar sum = 0.0;
						for (auto it = mmap.begin(); it != mmap.end(); ++it)
							weight_sum[it->first] += it->second;
						return true;
					});
					for (auto it : weight_sum)
						std::cout << "Sphere " << it.first << " weight sum: " << it.second << std::endl;
				}
				// ImGui::Checkbox("Auto split outside spheres", &p.auto_split_outside_spheres_);

				if (ImGui::Button("Update spheres"))
				{
					if (!p.running_)
					{
						std::lock_guard<std::mutex> lock(p.mutex_);
						update_spheres(p);
						compute_spheres_error(p);
						update_render_data(p);
					}
				}

				if (ImGui::Button("Compute clusters"))
				{
					if (!p.running_)
					{
						std::lock_guard<std::mutex> lock(p.mutex_);
						//compute_clusters(p);
						compute_membership_soft(p);
						materialize_top1_labels(p);
						compute_spheres_error(p);
						update_render_data(p);
					}
				}

				ImGui::Checkbox("Sphere correction step", &p.sphere_correction_);
				if (p.sphere_correction_)
				{
					ImGui::RadioButton("Always", (int*)&p.sphere_correction_mode_, CORRECT_ALWAYS);
					ImGui::SameLine();
					ImGui::RadioButton("On split", (int*)&p.sphere_correction_mode_, CORRECT_ON_SPLIT);
				}

				if (ImGui::Button("Correct spheres"))
				{
					if (!p.running_)
					{
						std::lock_guard<std::mutex> lock(p.mutex_);
						parallel_foreach_cell(*p.spheres_, [&](PVertex v) -> bool {
							correct_sphere(p, v);
							return true;
						});
						compute_spheres_error(p);
						update_render_data(p);
					}
				}

				ImGui::Checkbox("Slow down", &p.slow_down_);
				if (p.slow_down_)
					ImGui::SliderInt("Update rate", (int*)&p.update_rate_, 1, 100);
				ImGui::Separator();
				ImGui::RadioButton("Sphere Eculidean", (int*)&p.distance_mode_, SPHERE_EUCLIDEAN_DISTANCE);
				ImGui::SameLine();
				ImGui::RadioButton("Sphere Center", (int*)&p.distance_mode_, SPHERE_CENTER_DISTANCE);

				ImGui::RadioButton("Sphere Power", (int*)&p.distance_mode_, SPHERE_POWER_DISTANCE);
				ImGui::SameLine();
				ImGui::RadioButton("Sphere Power Squared", (int*)&p.distance_mode_, SPHERE_POWER_DISTANCE_SQUARED);
				if (!p.running_)
				{
					if (ImGui::Button("Start spheres update"))
						start_spheres_update(p);
				}
				else
				{
					if (ImGui::Button("Stop spheres update"))
						stop_spheres_update(p);
				}
				ImGui::Checkbox("Auto stop", &p.auto_stop_);

				ImGui::Separator();

				ImGui::Checkbox("Auto split", &p.auto_split_);
				if (p.auto_split_)
				{
					// ImGui::Checkbox("Auto simplify", &p.auto_simplify_);
					ImGui::RadioButton("Error threshold", (int*)&p.auto_split_mode_, ERROR_THRESHOLD);
					ImGui::SameLine();
					ImGui::RadioButton("Max nb sphere", (int*)&p.auto_split_mode_, MAX_NB_SPHERES);
					if (p.auto_split_mode_ == ERROR_THRESHOLD)
						ImGui::SliderFloat("Threshold", &p.auto_split_error_threshold_, 0.0f, 1.0f, "%.6f",
										   ImGuiSliderFlags_Logarithmic);
					else
						ImGui::InputScalar("Nb spheres", ImGuiDataType_U32, &p.auto_split_max_nb_spheres_);
				}

				if (ImGui::Checkbox("Error as color", &p.error_as_spheres_color_))
				{
					std::lock_guard<std::mutex> lock(p.mutex_);
					if (!p.running_)
						update_render_data(p);
				}
				if (ImGui::SliderFloat("Transparency", &p.spheres_transparency_, 0.0f, 1.0f))
				{
					std::lock_guard<std::mutex> lock(p.mutex_);
					if (!p.running_)
						update_render_data(p);
				}

				if (ImGui::Button("Split max error sphere"))
				{
					std::lock_guard<std::mutex> lock(p.mutex_);
					split_sphere(p, p.max_error_sphere_);
					//compute_clusters(p);
					compute_membership_soft(p);
					materialize_top1_labels(p);
					compute_spheres_error(p);
					if (!p.running_)
						update_render_data(p);
				}

				ImGui::Separator();

				ImGui::Text("Total error: %f", p.total_error_);
				ImGui::Text("Min error: %f", p.min_error_);
				ImGui::Text("Max error: %f", p.max_error_);

				ImGui::Separator();

				ImGui::Text("Pick the sphere under the mouse with I, split it with S, delete it with D");
				if (picked_sphere_.is_valid())
				{
					ImGui::Text("Picked sphere:");
					const Vec3& sp = value<Vec3>(*p.spheres_, p.spheres_position_, picked_sphere_);
					ImGui::Text("Index: %u", index_of(*p.spheres_, picked_sphere_));
					ImGui::Text("Center: (%f, %f, %f)", sp[0], sp[1], sp[2]);
					ImGui::Text("Radius: %f", value<Scalar>(*p.spheres_, p.spheres_radius_, picked_sphere_));
				}
			}
		}
	}

private:
	MeshProvider<SURFACE>* surface_provider_ = nullptr;
	MeshProvider<POINTS>* points_provider_ = nullptr;
	MeshProvider<NONMANIFOLD>* non_manifold_provider_ = nullptr;

	std::unordered_map<const SURFACE*, SurfaceParameters> surface_parameters_;

	SURFACE* selected_surface_ = nullptr;
	PVertex picked_sphere_;

	std::array<std::mutex, 43> spheres_mutex_;

	std::unordered_map<const SURFACE*, std::vector<std::shared_ptr<boost::synapse::connection>>> surface_connections_;
	std::shared_ptr<boost::synapse::connection> timer_connection_;
};

} // namespace ui

} // namespace cgogn

#endif // CGOGN_MODULE_VOLUMN_VMAS_H_
