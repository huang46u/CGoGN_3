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

#include <Eigen/Sparse>
#include <cgogn/geometry/algos/area.h>
#include <cgogn/geometry/algos/fitting.h>
#include <cgogn/geometry/algos/length.h>
#include <cgogn/geometry/algos/medial_axis.h>
#include <cgogn/geometry/functions/angle.h>
#include <cgogn/geometry/functions/distance.h>
#include <cgogn/geometry/types/line_quadric.h>
#include <cgogn/geometry/types/spherical_quadric.h>
#include <cgogn/geometry/types/vmas_delaunay.h>
#include <cgogn/geometry/types/volume_poisson_disk_sampler.h>
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

using geometry::Line_Quadric;
using geometry::Mat3;
using geometry::Mat4;
using geometry::Point_Type;
using geometry::Scalar;
using geometry::Spherical_Quadric;
using geometry::SQEM_CASE;
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
	using VolumePoissonDiskSampler = cgogn::geometry::VolumePoissonDiskSampler<POINTS>;
	using VMAS_Delaunay = cgogn::geometry::Vmas_Delaunay_Manager;

	enum AutoSplitMode : uint32
	{
		MAX_NB_SPHERES,
		ERROR_THRESHOLD
	};
	enum DistanceMode : uint32
	{
		SPHERE_EUCLIDEAN_DISTANCE,
		SPHERE_POWER_DISTANCE,
		LINE_QUADRIC_DISTANCE
	};

	enum CorrectionMode : uint32
	{
		CORRECT_ALWAYS,
		CORRECT_ON_SPLIT
	};

	const uint32 k = 20;

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
		acc::KDTree<3, uint32>* surface_kdt_ = nullptr;
		std::vector<SVertex> surface_kdt_vertices_;
		acc::KDTree<3, uint32>* samples_projections_kdt_ = nullptr;
		std::vector<PVertex> samples_projections_kdt_vertices_;

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
		std::shared_ptr<PAttribute<std::set<uint32>>> spheres_neighbor_clusters_ = nullptr;
		std::shared_ptr<PAttribute<bool>> spheres_do_not_split_ = nullptr;
		std::shared_ptr<PAttribute<Scalar>> spheres_error_ = nullptr;
		std::shared_ptr<PAttribute<Scalar>> spheres_error_not_normalized_ = nullptr;
		// PAttribute<Scalar>* selected_spheres_error_ = nullptr;

		// Volumn Samples
		POINTS* samples_;
		std::unique_ptr<VolumePoissonDiskSampler> volume_sampler_;
		std::shared_ptr<PAttribute<Vec3>> samples_position_ = nullptr;
		std::shared_ptr<PAttribute<Vec4>> samples_vertex_color_ = nullptr;
		std::shared_ptr<PAttribute<Vec3>> projected_samples_position_ = nullptr;
		std::shared_ptr<PAttribute<Vec3>> projected_samples_normal_ = nullptr;
		std::shared_ptr<PAttribute<Spherical_Quadric>> samples_quadric_ = nullptr;
		std::shared_ptr<PAttribute<Line_Quadric>> samples_line_quadric_ = nullptr;
		std::shared_ptr<PAttribute<Vec3>> medial_axis_position_ = nullptr;
		std::shared_ptr<PAttribute<Scalar>> medial_axis_radius_ = nullptr;
		std::shared_ptr<PAttribute<PVertex>> medial_axis_secondary_vertex_ = nullptr;
		std::shared_ptr<PAttribute<PVertex>> samples_vertex_sphere_ = nullptr;
		std::shared_ptr<PAttribute<Scalar>> samples_vertex_error_ = nullptr;
		std::shared_ptr<PAttribute<std::vector<PVertex>>> samples_projections_vertex_knn_ = nullptr;
		std::shared_ptr<PAttribute<Scalar>> samples_weight_ = nullptr;
		std::shared_ptr<PAttribute<uint8>> samples_poisson_depth_ = nullptr;
		std::vector<CellsSet<POINTS, PVertex>*> cluster_sets_;
		std::vector<CellsSet<POINTS, PVertex>*> last_cluster_sets_;
		std::unordered_map<uint32, uint32> cluster_components_;
		VMAS_Delaunay::ClusterComponentsSizes cluster_components_sizes_;
		VMAS_Delaunay::ClusterAdjacency cluster_adjacency_;
		VMAS_Delaunay::ClusterPairComponents cluster_pair_components_;
		bool cluster_cc_ready_ = false;

		uint32 samples_numbers_ = 100000;
		std::unique_ptr<Side_tester> inside_tester_;
		CGAL_SurfaceMesh cgal_surface_mesh_;
		std::unique_ptr<VMAS_Delaunay> vmas_delaunay_;

		NONMANIFOLD* skeleton_;
		std::shared_ptr<NMAttribute<Vec3>> skeleton_position_ = nullptr;

		float32 filter_radius_threshold_ = 0.0f;
		float32 filter_angle_threshold_ = 0.0f;

		float32 init_dilation_factor_ = 3.0f;

		// UpdateMethod update_method_ = SQEM;

		// bool auto_split_outside_spheres_ = false;

		bool sphere_correction_ = false;
		CorrectionMode sphere_correction_mode_ = CORRECT_ALWAYS;
		bool use_line_quadric_ = false;
		bool line_quadric_combined_ = false;
		bool auto_stop_ = false;
		bool auto_split_ = true;
		AutoSplitMode auto_split_mode_ = ERROR_THRESHOLD;
		DistanceMode distance_mode_ = LINE_QUADRIC_DISTANCE;

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

		uint32 iteration_count_ = 0;
		std::mutex mutex_;
		bool running_ = false;
		bool stopping_ = false;
		bool slow_down_ = true;
		bool local_sampling_ = false;
		uint32 update_rate_ = 20;
		// poisson disk sampling
		Scalar r_ = 0.01;
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

		for (auto& [surface, params] : surface_parameters_)
		{
			if (params.surface_bvh_)
			{
				delete params.surface_bvh_;
				params.surface_bvh_ = nullptr;
			}

			if (params.surface_kdt_)
			{
				delete params.surface_kdt_;
				params.surface_kdt_ = nullptr;
			}
			if (params.samples_projections_kdt_)
			{
				delete params.samples_projections_kdt_;
				params.samples_projections_kdt_ = nullptr;
			}
		}
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

	inline void on_accept_post(SurfaceParameters& p, PVertex& v, const Vec3& pos, uint8 depth)
	{
		uint32 vid = index_of(*p.samples_, v);
		std::pair<uint32, Vec3> bvh_res;
		p.surface_bvh_->closest_point(pos, &bvh_res);
		Vec3 closest_surface_position = bvh_res.second;

		(*p.samples_position_)[vid] = pos;
		(*p.samples_poisson_depth_)[vid] = depth;
		(*p.projected_samples_position_)[vid] = closest_surface_position;
		(*p.projected_samples_normal_)[vid] = (closest_surface_position - pos).normalized();

		// Insert into delaunay
		p.vmas_delaunay_->insert_sample(pos, vid, Point_Type::VOLUME_SAMPLE);
		p.vmas_delaunay_->insert_sample(closest_surface_position, vid, Point_Type::SURFACE_PROJECTION);
	}
	void poisson_disk_sampling(SurfaceParameters& p)
	{
		auto in_volume = [&](const Vec3& pos) -> bool { return is_inside(p, pos); };
		auto on_accept = [&](const Vec3& pos, PVertex& v, const uint8 depth) {
			cgogn_message_assert(is_inside(p, pos), "pos is not inside");
			on_accept_post(p, v, pos, depth);
		};
		p.volume_sampler_->sample_fill_at_depth(0, on_accept, in_volume);
		p.volume_sampler_->sample_fill_at_depth(1, on_accept, in_volume);
		p.volume_sampler_->sample_fill_at_depth(2, on_accept, in_volume);
		//p.volume_sampler_->sample_fill_at_depth(3, on_accept, in_volume);
		points_provider_->emit_connectivity_changed(*p.samples_);
	}

	void verify_conflicts(SurfaceParameters& p)
	{
		auto poisson_sample_radius_ = get_or_add_attribute<Scalar, PVertex>(*p.samples_, "poisson_sample_radius");
		foreach_cell(*p.samples_, [&](PVertex v) {
			const uint32 v_index = index_of(*p.samples_, v);
			const Vec3& p_pos = (*p.samples_position_)[v_index];
			Scalar R = (*poisson_sample_radius_)[v_index];
			auto adjacent_map = p.vmas_delaunay_->compute_adjacent_clusters();
			for (uint32 iv_index : adjacent_map[v_index])
			{
				if (iv_index == v_index)
					continue;
				const Vec3& q_pos = (*p.samples_position_)[iv_index];
				Scalar R_iv = (*poisson_sample_radius_)[iv_index];
				Scalar R_threshold = std::min(R, R_iv);
				if ((p_pos - q_pos).dot(p_pos - q_pos) < (R_threshold * R_threshold))
				{
					std::cout << "Conflict between " << v_index << " and " << iv_index << std::endl;
					std::cout << "Depth: " << (int)(*p.samples_poisson_depth_)[v_index] << " and "
							  << (int)(*p.samples_poisson_depth_)[iv_index] << std::endl;
					// radius
					std::cout << "Radius: " << R << " and " << R_iv << std::endl;
				}
			}
			return true;
		});
	}
	void poisson_disk_sampling_local(SurfaceParameters& p, PVertex sphere, uint32 target_count)
	{
		uint32 s_index = index_of(*p.spheres_, sphere);
		auto in_volume = [&](const Vec3& pos) -> bool { return is_inside(p, pos); };
		auto on_accept = [&](const Vec3& pos, PVertex& v, const uint8 depth) {
			cgogn_message_assert(is_inside(p, pos), "pos is not inside");
			uint32 vid = index_of(*p.samples_, v);
			std::pair<uint32, Vec3> bvh_res;
			p.surface_bvh_->closest_point(pos, &bvh_res);
			Vec3 closest_surface_position = bvh_res.second;

			SFace closest_face = p.surface_bvh_faces_[bvh_res.first];
			(*p.samples_position_)[vid] = pos;
			(*p.samples_poisson_depth_)[vid] = depth;
			(*p.projected_samples_position_)[vid] = closest_surface_position;
			(*p.projected_samples_normal_)[vid] = (closest_surface_position - pos).normalized();
			// compute quadric
			Line_Quadric& lq = (*p.samples_line_quadric_)[vid];
			lq.zero();

			Spherical_Quadric& q = (*p.samples_quadric_)[vid];
			q.clear();
			const Vec3& n = (*p.projected_samples_normal_)[vid];
			lq = Line_Quadric(pos, n);
			q = Spherical_Quadric(
				Vec4(closest_surface_position.x(), closest_surface_position.y(), closest_surface_position.z(), 0),
				Vec4(n.x(), n.y(), n.z(), 1));

			// update sphere assignment
			(*p.samples_vertex_sphere_)[vid] = sphere;
			// update color
			(*p.samples_vertex_color_)[vid] = (*p.spheres_color_)[s_index];
			// compute shrinking ball
			auto [center, radius] = geometry::shrinking_ball_center(closest_surface_position, n, p.surface_kdt_);

			if (!is_inside(p, center))
			{
				std::cout << "Warning: shrinking ball center is outside the volume." << std::endl;
				std::cout << "Vertex: " << vid << ", pos: " << pos.x() << ", " << pos.y() << ", " << pos.z()
						  << "Projected Normal:" << n.transpose() << std::endl;
				Vec3 flipped_normal = -n;
				auto result = geometry::shrinking_ball_center(closest_surface_position, flipped_normal, p.surface_kdt_);
				center = result.first;
				radius = result.second;
				if (!is_inside(p, center))
				{
					std::cout << "Error: shrinking ball center is still outside the volume after flipping normal."
							  << std::endl;
				}
			}
			(*p.medial_axis_position_)[vid] = center;
			(*p.medial_axis_radius_)[vid] = radius;
			// // add to sphere cluster
			(*p.spheres_cluster_)[s_index].push_back(v);
			// (*p.spheres_cluster_area_)[s_index] += (*p.samples_weight_)[vid];
		};

		auto in_cluster = [&](Vec3& pos) {
			if (!is_inside(p, pos))
				return false;

			std::pair<uint32, Vec3> bvh_res;
			p.surface_bvh_->closest_point(pos, &bvh_res);
			Vec3 closest_surface_position = bvh_res.second;
			SFace closest_face = p.surface_bvh_faces_[bvh_res.first];
			Vec3 normal = (closest_surface_position - pos).normalized();
			auto spheres_neighbor_clusters_ = (*p.spheres_neighbor_clusters_)[s_index];
			Scalar min_energy = (std::numeric_limits<Scalar>::max)();
			Scalar best_index = s_index;

			std::vector<PVertex> candidate_spheres;
			candidate_spheres.reserve(spheres_neighbor_clusters_.size() + 1);
			for (uint32 cid : spheres_neighbor_clusters_)
			{
				auto pv = of_index<PVertex>(*p.spheres_, cid);
				if (pv != PVertex())
					candidate_spheres.push_back(pv);
			}
			candidate_spheres.push_back(sphere);
			for (PVertex neighbor_sphere : candidate_spheres)
			{
				uint32 neighbor_index = index_of(*p.spheres_, neighbor_sphere);
				Scalar r = (*p.spheres_radius_)[neighbor_index];
				Vec3 center = (*p.spheres_position_)[neighbor_index];
				Scalar dist_SQEM = (closest_surface_position - center).dot(normal) - r;
				Line_Quadric lq = Line_Quadric(pos, normal);
				Scalar dist_line_quadric = lq.eval(center);
				Scalar energy = dist_SQEM * dist_SQEM + p.sqem_update_lambda_ * dist_line_quadric;
				if (energy < min_energy)
				{
					min_energy = energy;
					best_index = neighbor_index;
				}
			}
			return best_index == s_index;
		};
		auto& spheres_clusters_ = (*p.spheres_cluster_)[s_index];
		auto& center = (*p.spheres_position_)[s_index];
		auto& radius = (*p.spheres_radius_)[s_index];
		uint32 added = p.volume_sampler_->sample_cluster(center, radius, spheres_clusters_, on_accept, in_volume,
														 in_cluster, target_count);
		// recompute kd_tree
		if (added < target_count)
		{
			std::cout << "Sample cluster failed: target: " << target_count << ", but only " << added
					  << " samples added." << std::endl;
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

			const Vec3& pos = (*p.projected_samples_position_)[v_index];
			const Vec3& n = (*p.projected_samples_normal_)[v_index];

			Spherical_Quadric& q = (*p.samples_quadric_)[v_index];
			q.clear();
			q = Spherical_Quadric(Vec4(pos.x(), pos.y(), pos.z(), 0), Vec4(n.x(), n.y(), n.z(), 1));

			Line_Quadric& lq = (*p.samples_line_quadric_)[v_index];
			lq.zero();
			lq = Line_Quadric(pos, n);

			return true;
		});
	}

	template <class MESH>
	acc::KDTree<3, uint32>* construct_kd_tree(
		MESH& mesh, std::vector<typename mesh_traits<MESH>::Vertex>& vertices,
		std::shared_ptr<typename mesh_traits<MESH>::template Attribute<Vec3>> pos_attribute)
	{
		vertices.clear();
		uint32 nb_vertices = nb_cells<typename mesh_traits<MESH>::Vertex>(mesh);
		vertices.reserve(nb_vertices);
		std::vector<Vec3> position_vector;
		position_vector.reserve(nb_vertices);
		foreach_cell(mesh, [&](typename mesh_traits<MESH>::Vertex v) {
			vertices.push_back(v);
			position_vector.push_back(value<Vec3>(mesh, pos_attribute, v));
			return true;
		});
		auto* tree = new acc::KDTree<3, uint32>(position_vector);
		return tree;
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
				boost::synapse::connect<typename MeshProvider<SURFACE>::attribute_changed_t<Vec3>>(
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
		p.vmas_delaunay_ = std::make_unique<VMAS_Delaunay>();
		if (!p.samples_)
			p.samples_ = points_provider_->add_mesh(surface_provider_->mesh_name(s) + "_samples");
		// Volumn Samples
		p.samples_position_ = get_or_add_attribute<Vec3, PVertex>(*p.samples_, "position");
		p.samples_vertex_color_ = get_or_add_attribute<Vec4, PVertex>(*p.samples_, "color");
		p.projected_samples_position_ = get_or_add_attribute<Vec3, PVertex>(*p.samples_, "projected_position");
		p.projected_samples_normal_ = get_or_add_attribute<Vec3, PVertex>(*p.samples_, "projected_normal");
		p.samples_poisson_depth_ = get_or_add_attribute<uint8, PVertex>(*p.samples_, "poisson_sample_depth");
		p.samples_weight_ = get_or_add_attribute<Scalar, PVertex>(*p.samples_, "weight");
		p.volume_sampler_ = std::make_unique<VolumePoissonDiskSampler>(*p.samples_);
		// initialize volumn samples
		poisson_disk_sampling(p);

		// create KDTree for the projected surface vertices
		if (p.surface_kdt_)
			delete p.surface_kdt_;
		if (p.samples_projections_kdt_)
			delete p.samples_projections_kdt_;
		p.surface_kdt_ = construct_kd_tree(*p.surface_, p.surface_kdt_vertices_, p.surface_vertex_position_);
		p.samples_projections_kdt_ =
			construct_kd_tree(*p.samples_, p.samples_projections_kdt_vertices_, p.projected_samples_position_);

		// compute knn graph
		p.samples_projections_vertex_knn_ = get_or_add_attribute<std::vector<PVertex>, PVertex>(*p.samples_, "knn");
		foreach_cell(*p.samples_, [&](PVertex v) -> bool {
			uint32 v_index = index_of(*p.samples_, v);
			const Vec3& pos = (*p.projected_samples_position_)[v_index];
			std::vector<std::pair<uint32, Scalar>> k_res;
			p.samples_projections_kdt_->find_nns(pos, k, &k_res);
			(*p.samples_projections_vertex_knn_)[v_index].clear();
			(*p.samples_projections_vertex_knn_)[v_index].reserve(k_res.size());
			for (const auto& [idx, dist] : k_res)
			{
				PVertex nv = p.samples_projections_kdt_vertices_[idx];
				if (nv != v)
					(*p.samples_projections_vertex_knn_)[v_index].push_back(nv);
			}
			return true;
		});

		foreach_cell(*p.samples_, [&](PVertex v) -> bool {
			uint32 v_index = index_of(*p.samples_, v);
			Scalar sum = 0.0;
			const Vec3& pos = (*p.projected_samples_position_)[v_index];
			auto s_knn = (*p.samples_projections_vertex_knn_)[v_index];
			for (PVertex nv : s_knn)
			{
				auto pos_nv = value<Vec3>(*p.samples_,p.projected_samples_position_, nv);
				Scalar dist = (pos_nv - pos).norm();
				sum += dist;
			}
			Scalar area = 1.0;
			(*p.samples_weight_)[v_index] = area;
			if (area == 0)
			{
				std::cout << "Warning: sample " << v_index << " has zero weight." << std::endl;
				std::cout << std::setprecision(10) << area << std::endl;
			}
			return true;
		});

		// initialize SQEM quadrics
		p.samples_line_quadric_ = get_or_add_attribute<Line_Quadric, PVertex>(*p.samples_, "line_quadric");
		p.samples_quadric_ = get_or_add_attribute<Spherical_Quadric, PVertex>(*p.samples_, "quadric");
		compute_quadrics(p);

		// compute shrinking balls for the surface vertices

		p.medial_axis_position_ = get_or_add_attribute<Vec3, PVertex>(*p.samples_, "medial_axis_position");
		p.medial_axis_radius_ = get_or_add_attribute<Scalar, PVertex>(*p.samples_, "medial_axis_radius");
		p.medial_axis_secondary_vertex_ =
			get_or_add_attribute<PVertex, PVertex>(*p.samples_, "medial_axis_secondary_vertex_");

		parallel_foreach_cell(*p.samples_, [&](PVertex v) -> bool {
			uint32 v_index = index_of(*p.samples_, v);
			auto [c, r] = geometry::shrinking_ball_center((*p.projected_samples_position_)[v_index],
														  (*p.projected_samples_normal_)[v_index], p.surface_kdt_);
			(*p.medial_axis_position_)[v_index] = c;
			(*p.medial_axis_radius_)[v_index] = r;

			return true;
		});
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
			get_or_add_attribute<std::set<uint32>, PVertex>(*p.spheres_, "neighbor_clusters"); // neighbor clusters

		// create the skeleton mesh

		if (!p.skeleton_)
			p.skeleton_ = non_manifold_provider_->add_mesh(surface_provider_->mesh_name(s) + "_skeleton");
		p.skeleton_position_ = get_or_add_attribute<Vec3, NMVertex>(*p.skeleton_, "position");

		// if we already have spheres, we need to recompute the clusters and errors
		// (cleans out surface vertex sphere data)

		compute_clusters(p);

		
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

		compute_clusters(p);
		
		if (!p.running_)
			update_render_data(p);
	}

	void compute_clusters(SurfaceParameters& p)
	{
		// clean cluster affectation
		parallel_foreach_cell(*p.spheres_, [&](PVertex v) -> bool {
			uint32 v_index = index_of(*p.spheres_, v);
			(*p.spheres_cluster_)[v_index].clear();
			(*p.spheres_cluster_area_)[v_index] = 0.0;
			return true;
		});
		p.samples_vertex_sphere_->fill(PVertex());

		if (p.nb_spheres_ == 0)
			return;

		// auto start = std::chrono::high_resolution_clock::now();

		parallel_foreach_cell(*p.samples_, [&](PVertex v) -> bool {
			uint32 v_index = index_of(*p.samples_, v);
			Scalar a = value<Scalar>(*p.samples_, p.samples_weight_, v);
			const Vec3& vp = (*p.projected_samples_position_)[v_index];
			const Vec3& pos = (*p.samples_position_)[v_index];
			Scalar min_distance = std::numeric_limits<Scalar>::max();
			PVertex closest_sphere;
			uint32 closest_sphere_index;
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
					dist_other = (vp - center).norm() - radius;
					dist_other *= dist_other;
				}
				break;
				case SPHERE_POWER_DISTANCE: {
					dist_other = (pos - center).dot(pos - center) - radius * radius;
				}
				break;
				case LINE_QUADRIC_DISTANCE: {
					dist_other = (*p.samples_line_quadric_)[v_index].eval(center);
				}
				break;
				}
				dist_other *= a;
				Scalar dist = dist_sqem + p.sqem_clustering_lambda_ * dist_other;
				if (dist < min_distance)
				{
					min_distance = dist;
					closest_sphere = pv;
					closest_sphere_index = pv_index;
				}
				return true;
			});

			value<PVertex>(*p.samples_, p.samples_vertex_sphere_, v) = closest_sphere;
			p.vmas_delaunay_->set_cluster_id(v_index, closest_sphere_index);
			std::lock_guard<std::mutex> lock(spheres_mutex_[closest_sphere_index % spheres_mutex_.size()]);
			value<std::vector<PVertex>>(*p.spheres_, p.spheres_cluster_, closest_sphere).push_back(v);
			//value<Scalar>(*p.spheres_, p.spheres_cluster_area_, closest_sphere) += a;

			return true;
		});
		// 	break;
		// }
		// }

		//visulaize each cluster 
		update_samples_color(p);
		MeshData<POINTS>& md = points_provider_->mesh_data(*p.samples_);
		uint32 last_nb_clusters = p.cluster_sets_.size();
		p.last_cluster_sets_.clear();
		p.last_cluster_sets_.resize(last_nb_clusters, nullptr);
		for (uint32 pv_index = 0; pv_index < p.nb_spheres_; pv_index++)
		{
			if (pv_index >= last_nb_clusters)
			{
				p.cluster_sets_.push_back(&md.template get_or_add_cells_set<PVertex>
										  ("cluster_" + std::to_string(pv_index)));
			}
			else
			{
				// assign previous cluster

				p.last_cluster_sets_[pv_index] =
					&md.template get_or_add_cells_set<PVertex>("last_cluster_" + std::to_string(pv_index));
				p.last_cluster_sets_[pv_index]->clear();
				p.cluster_sets_[pv_index]->foreach_cell([&](PVertex v) -> bool {
					p.last_cluster_sets_[pv_index]->select(v);
					return true;
				});
				p.cluster_sets_[pv_index]->clear();
			}
			const std::vector<PVertex>& cluster = (*p.spheres_cluster_)[pv_index];
			for (PVertex sv : cluster)
			{
				p.cluster_sets_[pv_index]->select(sv);
			}
		}

		//points_provider_->emit_attribute_changed(*p.samples_, p.samples_vertex_color_.get());

		for (auto& cluster_set : p.cluster_sets_)
		{
			points_provider_->emit_cells_set_changed(*p.samples_, cluster_set);
		}
		for (auto& cluster_set : p.last_cluster_sets_)
		{
			points_provider_->emit_cells_set_changed(*p.samples_, cluster_set);
		}
		

		// auto end = std::chrono::high_resolution_clock::now();

		// std::cout << "Cluster computation time: " << std::chrono::duration<Scalar>(end - start).count() << "s"
		// 		  << std::endl;

		// local sampling
		if (p.local_sampling_ == true)
		{

			bool changed = false;
			foreach_cell(*p.spheres_, [&](PVertex v) -> bool {
				std::vector<PVertex>& cluster = value<std::vector<PVertex>>(*p.spheres_, p.spheres_cluster_, v);
				auto cluster_size = cluster.size();
				if (cluster_size < 20)
				{
					poisson_disk_sampling_local(p, v, 20 - cluster_size);
					changed = true;
				}
				return true;
			});

			if (changed)
			{
				// recompute kd_tree
				if (p.samples_projections_kdt_)
					delete p.samples_projections_kdt_;
				p.samples_projections_kdt_ =
					construct_kd_tree(*p.samples_, p.samples_projections_kdt_vertices_, p.samples_position_);

				// recompute knn graph
				p.samples_projections_vertex_knn_ =
					get_or_add_attribute<std::vector<PVertex>, PVertex>(*p.samples_, "knn");
				foreach_cell(*p.samples_, [&](PVertex v) -> bool {
					uint32 v_index = index_of(*p.samples_, v);
					const Vec3& pos = (*p.projected_samples_position_)[v_index];
					std::vector<std::pair<uint32, Scalar>> k_res;
					p.samples_projections_kdt_->find_nns(pos, k, &k_res);
					(*p.samples_projections_vertex_knn_)[v_index].clear();
					(*p.samples_projections_vertex_knn_)[v_index].reserve(k_res.size());
					for (const auto& [idx, dist] : k_res)
					{
						PVertex nv = p.samples_projections_kdt_vertices_[idx];
						if (nv != v)
							(*p.samples_projections_vertex_knn_)[v_index].push_back(nv);
					}
					return true;
				});
				// recompute weights
				parallel_foreach_cell(*p.samples_, [&](PVertex v) -> bool {
					uint32 v_index = index_of(*p.samples_, v);
					Scalar sum = 0.0;
					const Vec3& pos = (*p.projected_samples_position_)[v_index];
					for (PVertex nv : (*p.samples_projections_vertex_knn_)[v_index])
						sum += (value<Vec3>(*p.samples_, p.projected_samples_position_, nv) - pos).norm();
					(*p.samples_weight_)[v_index] = 1.0;
					
					return true;
				});
			}
		}
		// update sphere clusters area
		parallel_foreach_cell(*p.spheres_, [&](PVertex v) -> bool {
			uint32 v_index = index_of(*p.spheres_, v);
			(*p.spheres_cluster_area_)[v_index] = 0.0;
			const std::vector<PVertex>& cluster = (*p.spheres_cluster_)[v_index];
			for (PVertex sv : cluster)
			{
				uint32 sv_index = index_of(*p.samples_, sv);
				Scalar a = value<Scalar>(*p.samples_, p.samples_weight_, sv);
				(*p.spheres_cluster_area_)[v_index] += a;
			}
			
			return true;
		});
		
	}

	// recompute clusters using power distance only (no SQEM term)
	void compute_clusters_power(SurfaceParameters& p)
	{
		parallel_foreach_cell(*p.spheres_, [&](PVertex v) -> bool {
			uint32 v_index = index_of(*p.spheres_, v);
			(*p.spheres_cluster_)[v_index].clear();
			(*p.spheres_cluster_area_)[v_index] = 0.0;
			return true;
		});
		p.samples_vertex_sphere_->fill(PVertex());

		if (p.nb_spheres_ == 0)
			return;

		parallel_foreach_cell(*p.samples_, [&](PVertex v) -> bool {
			uint32 v_index = index_of(*p.samples_, v);
			Scalar a = value<Scalar>(*p.samples_, p.samples_weight_, v);
			const Vec3& pos = (*p.samples_position_)[v_index];
			Scalar min_distance = std::numeric_limits<Scalar>::max();
			PVertex closest_sphere;
			uint32 closest_sphere_index = 0;
			foreach_cell(*p.spheres_, [&](PVertex pv) {
				uint32 pv_index = index_of(*p.spheres_, pv);
				const Vec3& center = (*p.spheres_position_)[pv_index];
				Scalar radius = (*p.spheres_radius_)[pv_index];
				Scalar dist_other = (pos - center).dot(pos - center) - radius * radius;
				dist_other *= a;
				if (dist_other < min_distance)
				{
					min_distance = dist_other;
					closest_sphere = pv;
					closest_sphere_index = pv_index;
				}
				return true;
			});

			value<PVertex>(*p.samples_, p.samples_vertex_sphere_, v) = closest_sphere;
			p.vmas_delaunay_->set_cluster_id(v_index, closest_sphere_index);
			std::lock_guard<std::mutex> lock(spheres_mutex_[closest_sphere_index % spheres_mutex_.size()]);
			value<std::vector<PVertex>>(*p.spheres_, p.spheres_cluster_, closest_sphere).push_back(v);
			value<Scalar>(*p.spheres_, p.spheres_cluster_area_, closest_sphere) += a;
			return true;
		});

		update_samples_color(p);
		MeshData<POINTS>& md = points_provider_->mesh_data(*p.samples_);
		uint32 last_nb_clusters = p.cluster_sets_.size();
		p.last_cluster_sets_.clear();
		p.last_cluster_sets_.resize(last_nb_clusters, nullptr);
		for (uint32 pv_index = 0; pv_index < p.nb_spheres_; pv_index++)
		{
			if (pv_index >= last_nb_clusters)
			{
				p.cluster_sets_.push_back(
					&md.template get_or_add_cells_set<PVertex>("cluster_" + std::to_string(pv_index)));
			}
			else
			{
				p.last_cluster_sets_[pv_index] =
					&md.template get_or_add_cells_set<PVertex>("last_cluster_" + std::to_string(pv_index));
				p.last_cluster_sets_[pv_index]->clear();
				p.cluster_sets_[pv_index]->foreach_cell([&](PVertex v) -> bool {
					p.last_cluster_sets_[pv_index]->select(v);
					return true;
				});
				p.cluster_sets_[pv_index]->clear();
			}
			const std::vector<PVertex>& cluster = (*p.spheres_cluster_)[pv_index];
			for (PVertex sv : cluster)
				p.cluster_sets_[pv_index]->select(sv);
		}

		points_provider_->emit_attribute_changed(*p.samples_, p.samples_vertex_color_.get());

		for (auto& cluster_set : p.cluster_sets_)
			points_provider_->emit_cells_set_changed(*p.samples_, cluster_set);
		for (auto& cluster_set : p.last_cluster_sets_)
			points_provider_->emit_cells_set_changed(*p.samples_, cluster_set);
	}

	void compute_spheres_error(SurfaceParameters& p)
	{
		parallel_foreach_cell(*p.spheres_, [&](PVertex v) {
			uint32 v_index = index_of(*p.spheres_, v);

			const Vec3& center = (*p.spheres_position_)[v_index];
			Scalar radius = (*p.spheres_radius_)[v_index];
			const std::vector<PVertex>& cluster = (*p.spheres_cluster_)[v_index];

			Scalar cluster_error = 0.0;
			for (PVertex sv : cluster)
			{
				uint32 sv_index = index_of(*p.samples_, sv);

				Scalar a = value<Scalar>(*p.samples_, p.samples_weight_, sv);
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
				case SPHERE_POWER_DISTANCE: {
					dist_other =
						((*p.samples_position_)[sv_index] - center).dot((*p.samples_position_)[sv_index] - center) -
						radius * radius;
				}
				break;
				case LINE_QUADRIC_DISTANCE: {
					dist_other = (*p.samples_line_quadric_)[sv_index].eval(center);
				}
				break;
				}
				dist_other *= a;
				Scalar dist = dist_sqem + p.sqem_update_lambda_ * dist_other;
				(*p.samples_vertex_error_)[sv_index] = dist;
				cluster_error += dist;
			}
			if ((*p.spheres_cluster_area_)[v_index] > 0.0)
			{
				(*p.spheres_error_)[v_index] = cluster_error / (*p.spheres_cluster_area_)[v_index];
				(*p.spheres_error_not_normalized_)[v_index] = cluster_error;
			}
			else
			{
				std::cerr << "Warning: sphere with zero cluster area." << std::endl;
				(*p.spheres_error_)[v_index] = 0.0;
				(*p.spheres_error_not_normalized_)[v_index] = 0.0;
			}

			return true;
		});

		p.min_error_ = (std::numeric_limits<Scalar>::max)();
		p.max_error_ = (std::numeric_limits<Scalar>::lowest)();
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

	void update_samples_color(SurfaceParameters& p)
	{
		parallel_foreach_cell(*p.samples_, [&](PVertex v) -> bool {
			uint32 v_index = index_of(*p.samples_, v);
			PVertex sphere = (*p.samples_vertex_sphere_)[v_index];
			if (sphere.is_valid())
			{
				uint32 sphere_index = index_of(*p.spheres_, sphere);
				(*p.samples_vertex_color_)[v_index] = (*p.spheres_color_)[sphere_index];
			}
			else
			{
				(*p.samples_vertex_color_)[v_index] = Vec4(1.0, 1.0, 1.0, 1.0);
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
		if (cluster.size() == 0)
		{
			std::cout << "Warning: empty cluster for sphere " << sphere_index << std::endl;
			return;
		}
		Vec3 c = (*p.spheres_position_)[sphere_index];
		Scalar r = (*p.spheres_radius_)[sphere_index];
		// Verify if the SQEM is well conditioned
		Spherical_Quadric q;
		for (PVertex v : cluster)
		{
			uint32 v_index = index_of(*p.samples_, v);
			Scalar a = (*p.samples_weight_)[v_index];
			q += (*p.samples_quadric_)[v_index] * a;
		}
		Scalar radius = 0.0;
		SQEM_CASE sc = q.well_conditioned(radius);
		if (sc != SQEM_CASE::Case4_Degenerate)
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
					const Vec3& n = (*p.projected_samples_normal_)[v_index];
					// SQEM energy
					Eigen::Vector4d lhs = Eigen::Vector4d::Zero();
					Scalar rhs = 0.0;
					Vec4 n4 = Vec4(n.x(), n.y(), n.z(), 1.0);
					Scalar a = (*p.samples_weight_)[v_index];
					const Scalar w_sqem = std::sqrt(std::max(a, Scalar(0)));
					lhs += -n4 * w_sqem;
					rhs += -1.0 * ((pos - Vec3(s(0), s(1), s(2))).dot(n) - s(3)) * w_sqem;
					J.row(idx) = lhs;
					b(idx) = rhs;
					++idx;

					// distance energy
					Vec3 d = pos - Vec3(s(0), s(1), s(2));
					Scalar l = d.norm();

					const Scalar w_dist = w_sqem * p.sqem_update_lambda_;
					if (l > Scalar(1e-12))
					{
						J.row(idx) = Eigen::Vector4d(-(d[0] / l), -(d[1] / l), -(d[2] / l), -1.0) * w_dist;
						b(idx) = -(l - s(3)) * w_dist; // scale the row by the update lambda
					}
					else
					{
						J.row(idx).setZero();
						b(idx) = 0.0;
					}

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

			auto [center, radius] =
				geometry::shrinking_ball_center(closest_point_position, closest_point_dir, p.surface_kdt_);
			c = center;
			r = radius;
		}
		if (!is_inside(p, c))
		{
			std::cout << "Sphere " << sphere_index << " center is outside after optimization, applying shrinking ball"
					  << std::endl;
			std::pair<uint32, Vec3> bvh_res;
			p.surface_bvh_->closest_point(c, &bvh_res);
			Vec3 closest_point_position = bvh_res.second;
			Vec3 closest_point_dir = (closest_point_position - c).normalized();

			closest_point_dir = -closest_point_dir;

			auto [center, radius] =
				geometry::shrinking_ball_center(closest_point_position, closest_point_dir, p.surface_kdt_);
			c = center;
			r = radius;
		}

		(*p.spheres_position_)[sphere_index] = c;
		(*p.spheres_radius_)[sphere_index] = r;
	}

	void update_sphere_center_distance(SurfaceParameters& p, PVertex sphere)
	{
		uint32 sphere_index = index_of(*p.spheres_, sphere);

		const std::vector<PVertex>& cluster = (*p.spheres_cluster_)[sphere_index];
		if (cluster.size() == 0)
		{
			std::cout << "Warning: empty cluster for sphere " << sphere_index << std::endl;
			return;
		}
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
			Scalar weight = value<Scalar>(*p.samples_, p.samples_weight_, v);
			if (weight <= 0.0)
				std::cout << "Warning: sample with zero weight in sphere " << sphere_index << std::endl;
			q += (*p.samples_quadric_)[v_index] * weight;
			h += weight * (*p.samples_position_)[v_index];
			area += weight;
		}
		Scalar r_sqem = 0.0;
		// if well conditioned, compute optimal sphere parameters
		// Solve SQEM for optimal r
		//  Fix r and solve for c
		SQEM_CASE sc = q.well_conditioned(r_sqem);
		if (sc == SQEM_CASE::Case2_Line || sc == SQEM_CASE::Case3_Plane)
		{

			Mat3 As = 2 * Mat3::Identity() * area;
			Vec3 bs = 2 * h;
			Mat3 A = q._A.block<3, 3>(0, 0) + p.sqem_update_lambda_ * As;
			Vec3 b = (q._b.head<3>() + p.sqem_update_lambda_ * bs) - q._A.block<3, 1>(0, 3) * r_sqem;
			c = A.ldlt().solve(b);

			r = r_sqem;
		}
		else if (sc == SQEM_CASE::Case1_Full)
		{
			Vec4 s;
			q.optimized(s);
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
			if ((*p.inside_tester_)(Point_3(c.x(), c.y(), c.z())) == CGAL::ON_UNBOUNDED_SIDE)
				closest_point_dir = -closest_point_dir;

			auto [center, radius] =
				geometry::shrinking_ball_center(closest_point_position, closest_point_dir, p.surface_kdt_);
			c = center;
			r = radius;
			std::cout << "Case 4 (Degenerate)" << " for sphere " << sphere_index << std::endl;
		}
		if (!is_inside(p, c))
		{
			std::pair<uint32, Vec3> bvh_res;
			p.surface_bvh_->closest_point(c, &bvh_res);
			Vec3 closest_point_position = bvh_res.second;
			Vec3 closest_point_dir = (closest_point_position - c).normalized();

			closest_point_dir = -closest_point_dir;

			auto [center, radius] =
				geometry::shrinking_ball_center(closest_point_position, closest_point_dir, p.surface_kdt_);
			c = center;
			r = radius;
			std::cout << "Sphere " << sphere_index << " center is outside after optimization, applying shrinking ball"
					  << std::endl;
		}
		/*std::cout << "Updated sphere " << sphere_index << " center to (" << c.x() << ", " << c.y() << ", " << c.z()
				  << "), radius to " << r << std::endl;*/
		(*p.spheres_position_)[sphere_index] = c;
		(*p.spheres_radius_)[sphere_index] = r;
	}

	void update_sphere_line_quadric_distance(SurfaceParameters& p, PVertex sphere)
	{
		uint32 sphere_index = index_of(*p.spheres_, sphere);

		const std::vector<PVertex>& cluster = (*p.spheres_cluster_)[sphere_index];
		if (cluster.size() == 0)
		{
			std::cout << "Warning: empty cluster for sphere " << sphere_index << std::endl;
			return;
		}
		Vec3 c = (*p.spheres_position_)[sphere_index];
		Scalar r = (*p.spheres_radius_)[sphere_index];
		// Verify if the SQEM is well conditioned
		Spherical_Quadric q;
		Line_Quadric lq;

		Scalar area = 0.0;
		Vec3 h;
		h.setZero();
		for (PVertex v : cluster)
		{
			uint32 v_index = index_of(*p.samples_, v);
			Scalar weight = value<Scalar>(*p.samples_, p.samples_weight_, v);
			if (weight <= 0.0)
				std::cout << "Warning: sample with zero volume weight in sphere " << sphere_index << std::endl;
			q += (*p.samples_quadric_)[v_index] * weight;
			h += weight * (*p.samples_position_)[v_index];
			lq += (*p.samples_line_quadric_)[v_index] * weight;

			area += weight;
		}
		Scalar r_sqem = 0.0;
		// if well conditioned, compute optimal sphere parameters
		// Solve SQEM for optimal r
		//  Fix r and solve for c
		SQEM_CASE sc = q.well_conditioned(r_sqem);
		if (sc != SQEM_CASE::Case4_Degenerate)
		{
			/*Mat4 As = q._A;
			Vec4 bs = q._b;
			Mat4 Al = lq.get_quadric().matrix();
			Mat4 A = As + p.sqem_update_lambda_ * Al;
			Vec4 s = A.ldlt().solve(bs);*/

			Mat4 A = q._A;
			Vec4 b = q._b;

			Mat4 Ql = lq.get_quadric().matrix();
			Mat3 Al = Ql.block<3, 3>(0, 0);
			Vec3 bl = -Ql.block<3, 1>(0, 3);
			if (p.line_quadric_combined_)
			{
				Mat4 Al_ext = Mat4::Zero();
				Al_ext.block<3, 3>(0, 0) = Al;
				Vec4 bl_ext = Vec4::Zero();
				bl_ext.head<3>() = bl;
				Mat4 A_c = A + p.sqem_update_lambda_ * Al_ext;
				Vec4 b_c = b + p.sqem_update_lambda_ * bl_ext;
				Vec4 s = A_c.ldlt().solve(b_c);
				c = s.head<3>();
				r = s[3];
			}
			else
			{
				Mat3 As = A.block<3, 3>(0, 0);
				Vec3 bs = b.head<3>();
				Vec3 Asr = q._A.block<3, 1>(0, 3);
				Mat3 A = As + p.sqem_update_lambda_ * Al;
				Vec3 b = (bs + p.sqem_update_lambda_ * bl) - Asr * r_sqem;
				c = A.ldlt().solve(b);

				r = r_sqem;
			}
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
			if ((*p.inside_tester_)(Point_3(c.x(), c.y(), c.z())) == CGAL::ON_UNBOUNDED_SIDE)
				closest_point_dir = -closest_point_dir;

			auto [center, radius] =
				geometry::shrinking_ball_center(closest_point_position, closest_point_dir, p.surface_kdt_);
			c = center;
			r = radius;
			std::cout << "Case 4 (Degenerate)" << " for sphere " << sphere_index << std::endl;
		}
		if (!is_inside(p, c))
		{
			std::pair<uint32, Vec3> bvh_res;
			p.surface_bvh_->closest_point(c, &bvh_res);
			Vec3 closest_point_position = bvh_res.second;
			Vec3 closest_point_dir = (closest_point_position - c).normalized();

			closest_point_dir = -closest_point_dir;

			auto [center, radius] =
				geometry::shrinking_ball_center(closest_point_position, closest_point_dir, p.surface_kdt_);
			c = center;
			r = radius;
			std::cout << "Sphere " << sphere_index << " center is outside after optimization, applying shrinking ball"
					  << std::endl;
		}
		/*std::cout << "Updated sphere " << sphere_index << " center to (" << c.x() << ", " << c.y() << ", " << c.z()
				  << "), radius to " << r << std::endl;*/
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

		std::pair<uint32, Vec3> bvh_res;
		p.surface_bvh_->closest_point(c, &bvh_res);
		closest_point_position = bvh_res.second;
		closest_point_dir = (closest_point_position - c).normalized();

		const Vec3& closest_face_normal =
			value<Vec3>(*p.surface_, p.surface_face_normal_, p.surface_bvh_faces_[bvh_res.first]);
		// TODO: exterior detection is not reliable
		if (closest_point_dir.dot(closest_face_normal) <= 0.0)
			closest_point_dir = -closest_point_dir;

		auto [center, radius] =
			geometry::shrinking_ball_center(closest_point_position, closest_point_dir, p.surface_kdt_);

		c = center;
		r = radius;
	}

	void update_spheres(SurfaceParameters& p)
	{
		// auto start = std::chrono::high_resolution_clock::now();

		compute_clusters(p);
		// std::cout << "Computed clusters." << std::endl;
		parallel_foreach_cell(*p.spheres_, [&](PVertex v) -> bool {
			switch (p.distance_mode_)
			{
			case SPHERE_EUCLIDEAN_DISTANCE:
				update_sphere_euclidean_distance(p, v);
				break;
			case SPHERE_POWER_DISTANCE:
				update_sphere_center_distance(p, v);
				break;
			case LINE_QUADRIC_DISTANCE:
				update_sphere_line_quadric_distance(p, v);
				break;
			}
			if (p.sphere_correction_ && p.sphere_correction_mode_ == CORRECT_ALWAYS)
				correct_sphere(p, v);
			value<bool>(*p.spheres_, p.spheres_do_not_split_, v) = false;
			return true;
		});
		// std::cout << "Updated spheres." << std::endl;
		//  	break;
		//  }
		//  }

		compute_spheres_error(p); // compute spheres error
		std::cout << "Iteration " << p.iteration_count_ << ": min error = " << p.min_error_
				  << ", max error = " << p.max_error_ << ", total error = " << p.total_error_
				  << ", total error diff = " << p.total_error_diff_ << std::endl;
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
							for (uint32 neighbor_id : (*p.spheres_neighbor_clusters_)[s_index])
								(*p.spheres_do_not_split_)[neighbor_id] = true;
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
							for (uint32 neighbor_id : (*p.spheres_neighbor_clusters_)[s_index])
								(*p.spheres_do_not_split_)[neighbor_id] = true;
							split_sphere(p, sphere);
							--to_split_max;
						}
					}
				}
			}
			break;
			}
		}
		// std::cout << "Auto-split spheres." << std::endl;

		// auto end = std::chrono::high_resolution_clock::now();
		// std::cout << "Update spheres: " << std::chrono::duration<Scalar>(end - start).count() << "s" << std::endl;

		if (!p.running_)
			update_render_data(p);
		std::cout << "-------------------------------------" << std::endl;
	}

	void split_sphere(SurfaceParameters& p, PVertex v)
	{
		uint32 v_index = index_of(*p.spheres_, v);

		// find the vertex of the cluster with maximal error w.r.t. the sphere

		PVertex max_error_vertex;
		Scalar max_error_vertex_error = (std::numeric_limits<Scalar>::lowest)();
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
		auto adjacent_map = p.vmas_delaunay_->compute_adjacent_clusters_main_cc();
		p.cluster_adjacency_ = adjacent_map;
		p.cluster_cc_ready_ = false;
		foreach_cell(*p.spheres_, [&](PVertex sv) -> bool {
			uint32 sv_index = index_of(*p.spheres_, sv);
			(*p.spheres_neighbor_clusters_)[sv_index] = adjacent_map[sv_index];
			return true;
		});
	}

	void compute_cluster_cc_info(SurfaceParameters& p)
	{
		if (!p.vmas_delaunay_)
			return;
		p.cluster_components_sizes_.clear();
		p.cluster_components_ = p.vmas_delaunay_->compute_clusters_components(&p.cluster_components_sizes_);
		p.cluster_pair_components_ = p.vmas_delaunay_->compute_adjacent_clusters_components();
		if (p.cluster_adjacency_.empty())
			p.cluster_adjacency_ = p.vmas_delaunay_->compute_adjacent_clusters();
		p.cluster_cc_ready_ = true;
	}
	void compute_skeleton(SurfaceParameters& p, bool compute_neighbor_clusters_only = false)
	{
		if (p.nb_spheres_ == 0)
			return;
		compute_cluster_neighbour(p);
		// auto start = std::chrono::high_resolution_clock::now();

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
			const std::set<uint32>& neighbors_ids = (*p.spheres_neighbor_clusters_)[pv_index];
			for (uint32 n_id : neighbors_ids)
			{
				NMVertex nmv2 = (*spheres_skeleton_vertex_map)[n_id];
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
			uint32 pv_index = index_of(*p.spheres_, pv);
			NMVertex nmv1 = value<NMVertex>(*p.spheres_, spheres_skeleton_vertex_map, pv);
			const std::set<uint32>& n_pv = (*p.spheres_neighbor_clusters_)[pv_index];
			for (const uint32& ne1 : n_pv)
			{
				NMVertex nmv2 = (*spheres_skeleton_vertex_map)[ne1];
				const std::set<uint32>& ne_ne1 = (*p.spheres_neighbor_clusters_)[ne1];
				for (const uint32& ne2 : ne_ne1)
				{
					if (n_pv.find(ne2) != n_pv.end())
					{
						NMVertex nmv3 = (*spheres_skeleton_vertex_map)[ne2];
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

			parallel_foreach_cell(*p.samples_, [&](PVertex v) -> bool {
				uint32 v_index = index_of(*p.samples_, v);
				PVertex sphere = (*p.samples_vertex_sphere_)[v_index];
				if (sphere.is_valid())
						(*p.samples_vertex_color_)[v_index] = value<Vec4>(*p.spheres_, p.spheres_color_, sphere);
				return true;
			});
			points_provider_->emit_attribute_changed(*p.samples_, p.samples_position_.get());
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

			parallel_foreach_cell(*p.samples_, [&](PVertex v) -> bool {
				uint32 v_index = index_of(*p.samples_, v);
				PVertex sphere = (*p.samples_vertex_sphere_)[v_index];
				if (sphere.is_valid())
						(*p.samples_vertex_color_)[v_index] = value<Vec4>(*p.spheres_, p.spheres_color_, sphere);
				return true;
			});
			points_provider_->emit_attribute_changed(*p.samples_, p.samples_position_.get());
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
				compute_clusters(p);
				compute_spheres_error(p);
				if (!p.running_)
					update_render_data(p);
			}
			else if (key_code == GLFW_KEY_D && picked_sphere_.is_valid())
			{
				std::lock_guard<std::mutex> lock(p.mutex_);
				remove_sphere(p, picked_sphere_);
				compute_clusters(p);
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
						compute_clusters(p);
						compute_spheres_error(p);
						update_render_data(p);
					}
				}
				if (ImGui::Button("Verify Conflict"))
				{
					verify_conflicts(p);
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
				ImGui::RadioButton("Sphere Eculidean", (int*)&p.distance_mode_, SPHERE_EUCLIDEAN_DISTANCE);
				ImGui::SameLine();
				ImGui::RadioButton("Sphere Power", (int*)&p.distance_mode_, SPHERE_POWER_DISTANCE);
				ImGui::SameLine();
				ImGui::RadioButton("Line Quadric", (int*)&p.distance_mode_, LINE_QUADRIC_DISTANCE);
				if (p.distance_mode_ == LINE_QUADRIC_DISTANCE)
					ImGui::Checkbox("Line quadric combined", &p.line_quadric_combined_);

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
				if (ImGui::Button("Compute power distance clusters"))
				{
					std::lock_guard<std::mutex> lock(p.mutex_);
					compute_clusters_power(p);
					compute_cluster_neighbour(p);
					compute_skeleton(p);
					non_manifold_provider_->emit_connectivity_changed(*p.skeleton_);
					compute_spheres_error(p);
					if (!p.running_)
						update_render_data(p);
				}
				if (ImGui::Button("Compute cluster CCs"))
				{
					std::lock_guard<std::mutex> lock(p.mutex_);
					compute_cluster_cc_info(p);
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
					compute_clusters(p);
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

					// Connected components info for picked sphere's cluster (precomputed)
					if (p.vmas_delaunay_ && p.cluster_cc_ready_)
					{
						uint32 picked_idx = index_of(*p.spheres_, picked_sphere_);
						uint32 cid = picked_idx;
						uint32 cc = 0;
						if (auto it = p.cluster_components_.find(cid); it != p.cluster_components_.end())
							cc = it->second;
						ImGui::Text("Cluster %u CCs: %u", cid, cc);
						if (auto it = p.cluster_components_sizes_.find(cid); it != p.cluster_components_sizes_.end())
						{
							ImGui::TextUnformatted("CC sizes:");
							for (uint32 sz : it->second)
								ImGui::SameLine(), ImGui::Text("%u", sz);
						}

						if (auto it = p.cluster_adjacency_.find(cid); it != p.cluster_adjacency_.end())
						{
							for (uint32 ncid : it->second)
							{
								auto key = cid < ncid ? std::make_pair(cid, ncid) : std::make_pair(ncid, cid);
								uint32 comps = 0;
								if (auto itp = p.cluster_pair_components_.find(key);
									itp != p.cluster_pair_components_.end())
									comps = itp->second;
								ImGui::Text("With cluster %u, CCs: %u", ncid, comps);
							}
						}
					}
					else if (p.vmas_delaunay_)
					{
						ImGui::TextUnformatted("Cluster CCs not computed (use Compute cluster CCs).");
					}
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
