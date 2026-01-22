#ifndef CGOGN_MODULE_UDF_TRAINING_H_
#define CGOGN_MODULE_UDF_TRAINING_H_

#include <cgogn/ui/app.h>
#include <cgogn/ui/module.h>

#include <cgogn/core/ui_modules/mesh_provider.h>
#include <cgogn/geometry/algos/area.h>
#include <cgogn/geometry/algos/fitting.h>
#include <cgogn/geometry/algos/length.h>
#include <cgogn/geometry/algos/medial_axis.h>
#include <cgogn/geometry/functions/angle.h>
#include <cgogn/geometry/functions/distance.h>
#include <cgogn/geometry/functions/bounding_box.h>
#include <cgogn/geometry/types/line_quadric.h>
#include <cgogn/geometry/types/quadric.h>
#include <cgogn/geometry/types/ray_level_set_sampler.h>
#include <cgogn/geometry/types/spherical_quadric.h>
#include <cgogn/geometry/types/spatial_grid.h>
#include <cgogn/geometry/types/neural_field_forward.h>
#include <cgogn/geometry/types/vector_traits.h>
#include <cgogn/geometry/types/fast_winding_number_traits.h>
#include <cgogn/geometry/types/fast_winding_number.h>


#include <cgogn/rendering/ui_modules/point_cloud_render.h>
#include <cgogn/rendering/ui_modules/surface_render.h>

#include <Eigen/Dense>
#include <Eigen/Eigenvalues>
#include <Eigen/Sparse>

#include <libacc/bvh_tree.h>
#include <libacc/bvh_tree_spheres.h>
#include <libacc/kd_tree.h>

// import CGAL
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Polygon_mesh_processing/distance.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/poisson_eliminate.h>

#include <GLFW/glfw3.h>
#include <algorithm>
#include <array>
#include <filesystem>
#include <mutex>
#include <numeric>
#include <random>
#include <set>
#include <thread>
#include <torch/autograd.h>
#include <torch/script.h>
#include <torch/torch.h>
#include <unordered_map>

namespace cgogn
{

namespace ui
{

using geometry::Line_Quadric;
using geometry::Mat3;
using geometry::Mat4;
using geometry::Scalar;
using geometry::Spherical_Quadric;
using geometry::Quadric;
using geometry::SQEM_CASE;
using geometry::SpatialGrid;
using geometry::BatchUDFResult;
using geometry::RayLevelSetSampler;
using geometry::NeuralFieldForward;
using geometry::Vec3;
using geometry::Vec4;

template <typename SURFACE, typename POINTS, typename NONMANIFOLD>
class UDFTraining : public ViewModule
{
	using PVertex = typename mesh_traits<POINTS>::Vertex;
	using NMVertex = typename mesh_traits<NONMANIFOLD>::Vertex;
	using SVertex = typename mesh_traits<SURFACE>::Vertex;

	using SFace = typename mesh_traits<SURFACE>::Face;
	using NMFace = typename mesh_traits<NONMANIFOLD>::Face;
	using NMEdge = typename mesh_traits<NONMANIFOLD>::Edge;

	using PointTraits = geometry::FWN_Point_Traits<POINTS>;
	template <int ORDER>
	using PointFWN = geometry::Fast_Winding_Number<PointTraits, ORDER>;
	template <typename T>
	using PAttribute = typename mesh_traits<POINTS>::template Attribute<T>;
	template <typename T>
	using NMAttribute = typename mesh_traits<NONMANIFOLD>::template Attribute<T>;
	using K = CGAL::Exact_predicates_inexact_constructions_kernel;
	using CGAL_Mesh = CGAL::Surface_mesh<K::Point_3>;
	using Vector_3 = typename K::Vector_3;
	using Point_3 = typename K::Point_3;
	using NMFaceKey = std::array<uint32, 3>;

	enum AutoSplitMode : uint32
	{
		ERROR_THRESHOLD,
		MAX_NB_SPHERES
	};

	enum CorrectionMode : uint32
	{
		CORRECT_ALWAYS,
		CORRECT_ON_SPLIT
	};
	enum DistanceMode : uint32
	{
		SPHERE_EUCLIDEAN_DISTANCE,
		LINE_QUADRIC_DISTANCE,
		PURE_EUCLIDEAN_DISTANCE
	};
	enum InputMode : uint32
	{
		INPUT_POINT_CLOUD,
		INPUT_SURFACE_MESH,
		INPUT_NEURAL_UDF
	};

private:
	struct PointsParameters;


	struct Tet;
	struct PointsParameters
	{
		bool initialized_ = false;
		bool fitting_data_computed_ = false;
		InputMode input_mode_ = INPUT_POINT_CLOUD;
		// Input Points
		POINTS* points_ = nullptr;
		std::shared_ptr<PAttribute<Vec3>> position_ = nullptr;
		std::shared_ptr<PAttribute<Vec3>> normal_ = nullptr;			  // Computed on input for sampling
		std::shared_ptr<PAttribute<std::vector<PVertex>>> knn_ = nullptr; // For input normals
		std::vector<Vec3> input_position_backup_;
		std::vector<Vec3> input_normal_backup_;
		bool input_jitter_backup_valid_ = false;
		float32 input_jitter_pos_pct_ = 0.5f;

		// Nueral UDF
		bool neural_udf_loaded_ = false;
		torch::jit::Module neural_udf_model_;
		std::string neural_udf_model_path_ = "";

		// Ray Sampling
		std::unique_ptr<RayLevelSetSampler> ray_sampler_;

		// Sampling & Fitting Data
		POINTS* samples_mesh_ = nullptr;
		std::shared_ptr<PAttribute<Vec3>> samples_position_ = nullptr;
		std::shared_ptr<PAttribute<Vec3>> samples_normal_ = nullptr;
		std::shared_ptr<PAttribute<Scalar>> samples_area_ = nullptr;
		std::shared_ptr<PAttribute<std::vector<PVertex>>> samples_knn_ = nullptr;
		std::shared_ptr<PAttribute<Spherical_Quadric>> samples_quadric_ = nullptr;
		std::shared_ptr<PAttribute<Line_Quadric>> samples_line_quadric_ = nullptr;

		// Medial Axis (on samples)
		std::shared_ptr<PAttribute<Vec3>> samples_ma_position_ = nullptr;
		std::shared_ptr<PAttribute<Scalar>> samples_ma_radius_ = nullptr;
		std::shared_ptr<PAttribute<PVertex>> samples_ma_secondary_vertex_ = nullptr;

		// Clustering (on samples)
		std::shared_ptr<PAttribute<PVertex>> samples_sphere_ = nullptr; // Cluster ID for each sample
		std::shared_ptr<PAttribute<Scalar>> samples_error_ = nullptr;
		std::shared_ptr<PAttribute<Vec4>> samples_color_ = nullptr;
		std::shared_ptr<PAttribute<Vec4>> samples_normal_color_ = nullptr;
		std::vector<Vec3> samples_position_backup_;
		std::vector<Vec3> samples_normal_backup_;
		bool samples_jitter_backup_valid_ = false;
		float32 samples_jitter_pos_pct_ = 0.5f;
		float32 samples_jitter_normal_sigma_ = 0.02f;

		std::unique_ptr<acc::BVHTreeSpheres<uint32, Vec3>> samples_wn_bvh_;
		std::vector<Vec3> samples_wn_bvh_centers_;
		std::vector<Scalar> samples_wn_bvh_radii_;
		std::unique_ptr<PointFWN<3>> samples_winding_number_;
		std::vector<PVertex> samples_wn_vertices_;
		Scalar beta_ = Scalar(2.0);

		acc::KDTree<3, uint32>* samples_kdtree_ = nullptr; // KDTree of alpha-expanding samples
		std::vector<PVertex> samples_kdtree_vertices_;	   // Vertices of alpha-expanding samples in KDTree order

		acc::KDTree<3, uint32>* input_kdtree_ = nullptr; // KDTree of input points
		std::vector<PVertex> input_kdtree_vertices_;	 // Vertices of input points in KDTree order

		std::unique_ptr<SpatialGrid> samples_spatial_grid_ = nullptr;
		// Spheres
		POINTS* spheres_ = nullptr;
		uint32 nb_spheres_ = 0;
		std::shared_ptr<PAttribute<Vec3>> spheres_position_ = nullptr;
		std::shared_ptr<PAttribute<Scalar>> spheres_radius_ = nullptr;
		std::shared_ptr<PAttribute<Vec4>> spheres_color_ = nullptr;
		std::shared_ptr<PAttribute<std::vector<PVertex>>> spheres_cluster_ = nullptr; // Sample points in cluster
		std::shared_ptr<PAttribute<Scalar>> spheres_cluster_area_ = nullptr;
		std::shared_ptr<PAttribute<Vec4>> spheres_cluster_color_ = nullptr;
		std::shared_ptr<PAttribute<std::set<PVertex>>> spheres_neighbor_clusters_ = nullptr;
		std::shared_ptr<PAttribute<PVertex>> spheres_parent_ = nullptr;
		std::shared_ptr<PAttribute<bool>> spheres_do_not_split_ = nullptr;
		std::shared_ptr<PAttribute<Scalar>> spheres_error_ = nullptr;
		std::shared_ptr<PAttribute<Scalar>> spheres_error_not_normalized_ = nullptr;

		// GPU cluster cache
		bool cluster_gpu_dirty_ = true;
		int64_t cluster_gpu_samples_ = 0;
		int64_t cluster_gpu_spheres_ = 0;
		torch::Device cluster_gpu_device_ = torch::kCPU;

		torch::Tensor samples_pos_gpu_;
		torch::Tensor samples_area_gpu_;
		torch::Tensor spheres_center_gpu_;
		torch::Tensor spheres_radius_gpu_;
		torch::Tensor samples_sqem_A_gpu_;
		torch::Tensor samples_sqem_b_gpu_;
		torch::Tensor samples_sqem_c_gpu_;

		torch::Tensor samples_line_Q_gpu_;

		std::vector<PVertex> samples_gpu_order_;
		std::vector<PVertex> spheres_gpu_order_;

		// Skeleton
		NONMANIFOLD* skeleton_ = nullptr;
		std::shared_ptr<NMAttribute<Vec3>> skeleton_position_ = nullptr;
		std::shared_ptr<NMAttribute<Scalar>> skeleton_radius_ = nullptr;
		std::shared_ptr<NMAttribute<std::set<std::size_t>>> incident_tets_ = nullptr;

		std::shared_ptr<NMAttribute<Vec3>> skeleton_face_color_ = nullptr;
		std::shared_ptr<NMAttribute<Vec3>> skeleton_edge_color_ = nullptr;
		std::shared_ptr<NMAttribute<uint32>> edge_degree_ = nullptr;

		std::map<NMFaceKey, NMFace> skeleton_faces_map_;

		std::unordered_map<std::size_t, Tet> skeleton_tets_;

		float32 filter_radius_threshold_ = 0.0f;
		float32 init_dilation_factor_ = 3.0f;

		bool sphere_correction_ = false;
		CorrectionMode sphere_correction_mode_ = CORRECT_ALWAYS;
		DistanceMode distance_mode_ = SPHERE_EUCLIDEAN_DISTANCE;
		bool auto_stop_ = false;
		bool auto_split_ = false;
		AutoSplitMode auto_split_mode_ = ERROR_THRESHOLD;
		float32 auto_split_error_threshold_ = 0.00025f;
		uint32 auto_split_max_nb_spheres_ = 50;
		bool error_as_spheres_color_ = false;
		float32 spheres_transparency_ = 0.5f;
		float32 sqem_update_lambda_ = 0.20f;
		float32 sqem_clustering_lambda_ = 0.20f;

		// Filtering
		float32 target_radius_ = 0.1f;
		float32 radius_tolerance_ = 0.01f;

		// Sampling Parameters
		float alpha_ = 0.005f;
		float sample_radius_ = 0.0025f;
		int sample_iterations_ = 30; // Max attempts per point
		int knn_k_ = 10;
		int seed_ = 42;
		int cluster_min_points_ = 20;
		float grid_cell_size_ = 0.0025f;
		// Neural UDF Sampling
		int num_alpha_samples_ = 200000;
		int batch_size_ = 65532;	   // sample batch
		float tol_ = 1e-5; // convergence tolerance
		bool preprocss_sample_points_ = false;

		// Neural UDF ray sampling parameters
		float udf_bbox_expand_ = 0.1f;
		float udf_lipschitz_ = 4.0f;
		float udf_delta_enter_ = 0.003f;
		int udf_max_iterations_ = 3000;

		// State
		Scalar total_error_ = 0.0;
		Scalar total_error_not_normalized_ = 0.0;
		Scalar last_total_error_ = 0.0;
		Scalar total_error_diff_ = 0.0;
		Scalar min_error_ = 0.0;
		Scalar max_error_ = 0.0;
		PVertex max_error_sphere_ = PVertex();

		// Threading
		uint32 iteration_count_ = 0;
		std::mutex mutex_;
		bool running_ = false;
		bool stopping_ = false;
		bool slow_down_ = true;
		uint32 update_rate_ = 20;

		~PointsParameters()
		{
			if (samples_kdtree_)
				delete samples_kdtree_;
			if (input_kdtree_)
				delete input_kdtree_;
		}
	};

	struct Tet
	{
		std::size_t tet_id;
		NMFace faces[4];

		void print_tet_info(PointsParameters& p) const
		{
			std::cout << "Tet ID: " << tet_id << "\n";
			for (int i = 0; i < 4; ++i)
			{
				// print incidented edge incident faces
				NMFace f = faces[i];
				if (!f.is_valid())
					continue;
				uint32 idf = index_of(*p.skeleton_, f);
				uint32 nb_tets = (*p.incident_tets_)[idf].size();
				auto in_edges = incident_edges(*p.skeleton_, f);
				std::cout << " Face " << i << " incident to " << nb_tets << " tets. \n";
				for (const auto& edge : in_edges)
				{
					auto in_faces = incident_faces(*p.skeleton_, edge);
					uint32 nb_faces = in_faces.size();
					std::cout << "  Incident Edge (" << index_of(*p.skeleton_, edge) << ") with" << nb_faces
							  << " faces ";
				}

				std::cout << "\n";
			}
		}
	};

public:
	UDFTraining(const App& app) : ViewModule(app, "UDFTraining")
	{
	}

	~UDFTraining()
	{
	}

	void set_selected_surface(SURFACE& s)
	{
		selected_surface_ = &s;
		surface_bvh_dirty_ = true;
	} // Compatibility

	void set_selected_points(POINTS& p)
	{
		selected_points_ = &p;
		init_points_data(p);
	}

	void load_neural_udf_model(POINTS& points, const std::string& model_path)
	{
		PointsParameters& p = points_parameters_[&points];
		if (!std::filesystem::exists(model_path))
		{
			std::cout << "Neural UDF model file does not exist: " << model_path << std::endl;
			return;
		}
		try
		{
			std::cout << "Loading Neural UDF model from: " << model_path << std::endl;
			p.neural_udf_model_ = torch::jit::load(model_path, device_);
			p.neural_udf_model_.eval();
			p.neural_udf_loaded_ = true;
			p.neural_udf_model_path_ = model_path;
			p.input_mode_ = INPUT_NEURAL_UDF;
			std::cout << "Loaded neural UDF model from: " << model_path << std::endl;
		}
		catch (const c10::Error& e)
		{
			std::cerr << "Error loading Neural UDF model: " << e.what() << std::endl;
			p.neural_udf_loaded_ = false;
		}
	}

	NeuralFieldForward make_neural_field_forward(PointsParameters& p)
	{
		return NeuralFieldForward(&p.neural_udf_model_, p.neural_udf_loaded_, device_);
	}
	void load_alpha_samples_to_mesh(PointsParameters& p, size_t num_points)
	{
		if (!p.neural_udf_loaded_)
		{
			std::cerr << "Neural UDF model not loaded. Cannot sample alpha level set." << std::endl;
			return;
		}

		p.samples_spatial_grid_ = std::make_unique<SpatialGrid>(p.grid_cell_size_);

		std::cout << "Sampling " << num_points << " points on alpha=" << p.alpha_ << " level set..." << std::endl;

		// Sample points on alpha level set
		// std::vector<Vec3> sampled_points = sample_alpha_level_set(p, num_points*10);
		RayLevelSetSampler::Params ray_params;
		ray_params.bbox_expand = p.udf_bbox_expand_;
		ray_params.alpha = p.alpha_;
		ray_params.tol = p.tol_;
		ray_params.step_bound = Scalar(2.0);
		ray_params.batch_size = p.batch_size_;
		ray_params.max_iterations = p.udf_max_iterations_;
		ray_params.max_outer_iterations = 500;
		ray_params.seed = p.seed_;

		auto update_ray_sampler = [&](const torch::Device& device) {
			if (!p.ray_sampler_)
				p.ray_sampler_ = std::make_unique<RayLevelSetSampler>(ray_params, device);
			else
			{
				p.ray_sampler_->set_params(ray_params);
				p.ray_sampler_->set_device(device);
			}
		};

		NeuralFieldForward udf = make_neural_field_forward(p);
		Vec3 bbox_min(0, 0, 0);
		Vec3 bbox_max(1, 1, 1);
		auto bbox_valid = [&](const Vec3& min, const Vec3& max) {
			return min[0] <= max[0] && min[1] <= max[1] && min[2] <= max[2];
		};
		bool has_bbox = false;
		if (p.points_ && p.position_)
		{
			const uint32 count = nb_cells<PVertex>(*p.points_);
			if (count > 0)
			{
				auto bb = cgogn::geometry::bounding_box(*p.position_.get());
				if (bbox_valid(bb.first, bb.second))
				{
					bbox_min = bb.first;
					bbox_max = bb.second;
					has_bbox = true;
				}
			}
		}
		if (!has_bbox && selected_surface_)
		{
			auto s_pos = get_attribute<Vec3, SVertex>(*selected_surface_, "position");
			if (s_pos)
			{
				auto bb = cgogn::geometry::bounding_box(*s_pos.get());
				if (bbox_valid(bb.first, bb.second))
				{
					bbox_min = bb.first;
					bbox_max = bb.second;
					has_bbox = true;
				}
			}
		}
		update_ray_sampler(udf.device());
		bool used_sdf_filter = false;
		std::vector<Vec3> sampled_points = p.ray_sampler_->sample_alpha_level_set_rays(
			udf, num_points, p.samples_spatial_grid_.get(), p.grid_cell_size_, bbox_min, bbox_max, &used_sdf_filter);
		if (p.preprocss_sample_points_ && !used_sdf_filter)
			pre_process_sampling_points(p, sampled_points);
		// sampled_points = poisson_eliminate_points(sampled_points, num_points);
		if (sampled_points.empty())
		{
			std::cerr << "Failed to sample points on alpha level set." << std::endl;
			return;
		}

		std::cout << "Successfully sampled " << sampled_points.size() << " points." << std::endl;

		// Clear existing samples
		if (p.samples_mesh_)
			points_provider_->clear_mesh(*p.samples_mesh_);
		p.samples_jitter_backup_valid_ = false;
		p.samples_position_backup_.clear();
		p.samples_normal_backup_.clear();

		// Add sampled points to samples_mesh_
		for (const Vec3& pt : sampled_points)
		{
			PVertex v = add_vertex(*p.samples_mesh_);
			uint32 v_idx = index_of(*p.samples_mesh_, v);
			// std::cout << pt.transpose() << std::endl;
			(*p.samples_position_)[v_idx] = pt;
		}

		// Compute normals using gradient of UDF
		std::cout << "Computing normals from UDF gradients..." << std::endl;
		std::vector<Vec3> all_positions;
		all_positions.reserve(sampled_points.size());
		foreach_cell(*p.samples_mesh_, [&](PVertex v) {
			uint32 v_idx = index_of(*p.samples_mesh_, v);
			all_positions.push_back((*p.samples_position_)[v_idx]);
			return true;
		});

		auto apply_grad_normals = [&](const auto& grad_result) {
			if (grad_result.ok && grad_result.gradients.size() == all_positions.size())
			{
				uint32 idx = 0;
				foreach_cell(*p.samples_mesh_, [&](PVertex v) {
					uint32 v_idx = index_of(*p.samples_mesh_, v);
					Vec3 normal = grad_result.gradients[idx].normalized();
					(*p.samples_normal_)[v_idx] = normal;

					// Compute normal color for visualization
					(*p.samples_normal_color_)[v_idx] =
						Vec4((normal.x() + 1.0) * 0.5, (normal.y() + 1.0) * 0.5, (normal.z() + 1.0) * 0.5, 1.0);
					idx++;
					return true;
				});
				return true;
			}
			return false;
		};

		BatchUDFResult grad_result = udf.forward_batch_with_grad(all_positions);
		bool normals_ok = apply_grad_normals(grad_result);
		if (!normals_ok)
			std::cerr << "Failed to compute normals from UDF gradients." << std::endl;

		std::cout << "Building KDTree for sampled points..." << std::endl;
		build_kdtree(p);

		points_provider_->emit_connectivity_changed(*p.samples_mesh_);

		std::cout << "Alpha level set sampling complete. Ready for fitting." << std::endl;
	}

	void pre_process_sampling_points(PointsParameters& p, std::vector<Vec3>& points)
	{
		const Scalar tol = Scalar(1e-2);
		const Scalar max_dist = p.alpha_ + tol;
		std::cout << "Before pre-processing, " << points.size() << " points." << std::endl;
		if (p.input_kdtree_)
		{
			auto new_end = std::remove_if(points.begin(), points.end(), [&](const Vec3& pos) {
				std::pair<uint32, Scalar> knn_res;
				return !p.input_kdtree_->find_nn(pos, &knn_res, max_dist);
			});
			points.erase(new_end, points.end());
			return;
		}
		//Todo: the bvh shuld not be built here, to be moved
		build_surface_bvh();
		if (!surface_bvh_)
			return;

		auto new_end = std::remove_if(points.begin(), points.end(), [&](const Vec3& pos) {
			std::pair<uint32, Vec3> cp;
			return !surface_bvh_->closest_point(pos, &cp, max_dist);
		});
		points.erase(new_end, points.end());
		std::cout << "After pre-processing, " << points.size() << " points." << std::endl;
	}

	std::vector<Vec3> poisson_eliminate_points(const std::vector<Vec3>& points, size_t target_num)
	{
		std::vector<Point_3> cgal_in;
		cgal_in.reserve(points.size());
		for (const Vec3& p : points)
		{
			cgal_in.push_back(Point_3(p.x(), p.y(), p.z()));
		}
		target_num = std::min(target_num, points.size());
		std::vector<Point_3> cgal_out;
		cgal_out.reserve(target_num);
		CGAL::poisson_eliminate(cgal_in, target_num, std::back_inserter(cgal_out));
		std::vector<Vec3> out;
		out.reserve(cgal_out.size());
		for (const Point_3& p : cgal_out)
		{
			out.push_back(Vec3(p.x(), p.y(), p.z()));
		}
		return out;
	}

	std::pair<std::vector<Vec3>, std::vector<Vec3>> project_points_to_alpha_gpu_impl(PointsParameters& p,
																					 const std::vector<Vec3>& points)
	{
		const size_t nb_points = points.size();
		const int max_iters = 30;
		const Scalar tol = 1e-6f;
		NeuralFieldForward udf = make_neural_field_forward(p);

		try
		{
			torch::Tensor X_cpu =
				torch::empty({static_cast<int64_t>(nb_points), 3},
							 torch::TensorOptions().dtype(torch::kFloat32).device(torch::kCPU).pinned_memory(true));
			torch::Tensor active_cpu = torch::ones({static_cast<int64_t>(nb_points)},
												   torch::TensorOptions().dtype(torch::kBool).device(device_));
			auto x_acc = X_cpu.accessor<float, 2>();
			auto active_acc = active_cpu.accessor<bool, 1>();
			for (size_t i = 0; i < nb_points; ++i)
			{
				x_acc[i][0] = static_cast<float>(points[i].x());
				x_acc[i][1] = static_cast<float>(points[i].y());
				x_acc[i][2] = static_cast<float>(points[i].z());
			}
			torch::Tensor X = X_cpu.to(device_);
			torch::Tensor active = active_cpu.to(device_);

			for (int iter = 0; iter < max_iters; ++iter)
			{
				auto [values, grad] = udf.forward_values_grad_gpu(X);
				if (!values.defined() || !grad.defined())
					break;

				torch::Tensor f = values;
				if (f.dim() == 2 && f.size(1) == 1)
					f = f.squeeze(1);

				torch::Tensor residual = f - p.alpha_;
				torch::Tensor abs_residual = torch::abs(residual);

				//// DEBUG
				// float max_residual = abs_residual.max().item<float>();
				// float mean_residual = abs_residual.mean().item<float>();
				// int active_count = torch::sum(active).item<int>();

				// std::cout << "Iter " << iter << ": active=" << active_count << ", max_res=" << max_residual
				//		  << ", mean_res=" << mean_residual << std::endl;
				//// END DEBUG
				// Decide which points are still active

				active = active & (abs_residual > tol);

				torch::Tensor grad_norm2 = (grad * grad).sum(1, true);
				torch::Tensor step = residual.unsqueeze(1) / (grad_norm2 + 1e-12f);
				step = torch::where(active.unsqueeze(1), step, torch::zeros_like(step));

				X = X - step * grad;

				X = torch::clamp(X, 0.0f, 1.0f);
			}

			auto [values, grad] = udf.forward_values_grad_gpu(X);
			if (!values.defined() || !grad.defined())
				return {};

			torch::Tensor final_X_cpu = X.to(torch::kCPU);
			torch::Tensor final_grad_cpu = grad.to(torch::kCPU);

			auto final_X_acc = final_X_cpu.accessor<float, 2>();
			auto final_grad_acc = final_grad_cpu.accessor<float, 2>();

			std::vector<Vec3> projected_points;
			std::vector<Vec3> normals;
			projected_points.reserve(nb_points);
			normals.reserve(nb_points);
			for (size_t i = 0; i < nb_points; ++i)
			{
				projected_points.push_back(Vec3(final_X_acc[i][0], final_X_acc[i][1], final_X_acc[i][2]));
				normals.push_back(Vec3(final_grad_acc[i][0], final_grad_acc[i][1], final_grad_acc[i][2]).normalized());
			}
			return {projected_points, normals};
		}
		catch (const c10::Error& e)
		{
			std::cerr << "GPU projection to alpha level set failed: " << e.what() << std::endl;
			return {};
		}
	}

	std::pair<std::vector<Vec3>, std::vector<Vec3>> project_points_to_alpha_gpu(PointsParameters& p,
																				const std::vector<Vec3>& points)
	{
		return project_points_to_alpha_gpu_impl(p, points);
	}

	void test_batch_forward(PointsParameters& p)
	{
		// generate batch of random test points in [0,1]^3
		const size_t N = 1000;
		std::vector<Vec3> test_points(N);
		std::mt19937 gen(42);
		std::uniform_real_distribution<Scalar> dis(0.0, 1.0);
		for (size_t i = 0; i < N; ++i)
		{
			test_points[i] = Vec3(dis(gen), dis(gen), dis(gen));
		}
		NeuralFieldForward udf = make_neural_field_forward(p);
		BatchUDFResult result = udf.forward_batch_with_grad(test_points);
		if (result.ok)
		{
			std::cout << "Batch UDF evaluation successful. Sample results:" << std::endl;
			for (size_t i = 0; i < 5; ++i)
			{
				std::cout << "Point: " << test_points[i].transpose() << " UDF: " << result.values[i]
						  << " Grad: " << result.gradients[i].transpose() << std::endl;
			}
		}
		else
		{
			std::cerr << "Batch UDF evaluation failed." << std::endl;
		}
	}

public:
	// --- Surface Sampling (CGAL) ---
	void sample_surface_to_points(SURFACE& surface, POINTS& points, int num_samples)
	{
		points_provider_->clear_mesh(points);

		// 1. Convert CGoGN SURFACE to CGAL::Surface_mesh
		CGAL_Mesh cgal_mesh;
		//
		// auto pos = cgogn::get_attribute<Vec3, SVertex>(surface, "position");

		// std::unordered_map<uint32, CGAL_Mesh::Vertex_index> v_map;

		//// Add vertices
		// foreach_cell(surface, [&](SVertex v) {
		//	uint32 v_idx = index_of(surface, v);
		//	const Vec3& p = (*pos)[v_idx];
		//	v_map[v_idx] = cgal_mesh.add_vertex(Point_3(p[0], p[1], p[2]));
		//	return true;
		// });

		//// Add faces
		// foreach_cell(surface, [&](SFace f) {
		//	std::vector<CGAL_Mesh::Vertex_index> face_v;
		//	foreach_incident_vertex(surface, f, [&](SVertex v) {
		//		face_v.push_back(v_map[index_of(surface, v)]);
		//		return true;
		//	});
		//
		//	cgal_mesh.add_face(face_v);
		//	return true;
		// });
		std::string filename = surface_provider_->mesh_filename(surface);
		if (!filename.empty())
		{
			if (!CGAL::IO::read_polygon_mesh(filename, cgal_mesh) || cgal_mesh.is_empty())
			{
				std::cout << "Error loading CGAL surface mesh from file: " << filename << std::endl;
			}
		}
		normalize_surface_mesh(cgal_mesh);

		// 2. Sample mesh
		std::vector<Point_3> sampled_points;

		// Using simple random sampling on mesh
		CGAL::Polygon_mesh_processing::sample_triangle_mesh(
			cgal_mesh, std::back_inserter(sampled_points),
			CGAL::parameters::number_of_points_per_area_unit(num_samples));

		std::cout << "Sampled " << sampled_points.size() << " points from surface." << std::endl;

		// 3. Store in POINTS mesh
		auto p_pos = get_or_add_attribute<Vec3, PVertex>(points, "position");
		auto p_norm = get_or_add_attribute<Vec3, PVertex>(
			points, "normal"); // Need to compute normals if sampler doesn't give them

		// For now, reconstruct normals using input mesh or sampler?
		// The basic sampler might not give normals directly in the point vector.
		// We can use a location map if provided, but for now let's just add points.
		// NOTE: Ideally we want normals too. PMP::sample_triangle_mesh can take a property map for output,
		// or we can estimate them later. For UDF training, input normals are important.

		// Let's iterate and add points. We will re-compute normals using the source surface or simple estimation.
		// Since we have the CGAL mesh, we can use AABB tree to get normals for sampled points?
		// Or assume dense enough and use PCA later?
		// "compute_input_normals" is called later in the pipeline usually?
		// modify compute_point_cloud_normals to work on this input?

		for (const auto& pt : sampled_points)
		{
			PVertex v = add_vertex(points);
			uint32 v_idx = index_of(points, v);
			(*p_pos)[v_idx] = Vec3(pt.x(), pt.y(), pt.z());
		}

		// Compute normals for the new input point cloud
		// Since we just sampled from a surface, we can use the surface normals directly if we had a location map.
		// Alternatively, just re-use the generic compute_pca_normal or rely on external processing.
		// For robustness, let's just ensure the attribute exists.
		// If the user wants precise surface normals, we'd need to use a different sampler overload.
		points_provider_->emit_connectivity_changed(points);
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
			if (selected_points_)
			{
				PointsParameters& p = points_parameters_[selected_points_];
				update_render_data(p);
			}
		});
		// Initialize PyTorch device
		if (torch::cuda::is_available())
		{
			std::cout << "CUDA is available! Using GPU device 0." << std::endl;
			device_ = torch::Device(torch::kCUDA, 0);
		}
		else
		{
			std::cout << "CUDA is not available! Using the CPU." << std::endl;
			device_ = torch::kCPU;
		}
	}

private:
	void build_surface_bvh()
	{
		if (!selected_surface_ || !surface_provider_)
			return;
		if (surface_bvh_ && !surface_bvh_dirty_)
			return;

		auto s_pos = get_attribute<Vec3, SVertex>(*selected_surface_, "position");
		if (!s_pos)
		{
			surface_bvh_.reset();
			surface_bvh_dirty_ = false;
			return;
		}

		MeshData<SURFACE>& md = surface_provider_->mesh_data(*selected_surface_);
		uint32 nb_vertices = md.template nb_cells<SVertex>();
		uint32 nb_faces = md.template nb_cells<SFace>();
		if (nb_vertices == 0 || nb_faces == 0)
		{
			surface_bvh_.reset();
			surface_bvh_dirty_ = false;
			return;
		}

		auto bvh_vertex_index = get_or_add_attribute<uint32, SVertex>(*selected_surface_, "__bvh_vertex_index");

		surface_bvh_vertices_.clear();
		surface_bvh_vertices_.reserve(nb_vertices);
		surface_bvh_vertex_positions_.clear();
		surface_bvh_vertex_positions_.reserve(nb_vertices);

		uint32 idx = 0;
		foreach_cell(*selected_surface_, [&](SVertex v) -> bool {
			surface_bvh_vertices_.push_back(v);
			value<uint32>(*selected_surface_, bvh_vertex_index, v) = idx++;
			surface_bvh_vertex_positions_.push_back(value<Vec3>(*selected_surface_, s_pos, v));
			return true;
		});

		surface_bvh_faces_.clear();
		surface_bvh_faces_.reserve(nb_faces);
		std::vector<uint32> face_vertex_indices;
		face_vertex_indices.reserve(nb_faces * 3);
		foreach_cell(*selected_surface_, [&](SFace f) -> bool {
			surface_bvh_faces_.push_back(f);
			foreach_incident_vertex(*selected_surface_, f, [&](SVertex v) -> bool {
				face_vertex_indices.push_back(value<uint32>(*selected_surface_, bvh_vertex_index, v));
				return true;
			});
			return true;
		});

		surface_bvh_ = std::make_unique<acc::BVHTree<uint32, Vec3>>(face_vertex_indices, surface_bvh_vertex_positions_);

		remove_attribute<SVertex>(*selected_surface_, bvh_vertex_index);
		surface_bvh_dirty_ = false;
	}

	// --- Initialization ---
	void normalize_surface_mesh(CGAL_Mesh& mesh)
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
	void init_points_data(POINTS& m)
	{

		PointsParameters& p = points_parameters_[&m];
		if (p.initialized_)
			return;

		// Init Input Points

		p.points_ = &m;
		p.position_ = get_attribute<Vec3, PVertex>(m, "position");
		p.normal_ = get_or_add_attribute<Vec3, PVertex>(m, "normal");
		p.knn_ = get_or_add_attribute<std::vector<PVertex>, PVertex>(m, "knn");
		// Build KDTree for input points
		uint32 nb_vertices = nb_cells<PVertex>(*p.points_);
		if (nb_vertices > 0)
		{
			std::vector<Vec3> points;
			points.reserve(nb_vertices);
			p.input_kdtree_vertices_.reserve(nb_vertices);
			foreach_cell(*p.points_, [&](PVertex v) {
				uint32 v_idx = index_of(*p.points_, v);
				points.push_back((*p.position_)[v_idx]);
				p.input_kdtree_vertices_.push_back(v);
				return true;
			});
			p.input_kdtree_ = new acc::KDTree<3, uint32>(points);
		}

		// Init Samples Mesh
		std::string sample_name = points_provider_->mesh_name(*p.points_) + "_samples";
		if (!p.samples_mesh_)
			p.samples_mesh_ = points_provider_->has_mesh(sample_name) ? points_provider_->mesh(sample_name)
																	  : points_provider_->add_mesh(sample_name);
		else
			points_provider_->clear_mesh(*p.samples_mesh_);

		// Initialize attributes for samples
		p.samples_position_ = get_or_add_attribute<Vec3, PVertex>(*p.samples_mesh_, "position");
		p.samples_normal_ = get_or_add_attribute<Vec3, PVertex>(*p.samples_mesh_, "normal");
		p.samples_area_ = get_or_add_attribute<Scalar, PVertex>(*p.samples_mesh_, "area");
		p.samples_knn_ = get_or_add_attribute<std::vector<PVertex>, PVertex>(*p.samples_mesh_, "knn");
		p.samples_quadric_ = get_or_add_attribute<Spherical_Quadric, PVertex>(*p.samples_mesh_, "quadric");
		p.samples_line_quadric_ = get_or_add_attribute<Line_Quadric, PVertex>(*p.samples_mesh_, "line_quadric");
		p.samples_ma_position_ = get_or_add_attribute<Vec3, PVertex>(*p.samples_mesh_, "ma_position");
		p.samples_ma_radius_ = get_or_add_attribute<Scalar, PVertex>(*p.samples_mesh_, "ma_radius");
		p.samples_ma_secondary_vertex_ =
			get_or_add_attribute<PVertex, PVertex>(*p.samples_mesh_, "ma_secondary_vertex");
		p.samples_sphere_ = get_or_add_attribute<PVertex, PVertex>(*p.samples_mesh_, "sphere");
		p.samples_error_ = get_or_add_attribute<Scalar, PVertex>(*p.samples_mesh_, "error");
		p.samples_color_ = get_or_add_attribute<Vec4, PVertex>(*p.samples_mesh_, "color");
		p.samples_normal_color_ = get_or_add_attribute<Vec4, PVertex>(*p.samples_mesh_, "normal_color");

		// Init Spheres Mesh
		std::string sphere_name = points_provider_->mesh_name(m) + "_spheres";
		if (!p.spheres_)
			p.spheres_ = points_provider_->has_mesh(sphere_name) ? points_provider_->mesh(sphere_name)
																 : points_provider_->add_mesh(sphere_name);

		p.spheres_position_ = get_or_add_attribute<Vec3, PVertex>(*p.spheres_, "position");
		p.spheres_radius_ = get_or_add_attribute<Scalar, PVertex>(*p.spheres_, "radius");
		p.spheres_color_ = get_or_add_attribute<Vec4, PVertex>(*p.spheres_, "color");
		p.spheres_cluster_ = get_or_add_attribute<std::vector<PVertex>, PVertex>(*p.spheres_, "cluster");
		p.spheres_cluster_color_ = get_or_add_attribute<Vec4, PVertex>(*p.spheres_, "cluster_color");
		p.spheres_cluster_area_ = get_or_add_attribute<Scalar, PVertex>(*p.spheres_, "cluster_area");
		p.spheres_neighbor_clusters_ =
			get_or_add_attribute<std::set<PVertex>, PVertex>(*p.spheres_, "neighbor_clusters");
		p.spheres_parent_ = get_or_add_attribute<PVertex, PVertex>(*p.spheres_, "parent");
		p.spheres_do_not_split_ = get_or_add_attribute<bool, PVertex>(*p.spheres_, "do_not_split");
		p.spheres_error_ = get_or_add_attribute<Scalar, PVertex>(*p.spheres_, "error");
		p.spheres_error_not_normalized_ = get_or_add_attribute<Scalar, PVertex>(*p.spheres_, "error_not_normalized");

		// Init Skeleton Mesh
		std::string skel_name = points_provider_->mesh_name(m) + "_skeleton";
		if (!p.skeleton_)
			p.skeleton_ = non_manifold_provider_->has_mesh(skel_name) ? non_manifold_provider_->mesh(skel_name)
																	  : non_manifold_provider_->add_mesh(skel_name);
		p.skeleton_position_ = get_or_add_attribute<Vec3, NMVertex>(*p.skeleton_, "position");
		p.skeleton_radius_ = get_or_add_attribute<Scalar, NMVertex>(*p.skeleton_, "radius");
		p.incident_tets_ = get_or_add_attribute<std::set<std::size_t>, NMFace>(*p.skeleton_, "incident_tets");
		p.edge_degree_ = get_or_add_attribute<uint32, NMEdge>(*p.skeleton_, "degree"); // edge face degree
		p.skeleton_face_color_ = get_or_add_attribute<Vec3, NMFace>(*p.skeleton_, "color");
		p.skeleton_edge_color_ = get_or_add_attribute<Vec3, NMEdge>(*p.skeleton_, "color");
		p.samples_spatial_grid_ = std::make_unique<SpatialGrid>(p.grid_cell_size_);

		p.initialized_ = true;
	}

	void compute_fitting_data(PointsParameters& p)
	{
		if (p.fitting_data_computed_)
			return;
		if (!p.samples_mesh_)
		{
			std::cerr << "Error: No sampled points found. Please sample points first." << std::endl;
			return;
		}

		std::cout << "Building KDTree..." << std::endl;
		build_kdtree(p);
		std::cout << "Computing KNN and Area..." << std::endl;
		compute_samples_area(p); // Compute KNN and Area for samples
		std::cout << "Computing Winding Numbers..." << std::endl;
		compute_winding_numbers(p);
		std::cout << "Computing Quadrics..." << std::endl;
		compute_quadrics(p);
		std::cout << "Computing Initial Medial Axis..." << std::endl;
		compute_initial_medial_axis(p);

		std::cout << "Fitting Data Computed." << std::endl;

		p.fitting_data_computed_ = true;
	}

	void build_kdtree(PointsParameters& p)
	{
		if (p.samples_kdtree_)
			delete p.samples_kdtree_;

		std::vector<Vec3> points;
		p.samples_kdtree_vertices_.clear();
		points.reserve(nb_cells<PVertex>(*p.samples_mesh_));
		p.samples_kdtree_vertices_.reserve(points.size());

		foreach_cell(*p.samples_mesh_, [&](PVertex v) {
			uint32 idx = index_of(*p.samples_mesh_, v);
			points.push_back((*p.samples_position_)[idx]);
			p.samples_kdtree_vertices_.push_back(v);
			return true;
		});

		p.samples_kdtree_ = new acc::KDTree<3, uint32>(points);
	}

	// Generic PCA normal computation for any point cloud mesh
	template <typename MESH, typename POS_ATTR>
	Vec3 compute_pca_normal(const MESH& mesh, const POS_ATTR& position, const std::vector<uint32>& indices,
							const std::vector<typename mesh_traits<MESH>::Vertex>& kdtree_vertices)
	{
		using Vertex = typename mesh_traits<MESH>::Vertex;
		if (indices.size() < 3)
			return Vec3(0, 0, 1);

		Vec3 centroid(0, 0, 0);
		for (uint32 idx : indices)
		{
			Vertex v = kdtree_vertices[idx];
			uint32 v_idx = index_of(mesh, v);
			centroid += position[v_idx];
		}
		centroid /= Scalar(indices.size());

		Eigen::Matrix3d covariance = Eigen::Matrix3d::Zero();
		for (uint32 idx : indices)
		{
			Vertex v = kdtree_vertices[idx];
			uint32 v_idx = index_of(mesh, v);
			Vec3 diff = position[v_idx] - centroid;
			Eigen::Vector3d pe(diff[0], diff[1], diff[2]);
			covariance += pe * pe.transpose();
		}
		covariance /= Scalar(indices.size());

		Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> solver(covariance);
		Eigen::Vector3d normal = solver.eigenvectors().col(0);
		return Vec3(normal[0], normal[1], normal[2]).normalized();
	}

	void compute_input_normals(PointsParameters& p)
	{
		// Compute normals
		parallel_foreach_cell(*p.points_, [&](PVertex v) {
			uint32 v_idx = index_of(*p.points_, v);
			const Vec3& pt = (*p.position_)[v_idx];
			std::vector<std::pair<uint32, Scalar>> knn_res;
			p.input_kdtree_->find_nns(pt, p.knn_k_, &knn_res);

			std::vector<uint32> indices;
			for (auto& res : knn_res)
				indices.push_back(res.first);

			(*p.normal_)[v_idx] = compute_pca_normal(*p.points_, *p.position_, indices, p.input_kdtree_vertices_);
			return true;
		});
	}

	Vec3 compute_avg_normal(PointsParameters& p, PVertex v, acc::KDTree<3, uint32>& kdtree,
							const std::vector<PVertex>& vertices)
	{
		uint32 v_idx = index_of(*p.points_, v);
		const Vec3& pt = (*p.position_)[v_idx];
		const Vec3& n = (*p.normal_)[v_idx];

		std::vector<std::pair<uint32, Scalar>> knn_res;
		kdtree.find_nns(pt, p.knn_k_, &knn_res);

		Vec3 avg_n(0, 0, 0);
		for (auto& res : knn_res)
		{
			PVertex neighbor = vertices[res.first];
			uint32 n_idx = index_of(*p.points_, neighbor);
			Vec3 nn = (*p.normal_)[n_idx];
			if (nn.dot(n) < 0)
				nn = -nn;
			avg_n += nn;
		}
		return avg_n.normalized();
	}
	void compute_samples_area(PointsParameters& p)
	{
		// Compute KNN and Area on samples_mesh_
		if (!p.samples_mesh_ || !p.samples_kdtree_)
			return;

		parallel_foreach_cell(*p.samples_mesh_, [&](PVertex v) {
			uint32 v_idx = index_of(*p.samples_mesh_, v);
			const Vec3& pt = (*p.samples_position_)[v_idx];
			std::vector<std::pair<uint32, Scalar>> knn_res;
			p.samples_kdtree_->find_nns(pt, p.knn_k_, &knn_res);

			(*p.samples_knn_)[v_idx].clear();
			Scalar sum_dist = 0.0;
			for (auto& res : knn_res)
			{
				if (p.samples_kdtree_vertices_[res.first] != v)
				{
					(*p.samples_knn_)[v_idx].push_back(p.samples_kdtree_vertices_[res.first]);
					sum_dist += res.second;
				}
			}
			// Normals are already computed/oriented in sample_points
			(*p.samples_area_)[v_idx] = (sum_dist * sum_dist) / (2.0 * p.knn_k_); // Rough area estimate
			return true;
		});
	}

	void compute_winding_numbers(PointsParameters& p)
	{
		uint32 nb_samples = nb_cells<PVertex>(*p.samples_mesh_);
		p.samples_wn_bvh_centers_.clear();
		p.samples_wn_bvh_radii_.clear();
		p.samples_wn_vertices_.clear();
		p.samples_wn_bvh_centers_.reserve(nb_samples);
		p.samples_wn_bvh_radii_.reserve(nb_samples);
		p.samples_wn_vertices_.reserve(nb_samples);

		foreach_cell(*p.samples_mesh_, [&](PVertex v) {
			uint32 v_idx = index_of(*p.samples_mesh_, v);
			p.samples_wn_vertices_.push_back(v);
			p.samples_wn_bvh_centers_.push_back((*p.samples_position_)[v_idx]);
			Scalar area = (*p.samples_area_)[v_idx];
			Scalar radius = std::sqrt(area / M_PI);
			p.samples_wn_bvh_radii_.push_back(radius);
			return true;
		});

		p.samples_wn_bvh_ = std::make_unique<acc::BVHTreeSpheres<uint32, Vec3>>(p.samples_wn_bvh_centers_, p.samples_wn_bvh_radii_);

		PointTraits samples_wn_traits(*p.samples_mesh_, p.samples_wn_bvh_.get(), p.samples_wn_vertices_,
								 p.samples_position_.get(), p.samples_normal_.get(),
								 p.samples_area_.get());
		p.samples_winding_number_ = std::make_unique<PointFWN<3>>(*p.samples_wn_bvh_, samples_wn_traits, p.beta_);
	}

	void compute_quadrics(PointsParameters& p)
	{
		parallel_foreach_cell(*p.samples_mesh_, [&](PVertex v) {
			uint32 v_idx = index_of(*p.samples_mesh_, v);
			Spherical_Quadric& q = (*p.samples_quadric_)[v_idx];
			Line_Quadric& lq = (*p.samples_line_quadric_)[v_idx];
			q.clear();
			lq.clear();
			const Vec3& pos = (*p.samples_position_)[v_idx];
			const Vec3& n = (*p.samples_normal_)[v_idx];
			Scalar a = (*p.samples_area_)[v_idx] / (p.knn_k_ + 1.0);
			q += Spherical_Quadric(Vec4(pos.x(), pos.y(), pos.z(), 0), Vec4(n.x(), n.y(), n.z(), 1)) * a;
			lq += Line_Quadric(pos, n) * a;
			for (PVertex vn : (*p.samples_knn_)[v_idx])
			{
				uint32 vn_idx = index_of(*p.samples_mesh_, vn);
				const Vec3& pn = (*p.samples_position_)[vn_idx];
				const Vec3& nn = (*p.samples_normal_)[vn_idx];
				Scalar an = (*p.samples_area_)[vn_idx] / (p.knn_k_ + 1.0);
				q += Spherical_Quadric(Vec4(pn.x(), pn.y(), pn.z(), 0), Vec4(nn.x(), nn.y(), nn.z(), 1)) * an;
				lq += Line_Quadric(pn, nn) * an;
			}
			return true;
		});
	}

	// --- Sampling ---

	Vec3 random_sample_around(const Vec3& p, const Scalar radius, std::uniform_real_distribution<Scalar>& uni,
							  std::mt19937& rng)
	{
		Scalar u = uni(rng);
		Scalar v = uni(rng);
		Scalar w = uni(rng);

		Scalar R3 = radius * radius * radius;
		Scalar r = std::cbrt(R3 + u * (8 * R3 - R3)); // r = pow((R^3 + u(8R^3 - R^3)), 1/3)

		Scalar phi = v * 2.0 * M_PI;
		Scalar theta = std::acos(1.0 - 2.0 * w);

		Scalar x = r * std::sin(theta) * std::cos(phi);
		Scalar y = r * std::sin(theta) * std::sin(phi);
		Scalar z = r * std::cos(theta);

		return p + Vec3(x, y, z);
	}

	void backup_samples_state(PointsParameters& p)
	{
		const uint32 count = nb_cells<PVertex>(*p.samples_mesh_);
		p.samples_position_backup_.assign(count, Vec3(0, 0, 0));
		p.samples_normal_backup_.assign(count, Vec3(0, 0, 1));
		foreach_cell(*p.samples_mesh_, [&](PVertex v) {
			uint32 v_idx = index_of(*p.samples_mesh_, v);
			if (v_idx < count)
			{
				p.samples_position_backup_[v_idx] = (*p.samples_position_)[v_idx];
				p.samples_normal_backup_[v_idx] = (*p.samples_normal_)[v_idx];
			}
			return true;
		});
		p.samples_jitter_backup_valid_ = true;
	}

	void refresh_sample_normals_color(PointsParameters& p)
	{
		parallel_foreach_cell(*p.samples_mesh_, [&](PVertex v) -> bool {
			uint32 v_idx = index_of(*p.samples_mesh_, v);
			const Vec3& n = (*p.samples_normal_)[v_idx];
			(*p.samples_normal_color_)[v_idx] =
				Vec4((n.x() + 1.0) * 0.5, (n.y() + 1.0) * 0.5, (n.z() + 1.0) * 0.5, 1.0);
			return true;
		});
	}

	void rebuild_input_kdtree(PointsParameters& p)
	{
		if (p.input_kdtree_)
			delete p.input_kdtree_;
		p.input_kdtree_ = nullptr;
		p.input_kdtree_vertices_.clear();

		const uint32 nb_vertices = nb_cells<PVertex>(*p.points_);
		if (nb_vertices == 0)
			return;

		std::vector<Vec3> points;
		points.reserve(nb_vertices);
		p.input_kdtree_vertices_.reserve(nb_vertices);
		foreach_cell(*p.points_, [&](PVertex v) {
			uint32 v_idx = index_of(*p.points_, v);
			points.push_back((*p.position_)[v_idx]);
			p.input_kdtree_vertices_.push_back(v);
			return true;
		});
		p.input_kdtree_ = new acc::KDTree<3, uint32>(points);
	}

	void invalidate_samples_after_input_change(PointsParameters& p)
	{
		p.fitting_data_computed_ = false;
		p.samples_winding_number_.reset();
		p.samples_wn_bvh_.reset();
		p.cluster_gpu_dirty_ = true;
		p.samples_jitter_backup_valid_ = false;
		p.samples_position_backup_.clear();
		p.samples_normal_backup_.clear();
		if (p.samples_mesh_)
			points_provider_->clear_mesh(*p.samples_mesh_);
	}

	void backup_input_state(PointsParameters& p)
	{
		const uint32 count = nb_cells<PVertex>(*p.points_);
		p.input_position_backup_.assign(count, Vec3(0, 0, 0));
		p.input_normal_backup_.assign(count, Vec3(0, 0, 1));
		foreach_cell(*p.points_, [&](PVertex v) {
			uint32 v_idx = index_of(*p.points_, v);
			if (v_idx < count)
			{
				p.input_position_backup_[v_idx] = (*p.position_)[v_idx];
				p.input_normal_backup_[v_idx] = (*p.normal_)[v_idx];
			}
			return true;
		});
		p.input_jitter_backup_valid_ = true;
	}

	void apply_input_position_jitter(PointsParameters& p)
	{
		if (!p.points_ || !p.position_)
			return;
		const uint32 count = nb_cells<PVertex>(*p.points_);
		if (count == 0)
			return;
		auto [bb_min, bb_max] = cgogn::geometry::bounding_box(*p.position_.get());
		const Scalar diag = (bb_max - bb_min).norm();
		const Scalar pct = static_cast<Scalar>(p.input_jitter_pos_pct_) * Scalar(0.01);
		if (diag <= Scalar(0) || pct <= Scalar(0))
			return;
		if (!p.input_jitter_backup_valid_)
			backup_input_state(p);

		const Scalar sigma = diag * pct;
		std::mt19937 gen(p.seed_ + 4242);
		std::normal_distribution<Scalar> normal(Scalar(0), sigma);
		foreach_cell(*p.points_, [&](PVertex v) {
			uint32 v_idx = index_of(*p.points_, v);
			Vec3 pos = (*p.position_)[v_idx];
			pos += Vec3(normal(gen), normal(gen), normal(gen));
			(*p.position_)[v_idx] = pos;
			return true;
		});

		rebuild_input_kdtree(p);
		compute_input_normals(p);
		points_provider_->emit_attribute_changed(*p.points_, p.position_.get());
		points_provider_->emit_attribute_changed(*p.points_, p.normal_.get());
		invalidate_samples_after_input_change(p);
	}

	void restore_input_state(PointsParameters& p)
	{
		if (!p.points_ || !p.input_jitter_backup_valid_)
			return;
		const uint32 count = nb_cells<PVertex>(*p.points_);
		if (p.input_position_backup_.size() < count || p.input_normal_backup_.size() < count)
		{
			p.input_jitter_backup_valid_ = false;
			return;
		}
		foreach_cell(*p.points_, [&](PVertex v) {
			uint32 v_idx = index_of(*p.points_, v);
			(*p.position_)[v_idx] = p.input_position_backup_[v_idx];
			(*p.normal_)[v_idx] = p.input_normal_backup_[v_idx];
			return true;
		});

		rebuild_input_kdtree(p);
		points_provider_->emit_attribute_changed(*p.points_, p.position_.get());
		points_provider_->emit_attribute_changed(*p.points_, p.normal_.get());
		invalidate_samples_after_input_change(p);
	}

	void apply_samples_position_jitter(PointsParameters& p)
	{
		if (!p.samples_mesh_)
			return;
		const uint32 count = nb_cells<PVertex>(*p.samples_mesh_);
		if (count == 0)
			return;
		auto [bb_min, bb_max] = cgogn::geometry::bounding_box(*p.samples_position_.get());
		const Scalar diag = (bb_max - bb_min).norm();
		const Scalar pct = static_cast<Scalar>(p.samples_jitter_pos_pct_) * Scalar(0.01);
		if (diag <= Scalar(0) || pct <= Scalar(0))
			return;
		if (!p.samples_jitter_backup_valid_)
			backup_samples_state(p);

		const Scalar sigma = diag * pct;
		std::mt19937 gen(p.seed_ + 1337);
		std::normal_distribution<Scalar> normal(Scalar(0), sigma);
		foreach_cell(*p.samples_mesh_, [&](PVertex v) {
			uint32 v_idx = index_of(*p.samples_mesh_, v);
			Vec3 pos = (*p.samples_position_)[v_idx];
			pos += Vec3(normal(gen), normal(gen), normal(gen));
			(*p.samples_position_)[v_idx] = pos;
			return true;
		});

		p.fitting_data_computed_ = false;
		p.samples_winding_number_.reset();
		p.samples_wn_bvh_.reset();
		p.cluster_gpu_dirty_ = true;
		points_provider_->emit_attribute_changed(*p.samples_mesh_, p.samples_position_.get());
	}

	void apply_samples_normal_jitter(PointsParameters& p)
	{
		if (!p.samples_mesh_)
			return;
		const uint32 count = nb_cells<PVertex>(*p.samples_mesh_);
		if (count == 0)
			return;
		const Scalar sigma = static_cast<Scalar>(p.samples_jitter_normal_sigma_);
		if (sigma <= Scalar(0))
			return;
		if (!p.samples_jitter_backup_valid_)
			backup_samples_state(p);

		std::mt19937 gen(p.seed_ + 1337);
		std::normal_distribution<Scalar> normal(Scalar(0), sigma);
		foreach_cell(*p.samples_mesh_, [&](PVertex v) {
			uint32 v_idx = index_of(*p.samples_mesh_, v);
			Vec3 n = (*p.samples_normal_)[v_idx];
			n += Vec3(normal(gen), normal(gen), normal(gen));
			if (n.squaredNorm() > Scalar(0))
				n.normalize();
			else
				n = Vec3(0, 0, 1);
			(*p.samples_normal_)[v_idx] = n;
			return true;
		});

		refresh_sample_normals_color(p);
		p.fitting_data_computed_ = false;
		p.samples_winding_number_.reset();
		p.samples_wn_bvh_.reset();
		p.cluster_gpu_dirty_ = true;
		points_provider_->emit_attribute_changed(*p.samples_mesh_, p.samples_normal_.get());
		points_provider_->emit_attribute_changed(*p.samples_mesh_, p.samples_normal_color_.get());
	}

	void restore_samples_state(PointsParameters& p)
	{
		if (!p.samples_mesh_ || !p.samples_jitter_backup_valid_)
			return;
		const uint32 count = nb_cells<PVertex>(*p.samples_mesh_);
		if (p.samples_position_backup_.size() < count || p.samples_normal_backup_.size() < count)
		{
			p.samples_jitter_backup_valid_ = false;
			return;
		}
		foreach_cell(*p.samples_mesh_, [&](PVertex v) {
			uint32 v_idx = index_of(*p.samples_mesh_, v);
			(*p.samples_position_)[v_idx] = p.samples_position_backup_[v_idx];
			(*p.samples_normal_)[v_idx] = p.samples_normal_backup_[v_idx];
			return true;
		});

		refresh_sample_normals_color(p);
		p.fitting_data_computed_ = false;
		p.samples_winding_number_.reset();
		p.samples_wn_bvh_.reset();
		p.cluster_gpu_dirty_ = true;
		points_provider_->emit_attribute_changed(*p.samples_mesh_, p.samples_position_.get());
		points_provider_->emit_attribute_changed(*p.samples_mesh_, p.samples_normal_.get());
		points_provider_->emit_attribute_changed(*p.samples_mesh_, p.samples_normal_color_.get());
	}

	std::pair<Vec3, Vec3> project_and_normal(PointsParameters& p, const Vec3& sample_pos)
	{
		std::pair<uint32, Scalar> knn_res;
		p.input_kdtree_->find_nn(sample_pos, &knn_res);
		PVertex nn = p.input_kdtree_vertices_[knn_res.first];
		Vec3 nn_pos = (*p.position_)[index_of(*p.points_, nn)];
		Vec3 n = sample_pos - nn_pos;
		const Scalar eps = Scalar(1e-12);
		if (n.squaredNorm() < eps)
		{
			// Fallback to PCA normal when the sample coincides with its nearest neighbor.
			std::vector<std::pair<uint32, Scalar>> knn_res_all;
			p.input_kdtree_->find_nns(sample_pos, p.knn_k_, &knn_res_all);
			std::vector<uint32> indices;
			indices.reserve(knn_res_all.size());
			for (const auto& res : knn_res_all)
				indices.push_back(res.first);
			n = compute_pca_normal(*p.points_, *p.position_, indices, p.input_kdtree_vertices_);
		}
		if (n.squaredNorm() < eps)
			n = Vec3(0, 0, 1);
		else
			n.normalize();
		Vec3 query = nn_pos + n * p.alpha_;
		PVertex last_nn = PVertex();
		do
		{
			p.input_kdtree_->find_nn(query, &knn_res);
			nn = p.input_kdtree_vertices_[knn_res.first];
			nn_pos = (*p.position_)[index_of(*p.points_, nn)];
			Vec3 step_dir = query - nn_pos;
			if (step_dir.squaredNorm() < eps)
				step_dir = n;
			step_dir.normalize();
			query = nn_pos + step_dir * p.alpha_;
			last_nn = nn;
		} while (index_of(*p.points_, nn) != index_of(*p.points_, last_nn));

		return {query, n};
	}

	void sample_points(PointsParameters& p)
	{

		// Resampling from a point cloud invalidates cached fitting data and samples.
		p.fitting_data_computed_ = false;
		p.samples_winding_number_.reset();
		p.samples_wn_bvh_.reset();
		p.cluster_gpu_dirty_ = true;
		p.samples_jitter_backup_valid_ = false;
		p.samples_position_backup_.clear();
		p.samples_normal_backup_.clear();
		if (p.samples_mesh_)
			points_provider_->clear_mesh(*p.samples_mesh_);

		std::uniform_real_distribution<Scalar> uniform(0.0, 1.0);
		std::mt19937 gen(p.seed_);

		// Pick a random point from input as generator seed
		std::uniform_int_distribution<uint32> uniform_idx(0, uint32(p.input_kdtree_vertices_.size() - 1));
		uint32 rand_start_idx = uniform_idx(gen);
		PVertex start_seed_vertex = p.input_kdtree_vertices_[rand_start_idx];
		Vec3 generator = (*p.position_)[index_of(*p.points_, start_seed_vertex)];
		// Jitter the initial generator to avoid zero-length projection direction.
		Scalar jitter_radius =Scalar(0.001);
		if (jitter_radius > Scalar(0))
		{
			Vec3 jitter(uniform(gen) * Scalar(2.0) - Scalar(1.0),
						uniform(gen) * Scalar(2.0) - Scalar(1.0),
						uniform(gen) * Scalar(2.0) - Scalar(1.0));
			if (jitter.squaredNorm() > Scalar(0))
			{
				jitter.normalize();
				generator += jitter * jitter_radius;
			}
		}

		// --- Spatial Grid Optimization ---
		SpatialGrid grid(p.sample_radius_);

		PVertex start_vertex = add_vertex(*p.samples_mesh_);
		uint32 start_idx = index_of(*p.samples_mesh_, start_vertex);
		auto [pos, normal] = project_and_normal(p, generator);
		(*p.samples_position_)[start_idx] = pos;
		(*p.samples_normal_)[start_idx] = normal;

		grid.insert(pos, start_idx);

		std::vector<PVertex> active_list = {start_vertex};

		uint32 count = 1;

		while (!active_list.empty())
		{
			int rand_index = rand() % active_list.size();
			PVertex current_vertex = active_list[rand_index];
			uint32 current_idx = index_of(*p.samples_mesh_, current_vertex);
			Vec3 current_pos = (*p.samples_position_)[current_idx];
			bool found_new_sample = false;
			for (uint32 i = 0; i < p.sample_iterations_ && !found_new_sample; i++)
			{
				Vec3 sample_pos = random_sample_around(current_pos, p.sample_radius_, uniform, gen);
				auto [pos, normal] = project_and_normal(p, sample_pos);

				// Check using spatial grid (Poisson disk constraint)
				if (grid.is_valid_sample(pos, p.sample_radius_, *p.samples_position_))
				{
					// Verify distance to input cloud >= epsilon
					std::pair<uint32, Scalar> input_nn_res;
					p.input_kdtree_->find_nn(pos, &input_nn_res);
					if (input_nn_res.second < p.alpha_)
						continue; // Reject: too close to input surface

					PVertex new_vertex = add_vertex(*p.samples_mesh_);
					uint32 new_idx = index_of(*p.samples_mesh_, new_vertex);

					(*p.samples_position_)[new_idx] = pos;
					(*p.samples_normal_)[new_idx] = normal;

					grid.insert(pos, new_idx);

					active_list.push_back(new_vertex);
					count++;
					found_new_sample = true;
					if (count % 100 == 0)
						std::cout << "Sampled point " << count << "\r" << std::flush;
				}
			}
			if (!found_new_sample)
			{
				active_list[rand_index] = active_list.back();
				active_list.pop_back();
			}
		}

		std::cout << "Building Final KDTree..." << std::endl;
		build_kdtree(p);

		// Compute normal color
		parallel_foreach_cell(*p.samples_mesh_, [&](PVertex v) -> bool {
			uint32 v_idx = index_of(*p.samples_mesh_, v);
			const Vec3& n = (*p.samples_normal_)[v_idx];
			(*p.samples_normal_color_)[v_idx] =
				Vec4((n.x() + 1.0) * 0.5, (n.y() + 1.0) * 0.5, (n.z() + 1.0) * 0.5, 1.0);
			return true;
		});

		points_provider_->emit_connectivity_changed(*p.samples_mesh_);
	}

	// --- Shrinking Balls ---

	void compute_initial_medial_axis(PointsParameters& p)
	{
		uint32 total = nb_cells<PVertex>(*p.samples_mesh_);

		parallel_foreach_cell(*p.samples_mesh_, [&](PVertex v) {
			uint32 v_idx = index_of(*p.samples_mesh_, v);
			const Vec3& pt = (*p.samples_position_)[v_idx];
			const Vec3& n = (*p.samples_normal_)[v_idx];

			/*auto [c1, r1, q1] = cgogn::geometry::shrinking_ball_center<PVertex>(
				pt, n, p.samples_kdtree_, p.samples_kdtree_vertices_, p.alpha_ * 1.5
			);*/

			auto c = pt - n * p.alpha_;
			Scalar r = p.alpha_;
			auto q = c - n * r;
			std::pair<uint32, Scalar> knn_res;
			p.samples_kdtree_->find_nn(q, &knn_res);
			PVertex q1 = p.samples_kdtree_vertices_[knn_res.first];

			(*p.samples_ma_position_)[v_idx] = c;
			(*p.samples_ma_radius_)[v_idx] = r;
			(*p.samples_ma_secondary_vertex_)[v_idx] = *reinterpret_cast<PVertex*>(&q1);

			return true;
		});
	}

	void init_spheres(PointsParameters& p, uint32 max_nb_spheres)
	{
		points_provider_->clear_mesh(*p.spheres_);

		std::vector<PVertex> sorted_vertices;
		foreach_cell(*p.samples_mesh_, [&](PVertex v) {
			sorted_vertices.push_back(v);
			return true;
		});
		std::sort(sorted_vertices.begin(), sorted_vertices.end(), [&](PVertex a, PVertex b) {
			// sort candidate spheres by decreasing radius
			uint32 idx_a = index_of(*p.samples_mesh_, a);
			uint32 idx_b = index_of(*p.samples_mesh_, b);
			return (*p.samples_ma_radius_)[idx_a] > (*p.samples_ma_radius_)[idx_b];
		});

		auto covered = get_or_add_attribute<bool, PVertex>(*p.samples_mesh_, "__covered");
		covered->fill(false);

		p.nb_spheres_ = 0;

		for (PVertex v : sorted_vertices)
		{
			uint32 v_index = index_of(*p.samples_mesh_, v);

			if (p.nb_spheres_ >= max_nb_spheres)
				break;

			if ((*covered)[v_index])
				continue;

			const Vec3& vp = (*p.samples_ma_position_)[v_index];
			Scalar vr = (*p.samples_ma_radius_)[v_index];

			PVertex sphere = add_vertex(*p.spheres_);
			p.nb_spheres_++;
			uint32 sphere_index = index_of(*p.spheres_, sphere);

			(*p.spheres_position_)[sphere_index] = vp;
			(*p.spheres_radius_)[sphere_index] = vr;
			(*p.spheres_cluster_color_)[sphere_index] =
				Vec4(0.5 + 0.5 * (rand() % 256) / 256.0, 0.5 + 0.5 * (rand() % 256) / 256.0,
					 0.5 + 0.5 * (rand() % 256) / 256.0, 1.0);

			std::vector<PVertex> stack;
			stack.push_back(v);
			while (!stack.empty())
			{
				PVertex w = stack.back();
				stack.pop_back();
				uint32 w_idx = index_of(*p.samples_mesh_, w);
				(*covered)[w_idx] = true;

				// Use KNN for propagation on point cloud
				for (PVertex u : (*p.samples_knn_)[w_idx])
				{
					uint32 u_idx = index_of(*p.samples_mesh_, u);
					if (!(*covered)[u_idx] &&
						((*p.samples_position_)[u_idx] - vp).norm() < p.init_dilation_factor_ * vr)
						stack.push_back(u);
				}
			}

			// Also check secondary vertex logic if needed, but for now stick to KNN propagation
			PVertex secondary = (*p.samples_ma_secondary_vertex_)[v_index];
			if (secondary.is_valid())
			{
				stack.push_back(secondary);
				while (!stack.empty())
				{
					PVertex w = stack.back();
					stack.pop_back();
					uint32 w_idx = index_of(*p.samples_mesh_, w);
					(*covered)[w_idx] = true;

					for (PVertex u : (*p.samples_knn_)[w_idx])
					{
						uint32 u_idx = index_of(*p.samples_mesh_, u);
						if (!(*covered)[u_idx] &&
							((*p.samples_position_)[u_idx] - vp).norm() < p.init_dilation_factor_ * vr)
							stack.push_back(u);
					}
				}
			}
		}

		remove_attribute<PVertex>(*p.samples_mesh_, covered);

		compute_clusters(p);

		if (!p.running_)
			update_render_data(p);
	}

	void build_cluster_gpu_cache(PointsParameters& p)
	{
		const int64_t N = static_cast<int64_t>(nb_cells<PVertex>(*p.samples_mesh_));
		const int64_t M = static_cast<int64_t>(nb_cells<PVertex>(*p.spheres_));
		if (N == 0 || M == 0)
			return;

		const bool need_rebuild = p.cluster_gpu_dirty_ || !p.samples_pos_gpu_.defined() ||
								  !p.spheres_center_gpu_.defined() || p.cluster_gpu_samples_ != N ||
								  p.cluster_gpu_spheres_ != M || p.cluster_gpu_device_ != device_;

		if (!need_rebuild)
			return;

		p.cluster_gpu_samples_ = N;
		p.cluster_gpu_spheres_ = M;
		p.cluster_gpu_device_ = device_;
		p.cluster_gpu_dirty_ = false;

		auto cpu_opts = torch::TensorOptions().dtype(torch::kFloat32).device(torch::kCPU);
		if (device_.is_cuda())
			cpu_opts = cpu_opts.pinned_memory(true);

		// samples
		p.samples_gpu_order_.clear();
		p.samples_gpu_order_.reserve(N);

		torch::Tensor samples_pos_cpu = torch::empty({N, 3}, cpu_opts);
		torch::Tensor samples_area_cpu = torch::empty({N}, cpu_opts);
		torch::Tensor A_cpu = torch::empty({N, 4, 4}, cpu_opts);
		torch::Tensor b_cpu = torch::empty({N, 4}, cpu_opts);
		torch::Tensor c_cpu = torch::empty({N}, cpu_opts);
		torch::Tensor Q_cpu = torch::empty({N, 4, 4}, cpu_opts);

		auto pos_acc = samples_pos_cpu.accessor<float, 2>();
		auto area_acc = samples_area_cpu.accessor<float, 1>();
		auto A_acc = A_cpu.accessor<float, 3>();
		auto b_acc = b_cpu.accessor<float, 2>();
		auto c_acc = c_cpu.accessor<float, 1>();
		auto Q_acc = Q_cpu.accessor<float, 3>();

		int64_t i = 0;
		foreach_cell(*p.samples_mesh_, [&](PVertex v) {
			uint32 v_idx = index_of(*p.samples_mesh_, v);
			const Vec3& pt = (*p.samples_position_)[v_idx];

			// copy position
			pos_acc[i][0] = static_cast<float>(pt.x());
			pos_acc[i][1] = static_cast<float>(pt.y());
			pos_acc[i][2] = static_cast<float>(pt.z());
			// copy area
			area_acc[i] = static_cast<float>((*p.samples_area_)[v_idx]);
			// copy quadric

			const Spherical_Quadric& sq = (*p.samples_quadric_)[v_idx];
			const Quadric& lq = (*p.samples_line_quadric_)[v_idx].get_quadric();

			// Spherical quadric
			for (int r = 0; r < 4; ++r)
			{
				for (int c = 0; c < 4; ++c)
					A_acc[i][r][c] = static_cast<float>(sq._A(r, c));
				b_acc[i][r] = static_cast<float>(sq._b(r));
			}

			// Line quadric
			Mat4 q_mat = lq.matrix();
			for (int r = 0; r < 4; ++r)
				for (int c = 0; c < 4; ++c)
					Q_acc[i][r][c] = static_cast<float>(q_mat(r, c));

			p.samples_gpu_order_.push_back(v);
			++i;
			return true;
		});

		// spheres
		p.spheres_gpu_order_.clear();
		p.spheres_gpu_order_.reserve(M);

		torch::Tensor spheres_center_cpu = torch::empty({M, 3}, cpu_opts);
		torch::Tensor spheres_radius_cpu = torch::empty({M}, cpu_opts);

		auto cen_acc = spheres_center_cpu.accessor<float, 2>();
		auto rad_acc = spheres_radius_cpu.accessor<float, 1>();

		int64_t s = 0;
		foreach_cell(*p.spheres_, [&](PVertex v) {
			uint32 v_idx = index_of(*p.spheres_, v);
			const Vec3& c = (*p.spheres_position_)[v_idx];
			Scalar r = (*p.spheres_radius_)[v_idx];

			cen_acc[s][0] = static_cast<float>(c.x());
			cen_acc[s][1] = static_cast<float>(c.y());
			cen_acc[s][2] = static_cast<float>(c.z());
			rad_acc[s] = static_cast<float>(r);

			p.spheres_gpu_order_.push_back(v);
			++s;
			return true;
		});

		const bool nonblocking = device_.is_cuda();
		p.samples_pos_gpu_ = samples_pos_cpu.to(device_, nonblocking);
		p.samples_area_gpu_ = samples_area_cpu.to(device_, nonblocking);
		p.spheres_center_gpu_ = spheres_center_cpu.to(device_, nonblocking);
		p.spheres_radius_gpu_ = spheres_radius_cpu.to(device_, nonblocking);
		p.samples_sqem_A_gpu_ = A_cpu.to(device_, nonblocking);
		p.samples_sqem_b_gpu_ = b_cpu.to(device_, nonblocking);
		p.samples_sqem_c_gpu_ = c_cpu.to(device_, nonblocking);
		p.samples_line_Q_gpu_ = Q_cpu.to(device_, nonblocking);
	}

	torch::Tensor eval_sqem_batch(const torch::Tensor& A, // [N,4,4]
								  const torch::Tensor& b, // [N,4]
								  const torch::Tensor& c, // [N]
								  const torch::Tensor& v  // [B,3]
	)
	{
		torch::Tensor v_exp = v.unsqueeze(0).expand({A.size(0), v.size(0), 4}); // [N,B,4]
		torch::Tensor A_exp = A.unsqueeze(1);									// [N,1,4,4]

		torch::Tensor Av = torch::matmul(A_exp, v_exp.unsqueeze(-1)).squeeze(-1); // [N,B,4]

		torch::Tensor vTAv = (v_exp * Av).sum(-1);									 // [N,B]
		torch::Tensor b_exp = b.unsqueeze(1);										 // [N,1,4]
		torch::Tensor dist = 0.5f * vTAv - (b_exp * v_exp).sum(-1) + c.unsqueeze(1); // [N,B]
		return dist;
	}

	torch::Tensor eval_line_quadric_batch(const torch::Tensor& Q, // [N,4,4]
										  const torch::Tensor& v) // [B,4], v = [cx,cy,cz,1]
	{
		torch::Tensor v_exp = v.unsqueeze(0).expand({Q.size(0), v.size(0), 4});
		torch::Tensor Q_exp = Q.unsqueeze(1);
		//						 	  [N,1,4,4]	[N,B,4,1] -> [N,B,4,1] -> [N,B,4] 		
		torch::Tensor Qv = torch::matmul(Q_exp, v_exp.unsqueeze(-1)).squeeze(-1);
		torch::Tensor vTQv = (v_exp * Qv).sum(-1); // [N, B, 4] *[N,B, 4] -> [N,B]

		return vTQv;
	}

	void compute_clusters(PointsParameters& p)
	{
		// clean cluster affectation
		parallel_foreach_cell(*p.spheres_, [&](PVertex v) -> bool {
			uint32 v_index = index_of(*p.spheres_, v);
			(*p.spheres_cluster_)[v_index].clear();
			(*p.spheres_cluster_area_)[v_index] = 0.0;
			return true;
		});
		p.samples_sphere_->fill(PVertex());

		if (p.nb_spheres_ == 0)
			return;

		parallel_foreach_cell(*p.samples_mesh_, [&](PVertex v) -> bool {
			uint32 v_index = index_of(*p.samples_mesh_, v);
			// if (!(*p.medial_axis_selected_)[v_index])
			// 	return true;

			Scalar a = (*p.samples_area_)[v_index];

			const Vec3& vp = (*p.samples_position_)[v_index];
			Scalar min_distance = std::numeric_limits<Scalar>::max();
			PVertex closest_sphere;
			uint32 closest_sphere_index;

			foreach_cell(*p.spheres_, [&](PVertex pv) {
				uint32 pv_index = index_of(*p.spheres_, pv);

				const Vec3& center = (*p.spheres_position_)[pv_index];
				Scalar radius = (*p.spheres_radius_)[pv_index];

				Scalar dist = 0.0;
				Scalar dist_other = 0.0;
				switch (p.distance_mode_)
				{
				case SPHERE_EUCLIDEAN_DISTANCE: {
					Scalar dist_sqem =
						(*p.samples_quadric_)[v_index].eval(Vec4(center.x(), center.y(), center.z(), radius));
					dist_other = ((vp - center).norm() - radius);
					dist_other *= dist_other;
					dist_other *= a;
					dist = dist_sqem + p.sqem_clustering_lambda_ * dist_other;
				}
				break;

				case LINE_QUADRIC_DISTANCE: {
					Scalar dist_sqem =
						(*p.samples_quadric_)[v_index].eval(Vec4(center.x(), center.y(), center.z(), radius));
					dist_other = (*p.samples_line_quadric_)[v_index].eval(center);
					dist = dist_sqem + p.sqem_clustering_lambda_ * dist_other;
				}
				break;

				case PURE_EUCLIDEAN_DISTANCE: {
					dist_other = ((vp - center).norm() - radius);
					dist_other *= dist_other;
					dist = dist_other * a;
				}
				break;
				}
				if (dist < min_distance)
				{
					min_distance = dist;
					closest_sphere = pv;
					closest_sphere_index = pv_index;
				}
				return true;
			});

			(*p.samples_sphere_)[v_index] = closest_sphere;

			std::lock_guard<std::mutex> lock(spheres_mutex_[closest_sphere_index % spheres_mutex_.size()]);
			(*p.spheres_cluster_)[closest_sphere_index].push_back(v);
			(*p.spheres_cluster_area_)[closest_sphere_index] += a;

			return true;
		});
		// augment_insufficient_clusters(p);
	}

	void compute_clusters_gpu(PointsParameters& p)
	{
		if (!device_.is_cuda())
			return;

		if (!p.samples_mesh_ || !p.spheres_)
			return;

		const int64_t N = static_cast<int64_t>(nb_cells<PVertex>(*p.samples_mesh_));
		const int64_t M = static_cast<int64_t>(nb_cells<PVertex>(*p.spheres_));
		if (N == 0 || M == 0)
			return;

		build_cluster_gpu_cache(p);

		auto samples_pos = p.samples_pos_gpu_;	 // [N,3]
		auto samples_area = p.samples_area_gpu_; // [N]
		auto centers = p.spheres_center_gpu_;	 // [M,3]
		auto radius = p.spheres_radius_gpu_;	 // [M]

		auto fopts = torch::TensorOptions().dtype(torch::kFloat32).device(device_);
		auto iopts = torch::TensorOptions().dtype(torch::kInt64).device(device_);

		torch::Tensor min_dist = torch::full({N}, std::numeric_limits<float>::infinity(), fopts);
		torch::Tensor argmin = torch::full({N}, -1, iopts);

		const int64_t chunk = 128;
		for (int64_t k = 0; k < M; k += chunk)
		{
			const int64_t B = std::min(chunk, M - k);
			auto centers_c = centers.narrow(0, k, B); // [B,3]
			auto radius_c = radius.narrow(0, k, B);	  // [B]

			torch::Tensor dist;

			if(p.distance_mode_ == PURE_EUCLIDEAN_DISTANCE){
				
				torch::Tensor diff = samples_pos.unsqueeze(1) - centers_c.unsqueeze(0); // [N,B,3]
				torch::Tensor d = torch::sqrt(torch::sum(diff * diff, 2));				// [N,B]

				dist = d - radius_c.unsqueeze(0);
				dist = dist * dist;
				dist = dist * samples_area.unsqueeze(1);  

			}
			else{
				torch::Tensor v_sphere = torch::cat({centers_c, radius_c.unsqueeze(1)}, 1); // [B,4]
				torch::Tensor dist_sqem = eval_sqem_batch(
					p.samples_sqem_A_gpu_, p.samples_sqem_b_gpu_, p.samples_sqem_c_gpu_, v_sphere); // [N,B]
				torch::Tensor dist_other;
				if(p.distance_mode_ == SPHERE_EUCLIDEAN_DISTANCE){
					// euclidean part
					torch::Tensor diff = samples_pos.unsqueeze(1) - centers_c.unsqueeze(0); // [N,B,3]
					torch::Tensor d = torch::sqrt(torch::sum(diff * diff, 2));				// [N,B]

					dist_other = d - radius_c.unsqueeze(0);
					dist_other = dist_other * dist_other; // [N,B]
				}
				else{
					// line quadric part
					torch::Tensor v_sphere_line = torch::cat({centers_c, torch::ones({B,1}, fopts)}, 1); // [B,4]
					dist_other = eval_line_quadric_batch(
						p.samples_line_Q_gpu_, v_sphere_line); // [N,B]
				}
				dist = dist_sqem + p.sqem_clustering_lambda_ * dist_other;
			}

			auto min_pair = dist.min(1);
			torch::Tensor dist_min = std::get<0>(min_pair); // [N]
			torch::Tensor dist_idx = std::get<1>(min_pair); // [N] in [0,B)

			torch::Tensor global_idx = dist_idx + k;
			torch::Tensor mask = dist_min < min_dist;

			min_dist = torch::where(mask, dist_min, min_dist);
			argmin = torch::where(mask, global_idx, argmin);
		}

		// back to CPU, map results to PVertex
		torch::Tensor argmin_cpu = argmin.to(torch::kCPU, true);
		auto arg_acc = argmin_cpu.accessor<int64_t, 1>();

		// clear cluster
		parallel_foreach_cell(*p.spheres_, [&](PVertex v) -> bool {
			uint32 v_index = index_of(*p.spheres_, v);
			(*p.spheres_cluster_)[v_index].clear();
			(*p.spheres_cluster_area_)[v_index] = 0.0;
			return true;
		});
		p.samples_sphere_->fill(PVertex());

		for (int64_t i = 0; i < N; ++i)
		{
			const int64_t sid = arg_acc[i];
			if (sid < 0 || sid >= static_cast<int64_t>(p.spheres_gpu_order_.size()))
				continue;

			PVertex sv = p.samples_gpu_order_[i];
			PVertex sp = p.spheres_gpu_order_[sid];

			uint32 sv_idx = index_of(*p.samples_mesh_, sv);
			uint32 sp_idx = index_of(*p.spheres_, sp);

			(*p.samples_sphere_)[sv_idx] = sp;
			(*p.spheres_cluster_)[sp_idx].push_back(sv);
			(*p.spheres_cluster_area_)[sp_idx] += (*p.samples_area_)[sv_idx];
		}
	}

	void augment_insufficient_clusters(PointsParameters& p)
	{
		std::cout << "Augmenting insufficient clusters..." << std::endl;

		std::vector<PVertex> to_process;
		std::vector<uint32> to_process_indices;
		std::vector<uint32> needed;
		std::vector<uint32> samples_start;
		std::vector<uint32> samples_count;

		std::vector<Vec3> all_raw_samples;

		// Add points to insufficient clusters
		foreach_cell(*p.spheres_, [&](PVertex v) -> bool {
			uint32 v_index = index_of(*p.spheres_, v);
			uint32 cluster_size = (*p.spheres_cluster_)[v_index].size();
			if (cluster_size < p.cluster_min_points_)
			{
				const Vec3& center = (*p.spheres_position_)[v_index];
				Scalar radius = (*p.spheres_radius_)[v_index];
				uint32 needed_samples = p.cluster_min_points_ - cluster_size;

				std::mt19937 gen(p.seed_ + v_index);
				std::vector<Vec3> samples = sample_in_sphere_volume(center, radius, needed_samples * 3, gen);

				to_process.push_back(v);
				to_process_indices.push_back(v_index);
				needed.push_back(needed_samples);
				samples_start.push_back(uint32(all_raw_samples.size()));
				samples_count.push_back(uint32(samples.size()));

				all_raw_samples.insert(all_raw_samples.end(), samples.begin(), samples.end());
			}
			return true;
		});

		std::cout << "  Found " << to_process.size() << " clusters to augment" << std::endl;
		std::cout << "  Projecting " << all_raw_samples.size() << " samples..." << std::endl;

		auto [all_projected, all_normals] = project_points_to_alpha_gpu(p, all_raw_samples);

		std::cout << "  Projected: " << all_raw_samples.size() << " -> " << all_projected.size() << " points"
				  << std::endl;
		int total_added = 0;
		for (size_t i = 0; i < to_process.size(); ++i)
		{
			PVertex sphere = to_process[i];
			uint32 sphere_index = to_process_indices[i];
			uint32 needed_samples = needed[i];
			uint32 start_idx = samples_start[i];
			uint32 count = samples_count[i];

			std::vector<PVertex>& cluster = (*p.spheres_cluster_)[sphere_index];
			int added = 0;

			for (int i = start_idx; i < start_idx + count && added < needed_samples; ++i)
			{
				const Vec3& pos = all_projected[i];
				const Vec3& normal = all_normals[i];
				// Check if valid sample
				if (p.samples_spatial_grid_->is_valid_sample(pos, p.grid_cell_size_, *p.samples_position_))
				{
					PVertex new_vertex = add_vertex(*p.samples_mesh_);
					uint32 new_idx = index_of(*p.samples_mesh_, new_vertex);
					(*p.samples_position_)[new_idx] = pos;
					(*p.samples_normal_)[new_idx] = normal;
					(*p.samples_normal_color_)[new_idx] =
						Vec4((normal.x() + 1.0) * 0.5, (normal.y() + 1.0) * 0.5, (normal.z() + 1.0) * 0.5, 1.0);
					cluster.push_back(new_vertex);
					(*p.samples_sphere_)[new_idx] = sphere;
					p.samples_spatial_grid_->insert(pos, new_idx);
					added++;
					total_added++;
				}
			}
		}
		std::cout << "  Added a total of " << total_added << " samples to clusters." << std::endl;
		std::cout << "  Rebuilding KDTree..." << std::endl;
		build_kdtree(p);

		std::cout << "  Recomputing KNN and area..." << std::endl;
		compute_samples_area(p);

		std::cout << "  Recomputing quadrics..." << std::endl;
		compute_quadrics(p);

		points_provider_->emit_connectivity_changed(*p.samples_mesh_);
		std::cout << "  Cluster augmentation complete." << std::endl;
	}

	std::vector<Vec3> sample_in_sphere_volume(const Vec3& center, Scalar radius, int num_samples, std::mt19937& gen)
	{
		std::vector<Vec3> samples;
		samples.reserve(num_samples);

		std::uniform_real_distribution<Scalar> uni(-1.0, 1.0);

		for (int i = 0; i < num_samples; ++i)
		{
			Vec3 offset(uni(gen) * radius, uni(gen) * radius, uni(gen) * radius);
			samples.push_back(center + offset);
		}

		return samples;
	}
	// Power distance clustering: d_power(p, sphere) = |p - center|^2 - radius^2
	void compute_power_cluster(PointsParameters& p)
	{
		// clean cluster affectation
		parallel_foreach_cell(*p.spheres_, [&](PVertex v) -> bool {
			uint32 v_index = index_of(*p.spheres_, v);
			(*p.spheres_cluster_)[v_index].clear();
			(*p.spheres_cluster_area_)[v_index] = 0.0;
			return true;
		});
		p.samples_sphere_->fill(PVertex());

		if (p.nb_spheres_ == 0)
			return;

		parallel_foreach_cell(*p.samples_mesh_, [&](PVertex v) -> bool {
			uint32 v_index = index_of(*p.samples_mesh_, v);
			Scalar a = (*p.samples_area_)[v_index];
			const Vec3& vp = (*p.samples_position_)[v_index];

			Scalar min_power_distance = std::numeric_limits<Scalar>::max();
			PVertex closest_sphere;
			uint32 closest_sphere_index = 0;

			foreach_cell(*p.spheres_, [&](PVertex pv) {
				uint32 pv_index = index_of(*p.spheres_, pv);
				const Vec3& center = (*p.spheres_position_)[pv_index];
				Scalar radius = (*p.spheres_radius_)[pv_index];

				// Power distance: |p - center|^2 - radius^2
				Scalar dist_sq = (vp - center).squaredNorm();
				Scalar power_dist = dist_sq - radius * radius;

				if (power_dist < min_power_distance)
				{
					min_power_distance = power_dist;
					closest_sphere = pv;
					closest_sphere_index = pv_index;
				}
				return true;
			});

			(*p.samples_sphere_)[v_index] = closest_sphere;

			std::lock_guard<std::mutex> lock(spheres_mutex_[closest_sphere_index % spheres_mutex_.size()]);
			(*p.spheres_cluster_)[closest_sphere_index].push_back(v);
			(*p.spheres_cluster_area_)[closest_sphere_index] += a;

			return true;
		});

		//// remove small clusters
		// foreach_cell(*p.spheres_, [&](PVertex v) -> bool {
		//	uint32 v_idx = index_of(*p.spheres_, v);
		//	std::vector<PVertex>& cluster = (*p.spheres_cluster_)[v_idx];
		//	if (cluster.size() < 4)
		//	{
		//		for (PVertex sv : cluster) {
		//			uint32 sv_idx = index_of(*p.samples_mesh_, sv);
		//			(*p.samples_sphere_)[sv_idx] = PVertex();
		//		}
		//		remove_vertex(*p.spheres_, v);
		//		p.nb_spheres_--;
		//	}
		//	return true;
		// });
	}

	// Use power-based clusters to compute sphere neighbors and build skeleton
	void compute_skeleton_power(PointsParameters& p, bool only_neighbors = false)
	{
		// First compute power clusters
		compute_power_cluster(p);

		// Then build skeleton using the power-based cluster assignments
		compute_skeleton(p, only_neighbors);
	}

	void compute_spheres_error(PointsParameters& p)
	{
		parallel_foreach_cell(*p.spheres_, [&](PVertex v) {
			uint32 v_index = index_of(*p.spheres_, v);

			const Vec3& center = (*p.spheres_position_)[v_index];
			Scalar radius = (*p.spheres_radius_)[v_index];
			const std::vector<PVertex>& cluster = (*p.spheres_cluster_)[v_index];

			Scalar cluster_error = 0.0;
			for (PVertex sv : cluster)
			{
				uint32 sv_index = index_of(*p.samples_mesh_, sv);
				Vec3 vp = (*p.samples_position_)[sv_index];
				Scalar a = (*p.samples_area_)[sv_index];

				Scalar dist = 0.0;
				Scalar dist_other = 0.0;
				switch (p.distance_mode_)
				{
				case SPHERE_EUCLIDEAN_DISTANCE: {
					Scalar dist_sqem =
						(*p.samples_quadric_)[sv_index].eval(Vec4(center.x(), center.y(), center.z(), radius));
					dist_other = ((vp - center).norm() - radius);
					dist_other *= dist_other;
					dist_other *= a;
					dist = dist_sqem + p.sqem_clustering_lambda_ * dist_other;
				}
				break;

				case LINE_QUADRIC_DISTANCE: {
					Scalar dist_sqem =
						(*p.samples_quadric_)[sv_index].eval(Vec4(center.x(), center.y(), center.z(), radius));
					// Don't multiply by area here since line quadric already incorporates it
					dist_other = (*p.samples_line_quadric_)[sv_index].eval(center);
					dist = dist_sqem + p.sqem_clustering_lambda_ * dist_other;
				}
				break;

				case PURE_EUCLIDEAN_DISTANCE: {
					dist_other = ((vp - center).norm() - radius);
					dist_other *= dist_other;
					dist = dist_other * a;
				}
				break;
				}
				(*p.samples_error_)[sv_index] = dist;
				cluster_error += dist;
			}

			if ((*p.spheres_cluster_area_)[v_index] > 0)
				(*p.spheres_error_)[v_index] = cluster_error / (*p.spheres_cluster_area_)[v_index];
			else
				(*p.spheres_error_)[v_index] = 0.0;

			(*p.spheres_error_not_normalized_)[v_index] = cluster_error;

			return true;
		});

		p.min_error_ = std::numeric_limits<Scalar>::max();
		p.max_error_ = std::numeric_limits<Scalar>::min();
		p.max_error_sphere_ = PVertex();
		p.total_error_ = 0.0;
		p.total_error_not_normalized_ = 0.0;

		foreach_cell(*p.spheres_, [&](PVertex v) -> bool {
			uint32 v_idx = index_of(*p.spheres_, v);
			Scalar error = (*p.spheres_error_)[v_idx];
			Scalar error_not_normalized = (*p.spheres_error_not_normalized_)[v_idx];

			if (error < p.min_error_)
				p.min_error_ = error;
			if (error > p.max_error_)
			{
				p.max_error_ = error;
				p.max_error_sphere_ = v;
			}
			p.total_error_ += error;
			p.total_error_not_normalized_ += error_not_normalized;
			return true;
		});

		p.total_error_diff_ = std::abs(p.total_error_ - p.last_total_error_);
		p.last_total_error_ = p.total_error_;
		std::cout << "compute_spheres_error end" << std::endl;
	}

	void update_spheres_color(PointsParameters& p)
	{
		parallel_foreach_cell(*p.spheres_, [&](PVertex v) -> bool {
			uint32 v_index = index_of(*p.spheres_, v);
			if (p.error_as_spheres_color_)
				(*p.spheres_color_)[v_index] =
					color_map((*p.spheres_error_)[v_index], p.min_error_, p.max_error_, p.spheres_transparency_);
			else
			{
				const Vec4& c = (*p.spheres_cluster_color_)[v_index];
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

	void update_sphere_euclidean(PointsParameters& p, PVertex sphere)
	{
		uint32 sphere_index = index_of(*p.spheres_, sphere);

		const std::vector<PVertex>& cluster = (*p.spheres_cluster_)[sphere_index];
		if (cluster.empty())
			return;

		Vec3 c = (*p.spheres_position_)[sphere_index];

		Eigen::MatrixXd J(cluster.size(), 3);
		J.setZero();
		Eigen::VectorXd b(cluster.size());
		b.setZero();
		uint32 idx = 0;
		Eigen::VectorXd s(3);
		Scalar r_fixed = p.alpha_; // we only want optimize center, the radius should be fixed as alpha
		s << c[0], c[1], c[2];
		for (uint32 i = 0; i < 10; ++i)
		{
			idx = 0;
			for (PVertex v : cluster)
			{
				uint32 v_index = index_of(*p.samples_mesh_, v);
				const Vec3& pos = (*p.samples_position_)[v_index];

				Vec3 d = pos - Vec3(s(0), s(1), s(2));
				Scalar l = d.norm();

				Scalar a = std::sqrt((*p.samples_area_)[v_index]);
				J.row(idx) = (-d / l) * a;
				b(idx) = -(l - r_fixed) * a;

				++idx;
			};

			Eigen::Matrix3d H = J.transpose() * J;
			Eigen::Vector3d g = J.transpose() * b;
			Eigen::VectorXd delta_s = H.ldlt().solve(b);
			s += delta_s;
			if (delta_s.norm() < 1e-6) // stop early if converged
				break;
		}

		c = s;

		(*p.spheres_position_)[sphere_index] = c;
		(*p.spheres_radius_)[sphere_index] = r_fixed;
	}

	void update_sphere_sqem(PointsParameters& p, PVertex sphere)
	{
		uint32 sphere_index = index_of(*p.spheres_, sphere);

		const std::vector<PVertex>& cluster = (*p.spheres_cluster_)[sphere_index];
		if (cluster.empty())
			return;

		Vec3 c = (*p.spheres_position_)[sphere_index];
		Scalar r = (*p.spheres_radius_)[sphere_index];

		/*Spherical_Quadric q;
		for (PVertex v : cluster)
		{
			uint32 v_index = index_of(*p.samples_mesh_, v);
			Scalar a = (*p.samples_area_)[v_index];
			q += (*p.samples_quadric_)[v_index] * a;
		}
		Scalar radius = 0.0;
		SQEM_CASE sc = q.well_conditioned(radius);*/
		/*if (sc != SQEM_CASE::Case4_Degenerate)
		{*/
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
				uint32 v_index = index_of(*p.samples_mesh_, v);
				const Vec3& pos = (*p.samples_position_)[v_index];

				// SQEM energy
				Eigen::Vector4d lhs = Eigen::Vector4d::Zero();
				Scalar rhs = 0.0;
				const Vec3& n = (*p.samples_normal_)[v_index];
				Vec4 n4 = Vec4(n.x(), n.y(), n.z(), 1.0);
				Scalar a = sqrt((*p.samples_area_)[v_index] / (p.knn_k_ + 1.0));
				lhs += -n4 * a;
				rhs += -1.0 * ((pos - Vec3(s(0), s(1), s(2))).dot(n) - s(3)) * a;

				for (PVertex vn : (*p.samples_knn_)[v_index])
				{
					uint32 vn_index = index_of(*p.samples_mesh_, vn);
					const Vec3& pn = (*p.samples_position_)[vn_index];
					const Vec3& nn = (*p.samples_normal_)[vn_index];
					Vec4 nn4 = Vec4(nn.x(), nn.y(), nn.z(), 1.0);
					Scalar an = sqrt((*p.samples_area_)[vn_index] / (p.knn_k_ + 1.0));
					lhs += -nn4 * an;
					rhs += -1.0 * ((pn - Vec3(s(0), s(1), s(2))).dot(nn) - s(3)) * an;
				}
				J.row(idx) = lhs;
				b(idx) = rhs;
				++idx;

				// distance energy
				Vec3 d = pos - Vec3(s(0), s(1), s(2));
				Scalar l = d.norm();

				Scalar a_dist = sqrt((*p.samples_area_)[v_index]);
				J.row(idx) =
					Eigen::Vector4d(-(d[0] / l), -(d[1] / l), -(d[2] / l), -1.0) * a_dist * p.sqem_update_lambda_;
				b(idx) = -(l - s(3)) * a_dist * p.sqem_update_lambda_; // scale the row by the update lambda

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
		//}
		// else
		//{
		//	std::cout << "Sphere " << sphere_index << " is not well conditioned, using shrinking ball" << std::endl;
		//	// apply shrinking ball
		//	std::pair<uint32, Scalar> knn_res;
		//	p.samples_kdtree_->find_nn(c, &knn_res);
		//	PVertex nn = p.samples_kdtree_vertices_[knn_res.first];
		//	Vec3 closest_point_position = (*p.position_)[index_of(*p.points_, nn)];
		//	Vec3 closest_point_dir = (closest_point_position - c).normalized();
		//	auto [center, radius, q] = geometry::shrinking_ball_center(
		//		closest_point_position, closest_point_dir, p.samples_kdtree_, p.samples_kdtree_vertices_, p.alpha_*1.5);
		//	c = center;
		//	r = radius;
		//}

		(*p.spheres_position_)[sphere_index] = c;
		(*p.spheres_radius_)[sphere_index] = r;
	}

	void update_sphere_line_quadric_distance(PointsParameters& p, PVertex sphere)
	{
		uint32 sphere_index = index_of(*p.spheres_, sphere);

		const std::vector<PVertex>& cluster = (*p.spheres_cluster_)[sphere_index];
		Vec3 c = (*p.spheres_position_)[sphere_index];
		Scalar r = (*p.spheres_radius_)[sphere_index];
		Spherical_Quadric q;
		Line_Quadric lq;
		Scalar area = 0.0;
		Vec3 h;
		h.setZero();
		for (PVertex v : cluster)
		{
			uint32 v_index = index_of(*p.samples_mesh_, v);
			Scalar weight = value<Scalar>(*p.samples_mesh_, p.samples_area_, v);
			if (weight <= 0.0)
				std::cout << "Warning: sample with zero volume weight in sphere " << sphere_index << std::endl;
			q += (*p.samples_quadric_)[v_index] * weight;
			h += weight * (*p.samples_position_)[v_index];
			lq += (*p.samples_line_quadric_)[v_index] * weight;

			area += weight;
		}

		/*Mat4 A = q._A;
		Vec4 b = q._b;*/

		Mat4 Ql = lq.get_quadric().matrix();
		Mat3 Al = Ql.block<3, 3>(0, 0);
		Vec3 bl = -Ql.block<3, 1>(0, 3);

		/*Mat4 Al_ext = Mat4::Zero();
		Al_ext.block<3, 3>(0, 0) = Al;
		Vec4 bl_ext = Vec4::Zero();
		bl_ext.head<3>() = bl;
		Mat4 A_c = A + p.sqem_update_lambda_ * Al_ext;
		Vec4 b_c = b + p.sqem_update_lambda_ * bl_ext;
		Vec4 s = A_c.ldlt().solve(b_c);
		c = s.head<3>();
		r = s[3];*/
		Mat3 As = q._A.block<3, 3>(0, 0);
		Vec3 bs = q._b.head<3>();
		Vec3 Asr = q._A.block<3, 1>(0, 3);

		Mat3 A = As + p.sqem_update_lambda_ * Al;
		Vec3 b = (bs + p.sqem_update_lambda_ * bl) - Asr * p.alpha_;

		c = A.ldlt().solve(b);
		r = p.alpha_;
		(*p.spheres_position_)[sphere_index] = c;
		(*p.spheres_radius_)[sphere_index] = r;
	}

	void correct_sphere(PointsParameters& p, PVertex v)
	{
		uint32 v_index = index_of(*p.spheres_, v);

		Vec3& c = (*p.spheres_position_)[v_index];
		Scalar& r = (*p.spheres_radius_)[v_index];

		bool inside = p.samples_winding_number_->is_inside(c);

		std::pair<uint32, Scalar> k_res;
		p.samples_kdtree_->find_nn(c, &k_res);
		uint32 k_idx = index_of(*p.samples_mesh_, p.samples_kdtree_vertices_[k_res.first]);
		Vec3 closest_pos = (*p.samples_position_)[k_idx];
		Vec3 dir = (closest_pos - c).normalized();
		if (!inside)
			dir = -dir;

		c = closest_pos - dir * p.alpha_;
		r = p.alpha_;

		(*p.spheres_position_)[v_index] = c;
		(*p.spheres_radius_)[v_index] = r;
	}
	void update_spheres_batch_with_projection(PointsParameters& p)
	{
		compute_clusters(p);

		std::vector<PVertex> all_spheres;
		std::vector<uint32> sphere_sample_start;
		std::vector<uint32> sphere_sample_count;
		std::vector<Vec3> all_raw_samples;

		const int num_projection_samples = 30;

		foreach_cell(*p.spheres_, [&](PVertex sphere) -> bool {
			uint32 sphere_index = index_of(*p.spheres_, sphere);

			Vec3 c = (*p.spheres_position_)[sphere_index];
			Scalar r = (*p.spheres_radius_)[sphere_index];

			// ??
			std::mt19937 gen(p.seed_ + sphere_index);
			std::vector<Vec3> sphere_samples = sample_in_sphere_volume(c, r, num_projection_samples, gen);

			// ????
			all_spheres.push_back(sphere);
			sphere_sample_start.push_back(static_cast<uint32>(all_raw_samples.size()));
			sphere_sample_count.push_back(static_cast<uint32>(sphere_samples.size()));

			// ?????
			all_raw_samples.insert(all_raw_samples.end(), sphere_samples.begin(), sphere_samples.end());

			return true;
		});

		std::cout << "Collected " << all_raw_samples.size() << " samples from " << all_spheres.size() << " spheres"
				  << std::endl;

		// === 2. ????????? ===
		auto [all_projected, all_normals] = project_points_to_alpha_gpu(p, all_raw_samples);

		if (all_projected.empty())
		{
			std::cerr << "Projection failed, falling back to standard update" << std::endl;

			// ???????
			parallel_foreach_cell(*p.spheres_, [&](PVertex v) -> bool {
				if (p.use_sqem_term_)
					update_sphere_sqem(p, v);
				else
					update_sphere_euclidean(p, v);
				return true;
			});
			return;
		}

		// === 3. ??????UDF??????? ===
		std::cout << "Verifying projected points..." << std::endl;
		std::vector<bool> is_valid(all_projected.size(), false);
		int valid_count = 0;

		NeuralFieldForward udf = make_neural_field_forward(p);
		BatchUDFResult udf_result = udf.forward_batch(all_projected);
		if (udf_result.ok)
		{
			for (size_t i = 0; i < all_projected.size(); ++i)
			{
				Scalar error = std::abs(udf_result.values[i] - p.alpha_);
				is_valid[i] = (error < p.tol_);
				if (is_valid[i])
					valid_count++;
			}
		}

		std::cout << "Valid projected points: " << valid_count << " / " << all_projected.size() << std::endl;

		parallel_foreach_cell(*p.spheres_, [&](PVertex sphere) -> bool {
			uint32 sphere_index = index_of(*p.spheres_, sphere);

			auto it = std::find(all_spheres.begin(), all_spheres.end(), sphere);
			if (it == all_spheres.end())
				return true;

			size_t sphere_idx = std::distance(all_spheres.begin(), it);
			uint32 start = sphere_sample_start[sphere_idx];
			uint32 count = sphere_sample_count[sphere_idx];

			std::vector<Vec3> valid_projected;
			std::vector<Vec3> valid_normals;

			for (uint32 i = start; i < start + count; ++i)
			{
				if (i < is_valid.size() && is_valid[i])
				{
					valid_projected.push_back(all_projected[i]);
					valid_normals.push_back(all_normals[i]);
				}
			}

			if (valid_projected.empty())
				return true;

			Vec3 c = (*p.spheres_position_)[sphere_index];
			const Scalar r = p.alpha_;

			const int num_valid_proj = static_cast<int>(valid_projected.size());

			Eigen::MatrixXd J(num_valid_proj, 3);
			J.setZero();
			Eigen::VectorXd b(num_valid_proj);
			b.setZero();

			Eigen::Vector3d center(c[0], c[1], c[2]);

			for (uint32 iter = 0; iter < 10; ++iter)
			{
				int row = 0;

				// ??????|pos - center| = alpha
				for (size_t i = 0; i < valid_projected.size(); ++i)
				{
					const Vec3& pos = valid_projected[i];

					Vec3 d = pos - Vec3(center[0], center[1], center[2]);
					Scalar l = d.norm();

					if (l < 1e-10) // ????
						continue;

					Scalar weight = 1.0;

					// **??????? - ????alpha**
					Scalar residual = l - r;

					// **Jacobian?d(residual)/d(center) = -d/|d|**
					J.row(row) = Eigen::Vector3d(-(d[0] / l), -(d[1] / l), -(d[2] / l)) * weight;
					b(row) = -residual * weight;
					++row;
				}

				// ???????????
				if (row < 3)
					break;

				// ??????
				J.conservativeResize(row, 3);
				b.conservativeResize(row);

				// Solve: (J^T * J) * delta_center = J^T * b
				Eigen::LDLT<Eigen::MatrixXd> solver(J.transpose() * J);
				Eigen::Vector3d delta_center = solver.solve(J.transpose() * b);

				center += delta_center;

				// ????
				if (delta_center.norm() < 1e-6)
					break;
			}

			// Update sphere center (radius remains alpha)
			(*p.spheres_position_)[sphere_index] = Vec3(center[0], center[1], center[2]);
			(*p.spheres_radius_)[sphere_index] = p.alpha_; // **???????alpha**

			return true;
		});

		// Sphere correction
		if (p.sphere_correction_ && p.sphere_correction_mode_ == CORRECT_ALWAYS)
		{
			parallel_foreach_cell(*p.spheres_, [&](PVertex v) -> bool {
				correct_sphere(p, v);
				return true;
			});
		}

		parallel_foreach_cell(*p.spheres_, [&](PVertex v) -> bool {
			(*p.spheres_do_not_split_)[index_of(*p.spheres_, v)] = false;
			return true;
		});
	}
	void update_spheres(PointsParameters& p)
	{
		compute_clusters(p);

		//compute_clusters_gpu(p);

		parallel_foreach_cell(*p.spheres_, [&](PVertex v) -> bool {
			switch (p.distance_mode_)
			{
			case SPHERE_EUCLIDEAN_DISTANCE: {
				update_sphere_sqem(p, v);
			}
			break;
			case LINE_QUADRIC_DISTANCE: {
				update_sphere_line_quadric_distance(p, v);
			}
			break;
			case PURE_EUCLIDEAN_DISTANCE: {
				update_sphere_euclidean(p, v);
			}
			break;
			}
			if (p.sphere_correction_ && p.sphere_correction_mode_ == CORRECT_ALWAYS)
				correct_sphere(p, v);

			(*p.spheres_do_not_split_)[index_of(*p.spheres_, v)] = false;
			return true;
		});
		// update_spheres_batch_with_projection(p);

		compute_spheres_error(p);

		if (p.auto_split_ && (p.total_error_diff_ < 1e-5 || p.iteration_count_ % 10 == 0))
		{
			switch (p.auto_split_mode_)
			{
			case ERROR_THRESHOLD: {
				if (p.max_error_ > p.auto_split_error_threshold_)
				{
					compute_skeleton(p, true);

					if (p.sphere_correction_ && p.sphere_correction_mode_ == CORRECT_ON_SPLIT)
					{
						parallel_foreach_cell(*p.spheres_, [&](PVertex v) -> bool {
							correct_sphere(p, v);
							return true;
						});
					}

					std::vector<PVertex> sorted_spheres;
					sorted_spheres.reserve(p.nb_spheres_);
					foreach_cell(*p.spheres_, [&](PVertex v) -> bool {
						sorted_spheres.push_back(v);
						return true;
					});
					std::sort(sorted_spheres.begin(), sorted_spheres.end(), [&](PVertex a, PVertex b) {
						return (*p.spheres_error_)[index_of(*p.spheres_, a)] >
							   (*p.spheres_error_)[index_of(*p.spheres_, b)];
					});

					uint32 to_split_max = std::min(uint32(std::ceil(p.nb_spheres_ * 0.2)), 10u);
					for (PVertex sphere : sorted_spheres)
					{
						uint32 s_index = index_of(*p.spheres_, sphere);
						Scalar error = (*p.spheres_error_)[s_index];
						if (error < p.auto_split_error_threshold_)
							break;
						if (to_split_max > 0 && !(*p.spheres_do_not_split_)[s_index])
						{
							for (PVertex neighbor : (*p.spheres_neighbor_clusters_)[s_index])
								(*p.spheres_do_not_split_)[index_of(*p.spheres_, neighbor)] = true;
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
					compute_skeleton(p, true);

					if (p.sphere_correction_ && p.sphere_correction_mode_ == CORRECT_ON_SPLIT)
					{
						parallel_foreach_cell(*p.spheres_, [&](PVertex v) -> bool {
							correct_sphere(p, v);
							return true;
						});
					}

					std::vector<PVertex> sorted_spheres;
					sorted_spheres.reserve(p.nb_spheres_);
					foreach_cell(*p.spheres_, [&](PVertex v) -> bool {
						sorted_spheres.push_back(v);
						return true;
					});
					std::sort(sorted_spheres.begin(), sorted_spheres.end(), [&](PVertex a, PVertex b) {
						return (*p.spheres_error_)[index_of(*p.spheres_, a)] >
							   (*p.spheres_error_)[index_of(*p.spheres_, b)];
					});

					// uint32 to_split_max = std::min(uint32(std::ceil(p.nb_spheres_ * 0.2)), 100u);
					uint32 to_split_max = std::max(0.5 * p.nb_spheres_, 1.0);
					for (PVertex sphere : sorted_spheres)
					{
						uint32 s_index = index_of(*p.spheres_, sphere);
						if (p.auto_split_max_nb_spheres_ - p.nb_spheres_ <= 0)
							break;
						if (to_split_max > 0 && !(*p.spheres_do_not_split_)[s_index])
						{
							for (PVertex neighbor : (*p.spheres_neighbor_clusters_)[s_index])
								(*p.spheres_do_not_split_)[index_of(*p.spheres_, neighbor)] = true;
							split_sphere(p, sphere);
							--to_split_max;
						}
					}
				}
			}
			break;
			}
		}

		if (!p.running_)
			update_render_data(p);
	}

	void remove_sphere(PointsParameters& p, PVertex v)
	{
		const std::vector<PVertex>& cluster = (*p.spheres_cluster_)[index_of(*p.spheres_, v)];
		for (PVertex s : cluster)
		{
			uint32 s_idx = index_of(*p.samples_mesh_, s);
			(*p.samples_sphere_)[s_idx] = PVertex();
		}
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

	void compute_skeleton(PointsParameters& p, bool only_neighbors = false)
	{
		// clear graph
		parallel_foreach_cell(*p.spheres_, [&](PVertex v) -> bool {
			uint32 v_index = index_of(*p.spheres_, v);
			(*p.spheres_neighbor_clusters_)[v_index].clear();
			return true;
		});

		foreach_cell(*p.samples_mesh_, [&](PVertex v) -> bool {
			uint32 v_index = index_of(*p.samples_mesh_, v);
			PVertex v_sphere = (*p.samples_sphere_)[v_index];
			for (PVertex w : (*p.samples_knn_)[v_index])
			{
				PVertex w_sphere = (*p.samples_sphere_)[index_of(*p.samples_mesh_, w)];
				if (v_sphere.is_valid() && w_sphere.is_valid() && v_sphere != w_sphere)
				{
					uint32 v_index = index_of(*p.spheres_, v_sphere);
					uint32 w_index = index_of(*p.spheres_, w_sphere);
					(*p.spheres_neighbor_clusters_)[v_index].insert(w_sphere);
					(*p.spheres_neighbor_clusters_)[w_index].insert(v_sphere);
				}
			}
			return true;
		});

		if (only_neighbors)
			return;

		clear(*p.skeleton_);
		p.skeleton_faces_map_.clear();
		auto get_face_key = [](uint32 i1, uint32 i2, uint32 i3) -> NMFaceKey {
			std::array<uint32, 3> key = {i1, i2, i3};
			std::sort(key.begin(), key.end());
			return key;
		};
		auto spheres_skeleton_vertex_map =
			add_attribute<NMVertex, PVertex>(*p.spheres_, "__spheres_skeleton_vertex_map");
		std::vector<std::array<uint32, 4>> raw_tets;
		foreach_cell(*p.spheres_, [&](PVertex pv) -> bool {
			uint32 pv_index = index_of(*p.spheres_, pv);
			NMVertex nmv = add_vertex(*p.skeleton_);
			(*p.skeleton_position_)[index_of(*p.skeleton_, nmv)] = (*p.spheres_position_)[pv_index];
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
				uint32 n_index = index_of(*p.spheres_, neighbor);
				NMVertex nmv2 = (*spheres_skeleton_vertex_map)[n_index];
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
			uint32 idx1 = index_of(*p.spheres_, pv);
			NMVertex nmv1 = (*spheres_skeleton_vertex_map)[idx1];
			const std::set<PVertex>& n_pv = (*p.spheres_neighbor_clusters_)[idx1];
			for (const PVertex& ne1 : n_pv)
			{
				uint32 idx2 = index_of(*p.spheres_, ne1);
				if (idx1 >= idx2)
					continue;
				NMVertex nmv2 = (*spheres_skeleton_vertex_map)[idx2];
				const std::set<PVertex>& ne_ne1 = (*p.spheres_neighbor_clusters_)[idx2];
				for (const PVertex& ne2 : ne_ne1)
				{
					if (n_pv.find(ne2) == n_pv.end())
						continue;
					uint32 idx3 = index_of(*p.spheres_, ne2);
					if (idx2 >= idx3)
						continue;

					NMVertex nmv3 = (*spheres_skeleton_vertex_map)[idx3];

					std::vector<NMEdge> edges;
					edges.reserve(3);
					edges.push_back(edge_indices[{index_of(*p.skeleton_, nmv1), index_of(*p.skeleton_, nmv2)}]);
					edges.push_back(edge_indices[{index_of(*p.skeleton_, nmv2), index_of(*p.skeleton_, nmv3)}]);
					edges.push_back(edge_indices[{index_of(*p.skeleton_, nmv1), index_of(*p.skeleton_, nmv3)}]);
					NMFace new_face = add_face(*p.skeleton_, edges);

					p.skeleton_faces_map_[get_face_key(idx1, idx2, idx3)] = new_face;

					const std::set<PVertex>& ne_ne2 = (*p.spheres_neighbor_clusters_)[idx3];

					for (const PVertex& ne3 : ne_ne2)
					{
						uint32 idx4 = index_of(*p.spheres_, ne3);
						if (idx3 >= idx4)
							continue;
						bool connected_v1 = (n_pv.find(ne3) != n_pv.end());
						bool connected_v2 = (ne_ne1.find(ne3) != ne_ne1.end());
						if (connected_v1 && connected_v2)
						{
							raw_tets.push_back({idx1, idx2, idx3, idx4});
						}
					}
				}
			}
			return true;
		});

		// Resolve Tets
		p.skeleton_tets_.clear();
		p.skeleton_tets_.reserve(raw_tets.size());
		foreach_cell(*p.skeleton_, [&](NMFace f) -> bool {
			value<std::set<std::size_t>>(*p.skeleton_, p.incident_tets_, f).clear();
			value<Vec3>(*p.skeleton_, p.skeleton_face_color_, f) = Vec3(0.0, 0.0, 0.0);
			return true;
		});
		std::size_t tet_index = 0;
		for (const auto& rt : raw_tets)
		{
			Tet new_tet;
			uint32 v[4] = {rt[0], rt[1], rt[2], rt[3]};

			NMFaceKey keys[4] = {get_face_key(v[1], v[2], v[3]), get_face_key(v[0], v[2], v[3]),
								 get_face_key(v[0], v[1], v[3]), get_face_key(v[0], v[1], v[2])};

			for (uint32 i = 0; i < 4; ++i)
			{
				auto it = p.skeleton_faces_map_.find(keys[i]);
				if (it != p.skeleton_faces_map_.end())
				{
					NMFace f = it->second;
					new_tet.faces[i] = f;

					value<std::set<size_t>>(*p.skeleton_, p.incident_tets_, f).insert(tet_index);

					value<Vec3>(*p.skeleton_, p.skeleton_face_color_, f) = Vec3(0.8, 0.5, 0.5);
				}
			}
			new_tet.tet_id = tet_index;
			p.skeleton_tets_.insert({tet_index, new_tet});
			tet_index++;
		}
		parallel_foreach_cell(*p.skeleton_, [&](NMFace f) -> bool {
			if (value<std::set<size_t>>(*p.skeleton_, p.incident_tets_, f).size() == 1)
			{
				// std::cout << "Find simple face" << std::endl;
				value<Vec3>(*p.skeleton_, p.skeleton_face_color_, f) = Vec3(1.0, 0.0, 0.0);
			}
			return true;
		});

		compute_edge_degree(p);
		// std::cout << "Found " << p.skeleton_tets_.size() << " tets in the skeleton." << std::endl;
		remove_attribute<PVertex>(*p.spheres_, spheres_skeleton_vertex_map);

		non_manifold_provider_->emit_attribute_changed(*p.skeleton_, p.skeleton_face_color_.get());
		non_manifold_provider_->emit_attribute_changed(*p.skeleton_, p.skeleton_edge_color_.get());
	}

protected:
	void split_sphere(PointsParameters& p, PVertex sphere)
	{
		if (!sphere.is_valid())
			return;
		uint32 s_index = index_of(*p.spheres_, sphere);
		Vec3 c = (*p.spheres_position_)[s_index];
		Scalar r = (*p.spheres_radius_)[s_index];

		// find the point in the cluster with the max error
		const std::vector<PVertex>& cluster = (*p.spheres_cluster_)[s_index];
		if (cluster.empty())
			return;

		Scalar max_err = -1.0;
		PVertex max_err_v;
		for (PVertex v : cluster)
		{
			uint32 v_idx = index_of(*p.samples_mesh_, v);
			Scalar err = (*p.samples_error_)[v_idx];
			if (err > max_err)
			{
				max_err = err;
				max_err_v = v;
			}
		}

		if (!max_err_v.is_valid())
			return;

		uint32 max_err_v_idx = index_of(*p.samples_mesh_, max_err_v);

		PVertex new_sphere = add_vertex(*p.spheres_);
		uint32 new_s_index = index_of(*p.spheres_, new_sphere);

		(*p.spheres_position_)[new_s_index] = (*p.samples_ma_position_)[max_err_v_idx];
		(*p.spheres_radius_)[new_s_index] = (*p.samples_ma_radius_)[max_err_v_idx];
		(*p.spheres_parent_)[new_s_index] = sphere;
		(*p.spheres_cluster_color_)[new_s_index] =
			Vec4(0.5 + 0.5 * (rand() % 256) / 256.0, 0.5 + 0.5 * (rand() % 256) / 256.0,
				 0.5 + 0.5 * (rand() % 256) / 256.0, 1.0);

		p.nb_spheres_++;
	}

	//------------------------------//
	//-----Topology correction------//
	//------------------------------//

	void compute_edge_degree(PointsParameters& p)
	{
		parallel_foreach_cell(*p.skeleton_, [&](NMEdge e) {
			auto in_face = incident_faces(*p.skeleton_, e);

			if (in_face.size() == 2)
				value<Vec3>(*p.skeleton_, p.skeleton_edge_color_, e) = Vec3(0.0, 1.0, 0.0);

			return true;
		});
	}

	bool is_simple_face_3d(PointsParameters& p, NMFace& f)
	{
		if (!f.is_valid())
			return false;
		uint32 idf = index_of(*p.skeleton_, f);
		if ((*p.incident_tets_)[idf].size() != 1)
			return false;
		auto edges = incident_edges(*p.skeleton_, f);

		for (NMEdge e : edges)
		{
			auto in_faces = incident_faces(*p.skeleton_, e);
			return in_faces.size() <= 2;
		}
		return false;
	}

	bool is_simple_face_2d(PointsParameters& p, NMFace& f)
	{
		if (!f.is_valid())
			return false;
		uint32 idf = index_of(*p.skeleton_, f);
		auto edges = incident_edges(*p.skeleton_, f);
		for (NMEdge e : edges)
		{
			auto in_faces = incident_faces(*p.skeleton_, e);
			if (in_faces.size() == 1)
				return true;
		}
		return false;
	}

	bool is_simple_tet(PointsParameters& p, Tet& t)
	{
		for (uint32 i = 0; i < 4; ++i)
		{
			NMFace f = t.faces[i];
			if (!f.is_valid())
			{
				continue;
			}
			uint32 idf = index_of(*p.skeleton_, f);
			if (is_simple_face_3d(p, f))
			{
				return true;
			}
		}
		return false;
	}

	// find a face of tet t that has at least one edge with degree 2
	NMFace find_simple_tet_face(PointsParameters& p, Tet& t)
	{

		NMEdge e;
		for (uint32 i = 0; i < 4; ++i)
		{
			NMFace f = t.faces[i];
			if (!f.is_valid())
				continue;
			if (is_simple_face_3d(p, f))
			{
				return f;
			}
		}
		return NMFace();
	}

	void skeleton_post_pocessing(PointsParameters& p)
	{
		// compute_edge_degree(p);
		std::queue<std::size_t> Q_tet;
		std::queue<NMFace> Q_face;
		std::vector<bool> visited_tet(p.skeleton_tets_.size(), false);
		std::unordered_map<uint32, std::size_t> face_id_map;
		std::size_t face_count = 0;
		foreach_cell(*p.skeleton_, [&](NMFace f) -> bool {
			face_id_map[index_of(*p.skeleton_, f)] = face_count++;
			return true;
		});
		auto get_face_id = [&](NMFace f) -> uint32 {
			uint32 idf = index_of(*p.skeleton_, f);
			return face_id_map[idf];
		};
		std::vector<bool> visited_face(face_count, false);

		// Remove simple tet/face then remove simple face/edge iteratively

		// push all simple tets in queue
		for (auto& [id, tet] : p.skeleton_tets_)
		{
			if (is_simple_tet(p, tet))
			{
				auto tet_id = tet.tet_id;
				visited_tet[tet_id] = true;

				Q_tet.push(tet_id);
			}
		}
		uint32 removed_face = 0;
		uint32 removed_edge = 0;
		while (!Q_tet.empty() || !Q_face.empty())
		{
			if (!Q_tet.empty())
			{
				std::size_t current_tet_id = Q_tet.front();

				Q_tet.pop();

				auto it = p.skeleton_tets_.find(current_tet_id);
				if (it == p.skeleton_tets_.end())
					continue;
				Tet& current_tet = it->second;
				auto faces = current_tet.faces;

				// current_tet.print_tet_info(p);

				NMFace face_to_delete = find_simple_tet_face(p, current_tet);
				std::size_t tet_id = current_tet.tet_id;

				if (!face_to_delete.is_valid()) // should not happen
				{
					std::cout << "Error: It is not a simple tet" << std::endl;
					current_tet.print_tet_info(p);
					continue;
				}
				auto in_edges = incident_edges(*p.skeleton_, face_to_delete);
				remove_face(*p.skeleton_, face_to_delete);
				removed_face++;
				p.skeleton_tets_.erase(tet_id);
				// update incident tet faces
				for (uint32 i = 0; i < 4; ++i)
				{
					NMFace f = faces[i];
					if (!f.is_valid())
						continue;
					uint32 idf = index_of(*p.skeleton_, f);
					(*p.incident_tets_)[idf].erase(tet_id);
					if ((*p.incident_tets_)[idf].size() == 1)
					{
						// push the other tet to the queue, now it is a simple tet
						size_t other_tet_id = *(*p.incident_tets_)[idf].begin();
						if (is_simple_tet(p, p.skeleton_tets_[other_tet_id]) && visited_tet[other_tet_id] == false)
						{
							visited_tet[other_tet_id] = true;
							Q_tet.push(other_tet_id);
						}
					}
				}

				// update edge degrees
				for (NMEdge e : in_edges)
				{
					uint32 ide = index_of(*p.skeleton_, e);
					auto in_faces = incident_faces(*p.skeleton_, e);
					auto degree = in_faces.size();
					if (degree == 0) // should not happen
					{
						std::cout << "Error: edge degree is already zero" << std::endl;
						remove_edge(*p.skeleton_, e);
						continue;
					}

					for (NMFace f : in_faces)
					{
						uint32 idf = index_of(*p.skeleton_, f);
						std::size_t id_vector = get_face_id(f);
						if (is_simple_face_2d(p, f) && visited_face[idf] == false)
						{
							visited_face[idf] = true;
							Q_face.push(f);
						}
					}
				}
			}
			else
			{
				NMFace current_face = Q_face.front();
				Q_face.pop();
				if (!current_face.is_valid())
				{
					std::cout << "Warning: face to delete is not valid. This should not happen" << std::endl;
				}
				if (!is_simple_face_2d(p, current_face))
				{
					std::cout << "Warning: face to delete is not simple. This should not happen" << std::endl;
					continue;
				}
				auto in_edges = incident_edges(*p.skeleton_, current_face);
				auto in_tets = (*p.incident_tets_)[index_of(*p.skeleton_, current_face)];
				remove_face(*p.skeleton_, current_face);
				removed_face++;
				for (NMEdge e : in_edges)
				{
					uint32 ide = index_of(*p.skeleton_, e);
					auto in_faces = incident_faces(*p.skeleton_, e);
					auto degre = in_faces.size();
					/* auto in_faces = incident_faces(*p.skeleton_, e);
				   for (NMFace f : in_faces)
				   {
					   if (is_simple_face_2d(p, f) &&f != current_face)
						   Q_face.push(f);
				   }*/
					if (degre == 0)
					{
						remove_edge(*p.skeleton_, e);
						removed_edge++;
					}
				}
			}
		}
		std::cout << "Removed " << removed_face << " faces in skeleton post-processing." << std::endl;
		std::cout << p.skeleton_tets_.size() << " tets remain after post-processing." << std::endl;

		foreach_cell(*p.skeleton_, [&](NMEdge e) {
			auto in_face = incident_faces(*p.skeleton_, e);
			if (in_face.size() == 1)
				value<Vec3>(*p.skeleton_, p.skeleton_edge_color_, e) = Vec3(0.0, 0.0, 1.0);
			return true;
		});
		parallel_foreach_cell(*p.skeleton_, [&](NMFace f) -> bool {
			auto in_tets = (*p.incident_tets_)[index_of(*p.skeleton_, f)];
			if (in_tets.size() == 1)
			{
				// std::cout << "Find simple face" << std::endl;
				value<Vec3>(*p.skeleton_, p.skeleton_face_color_, f) = Vec3(1.0, 0.0, 0.0);
			}
			else if (in_tets.size() > 1)
			{
				value<Vec3>(*p.skeleton_, p.skeleton_face_color_, f) = Vec3(0.8, 0.5, 0.5);
			}
			else
			{
				value<Vec3>(*p.skeleton_, p.skeleton_face_color_, f) = Vec3(0.0, 0.0, 0.0);
			}
			return true;
		});
		non_manifold_provider_->emit_connectivity_changed(*p.skeleton_);
		non_manifold_provider_->emit_attribute_changed(*p.skeleton_, p.skeleton_edge_color_.get());
		non_manifold_provider_->emit_attribute_changed(*p.skeleton_, p.skeleton_face_color_.get());
	}

protected:
	void update_render_data(PointsParameters& p)
	{
		if (p.running_)
		{
			std::lock_guard<std::mutex> lock(p.mutex_);
			points_provider_->emit_connectivity_changed(*p.spheres_);
			points_provider_->emit_attribute_changed(*p.spheres_, p.spheres_position_.get());
			points_provider_->emit_attribute_changed(*p.spheres_, p.spheres_radius_.get());

			update_spheres_color(p);
			points_provider_->emit_attribute_changed(*p.spheres_, p.spheres_color_.get());

			parallel_foreach_cell(*p.samples_mesh_, [&](PVertex v) -> bool {
				uint32 v_index = index_of(*p.samples_mesh_, v);
				PVertex sphere = (*p.samples_sphere_)[v_index];
				if (sphere.is_valid())
				{
					Vec4 c = value<Vec4>(*p.spheres_, p.spheres_cluster_color_, sphere);
					c[3] = p.spheres_transparency_;
					(*p.samples_color_)[v_index] = c;
				}
				return true;
			});
			points_provider_->emit_attribute_changed(*p.samples_mesh_, p.samples_color_.get());
			points_provider_->emit_attribute_changed(*p.samples_mesh_, p.samples_normal_color_.get());

			compute_skeleton(p);
		}
		else
		{
			points_provider_->emit_connectivity_changed(*p.spheres_);
			points_provider_->emit_attribute_changed(*p.spheres_, p.spheres_position_.get());
			points_provider_->emit_attribute_changed(*p.spheres_, p.spheres_radius_.get());

			update_spheres_color(p);
			points_provider_->emit_attribute_changed(*p.spheres_, p.spheres_color_.get());

			parallel_foreach_cell(*p.samples_mesh_, [&](PVertex v) -> bool {
				uint32 v_index = index_of(*p.samples_mesh_, v);
				PVertex sphere = (*p.samples_sphere_)[v_index];
				if (sphere.is_valid())
				{
					Vec4 c = value<Vec4>(*p.spheres_, p.spheres_cluster_color_, sphere);
					c[3] = p.spheres_transparency_;
					(*p.samples_color_)[v_index] = c;
				}
				return true;
			});
			points_provider_->emit_attribute_changed(*p.samples_mesh_, p.samples_color_.get());
			points_provider_->emit_attribute_changed(*p.samples_mesh_, p.samples_normal_color_.get());

			compute_skeleton(p);
		}

		non_manifold_provider_->emit_connectivity_changed(*p.skeleton_);
		non_manifold_provider_->emit_attribute_changed(*p.skeleton_, p.skeleton_position_.get());
	}

	void start_spheres_update(PointsParameters& p)
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

				std::cout << "Iteration: " << p.iteration_count_ << " | Spheres: " << p.nb_spheres_
						  << " | Error: " << p.total_error_ << " | Diff: " << p.total_error_diff_ << std::endl;

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

	void stop_spheres_update(PointsParameters& p)
	{
		p.stopping_ = true;
	}

	void key_press_event(View* view, int32 key_code) override
	{
		if (!selected_points_)
			return;
		PointsParameters& p = points_parameters_[selected_points_];

		if (key_code == GLFW_KEY_G && view->control_pressed())
		{
			if (p.running_)
				stop_spheres_update(p);
			return;
		}

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
					picked_sphere_center = (*p.spheres_position_)[index_of(*p.spheres_, picked_sphere_)];
					return true;
				}
				const Vec3& sp = (*p.spheres_position_)[index_of(*p.spheres_, v)];
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

protected:
	void left_panel() override
	{
		if (!points_provider_)
		{
			ImGui::Text("Point Mesh Provider not linked.");
			return;
		}

		// Input Point Cloud Selection
		if (ImGui::BeginCombo("Input Point Cloud",
							  selected_points_ ? points_provider_->mesh_name(*selected_points_).c_str() : "None"))
		{
			points_provider_->foreach_mesh([&](POINTS& m, const std::string& name) {
				bool is_selected = (&m == selected_points_);
				if (ImGui::Selectable(name.c_str(), is_selected))
				{
					selected_points_ = &m;
					init_points_data(*selected_points_);
				}
				if (is_selected)
					ImGui::SetItemDefaultFocus();
			});
			ImGui::EndCombo();
		}

		if (!selected_points_)
			return;
		PointsParameters& p = points_parameters_[selected_points_];

		ImGui::Separator();
		if (ImGui::CollapsingHeader("Neural UDF", ImGuiTreeNodeFlags_DefaultOpen))
		{
			if (p.neural_udf_loaded_)
			{
				ImGui::TextColored(ImVec4(0, 1, 0, 1), "Model loaded: %s", p.neural_udf_model_path_.c_str());

				ImGui::Separator();
				ImGui::Text("Alpha Level Set Sampling Settings");

				if (ImGui::Button("Test forward model"))
				{
					test_batch_forward(p);
				}
				ImGui::SameLine();
				ImGui::TextColored(ImVec4(1, 1, 0, 1), "Test batch evaluation");
			}
		}
		// Sampling
		if (ImGui::CollapsingHeader("Sampling", ImGuiTreeNodeFlags_DefaultOpen))
		{
			ImGui::Checkbox("Preprocess sample points", &p.preprocss_sample_points_);
			if (p.input_mode_ == INPUT_NEURAL_UDF)
			{
				ImGui::InputFloat("Alpha", &p.alpha_, 0.001f, 0.1f, "%.4f");
				ImGui::InputInt("Num Samples", &p.num_alpha_samples_, 1000, 10000);
				ImGui::InputFloat("Grid Cell Size", &p.grid_cell_size_, 0.001f, 0.01f, "%.4f");
				ImGui::InputInt("Batch Size", &p.batch_size_, 256, 1024);
				ImGui::InputInt("Max Iterations", &p.udf_max_iterations_, 1000, 8000);
				ImGui::InputFloat("Tolerance", &p.tol_, 0.0f, 0.0f, "%.6f");

				if (ImGui::Button("Sample UDF"))
				{
					load_alpha_samples_to_mesh(p, static_cast<size_t>(p.num_alpha_samples_));
					p.fitting_data_computed_ = false;
				}
				ImGui::SameLine();
				if (ImGui::Button("Clear Samples"))
				{
					if (p.samples_mesh_)
						points_provider_->clear_mesh(*p.samples_mesh_);
					p.fitting_data_computed_ = false;
					p.samples_jitter_backup_valid_ = false;
					p.samples_position_backup_.clear();
					p.samples_normal_backup_.clear();
				}

				if (p.samples_mesh_)
					ImGui::Text("Number of samples: %zu", nb_cells<PVertex>(*p.samples_mesh_));
			}
			else
			{
				ImGui::InputFloat("Alpha", &p.alpha_, 0.001f, 0.1f, "%.4f");
				ImGui::InputFloat("Sample Radius", &p.sample_radius_, 0.0001f, 0.01f, "%.4f");
				ImGui::SliderInt("Sample Iterations", &p.sample_iterations_, 10, 100);
				ImGui::SliderInt("KNN for Normal", &p.knn_k_, 3, 50);

				if (ImGui::Button("Sample Points"))
					sample_points(p);
				ImGui::SameLine();
				if (ImGui::Button("Clear Samples"))
				{
					if (p.samples_mesh_)
						points_provider_->clear_mesh(*p.samples_mesh_);
					p.fitting_data_computed_ = false;
					p.samples_jitter_backup_valid_ = false;
					p.samples_position_backup_.clear();
					p.samples_normal_backup_.clear();
				}
			}
		}

		bool has_samples = p.samples_mesh_ && nb_cells<PVertex>(*p.samples_mesh_) > 0;

		ImGui::Separator();
		ImGui::Text("Input Point Cloud Noise");
		ImGui::SliderFloat("Input Pos Noise (%)", &p.input_jitter_pos_pct_, 0.0f, 5.0f, "%.3f");
		if (ImGui::Button("Jitter Input Positions"))
		{
			std::lock_guard<std::mutex> lock(p.mutex_);
			apply_input_position_jitter(p);
		}
		ImGui::SameLine();
		if (ImGui::Button("Restore Input Positions"))
		{
			std::lock_guard<std::mutex> lock(p.mutex_);
			restore_input_state(p);
		}

		if (has_samples)
		{
			ImGui::Separator();
			ImGui::Text("Sample Point Noise");
			ImGui::SliderFloat("Sample Pos Noise (%)", &p.samples_jitter_pos_pct_, 0.0f, 5.0f, "%.3f");
			if (ImGui::Button("Jitter Sample Positions"))
			{
				std::lock_guard<std::mutex> lock(p.mutex_);
				apply_samples_position_jitter(p);
			}
			ImGui::SameLine();
			if (ImGui::Button("Restore Sample Positions"))
			{
				std::lock_guard<std::mutex> lock(p.mutex_);
				restore_samples_state(p);
			}

			ImGui::SliderFloat("Sample Normal Sigma", &p.samples_jitter_normal_sigma_, 0.0f, 0.2f, "%.4f");
			if (ImGui::Button("Jitter Sample Normals"))
			{
				std::lock_guard<std::mutex> lock(p.mutex_);
				apply_samples_normal_jitter(p);
			}
			ImGui::SameLine();
			if (ImGui::Button("Restore Sample Normals"))
			{
				std::lock_guard<std::mutex> lock(p.mutex_);
				restore_samples_state(p);
			}
		}

		if (!has_samples)
		{
			ImGui::TextColored(ImVec4(1, 1, 0, 1), "Please sample points first.");
		}
		else
		{
			static uint32 init_max_nb_spheres = 1;
			if (ImGui::Button("Compute Fitting Data"))
			{
				compute_fitting_data(p);
				{
					std::lock_guard<std::mutex> lock(p.mutex_);
					init_spheres(p, init_max_nb_spheres);
				}
				update_render_data(p);
			}

			if (p.fitting_data_computed_)
			{
				// Sphere Fitting
				if (ImGui::CollapsingHeader("Sphere Fitting", ImGuiTreeNodeFlags_DefaultOpen))
				{
					ImGui::SliderFloat("Init dilation factor", &p.init_dilation_factor_, 1.0, 4.0);
					ImGui::InputScalar("Init nb spheres", ImGuiDataType_U32, &init_max_nb_spheres);
					if (ImGui::Button("Init spheres"))
					{
						std::lock_guard<std::mutex> lock(p.mutex_);
						init_spheres(p, init_max_nb_spheres);
						update_render_data(p);
					}

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

					ImGui::RadioButton("Sphere Euclidean", (int*)&p.distance_mode_, SPHERE_EUCLIDEAN_DISTANCE);
					ImGui::SameLine();
					ImGui::RadioButton("Pure Euclidean", (int*)&p.distance_mode_, PURE_EUCLIDEAN_DISTANCE);
					ImGui::SameLine();
					ImGui::RadioButton("Line Quadric", (int*)&p.distance_mode_, LINE_QUADRIC_DISTANCE);

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
					ImGui::SameLine();
					if (ImGui::Button("Power Cluster"))
					{
						if (!p.running_)
						{
							std::lock_guard<std::mutex> lock(p.mutex_);
							compute_power_cluster(p);
							compute_spheres_error(p);
							update_render_data(p);
						}
					}

					if (ImGui::Button("Build Skeleton"))
					{
						if (!p.running_)
						{
							std::lock_guard<std::mutex> lock(p.mutex_);
							compute_skeleton(p);
							update_render_data(p);
							non_manifold_provider_->emit_connectivity_changed(*p.skeleton_);
						}
					}
					ImGui::SameLine();
					if (ImGui::Button("Power Skeleton"))
					{
						if (!p.running_)
						{
							std::lock_guard<std::mutex> lock(p.mutex_);
							compute_skeleton_power(p);
							update_render_data(p);
							non_manifold_provider_->emit_connectivity_changed(*p.skeleton_);
						}
					}
					if (ImGui::Button("Topology fix"))
					{
						skeleton_post_pocessing(p);
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
						compute_clusters(p);
						compute_spheres_error(p);
						if (!p.running_)
							update_render_data(p);
					}

					ImGui::Separator();

					ImGui::Text("Total error: %f", p.total_error_ / p.nb_spheres_);
					ImGui::Text("Min error: %f", p.min_error_);
					ImGui::Text("Max error: %f", p.max_error_);

					ImGui::Separator();

					ImGui::Text("Pick the sphere under the mouse with I, split it with S, delete it with D");
					if (picked_sphere_.is_valid())
					{
						ImGui::Text("Picked sphere:");
						const Vec3& sp = (*p.spheres_position_)[index_of(*p.spheres_, picked_sphere_)];
						ImGui::Text("Center: (%f, %f, %f)", sp[0], sp[1], sp[2]);
						ImGui::Text("Radius: %f", (*p.spheres_radius_)[index_of(*p.spheres_, picked_sphere_)]);
					}
				}
			}
		}
	}

private:
	MeshProvider<SURFACE>* surface_provider_ = nullptr;
	MeshProvider<POINTS>* points_provider_ = nullptr;
	PointCloudRender<POINTS>* pcr_ = nullptr;
	MeshProvider<NONMANIFOLD>* non_manifold_provider_ = nullptr;

	SURFACE* selected_surface_ = nullptr;
	std::unique_ptr<acc::BVHTree<uint32, Vec3>> surface_bvh_;
	std::vector<SFace> surface_bvh_faces_;
	std::vector<SVertex> surface_bvh_vertices_;
	std::vector<Vec3> surface_bvh_vertex_positions_;
	bool surface_bvh_dirty_ = false;

	POINTS* selected_points_ = nullptr;
	std::map<POINTS*, PointsParameters> points_parameters_;
	PVertex picked_sphere_;
	std::shared_ptr<boost::synapse::connection> timer_connection_;
	std::array<std::mutex, 43> spheres_mutex_;

	torch::Device device_ = torch::kCPU;
};

} // namespace ui

} // namespace cgogn

#endif // CGOGN_MODULE_UDF_TRAINING_H_
