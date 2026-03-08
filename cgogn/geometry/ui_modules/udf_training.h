#ifndef CGOGN_MODULE_UDF_TRAINING_H_
#define CGOGN_MODULE_UDF_TRAINING_H_

#include <cgogn/ui/app.h>
#include <cgogn/ui/module.h>

#include <cgogn/core/ui_modules/mesh_provider.h>
#include <cgogn/geometry/algos/area.h>
#include <cgogn/geometry/algos/fitting.h>
#include <cgogn/geometry/algos/length.h>
#include <cgogn/geometry/algos/medial_axis.h>
#include <cgogn/geometry/algos/normal.h>
#include <cgogn/geometry/functions/angle.h>
#include <cgogn/geometry/functions/distance.h>
#include <cgogn/geometry/functions/bounding_box.h>
#include <cgogn/geometry/functions/normal.h>
#include <cgogn/geometry/types/line_quadric.h>
#include <cgogn/geometry/types/quadric.h>
#include <cgogn/geometry/types/ray_level_set_sampler.h>
#include <cgogn/geometry/types/ray_level_set_sampler_traits.h>
#include <cgogn/geometry/types/spherical_quadric.h>
#include <cgogn/geometry/types/spatial_grid.h>
#include <cgogn/geometry/types/neural_field_forward.h>
#include <cgogn/geometry/types/vector_traits.h>
#include <cgogn/geometry/types/fast_winding_number_traits.h>
#include <cgogn/geometry/types/fast_winding_number.h>


#include <cgogn/rendering/ui_modules/point_cloud_render.h>
#include <cgogn/rendering/ui_modules/surface_render.h>

#include <atomic>
#include <Eigen/Dense>
#include <Eigen/Eigenvalues>
#include <Eigen/Sparse>

#include <libacc/bvh_tree.h>
#include <libacc/bvh_tree_spheres.h>
#include <libacc/kd_tree.h>

// import CGAL
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/poisson_eliminate.h>

#include <GLFW/glfw3.h>
#include <algorithm>
#include <array>
#include <cmath>
#include <filesystem>
#include <limits>
#include <mutex>
#include <numeric>
#include <queue>
#include <random>
#include <set>
#include <thread>
#include <torch/autograd.h>
#include <torch/script.h>
#include <torch/torch.h>
#include <type_traits>
#include <unordered_map>
#include <unordered_set>

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
using geometry::NeuralFieldForward;
using geometry::Vec3;
using geometry::Vec4;

template <typename SURFACE, typename POINTS, typename NONMANIFOLD, typename RaySamplerTag>
class UDFTraining : public ViewModule
{
public:	
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
	template <typename T>
	using SAttribute = typename mesh_traits<SURFACE>::template Attribute<T>;
	using K = CGAL::Exact_predicates_inexact_constructions_kernel;
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
		LINE_QUADRIC_DISTANCE,
		LINE_QUADRIC_DISTANCE_FREE_RADIUS
	};
	enum InputMode : uint32
	{
		INPUT_POINT_CLOUD,
		INPUT_SURFACE_MESH,
		INPUT_NEURAL_UDF
	};
	enum NeuralModelType : uint32
	{
		NEURAL_MODEL_UDF,
		NEURAL_MODEL_MF
	};

	enum ConnectivityMode : uint32
	{
		MODE_A_MIXED_KNN = 0,
		MODE_B_INSIDE_POWER = 1
	};
	enum InitialMAMode : uint32
	{
		INITIAL_MA_AUTO,
		INITIAL_MA_DISPLACEMENT,
		INITIAL_MA_SHRINKING_BALL
	};

private:
	struct PointsParameters;

	using RaySamplerConfig = geometry::RaySamplerTraits<RaySamplerTag, SURFACE, POINTS>;
	using RaySamplerTraitsType = typename RaySamplerConfig::SamplerTraits;
	using RaySampler = geometry::RayLevelSetSampler<RaySamplerTraitsType>;
	using RaySamplerParams = typename RaySampler::Params;


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
		NeuralModelType neural_model_type_ = NEURAL_MODEL_UDF;
		InitialMAMode initial_ma_mode_override_ = INITIAL_MA_AUTO;
		bool udf_input_normalized_ = false;
		const void* udf_normalized_source_ = nullptr;

		// Ray Sampling
		std::unique_ptr<RaySampler> ray_sampler_;

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
		std::shared_ptr<PAttribute<Vec4>> samples_knn_color_ = nullptr;
		bool show_knn_hover_ = false;
		bool knn_hover_locked_ = false;
		PVertex hovered_sample_;
		std::vector<uint32> hovered_color_indices_;
		std::vector<Vec4> hovered_color_backup_;
		bool filter_positive_projection_ = false;
		std::vector<uint8_t> last_projection_keep_mask_;
		std::vector<Vec3> samples_position_backup_;
		std::vector<Vec3> samples_normal_backup_;
		bool samples_jitter_backup_valid_ = false;
		float32 samples_jitter_pos_pct_ = 0.5f;
		float32 samples_jitter_normal_sigma_ = 0.02f;

		// Alpha-inside samples (inside space, plus cached projection to alpha level set)
		POINTS* alpha_inside_mesh_ = nullptr;
		std::shared_ptr<PAttribute<Vec3>> alpha_inside_position_ = nullptr;
		std::shared_ptr<PAttribute<Vec4>> alpha_inside_color_ = nullptr;
		std::shared_ptr<PAttribute<PVertex>> alpha_inside_sphere_ = nullptr;
		std::shared_ptr<PAttribute<std::vector<PVertex>>> alpha_inside_knn_ = nullptr;
		std::shared_ptr<PAttribute<Vec3>> alpha_inside_projected_position_ = nullptr;
		std::shared_ptr<PAttribute<Vec3>> alpha_inside_projected_normal_ = nullptr;
		acc::KDTree<3, uint32>* alpha_inside_kdtree_ = nullptr;
		std::vector<PVertex> alpha_inside_kdtree_vertices_;

		std::unique_ptr<acc::BVHTreeSpheres<uint32, Vec3>> samples_wn_bvh_;
		std::vector<Vec3> samples_wn_bvh_centers_;
		std::vector<Scalar> samples_wn_bvh_radii_;
		std::unique_ptr<PointFWN<3>> samples_winding_number_;
		std::vector<PVertex> samples_wn_vertices_;
		std::unique_ptr<acc::BVHTreeSpheres<uint32, Vec3>> opt_wn_bvh_;
		std::vector<Vec3> opt_wn_bvh_centers_;
		std::vector<Scalar> opt_wn_bvh_radii_;
		std::unique_ptr<PointFWN<3>> opt_winding_number_;
		std::vector<PVertex> opt_wn_vertices_;
		Scalar beta_ = Scalar(2.0);

		acc::KDTree<3, uint32>* samples_kdtree_ = nullptr; // KDTree of alpha-expanding samples
		std::vector<PVertex> samples_kdtree_vertices_;	   // Vertices of alpha-expanding samples in KDTree order
		acc::KDTree<3, uint32>* samples_ma_kdtree_ = nullptr; // KDTree of sample medial-axis positions
		std::vector<PVertex> samples_ma_kdtree_vertices_;	   // Vertices in MA KDTree order
		acc::KDTree<3, uint32>* opt_kdtree_ = nullptr;		 // KDTree of mixed optimization samples
		std::vector<PVertex> opt_kdtree_vertices_;			   // Vertices in mixed KDTree order
		acc::KDTree<3, uint32>* opt_ma_kdtree_ = nullptr;	   // KDTree of mixed medial-axis positions
		std::vector<PVertex> opt_ma_kdtree_vertices_;		   // Vertices in mixed MA KDTree order

		// Mixed optimization cloud (surface samples + projected inside samples)
		POINTS* opt_mesh_ = nullptr;
		std::shared_ptr<PAttribute<Vec3>> opt_position_ = nullptr;
		std::shared_ptr<PAttribute<Vec3>> opt_normal_ = nullptr;
		std::shared_ptr<PAttribute<Scalar>> opt_area_ = nullptr;
		std::shared_ptr<PAttribute<std::vector<PVertex>>> opt_knn_ = nullptr;
		std::shared_ptr<PAttribute<Spherical_Quadric>> opt_quadric_ = nullptr;
		std::shared_ptr<PAttribute<Line_Quadric>> opt_line_quadric_ = nullptr;
		std::shared_ptr<PAttribute<PVertex>> opt_sphere_ = nullptr;
		std::shared_ptr<PAttribute<Scalar>> opt_error_ = nullptr;
		std::shared_ptr<PAttribute<Vec4>> opt_color_ = nullptr;
		std::shared_ptr<PAttribute<Vec4>> opt_normal_color_ = nullptr;
		std::shared_ptr<PAttribute<Vec3>> opt_ma_position_ = nullptr;
		std::shared_ptr<PAttribute<Scalar>> opt_ma_radius_ = nullptr;
		std::shared_ptr<PAttribute<PVertex>> opt_ma_secondary_vertex_ = nullptr;

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
		std::shared_ptr<PAttribute<Vec4>> spheres_topup_zero_inside_color_ = nullptr;
		std::shared_ptr<PAttribute<std::set<PVertex>>> spheres_neighbor_clusters_ = nullptr;
		std::shared_ptr<PAttribute<PVertex>> spheres_parent_ = nullptr;
		std::shared_ptr<PAttribute<bool>> spheres_do_not_split_ = nullptr;
		std::shared_ptr<PAttribute<Scalar>> spheres_error_ = nullptr;
		std::shared_ptr<PAttribute<Scalar>> spheres_error_not_normalized_ = nullptr;

		// Skeleton
		NONMANIFOLD* skeleton_ = nullptr;
		std::shared_ptr<NMAttribute<Vec3>> skeleton_position_ = nullptr;
		std::shared_ptr<NMAttribute<Scalar>> skeleton_radius_ = nullptr;
		std::shared_ptr<NMAttribute<std::set<std::size_t>>> incident_tets_ = nullptr;

		std::shared_ptr<NMAttribute<Vec3>> skeleton_face_color_ = nullptr;
		std::shared_ptr<NMAttribute<Vec3>> skeleton_face_udf_color_ = nullptr;
		std::shared_ptr<NMAttribute<Vec3>> skeleton_face_boundary_tet_color_ = nullptr;
		std::shared_ptr<NMAttribute<Vec3>> skeleton_edge_color_ = nullptr;
		std::shared_ptr<NMAttribute<Vec3>> skeleton_edge_udf_score_color_ = nullptr;
		std::shared_ptr<NMAttribute<Vec3>> skeleton_edge_boundary_tet_color_ = nullptr;
		std::shared_ptr<NMAttribute<Vec3>> skeleton_edge_non_manifold_color_ = nullptr;
		std::shared_ptr<NMAttribute<uint32>> edge_degree_ = nullptr;

		std::map<NMFaceKey, NMFace> skeleton_faces_map_;

		std::unordered_map<std::size_t, Tet> skeleton_tets_;
		bool topology_stage_snapshot_valid_ = false;
		NONMANIFOLD* topology_stage_snapshot_mesh_ = nullptr;
		std::map<NMFaceKey, NMFace> topology_stage_snapshot_faces_map_;
		std::unordered_map<std::size_t, Tet> topology_stage_snapshot_tets_;

		float32 filter_radius_threshold_ = 0.0f;

		bool sphere_correction_ = false;
		CorrectionMode sphere_correction_mode_ = CORRECT_ALWAYS;
		bool lock_skeleton_connectivity_ = false;
		DistanceMode distance_mode_ = LINE_QUADRIC_DISTANCE;
		bool use_local_clusters_ = false;
		uint32 local_cluster_connectivity_refresh_interval_ = 10;
		ConnectivityMode connectivity_mode_ = MODE_A_MIXED_KNN;
		bool auto_stop_ = false;
		bool auto_split_ = false;
		AutoSplitMode auto_split_mode_ = ERROR_THRESHOLD;
		float32 auto_split_error_threshold_ = 0.00025f;
		uint32 auto_split_max_nb_spheres_ = 500;
		float32 auto_split_ratio_ = 0.2f;
		uint32 auto_split_max_per_iter_error_ = 10;
		uint32 auto_split_max_per_iter_max_ = 100;
		bool error_as_spheres_color_ = false;
		float32 spheres_transparency_ = 0.5f;
		float32 sqem_update_lambda_ = 0.20f;
		float32 sqem_clustering_lambda_ = 0.20f;
		float32 sqem_fix_radius_scale_ = 1.0f;
		bool udf_center_enabled_ = false;
		float32 udf_center_lambda_ = 0.10f;

		// Filtering
		float32 target_radius_ = 0.1f;
		float32 radius_tolerance_ = 0.01f;
		float32 skeleton_edge_udf_zero_tol_ = 1e-4f;
		float32 skeleton_face_udf_zero_tol_ = 1e-4f;
		bool topology_edge_stage_diffuse_ = true;
		bool skeleton_face_score_normalize_by_area_ = false;
		float32 init_dilation_constant_ = 0.001f;

		// Sampling Parameters
		float alpha_ = 0.005f;
		float sample_radius_ = 0.0025f;
		int sample_iterations_ = 30; // Max attempts per point
		int knn_k_ = 10;
		int alpha_inside_knn_k_ = 10;
		int seed_ = 42;
		int cluster_min_points_ = 20;
		float grid_cell_size_ = 0.0025f;
		// Neural UDF Sampling
		int num_alpha_samples_ = 200000;
		int num_alpha_inside_samples_ = 200000;
		int poisson_eliminate_target_samples_ = 50000;
		int batch_size_ = 1310640;	   // NN evaluation batch
		int ray_sampler_batch_size_ = 4096; // Rays per sampling iteration
		float tol_ = 1e-5f; // convergence tolerance
		float alpha_inside_poisson_radius_ = 0.0005f;
		int alpha_inside_min_points_per_sphere_ = 100;
		float alpha_inside_topup_bbox_expand_ = 0.001f;
		bool topup_inside_before_reassign_in_loop_ = true;
		bool inside_skeleton_build_topup_first_ = false;
		std::vector<PVertex> spheres_pending_delete_;
		uint32 pending_delete_count_ = 0;
		bool force_full_refresh_on_stop_ = false;
		std::string stop_reason_;
		uint32 last_surface_count_ = 0;
		uint32 last_inside_count_ = 0;
		uint32 last_mixed_count_ = 0;
		bool mixed_cloud_dirty_ = true;
		std::string last_update_stage_ = "idle";

		// Neural UDF ray sampling parameters
		float udf_bbox_expand_ = 0.05f;
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
		bool manual_stop_requested_ = false;
		bool pending_full_refresh_after_auto_stop_ = false;
		bool slow_down_ = true;
		uint32 update_rate_ = 20;

		~PointsParameters()
		{
			if (samples_kdtree_)
				delete samples_kdtree_;
			if (samples_ma_kdtree_)
				delete samples_ma_kdtree_;
			if (opt_kdtree_)
				delete opt_kdtree_;
			if (opt_ma_kdtree_)
				delete opt_ma_kdtree_;
			if (alpha_inside_kdtree_)
				delete alpha_inside_kdtree_;
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
		if (selected_points_)
		{
			PointsParameters& p = points_parameters_[selected_points_];
			p.input_mode_ = ray_sampler_input_mode();
		}
	} // Compatibility

	void set_selected_points(POINTS& p)
	{
		selected_points_ = &p;
		init_points_data(p);
		PointsParameters& params = points_parameters_[selected_points_];
		params.input_mode_ = ray_sampler_input_mode();
	}

	void load_neural_udf_model(POINTS& points, const std::string& model_path, NeuralModelType model_type)
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
			p.neural_model_type_ = model_type;
			apply_sqem_defaults_by_model_type(p);
			p.udf_input_normalized_ = false;
			p.udf_normalized_source_ = nullptr;
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

	void apply_sqem_defaults_by_model_type(PointsParameters& p)
	{
		if (p.neural_model_type_ == NEURAL_MODEL_MF)
		{
			p.sqem_fix_radius_scale_ = 2.0f;
			p.sqem_update_lambda_ = 2.0f;
			p.sqem_clustering_lambda_ = 2.0f;
		}
		else
		{
			p.sqem_fix_radius_scale_ = 1.0f;
			p.sqem_update_lambda_ = 0.2f;
			p.sqem_clustering_lambda_ = 0.2f;
		}
	}
	
	template <typename Tag = RaySamplerTag>
	void load_alpha_samples_to_mesh(PointsParameters& p, size_t num_points)
	{
		load_alpha_samples_to_mesh_impl(p, num_points, Tag{});
	}

	RaySamplerParams make_ray_params(const PointsParameters& p) const
	{
		RaySamplerParams params;
		params.bbox_expand = p.udf_bbox_expand_;
		params.alpha = p.alpha_;
		params.tol = p.tol_;
		params.step_bound = Scalar(2.0);
		params.batch_size = std::max(1, p.ray_sampler_batch_size_);
		params.max_iterations = p.udf_max_iterations_;
		params.max_outer_iterations = 10;
		params.seed = p.seed_;
		return params;
	}

	std::pair<Vec3, Vec3> model_domain_bbox(const PointsParameters& p) const
	{
		// Projection domain should follow the neural model's native normalization.
		if (p.input_mode_ == INPUT_NEURAL_UDF)
		{
			switch (p.neural_model_type_)
			{
			case NEURAL_MODEL_UDF:
				return {Vec3(-0.5, -0.5, -0.5), Vec3(0.5, 0.5, 0.5)};
			case NEURAL_MODEL_MF:
			default:
				return {Vec3(0, 0, 0), Vec3(1, 1, 1)};
			}
		}
		return {Vec3(0, 0, 0), Vec3(1, 1, 1)};
	}

	std::pair<Vec3, Vec3> compute_sampling_bbox(PointsParameters& p) const
	{
		Vec3 bbox_min(0, 0, 0);
		Vec3 bbox_max(1, 1, 1);
		if (p.input_mode_ == INPUT_NEURAL_UDF)
		{
			if (p.neural_model_type_ == NEURAL_MODEL_UDF)
			{
				bbox_min = Vec3(-0.5, -0.5, -0.5);
				bbox_max = Vec3(0.5, 0.5, 0.5);
			}
			else
			{
				bbox_min = Vec3(0, 0, 0);
				bbox_max = Vec3(1, 1, 1);
			}
		}
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
		const Scalar expand = Scalar(0.1);
		bbox_min -= Vec3(expand, expand, expand);
		bbox_max += Vec3(expand, expand, expand);
		return {bbox_min, bbox_max};
	}

	void normalize_input_for_udf_model(PointsParameters& p)
	{
		const void* source = nullptr;
		if (p.points_ && p.position_)
			source = p.points_;
		else if (selected_surface_)
			source = selected_surface_;
		else
			return;

		if (p.udf_input_normalized_ && p.udf_normalized_source_ == source)
			return;

		if (p.points_ && p.position_ && source == p.points_)
		{
			geometry::normalize_centered(*p.position_);
			rebuild_input_kdtree(p);
			compute_input_normals(p);
			points_provider_->emit_attribute_changed(*p.points_, p.position_.get());
			points_provider_->emit_attribute_changed(*p.points_, p.normal_.get());
			points_provider_->set_mesh_bb_vertex_position(*p.points_, p.position_);
			invalidate_samples_after_input_change(p);
		}
		else if (selected_surface_ && surface_provider_)
		{
			auto s_pos = get_attribute<Vec3, SVertex>(*selected_surface_, "position");
			if (s_pos)
			{
				geometry::normalize_centered(*s_pos.get());
				surface_bvh_dirty_ = true;
				surface_provider_->set_mesh_bb_vertex_position(*selected_surface_, s_pos);
				surface_provider_->emit_attribute_changed(*selected_surface_, s_pos.get());
			}
		}

		p.udf_input_normalized_ = true;
		p.udf_normalized_source_ = source;
	}

	void build_alpha_inside_kdtree(PointsParameters& p)
	{
		if (p.alpha_inside_kdtree_)
		{
			delete p.alpha_inside_kdtree_;
			p.alpha_inside_kdtree_ = nullptr;
		}
		p.alpha_inside_kdtree_vertices_.clear();
		if (!p.alpha_inside_mesh_ || !p.alpha_inside_position_)
			return;
		const uint32 n = nb_cells<PVertex>(*p.alpha_inside_mesh_);
		if (n == 0)
			return;
		std::vector<Vec3> points;
		points.reserve(n);
		p.alpha_inside_kdtree_vertices_.reserve(n);
		foreach_cell(*p.alpha_inside_mesh_, [&](PVertex v) -> bool {
			const uint32 vid = index_of(*p.alpha_inside_mesh_, v);
			points.push_back((*p.alpha_inside_position_)[vid]);
			p.alpha_inside_kdtree_vertices_.push_back(v);
			return true;
		});
		p.alpha_inside_kdtree_ = new acc::KDTree<3, uint32>(points);
	}

	void compute_alpha_inside_knn_graph(PointsParameters& p)
	{
		if (!p.alpha_inside_mesh_ || !p.alpha_inside_kdtree_ || !p.alpha_inside_position_ || !p.alpha_inside_knn_)
			return;
		const int k = std::max(1, p.alpha_inside_knn_k_);
		parallel_foreach_cell(*p.alpha_inside_mesh_, [&](PVertex v) -> bool {
			const uint32 vid = index_of(*p.alpha_inside_mesh_, v);
			const Vec3& pt = (*p.alpha_inside_position_)[vid];
			std::vector<std::pair<uint32, Scalar>> knn_res;
			p.alpha_inside_kdtree_->find_nns(pt, k + 1, &knn_res);
			auto& knn = (*p.alpha_inside_knn_)[vid];
			knn.clear();
			knn.reserve(static_cast<size_t>(k));
			for (const auto& r : knn_res)
			{
				PVertex nb = p.alpha_inside_kdtree_vertices_[r.first];
				if (nb == v)
					continue;
				knn.push_back(nb);
				if (static_cast<int>(knn.size()) >= k)
					break;
			}
			return true;
		});
	}

	void mark_mixed_optimization_cloud_dirty(PointsParameters& p)
	{
		p.mixed_cloud_dirty_ = true;
	}

	void ensure_mixed_optimization_cloud_ready(PointsParameters& p)
	{
		if (!p.mixed_cloud_dirty_)
			return;
		rebuild_mixed_optimization_cloud(p);
	}
	void clear_alpha_inside_samples(PointsParameters& p)
	{
		mark_mixed_optimization_cloud_dirty(p);
		if (p.alpha_inside_mesh_)
			points_provider_->clear_mesh(*p.alpha_inside_mesh_);
		if (p.alpha_inside_kdtree_)
		{
			delete p.alpha_inside_kdtree_;
			p.alpha_inside_kdtree_ = nullptr;
		}
		p.alpha_inside_kdtree_vertices_.clear();
		p.last_inside_count_ = 0;
		if (p.alpha_inside_mesh_)
			points_provider_->emit_connectivity_changed(*p.alpha_inside_mesh_);
	}

	void gather_alpha_inside_state(PointsParameters& p, std::vector<Vec3>& inside_points, std::vector<Vec3>& projected,
								   std::vector<Vec3>& projected_normals, std::vector<PVertex>& owners)
	{
		inside_points.clear();
		projected.clear();
		projected_normals.clear();
		owners.clear();
		if (!p.alpha_inside_mesh_ || !p.alpha_inside_position_)
			return;
		const uint32 n = nb_cells<PVertex>(*p.alpha_inside_mesh_);
		inside_points.reserve(n);
		projected.reserve(n);
		projected_normals.reserve(n);
		owners.reserve(n);
		foreach_cell(*p.alpha_inside_mesh_, [&](PVertex v) -> bool {
			const uint32 vid = index_of(*p.alpha_inside_mesh_, v);
			inside_points.push_back((*p.alpha_inside_position_)[vid]);
			projected.push_back(p.alpha_inside_projected_position_ ? (*p.alpha_inside_projected_position_)[vid]
																  : (*p.alpha_inside_position_)[vid]);
			projected_normals.push_back(p.alpha_inside_projected_normal_ ? (*p.alpha_inside_projected_normal_)[vid]
																		: Vec3(0, 0, 1));
			owners.push_back(p.alpha_inside_sphere_ ? (*p.alpha_inside_sphere_)[vid] : PVertex());
			return true;
		});
	}

	void write_alpha_inside_samples(PointsParameters& p, const std::vector<Vec3>& inside_points,
									const std::vector<Vec3>& projected_points,
									const std::vector<Vec3>& projected_normals,
									const std::vector<PVertex>* inside_owners = nullptr)
	{
		if (!p.alpha_inside_mesh_ || !p.alpha_inside_position_ || !p.alpha_inside_projected_position_ ||
			!p.alpha_inside_projected_normal_)
			return;
		mark_mixed_optimization_cloud_dirty(p);
		points_provider_->clear_mesh(*p.alpha_inside_mesh_);
		const size_t n = inside_points.size();
		for (size_t i = 0; i < n; ++i)
		{
			PVertex v = add_vertex(*p.alpha_inside_mesh_);
			const uint32 vid = index_of(*p.alpha_inside_mesh_, v);
			(*p.alpha_inside_position_)[vid] = inside_points[i];
			(*p.alpha_inside_projected_position_)[vid] =
				(i < projected_points.size()) ? projected_points[i] : inside_points[i];
			Vec3 nrm = (i < projected_normals.size()) ? projected_normals[i] : Vec3(0, 0, 1);
			if (!nrm.allFinite() || nrm.squaredNorm() < Scalar(1e-12))
				nrm = Vec3(0, 0, 1);
			else
				nrm.normalize();
			(*p.alpha_inside_projected_normal_)[vid] = nrm;

			PVertex owner;
			if (inside_owners && i < inside_owners->size())
				owner = (*inside_owners)[i];
			if (p.alpha_inside_sphere_)
				(*p.alpha_inside_sphere_)[vid] = owner;
			if (p.alpha_inside_knn_)
				(*p.alpha_inside_knn_)[vid].clear();
			if (p.alpha_inside_color_)
			{
				Vec4 c(0.15f, 0.55f, 0.95f, 1.0f);
				if (owner.is_valid() && p.spheres_ && p.spheres_cluster_color_)
				{
					const uint32 sid = index_of(*p.spheres_, owner);
					if (sid != INVALID_INDEX && sid < nb_cells<PVertex>(*p.spheres_))
						c = (*p.spheres_cluster_color_)[sid];
				}
				(*p.alpha_inside_color_)[vid] = c;
			}
		}
		build_alpha_inside_kdtree(p);
		compute_alpha_inside_knn_graph(p);
		p.last_inside_count_ = static_cast<uint32>(n);
		points_provider_->emit_connectivity_changed(*p.alpha_inside_mesh_);
		points_provider_->emit_attribute_changed(*p.alpha_inside_mesh_, p.alpha_inside_position_.get());
		points_provider_->emit_attribute_changed(*p.alpha_inside_mesh_, p.alpha_inside_projected_position_.get());
		points_provider_->emit_attribute_changed(*p.alpha_inside_mesh_, p.alpha_inside_projected_normal_.get());
		if (p.alpha_inside_color_)
			points_provider_->emit_attribute_changed(*p.alpha_inside_mesh_, p.alpha_inside_color_.get());
	}

	void filter_alpha_inside_points_by_input_distance(PointsParameters& p, std::vector<Vec3>& inside_points,
													  std::vector<Vec3>* projected_points = nullptr,
													  std::vector<Vec3>* projected_normals = nullptr,
													  std::vector<PVertex>* owners = nullptr)
	{
		if (inside_points.empty())
			return;
		const Scalar max_dist = std::max<Scalar>(Scalar(0), Scalar(p.alpha_));
		if (max_dist <= Scalar(0))
			return;
		std::vector<uint8_t> keep(inside_points.size(), uint8_t(1));
		if (p.input_kdtree_)
		{
			for (size_t i = 0; i < inside_points.size(); ++i)
			{
				std::pair<uint32, Scalar> nn;
				if (!p.input_kdtree_->find_nn(inside_points[i], &nn, max_dist))
					keep[i] = uint8_t(0);
			}
		}
		else
		{
			build_surface_bvh();
			if (!surface_bvh_)
				return;
			for (size_t i = 0; i < inside_points.size(); ++i)
			{
				std::pair<uint32, Vec3> cp;
				if (!surface_bvh_->closest_point(inside_points[i], &cp, max_dist))
					keep[i] = uint8_t(0);
			}
		}

		auto compact_vec3 = [&](std::vector<Vec3>& v) {
			std::vector<Vec3> out;
			out.reserve(v.size());
			for (size_t i = 0; i < v.size() && i < keep.size(); ++i)
				if (keep[i])
					out.push_back(v[i]);
			v.swap(out);
		};
		compact_vec3(inside_points);
		if (projected_points)
			compact_vec3(*projected_points);
		if (projected_normals)
			compact_vec3(*projected_normals);
		if (owners)
		{
			std::vector<PVertex> out;
			out.reserve(owners->size());
			for (size_t i = 0; i < owners->size() && i < keep.size(); ++i)
				if (keep[i])
					out.push_back((*owners)[i]);
			owners->swap(out);
		}
	}

	bool filter_alpha_inside_projected_by_model(PointsParameters& p, std::vector<Vec3>& inside_points,
												std::vector<Vec3>& projected_points,
												std::vector<Vec3>& projected_normals,
												std::vector<PVertex>* owners = nullptr)
	{
		if (inside_points.empty())
			return true;
		if (projected_points.size() != inside_points.size() || projected_normals.size() != inside_points.size())
		{
			std::cerr << "[AlphaInsideProjectionFilter] skipped: size mismatch inside=" << inside_points.size()
					  << " projected=" << projected_points.size() << " normals=" << projected_normals.size() << std::endl;
			return false;
		}
		if (p.input_mode_ != INPUT_NEURAL_UDF)
			return true;
		if (!p.neural_udf_loaded_)
		{
			std::cerr << "[AlphaInsideProjectionFilter] skipped: neural model not loaded." << std::endl;
			return false;
		}

		NeuralFieldForward udf = make_neural_field_forward(p);
		if (!udf.is_loaded())
		{
			std::cerr << "[AlphaInsideProjectionFilter] skipped: neural model not ready." << std::endl;
			return false;
		}

		const size_t total = projected_points.size();
		const size_t batch = std::max<size_t>(1, static_cast<size_t>(p.batch_size_));
		const Scalar proj_tol = Scalar(1e-5);
		std::vector<uint8_t> keep(total, uint8_t(1));
		uint32 reject_value = 0;
		uint32 reject_sdf = 0;

		for (size_t offset = 0; offset < total; offset += batch)
		{
			const size_t count = std::min(batch, total - offset);
			auto cpu_opts = torch::TensorOptions().dtype(torch::kFloat32).device(torch::kCPU);
			if (device_.is_cuda())
				cpu_opts = cpu_opts.pinned_memory(true);
			torch::Tensor pts_cpu = torch::empty({static_cast<int64_t>(count), 3}, cpu_opts);
			auto pts_acc = pts_cpu.accessor<float, 2>();
			for (size_t i = 0; i < count; ++i)
			{
				const Vec3& pnt = projected_points[offset + i];
				pts_acc[static_cast<long>(i)][0] = static_cast<float>(pnt.x());
				pts_acc[static_cast<long>(i)][1] = static_cast<float>(pnt.y());
				pts_acc[static_cast<long>(i)][2] = static_cast<float>(pnt.z());
			}

			if (p.neural_model_type_ == NEURAL_MODEL_MF)
			{
				auto [values_t, sdf_t] = udf.forward_values_sdf_gpu(pts_cpu);
				if (!values_t.defined() || !sdf_t.defined())
				{
					std::cerr << "[AlphaInsideProjectionFilter][MF] forward_values_sdf_gpu failed." << std::endl;
					return false;
				}
				if (values_t.dim() == 2 && values_t.size(1) == 1)
					values_t = values_t.squeeze(1);
				if (sdf_t.dim() == 2 && sdf_t.size(1) == 1)
					sdf_t = sdf_t.squeeze(1);
				torch::Tensor values_cpu = values_t.to(torch::kCPU).contiguous();
				torch::Tensor sdf_cpu = sdf_t.to(torch::kCPU).contiguous();
				auto values_acc = values_cpu.accessor<float, 1>();
				auto sdf_acc = sdf_cpu.accessor<float, 1>();
				for (size_t i = 0; i < count; ++i)
				{
					const bool value_ok = static_cast<Scalar>(values_acc[static_cast<long>(i)]) <= p.alpha_;
					const bool sdf_ok = static_cast<Scalar>(sdf_acc[static_cast<long>(i)]) < Scalar(0);
					if (!value_ok || !sdf_ok)
					{
						keep[offset + i] = uint8_t(0);
						if (!value_ok)
							++reject_value;
						if (!sdf_ok)
							++reject_sdf;
					}
				}
			}
			else
			{
				torch::Tensor values_t = udf.forward_values_gpu(pts_cpu);
				if (!values_t.defined())
				{
					std::cerr << "[AlphaInsideProjectionFilter][UDF] forward_values_gpu failed." << std::endl;
					return false;
				}
				if (values_t.dim() == 2 && values_t.size(1) == 1)
					values_t = values_t.squeeze(1);
				torch::Tensor values_cpu = values_t.to(torch::kCPU).contiguous();
				auto values_acc = values_cpu.accessor<float, 1>();
				for (size_t i = 0; i < count; ++i)
				{
					const Scalar v = static_cast<Scalar>(values_acc[static_cast<long>(i)]);
					if (std::abs(v - p.alpha_) > proj_tol)
					{
						keep[offset + i] = uint8_t(0);
						++reject_value;
					}
				}
			}
		}

		const size_t before = total;
		auto compact_vec3 = [&](std::vector<Vec3>& vec) {
			std::vector<Vec3> out;
			out.reserve(vec.size());
			for (size_t i = 0; i < vec.size() && i < keep.size(); ++i)
				if (keep[i])
					out.push_back(vec[i]);
			vec.swap(out);
		};
		compact_vec3(inside_points);
		compact_vec3(projected_points);
		compact_vec3(projected_normals);
		if (owners)
		{
			std::vector<PVertex> out;
			out.reserve(owners->size());
			for (size_t i = 0; i < owners->size() && i < keep.size(); ++i)
				if (keep[i])
					out.push_back((*owners)[i]);
			owners->swap(out);
		}
		const size_t kept = projected_points.size();
		if (p.neural_model_type_ == NEURAL_MODEL_MF)
		{
			std::cout << "[AlphaInsideProjectionFilter][MF] " << before << " -> " << kept
					  << " kept, reject_value=" << reject_value << " reject_sdf=" << reject_sdf << std::endl;
		}
		else
		{
			std::cout << "[AlphaInsideProjectionFilter][UDF] " << before << " -> " << kept
					  << " kept, reject_abs_udf=" << reject_value << " tol=" << proj_tol << std::endl;
		}
		return true;
	}

	struct PowerSphereSnapshot
	{
		std::vector<PVertex> spheres;
		std::vector<uint32> sphere_indices;
		std::vector<Vec3> centers;
		std::vector<Scalar> radii;
	};

	bool build_power_sphere_snapshot(PointsParameters& p, PowerSphereSnapshot& snap)
	{
		snap.spheres.clear();
		snap.sphere_indices.clear();
		snap.centers.clear();
		snap.radii.clear();
		if (!p.spheres_ || !p.spheres_position_ || !p.spheres_radius_)
			return false;
		foreach_cell(*p.spheres_, [&](PVertex s) -> bool {
			const uint32 sid = index_of(*p.spheres_, s);
			snap.spheres.push_back(s);
			snap.sphere_indices.push_back(sid);
			snap.centers.push_back((*p.spheres_position_)[sid]);
			snap.radii.push_back((*p.spheres_radius_)[sid]);
			return true;
		});
		return !snap.spheres.empty();
	}

	int closest_power_sphere(const PowerSphereSnapshot& snap, const Vec3& pnt)
	{
		if (snap.spheres.empty())
			return -1;
		Scalar min_pd = std::numeric_limits<Scalar>::max();
		int best = -1;
		for (size_t i = 0; i < snap.spheres.size(); ++i)
		{
			const Vec3 d = pnt - snap.centers[i];
			const Scalar pd = d.squaredNorm() - snap.radii[i] * snap.radii[i];
			if (pd < min_pd)
			{
				min_pd = pd;
				best = static_cast<int>(i);
			}
		}
		return best;
	}

	void assign_alpha_inside_sphere_power(PointsParameters& p)
	{
		if (!p.alpha_inside_mesh_ || !p.alpha_inside_position_ || !p.alpha_inside_sphere_)
			return;
		PowerSphereSnapshot snap;
		if (!build_power_sphere_snapshot(p, snap))
		{
			p.alpha_inside_sphere_->fill(PVertex());
			return;
		}
		parallel_foreach_cell(*p.alpha_inside_mesh_, [&](PVertex v) -> bool {
			const uint32 vid = index_of(*p.alpha_inside_mesh_, v);
			const int best = closest_power_sphere(snap, (*p.alpha_inside_position_)[vid]);
			PVertex owner = (best >= 0) ? snap.spheres[static_cast<size_t>(best)] : PVertex();
			(*p.alpha_inside_sphere_)[vid] = owner;
			if (p.alpha_inside_color_)
			{
				Vec4 c(0.15f, 0.55f, 0.95f, 1.0f);
				if (owner.is_valid() && p.spheres_cluster_color_)
				{
					const uint32 sid = index_of(*p.spheres_, owner);
					if (sid != INVALID_INDEX && sid < nb_cells<PVertex>(*p.spheres_))
						c = (*p.spheres_cluster_color_)[sid];
				}
				(*p.alpha_inside_color_)[vid] = c;
			}
			return true;
		});
		if (p.alpha_inside_color_)
			points_provider_->emit_attribute_changed(*p.alpha_inside_mesh_, p.alpha_inside_color_.get());
	}

	void set_alpha_inside_samples_with_projection(PointsParameters& p, std::vector<Vec3> inside_points)
	{
		filter_alpha_inside_points_by_input_distance(p, inside_points);
		std::vector<Vec3> projected;
		std::vector<Vec3> normals;
		if (!inside_points.empty())
		{
			auto proj = project_points_to_alpha_gpu(p, inside_points);
			projected = std::move(proj.first);
			normals = std::move(proj.second);
			if (projected.size() != inside_points.size() || normals.size() != inside_points.size())
			{
				projected.resize(inside_points.size());
				normals.resize(inside_points.size(), Vec3(0, 0, 1));
				for (size_t i = 0; i < inside_points.size(); ++i)
					if (!projected[i].allFinite())
						projected[i] = inside_points[i];
			}
			filter_alpha_inside_projected_by_model(p, inside_points, projected, normals, nullptr);
		}
		write_alpha_inside_samples(p, inside_points, projected, normals, nullptr);
		if (p.spheres_ && nb_cells<PVertex>(*p.spheres_) > 0)
			assign_alpha_inside_sphere_power(p);
	}

	void ensure_alpha_inside_min_points_per_sphere(PointsParameters& p, std::vector<Vec3>& inside_points,
												   std::vector<PVertex>& inside_owners)
	{
		const int min_points = std::max(0, p.alpha_inside_min_points_per_sphere_);
		if (min_points <= 0)
		{
			std::cout << "[AlphaInsideTopup] skipped: min_per_sphere<=0" << std::endl;
			return;
		}
		PowerSphereSnapshot snap;
		if (!build_power_sphere_snapshot(p, snap))
		{
			std::cout << "[AlphaInsideTopup] aborted: no valid spheres." << std::endl;
			return;
		}
		if (snap.spheres.empty())
			return;
		std::unordered_map<uint32, std::size_t> sphere_to_snap;
		sphere_to_snap.reserve(snap.sphere_indices.size());
		for (std::size_t i = 0; i < snap.sphere_indices.size(); ++i)
			sphere_to_snap.emplace(snap.sphere_indices[i], i);

		inside_owners.resize(inside_points.size());
		for (std::size_t i = 0; i < inside_points.size(); ++i)
		{
			const int owner = closest_power_sphere(snap, inside_points[i]);
			inside_owners[i] = (owner >= 0) ? snap.spheres[static_cast<std::size_t>(owner)] : PVertex();
		}

		auto build_cluster_boxes = [&](std::vector<Vec3>& box_min, std::vector<Vec3>& box_max,
									   std::vector<int>& cluster_sizes, Scalar expand, int& valid_boxes) {
			box_min.assign(snap.spheres.size(), Vec3(0, 0, 0));
			box_max.assign(snap.spheres.size(), Vec3(0, 0, 0));
			cluster_sizes.assign(snap.spheres.size(), 0);
			std::vector<bool> has_box(snap.spheres.size(), false);
			valid_boxes = 0;

			for (std::size_t pid = 0; pid < inside_points.size(); ++pid)
			{
				std::size_t i = static_cast<std::size_t>(-1);
				if (pid < inside_owners.size() && inside_owners[pid].is_valid() && p.spheres_)
				{
					const uint32 sid = index_of(*p.spheres_, inside_owners[pid]);
					auto it = sphere_to_snap.find(sid);
					if (it != sphere_to_snap.end())
						i = it->second;
				}
				else
				{
					const int owner = closest_power_sphere(snap, inside_points[pid]);
					if (owner >= 0)
						i = static_cast<std::size_t>(owner);
				}
				if (i == static_cast<std::size_t>(-1))
					continue;
				const Vec3& pt = inside_points[pid];
				if (!has_box[i])
				{
					box_min[i] = pt;
					box_max[i] = pt;
					has_box[i] = true;
				}
				else
				{
					box_min[i] = box_min[i].cwiseMin(pt);
					box_max[i] = box_max[i].cwiseMax(pt);
				}
				++cluster_sizes[i];
			}

			if (p.samples_mesh_ && p.samples_position_)
			{
				foreach_cell(*p.samples_mesh_, [&](PVertex v) -> bool {
					const uint32 vid = index_of(*p.samples_mesh_, v);
					if (vid == INVALID_INDEX)
						return true;
					const Vec3& pos = (*p.samples_position_)[vid];
					const int owner = closest_power_sphere(snap, pos);
					if (owner < 0)
						return true;
					const size_t i = static_cast<size_t>(owner);
					if (!has_box[i])
					{
						box_min[i] = pos;
						box_max[i] = pos;
						has_box[i] = true;
					}
					else
					{
						box_min[i] = box_min[i].cwiseMin(pos);
						box_max[i] = box_max[i].cwiseMax(pos);
					}
					++cluster_sizes[i];
					return true;
				});
			}

			for (size_t i = 0; i < snap.spheres.size(); ++i)
			{
				if (!has_box[i])
					continue;
				const Vec3 e(expand, expand, expand);
				box_min[i] = box_min[i] - e;
				box_max[i] = box_max[i] + e;
				++valid_boxes;
			}
		};

		std::vector<int> owned_counts(snap.spheres.size(), 0);
		for (std::size_t pid = 0; pid < inside_points.size(); ++pid)
		{
			if (pid < inside_owners.size() && inside_owners[pid].is_valid() && p.spheres_)
			{
				const uint32 sid = index_of(*p.spheres_, inside_owners[pid]);
				auto it = sphere_to_snap.find(sid);
				if (it != sphere_to_snap.end())
					++owned_counts[it->second];
				continue;
			}
			const int owner = closest_power_sphere(snap, inside_points[pid]);
			if (owner >= 0)
				++owned_counts[static_cast<size_t>(owner)];
		}
		std::vector<int> deficits(snap.spheres.size(), 0);
		std::vector<int> initial_deficits(snap.spheres.size(), 0);
		for (size_t i = 0; i < snap.spheres.size(); ++i)
		{
			const int need = std::max(0, min_points - owned_counts[i]);
			deficits[i] = need;
			initial_deficits[i] = need;
		}

		int total_need = 0;
		for (int d : deficits)
			if (d > 0)
				total_need += d;
		if (total_need <= 0)
		{
			std::cout << "[AlphaInsideTopup] already satisfied: spheres=" << snap.spheres.size()
					  << " min_per_sphere=" << min_points << " inside_points=" << inside_points.size() << std::endl;
			return;
		}

		const Scalar poisson_r = (p.alpha_inside_poisson_radius_ > Scalar(0)) ? Scalar(p.alpha_inside_poisson_radius_)
																				: Scalar(p.grid_cell_size_);
		const Scalar expand = (p.alpha_inside_topup_bbox_expand_ > 0.0f) ? Scalar(p.alpha_inside_topup_bbox_expand_)
																		 : std::max(poisson_r, Scalar(p.grid_cell_size_));
		std::vector<Vec3> cluster_box_min;
		std::vector<Vec3> cluster_box_max;
		std::vector<int> cluster_sizes;
		int valid_cluster_boxes = 0;
		build_cluster_boxes(cluster_box_min, cluster_box_max, cluster_sizes, expand, valid_cluster_boxes);
		std::vector<int> raw_cluster_sizes = cluster_sizes;
		std::vector<uint8_t> zero_start_sphere(snap.spheres.size(), uint8_t(0));
		std::vector<uint8_t> zero_start_bbox_fallback(snap.spheres.size(), uint8_t(0));
		auto sphere_bbox_with_expand = [&](size_t i) {
			const Scalar base_r = std::max(Scalar(0), snap.radii[i]);
			const Scalar half = std::max(base_r, std::max(poisson_r, Scalar(p.grid_cell_size_))) + expand;
			const Vec3 e(half, half, half);
			return std::make_pair(snap.centers[i] - e, snap.centers[i] + e);
		};
		for (size_t i = 0; i < snap.spheres.size(); ++i)
		{
			if (owned_counts[i] != 0 || deficits[i] <= 0)
				continue;
			zero_start_sphere[i] = uint8_t(1);
			zero_start_bbox_fallback[i] = uint8_t(1);
			auto bb = sphere_bbox_with_expand(i);
			cluster_box_min[i] = bb.first;
			cluster_box_max[i] = bb.second;
			if (cluster_sizes[i] <= 0)
				cluster_sizes[i] = 1;
		}
		if (valid_cluster_boxes == 0)
		{
			std::cout << "[AlphaInsideTopup] aborted: no valid cluster bbox found. "
					  << "Please compute clusters first." << std::endl;
			return;
		}

		const bool use_udf_value_filter =
			(p.input_mode_ == INPUT_NEURAL_UDF && p.neural_model_type_ == NEURAL_MODEL_UDF);
		const bool skip_input_filter_for_mf =
			(p.input_mode_ == INPUT_NEURAL_UDF && p.neural_model_type_ == NEURAL_MODEL_MF);
		if (use_udf_value_filter && !p.neural_udf_loaded_)
		{
			std::cout << "[AlphaInsideTopup] aborted: neural UDF model not loaded for udf(x)<alpha input filter."
					  << std::endl;
			return;
		}
		const Scalar max_input_dist = std::max(Scalar(0), Scalar(p.alpha_));
		const bool has_input_kdtree = (p.input_kdtree_ != nullptr && !p.input_kdtree_vertices_.empty());
		if (!use_udf_value_filter && !skip_input_filter_for_mf && !has_input_kdtree)
		{
			build_surface_bvh();
			if (!surface_bvh_)
			{
				std::cout << "[AlphaInsideTopup] aborted: no input distance source." << std::endl;
				return;
			}
		}

		auto pass_input_filter = [&](const Vec3& pos) -> bool {
			if (use_udf_value_filter)
			{
				Scalar udf_value = Scalar(0);
				Vec3 grad(0, 0, 0);
				if (!eval_udf_and_grad(p, pos, udf_value, grad))
					return false;
				return udf_value < p.alpha_;
			}
			if (skip_input_filter_for_mf)
				return true;
			if (has_input_kdtree)
			{
				std::pair<uint32, Scalar> nn;
				return p.input_kdtree_->find_nn(pos, &nn, max_input_dist);
			}
			std::pair<uint32, Vec3> cp;
			return surface_bvh_->closest_point(pos, &cp, max_input_dist);
		};

		const char* input_filter_desc =
			use_udf_value_filter ? "udf(x)<alpha"
								 : (skip_input_filter_for_mf ? "deferred(MF-final-filter)" : "surface_distance<=alpha");
		std::cout << "[AlphaInsideTopup] begin spheres=" << snap.spheres.size() << " min_per_sphere=" << min_points
				  << " poisson_r=" << poisson_r << " bbox_expand=" << expand
				  << " initial_inside_points=" << inside_points.size() << " valid_cluster_boxes=" << valid_cluster_boxes
				  << " total_need=" << total_need
				  << " input_filter=" << input_filter_desc << std::endl;
		const int udf_input_filter_batch_size =
			std::max(1, p.batch_size_);
		constexpr std::size_t kTopupLogStride = 200;
		for (size_t i = 0; i < snap.spheres.size(); ++i)
		{
			if ((i % kTopupLogStride) != 0)
				continue;
			std::cout << "[AlphaInsideTopup][SphereInit] sphere_rank=" << i
					  << " sphere_mesh_index=" << snap.sphere_indices[i]
					  << " existing=" << owned_counts[i] << " need=" << initial_deficits[i]
					  << " cluster_points_raw=" << raw_cluster_sizes[i]
					  << " cluster_points_effective=" << cluster_sizes[i]
					  << " bbox_source=" << (zero_start_bbox_fallback[i] ? "sphere_fallback" : "cluster_box")
					  << " has_bbox=" << (cluster_sizes[i] > 0 ? "yes" : "no") << std::endl;
        }

		SpatialGrid inside_grid(poisson_r);
		for (size_t i = 0; i < inside_points.size(); ++i)
			inside_grid.insert(inside_points[i], static_cast<uint32>(i));
		std::mt19937 gen(static_cast<uint32>(p.seed_) + 24681357u);

		int total_added = 0;
		std::vector<int> attempts_per_sphere(snap.spheres.size(), 0);
		std::vector<int> accepted_per_sphere(snap.spheres.size(), 0);
		std::vector<int> zero_attempts_per_sphere(snap.spheres.size(), 0);
		std::vector<int> zero_reject_input_per_sphere(snap.spheres.size(), 0);
		std::vector<int> zero_reject_owner_per_sphere(snap.spheres.size(), 0);
		std::vector<int> zero_reject_poisson_per_sphere(snap.spheres.size(), 0);
		std::vector<int> zero_added_per_sphere(snap.spheres.size(), 0);

		int zero_start_count = 0;
		for (uint8_t f : zero_start_sphere)
			zero_start_count += (f != 0) ? 1 : 0;
		if (zero_start_count > 0)
		{
			const int max_zero_prepass_attempts = 2000000;
			const int max_zero_per_sphere_attempts = 200000;
			int zero_attempted = 0;
			int zero_added = 0;
			int zero_reject_input = 0;
			int zero_reject_owner = 0;
			int zero_reject_poisson = 0;
			int zero_unresolved = 0;
			std::vector<Vec3> zero_batch_points;
			std::vector<size_t> zero_batch_spheres;
			zero_batch_points.reserve(static_cast<size_t>(udf_input_filter_batch_size));
			zero_batch_spheres.reserve(static_cast<size_t>(udf_input_filter_batch_size));

			auto flush_zero_batch = [&]() {
				if (zero_batch_points.empty())
					return;
				std::vector<Scalar> udf_values;
				const bool batch_ok = use_udf_value_filter ? eval_udf_values(p, zero_batch_points, udf_values) : true;
				for (size_t bi = 0; bi < zero_batch_points.size() && total_need > 0; ++bi)
				{
					const size_t si = zero_batch_spheres[bi];
					if (si >= deficits.size() || deficits[si] <= 0)
						continue;
					const Vec3& pt = zero_batch_points[bi];
					bool pass = false;
					if (use_udf_value_filter)
					{
						if (batch_ok && bi < udf_values.size())
						{
							const Scalar v = udf_values[bi];
							pass = std::isfinite(v) && (v < p.alpha_);
						}
						else
						{
							pass = pass_input_filter(pt);
						}
					}
					else
					{
						pass = pass_input_filter(pt);
					}
					if (!pass)
					{
						++zero_reject_input;
						++zero_reject_input_per_sphere[si];
						continue;
					}
					const int owner = closest_power_sphere(snap, pt);
					if (owner != static_cast<int>(si))
					{
						++zero_reject_owner;
						++zero_reject_owner_per_sphere[si];
						continue;
					}
					if (!inside_grid.is_valid_sample(pt, poisson_r, inside_points))
					{
						++zero_reject_poisson;
						++zero_reject_poisson_per_sphere[si];
						continue;
					}
					inside_points.push_back(pt);
					inside_owners.push_back(snap.spheres[si]);
					inside_grid.insert(pt, static_cast<uint32>(inside_points.size() - 1));
					--deficits[si];
					--total_need;
					++zero_added;
					++total_added;
					++accepted_per_sphere[si];
					++zero_added_per_sphere[si];
				}
				zero_batch_points.clear();
				zero_batch_spheres.clear();
			};

			for (size_t i = 0; i < snap.spheres.size() && total_need > 0; ++i)
			{
				if (!zero_start_sphere[i] || deficits[i] <= 0)
					continue;
				const int per_sphere_budget =
					std::min(max_zero_per_sphere_attempts, std::max(5000, deficits[i] * 2000));
				int attempted_i = 0;
				std::uniform_real_distribution<Scalar> ux(cluster_box_min[i].x(), cluster_box_max[i].x());
				std::uniform_real_distribution<Scalar> uy(cluster_box_min[i].y(), cluster_box_max[i].y());
				std::uniform_real_distribution<Scalar> uz(cluster_box_min[i].z(), cluster_box_max[i].z());
				while (attempted_i < per_sphere_budget && deficits[i] > 0 && total_need > 0 &&
					   zero_attempted < max_zero_prepass_attempts)
				{
					Vec3 pt(ux(gen), uy(gen), uz(gen));
					++attempted_i;
					++zero_attempted;
					++attempts_per_sphere[i];
					++zero_attempts_per_sphere[i];
					if ((zero_attempted % 200000) == 0)
					{
						std::cout << "[AlphaInsideTopup][ZeroStartPrepass] progress attempted=" << zero_attempted
								  << " added=" << zero_added << " remaining_need=" << total_need << std::endl;
					}
					if (use_udf_value_filter)
					{
						zero_batch_points.push_back(pt);
						zero_batch_spheres.push_back(i);
						if (static_cast<int>(zero_batch_points.size()) >= udf_input_filter_batch_size)
							flush_zero_batch();
						continue;
					}
					if (!pass_input_filter(pt))
					{
						++zero_reject_input;
						++zero_reject_input_per_sphere[i];
						continue;
					}
					const int owner = closest_power_sphere(snap, pt);
					if (owner != static_cast<int>(i))
					{
						++zero_reject_owner;
						++zero_reject_owner_per_sphere[i];
						continue;
					}
					if (!inside_grid.is_valid_sample(pt, poisson_r, inside_points))
					{
						++zero_reject_poisson;
						++zero_reject_poisson_per_sphere[i];
						continue;
					}
					inside_points.push_back(pt);
					inside_owners.push_back(snap.spheres[i]);
					inside_grid.insert(pt, static_cast<uint32>(inside_points.size() - 1));
					--deficits[i];
					--total_need;
					++zero_added;
					++total_added;
					++accepted_per_sphere[i];
					++zero_added_per_sphere[i];
				}
				if (deficits[i] > 0)
					++zero_unresolved;
				if (zero_attempted >= max_zero_prepass_attempts)
					break;
			}
			flush_zero_batch();
			zero_unresolved = 0;
			for (size_t i = 0; i < snap.spheres.size(); ++i)
				if (zero_start_sphere[i] && deficits[i] > 0)
					++zero_unresolved;
			std::cout << "[AlphaInsideTopup][ZeroStartPrepass] spheres=" << zero_start_count
					  << " attempted=" << zero_attempted << " added=" << zero_added
					  << " reject_input=" << zero_reject_input
					  << " reject_owner=" << zero_reject_owner
					  << " reject_poisson=" << zero_reject_poisson
					  << " unresolved_spheres=" << zero_unresolved
					  << " remaining_need=" << total_need << std::endl;
			for (size_t i = 0; i < snap.spheres.size(); ++i)
			{
				if (!zero_start_sphere[i])
					continue;
				const bool unresolved = deficits[i] > 0;
				const bool failed_all = zero_added_per_sphere[i] == 0;
				if (!unresolved && !failed_all)
					continue;
				const Vec3 bbox_extent = cluster_box_max[i] - cluster_box_min[i];
				const Scalar bbox_diag = bbox_extent.norm();
				const Scalar bbox_volume = bbox_extent.x() * bbox_extent.y() * bbox_extent.z();
				const int post_input_candidates = std::max(0, zero_attempts_per_sphere[i] - zero_reject_input_per_sphere[i]);
				const int owner_pass = std::max(0, post_input_candidates - zero_reject_owner_per_sphere[i]);
				const bool prepass_reached = zero_attempts_per_sphere[i] > 0;
				const bool skipped_due_budget = (!prepass_reached && zero_attempted >= max_zero_prepass_attempts);
				std::cout << "[AlphaInsideTopup][ZeroStartPrepass][Sphere] sphere_rank=" << i
						  << " sphere_mesh_index=" << snap.sphere_indices[i]
						  << " existing_start=" << owned_counts[i]
						  << " need_start=" << initial_deficits[i]
						  << " cluster_points_raw=" << raw_cluster_sizes[i]
						  << " cluster_points_effective=" << cluster_sizes[i]
						  << " bbox_source=" << (zero_start_bbox_fallback[i] ? "sphere_fallback" : "cluster_box")
						  << " prepass_reached=" << (prepass_reached ? "yes" : "no")
						  << " skipped_due_budget=" << (skipped_due_budget ? "yes" : "no")
						  << " attempted=" << zero_attempts_per_sphere[i]
						  << " added=" << zero_added_per_sphere[i]
						  << " reject_input=" << zero_reject_input_per_sphere[i]
						  << " reject_owner=" << zero_reject_owner_per_sphere[i]
						  << " reject_poisson=" << zero_reject_poisson_per_sphere[i]
						  << " post_input_candidates=" << post_input_candidates
						  << " owner_pass=" << owner_pass
						  << " need_after_prepass=" << deficits[i] << std::endl;
			}
		}
		int stalled_passes = 0;
		const int max_passes = 10;
		const int max_attempts_per_pass = 1000000;
		for (int pass = 0; pass < max_passes && total_need > 0; ++pass)
		{
			std::vector<size_t> active_spheres;
			active_spheres.reserve(snap.spheres.size());
			int skipped_no_bbox = 0;
			for (size_t i = 0; i < snap.spheres.size(); ++i)
			{
				if (deficits[i] <= 0)
					continue;
				if (cluster_sizes[i] <= 0)
				{
					++skipped_no_bbox;
					continue;
				}
				active_spheres.push_back(i);
			}
			if (active_spheres.empty())
			{
				std::cout << "[AlphaInsideTopup][Pass " << (pass + 1)
						  << "] no active spheres (all missing bbox or already satisfied), remaining_need=" << total_need
						  << " skipped_no_bbox=" << skipped_no_bbox << std::endl;
				break;
			}
			std::shuffle(active_spheres.begin(), active_spheres.end(), gen);

			int added_this_pass = 0;
			int attempted_this_pass = 0;
			int rejected_input_pass = 0;
			int rejected_owner_pass = 0;
			int rejected_poisson_pass = 0;
			std::vector<int> attempted_per_sphere_pass(snap.spheres.size(), 0);
			std::vector<int> accepted_per_sphere_pass(snap.spheres.size(), 0);
			std::vector<Vec3> pass_batch_points;
			std::vector<size_t> pass_batch_spheres;
			pass_batch_points.reserve(static_cast<size_t>(udf_input_filter_batch_size));
			pass_batch_spheres.reserve(static_cast<size_t>(udf_input_filter_batch_size));
			auto flush_pass_batch = [&]() {
				if (pass_batch_points.empty())
					return;
				std::vector<Scalar> udf_values;
				const bool batch_ok = use_udf_value_filter ? eval_udf_values(p, pass_batch_points, udf_values) : true;
				for (size_t bi = 0; bi < pass_batch_points.size() && total_need > 0; ++bi)
				{
					const size_t si = pass_batch_spheres[bi];
					if (si >= deficits.size() || deficits[si] <= 0)
						continue;
					const Vec3& pt = pass_batch_points[bi];
					bool pass_ok = false;
					if (use_udf_value_filter)
					{
						if (batch_ok && bi < udf_values.size())
						{
							const Scalar v = udf_values[bi];
							pass_ok = std::isfinite(v) && (v < p.alpha_);
						}
						else
						{
							pass_ok = pass_input_filter(pt);
						}
					}
					else
					{
						pass_ok = pass_input_filter(pt);
					}
					if (!pass_ok)
					{
						++rejected_input_pass;
						continue;
					}
					const int owner = closest_power_sphere(snap, pt);
					if (owner != static_cast<int>(si))
					{
						++rejected_owner_pass;
						continue;
					}
					if (!inside_grid.is_valid_sample(pt, poisson_r, inside_points))
					{
						++rejected_poisson_pass;
						continue;
					}
					inside_points.push_back(pt);
					inside_owners.push_back(snap.spheres[si]);
					inside_grid.insert(pt, static_cast<uint32>(inside_points.size() - 1));
					--deficits[si];
					--total_need;
					++added_this_pass;
					++total_added;
					++accepted_per_sphere_pass[si];
					++accepted_per_sphere[si];
				}
				pass_batch_points.clear();
				pass_batch_spheres.clear();
			};

			const int fair_budget = std::max(16, max_attempts_per_pass / std::max(1, static_cast<int>(active_spheres.size())));
			for (size_t idx = 0; idx < active_spheres.size() && total_need > 0; ++idx)
			{
				const size_t i = active_spheres[idx];
				if (deficits[i] <= 0)
					continue;
				const int pass_budget_left = max_attempts_per_pass - attempted_this_pass;
				if (pass_budget_left <= 0)
					break;

				const int need_i = deficits[i];
				const int base_target = std::max(fair_budget, std::min(need_i * 20, fair_budget * 4));
				const int attempts_target = std::min(pass_budget_left, base_target);
				std::uniform_real_distribution<Scalar> ux(cluster_box_min[i].x(), cluster_box_max[i].x());
				std::uniform_real_distribution<Scalar> uy(cluster_box_min[i].y(), cluster_box_max[i].y());
				std::uniform_real_distribution<Scalar> uz(cluster_box_min[i].z(), cluster_box_max[i].z());

				for (int a = 0; a < attempts_target && deficits[i] > 0 && total_need > 0; ++a)
				{
					Vec3 pt(ux(gen), uy(gen), uz(gen));
					++attempted_this_pass;
					++attempted_per_sphere_pass[i];
					++attempts_per_sphere[i];
					++zero_attempts_per_sphere[i];

					if (use_udf_value_filter)
					{
						pass_batch_points.push_back(pt);
						pass_batch_spheres.push_back(i);
						if (static_cast<int>(pass_batch_points.size()) >= udf_input_filter_batch_size)
							flush_pass_batch();
						continue;
					}

					if (!pass_input_filter(pt))
					{
						++rejected_input_pass;
						continue;
					}

					const int owner = closest_power_sphere(snap, pt);
					if (owner != static_cast<int>(i))
					{
						++rejected_owner_pass;
						continue;
					}

					if (!inside_grid.is_valid_sample(pt, poisson_r, inside_points))
					{
						++rejected_poisson_pass;
						continue;
					}

					inside_points.push_back(pt);
					inside_owners.push_back(snap.spheres[i]);
					inside_grid.insert(pt, static_cast<uint32>(inside_points.size() - 1));
					--deficits[i];
					--total_need;
					++added_this_pass;
					++total_added;
					++accepted_per_sphere_pass[i];
					++accepted_per_sphere[i];
				}
				if (attempted_this_pass >= max_attempts_per_pass)
					continue;
			}
			flush_pass_batch();
			std::cout << "[AlphaInsideTopup][Pass " << (pass + 1) << "] active_spheres=" << active_spheres.size()
					  << " fair_budget=" << fair_budget << " skipped_no_bbox=" << skipped_no_bbox
					  << " attempted=" << attempted_this_pass << " added=" << added_this_pass
					  << " reject_input=" << rejected_input_pass
					  << " reject_owner=" << rejected_owner_pass
					  << " reject_poisson=" << rejected_poisson_pass
					  << " remaining_need=" << total_need << std::endl;

			if (added_this_pass == 0)
			{
				++stalled_passes;
				if (stalled_passes >= 2)
					break;
			}
			else
			{
				stalled_passes = 0;
			}
		}
		int remaining_need = 0;
		for (int d : deficits)
			if (d > 0)
				remaining_need += d;
		for (size_t i = 0; i < snap.spheres.size(); ++i)
		{
			if (deficits[i] <= 0)
				continue;
			std::cout << "[AlphaInsideTopup][SphereFinal] sphere_rank=" << i
					  << " sphere_mesh_index=" << snap.sphere_indices[i]
					  << " existing_start=" << owned_counts[i] << " need_start=" << initial_deficits[i]
                      << " attempted_total=" << attempts_per_sphere[i]
					  << " added_total=" << accepted_per_sphere[i]
					  << " need_final=" << deficits[i] << std::endl;
		}

		std::cout << "[AlphaInsideTopup] min_per_sphere=" << min_points << " added=" << total_added
				  << " remaining_need=" << remaining_need << " final_inside_points=" << inside_points.size()
				  << std::endl;
	}

	void topup_alpha_inside_samples_from_current_mesh(PointsParameters& p, const char* source)
	{
		std::cout << "[AlphaInsideTopupSource] begin source=" << (source ? source : "unknown") << std::endl;
		topup_alpha_inside_samples_from_current_mesh(p);
		std::cout << "[AlphaInsideTopupSource] end source=" << (source ? source : "unknown") << std::endl;
	}

	void topup_alpha_inside_samples_from_current_mesh(PointsParameters& p)
	{
		if (!p.alpha_inside_mesh_ || !p.alpha_inside_position_)
			return;
		if (!p.spheres_ || nb_cells<PVertex>(*p.spheres_) == 0)
		{
			std::cout << "[AlphaInsideTopup] skipped: spheres not initialized." << std::endl;
			return;
		}
		const uint32 inside_before =
			(p.alpha_inside_mesh_ != nullptr) ? static_cast<uint32>(nb_cells<PVertex>(*p.alpha_inside_mesh_)) : uint32(0);
		std::vector<Vec3> inside_points, projected_points, projected_normals;
		std::vector<PVertex> owners;
		assign_alpha_inside_sphere_power(p);
		gather_alpha_inside_state(p, inside_points, projected_points, projected_normals, owners);
		const size_t gathered_count = inside_points.size();
		filter_alpha_inside_points_by_input_distance(p, inside_points, &projected_points, &projected_normals, &owners);
		const size_t filtered_count = inside_points.size();
		const size_t old_size = inside_points.size();
		std::cout << "[AlphaInsideTopup] wrapper begin inside_before_mesh=" << inside_before
				  << " gathered=" << gathered_count << " after_filter=" << filtered_count << std::endl;
		ensure_alpha_inside_min_points_per_sphere(p, inside_points, owners);
		const size_t added_inside = (inside_points.size() > old_size) ? (inside_points.size() - old_size) : size_t(0);
		if (inside_points.size() > old_size)
		{
			std::vector<Vec3> new_points;
			new_points.reserve(inside_points.size() - old_size);
			for (size_t i = old_size; i < inside_points.size(); ++i)
				new_points.push_back(inside_points[i]);
			auto projected_new = project_points_to_alpha_gpu(p, new_points);
			auto& new_proj = projected_new.first;
			auto& new_normals = projected_new.second;
			if (new_proj.size() != new_points.size() || new_normals.size() != new_points.size())
			{
				std::cout << "[AlphaInsideTopup] projection size mismatch: new=" << new_points.size()
						  << " proj=" << new_proj.size() << " normal=" << new_normals.size() << " (fallback resize)"
						  << std::endl;
				new_proj.resize(new_points.size());
				new_normals.resize(new_points.size(), Vec3(0, 0, 1));
				for (size_t i = 0; i < new_points.size(); ++i)
					if (!new_proj[i].allFinite())
						new_proj[i] = new_points[i];
			}
			projected_points.insert(projected_points.end(), new_proj.begin(), new_proj.end());
			projected_normals.insert(projected_normals.end(), new_normals.begin(), new_normals.end());
		}
		filter_alpha_inside_projected_by_model(p, inside_points, projected_points, projected_normals, &owners);
		write_alpha_inside_samples(p, inside_points, projected_points, projected_normals, &owners);
		update_spheres_topup_zero_inside_color(p);
		if (p.spheres_topup_zero_inside_color_)
			points_provider_->emit_attribute_changed(*p.spheres_, p.spheres_topup_zero_inside_color_.get());
		const uint32 inside_after =
			(p.alpha_inside_mesh_ != nullptr) ? static_cast<uint32>(nb_cells<PVertex>(*p.alpha_inside_mesh_)) : uint32(0);
		std::cout << "[AlphaInsideTopup] wrapper end added=" << added_inside << " mesh_inside_after=" << inside_after
				  << " projected_inside=" << projected_points.size() << " owners=" << owners.size() << std::endl;
	}
	void recompute_opt_knn_area_quadrics(PointsParameters& p)
	{
		if (!p.opt_mesh_ || !p.opt_position_ || !p.opt_normal_ || !p.opt_area_ || !p.opt_knn_ || !p.opt_quadric_ ||
			!p.opt_line_quadric_)
			return;
		const uint32 n = nb_cells<PVertex>(*p.opt_mesh_);
		if (n == 0)
			return;
		std::vector<Vec3> points;
		std::vector<PVertex> vertices;
		points.reserve(n);
		vertices.reserve(n);
		foreach_cell(*p.opt_mesh_, [&](PVertex v) -> bool {
			const uint32 vid = index_of(*p.opt_mesh_, v);
			points.push_back((*p.opt_position_)[vid]);
			vertices.push_back(v);
			return true;
		});
		acc::KDTree<3, uint32> kdtree(points);
		parallel_foreach_cell(*p.opt_mesh_, [&](PVertex v) -> bool {
			const uint32 vid = index_of(*p.opt_mesh_, v);
			const Vec3& pt = (*p.opt_position_)[vid];
			std::vector<std::pair<uint32, Scalar>> knn_res;
			kdtree.find_nns(pt, p.knn_k_ + 10, &knn_res);
			(*p.opt_knn_)[vid].clear();
			Scalar sum_dist = 0.0;
			Vec3 nrm = (*p.opt_normal_)[vid];
			if (nrm.squaredNorm() < Scalar(1e-12))
				nrm = Vec3(0, 0, 1);
			else
				nrm.normalize();
			int kept = 0;
			for (const auto& res : knn_res)
			{
				PVertex nb = vertices[res.first];
				if (nb == v)
					continue;
				const uint32 nb_idx = index_of(*p.opt_mesh_, nb);
				const Vec3& q = (*p.opt_position_)[nb_idx];
				if (std::abs((q - pt).dot(nrm)) > p.alpha_)
					continue;
				(*p.opt_knn_)[vid].push_back(nb);
				sum_dist += res.second;
				++kept;
				if (kept >= p.knn_k_ + 1)
					break;
			}
			if (kept == 0)
			{
				for (const auto& res : knn_res)
				{
					PVertex nb = vertices[res.first];
					if (nb == v)
						continue;
					(*p.opt_knn_)[vid].push_back(nb);
					sum_dist += res.second;
				}
			}
			(*p.opt_area_)[vid] = (sum_dist * sum_dist) / (2.0 * p.knn_k_);
			return true;
		});

		parallel_foreach_cell(*p.opt_mesh_, [&](PVertex v) -> bool {
			const uint32 vid = index_of(*p.opt_mesh_, v);
			Spherical_Quadric& q = (*p.opt_quadric_)[vid];
			Line_Quadric& lq = (*p.opt_line_quadric_)[vid];
			q.clear();
			lq.clear();
			const Vec3& pos = (*p.opt_position_)[vid];
			Vec3 nrm = (*p.opt_normal_)[vid];
			if (nrm.squaredNorm() < Scalar(1e-12))
				nrm = Vec3(0, 0, 1);
			else
				nrm.normalize();
			Scalar a = (*p.opt_area_)[vid] / (p.knn_k_ + 1.0);
			q += Spherical_Quadric(Vec4(pos.x(), pos.y(), pos.z(), 0), Vec4(nrm.x(), nrm.y(), nrm.z(), 1)) * a;
			lq += Line_Quadric(pos, nrm) * a;
			for (PVertex vn : (*p.opt_knn_)[vid])
			{
				const uint32 vn_idx = index_of(*p.opt_mesh_, vn);
				const Vec3& pn = (*p.opt_position_)[vn_idx];
				Vec3 nn = (*p.opt_normal_)[vn_idx];
				if (nn.squaredNorm() < Scalar(1e-12))
					nn = Vec3(0, 0, 1);
				else
					nn.normalize();
				Scalar an = (*p.opt_area_)[vn_idx] / (p.knn_k_ + 1.0);
				q += Spherical_Quadric(Vec4(pn.x(), pn.y(), pn.z(), 0), Vec4(nn.x(), nn.y(), nn.z(), 1)) * an;
				lq += Line_Quadric(pn, nn) * an;
			}
			return true;
		});
	}

	void rebuild_mixed_optimization_cloud(PointsParameters& p)
	{
		if (!p.opt_mesh_ || !p.opt_position_ || !p.opt_normal_ || !p.samples_mesh_ || !p.samples_position_ ||
			!p.samples_normal_)
			return;
		if (p.opt_kdtree_)
		{
			delete p.opt_kdtree_;
			p.opt_kdtree_ = nullptr;
		}
		if (p.opt_ma_kdtree_)
		{
			delete p.opt_ma_kdtree_;
			p.opt_ma_kdtree_ = nullptr;
		}
		p.opt_kdtree_vertices_.clear();
		p.opt_ma_kdtree_vertices_.clear();
		p.opt_winding_number_.reset();
		p.opt_wn_bvh_.reset();
		p.opt_wn_vertices_.clear();
		p.opt_wn_bvh_centers_.clear();
		p.opt_wn_bvh_radii_.clear();

		points_provider_->clear_mesh(*p.opt_mesh_);

		uint32 surface_count = 0;
		foreach_cell(*p.samples_mesh_, [&](PVertex v) -> bool {
			const uint32 vid = index_of(*p.samples_mesh_, v);
			PVertex ov = add_vertex(*p.opt_mesh_);
			const uint32 oid = index_of(*p.opt_mesh_, ov);
			(*p.opt_position_)[oid] = (*p.samples_position_)[vid];
			(*p.opt_normal_)[oid] = (*p.samples_normal_)[vid];
			if (p.opt_normal_color_)
				(*p.opt_normal_color_)[oid] =
					Vec4(((*p.opt_normal_)[oid].x() + 1.0) * 0.5, ((*p.opt_normal_)[oid].y() + 1.0) * 0.5,
						 (((*p.opt_normal_)[oid].z() + 1.0) * 0.5), 1.0);
			if (p.opt_color_)
				(*p.opt_color_)[oid] = Vec4(0.2f, 0.8f, 0.2f, 1.0f);
			++surface_count;
			return true;
		});

		uint32 inside_count = 0;
		if (p.alpha_inside_mesh_ && p.alpha_inside_projected_position_ && p.alpha_inside_projected_normal_)
		{
			foreach_cell(*p.alpha_inside_mesh_, [&](PVertex v) -> bool {
				const uint32 vid = index_of(*p.alpha_inside_mesh_, v);
				PVertex ov = add_vertex(*p.opt_mesh_);
				const uint32 oid = index_of(*p.opt_mesh_, ov);
				(*p.opt_position_)[oid] = (*p.alpha_inside_projected_position_)[vid];
				Vec3 nrm = (*p.alpha_inside_projected_normal_)[vid];
				if (!nrm.allFinite() || nrm.squaredNorm() < Scalar(1e-12))
					nrm = Vec3(0, 0, 1);
				else
					nrm.normalize();
				(*p.opt_normal_)[oid] = nrm;
				if (p.opt_normal_color_)
					(*p.opt_normal_color_)[oid] = Vec4((nrm.x() + 1.0) * 0.5, (nrm.y() + 1.0) * 0.5,
									 (nrm.z() + 1.0) * 0.5, 1.0);
				if (p.opt_color_)
					(*p.opt_color_)[oid] = Vec4(0.95f, 0.45f, 0.15f, 1.0f);
				++inside_count;
				return true;
			});
		}

		recompute_opt_knn_area_quadrics(p);
		p.mixed_cloud_dirty_ = false;
		p.last_surface_count_ = surface_count;
		p.last_inside_count_ = inside_count;
		p.last_mixed_count_ = surface_count + inside_count;
		points_provider_->set_mesh_bb_vertex_position(*p.opt_mesh_, p.opt_position_);
		points_provider_->emit_connectivity_changed(*p.opt_mesh_);
		points_provider_->emit_attribute_changed(*p.opt_mesh_, p.opt_position_.get());
		points_provider_->emit_attribute_changed(*p.opt_mesh_, p.opt_normal_.get());
		if (p.opt_normal_color_)
			points_provider_->emit_attribute_changed(*p.opt_mesh_, p.opt_normal_color_.get());
		if (p.opt_color_)
			points_provider_->emit_attribute_changed(*p.opt_mesh_, p.opt_color_.get());
	}

	void load_alpha_samples_to_mesh_impl(PointsParameters& p, size_t num_points, geometry::RaySamplerNeural)
	{
		if (!p.neural_udf_loaded_)
		{
			std::cerr << "Neural UDF model not loaded. Cannot sample alpha level set." << std::endl;
			return;
		}

		const bool use_spatial_grid = (p.grid_cell_size_ > Scalar(0));
		if (use_spatial_grid)
			p.samples_spatial_grid_ = std::make_unique<SpatialGrid>(p.grid_cell_size_);
		else
		{
			p.samples_spatial_grid_.reset();
			std::cout << "Grid Cell Size <= 0: sampling without SpatialGrid deduplication." << std::endl;
		}
		std::cout << "Sampling " << num_points << " points on alpha=" << p.alpha_ << " level set..." << std::endl;

		RaySamplerParams ray_params = make_ray_params(p);
		NeuralFieldForward udf = make_neural_field_forward(p);
		auto [bbox_min, bbox_max] = compute_sampling_bbox(p);
		if (!p.ray_sampler_)
			p.ray_sampler_ = std::make_unique<RaySampler>(ray_params, udf.device());
		else
		{
			p.ray_sampler_->set_params(ray_params);
			p.ray_sampler_->set_device(udf.device());
		}
		std::cout << "Neural sampling bbox: min(" << bbox_min.transpose() << "), max(" << bbox_max.transpose() << ")"
				  << std::endl;
		bool used_sdf_filter = false;
		auto traits = RaySamplerConfig::make(udf);
		const size_t inside_target =
			p.num_alpha_inside_samples_ > 0 ? static_cast<size_t>(p.num_alpha_inside_samples_) : size_t(0);
		SpatialGrid inside_grid(
			(p.alpha_inside_poisson_radius_ > Scalar(0)) ? Scalar(p.alpha_inside_poisson_radius_) : Scalar(p.grid_cell_size_));
		auto sample_result = p.ray_sampler_->sample_alpha_level_set_rays_with_inside(
			traits, num_points, inside_target, p.samples_spatial_grid_.get(), &inside_grid, p.grid_cell_size_, bbox_min,
			bbox_max, &used_sdf_filter);
		std::vector<Vec3> sampled_points = std::move(sample_result.surface_samples);
		std::vector<Vec3> inside_points = std::move(sample_result.inside_samples);
		set_alpha_inside_samples_with_projection(p, std::move(inside_points));

		if (sampled_points.empty())
		{
			std::cerr << "Failed to sample points on alpha level set." << std::endl;
			return;
		}

		std::cout << "Successfully sampled " << sampled_points.size() << " surface points and "
				  << nb_cells<PVertex>(*p.alpha_inside_mesh_) << " inside points." << std::endl;
		{
			const size_t check_n = std::min<size_t>(30, sampled_points.size());
			if (check_n > 0)
			{
				std::vector<size_t> indices(sampled_points.size());
				std::iota(indices.begin(), indices.end(), size_t(0));
				std::mt19937 gen(p.seed_ + 1337);
				std::shuffle(indices.begin(), indices.end(), gen);
				std::vector<Vec3> check_points;
				check_points.reserve(check_n);
				for (size_t i = 0; i < check_n; ++i)
					check_points.push_back(sampled_points[indices[i]]);

				BatchUDFResult check_res = udf.forward_batch(check_points);
				if (check_res.ok && check_res.values.size() == check_n)
				{
					std::cout << "UDF check (random " << check_n << "):" << std::endl;
					for (size_t i = 0; i < check_n; ++i)
					{
						const Scalar v = check_res.values[i];
						std::cout << "  " << i << ": udf=" << v << " |udf-alpha|=" << std::abs(v - p.alpha_)
								  << std::endl;
					}
				}
				else
				{
					std::cerr << "UDF check failed (forward_batch)." << std::endl;
				}
			}
		}
		if (p.samples_mesh_)
			points_provider_->clear_mesh(*p.samples_mesh_);
		p.samples_jitter_backup_valid_ = false;
		p.samples_position_backup_.clear();
		p.samples_normal_backup_.clear();
		for (const Vec3& pt : sampled_points)
		{
			PVertex v = add_vertex(*p.samples_mesh_);
			uint32 v_idx = index_of(*p.samples_mesh_, v);
			(*p.samples_position_)[v_idx] = pt;
			if (p.samples_knn_color_)
				(*p.samples_knn_color_)[v_idx] = Vec4(0.0, 0.0, 0.0, 1.0);
		}

		std::cout << "Computing normals from UDF gradients..." << std::endl;
		std::vector<Vec3> all_positions;
		all_positions.reserve(sampled_points.size());
		foreach_cell(*p.samples_mesh_, [&](PVertex v) {
			uint32 v_idx = index_of(*p.samples_mesh_, v);
			all_positions.push_back((*p.samples_position_)[v_idx]);
			return true;
		});

		bool normals_ok = true;
		std::vector<Vec3> gradients(all_positions.size(), Vec3(0, 0, 1));
		const size_t grad_batch_cap = 131464;
		const size_t grad_batch =
			std::max<size_t>(1, std::min<size_t>(static_cast<size_t>(p.batch_size_), grad_batch_cap));
		for (size_t offset = 0; offset < all_positions.size(); offset += grad_batch)
		{
			const size_t count = std::min(grad_batch, all_positions.size() - offset);
			auto cpu_opts = torch::TensorOptions().dtype(torch::kFloat32).device(torch::kCPU);
			if (device_.is_cuda())
				cpu_opts = cpu_opts.pinned_memory(true);
			torch::Tensor pts_cpu = torch::empty({static_cast<int64_t>(count), 3}, cpu_opts);
			auto pts_acc = pts_cpu.accessor<float, 2>();
			for (size_t i = 0; i < count; ++i)
			{
				const Vec3& pnt = all_positions[offset + i];
				pts_acc[static_cast<long>(i)][0] = static_cast<float>(pnt.x());
				pts_acc[static_cast<long>(i)][1] = static_cast<float>(pnt.y());
				pts_acc[static_cast<long>(i)][2] = static_cast<float>(pnt.z());
			}

			auto [values_t, grad_t] = udf.forward_values_grad_gpu(pts_cpu);
			if (!values_t.defined() || !grad_t.defined())
			{
				normals_ok = false;
				break;
			}
			if (grad_t.dim() != 2 || grad_t.size(0) != static_cast<long>(count) || grad_t.size(1) != 3)
			{
				normals_ok = false;
				break;
			}

			torch::Tensor grad_cpu = grad_t.to(torch::kCPU).contiguous();
			auto grad_acc = grad_cpu.accessor<float, 2>();
			for (size_t i = 0; i < count; ++i)
			{
				Vec3 normal(static_cast<Scalar>(grad_acc[static_cast<long>(i)][0]),
							static_cast<Scalar>(grad_acc[static_cast<long>(i)][1]),
							static_cast<Scalar>(grad_acc[static_cast<long>(i)][2]));
				if (normal.squaredNorm() > Scalar(1e-12))
					normal.normalize();
				else
					normal = Vec3(0, 0, 1);
				gradients[offset + i] = normal;
			}
		}

		if (normals_ok)
		{
			uint32 idx = 0;
			foreach_cell(*p.samples_mesh_, [&](PVertex v) {
				uint32 v_idx = index_of(*p.samples_mesh_, v);
				const Vec3& normal = gradients[idx];
				(*p.samples_normal_)[v_idx] = normal;
				(*p.samples_normal_color_)[v_idx] =
					Vec4((normal.x() + 1.0) * 0.5, (normal.y() + 1.0) * 0.5, (normal.z() + 1.0) * 0.5, 1.0);
				idx++;
				return true;
			});
		}
		if (!normals_ok)
			std::cerr << "Failed to compute normals from UDF gradients." << std::endl;

		std::cout << "Building KDTree for sampled points..." << std::endl;
		build_kdtree(p);
		points_provider_->emit_connectivity_changed(*p.samples_mesh_);
		mark_mixed_optimization_cloud_dirty(p);
		assign_alpha_inside_sphere_power(p);
		ensure_mixed_optimization_cloud_ready(p);
		std::cout << "Alpha level set sampling complete. Ready for fitting." << std::endl;
	}

	void load_alpha_samples_to_mesh_impl(PointsParameters& p, size_t num_points, geometry::RaySamplerSurface)
	{
		build_surface_bvh();
		if (!surface_bvh_)
		{
			std::cerr << "Surface BVH not available. Cannot sample alpha level set." << std::endl;
			return;
		}

		auto s_pos = get_attribute<Vec3, SVertex>(*selected_surface_, "position");
		if (!s_pos)
		{
			std::cerr << "Surface position attribute not available. Cannot sample alpha level set." << std::endl;
			return;
		}

		const bool use_spatial_grid = (p.grid_cell_size_ > Scalar(0));
		if (use_spatial_grid)
			p.samples_spatial_grid_ = std::make_unique<SpatialGrid>(p.grid_cell_size_);
		else
		{
			p.samples_spatial_grid_.reset();
			std::cout << "Grid Cell Size <= 0: sampling without SpatialGrid deduplication." << std::endl;
		}
		std::cout << "Sampling " << num_points << " points on alpha=" << p.alpha_ << " level set..." << std::endl;

		RaySamplerParams ray_params = make_ray_params(p);
		if (!p.ray_sampler_)
			p.ray_sampler_ = std::make_unique<RaySampler>(ray_params, torch::kCPU);
		else
		{
			p.ray_sampler_->set_params(ray_params);
			p.ray_sampler_->set_device(torch::kCPU);
		}
		auto [bbox_min, bbox_max] = compute_sampling_bbox(p);
		std::cout << "Surface sampling bbox: min(" << bbox_min.transpose() << "), max(" << bbox_max.transpose() << ")"
				  << std::endl;

		auto traits =
			RaySamplerConfig::make(*selected_surface_, s_pos.get(), surface_bvh_.get(), &surface_bvh_faces_);
		const size_t inside_target =
			p.num_alpha_inside_samples_ > 0 ? static_cast<size_t>(p.num_alpha_inside_samples_) : size_t(0);
		SpatialGrid inside_grid(
			(p.alpha_inside_poisson_radius_ > Scalar(0)) ? Scalar(p.alpha_inside_poisson_radius_) : Scalar(p.grid_cell_size_));
		auto sample_result = p.ray_sampler_->sample_alpha_level_set_rays_with_inside(
			traits, num_points, inside_target, p.samples_spatial_grid_.get(), &inside_grid, p.grid_cell_size_, bbox_min,
			bbox_max);
		std::vector<Vec3> sampled_points = std::move(sample_result.surface_samples);
		std::vector<Vec3> inside_points = std::move(sample_result.inside_samples);
		set_alpha_inside_samples_with_projection(p, std::move(inside_points));

		if (sampled_points.empty())
		{
			std::cerr << "Failed to sample points on alpha level set." << std::endl;
			return;
		}

		std::cout << "Successfully sampled " << sampled_points.size() << " surface points and "
				  << nb_cells<PVertex>(*p.alpha_inside_mesh_) << " inside points." << std::endl;
		if (p.samples_mesh_)
			points_provider_->clear_mesh(*p.samples_mesh_);
		p.samples_jitter_backup_valid_ = false;
		p.samples_position_backup_.clear();
		p.samples_normal_backup_.clear();
		for (const Vec3& pt : sampled_points)
		{
			PVertex v = add_vertex(*p.samples_mesh_);
			uint32 v_idx = index_of(*p.samples_mesh_, v);
			(*p.samples_position_)[v_idx] = pt;
			if (p.samples_knn_color_)
				(*p.samples_knn_color_)[v_idx] = Vec4(0.0, 0.0, 0.0, 1.0);
		}

		std::cout << "Computing normals from surface mesh..." << std::endl;
		foreach_cell(*p.samples_mesh_, [&](PVertex v) {
			uint32 v_idx = index_of(*p.samples_mesh_, v);
			const Vec3& pos = (*p.samples_position_)[v_idx];
			std::pair<uint32, Vec3> cp;
			surface_bvh_->closest_point(pos, &cp);
			SFace face = surface_bvh_faces_[cp.first];
			Vec3 n = geometry::normal(*selected_surface_, face, s_pos.get());
			if (surface_vertex_normal_)
			{
				std::array<SVertex, 3> vertices;
				uint32 vi = 0;
				foreach_incident_vertex(*selected_surface_, face, [&](SVertex sv) -> bool {
					if (vi < vertices.size())
						vertices[vi++] = sv;
					return true;
				});
				if (vi == vertices.size())
				{
					const Vec3& p0 = value<Vec3>(*selected_surface_, s_pos, vertices[0]);
					const Vec3& p1 = value<Vec3>(*selected_surface_, s_pos, vertices[1]);
					const Vec3& p2 = value<Vec3>(*selected_surface_, s_pos, vertices[2]);
					Scalar u = 0.0, v_bary = 0.0, w = 0.0;
					cgogn::geometry::closest_point_in_triangle(cp.second, p0, p1, p2, u, v_bary, w);
					const Vec3& n0 = value<Vec3>(*selected_surface_, surface_vertex_normal_, vertices[0]);
					const Vec3& n1 = value<Vec3>(*selected_surface_, surface_vertex_normal_, vertices[1]);
					const Vec3& n2 = value<Vec3>(*selected_surface_, surface_vertex_normal_, vertices[2]);
					n = u * n0 + v_bary * n1 + w * n2;
				}
			}
			if (n.squaredNorm() < Scalar(1e-12))
				n = Vec3(0, 0, 1);
			else
				n.normalize();
			Vec3 to_sample = pos - cp.second;
			if (to_sample.squaredNorm() > Scalar(1e-12))
			{
				if (to_sample.dot(n) < Scalar(0))
					n = -n;
			}
			(*p.samples_normal_)[v_idx] = n;
			(*p.samples_normal_color_)[v_idx] =
				Vec4((n.x() + 1.0) * 0.5, (n.y() + 1.0) * 0.5, (n.z() + 1.0) * 0.5, 1.0);
	
			return true;
		});

		std::cout << "Building KDTree for sampled points..." << std::endl;
		build_kdtree(p);
		points_provider_->emit_connectivity_changed(*p.samples_mesh_);
		mark_mixed_optimization_cloud_dirty(p);
		assign_alpha_inside_sphere_power(p);
		ensure_mixed_optimization_cloud_ready(p);
		std::cout << "Alpha level set sampling complete. Ready for fitting." << std::endl;
	}

	void load_alpha_samples_to_mesh_impl(PointsParameters& p, size_t num_points, geometry::RaySamplerPointCloud)
	{
		if (!p.input_kdtree_)
		{
			std::cerr << "Input point cloud KDTree not available. Cannot sample alpha level set." << std::endl;
			return;
		}
		if (!p.normal_ || !p.knn_)
		{
			std::cerr << "Input point cloud normals/KNN not available. Cannot sample alpha level set." << std::endl;
			return;
		}

		const bool use_spatial_grid = (p.grid_cell_size_ > Scalar(0));
		if (use_spatial_grid)
			p.samples_spatial_grid_ = std::make_unique<SpatialGrid>(p.grid_cell_size_);
		else
		{
			p.samples_spatial_grid_.reset();
			std::cout << "Grid Cell Size <= 0: sampling without SpatialGrid deduplication." << std::endl;
		}
		std::cout << "Sampling " << num_points << " points on alpha=" << p.alpha_ << " level set..." << std::endl;

		RaySamplerParams ray_params = make_ray_params(p);
		if (!p.ray_sampler_)
			p.ray_sampler_ = std::make_unique<RaySampler>(ray_params, torch::kCPU);
		else
		{
			p.ray_sampler_->set_params(ray_params);
			p.ray_sampler_->set_device(torch::kCPU);
		}
		auto [bbox_min, bbox_max] = compute_sampling_bbox(p);
		std::cout << "Point cloud sampling bbox: min(" << bbox_min.transpose() << "), max(" << bbox_max.transpose()
				  << ")" << std::endl;

		auto traits = RaySamplerConfig::make(*p.points_, p.position_.get(), p.normal_.get(), p.knn_.get(),
											 p.input_kdtree_, &p.input_kdtree_vertices_);
		const size_t inside_target =
			p.num_alpha_inside_samples_ > 0 ? static_cast<size_t>(p.num_alpha_inside_samples_) : size_t(0);
		SpatialGrid inside_grid(
			(p.alpha_inside_poisson_radius_ > Scalar(0)) ? Scalar(p.alpha_inside_poisson_radius_) : Scalar(p.grid_cell_size_));
		auto sample_result = p.ray_sampler_->sample_alpha_level_set_rays_with_inside(
			traits, num_points, inside_target, p.samples_spatial_grid_.get(), &inside_grid, p.grid_cell_size_, bbox_min,
			bbox_max);
		std::vector<Vec3> sampled_points = std::move(sample_result.surface_samples);
		std::vector<Vec3> inside_points = std::move(sample_result.inside_samples);
		set_alpha_inside_samples_with_projection(p, std::move(inside_points));

		if (sampled_points.empty())
		{
			std::cerr << "Failed to sample points on alpha level set." << std::endl;
			return;
		}

		std::cout << "Successfully sampled " << sampled_points.size() << " surface points and "
				  << nb_cells<PVertex>(*p.alpha_inside_mesh_) << " inside points." << std::endl;
		if (p.samples_mesh_)
			points_provider_->clear_mesh(*p.samples_mesh_);
		p.samples_jitter_backup_valid_ = false;
		p.samples_position_backup_.clear();
		p.samples_normal_backup_.clear();
		for (const Vec3& pt : sampled_points)
		{
			PVertex v = add_vertex(*p.samples_mesh_);
			uint32 v_idx = index_of(*p.samples_mesh_, v);
			(*p.samples_position_)[v_idx] = pt;
			if (p.samples_knn_color_)
				(*p.samples_knn_color_)[v_idx] = Vec4(0.0, 0.0, 0.0, 1.0);
		}

		std::cout << "Computing normals from input point cloud..." << std::endl;
		const Scalar eps = Scalar(1e-12);
		foreach_cell(*p.samples_mesh_, [&](PVertex v) {
			uint32 v_idx = index_of(*p.samples_mesh_, v);
			const Vec3& pos = (*p.samples_position_)[v_idx];
			std::pair<uint32, Scalar> knn_res;
			p.input_kdtree_->find_nn(pos, &knn_res);
			
			uint32 idx = knn_res.first;
			PVertex vn = p.input_kdtree_vertices_[idx];
			uint32 vn_idx = index_of(*p.points_, vn);
			Vec3 n = pos - (*p.position_)[vn_idx];
			n.normalize();
			(*p.samples_normal_)[v_idx] = n;
			(*p.samples_normal_color_)[v_idx] =
				Vec4((n.x() + 1.0) * 0.5, (n.y() + 1.0) * 0.5, (n.z() + 1.0) * 0.5, 1.0);
		
			return true;
		});

		std::cout << "Building KDTree for sampled points..." << std::endl;
		build_kdtree(p);
		points_provider_->emit_connectivity_changed(*p.samples_mesh_);
		mark_mixed_optimization_cloud_dirty(p);
		assign_alpha_inside_sphere_power(p);
		ensure_mixed_optimization_cloud_ready(p);
		std::cout << "Alpha level set sampling complete. Ready for fitting." << std::endl;
	}

	void pre_process_sampling_points_kdtree(PointsParameters& p, std::vector<Vec3>& points)
	{
		const Scalar tol = Scalar(1e-2);
		const Scalar max_dist = p.alpha_ + tol;
		std::cout << "Before pre-processing, " << points.size() << " points." << std::endl;
		if (!p.input_kdtree_)
			return;
		auto new_end = std::remove_if(points.begin(), points.end(), [&](const Vec3& pos) {
			std::pair<uint32, Scalar> knn_res;
			return !p.input_kdtree_->find_nn(pos, &knn_res, max_dist);
		});
		points.erase(new_end, points.end());
		std::cout << "After pre-processing, " << points.size() << " points." << std::endl;
	}

	void pre_process_sampling_points_bvh(PointsParameters& p, std::vector<Vec3>& points)
	{
		const Scalar tol = Scalar(1e-2);
		const Scalar max_dist = p.alpha_ + tol;
		std::cout << "Before pre-processing, " << points.size() << " points." << std::endl;
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

	void recompute_samples_normals_from_current_input(PointsParameters& p)
	{
		if (!p.samples_mesh_ || !p.samples_position_ || !p.samples_normal_)
			return;

		const Scalar eps = Scalar(1e-12);
		bool normals_ok = false;

		if (p.input_mode_ == INPUT_NEURAL_UDF && p.neural_udf_loaded_)
		{
			std::vector<Vec3> all_positions;
			all_positions.reserve(nb_cells<PVertex>(*p.samples_mesh_));
			foreach_cell(*p.samples_mesh_, [&](PVertex v) {
				uint32 v_idx = index_of(*p.samples_mesh_, v);
				all_positions.push_back((*p.samples_position_)[v_idx]);
				return true;
			});
			NeuralFieldForward udf = make_neural_field_forward(p);
			std::vector<Vec3> gradients(all_positions.size(), Vec3(0, 0, 1));
			const size_t grad_batch_cap = 131464;
			const size_t grad_batch =
				std::max<size_t>(1, std::min<size_t>(static_cast<size_t>(p.batch_size_), grad_batch_cap));
			bool grad_ok = true;
			for (size_t offset = 0; offset < all_positions.size(); offset += grad_batch)
			{
				const size_t count = std::min(grad_batch, all_positions.size() - offset);
				auto cpu_opts = torch::TensorOptions().dtype(torch::kFloat32).device(torch::kCPU);
				if (device_.is_cuda())
					cpu_opts = cpu_opts.pinned_memory(true);
				torch::Tensor pts_cpu = torch::empty({static_cast<int64_t>(count), 3}, cpu_opts);
				auto pts_acc = pts_cpu.accessor<float, 2>();
				for (size_t i = 0; i < count; ++i)
				{
					const Vec3& pnt = all_positions[offset + i];
					pts_acc[static_cast<long>(i)][0] = static_cast<float>(pnt.x());
					pts_acc[static_cast<long>(i)][1] = static_cast<float>(pnt.y());
					pts_acc[static_cast<long>(i)][2] = static_cast<float>(pnt.z());
				}

				auto [values_t, grad_t] = udf.forward_values_grad_gpu(pts_cpu);
				if (!values_t.defined() || !grad_t.defined())
				{
					grad_ok = false;
					break;
				}
				if (grad_t.dim() != 2 || grad_t.size(0) != static_cast<long>(count) || grad_t.size(1) != 3)
				{
					grad_ok = false;
					break;
				}

				torch::Tensor grad_cpu = grad_t.to(torch::kCPU).contiguous();
				auto grad_acc = grad_cpu.accessor<float, 2>();
				for (size_t i = 0; i < count; ++i)
				{
					Vec3 normal(static_cast<Scalar>(grad_acc[static_cast<long>(i)][0]),
								static_cast<Scalar>(grad_acc[static_cast<long>(i)][1]),
								static_cast<Scalar>(grad_acc[static_cast<long>(i)][2]));
					if (normal.squaredNorm() > eps)
						normal.normalize();
					else
						normal = Vec3(0, 0, 1);
					gradients[offset + i] = normal;
				}
			}

			if (grad_ok)
			{
				uint32 idx = 0;
				foreach_cell(*p.samples_mesh_, [&](PVertex v) {
					uint32 v_idx = index_of(*p.samples_mesh_, v);
					const Vec3& normal = gradients[idx];
					(*p.samples_normal_)[v_idx] = normal;
					idx++;
					return true;
				});
				normals_ok = true;
			}
		}

		if (!normals_ok && p.input_mode_ == INPUT_SURFACE_MESH && selected_surface_)
		{
			auto s_pos = get_attribute<Vec3, SVertex>(*selected_surface_, "position");
			build_surface_bvh();
			if (s_pos && surface_bvh_)
			{
				foreach_cell(*p.samples_mesh_, [&](PVertex v) {
					uint32 v_idx = index_of(*p.samples_mesh_, v);
					const Vec3& pos = (*p.samples_position_)[v_idx];
					std::pair<uint32, Vec3> cp;
					Vec3 n(0, 0, 1);
					if (surface_bvh_->closest_point(pos, &cp))
					{
						SFace face = surface_bvh_faces_[cp.first];
						n = geometry::normal(*selected_surface_, face, s_pos.get());
						if (surface_vertex_normal_)
						{
							std::array<SVertex, 3> vertices;
							uint32 vi = 0;
							foreach_incident_vertex(*selected_surface_, face, [&](SVertex sv) -> bool {
								if (vi < vertices.size())
									vertices[vi++] = sv;
								return true;
							});
							if (vi == vertices.size())
							{
								const Vec3& p0 = value<Vec3>(*selected_surface_, s_pos, vertices[0]);
								const Vec3& p1 = value<Vec3>(*selected_surface_, s_pos, vertices[1]);
								const Vec3& p2 = value<Vec3>(*selected_surface_, s_pos, vertices[2]);
								Scalar u = 0.0, v_bary = 0.0, w = 0.0;
								cgogn::geometry::closest_point_in_triangle(cp.second, p0, p1, p2, u, v_bary, w);
								const Vec3& n0 = value<Vec3>(*selected_surface_, surface_vertex_normal_, vertices[0]);
								const Vec3& n1 = value<Vec3>(*selected_surface_, surface_vertex_normal_, vertices[1]);
								const Vec3& n2 = value<Vec3>(*selected_surface_, surface_vertex_normal_, vertices[2]);
								n = u * n0 + v_bary * n1 + w * n2;
							}
						}
						if (n.squaredNorm() < eps)
							n = Vec3(0, 0, 1);
						else
							n.normalize();
						Vec3 to_sample = pos - cp.second;
						if (to_sample.squaredNorm() > eps && to_sample.dot(n) < Scalar(0))
							n = -n;
					}
					(*p.samples_normal_)[v_idx] = n;
					return true;
				});
				normals_ok = true;
			}
		}

		if (!normals_ok && p.input_kdtree_ && p.points_ && p.position_)
		{
			foreach_cell(*p.samples_mesh_, [&](PVertex v) {
				uint32 v_idx = index_of(*p.samples_mesh_, v);
				const Vec3& pos = (*p.samples_position_)[v_idx];
				std::pair<uint32, Scalar> knn_res;
				Vec3 n(0, 0, 1);
				if (p.input_kdtree_->find_nn(pos, &knn_res))
				{
					uint32 idx = knn_res.first;
					PVertex vn = p.input_kdtree_vertices_[idx];
					uint32 vn_idx = index_of(*p.points_, vn);
					n = pos - (*p.position_)[vn_idx];
					if (n.squaredNorm() < eps && p.normal_)
						n = (*p.normal_)[vn_idx];
					if (n.squaredNorm() < eps)
						n = Vec3(0, 0, 1);
					else
						n.normalize();
				}
				(*p.samples_normal_)[v_idx] = n;
				return true;
			});
			normals_ok = true;
		}

		if (!normals_ok)
		{
			foreach_cell(*p.samples_mesh_, [&](PVertex v) {
				uint32 v_idx = index_of(*p.samples_mesh_, v);
				(*p.samples_normal_)[v_idx] = Vec3(0, 0, 1);
				return true;
			});
		}

		refresh_sample_normals_color(p);
	}

	void apply_sampling_preprocess_filtering(PointsParameters& p)
	{
		if (!p.samples_mesh_ || !p.samples_position_)
			return;
		if (p.running_)
		{
			std::cerr << "Stop spheres update before filtering sampled points." << std::endl;
			return;
		}

		const uint32 count = nb_cells<PVertex>(*p.samples_mesh_);
		if (count == 0)
		{
			std::cout << "No sampled points to filter." << std::endl;
			return;
		}

		std::vector<Vec3> filtered_points;
		filtered_points.reserve(count);
		foreach_cell(*p.samples_mesh_, [&](PVertex v) {
			uint32 v_idx = index_of(*p.samples_mesh_, v);
			filtered_points.push_back((*p.samples_position_)[v_idx]);
			return true;
		});

		const size_t before = filtered_points.size();
		if (p.input_kdtree_)
			pre_process_sampling_points_kdtree(p, filtered_points);
		else
			pre_process_sampling_points_bvh(p, filtered_points);
		const size_t after = filtered_points.size();
		const bool sample_changed = (after != before);

		std::vector<Vec3> inside_points;
		std::vector<Vec3> inside_projected_points;
		std::vector<Vec3> inside_projected_normals;
		std::vector<PVertex> inside_owners;
		gather_alpha_inside_state(
			p, inside_points, inside_projected_points, inside_projected_normals, inside_owners);
		const size_t inside_before = inside_points.size();
		filter_alpha_inside_points_by_input_distance(
			p, inside_points, &inside_projected_points, &inside_projected_normals, &inside_owners);
		filter_alpha_inside_projected_by_model(
			p, inside_points, inside_projected_points, inside_projected_normals, &inside_owners);
		const size_t inside_after = inside_points.size();
		const bool inside_changed = (inside_after != inside_before);

		if (!sample_changed && !inside_changed)
		{
			std::cout << "Sampling filtering applied: no points removed (samples " << before << " -> " << after
					  << ", inside " << inside_before << " -> " << inside_after << ")." << std::endl;
			return;
		}

		if (sample_changed)
		{
			clear_knn_hover(p);
			p.knn_hover_locked_ = false;
			p.fitting_data_computed_ = false;
			p.samples_winding_number_.reset();
			p.samples_wn_bvh_.reset();
			p.opt_winding_number_.reset();
			p.opt_wn_bvh_.reset();
			p.samples_jitter_backup_valid_ = false;
			p.samples_position_backup_.clear();
			p.samples_normal_backup_.clear();
			points_provider_->clear_mesh(*p.samples_mesh_);

			if (filtered_points.empty())
			{
				if (p.samples_kdtree_)
				{
					delete p.samples_kdtree_;
					p.samples_kdtree_ = nullptr;
				}
				if (p.samples_ma_kdtree_)
				{
					delete p.samples_ma_kdtree_;
					p.samples_ma_kdtree_ = nullptr;
				}
				p.samples_kdtree_vertices_.clear();
				p.samples_ma_kdtree_vertices_.clear();
				points_provider_->emit_connectivity_changed(*p.samples_mesh_);
			}
			else
			{
				for (const Vec3& pt : filtered_points)
				{
					PVertex v = add_vertex(*p.samples_mesh_);
					uint32 v_idx = index_of(*p.samples_mesh_, v);
					(*p.samples_position_)[v_idx] = pt;
					if (p.samples_color_)
						(*p.samples_color_)[v_idx] = Vec4(0.0, 0.0, 0.0, 1.0);
					if (p.samples_knn_color_)
						(*p.samples_knn_color_)[v_idx] = Vec4(0.0, 0.0, 0.0, 1.0);
				}

				recompute_samples_normals_from_current_input(p);
				build_kdtree(p);
				points_provider_->emit_connectivity_changed(*p.samples_mesh_);
				points_provider_->emit_attribute_changed(*p.samples_mesh_, p.samples_position_.get());
				points_provider_->emit_attribute_changed(*p.samples_mesh_, p.samples_normal_.get());
				points_provider_->emit_attribute_changed(*p.samples_mesh_, p.samples_normal_color_.get());
				if (p.samples_knn_color_)
					points_provider_->emit_attribute_changed(*p.samples_mesh_, p.samples_knn_color_.get());
				if (p.samples_color_)
					points_provider_->emit_attribute_changed(*p.samples_mesh_, p.samples_color_.get());
			}
		}

		if (inside_changed)
			write_alpha_inside_samples(
				p, inside_points, inside_projected_points, inside_projected_normals, &inside_owners);

		mark_mixed_optimization_cloud_dirty(p);
		assign_alpha_inside_sphere_power(p);
		ensure_mixed_optimization_cloud_ready(p);

		std::cout << "Sampling filtering applied: samples " << before << " -> " << after << ", inside "
				  << inside_before << " -> " << inside_after << "." << std::endl;
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

	void apply_poisson_eliminate_samples(PointsParameters& p, size_t target_num)
	{
		if (!p.samples_mesh_ || !p.samples_position_)
			return;
		if (p.running_)
		{
			std::cerr << "Stop spheres update before Poisson eliminate on sampled points." << std::endl;
			return;
		}

		const uint32 count = nb_cells<PVertex>(*p.samples_mesh_);
		if (count == 0)
		{
			std::cout << "No sampled points to downsample." << std::endl;
			return;
		}

		std::vector<Vec3> input_points;
		input_points.reserve(count);
		foreach_cell(*p.samples_mesh_, [&](PVertex v) {
			const uint32 v_idx = index_of(*p.samples_mesh_, v);
			input_points.push_back((*p.samples_position_)[v_idx]);
			return true;
		});

		const size_t before = input_points.size();
		const size_t clamped_target = std::max<size_t>(1, std::min(target_num, before));
		if (clamped_target >= before)
		{
			std::cout << "Poisson eliminate skipped: target >= current (" << clamped_target << " >= " << before
					  << ")." << std::endl;
			return;
		}

		std::vector<Vec3> reduced_points = poisson_eliminate_points(input_points, clamped_target);
		if (reduced_points.empty())
		{
			std::cerr << "Poisson eliminate failed: no points generated." << std::endl;
			return;
		}

		clear_knn_hover(p);
		p.knn_hover_locked_ = false;
		p.fitting_data_computed_ = false;
		p.samples_winding_number_.reset();
		p.samples_wn_bvh_.reset();
		p.opt_winding_number_.reset();
		p.opt_wn_bvh_.reset();
		p.samples_jitter_backup_valid_ = false;
		p.samples_position_backup_.clear();
		p.samples_normal_backup_.clear();
		points_provider_->clear_mesh(*p.samples_mesh_);

		for (const Vec3& pt : reduced_points)
		{
			PVertex v = add_vertex(*p.samples_mesh_);
			const uint32 v_idx = index_of(*p.samples_mesh_, v);
			(*p.samples_position_)[v_idx] = pt;
			if (p.samples_color_)
				(*p.samples_color_)[v_idx] = Vec4(0.0, 0.0, 0.0, 1.0);
			if (p.samples_knn_color_)
				(*p.samples_knn_color_)[v_idx] = Vec4(0.0, 0.0, 0.0, 1.0);
		}

		recompute_samples_normals_from_current_input(p);
		build_kdtree(p);
		points_provider_->emit_connectivity_changed(*p.samples_mesh_);
		points_provider_->emit_attribute_changed(*p.samples_mesh_, p.samples_position_.get());
		points_provider_->emit_attribute_changed(*p.samples_mesh_, p.samples_normal_.get());
		points_provider_->emit_attribute_changed(*p.samples_mesh_, p.samples_normal_color_.get());

		if (p.samples_knn_color_)
			points_provider_->emit_attribute_changed(*p.samples_mesh_, p.samples_knn_color_.get());
		if (p.samples_color_)
			points_provider_->emit_attribute_changed(*p.samples_mesh_, p.samples_color_.get());
		mark_mixed_optimization_cloud_dirty(p);

		std::cout << "Poisson eliminate applied: " << before << " -> " << reduced_points.size()
				  << " points (target=" << clamped_target << ")." << std::endl;
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
			const auto [bbox_min, bbox_max] = model_domain_bbox(p);
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
			torch::Tensor bbox_min_t = torch::tensor(
				{static_cast<float>(bbox_min.x()), static_cast<float>(bbox_min.y()), static_cast<float>(bbox_min.z())},
				torch::TensorOptions().dtype(torch::kFloat32).device(device_));
			torch::Tensor bbox_max_t = torch::tensor(
				{static_cast<float>(bbox_max.x()), static_cast<float>(bbox_max.y()), static_cast<float>(bbox_max.z())},
				torch::TensorOptions().dtype(torch::kFloat32).device(device_));
			bbox_min_t = bbox_min_t.view({1, 3});
			bbox_max_t = bbox_max_t.view({1, 3});

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
				X = torch::minimum(torch::maximum(X, bbox_min_t), bbox_max_t);
			}

			auto [values, grad] = udf.forward_values_grad_gpu(X);
			if (!values.defined() || !grad.defined())
				return {};

			torch::Tensor final_X_cpu = X.to(torch::kCPU);
			torch::Tensor final_grad_cpu = grad.to(torch::kCPU);

			torch::Tensor sdf_cpu;
			{
				auto [values_sdf, sdf] = udf.forward_values_sdf_gpu(X);
				if (sdf.defined())
				{
					if (sdf.dim() == 2 && sdf.size(1) == 1)
						sdf = sdf.squeeze(1);
					sdf_cpu = sdf.to(torch::kCPU);
				}
			}

			auto final_X_acc = final_X_cpu.accessor<float, 2>();
			auto final_grad_acc = final_grad_cpu.accessor<float, 2>();

			std::vector<Vec3> projected_points;
			std::vector<Vec3> normals;
			projected_points.reserve(nb_points);
			normals.reserve(nb_points);

			p.last_projection_keep_mask_.assign(nb_points, 1);
			if (sdf_cpu.defined() && sdf_cpu.dim() == 1 &&
				sdf_cpu.size(0) == static_cast<long>(nb_points))
			{
				auto sdf_acc = sdf_cpu.accessor<float, 1>();
				for (size_t i = 0; i < nb_points; ++i)
				{
					const bool keep = (sdf_acc[i] <= 0.0f);
					p.last_projection_keep_mask_[i] = keep ? uint8_t(1) : uint8_t(0);
					if (p.filter_positive_projection_ && !keep)
						continue;
					projected_points.push_back(Vec3(final_X_acc[i][0], final_X_acc[i][1], final_X_acc[i][2]));
					normals.push_back(
						Vec3(final_grad_acc[i][0], final_grad_acc[i][1], final_grad_acc[i][2]).normalized());
				}
			}
			else
			{
				for (size_t i = 0; i < nb_points; ++i)
				{
					projected_points.push_back(Vec3(final_X_acc[i][0], final_X_acc[i][1], final_X_acc[i][2]));
					normals.push_back(
						Vec3(final_grad_acc[i][0], final_grad_acc[i][1], final_grad_acc[i][2]).normalized());
				}
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

protected:
	void init() override
	{
		surface_provider_ = static_cast<ui::MeshProvider<SURFACE>*>(
			app_.module("MeshProvider (" + std::string{mesh_traits<SURFACE>::name} + ")"));

		points_provider_ = static_cast<ui::MeshProvider<POINTS>*>(
			app_.module("MeshProvider (" + std::string{mesh_traits<POINTS>::name} + ")"));

		non_manifold_provider_ = static_cast<ui::MeshProvider<NONMANIFOLD>*>(
			app_.module("MeshProvider (" + std::string{mesh_traits<NONMANIFOLD>::name} + ")"));

		pcr_ = static_cast<PointCloudRender<POINTS>*>(
			app_.module("PointCloudRender (" + std::string{mesh_traits<POINTS>::name} + ")"));

		timer_connection_ = boost::synapse::connect<App::timer_tick>(&app_, [this]() {
			if (selected_points_)
			{
				PointsParameters& p = points_parameters_[selected_points_];
				if (p.running_)
					update_render_data(p, true, false);
				else if (p.pending_full_refresh_after_auto_stop_)
				{
					update_render_data(p);
					p.pending_full_refresh_after_auto_stop_ = false;
				}
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
	static constexpr InputMode ray_sampler_input_mode()
	{
		if constexpr (std::is_same_v<RaySamplerTag, geometry::RaySamplerNeural>)
			return INPUT_NEURAL_UDF;
		else if constexpr (std::is_same_v<RaySamplerTag, geometry::RaySamplerSurface>)
			return INPUT_SURFACE_MESH;
		else
			return INPUT_POINT_CLOUD;
	}

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

		surface_vertex_normal_ = get_or_add_attribute<Vec3, SVertex>(*selected_surface_, "normal");
		geometry::compute_normal<SVertex>(*selected_surface_, s_pos.get(), surface_vertex_normal_.get());
		surface_provider_->emit_attribute_changed(*selected_surface_, surface_vertex_normal_.get());
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
			compute_input_normals(p);
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
		p.samples_knn_color_ = get_or_add_attribute<Vec4, PVertex>(*p.samples_mesh_, "knn_color");

		// Init Alpha-Inside Mesh
		std::string inside_name = points_provider_->mesh_name(*p.points_) + "_alpha_inside";
		if (!p.alpha_inside_mesh_)
			p.alpha_inside_mesh_ = points_provider_->has_mesh(inside_name) ? points_provider_->mesh(inside_name)
																			: points_provider_->add_mesh(inside_name);
		else
			points_provider_->clear_mesh(*p.alpha_inside_mesh_);
		p.alpha_inside_position_ = get_or_add_attribute<Vec3, PVertex>(*p.alpha_inside_mesh_, "position");
		p.alpha_inside_color_ = get_or_add_attribute<Vec4, PVertex>(*p.alpha_inside_mesh_, "color");
		p.alpha_inside_sphere_ = get_or_add_attribute<PVertex, PVertex>(*p.alpha_inside_mesh_, "sphere");
		p.alpha_inside_knn_ = get_or_add_attribute<std::vector<PVertex>, PVertex>(*p.alpha_inside_mesh_, "knn");
		p.alpha_inside_projected_position_ =
			get_or_add_attribute<Vec3, PVertex>(*p.alpha_inside_mesh_, "projected_position");
		p.alpha_inside_projected_normal_ =
			get_or_add_attribute<Vec3, PVertex>(*p.alpha_inside_mesh_, "projected_normal");

		// Init Mixed Optimization Mesh
		std::string opt_name = points_provider_->mesh_name(*p.points_) + "_opt";
		if (!p.opt_mesh_)
			p.opt_mesh_ = points_provider_->has_mesh(opt_name) ? points_provider_->mesh(opt_name)
															   : points_provider_->add_mesh(opt_name);
		else
			points_provider_->clear_mesh(*p.opt_mesh_);
		p.opt_position_ = get_or_add_attribute<Vec3, PVertex>(*p.opt_mesh_, "position");
		p.opt_normal_ = get_or_add_attribute<Vec3, PVertex>(*p.opt_mesh_, "normal");
		p.opt_area_ = get_or_add_attribute<Scalar, PVertex>(*p.opt_mesh_, "area");
		p.opt_knn_ = get_or_add_attribute<std::vector<PVertex>, PVertex>(*p.opt_mesh_, "knn");
		p.opt_quadric_ = get_or_add_attribute<Spherical_Quadric, PVertex>(*p.opt_mesh_, "quadric");
		p.opt_line_quadric_ = get_or_add_attribute<Line_Quadric, PVertex>(*p.opt_mesh_, "line_quadric");
		p.opt_sphere_ = get_or_add_attribute<PVertex, PVertex>(*p.opt_mesh_, "sphere");
		p.opt_error_ = get_or_add_attribute<Scalar, PVertex>(*p.opt_mesh_, "error");
		p.opt_color_ = get_or_add_attribute<Vec4, PVertex>(*p.opt_mesh_, "color");
		p.opt_normal_color_ = get_or_add_attribute<Vec4, PVertex>(*p.opt_mesh_, "normal_color");
		p.opt_ma_position_ = get_or_add_attribute<Vec3, PVertex>(*p.opt_mesh_, "ma_position");
		p.opt_ma_radius_ = get_or_add_attribute<Scalar, PVertex>(*p.opt_mesh_, "ma_radius");
		p.opt_ma_secondary_vertex_ = get_or_add_attribute<PVertex, PVertex>(*p.opt_mesh_, "ma_secondary_vertex");

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
		p.spheres_topup_zero_inside_color_ = get_or_add_attribute<Vec4, PVertex>(*p.spheres_, "topup_zero_inside_color");
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
		p.skeleton_face_udf_color_ = get_or_add_attribute<Vec3, NMFace>(*p.skeleton_, "udf_color");
		p.skeleton_face_boundary_tet_color_ = get_or_add_attribute<Vec3, NMFace>(*p.skeleton_, "boundary_tet_color");
		p.skeleton_edge_color_ = get_or_add_attribute<Vec3, NMEdge>(*p.skeleton_, "color");
		p.skeleton_edge_udf_score_color_ = get_or_add_attribute<Vec3, NMEdge>(*p.skeleton_, "edge_udf_score_color");
		p.skeleton_edge_boundary_tet_color_ =
			get_or_add_attribute<Vec3, NMEdge>(*p.skeleton_, "boundary_tet_best_edge_color");
		p.skeleton_edge_non_manifold_color_ =
			get_or_add_attribute<Vec3, NMEdge>(*p.skeleton_, "non_manifold_edge_color");
		if (p.grid_cell_size_ > Scalar(0))
			p.samples_spatial_grid_ = std::make_unique<SpatialGrid>(p.grid_cell_size_);
		else
			p.samples_spatial_grid_.reset();

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

		// Build mixed cloud first so fitting can run directly on (surface + projected-inside).
		ensure_mixed_optimization_cloud_ready(p);
		const bool use_mixed_fitting = p.opt_mesh_ && p.opt_position_ && p.opt_normal_ && p.opt_area_ && p.opt_knn_ &&
									   p.opt_quadric_ && p.opt_line_quadric_ && p.opt_sphere_ && p.opt_error_ &&
									   p.opt_ma_position_ && p.opt_ma_radius_ && p.opt_ma_secondary_vertex_ &&
									   (nb_cells<PVertex>(*p.opt_mesh_) > 0);
		POINTS* fit_mesh = use_mixed_fitting ? p.opt_mesh_ : p.samples_mesh_;
		auto fit_position = use_mixed_fitting ? p.opt_position_ : p.samples_position_;
		auto fit_ma_position = use_mixed_fitting ? p.opt_ma_position_ : p.samples_ma_position_;

		if (use_mixed_fitting)
		{
			std::cout << "Compute fitting data on mixed cloud (surface + projected inside): "
					  << nb_cells<PVertex>(*p.opt_mesh_) << " points." << std::endl;
		}

		std::cout << "Building KDTree..." << std::endl;
		build_kdtree(p, use_mixed_fitting);
		std::cout << "Recomputing Normals (PCA)..." << std::endl;
		recompute_samples_normals_pca(p, use_mixed_fitting);
		std::cout << "Computing KNN and Area..." << std::endl;
		compute_samples_area(p, use_mixed_fitting);
		std::cout << "Computing Winding Numbers..." << std::endl;
		compute_winding_numbers(p, use_mixed_fitting);
		std::cout << "Computing Quadrics..." << std::endl;
		compute_quadrics(p, use_mixed_fitting);
		std::cout << "Computing Initial Medial Axis..." << std::endl;
		compute_initial_medial_axis(p, use_mixed_fitting);
		std::cout << "Initial Medial Axis computed." << std::endl;

		if (p.neural_udf_loaded_ && fit_ma_position && fit_mesh && fit_position)
		{
			NeuralFieldForward udf = make_neural_field_forward(p);
			if (udf.is_loaded())
			{
				std::vector<Vec3> medial_positions;
				medial_positions.reserve(nb_cells<PVertex>(*fit_mesh));
				foreach_cell(*fit_mesh, [&](PVertex v) {
					uint32 v_idx = index_of(*fit_mesh, v);
					medial_positions.push_back((*fit_ma_position)[v_idx]);
					return true;
				});

				const size_t check_n = std::min<size_t>(30, medial_positions.size());
				if (check_n > 0)
				{
					std::vector<size_t> indices(medial_positions.size());
					std::iota(indices.begin(), indices.end(), size_t(0));
					std::mt19937 gen(p.seed_ + 2024);
					std::shuffle(indices.begin(), indices.end(), gen);

					std::vector<Vec3> check_medial_points;
					check_medial_points.reserve(check_n);
					std::vector<Vec3> check_sample_points;
					check_sample_points.reserve(check_n);
					for (size_t i = 0; i < check_n; ++i)
					{
						size_t idx = indices[i];
						check_medial_points.push_back(medial_positions[idx]);
						check_sample_points.push_back((*fit_position)[static_cast<uint32>(idx)]);
					}

					BatchUDFResult medial_res = udf.forward_batch_with_grad(check_medial_points);
					BatchUDFResult sample_res = udf.forward_batch_with_grad(check_sample_points);
					if (medial_res.ok && sample_res.ok && medial_res.values.size() == check_n &&
						sample_res.values.size() == check_n && medial_res.gradients.size() == check_n &&
						sample_res.gradients.size() == check_n)
					{
						std::cout << "Medial UDF check (random " << check_n << "):" << std::endl;
						for (size_t i = 0; i < check_n; ++i)
						{
							const Scalar medial_udf = medial_res.values[i];
							const Scalar sample_udf = sample_res.values[i];
							const Scalar medial_grad_norm = medial_res.gradients[i].norm();
							const Scalar sample_grad_norm = sample_res.gradients[i].norm();
							std::cout << "  " << i << ": medial_udf=" << medial_udf << " | sample_udf=" << sample_udf
									  << " | sample_udf-alpha=" << (sample_udf - p.alpha_)
									  << " | medial_grad_norm=" << medial_grad_norm
									  << " | sample_grad_norm=" << sample_grad_norm << std::endl;
						}
					}
					else
					{
						std::cerr << "Medial UDF check failed (forward_batch)." << std::endl;
					}
				}
			}
		}

		std::cout << "Fitting Data Computed." << std::endl;

		p.fitting_data_computed_ = true;
		assign_alpha_inside_sphere_power(p);
		// Keep mixed cloud as fitting result source when we already computed fitting directly on it.
		if (!use_mixed_fitting)
			ensure_mixed_optimization_cloud_ready(p);
	}

	void compute_fitting_data_with_initial_ma_mode(PointsParameters& p, InitialMAMode mode, bool force_recompute)
	{
		const InitialMAMode prev_mode = p.initial_ma_mode_override_;
		p.initial_ma_mode_override_ = mode;
		if (force_recompute)
			p.fitting_data_computed_ = false;
		compute_fitting_data(p);
		p.initial_ma_mode_override_ = prev_mode;
	}
	void recompute_sample_knn_graph(PointsParameters& p, int requested_k)
	{
		if (!p.samples_mesh_ || !p.samples_position_ || !p.samples_normal_ || !p.samples_knn_ || !p.samples_area_)
		{
			std::cerr << "[KNNRebuild] samples are not initialized." << std::endl;
			return;
		}

		p.knn_k_ = std::max(1, requested_k);

		std::cout << "[KNNRebuild] rebuilding sample KNN with k=" << p.knn_k_ << std::endl;
		build_kdtree(p);
		compute_samples_area(p); // updates samples_knn_ and samples_area_

		if (p.fitting_data_computed_)
		{
			// Keep fitting-state attributes coherent after KNN changes.
			compute_winding_numbers(p);
			compute_quadrics(p);
			compute_initial_medial_axis(p);
			std::cout << "[KNNRebuild] fitting-dependent attributes refreshed." << std::endl;
		}

		clear_knn_hover(p);
		if (p.samples_knn_color_)
			points_provider_->emit_attribute_changed(*p.samples_mesh_, p.samples_knn_color_.get());

		std::cout << "[KNNRebuild] done samples=" << nb_cells<PVertex>(*p.samples_mesh_)
				  << " k=" << p.knn_k_ << std::endl;
	}

	void build_kdtree(PointsParameters& p, bool use_mixed_cloud = false)
	{
		const bool fit_mixed = use_mixed_cloud && p.opt_mesh_ && p.opt_position_ && (nb_cells<PVertex>(*p.opt_mesh_) > 0);
		POINTS* fit_mesh = p.samples_mesh_;
		std::shared_ptr<PAttribute<Vec3>> fit_position = p.samples_position_;
		std::shared_ptr<PAttribute<Vec3>> fit_ma_position = p.samples_ma_position_;
		std::shared_ptr<PAttribute<Scalar>> fit_ma_radius = p.samples_ma_radius_;
		if (fit_mixed)
		{
			fit_mesh = p.opt_mesh_;
			fit_position = p.opt_position_;
			fit_ma_position = p.opt_ma_position_;
			fit_ma_radius = p.opt_ma_radius_;
		}

		acc::KDTree<3, uint32>*& fit_kdtree = fit_mixed ? p.opt_kdtree_ : p.samples_kdtree_;
		acc::KDTree<3, uint32>*& fit_ma_kdtree = fit_mixed ? p.opt_ma_kdtree_ : p.samples_ma_kdtree_;
		std::vector<PVertex>& fit_kdtree_vertices = fit_mixed ? p.opt_kdtree_vertices_ : p.samples_kdtree_vertices_;
		std::vector<PVertex>& fit_ma_kdtree_vertices = fit_mixed ? p.opt_ma_kdtree_vertices_ : p.samples_ma_kdtree_vertices_;
		if (fit_kdtree)
			delete fit_kdtree;
		if (fit_ma_kdtree)
			delete fit_ma_kdtree;

		if (!fit_mesh || !fit_position)
		{
			fit_kdtree = nullptr;
			fit_ma_kdtree = nullptr;
			fit_kdtree_vertices.clear();
			fit_ma_kdtree_vertices.clear();
			return;
		}

		std::vector<Vec3> points;
		const uint32 sample_count = nb_cells<PVertex>(*fit_mesh);
		fit_kdtree_vertices.clear();
		fit_ma_kdtree_vertices.clear();
		points.reserve(sample_count);
		fit_kdtree_vertices.reserve(sample_count);
		fit_ma_kdtree_vertices.reserve(sample_count);
		std::vector<Vec3> points_ma;
		points_ma.reserve(sample_count);

		foreach_cell(*fit_mesh, [&](PVertex v) {
			uint32 idx = index_of(*fit_mesh, v);
			points.push_back((*fit_position)[idx]);
			fit_kdtree_vertices.push_back(v);
			if (fit_ma_position)
			{
				const Vec3& ma_pos = (*fit_ma_position)[idx];
				const bool ma_finite = ma_pos.allFinite();
				const bool ma_radius_ok = (!fit_ma_radius) || ((*fit_ma_radius)[idx] > Scalar(0));
				if (ma_finite && ma_radius_ok)
				{
					points_ma.push_back(ma_pos);
					fit_ma_kdtree_vertices.push_back(v);
				}
			}
			return true;
		});

		fit_kdtree = points.empty() ? nullptr : new acc::KDTree<3, uint32>(points);
		fit_ma_kdtree = points_ma.empty() ? nullptr : new acc::KDTree<3, uint32>(points_ma);
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
		if (!p.input_kdtree_ || !p.points_ || !p.position_ || !p.normal_ || !p.knn_)
			return;

		// Compute normals
		parallel_foreach_cell(*p.points_, [&](PVertex v) {
			uint32 v_idx = index_of(*p.points_, v);
			const Vec3& pt = (*p.position_)[v_idx];
			std::vector<std::pair<uint32, Scalar>> knn_res;
			p.input_kdtree_->find_nns(pt, p.knn_k_+1, &knn_res);

			std::vector<uint32> indices;
			(*p.knn_)[v_idx].clear();
			for (auto& res : knn_res)
			{
				indices.push_back(res.first);
				(*p.knn_)[v_idx].push_back(p.input_kdtree_vertices_[res.first]);
			}

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
		kdtree.find_nns(pt, p.knn_k_+1, &knn_res);

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
	void compute_samples_area(PointsParameters& p, bool use_mixed_cloud = false)
	{
		const bool fit_mixed =
			use_mixed_cloud && p.opt_mesh_ && p.opt_position_ && p.opt_normal_ && p.opt_knn_ && p.opt_area_ &&
			(nb_cells<PVertex>(*p.opt_mesh_) > 0);
		POINTS* fit_mesh = p.samples_mesh_;
		std::shared_ptr<PAttribute<Vec3>> fit_position = p.samples_position_;
		std::shared_ptr<PAttribute<Vec3>> fit_normal = p.samples_normal_;
		std::shared_ptr<PAttribute<std::vector<PVertex>>> fit_knn = p.samples_knn_;
		std::shared_ptr<PAttribute<Scalar>> fit_area = p.samples_area_;
		if (fit_mixed)
		{
			fit_mesh = p.opt_mesh_;
			fit_position = p.opt_position_;
			fit_normal = p.opt_normal_;
			fit_knn = p.opt_knn_;
			fit_area = p.opt_area_;
		}

		acc::KDTree<3, uint32>* fit_kdtree = fit_mixed ? p.opt_kdtree_ : p.samples_kdtree_;
		const std::vector<PVertex>& fit_kdtree_vertices =
			fit_mixed ? p.opt_kdtree_vertices_ : p.samples_kdtree_vertices_;
		if (!fit_mesh || !fit_position || !fit_normal || !fit_knn || !fit_area || !fit_kdtree)
			return;

		parallel_foreach_cell(*fit_mesh, [&](PVertex v) {
			uint32 v_idx = index_of(*fit_mesh, v);
			const Vec3& pt = (*fit_position)[v_idx];
			std::vector<std::pair<uint32, Scalar>> knn_res;
			fit_kdtree->find_nns(pt, p.knn_k_ + 10, &knn_res);

			(*fit_knn)[v_idx].clear();
			Scalar sum_dist = 0.0;
			const Scalar eps = Scalar(1e-12);
			float band = p.alpha_;
			Vec3 n = (*fit_normal)[v_idx];
			const Scalar n2 = n.squaredNorm();
			if (n2 > eps)
				n /= std::sqrt(n2);
			else
				n = Vec3(0, 0, 1);
			int kept = 0;
			for (auto& res : knn_res)
			{
				if (fit_kdtree_vertices[res.first] != v)
				{
					PVertex nb = fit_kdtree_vertices[res.first];
					uint32 nb_idx = index_of(*fit_mesh, nb);
					const Vec3& q = (*fit_position)[nb_idx];
					const Scalar dn = (q - pt).dot(n);
					if (std::abs(dn) > band)
						continue;
					(*fit_knn)[v_idx].push_back(nb);
					sum_dist += res.second;
					++kept;
					if (kept >= p.knn_k_ +1)
					{
						break;
					}
				}
			}
			if (kept == 0) // fall back
			{
				for (auto& res : knn_res)
				{
					if (fit_kdtree_vertices[res.first] != v)
					{
						PVertex nb = fit_kdtree_vertices[res.first];
						(*fit_knn)[v_idx].push_back(nb);
						sum_dist += res.second;
					}
				}
			}
			// Normals are already computed/oriented in sample_points
			(*fit_area)[v_idx] = (sum_dist * sum_dist) / (2.0 * p.knn_k_); // Rough area estimate
			return true;
		});
		// MA positions changed; keep MA-KDTree in sync for MF topology scoring.
		build_kdtree(p, fit_mixed);
	}

	void compute_winding_numbers(PointsParameters& p, bool use_mixed_cloud = false)
	{
		const bool fit_mixed = use_mixed_cloud && p.opt_mesh_ && p.opt_position_ && p.opt_normal_ && p.opt_area_ &&
								(nb_cells<PVertex>(*p.opt_mesh_) > 0);
		POINTS* fit_mesh = p.samples_mesh_;
		std::shared_ptr<PAttribute<Vec3>> fit_position = p.samples_position_;
		std::shared_ptr<PAttribute<Vec3>> fit_normal = p.samples_normal_;
		std::shared_ptr<PAttribute<Scalar>> fit_area = p.samples_area_;
		if (fit_mixed)
		{
			fit_mesh = p.opt_mesh_;
			fit_position = p.opt_position_;
			fit_normal = p.opt_normal_;
			fit_area = p.opt_area_;
		}
		if (!fit_mesh || !fit_position || !fit_normal || !fit_area)
			return;

		std::vector<Vec3>& wn_centers = fit_mixed ? p.opt_wn_bvh_centers_ : p.samples_wn_bvh_centers_;
		std::vector<Scalar>& wn_radii = fit_mixed ? p.opt_wn_bvh_radii_ : p.samples_wn_bvh_radii_;
		std::vector<PVertex>& wn_vertices = fit_mixed ? p.opt_wn_vertices_ : p.samples_wn_vertices_;
		std::unique_ptr<acc::BVHTreeSpheres<uint32, Vec3>>& wn_bvh = fit_mixed ? p.opt_wn_bvh_ : p.samples_wn_bvh_;
		std::unique_ptr<PointFWN<3>>& winding_number = fit_mixed ? p.opt_winding_number_ : p.samples_winding_number_;
		uint32 nb_samples = nb_cells<PVertex>(*fit_mesh);
		wn_centers.clear();
		wn_radii.clear();
		wn_vertices.clear();
		wn_centers.reserve(nb_samples);
		wn_radii.reserve(nb_samples);
		wn_vertices.reserve(nb_samples);

		foreach_cell(*fit_mesh, [&](PVertex v) {
			uint32 v_idx = index_of(*fit_mesh, v);
			wn_vertices.push_back(v);
			wn_centers.push_back((*fit_position)[v_idx]);
			Scalar area = (*fit_area)[v_idx];
			Scalar radius = std::sqrt(area / M_PI);
			wn_radii.push_back(radius);
			return true;
		});

		wn_bvh = std::make_unique<acc::BVHTreeSpheres<uint32, Vec3>>(wn_centers, wn_radii);

		PointTraits samples_wn_traits(*fit_mesh, wn_bvh.get(), wn_vertices, fit_position.get(), fit_normal.get(),
									  fit_area.get());
		winding_number = std::make_unique<PointFWN<3>>(*wn_bvh, samples_wn_traits, p.beta_);
	}

	void compute_quadrics(PointsParameters& p, bool use_mixed_cloud = false)
	{
		POINTS* fit_mesh = p.samples_mesh_;
		std::shared_ptr<PAttribute<Vec3>> fit_position = p.samples_position_;
		std::shared_ptr<PAttribute<Vec3>> fit_normal = p.samples_normal_;
		std::shared_ptr<PAttribute<Scalar>> fit_area = p.samples_area_;
		std::shared_ptr<PAttribute<std::vector<PVertex>>> fit_knn = p.samples_knn_;
		std::shared_ptr<PAttribute<Spherical_Quadric>> fit_quadric = p.samples_quadric_;
		std::shared_ptr<PAttribute<Line_Quadric>> fit_line_quadric = p.samples_line_quadric_;
		if (use_mixed_cloud && p.opt_mesh_ && p.opt_position_ && p.opt_normal_ && p.opt_area_ && p.opt_knn_ &&
			p.opt_quadric_ && p.opt_line_quadric_ && nb_cells<PVertex>(*p.opt_mesh_) > 0)
		{
			fit_mesh = p.opt_mesh_;
			fit_position = p.opt_position_;
			fit_normal = p.opt_normal_;
			fit_area = p.opt_area_;
			fit_knn = p.opt_knn_;
			fit_quadric = p.opt_quadric_;
			fit_line_quadric = p.opt_line_quadric_;
		}
		if (!fit_mesh || !fit_position || !fit_normal || !fit_area || !fit_knn || !fit_quadric || !fit_line_quadric)
			return;

		parallel_foreach_cell(*fit_mesh, [&](PVertex v) {
			uint32 v_idx = index_of(*fit_mesh, v);
			Spherical_Quadric& q = (*fit_quadric)[v_idx];
			Line_Quadric& lq = (*fit_line_quadric)[v_idx];
			q.clear();
			lq.clear();
			const Vec3& pos = (*fit_position)[v_idx];
			const Vec3& n = (*fit_normal)[v_idx];
			Scalar a = (*fit_area)[v_idx] / (p.knn_k_ + 1.0);
			q += Spherical_Quadric(Vec4(pos.x(), pos.y(), pos.z(), 0), Vec4(n.x(), n.y(), n.z(), 1)) * a;
			lq += Line_Quadric(pos, n) * a;
			for (PVertex vn : (*fit_knn)[v_idx])
			{
				uint32 vn_idx = index_of(*fit_mesh, vn);
				const Vec3& pn = (*fit_position)[vn_idx];
				const Vec3& nn = (*fit_normal)[vn_idx];
				Scalar an = (*fit_area)[vn_idx] / (p.knn_k_ + 1.0);
				q += Spherical_Quadric(Vec4(pn.x(), pn.y(), pn.z(), 0), Vec4(nn.x(), nn.y(), nn.z(), 1)) * an;
				lq += Line_Quadric(pn, nn) * an;
			}
			return true;
		});
	}

	// --- Sampling ---

	void recompute_samples_normals_pca(PointsParameters& p, bool use_mixed_cloud = false)
	{
		const bool fit_mixed = use_mixed_cloud && p.opt_mesh_ && p.opt_position_ && p.opt_normal_ && p.opt_normal_color_ &&
								(nb_cells<PVertex>(*p.opt_mesh_) > 0);
		POINTS* fit_mesh = p.samples_mesh_;
		std::shared_ptr<PAttribute<Vec3>> fit_position = p.samples_position_;
		std::shared_ptr<PAttribute<Vec3>> fit_normal = p.samples_normal_;
		std::shared_ptr<PAttribute<Vec4>> fit_normal_color = p.samples_normal_color_;
		if (fit_mixed)
		{
			fit_mesh = p.opt_mesh_;
			fit_position = p.opt_position_;
			fit_normal = p.opt_normal_;
			fit_normal_color = p.opt_normal_color_;
		}
		acc::KDTree<3, uint32>* fit_kdtree = fit_mixed ? p.opt_kdtree_ : p.samples_kdtree_;
		const std::vector<PVertex>& fit_kdtree_vertices =
			fit_mixed ? p.opt_kdtree_vertices_ : p.samples_kdtree_vertices_;
		if (!fit_mesh || !fit_kdtree || !fit_position || !fit_normal)
			return;

		const int k = std::max(3, p.knn_k_);
		const Scalar eps = Scalar(1e-12);

		parallel_foreach_cell(*fit_mesh, [&](PVertex v) {
			uint32 v_idx = index_of(*fit_mesh, v);
			const Vec3& center = (*fit_position)[v_idx];

			std::vector<std::pair<uint32, Scalar>> knn_res;
			fit_kdtree->find_nns(center, k + 1, &knn_res);

			std::vector<Vec3> neighbors;
			neighbors.reserve(k);
			for (const auto& res : knn_res)
			{
				PVertex nb = fit_kdtree_vertices[res.first];
				if (nb == v)
					continue;
				uint32 nb_idx = index_of(*fit_mesh, nb);
				neighbors.push_back((*fit_position)[nb_idx]);
				if (static_cast<int>(neighbors.size()) >= k)
					break;
			}
			if (neighbors.size() < 3)
				return true;

			Vec3 mean(0, 0, 0);
			for (const Vec3& q : neighbors)
				mean += q;
			mean /= Scalar(neighbors.size());

			Eigen::Matrix<Scalar, 3, 3> cov = Eigen::Matrix<Scalar, 3, 3>::Zero();
			for (const Vec3& q : neighbors)
			{
				Vec3 d = q - mean;
				cov(0, 0) += d.x() * d.x();
				cov(0, 1) += d.x() * d.y();
				cov(0, 2) += d.x() * d.z();
				cov(1, 1) += d.y() * d.y();
				cov(1, 2) += d.y() * d.z();
				cov(2, 2) += d.z() * d.z();
			}
			cov(1, 0) = cov(0, 1);
			cov(2, 0) = cov(0, 2);
			cov(2, 1) = cov(1, 2);

			Eigen::SelfAdjointEigenSolver<Eigen::Matrix<Scalar, 3, 3>> solver(cov);
			if (solver.info() != Eigen::Success)
				return true;

			Eigen::Matrix<Scalar, 3, 1> ev = solver.eigenvectors().col(0);
			Vec3 n(ev(0), ev(1), ev(2));
			Scalar n2 = n.squaredNorm();
			if (n2 < eps)
				return true;

			Vec3 n0 = (*fit_normal)[v_idx];
			if (n0.squaredNorm() > eps && n.dot(n0) < Scalar(0))
				n = -n;
			n.normalize();
			(*fit_normal)[v_idx] = n;
			return true;
		});

		if (fit_normal_color)
		{
			parallel_foreach_cell(*fit_mesh, [&](PVertex v) -> bool {
				uint32 v_idx = index_of(*fit_mesh, v);
				const Vec3& n = (*fit_normal)[v_idx];
				(*fit_normal_color)[v_idx] = Vec4((n.x() + 1.0) * 0.5, (n.y() + 1.0) * 0.5, (n.z() + 1.0) * 0.5, 1.0);
				return true;
			});
		}
		points_provider_->emit_attribute_changed(*fit_mesh, fit_normal.get());
		if (fit_normal_color)
			points_provider_->emit_attribute_changed(*fit_mesh, fit_normal_color.get());
	}
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

	Vec3 random_sample_in_sphere(const Vec3& p, const Scalar radius, std::uniform_real_distribution<Scalar>& uni,
								 std::mt19937& rng)
	{
		Scalar u = uni(rng);
		Scalar v = uni(rng);
		Scalar w = uni(rng);

		Scalar r = radius * std::cbrt(u);
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
		p.opt_winding_number_.reset();
		p.opt_wn_bvh_.reset();
		p.samples_jitter_backup_valid_ = false;
		p.samples_position_backup_.clear();
		p.samples_normal_backup_.clear();
		if (p.samples_mesh_)
			points_provider_->clear_mesh(*p.samples_mesh_);
		clear_alpha_inside_samples(p);
		if (p.opt_mesh_)
			points_provider_->clear_mesh(*p.opt_mesh_);
		p.last_surface_count_ = 0;
		p.last_inside_count_ = 0;
		p.last_mixed_count_ = 0;
		p.mixed_cloud_dirty_ = true;
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
		p.opt_winding_number_.reset();
		p.opt_wn_bvh_.reset();
		mark_mixed_optimization_cloud_dirty(p);
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
		p.opt_winding_number_.reset();
		p.opt_wn_bvh_.reset();
		mark_mixed_optimization_cloud_dirty(p);
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
		p.opt_winding_number_.reset();
		p.opt_wn_bvh_.reset();
		mark_mixed_optimization_cloud_dirty(p);
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
		p.opt_winding_number_.reset();
		p.opt_wn_bvh_.reset();
		p.samples_jitter_backup_valid_ = false;
		p.samples_position_backup_.clear();
		p.samples_normal_backup_.clear();
		mark_mixed_optimization_cloud_dirty(p);
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
	// --- Initial Medial Axis ---

	void compute_initial_medial_axis(PointsParameters& p, bool use_mixed_cloud = false)
	{
		const bool fit_mixed =
			use_mixed_cloud && p.opt_mesh_ && p.opt_position_ && p.opt_normal_ && p.opt_ma_position_ &&
			p.opt_ma_radius_ && p.opt_ma_secondary_vertex_ && (nb_cells<PVertex>(*p.opt_mesh_) > 0);
		POINTS* fit_mesh = p.samples_mesh_;
		std::shared_ptr<PAttribute<Vec3>> fit_position = p.samples_position_;
		std::shared_ptr<PAttribute<Vec3>> fit_normal = p.samples_normal_;
		std::shared_ptr<PAttribute<Vec3>> fit_ma_position = p.samples_ma_position_;
		std::shared_ptr<PAttribute<Scalar>> fit_ma_radius = p.samples_ma_radius_;
		std::shared_ptr<PAttribute<PVertex>> fit_ma_secondary = p.samples_ma_secondary_vertex_;
		if (fit_mixed)
		{
			fit_mesh = p.opt_mesh_;
			fit_position = p.opt_position_;
			fit_normal = p.opt_normal_;
			fit_ma_position = p.opt_ma_position_;
			fit_ma_radius = p.opt_ma_radius_;
			fit_ma_secondary = p.opt_ma_secondary_vertex_;
		}
		acc::KDTree<3, uint32>* fit_kdtree = fit_mixed ? p.opt_kdtree_ : p.samples_kdtree_;
		const std::vector<PVertex>& fit_kdtree_vertices =
			fit_mixed ? p.opt_kdtree_vertices_ : p.samples_kdtree_vertices_;

		if (!fit_mesh || !fit_kdtree || !fit_position || !fit_normal || !fit_ma_position || !fit_ma_radius ||
			!fit_ma_secondary)
			return;

		const Scalar fallback_radius = p.alpha_;
		const Scalar min_norm = Scalar(1e-12);
		const bool neural_input = (p.input_mode_ == INPUT_NEURAL_UDF && p.neural_udf_loaded_);
		const bool mf_model = neural_input && (p.neural_model_type_ == NEURAL_MODEL_MF);
		const bool udf_model = neural_input && (p.neural_model_type_ == NEURAL_MODEL_UDF);
		const bool force_displacement = (p.initial_ma_mode_override_ == INITIAL_MA_DISPLACEMENT);
		const bool force_shrinking_ball = (p.initial_ma_mode_override_ == INITIAL_MA_SHRINKING_BALL);
		const bool use_displacement = force_displacement || (!force_shrinking_ball && udf_model);

		// UDF model (or forced displacement mode): direct MA initialization.
		if (use_displacement)
		{
			parallel_foreach_cell(*fit_mesh, [&](PVertex v) -> bool {
				if (!v.is_valid())
					return true;
				const uint32 v_idx = index_of(*fit_mesh, v);
				if (v_idx == INVALID_INDEX)
					return true;

				const Vec3& pt = (*fit_position)[v_idx];
				Vec3 n = (*fit_normal)[v_idx];
				if (n.squaredNorm() > min_norm)
					n.normalize();
				else
					n = Vec3(0, 0, 1);

				const Vec3 c = pt - n * fallback_radius;
				PVertex secondary;
				if (fit_kdtree && !fit_kdtree_vertices.empty())
				{
					const Vec3 q = c - n * fallback_radius;
					std::pair<uint32, Scalar> knn_res;
					fit_kdtree->find_nn(q, &knn_res);
					if (knn_res.first < fit_kdtree_vertices.size())
						secondary = fit_kdtree_vertices[knn_res.first];
				}

				(*fit_ma_position)[v_idx] = c;
				(*fit_ma_radius)[v_idx] = fallback_radius;
				(*fit_ma_secondary)[v_idx] = secondary;
				return true;
			});

			// Keep MA-KDTree in sync for MF/UDF topology scoring paths.
			build_kdtree(p, fit_mixed);
			return;
		}

		// MF model (and non-neural fallback): shrinking-ball based MA initialization.
		const Scalar initial_radius = std::max<Scalar>(fallback_radius * Scalar(10), Scalar(0));

		auto run_shrinking_ball_for_vertex = [&](PVertex v) -> bool {
			if (!v.is_valid())
				return false;
			const uint32 v_idx = index_of(*fit_mesh, v);
			if (v_idx == INVALID_INDEX)
				return false;
			const Vec3& pt = (*fit_position)[v_idx];
			Vec3 n = (*fit_normal)[v_idx];

			if (n.squaredNorm() > min_norm)
				n.normalize();

			Vec3 c = pt - n * fallback_radius;
			Scalar r = fallback_radius;
			PVertex secondary;

			if (n.squaredNorm() > min_norm)
			{
				auto [c1, r1, q1] = geometry::shrinking_ball_center<PVertex>(
					pt, n, fit_kdtree, fit_kdtree_vertices, initial_radius);

				if (std::isfinite(static_cast<double>(r1)) && r1 > Scalar(0))
				{
					c = c1;
					r = r1;
					secondary = q1;
				}
			}

			if (!secondary.is_valid() && fit_kdtree && !fit_kdtree_vertices.empty())
			{
				const Vec3 q = c - n * r;
				std::pair<uint32, Scalar> knn_res;
				fit_kdtree->find_nn(q, &knn_res);
				if (knn_res.first < fit_kdtree_vertices.size())
					secondary = fit_kdtree_vertices[knn_res.first];
			}

			(*fit_ma_position)[v_idx] = c;
			const Scalar expected_radius = p.alpha_ * p.sqem_fix_radius_scale_;
			(*fit_ma_radius)[v_idx] = (r > expected_radius) ? r : expected_radius;
			(*fit_ma_secondary)[v_idx] = secondary;
			return true;
		};

		auto run_shrinking_ball_for_all = [&]() {
			parallel_foreach_cell(*fit_mesh, [&](PVertex v) -> bool {
				run_shrinking_ball_for_vertex(v);
				return true;
			});
		};

		run_shrinking_ball_for_all();

		// Post-process on shrinking-ball centers:
		// - UDF model: retry if udf(center) > alpha; still > alpha -> delete sample.
		// - MF model: retry if sdf(center) > 0; still sdf > 0 -> delete sample.
		uint32 flipped_normals = 0;
		uint32 flip_triggered_points = 0;
		uint32 deleted_samples = 0;
		uint32 suppressed_deletes = 0;
		if (p.neural_udf_loaded_)
		{
			auto eval_mf_values_sdf = [&](const std::vector<Vec3>& query_points, std::vector<Scalar>& out_values,
										  std::vector<Scalar>& out_sdf) -> bool {
				out_values.clear();
				out_sdf.clear();
				out_values.resize(query_points.size(), Scalar(0));
				out_sdf.resize(query_points.size(), Scalar(0));
				if (query_points.empty())
					return true;
				NeuralFieldForward udf = make_neural_field_forward(p);
				if (!udf.is_loaded())
					return false;
				const size_t batch = std::max<size_t>(1, static_cast<size_t>(p.batch_size_));
				for (size_t offset = 0; offset < query_points.size(); offset += batch)
				{
					const size_t count = std::min(batch, query_points.size() - offset);
					auto cpu_opts = torch::TensorOptions().dtype(torch::kFloat32).device(torch::kCPU);
					if (device_.is_cuda())
						cpu_opts = cpu_opts.pinned_memory(true);
					torch::Tensor points_cpu = torch::empty({static_cast<long>(count), 3}, cpu_opts);
					auto points_acc = points_cpu.accessor<float, 2>();
					for (size_t i = 0; i < count; ++i)
					{
						const Vec3& pnt = query_points[offset + i];
						points_acc[(long)i][0] = static_cast<float>(pnt.x());
						points_acc[(long)i][1] = static_cast<float>(pnt.y());
						points_acc[(long)i][2] = static_cast<float>(pnt.z());
					}
					auto [values_t, sdf_t] = udf.forward_values_sdf_gpu(points_cpu);
					if (!values_t.defined() || !sdf_t.defined() || values_t.numel() != static_cast<long>(count) ||
						sdf_t.numel() != static_cast<long>(count))
						return false;
					if (values_t.dim() == 2 && values_t.size(1) == 1)
						values_t = values_t.squeeze(1);
					if (sdf_t.dim() == 2 && sdf_t.size(1) == 1)
						sdf_t = sdf_t.squeeze(1);
					torch::Tensor values_cpu = values_t.to(torch::kCPU).contiguous();
					torch::Tensor sdf_cpu = sdf_t.to(torch::kCPU).contiguous();
					auto values_acc = values_cpu.accessor<float, 1>();
					auto sdf_acc = sdf_cpu.accessor<float, 1>();
					for (size_t i = 0; i < count; ++i)
					{
						out_values[offset + i] = static_cast<Scalar>(values_acc[(long)i]);
						out_sdf[offset + i] = static_cast<Scalar>(sdf_acc[(long)i]);
					}
				}
				return true;
			};

			std::vector<PVertex> vertices;
			std::vector<Vec3> centers;
			vertices.reserve(nb_cells<PVertex>(*fit_mesh));
			centers.reserve(nb_cells<PVertex>(*fit_mesh));
			foreach_cell(*fit_mesh, [&](PVertex v) -> bool {
				const uint32 vid = index_of(*fit_mesh, v);
				if (vid == INVALID_INDEX)
					return true;
				vertices.push_back(v);
				centers.push_back((*fit_ma_position)[vid]);
				return true;
			});

			std::vector<Scalar> score_values;
			std::vector<Scalar> sdf_values;
			bool eval_ok = false;
			if (mf_model)
				eval_ok = eval_mf_values_sdf(centers, score_values, sdf_values);
			else
				eval_ok = eval_udf_values(p, centers, score_values);

			if (eval_ok && score_values.size() == centers.size())
			{
				if (mf_model)
				{
					const size_t debug_n = std::min<size_t>(20, score_values.size());
					for (size_t i = 0; i < debug_n; ++i)
					{
						const Scalar sdf_i = (i < sdf_values.size()) ? sdf_values[i] : Scalar(0);
						std::cout << "[MAFlipPrune][MF] center_udf#" << i << " udf=" << score_values[i]
								  << " sdf=" << sdf_i << std::endl;
					}
				}

				std::vector<PVertex> need_retry;
				need_retry.reserve(vertices.size() / 8 + 1);
				for (size_t i = 0; i < vertices.size(); ++i)
				{
					const bool by_udf = (score_values[i] > p.alpha_);
					const bool by_mf_sdf = mf_model && i < sdf_values.size() && (sdf_values[i] > Scalar(0));
					const bool trigger_retry = mf_model ? (by_mf_sdf || by_udf) : by_udf;
					if (trigger_retry)
						need_retry.push_back(vertices[i]);
				}
				flip_triggered_points += static_cast<uint32>(need_retry.size());

				for (PVertex v : need_retry)
				{
					const uint32 vid = index_of(*fit_mesh, v);
					if (vid == INVALID_INDEX)
						continue;
					Vec3 n = (*fit_normal)[vid];
					if (n.squaredNorm() > min_norm)
					{
						n.normalize();
						(*fit_normal)[vid] = -n;
						++flipped_normals;
					}
					run_shrinking_ball_for_vertex(v);
				}

				std::vector<Vec3> retry_centers;
				retry_centers.reserve(need_retry.size());
				for (PVertex v : need_retry)
				{
					const uint32 vid = index_of(*fit_mesh, v);
					if (vid == INVALID_INDEX)
					{
						retry_centers.push_back(Vec3(0, 0, 0));
						continue;
					}
					retry_centers.push_back((*fit_ma_position)[vid]);
				}

				std::vector<Scalar> retry_score;
				std::vector<Scalar> retry_sdf;
				bool retry_eval_ok = false;
				if (mf_model)
					retry_eval_ok = eval_mf_values_sdf(retry_centers, retry_score, retry_sdf);
				else
					retry_eval_ok = eval_udf_values(p, retry_centers, retry_score);

				if (retry_eval_ok && retry_score.size() == retry_centers.size())
				{
					std::unordered_set<uint32> deleted_ids;
					deleted_ids.reserve(need_retry.size());
					for (size_t i = 0; i < need_retry.size(); ++i)
					{
						bool should_delete = false;
						if (mf_model)
						{
							if (i < retry_sdf.size() && retry_sdf[i] > Scalar(0))
								should_delete = true;
						}
						else
						{
							if (retry_score[i] > p.alpha_)
								should_delete = true;
						}
						if (!should_delete)
							continue;
						const uint32 vid = index_of(*fit_mesh, need_retry[i]);
						if (vid == INVALID_INDEX || !deleted_ids.insert(vid).second)
							continue;
						if (fit_mixed)
						{
							// Stability path for mixed cloud: avoid in-loop vertex deletion.
							Vec3 n = (*fit_normal)[vid];
							if (n.squaredNorm() > min_norm)
								n.normalize();
							else
								n = Vec3(0, 0, 1);
							(*fit_ma_position)[vid] = (*fit_position)[vid] - n * fallback_radius;
							(*fit_ma_radius)[vid] =
								std::max<Scalar>(fallback_radius, p.alpha_ * p.sqem_fix_radius_scale_);
							(*fit_ma_secondary)[vid] = PVertex();
							++suppressed_deletes;
						}
						else
						{
							remove_vertex(*fit_mesh, need_retry[i]);
							++deleted_samples;
						}
					}
				}
			}
		}

		if (deleted_samples > 0)
		{
			std::cout << "[MAFlipPrune] flip_triggered_points=" << flip_triggered_points
					  << " flipped_normals=" << flipped_normals
					  << " deleted_samples=" << deleted_samples
					  << " suppressed_deletes=" << suppressed_deletes
					  << " remaining_samples=" << nb_cells<PVertex>(*fit_mesh) << std::endl;

			// Keep fitting-dependent structures coherent after sample removals.
			build_kdtree(p, fit_mixed);
			if (nb_cells<PVertex>(*fit_mesh) > 0)
			{
				recompute_samples_normals_pca(p, fit_mixed);
				compute_samples_area(p, fit_mixed);
				compute_winding_numbers(p, fit_mixed);
				compute_quadrics(p, fit_mixed);
				run_shrinking_ball_for_all();
			}
		}
		else if (flip_triggered_points > 0 || flipped_normals > 0)
		{
			std::cout << "[MAFlipPrune] flip_triggered_points=" << flip_triggered_points
					  << " flipped_normals=" << flipped_normals << " deleted_samples=0"
					  << " suppressed_deletes=" << suppressed_deletes << std::endl;
		}

		// MA positions changed; keep MA-KDTree in sync for MF topology scoring.
		build_kdtree(p, fit_mixed);
	}

	// MF post-process: search along opposite normal direction for minimal mf-abs(sdf)

	void init_spheres(PointsParameters& p, uint32 max_nb_spheres)
	{
		points_provider_->clear_mesh(*p.spheres_);
		init_spheres_from_samples(p, max_nb_spheres);

		if (!p.running_)
			update_render_data(p);
	}

	void init_spheres_from_samples(PointsParameters& p, uint32 max_nb_spheres)
	{
		SphereFitData data;
		if (!get_sphere_fit_data(p, data))
			return;
		if (!data.mesh || !data.position || !data.knn || !data.ma_position || !data.ma_radius)
			return;

		std::vector<PVertex> sorted_vertices;
		foreach_cell(*data.mesh, [&](PVertex v) {
			sorted_vertices.push_back(v);
			return true;
		});
		std::sort(sorted_vertices.begin(), sorted_vertices.end(), [&](PVertex a, PVertex b) {
			// sort candidate spheres by decreasing radius
			uint32 idx_a = index_of(*data.mesh, a);
			uint32 idx_b = index_of(*data.mesh, b);
			return (*data.ma_radius)[idx_a] > (*data.ma_radius)[idx_b];
		});

		auto covered = get_or_add_attribute<bool, PVertex>(*data.mesh, "__covered");
		covered->fill(false);

		p.nb_spheres_ = 0;

		for (PVertex v : sorted_vertices)
		{
			uint32 v_index = index_of(*data.mesh, v);

			if (p.nb_spheres_ >= max_nb_spheres)
				break;

			if ((*covered)[v_index])
				continue;

			const Vec3& vp = (*data.ma_position)[v_index];
			Scalar vr = (*data.ma_radius)[v_index];
			const Scalar dilation_radius = std::max<Scalar>(vr + Scalar(p.init_dilation_constant_), Scalar(0));
			const Scalar dilation_radius_sq = dilation_radius * dilation_radius;

			PVertex sphere = add_vertex(*p.spheres_);
			p.nb_spheres_++;
			uint32 sphere_index = index_of(*p.spheres_, sphere);

			(*p.spheres_position_)[sphere_index] = vp;
			(*p.spheres_radius_)[sphere_index] = vr;
			(*p.spheres_cluster_color_)[sphere_index] =
				Vec4(0.5 + 0.5 * (rand() % 256) / 256.0, 0.5 + 0.5 * (rand() % 256) / 256.0,
					 0.5 + 0.5 * (rand() % 256) / 256.0, 1.0);

			auto flood_cover = [&](PVertex seed) {
				if (!seed.is_valid())
					return;
				uint32 seed_idx = index_of(*data.mesh, seed);
				if (seed_idx == INVALID_INDEX || (*covered)[seed_idx])
					return;

				std::vector<PVertex> stack;
				stack.reserve(128);
				// Mark when enqueued to avoid duplicate pushes through overlapping KNN neighborhoods.
				(*covered)[seed_idx] = true;
				stack.push_back(seed);
				while (!stack.empty())
				{
					PVertex w = stack.back();
					stack.pop_back();
					uint32 w_idx = index_of(*data.mesh, w);

					// Use KNN for propagation on point cloud
					for (PVertex u : (*data.knn)[w_idx])
					{
						uint32 u_idx = index_of(*data.mesh, u);
						if (!(*covered)[u_idx] &&
							((*data.position)[u_idx] - vp).squaredNorm() < dilation_radius_sq)
						{
							(*covered)[u_idx] = true;
							stack.push_back(u);
						}
					}
				}
			};

			flood_cover(v);

			// Secondary seed can reach another lobe; skip if already covered to avoid redundant traversal.
			PVertex secondary;
			if (data.ma_secondary)
				secondary = (*data.ma_secondary)[v_index];
			if (secondary.is_valid())
			{
				flood_cover(secondary);
			}
		}

		remove_attribute<PVertex>(*data.mesh, covered);

		compute_clusters_full(p);
	}

	struct SphereFitData
	{
		POINTS* mesh = nullptr;
		std::shared_ptr<PAttribute<Vec3>> position = nullptr;
		std::shared_ptr<PAttribute<Vec3>> normal = nullptr;
		std::shared_ptr<PAttribute<Scalar>> area = nullptr;
		std::shared_ptr<PAttribute<std::vector<PVertex>>> knn = nullptr;
		std::shared_ptr<PAttribute<Spherical_Quadric>> quadric = nullptr;
		std::shared_ptr<PAttribute<Line_Quadric>> line_quadric = nullptr;
		std::shared_ptr<PAttribute<PVertex>> sphere = nullptr;
		std::shared_ptr<PAttribute<Scalar>> error = nullptr;
		std::shared_ptr<PAttribute<Vec3>> ma_position = nullptr;
		std::shared_ptr<PAttribute<Scalar>> ma_radius = nullptr;
		std::shared_ptr<PAttribute<PVertex>> ma_secondary = nullptr;
	};

	bool get_sphere_fit_data(PointsParameters& p, SphereFitData& data)
	{
		const bool use_opt = !p.mixed_cloud_dirty_ && p.opt_mesh_ && p.opt_position_ && p.opt_area_ && p.opt_knn_ && p.opt_quadric_ &&
							 p.opt_line_quadric_ && p.opt_sphere_ && (nb_cells<PVertex>(*p.opt_mesh_) > 0);
		if (use_opt)
		{
			data.mesh = p.opt_mesh_;
			data.position = p.opt_position_;
			data.normal = p.opt_normal_;
			data.area = p.opt_area_;
			data.knn = p.opt_knn_;
			data.quadric = p.opt_quadric_;
			data.line_quadric = p.opt_line_quadric_;
			data.sphere = p.opt_sphere_;
			data.error = p.opt_error_;
			data.ma_position = p.opt_ma_position_;
			data.ma_radius = p.opt_ma_radius_;
			data.ma_secondary = p.opt_ma_secondary_vertex_;
		}
		else
		{
			data.mesh = p.samples_mesh_;
			data.position = p.samples_position_;
			data.normal = p.samples_normal_;
			data.area = p.samples_area_;
			data.knn = p.samples_knn_;
			data.quadric = p.samples_quadric_;
			data.line_quadric = p.samples_line_quadric_;
			data.sphere = p.samples_sphere_;
			data.error = p.samples_error_;
			data.ma_position = p.samples_ma_position_;
			data.ma_radius = p.samples_ma_radius_;
			data.ma_secondary = p.samples_ma_secondary_vertex_;
		}

		if (!data.mesh || !data.position || !data.area || !data.quadric || !data.line_quadric || !data.sphere)
			return false;
		return true;
	}

	void compute_clusters_full(PointsParameters& p)
	{
		SphereFitData data;
		if (!get_sphere_fit_data(p, data))
			return;
		// clean cluster affectation
		parallel_foreach_cell(*p.spheres_, [&](PVertex v) -> bool {
			uint32 v_index = index_of(*p.spheres_, v);
			(*p.spheres_cluster_)[v_index].clear();
			(*p.spheres_cluster_area_)[v_index] = 0.0;
			return true;
		});
		data.sphere->fill(PVertex());

		if (p.nb_spheres_ == 0)
			return;

		parallel_foreach_cell(*data.mesh, [&](PVertex v) -> bool {
			uint32 v_index = index_of(*data.mesh, v);
			// if (!(*p.medial_axis_selected_)[v_index])
			// 	return true;

			Scalar a = (*data.area)[v_index];

			Scalar min_distance = std::numeric_limits<Scalar>::max();
			PVertex closest_sphere;
			uint32 closest_sphere_index;

			foreach_cell(*p.spheres_, [&](PVertex pv) {
				uint32 pv_index = index_of(*p.spheres_, pv);

				const Vec3& center = (*p.spheres_position_)[pv_index];
				Scalar radius = (*p.spheres_radius_)[pv_index];

				Scalar dist_sqem = (*data.quadric)[v_index].eval(Vec4(center.x(), center.y(), center.z(), radius));
				Scalar dist_other = (*data.line_quadric)[v_index].eval(center);
				Scalar dist = dist_sqem + p.sqem_clustering_lambda_ * dist_other;
				if (dist < min_distance)
				{
					min_distance = dist;
					closest_sphere = pv;
					closest_sphere_index = pv_index;
				}
				return true;
			});

			(*data.sphere)[v_index] = closest_sphere;

			std::lock_guard<std::mutex> lock(spheres_mutex_[closest_sphere_index % spheres_mutex_.size()]);
			(*p.spheres_cluster_)[closest_sphere_index].push_back(v);
			(*p.spheres_cluster_area_)[closest_sphere_index] += a;

			return true;
		});
		// augment_insufficient_clusters(p);
	}

	void compute_clusters_local(PointsParameters& p)
	{
		SphereFitData data;
		if (!get_sphere_fit_data(p, data))
			return;
		parallel_foreach_cell(*p.spheres_, [&](PVertex v) -> bool {
			uint32 v_index = index_of(*p.spheres_, v);
			(*p.spheres_cluster_)[v_index].clear();
			(*p.spheres_cluster_area_)[v_index] = 0.0;
			return true;
		});

		if (p.nb_spheres_ == 0)
			return;

		std::atomic<uint32> fallback_global_search_points(0);
		parallel_foreach_cell(*data.mesh, [&](PVertex v) -> bool {
			uint32 v_index = index_of(*data.mesh, v);

			Scalar a = (*data.area)[v_index];
			PVertex cluster_sphere = (*data.sphere)[v_index];
			const uint32 cs_index = cluster_sphere.is_valid() ? index_of(*p.spheres_, cluster_sphere) : INVALID_INDEX;
			const bool use_global_search = (!cluster_sphere.is_valid() || cs_index == INVALID_INDEX);

			Scalar min_distance = std::numeric_limits<Scalar>::max();
			PVertex closest_sphere;
			uint32 closest_sphere_index = INVALID_INDEX;

			auto eval_candidate = [&](PVertex pv) {
				const uint32 pv_index = index_of(*p.spheres_, pv);
				if (pv_index == INVALID_INDEX)
					return;
				const Vec3& center = (*p.spheres_position_)[pv_index];
				const Scalar radius = (*p.spheres_radius_)[pv_index];
				const Scalar dist_sqem = (*data.quadric)[v_index].eval(Vec4(center.x(), center.y(), center.z(), radius));
				const Scalar dist_other = (*data.line_quadric)[v_index].eval(center);
				const Scalar dist = dist_sqem + p.sqem_clustering_lambda_ * dist_other;
				if (dist < min_distance)
				{
					min_distance = dist;
					closest_sphere = pv;
					closest_sphere_index = pv_index;
				}
			};

			if (use_global_search)
			{
				fallback_global_search_points.fetch_add(1, std::memory_order_relaxed);
				foreach_cell(*p.spheres_, [&](PVertex pv) -> bool {
					eval_candidate(pv);
					return true;
				});
			}
			else
			{
				auto neighbors_spheres = (*p.spheres_neighbor_clusters_)[cs_index];
				neighbors_spheres.insert(cluster_sphere);
				for (PVertex pv : neighbors_spheres)
					eval_candidate(pv);
			}

			if (!closest_sphere.is_valid() || closest_sphere_index == INVALID_INDEX)
				return true;

			(*data.sphere)[v_index] = closest_sphere;
			std::lock_guard<std::mutex> lock(spheres_mutex_[closest_sphere_index % spheres_mutex_.size()]);
			(*p.spheres_cluster_)[closest_sphere_index].push_back(v);
			(*p.spheres_cluster_area_)[closest_sphere_index] += a;
			return true;
		});

		const uint32 fallback_points = fallback_global_search_points.load(std::memory_order_relaxed);
		if (fallback_points > 0)
		{
			std::cout << "[ComputeClustersLocal] fallback_global_search_points=" << fallback_points
					  << " total_points=" << nb_cells<PVertex>(*data.mesh) << std::endl;
		}
	}
	void compute_clusters(PointsParameters& p)
	{
		if (p.use_local_clusters_)
			compute_clusters_local(p);
		else
			compute_clusters_full(p);
		prune_empty_clusters(p);
	}

	void prune_empty_clusters(PointsParameters& p)
	{
		if (!p.spheres_ || p.nb_spheres_ == 0)
			return;
		std::vector<PVertex> to_remove;
		to_remove.reserve(p.nb_spheres_);
		foreach_cell(*p.spheres_, [&](PVertex v) -> bool {
			uint32 v_index = index_of(*p.spheres_, v);
			if ((*p.spheres_cluster_)[v_index].empty())
				to_remove.push_back(v);
			return true;
		});
		if (to_remove.empty())
			return;
		for (PVertex v : to_remove)
		{
			if (v.is_valid())
				remove_sphere(p, v);
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
				const bool accept_without_dedup = !p.samples_spatial_grid_ || (p.grid_cell_size_ <= Scalar(0));
				if (accept_without_dedup ||
					p.samples_spatial_grid_->is_valid_sample(pos, p.grid_cell_size_, *p.samples_position_))
				{
					PVertex new_vertex = add_vertex(*p.samples_mesh_);
					uint32 new_idx = index_of(*p.samples_mesh_, new_vertex);
					(*p.samples_position_)[new_idx] = pos;
					(*p.samples_normal_)[new_idx] = normal;
					(*p.samples_normal_color_)[new_idx] =
						Vec4((normal.x() + 1.0) * 0.5, (normal.y() + 1.0) * 0.5, (normal.z() + 1.0) * 0.5, 1.0);
					cluster.push_back(new_vertex);
					(*p.samples_sphere_)[new_idx] = sphere;
					if (p.samples_spatial_grid_)
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
		SphereFitData data;
		if (!get_sphere_fit_data(p, data))
			return;
		// clean cluster affectation
		parallel_foreach_cell(*p.spheres_, [&](PVertex v) -> bool {
			uint32 v_index = index_of(*p.spheres_, v);
			(*p.spheres_cluster_)[v_index].clear();
			(*p.spheres_cluster_area_)[v_index] = 0.0;
			return true;
		});
		data.sphere->fill(PVertex());

		if (p.nb_spheres_ == 0)
			return;

		parallel_foreach_cell(*data.mesh, [&](PVertex v) -> bool {
			uint32 v_index = index_of(*data.mesh, v);
			Scalar a = (*data.area)[v_index];
			const Vec3& vp = (*data.position)[v_index];

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

			(*data.sphere)[v_index] = closest_sphere;

			std::lock_guard<std::mutex> lock(spheres_mutex_[closest_sphere_index % spheres_mutex_.size()]);
			(*p.spheres_cluster_)[closest_sphere_index].push_back(v);
			(*p.spheres_cluster_area_)[closest_sphere_index] += a;

			return true;
		});
		prune_empty_clusters(p);

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

	// Power distance clustering for alpha-inside points (uses original alpha-inside positions).
	// Power distance clustering for alpha-inside points into main clusters (uses original alpha-inside positions).
	// Use power-based clusters to compute sphere neighbors and build skeleton.
	void compute_skeleton_power(PointsParameters& p, bool only_neighbors = false)
	{
		// First compute power clusters.
		compute_power_cluster(p);
		// Then build skeleton using the current cluster assignments.
		compute_skeleton(p, only_neighbors);
	}

	void compute_spheres_error(PointsParameters& p)
	{
		SphereFitData data;
		if (!get_sphere_fit_data(p, data))
			return;
		parallel_foreach_cell(*p.spheres_, [&](PVertex v) {
			uint32 v_index = index_of(*p.spheres_, v);

			const Vec3& center = (*p.spheres_position_)[v_index];
			Scalar radius = (*p.spheres_radius_)[v_index];
			const std::vector<PVertex>& cluster = (*p.spheres_cluster_)[v_index];

			Scalar cluster_error = 0.0;
			for (PVertex sv : cluster)
			{
				uint32 sv_index = index_of(*data.mesh, sv);
				Scalar dist_sqem =
					(*data.quadric)[sv_index].eval(Vec4(center.x(), center.y(), center.z(), radius));
				// Don't multiply by area here since line quadric already incorporates it
				Scalar dist_other = (*data.line_quadric)[sv_index].eval(center);
				Scalar dist = dist_sqem + p.sqem_clustering_lambda_ * dist_other;
				if (data.error)
					(*data.error)[sv_index] = dist;
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
		
	}

	void update_opt_mesh_colors(PointsParameters& p)
	{
		if (!p.opt_mesh_ || !p.opt_color_ || !p.opt_sphere_)
			return;
		parallel_foreach_cell(*p.opt_mesh_, [&](PVertex v) -> bool {
			const uint32 v_index = index_of(*p.opt_mesh_, v);
			Vec4 c(0.35f, 0.35f, 0.35f, 1.0f);
			const PVertex sphere = (*p.opt_sphere_)[v_index];
			if (sphere.is_valid() && p.spheres_ && p.spheres_cluster_color_)
			{
				const uint32 sid = index_of(*p.spheres_, sphere);
				if (sid != INVALID_INDEX && sid < nb_cells<PVertex>(*p.spheres_))
					c = (*p.spheres_cluster_color_)[sid];
			}
			(*p.opt_color_)[v_index] = c;
			return true;
		});
	}
	void update_spheres_topup_zero_inside_color(PointsParameters& p)
	{
		if (!p.spheres_ || !p.spheres_topup_zero_inside_color_)
			return;
		std::vector<uint32> inside_counts(nb_cells<PVertex>(*p.spheres_), uint32(0));
		if (p.alpha_inside_mesh_ && p.alpha_inside_sphere_)
		{
			foreach_cell(*p.alpha_inside_mesh_, [&](PVertex v) -> bool {
				const uint32 vid = index_of(*p.alpha_inside_mesh_, v);
				const PVertex owner = (*p.alpha_inside_sphere_)[vid];
				if (!owner.is_valid())
					return true;
				const uint32 sid = index_of(*p.spheres_, owner);
				if (sid != INVALID_INDEX && sid < inside_counts.size())
					++inside_counts[sid];
				return true;
			});
		}
		parallel_foreach_cell(*p.spheres_, [&](PVertex v) -> bool {
			const uint32 sid = index_of(*p.spheres_, v);
			const bool zero_inside = (sid >= inside_counts.size()) || (inside_counts[sid] == 0);
			(*p.spheres_topup_zero_inside_color_)[sid] = zero_inside
				? Vec4(1.0f, 0.15f, 0.15f, p.spheres_transparency_)
				: Vec4(0.20f, 0.20f, 0.20f, p.spheres_transparency_);
			return true;
		});
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

	bool eval_udf_and_grad(PointsParameters& p, const Vec3& query_point, Scalar& value, Vec3& grad)
	{
		if (!p.neural_udf_loaded_)
			return false;
		NeuralFieldForward udf = make_neural_field_forward(p);
		if (!udf.is_loaded())
			return false;
		auto res = udf.forward_point_with_grad(query_point);
		value = res.first;
		grad = res.second;
		if (!std::isfinite(value) || !grad.allFinite())
			return false;
		return true;
	}

	bool eval_udf_values(PointsParameters& p, const std::vector<Vec3>& points, std::vector<Scalar>& out_values)
	{
		out_values.clear();
		out_values.resize(points.size(), Scalar(0));
		if (points.empty())
			return true;
		if (!p.neural_udf_loaded_)
			return false;
		NeuralFieldForward udf = make_neural_field_forward(p);
		if (!udf.is_loaded())
			return false;

		const size_t batch = std::max<size_t>(1, static_cast<size_t>(p.batch_size_));
		for (size_t offset = 0; offset < points.size(); offset += batch)
		{
			const size_t count = std::min(batch, points.size() - offset);
			auto cpu_opts = torch::TensorOptions().dtype(torch::kFloat32).device(torch::kCPU);
			if (device_.is_cuda())
				cpu_opts = cpu_opts.pinned_memory(true);

			torch::Tensor pts_cpu = torch::empty({static_cast<int64_t>(count), 3}, cpu_opts);
			auto pts_acc = pts_cpu.accessor<float, 2>();
			for (size_t i = 0; i < count; ++i)
			{
				const Vec3& pnt = points[offset + i];
				pts_acc[static_cast<long>(i)][0] = static_cast<float>(pnt.x());
				pts_acc[static_cast<long>(i)][1] = static_cast<float>(pnt.y());
				pts_acc[static_cast<long>(i)][2] = static_cast<float>(pnt.z());
			}

			torch::Tensor pts = pts_cpu.to(device_);
			torch::Tensor values = udf.forward_values_gpu(pts);
			if (!values.defined())
				return false;
			if (values.dim() == 2 && values.size(1) == 1)
				values = values.squeeze(1);
			torch::Tensor values_cpu = values.to(torch::kCPU).contiguous();
			auto values_acc = values_cpu.accessor<float, 1>();
			for (size_t i = 0; i < count; ++i)
				out_values[offset + i] = static_cast<Scalar>(values_acc[static_cast<long>(i)]);
		}
		return true;
	}

	bool eval_topology_score_values(PointsParameters& p, const std::vector<Vec3>& points,
									std::vector<Scalar>& out_values)
	{
		out_values.clear();
		out_values.resize(points.size(), Scalar(0));
		if (points.empty())
			return true;

		const bool use_mf_ma_distance =
			(p.input_mode_ == INPUT_NEURAL_UDF && p.neural_model_type_ == NEURAL_MODEL_MF);
		if (!use_mf_ma_distance)
			return eval_udf_values(p, points, out_values);

		SphereFitData data;
		if (!get_sphere_fit_data(p, data))
		{
			std::cerr << "[TopologyScore][MF] missing sample/ma data for MF score evaluation." << std::endl;
			return false;
		}
		const bool fit_mixed = (data.mesh == p.opt_mesh_);
		acc::KDTree<3, uint32>* fit_ma_kdtree = fit_mixed ? p.opt_ma_kdtree_ : p.samples_ma_kdtree_;
		const std::vector<PVertex>& fit_ma_vertices =
			fit_mixed ? p.opt_ma_kdtree_vertices_ : p.samples_ma_kdtree_vertices_;
		if (!fit_ma_kdtree || fit_ma_vertices.empty())
			build_kdtree(p, fit_mixed);
		fit_ma_kdtree = fit_mixed ? p.opt_ma_kdtree_ : p.samples_ma_kdtree_;
		const std::vector<PVertex>& fit_ma_vertices_ref =
			fit_mixed ? p.opt_ma_kdtree_vertices_ : p.samples_ma_kdtree_vertices_;
		if (!fit_ma_kdtree || fit_ma_vertices_ref.empty() || !data.ma_position || !data.position)
		{
			std::cerr << "[TopologyScore][MF] MA KDTree unavailable (no valid ma_position / ma_radius). "
					  << "Run fitting data computation first." << std::endl;
			return false;
		}

		for (size_t i = 0; i < points.size(); ++i)
		{
			std::pair<uint32, Scalar> knn_res;
			fit_ma_kdtree->find_nn(points[i], &knn_res);
			const uint32 nn_idx = knn_res.first;
			if (nn_idx >= fit_ma_vertices_ref.size())
				continue;
			const PVertex nearest_sample = fit_ma_vertices_ref[nn_idx];
			const uint32 sid = index_of(*data.mesh, nearest_sample);
			if (sid == INVALID_INDEX)
				continue;
			const Vec3& ma_pos = (*data.ma_position)[sid];
			out_values[i] = (points[i] - ma_pos).norm();
		}
		return true;
	}

	void update_sphere_line_quadric_distance_fix_radius(PointsParameters& p, PVertex sphere)
	{
		SphereFitData data;
		if (!get_sphere_fit_data(p, data))
			return;
		uint32 sphere_index = index_of(*p.spheres_, sphere);

		const std::vector<PVertex>& cluster = (*p.spheres_cluster_)[sphere_index];
		if (cluster.empty())
			return;
		Vec3 c = (*p.spheres_position_)[sphere_index];
		Scalar r = (*p.spheres_radius_)[sphere_index];
		Spherical_Quadric q;
		Line_Quadric lq;
		Scalar area = 0.0;
		Vec3 h;
		h.setZero();
		for (PVertex v : cluster)
		{
			uint32 v_index = index_of(*data.mesh, v);
			Scalar weight = value<Scalar>(*data.mesh, data.area, v);
			if (weight <= 0.0)
				std::cout << "Warning: sample with zero volume weight in sphere " << sphere_index << std::endl;
			q += (*data.quadric)[v_index] * weight;
			h += weight * (*data.position)[v_index];
			lq += (*data.line_quadric)[v_index] * weight;

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

		if (p.udf_center_enabled_)
		{
			Scalar f0;
			Vec3 g;
			if (eval_udf_and_grad(p, c, f0, g))
			{
				const Scalar g2 = g.squaredNorm();
				const Scalar eps = Scalar(1e-12);
				if (g2 > eps && area > Scalar(0))
				{
					const Scalar mu = Scalar(p.udf_center_lambda_) /** area*/;
					const Scalar t = g.dot(c) - f0;
					A.noalias() += mu * (g * g.transpose());
					b.noalias() += mu * t * g;
				}
			}
		}

		c = A.ldlt().solve(b);
		r = p.alpha_ * p.sqem_fix_radius_scale_;
		(*p.spheres_position_)[sphere_index] = c;
		(*p.spheres_radius_)[sphere_index] = r;
	}

	void update_sphere_line_quadric_distance_free_radius(PointsParameters& p, PVertex sphere)
	{
		SphereFitData data;
		if (!get_sphere_fit_data(p, data))
			return;
		uint32 sphere_index = index_of(*p.spheres_, sphere);

		const std::vector<PVertex>& cluster = (*p.spheres_cluster_)[sphere_index];
		if (cluster.empty())
			return;
		Spherical_Quadric q;
		Line_Quadric lq;
		Scalar weight_sum = Scalar(0);
		for (PVertex v : cluster)
		{
			uint32 v_index = index_of(*data.mesh, v);
			Scalar weight = value<Scalar>(*data.mesh, data.area, v);
			if (weight <= 0.0)
				continue;
			q += (*data.quadric)[v_index] * weight;
			lq += (*data.line_quadric)[v_index] * weight;
			weight_sum += weight;
		}
		if (weight_sum <= Scalar(0))
			return;

		Mat4 Ql = lq.get_quadric().matrix();
		Mat3 Al = Ql.block<3, 3>(0, 0);
		Vec3 bl = -Ql.block<3, 1>(0, 3);

		Mat4 Al_ext = Mat4::Zero();
		Al_ext.block<3, 3>(0, 0) = Al;
		Vec4 bl_ext = Vec4::Zero();
		bl_ext.head<3>() = bl;

		Mat4 A = q._A + p.sqem_update_lambda_ * Al_ext;
		Vec4 b = q._b + p.sqem_update_lambda_ * bl_ext;
		if (p.udf_center_enabled_)
		{
			const Vec3 c0 = (*p.spheres_position_)[sphere_index];
			Scalar f0;
			Vec3 g;
			if (eval_udf_and_grad(p, c0, f0, g))
			{
				const Scalar g2 = g.squaredNorm();
				const Scalar eps = Scalar(1e-12);
				if (g2 > eps)
				{
					const Scalar mu = Scalar(p.udf_center_lambda_) /** weight_sum*/;
					const Scalar t = g.dot(c0) - f0;
					A.block<3, 3>(0, 0).noalias() += mu * (g * g.transpose());
					b.head<3>().noalias() += mu * t * g;
				}
			}
		}
		Vec4 s = A.completeOrthogonalDecomposition().solve(b);
		if (!s.allFinite())
			return;

		if (s[3] > Scalar(0) && s[3] <= p.alpha_)
		{
			(*p.spheres_position_)[sphere_index] = s.head<3>();
			Scalar expected_radius = p.alpha_ * p.sqem_fix_radius_scale_;
			(*p.spheres_radius_)[sphere_index] = s[3] > expected_radius ? s[3] : expected_radius;
			return;
		}

		update_sphere_line_quadric_distance_fix_radius(p, sphere);
	}
	void correct_sphere(PointsParameters& p, PVertex v)
	{
		SphereFitData data;
		if (!get_sphere_fit_data(p, data))
			return;
		const bool fit_mixed = (data.mesh == p.opt_mesh_);
		std::unique_ptr<PointFWN<3>>& winding_number = fit_mixed ? p.opt_winding_number_ : p.samples_winding_number_;
		acc::KDTree<3, uint32>* fit_kdtree = fit_mixed ? p.opt_kdtree_ : p.samples_kdtree_;
		const std::vector<PVertex>& fit_kdtree_vertices =
			fit_mixed ? p.opt_kdtree_vertices_ : p.samples_kdtree_vertices_;
		if (!winding_number || !fit_kdtree || fit_kdtree_vertices.empty() || !data.position)
			return;

		uint32 v_index = index_of(*p.spheres_, v);

		Vec3& c = (*p.spheres_position_)[v_index];
		Scalar& r = (*p.spheres_radius_)[v_index];

		bool inside = winding_number->is_inside(c);

		std::pair<uint32, Scalar> k_res;
		fit_kdtree->find_nn(c, &k_res);
		if (k_res.first >= fit_kdtree_vertices.size())
			return;
		uint32 k_idx = index_of(*data.mesh, fit_kdtree_vertices[k_res.first]);
		if (k_idx == INVALID_INDEX)
			return;
		Vec3 closest_pos = (*data.position)[k_idx];
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
				if (p.distance_mode_ == LINE_QUADRIC_DISTANCE_FREE_RADIUS)
					update_sphere_line_quadric_distance_free_radius(p, v);
				else
					update_sphere_line_quadric_distance_fix_radius(p, v);
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

	void mark_empty_cluster_spheres_for_deletion(PointsParameters& p)
	{
		p.spheres_pending_delete_.clear();
		p.pending_delete_count_ = 0;
		if (!p.spheres_)
			return;
		const uint32 ns = nb_cells<PVertex>(*p.spheres_);
		if (ns == 0)
			return;
		if (p.connectivity_mode_ == MODE_A_MIXED_KNN)
		{
			foreach_cell(*p.spheres_, [&](PVertex s) -> bool {
				const uint32 sid = index_of(*p.spheres_, s);
				if ((*p.spheres_cluster_)[sid].empty())
					p.spheres_pending_delete_.push_back(s);
				return true;
			});
			p.pending_delete_count_ = static_cast<uint32>(p.spheres_pending_delete_.size());
			return;
		}
		if (!p.alpha_inside_mesh_ || !p.alpha_inside_sphere_)
			return;
		std::vector<uint32> inside_count(ns, 0);
		foreach_cell(*p.alpha_inside_mesh_, [&](PVertex v) -> bool {
			const uint32 vid = index_of(*p.alpha_inside_mesh_, v);
			const PVertex owner = (*p.alpha_inside_sphere_)[vid];
			if (!owner.is_valid())
				return true;
			const uint32 sid = index_of(*p.spheres_, owner);
			if (sid != INVALID_INDEX && sid < ns)
				++inside_count[sid];
			return true;
		});
		foreach_cell(*p.spheres_, [&](PVertex s) -> bool {
			const uint32 sid = index_of(*p.spheres_, s);
			if (inside_count[sid] == 0)
				p.spheres_pending_delete_.push_back(s);
			return true;
		});
		p.pending_delete_count_ = static_cast<uint32>(p.spheres_pending_delete_.size());
	}
	void apply_pending_sphere_deletions(PointsParameters& p)
	{
		if (p.spheres_pending_delete_.empty())
			return;
		for (PVertex s : p.spheres_pending_delete_)
		{
			if (!s.is_valid())
				continue;
			const uint32 sid = index_of(*p.spheres_, s);
			if (sid == INVALID_INDEX)
				continue;
			remove_sphere(p, s);
		}
		p.spheres_pending_delete_.clear();
		p.pending_delete_count_ = 0;
	}

	void prepare_hybrid_iteration_data(PointsParameters& p, bool run_topup)
	{
		if (run_topup)
			topup_alpha_inside_samples_from_current_mesh(p, "hybrid_prepare");
		else if (p.connectivity_mode_ == MODE_B_INSIDE_POWER)
			assign_alpha_inside_sphere_power(p);
		ensure_mixed_optimization_cloud_ready(p);
	}

	bool should_refresh_local_connectivity(const PointsParameters& p)
	{
		if (!p.use_local_clusters_ || p.lock_skeleton_connectivity_)
			return false;
		const uint32 interval = std::max<uint32>(1, p.local_cluster_connectivity_refresh_interval_);
		return (p.iteration_count_ % interval) == 0;
	}

	void update_spheres(PointsParameters& p)
	{
		auto set_stage = [&](const char* stage) {
			p.last_update_stage_ = stage;
			std::cout << "[UpdateSpheresStage] iter=" << (p.iteration_count_ + 1) << " stage=" << stage
					  << " spheres=" << p.nb_spheres_ << std::endl;
		};
		const uint32 inside_before =
			(p.alpha_inside_mesh_ != nullptr) ? static_cast<uint32>(nb_cells<PVertex>(*p.alpha_inside_mesh_)) : uint32(0);
		set_stage("prepare_hybrid_iteration_data");
		prepare_hybrid_iteration_data(p, p.topup_inside_before_reassign_in_loop_);
		const uint32 inside_after_prepare =
			(p.alpha_inside_mesh_ != nullptr) ? static_cast<uint32>(nb_cells<PVertex>(*p.alpha_inside_mesh_)) : uint32(0);
		const uint32 topup_added =
			(inside_after_prepare > inside_before) ? (inside_after_prepare - inside_before) : uint32(0);
		if (p.nb_spheres_ == 0)
		{
			p.stop_reason_ = "No spheres available after hybrid iteration preparation.";
			std::cerr << "[UpdateSpheresStage] " << p.stop_reason_ << std::endl;
			p.force_full_refresh_on_stop_ = true;
			p.stopping_ = true;
			return;
		}
		if (should_refresh_local_connectivity(p))
		{
			set_stage("refresh_connectivity_interval");
			compute_skeleton_current_mode(p, true);
		}

		set_stage("compute_clusters");
		compute_clusters(p);
		if (p.nb_spheres_ == 0)
		{
			p.stop_reason_ = "All spheres were removed during cluster computation.";
			std::cerr << "[UpdateSpheresStage] " << p.stop_reason_ << std::endl;
			p.force_full_refresh_on_stop_ = true;
			p.stopping_ = true;
			return;
		}

		set_stage("update_sphere_parameters");
		parallel_foreach_cell(*p.spheres_, [&](PVertex v) -> bool {
			if (p.distance_mode_ == LINE_QUADRIC_DISTANCE_FREE_RADIUS)
				update_sphere_line_quadric_distance_free_radius(p, v);
			else
				update_sphere_line_quadric_distance_fix_radius(p, v);
			return true;
		});

		if (p.sphere_correction_ && p.sphere_correction_mode_ == CORRECT_ALWAYS)
		{
			set_stage("correct_spheres");
			parallel_foreach_cell(*p.spheres_, [&](PVertex v) -> bool {
				correct_sphere(p, v);
				return true;
			});
		}

		parallel_foreach_cell(*p.spheres_, [&](PVertex v) -> bool {
			(*p.spheres_do_not_split_)[index_of(*p.spheres_, v)] = false;
			return true;
		});

		set_stage("compute_spheres_error");
		compute_spheres_error(p);

		if (p.auto_split_ && (p.total_error_diff_ < 1e-5 || p.iteration_count_ % 10 == 0))
		{
			switch (p.auto_split_mode_)
			{
			case ERROR_THRESHOLD: {
				if (p.max_error_ > p.auto_split_error_threshold_)
				{
					if (!p.lock_skeleton_connectivity_)
					{
						set_stage("compute_skeleton_for_split");
						compute_skeleton_current_mode(p, true);
					}

					if (p.sphere_correction_ && p.sphere_correction_mode_ == CORRECT_ON_SPLIT)
					{
						set_stage("correct_spheres_on_split");
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

					uint32 to_split_max = std::min(uint32(std::ceil(p.nb_spheres_ * p.auto_split_ratio_)),
												  p.auto_split_max_per_iter_error_);
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
					if (!p.lock_skeleton_connectivity_)
					{
						set_stage("compute_skeleton_for_split");
						compute_skeleton_current_mode(p, true);
					}

					if (p.sphere_correction_ && p.sphere_correction_mode_ == CORRECT_ON_SPLIT)
					{
						set_stage("correct_spheres_on_split");
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

					uint32 to_split_max = std::min(uint32(std::ceil(p.nb_spheres_ * p.auto_split_ratio_)),
												  p.auto_split_max_per_iter_max_);
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

		set_stage("cleanup_empty_clusters");
		if (p.connectivity_mode_ == MODE_B_INSIDE_POWER)
			assign_alpha_inside_sphere_power(p);
		mark_empty_cluster_spheres_for_deletion(p);
		const uint32 marked_for_delete = p.pending_delete_count_;
		apply_pending_sphere_deletions(p);

		std::cout << "[HybridLoop] iter=" << (p.iteration_count_ + 1)
				  << " topup=" << (p.topup_inside_before_reassign_in_loop_ ? "on" : "off")
				  << " topup_added=" << topup_added << " surface=" << p.last_surface_count_
				  << " inside=" << p.last_inside_count_ << " mixed=" << p.last_mixed_count_
				  << " marked_delete=" << marked_for_delete << " spheres=" << p.nb_spheres_ << std::endl;

		if (p.nb_spheres_ == 0)
		{
			p.stop_reason_ = "No spheres left after empty-cluster cleanup.";
			std::cerr << "[HybridLoop] " << p.stop_reason_ << std::endl;
			p.force_full_refresh_on_stop_ = true;
			p.stopping_ = true;
		}
		set_stage("done");

		if (!p.running_)
			update_render_data(p);
	}

	void remove_sphere(PointsParameters& p, PVertex v)
	{
		SphereFitData data;
		if (!get_sphere_fit_data(p, data))
			return;
		const uint32 v_index = index_of(*p.spheres_, v);
		const std::vector<PVertex>& cluster = (*p.spheres_cluster_)[v_index];
		for (PVertex s : cluster)
		{
			uint32 s_idx = index_of(*data.mesh, s);
			(*data.sphere)[s_idx] = PVertex();
		}
		if (p.alpha_inside_mesh_ && p.alpha_inside_sphere_)
		{
			parallel_foreach_cell(*p.alpha_inside_mesh_, [&](PVertex iv) -> bool {
				const uint32 iid = index_of(*p.alpha_inside_mesh_, iv);
				if ((*p.alpha_inside_sphere_)[iid] == v)
					(*p.alpha_inside_sphere_)[iid] = PVertex();
				return true;
			});
		}

		// Keep neighbor lists consistent by removing the deleted sphere from adjacent sets.
		if (p.spheres_neighbor_clusters_)
		{
			const std::set<PVertex>& neighbors = (*p.spheres_neighbor_clusters_)[v_index];
			for (PVertex neighbor : neighbors)
			{
				if (!neighbor.is_valid())
					continue;
				const uint32 d_idx = neighbor.dart_.index_;
				if (d_idx >= p.spheres_->darts_.maximum_index())
					continue;
				const uint32 n_index = index_of(*p.spheres_, neighbor);
				if (n_index == INVALID_INDEX)
					continue;
				(*p.spheres_neighbor_clusters_)[n_index].erase(v);
			}
			(*p.spheres_neighbor_clusters_)[v_index].clear();
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

	void compute_neighbors_from_alpha_inside_power(PointsParameters& p)
	{
		if (!p.spheres_ || !p.spheres_neighbor_clusters_)
			return;
		parallel_foreach_cell(*p.spheres_, [&](PVertex v) -> bool {
			const uint32 sid = index_of(*p.spheres_, v);
			(*p.spheres_neighbor_clusters_)[sid].clear();
			return true;
		});
		if (!p.alpha_inside_mesh_ || !p.alpha_inside_position_ || !p.alpha_inside_sphere_ || !p.alpha_inside_knn_)
			return;
		if (nb_cells<PVertex>(*p.alpha_inside_mesh_) == 0)
			return;
		if (!p.alpha_inside_kdtree_)
			build_alpha_inside_kdtree(p);
		compute_alpha_inside_knn_graph(p);
		assign_alpha_inside_sphere_power(p);
		foreach_cell(*p.alpha_inside_mesh_, [&](PVertex v) -> bool {
			const uint32 vid = index_of(*p.alpha_inside_mesh_, v);
			const PVertex vs = (*p.alpha_inside_sphere_)[vid];
			if (!vs.is_valid())
				return true;
			for (PVertex w : (*p.alpha_inside_knn_)[vid])
			{
				const uint32 wid = index_of(*p.alpha_inside_mesh_, w);
				if (wid == INVALID_INDEX)
					continue;
				const PVertex ws = (*p.alpha_inside_sphere_)[wid];
				if (!ws.is_valid() || ws == vs)
					continue;
				const uint32 sv = index_of(*p.spheres_, vs);
				const uint32 sw = index_of(*p.spheres_, ws);
				if (sv == INVALID_INDEX || sw == INVALID_INDEX)
					continue;
				(*p.spheres_neighbor_clusters_)[sv].insert(ws);
				(*p.spheres_neighbor_clusters_)[sw].insert(vs);
			}
			return true;
		});
	}

	void build_skeleton_from_neighbor_graph(PointsParameters& p)
	{
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
							raw_tets.push_back({idx1, idx2, idx3, idx4});
					}
				}
			}
			return true;
		});

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
			++tet_index;
		}
		parallel_foreach_cell(*p.skeleton_, [&](NMFace f) -> bool {
			if (value<std::set<size_t>>(*p.skeleton_, p.incident_tets_, f).size() == 1)
				value<Vec3>(*p.skeleton_, p.skeleton_face_color_, f) = Vec3(1.0, 0.0, 0.0);
			return true;
		});

		compute_edge_degree(p);
		remove_attribute<PVertex>(*p.spheres_, spheres_skeleton_vertex_map);
	}

	void compute_skeleton_inside_connectivity(
		PointsParameters& p, bool only_neighbors = false, bool run_topup_before_build = false)
	{
		if (run_topup_before_build)
		{
			topup_alpha_inside_samples_from_current_mesh(p, "skeleton_build");
			ensure_mixed_optimization_cloud_ready(p);
		}
		compute_neighbors_from_alpha_inside_power(p);
		if (!only_neighbors)
			build_skeleton_from_neighbor_graph(p);
	}

	void compute_skeleton_current_mode(PointsParameters& p, bool only_neighbors = false)
	{
		if (p.connectivity_mode_ == MODE_B_INSIDE_POWER)
		{
			compute_skeleton_inside_connectivity(p, only_neighbors, false);
			return;
		}
		compute_skeleton(p, only_neighbors);
	}
	void compute_skeleton(PointsParameters& p, bool only_neighbors = false)
	{
		SphereFitData data;
		if (!get_sphere_fit_data(p, data))
			return;
		auto knn_attr = data.knn;
		if (!data.mesh || !knn_attr || !data.sphere)
			return;
		// clear graph
		parallel_foreach_cell(*p.spheres_, [&](PVertex v) -> bool {
			uint32 v_index = index_of(*p.spheres_, v);
			(*p.spheres_neighbor_clusters_)[v_index].clear();
			return true;
		});

		foreach_cell(*data.mesh, [&](PVertex v) -> bool {
			uint32 v_index = index_of(*data.mesh, v);
			PVertex v_sphere = (*data.sphere)[v_index];
			for (PVertex w : (*knn_attr)[v_index])
			{
				PVertex w_sphere = (*data.sphere)[index_of(*data.mesh, w)];
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
		non_manifold_provider_->emit_attribute_changed(*p.skeleton_, p.skeleton_edge_non_manifold_color_.get());
		invalidate_topology_stage_snapshot(p);
	}

	void skeleton_geometry_filter(PointsParameters& p)
	{
		if (!p.neural_udf_loaded_)
		{
			std::cerr << "Geometry filter requires a loaded neural UDF model." << std::endl;
			return;
		}

		compute_skeleton_current_mode(p, true);
		if (!p.spheres_ || !p.spheres_position_ || !p.spheres_neighbor_clusters_ || !p.skeleton_)
			return;

		auto get_face_key = [](uint32 i1, uint32 i2, uint32 i3) -> NMFaceKey {
			std::array<uint32, 3> key = {i1, i2, i3};
			std::sort(key.begin(), key.end());
			return key;
		};
		auto edge_key = [](uint32 a, uint32 b) -> std::pair<uint32, uint32> {
			return {std::min(a, b), std::max(a, b)};
		};
		std::set<std::pair<uint32, uint32>> edge_candidates_set;
		std::vector<std::pair<uint32, uint32>> edge_candidates;
		std::vector<Vec3> edge_midpoints;
		foreach_cell(*p.spheres_, [&](PVertex pv) -> bool {
			const uint32 idx1 = index_of(*p.spheres_, pv);
			const std::set<PVertex>& neighbors = (*p.spheres_neighbor_clusters_)[idx1];
			for (PVertex neighbor : neighbors)
			{
				const uint32 idx2 = index_of(*p.spheres_, neighbor);
				if (idx2 == INVALID_INDEX || idx1 == idx2)
					continue;
				const auto key = edge_key(idx1, idx2);
				if (!edge_candidates_set.insert(key).second)
					continue;
				edge_candidates.push_back(key);
				const Vec3& c1 = (*p.spheres_position_)[key.first];
				const Vec3& c2 = (*p.spheres_position_)[key.second];
				edge_midpoints.push_back((c1 + c2) * Scalar(0.5));
			}
			return true;
		});

		std::vector<Scalar> edge_udf_values;
		if (!edge_midpoints.empty() && !eval_udf_values(p, edge_midpoints, edge_udf_values))
		{
			std::cerr << "Geometry filter failed: unable to evaluate UDF on edge midpoints." << std::endl;
			return;
		}

		std::set<std::pair<uint32, uint32>> kept_edges;
		for (size_t i = 0; i < edge_candidates.size(); ++i)
		{
			if (std::abs(edge_udf_values[i]) <= Scalar(p.skeleton_edge_udf_zero_tol_))
				kept_edges.insert(edge_candidates[i]);
		}

		std::set<NMFaceKey> face_seen;
		std::vector<std::array<uint32, 3>> face_candidates;
		std::vector<Vec3> face_centers;
		std::vector<std::array<uint32, 4>> raw_tets;
		std::set<std::array<uint32, 4>> raw_tet_seen;

		foreach_cell(*p.spheres_, [&](PVertex pv) -> bool {
			const uint32 idx1 = index_of(*p.spheres_, pv);
			const std::set<PVertex>& n_pv = (*p.spheres_neighbor_clusters_)[idx1];
			for (const PVertex& ne1 : n_pv)
			{
				const uint32 idx2 = index_of(*p.spheres_, ne1);
				if (idx2 == INVALID_INDEX || idx1 >= idx2)
					continue;
				const std::set<PVertex>& ne_ne1 = (*p.spheres_neighbor_clusters_)[idx2];
				for (const PVertex& ne2 : ne_ne1)
				{
					if (n_pv.find(ne2) == n_pv.end())
						continue;
					const uint32 idx3 = index_of(*p.spheres_, ne2);
					if (idx3 == INVALID_INDEX || idx2 >= idx3)
						continue;

					const NMFaceKey fkey = get_face_key(idx1, idx2, idx3);
					if (face_seen.insert(fkey).second)
					{
						const bool e12 = kept_edges.find(edge_key(idx1, idx2)) != kept_edges.end();
						const bool e13 = kept_edges.find(edge_key(idx1, idx3)) != kept_edges.end();
						const bool e23 = kept_edges.find(edge_key(idx2, idx3)) != kept_edges.end();
						if (e12 && e13 && e23)
						{
							face_candidates.push_back({idx1, idx2, idx3});
							const Vec3& c1 = (*p.spheres_position_)[idx1];
							const Vec3& c2 = (*p.spheres_position_)[idx2];
							const Vec3& c3 = (*p.spheres_position_)[idx3];
							face_centers.push_back((c1 + c2 + c3) / Scalar(3.0));
						}
					}

					const std::set<PVertex>& ne_ne2 = (*p.spheres_neighbor_clusters_)[idx3];
					for (const PVertex& ne3 : ne_ne2)
					{
						const uint32 idx4 = index_of(*p.spheres_, ne3);
						if (idx4 == INVALID_INDEX || idx3 >= idx4)
							continue;
						const bool connected_v1 = (n_pv.find(ne3) != n_pv.end());
						const bool connected_v2 = (ne_ne1.find(ne3) != ne_ne1.end());
						if (!connected_v1 || !connected_v2)
							continue;

						std::array<uint32, 4> tet = {idx1, idx2, idx3, idx4};
						std::sort(tet.begin(), tet.end());
						if (raw_tet_seen.insert(tet).second)
							raw_tets.push_back(tet);
					}
				}
			}
			return true;
		});

		std::vector<Scalar> face_udf_values;
		if (!face_centers.empty() && !eval_udf_values(p, face_centers, face_udf_values))
		{
			std::cerr << "Geometry filter failed: unable to evaluate UDF on face centers." << std::endl;
			return;
		}

		std::set<NMFaceKey> kept_faces;
		for (size_t i = 0; i < face_candidates.size(); ++i)
		{
			if (std::abs(face_udf_values[i]) <= Scalar(p.skeleton_face_udf_zero_tol_))
				kept_faces.insert(get_face_key(face_candidates[i][0], face_candidates[i][1], face_candidates[i][2]));
		}

		clear(*p.skeleton_);
		p.skeleton_faces_map_.clear();
		auto spheres_skeleton_vertex_map =
			add_attribute<NMVertex, PVertex>(*p.spheres_, "__spheres_skeleton_vertex_map");
		foreach_cell(*p.spheres_, [&](PVertex pv) -> bool {
			const uint32 pv_index = index_of(*p.spheres_, pv);
			const NMVertex nmv = add_vertex(*p.skeleton_);
			(*p.skeleton_position_)[index_of(*p.skeleton_, nmv)] = (*p.spheres_position_)[pv_index];
			(*spheres_skeleton_vertex_map)[pv_index] = nmv;
			return true;
		});

		std::unordered_map<std::pair<uint32, uint32>, NMEdge, edge_hash, edge_equal> edge_indices;
		for (const auto& ekey : kept_edges)
		{
			const NMVertex v1 = (*spheres_skeleton_vertex_map)[ekey.first];
			const NMVertex v2 = (*spheres_skeleton_vertex_map)[ekey.second];
			if (!v1.is_valid() || !v2.is_valid() || v1 == v2)
				continue;
			const NMEdge e = add_edge(*p.skeleton_, v1, v2);
			edge_indices[{index_of(*p.skeleton_, v1), index_of(*p.skeleton_, v2)}] = e;
		}

		auto find_edge = [&](uint32 a, uint32 b, NMEdge& out) -> bool {
			auto it = edge_indices.find({a, b});
			if (it == edge_indices.end() || !it->second.is_valid())
				return false;
			out = it->second;
			return true;
		};

		for (const auto& f : face_candidates)
		{
			const NMFaceKey fkey = get_face_key(f[0], f[1], f[2]);
			if (kept_faces.find(fkey) == kept_faces.end())
				continue;

			const NMVertex v1 = (*spheres_skeleton_vertex_map)[f[0]];
			const NMVertex v2 = (*spheres_skeleton_vertex_map)[f[1]];
			const NMVertex v3 = (*spheres_skeleton_vertex_map)[f[2]];
			if (!v1.is_valid() || !v2.is_valid() || !v3.is_valid())
				continue;

			const uint32 i1 = index_of(*p.skeleton_, v1);
			const uint32 i2 = index_of(*p.skeleton_, v2);
			const uint32 i3 = index_of(*p.skeleton_, v3);
			NMEdge e12, e23, e13;
			if (!find_edge(i1, i2, e12) || !find_edge(i2, i3, e23) || !find_edge(i1, i3, e13))
				continue;

			std::vector<NMEdge> fedges;
			fedges.reserve(3);
			fedges.push_back(e12);
			fedges.push_back(e23);
			fedges.push_back(e13);
			NMFace new_face = add_face(*p.skeleton_, fedges);
			p.skeleton_faces_map_[fkey] = new_face;
		}

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
			const uint32 v[4] = {rt[0], rt[1], rt[2], rt[3]};
			const NMFaceKey keys[4] = {get_face_key(v[1], v[2], v[3]), get_face_key(v[0], v[2], v[3]),
									   get_face_key(v[0], v[1], v[3]), get_face_key(v[0], v[1], v[2])};

			std::array<NMFace, 4> tet_faces = {NMFace(), NMFace(), NMFace(), NMFace()};
			bool complete_tet = true;
			for (uint32 i = 0; i < 4; ++i)
			{
				auto it = p.skeleton_faces_map_.find(keys[i]);
				if (it == p.skeleton_faces_map_.end())
				{
					complete_tet = false;
					break;
				}
				tet_faces[i] = it->second;
			}
			if (!complete_tet)
				continue;

			Tet new_tet;
			for (uint32 i = 0; i < 4; ++i)
			{
				new_tet.faces[i] = tet_faces[i];
				value<std::set<size_t>>(*p.skeleton_, p.incident_tets_, tet_faces[i]).insert(tet_index);
				value<Vec3>(*p.skeleton_, p.skeleton_face_color_, tet_faces[i]) = Vec3(0.8, 0.5, 0.5);
			}
			new_tet.tet_id = tet_index;
			p.skeleton_tets_.insert({tet_index, new_tet});
			++tet_index;
		}

		parallel_foreach_cell(*p.skeleton_, [&](NMFace f) -> bool {
			const auto in_tets = value<std::set<size_t>>(*p.skeleton_, p.incident_tets_, f);
			if (in_tets.size() == 1)
				value<Vec3>(*p.skeleton_, p.skeleton_face_color_, f) = Vec3(1.0, 0.0, 0.0);
			else if (in_tets.empty())
				value<Vec3>(*p.skeleton_, p.skeleton_face_color_, f) = Vec3(0.0, 0.0, 0.0);
			return true;
		});

		parallel_foreach_cell(*p.skeleton_, [&](NMEdge e) -> bool {
			value<Vec3>(*p.skeleton_, p.skeleton_edge_color_, e) = Vec3(0.0, 0.0, 0.0);
			return true;
		});
		compute_edge_degree(p);
		remove_attribute<PVertex>(*p.spheres_, spheres_skeleton_vertex_map);

		std::cout << "Geometry filter: edges " << edge_candidates.size() << " -> " << kept_edges.size()
				  << ", faces " << face_candidates.size() << " -> " << kept_faces.size() << std::endl;

		non_manifold_provider_->emit_connectivity_changed(*p.skeleton_);
		non_manifold_provider_->emit_attribute_changed(*p.skeleton_, p.skeleton_edge_color_.get());
		non_manifold_provider_->emit_attribute_changed(*p.skeleton_, p.skeleton_edge_non_manifold_color_.get());
		non_manifold_provider_->emit_attribute_changed(*p.skeleton_, p.skeleton_face_color_.get());
		invalidate_topology_stage_snapshot(p);
	}

	void inherit_sphere_neighbors(PointsParameters& p, PVertex parent, PVertex child)
	{
		if (!parent.is_valid() || !child.is_valid())
			return;
		if (p.nb_spheres_ == 2)
		{
			uint32 parent_index = index_of(*p.spheres_, parent);
			uint32 child_index = index_of(*p.spheres_, child);
			(*p.spheres_neighbor_clusters_)[parent_index].clear();
			(*p.spheres_neighbor_clusters_)[child_index].clear();
			(*p.spheres_neighbor_clusters_)[parent_index].insert(child);
			(*p.spheres_neighbor_clusters_)[child_index].insert(parent);
			return;
		}
		uint32 parent_index = index_of(*p.spheres_, parent);
		uint32 child_index = index_of(*p.spheres_, child);

		(*p.spheres_neighbor_clusters_)[child_index].clear();
		(*p.spheres_neighbor_clusters_)[child_index].insert(parent);
		for (PVertex neighbor : (*p.spheres_neighbor_clusters_)[parent_index])
		{
			if (!neighbor.is_valid() || neighbor == child)
				continue;
			(*p.spheres_neighbor_clusters_)[child_index].insert(neighbor);
			uint32 n_index = index_of(*p.spheres_, neighbor);
			(*p.spheres_neighbor_clusters_)[n_index].insert(child);
		}
		(*p.spheres_neighbor_clusters_)[parent_index].insert(child);
	}

	void recompute_clusters_local_neighborhood(PointsParameters& p, PVertex center_sphere)
	{
		if (!center_sphere.is_valid() || p.nb_spheres_ == 0)
			return;

		SphereFitData data;
		if (!get_sphere_fit_data(p, data))
			return;
		const uint32 nb_samples = nb_cells<PVertex>(*data.mesh);
		if (nb_samples == 0)
			return;

		uint32 center_index = index_of(*p.spheres_, center_sphere);

		std::vector<PVertex> candidate_spheres;
		candidate_spheres.push_back(center_sphere);
		for (PVertex neighbor : (*p.spheres_neighbor_clusters_)[center_index])
		{
			if (neighbor.is_valid())
				candidate_spheres.push_back(neighbor);
		}

		std::vector<PVertex> samples;
		samples.reserve(nb_samples);
		std::vector<uint8_t> marked(nb_samples, 0);
		for (PVertex sphere : candidate_spheres)
		{
			uint32 s_index = index_of(*p.spheres_, sphere);
			for (PVertex v : (*p.spheres_cluster_)[s_index])
			{
				uint32 v_idx = index_of(*data.mesh, v);
				if (v_idx >= nb_samples || marked[v_idx])
					continue;
				marked[v_idx] = 1;
				samples.push_back(v);
			}
		}

		for (PVertex sphere : candidate_spheres)
		{
			uint32 s_index = index_of(*p.spheres_, sphere);
			(*p.spheres_cluster_)[s_index].clear();
			(*p.spheres_cluster_area_)[s_index] = 0.0;
		}

		if (samples.empty())
			return;

		auto eval_distance = [&](uint32 v_index, uint32 sphere_index) -> Scalar {
			const Vec3& center = (*p.spheres_position_)[sphere_index];
			Scalar radius = (*p.spheres_radius_)[sphere_index];
			Scalar dist_sqem = (*data.quadric)[v_index].eval(Vec4(center.x(), center.y(), center.z(), radius));
			Scalar dist_other = (*data.line_quadric)[v_index].eval(center);
			Scalar dist = dist_sqem + p.sqem_clustering_lambda_ * dist_other;
			return dist;
		};

		for (PVertex v : samples)
		{
			uint32 v_index = index_of(*data.mesh, v);
			Scalar a = (*data.area)[v_index];
			Scalar min_distance = std::numeric_limits<Scalar>::max();
			PVertex closest_sphere;
			uint32 closest_sphere_index = 0;

			for (PVertex sphere : candidate_spheres)
			{
				uint32 s_index = index_of(*p.spheres_, sphere);
				Scalar dist = eval_distance(v_index, s_index);
				if (dist < min_distance)
				{
					min_distance = dist;
					closest_sphere = sphere;
					closest_sphere_index = s_index;
				}
			}

			if (!closest_sphere.is_valid())
				continue;
			(*data.sphere)[v_index] = closest_sphere;
			(*p.spheres_cluster_)[closest_sphere_index].push_back(v);
			(*p.spheres_cluster_area_)[closest_sphere_index] += a;
		}
	}

protected:
	void split_sphere(PointsParameters& p, PVertex sphere)
	{
		if (!sphere.is_valid())
			return;
		SphereFitData data;
		if (!get_sphere_fit_data(p, data))
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
			uint32 v_idx = index_of(*data.mesh, v);
			Scalar err = data.error ? (*data.error)[v_idx] : Scalar(0);
			if (err > max_err)
			{
				max_err = err;
				max_err_v = v;
			}
		}

		if (!max_err_v.is_valid())
			return;

		uint32 max_err_v_idx = index_of(*data.mesh, max_err_v);

		PVertex new_sphere = add_vertex(*p.spheres_);
		uint32 new_s_index = index_of(*p.spheres_, new_sphere);

		if (data.ma_position && data.ma_radius)
		{
			(*p.spheres_position_)[new_s_index] = (*data.ma_position)[max_err_v_idx];
			(*p.spheres_radius_)[new_s_index] = (*data.ma_radius)[max_err_v_idx];
		}
		else
		{
			(*p.spheres_position_)[new_s_index] = c;
			(*p.spheres_radius_)[new_s_index] = r;
		}
		(*p.spheres_parent_)[new_s_index] = sphere;
		(*p.spheres_cluster_color_)[new_s_index] =
			Vec4(0.5 + 0.5 * (rand() % 256) / 256.0, 0.5 + 0.5 * (rand() % 256) / 256.0,
				 0.5 + 0.5 * (rand() % 256) / 256.0, 1.0);

		p.nb_spheres_++;
		inherit_sphere_neighbors(p, sphere, new_sphere);
		recompute_clusters_local_neighborhood(p, new_sphere);
	}

	//------------------------------//
	//-----Topology correction------//
	//------------------------------//

	void compute_edge_degree(PointsParameters& p)
	{
		parallel_foreach_cell(*p.skeleton_, [&](NMEdge e) {
			auto in_face = incident_faces(*p.skeleton_, e);
			const uint32 ide = index_of(*p.skeleton_, e);
			const uint32 deg = static_cast<uint32>(in_face.size());

			if (ide != INVALID_INDEX && p.edge_degree_)
				(*p.edge_degree_)[ide] = deg;

			if (deg == 2)
				value<Vec3>(*p.skeleton_, p.skeleton_edge_color_, e) = Vec3(0.0, 1.0, 0.0);
			if (p.skeleton_edge_non_manifold_color_)
			{
				// Bright color for non-manifold edges, dim color for all others.
				value<Vec3>(*p.skeleton_, p.skeleton_edge_non_manifold_color_, e) =
					(deg > 2) ? Vec3(1.0, 0.2, 0.1) : Vec3(0.15, 0.15, 0.15);
			}

			return true;
		});
	}

	bool is_simple_face_3d(PointsParameters& p, const NMFace& f)
	{
		if (!f.is_valid())
			return false;
		uint32 idf = index_of(*p.skeleton_, f);
		if (idf == INVALID_INDEX)
			return false;
		return (*p.incident_tets_)[idf].size() == 1;
	}

	bool violates_face_removal_edge_guard(PointsParameters& p, const NMFace& f,
										  uint32 allow_degree0_edge_id = INVALID_INDEX)
	{
		if (!f.is_valid())
			return true;
		uint32 new_boundary_edge_count = 0;
		for (NMEdge e : incident_edges(*p.skeleton_, f))
		{
			if (!e.is_valid())
				continue;
			const uint32 ide = index_of(*p.skeleton_, e);
			if (ide == INVALID_INDEX)
				return true;
			const size_t deg = incident_faces(*p.skeleton_, e).size();
			// Deleting f would make this edge isolated.
			if (deg <= 1 && ide != allow_degree0_edge_id)
				return true;
			// Deleting f would turn this edge into boundary (degree 1).
			if (deg == 2)
			{
				++new_boundary_edge_count;
				if (new_boundary_edge_count >= 2)
					return true;
			}
		}
		return false;
	}

	uint32 remove_orphan_edges_from_removed_face_edges(PointsParameters& p, const std::vector<NMEdge>& face_edges)
	{
		uint32 removed_edges = 0;
		for (NMEdge e : face_edges)
		{
			if (!e.is_valid())
				continue;
			if (incident_faces(*p.skeleton_, e).empty())
			{
				remove_edge(*p.skeleton_, e);
				++removed_edges;
			}
		}
		return removed_edges;
	}

	bool compute_skeleton_face_udf_integrals_adaptive(PointsParameters& p, std::unordered_map<uint32, Scalar>& face_score_cache,
		uint32 max_adaptive_depth = 2, Scalar adaptive_abs_range_tol = Scalar(1e-4),
		Scalar adaptive_rel_range_tol = Scalar(0.35), const char* log_prefix = "[FaceUDF]")
	{
		if (!p.skeleton_ || !p.skeleton_position_ || !p.incident_tets_)
			return false;

		struct AdaptiveTriangleTask
		{
			uint32 face_id = INVALID_INDEX;
			Vec3 a;
			Vec3 b;
			Vec3 c;
			uint32 depth = 0;
		};

		static const std::array<std::array<Scalar, 3>, 7> dunavant7_bary = {{
			{{Scalar(1.0 / 3.0), Scalar(1.0 / 3.0), Scalar(1.0 / 3.0)}},
			{{Scalar(0.470142064105115), Scalar(0.470142064105115), Scalar(0.059715871789770)}},
			{{Scalar(0.470142064105115), Scalar(0.059715871789770), Scalar(0.470142064105115)}},
			{{Scalar(0.059715871789770), Scalar(0.470142064105115), Scalar(0.470142064105115)}},
			{{Scalar(0.101286507323456), Scalar(0.101286507323456), Scalar(0.797426985353087)}},
			{{Scalar(0.101286507323456), Scalar(0.797426985353087), Scalar(0.101286507323456)}},
			{{Scalar(0.797426985353087), Scalar(0.101286507323456), Scalar(0.101286507323456)}},
		}};
		static const std::array<Scalar, 7> dunavant7_w = {
			Scalar(0.225000000000000), Scalar(0.132394152788506), Scalar(0.132394152788506),
			Scalar(0.132394152788506), Scalar(0.125939180544827), Scalar(0.125939180544827),
			Scalar(0.125939180544827)};

		auto subdivide_triangle = [&](const AdaptiveTriangleTask& tri, std::vector<AdaptiveTriangleTask>& out) {
			const Vec3 ab = (tri.a + tri.b) * Scalar(0.5);
			const Vec3 bc = (tri.b + tri.c) * Scalar(0.5);
			const Vec3 ca = (tri.c + tri.a) * Scalar(0.5);
			const uint32 next_depth = tri.depth + 1;
			out.push_back({tri.face_id, tri.a, ab, ca, next_depth});
			out.push_back({tri.face_id, ab, tri.b, bc, next_depth});
			out.push_back({tri.face_id, ca, bc, tri.c, next_depth});
			out.push_back({tri.face_id, ab, bc, ca, next_depth});
		};

		face_score_cache.clear();
		face_score_cache.reserve(nb_cells<NMFace>(*p.skeleton_));
		std::vector<AdaptiveTriangleTask> pending_tris;
		pending_tris.reserve(nb_cells<NMFace>(*p.skeleton_) * 2);
		foreach_cell(*p.skeleton_, [&](NMFace f) -> bool {
			const uint32 idf = index_of(*p.skeleton_, f);
			if (idf == INVALID_INDEX)
				return true;
			const auto& in_tets = (*p.incident_tets_)[idf];
			// Only score faces that are currently adjacent to at least one tet.
			if (in_tets.empty())
				return true;
			face_score_cache[idf] = Scalar(0);
			if (!f.is_valid())
				return true;
			const std::vector<NMVertex> vertices = incident_vertices(*p.skeleton_, f);
			if (vertices.size() < 3)
				return true;
			const Vec3 p0 = value<Vec3>(*p.skeleton_, p.skeleton_position_, vertices[0]);
			for (uint32 i = 1; i + 1 < vertices.size(); ++i)
			{
				const Vec3 p1 = value<Vec3>(*p.skeleton_, p.skeleton_position_, vertices[i]);
				const Vec3 p2 = value<Vec3>(*p.skeleton_, p.skeleton_position_, vertices[i + 1]);
				if (!(geometry::area(p0, p1, p2) > Scalar(0)))
					continue;
				pending_tris.push_back({idf, p0, p1, p2, 0});
			}
			return true;
		});

		while (!pending_tris.empty())
		{
			std::vector<AdaptiveTriangleTask> next_tris;
			next_tris.reserve(pending_tris.size() * 2);

			struct TriBatchMeta
			{
				size_t tri_idx = 0;
				size_t sample_offset = 0;
				Scalar area = Scalar(0);
			};
			std::vector<TriBatchMeta> tri_batch;
			tri_batch.reserve(pending_tris.size());

			std::vector<Vec3> sample_points;
			sample_points.reserve(pending_tris.size() * 7);

			for (size_t tri_idx = 0; tri_idx < pending_tris.size(); ++tri_idx)
			{
				const AdaptiveTriangleTask& tri = pending_tris[tri_idx];
				const Scalar tri_area = geometry::area(tri.a, tri.b, tri.c);
				if (!(tri_area > Scalar(0)))
					continue;
				const size_t sample_offset = sample_points.size();
				for (uint32 k = 0; k < 7; ++k)
				{
					const Scalar l0 = dunavant7_bary[k][0];
					const Scalar l1 = dunavant7_bary[k][1];
					const Scalar l2 = dunavant7_bary[k][2];
					sample_points.push_back(tri.a * l0 + tri.b * l1 + tri.c * l2);
				}
				tri_batch.push_back({tri_idx, sample_offset, tri_area});
			}

			if (sample_points.empty())
			{
				pending_tris.clear();
				break;
			}

			std::vector<Scalar> score_values;
			if (!eval_topology_score_values(p, sample_points, score_values))
			{
				std::cerr << log_prefix << " face score field evaluation failed." << std::endl;
				return false;
			}

			for (const TriBatchMeta& tri_eval : tri_batch)
			{
				const AdaptiveTriangleTask& tri = pending_tris[tri_eval.tri_idx];
				Scalar weighted_mean = Scalar(0);
				Scalar min_udf = std::numeric_limits<Scalar>::max();
				Scalar max_udf = std::numeric_limits<Scalar>::lowest();
				for (uint32 k = 0; k < 7; ++k)
				{
					const Scalar u = score_values[tri_eval.sample_offset + k];
					weighted_mean += dunavant7_w[k] * u;
					min_udf = std::min(min_udf, u);
					max_udf = std::max(max_udf, u);
				}

				const Scalar tri_integral = tri_eval.area * weighted_mean;
				const Scalar udf_range = max_udf - min_udf;
				const Scalar refine_threshold =
					std::max(adaptive_abs_range_tol, adaptive_rel_range_tol * std::max(weighted_mean, Scalar(0)));
				const bool should_refine = (tri.depth < max_adaptive_depth) && (udf_range > refine_threshold);

				if (should_refine)
				{
					subdivide_triangle(tri, next_tris);
				}
				else
				{
					face_score_cache[tri.face_id] += tri_integral;
				}
			}
			pending_tris.swap(next_tris);
		}

		return true;
	}

	bool normalize_face_scores_by_area(PointsParameters& p, std::unordered_map<uint32, Scalar>& face_scores,
									   const char* log_prefix)
	{
		(void)log_prefix;
		if (!p.skeleton_ || !p.skeleton_position_)
			return false;
		foreach_cell(*p.skeleton_, [&](NMFace f) -> bool {
			const uint32 idf = index_of(*p.skeleton_, f);
			if (idf == INVALID_INDEX)
				return true;
			auto it_score = face_scores.find(idf);
			if (it_score == face_scores.end())
				return true;

			const std::vector<NMVertex> vertices = incident_vertices(*p.skeleton_, f);
			if (vertices.size() < 3)
			{
				it_score->second = Scalar(0);
				return true;
			}

			const Vec3 p0 = value<Vec3>(*p.skeleton_, p.skeleton_position_, vertices[0]);
			Scalar area = Scalar(0);
			for (uint32 i = 1; i + 1 < vertices.size(); ++i)
			{
				const Vec3 p1 = value<Vec3>(*p.skeleton_, p.skeleton_position_, vertices[i]);
				const Vec3 p2 = value<Vec3>(*p.skeleton_, p.skeleton_position_, vertices[i + 1]);
				area += geometry::area(p0, p1, p2);
			}

			if (area > Scalar(0))
				it_score->second /= area;
			else
			{
				it_score->second = Scalar(0);
			}
			return true;
		});
		return true;
	}

	bool compute_skeleton_face_scores(PointsParameters& p, std::unordered_map<uint32, Scalar>& face_scores,
									  bool normalize_by_area, const char* log_prefix)
	{
		if (!compute_skeleton_face_udf_integrals_adaptive(
				p, face_scores, 2, Scalar(1e-4), Scalar(0.35), log_prefix))
			return false;
		if (normalize_by_area)
		{
			if (!normalize_face_scores_by_area(p, face_scores, log_prefix))
				return false;
		}
		return true;
	}

	bool compute_skeleton_edge_scores_gauss3_normalized(PointsParameters& p,
														std::unordered_map<uint32, Scalar>& edge_scores,
														const char* log_prefix = "[EdgeUDF]")
	{
		if (!p.skeleton_ || !p.skeleton_position_)
			return false;

		static const Scalar gauss3_xi[3] = {Scalar(-0.7745966692414834), Scalar(0.0), Scalar(0.7745966692414834)};
		static const Scalar gauss3_w_ref[3] = {Scalar(0.5555555555555556), Scalar(0.8888888888888888),
											   Scalar(0.5555555555555556)};
		static const Scalar gauss3_w_01[3] = {Scalar(0.5) * gauss3_w_ref[0], Scalar(0.5) * gauss3_w_ref[1],
											  Scalar(0.5) * gauss3_w_ref[2]};

		struct EdgeBatchMeta
		{
			uint32 edge_id = INVALID_INDEX;
			size_t sample_offset = 0;
		};

		edge_scores.clear();
		edge_scores.reserve(nb_cells<NMEdge>(*p.skeleton_));

		// Only score edges that belong to current tet faces.
		std::unordered_set<uint32> tet_edge_ids;
		tet_edge_ids.reserve(nb_cells<NMEdge>(*p.skeleton_));
		for (const auto& kv : p.skeleton_tets_)
		{
			const Tet& tet = kv.second;
			for (uint32 i = 0; i < 4; ++i)
			{
				const NMFace f = tet.faces[i];
				if (!f.is_valid())
					continue;
				for (NMEdge e : incident_edges(*p.skeleton_, f))
				{
					if (!e.is_valid())
						continue;
					const uint32 ide = index_of(*p.skeleton_, e);
					if (ide != INVALID_INDEX)
						tet_edge_ids.insert(ide);
				}
			}
		}

		std::vector<EdgeBatchMeta> edge_batch;
		edge_batch.reserve(tet_edge_ids.size());
		std::vector<Vec3> sample_points;
		sample_points.reserve(tet_edge_ids.size() * 3);

		foreach_cell(*p.skeleton_, [&](NMEdge e) -> bool {
			if (!e.is_valid())
				return true;
			const uint32 ide = index_of(*p.skeleton_, e);
			if (ide == INVALID_INDEX)
				return true;
			if (tet_edge_ids.find(ide) == tet_edge_ids.end())
				return true;
			edge_scores[ide] = Scalar(0);

			const std::vector<NMVertex> verts = incident_vertices(*p.skeleton_, e);
			if (verts.size() != 2)
				return true;

			const Vec3 p0 = value<Vec3>(*p.skeleton_, p.skeleton_position_, verts[0]);
			const Vec3 p1 = value<Vec3>(*p.skeleton_, p.skeleton_position_, verts[1]);
			const Scalar edge_len = (p1 - p0).norm();
			if (!(edge_len > Scalar(0)))
				return true;

			const size_t sample_offset = sample_points.size();
			for (uint32 k = 0; k < 3; ++k)
			{
				const Scalar t = Scalar(0.5) * (gauss3_xi[k] + Scalar(1.0));
				sample_points.push_back(p0 * (Scalar(1.0) - t) + p1 * t);
			}
			edge_batch.push_back({ide, sample_offset});
			return true;
		});

		if (sample_points.empty())
			return true;

		std::vector<Scalar> score_values;
		if (!eval_topology_score_values(p, sample_points, score_values))
		{
			std::cerr << log_prefix << " edge score field evaluation failed." << std::endl;
			return false;
		}

		for (const EdgeBatchMeta& edge_eval : edge_batch)
		{
			Scalar avg_udf = Scalar(0);
			for (uint32 k = 0; k < 3; ++k)
				avg_udf += gauss3_w_01[k] * score_values[edge_eval.sample_offset + k];
			edge_scores[edge_eval.edge_id] = avg_udf;
		}
		std::cout << log_prefix << " scored_tet_edges=" << edge_scores.size()
				  << " total_edges=" << nb_cells<NMEdge>(*p.skeleton_) << std::endl;
		return true;
	}

	void visualize_skeleton_edge_udf_scores(PointsParameters& p)
	{
		if (!p.neural_udf_loaded_)
		{
			std::cerr << "Edge score visualization requires a loaded neural UDF model." << std::endl;
			return;
		}
		if (!p.skeleton_ || !p.skeleton_edge_udf_score_color_)
		{
			std::cerr << "Skeleton or edge score color attribute is not initialized." << std::endl;
			return;
		}

		std::unordered_map<uint32, Scalar> edge_scores;
		if (!compute_skeleton_edge_scores_gauss3_normalized(p, edge_scores, "[EdgeUDFColormap]"))
		{
			std::cerr << "Failed to compute edge UDF scores." << std::endl;
			return;
		}

		Scalar e_min = std::numeric_limits<Scalar>::max();
		Scalar e_max = std::numeric_limits<Scalar>::lowest();
		size_t edge_count = 0;
		size_t deg2_edge_count = 0;
		foreach_cell(*p.skeleton_, [&](NMEdge e) -> bool {
			const uint32 ide = index_of(*p.skeleton_, e);
			auto it = edge_scores.find(ide);
			if (it == edge_scores.end())
				return true;
			const Scalar v = it->second;
			e_min = std::min(e_min, v);
			e_max = std::max(e_max, v);
			++edge_count;
			if (incident_faces(*p.skeleton_, e).size() == 2)
				++deg2_edge_count;
			return true;
		});

		const bool has_range = (edge_count > 0) && (e_max > e_min);
		// Non-tet edges (not scored) stay in a low-visibility fixed color.
		const Vec3 fallback_color(0.08, 0.08, 0.08);
		foreach_cell(*p.skeleton_, [&](NMEdge e) -> bool {
			const uint32 ide = index_of(*p.skeleton_, e);
			Vec3 color = fallback_color;
			auto it = edge_scores.find(ide);
			if (it != edge_scores.end())
			{
				const Scalar v = it->second;
				if (has_range)
				{
					const Scalar vc = std::clamp(v, e_min, e_max);
					const Vec4 c = color_map(vc, e_min, e_max, 1.0f);
					color = Vec3(c[0], c[1], c[2]);
				}
				else
				{
					const Vec4 c = color_map(Scalar(0.5), Scalar(0), Scalar(1), 1.0f);
					color = Vec3(c[0], c[1], c[2]);
				}
			}
			value<Vec3>(*p.skeleton_, p.skeleton_edge_udf_score_color_, e) = color;
			return true;
		});

		non_manifold_provider_->emit_attribute_changed(*p.skeleton_, p.skeleton_edge_udf_score_color_.get());
		std::cout << "[EdgeUDFColormap] edges=" << edge_count << " deg2_edges=" << deg2_edge_count;
		if (edge_count > 0)
			std::cout << " clamp_min=" << e_min << " clamp_max=" << e_max;
		std::cout << " quadrature=gauss3 normalized_by_length=true" << std::endl;
	}

	enum class EdgeTetDeleteMode
	{
		SimpleTet,
		NonSimpleTet
	};

	void run_edge_score_tet_mode_topology_fix(PointsParameters& p, EdgeTetDeleteMode mode, bool single_step)
	{
		const char* log_tag = (mode == EdgeTetDeleteMode::SimpleTet) ? "[EdgeTetSimple]" : "[EdgeTetNonSimple]";
		if (!p.neural_udf_loaded_)
		{
			std::cerr << log_tag << " requires a loaded neural UDF model." << std::endl;
			return;
		}
		if (!p.skeleton_ || !p.incident_tets_)
		{
			std::cerr << log_tag << " requires a built skeleton." << std::endl;
			return;
		}

		std::unordered_map<uint32, Scalar> edge_score_cache;
		if (!compute_skeleton_edge_scores_gauss3_normalized(p, edge_score_cache, log_tag))
		{
			std::cerr << log_tag << " failed to compute edge scores." << std::endl;
			return;
		}
		std::unordered_map<uint32, Scalar> face_score_cache;
		if (!compute_skeleton_face_scores(
				p, face_score_cache, p.skeleton_face_score_normalize_by_area_, log_tag))
		{
			std::cerr << log_tag << " failed to compute face scores." << std::endl;
			return;
		}

		std::unordered_map<uint32, std::vector<std::size_t>> face_owner_tets;
		face_owner_tets.reserve(nb_cells<NMFace>(*p.skeleton_));
		std::unordered_map<std::size_t, uint32> deleted_face_count_per_tet;
		deleted_face_count_per_tet.reserve(p.skeleton_tets_.size());
		for (const auto& kv : p.skeleton_tets_)
		{
			const std::size_t tet_id = kv.first;
			const Tet& tet = kv.second;
			uint32 initial_deleted_faces = 0;
			for (uint32 i = 0; i < 4; ++i)
			{
				const NMFace f = tet.faces[i];
				if (!f.is_valid())
				{
					++initial_deleted_faces;
					continue;
				}
				const uint32 idf = index_of(*p.skeleton_, f);
				if (idf == INVALID_INDEX)
				{
					++initial_deleted_faces;
					continue;
				}
				face_owner_tets[idf].push_back(tet_id);
			}
			deleted_face_count_per_tet[tet_id] = initial_deleted_faces;
		}

		auto get_face_score = [&](uint32 idf) -> Scalar {
			auto it = face_score_cache.find(idf);
			return (it != face_score_cache.end()) ? it->second : Scalar(0);
		};
		auto has_face_budget = [&](uint32 face_id) -> bool {
			auto it_owner = face_owner_tets.find(face_id);
			if (it_owner == face_owner_tets.end())
				return false;
			for (std::size_t tet_id : it_owner->second)
			{
				if (p.skeleton_tets_.find(tet_id) == p.skeleton_tets_.end())
					continue;
				auto it_count = deleted_face_count_per_tet.find(tet_id);
				if (it_count != deleted_face_count_per_tet.end() && it_count->second >= 2)
					return false;
			}
			return true;
		};
		auto account_face_deletion_budget = [&](uint32 face_id) {
			auto it_owner = face_owner_tets.find(face_id);
			if (it_owner == face_owner_tets.end())
				return;
			for (std::size_t tet_id : it_owner->second)
			{
				if (p.skeleton_tets_.find(tet_id) == p.skeleton_tets_.end())
					continue;
				auto it_count = deleted_face_count_per_tet.find(tet_id);
				if (it_count != deleted_face_count_per_tet.end())
					++(it_count->second);
			}
		};
		auto face_belongs_to_tet = [&](const NMFace& f, std::size_t tet_id) -> bool {
			if (!f.is_valid())
				return false;
			const uint32 idf = index_of(*p.skeleton_, f);
			if (idf == INVALID_INDEX)
				return false;
			const auto& in_tets = (*p.incident_tets_)[idf];
			return in_tets.find(tet_id) != in_tets.end();
		};
		auto face_has_degree1_edge = [&](const NMFace& f) -> bool {
			if (!f.is_valid())
				return false;
			for (NMEdge e : incident_edges(*p.skeleton_, f))
			{
				if (!e.is_valid())
					continue;
				if (incident_faces(*p.skeleton_, e).size() == 1)
					return true;
			}
			return false;
		};
		auto face_has_edge_with_tet_face_count_gt2 =
			[&](const NMFace& f, const std::unordered_set<std::size_t>* ignored_tets) -> bool {
			if (!f.is_valid())
				return false;
			for (NMEdge e : incident_edges(*p.skeleton_, f))
			{
				if (!e.is_valid())
					continue;
				size_t tet_face_count = 0;
				for (NMFace ef : incident_faces(*p.skeleton_, e))
				{
					if (!ef.is_valid())
						continue;
					const uint32 idef = index_of(*p.skeleton_, ef);
					if (idef == INVALID_INDEX)
						continue;
					bool has_alive_owner = false;
					for (std::size_t tet_id : (*p.incident_tets_)[idef])
					{
						if (p.skeleton_tets_.find(tet_id) == p.skeleton_tets_.end())
							continue;
						if (ignored_tets && ignored_tets->find(tet_id) != ignored_tets->end())
							continue;
						has_alive_owner = true;
						break;
					}
					if (has_alive_owner)
						++tet_face_count;
				}
				if (tet_face_count > 2)
					return true;
			}
			return false;
		};
		auto tet_has_simple_face = [&](std::size_t tet_id) -> bool {
			auto it_tet = p.skeleton_tets_.find(tet_id);
			if (it_tet == p.skeleton_tets_.end())
				return false;
			const Tet& tet = it_tet->second;
			for (uint32 i = 0; i < 4; ++i)
			{
				const NMFace f = tet.faces[i];
				if (!f.is_valid())
					continue;
				const uint32 idf = index_of(*p.skeleton_, f);
				if (idf == INVALID_INDEX)
					continue;
				const auto& in_tets = (*p.incident_tets_)[idf];
				if (in_tets.find(tet_id) == in_tets.end())
					continue;
				if (in_tets.size() == 1)
					return true;
			}
			return false;
		};
		auto pick_tet_best_edge = [&](std::size_t tet_id, NMEdge& out_edge, uint32& out_edge_id, Scalar& out_score) -> bool {
			auto it_tet = p.skeleton_tets_.find(tet_id);
			if (it_tet == p.skeleton_tets_.end())
				return false;
			const Tet& tet = it_tet->second;
			std::unordered_set<uint32> seen_edges;
			seen_edges.reserve(8);
			bool found = false;
			for (uint32 i = 0; i < 4; ++i)
			{
				const NMFace f = tet.faces[i];
				if (!face_belongs_to_tet(f, tet_id))
					continue;
				for (NMEdge e : incident_edges(*p.skeleton_, f))
				{
					if (!e.is_valid())
						continue;
					const uint32 ide = index_of(*p.skeleton_, e);
					if (ide == INVALID_INDEX)
						continue;
					if (!seen_edges.insert(ide).second)
						continue;
					const auto it_score = edge_score_cache.find(ide);
					const Scalar score = (it_score != edge_score_cache.end()) ? it_score->second : Scalar(0);
					if (!found || score > out_score)
					{
						found = true;
						out_edge = e;
						out_edge_id = ide;
						out_score = score;
					}
				}
			}
			return found;
		};
		auto pick_first_face = [&](std::size_t tet_id, const NMEdge& e, NMFace& out_face, uint32& out_face_id,
								  Scalar& out_face_score) -> bool {
			const auto in_faces = incident_faces(*p.skeleton_, e);
			bool found = false;
			for (const NMFace& f : in_faces)
			{
				if (!face_belongs_to_tet(f, tet_id))
					continue;
				const uint32 idf = index_of(*p.skeleton_, f);
				if (idf == INVALID_INDEX || !has_face_budget(idf))
					continue;
				if (mode == EdgeTetDeleteMode::SimpleTet && (*p.incident_tets_)[idf].size() != 1)
					continue;
				const Scalar fs = get_face_score(idf);
				if (!found || fs > out_face_score)
				{
					found = true;
					out_face = f;
					out_face_id = idf;
					out_face_score = fs;
				}
			}
			return found;
		};
		auto remove_faces_collect_tets =
			[&](const std::vector<NMFace>& faces_to_remove, std::unordered_set<std::size_t>& out_tets_to_erase,
				uint32& out_removed_faces, uint32& out_removed_edges) -> bool {
			out_tets_to_erase.clear();
			std::vector<std::pair<uint32, NMFace>> unique_faces;
			unique_faces.reserve(faces_to_remove.size());
			std::unordered_set<uint32> seen_face_ids;
			seen_face_ids.reserve(faces_to_remove.size() * 2 + 1);
			for (const NMFace& f : faces_to_remove)
			{
				if (!f.is_valid())
					continue;
				const uint32 idf = index_of(*p.skeleton_, f);
				if (idf == INVALID_INDEX)
					continue;
				if (!seen_face_ids.insert(idf).second)
					continue;
				unique_faces.push_back({idf, f});
			}
			if (unique_faces.empty())
				return false;

			std::unordered_set<std::size_t> tets_to_remove_local;
			tets_to_remove_local.reserve(16);
			for (const auto& fpair : unique_faces)
			{
				const uint32 idf = fpair.first;
				for (std::size_t tet_id : (*p.incident_tets_)[idf])
					if (p.skeleton_tets_.find(tet_id) != p.skeleton_tets_.end())
						tets_to_remove_local.insert(tet_id);
			}
			if (tets_to_remove_local.empty())
				return false;

			std::vector<NMEdge> affected_edges;
			affected_edges.reserve(unique_faces.size() * 3);
			for (const auto& fpair : unique_faces)
			{
				const NMFace f = fpair.second;
				for (NMEdge e : incident_edges(*p.skeleton_, f))
					affected_edges.push_back(e);
			}

			for (const auto& fpair : unique_faces)
			{
				const NMFace f = fpair.second;
				if (!f.is_valid())
					continue;
				const uint32 idf_now = index_of(*p.skeleton_, f);
				if (idf_now == INVALID_INDEX)
					continue;
				account_face_deletion_budget(idf_now);
				remove_face(*p.skeleton_, f);
				++out_removed_faces;
			}

			out_removed_edges += remove_orphan_edges_from_removed_face_edges(p, affected_edges);
			out_tets_to_erase = std::move(tets_to_remove_local);
			return true;
		};
		auto finalize_erase_tets = [&](const std::unordered_set<std::size_t>& tets_to_erase, uint32& out_removed_tets) {
			for (std::size_t tet_id : tets_to_erase)
			{
				auto it_tet = p.skeleton_tets_.find(tet_id);
				if (it_tet == p.skeleton_tets_.end())
					continue;
				const Tet old_tet = it_tet->second;
				for (uint32 i = 0; i < 4; ++i)
				{
					const NMFace f = old_tet.faces[i];
					if (!f.is_valid())
						continue;
					const uint32 idf = index_of(*p.skeleton_, f);
					if (idf == INVALID_INDEX)
						continue;
					(*p.incident_tets_)[idf].erase(tet_id);
				}
				p.skeleton_tets_.erase(it_tet);
				++out_removed_tets;
			}
		};
		struct StepCandidate
		{
			std::size_t tet_id = 0;
			NMEdge edge;
			uint32 edge_id = INVALID_INDEX;
			Scalar edge_score = Scalar(0);
			NMFace first_face;
			uint32 first_face_id = INVALID_INDEX;
			Scalar first_face_score = Scalar(0);
		};
		struct FaceChoice
		{
			NMFace face;
			uint32 face_id = INVALID_INDEX;
			Scalar score = Scalar(0);
		};

		uint32 removed_faces = 0;
		uint32 removed_edges = 0;
		uint32 removed_tets = 0;
		uint32 processed_steps = 0;
		uint32 skipped_no_candidate = 0;

		while (true)
		{
			StepCandidate best;
			bool has_best = false;
			for (const auto& kv : p.skeleton_tets_)
			{
				const std::size_t tet_id = kv.first;
				const bool is_simple_tet = tet_has_simple_face(tet_id);
				if (mode == EdgeTetDeleteMode::SimpleTet && !is_simple_tet)
					continue;
				if (mode == EdgeTetDeleteMode::NonSimpleTet && is_simple_tet)
					continue;

				NMEdge best_edge_local;
				uint32 best_edge_local_id = INVALID_INDEX;
				Scalar best_edge_local_score = Scalar(0);
				if (!pick_tet_best_edge(tet_id, best_edge_local, best_edge_local_id, best_edge_local_score))
					continue;

				NMFace first_face_local;
				uint32 first_face_local_id = INVALID_INDEX;
				Scalar first_face_local_score = Scalar(0);
				if (!pick_first_face(tet_id, best_edge_local, first_face_local, first_face_local_id, first_face_local_score))
					continue;

				const Scalar edge_eps = Scalar(1e-12);
				const Scalar face_eps = Scalar(1e-12);
				const bool better_edge = (!has_best) || (best_edge_local_score > best.edge_score + edge_eps);
				const bool tie_edge = has_best && (std::abs(best_edge_local_score - best.edge_score) <= edge_eps);
				const bool better_face_on_tie = tie_edge && (first_face_local_score > best.first_face_score + face_eps);
				const bool tie_face = tie_edge && (std::abs(first_face_local_score - best.first_face_score) <= face_eps);
				const bool better_id_on_full_tie = tie_face && (tet_id < best.tet_id);
				if (better_edge || better_face_on_tie || better_id_on_full_tie)
				{
					has_best = true;
					best.tet_id = tet_id;
					best.edge = best_edge_local;
					best.edge_id = best_edge_local_id;
					best.edge_score = best_edge_local_score;
					best.first_face = first_face_local;
					best.first_face_id = first_face_local_id;
					best.first_face_score = first_face_local_score;
				}
			}

			if (!has_best)
			{
				++skipped_no_candidate;
				break;
			}

			std::unordered_set<std::size_t> first_removed_tets;
			for (std::size_t tet_id : (*p.incident_tets_)[best.first_face_id])
				if (p.skeleton_tets_.find(tet_id) != p.skeleton_tets_.end())
					first_removed_tets.insert(tet_id);

			uint32 step_removed_faces = 0;
			uint32 step_removed_edges = 0;
			uint32 step_removed_tets = 0;
			std::unordered_set<std::size_t> step_tets_to_erase;
			std::unordered_set<std::size_t> first_tets_to_erase;
			if (!remove_faces_collect_tets({best.first_face}, first_tets_to_erase, step_removed_faces, step_removed_edges))
			{
				std::cout << log_tag << " skip tet=" << best.tet_id << " edge=" << best.edge_id
						  << " first_face=" << best.first_face_id << " reason=first_face_remove_failed" << std::endl;
				continue;
			}
			step_tets_to_erase.insert(first_tets_to_erase.begin(), first_tets_to_erase.end());

			std::vector<FaceChoice> created_deg1_faces;
			created_deg1_faces.reserve(16);
			std::unordered_set<uint32> seen_created_ids;
			seen_created_ids.reserve(32);
			uint32 second_scan_total = 0;
			uint32 second_scan_skip_first_face = 0;
			uint32 second_scan_skip_invalid_index = 0;
			uint32 second_scan_skip_duplicate = 0;
			uint32 second_scan_skip_no_budget = 0;
			uint32 second_scan_skip_not_deg1 = 0;
			uint32 second_scan_skip_edge_tet_face_gt2 = 0;
			for (std::size_t tet_id : first_removed_tets)
			{
				auto it_tet = p.skeleton_tets_.find(tet_id);
				if (it_tet == p.skeleton_tets_.end())
					continue;
				const Tet& tet = it_tet->second;
				for (uint32 i = 0; i < 4; ++i)
				{
					const NMFace f = tet.faces[i];
					const uint32 idf = index_of(*p.skeleton_, f);
					if (idf == INVALID_INDEX)
					{
						++second_scan_skip_invalid_index;
						continue;
					}
					if (idf == best.first_face_id)
					{
						++second_scan_skip_first_face;
						continue;
					}
					++second_scan_total;
					if (!seen_created_ids.insert(idf).second)
					{
						++second_scan_skip_duplicate;
						continue;
					}
					if (!has_face_budget(idf))
					{
						++second_scan_skip_no_budget;
						continue;
					}
					if (!face_has_degree1_edge(f))
					{
						++second_scan_skip_not_deg1;
						continue;
					}
					if (face_has_edge_with_tet_face_count_gt2(f, &step_tets_to_erase))
					{
						++second_scan_skip_edge_tet_face_gt2;
						continue;
					}
					created_deg1_faces.push_back({f, idf, get_face_score(idf)});
				}
			}
			std::sort(created_deg1_faces.begin(), created_deg1_faces.end(),
					  [](const FaceChoice& a, const FaceChoice& b) { return a.score > b.score; });
			const uint32 second_candidates_before_mode_cap = static_cast<uint32>(created_deg1_faces.size());
			uint32 second_candidates_trimmed_by_simple_mode = 0;
			if (mode == EdgeTetDeleteMode::SimpleTet && created_deg1_faces.size() > 1)
			{
				second_candidates_trimmed_by_simple_mode = static_cast<uint32>(created_deg1_faces.size() - 1);
				created_deg1_faces.resize(1);
			}

			uint32 second_faces_removed = 0;
			uint32 second_exec_skip_invalid_face = 0;
			uint32 second_exec_skip_invalid_index = 0;
			uint32 second_exec_skip_no_budget = 0;
			uint32 second_exec_skip_not_deg1 = 0;
			uint32 second_exec_skip_edge_tet_face_gt2 = 0;
			uint32 second_exec_remove_failed = 0;
			for (const FaceChoice& fc : created_deg1_faces)
			{
				if (!fc.face.is_valid())
				{
					++second_exec_skip_invalid_face;
					continue;
				}
				const uint32 idf_now = index_of(*p.skeleton_, fc.face);
				if (idf_now == INVALID_INDEX)
				{
					++second_exec_skip_invalid_index;
					continue;
				}
				if (!has_face_budget(idf_now))
				{
					++second_exec_skip_no_budget;
					continue;
				}
				if (!face_has_degree1_edge(fc.face))
				{
					++second_exec_skip_not_deg1;
					continue;
				}
				if (face_has_edge_with_tet_face_count_gt2(fc.face, &step_tets_to_erase))
				{
					++second_exec_skip_edge_tet_face_gt2;
					continue;
				}
				std::unordered_set<std::size_t> second_tets_to_erase;
				uint32 rf = 0, re = 0;
				if (remove_faces_collect_tets({fc.face}, second_tets_to_erase, rf, re))
				{
					step_removed_faces += rf;
					step_removed_edges += re;
					step_tets_to_erase.insert(second_tets_to_erase.begin(), second_tets_to_erase.end());
					++second_faces_removed;
				}
				else
				{
					++second_exec_remove_failed;
				}
			}
			finalize_erase_tets(step_tets_to_erase, step_removed_tets);

			removed_faces += step_removed_faces;
			removed_edges += step_removed_edges;
			removed_tets += step_removed_tets;
			++processed_steps;
			std::cout << log_tag << " remove#" << processed_steps
					  << " tet=" << best.tet_id
					  << " edge=" << best.edge_id
					  << " edge_score=" << best.edge_score
					  << " first_face=" << best.first_face_id
					  << " first_face_score=" << best.first_face_score
					  << " second_faces_removed=" << second_faces_removed
					  << " second_scan_total=" << second_scan_total
					  << " second_scan_candidates=" << second_candidates_before_mode_cap
					  << " second_scan_trimmed_simple=" << second_candidates_trimmed_by_simple_mode
					  << " second_scan_skip_first_face=" << second_scan_skip_first_face
					  << " second_scan_skip_invalid_index=" << second_scan_skip_invalid_index
					  << " second_scan_skip_duplicate=" << second_scan_skip_duplicate
					  << " second_scan_skip_no_budget=" << second_scan_skip_no_budget
					  << " second_scan_skip_not_deg1=" << second_scan_skip_not_deg1
					  << " second_scan_skip_edge_tet_face_gt2=" << second_scan_skip_edge_tet_face_gt2
					  << " second_exec_skip_invalid_face=" << second_exec_skip_invalid_face
					  << " second_exec_skip_invalid_index=" << second_exec_skip_invalid_index
					  << " second_exec_skip_no_budget=" << second_exec_skip_no_budget
					  << " second_exec_skip_not_deg1=" << second_exec_skip_not_deg1
					  << " second_exec_skip_edge_tet_face_gt2=" << second_exec_skip_edge_tet_face_gt2
					  << " second_exec_remove_failed=" << second_exec_remove_failed
					  << " removed_faces=" << step_removed_faces
					  << " removed_edges=" << step_removed_edges
					  << " removed_tets=" << step_removed_tets
					  << " remaining_tets=" << p.skeleton_tets_.size() << std::endl;

			if (single_step)
				break;
		}

		if (single_step && processed_steps == 0)
			std::cout << log_tag << " no removable tet in current mode." << std::endl;
		std::cout << log_tag << " done"
				  << " steps=" << processed_steps
				  << " removed_faces=" << removed_faces
				  << " removed_edges=" << removed_edges
				  << " removed_tets=" << removed_tets
				  << " skipped_no_candidate=" << skipped_no_candidate
				  << " remaining_tets=" << p.skeleton_tets_.size() << std::endl;

		refresh_skeleton_topology_colors(p);
		mark_boundary_tets_color(p);
		visualize_skeleton_edge_udf_scores(p);
	}

	void run_edge_score_simple_tet_topology_fix(PointsParameters& p)
	{
		run_edge_score_tet_mode_topology_fix(p, EdgeTetDeleteMode::SimpleTet, false);
	}

	void run_edge_score_simple_tet_topology_fix_single_step(PointsParameters& p)
	{
		run_edge_score_tet_mode_topology_fix(p, EdgeTetDeleteMode::SimpleTet, true);
	}

	void run_edge_score_nonsimple_tet_topology_fix(PointsParameters& p)
	{
		run_edge_score_tet_mode_topology_fix(p, EdgeTetDeleteMode::NonSimpleTet, false);
	}

	void run_edge_score_nonsimple_tet_topology_fix_single_step(PointsParameters& p)
	{
		run_edge_score_tet_mode_topology_fix(p, EdgeTetDeleteMode::NonSimpleTet, true);
	}

	void visualize_skeleton_tet_face_udf(PointsParameters& p, bool normalize_by_area = false)
	{
		if (!p.neural_udf_loaded_)
		{
			std::cerr << "UDF face visualization requires a loaded neural UDF model." << std::endl;
			return;
		}
		if (!p.skeleton_ || !p.incident_tets_ || !p.skeleton_face_udf_color_)
		{
			std::cerr << "Skeleton or face attributes are not initialized." << std::endl;
			return;
		}

		std::unordered_map<uint32, Scalar> face_udf_integrals;
		if (!compute_skeleton_face_scores(p, face_udf_integrals, normalize_by_area, "[UDFColormap]"))
		{
			std::cerr << "Failed to compute face UDF scores for visualization." << std::endl;
			return;
		}

		const Vec3 non_tet_face_color(0.2, 0.2, 0.2);
		Scalar tri_min = std::numeric_limits<Scalar>::max();
		Scalar tri_max = std::numeric_limits<Scalar>::lowest();
		size_t tet_face_count = 0;
		size_t non_tet_face_count = 0;

		foreach_cell(*p.skeleton_, [&](NMFace f) -> bool {
			const uint32 idf = index_of(*p.skeleton_, f);
			if (idf == INVALID_INDEX)
				return true;
			const auto& in_tets = (*p.incident_tets_)[idf];
			if (in_tets.empty())
			{
				++non_tet_face_count;
				return true;
			}
			++tet_face_count;
			const auto it = face_udf_integrals.find(idf);
			const Scalar v = (it != face_udf_integrals.end()) ? it->second : Scalar(0);
			tri_min = std::min(tri_min, v);
			tri_max = std::max(tri_max, v);
			return true;
		});

		const bool has_tet_faces = (tet_face_count > 0);
		const bool has_range = has_tet_faces && (tri_max > tri_min);
		foreach_cell(*p.skeleton_, [&](NMFace f) -> bool {
			const uint32 idf = index_of(*p.skeleton_, f);
			const auto& in_tets = (*p.incident_tets_)[idf];
			Vec3 color = non_tet_face_color;
			if (!in_tets.empty())
			{
				const auto it = face_udf_integrals.find(idf);
				const Scalar v = (it != face_udf_integrals.end()) ? it->second : Scalar(0);
				if (has_range)
				{
					const Scalar vc = std::clamp(v, tri_min, tri_max);
					const Vec4 c = color_map(vc, tri_min, tri_max, 1.0f);
					color = Vec3(c[0], c[1], c[2]);
				}
				else
				{
					const Vec4 c = color_map(Scalar(0.5), Scalar(0), Scalar(1), 1.0f);
					color = Vec3(c[0], c[1], c[2]);
				}
			}
			value<Vec3>(*p.skeleton_, p.skeleton_face_udf_color_, f) = color;
			return true;
		});

		non_manifold_provider_->emit_attribute_changed(*p.skeleton_, p.skeleton_face_udf_color_.get());
		std::cout << "[UDFColormap] tet_faces=" << tet_face_count << " non_tet_faces=" << non_tet_face_count;
		if (has_tet_faces)
			std::cout << " clamp_min=" << tri_min << " clamp_max=" << tri_max;
		std::cout << " normalized_by_area=" << (normalize_by_area ? "true" : "false");
		std::cout << std::endl;
	}

	struct K5DetectionResult
	{
		std::vector<std::array<uint32, 5>> cliques;
		std::unordered_set<uint32> face_ids;
		std::unordered_set<uint32> edge_ids;
	};

	K5DetectionResult detect_k5_cells(PointsParameters& p)
	{
		K5DetectionResult result;
		if (!p.skeleton_)
			return result;

		auto make_edge_key = [](uint32 a, uint32 b) -> uint64 {
			if (a > b)
				std::swap(a, b);
			return (uint64(a) << 32) | uint64(b);
		};
		auto make_face_key = [](uint32 a, uint32 b, uint32 c) -> std::array<uint32, 3> {
			std::array<uint32, 3> key = {a, b, c};
			std::sort(key.begin(), key.end());
			return key;
		};

		std::vector<uint32> vertex_ids;
		vertex_ids.reserve(nb_cells<NMVertex>(*p.skeleton_));
		std::unordered_map<uint32, std::unordered_set<uint32>> adjacency;
		adjacency.reserve(nb_cells<NMVertex>(*p.skeleton_));
		std::unordered_map<uint64, uint32> edge_pair_to_id;
		edge_pair_to_id.reserve(nb_cells<NMEdge>(*p.skeleton_) * 2 + 1);

		foreach_cell(*p.skeleton_, [&](NMVertex v) -> bool {
			const uint32 idv = index_of(*p.skeleton_, v);
			if (idv != INVALID_INDEX)
				vertex_ids.push_back(idv);
			return true;
		});
		std::sort(vertex_ids.begin(), vertex_ids.end());
		vertex_ids.erase(std::unique(vertex_ids.begin(), vertex_ids.end()), vertex_ids.end());

		foreach_cell(*p.skeleton_, [&](NMEdge e) -> bool {
			const uint32 ide = index_of(*p.skeleton_, e);
			if (ide == INVALID_INDEX)
				return true;
			const std::vector<NMVertex> vv = incident_vertices(*p.skeleton_, e);
			if (vv.size() != 2)
				return true;
			const uint32 a = index_of(*p.skeleton_, vv[0]);
			const uint32 b = index_of(*p.skeleton_, vv[1]);
			if (a == INVALID_INDEX || b == INVALID_INDEX || a == b)
				return true;
			adjacency[a].insert(b);
			adjacency[b].insert(a);
			edge_pair_to_id[make_edge_key(a, b)] = ide;
			return true;
		});

		if (vertex_ids.size() < 5)
			return result;

		auto has_edge = [&](uint32 a, uint32 b) -> bool {
			auto it = adjacency.find(a);
			if (it == adjacency.end())
				return false;
			return it->second.find(b) != it->second.end();
		};

		const size_t n = vertex_ids.size();
		for (size_t ia = 0; ia < n; ++ia)
		{
			const uint32 a = vertex_ids[ia];
			for (size_t ib = ia + 1; ib < n; ++ib)
			{
				const uint32 b = vertex_ids[ib];
				if (!has_edge(a, b))
					continue;
				for (size_t ic = ib + 1; ic < n; ++ic)
				{
					const uint32 c = vertex_ids[ic];
					if (!has_edge(a, c) || !has_edge(b, c))
						continue;
					for (size_t id = ic + 1; id < n; ++id)
					{
						const uint32 d = vertex_ids[id];
						if (!has_edge(a, d) || !has_edge(b, d) || !has_edge(c, d))
							continue;
						for (size_t ie = id + 1; ie < n; ++ie)
						{
							const uint32 e = vertex_ids[ie];
							if (!has_edge(a, e) || !has_edge(b, e) || !has_edge(c, e) || !has_edge(d, e))
								continue;

							const std::array<uint32, 5> clique = {a, b, c, d, e};
							result.cliques.push_back(clique);

							for (uint32 i = 0; i < 5; ++i)
							{
								for (uint32 j = i + 1; j < 5; ++j)
								{
									const auto it_eid = edge_pair_to_id.find(make_edge_key(clique[i], clique[j]));
									if (it_eid != edge_pair_to_id.end())
										result.edge_ids.insert(it_eid->second);
								}
							}
						}
					}
				}
			}
		}

		std::map<std::array<uint32, 3>, uint32> face_key_to_id;
		face_key_to_id.clear();
		foreach_cell(*p.skeleton_, [&](NMFace f) -> bool {
			const uint32 idf = index_of(*p.skeleton_, f);
			if (idf == INVALID_INDEX)
				return true;
			const std::vector<NMVertex> fv = incident_vertices(*p.skeleton_, f);
			if (fv.size() != 3)
				return true;
			const uint32 a = index_of(*p.skeleton_, fv[0]);
			const uint32 b = index_of(*p.skeleton_, fv[1]);
			const uint32 c = index_of(*p.skeleton_, fv[2]);
			if (a == INVALID_INDEX || b == INVALID_INDEX || c == INVALID_INDEX)
				return true;
			face_key_to_id[make_face_key(a, b, c)] = idf;
			return true;
		});

		for (const auto& clique : result.cliques)
		{
			for (uint32 i = 0; i < 5; ++i)
			{
				for (uint32 j = i + 1; j < 5; ++j)
				{
					for (uint32 k = j + 1; k < 5; ++k)
					{
						const auto it_f = face_key_to_id.find(make_face_key(clique[i], clique[j], clique[k]));
						if (it_f != face_key_to_id.end())
							result.face_ids.insert(it_f->second);
					}
				}
			}
		}

		return result;
	}

	void mark_boundary_tets_color(PointsParameters& p)
	{
		if (!p.skeleton_ || !p.skeleton_face_boundary_tet_color_ || !p.skeleton_edge_boundary_tet_color_)
		{
			std::cerr << "Boundary tet coloring requires skeleton face/edge color attributes." << std::endl;
			return;
		}

		const Vec3 boundary_face_color(1.0, 0.85, 0.1);
		const Vec3 k5_face_color(0.05, 0.9, 0.95);
		const Vec3 non_boundary_face_color(0.12, 0.12, 0.12);
		const Vec3 boundary_best_edge_color(1.0, 0.45, 0.05);
		const Vec3 k5_edge_color(0.0, 0.85, 0.85);
		const Vec3 overlap_edge_color(0.95, 0.95, 0.2);
		const Vec3 non_boundary_best_edge_color(0.1, 0.1, 0.1);
		std::unordered_set<uint32> boundary_face_ids;
		boundary_face_ids.reserve(nb_cells<NMFace>(*p.skeleton_));
		std::unordered_set<uint32> boundary_best_edge_ids;
		boundary_best_edge_ids.reserve(nb_cells<NMEdge>(*p.skeleton_));

		uint32 boundary_tet_count = 0;
		for (const auto& kv : p.skeleton_tets_)
		{
			const Tet& tet = kv.second;
			bool is_boundary_tet = false;
			std::unordered_set<uint32> tet_face_ids;
			tet_face_ids.reserve(4);
			for (uint32 i = 0; i < 4; ++i)
			{
				const NMFace f = tet.faces[i];
				const uint32 idf = index_of(*p.skeleton_, f);

				tet_face_ids.insert(idf);

				uint32 deg2_edge_count = 0;
				for (NMEdge e : incident_edges(*p.skeleton_, f))
				{
					if (incident_faces(*p.skeleton_, e).size() == 2)
						++deg2_edge_count;
				}
				if (deg2_edge_count >= 2)
				{
					is_boundary_tet = true;
					break;
				}
			}
			if (!is_boundary_tet)
				continue;

			struct BoundaryEdgeChoice
			{
				NMEdge edge;
				Scalar length = Scalar(-1);
			};
			BoundaryEdgeChoice best_choice;
			std::unordered_set<uint32> seen_edges;
			seen_edges.reserve(8);
			for (uint32 i = 0; i < 4; ++i)
			{
				const NMFace f = tet.faces[i];
				for (NMEdge e : incident_edges(*p.skeleton_, f))
				{
					const uint32 ide = index_of(*p.skeleton_, e);
					if (!seen_edges.insert(ide).second)
						continue;
					const auto in_faces = incident_faces(*p.skeleton_, e);
					if (in_faces.size() != 2)
						continue;
					const NMFace f0 = in_faces[0];
					const NMFace f1 = in_faces[1];
					const uint32 idf0 = index_of(*p.skeleton_, f0);
					const uint32 idf1 = index_of(*p.skeleton_, f1);
					const std::vector<NMVertex> edge_vertices = incident_vertices(*p.skeleton_, e);
					const Vec3 p0 = value<Vec3>(*p.skeleton_, p.skeleton_position_, edge_vertices[0]);
					const Vec3 p1 = value<Vec3>(*p.skeleton_, p.skeleton_position_, edge_vertices[1]);
					const Scalar edge_len = (p1 - p0).norm();
					if (edge_len > best_choice.length)
					{
						best_choice.edge = e;
						best_choice.length = edge_len;
					}
				}
			}
			if (best_choice.length >= Scalar(0))
			{
				const uint32 ide = index_of(*p.skeleton_, best_choice.edge);
				boundary_best_edge_ids.insert(ide);
			}

			++boundary_tet_count;
			for (uint32 i = 0; i < 4; ++i)
			{
				const NMFace f = tet.faces[i];

				const uint32 idf = index_of(*p.skeleton_, f);
				boundary_face_ids.insert(idf);
			}
		}

		const K5DetectionResult k5_info = detect_k5_cells(p);

		foreach_cell(*p.skeleton_, [&](NMFace f) -> bool {
			const uint32 idf = index_of(*p.skeleton_, f);
			const bool is_boundary_face = boundary_face_ids.find(idf) != boundary_face_ids.end();
			const bool is_k5_face = k5_info.face_ids.find(idf) != k5_info.face_ids.end();
			Vec3 c = non_boundary_face_color;
			if (is_k5_face)
				c = k5_face_color;
			else if (is_boundary_face)
				c = boundary_face_color;
			value<Vec3>(*p.skeleton_, p.skeleton_face_boundary_tet_color_, f) = c;
			return true;
		});
		foreach_cell(*p.skeleton_, [&](NMEdge e) -> bool {
			const uint32 ide = index_of(*p.skeleton_, e);
			const bool is_best_edge = boundary_best_edge_ids.find(ide) != boundary_best_edge_ids.end();
			const bool is_k5_edge = k5_info.edge_ids.find(ide) != k5_info.edge_ids.end();
			Vec3 c = non_boundary_best_edge_color;
			if (is_best_edge && is_k5_edge)
				c = overlap_edge_color;
			else if (is_best_edge)
				c = boundary_best_edge_color;
			else if (is_k5_edge)
				c = k5_edge_color;
			value<Vec3>(*p.skeleton_, p.skeleton_edge_boundary_tet_color_, e) = c;
			return true;
		});

		non_manifold_provider_->emit_attribute_changed(*p.skeleton_, p.skeleton_face_boundary_tet_color_.get());
		non_manifold_provider_->emit_attribute_changed(*p.skeleton_, p.skeleton_edge_boundary_tet_color_.get());

		std::cout << "[BoundaryTetMask] boundary_tets=" << boundary_tet_count
				  << " boundary_faces=" << boundary_face_ids.size()
				  << " boundary_best_edges=" << boundary_best_edge_ids.size()
				  << " k5_count=" << k5_info.cliques.size()
				  << " k5_faces=" << k5_info.face_ids.size()
				  << " k5_edges=" << k5_info.edge_ids.size()
				  << " total_faces=" << nb_cells<NMFace>(*p.skeleton_) << std::endl;
	}

	void run_k5_face_deletion(PointsParameters& p)
	{
		if (!p.skeleton_ || !p.incident_tets_)
		{
			std::cerr << "[K5Delete] requires a built skeleton." << std::endl;
			return;
		}
		if (!p.neural_udf_loaded_)
		{
			std::cerr << "[K5Delete] requires a loaded neural UDF model for face score ranking." << std::endl;
			return;
		}

		std::unordered_map<uint32, Scalar> face_score_cache;
		if (!compute_skeleton_face_scores(
				p, face_score_cache, p.skeleton_face_score_normalize_by_area_, "[K5Delete]"))
		{
			std::cerr << "[K5Delete] failed to compute face scores." << std::endl;
			return;
		}

		const K5DetectionResult k5_info = detect_k5_cells(p);
		if (k5_info.cliques.empty())
		{
			std::cout << "[K5Delete] no K5 detected." << std::endl;
			return;
		}

		auto make_face_key = [](uint32 a, uint32 b, uint32 c) -> std::array<uint32, 3> {
			std::array<uint32, 3> key = {a, b, c};
			std::sort(key.begin(), key.end());
			return key;
		};

		std::map<std::array<uint32, 3>, NMFace> face_key_to_face;
		foreach_cell(*p.skeleton_, [&](NMFace f) -> bool {
			const std::vector<NMVertex> fv = incident_vertices(*p.skeleton_, f);
			if (fv.size() != 3)
				return true;
			const uint32 a = index_of(*p.skeleton_, fv[0]);
			const uint32 b = index_of(*p.skeleton_, fv[1]);
			const uint32 c = index_of(*p.skeleton_, fv[2]);
			if (a == INVALID_INDEX || b == INVALID_INDEX || c == INVALID_INDEX)
				return true;
			face_key_to_face[make_face_key(a, b, c)] = f;
			return true;
		});

		uint32 processed_k5 = 0;
		uint32 removed_faces = 0;
		uint32 removed_edges = 0;
		uint32 removed_tets = 0;

		for (const auto& clique : k5_info.cliques)
		{
			std::unordered_set<uint32> clique_face_ids;
			clique_face_ids.reserve(16);
			std::vector<NMFace> clique_faces;
			clique_faces.reserve(10);
			for (uint32 i = 0; i < 5; ++i)
			{
				for (uint32 j = i + 1; j < 5; ++j)
				{
					for (uint32 k = j + 1; k < 5; ++k)
					{
						const auto it_f = face_key_to_face.find(make_face_key(clique[i], clique[j], clique[k]));
						if (it_f == face_key_to_face.end())
							continue;
						const NMFace f = it_f->second;
						const uint32 idf = index_of(*p.skeleton_, f);
						if (idf == INVALID_INDEX)
							continue;
						clique_face_ids.insert(idf);
						clique_faces.push_back(f);
					}
				}
			}
			if (clique_faces.empty())
				continue;

			struct FaceCand
			{
				NMFace f;
				uint32 idf = INVALID_INDEX;
				Scalar score = Scalar(0);
				std::vector<NMEdge> deg2_edges;
			};

			bool has_best = false;
			FaceCand best;
			for (NMFace f : clique_faces)
			{
				if (!f.is_valid())
					continue;
				const uint32 idf = index_of(*p.skeleton_, f);
				if (idf == INVALID_INDEX)
					continue;
				std::vector<NMEdge> deg2_edges;
				deg2_edges.reserve(3);
				for (NMEdge e : incident_edges(*p.skeleton_, f))
				{
					if (!e.is_valid())
						continue;
					if (incident_faces(*p.skeleton_, e).size() == 2)
						deg2_edges.push_back(e);
				}
				if (deg2_edges.size() < 2)
					continue;
				const auto it_score = face_score_cache.find(idf);
				const Scalar score = (it_score != face_score_cache.end()) ? it_score->second : Scalar(0);
				if (!has_best || score > best.score)
				{
					has_best = true;
					best.f = f;
					best.idf = idf;
					best.score = score;
					best.deg2_edges = std::move(deg2_edges);
				}
			}
			if (!has_best)
				continue;

			struct NeighborFace
			{
				NMFace f;
				uint32 idf = INVALID_INDEX;
				Scalar score = Scalar(0);
			};
			std::unordered_map<uint32, NeighborFace> neighbor_faces;
			neighbor_faces.reserve(8);
			for (NMEdge e : best.deg2_edges)
			{
				if (!e.is_valid())
					continue;
				const auto in_faces = incident_faces(*p.skeleton_, e);
				if (in_faces.size() != 2)
					continue;
				NMFace nf = in_faces[0];
				if (nf == best.f)
					nf = in_faces[1];
				if (!nf.is_valid() || nf == best.f)
					continue;
				const uint32 idn = index_of(*p.skeleton_, nf);
				if (idn == INVALID_INDEX || idn == best.idf)
					continue;
				if (clique_face_ids.find(idn) == clique_face_ids.end())
					continue;
				const auto it_score = face_score_cache.find(idn);
				const Scalar score = (it_score != face_score_cache.end()) ? it_score->second : Scalar(0);
				auto itn = neighbor_faces.find(idn);
				if (itn == neighbor_faces.end() || score > itn->second.score)
					neighbor_faces[idn] = {nf, idn, score};
			}

			std::vector<NeighborFace> neighbors;
			neighbors.reserve(neighbor_faces.size());
			for (const auto& kvn : neighbor_faces)
				neighbors.push_back(kvn.second);
			std::sort(neighbors.begin(), neighbors.end(),
					  [](const NeighborFace& a, const NeighborFace& b) { return a.score > b.score; });
			if (neighbors.size() < 2)
				continue;

			std::vector<NMFace> faces_to_remove = {best.f, neighbors[0].f, neighbors[1].f};
			std::unordered_set<uint32> unique_face_ids;
			unique_face_ids.reserve(4);
			std::unordered_set<std::size_t> tets_to_erase;
			std::vector<NMEdge> affected_edges;
			affected_edges.reserve(9);
			for (NMFace f : faces_to_remove)
			{
				if (!f.is_valid())
					continue;
				const uint32 idf = index_of(*p.skeleton_, f);
				if (idf == INVALID_INDEX)
					continue;
				if (!unique_face_ids.insert(idf).second)
					continue;
				for (std::size_t tet_id : (*p.incident_tets_)[idf])
					if (p.skeleton_tets_.find(tet_id) != p.skeleton_tets_.end())
						tets_to_erase.insert(tet_id);
				for (NMEdge e : incident_edges(*p.skeleton_, f))
					affected_edges.push_back(e);
			}
			if (unique_face_ids.size() < 3)
				continue;

			uint32 removed_faces_this_k5 = 0;
			for (NMFace f : faces_to_remove)
			{
				if (!f.is_valid())
					continue;
				const uint32 idf = index_of(*p.skeleton_, f);
				if (idf == INVALID_INDEX)
					continue;
				remove_face(*p.skeleton_, f);
				++removed_faces_this_k5;
			}
			if (removed_faces_this_k5 == 0)
				continue;

			const uint32 removed_edges_this_k5 = remove_orphan_edges_from_removed_face_edges(p, affected_edges);
			uint32 removed_tets_this_k5 = 0;
			for (std::size_t tet_id : tets_to_erase)
			{
				auto it_tet = p.skeleton_tets_.find(tet_id);
				if (it_tet == p.skeleton_tets_.end())
					continue;
				const Tet old_tet = it_tet->second;
				for (uint32 i = 0; i < 4; ++i)
				{
					const NMFace tf = old_tet.faces[i];
					if (!tf.is_valid())
						continue;
					const uint32 idtf = index_of(*p.skeleton_, tf);
					if (idtf == INVALID_INDEX)
						continue;
					(*p.incident_tets_)[idtf].erase(tet_id);
				}
				p.skeleton_tets_.erase(it_tet);
				++removed_tets_this_k5;
			}

			++processed_k5;
			removed_faces += removed_faces_this_k5;
			removed_edges += removed_edges_this_k5;
			removed_tets += removed_tets_this_k5;

			std::cout << "[K5Delete] remove#" << processed_k5
					  << " seed_face=" << best.idf
					  << " seed_score=" << best.score
					  << " removed_faces=" << removed_faces_this_k5
					  << " removed_edges=" << removed_edges_this_k5
					  << " removed_tets=" << removed_tets_this_k5
					  << " remaining_tets=" << p.skeleton_tets_.size() << std::endl;
		}

		std::cout << "[K5Delete] done"
				  << " detected_k5=" << k5_info.cliques.size()
				  << " processed_k5=" << processed_k5
				  << " removed_faces=" << removed_faces
				  << " removed_edges=" << removed_edges
				  << " removed_tets=" << removed_tets
				  << " remaining_tets=" << p.skeleton_tets_.size() << std::endl;

		refresh_skeleton_topology_colors(p);
		mark_boundary_tets_color(p);
	}

	struct BoundaryTetPrepassStats
	{
		uint32 removed_faces = 0;
		uint32 removed_tets = 0;
		uint32 removed_edges = 0;
		uint32 processed_tets = 0;
		uint32 passes = 0;
	};

	BoundaryTetPrepassStats run_boundary_tet_face_deletion(
		PointsParameters& p, const char* log_prefix = "[TopologyFilter]")
	{
		BoundaryTetPrepassStats stats;

		struct BoundaryTetCandidate
		{
			std::size_t tet_id = std::size_t(-1);
			NMEdge edge;
			NMFace f0;
			NMFace f1;
		};

		std::vector<std::size_t> tet_ids;
		tet_ids.reserve(p.skeleton_tets_.size());
		for (const auto& kv : p.skeleton_tets_)
			tet_ids.push_back(kv.first);

		std::vector<BoundaryTetCandidate> candidates;
		candidates.reserve(tet_ids.size());

		// Phase 1: scan all current tets once and collect boundary-tet deletion candidates.
		for (std::size_t tet_id : tet_ids)
		{
			auto it_tet = p.skeleton_tets_.find(tet_id);
			if (it_tet == p.skeleton_tets_.end())
				continue;
			const Tet tet = it_tet->second;

			bool is_boundary_tet = false;
			std::unordered_set<uint32> tet_face_ids;
			tet_face_ids.reserve(4);
			for (uint32 i = 0; i < 4; ++i)
			{
				const NMFace f = tet.faces[i];
				const uint32 idf = index_of(*p.skeleton_, f);
				tet_face_ids.insert(idf);

				uint32 deg2_edge_count = 0;
				for (NMEdge e : incident_edges(*p.skeleton_, f))
				{
					if (incident_faces(*p.skeleton_, e).size() == 2)
						++deg2_edge_count;
				}
				if (deg2_edge_count >= 2)
				{
					is_boundary_tet = true;
					break;
				}
			}
			if (!is_boundary_tet || tet_face_ids.empty())
				continue;

			struct BoundaryEdgeChoice
			{
				NMEdge edge;
				NMFace f0;
				NMFace f1;
				Scalar length = Scalar(-1);
			};
			BoundaryEdgeChoice best_choice;
			std::unordered_set<uint32> seen_edges;
			seen_edges.reserve(8);

			for (uint32 i = 0; i < 4; ++i)
			{
				const NMFace f = tet.faces[i];
				for (NMEdge e : incident_edges(*p.skeleton_, f))
				{
					const uint32 ide = index_of(*p.skeleton_, e);
					if (!seen_edges.insert(ide).second)
						continue;

					const auto in_faces = incident_faces(*p.skeleton_, e);
					if (in_faces.size() != 2)
						continue;
					const NMFace f0 = in_faces[0];
					const NMFace f1 = in_faces[1];
					const uint32 idf0 = index_of(*p.skeleton_, f0);
					const uint32 idf1 = index_of(*p.skeleton_, f1);
					const std::vector<NMVertex> edge_vertices = incident_vertices(*p.skeleton_, e);
					if (edge_vertices.size() != 2)
						continue;
					const Vec3 p0 = value<Vec3>(*p.skeleton_, p.skeleton_position_, edge_vertices[0]);
					const Vec3 p1 = value<Vec3>(*p.skeleton_, p.skeleton_position_, edge_vertices[1]);
					const Scalar edge_len = (p1 - p0).norm();
					if (edge_len > best_choice.length)
					{
						best_choice.edge = e;
						best_choice.f0 = f0;
						best_choice.f1 = f1;
						best_choice.length = edge_len;
					}
				}
			}

			if (best_choice.length >= Scalar(0))
				candidates.push_back({tet_id, best_choice.edge, best_choice.f0, best_choice.f1});
		}

		// Phase 2: delete using precomputed candidates only (no boundary re-check during deletion).
		stats.passes = 1;
		for (const BoundaryTetCandidate& cand : candidates)
		{
			auto it_remove_tet = p.skeleton_tets_.find(cand.tet_id);
			if (it_remove_tet == p.skeleton_tets_.end())
				continue;
			const Tet old_tet = it_remove_tet->second;

			const std::array<NMFace, 2> faces_to_remove = {cand.f0, cand.f1};
			for (const NMFace& f_del : faces_to_remove)
			{
				const uint32 idf_del = index_of(*p.skeleton_, f_del);
				remove_face(*p.skeleton_, f_del);
				++stats.removed_faces;
			}
			if (cand.edge.is_valid())
			{
				const auto edge_in_faces = incident_faces(*p.skeleton_, cand.edge);
				if (edge_in_faces.empty())
				{
					remove_edge(*p.skeleton_, cand.edge);
					++stats.removed_edges;
				}
			}
			p.skeleton_tets_.erase(it_remove_tet);
			++stats.removed_tets;
			for (uint32 i = 0; i < 4; ++i)
			{
				const NMFace f = old_tet.faces[i];
				const uint32 idf = index_of(*p.skeleton_, f);
				(*p.incident_tets_)[idf].erase(cand.tet_id);
			}
			++stats.processed_tets;
		}

		std::cout << log_prefix << " boundary_prepass_passes=" << stats.passes
				  << " processed_boundary_tets=" << stats.processed_tets
				  << " removed_faces=" << stats.removed_faces
				  << " removed_edges=" << stats.removed_edges
				  << " removed_tets=" << stats.removed_tets << std::endl;
		return stats;
	}

	void refresh_skeleton_topology_colors(PointsParameters& p)
	{
		compute_edge_degree(p);
		foreach_cell(*p.skeleton_, [&](NMEdge e) {
			auto in_face = incident_faces(*p.skeleton_, e);
			if (in_face.size() == 1)
				value<Vec3>(*p.skeleton_, p.skeleton_edge_color_, e) = Vec3(0.0, 0.0, 1.0);
			return true;
		});
		parallel_foreach_cell(*p.skeleton_, [&](NMFace f) -> bool {
			auto in_tets = (*p.incident_tets_)[index_of(*p.skeleton_, f)];
			if (in_tets.size() == 1)
				value<Vec3>(*p.skeleton_, p.skeleton_face_color_, f) = Vec3(1.0, 0.0, 0.0);
			else if (in_tets.size() > 1)
				value<Vec3>(*p.skeleton_, p.skeleton_face_color_, f) = Vec3(0.8, 0.5, 0.5);
			else
				value<Vec3>(*p.skeleton_, p.skeleton_face_color_, f) = Vec3(0.0, 0.0, 0.0);
			return true;
		});
		non_manifold_provider_->emit_connectivity_changed(*p.skeleton_);
		non_manifold_provider_->emit_attribute_changed(*p.skeleton_, p.skeleton_edge_color_.get());
		non_manifold_provider_->emit_attribute_changed(*p.skeleton_, p.skeleton_edge_non_manifold_color_.get());
		non_manifold_provider_->emit_attribute_changed(*p.skeleton_, p.skeleton_face_color_.get());
	}

	void refresh_skeleton_attribute_handles(PointsParameters& p)
	{
		p.skeleton_position_ = get_or_add_attribute<Vec3, NMVertex>(*p.skeleton_, "position");
		p.skeleton_radius_ = get_or_add_attribute<Scalar, NMVertex>(*p.skeleton_, "radius");
		p.incident_tets_ = get_or_add_attribute<std::set<std::size_t>, NMFace>(*p.skeleton_, "incident_tets");
		p.edge_degree_ = get_or_add_attribute<uint32, NMEdge>(*p.skeleton_, "degree");
		p.skeleton_face_color_ = get_or_add_attribute<Vec3, NMFace>(*p.skeleton_, "color");
		p.skeleton_face_udf_color_ = get_or_add_attribute<Vec3, NMFace>(*p.skeleton_, "udf_color");
		p.skeleton_face_boundary_tet_color_ = get_or_add_attribute<Vec3, NMFace>(*p.skeleton_, "boundary_tet_color");
		p.skeleton_edge_color_ = get_or_add_attribute<Vec3, NMEdge>(*p.skeleton_, "color");
		p.skeleton_edge_udf_score_color_ = get_or_add_attribute<Vec3, NMEdge>(*p.skeleton_, "edge_udf_score_color");
		p.skeleton_edge_boundary_tet_color_ =
			get_or_add_attribute<Vec3, NMEdge>(*p.skeleton_, "boundary_tet_best_edge_color");
		p.skeleton_edge_non_manifold_color_ =
			get_or_add_attribute<Vec3, NMEdge>(*p.skeleton_, "non_manifold_edge_color");
	}

	void invalidate_topology_stage_snapshot(PointsParameters& p)
	{
		p.topology_stage_snapshot_valid_ = false;
		p.topology_stage_snapshot_faces_map_.clear();
		p.topology_stage_snapshot_tets_.clear();
	}

	bool capture_topology_stage_snapshot(PointsParameters& p)
	{
		if (!p.skeleton_)
		{
			std::cerr << "[TopologyStage] cannot capture snapshot: skeleton is null." << std::endl;
			return false;
		}

		if (!p.topology_stage_snapshot_mesh_)
		{
			const std::string snapshot_name = non_manifold_provider_->mesh_name(*p.skeleton_) + "_topology_stage_snapshot";
			p.topology_stage_snapshot_mesh_ = non_manifold_provider_->has_mesh(snapshot_name)
												  ? non_manifold_provider_->mesh(snapshot_name)
												  : non_manifold_provider_->add_mesh(snapshot_name);
		}
		non_manifold_provider_->copy_mesh(*p.topology_stage_snapshot_mesh_, *p.skeleton_);
		p.topology_stage_snapshot_faces_map_ = p.skeleton_faces_map_;
		p.topology_stage_snapshot_tets_ = p.skeleton_tets_;
		p.topology_stage_snapshot_valid_ = true;
		std::cout << "[TopologyStage] snapshot captured: faces=" << nb_cells<NMFace>(*p.skeleton_)
				  << " tets=" << p.skeleton_tets_.size() << std::endl;
		return true;
	}

	bool restore_topology_stage_snapshot(PointsParameters& p)
	{
		if (!p.topology_stage_snapshot_valid_ || !p.topology_stage_snapshot_mesh_)
		{
			std::cerr << "[TopologyStage] snapshot not available." << std::endl;
			return false;
		}
		non_manifold_provider_->copy_mesh(*p.skeleton_, *p.topology_stage_snapshot_mesh_);
		refresh_skeleton_attribute_handles(p);
		p.skeleton_faces_map_ = p.topology_stage_snapshot_faces_map_;
		p.skeleton_tets_ = p.topology_stage_snapshot_tets_;
		return true;
	}

	void boundary_tet_face_delete_only(PointsParameters& p)
	{
		if (!p.skeleton_ || !p.incident_tets_)
		{
			std::cerr << "Boundary tet delete requires a built skeleton." << std::endl;
			return;
		}
		run_boundary_tet_face_deletion(p, "[BoundaryTetDelete]");
		refresh_skeleton_topology_colors(p);
		mark_boundary_tets_color(p);
	}

	std::unordered_set<uint32> collect_current_tet_face_id_whitelist(PointsParameters& p)
	{
		std::unordered_set<uint32> tet_face_ids;
		tet_face_ids.reserve(nb_cells<NMFace>(*p.skeleton_));
		for (const auto& kv : p.skeleton_tets_)
		{
			const Tet& tet = kv.second;
			for (uint32 i = 0; i < 4; ++i)
			{
				const NMFace f = tet.faces[i];
				if (!f.is_valid())
					continue;
				const uint32 idf = index_of(*p.skeleton_, f);
				if (idf != INVALID_INDEX)
					tet_face_ids.insert(idf);
			}
		}
		return tet_face_ids;
	}

	void prune_deg_faces_from_whitelist_and_orphan_edges(
		PointsParameters& p, const std::unordered_set<uint32>& tet_face_whitelist, const char* log_prefix = "[DegFacePrune]")
	{
		if (!p.skeleton_ || !p.incident_tets_)
		{
			std::cerr << log_prefix << " requires a built skeleton." << std::endl;
			return;
		}

		auto is_deg_face = [&](const NMFace& f) -> bool {
			if (!f.is_valid())
				return false;
			const uint32 idf = index_of(*p.skeleton_, f);
			if (idf == INVALID_INDEX)
				return false;
			if (tet_face_whitelist.find(idf) == tet_face_whitelist.end())
				return false;

			uint32 deg1_count = 0;
			uint32 gt2_count = 0;
			for (NMEdge e : incident_edges(*p.skeleton_, f))
			{
				if (!e.is_valid())
					continue;
				const size_t edge_degree = incident_faces(*p.skeleton_, e).size();
				if (edge_degree == 1)
					++deg1_count;
				else if (edge_degree > 2)
					++gt2_count;
			}
			return (deg1_count == 1 && gt2_count == 2);
		};

		uint32 total_removed_faces = 0;
		uint32 total_removed_edges = 0;
		uint32 total_removed_tets = 0;
		uint32 total_passes = 0;

		while (true)
		{
			std::vector<NMFace> faces_to_remove;
			std::vector<NMEdge> affected_edges;
			std::unordered_set<std::size_t> tets_to_erase;

			foreach_cell(*p.skeleton_, [&](NMFace f) -> bool {
				if (!is_deg_face(f))
					return true;
				const uint32 idf = index_of(*p.skeleton_, f);
				for (std::size_t tet_id : (*p.incident_tets_)[idf])
					if (p.skeleton_tets_.find(tet_id) != p.skeleton_tets_.end())
						tets_to_erase.insert(tet_id);
				faces_to_remove.push_back(f);
				for (NMEdge e : incident_edges(*p.skeleton_, f))
					affected_edges.push_back(e);
				return true;
			});

			if (faces_to_remove.empty())
				break;

			++total_passes;
			uint32 removed_faces_this_pass = 0;
			for (const NMFace& f : faces_to_remove)
			{
				if (!f.is_valid())
					continue;
				const uint32 idf = index_of(*p.skeleton_, f);
				if (idf == INVALID_INDEX)
					continue;
				remove_face(*p.skeleton_, f);
				++removed_faces_this_pass;
			}

			uint32 removed_edges_this_pass = remove_orphan_edges_from_removed_face_edges(p, affected_edges);
			uint32 removed_tets_this_pass = 0;
			for (std::size_t tet_id : tets_to_erase)
			{
				auto it_tet = p.skeleton_tets_.find(tet_id);
				if (it_tet == p.skeleton_tets_.end())
					continue;
				const Tet old_tet = it_tet->second;
				for (uint32 i = 0; i < 4; ++i)
				{
					const NMFace f = old_tet.faces[i];
					if (!f.is_valid())
						continue;
					const uint32 idf = index_of(*p.skeleton_, f);
					if (idf == INVALID_INDEX)
						continue;
					(*p.incident_tets_)[idf].erase(tet_id);
				}
				p.skeleton_tets_.erase(it_tet);
				++removed_tets_this_pass;
			}

			total_removed_faces += removed_faces_this_pass;
			total_removed_edges += removed_edges_this_pass;
			total_removed_tets += removed_tets_this_pass;

			std::cout << log_prefix << " pass=" << total_passes
					  << " removed_deg_faces=" << removed_faces_this_pass
					  << " removed_orphan_edges=" << removed_edges_this_pass
					  << " removed_tets=" << removed_tets_this_pass
					  << " remaining_tets=" << p.skeleton_tets_.size() << std::endl;
		}

		std::cout << log_prefix << " done"
				  << " passes=" << total_passes
				  << " removed_deg_faces=" << total_removed_faces
				  << " removed_orphan_edges=" << total_removed_edges
				  << " removed_tets=" << total_removed_tets
				  << " remaining_tets=" << p.skeleton_tets_.size() << std::endl;
	}

	void run_complete_topology_fix_pipeline(PointsParameters& p)
	{
		if (!p.skeleton_ || !p.incident_tets_)
		{
			std::cerr << "[TopologyFull] requires a built skeleton." << std::endl;
			return;
		}

		const std::unordered_set<uint32> tet_face_whitelist = collect_current_tet_face_id_whitelist(p);
		std::cout << "[TopologyFull] start"
				  << " whitelist_tet_faces=" << tet_face_whitelist.size()
				  << " initial_tets=" << p.skeleton_tets_.size() << std::endl;

		// 1) Boundary tet mask -> boundary tet delete
		mark_boundary_tets_color(p);
		run_boundary_tet_face_deletion(p, "[TopologyFull]");
		refresh_skeleton_topology_colors(p);

		// 2) Edge score simple tet delete
		run_edge_score_simple_tet_topology_fix(p);
		// 3) Edge score non-simple tet delete
		run_edge_score_nonsimple_tet_topology_fix(p);

		// 4) Boundary tet mask -> boundary tet delete
		mark_boundary_tets_color(p);
		run_boundary_tet_face_deletion(p, "[TopologyFull]");
		refresh_skeleton_topology_colors(p);

		// 5) Edge score simple tet delete
		run_edge_score_simple_tet_topology_fix(p);

		// 6) Deg-face prune + orphan-edge cleanup (faces that used to be tet faces only)
		prune_deg_faces_from_whitelist_and_orphan_edges(p, tet_face_whitelist, "[TopologyFullDeg]");

		refresh_skeleton_topology_colors(p);
		mark_boundary_tets_color(p);
		visualize_skeleton_edge_udf_scores(p);
		std::cout << "[TopologyFull] done remaining_tets=" << p.skeleton_tets_.size() << std::endl;
	}

	void run_topology_stage_filter_from_snapshot(PointsParameters& p, bool run_edge_stage)
	{
		if (!p.skeleton_ || !p.incident_tets_)
		{
			std::cerr << "Topology stage filter requires a built skeleton." << std::endl;
			return;
		}
		if (!p.topology_stage_snapshot_valid_)
		{
			if (!capture_topology_stage_snapshot(p))
				return;
		}
		if (!restore_topology_stage_snapshot(p))
			return;

		std::cout << "[TopologyStage] run=" << (run_edge_stage ? "edge(face+edge)" : "face-only")
				  << " from_snapshot=true"
				  << " edge_diffuse=" << (p.topology_edge_stage_diffuse_ ? "on" : "off") << std::endl;
		skeleton_post_pocessing(p, run_edge_stage, p.topology_edge_stage_diffuse_);
	}

	struct FaceStageCandidate
	{
		NMFace face;
		uint32 face_id = INVALID_INDEX;
		Scalar score = Scalar(0);
	};

	bool pick_next_face_stage_candidate(PointsParameters& p, const std::unordered_map<uint32, Scalar>& face_score_cache,
										FaceStageCandidate& out_candidate)
	{
		if (!p.skeleton_ || !p.incident_tets_)
			return false;

		// Same tet-face cap model as main topology flow.
		std::unordered_map<uint32, std::vector<std::size_t>> face_owner_tets;
		face_owner_tets.reserve(nb_cells<NMFace>(*p.skeleton_));
		std::unordered_map<std::size_t, uint32> deleted_face_count_per_tet;
		deleted_face_count_per_tet.reserve(p.skeleton_tets_.size());
		for (const auto& kv : p.skeleton_tets_)
		{
			const std::size_t tet_id = kv.first;
			const Tet& tet = kv.second;
			uint32 initial_deleted_faces = 0;
			for (uint32 i = 0; i < 4; ++i)
			{
				const NMFace f = tet.faces[i];
				if (!f.is_valid())
				{
					++initial_deleted_faces;
					continue;
				}
				const uint32 idf = index_of(*p.skeleton_, f);
				if (idf == INVALID_INDEX)
				{
					++initial_deleted_faces;
					continue;
				}
				face_owner_tets[idf].push_back(tet_id);
			}
			deleted_face_count_per_tet[tet_id] = initial_deleted_faces;
		}
		auto is_blocked_by_tet_face_cap = [&](uint32 face_id) -> bool {
			auto owner_it = face_owner_tets.find(face_id);
			if (owner_it == face_owner_tets.end())
				return false;
			for (std::size_t owner_tet : owner_it->second)
			{
				if (p.skeleton_tets_.find(owner_tet) == p.skeleton_tets_.end())
					continue;
				auto it_count = deleted_face_count_per_tet.find(owner_tet);
				if (it_count != deleted_face_count_per_tet.end() && it_count->second >= 2)
					return true;
			}
			return false;
		};

		struct FaceQueueItem
		{
			Scalar score = Scalar(0);
			NMFace face;
			uint32 face_id = INVALID_INDEX;
		};
		auto cmp = [](const FaceQueueItem& a, const FaceQueueItem& b) { return a.score < b.score; };
		std::priority_queue<FaceQueueItem, std::vector<FaceQueueItem>, decltype(cmp)> simple_queue(cmp);
		std::unordered_set<uint32> simple_face_ids;
		simple_face_ids.reserve(nb_cells<NMFace>(*p.skeleton_));

		auto push_face_by_type = [&](NMFace f) {
			if (!f.is_valid())
				return;
			const uint32 idf = index_of(*p.skeleton_, f);
			if (idf == INVALID_INDEX)
				return;
			const size_t tet_count = (*p.incident_tets_)[idf].size();
			if (tet_count == 0)
			{
				simple_face_ids.erase(idf);
				return;
			}
			auto it = face_score_cache.find(idf);
			const Scalar score = (it != face_score_cache.end()) ? it->second : Scalar(0);
			if (!is_simple_face_3d(p, f))
			{
				simple_face_ids.erase(idf);
				return;
			}
			if (simple_face_ids.insert(idf).second)
				simple_queue.push({score, f, idf});
		};

		foreach_cell(*p.skeleton_, [&](NMFace f) -> bool {
			push_face_by_type(f);
			return true;
		});

		auto pop_next_valid = [&](FaceQueueItem& out) -> bool {
			while (!simple_queue.empty())
			{
				FaceQueueItem item = simple_queue.top();
				simple_queue.pop();

				const auto itid = simple_face_ids.find(item.face_id);
				if (itid == simple_face_ids.end())
					continue;
				simple_face_ids.erase(itid);

				if (!item.face.is_valid())
					continue;
				const uint32 idf = index_of(*p.skeleton_, item.face);
				if (idf == INVALID_INDEX || idf != item.face_id)
					continue;
				if (!is_simple_face_3d(p, item.face))
					continue;
				out = item;
				return true;
			}
			return false;
		};

		while (true)
		{
			FaceQueueItem item;
			if (!pop_next_valid(item))
			{
				if (simple_face_ids.empty())
					return false;
				continue;
			}

			const uint32 idf = item.face_id;
			if (is_blocked_by_tet_face_cap(idf))
				continue;
			if (violates_face_removal_edge_guard(p, item.face))
				continue;

			std::set<std::size_t> in_tets = (*p.incident_tets_)[idf];
			if (in_tets.empty() || !is_simple_face_3d(p, item.face))
				continue;

			out_candidate.face = item.face;
			out_candidate.face_id = idf;
			out_candidate.score = item.score;
			return true;
		}
	}

	void run_face_stage_single_step_and_highlight_next(PointsParameters& p)
	{
		if (!p.skeleton_ || !p.incident_tets_)
		{
			std::cerr << "Face-stage single step requires a built skeleton." << std::endl;
			return;
		}
		if (!p.neural_udf_loaded_)
		{
			std::cerr << "Face-stage single step requires a loaded neural UDF model." << std::endl;
			return;
		}

		std::unordered_map<uint32, Scalar> face_score_cache;
		if (!compute_skeleton_face_scores(
				p, face_score_cache, p.skeleton_face_score_normalize_by_area_, "[TopologyStep]"))
		{
			std::cerr << "[TopologyStep] failed to compute face scores." << std::endl;
			return;
		}

		FaceStageCandidate current;
		if (!pick_next_face_stage_candidate(p, face_score_cache, current))
		{
			refresh_skeleton_topology_colors(p);
			mark_boundary_tets_color(p);
			std::cout << "[TopologyStep] no removable face candidate." << std::endl;
			return;
		}

		const std::set<std::size_t> in_tets = (*p.incident_tets_)[current.face_id];
		const std::vector<NMEdge> in_edges = incident_edges(*p.skeleton_, current.face);
		remove_face(*p.skeleton_, current.face);

		uint32 removed_tets_this_step = 0;
		for (std::size_t tet_id : in_tets)
		{
			auto it_tet = p.skeleton_tets_.find(tet_id);
			if (it_tet != p.skeleton_tets_.end())
			{
				const Tet old_tet = it_tet->second;
				p.skeleton_tets_.erase(it_tet);
				++removed_tets_this_step;
				for (uint32 i = 0; i < 4; ++i)
				{
					const NMFace f = old_tet.faces[i];
					if (!f.is_valid() || f == current.face)
						continue;
					const uint32 idf = index_of(*p.skeleton_, f);
					if (idf == INVALID_INDEX)
						continue;
					(*p.incident_tets_)[idf].erase(tet_id);
				}
			}
			else
			{
				foreach_cell(*p.skeleton_, [&](NMFace f) -> bool {
					const uint32 idf = index_of(*p.skeleton_, f);
					if (idf != INVALID_INDEX)
						(*p.incident_tets_)[idf].erase(tet_id);
					return true;
				});
			}
		}

		uint32 removed_deg1_edges = 0;
		for (NMEdge e : in_edges)
		{
			if (!e.is_valid())
				continue;
			const auto in_faces_after = incident_faces(*p.skeleton_, e);
			if (in_faces_after.size() == 1)
			{
				remove_edge(*p.skeleton_, e);
				++removed_deg1_edges;
			}
		}

		refresh_skeleton_topology_colors(p);
		mark_boundary_tets_color(p);

		FaceStageCandidate next_candidate;
		bool has_next = false;
		has_next = pick_next_face_stage_candidate(p, face_score_cache, next_candidate);

		if (has_next && next_candidate.face.is_valid() && p.skeleton_face_color_)
		{
			value<Vec3>(*p.skeleton_, p.skeleton_face_color_, next_candidate.face) = Vec3(1.0, 1.0, 0.0);
			non_manifold_provider_->emit_attribute_changed(*p.skeleton_, p.skeleton_face_color_.get());
		}

		std::cout << "[TopologyStep] removed_face=" << current.face_id
				  << " type=simple"
				  << " score=" << current.score
				  << " removed_tets=" << removed_tets_this_step
				  << " removed_deg1_edges=" << removed_deg1_edges << std::endl;
		if (has_next)
			std::cout << "[TopologyStep] next_face=" << next_candidate.face_id
					  << " next_type=simple"
					  << " next_score=" << next_candidate.score << " highlighted_color=(1,1,0)" << std::endl;
		else
			std::cout << "[TopologyStep] next_face=none" << std::endl;
	}

	void skeleton_post_pocessing(PointsParameters& p, bool run_edge_stage = true, bool edge_stage_diffuse = true)
	{
		(void)edge_stage_diffuse;
		if (!p.neural_udf_loaded_)
		{
			std::cerr << "Topology filter requires a loaded neural UDF model." << std::endl;
			return;
		}
		std::unordered_map<uint32, Scalar> face_score_cache;
		if (!compute_skeleton_face_scores(
				p, face_score_cache, p.skeleton_face_score_normalize_by_area_, "[TopologyFilter]"))
		{
			std::cerr << "[TopologyFilter] Face scores may be incomplete due to UDF evaluation failure." << std::endl;
		}

		const BoundaryTetPrepassStats boundary_prepass = run_boundary_tet_face_deletion(p, "[TopologyFilter]");
		std::unordered_set<uint32> initial_tet_face_ids;
		initial_tet_face_ids.reserve(nb_cells<NMFace>(*p.skeleton_));
		for (const auto& kv : p.skeleton_tets_)
		{
			const Tet& tet = kv.second;
			for (uint32 i = 0; i < 4; ++i)
			{
				const NMFace f = tet.faces[i];
				const uint32 idf = index_of(*p.skeleton_, f);
				if (idf != INVALID_INDEX)
					initial_tet_face_ids.insert(idf);
			}
		}

		struct FaceQueueItem
		{
			Scalar score = Scalar(0);
			NMFace face;
			uint32 face_id = INVALID_INDEX;
		};
		auto cmp = [](const FaceQueueItem& a, const FaceQueueItem& b) { return a.score < b.score; };
		std::priority_queue<FaceQueueItem, std::vector<FaceQueueItem>, decltype(cmp)> simple_queue(cmp);
		std::unordered_set<uint32> simple_face_ids;
		simple_face_ids.reserve(nb_cells<NMFace>(*p.skeleton_));
		size_t total_pushes_simple = 0;

		auto push_face_by_type = [&](NMFace f) {
			const uint32 idf = index_of(*p.skeleton_, f);
			const size_t tet_count = (*p.incident_tets_)[idf].size();
			if (tet_count == 0)
			{
				simple_face_ids.erase(idf);
				return;
			}
			auto it = face_score_cache.find(idf);
			const Scalar score = (it != face_score_cache.end()) ? it->second : Scalar(0);
			if (it == face_score_cache.end())
				face_score_cache[idf] = Scalar(0);

			if (!is_simple_face_3d(p, f))
			{
				simple_face_ids.erase(idf);
				return;
			}
			if (simple_face_ids.insert(idf).second)
			{
				simple_queue.push({score, f, idf});
				++total_pushes_simple;
			}
		};

		size_t initial_simple_faces = 0;
		foreach_cell(*p.skeleton_, [&](NMFace f) -> bool {
			if (is_simple_face_3d(p, f))
				++initial_simple_faces;
			push_face_by_type(f);
			return true;
		});
		std::cout << "[TopologyFilter] initial_faces=" << nb_cells<NMFace>(*p.skeleton_)
				  << " initial_tets=" << p.skeleton_tets_.size()
				  << " initial_simple_faces=" << initial_simple_faces
				  << " initial_simple_queue_size=" << simple_queue.size() << std::endl;
		std::cout << "[TopologyFilter] initial_tet_faces=" << initial_tet_face_ids.size() << std::endl;

		// Global budget (face+edge stage): each tet can lose at most two faces.
		std::unordered_map<uint32, std::vector<std::size_t>> face_owner_tets;
		face_owner_tets.reserve(nb_cells<NMFace>(*p.skeleton_));
		std::unordered_map<std::size_t, uint32> deleted_face_count_per_tet;
		deleted_face_count_per_tet.reserve(p.skeleton_tets_.size());
		for (const auto& kv : p.skeleton_tets_)
		{
			const std::size_t tet_id = kv.first;
			const Tet& tet = kv.second;
			uint32 initial_deleted_faces = 0;
			for (uint32 i = 0; i < 4; ++i)
			{
				const NMFace f = tet.faces[i];
				if (!f.is_valid())
				{
					++initial_deleted_faces;
					continue;
				}
				const uint32 idf = index_of(*p.skeleton_, f);
				if (idf == INVALID_INDEX)
				{
					++initial_deleted_faces;
					continue;
				}
				face_owner_tets[idf].push_back(tet_id);
			}
			// "From the beginning" budget: include faces already missing after boundary prepass.
			deleted_face_count_per_tet[tet_id] = initial_deleted_faces;
		}

		uint32 removed_face = boundary_prepass.removed_faces;
		uint32 removed_simple_face = 0;
		uint32 removed_tets = boundary_prepass.removed_tets;
		uint32 removed_edges = boundary_prepass.removed_edges;
		size_t pop_count = 0;
		size_t skipped_invalid = 0;
		size_t skipped_not_simple = 0;
		size_t skipped_not_single_tet = 0;
		size_t skipped_tet_face_cap = 0;
		size_t skipped_edge_guard = 0;
		size_t total_edge_followup_pushes = 0;
		size_t total_degree1_edges_removed_after_face = 0;
		auto is_blocked_by_tet_face_cap = [&](uint32 face_id) -> bool {
			auto owner_it = face_owner_tets.find(face_id);
			if (owner_it == face_owner_tets.end())
				return false;
			for (std::size_t owner_tet : owner_it->second)
			{
				if (p.skeleton_tets_.find(owner_tet) == p.skeleton_tets_.end())
					continue;
				auto it_count = deleted_face_count_per_tet.find(owner_tet);
				if (it_count != deleted_face_count_per_tet.end() && it_count->second >= 2)
					return true;
			}
			return false;
		};
		auto account_face_deletion_to_owner_tets = [&](uint32 face_id) {
			auto owner_it = face_owner_tets.find(face_id);
			if (owner_it == face_owner_tets.end())
				return;
			for (std::size_t owner_tet : owner_it->second)
			{
				if (p.skeleton_tets_.find(owner_tet) == p.skeleton_tets_.end())
					continue;
				auto it_count = deleted_face_count_per_tet.find(owner_tet);
				if (it_count != deleted_face_count_per_tet.end())
					++(it_count->second);
			}
		};
		auto remove_incident_tets_from_face = [&](const std::set<std::size_t>& in_tets, const NMFace& removed_face,
												  bool update_face_queues) -> uint32 {
			uint32 removed_tets_local = 0;
			for (std::size_t tet_id : in_tets)
			{
				auto it_tet = p.skeleton_tets_.find(tet_id);
				if (it_tet != p.skeleton_tets_.end())
				{
					const Tet old_tet = it_tet->second;
					p.skeleton_tets_.erase(it_tet);
					++removed_tets_local;
					for (uint32 i = 0; i < 4; ++i)
					{
						const NMFace f = old_tet.faces[i];
						if (!f.is_valid() || f == removed_face)
							continue;
						const uint32 idf = index_of(*p.skeleton_, f);
						if (idf == INVALID_INDEX)
							continue;
						(*p.incident_tets_)[idf].erase(tet_id);
						if (update_face_queues)
							push_face_by_type(f);
					}
				}
				else
				{
					foreach_cell(*p.skeleton_, [&](NMFace f) -> bool {
						const uint32 idf = index_of(*p.skeleton_, f);
						if (idf != INVALID_INDEX)
						{
							(*p.incident_tets_)[idf].erase(tet_id);
							if (update_face_queues)
								push_face_by_type(f);
						}
						return true;
					});
				}
			}
			return removed_tets_local;
		};
		auto pop_next_valid = [&](FaceQueueItem& out) -> bool {
			while (!simple_queue.empty())
			{
				FaceQueueItem item = simple_queue.top();
				simple_queue.pop();
				++pop_count;

				const auto itid = simple_face_ids.find(item.face_id);
				if (itid == simple_face_ids.end())
					continue;
				simple_face_ids.erase(itid);

				if (!item.face.is_valid())
				{
					++skipped_invalid;
					continue;
				}
				const uint32 idf = index_of(*p.skeleton_, item.face);
				if (idf == INVALID_INDEX || idf != item.face_id)
				{
					++skipped_invalid;
					continue;
				}
				if (!is_simple_face_3d(p, item.face))
				{
					++skipped_not_simple;
					continue;
				}
				out = item;
				return true;
			}
			return false;
		};

		while (true)
		{
			FaceQueueItem current_item;
			if (!pop_next_valid(current_item))
			{
				if (simple_face_ids.empty())
					break;
				continue;
			}

			const NMFace current_face = current_item.face;
			const uint32 idf_current = index_of(*p.skeleton_, current_face);
			if (is_blocked_by_tet_face_cap(idf_current))
			{
				++skipped_tet_face_cap;
				continue;
			}
			if (violates_face_removal_edge_guard(p, current_face))
			{
				++skipped_edge_guard;
				continue;
			}
			std::set<std::size_t> in_tets = (*p.incident_tets_)[idf_current];
			if (!is_simple_face_3d(p, current_face))
			{
				++skipped_not_single_tet;
				continue;
			}
			if (in_tets.empty())
			{
				++skipped_not_single_tet;
				continue;
			}
			const std::vector<NMEdge> in_edges = incident_edges(*p.skeleton_, current_face);
			remove_face(*p.skeleton_, current_face);
			++removed_face;
			++removed_simple_face;
			simple_face_ids.erase(idf_current);
			account_face_deletion_to_owner_tets(idf_current);

			size_t newly_pushed = 0;
			size_t removed_deg1_edges_this_step = 0;
			for (NMEdge e : in_edges)
			{
				if (!e.is_valid())
					continue;

				const auto in_faces_after = incident_faces(*p.skeleton_, e);
				if (in_faces_after.size() == 1)
				{
					const NMFace next_face = in_faces_after[0];
					if (next_face.is_valid())
					{
						const uint32 idf_next = index_of(*p.skeleton_, next_face);
						if (idf_next != INVALID_INDEX)
						{
							if (initial_tet_face_ids.find(idf_next) != initial_tet_face_ids.end())
							{
								const size_t before_simple = simple_face_ids.size();
								push_face_by_type(next_face);
								if (simple_face_ids.size() > before_simple)
									++newly_pushed;
							}
						}
					}
					remove_edge(*p.skeleton_, e);
					++removed_edges;
					++removed_deg1_edges_this_step;
				}
			}
			uint32 removed_tets_this_step = remove_incident_tets_from_face(in_tets, current_face, true);
			removed_tets += removed_tets_this_step;
			total_edge_followup_pushes += newly_pushed;
			total_degree1_edges_removed_after_face += removed_deg1_edges_this_step;
			std::cout << "[TopologyFilter] remove#" << removed_face << " type=simple"
					  << " face=" << idf_current
					  << " score=" << current_item.score << " incident_tets_removed=" << removed_tets_this_step
					  << " newly_pushed=" << newly_pushed
					  << " removed_deg1_edges=" << removed_deg1_edges_this_step
					  << " simple_queue_size=" << simple_queue.size()
					  << " remaining_tets=" << p.skeleton_tets_.size() << std::endl;
		}

		if (run_edge_stage)
			std::cout << "[TopologyFilter] edge_stage_merged=true (face deletion now handles edge-driven propagation)"
					  << std::endl;
		else
			std::cout << "[TopologyFilter] unified_face_flow_only=true" << std::endl;

		std::cout << "Removed " << removed_face << " faces in skeleton post-processing." << std::endl;
		std::cout << "[TopologyFilter] pop_count=" << pop_count << " skipped_invalid=" << skipped_invalid
				  << " skipped_not_simple=" << skipped_not_simple
				  << " skipped_not_single_tet=" << skipped_not_single_tet
				  << " skipped_tet_face_cap=" << skipped_tet_face_cap
				  << " skipped_edge_guard=" << skipped_edge_guard
				  << " boundary_prepass_processed_tets=" << boundary_prepass.processed_tets
				  << " boundary_prepass_removed_faces=" << boundary_prepass.removed_faces
				  << " boundary_prepass_removed_edges=" << boundary_prepass.removed_edges
				  << " boundary_prepass_removed_tets=" << boundary_prepass.removed_tets
				  << " total_edge_followup_pushes=" << total_edge_followup_pushes
				  << " total_removed_deg1_edges_after_face=" << total_degree1_edges_removed_after_face
				  << " removed_edges=" << removed_edges
				  << " removed_simple_faces=" << removed_simple_face
				  << " removed_tets=" << removed_tets
				  << " total_pushes_simple=" << total_pushes_simple << std::endl;
		std::cout << p.skeleton_tets_.size() << " tets remain after post-processing." << std::endl;

		refresh_skeleton_topology_colors(p);
		mark_boundary_tets_color(p);
	}

protected:
	void update_render_data(PointsParameters& p, bool non_blocking_running_lock = false, bool full_refresh = true)
	{
		std::unique_lock<std::mutex> lock(p.mutex_, std::defer_lock);
		if (p.running_)
		{
			if (non_blocking_running_lock)
			{
				if (!lock.try_lock())
					return;
			}
			else
			{
				lock.lock();
			}
		}

		points_provider_->emit_connectivity_changed(*p.spheres_);
		points_provider_->emit_attribute_changed(*p.spheres_, p.spheres_position_.get());
		points_provider_->emit_attribute_changed(*p.spheres_, p.spheres_radius_.get());

		update_spheres_color(p);
		update_spheres_topup_zero_inside_color(p);
		points_provider_->emit_attribute_changed(*p.spheres_, p.spheres_color_.get());
		if (p.spheres_topup_zero_inside_color_)
			points_provider_->emit_attribute_changed(*p.spheres_, p.spheres_topup_zero_inside_color_.get());

		if (!full_refresh)
			return;

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

		update_opt_mesh_colors(p);
		if (p.opt_mesh_)
		{
			points_provider_->emit_attribute_changed(*p.opt_mesh_, p.opt_color_.get());
			points_provider_->emit_attribute_changed(*p.opt_mesh_, p.opt_normal_.get());
			if (p.opt_normal_color_)
				points_provider_->emit_attribute_changed(*p.opt_mesh_, p.opt_normal_color_.get());
		}

		compute_skeleton_current_mode(p);
		non_manifold_provider_->emit_connectivity_changed(*p.skeleton_);
		non_manifold_provider_->emit_attribute_changed(*p.skeleton_, p.skeleton_position_.get());
	}

	void start_spheres_update(PointsParameters& p)
	{
		const Scalar convergence_eps = Scalar(1e-10);
		const uint32 max_post_convergence_iterations = 10;
		const uint32 max_iterations_without_autosplit = 300;
		const uint32 max_iterations_after_reaching_max_spheres = 100;
		p.running_ = true;
		p.iteration_count_ = 0;
		p.total_error_diff_ = 0.0;
		p.last_total_error_ = std::numeric_limits<Scalar>::max();
		p.manual_stop_requested_ = false;
		p.pending_full_refresh_after_auto_stop_ = false;
		p.force_full_refresh_on_stop_ = false;
		p.stop_reason_.clear();
		p.last_update_stage_ = "thread_start";

		launch_thread([&, convergence_eps, max_post_convergence_iterations, max_iterations_without_autosplit,
						  max_iterations_after_reaching_max_spheres]() {
			try
			{
				bool convergence_reached = false;
				bool stopped_by_automatic_condition = false;
				uint32 post_convergence_iterations = 0;
				bool target_reached_reported = false;
				bool max_spheres_reached_once = false;
				uint32 post_max_spheres_iterations = 0;
				auto start = std::chrono::high_resolution_clock::now();
				while (true)
				{
					{
						std::lock_guard<std::mutex> lock(p.mutex_);
						p.last_update_stage_ = "update_spheres";
						update_spheres(p);
						p.iteration_count_++;
					}
					if (p.slow_down_)
						std::this_thread::sleep_for(std::chrono::microseconds(1000000 / p.update_rate_));
					else
						std::this_thread::yield();

					if (p.auto_stop_)
					{
						const bool converged = (p.total_error_diff_ < convergence_eps);
						if (converged)
						{
							if (!convergence_reached)
							{
								convergence_reached = true;
								post_convergence_iterations = 0;
								std::cout << "Auto stop: error converged (Diff < " << convergence_eps
										  << "), start post-convergence countdown (" << max_post_convergence_iterations
										  << ")." << std::endl;
							}
							else
							{
								++post_convergence_iterations;
							}

							bool reached_target = false;
							switch (p.auto_split_mode_)
							{
							case ERROR_THRESHOLD:
								if (p.max_error_ < p.auto_split_error_threshold_)
									reached_target = true;
								break;
							case MAX_NB_SPHERES:
								if (p.nb_spheres_ >= p.auto_split_max_nb_spheres_)
									reached_target = true;
								break;
							}

							if (reached_target)
							{
								if (!target_reached_reported)
								{
									std::cout << "Auto stop: target reached under current auto-split mode." << std::endl;
									target_reached_reported = true;
								}
							}
							if (post_convergence_iterations >= max_post_convergence_iterations)
							{
								std::cout << "Auto stop: reached max post-convergence iterations ("
										  << max_post_convergence_iterations << ")." << std::endl;
								stopped_by_automatic_condition = true;
								p.stopping_ = true;
							}
						}
						else if (convergence_reached)
						{
							convergence_reached = false;
							post_convergence_iterations = 0;
							target_reached_reported = false;
						}
					}

					if (p.auto_split_)
					{
						if (p.nb_spheres_ >= p.auto_split_max_nb_spheres_)
						{
							if (!max_spheres_reached_once)
							{
								max_spheres_reached_once = true;
								post_max_spheres_iterations = 0;
								std::cout << "Auto split stop: reached max spheres (" << p.auto_split_max_nb_spheres_
										  << "), start post-max countdown (" << max_iterations_after_reaching_max_spheres
										  << ")." << std::endl;
							}
							else
							{
								++post_max_spheres_iterations;
							}

							if (post_max_spheres_iterations >= max_iterations_after_reaching_max_spheres)
							{
								std::cout << "Auto split stop: reached max post-max-sphere iterations ("
										  << max_iterations_after_reaching_max_spheres << ")." << std::endl;
								stopped_by_automatic_condition = true;
								p.stopping_ = true;
							}
						}
					}
					else if (p.iteration_count_ >= max_iterations_without_autosplit)
					{
						std::cout << "Stop: reached max iterations without auto split ("
								  << max_iterations_without_autosplit << ")." << std::endl;
						stopped_by_automatic_condition = true;
						p.stopping_ = true;
					}

					std::cout << "Iteration: " << p.iteration_count_ << " | Spheres: " << p.nb_spheres_
							  << " | Error: " << p.total_error_ << " | Diff: " << p.total_error_diff_ << std::endl;

					if (p.stopping_)
					{
						if (!p.stop_reason_.empty())
							std::cerr << "Stop reason: " << p.stop_reason_ << std::endl;
						const bool should_full_refresh = p.force_full_refresh_on_stop_ ||
														 (stopped_by_automatic_condition && !p.manual_stop_requested_);
						p.stopping_ = false;
						p.running_ = false;
						p.manual_stop_requested_ = false;
						p.pending_full_refresh_after_auto_stop_ = should_full_refresh;
						p.force_full_refresh_on_stop_ = false;
						break;
					}
				}
				auto end = std::chrono::high_resolution_clock::now();
				std::cout << "Sphere optimizations time: " << std::chrono::duration<Scalar>(end - start).count() << "s"
						  << std::endl;
				std::cout << "Nb iterations: " << p.iteration_count_ << std::endl;
			}
			catch (const std::exception& e)
			{
				std::lock_guard<std::mutex> lock(p.mutex_);
				p.stop_reason_ =
					std::string("Exception in sphere update thread at stage '") + p.last_update_stage_ + "': " + e.what();
				p.running_ = false;
				p.stopping_ = false;
				p.manual_stop_requested_ = false;
				p.pending_full_refresh_after_auto_stop_ = false;
				p.force_full_refresh_on_stop_ = false;
				std::cerr << "[SphereUpdateThread] " << p.stop_reason_ << std::endl;
			}
			catch (...)
			{
				std::lock_guard<std::mutex> lock(p.mutex_);
				p.stop_reason_ =
					std::string("Unknown exception in sphere update thread at stage '") + p.last_update_stage_ + "'.";
				p.running_ = false;
				p.stopping_ = false;
				p.manual_stop_requested_ = false;
				p.pending_full_refresh_after_auto_stop_ = false;
				p.force_full_refresh_on_stop_ = false;
				std::cerr << "[SphereUpdateThread] " << p.stop_reason_ << std::endl;
			}
		});

		app_.start_timer(100, [&]() -> bool { return !p.running_ && !p.pending_full_refresh_after_auto_stop_; });
	}

	void stop_spheres_update(PointsParameters& p)
	{
		p.manual_stop_requested_ = true;
		p.stopping_ = true;
	}

	void key_press_event(View* view, int32 key_code) override
	{
		if (!selected_points_)
			return;
		PointsParameters& p = points_parameters_[selected_points_];

		if (key_code == GLFW_KEY_S && p.show_knn_hover_)
		{
			if (p.knn_hover_locked_)
			{
				clear_knn_hover(p);
				p.knn_hover_locked_ = false;
			}
			else
			{
				update_knn_hover(view, p, view->mouse_x(), view->mouse_y(), true);
				if (p.hovered_sample_.is_valid())
					p.knn_hover_locked_ = true;
			}
			view->request_update();
			return;
		}

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

	void mouse_move_event(View* view, int32 x, int32 y) override
	{
		if (!selected_points_)
			return;
		// No hover update: highlight is toggled by key press.
	}

	void clear_knn_hover(PointsParameters& p)
	{
		if (!p.samples_mesh_ || !p.samples_knn_color_)
			return;
		if (p.hovered_color_indices_.empty())
			return;
		for (size_t i = 0; i < p.hovered_color_indices_.size(); ++i)
		{
			uint32 v_idx = p.hovered_color_indices_[i];
			if (v_idx < nb_cells<PVertex>(*p.samples_mesh_))
				(*p.samples_knn_color_)[v_idx] = p.hovered_color_backup_[i];
		}
		p.hovered_color_indices_.clear();
		p.hovered_color_backup_.clear();
		p.hovered_sample_ = PVertex();
		p.knn_hover_locked_ = false;
		points_provider_->emit_attribute_changed(*p.samples_mesh_, p.samples_knn_color_.get());
	}

	void update_knn_hover(View* view, PointsParameters& p, int32 x, int32 y, bool force)
	{
		if (!p.show_knn_hover_)
			return;
		if (!force && !p.knn_hover_locked_)
			return;
		if (!p.samples_mesh_ || !p.samples_kdtree_ || !p.samples_position_ || !p.samples_knn_ ||
			!p.samples_knn_color_)
		{
			clear_knn_hover(p);
			return;
		}

		rendering::GLVec3d near_ = view->unproject(x, y, 0.0);
		rendering::GLVec3d far_d = view->unproject(x, y, 1.0);
		Vec3 A{near_.x(), near_.y(), near_.z()};
		Vec3 B{far_d.x(), far_d.y(), far_d.z()};
		Vec3 D = (B - A);
		const Scalar d2 = D.squaredNorm();
		if (d2 <= Scalar(1e-12))
		{
			clear_knn_hover(p);
			return;
		}

		Scalar best_d2 = std::numeric_limits<Scalar>::max();
		PVertex best_v;
		foreach_cell(*p.samples_mesh_, [&](PVertex v) -> bool {
			uint32 v_idx = index_of(*p.samples_mesh_, v);
			const Vec3& pos = (*p.samples_position_)[v_idx];
			Scalar dist2 = geometry::squared_distance_line_point(A, B, pos);
			if (dist2 < best_d2)
			{
				best_d2 = dist2;
				best_v = v;
			}
			return true;
		});

		if (!best_v.is_valid())
		{
			clear_knn_hover(p);
			return;
		}

		if (p.hovered_sample_.is_valid() && best_v == p.hovered_sample_)
			return;

		clear_knn_hover(p);

		uint32 best_idx = index_of(*p.samples_mesh_, best_v);
		p.hovered_sample_ = best_v;

		p.hovered_color_indices_.clear();
		p.hovered_color_backup_.clear();

		auto backup_color = [&](uint32 idx) {
			p.hovered_color_indices_.push_back(idx);
			p.hovered_color_backup_.push_back((*p.samples_knn_color_)[idx]);
		};

		backup_color(best_idx);
		(*p.samples_knn_color_)[best_idx] = Vec4(1.0, 0.2, 0.2, 1.0);

		const auto& neighbors = (*p.samples_knn_)[best_idx];
		for (PVertex vn : neighbors)
		{
			uint32 vn_idx = index_of(*p.samples_mesh_, vn);
			backup_color(vn_idx);
			(*p.samples_knn_color_)[vn_idx] = Vec4(1.0, 1.0, 0.2, 1.0);
		}

		points_provider_->emit_attribute_changed(*p.samples_mesh_, p.samples_knn_color_.get());
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
			ImGui::InputFloat("Alpha", &p.alpha_, 0.001f, 0.1f, "%.4f");
			ImGui::InputInt("Num Samples", &p.num_alpha_samples_, 1000, 10000);
			ImGui::InputInt("Num Inside Samples", &p.num_alpha_inside_samples_, 1000, 10000);
			ImGui::InputInt("Inside KNN K", &p.alpha_inside_knn_k_, 1, 10);
			ImGui::InputInt("Inside Min / Sphere", &p.alpha_inside_min_points_per_sphere_, 1, 50);
			ImGui::InputFloat("Inside Poisson Radius", &p.alpha_inside_poisson_radius_, 0.0001f, 0.001f, "%.6f");
			ImGui::InputFloat("Inside Topup BBox Expand", &p.alpha_inside_topup_bbox_expand_, 0.0001f, 0.001f, "%.6f");
			ImGui::InputFloat("Grid Cell Size", &p.grid_cell_size_, 0.001f, 0.01f, "%.4f");
			ImGui::InputInt("Eval Batch Size", &p.batch_size_, 256, 1024);
			ImGui::InputInt("Ray Batch Size", &p.ray_sampler_batch_size_, 256, 2048);
			ImGui::InputInt("Max Iterations", &p.udf_max_iterations_, 1000, 8000);
			ImGui::InputFloat("Tolerance", &p.tol_, 0.0f, 0.0f, "%.6f");
			ImGui::InputInt("KNN K", &p.knn_k_, 1, 5);
			ImGui::InputInt("Poisson Target", &p.poisson_eliminate_target_samples_, 1000, 10000);
			if (ImGui::Button("Recompute Sample KNN"))
			{
				std::lock_guard<std::mutex> lock(p.mutex_);
				recompute_sample_knn_graph(p, p.knn_k_);
			}
			ImGui::SameLine();
			if (ImGui::Button("Poisson Eliminate Samples"))
			{
				std::lock_guard<std::mutex> lock(p.mutex_);
				const size_t target = static_cast<size_t>(std::max(1, p.poisson_eliminate_target_samples_));
				apply_poisson_eliminate_samples(p, target);
			}
			if (ImGui::Button("Apply Sampling Filtering"))
			{
				std::lock_guard<std::mutex> lock(p.mutex_);
				apply_sampling_preprocess_filtering(p);
			}
			if (ImGui::Checkbox("Hover KNN", &p.show_knn_hover_))
			{
				if (!p.show_knn_hover_)
				{
					clear_knn_hover(p);
					p.knn_hover_locked_ = false;
					if (pcr_ && p.samples_mesh_)
						pcr_->set_vertex_color(*app_.current_view(), *p.samples_mesh_, p.samples_color_);
				}
				else
				{
					if (pcr_ && p.samples_mesh_)
						pcr_->set_vertex_color(*app_.current_view(), *p.samples_mesh_, p.samples_knn_color_);
				}
			}

			if (ImGui::Button("Sample UDF"))
			{
				load_alpha_samples_to_mesh(p, static_cast<size_t>(p.num_alpha_samples_));
				p.fitting_data_computed_ = false;
			}
			ImGui::SameLine();
			if (ImGui::Button("Top-up Inside"))
			{
				std::lock_guard<std::mutex> lock(p.mutex_);
				topup_alpha_inside_samples_from_current_mesh(p, "manual_button");
				ensure_mixed_optimization_cloud_ready(p);
			}
			ImGui::SameLine();
			if (ImGui::Button("Clear Samples"))
			{
				if (p.samples_mesh_)
					points_provider_->clear_mesh(*p.samples_mesh_);
				clear_alpha_inside_samples(p);
				if (p.opt_mesh_)
					points_provider_->clear_mesh(*p.opt_mesh_);
				p.fitting_data_computed_ = false;
				p.samples_jitter_backup_valid_ = false;
				p.samples_position_backup_.clear();
				p.samples_normal_backup_.clear();
				p.last_surface_count_ = 0;
				p.last_inside_count_ = 0;
				p.last_mixed_count_ = 0;
			}

			if (p.samples_mesh_)
				ImGui::Text("Number of samples: %zu", nb_cells<PVertex>(*p.samples_mesh_));
			if (p.alpha_inside_mesh_)
				ImGui::Text("Inside samples: %zu", nb_cells<PVertex>(*p.alpha_inside_mesh_));
			if (p.opt_mesh_)
				ImGui::Text("Mixed optimization points: %zu", nb_cells<PVertex>(*p.opt_mesh_));
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
			ImGui::Separator();
			static uint32 init_max_nb_spheres = 1;
			if (ImGui::Button("Compute Fitting Data (MA: Displacement)"))
			{
				compute_fitting_data_with_initial_ma_mode(p, INITIAL_MA_DISPLACEMENT, true);
				{
					std::lock_guard<std::mutex> lock(p.mutex_);
					init_spheres(p, init_max_nb_spheres);
				}
				update_render_data(p);
			}
			if (ImGui::Button("Compute Fitting Data (MA: Shrinking Ball)"))
			{
				compute_fitting_data_with_initial_ma_mode(p, INITIAL_MA_SHRINKING_BALL, true);
				{
					std::lock_guard<std::mutex> lock(p.mutex_);
					init_spheres(p, init_max_nb_spheres);
				}
				update_render_data(p);
			}
			const bool sphere_fit_ready = p.fitting_data_computed_;
			if (sphere_fit_ready)
			{
				// Sphere Fitting
				if (ImGui::CollapsingHeader("Sphere Fitting", ImGuiTreeNodeFlags_DefaultOpen))
				{
					ImGui::SliderFloat("Init dilation constant", &p.init_dilation_constant_, 0.001f, 0.01f, "%.4f");
					ImGui::InputScalar("Init nb spheres", ImGuiDataType_U32, &init_max_nb_spheres);
					if (ImGui::Button("Init spheres"))
					{
						std::lock_guard<std::mutex> lock(p.mutex_);
						init_spheres(p, init_max_nb_spheres);
						update_render_data(p);
					}

					static bool sync_lambda = true;
					ImGui::Checkbox("Sync lambda", &sync_lambda);
					if (ImGui::SliderFloat("update lambda", &p.sqem_update_lambda_, 0.0f, 4.0f, "%.6f"))
					{
						if (sync_lambda)
							p.sqem_clustering_lambda_ = p.sqem_update_lambda_;
					}
					if (ImGui::SliderFloat("clustering lambda", &p.sqem_clustering_lambda_, 0.0f, 4.0f, "%.6f"))
					{
						if (sync_lambda)
							p.sqem_update_lambda_ = p.sqem_clustering_lambda_;
					}
					const bool fix_r_mode = (p.distance_mode_ == LINE_QUADRIC_DISTANCE);
					if (!fix_r_mode)
						ImGui::BeginDisabled();
					ImGui::SliderFloat("fix radius scale", &p.sqem_fix_radius_scale_, 1.0f, 5.0f, "%.3f");
					if (!fix_r_mode)
						ImGui::EndDisabled();

					if (p.neural_udf_loaded_)
					{
						ImGui::Checkbox("UDF center term", &p.udf_center_enabled_);
						if (p.udf_center_enabled_)
							ImGui::SliderFloat("UDF lambda", &p.udf_center_lambda_, 0.0f, 2.0f, "%.6f");
					}
					else
					{
						ImGui::TextColored(ImVec4(1, 1, 0, 1), "UDF center term requires loaded model.");
					}

					ImGui::RadioButton("Line Quadric (fix r)", (int*)&p.distance_mode_, LINE_QUADRIC_DISTANCE);
					ImGui::SameLine();
					ImGui::RadioButton("Line Quadric (free r)", (int*)&p.distance_mode_, LINE_QUADRIC_DISTANCE_FREE_RADIUS);
					if (ImGui::Button(p.lock_skeleton_connectivity_ ? "Skeleton: Locked" : "Skeleton: Unlocked"))
						p.lock_skeleton_connectivity_ = !p.lock_skeleton_connectivity_;
					int connectivity_mode_ui = static_cast<int>(p.connectivity_mode_);
					const char* connectivity_labels[] = {"Mode A: Mixed KNN", "Mode B: Inside Power"};
					if (ImGui::Combo("Connectivity Mode", &connectivity_mode_ui, connectivity_labels, 2))
						p.connectivity_mode_ = static_cast<ConnectivityMode>(connectivity_mode_ui);
					ImGui::Checkbox("Top-up in training loop", &p.topup_inside_before_reassign_in_loop_);
					ImGui::Checkbox("Inside skeleton top-up first", &p.inside_skeleton_build_topup_first_);
					const uint32 ui_surface_count =
						p.samples_mesh_ ? static_cast<uint32>(nb_cells<PVertex>(*p.samples_mesh_)) : uint32(0);
					const uint32 ui_inside_count =
						p.alpha_inside_mesh_ ? static_cast<uint32>(nb_cells<PVertex>(*p.alpha_inside_mesh_)) : uint32(0);
					const uint32 ui_mixed_count =
						p.opt_mesh_ ? static_cast<uint32>(nb_cells<PVertex>(*p.opt_mesh_)) : uint32(0);
					ImGui::Text("surface=%u inside=%u mixed=%u pending_delete=%u", ui_surface_count, ui_inside_count,
								ui_mixed_count, p.pending_delete_count_);

					if (ImGui::Button("Update spheres"))
					{
						if (!p.running_)
						{
							std::lock_guard<std::mutex> lock(p.mutex_);
							update_spheres(p);
						}
					}

					if (ImGui::Button("Compute clusters"))
					{
						if (!p.running_)
						{
							std::lock_guard<std::mutex> lock(p.mutex_);
							prepare_hybrid_iteration_data(p, p.topup_inside_before_reassign_in_loop_);
							compute_clusters(p);
							compute_spheres_error(p);
							if (!p.lock_skeleton_connectivity_)
								compute_skeleton_current_mode(p, true);
							update_render_data(p);
						}
					}
					ImGui::SameLine();
					if (ImGui::Button(p.use_local_clusters_ ? "Cluster Mode: Neighbor" : "Cluster Mode: Global"))
						p.use_local_clusters_ = !p.use_local_clusters_;
					if (p.use_local_clusters_)
					{
						ImGui::InputScalar("Neighbor refresh/iter", ImGuiDataType_U32, &p.local_cluster_connectivity_refresh_interval_);
						if (p.local_cluster_connectivity_refresh_interval_ < 1)
							p.local_cluster_connectivity_refresh_interval_ = 1;
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

					if (ImGui::Button("Build Skeleton (Current Mode)"))
					{
						if (!p.running_)
						{
							std::lock_guard<std::mutex> lock(p.mutex_);
							compute_skeleton_current_mode(p);
							update_render_data(p);
							non_manifold_provider_->emit_connectivity_changed(*p.skeleton_);
						}
					}
					ImGui::SameLine();
					if (ImGui::Button("Inside Connectivity Skeleton"))
					{
						if (!p.running_)
						{
							std::lock_guard<std::mutex> lock(p.mutex_);
							compute_skeleton_inside_connectivity(p, false, p.inside_skeleton_build_topup_first_);
							update_render_data(p);
							non_manifold_provider_->emit_connectivity_changed(*p.skeleton_);
						}
					}
					ImGui::InputFloat("Edge UDF |0| tol", &p.skeleton_edge_udf_zero_tol_, 0.0f, 0.0f, "%.6f");
					ImGui::InputFloat("Face UDF |0| tol", &p.skeleton_face_udf_zero_tol_, 0.0f, 0.0f, "%.6f");
					ImGui::Checkbox("Face score normalize(area)", &p.skeleton_face_score_normalize_by_area_);
					if (ImGui::Button("Geometry filter"))
					{
						if (!p.running_)
						{
							std::lock_guard<std::mutex> lock(p.mutex_);
							skeleton_geometry_filter(p);
						}
					}
					ImGui::SameLine();
					if (ImGui::Button("Face stage filter"))
					{
						if (!p.running_)
						{
							std::lock_guard<std::mutex> lock(p.mutex_);
							run_topology_stage_filter_from_snapshot(p, false);
						}
					}
					ImGui::SameLine();
					if (ImGui::Button("Face stage step (highlight next)"))
					{
						if (!p.running_)
						{
							std::lock_guard<std::mutex> lock(p.mutex_);
							run_face_stage_single_step_and_highlight_next(p);
						}
					}
					if (ImGui::Button("Reset stage base (current skeleton)"))
					{
						if (!p.running_)
						{
							std::lock_guard<std::mutex> lock(p.mutex_);
							capture_topology_stage_snapshot(p);
						}
					}
					
					if (ImGui::Button("Boundary tet delete"))
					{
						if (!p.running_)
						{
							std::lock_guard<std::mutex> lock(p.mutex_);
							boundary_tet_face_delete_only(p);
						}
					}
					ImGui::SameLine();
					if (ImGui::Button("Tet face UDF colormap (Integral)"))
					{
						if (!p.running_)
						{
							std::lock_guard<std::mutex> lock(p.mutex_);
							p.skeleton_face_score_normalize_by_area_ = false;
							visualize_skeleton_tet_face_udf(p, false);
						}
					}
					ImGui::SameLine();
					if (ImGui::Button("Tet face UDF colormap (Normalized)"))
					{
						if (!p.running_)
						{
							std::lock_guard<std::mutex> lock(p.mutex_);
							p.skeleton_face_score_normalize_by_area_ = true;
							visualize_skeleton_tet_face_udf(p, true);
						}
					}
					if (ImGui::Button("Edge UDF score colormap (Gauss3)"))
					{
						if (!p.running_)
						{
							std::lock_guard<std::mutex> lock(p.mutex_);
							visualize_skeleton_edge_udf_scores(p);
						}
					}
					ImGui::SameLine();
					if (ImGui::Button("Edge score simple tet delete"))
					{
						if (!p.running_)
						{
							std::lock_guard<std::mutex> lock(p.mutex_);
							run_edge_score_simple_tet_topology_fix(p);
						}
					}
					ImGui::SameLine();
					if (ImGui::Button("Edge score simple tet step"))
					{
						if (!p.running_)
						{
							std::lock_guard<std::mutex> lock(p.mutex_);
							run_edge_score_simple_tet_topology_fix_single_step(p);
						}
					}
					if (ImGui::Button("Edge score non-simple tet delete"))
					{
						if (!p.running_)
						{
							std::lock_guard<std::mutex> lock(p.mutex_);
							run_edge_score_nonsimple_tet_topology_fix(p);
						}
					}
					ImGui::SameLine();
					if (ImGui::Button("Edge score non-simple tet step"))
					{
						if (!p.running_)
						{
							std::lock_guard<std::mutex> lock(p.mutex_);
							run_edge_score_nonsimple_tet_topology_fix_single_step(p);
						}
					}
					
					if (ImGui::Button("Boundary tet mask"))
					{
						if (!p.running_)
						{
							std::lock_guard<std::mutex> lock(p.mutex_);
							mark_boundary_tets_color(p);
						}
					}
					ImGui::SameLine();
					if (ImGui::Button("K5 delete"))
					{
						if (!p.running_)
						{
							std::lock_guard<std::mutex> lock(p.mutex_);
							run_k5_face_deletion(p);
						}
					}
					ImGui::SameLine();
					if (ImGui::Button("Topology fix full pipeline"))
					{
						if (!p.running_)
						{
							std::lock_guard<std::mutex> lock(p.mutex_);
							run_complete_topology_fix_pipeline(p);
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
						{
							ImGui::SliderFloat("Threshold", &p.auto_split_error_threshold_, 0.0f, 1.0f, "%.6f",
											   ImGuiSliderFlags_Logarithmic);
							ImGui::SliderFloat("Split ratio", &p.auto_split_ratio_, 0.0f, 1.0f, "%.3f");
							ImGui::InputScalar("Max split/iter", ImGuiDataType_U32, &p.auto_split_max_per_iter_error_);
						}
						else
						{
							ImGui::InputScalar("Nb spheres", ImGuiDataType_U32, &p.auto_split_max_nb_spheres_);
							ImGui::SliderFloat("Split ratio", &p.auto_split_ratio_, 0.0f, 1.0f, "%.3f");
							ImGui::InputScalar("Max split/iter", ImGuiDataType_U32, &p.auto_split_max_per_iter_max_);
						}
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
			else
			{
				ImGui::TextColored(ImVec4(1, 1, 0, 1),
								   "Compute fitting data to enable sphere fitting.");
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
	std::shared_ptr<SAttribute<Vec3>> surface_vertex_normal_ = nullptr;
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














