#ifndef CGOGN_MODULE_UDF_TRAINING_H_
#define CGOGN_MODULE_UDF_TRAINING_H_

#include <cgogn/ui/app.h>
#include <cgogn/ui/module.h>

#include <cgogn/core/ui_modules/mesh_provider.h>
#include <cgogn/geometry/algos/area.h>
#include <cgogn/geometry/algos/fitting.h>
#include <cgogn/geometry/algos/length.h>
#include <cgogn/geometry/algos/medial_axis.h>
#include <cgogn/geometry/algos/udf/alpha_projection.h>
#include <cgogn/geometry/algos/udf/alpha_sampling.h>
#include <cgogn/geometry/algos/udf/neural_alpha_projection.h>
#include <cgogn/geometry/algos/udf/neural_field_query.h>
#include <cgogn/geometry/algos/udf/spatial_query.h>
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
#include <cgogn/io/point/export_options.h>
#include <cgogn/io/surface/export_options.h>
#include <cgogn/io/surface/ply.h>
#include <cgogn/ui/portable-file-dialogs.h>


#include <cgogn/rendering/ui_modules/point_cloud_render.h>
#include <cgogn/rendering/ui_modules/surface_render.h>

#include <Eigen/Dense>
#include <Eigen/Eigenvalues>
#include <Eigen/Sparse>

#include <libacc/bvh_tree.h>
#include <libacc/kd_tree.h>

#include <GLFW/glfw3.h>
#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <filesystem>
#include <limits>
#include <mutex>
#include <numeric>
#include <queue>
#include <random>
#include <set>
#include <string>
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
using geometry::AlphaProjectionInfo;
using geometry::Mat3;
using geometry::Mat4;
using geometry::Scalar;
using geometry::Spherical_Quadric;
using geometry::Quadric;
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

	template <typename T>
	using PAttribute = typename mesh_traits<POINTS>::template Attribute<T>;
	template <typename T>
	using NMAttribute = typename mesh_traits<NONMANIFOLD>::template Attribute<T>;
	template <typename T>
	using SAttribute = typename mesh_traits<SURFACE>::template Attribute<T>;
	using NMFaceKey = std::array<uint32, 3>;

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
	enum OutputVerbosity : uint32
	{
		OUTPUT_MUTE,
		OUTPUT_NORMAL,
		OUTPUT_VERBOSE
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

		// Nueral UDF
		bool neural_udf_loaded_ = false;
		torch::jit::Module neural_udf_model_;
		std::string neural_udf_model_path_ = "";
		NeuralModelType neural_model_type_ = NEURAL_MODEL_UDF;
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

		acc::KDTree<3, uint32>* samples_kdtree_ = nullptr; // KDTree of alpha-expanding samples
		std::vector<PVertex> samples_kdtree_vertices_;	   // Vertices of alpha-expanding samples in KDTree order
		acc::KDTree<3, uint32>* samples_ma_kdtree_ = nullptr; // KDTree of sample medial-axis positions
		std::vector<PVertex> samples_ma_kdtree_vertices_;	   // Vertices in MA KDTree order
		acc::KDTree<3, uint32>* input_kdtree_ = nullptr; // KDTree of input points
		std::vector<PVertex> input_kdtree_vertices_;	 // Vertices of input points in KDTree order
		std::unique_ptr<SpatialGrid> samples_spatial_grid_ = nullptr;
		geometry::NeuralProjectionWorkspace neural_projection_workspace_;
		bool filter_positive_projection_ = false;
		bool ma_flip_prune_enabled_ = false;
		Scalar ma_flip_prune_alpha_factor_ = Scalar(1);
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
		std::shared_ptr<PAttribute<Scalar>> spheres_error_ = nullptr;
		std::shared_ptr<PAttribute<Scalar>> spheres_error_not_normalized_ = nullptr;
		std::shared_ptr<PAttribute<NMVertex>> spheres_skeleton_vertex_ = nullptr;

		// Skeleton
		NONMANIFOLD* skeleton_ = nullptr;
		std::shared_ptr<NMAttribute<Vec3>> skeleton_position_ = nullptr;
		std::shared_ptr<NMAttribute<Scalar>> skeleton_radius_ = nullptr;
		std::shared_ptr<NMAttribute<std::set<std::size_t>>> incident_tets_ = nullptr;

		std::shared_ptr<NMAttribute<Vec3>> skeleton_face_color_ = nullptr;
		std::shared_ptr<NMAttribute<uint32>> skeleton_face_component_id_ = nullptr;
		std::shared_ptr<NMAttribute<Vec3>> skeleton_face_component_color_ = nullptr;
		std::shared_ptr<NMAttribute<PVertex>> skeleton_source_sphere_ = nullptr;
		std::shared_ptr<NMAttribute<Vec3>> skeleton_edge_color_ = nullptr;
		std::shared_ptr<NMAttribute<Vec3>> skeleton_edge_non_manifold_color_ = nullptr;
		std::shared_ptr<NMAttribute<uint32>> edge_degree_ = nullptr;


		std::unordered_map<std::size_t, Tet> skeleton_tets_;

		float32 filter_radius_threshold_ = 0.0f;

		bool error_as_spheres_color_ = false;
		float32 spheres_transparency_ = 0.5f;
		float32 sqem_update_lambda_line_plane_ = 0.20f;
		bool export_samples_mesh_selected_ = true;
		bool export_samples_mesh_normal_color_selected_ = false;
		bool export_samples_spheres_selected_ = true;
		bool export_skeleton_selected_ = true;

		// Filtering
		float32 target_radius_ = 0.1f;
		float32 radius_tolerance_ = 0.01f;
		float32 init_dilation_constant_ = 0.001f;
		// Sampling Parameters
		float alpha_ = 0.005f;
		float sample_radius_ = 0.0025f;
		int knn_k_ = 10;
		int seed_ = 42;
		float sample_grid_cell_size_ = 0.0025f;
		int batch_size_ = 1310640;	   // NN evaluation batch
		int ray_sampler_batch_size_ = 8192; // Rays per sampling iteration
		float tol_ = 1e-6f; // convergence tolerance
		bool recompute_sample_normals_after_sampling_ = false;

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
		OutputVerbosity output_verbosity_ = OUTPUT_NORMAL;

		// Threading
		uint32 iteration_count_ = 0;
		std::mutex mutex_;
		bool running_ = false;
		bool stopping_ = false;
		bool preview_render_during_sphere_update_ = true;
		bool pending_full_refresh_after_stop_ = false;
		bool slow_down_ = true;
		uint32 update_rate_ = 20;

		~PointsParameters()
		{
			if (samples_kdtree_)
				delete samples_kdtree_;
			if (samples_ma_kdtree_)
				delete samples_ma_kdtree_;
			if (input_kdtree_)
				delete input_kdtree_;
		}
	};

	struct Tet
	{
		NMFace faces[4];
	};

public:
	UDFTraining(const App& app) : ViewModule(app, "UDFTraining")
	{
	}

	~UDFTraining()
	{
	}

	struct HeadlessBenchmarkOptions
	{
		bool verbose_ = true;
		OutputVerbosity output_verbosity_ = OUTPUT_NORMAL;
		bool ma_flip_prune_enabled_ = true;
		float ma_flip_prune_alpha_factor_ = 1.0f;
		float32 filter_radius_threshold_ = 0.0f;
		float32 sqem_update_lambda_line_plane_ = 0.20f;
		float32 init_dilation_constant_ = 0.001f;
		float alpha_ = 0.005f;
		float sample_radius_ = 0.0025f;
		int knn_k_ = 10;
		int seed_ = 42;
		int batch_size_ = 1310640;
		int ray_sampler_batch_size_ = 4096;
		float tol_ = 1e-6f;
		int udf_max_iterations_ = 3000;
		bool recompute_sample_normals_after_sampling_ = false;
	};

	struct HeadlessOptimizationStats
	{
		float64 optimization_total_ms_ = 0.0;
		float64 cluster_total_ms_ = 0.0;
		float64 sphere_update_total_ms_ = 0.0;
		float64 error_total_ms_ = 0.0;
		float64 average_iteration_ms_ = 0.0;
		uint32 optimization_iterations_ = 0;
	};

	struct HeadlessCounts
	{
		uint32 input_vertices_ = 0;
		uint32 input_points_ = 0;
		uint32 sample_points_ = 0;
		uint32 final_spheres_ = 0;
		uint32 skeleton_vertices_ = 0;
		uint32 skeleton_edges_ = 0;
		uint32 skeleton_faces_ = 0;
		uint32 optimization_iterations_ = 0;
	};

	static int output_verbosity_index(OutputVerbosity value)
	{
		switch (value)
		{
		case OUTPUT_MUTE:
			return 0;
		case OUTPUT_VERBOSE:
			return 2;
		case OUTPUT_NORMAL:
		default:
			return 1;
		}
	}

	static OutputVerbosity output_verbosity_from_index(int index)
	{
		switch (index)
		{
		case 0:
			return OUTPUT_MUTE;
		case 2:
			return OUTPUT_VERBOSE;
		case 1:
		default:
			return OUTPUT_NORMAL;
		}
	}

	bool is_basic_logging_enabled(const PointsParameters& p) const
	{
		return p.output_verbosity_ != OUTPUT_MUTE;
	}

	bool is_verbose_logging_enabled(const PointsParameters& p) const
	{
		return p.output_verbosity_ == OUTPUT_VERBOSE;
	}

	bool is_basic_logging_enabled(OutputVerbosity output_verbosity) const
	{
		return output_verbosity != OUTPUT_MUTE;
	}

	bool is_verbose_logging_enabled(OutputVerbosity output_verbosity) const
	{
		return output_verbosity == OUTPUT_VERBOSE;
	}

	template <typename... Args>
	void log_basic(const PointsParameters& p, Args&&... args) const
	{
		if (!is_basic_logging_enabled(p))
			return;
		(std::cout << ... << std::forward<Args>(args));
	}

	template <typename... Args>
	void log_basic(OutputVerbosity output_verbosity, Args&&... args) const
	{
		if (!is_basic_logging_enabled(output_verbosity))
			return;
		(std::cout << ... << std::forward<Args>(args));
	}

	template <typename... Args>
	void log_basic(const PointsParameters& p, Args&&... args, std::ostream& (*manip)(std::ostream&)) const
	{
		if (!is_basic_logging_enabled(p))
			return;
		(std::cout << ... << std::forward<Args>(args));
		manip(std::cout);
	}

	template <typename... Args>
	void log_verbose(const PointsParameters& p, Args&&... args) const
	{
		if (!is_verbose_logging_enabled(p))
			return;
		(std::cout << ... << std::forward<Args>(args));
	}

	template <typename... Args>
	void log_verbose(OutputVerbosity output_verbosity, Args&&... args) const
	{
		if (!is_verbose_logging_enabled(output_verbosity))
			return;
		(std::cout << ... << std::forward<Args>(args));
	}

	template <typename... Args>
	void log_verbose(const PointsParameters& p, Args&&... args, std::ostream& (*manip)(std::ostream&)) const
	{
		if (!is_verbose_logging_enabled(p))
			return;
		(std::cout << ... << std::forward<Args>(args));
		manip(std::cout);
	}

	template <typename... Args>
	void log_error(const PointsParameters& p, Args&&... args) const
	{
		if (!is_basic_logging_enabled(p))
			return;
		(std::cerr << ... << std::forward<Args>(args));
	}

	template <typename... Args>
	void log_error(OutputVerbosity output_verbosity, Args&&... args) const
	{
		if (!is_basic_logging_enabled(output_verbosity))
			return;
		(std::cerr << ... << std::forward<Args>(args));
	}

	template <typename... Args>
	void log_error(const PointsParameters& p, Args&&... args, std::ostream& (*manip)(std::ostream&)) const
	{
		if (!is_basic_logging_enabled(p))
			return;
		(std::cerr << ... << std::forward<Args>(args));
		manip(std::cerr);
	}

	void headless_prepare_points(POINTS& points)
	{
		init_points_data(points);
	}

	void apply_headless_benchmark_options_prepared(POINTS& points, const HeadlessBenchmarkOptions& options)
	{
		apply_headless_benchmark_options_prepared(points_parameters_[&points], options);
	}

	void apply_headless_benchmark_options_prepared(PointsParameters& p, const HeadlessBenchmarkOptions& options)
	{
		p.output_verbosity_ = options.output_verbosity_;
		if (options.verbose_)
			p.output_verbosity_ = OUTPUT_VERBOSE;
		p.ma_flip_prune_enabled_ = options.ma_flip_prune_enabled_;
		p.ma_flip_prune_alpha_factor_ = options.ma_flip_prune_alpha_factor_;
		p.filter_radius_threshold_ = options.filter_radius_threshold_;
		p.sqem_update_lambda_line_plane_ = options.sqem_update_lambda_line_plane_;
		p.init_dilation_constant_ = options.init_dilation_constant_;
		p.alpha_ = options.alpha_;
		p.sample_radius_ = options.sample_radius_;
		p.knn_k_ = options.knn_k_;
		p.seed_ = options.seed_;
		p.sample_grid_cell_size_ = std::max(p.sample_radius_ / std::sqrt(Scalar(3.0)), Scalar(1e-8));
		p.batch_size_ = options.batch_size_;
		p.ray_sampler_batch_size_ = options.ray_sampler_batch_size_;
		p.tol_ = options.tol_;
		p.udf_max_iterations_ = options.udf_max_iterations_;
		p.recompute_sample_normals_after_sampling_ = options.recompute_sample_normals_after_sampling_;
		if (p.sample_grid_cell_size_ > Scalar(0))
			p.samples_spatial_grid_ = std::make_unique<SpatialGrid>(p.sample_grid_cell_size_);
		else
			p.samples_spatial_grid_.reset();
	}

	void headless_sample_alpha_level_set_prepared(PointsParameters& p, const HeadlessBenchmarkOptions& options)
	{
		(void)options;
		load_alpha_samples_to_mesh(p);
	}

	void headless_sample_alpha_level_set_prepared(POINTS& points, const HeadlessBenchmarkOptions& options)
	{
		headless_sample_alpha_level_set_prepared(points_parameters_[&points], options);
	}

	void headless_apply_sampling_filtering_prepared(PointsParameters& p)
	{
		apply_sampling_preprocess_filtering(p);
	}

	void headless_apply_sampling_filtering_prepared(POINTS& points)
	{
		headless_apply_sampling_filtering_prepared(points_parameters_[&points]);
	}

	void headless_build_sample_kdtree_and_normals_prepared(PointsParameters& p)
	{
		log_basic(p, "Headless sample normal recompute after sampling: ",
				  (p.recompute_sample_normals_after_sampling_ ? "true" : "false"), '\n');
		build_kdtree(p);
		if (p.recompute_sample_normals_after_sampling_)
		{
			log_basic(p, "Recomputing sampled normals with PCA after sampling...", '\n');
			recompute_samples_normals_pca(p);
		}
		else
		{
			log_basic(p, "Keeping sampled normals from sampling output without PCA recompute.", '\n');
		}
	}

	void headless_build_sample_kdtree_and_normals_prepared(POINTS& points)
	{
		headless_build_sample_kdtree_and_normals_prepared(points_parameters_[&points]);
	}

	void headless_compute_initial_medial_axis_prepared(PointsParameters& p)
	{
		compute_initial_medial_axis(p);
		p.fitting_data_computed_ = true;
	}

	void headless_compute_initial_medial_axis_prepared(POINTS& points)
	{
		headless_compute_initial_medial_axis_prepared(points_parameters_[&points]);
	}

	void headless_compute_fitting_primitives_prepared(PointsParameters& p)
	{
		build_kdtree(p);
		compute_samples_area(p);
		compute_quadrics(p);
	}

	void headless_compute_fitting_primitives_prepared(POINTS& points)
	{
		headless_compute_fitting_primitives_prepared(points_parameters_[&points]);
	}

	void headless_init_spheres_prepared(POINTS& points)
	{
		auto& p = points_parameters_[&points];
		clear(*p.spheres_);
		init_spheres_from_samples(p);
		if (p.samples_sphere_)
			p.samples_sphere_->fill(PVertex());
		compute_clusters_local(p);
	}

	HeadlessOptimizationStats headless_optimize_spheres_prepared(PointsParameters& p, bool verbose = false)
	{
		const OutputVerbosity output_verbosity = verbose ? OUTPUT_VERBOSE : p.output_verbosity_;
		constexpr uint32 max_iterations = 150;
		const Scalar convergence_eps = Scalar(1e-10);
		const uint32 max_post_convergence_iterations = 10;
		HeadlessOptimizationStats stats;
		p.running_ = true;
		p.stopping_ = false;
		p.iteration_count_ = 0;
		p.total_error_diff_ = 0.0;
		p.last_total_error_ = std::numeric_limits<Scalar>::max();
		p.pending_full_refresh_after_stop_ = false;

		bool convergence_reached = false;
		uint32 post_convergence_iterations = 0;
		auto optimization_start = std::chrono::high_resolution_clock::now();
		while (p.iteration_count_ < max_iterations)
		{
			update_spheres(p, &stats);
			++p.iteration_count_;
			log_basic(output_verbosity, "[HeadlessOptimize] iteration=", p.iteration_count_, " spheres=",
					  p.nb_spheres_, " error=", p.total_error_, " diff=", p.total_error_diff_, '\n');

			if (p.total_error_diff_ < convergence_eps)
			{
				if (!convergence_reached)
				{
					convergence_reached = true;
					post_convergence_iterations = 0;
				}
				else
				{
					++post_convergence_iterations;
				}
				if (post_convergence_iterations >= max_post_convergence_iterations)
					break;
			}
			else if (convergence_reached)
			{
				convergence_reached = false;
				post_convergence_iterations = 0;
			}
		}

		auto optimization_end = std::chrono::high_resolution_clock::now();
		stats.optimization_total_ms_ =
			std::chrono::duration<float64, std::milli>(optimization_end - optimization_start).count();
		stats.optimization_iterations_ = p.iteration_count_;
		stats.average_iteration_ms_ = (p.iteration_count_ > 0)
										 ? stats.optimization_total_ms_ / static_cast<float64>(p.iteration_count_)
										 : 0.0;
		p.running_ = false;
		p.stopping_ = false;
		p.pending_full_refresh_after_stop_ = false;
		return stats;
	}

	HeadlessOptimizationStats headless_optimize_spheres_prepared(POINTS& points, bool verbose = false)
	{
		return headless_optimize_spheres_prepared(points_parameters_[&points], verbose);
	}
	void headless_build_skeleton_prepared(PointsParameters& p)
	{
		compute_skeleton(p);
	}

	void headless_build_skeleton_prepared(POINTS& points)
	{
		headless_build_skeleton_prepared(points_parameters_[&points]);
	}

	void headless_run_topology_fix_prepared(PointsParameters& p)
	{
		run_topology_fix_pipeline(p);
	}

	void headless_run_topology_fix_prepared(POINTS& points)
	{
		headless_run_topology_fix_prepared(points_parameters_[&points]);
	}
	void headless_run_completion_residual_prune_prepared(PointsParameters& p)
	{
		run_residual_sheet_prune(p, "[ResidualSheetPrune]");
	}

	void headless_run_completion_residual_prune_prepared(POINTS& points)
	{
		headless_run_completion_residual_prune_prepared(points_parameters_[&points]);
	}

	HeadlessCounts headless_collect_counts(const POINTS& points) const
	{
		HeadlessCounts counts;
		auto it = points_parameters_.find(const_cast<POINTS*>(&points));
		if (it == points_parameters_.end())
			return counts;
		const PointsParameters& p = it->second;
		counts.input_points_ = p.points_ ? nb_cells<PVertex>(*p.points_) : 0;
		counts.input_vertices_ = (selected_surface_ ? nb_cells<SVertex>(*selected_surface_) : counts.input_points_);
		counts.sample_points_ = p.samples_mesh_ ? nb_cells<PVertex>(*p.samples_mesh_) : 0;
		counts.final_spheres_ = p.spheres_ ? nb_cells<PVertex>(*p.spheres_) : 0;
		counts.skeleton_vertices_ = p.skeleton_ ? nb_cells<NMVertex>(*p.skeleton_) : 0;
		counts.skeleton_edges_ = p.skeleton_ ? nb_cells<NMEdge>(*p.skeleton_) : 0;
		counts.skeleton_faces_ = p.skeleton_ ? nb_cells<NMFace>(*p.skeleton_) : 0;
		counts.optimization_iterations_ = p.iteration_count_;
		return counts;
	}

	std::filesystem::path default_training_export_directory(const PointsParameters& p) const
	{
		if (points_provider_ && p.points_)
		{
			const std::string source_filename = points_provider_->mesh_filename(*p.points_);
			if (!source_filename.empty())
			{
				const std::filesystem::path parent = std::filesystem::path(source_filename).parent_path();
				if (!parent.empty())
					return parent;
			}
		}
		return std::filesystem::current_path();
	}

	std::string training_export_basename(const PointsParameters& p) const
	{
		if (points_provider_ && p.points_)
		{
			const std::string mesh_name = points_provider_->mesh_name(*p.points_);
			if (!mesh_name.empty())
				return std::filesystem::path(mesh_name).stem().string();
		}
		return "udf_training";
	}

	bool export_samples_mesh_cluster_ply(PointsParameters& p, const std::string& filename)
	{
		if (!points_provider_ || !p.samples_mesh_ || !p.samples_position_ || !p.samples_sphere_ || !p.spheres_ ||
			!p.spheres_cluster_color_)
			return false;

		auto export_cluster_color = get_attribute<Vec4, PVertex>(*p.samples_mesh_, "cluster_color");
		const bool created_export_cluster_color = !export_cluster_color;
		if (!export_cluster_color)
			export_cluster_color = add_attribute<Vec4, PVertex>(*p.samples_mesh_, "cluster_color");
		if (!export_cluster_color)
			return false;
		foreach_cell(*p.samples_mesh_, [&](PVertex sample) -> bool {
			const uint32 sample_index = index_of(*p.samples_mesh_, sample);
			Vec4 color(0.0, 0.0, 0.0, 1.0);
			if (sample_index != INVALID_INDEX)
			{
				const PVertex sphere = (*p.samples_sphere_)[sample_index];
				if (sphere.is_valid())
					color = value<Vec4>(*p.spheres_, p.spheres_cluster_color_, sphere);
				(*export_cluster_color)[sample_index] = color;
			}
			return true;
		});

		io::PointExportAttributeSelection<POINTS> export_attributes;
		export_attributes.vertex_color_attribute = export_cluster_color;
		points_provider_->save_points_ply_to_file(*p.samples_mesh_, p.samples_position_.get(), filename, export_attributes);
		if (created_export_cluster_color)
			remove_attribute<PVertex>(*p.samples_mesh_, export_cluster_color);
		return true;
	}

	bool export_samples_spheres_ply(const PointsParameters& p, const std::string& filename) const
	{
		if (!points_provider_ || !p.spheres_ || !p.spheres_position_ || !p.spheres_radius_ || !p.spheres_cluster_color_)
			return false;

		io::PointExportAttributeSelection<POINTS> export_attributes;
		export_attributes.vertex_attributes.push_back(p.spheres_radius_);
		export_attributes.vertex_color_attribute = p.spheres_cluster_color_;
		points_provider_->save_points_ply_to_file(*p.spheres_, p.spheres_position_.get(), filename, export_attributes);
		return true;
	}

	bool export_samples_mesh_normal_color_ply(const PointsParameters& p, const std::string& filename) const
	{
		if (!points_provider_ || !p.samples_mesh_ || !p.samples_position_ || !p.samples_normal_color_)
			return false;

		io::PointExportAttributeSelection<POINTS> export_attributes;
		export_attributes.vertex_color_attribute = p.samples_normal_color_;
		points_provider_->save_points_ply_to_file(*p.samples_mesh_, p.samples_position_.get(), filename, export_attributes);
		return true;
	}

	bool export_skeleton_mesh_ply(PointsParameters& p, const std::string& filename)
	{
		if (!p.skeleton_ || !p.skeleton_position_ || !p.skeleton_radius_ || !p.spheres_ || !p.spheres_radius_)
			return false;

		if (!std::filesystem::path(filename).parent_path().empty())
			std::filesystem::create_directories(std::filesystem::path(filename).parent_path());

		foreach_cell(*p.skeleton_, [&](NMVertex v) -> bool {
			const uint32 v_index = index_of(*p.skeleton_, v);
			if (v_index < nb_cells<PVertex>(*p.spheres_))
				(*p.skeleton_radius_)[v_index] = (*p.spheres_radius_)[v_index];
			return true;
		});

		io::SurfaceExportAttributeSelection<NONMANIFOLD> export_attributes;
		export_attributes.vertex_attributes.push_back(p.skeleton_radius_);
		if (p.skeleton_face_color_)
			export_attributes.face_attributes.push_back(p.skeleton_face_color_);
		if (non_manifold_provider_)
			non_manifold_provider_->save_surface_ply_to_file(
				*p.skeleton_, p.skeleton_position_.get(), filename, export_attributes);
		else
			io::export_PLY(*p.skeleton_, p.skeleton_position_.get(), filename, &export_attributes);
		return true;
	}

	bool export_training_ply_bundle(PointsParameters& p, const std::filesystem::path& output_directory)
	{
		const std::string basename = training_export_basename(p);
		const std::filesystem::path samples_mesh_path = output_directory / (basename + "_samples_mesh.ply");
		const std::filesystem::path samples_mesh_normal_color_path =
			output_directory / (basename + "_samples_mesh_normal_color.ply");
		const std::filesystem::path samples_spheres_path = output_directory / (basename + "_samples_spheres.ply");
		const std::filesystem::path skeleton_path = output_directory / (basename + "_skeleton.ply");

		const bool any_selected =
			p.export_samples_mesh_selected_ || p.export_samples_mesh_normal_color_selected_ ||
			p.export_samples_spheres_selected_ || p.export_skeleton_selected_;
		const bool samples_ok =
			!p.export_samples_mesh_selected_ || export_samples_mesh_cluster_ply(p, samples_mesh_path.string());
		const bool samples_normal_color_ok =
			!p.export_samples_mesh_normal_color_selected_ ||
			export_samples_mesh_normal_color_ply(p, samples_mesh_normal_color_path.string());
		const bool spheres_ok =
			!p.export_samples_spheres_selected_ || export_samples_spheres_ply(p, samples_spheres_path.string());
		const bool skeleton_ok = !p.export_skeleton_selected_ || export_skeleton_mesh_ply(p, skeleton_path.string());
		const std::string samples_status =
			p.export_samples_mesh_selected_ ? (samples_ok ? samples_mesh_path.string() : "FAILED") : "SKIPPED";
		const std::string samples_normal_color_status = p.export_samples_mesh_normal_color_selected_
														 ? (samples_normal_color_ok ? samples_mesh_normal_color_path.string()
																					: "FAILED")
														 : "SKIPPED";
		const std::string spheres_status =
			p.export_samples_spheres_selected_ ? (spheres_ok ? samples_spheres_path.string() : "FAILED") : "SKIPPED";
		const std::string skeleton_status =
			p.export_skeleton_selected_ ? (skeleton_ok ? skeleton_path.string() : "FAILED") : "SKIPPED";

		std::cout << "[UDFExport] samples_mesh=" << samples_status
				  << " samples_mesh_normal_color=" << samples_normal_color_status
				  << " samples_spheres=" << spheres_status << " skeleton=" << skeleton_status << std::endl;
		return any_selected && samples_ok && samples_normal_color_ok && spheres_ok && skeleton_ok;
	}

	void headless_export_skeleton_ply_prepared(PointsParameters& p, const std::string& filename,
											   bool save_face_components = false)
	{
		if (!p.skeleton_ || !p.skeleton_position_ || !p.skeleton_radius_ || !p.spheres_ || !p.spheres_radius_)
			return;

		if (!std::filesystem::path(filename).parent_path().empty())
			std::filesystem::create_directories(std::filesystem::path(filename).parent_path());

		foreach_cell(*p.skeleton_, [&](NMVertex v) -> bool {
			const uint32 v_index = index_of(*p.skeleton_, v);
			if (v_index < nb_cells<PVertex>(*p.spheres_))
				(*p.skeleton_radius_)[v_index] = (*p.spheres_radius_)[v_index];
			return true;
		});

		io::SurfaceExportAttributeSelection<NONMANIFOLD> export_attributes;
		export_attributes.vertex_attributes.push_back(p.skeleton_radius_);
		if (save_face_components)
		{
			prune_fully_non_manifold_triangles(p, "[FaceComponentsPreDelete]");
			if (compute_skeleton_face_components_union_find(p, "[FaceComponentsExport]"))
				colorize_skeleton_face_components(p, "[FaceComponentsExport]");
			if (p.skeleton_face_component_color_)
				export_attributes.face_attributes.push_back(p.skeleton_face_component_color_);
		}
		io::export_PLY(*p.skeleton_, p.skeleton_position_.get(), filename, &export_attributes);
	}

	void headless_export_skeleton_ply_prepared(POINTS& points, const std::string& filename,
											   bool save_face_components = false)
	{
		headless_export_skeleton_ply_prepared(points_parameters_[&points], filename, save_face_components);
	}

	void reset_headless_state()
	{
		if (surface_provider_ && selected_surface_)
			surface_provider_->clear_mesh(*selected_surface_);

		for (auto& [points, p] : points_parameters_)
		{
			(void)points;
			if (points_provider_)
			{
				if (p.points_)
					points_provider_->clear_mesh(*p.points_);
				if (p.samples_mesh_ && p.samples_mesh_ != p.points_)
					points_provider_->clear_mesh(*p.samples_mesh_);
				if (p.spheres_ && p.spheres_ != p.points_ && p.spheres_ != p.samples_mesh_)
					points_provider_->clear_mesh(*p.spheres_);
			}
			if (non_manifold_provider_ && p.skeleton_)
				non_manifold_provider_->clear_mesh(*p.skeleton_);
		}

		points_parameters_.clear();
		selected_surface_ = nullptr;
		surface_bvh_.reset();
		surface_bvh_faces_.clear();
		surface_bvh_vertices_.clear();
		surface_bvh_vertex_positions_.clear();
		surface_vertex_normal_ = nullptr;
		surface_bvh_dirty_ = false;

		selected_points_ = nullptr;
		picked_sphere_ = PVertex();
		timer_connection_.reset();
	}

	void set_selected_surface(SURFACE& s)
	{
		if (!surface_provider_)
			throw std::runtime_error("UDFTraining surface provider is not initialized before set_selected_surface().");
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
		if (!points_provider_)
			throw std::runtime_error("UDFTraining points provider is not initialized before set_selected_points().");
		if (!non_manifold_provider_)
			throw std::runtime_error("UDFTraining non-manifold provider is not initialized before set_selected_points().");
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
			log_error(p, "Neural UDF model file does not exist: ", model_path, '\n');
			return;
		}
		try
		{
			log_basic(p, "Loading Neural UDF model from: ", model_path, '\n');
			p.neural_udf_model_ = torch::jit::load(model_path, device_);
			p.neural_udf_model_.eval();
			p.neural_udf_loaded_ = true;
			p.neural_udf_model_path_ = model_path;
			p.neural_model_type_ = model_type;
			p.udf_input_normalized_ = false;
			p.udf_normalized_source_ = nullptr;
			p.input_mode_ = INPUT_NEURAL_UDF;
			log_basic(p, "Loaded neural UDF model from: ", model_path, '\n');
		}
		catch (const c10::Error& e)
		{
			log_error(p, "Error loading Neural UDF model: ", e.what(), '\n');
			p.neural_udf_loaded_ = false;
		}
	}

	bool try_eval_debug_udf_value(POINTS& points, const Vec3& query_point, Scalar& out_value)
	{
		PointsParameters* params = nullptr;
		auto it = points_parameters_.find(&points);
		if (it != points_parameters_.end())
		{
			params = &it->second;
		}
		else
		{
			for (auto& [input_points, p] : points_parameters_)
			{
				(void)input_points;
				if (p.points_ == &points || p.samples_mesh_ == &points || p.spheres_ == &points)
				{
					params = &p;
					break;
				}
			}
		}
		if (!params)
			return false;
		PointsParameters& p = *params;
		if (!p.neural_udf_loaded_)
			return false;
		std::vector<Scalar> values;
		if (!eval_udf_values(p, std::vector<Vec3>{query_point}, values) || values.empty())
			return false;
		out_value = values[0];
		return true;
	}

	NeuralFieldForward make_neural_field_forward(PointsParameters& p)
	{
		return NeuralFieldForward(&p.neural_udf_model_, p.neural_udf_loaded_, device_);
	}

	bool is_valid_sphere_for_clustering(const PointsParameters& p, uint32 sphere_index) const
	{
		if (!p.spheres_position_ || !p.spheres_radius_ || sphere_index == INVALID_INDEX)
			return false;
		const Vec3& center = (*p.spheres_position_)[sphere_index];
		const Scalar radius = (*p.spheres_radius_)[sphere_index];
		return center.allFinite() && std::isfinite(static_cast<double>(radius)) && radius > Scalar(0);
	}

	template <typename Tag = RaySamplerTag>
	void load_alpha_samples_to_mesh(PointsParameters& p)
	{
		load_alpha_samples_to_mesh_impl(p, Tag{});
	}

	void finalize_sample_mesh_after_sampling(PointsParameters& p)
	{
		log_basic(p, "Building KDTree for sampled points...", '\n');
		build_kdtree(p);
		points_provider_->emit_connectivity_changed(*p.samples_mesh_);
		log_basic(p, "Alpha level set sampling complete. Ready for fitting.", '\n');
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

	bool compute_neural_gradients_gpu(PointsParameters& p, const std::vector<Vec3>& points, std::vector<Vec3>& normals)
	{
		NeuralFieldForward udf = make_neural_field_forward(p);
		return geometry::evaluate_udf_normals(udf, points, static_cast<size_t>(std::max(1, p.batch_size_)),
									  p.neural_projection_workspace_, normals);
	}
	bool project_points_to_alpha_impl(PointsParameters& p, const std::vector<Vec3>& points, const Vec3& bbox_min,
									  const Vec3& bbox_max, std::vector<Vec3>& projected_points,
									  std::vector<Vec3>& normals, std::vector<uint8_t>& keep_mask,
									  geometry::RaySamplerNeural, int max_iters)
	{
		geometry::AlphaProjectionResult result;
		NeuralFieldForward udf = make_neural_field_forward(p);
		const bool ok = geometry::project_points_to_alpha(
			udf, points, geometry::AlphaProjectionParameters{p.alpha_, p.tol_, bbox_min, bbox_max, max_iters},
			p.filter_positive_projection_, p.neural_projection_workspace_, result);
		projected_points = std::move(result.projected_points);
		normals = std::move(result.normals);
		keep_mask = std::move(result.keep_mask);
		return ok;
	}
	bool query_surface_projection_info(const Vec3& sample_pos, const SAttribute<Vec3>* surface_position,
									   AlphaProjectionInfo& info)
	{
		info = AlphaProjectionInfo{};
		if (!selected_surface_ || !surface_position || !surface_bvh_)
			return false;
		uint32 primitive_index;
		if (!geometry::query_surface_closest_point(
				*surface_bvh_, sample_pos, primitive_index, info.closest_point, info.distance))
			return false;
		info.normal = sample_pos - info.closest_point;
		if (info.normal.squaredNorm() < Scalar(1e-12))
		{
			if (primitive_index < surface_bvh_faces_.size())
				info.normal = geometry::normal(*selected_surface_, surface_bvh_faces_[primitive_index], surface_position);
		}
		if (info.normal.squaredNorm() < Scalar(1e-12))
			info.normal = Vec3(0, 0, 1);
		else
			info.normal.normalize();
		info.valid = info.closest_point.allFinite() && info.normal.allFinite() &&
					 std::isfinite(static_cast<double>(info.distance));
		return info.valid;
	}

	bool query_point_cloud_projection_info(PointsParameters& p, const Vec3& sample_pos, AlphaProjectionInfo& info)
	{
		info = AlphaProjectionInfo{};
		if (!p.input_kdtree_ || !p.points_ || !p.position_ || p.input_kdtree_vertices_.empty())
			return false;
		uint32 point_index;
		if (!geometry::query_point_cloud_nearest_point(
				*p.input_kdtree_, sample_pos, point_index, info.distance))
			return false;
		if (point_index >= p.input_kdtree_vertices_.size())
			return false;
		const PVertex nearest_vertex = p.input_kdtree_vertices_[point_index];
		const uint32 nearest_index = index_of(*p.points_, nearest_vertex);
		info.closest_point = (*p.position_)[nearest_index];
		info.normal = sample_pos - info.closest_point;
		if (info.normal.squaredNorm() < Scalar(1e-12) && p.normal_)
			info.normal = (*p.normal_)[nearest_index];
		if (info.normal.squaredNorm() < Scalar(1e-12))
			info.normal = Vec3(0, 0, 1);
		else
			info.normal.normalize();
		info.valid = info.closest_point.allFinite() && info.normal.allFinite() &&
					 std::isfinite(static_cast<double>(info.distance));
		return info.valid;
	}
	bool project_points_to_alpha_impl(PointsParameters& p, const std::vector<Vec3>& points, const Vec3& bbox_min,
									  const Vec3& bbox_max, std::vector<Vec3>& projected_points,
									  std::vector<Vec3>& normals, std::vector<uint8_t>& keep_mask,
									  geometry::RaySamplerSurface, int max_iters)
	{
		auto surface_position = get_attribute<Vec3, SVertex>(*selected_surface_, "position");
		if (!surface_position)
			return false;
		geometry::AlphaProjectionResult result;
		const bool ok = geometry::project_points_to_alpha(
			points, geometry::AlphaProjectionParameters{p.alpha_, p.tol_, bbox_min, bbox_max, max_iters},
			[&](const Vec3& sample_pos, geometry::AlphaProjectionInfo& info) {
				return query_surface_projection_info(sample_pos, surface_position.get(), info);
			},
			result);
		projected_points = std::move(result.projected_points);
		normals = std::move(result.normals);
		keep_mask = std::move(result.keep_mask);
		return ok;
	}
	bool project_points_to_alpha_impl(PointsParameters& p, const std::vector<Vec3>& points, const Vec3& bbox_min,
									  const Vec3& bbox_max, std::vector<Vec3>& projected_points,
									  std::vector<Vec3>& normals, std::vector<uint8_t>& keep_mask,
									  geometry::RaySamplerPointCloud, int max_iters)
	{
		geometry::AlphaProjectionResult result;
		const bool ok = geometry::project_points_to_alpha(
			points, geometry::AlphaProjectionParameters{p.alpha_, p.tol_, bbox_min, bbox_max, max_iters},
			[&](const Vec3& sample_pos, geometry::AlphaProjectionInfo& info) {
				return query_point_cloud_projection_info(p, sample_pos, info);
			},
			result);
		projected_points = std::move(result.projected_points);
		normals = std::move(result.normals);
		keep_mask = std::move(result.keep_mask);
		return ok;
	}
	bool compute_alpha_point_normals_impl(PointsParameters& p, const std::vector<Vec3>& points,
										  std::vector<Vec3>& normals, std::vector<uint8_t>& keep_mask,
										  geometry::RaySamplerNeural)
	{
		keep_mask.assign(points.size(), uint8_t(0));
		if (!compute_neural_gradients_gpu(p, points, normals))
			return false;
		if (normals.size() != points.size())
			return false;
		for (size_t i = 0; i < points.size(); ++i)
		{
			if (!points[i].allFinite() || !normals[i].allFinite())
				continue;
			keep_mask[i] = uint8_t(1);
		}
		return true;
	}

	bool compute_alpha_point_normals_impl(PointsParameters& p, const std::vector<Vec3>& points,
										  std::vector<Vec3>& normals, std::vector<uint8_t>& keep_mask,
										  geometry::RaySamplerSurface)
	{
		auto surface_position = get_attribute<Vec3, SVertex>(*selected_surface_, "position");
		if (!surface_position)
			return false;
		keep_mask.assign(points.size(), uint8_t(0));
		normals.assign(points.size(), Vec3(0, 0, 1));
		for (size_t i = 0; i < points.size(); ++i)
		{
			AlphaProjectionInfo info;
			if (!query_surface_projection_info(points[i], surface_position.get(), info) || !info.valid)
				continue;
			if (!points[i].allFinite() || !info.normal.allFinite())
				continue;
			Vec3 n = info.normal;
			if (n.squaredNorm() < Scalar(1e-12))
				n = Vec3(0, 0, 1);
			else
				n.normalize();
			normals[i] = n;
			keep_mask[i] = uint8_t(1);
		}
		return true;
	}

	bool compute_alpha_point_normals_impl(PointsParameters& p, const std::vector<Vec3>& points,
										  std::vector<Vec3>& normals, std::vector<uint8_t>& keep_mask,
										  geometry::RaySamplerPointCloud)
	{
		keep_mask.assign(points.size(), uint8_t(0));
		normals.assign(points.size(), Vec3(0, 0, 1));
		for (size_t i = 0; i < points.size(); ++i)
		{
			AlphaProjectionInfo info;
			if (!query_point_cloud_projection_info(p, points[i], info) || !info.valid)
				continue;
			if (!points[i].allFinite() || !info.normal.allFinite())
				continue;
			Vec3 n = info.normal;
			if (n.squaredNorm() < Scalar(1e-12))
				n = Vec3(0, 0, 1);
			else
				n.normalize();
			normals[i] = n;
			keep_mask[i] = uint8_t(1);
		}
		return true;
	}

	void materialize_samples_to_mesh(PointsParameters& p, const std::vector<Vec3>& positions,
									 const std::vector<Vec3>& normals)
	{
		if (p.samples_mesh_)
			points_provider_->clear_mesh(*p.samples_mesh_);
		for (size_t i = 0; i < positions.size(); ++i)
		{
			PVertex v = add_vertex(*p.samples_mesh_);
			uint32 v_idx = index_of(*p.samples_mesh_, v);
			(*p.samples_position_)[v_idx] = positions[i];
			const Vec3& n = normals[i];
			(*p.samples_normal_)[v_idx] = n;
			if (p.samples_color_)
				(*p.samples_color_)[v_idx] = Vec4(0.0, 0.0, 0.0, 1.0);
			if (p.samples_normal_color_)
				(*p.samples_normal_color_)[v_idx] =
					Vec4((n.x() + 1.0) * 0.5, (n.y() + 1.0) * 0.5, (n.z() + 1.0) * 0.5, 1.0);
		}
	}

	template <typename Traits, typename Tag>
	void load_alpha_samples_to_mesh_bridson_impl(PointsParameters& p, Traits& traits, const torch::Device& device,
											 const char* mode_label, Tag tag)
	{
		if (p.sample_radius_ <= Scalar(0))
		{
			log_error(p, "Sample radius must be positive for alpha level set sampling.", '\n');
			return;
		}
		if (!traits.is_ready())
		{
			log_error(p, "Sampling traits are not ready. Cannot sample alpha level set.", '\n');
			return;
		}
		RaySamplerParams ray_params = make_ray_params(p);
		if (!p.ray_sampler_)
			p.ray_sampler_ = std::make_unique<RaySampler>(ray_params, device);
		else
		{
			p.ray_sampler_->set_params(ray_params);
			p.ray_sampler_->set_device(device);
		}
		auto [bbox_min, bbox_max] = compute_sampling_bbox(p);
		geometry::AlphaSamplingParameters parameters;
		parameters.sample_radius = p.sample_radius_;
		parameters.seed = p.seed_;
		parameters.ray_sampler_batch_size = static_cast<size_t>(std::max(1, p.ray_sampler_batch_size_));
		parameters.bbox_min = bbox_min;
		parameters.bbox_max = bbox_max;
		geometry::AlphaSamplingResult result = geometry::sample_alpha_level_set(
			traits, *p.ray_sampler_, parameters,
			[&](const std::vector<Vec3>& points, std::vector<Vec3>& normals, std::vector<uint8_t>& keep_mask) {
				return compute_alpha_point_normals_impl(p, points, normals, keep_mask, tag);
			},
			[&](const std::vector<Vec3>& points, std::vector<Vec3>& projected_points, std::vector<Vec3>& normals,
				std::vector<uint8_t>& keep_mask) {
				return project_points_to_alpha_impl(p, points, bbox_min, bbox_max, projected_points, normals, keep_mask, tag, 1);
			});
		p.samples_spatial_grid_ = std::move(result.spatial_grid);
		if (!result.success)
			log_error(p, "Alpha level set sampling failed for ", mode_label, ".", '\n');
		if (result.positions.empty())
		{
			if (p.samples_mesh_)
				points_provider_->clear_mesh(*p.samples_mesh_);
			points_provider_->emit_connectivity_changed(*p.samples_mesh_);
			if (result.success)
				log_error(p, "Failed to generate any alpha level set samples.", '\n');
			return;
		}
		log_basic(p, "[AlphaSampling:", mode_label, "] accepted=", result.positions.size(), '\n');
		materialize_samples_to_mesh(p, result.positions, result.normals);
		finalize_sample_mesh_after_sampling(p);
	}

	void load_alpha_samples_to_mesh_impl(PointsParameters& p, geometry::RaySamplerNeural)
	{
		if (!p.neural_udf_loaded_)
		{
			log_error(p, "Neural UDF model not loaded. Cannot sample alpha level set.", '\n');
			return;
		}

		NeuralFieldForward udf = make_neural_field_forward(p);
		auto traits = RaySamplerConfig::make(udf);
		const bool prev_filter_positive_projection = p.filter_positive_projection_;
		if (p.neural_model_type_ == NEURAL_MODEL_MF)
			p.filter_positive_projection_ = true;
		load_alpha_samples_to_mesh_bridson_impl(p, traits, udf.device(), "neural_udf", geometry::RaySamplerNeural{});
		p.filter_positive_projection_ = prev_filter_positive_projection;
	}

	void load_alpha_samples_to_mesh_impl(PointsParameters& p, geometry::RaySamplerSurface)
	{
		build_surface_bvh();
		if (!surface_bvh_)
		{
			log_error(p, "Surface BVH not available. Cannot sample alpha level set.", '\n');
			return;
		}

		auto s_pos = get_attribute<Vec3, SVertex>(*selected_surface_, "position");
		if (!s_pos)
		{
			log_error(p, "Surface position attribute not available. Cannot sample alpha level set.", '\n');
			return;
		}

		auto traits = RaySamplerConfig::make(*selected_surface_, s_pos.get(), surface_bvh_.get(), &surface_bvh_faces_);
		torch::Device cpu_device(torch::kCPU);
		load_alpha_samples_to_mesh_bridson_impl(p, traits, cpu_device, "surface_mesh", geometry::RaySamplerSurface{});
	}

	void load_alpha_samples_to_mesh_impl(PointsParameters& p, geometry::RaySamplerPointCloud)
	{
		if (!p.input_kdtree_)
		{
			log_error(p, "Input point cloud KDTree not available. Cannot sample alpha level set.", '\n');
			return;
		}
		if (!p.normal_ || !p.knn_)
		{
			log_error(p, "Input point cloud normals/KNN not available. Cannot sample alpha level set.", '\n');
			return;
		}

		auto traits = RaySamplerConfig::make(*p.points_, p.position_.get(), p.normal_.get(), p.knn_.get(),
											 p.input_kdtree_, &p.input_kdtree_vertices_);
		torch::Device cpu_device(torch::kCPU);
		load_alpha_samples_to_mesh_bridson_impl(p, traits, cpu_device, "point_cloud",
												geometry::RaySamplerPointCloud{});
	}

	void pre_process_sampling_points_kdtree(PointsParameters& p, std::vector<Vec3>& points)
	{
		const Scalar tol = Scalar(1e-2);
		const Scalar max_dist = p.alpha_ + tol;
		log_basic(p, "Before pre-processing, ", points.size(), " points.", '\n');
		if (!p.input_kdtree_)
			return;
		auto new_end = std::remove_if(points.begin(), points.end(), [&](const Vec3& pos) {
			std::pair<uint32, Scalar> knn_res;
			return !p.input_kdtree_->find_nn(pos, &knn_res, max_dist);
		});
		points.erase(new_end, points.end());
		log_basic(p, "After pre-processing, ", points.size(), " points.", '\n');
	}

	void pre_process_sampling_points_bvh(PointsParameters& p, std::vector<Vec3>& points)
	{
		const Scalar tol = Scalar(1e-2);
		const Scalar max_dist = p.alpha_ + tol;
		log_basic(p, "Before pre-processing, ", points.size(), " points.", '\n');
		//Todo: the bvh shuld not be built here, to be moved
		build_surface_bvh();
		if (!surface_bvh_)
			return;

		auto new_end = std::remove_if(points.begin(), points.end(), [&](const Vec3& pos) {
			std::pair<uint32, Vec3> cp;
			return !surface_bvh_->closest_point(pos, &cp, max_dist);
		});
		points.erase(new_end, points.end());
		log_basic(p, "After pre-processing, ", points.size(), " points.", '\n');
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
						n = pos - cp.second;
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
			log_error(p, "Stop spheres update before filtering sampled points.", '\n');
			return;
		}

		const uint32 count = nb_cells<PVertex>(*p.samples_mesh_);
		if (count == 0)
		{
			log_basic(p, "No sampled points to filter.", '\n');
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
		if (after == before)
		{
			log_basic(p, "Sampling filtering applied: no points removed (", before, " -> ", after, ").", '\n');
			return;
		}

		p.fitting_data_computed_ = false;
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
			log_basic(p, "Sampling filtering applied: ", before, " -> 0 points.", '\n');
			return;
		}

		for (const Vec3& pt : filtered_points)
		{
			PVertex v = add_vertex(*p.samples_mesh_);
			uint32 v_idx = index_of(*p.samples_mesh_, v);
			(*p.samples_position_)[v_idx] = pt;
			if (p.samples_color_)
				(*p.samples_color_)[v_idx] = Vec4(0.0, 0.0, 0.0, 1.0);
		}

		recompute_samples_normals_from_current_input(p);
		build_kdtree(p);
		points_provider_->emit_connectivity_changed(*p.samples_mesh_);
		points_provider_->emit_attribute_changed(*p.samples_mesh_, p.samples_position_.get());
		points_provider_->emit_attribute_changed(*p.samples_mesh_, p.samples_normal_.get());
		points_provider_->emit_attribute_changed(*p.samples_mesh_, p.samples_normal_color_.get());
		if (p.samples_color_)
			points_provider_->emit_attribute_changed(*p.samples_mesh_, p.samples_color_.get());
		log_basic(p, "Sampling filtering applied: ", before, " -> ", after, " points.", '\n');
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
			log_verbose(p, "Batch UDF evaluation successful. Sample results:", '\n');
			for (size_t i = 0; i < 5; ++i)
			{
				log_verbose(p, "Point: ", test_points[i].transpose(), " UDF: ", result.values[i],
						   " Grad: ", result.gradients[i].transpose(), '\n');
			}
		}
		else
		{
			log_error(p, "Batch UDF evaluation failed.", '\n');
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
		skeleton_render_ = static_cast<SurfaceRender<NONMANIFOLD>*>(
			app_.module("SurfaceRender (" + std::string{mesh_traits<NONMANIFOLD>::name} + ")"));

		timer_connection_ = boost::synapse::connect<App::timer_tick>(&app_, [this]() {
			if (selected_points_)
			{
				PointsParameters& p = points_parameters_[selected_points_];
				if (p.running_ && p.preview_render_during_sphere_update_)
				{
					update_render_data(p, false, true, true);
					request_linked_views_update();
				}
				else if (p.pending_full_refresh_after_stop_)
				{
					update_render_data(p, false, true, true);
					request_linked_views_update();
					p.pending_full_refresh_after_stop_ = false;
				}
			}
		});
		// Initialize PyTorch device
		const OutputVerbosity startup_output_verbosity =
			selected_points_ ? points_parameters_[selected_points_].output_verbosity_ : OUTPUT_NORMAL;
		if (torch::cuda::is_available())
		{
			log_basic(startup_output_verbosity, "CUDA is available! Using GPU device 0.", '\n');
			device_ = torch::Device(torch::kCUDA, 0);
		}
		else
		{
			log_basic(startup_output_verbosity, "CUDA is not available! Using the CPU.", '\n');
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
		p.spheres_error_ = get_or_add_attribute<Scalar, PVertex>(*p.spheres_, "error");
		p.spheres_error_not_normalized_ = get_or_add_attribute<Scalar, PVertex>(*p.spheres_, "error_not_normalized");
		p.spheres_skeleton_vertex_ = get_or_add_attribute<NMVertex, PVertex>(*p.spheres_, "skeleton_vertex");

		// Init Skeleton Mesh
		std::string skel_name = points_provider_->mesh_name(m) + "_skeleton";
		if (!p.skeleton_)
			p.skeleton_ = non_manifold_provider_->has_mesh(skel_name) ? non_manifold_provider_->mesh(skel_name)
															  : non_manifold_provider_->add_mesh(skel_name);
		p.skeleton_position_ = get_or_add_attribute<Vec3, NMVertex>(*p.skeleton_, "position");
		p.skeleton_radius_ = get_or_add_attribute<Scalar, NMVertex>(*p.skeleton_, "radius");
		p.skeleton_source_sphere_ = get_or_add_attribute<PVertex, NMVertex>(*p.skeleton_, "source_sphere");
		p.incident_tets_ = get_or_add_attribute<std::set<std::size_t>, NMFace>(*p.skeleton_, "incident_tets");
		p.edge_degree_ = get_or_add_attribute<uint32, NMEdge>(*p.skeleton_, "degree");
		p.skeleton_face_color_ = get_or_add_attribute<Vec3, NMFace>(*p.skeleton_, "color");
		p.skeleton_face_component_id_ = get_or_add_attribute<uint32, NMFace>(*p.skeleton_, "face_component_id");
		p.skeleton_face_component_color_ = get_or_add_attribute<Vec3, NMFace>(*p.skeleton_, "face_component_color");
		p.skeleton_edge_color_ = get_or_add_attribute<Vec3, NMEdge>(*p.skeleton_, "color");
		p.skeleton_edge_non_manifold_color_ =
			get_or_add_attribute<Vec3, NMEdge>(*p.skeleton_, "non_manifold_edge_color");
		if (p.sample_grid_cell_size_ > Scalar(0))
			p.samples_spatial_grid_ = std::make_unique<SpatialGrid>(p.sample_grid_cell_size_);
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
			log_error(p, "Error: No sampled points found. Please sample points first.", '\n');
			return;
		}

		log_basic(p, "Building KDTree...", '\n');
		build_kdtree(p);
		if (p.recompute_sample_normals_after_sampling_)
		{
			log_basic(p, "Recomputing Normals (PCA)...", '\n');
			recompute_samples_normals_pca(p);
		}
		log_basic(p, "Computing Initial Medial Axis...", '\n');
		compute_initial_medial_axis(p);
		build_kdtree(p);
		log_basic(p, "Computing KNN and Area...", '\n');
		compute_samples_area(p); // Compute KNN and Area for samples
		log_basic(p, "Computing Quadrics...", '\n');
		compute_quadrics(p);

		log_basic(p, "Fitting Data Computed.", '\n');

		p.fitting_data_computed_ = true;
	}

	void build_kdtree(PointsParameters& p)
	{
		if (p.samples_kdtree_)
			delete p.samples_kdtree_;
		if (p.samples_ma_kdtree_)
			delete p.samples_ma_kdtree_;

		std::vector<Vec3> points;
		const uint32 sample_count = nb_cells<PVertex>(*p.samples_mesh_);
		p.samples_kdtree_vertices_.clear();
		p.samples_ma_kdtree_vertices_.clear();
		points.reserve(sample_count);
		p.samples_kdtree_vertices_.reserve(sample_count);
		p.samples_ma_kdtree_vertices_.reserve(sample_count);
		std::vector<Vec3> points_ma;
		points_ma.reserve(sample_count);

		foreach_cell(*p.samples_mesh_, [&](PVertex v) {
			uint32 idx = index_of(*p.samples_mesh_, v);
			points.push_back((*p.samples_position_)[idx]);
			p.samples_kdtree_vertices_.push_back(v);
			if (p.samples_ma_position_)
			{
				const Vec3& ma_pos = (*p.samples_ma_position_)[idx];
				const bool ma_finite = ma_pos.allFinite();
				const bool ma_radius_ok = (!p.samples_ma_radius_) || ((*p.samples_ma_radius_)[idx] > Scalar(0));
				if (ma_finite && ma_radius_ok)
				{
					points_ma.push_back(ma_pos);
					p.samples_ma_kdtree_vertices_.push_back(v);
				}
			}
			return true;
		});

		p.samples_kdtree_ = points.empty() ? nullptr : new acc::KDTree<3, uint32>(points);
		p.samples_ma_kdtree_ = points_ma.empty() ? nullptr : new acc::KDTree<3, uint32>(points_ma);
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
			centroid += position[index_of(mesh, kdtree_vertices[idx])];
		centroid /= Scalar(indices.size());

		Eigen::Matrix3d covariance = Eigen::Matrix3d::Zero();
		for (uint32 idx : indices)
		{
			const Vertex v = kdtree_vertices[idx];
			const Vec3 diff = position[index_of(mesh, v)] - centroid;
			const Eigen::Vector3d point(diff[0], diff[1], diff[2]);
			covariance += point * point.transpose();
		}
		covariance /= Scalar(indices.size());
		const Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> solver(covariance);
		const Eigen::Vector3d normal = solver.eigenvectors().col(0);
		return Vec3(normal[0], normal[1], normal[2]).normalized();
	}

	void compute_input_normals(PointsParameters& p)
	{
		if (!p.input_kdtree_ || !p.points_ || !p.position_ || !p.normal_ || !p.knn_)
			return;
		parallel_foreach_cell(*p.points_, [&](PVertex v) -> bool {
			const uint32 v_idx = index_of(*p.points_, v);
			const Vec3& pt = (*p.position_)[v_idx];
			std::vector<std::pair<uint32, Scalar>> knn_res;
			p.input_kdtree_->find_nns(pt, p.knn_k_ + 1, &knn_res);
			std::vector<uint32> indices;
			(*p.knn_)[v_idx].clear();
			for (const auto& res : knn_res)
			{
				indices.push_back(res.first);
				(*p.knn_)[v_idx].push_back(p.input_kdtree_vertices_[res.first]);
			}
			(*p.normal_)[v_idx] = compute_pca_normal(*p.points_, *p.position_, indices, p.input_kdtree_vertices_);
			return true;
		});
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
			p.samples_kdtree_->find_nns(pt, p.knn_k_ + 10, &knn_res);

			(*p.samples_knn_)[v_idx].clear();
			Scalar sum_dist = 0.0;
			const Scalar eps = Scalar(1e-12);
			const Scalar band = p.alpha_;
			Vec3 n = (*p.samples_normal_)[v_idx];
			const Scalar n2 = n.squaredNorm();
			if (n2 > eps)
				n /= std::sqrt(n2);
			else
				n = Vec3(0, 0, 1);
			int kept = 0;
			for (auto& res : knn_res)
			{
				if (p.samples_kdtree_vertices_[res.first] != v)
				{
					PVertex nb = p.samples_kdtree_vertices_[res.first];
					uint32 nb_idx = index_of(*p.samples_mesh_, nb);
					const Vec3& q = (*p.samples_position_)[nb_idx];
					const Scalar dn = (q - pt).dot(n);
					if (std::abs(dn) > band)
						continue;
					(*p.samples_knn_)[v_idx].push_back(nb);
					sum_dist += res.second;
					++kept;
					if (kept >= p.knn_k_ + 1)
						break;
				}
			}
			if (kept == 0) // fall back
			{
				for (auto& res : knn_res)
				{
					if (p.samples_kdtree_vertices_[res.first] != v)
					{
						PVertex nb = p.samples_kdtree_vertices_[res.first];
						(*p.samples_knn_)[v_idx].push_back(nb);
						sum_dist += res.second;
					}
				}
			}
			// Normals are already computed/oriented in sample_points
			(*p.samples_area_)[v_idx] = (sum_dist * sum_dist) / (2.0 * p.knn_k_); // Rough area estimate
			return true;
		});
	}

	void compute_quadrics(PointsParameters& p)
	{
		if (!p.samples_mesh_)
			return;
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

	bool recompute_sample_normal_pca_for_vertex(PointsParameters& p, PVertex v)
	{
		if (!p.samples_mesh_ || !p.samples_kdtree_ || !p.samples_position_ || !p.samples_normal_ || !v.is_valid())
			return false;

		const int k = std::max(3, p.knn_k_);
		const Scalar eps = Scalar(1e-12);
		const uint32 v_idx = index_of(*p.samples_mesh_, v);
		if (v_idx == INVALID_INDEX)
			return false;
		const Vec3& center = (*p.samples_position_)[v_idx];

		std::vector<std::pair<uint32, Scalar>> knn_res;
		p.samples_kdtree_->find_nns(center, k + 1, &knn_res);

		std::vector<Vec3> neighbors;
		neighbors.reserve(k);
		for (const auto& res : knn_res)
		{
			PVertex nb = p.samples_kdtree_vertices_[res.first];
			if (nb == v)
				continue;
			uint32 nb_idx = index_of(*p.samples_mesh_, nb);
			neighbors.push_back((*p.samples_position_)[nb_idx]);
			if (static_cast<int>(neighbors.size()) >= k)
				break;
		}
		if (neighbors.size() < 3)
			return false;

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
			return false;

		Eigen::Matrix<Scalar, 3, 1> ev = solver.eigenvectors().col(0);
		Vec3 n(ev(0), ev(1), ev(2));
		Scalar n2 = n.squaredNorm();
		if (n2 < eps)
			return false;

		Vec3 n0 = (*p.samples_normal_)[v_idx];
		if (n0.squaredNorm() > eps && n.dot(n0) < Scalar(0))
			n = -n;
		n.normalize();
		(*p.samples_normal_)[v_idx] = n;
		return true;
	}

	void recompute_samples_normals_pca(PointsParameters& p)
	{
		if (!p.samples_mesh_ || !p.samples_kdtree_ || !p.samples_position_ || !p.samples_normal_)
			return;
		parallel_foreach_cell(*p.samples_mesh_, [&](PVertex v) -> bool {
			recompute_sample_normal_pca_for_vertex(p, v);
			return true;
		});
		refresh_sample_normals_color(p);
		points_provider_->emit_attribute_changed(*p.samples_mesh_, p.samples_normal_.get());
		points_provider_->emit_attribute_changed(*p.samples_mesh_, p.samples_normal_color_.get());
	}

	void recompute_samples_normals_pca_for_vertices(PointsParameters& p, const std::vector<PVertex>& vertices)
	{
		if (!p.samples_mesh_ || !p.samples_kdtree_ || !p.samples_position_ || !p.samples_normal_ || vertices.empty())
			return;
		std::unordered_set<uint32> visited_ids;
		std::vector<PVertex> updated_vertices;
		updated_vertices.reserve(vertices.size());
		for (PVertex v : vertices)
		{
			if (!v.is_valid())
				continue;
			const uint32 v_idx = index_of(*p.samples_mesh_, v);
			if (v_idx == INVALID_INDEX || !visited_ids.insert(v_idx).second)
				continue;
			if (recompute_sample_normal_pca_for_vertex(p, v))
				updated_vertices.push_back(v);
		}
		if (updated_vertices.empty())
			return;
		for (PVertex v : updated_vertices)
		{
			const uint32 v_idx = index_of(*p.samples_mesh_, v);
			const Vec3& n = (*p.samples_normal_)[v_idx];
			(*p.samples_normal_color_)[v_idx] =
				Vec4((n.x() + 1.0) * 0.5, (n.y() + 1.0) * 0.5, (n.z() + 1.0) * 0.5, 1.0);
		}
		points_provider_->emit_attribute_changed(*p.samples_mesh_, p.samples_normal_.get());
		points_provider_->emit_attribute_changed(*p.samples_mesh_, p.samples_normal_color_.get());
	}

	Vec3 random_sample_in_sphere(const Vec3& p, const Scalar radius,
							 std::uniform_real_distribution<Scalar>& uni, std::mt19937& rng)
	{
		const Scalar u = uni(rng);
		const Scalar v = uni(rng);
		const Scalar w = uni(rng);
		const Scalar r = radius * std::cbrt(u);
		const Scalar phi = v * 2.0 * M_PI;
		const Scalar theta = std::acos(1.0 - 2.0 * w);
		return p + Vec3(r * std::sin(theta) * std::cos(phi), r * std::sin(theta) * std::sin(phi),
					r * std::cos(theta));
	}

	void refresh_sample_normals_color(PointsParameters& p)
	{
		parallel_foreach_cell(*p.samples_mesh_, [&](PVertex v) -> bool {
			const uint32 v_idx = index_of(*p.samples_mesh_, v);
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
		if (p.samples_mesh_)
			points_provider_->clear_mesh(*p.samples_mesh_);
	}

	// --- Initial Medial Axis ---

	void compute_initial_medial_axis(PointsParameters& p)
	{
		if (!p.samples_mesh_ || !p.samples_kdtree_ || !p.samples_position_ || !p.samples_normal_ ||
			!p.samples_ma_position_ || !p.samples_ma_radius_ || !p.samples_ma_secondary_vertex_)
			return;

		const Scalar fallback_radius = p.alpha_;
		const Scalar min_norm = Scalar(1e-12);
		const bool neural_input = (p.input_mode_ == INPUT_NEURAL_UDF && p.neural_udf_loaded_);
		const bool mf_model = neural_input && (p.neural_model_type_ == NEURAL_MODEL_MF);
		const bool udf_model = neural_input && (p.neural_model_type_ == NEURAL_MODEL_UDF);
		// All input modes use shrinking-ball medial-axis initialization.
		const Scalar initial_radius = std::max<Scalar>(fallback_radius * Scalar(10), Scalar(0));

		auto run_shrinking_ball_for_vertex = [&](PVertex v) -> bool {
			if (!v.is_valid())
				return false;
			const uint32 v_idx = index_of(*p.samples_mesh_, v);
			if (v_idx == INVALID_INDEX)
				return false;
			const Vec3& pt = (*p.samples_position_)[v_idx];
			Vec3 n = (*p.samples_normal_)[v_idx];
			const bool normal_finite = n.allFinite();
			const Scalar n2 = normal_finite ? n.squaredNorm() : Scalar(0);
			const bool can_run_shrinking_ball = normal_finite && (n2 > min_norm);
			if (can_run_shrinking_ball)
				n.normalize();
			else
				n = Vec3(0, 0, 1);

			Vec3 c = pt - n * fallback_radius;
			Scalar r = fallback_radius;
			PVertex secondary;

			if (can_run_shrinking_ball)
			{
				auto [c1, r1, q1] = geometry::shrinking_ball_center<PVertex>(
					pt, n, p.samples_kdtree_, p.samples_kdtree_vertices_, initial_radius);

				if (c1.allFinite() && std::isfinite(static_cast<double>(r1)) && r1 > Scalar(0))
				{
					c = c1;
					r = r1;
					secondary = q1;
				}
			}

			if (!c.allFinite())
				c = pt - n * fallback_radius;
			if (!std::isfinite(static_cast<double>(r)) || r <= Scalar(0))
				r = fallback_radius;

			if (!secondary.is_valid() && p.samples_kdtree_ && !p.samples_kdtree_vertices_.empty())
			{
				const Vec3 q = c - n * r;
				std::pair<uint32, Scalar> knn_res;
				p.samples_kdtree_->find_nn(q, &knn_res);
				if (knn_res.first < p.samples_kdtree_vertices_.size())
					secondary = p.samples_kdtree_vertices_[knn_res.first];
			}

			(*p.samples_ma_position_)[v_idx] = c;
			const Scalar expected_radius = p.alpha_;
			(*p.samples_ma_radius_)[v_idx] = (r > expected_radius) ? r : expected_radius;
			(*p.samples_ma_secondary_vertex_)[v_idx] = secondary;
			return true;
		};

		auto run_shrinking_ball_for_all = [&]() {
			parallel_foreach_cell(*p.samples_mesh_, [&](PVertex v) -> bool {
				run_shrinking_ball_for_vertex(v);
				return true;
			});
		};

		run_shrinking_ball_for_all();

		// Optional post-process on shrinking-ball centers:
		// - UDF model: retry/delete only if udf(center) > 1.2 * alpha.
		// - Surface mesh / point cloud: retry if distance(center) > alpha; still > alpha -> delete sample.
		// - MF model: retry if sdf(center) > 0 or udf(center) > alpha; after flip, delete if sdf(center) > 0
		//   or udf(center) > alpha.
		if (p.ma_flip_prune_enabled_)
		{
			uint32 flipped_normals = 0;
			uint32 flip_triggered_points = 0;
			uint32 deleted_samples = 0;
			std::vector<uint32> flipped_vertex_indices;
			std::function<uint32(const std::vector<PVertex>&)> prune_vertices_by_current_centers;
			const Scalar udf_center_prune_threshold = udf_model ? (Scalar(1.2) * p.alpha_) : p.alpha_;
			const Scalar mf_center_prune_threshold =
				p.alpha_ * std::max(Scalar(1), std::min(Scalar(5), p.ma_flip_prune_alpha_factor_));
			const bool supports_center_retry_prune =
				mf_model || p.input_mode_ == INPUT_SURFACE_MESH || p.input_mode_ == INPUT_POINT_CLOUD ||
				(p.input_mode_ == INPUT_NEURAL_UDF && p.neural_udf_loaded_);
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

			if (supports_center_retry_prune)
			{
				auto eval_center_scores = [&](const std::vector<Vec3>& query_points, std::vector<Scalar>& out_values,
										  std::vector<Scalar>& out_sdf) -> bool {
					out_sdf.clear();
					if (mf_model)
						return eval_mf_values_sdf(query_points, out_values, out_sdf);
					return eval_topology_score_values(p, query_points, out_values);
				};

				prune_vertices_by_current_centers = [&](const std::vector<PVertex>& candidate_vertices) -> uint32 {
					std::vector<PVertex> active_vertices;
					std::vector<Vec3> active_centers;
					active_vertices.reserve(candidate_vertices.size());
					active_centers.reserve(candidate_vertices.size());
					for (PVertex v : candidate_vertices)
					{
						const uint32 vid = index_of(*p.samples_mesh_, v);
						if (vid == INVALID_INDEX)
							continue;
						active_vertices.push_back(v);
						active_centers.push_back((*p.samples_ma_position_)[vid]);
					}
					if (active_vertices.empty())
						return 0;

					std::vector<Scalar> current_score;
					std::vector<Scalar> current_sdf;
					const bool eval_ok = eval_center_scores(active_centers, current_score, current_sdf);
					if (!eval_ok || current_score.size() != active_centers.size())
						return 0;

					std::unordered_set<uint32> deleted_ids;
					deleted_ids.reserve(active_vertices.size());
					uint32 removed_count = 0;
					for (size_t i = 0; i < active_vertices.size(); ++i)
					{
						bool should_delete = false;
						if (mf_model)
						{
							if ((i < current_sdf.size() && current_sdf[i] > Scalar(0)) ||
								current_score[i] > mf_center_prune_threshold)
								should_delete = true;
						}
						else
						{
							if (current_score[i] > udf_center_prune_threshold)
								should_delete = true;
						}
						if (!should_delete)
							continue;

						const uint32 vid = index_of(*p.samples_mesh_, active_vertices[i]);
						if (vid == INVALID_INDEX || !deleted_ids.insert(vid).second)
							continue;
						remove_vertex(*p.samples_mesh_, active_vertices[i]);
						++removed_count;
					}
					return removed_count;
				};

				std::vector<PVertex> vertices;
				std::vector<Vec3> centers;
				vertices.reserve(nb_cells<PVertex>(*p.samples_mesh_));
				centers.reserve(nb_cells<PVertex>(*p.samples_mesh_));
				foreach_cell(*p.samples_mesh_, [&](PVertex v) -> bool {
					const uint32 vid = index_of(*p.samples_mesh_, v);
					if (vid == INVALID_INDEX)
						return true;
					vertices.push_back(v);
					centers.push_back((*p.samples_ma_position_)[vid]);
					return true;
				});

				std::vector<Scalar> score_values;
				std::vector<Scalar> sdf_values;
				const bool eval_ok = eval_center_scores(centers, score_values, sdf_values);

				if (eval_ok && score_values.size() == centers.size())
				{
					std::vector<PVertex> need_retry;
					need_retry.reserve(vertices.size() / 8 + 1);
					for (size_t i = 0; i < vertices.size(); ++i)
					{
						const bool by_udf =
							mf_model ? (score_values[i] > mf_center_prune_threshold)
									 : (score_values[i] > udf_center_prune_threshold);
						const bool by_mf_sdf = mf_model && i < sdf_values.size() && (sdf_values[i] > Scalar(0));
						const bool trigger_retry = mf_model ? (by_mf_sdf || by_udf) : by_udf;
						if (trigger_retry)
							need_retry.push_back(vertices[i]);
					}
					flip_triggered_points += static_cast<uint32>(need_retry.size());
					flipped_vertex_indices.reserve(need_retry.size());

					for (PVertex v : need_retry)
					{
						const uint32 vid = index_of(*p.samples_mesh_, v);
						if (vid == INVALID_INDEX)
							continue;
						Vec3 n = (*p.samples_normal_)[vid];
						if (n.squaredNorm() > min_norm)
						{
							n.normalize();
							(*p.samples_normal_)[vid] = -n;
							flipped_vertex_indices.push_back(vid);
							++flipped_normals;
						}
						run_shrinking_ball_for_vertex(v);
					}

					deleted_samples += prune_vertices_by_current_centers(need_retry);
				}
			}

			if (deleted_samples > 0)
			{
				// Sample connectivity changed; refresh the geometry needed to stabilize MA first.
				build_kdtree(p);
				if (nb_cells<PVertex>(*p.samples_mesh_) > 0 && !flipped_vertex_indices.empty())
				{
					std::vector<PVertex> surviving_flipped_vertices;
					surviving_flipped_vertices.reserve(flipped_vertex_indices.size());
					for (uint32 vid : flipped_vertex_indices)
					{
						PVertex v = of_index<PVertex>(*p.samples_mesh_, vid);
						if (v.is_valid())
							surviving_flipped_vertices.push_back(v);
					}

					recompute_samples_normals_pca_for_vertices(p, surviving_flipped_vertices);
					for (PVertex v : surviving_flipped_vertices)
						run_shrinking_ball_for_vertex(v);
					if (prune_vertices_by_current_centers)
						deleted_samples += prune_vertices_by_current_centers(surviving_flipped_vertices);
				}
				points_provider_->emit_connectivity_changed(*p.samples_mesh_);
				log_basic(p, "[MAFlipPrune] flip_triggered_points=", flip_triggered_points, " flipped_normals=",
						  flipped_normals, " deleted_samples=", deleted_samples, " remaining_samples=",
						  nb_cells<PVertex>(*p.samples_mesh_), '\n');
			}
			else if (flip_triggered_points > 0 || flipped_normals > 0)
			{
				log_basic(p, "[MAFlipPrune] flip_triggered_points=", flip_triggered_points, " flipped_normals=",
						  flipped_normals, " deleted_samples=0", '\n');
			}
		}

		// MA positions changed; keep MA-KDTree in sync for MF topology scoring.
		build_kdtree(p);
	}

	// MF post-process: search along opposite normal direction for minimal mf-abs(sdf)

	void init_spheres(PointsParameters& p)
	{
		points_provider_->clear_mesh(*p.spheres_);
		init_spheres_from_samples(p);
		if (p.samples_sphere_)
			p.samples_sphere_->fill(PVertex());
		compute_clusters_local(p);

		if (!p.running_)
		{
			set_post_init_sphere_render_state(p);
		}
	}

	void init_spheres_from_samples(PointsParameters& p)
	{
		if (!p.samples_mesh_ || !p.samples_position_ || !p.samples_knn_ || !p.samples_ma_position_ ||
			!p.samples_ma_radius_ || !p.samples_ma_secondary_vertex_)
			return;
		constexpr uint32 min_cover_points = 10;
		constexpr uint32 max_nb_spheres = 100000;
		std::vector<PVertex> sorted_vertices;
		uint32 max_sample_index = 0;
		foreach_cell(*p.samples_mesh_, [&](PVertex v) {
			sorted_vertices.push_back(v);
			const uint32 idx = index_of(*p.samples_mesh_, v);
			if (idx != INVALID_INDEX)
				max_sample_index = std::max(max_sample_index, idx);
			return true;
		});
		std::sort(sorted_vertices.begin(), sorted_vertices.end(), [&](PVertex a, PVertex b) {
			// sort candidate spheres by decreasing radius
			uint32 idx_a = index_of(*p.samples_mesh_, a);
			uint32 idx_b = index_of(*p.samples_mesh_, b);
			const Scalar ra = (idx_a != INVALID_INDEX) ? (*p.samples_ma_radius_)[idx_a] : Scalar(-1);
			const Scalar rb = (idx_b != INVALID_INDEX) ? (*p.samples_ma_radius_)[idx_b] : Scalar(-1);
			const Scalar safe_ra = std::isfinite(static_cast<double>(ra)) ? ra : Scalar(-1);
			const Scalar safe_rb = std::isfinite(static_cast<double>(rb)) ? rb : Scalar(-1);
			return safe_ra > safe_rb;
		});
		auto covered = get_or_add_attribute<bool, PVertex>(*p.samples_mesh_, "__covered");
		covered->fill(false);
		const std::size_t candidate_marks_size =
			sorted_vertices.empty() ? std::size_t(0) : (static_cast<std::size_t>(max_sample_index) + 1);
		std::vector<uint32> candidate_marks(candidate_marks_size, 0);
		uint32 candidate_mark_token = 1;

		p.nb_spheres_ = 0;
		uint32 skipped_invalid_ma_seeds = 0;

		for (PVertex v : sorted_vertices)
		{
			uint32 v_index = index_of(*p.samples_mesh_, v);
			if (v_index == INVALID_INDEX)
				continue;

			if (p.nb_spheres_ >= max_nb_spheres)
				break;

			if ((*covered)[v_index])
				continue;

			const Vec3& vp = (*p.samples_ma_position_)[v_index];
			Scalar vr = (*p.samples_ma_radius_)[v_index];
			if (!vp.allFinite() || !std::isfinite(static_cast<double>(vr)) || vr <= Scalar(0))
			{
				++skipped_invalid_ma_seeds;
				continue;
			}
			const Scalar dilation_radius = std::max<Scalar>(vr + Scalar(p.init_dilation_constant_), Scalar(0));
			const Scalar dilation_radius_sq = dilation_radius * dilation_radius;
			if (candidate_mark_token == std::numeric_limits<uint32>::max())
			{
				std::fill(candidate_marks.begin(), candidate_marks.end(), 0u);
				candidate_mark_token = 1;
			}
			const uint32 current_mark = candidate_mark_token++;
			std::vector<PVertex> candidate_cover;
			candidate_cover.reserve(128);

			auto flood_cover = [&](PVertex seed) {
				if (!seed.is_valid())
					return;
				uint32 seed_idx = index_of(*p.samples_mesh_, seed);
				if (seed_idx == INVALID_INDEX || seed_idx >= candidate_marks.size() || (*covered)[seed_idx] ||
					candidate_marks[seed_idx] == current_mark)
					return;

				std::vector<PVertex> stack;
				stack.reserve(128);
				candidate_marks[seed_idx] = current_mark;
				candidate_cover.push_back(seed);
				stack.push_back(seed);
				while (!stack.empty())
				{
					PVertex w = stack.back();
					stack.pop_back();
					uint32 w_idx = index_of(*p.samples_mesh_, w);
					if (w_idx == INVALID_INDEX)
						continue;

					for (PVertex u : (*p.samples_knn_)[w_idx])
					{
						uint32 u_idx = index_of(*p.samples_mesh_, u);
						if (u_idx == INVALID_INDEX || u_idx >= candidate_marks.size())
							continue;
						if (!(*covered)[u_idx] && candidate_marks[u_idx] != current_mark &&
							((*p.samples_position_)[u_idx] - vp).squaredNorm() < dilation_radius_sq)
						{
							candidate_marks[u_idx] = current_mark;
							candidate_cover.push_back(u);
							stack.push_back(u);
						}
					}
				}
			};

			flood_cover(v);

			// Secondary seed can reach another lobe; skip if already covered to avoid redundant traversal.
			PVertex secondary = (*p.samples_ma_secondary_vertex_)[v_index];
			if (secondary.is_valid())
			{
				flood_cover(secondary);
			}

			for (PVertex covered_vertex : candidate_cover)
			{
				uint32 covered_index = index_of(*p.samples_mesh_, covered_vertex);
				if (covered_index == INVALID_INDEX)
					continue;
				(*covered)[covered_index] = true;
			}

			const bool keep_sphere =
				(candidate_cover.size() >= min_cover_points) || (p.nb_spheres_ == 0 && !candidate_cover.empty());
			if (!keep_sphere)
				continue;

			PVertex sphere = add_vertex(*p.spheres_);
			p.nb_spheres_++;
			uint32 sphere_index = index_of(*p.spheres_, sphere);

			(*p.spheres_position_)[sphere_index] = vp;
			(*p.spheres_radius_)[sphere_index] = vr;
			(*p.spheres_cluster_color_)[sphere_index] =
				Vec4(0.5 + 0.5 * (rand() % 256) / 256.0, 0.5 + 0.5 * (rand() % 256) / 256.0,
					 0.5 + 0.5 * (rand() % 256) / 256.0, 1.0);
		}
		remove_attribute<PVertex>(*p.samples_mesh_, covered);
		if (skipped_invalid_ma_seeds > 0)
			log_basic(p, "[InitSpheres] skipped_invalid_ma_seeds=", skipped_invalid_ma_seeds, '\n');

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
	};

	bool get_sphere_fit_data(PointsParameters& p, SphereFitData& data)
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

		if (!data.mesh || !data.position || !data.area || !data.quadric || !data.line_quadric || !data.sphere)
			return false;
		return true;
	}

	void compute_clusters_local(PointsParameters& p)
	{
		if (p.nb_spheres_ == 0)
		{
			if (p.samples_sphere_)
				p.samples_sphere_->fill(PVertex());
			return;
		}
		SphereFitData data;
		if (!get_sphere_fit_data(p, data))
			return;
		std::atomic<uint64> invalid_sphere_candidates(0);
		std::atomic<uint64> nonfinite_distance_candidates(0);
		std::atomic<uint64> unassigned_samples(0);
		parallel_foreach_cell(*p.spheres_, [&](PVertex v) -> bool {
			uint32 v_index = index_of(*p.spheres_, v);
			(*p.spheres_cluster_)[v_index].clear();
			(*p.spheres_cluster_area_)[v_index] = 0.0;
			return true;
		});

		parallel_foreach_cell(*data.mesh, [&](PVertex v) -> bool {
			uint32 v_index = index_of(*data.mesh, v);

			Scalar a = (*data.area)[v_index];
			PVertex cluster_sphere = (*data.sphere)[v_index];
			const uint32 cs_index = cluster_sphere.is_valid() ? index_of(*p.spheres_, cluster_sphere) : INVALID_INDEX;

			Scalar min_distance = std::numeric_limits<Scalar>::max();
			PVertex closest_sphere;
			uint32 closest_sphere_index = INVALID_INDEX;

			auto evaluate_candidate = [&](PVertex pv) {
				uint32 pv_index = index_of(*p.spheres_, pv);
				if (!is_valid_sphere_for_clustering(p, pv_index))
				{
					++invalid_sphere_candidates;
					return;
				}
				const Vec3& center = (*p.spheres_position_)[pv_index];
				Scalar radius = (*p.spheres_radius_)[pv_index];

				Scalar dist_sqem = (*data.quadric)[v_index].eval(Vec4(center.x(), center.y(), center.z(), radius));
				Scalar dist_other = (*data.line_quadric)[v_index].eval(center);
				Scalar dist = dist_sqem + Scalar(p.sqem_update_lambda_line_plane_) * dist_other;
				if (!std::isfinite(static_cast<double>(dist)))
				{
					++nonfinite_distance_candidates;
					return;
				}
				if (dist < min_distance)
				{
					min_distance = dist;
					closest_sphere = pv;
					closest_sphere_index = pv_index;
				}
			};

			if (cs_index != INVALID_INDEX && is_valid_sphere_for_clustering(p, cs_index))
			{
				const std::set<PVertex>& neighbor_spheres = (*p.spheres_neighbor_clusters_)[cs_index];
				bool owner_evaluated = false;
				for (PVertex pv : neighbor_spheres)
				{
					if (!owner_evaluated && cluster_sphere < pv)
					{
						evaluate_candidate(cluster_sphere);
						owner_evaluated = true;
					}
					if (pv == cluster_sphere)
						owner_evaluated = true;
					evaluate_candidate(pv);
				}
				if (!owner_evaluated)
					evaluate_candidate(cluster_sphere);
			}
			else
			{
				// Invalid or uninitialized ownership requires a global search.
				foreach_cell(*p.spheres_, [&](PVertex pv) -> bool {
					evaluate_candidate(pv);
					return true;
				});
			}

			(*data.sphere)[v_index] = closest_sphere;
			if (!closest_sphere.is_valid() || closest_sphere_index == INVALID_INDEX)
			{
				++unassigned_samples;
				return true;
			}

			std::lock_guard<std::mutex> lock(spheres_mutex_[closest_sphere_index % spheres_mutex_.size()]);
			(*p.spheres_cluster_)[closest_sphere_index].push_back(v);
			(*p.spheres_cluster_area_)[closest_sphere_index] += a;

			return true;
		});
		if (invalid_sphere_candidates.load() > 0 || nonfinite_distance_candidates.load() > 0 || unassigned_samples.load() > 0)
		{
			log_error(p, "[ClusterLocalGuard] invalid_sphere_candidates=", invalid_sphere_candidates.load(),
					  " nonfinite_distance_candidates=", nonfinite_distance_candidates.load(),
					  " unassigned_samples=", unassigned_samples.load(), '\n');
		}
	}

	void compute_clusters(PointsParameters& p)
	{
		compute_clusters_local(p);
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
				Scalar dist = dist_sqem + Scalar(p.sqem_update_lambda_line_plane_) * dist_other;
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
	bool eval_udf_values(PointsParameters& p, const std::vector<Vec3>& points, std::vector<Scalar>& out_values)
	{
		return geometry::evaluate_udf_values(make_neural_field_forward(p), points,
											std::max<size_t>(1, static_cast<size_t>(p.batch_size_)), out_values);
	}
	bool ensure_topology_score_backend(PointsParameters& p, const char* log_prefix = "[TopologyScore]")
	{
		const bool use_mf_ma_distance =
			(p.input_mode_ == INPUT_NEURAL_UDF && p.neural_model_type_ == NEURAL_MODEL_MF);
		if (use_mf_ma_distance)
		{
			if (!p.samples_mesh_ || !p.samples_ma_position_)
			{
				log_error(p, log_prefix, " missing sample/ma data for MF score evaluation.", '\n');
				return false;
			}
			if (!p.samples_ma_kdtree_ || p.samples_ma_kdtree_vertices_.empty())
				build_kdtree(p);
			if (!p.samples_ma_kdtree_ || p.samples_ma_kdtree_vertices_.empty())
			{
				log_error(p, log_prefix, " MA KDTree unavailable (no valid ma_position / ma_radius). ",
						  "Run fitting data computation first.", '\n');
				return false;
			}
			return true;
		}

		if (p.input_mode_ == INPUT_NEURAL_UDF)
		{
			if (!p.neural_udf_loaded_)
			{
				log_error(p, log_prefix, " requires a loaded neural UDF model.", '\n');
				return false;
			}
			return true;
		}

		if (p.input_mode_ == INPUT_SURFACE_MESH)
		{
			build_surface_bvh();
			if (!surface_bvh_)
			{
				log_error(p, log_prefix, " requires a valid surface BVH for mesh distance evaluation.", '\n');
				return false;
			}
			return true;
		}

		if (p.input_mode_ == INPUT_POINT_CLOUD)
		{
			if ((!p.input_kdtree_ || p.input_kdtree_vertices_.empty()) && p.points_ && p.position_)
				rebuild_input_kdtree(p);
			if (!p.input_kdtree_ || p.input_kdtree_vertices_.empty() || !p.points_ || !p.position_)
			{
				log_error(p, log_prefix, " requires a valid input KDTree for point-cloud distance evaluation.", '\n');
				return false;
			}
			return true;
		}

		log_error(p, log_prefix, " no supported topology-score backend for current input mode.", '\n');
		return false;
	}

	bool eval_topology_score_values(PointsParameters& p, const std::vector<Vec3>& points,
									std::vector<Scalar>& out_values)
	{
		out_values.clear();
		out_values.resize(points.size(), Scalar(0));
		if (points.empty())
			return true;
		if (!ensure_topology_score_backend(p))
			return false;

		const bool use_mf_ma_distance =
			(p.input_mode_ == INPUT_NEURAL_UDF && p.neural_model_type_ == NEURAL_MODEL_MF);
		if (!use_mf_ma_distance)
		{
			if (p.input_mode_ == INPUT_NEURAL_UDF)
				return eval_udf_values(p, points, out_values);
			if (p.input_mode_ == INPUT_SURFACE_MESH)
			{
				for (size_t i = 0; i < points.size(); ++i)
				{
					uint32 primitive_index;
					Vec3 closest_point;
					Scalar distance;
					if (!geometry::query_surface_closest_point(
							*surface_bvh_, points[i], primitive_index, closest_point, distance))
						continue;
					(void)primitive_index;
					out_values[i] = distance;
				}
				return true;
			}
			if (p.input_mode_ == INPUT_POINT_CLOUD)
			{
				for (size_t i = 0; i < points.size(); ++i)
				{
					uint32 point_index;
					Scalar distance;
					if (!geometry::query_point_cloud_nearest_point(
							*p.input_kdtree_, points[i], point_index, distance))
						continue;
					(void)distance;
					const uint32 nn_idx = point_index;
					if (nn_idx >= p.input_kdtree_vertices_.size())
						continue;
					const PVertex nn = p.input_kdtree_vertices_[nn_idx];
					const uint32 vid = index_of(*p.points_, nn);
					if (vid == INVALID_INDEX)
						continue;
					out_values[i] = (points[i] - (*p.position_)[vid]).norm();
				}
				return true;
			}
			log_error(p, "[TopologyScore] unsupported backend during score evaluation.", '\n');
			return false;
		}

		for (size_t i = 0; i < points.size(); ++i)
		{
			std::pair<uint32, Scalar> knn_res;
			if (!p.samples_ma_kdtree_->find_nn(points[i], &knn_res))
				continue;
			const uint32 nn_idx = knn_res.first;
			if (nn_idx >= p.samples_ma_kdtree_vertices_.size())
				continue;
			const PVertex nearest_sample = p.samples_ma_kdtree_vertices_[nn_idx];
			const uint32 sid = index_of(*p.samples_mesh_, nearest_sample);
			if (sid == INVALID_INDEX)
				continue;
			const Vec3& ma_pos = (*p.samples_ma_position_)[sid];
			out_values[i] = (points[i] - ma_pos).norm();
		}
		return true;
	}

	bool try_get_nearest_sample_ma_radius(PointsParameters& p, const Vec3& query, Scalar& out_radius)
	{
		if (!p.samples_mesh_ || !p.samples_ma_radius_)
			return false;
		if (!p.samples_kdtree_ || p.samples_kdtree_vertices_.empty())
			build_kdtree(p);
		if (!p.samples_kdtree_ || p.samples_kdtree_vertices_.empty())
			return false;

		std::pair<uint32, Scalar> knn_res;
		p.samples_kdtree_->find_nn(query, &knn_res);
		if (knn_res.first >= p.samples_kdtree_vertices_.size())
			return false;

		const uint32 sample_index = index_of(*p.samples_mesh_, p.samples_kdtree_vertices_[knn_res.first]);
		if (sample_index == INVALID_INDEX)
			return false;

		const Scalar ma_radius = (*p.samples_ma_radius_)[sample_index];
		if (!std::isfinite(ma_radius) || ma_radius <= Scalar(0))
			return false;

		out_radius = ma_radius;
		return true;
	}

	void update_sphere_line_quadric_distance_fix_current_radius(PointsParameters& p, PVertex sphere, Scalar fixed_radius)
	{
		if (!std::isfinite(fixed_radius) || fixed_radius <= Scalar(0))
			return;

		SphereFitData data;
		if (!get_sphere_fit_data(p, data))
			return;
		uint32 sphere_index = index_of(*p.spheres_, sphere);

		const std::vector<PVertex>& cluster = (*p.spheres_cluster_)[sphere_index];
		if (cluster.empty())
			return;
		Vec3 c = (*p.spheres_position_)[sphere_index];
		Spherical_Quadric q;
		Line_Quadric lq;
		for (PVertex v : cluster)
		{
			uint32 v_index = index_of(*data.mesh, v);
			Scalar weight = value<Scalar>(*data.mesh, data.area, v);
			if (weight <= Scalar(0))
			{
				log_verbose(p, "Warning: sample with zero volume weight in sphere ", sphere_index, '\n');
				continue;
			}
			q += (*data.quadric)[v_index] * weight;
			lq += (*data.line_quadric)[v_index] * weight;
		}

		Mat4 Ql = lq.get_quadric().matrix();
		Mat3 Al = Ql.block<3, 3>(0, 0);
		Vec3 bl = -Ql.block<3, 1>(0, 3);

		Mat3 As = q._A.block<3, 3>(0, 0);
		Vec3 bs = q._b.head<3>();
		Vec3 Asr = q._A.block<3, 1>(0, 3);

		const Scalar update_lambda = Scalar(p.sqem_update_lambda_line_plane_);
		Mat3 A = As + update_lambda * Al;
		Vec3 b = (bs + update_lambda * bl) - Asr * fixed_radius;

		c = A.ldlt().solve(b);
		if (!c.allFinite())
			return;

		(*p.spheres_position_)[sphere_index] = c;
		(*p.spheres_radius_)[sphere_index] = fixed_radius;
	}

	void update_sphere_line_quadric_distance_free_radius(PointsParameters& p, PVertex sphere)
	{
		SphereFitData data;
		if (!get_sphere_fit_data(p, data))
			return;
		uint32 sphere_index = index_of(*p.spheres_, sphere);

		const std::vector<PVertex>& cluster = (*p.spheres_cluster_)[sphere_index];
		const Scalar radius = (*p.spheres_radius_)[sphere_index];
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

		const Scalar update_lambda = Scalar(p.sqem_update_lambda_line_plane_);
		Mat4 A = q._A + update_lambda * Al_ext;
		Vec4 b = q._b + update_lambda * bl_ext;
		Vec4 s = A.completeOrthogonalDecomposition().solve(b);
		if (!s.allFinite())
			return;
		Scalar nearest_ma_radius = Scalar(0);
		const bool radius_too_large =
			try_get_nearest_sample_ma_radius(p, s.head<3>(), nearest_ma_radius) &&
			(s[3] > nearest_ma_radius * Scalar(1.5));
		const bool radius_non_positive = (s[3] <= Scalar(0));
		if (radius_too_large || radius_non_positive)
		{
			update_sphere_line_quadric_distance_fix_current_radius(p, sphere, radius);
			return;
		}

		(*p.spheres_position_)[sphere_index] = s.head<3>();
		(*p.spheres_radius_)[sphere_index] = s[3];
	}
	bool should_refresh_local_connectivity(const PointsParameters& p)
	{
		return (p.iteration_count_ % 10) == 0;
	}


	void update_spheres(PointsParameters& p, HeadlessOptimizationStats* stats = nullptr)
	{
		auto cluster_start = std::chrono::high_resolution_clock::now();
		if (should_refresh_local_connectivity(p))
			compute_sphere_neighbors(p);
		compute_clusters(p);
		auto cluster_end = std::chrono::high_resolution_clock::now();
		if (stats)
			stats->cluster_total_ms_ +=
				std::chrono::duration<float64, std::milli>(cluster_end - cluster_start).count();

		auto update_start = std::chrono::high_resolution_clock::now();
		parallel_foreach_cell(*p.spheres_, [&](PVertex v) -> bool {
			update_sphere_line_quadric_distance_free_radius(p, v);
			return true;
		});
		auto update_end = std::chrono::high_resolution_clock::now();
		if (stats)
			stats->sphere_update_total_ms_ +=
				std::chrono::duration<float64, std::milli>(update_end - update_start).count();

		auto error_start = std::chrono::high_resolution_clock::now();
		compute_spheres_error(p);
		auto error_end = std::chrono::high_resolution_clock::now();
		if (stats)
			stats->error_total_ms_ +=
				std::chrono::duration<float64, std::milli>(error_end - error_start).count();
	}

	void remove_sphere(PointsParameters& p, PVertex v)
	{
		uint32 v_index = INVALID_INDEX;
		if (!is_live_sphere_vertex(p, v, v_index))
			return;
		SphereFitData data;
		if (!get_sphere_fit_data(p, data))
			return;
		const std::vector<PVertex> removed_cluster = (*p.spheres_cluster_)[v_index];
		std::set<PVertex> removed_neighbors;
		if (p.spheres_neighbor_clusters_)
			removed_neighbors = (*p.spheres_neighbor_clusters_)[v_index];
		for (PVertex s : removed_cluster)
		{
			uint32 s_idx = index_of(*data.mesh, s);
			(*data.sphere)[s_idx] = PVertex();
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

		if (p.spheres_skeleton_vertex_)
		{
			const NMVertex linked_skeleton_vertex = (*p.spheres_skeleton_vertex_)[v_index];
			if (linked_skeleton_vertex.is_valid() && p.skeleton_ && p.skeleton_source_sphere_)
				value<PVertex>(*p.skeleton_, p.skeleton_source_sphere_, linked_skeleton_vertex) = PVertex();
			(*p.spheres_skeleton_vertex_)[v_index] = NMVertex();
		}

		remove_vertex(*p.spheres_, v);
		p.nb_spheres_--;
		redistribute_removed_sphere_cluster(p, removed_cluster, removed_neighbors);
	}

	bool remove_skeleton_vertex_and_linked_sphere(PointsParameters& p, const NMVertex& v)
	{
		if (!p.skeleton_ || !v.is_valid())
			return false;
		const uint32 skeleton_vertex_id = index_of(*p.skeleton_, v);
		if (skeleton_vertex_id == INVALID_INDEX)
			return false;

		PVertex linked_sphere;
		if (p.skeleton_source_sphere_)
			linked_sphere = value<PVertex>(*p.skeleton_, p.skeleton_source_sphere_, v);

		if (linked_sphere.is_valid() && p.spheres_)
		{
			const uint32 sphere_index = index_of(*p.spheres_, linked_sphere);
			if (sphere_index != INVALID_INDEX && p.spheres_skeleton_vertex_)
				(*p.spheres_skeleton_vertex_)[sphere_index] = NMVertex();
			remove_sphere(p, linked_sphere);
		}

		if (p.skeleton_source_sphere_)
			value<PVertex>(*p.skeleton_, p.skeleton_source_sphere_, v) = PVertex();
		remove_vertex(*p.skeleton_, v);
		return true;
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

	void compute_sphere_neighbors(PointsParameters& p)
	{
		SphereFitData data;
		if (!get_sphere_fit_data(p, data))
			return;
		auto knn_attr = data.knn;
		if (!data.mesh || !data.sphere)
			return;
		if (!knn_attr)
			return;
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

	}

	void compute_skeleton(PointsParameters& p)
	{
		compute_sphere_neighbors(p);
		SphereFitData data;
		if (!get_sphere_fit_data(p, data))
			return;
		auto knn_attr = data.knn;
		if (!data.mesh || !data.sphere || !knn_attr)
			return;
		clear(*p.skeleton_);
		std::map<NMFaceKey, NMFace> skeleton_faces_map;
		auto get_face_key = [](uint32 i1, uint32 i2, uint32 i3) -> NMFaceKey {
			std::array<uint32, 3> key = {i1, i2, i3};
			std::sort(key.begin(), key.end());
			return key;
		};
		std::vector<std::array<uint32, 4>> raw_tets;
		if (p.spheres_skeleton_vertex_)
		{
			foreach_cell(*p.spheres_, [&](PVertex pv) -> bool {
				const uint32 sphere_index = index_of(*p.spheres_, pv);
				if (sphere_index != INVALID_INDEX)
					(*p.spheres_skeleton_vertex_)[sphere_index] = NMVertex();
				return true;
			});
		}
		foreach_cell(*p.spheres_, [&](PVertex pv) -> bool {
			uint32 pv_index = index_of(*p.spheres_, pv);
			NMVertex nmv = add_vertex(*p.skeleton_);
			(*p.skeleton_position_)[index_of(*p.skeleton_, nmv)] = (*p.spheres_position_)[pv_index];
			if (p.spheres_skeleton_vertex_)
				(*p.spheres_skeleton_vertex_)[pv_index] = nmv;
			if (p.skeleton_source_sphere_)
				value<PVertex>(*p.skeleton_, p.skeleton_source_sphere_, nmv) = pv;
			return true;
		});

		std::unordered_map<std::pair<uint32, uint32>, NMEdge, edge_hash, edge_equal> edge_indices;

		foreach_cell(*p.spheres_, [&](PVertex pv) -> bool {
			uint32 pv_index = index_of(*p.spheres_, pv);
			NMVertex nmv1 = p.spheres_skeleton_vertex_ ? (*p.spheres_skeleton_vertex_)[pv_index] : NMVertex();
			const std::set<PVertex>& neighbors = (*p.spheres_neighbor_clusters_)[pv_index];
			for (PVertex neighbor : neighbors)
			{
				uint32 n_index = index_of(*p.spheres_, neighbor);
				NMVertex nmv2 = p.spheres_skeleton_vertex_ ? (*p.spheres_skeleton_vertex_)[n_index] : NMVertex();
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
			NMVertex nmv1 = p.spheres_skeleton_vertex_ ? (*p.spheres_skeleton_vertex_)[idx1] : NMVertex();
			const std::set<PVertex>& n_pv = (*p.spheres_neighbor_clusters_)[idx1];
			for (const PVertex& ne1 : n_pv)
			{
				uint32 idx2 = index_of(*p.spheres_, ne1);
				if (idx1 >= idx2)
					continue;
				NMVertex nmv2 = p.spheres_skeleton_vertex_ ? (*p.spheres_skeleton_vertex_)[idx2] : NMVertex();
				const std::set<PVertex>& ne_ne1 = (*p.spheres_neighbor_clusters_)[idx2];
				for (const PVertex& ne2 : ne_ne1)
				{
					if (n_pv.find(ne2) == n_pv.end())
						continue;
					uint32 idx3 = index_of(*p.spheres_, ne2);
					if (idx2 >= idx3)
						continue;

					NMVertex nmv3 = p.spheres_skeleton_vertex_ ? (*p.spheres_skeleton_vertex_)[idx3] : NMVertex();

					std::vector<NMEdge> edges;
					edges.reserve(3);
					edges.push_back(edge_indices[{index_of(*p.skeleton_, nmv1), index_of(*p.skeleton_, nmv2)}]);
					edges.push_back(edge_indices[{index_of(*p.skeleton_, nmv2), index_of(*p.skeleton_, nmv3)}]);
					edges.push_back(edge_indices[{index_of(*p.skeleton_, nmv1), index_of(*p.skeleton_, nmv3)}]);
					NMFace new_face = add_face(*p.skeleton_, edges);

					skeleton_faces_map[get_face_key(idx1, idx2, idx3)] = new_face;

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
				auto it = skeleton_faces_map.find(keys[i]);
				if (it != skeleton_faces_map.end())
				{
					NMFace f = it->second;
					new_tet.faces[i] = f;

					value<std::set<size_t>>(*p.skeleton_, p.incident_tets_, f).insert(tet_index);

					value<Vec3>(*p.skeleton_, p.skeleton_face_color_, f) = Vec3(0.8, 0.5, 0.5);
				}
			}
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
		non_manifold_provider_->emit_attribute_changed(*p.skeleton_, p.skeleton_face_color_.get());
		non_manifold_provider_->emit_attribute_changed(*p.skeleton_, p.skeleton_edge_color_.get());
		non_manifold_provider_->emit_attribute_changed(*p.skeleton_, p.skeleton_edge_non_manifold_color_.get());
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
		if (!center_sphere.is_valid())
			return;

		std::vector<PVertex> center_spheres = {center_sphere};
		recompute_clusters_local_neighborhoods(p, center_spheres);
	}

	void recompute_clusters_local_neighborhoods(PointsParameters& p, const std::vector<PVertex>& center_spheres)
	{
		if (center_spheres.empty() || p.nb_spheres_ == 0)
			return;

		SphereFitData data;
		if (!get_sphere_fit_data(p, data))
			return;
		const uint32 nb_samples = nb_cells<PVertex>(*data.mesh);
		if (nb_samples == 0)
			return;

		std::vector<PVertex> candidate_spheres;
		std::unordered_set<uint32> candidate_indices;
		candidate_spheres.reserve(center_spheres.size() * 8);

		auto try_add_candidate = [&](PVertex sphere) {
			if (!sphere.is_valid())
				return;
			const uint32 s_index = index_of(*p.spheres_, sphere);
			if (s_index == INVALID_INDEX)
				return;
			if (candidate_indices.insert(s_index).second)
				candidate_spheres.push_back(sphere);
		};

		for (PVertex center_sphere : center_spheres)
		{
			if (!center_sphere.is_valid())
				continue;
			const uint32 center_index = index_of(*p.spheres_, center_sphere);
			if (center_index == INVALID_INDEX)
				continue;

			try_add_candidate(center_sphere);
			for (PVertex neighbor : (*p.spheres_neighbor_clusters_)[center_index])
				try_add_candidate(neighbor);
		}
		if (candidate_spheres.empty())
			return;

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
			Scalar dist = dist_sqem + Scalar(p.sqem_update_lambda_line_plane_) * dist_other;
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

	void redistribute_removed_sphere_cluster(PointsParameters& p, const std::vector<PVertex>& removed_cluster,
											 const std::set<PVertex>& removed_neighbors)
	{
		SphereFitData data;
		if (!get_sphere_fit_data(p, data))
			return;
		if (removed_cluster.empty() || p.nb_spheres_ == 0)
		{
			return;
		}

		std::vector<PVertex> candidate_spheres;
		candidate_spheres.reserve(removed_neighbors.size() + 4);
		std::unordered_set<uint32> candidate_indices;
		auto try_add_candidate = [&](PVertex sphere) {
			if (!sphere.is_valid())
				return;
			const uint32 s_index = index_of(*p.spheres_, sphere);
			if (s_index == INVALID_INDEX)
				return;
			if (candidate_indices.insert(s_index).second)
				candidate_spheres.push_back(sphere);
		};

		for (PVertex neighbor : removed_neighbors)
			try_add_candidate(neighbor);

		if (candidate_spheres.empty())
		{
			foreach_cell(*p.spheres_, [&](PVertex sphere) -> bool {
				try_add_candidate(sphere);
				return true;
			});
		}

		if (candidate_spheres.empty())
		{
			return;
		}
		auto eval_distance = [&](uint32 sample_index, uint32 sphere_index) -> Scalar {
			const Vec3& center = (*p.spheres_position_)[sphere_index];
			const Scalar radius = (*p.spheres_radius_)[sphere_index];
			const Scalar dist_sqem =
				(*data.quadric)[sample_index].eval(Vec4(center.x(), center.y(), center.z(), radius));
			const Scalar dist_other = (*data.line_quadric)[sample_index].eval(center);
			return dist_sqem + Scalar(p.sqem_update_lambda_line_plane_) * dist_other;
		};

		for (PVertex sample : removed_cluster)
		{
			const uint32 sample_index = index_of(*data.mesh, sample);
			if (sample_index == INVALID_INDEX)
				continue;

			Scalar min_distance = std::numeric_limits<Scalar>::max();
			PVertex closest_sphere;
			uint32 closest_sphere_index = INVALID_INDEX;
			for (PVertex sphere : candidate_spheres)
			{
				const uint32 sphere_index = index_of(*p.spheres_, sphere);
				if (sphere_index == INVALID_INDEX)
					continue;
				const Scalar dist = eval_distance(sample_index, sphere_index);
				if (dist < min_distance)
				{
					min_distance = dist;
					closest_sphere = sphere;
					closest_sphere_index = sphere_index;
				}
			}

			if (!closest_sphere.is_valid() || closest_sphere_index == INVALID_INDEX)
				continue;

			(*data.sphere)[sample_index] = closest_sphere;
			(*p.spheres_cluster_)[closest_sphere_index].push_back(sample);
			(*p.spheres_cluster_area_)[closest_sphere_index] += (*data.area)[sample_index];
		}
	}

protected:
	PVertex split_sphere(PointsParameters& p, PVertex sphere)
	{
		if (!sphere.is_valid())
			return PVertex();
		SphereFitData data;
		if (!get_sphere_fit_data(p, data))
			return PVertex();
		uint32 s_index = index_of(*p.spheres_, sphere);

		// find the point in the cluster with the max error
		const std::vector<PVertex>& cluster = (*p.spheres_cluster_)[s_index];
		if (cluster.empty())
			return PVertex();

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
			return PVertex();

		uint32 max_err_v_idx = index_of(*data.mesh, max_err_v);

		PVertex new_sphere = add_vertex(*p.spheres_);
		uint32 new_s_index = index_of(*p.spheres_, new_sphere);

		(*p.spheres_position_)[new_s_index] = (*p.samples_ma_position_)[max_err_v_idx];
		(*p.spheres_radius_)[new_s_index] = (*p.samples_ma_radius_)[max_err_v_idx];
		(*p.spheres_parent_)[new_s_index] = sphere;
		(*p.spheres_cluster_color_)[new_s_index] =
			Vec4(0.5 + 0.5 * (rand() % 256) / 256.0, 0.5 + 0.5 * (rand() % 256) / 256.0,
				 0.5 + 0.5 * (rand() % 256) / 256.0, 1.0);

		p.nb_spheres_++;
		inherit_sphere_neighbors(p, sphere, new_sphere);
		recompute_clusters_local_neighborhood(p, new_sphere);
		return new_sphere;
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

	void update_edge_degree_for_edge_ids(PointsParameters& p, const std::unordered_set<uint32>& edge_ids)
	{
		if (!p.skeleton_ || !p.edge_degree_)
			return;
		for (uint32 edge_id : edge_ids)
		{
			if (edge_id == INVALID_INDEX)
				continue;
			const NMEdge e = of_index<NMEdge>(*p.skeleton_, edge_id);
			if (!e.is_valid())
				continue;
			if (index_of(*p.skeleton_, e) != edge_id)
				continue;

			const uint32 deg = static_cast<uint32>(incident_faces(*p.skeleton_, e).size());
			(*p.edge_degree_)[edge_id] = deg;
			if (deg == 2)
				value<Vec3>(*p.skeleton_, p.skeleton_edge_color_, e) = Vec3(0.0, 1.0, 0.0);
			if (p.skeleton_edge_non_manifold_color_)
			{
				value<Vec3>(*p.skeleton_, p.skeleton_edge_non_manifold_color_, e) =
					(deg > 2) ? Vec3(1.0, 0.2, 0.1) : Vec3(0.15, 0.15, 0.15);
			}
		}
	}

	struct FaceComponentUnionFind
	{
		std::vector<uint32> parent_;
		std::vector<uint32> rank_;

		FaceComponentUnionFind() = default;

		explicit FaceComponentUnionFind(size_t n)
		{
			reset(n);
		}

		void reset(size_t n)
		{
			parent_.resize(n);
			rank_.assign(n, 0);
			for (uint32 i = 0; i < static_cast<uint32>(n); ++i)
				parent_[i] = i;
		}

		uint32 find(uint32 x)
		{
			uint32 root = x;
			while (parent_[root] != root)
				root = parent_[root];
			while (parent_[x] != x)
			{
				const uint32 next = parent_[x];
				parent_[x] = root;
				x = next;
			}
			return root;
		}

		void unite(uint32 a, uint32 b)
		{
			uint32 ra = find(a);
			uint32 rb = find(b);
			if (ra == rb)
				return;
			if (rank_[ra] < rank_[rb])
				std::swap(ra, rb);
			parent_[rb] = ra;
			if (rank_[ra] == rank_[rb])
				++rank_[ra];
		}
	};

	Vec3 hsv_to_rgb(Scalar h, Scalar s, Scalar v) const
	{
		const Scalar hh = h - std::floor(h);
		if (!(s > Scalar(0)))
			return Vec3(v, v, v);
		const Scalar scaled = hh * Scalar(6);
		const int sector = static_cast<int>(std::floor(scaled));
		const Scalar f = scaled - Scalar(sector);
		const Scalar p = v * (Scalar(1) - s);
		const Scalar q = v * (Scalar(1) - s * f);
		const Scalar t = v * (Scalar(1) - s * (Scalar(1) - f));

		switch (sector % 6)
		{
		case 0:
			return Vec3(v, t, p);
		case 1:
			return Vec3(q, v, p);
		case 2:
			return Vec3(p, v, t);
		case 3:
			return Vec3(p, q, v);
		case 4:
			return Vec3(t, p, v);
		default:
			return Vec3(v, p, q);
		}
	}

	Vec3 component_palette_color(uint32 component_id) const
	{
		static constexpr double golden_ratio_conjugate = 0.6180339887498949;
		const Scalar hue = Scalar(std::fmod(0.17 + golden_ratio_conjugate * static_cast<double>(component_id), 1.0));
		const Scalar saturation = (component_id % 3 == 0) ? Scalar(0.82) : (component_id % 3 == 1) ? Scalar(0.72)
																									  : Scalar(0.90);
		const Scalar value = (component_id % 4 == 0) ? Scalar(0.95) : (component_id % 4 == 1) ? Scalar(0.88)
																							 : (component_id % 4 == 2) ? Scalar(0.80)
																													 : Scalar(0.92);
		return hsv_to_rgb(hue, saturation, value);
	}

	bool prune_fully_non_manifold_triangles(
		PointsParameters& p, const char* log_prefix = "[FaceComponentsPreDelete]")
	{
		if (!p.skeleton_)
			return false;

		const bool collect_verbose_stats = is_verbose_logging_enabled(p);
		uint32 total_removed_faces = 0;
		uint32 total_removed_edges = 0;
		uint32 total_removed_vertices = 0;
		uint32 total_passes = 0;
		bool removed_any_faces = false;

		for (;;)
		{
			std::vector<uint32> face_ids_to_remove;
			foreach_cell(*p.skeleton_, [&](NMFace f) -> bool {
				if (!f.is_valid())
					return true;
				const uint32 face_id = index_of(*p.skeleton_, f);
				if (face_id == INVALID_INDEX)
					return true;

				const std::vector<NMEdge> face_edges = incident_edges(*p.skeleton_, f);
				if (face_edges.size() != 3)
					return true;

				for (const NMEdge& e : face_edges)
				{
					if (!e.is_valid() || incident_faces(*p.skeleton_, e).size() <= 2)
						return true;
				}

				face_ids_to_remove.push_back(face_id);
				return true;
			});

			if (face_ids_to_remove.empty())
				break;

			if (collect_verbose_stats)
				++total_passes;
			uint32 removed_faces_this_pass = 0;
			uint32 removed_edges_this_pass = 0;
			uint32 removed_vertices_this_pass = 0;
			for (uint32 face_id : face_ids_to_remove)
			{
				const NMFace f = of_index<NMFace>(*p.skeleton_, face_id);
				if (!f.is_valid())
					continue;

				const std::vector<NMEdge> affected_edges = incident_edges(*p.skeleton_, f);
				remove_face(*p.skeleton_, f);
				removed_any_faces = true;
				if (collect_verbose_stats)
					++removed_faces_this_pass;

				std::vector<NMVertex> candidate_vertices;
				candidate_vertices.reserve(affected_edges.size() * 2);
				for (const NMEdge& e : affected_edges)
				{
					if (!e.is_valid())
						continue;
					for (const NMVertex& v : incident_vertices(*p.skeleton_, e))
						candidate_vertices.push_back(v);
					if (incident_faces(*p.skeleton_, e).empty())
					{
						remove_edge(*p.skeleton_, e);
						if (collect_verbose_stats)
							++removed_edges_this_pass;
					}
				}
				for (const NMVertex& v : candidate_vertices)
				{
					if (!v.is_valid())
						continue;
					const uint32 vertex_id = index_of(*p.skeleton_, v);
					if (vertex_id == INVALID_INDEX || !incident_edges(*p.skeleton_, v).empty())
						continue;
					if (remove_skeleton_vertex_and_linked_sphere(p, v))
					{
						if (collect_verbose_stats)
							++removed_vertices_this_pass;
					}
				}
			}

			if (collect_verbose_stats)
			{
				total_removed_faces += removed_faces_this_pass;
				total_removed_edges += removed_edges_this_pass;
				total_removed_vertices += removed_vertices_this_pass;
			}

		}

		if (p.edge_degree_ && removed_any_faces)
			compute_edge_degree(p);

		log_verbose(p, log_prefix, " total_passes=", total_passes, " total_removed_faces=", total_removed_faces,
					" total_removed_edges=", total_removed_edges, " total_removed_vertices=",
					total_removed_vertices, '\n');
		return true;
	}

	bool compute_skeleton_face_components_union_find(PointsParameters& p, const char* log_prefix = "[FaceComponentsUF]")
	{
		if (!p.skeleton_ || !p.skeleton_face_component_id_)
		{
			log_error(p, log_prefix, " requires skeleton and face component id attribute.", '\n');
			return false;
		}

		std::vector<NMFace> faces;
		faces.reserve(nb_cells<NMFace>(*p.skeleton_));
		std::unordered_map<uint32, uint32> face_id_to_uf_index;
		face_id_to_uf_index.reserve(nb_cells<NMFace>(*p.skeleton_));
		foreach_cell(*p.skeleton_, [&](NMFace f) -> bool {
			if (!f.is_valid())
				return true;
			const uint32 face_id = index_of(*p.skeleton_, f);
			if (face_id == INVALID_INDEX)
				return true;
			(*p.skeleton_face_component_id_)[face_id] = INVALID_INDEX;
			face_id_to_uf_index.emplace(face_id, static_cast<uint32>(faces.size()));
			faces.push_back(f);
			return true;
		});

		FaceComponentUnionFind uf(faces.size());
		foreach_cell(*p.skeleton_, [&](NMEdge e) -> bool {
			if (!e.is_valid())
				return true;
			const auto incident = incident_faces(*p.skeleton_, e);
			if (incident.size() != 2)
				return true;
			const uint32 f0 = index_of(*p.skeleton_, incident[0]);
			const uint32 f1 = index_of(*p.skeleton_, incident[1]);
			const auto i0 = face_id_to_uf_index.find(f0);
			const auto i1 = face_id_to_uf_index.find(f1);
			if (i0 != face_id_to_uf_index.end() && i1 != face_id_to_uf_index.end())
				uf.unite(i0->second, i1->second);
			return true;
		});

		std::unordered_map<uint32, uint32> root_to_component_id;
		root_to_component_id.reserve(faces.size());
		for (uint32 i = 0; i < static_cast<uint32>(faces.size()); ++i)
		{
			const uint32 root = uf.find(i);
			const auto [it, inserted] =
				root_to_component_id.emplace(root, static_cast<uint32>(root_to_component_id.size()));
			(void)inserted;
			const uint32 face_id = index_of(*p.skeleton_, faces[i]);
			if (face_id != INVALID_INDEX)
				(*p.skeleton_face_component_id_)[face_id] = it->second;
		}
		log_verbose(p, log_prefix, " faces=", faces.size(), " components=", root_to_component_id.size(), '\n');
		return true;
	}

	bool colorize_skeleton_face_components(PointsParameters& p, const char* log_prefix = "[FaceComponentsColor]")
	{
		if (!p.skeleton_ || !p.skeleton_face_component_id_ || !p.skeleton_face_component_color_)
		{
			log_error(p, log_prefix, " requires skeleton, component id and component color attributes.", '\n');
			return false;
		}

		const bool collect_verbose_stats = is_verbose_logging_enabled(p);
		uint32 colored_faces = 0;
		std::unordered_map<uint32, Vec3> component_colors;
		component_colors.reserve(nb_cells<NMFace>(*p.skeleton_));

		foreach_cell(*p.skeleton_, [&](NMFace f) -> bool {
			if (!f.is_valid())
				return true;
			const uint32 idf = index_of(*p.skeleton_, f);
			if (idf == INVALID_INDEX)
				return true;

			const uint32 component_id = (*p.skeleton_face_component_id_)[idf];
			Vec3 color(0.0, 0.0, 0.0);
			if (component_id != INVALID_INDEX)
			{
				auto it = component_colors.find(component_id);
				if (it == component_colors.end())
					it = component_colors.emplace(component_id, component_palette_color(component_id)).first;
				color = it->second;
				if (collect_verbose_stats)
					++colored_faces;
			}
			(*p.skeleton_face_component_color_)[idf] = color;
			return true;
		});

		log_verbose(p, log_prefix, " colored_faces=", colored_faces, " components=", component_colors.size(), '\n');
		return true;
	}

	void run_residual_sheet_prune(
		PointsParameters& p, const char* log_prefix = "[ResidualSheetPrune]")
	{
		if (!p.skeleton_ || !p.skeleton_face_component_id_ || !p.incident_tets_)
		{
			log_error(p, log_prefix, " requires skeleton/component-id/incident_tets.", '\n');
			return;
		}

		compute_edge_degree(p);
		prune_fully_non_manifold_triangles(p, "[FaceComponentsPreDelete]");
		if (!compute_skeleton_face_components_union_find(p, "[FaceComponentsUF-Residual]"))
		{
			log_error(p, log_prefix, " failed to rebuild face components.", '\n');
			return;
		}

		std::unordered_map<uint32, std::vector<uint32>> label_to_faces;
		label_to_faces.reserve(nb_cells<NMFace>(*p.skeleton_));
		std::unordered_set<uint32> all_face_ids;
		all_face_ids.reserve(nb_cells<NMFace>(*p.skeleton_));
		foreach_cell(*p.skeleton_, [&](NMFace f) -> bool {
			if (!f.is_valid())
				return true;
			const uint32 face_id = index_of(*p.skeleton_, f);
			if (face_id == INVALID_INDEX)
				return true;
			const uint32 label = (*p.skeleton_face_component_id_)[face_id];
			if (label == INVALID_INDEX)
				return true;
			label_to_faces[label].push_back(face_id);
			all_face_ids.insert(face_id);
			return true;
		});
		if (all_face_ids.empty())
		{
			log_basic(p, log_prefix, " no valid faces.", '\n');
			return;
		}

		std::vector<uint32> sorted_sheet_labels;
		sorted_sheet_labels.reserve(label_to_faces.size());
		for (const auto& kv : label_to_faces)
			sorted_sheet_labels.push_back(kv.first);
		std::sort(sorted_sheet_labels.begin(), sorted_sheet_labels.end());

		const size_t forced_delete_max_face_count = 2;
		const bool collect_basic_stats = is_basic_logging_enabled(p);
		const bool collect_verbose_stats = is_verbose_logging_enabled(p);
		auto sheet_has_vertex_connected_to_isolated_edge = [&](const std::vector<uint32>& face_ids) -> bool {
			std::unordered_set<uint32> visited_vertex_ids;
			visited_vertex_ids.reserve(face_ids.size() * 3);
			for (uint32 face_id : face_ids)
			{
				if (face_id == INVALID_INDEX)
					continue;
				const NMFace f = of_index<NMFace>(*p.skeleton_, face_id);
				if (!f.is_valid() || index_of(*p.skeleton_, f) != face_id)
					continue;
				const std::vector<NMVertex> vertices = incident_vertices(*p.skeleton_, f);
				for (const NMVertex& v : vertices)
				{
					if (!v.is_valid())
						continue;
					const uint32 vertex_id = index_of(*p.skeleton_, v);
					if (vertex_id == INVALID_INDEX || !visited_vertex_ids.insert(vertex_id).second)
						continue;
					for (const NMEdge& e : incident_edges(*p.skeleton_, v))
					{
						if (!e.is_valid())
							continue;
						if (incident_faces(*p.skeleton_, e).empty())
							return true;
					}
				}
			}
			return false;
		};

		struct ResidualSheetScoreEntry
		{
			uint32 sheet_label = INVALID_INDEX;
			size_t face_count = 0;
			bool delete_sheet = false;
			bool delete_by_small_sheet = false;
			bool keep_by_isolated_edge_guard = false;
		};

		uint32 deleted_sheet_count = 0;
		std::unordered_set<uint32> face_ids_to_delete;
		std::vector<ResidualSheetScoreEntry> residual_sheet_scores;
		if (collect_verbose_stats)
			residual_sheet_scores.reserve(sorted_sheet_labels.size());
		for (uint32 sheet_label : sorted_sheet_labels)
		{
			const auto it_faces = label_to_faces.find(sheet_label);
			if (it_faces == label_to_faces.end() || it_faces->second.empty())
				continue;

			const size_t face_count = it_faces->second.size();
			const bool delete_by_small_sheet = (face_count <= forced_delete_max_face_count);
			const bool keep_by_isolated_edge_guard =
				delete_by_small_sheet && sheet_has_vertex_connected_to_isolated_edge(it_faces->second);
			const bool delete_sheet = delete_by_small_sheet && !keep_by_isolated_edge_guard;
			if (collect_verbose_stats)
			{
				residual_sheet_scores.push_back(ResidualSheetScoreEntry{
					sheet_label,
					face_count,
					delete_sheet,
					delete_by_small_sheet,
					keep_by_isolated_edge_guard});
			}

			if (!delete_sheet)
				continue;

			if (collect_basic_stats)
				++deleted_sheet_count;
			face_ids_to_delete.reserve(face_ids_to_delete.size() + it_faces->second.size());
			for (uint32 face_id : it_faces->second)
				face_ids_to_delete.insert(face_id);
		}

		if (collect_verbose_stats)
		{
			log_verbose(p, log_prefix, " current_residual_sheet_scores", " evaluated_sheets=",
						residual_sheet_scores.size(), " delete_face_count_leq=", forced_delete_max_face_count, '\n');
			for (const ResidualSheetScoreEntry& entry : residual_sheet_scores)
			{
				const char* action =
					entry.delete_sheet ? "delete_small_sheet" :
					(entry.keep_by_isolated_edge_guard ? "keep_isolated_edge_guard" : "keep");
				log_verbose(p, log_prefix, " residual_sheet_score", " sheet=", entry.sheet_label,
							" face_count=", entry.face_count, " small_sheet_delete=",
							(entry.delete_by_small_sheet ? "true" : "false"), " isolated_edge_guard=",
							(entry.keep_by_isolated_edge_guard ? "true" : "false"), " action=", action, '\n');
			}
		}

		if (face_ids_to_delete.empty())
		{
			if (collect_verbose_stats)
			{
				log_verbose(p, log_prefix, " done", " deleted_sheets=0", " delete_face_count_leq=",
							forced_delete_max_face_count, '\n');
			}
			else
			{
				log_basic(p, log_prefix, " done", " deleted_sheets=0", " delete_face_count_leq=",
						  forced_delete_max_face_count, '\n');
			}
			compute_edge_degree(p);
			prune_fully_non_manifold_triangles(p, "[FaceComponentsPreDelete]");
			if (compute_skeleton_face_components_union_find(p, "[FaceComponentsUF-Residual-Final]"))
				colorize_skeleton_face_components(p, "[FaceComponentsColor-Residual-Final]");
			return;
		}

		SkeletonFaceDeletionStats deletion_stats;
		if (!delete_skeleton_face_id_set(p, face_ids_to_delete, deletion_stats))
		{
			log_basic(p, log_prefix, " selected ", deleted_sheet_count,
					  " sheet(s) but no face was actually removed.", '\n');
			return;
		}

		refresh_skeleton_topology_colors(p);
		compute_edge_degree(p);
		prune_fully_non_manifold_triangles(p, "[FaceComponentsPreDelete]");
		if (compute_skeleton_face_components_union_find(p, "[FaceComponentsUF-Residual-Final]"))
			colorize_skeleton_face_components(p, "[FaceComponentsColor-Residual-Final]");

		if (collect_verbose_stats)
		{
			log_verbose(p, log_prefix, " done", " deleted_sheets=", deleted_sheet_count, " delete_face_count_leq=",
						forced_delete_max_face_count, " removed_faces=",
						deletion_stats.removed_faces, " removed_edges=",
						deletion_stats.removed_edges, " removed_vertices=", deletion_stats.removed_vertices,
						" removed_tets=", deletion_stats.removed_tets, '\n');
		}
		else
		{
			log_basic(p, log_prefix, " done", " deleted_sheets=", deleted_sheet_count, " delete_face_count_leq=",
					  forced_delete_max_face_count, " removed_faces=",
					  deletion_stats.removed_faces, " removed_edges=", deletion_stats.removed_edges,
					  " removed_vertices=", deletion_stats.removed_vertices, " removed_tets=",
					  deletion_stats.removed_tets, '\n');
		}
	}
	uint32 remove_orphan_edges_from_removed_face_edges(PointsParameters& p, const std::vector<NMEdge>& face_edges,
													   uint32* out_removed_vertices = nullptr)
	{
		uint32 removed_edges = 0;
		std::vector<NMVertex> candidate_vertices;
		candidate_vertices.reserve(face_edges.size() * 2);
		for (NMEdge e : face_edges)
		{
			if (!e.is_valid())
				continue;
			const std::vector<NMVertex> edge_vertices = incident_vertices(*p.skeleton_, e);
			for (const NMVertex& v : edge_vertices)
				candidate_vertices.push_back(v);
			if (incident_faces(*p.skeleton_, e).empty())
			{
				remove_edge(*p.skeleton_, e);
				++removed_edges;
			}
		}
		uint32 removed_vertices = 0;
		for (const NMVertex& v : candidate_vertices)
		{
			if (!v.is_valid())
				continue;
			const uint32 idv = index_of(*p.skeleton_, v);
			if (idv == INVALID_INDEX)
				continue;
			if (!incident_edges(*p.skeleton_, v).empty())
				continue;
			if (remove_skeleton_vertex_and_linked_sphere(p, v))
				++removed_vertices;
		}
		if (out_removed_vertices)
			*out_removed_vertices = removed_vertices;
		return removed_edges;
	}

	std::unordered_set<uint32> collect_current_orphan_edge_ids(PointsParameters& p)
	{
		std::unordered_set<uint32> orphan_edge_ids;
		if (!p.skeleton_)
			return orphan_edge_ids;
		orphan_edge_ids.reserve(nb_cells<NMEdge>(*p.skeleton_));
		foreach_cell(*p.skeleton_, [&](NMEdge e) -> bool {
			if (!e.is_valid())
				return true;
			const uint32 ide = index_of(*p.skeleton_, e);
			if (ide == INVALID_INDEX)
				return true;
			if (incident_faces(*p.skeleton_, e).empty())
				orphan_edge_ids.insert(ide);
			return true;
		});
		return orphan_edge_ids;
	}

	void remove_global_orphan_skeleton_elements(PointsParameters& p, uint32* out_removed_edges = nullptr,
												uint32* out_removed_vertices = nullptr,
												const std::unordered_set<uint32>* preserved_orphan_edge_ids = nullptr)
	{
		uint32 removed_edges = 0;
		uint32 removed_vertices = 0;
		if (!p.skeleton_)
		{
			if (out_removed_edges)
				*out_removed_edges = 0;
			if (out_removed_vertices)
				*out_removed_vertices = 0;
			return;
		}

		std::vector<NMEdge> orphan_edges;
		orphan_edges.reserve(nb_cells<NMEdge>(*p.skeleton_));
		foreach_cell(*p.skeleton_, [&](NMEdge e) -> bool {
			if (!e.is_valid())
				return true;
			if (!incident_faces(*p.skeleton_, e).empty())
				return true;
			const uint32 ide = index_of(*p.skeleton_, e);
			if (ide == INVALID_INDEX)
				return true;
			if (preserved_orphan_edge_ids && preserved_orphan_edge_ids->find(ide) != preserved_orphan_edge_ids->end())
				return true;
			orphan_edges.push_back(e);
			return true;
		});

		std::vector<NMVertex> candidate_vertices;
		candidate_vertices.reserve(orphan_edges.size() * 2);
		for (const NMEdge& e : orphan_edges)
		{
			if (!e.is_valid())
				continue;
			for (const NMVertex& v : incident_vertices(*p.skeleton_, e))
				candidate_vertices.push_back(v);
			remove_edge(*p.skeleton_, e);
			++removed_edges;
		}

		foreach_cell(*p.skeleton_, [&](NMVertex v) -> bool {
			if (!v.is_valid())
				return true;
			candidate_vertices.push_back(v);
			return true;
		});

		std::unordered_set<uint32> seen_vertex_ids;
		seen_vertex_ids.reserve(candidate_vertices.size());

		for (const NMVertex& v : candidate_vertices)
		{
			if (!v.is_valid())
				continue;
			const uint32 idv = index_of(*p.skeleton_, v);
			if (idv == INVALID_INDEX)
				continue;
			if (!seen_vertex_ids.insert(idv).second)
				continue;
			if (!incident_edges(*p.skeleton_, v).empty())
				continue;
			if (remove_skeleton_vertex_and_linked_sphere(p, v))
				++removed_vertices;
		}

		if (out_removed_edges)
			*out_removed_edges = removed_edges;
		if (out_removed_vertices)
			*out_removed_vertices = removed_vertices;
	}

	struct SkeletonFaceDeletionStats
	{
		uint32 removed_faces = 0;
		uint32 removed_edges = 0;
		uint32 removed_vertices = 0;
		uint32 removed_tets = 0;
	};

	struct SkeletonFaceDeletionOptions
	{
		bool remove_incident_orphan_edges_immediately = true;
		bool remove_global_orphan_elements_after_batch = true;
	};

	bool delete_skeleton_face_id_set(PointsParameters& p, const std::unordered_set<uint32>& face_ids_to_delete,
									 SkeletonFaceDeletionStats& out_stats,
									 const SkeletonFaceDeletionOptions& options = {})
	{
		out_stats = SkeletonFaceDeletionStats{};
		if (!p.skeleton_ || !p.incident_tets_ || face_ids_to_delete.empty())
			return false;
		const std::unordered_set<uint32> preexisting_orphan_edge_ids =
			options.remove_global_orphan_elements_after_batch ? collect_current_orphan_edge_ids(p) : std::unordered_set<uint32>{};

		for (uint32 face_id : face_ids_to_delete)
		{
			const NMFace face = of_index<NMFace>(*p.skeleton_, face_id);
			if (!face.is_valid())
				continue;
			const uint32 current_face_id = index_of(*p.skeleton_, face);
			if (current_face_id == INVALID_INDEX || current_face_id != face_id)
				continue;

			bool face_still_attached = false;
			for (NMEdge e : incident_edges(*p.skeleton_, face))
			{
				if (!e.is_valid())
					continue;
				for (NMFace ef : incident_faces(*p.skeleton_, e))
				{
					if (ef == face)
					{
						face_still_attached = true;
						break;
					}
				}
				if (face_still_attached)
					break;
			}
			if (!face_still_attached)
				continue;

			const std::set<std::size_t> in_tets = (*p.incident_tets_)[face_id];
			std::vector<NMEdge> affected_edges = incident_edges(*p.skeleton_, face);
			remove_face(*p.skeleton_, face);
			++out_stats.removed_faces;

			if (options.remove_incident_orphan_edges_immediately)
			{
				uint32 removed_vertices_this_face = 0;
				out_stats.removed_edges +=
					remove_orphan_edges_from_removed_face_edges(p, affected_edges, &removed_vertices_this_face);
				out_stats.removed_vertices += removed_vertices_this_face;
			}

			for (std::size_t tet_id : in_tets)
			{
				auto it_tet = p.skeleton_tets_.find(tet_id);
				if (it_tet == p.skeleton_tets_.end())
					continue;
				const Tet old_tet = it_tet->second;
				for (uint32 i = 0; i < 4; ++i)
				{
					const NMFace f = old_tet.faces[i];
					if (!f.is_valid() || f == face)
						continue;
					const uint32 idf = index_of(*p.skeleton_, f);
					if (idf != INVALID_INDEX)
						(*p.incident_tets_)[idf].erase(tet_id);
				}
				p.skeleton_tets_.erase(it_tet);
				++out_stats.removed_tets;
			}
		}

		if (out_stats.removed_faces == 0)
			return false;

		if (options.remove_global_orphan_elements_after_batch)
		{
			uint32 global_orphan_edges = 0;
			uint32 global_orphan_vertices = 0;
			remove_global_orphan_skeleton_elements(
				p, &global_orphan_edges, &global_orphan_vertices, &preexisting_orphan_edge_ids);
			out_stats.removed_edges += global_orphan_edges;
			out_stats.removed_vertices += global_orphan_vertices;
		}
		return true;
	}

	bool compute_skeleton_face_scores(PointsParameters& p, const std::unordered_set<uint32>& face_ids,
								  std::unordered_map<uint32, Scalar>& face_scores, bool normalize_by_area,
								  const char* log_prefix = "[FaceUDF]")
	{
		if (!p.skeleton_ || !p.skeleton_position_)
			return false;
		if (face_ids.empty())
		{
			face_scores.clear();
			return true;
		}

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

		face_scores.clear();
		face_scores.reserve(face_ids.size());
		std::vector<AdaptiveTriangleTask> pending_tris;
		pending_tris.reserve(face_ids.size() * 2);
		foreach_cell(*p.skeleton_, [&](NMFace f) -> bool {
			const uint32 idf = index_of(*p.skeleton_, f);
			if (idf == INVALID_INDEX)
				return true;
			if (face_ids.find(idf) == face_ids.end())
				return true;
			face_scores[idf] = Scalar(0);
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
				log_error(p, log_prefix, " face score field evaluation failed.", '\n');
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
					std::max(Scalar(1e-4), Scalar(0.35) * std::max(weighted_mean, Scalar(0)));
				const bool should_refine = (tri.depth < 2) && (udf_range > refine_threshold);

				if (should_refine)
				{
					subdivide_triangle(tri, next_tris);
				}
				else
				{
					face_scores[tri.face_id] += tri_integral;
				}
			}
			pending_tris.swap(next_tris);
		}

		if (normalize_by_area)
		{
			foreach_cell(*p.skeleton_, [&](NMFace f) -> bool {
				const uint32 face_id = index_of(*p.skeleton_, f);
				auto it_score = face_scores.find(face_id);
				if (face_id == INVALID_INDEX || it_score == face_scores.end())
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
				it_score->second = (area > Scalar(0)) ? (it_score->second / area) : Scalar(0);
				return true;
			});
		}
		return true;
	}

	std::unordered_set<uint32> collect_current_tet_edge_ids(PointsParameters& p)
	{
		std::unordered_set<uint32> tet_edge_ids;
		if (!p.skeleton_)
			return tet_edge_ids;
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
		return tet_edge_ids;
	}

	bool compute_skeleton_edge_scores_gauss3_normalized(PointsParameters& p, const std::unordered_set<uint32>& edge_ids,
										std::unordered_map<uint32, Scalar>& edge_scores,
										const char* log_prefix = "[EdgeUDF]")
	{
		if (!p.skeleton_ || !p.skeleton_position_)
			return false;
		edge_scores.clear();
		if (edge_ids.empty())
			return true;

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

		edge_scores.reserve(edge_ids.size());
		std::vector<EdgeBatchMeta> edge_batch;
		edge_batch.reserve(edge_ids.size());
		std::vector<Vec3> sample_points;
		sample_points.reserve(edge_ids.size() * 3);

		foreach_cell(*p.skeleton_, [&](NMEdge e) -> bool {
			if (!e.is_valid())
				return true;
			const uint32 ide = index_of(*p.skeleton_, e);
			if (ide == INVALID_INDEX)
				return true;
			if (edge_ids.find(ide) == edge_ids.end())
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
			log_error(p, log_prefix, " edge score field evaluation failed.", '\n');
			return false;
		}

		for (const EdgeBatchMeta& edge_eval : edge_batch)
		{
			Scalar avg_udf = Scalar(0);
			for (uint32 k = 0; k < 3; ++k)
				avg_udf += gauss3_w_01[k] * score_values[edge_eval.sample_offset + k];
			edge_scores[edge_eval.edge_id] = avg_udf;
		}
		return true;
	}

	struct TopologyFixScoreCache
	{
		std::unordered_map<uint32, Scalar> edge_scores;
		std::unordered_map<uint32, Scalar> face_scores;
		bool initialized = false;
	};

	bool initialize_topology_fix_score_cache(PointsParameters& p, TopologyFixScoreCache& cache,
												 const char* log_prefix = "[TopologyFixScore]")
	{
		if (cache.initialized)
			return true;
		const std::unordered_set<uint32> tet_face_ids = collect_current_tet_face_id_whitelist(p);
		const std::unordered_set<uint32> tet_edge_ids = collect_current_tet_edge_ids(p);
		if (!compute_skeleton_edge_scores_gauss3_normalized(p, tet_edge_ids, cache.edge_scores, log_prefix))
		{
			log_error(p, log_prefix, " failed to compute edge scores.", '\n');
			return false;
		}
		if (!compute_skeleton_face_scores(p, tet_face_ids, cache.face_scores, false, log_prefix))
		{
			log_error(p, log_prefix, " failed to compute face scores.", '\n');
			return false;
		}
		cache.initialized = true;
		return true;
	}

	void erase_topology_fix_scores_for_removed_faces_and_orphan_edges(
		PointsParameters& p, TopologyFixScoreCache* cache, const std::unordered_set<uint32>& removed_face_ids,
		const std::vector<NMEdge>& affected_edges)
	{
		if (!cache)
			return;
		for (uint32 face_id : removed_face_ids)
			cache->face_scores.erase(face_id);

		std::unordered_set<uint32> seen_edge_ids;
		seen_edge_ids.reserve(affected_edges.size() * 2 + 1);
		for (NMEdge e : affected_edges)
		{
			if (!e.is_valid())
				continue;
			const uint32 ide = index_of(*p.skeleton_, e);
			if (ide == INVALID_INDEX || !seen_edge_ids.insert(ide).second)
				continue;
			if (incident_faces(*p.skeleton_, e).empty())
				cache->edge_scores.erase(ide);
		}
	}

	enum class EdgeTetDeleteMode
	{
		SimpleTet,
		NonSimpleTet
	};

	struct EdgeTetModeRunStats
	{
		uint32 steps = 0;
		uint32 removed_faces = 0;
		uint32 removed_edges = 0;
		uint32 removed_tets = 0;
		uint32 remaining_tets = 0;
	};

	EdgeTetModeRunStats run_edge_score_tet_mode_topology_fix(
		PointsParameters& p, EdgeTetDeleteMode mode, TopologyFixScoreCache& score_cache)
	{
		EdgeTetModeRunStats stats;
		stats.remaining_tets = static_cast<uint32>(p.skeleton_tets_.size());
		const char* log_tag = (mode == EdgeTetDeleteMode::SimpleTet) ? "[EdgeTetSimple]" : "[EdgeTetNonSimple]";
		if (!ensure_topology_score_backend(p, log_tag))
		{
			return stats;
		}
		if (!p.skeleton_ || !p.incident_tets_)
		{
			log_error(p, log_tag, " requires a built skeleton.", '\n');
			return stats;
		}

		if (!initialize_topology_fix_score_cache(p, score_cache, log_tag))
			return stats;
		auto& edge_score_cache = score_cache.edge_scores;
		auto& face_score_cache = score_cache.face_scores;

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
		auto pick_tet_best_edge_and_face = [&](std::size_t tet_id, NMEdge& out_edge, uint32& out_edge_id,
											  Scalar& out_edge_score, NMFace& out_face, uint32& out_face_id,
											  Scalar& out_face_score) -> bool {
			auto it_tet = p.skeleton_tets_.find(tet_id);
			if (it_tet == p.skeleton_tets_.end())
				return false;
			const Tet& tet = it_tet->second;
			std::unordered_set<uint32> seen_edges;
			seen_edges.reserve(8);
			bool found = false;
			const Scalar edge_eps = Scalar(1e-12);
			const Scalar face_eps = Scalar(1e-12);
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
					NMFace best_face_on_edge;
					uint32 best_face_on_edge_id = INVALID_INDEX;
					Scalar best_face_on_edge_score = Scalar(0);
					bool edge_has_candidate = false;
					for (const NMFace& ef : incident_faces(*p.skeleton_, e))
					{
						if (!face_belongs_to_tet(ef, tet_id))
							continue;
						const uint32 idf = index_of(*p.skeleton_, ef);
						if (idf == INVALID_INDEX || !has_face_budget(idf))
							continue;
						if (mode == EdgeTetDeleteMode::SimpleTet && (*p.incident_tets_)[idf].size() != 1)
							continue;
						const Scalar face_score = get_face_score(idf);
						if (!edge_has_candidate || face_score > best_face_on_edge_score + face_eps)
						{
							edge_has_candidate = true;
							best_face_on_edge = ef;
							best_face_on_edge_id = idf;
							best_face_on_edge_score = face_score;
						}
					}
					if (!edge_has_candidate)
						continue;
					const auto it_score = edge_score_cache.find(ide);
					const Scalar score = (it_score != edge_score_cache.end()) ? it_score->second : Scalar(0);
					const bool better_edge = (!found) || (score > out_edge_score + edge_eps);
					const bool tie_edge = found && (std::abs(score - out_edge_score) <= edge_eps);
					const bool better_face_on_tie = tie_edge && (best_face_on_edge_score > out_face_score + face_eps);
					if (better_edge || better_face_on_tie)
					{
						found = true;
						out_edge = e;
						out_edge_id = ide;
						out_edge_score = score;
						out_face = best_face_on_edge;
						out_face_id = best_face_on_edge_id;
						out_face_score = best_face_on_edge_score;
					}
				}
			}
			return found;
		};
		auto remove_faces_collect_tets =
			[&](const std::vector<NMFace>& faces_to_remove, std::unordered_set<std::size_t>& out_tets_to_erase,
				uint32& out_removed_faces, uint32& out_removed_edges,
				std::vector<NMEdge>* out_affected_edges = nullptr) -> bool {
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
			if (out_affected_edges)
				out_affected_edges->insert(out_affected_edges->end(), affected_edges.begin(), affected_edges.end());

			std::unordered_set<uint32> removed_face_ids;
			removed_face_ids.reserve(unique_faces.size() * 2 + 1);
			for (const auto& fpair : unique_faces)
			{
				const NMFace f = fpair.second;
				if (!f.is_valid())
					continue;
				const uint32 idf_now = index_of(*p.skeleton_, f);
				if (idf_now == INVALID_INDEX)
					continue;
				account_face_deletion_budget(idf_now);
				removed_face_ids.insert(idf_now);
				remove_face(*p.skeleton_, f);
				++out_removed_faces;
			}

			erase_topology_fix_scores_for_removed_faces_and_orphan_edges(p, &score_cache, removed_face_ids, affected_edges);
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
		struct TetCandidateQueueEntry
		{
			std::size_t tet_id = 0;
			uint32 edge_id = INVALID_INDEX;
			uint32 face_id = INVALID_INDEX;
			Scalar edge_score = Scalar(0);
			Scalar face_score = Scalar(0);
			uint32 version = 0;
		};
		struct TetCandidateQueueCompare
		{
			bool operator()(const TetCandidateQueueEntry& a, const TetCandidateQueueEntry& b) const
			{
				if (a.edge_score != b.edge_score)
					return a.edge_score < b.edge_score;
				if (a.face_score != b.face_score)
					return a.face_score < b.face_score;
				return a.tet_id > b.tet_id;
			}
		};
		using TetCandidateQueue =
			std::priority_queue<TetCandidateQueueEntry, std::vector<TetCandidateQueueEntry>, TetCandidateQueueCompare>;

		uint32 removed_faces = 0;
		uint32 removed_edges = 0;
		uint32 removed_tets = 0;
		uint32 processed_steps = 0;
		std::unordered_map<std::size_t, uint32> tet_versions;
		tet_versions.reserve(p.skeleton_tets_.size() * 2 + 1);
		TetCandidateQueue candidate_queue;

		auto build_step_candidate = [&](std::size_t tet_id, StepCandidate& out_candidate) -> bool {
			const bool is_simple_tet = tet_has_simple_face(tet_id);
			if (mode == EdgeTetDeleteMode::SimpleTet && !is_simple_tet)
				return false;
			if (mode == EdgeTetDeleteMode::NonSimpleTet && is_simple_tet)
				return false;

			NMEdge best_edge_local;
			uint32 best_edge_local_id = INVALID_INDEX;
			Scalar best_edge_local_score = Scalar(0);
			NMFace first_face_local;
			uint32 first_face_local_id = INVALID_INDEX;
			Scalar first_face_local_score = Scalar(0);
			if (!pick_tet_best_edge_and_face(
					tet_id, best_edge_local, best_edge_local_id, best_edge_local_score, first_face_local,
					first_face_local_id, first_face_local_score))
				return false;

			out_candidate.tet_id = tet_id;
			out_candidate.edge = best_edge_local;
			out_candidate.edge_id = best_edge_local_id;
			out_candidate.edge_score = best_edge_local_score;
			out_candidate.first_face = first_face_local;
			out_candidate.first_face_id = first_face_local_id;
			out_candidate.first_face_score = first_face_local_score;
			return true;
		};
		auto queue_tet_candidate = [&](std::size_t tet_id, bool bump_version) {
			if (p.skeleton_tets_.find(tet_id) == p.skeleton_tets_.end())
			{
				tet_versions.erase(tet_id);
				return;
			}
			uint32& version = tet_versions[tet_id];
			if (bump_version)
				++version;
			StepCandidate candidate;
			if (!build_step_candidate(tet_id, candidate))
				return;
			candidate_queue.push(
				{candidate.tet_id, candidate.edge_id, candidate.first_face_id, candidate.edge_score,
				 candidate.first_face_score, version});
		};
		auto collect_dirty_tets_from_edges =
			[&](const std::vector<NMEdge>& affected_edges, std::unordered_set<std::size_t>& dirty_tets) {
				dirty_tets.clear();
				for (NMEdge e : affected_edges)
				{
					if (!e.is_valid())
						continue;
					for (NMFace f : incident_faces(*p.skeleton_, e))
					{
						if (!f.is_valid())
							continue;
						const uint32 idf = index_of(*p.skeleton_, f);
						if (idf == INVALID_INDEX)
							continue;
						for (std::size_t tet_id : (*p.incident_tets_)[idf])
							if (p.skeleton_tets_.find(tet_id) != p.skeleton_tets_.end())
								dirty_tets.insert(tet_id);
					}
				}
			};

		for (const auto& kv : p.skeleton_tets_)
		{
			tet_versions[kv.first] = 0;
			queue_tet_candidate(kv.first, false);
		}

		while (true)
		{
			StepCandidate best;
			bool has_best = false;
			const Scalar edge_eps = Scalar(1e-12);
			const Scalar face_eps = Scalar(1e-12);
			while (!candidate_queue.empty())
			{
				const TetCandidateQueueEntry entry = candidate_queue.top();
				candidate_queue.pop();

				auto it_version = tet_versions.find(entry.tet_id);
				if (it_version == tet_versions.end() || it_version->second != entry.version)
					continue;
				if (p.skeleton_tets_.find(entry.tet_id) == p.skeleton_tets_.end())
				{
					tet_versions.erase(entry.tet_id);
					continue;
				}

				StepCandidate refreshed;
				if (!build_step_candidate(entry.tet_id, refreshed))
				{
					++(it_version->second);
					continue;
				}

				const bool matches_cached =
					(refreshed.edge_id == entry.edge_id) && (refreshed.first_face_id == entry.face_id) &&
					(std::abs(refreshed.edge_score - entry.edge_score) <= edge_eps) &&
					(std::abs(refreshed.first_face_score - entry.face_score) <= face_eps);
				if (!matches_cached)
				{
					++(it_version->second);
					candidate_queue.push(
						{refreshed.tet_id, refreshed.edge_id, refreshed.first_face_id, refreshed.edge_score,
						 refreshed.first_face_score, it_version->second});
					continue;
				}

				best = refreshed;
				has_best = true;
				break;
			}

			if (!has_best)
			{
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
			std::vector<NMEdge> step_affected_edges;
			if (!remove_faces_collect_tets(
					{best.first_face}, first_tets_to_erase, step_removed_faces, step_removed_edges, &step_affected_edges))
			{
				queue_tet_candidate(best.tet_id, true);
				continue;
			}
			step_tets_to_erase.insert(first_tets_to_erase.begin(), first_tets_to_erase.end());

			std::vector<FaceChoice> created_deg1_faces;
			created_deg1_faces.reserve(16);
			std::unordered_set<uint32> seen_created_ids;
			seen_created_ids.reserve(32);
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
						continue;
					}
					if (idf == best.first_face_id)
						continue;
					if (!seen_created_ids.insert(idf).second)
						continue;
					if (!has_face_budget(idf))
						continue;
					if (!face_has_degree1_edge(f))
						continue;
					if (face_has_edge_with_tet_face_count_gt2(f, &step_tets_to_erase))
						continue;
					created_deg1_faces.push_back({f, idf, get_face_score(idf)});
				}
			}
			std::sort(created_deg1_faces.begin(), created_deg1_faces.end(),
					  [](const FaceChoice& a, const FaceChoice& b) { return a.score > b.score; });
			if (mode == EdgeTetDeleteMode::SimpleTet && created_deg1_faces.size() > 1)
				created_deg1_faces.resize(1);

			for (const FaceChoice& fc : created_deg1_faces)
			{
				if (!fc.face.is_valid())
					continue;
				const uint32 idf_now = index_of(*p.skeleton_, fc.face);
				if (idf_now == INVALID_INDEX)
					continue;
				if (!has_face_budget(idf_now))
					continue;
				if (!face_has_degree1_edge(fc.face))
					continue;
				if (face_has_edge_with_tet_face_count_gt2(fc.face, &step_tets_to_erase))
					continue;
				std::unordered_set<std::size_t> second_tets_to_erase;
				uint32 rf = 0, re = 0;
				if (remove_faces_collect_tets(
						{fc.face}, second_tets_to_erase, rf, re, &step_affected_edges))
				{
					step_removed_faces += rf;
					step_removed_edges += re;
					step_tets_to_erase.insert(second_tets_to_erase.begin(), second_tets_to_erase.end());
				}
			}
			finalize_erase_tets(step_tets_to_erase, step_removed_tets);
			for (std::size_t tet_id : step_tets_to_erase)
				tet_versions.erase(tet_id);

			std::unordered_set<std::size_t> dirty_tets;
			collect_dirty_tets_from_edges(step_affected_edges, dirty_tets);
			for (std::size_t tet_id : dirty_tets)
				queue_tet_candidate(tet_id, true);

			removed_faces += step_removed_faces;
			removed_edges += step_removed_edges;
			removed_tets += step_removed_tets;
			++processed_steps;
		}

		stats.steps = processed_steps;
		stats.removed_faces = removed_faces;
		stats.removed_edges = removed_edges;
		stats.removed_tets = removed_tets;
		stats.remaining_tets = static_cast<uint32>(p.skeleton_tets_.size());
		log_basic(p, log_tag, " done", " steps=", stats.steps, " removed_faces=", stats.removed_faces,
				  " removed_edges=", stats.removed_edges, " removed_tets=", stats.removed_tets,
				  " remaining_tets=", stats.remaining_tets,
				  '\n');

		return stats;
	}

	struct DenseFiveVertexDetectionResult
	{
		std::vector<std::array<uint32, 5>> independent_exact_k5_regions;
		std::vector<std::array<uint32, 5>> independent_k5_minus_1_regions;
	};

	DenseFiveVertexDetectionResult detect_k5_and_k5_minus_1_cells(PointsParameters& p)
	{
		DenseFiveVertexDetectionResult result;
		if (!p.skeleton_)
			return result;

		auto make_edge_key = [](uint32 a, uint32 b) -> uint64 {
			if (a > b)
				std::swap(a, b);
			return (uint64(a) << 32) | uint64(b);
		};
		std::vector<std::array<uint32, 5>> exact_k5_regions;
		std::vector<std::array<uint32, 5>> k5_minus_1_regions;

		std::unordered_map<uint32, std::unordered_set<uint32>> adjacency;
		adjacency.reserve(nb_cells<NMVertex>(*p.skeleton_));
		std::unordered_map<uint64, uint32> edge_pair_to_id;
		edge_pair_to_id.reserve(nb_cells<NMEdge>(*p.skeleton_) * 2 + 1);

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

		std::set<std::array<uint32, 5>> seen_region_keys;
		for (const auto& tet_entry : p.skeleton_tets_)
		{
			const Tet& tet = tet_entry.second;
			std::array<uint32, 4> base_vertices = {INVALID_INDEX, INVALID_INDEX, INVALID_INDEX, INVALID_INDEX};
			uint32 base_vertex_count = 0;
			for (uint32 face_index = 0; face_index < 4; ++face_index)
			{
				const NMFace face = tet.faces[face_index];
				if (!face.is_valid())
					continue;
				for (NMVertex vertex : incident_vertices(*p.skeleton_, face))
				{
					const uint32 vertex_id = index_of(*p.skeleton_, vertex);
					if (vertex_id == INVALID_INDEX ||
						std::find(base_vertices.begin(), base_vertices.begin() + base_vertex_count, vertex_id) !=
							base_vertices.begin() + base_vertex_count)
						continue;
					if (base_vertex_count < 4)
						base_vertices[base_vertex_count++] = vertex_id;
				}
			}
			if (base_vertex_count != 4)
				continue;
			std::sort(base_vertices.begin(), base_vertices.end());

			std::unordered_map<uint32, uint32> fifth_vertex_connection_count;
			for (uint32 base_vertex : base_vertices)
			{
				const auto adjacency_it = adjacency.find(base_vertex);
				if (adjacency_it == adjacency.end())
					continue;
				for (uint32 neighbor : adjacency_it->second)
				{
					if (!std::binary_search(base_vertices.begin(), base_vertices.end(), neighbor))
						++fifth_vertex_connection_count[neighbor];
				}
			}

			for (const auto& candidate : fifth_vertex_connection_count)
			{
				if (candidate.second < 3)
					continue;

				std::array<uint32, 5> region = {
					base_vertices[0], base_vertices[1], base_vertices[2], base_vertices[3], candidate.first};
				std::sort(region.begin(), region.end());
				if (!seen_region_keys.insert(region).second)
					continue;

				const bool is_exact_k5 = candidate.second == 4;
				if (is_exact_k5)
					exact_k5_regions.push_back(region);
				else
					k5_minus_1_regions.push_back(region);
			}
		}

		std::unordered_map<uint64, uint32> edge_region_count;
		edge_region_count.reserve((exact_k5_regions.size() + k5_minus_1_regions.size()) * 10);
		auto count_region_edges = [&](const std::array<uint32, 5>& region) {
			for (uint32 i = 0; i < 5; ++i)
			{
				for (uint32 j = i + 1; j < 5; ++j)
				{
					const uint64 edge_key = make_edge_key(region[i], region[j]);
					if (edge_pair_to_id.find(edge_key) != edge_pair_to_id.end())
						++edge_region_count[edge_key];
				}
			}
		};
		for (const auto& region : exact_k5_regions)
			count_region_edges(region);
		for (const auto& region : k5_minus_1_regions)
			count_region_edges(region);

		auto append_independent_region = [&](const std::array<uint32, 5>& region, bool is_exact_k5) {
			for (uint32 i = 0; i < 5; ++i)
			{
				for (uint32 j = i + 1; j < 5; ++j)
				{
					const uint64 edge_key = make_edge_key(region[i], region[j]);
					const auto count_it = edge_region_count.find(edge_key);
					if (count_it != edge_region_count.end() && count_it->second > 1)
						return;
				}
			}

			if (is_exact_k5)
				result.independent_exact_k5_regions.push_back(region);
			else
				result.independent_k5_minus_1_regions.push_back(region);
		};
		for (const auto& region : exact_k5_regions)
			append_independent_region(region, true);
		for (const auto& region : k5_minus_1_regions)
			append_independent_region(region, false);

		return result;
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
		PointsParameters& p, const char* log_prefix = "[TopologyFilter]",
		TopologyFixScoreCache* score_cache = nullptr)
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
			std::vector<NMEdge> affected_edges;
			affected_edges.reserve(6);
			for (const NMFace& f_del : faces_to_remove)
			{
				if (!f_del.is_valid())
					continue;
				for (NMEdge e : incident_edges(*p.skeleton_, f_del))
					affected_edges.push_back(e);
			}
			std::unordered_set<uint32> removed_face_ids;
			removed_face_ids.reserve(faces_to_remove.size() * 2 + 1);
			for (const NMFace& f_del : faces_to_remove)
			{
				if (const uint32 idf_del = index_of(*p.skeleton_, f_del); idf_del != INVALID_INDEX)
					removed_face_ids.insert(idf_del);
				remove_face(*p.skeleton_, f_del);
				++stats.removed_faces;
			}
			erase_topology_fix_scores_for_removed_faces_and_orphan_edges(p, score_cache, removed_face_ids, affected_edges);
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

		log_basic(p, log_prefix, " boundary_prepass_passes=", stats.passes, " processed_boundary_tets=",
				  stats.processed_tets, " removed_faces=", stats.removed_faces, " removed_edges=",
				  stats.removed_edges, " removed_tets=", stats.removed_tets, '\n');
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
		p.skeleton_source_sphere_ = get_or_add_attribute<PVertex, NMVertex>(*p.skeleton_, "source_sphere");
		p.incident_tets_ = get_or_add_attribute<std::set<std::size_t>, NMFace>(*p.skeleton_, "incident_tets");
		p.edge_degree_ = get_or_add_attribute<uint32, NMEdge>(*p.skeleton_, "degree");
		p.skeleton_face_color_ = get_or_add_attribute<Vec3, NMFace>(*p.skeleton_, "color");
		p.skeleton_face_component_id_ = get_or_add_attribute<uint32, NMFace>(*p.skeleton_, "face_component_id");
		p.skeleton_face_component_color_ = get_or_add_attribute<Vec3, NMFace>(*p.skeleton_, "face_component_color");
		p.skeleton_edge_color_ = get_or_add_attribute<Vec3, NMEdge>(*p.skeleton_, "color");
		p.skeleton_edge_non_manifold_color_ =
			get_or_add_attribute<Vec3, NMEdge>(*p.skeleton_, "non_manifold_edge_color");
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

	struct IndependentDenseFiveVertexRepairSeed
	{
		std::array<PVertex, 5> spheres;
	};

	struct IndependentDenseFiveVertexRepairStats
	{
		uint32 snapshotted_regions = 0;
		uint32 regions_with_isolated_nm_edges = 0;
		uint32 repaired_regions = 0;
		uint32 added_spheres = 0;
		uint32 added_faces = 0;
		uint32 skipped_missing_vertices = 0;
		uint32 skipped_adjacent_nm_edges = 0;
		uint32 skipped_unsafe_external_faces = 0;
		uint32 skipped_degenerate_geometry = 0;
	};

	bool is_live_sphere_vertex(PointsParameters& p, const PVertex& sphere, uint32& sphere_index)
	{
		sphere_index = INVALID_INDEX;
		if (!p.spheres_ || !sphere.is_valid() || sphere.dart_.index_ >= p.spheres_->darts_.maximum_index())
			return false;
		sphere_index = index_of(*p.spheres_, sphere);
		return sphere_index != INVALID_INDEX;
	}

	bool is_live_skeleton_vertex(PointsParameters& p, const NMVertex& vertex, uint32& vertex_index)
	{
		vertex_index = INVALID_INDEX;
		if (!p.skeleton_ || !vertex.is_valid())
			return false;
		vertex_index = index_of(*p.skeleton_, vertex);
		return vertex_index != INVALID_INDEX;
	}

	std::vector<IndependentDenseFiveVertexRepairSeed> collect_independent_dense_five_vertex_repair_seeds(
		PointsParameters& p, const DenseFiveVertexDetectionResult& detection)
	{
		std::vector<IndependentDenseFiveVertexRepairSeed> seeds;
		seeds.reserve(detection.independent_exact_k5_regions.size() +
					  detection.independent_k5_minus_1_regions.size());
		if (!p.skeleton_ || !p.spheres_ || !p.skeleton_source_sphere_)
			return seeds;

		auto append_regions = [&](const std::vector<std::array<uint32, 5>>& regions) {
			for (const auto& region : regions)
			{
				IndependentDenseFiveVertexRepairSeed seed;
				bool valid = true;
				for (uint32 i = 0; i < 5; ++i)
				{
					const NMVertex vertex = of_index<NMVertex>(*p.skeleton_, region[i]);
					uint32 current_vertex_index = INVALID_INDEX;
					if (!is_live_skeleton_vertex(p, vertex, current_vertex_index) || current_vertex_index != region[i])
					{
						valid = false;
						break;
					}
					const PVertex sphere = value<PVertex>(*p.skeleton_, p.skeleton_source_sphere_, vertex);
					uint32 sphere_index = INVALID_INDEX;
					if (!is_live_sphere_vertex(p, sphere, sphere_index))
					{
						valid = false;
						break;
					}
					seed.spheres[i] = sphere;
				}
				if (valid)
					seeds.push_back(seed);
			}
		};
		append_regions(detection.independent_exact_k5_regions);
		append_regions(detection.independent_k5_minus_1_regions);
		return seeds;
	}

	bool find_skeleton_edge_between(PointsParameters& p, const NMVertex& a, const NMVertex& b, NMEdge& out_edge)
	{
		out_edge = NMEdge();
		if (!p.skeleton_ || !a.is_valid() || !b.is_valid() || a == b)
			return false;
		for (NMEdge edge : incident_edges(*p.skeleton_, a))
		{
			if (!edge.is_valid())
				continue;
			const std::vector<NMVertex> edge_vertices = incident_vertices(*p.skeleton_, edge);
			if (edge_vertices.size() == 2 &&
				((edge_vertices[0] == a && edge_vertices[1] == b) ||
				 (edge_vertices[0] == b && edge_vertices[1] == a)))
			{
				out_edge = edge;
				return true;
			}
		}
		return false;
	}

	bool face_is_inside_vertex_set(PointsParameters& p, const NMFace& face,
								   const std::unordered_set<uint32>& vertex_ids)
	{
		if (!p.skeleton_ || !face.is_valid())
			return false;
		const std::vector<NMVertex> vertices = incident_vertices(*p.skeleton_, face);
		if (vertices.size() != 3)
			return false;
		for (NMVertex vertex : vertices)
		{
			const uint32 vertex_id = index_of(*p.skeleton_, vertex);
			if (vertex_id == INVALID_INDEX || vertex_ids.find(vertex_id) == vertex_ids.end())
				return false;
		}
		return true;
	}

	enum class IndependentDenseFiveVertexRepairResult
	{
		Repaired,
		MissingVertices,
		NoNonManifoldEdges,
		AdjacentNonManifoldEdges,
		UnsafeExternalFaces,
		DegenerateGeometry
	};

	IndependentDenseFiveVertexRepairResult repair_independent_dense_five_vertex_region(
		PointsParameters& p, const IndependentDenseFiveVertexRepairSeed& seed, uint32& out_added_faces)
	{
		out_added_faces = 0;
		if (!p.spheres_ || !p.spheres_position_ || !p.spheres_radius_ || !p.spheres_neighbor_clusters_ ||
			!p.spheres_skeleton_vertex_ || !p.skeleton_ || !p.skeleton_position_ || !p.skeleton_radius_ ||
			!p.skeleton_source_sphere_ || !p.incident_tets_)
			return IndependentDenseFiveVertexRepairResult::MissingVertices;

		std::array<uint32, 5> sphere_ids;
		std::array<NMVertex, 5> vertices;
		std::array<uint32, 5> vertex_ids;
		std::array<Vec3, 5> positions;
		Scalar average_radius = Scalar(0);
		Vec3 centroid = Vec3::Zero();
		std::unordered_set<uint32> region_vertex_ids;
		region_vertex_ids.reserve(8);
		for (uint32 i = 0; i < 5; ++i)
		{
			if (!is_live_sphere_vertex(p, seed.spheres[i], sphere_ids[i]))
				return IndependentDenseFiveVertexRepairResult::MissingVertices;
			vertices[i] = (*p.spheres_skeleton_vertex_)[sphere_ids[i]];
			if (!is_live_skeleton_vertex(p, vertices[i], vertex_ids[i]))
				return IndependentDenseFiveVertexRepairResult::MissingVertices;
			if (value<PVertex>(*p.skeleton_, p.skeleton_source_sphere_, vertices[i]) != seed.spheres[i])
				return IndependentDenseFiveVertexRepairResult::MissingVertices;
			positions[i] = (*p.spheres_position_)[sphere_ids[i]];
			centroid += positions[i];
			average_radius += (*p.spheres_radius_)[sphere_ids[i]];
			region_vertex_ids.insert(vertex_ids[i]);
		}
		if (region_vertex_ids.size() != 5)
			return IndependentDenseFiveVertexRepairResult::MissingVertices;
		centroid /= Scalar(5);
		average_radius /= Scalar(5);
		if (!centroid.allFinite() || !std::isfinite(static_cast<double>(average_radius)) || average_radius <= Scalar(0))
			return IndependentDenseFiveVertexRepairResult::DegenerateGeometry;

		auto make_edge_key = [](uint32 a, uint32 b) -> uint64 {
			if (a > b)
				std::swap(a, b);
			return (uint64(a) << 32) | uint64(b);
		};
		std::unordered_map<uint64, NMEdge> region_edges;
		region_edges.reserve(16);
		std::unordered_set<uint32> region_nm_edge_ids;
		region_nm_edge_ids.reserve(16);
		for (uint32 i = 0; i < 5; ++i)
		{
			for (uint32 j = i + 1; j < 5; ++j)
			{
				NMEdge edge;
				if (!find_skeleton_edge_between(p, vertices[i], vertices[j], edge))
					continue;
				region_edges[make_edge_key(vertex_ids[i], vertex_ids[j])] = edge;
				if (incident_faces(*p.skeleton_, edge).size() > 2)
				{
					const uint32 edge_id = index_of(*p.skeleton_, edge);
					if (edge_id != INVALID_INDEX)
						region_nm_edge_ids.insert(edge_id);
				}
			}
		}
		if (region_nm_edge_ids.empty())
			return IndependentDenseFiveVertexRepairResult::NoNonManifoldEdges;

		for (uint32 nm_edge_id : region_nm_edge_ids)
		{
			const NMEdge nm_edge = of_index<NMEdge>(*p.skeleton_, nm_edge_id);
			if (!nm_edge.is_valid())
				continue;
			for (NMVertex endpoint : incident_vertices(*p.skeleton_, nm_edge))
			{
				for (NMEdge adjacent_edge : incident_edges(*p.skeleton_, endpoint))
				{
					if (!adjacent_edge.is_valid() || incident_faces(*p.skeleton_, adjacent_edge).size() <= 2)
						continue;
					const uint32 adjacent_edge_id = index_of(*p.skeleton_, adjacent_edge);
					if (adjacent_edge_id != INVALID_INDEX &&
						region_nm_edge_ids.find(adjacent_edge_id) == region_nm_edge_ids.end())
						return IndependentDenseFiveVertexRepairResult::AdjacentNonManifoldEdges;
				}
			}
		}

		Eigen::Matrix<Scalar, 3, 3> covariance = Eigen::Matrix<Scalar, 3, 3>::Zero();
		for (const Vec3& position : positions)
		{
			const Vec3 delta = position - centroid;
			covariance += delta * delta.transpose();
		}
		Eigen::SelfAdjointEigenSolver<Eigen::Matrix<Scalar, 3, 3>> eigen_solver(covariance);
		if (eigen_solver.info() != Eigen::Success ||
			eigen_solver.eigenvalues()[1] <= std::numeric_limits<Scalar>::epsilon())
			return IndependentDenseFiveVertexRepairResult::DegenerateGeometry;
		const Vec3 axis_u = eigen_solver.eigenvectors().col(2);
		const Vec3 axis_v = eigen_solver.eigenvectors().col(1);
		std::array<uint32, 5> cyclic_order = {0, 1, 2, 3, 4};
		std::sort(cyclic_order.begin(), cyclic_order.end(), [&](uint32 lhs, uint32 rhs) {
			const Vec3 lhs_delta = positions[lhs] - centroid;
			const Vec3 rhs_delta = positions[rhs] - centroid;
			const Scalar lhs_angle = std::atan2(lhs_delta.dot(axis_v), lhs_delta.dot(axis_u));
			const Scalar rhs_angle = std::atan2(rhs_delta.dot(axis_v), rhs_delta.dot(axis_u));
			return lhs_angle < rhs_angle;
		});
		Scalar max_radial_squared = Scalar(0);
		for (const Vec3& position : positions)
			max_radial_squared = std::max(max_radial_squared, (position - centroid).squaredNorm());
		const Scalar fan_area_epsilon =
			std::max(std::numeric_limits<Scalar>::epsilon(), max_radial_squared * Scalar(1e-10));
		for (uint32 i = 0; i < 5; ++i)
		{
			const Vec3 radial_a = positions[cyclic_order[i]] - centroid;
			const Vec3 radial_b = positions[cyclic_order[(i + 1) % 5]] - centroid;
			if (radial_a.cross(radial_b).norm() <= fan_area_epsilon)
				return IndependentDenseFiveVertexRepairResult::DegenerateGeometry;
		}

		std::unordered_set<uint64> cycle_edge_keys;
		cycle_edge_keys.reserve(8);
		for (uint32 i = 0; i < 5; ++i)
		{
			const uint32 a = cyclic_order[i];
			const uint32 b = cyclic_order[(i + 1) % 5];
			cycle_edge_keys.insert(make_edge_key(vertex_ids[a], vertex_ids[b]));
		}

		std::unordered_set<uint32> internal_face_ids;
		internal_face_ids.reserve(16);
		foreach_cell(*p.skeleton_, [&](NMFace face) -> bool {
			if (face_is_inside_vertex_set(p, face, region_vertex_ids))
			{
				const uint32 face_id = index_of(*p.skeleton_, face);
				if (face_id != INVALID_INDEX)
					internal_face_ids.insert(face_id);
			}
			return true;
		});
		if (internal_face_ids.empty())
			return IndependentDenseFiveVertexRepairResult::UnsafeExternalFaces;

		for (const auto& edge_entry : region_edges)
		{
			uint32 external_face_count = 0;
			for (NMFace face : incident_faces(*p.skeleton_, edge_entry.second))
			{
				const uint32 face_id = index_of(*p.skeleton_, face);
				if (face_id == INVALID_INDEX || internal_face_ids.find(face_id) == internal_face_ids.end())
					++external_face_count;
			}
			if (cycle_edge_keys.find(edge_entry.first) == cycle_edge_keys.end())
			{
				if (external_face_count != 0)
					return IndependentDenseFiveVertexRepairResult::UnsafeExternalFaces;
			}
			else if (external_face_count > 1)
			{
				return IndependentDenseFiveVertexRepairResult::UnsafeExternalFaces;
			}
		}

		SkeletonFaceDeletionStats deletion_stats;
		SkeletonFaceDeletionOptions deletion_options;
		deletion_options.remove_incident_orphan_edges_immediately = false;
		deletion_options.remove_global_orphan_elements_after_batch = false;
		if (!delete_skeleton_face_id_set(p, internal_face_ids, deletion_stats, deletion_options))
			return IndependentDenseFiveVertexRepairResult::UnsafeExternalFaces;

		for (const auto& edge_entry : region_edges)
		{
			if (cycle_edge_keys.find(edge_entry.first) != cycle_edge_keys.end())
				continue;
			if (edge_entry.second.is_valid() && incident_faces(*p.skeleton_, edge_entry.second).empty())
				remove_edge(*p.skeleton_, edge_entry.second);
		}

		for (uint32 i = 0; i < 5; ++i)
		{
			for (uint32 j = i + 1; j < 5; ++j)
			{
				const uint64 edge_key = make_edge_key(vertex_ids[i], vertex_ids[j]);
				const bool is_cycle_edge = cycle_edge_keys.find(edge_key) != cycle_edge_keys.end();
				if (is_cycle_edge)
				{
					(*p.spheres_neighbor_clusters_)[sphere_ids[i]].insert(seed.spheres[j]);
					(*p.spheres_neighbor_clusters_)[sphere_ids[j]].insert(seed.spheres[i]);
				}
				else
				{
					(*p.spheres_neighbor_clusters_)[sphere_ids[i]].erase(seed.spheres[j]);
					(*p.spheres_neighbor_clusters_)[sphere_ids[j]].erase(seed.spheres[i]);
				}
			}
		}

		const PVertex center_sphere = add_vertex(*p.spheres_);
		const uint32 center_sphere_id = index_of(*p.spheres_, center_sphere);
		(*p.spheres_position_)[center_sphere_id] = centroid;
		(*p.spheres_radius_)[center_sphere_id] = average_radius;
		if (p.spheres_color_)
			(*p.spheres_color_)[center_sphere_id] = Vec4(0.95, 0.25, 0.15, 1.0);
		if (p.spheres_cluster_)
			(*p.spheres_cluster_)[center_sphere_id].clear();
		if (p.spheres_cluster_area_)
			(*p.spheres_cluster_area_)[center_sphere_id] = Scalar(0);
		if (p.spheres_cluster_color_)
			(*p.spheres_cluster_color_)[center_sphere_id] = Vec4(0.95, 0.25, 0.15, 1.0);
		if (p.spheres_parent_)
			(*p.spheres_parent_)[center_sphere_id] = PVertex();
		if (p.spheres_error_)
			(*p.spheres_error_)[center_sphere_id] = Scalar(0);
		if (p.spheres_error_not_normalized_)
			(*p.spheres_error_not_normalized_)[center_sphere_id] = Scalar(0);
		(*p.spheres_neighbor_clusters_)[center_sphere_id].clear();
		++p.nb_spheres_;

		const NMVertex center_vertex = add_vertex(*p.skeleton_);
		value<Vec3>(*p.skeleton_, p.skeleton_position_, center_vertex) = centroid;
		value<Scalar>(*p.skeleton_, p.skeleton_radius_, center_vertex) = average_radius;
		value<PVertex>(*p.skeleton_, p.skeleton_source_sphere_, center_vertex) = center_sphere;
		(*p.spheres_skeleton_vertex_)[center_sphere_id] = center_vertex;

		std::array<NMEdge, 5> cycle_edges;
		std::array<NMEdge, 5> spoke_edges;
		for (uint32 i = 0; i < 5; ++i)
		{
			const uint32 a = cyclic_order[i];
			const uint32 b = cyclic_order[(i + 1) % 5];
			if (!find_skeleton_edge_between(p, vertices[a], vertices[b], cycle_edges[i]))
				cycle_edges[i] = add_edge(*p.skeleton_, vertices[a], vertices[b]);
			spoke_edges[i] = add_edge(*p.skeleton_, center_vertex, vertices[a]);
			(*p.spheres_neighbor_clusters_)[center_sphere_id].insert(seed.spheres[a]);
			(*p.spheres_neighbor_clusters_)[sphere_ids[a]].insert(center_sphere);
		}

		for (uint32 i = 0; i < 5; ++i)
		{
			std::vector<NMEdge> face_edges = {cycle_edges[i], spoke_edges[i], spoke_edges[(i + 1) % 5]};
			const NMFace new_face = add_face(*p.skeleton_, face_edges);
			value<std::set<std::size_t>>(*p.skeleton_, p.incident_tets_, new_face).clear();
			if (p.skeleton_face_color_)
				value<Vec3>(*p.skeleton_, p.skeleton_face_color_, new_face) = Vec3(0.0, 0.0, 0.0);
			++out_added_faces;
		}

		return IndependentDenseFiveVertexRepairResult::Repaired;
	}

	IndependentDenseFiveVertexRepairStats repair_isolated_independent_dense_five_vertex_regions(
		PointsParameters& p, const std::vector<IndependentDenseFiveVertexRepairSeed>& seeds)
	{
		IndependentDenseFiveVertexRepairStats stats;
		stats.snapshotted_regions = static_cast<uint32>(seeds.size());
		std::vector<IndependentDenseFiveVertexRepairSeed> live_seeds;
		live_seeds.reserve(seeds.size());
		for (const auto& seed : seeds)
		{
			bool all_vertices_live = true;
			std::unordered_set<uint32> unique_skeleton_vertex_ids;
			unique_skeleton_vertex_ids.reserve(8);
			for (const PVertex& sphere : seed.spheres)
			{
				uint32 sphere_id = INVALID_INDEX;
				if (!is_live_sphere_vertex(p, sphere, sphere_id))
				{
					all_vertices_live = false;
					break;
				}
				const NMVertex vertex = (*p.spheres_skeleton_vertex_)[sphere_id];
				uint32 vertex_id = INVALID_INDEX;
				if (!is_live_skeleton_vertex(p, vertex, vertex_id) ||
					value<PVertex>(*p.skeleton_, p.skeleton_source_sphere_, vertex) != sphere)
				{
					all_vertices_live = false;
					break;
				}
				unique_skeleton_vertex_ids.insert(vertex_id);
			}
			if (all_vertices_live && unique_skeleton_vertex_ids.size() == 5)
				live_seeds.push_back(seed);
			else
				++stats.skipped_missing_vertices;
		}

		for (const auto& seed : live_seeds)
		{
			uint32 added_faces = 0;
			const IndependentDenseFiveVertexRepairResult result =
				repair_independent_dense_five_vertex_region(p, seed, added_faces);
			switch (result)
			{
				case IndependentDenseFiveVertexRepairResult::Repaired:
					++stats.regions_with_isolated_nm_edges;
					++stats.repaired_regions;
					++stats.added_spheres;
					stats.added_faces += added_faces;
					break;
				case IndependentDenseFiveVertexRepairResult::NoNonManifoldEdges:
					break;
				case IndependentDenseFiveVertexRepairResult::MissingVertices:
					++stats.skipped_missing_vertices;
					break;
				case IndependentDenseFiveVertexRepairResult::AdjacentNonManifoldEdges:
					++stats.skipped_adjacent_nm_edges;
					break;
				case IndependentDenseFiveVertexRepairResult::UnsafeExternalFaces:
					++stats.regions_with_isolated_nm_edges;
					++stats.skipped_unsafe_external_faces;
					break;
				case IndependentDenseFiveVertexRepairResult::DegenerateGeometry:
					++stats.regions_with_isolated_nm_edges;
					++stats.skipped_degenerate_geometry;
					break;
			}
		}

		if (stats.repaired_regions > 0)
		{
			compute_edge_degree(p);
			non_manifold_provider_->emit_connectivity_changed(*p.skeleton_);
			non_manifold_provider_->emit_attribute_changed(*p.skeleton_, p.skeleton_position_.get());
			non_manifold_provider_->emit_attribute_changed(*p.skeleton_, p.skeleton_radius_.get());
			points_provider_->emit_connectivity_changed(*p.spheres_);
			points_provider_->emit_attribute_changed(*p.spheres_, p.spheres_position_.get());
			points_provider_->emit_attribute_changed(*p.spheres_, p.spheres_radius_.get());
			if (p.spheres_color_)
				points_provider_->emit_attribute_changed(*p.spheres_, p.spheres_color_.get());
		}
		return stats;
	}

	void run_topology_fix_pipeline(PointsParameters& p)
	{
		if (!p.skeleton_ || !p.incident_tets_)
		{
			log_error(p, "[TopologyFull] requires a built skeleton.", '\n');
			return;
		}

		const DenseFiveVertexDetectionResult dense_five_vertex_detection = detect_k5_and_k5_minus_1_cells(p);
		const std::vector<IndependentDenseFiveVertexRepairSeed> independent_dense_repair_seeds =
			collect_independent_dense_five_vertex_repair_seeds(p, dense_five_vertex_detection);

		log_basic(p, "[TopologyFull] start initial_tets=",
				  p.skeleton_tets_.size(), " independent_k5=",
				  dense_five_vertex_detection.independent_exact_k5_regions.size(), " independent_k5_minus_1=",
				  dense_five_vertex_detection.independent_k5_minus_1_regions.size(), " repair_seeds=",
				  independent_dense_repair_seeds.size(), '\n');
		TopologyFixScoreCache score_cache;
		if (!initialize_topology_fix_score_cache(p, score_cache, "[TopologyFull]"))
			return;

		uint32 round = 0;
		while (!p.skeleton_tets_.empty())
		{
			++round;
			const uint32 round_start_tets = static_cast<uint32>(p.skeleton_tets_.size());
			log_basic(p, "[TopologyFull] round=", round, " start_tets=", round_start_tets, '\n');

			const BoundaryTetPrepassStats boundary_stats =
				run_boundary_tet_face_deletion(p, "[TopologyFull]", &score_cache);

			EdgeTetModeRunStats simple_stats;
			if (!p.skeleton_tets_.empty())
				simple_stats = run_edge_score_tet_mode_topology_fix(p, EdgeTetDeleteMode::SimpleTet, score_cache);

			EdgeTetModeRunStats nonsimple_stats;
			if (!p.skeleton_tets_.empty())
				nonsimple_stats = run_edge_score_tet_mode_topology_fix(p, EdgeTetDeleteMode::NonSimpleTet, score_cache);

			const uint32 round_end_tets = static_cast<uint32>(p.skeleton_tets_.size());
			const uint32 round_removed_tets = round_start_tets - round_end_tets;
			log_basic(p, "[TopologyFull] round=", round, " done boundary_removed_tets=", boundary_stats.removed_tets,
				  " simple_removed_tets=", simple_stats.removed_tets, " nonsimple_removed_tets=",
				  nonsimple_stats.removed_tets, " total_removed_tets=", round_removed_tets, " remaining_tets=",
				  round_end_tets, '\n');

			if (round_removed_tets == 0)
			{
				log_basic(p, "[TopologyFull] round=", round, " stop reason=no_progress", '\n');
				break;
			}
		}

		const IndependentDenseFiveVertexRepairStats dense_repair_stats =
			repair_isolated_independent_dense_five_vertex_regions(p, independent_dense_repair_seeds);
		refresh_skeleton_topology_colors(p);
		log_basic(p, "[TopologyFull] done remaining_tets=", p.skeleton_tets_.size(), " repair_snapshots=",
				  dense_repair_stats.snapshotted_regions, " isolated_nm_regions=",
				  dense_repair_stats.regions_with_isolated_nm_edges, " repaired_regions=",
				  dense_repair_stats.repaired_regions, " added_spheres=", dense_repair_stats.added_spheres,
				  " added_faces=", dense_repair_stats.added_faces, " skipped_missing_vertices=",
				  dense_repair_stats.skipped_missing_vertices, " skipped_adjacent_nm_edges=",
				  dense_repair_stats.skipped_adjacent_nm_edges, " skipped_unsafe_external_faces=",
				  dense_repair_stats.skipped_unsafe_external_faces, " skipped_degenerate_geometry=",
				  dense_repair_stats.skipped_degenerate_geometry, '\n');
	}

protected:
	void update_render_data(PointsParameters& p, bool non_blocking_running_lock = false, bool full_refresh = true,
						bool rebuild_skeleton = false)
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
		points_provider_->emit_attribute_changed(*p.spheres_, p.spheres_color_.get());

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

		if (rebuild_skeleton)
			compute_skeleton(p);
		non_manifold_provider_->emit_connectivity_changed(*p.skeleton_);
		non_manifold_provider_->emit_attribute_changed(*p.skeleton_, p.skeleton_position_.get());
	}

	void request_linked_views_update()
	{
		for (View* v : linked_views_)
			v->request_update();
	}

	void set_post_init_sphere_render_state(PointsParameters& p)
	{
		if (p.points_ && pcr_ && p.position_)
		{
			for (View* v : linked_views_)
			{
				pcr_->set_vertex_position(*v, *p.points_, p.position_);
				pcr_->set_render_vertices(*v, *p.points_, false);
			}
		}

		if (p.skeleton_ && skeleton_render_ && p.skeleton_position_)
		{
			for (View* v : linked_views_)
			{
				skeleton_render_->set_vertex_position(*v, *p.skeleton_, p.skeleton_position_);
				skeleton_render_->set_render_vertices(*v, *p.skeleton_, true);
				skeleton_render_->set_render_edges(*v, *p.skeleton_, true);
				skeleton_render_->set_render_faces(*v, *p.skeleton_, true);
			}
		}
	}

	void start_spheres_update(PointsParameters& p)
	{
		p.running_ = true;
		p.stopping_ = false;
		p.iteration_count_ = 0;
		p.total_error_diff_ = 0.0;
		p.last_total_error_ = std::numeric_limits<Scalar>::max();
		p.pending_full_refresh_after_stop_ = false;

		launch_thread([this, &p]() {
			constexpr Scalar convergence_eps = Scalar(1e-10);
			constexpr uint32 max_post_convergence_iterations = 10;
			constexpr uint32 max_iterations = 150;
			bool convergence_reached = false;
			uint32 post_convergence_iterations = 0;
			auto start = std::chrono::high_resolution_clock::now();
			while (p.iteration_count_ < max_iterations)
			{
				{
					std::lock_guard<std::mutex> lock(p.mutex_);
					if (!p.stopping_)
					{
						log_basic(p, "Start Sphere update", '\n');
						update_spheres(p);
						++p.iteration_count_;
					}
				}
				if (p.stopping_)
					break;
				if (p.slow_down_)
					std::this_thread::sleep_for(std::chrono::microseconds(1000000 / p.update_rate_));
				else
					std::this_thread::yield();

				if (p.total_error_diff_ < convergence_eps)
				{
					if (!convergence_reached)
					{
						convergence_reached = true;
						post_convergence_iterations = 0;
						log_basic(p, "Auto stop: error converged (Diff < ", convergence_eps,
								  "), start post-convergence countdown (", max_post_convergence_iterations, ").", '\n');
					}
					else
					{
						++post_convergence_iterations;
					}
					if (post_convergence_iterations >= max_post_convergence_iterations)
					{
						log_basic(p, "Auto stop: reached max post-convergence iterations (",
								  max_post_convergence_iterations, ").", '\n');
						p.stopping_ = true;
					}
				}
				else if (convergence_reached)
				{
					convergence_reached = false;
					post_convergence_iterations = 0;
				}

				if (p.iteration_count_ >= max_iterations)
				{
					log_basic(p, "Stop: reached max iterations (", max_iterations, ").", '\n');
					p.stopping_ = true;
				}
				log_basic(p, "Iteration: ", p.iteration_count_, " | Spheres: ", p.nb_spheres_, " | Error: ",
						  p.total_error_, " | Diff: ", p.total_error_diff_, '\n');
				if (p.stopping_)
					break;
			}
			p.stopping_ = false;
			p.running_ = false;
			p.pending_full_refresh_after_stop_ = true;
			auto end = std::chrono::high_resolution_clock::now();
			log_basic(p, "Sphere optimizations time: ", std::chrono::duration<Scalar>(end - start).count(), "s", '\n');
			log_basic(p, "Nb iterations: ", p.iteration_count_, '\n');
		});

		app_.start_timer(100, [&]() -> bool { return !p.running_ && !p.pending_full_refresh_after_stop_; });
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

			picked_sphere_ = PVertex();
			Vec3 picked_sphere_center = Vec3::Zero();
			bool has_picked_sphere_center = false;
			foreach_cell(*p.spheres_, [&](PVertex v) -> bool {
				if (!has_picked_sphere_center)
				{
					picked_sphere_ = v;
					picked_sphere_center = (*p.spheres_position_)[index_of(*p.spheres_, picked_sphere_)];
					has_picked_sphere_center = true;
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
					update_render_data(p, false, true, true);
			}
			else if (key_code == GLFW_KEY_D && picked_sphere_.is_valid())
			{
				std::lock_guard<std::mutex> lock(p.mutex_);
				remove_sphere(p, picked_sphere_);
				compute_spheres_error(p);
				picked_sphere_ = PVertex();
				if (!p.running_)
					update_render_data(p, false, true, true);
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
		int output_verbosity = output_verbosity_index(p.output_verbosity_);
		if (ImGui::Combo("Output", &output_verbosity, "Mute\0Normal\0Verbose\0"))
			p.output_verbosity_ = output_verbosity_from_index(output_verbosity);

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
			ImGui::InputFloat("Sample Radius", &p.sample_radius_, 0.001f, 0.01f, "%.4f");
			ImGui::InputInt("Eval Batch Size", &p.batch_size_, 256, 1024);
			ImGui::InputInt("Ray Batch Size", &p.ray_sampler_batch_size_, 256, 2048);
			ImGui::InputInt("Max Iterations", &p.udf_max_iterations_, 1000, 8000);
			ImGui::InputFloat("Tolerance", &p.tol_, 0.0f, 0.0f, "%.6f");
			ImGui::InputInt("KNN K", &p.knn_k_, 1, 5);
			ImGui::Checkbox("Recompute Normals In Fitting Data", &p.recompute_sample_normals_after_sampling_);
			if (ImGui::Button("Apply Sampling Filtering"))
			{
				std::lock_guard<std::mutex> lock(p.mutex_);
				apply_sampling_preprocess_filtering(p);
			}
			if (ImGui::Button("Sample UDF"))
			{
				load_alpha_samples_to_mesh(p);
				p.fitting_data_computed_ = false;
			}
			ImGui::SameLine();
			if (ImGui::Button("Clear Samples"))
			{
				if (p.samples_mesh_)
					points_provider_->clear_mesh(*p.samples_mesh_);
				p.fitting_data_computed_ = false;
					}

			if (p.samples_mesh_)
				ImGui::Text("Number of samples: %zu", nb_cells<PVertex>(*p.samples_mesh_));
		}

		bool has_samples = p.samples_mesh_ && nb_cells<PVertex>(*p.samples_mesh_) > 0;

		if (!has_samples)
		{
			ImGui::TextColored(ImVec4(1, 1, 0, 1), "Please sample points first.");
		}
		else
		{
			ImGui::Separator();
			ImGui::Checkbox("Enable MAFlipPrune", &p.ma_flip_prune_enabled_);
			float ma_flip_prune_alpha_factor = static_cast<float>(p.ma_flip_prune_alpha_factor_);
			if (ImGui::SliderFloat("MAFlipPrune MF Alpha Factor", &ma_flip_prune_alpha_factor, 1.0f, 5.0f, "%.2f"))
				p.ma_flip_prune_alpha_factor_ = Scalar(ma_flip_prune_alpha_factor);
			if (ImGui::Button("Compute Fitting Data"))
			{
				p.fitting_data_computed_ = false;
				compute_fitting_data(p);
				update_render_data(p);
			}
			const bool sphere_fit_ready = p.fitting_data_computed_;
			if (sphere_fit_ready)
			{
				// Sphere Fitting
				if (ImGui::CollapsingHeader("Sphere Fitting", ImGuiTreeNodeFlags_DefaultOpen))
				{
					ImGui::SliderFloat("Init dilation constant", &p.init_dilation_constant_, 0.001f, 0.01f, "%.4f");
					if (ImGui::Button("Init spheres"))
					{
						std::lock_guard<std::mutex> lock(p.mutex_);
						init_spheres(p);
						update_render_data(p, false, true, true);
					}
					ImGui::SliderFloat("Update lambda", &p.sqem_update_lambda_line_plane_, 0.0f, 4.0f, "%.6f");
					if (ImGui::Button("Update spheres"))
					{
						if (!p.running_)
						{
							std::lock_guard<std::mutex> lock(p.mutex_);
							update_spheres(p);
							update_render_data(p, false, true, true);
						}
					}


					if (ImGui::Button("Build Skeleton"))
					{
						if (!p.running_)
						{
							std::lock_guard<std::mutex> lock(p.mutex_);
							compute_skeleton(p);
							update_render_data(p, false, true, false);
						}
					}
					const bool samples_mesh_export_available =
						p.samples_mesh_ && p.samples_position_ && nb_cells<PVertex>(*p.samples_mesh_) > 0;
					const bool samples_mesh_normal_color_export_available =
						p.samples_mesh_ && p.samples_position_ && p.samples_normal_color_ &&
						nb_cells<PVertex>(*p.samples_mesh_) > 0;
					const bool samples_spheres_export_available =
						p.spheres_ && p.spheres_position_ && p.spheres_radius_ && p.spheres_cluster_color_ &&
						nb_cells<PVertex>(*p.spheres_) > 0;
					const bool skeleton_export_available =
						p.skeleton_ && p.skeleton_position_ && nb_cells<NMVertex>(*p.skeleton_) > 0;
					const bool all_selected_exports_available =
						(!p.export_samples_mesh_selected_ || samples_mesh_export_available) &&
						(!p.export_samples_mesh_normal_color_selected_ || samples_mesh_normal_color_export_available) &&
						(!p.export_samples_spheres_selected_ || samples_spheres_export_available) &&
						(!p.export_skeleton_selected_ || skeleton_export_available);
					const bool any_export_selected =
						p.export_samples_mesh_selected_ || p.export_samples_mesh_normal_color_selected_ ||
						p.export_samples_spheres_selected_ || p.export_skeleton_selected_;
					const bool can_export_training_bundle =
						!p.running_ && any_export_selected && all_selected_exports_available;
					ImGui::Checkbox("Save samples_mesh", &p.export_samples_mesh_selected_);
					ImGui::SameLine();
					ImGui::TextDisabled(samples_mesh_export_available ? "ready" : "missing");
					ImGui::Checkbox("Save sample mesh + normal color", &p.export_samples_mesh_normal_color_selected_);
					ImGui::SameLine();
					ImGui::TextDisabled(samples_mesh_normal_color_export_available ? "ready" : "missing");
					ImGui::Checkbox("Save samples_spheres", &p.export_samples_spheres_selected_);
					ImGui::SameLine();
					ImGui::TextDisabled(samples_spheres_export_available ? "ready" : "missing");
					ImGui::Checkbox("Save skeleton", &p.export_skeleton_selected_);
					ImGui::SameLine();
					ImGui::TextDisabled(skeleton_export_available ? "ready" : "missing");
					if (!can_export_training_bundle)
						ImGui::BeginDisabled();
					if (ImGui::Button("Export selected PLY"))
					{
						const std::string export_directory =
							pfd::select_folder("Select export folder", default_training_export_directory(p).string())
								.result();
						if (!export_directory.empty())
						{
							std::lock_guard<std::mutex> lock(p.mutex_);
							export_training_ply_bundle(p, std::filesystem::path(export_directory));
						}
					}
					if (!can_export_training_bundle)
						ImGui::EndDisabled();
					ImGui::SameLine();
					ImGui::TextDisabled("Select one or more ready targets.");
					if (ImGui::Button("Face components (UF)"))
					{
						if (!p.running_)
						{
							std::lock_guard<std::mutex> lock(p.mutex_);
							prune_fully_non_manifold_triangles(p, "[FaceComponentsPreDelete]");
							if (compute_skeleton_face_components_union_find(p) &&
								colorize_skeleton_face_components(p) && non_manifold_provider_ && p.skeleton_ &&
								p.skeleton_face_component_color_)
							{
								non_manifold_provider_->emit_attribute_changed(*p.skeleton_,
																p.skeleton_face_component_color_.get());
							}
						}
					}
					ImGui::SameLine();
					if (ImGui::Button("Topology fix full pipeline"))
					{
						if (!p.running_)
						{
							std::lock_guard<std::mutex> lock(p.mutex_);
							run_topology_fix_pipeline(p);
						}
					}
					ImGui::SameLine();

					ImGui::Checkbox("Slow down", &p.slow_down_);
					if (p.slow_down_)
						ImGui::SliderInt("Update rate", (int*)&p.update_rate_, 1, 100);
					ImGui::Checkbox("Preview during update", &p.preview_render_during_sphere_update_);
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
					ImGui::Separator();

					if (ImGui::Checkbox("Error as color", &p.error_as_spheres_color_))
					{
						std::lock_guard<std::mutex> lock(p.mutex_);
						if (!p.running_)
							update_render_data(p, false, false);
					}
					if (ImGui::SliderFloat("Transparency", &p.spheres_transparency_, 0.0f, 1.0f))
					{
						std::lock_guard<std::mutex> lock(p.mutex_);
						if (!p.running_)
							update_render_data(p, false, true, false);
					}

					if (ImGui::Button("Split max error sphere"))
					{
						std::lock_guard<std::mutex> lock(p.mutex_);
						split_sphere(p, p.max_error_sphere_);
						compute_clusters(p);
						compute_spheres_error(p);
						if (!p.running_)
							update_render_data(p, false, true, true);
					}

					ImGui::Separator();

					ImGui::Text("Total error: %f", p.total_error_ / p.nb_spheres_);
					ImGui::Text("Min error: %f", p.min_error_);
					ImGui::Text("Max error: %f", p.max_error_);

					ImGui::Separator();

					ImGui::Text("Pick the sphere under the mouse with I, split it with S, delete it with D");
					uint32 picked_index = INVALID_INDEX;
					if (picked_sphere_.is_valid() && !is_live_sphere_vertex(p, picked_sphere_, picked_index))
						picked_sphere_ = PVertex();
					if (picked_sphere_.is_valid())
					{
						ImGui::Text("Picked sphere:");
						const Vec3& sp = (*p.spheres_position_)[picked_index];
						ImGui::Text("Index: %u", picked_index);
						ImGui::Text("Center: (%f, %f, %f)", sp[0], sp[1], sp[2]);
						ImGui::Text("Radius: %f", (*p.spheres_radius_)[picked_index]);
						if (p.spheres_error_)
							ImGui::Text("Error: %f", (*p.spheres_error_)[picked_index]);
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
	SurfaceRender<NONMANIFOLD>* skeleton_render_ = nullptr;

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
