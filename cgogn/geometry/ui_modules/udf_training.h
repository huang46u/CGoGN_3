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
#include <cgogn/geometry/algos/udf/sample_processing.h>
#include <cgogn/geometry/algos/udf/spatial_query.h>
#include <cgogn/geometry/algos/udf/spheres_optimizer.h>
#include <cgogn/geometry/algos/udf/skeleton_topology.h>
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

#include <libacc/bvh_tree.h>
#include <libacc/kd_tree.h>

#include <GLFW/glfw3.h>
#include <algorithm>
#include <chrono>
#include <cmath>
#include <filesystem>
#include <limits>
#include <mutex>
#include <numeric>
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
	using SpheresOptimizerType = geometry::SpheresOptimizer<POINTS>;
	using SpheresOptimizerData = typename SpheresOptimizerType::Data;
	using SpheresOptimizerMetrics = typename SpheresOptimizerType::Metrics;
	using SpheresOptimizerStatus = typename SpheresOptimizerType::Status;
	using SkeletonTopologyType = geometry::SkeletonTopology<POINTS, NONMANIFOLD>;
	using SkeletonTopologyData = typename SkeletonTopologyType::Data;
	using SkeletonTopologyMetrics = typename SkeletonTopologyType::Metrics;
	using SkeletonTopologyStatus = typename SkeletonTopologyType::Status;

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
		std::shared_ptr<PAttribute<Scalar>> spheres_error_ = nullptr;
		std::shared_ptr<PAttribute<Scalar>> spheres_error_not_normalized_ = nullptr;
		std::shared_ptr<PAttribute<NMVertex>> spheres_skeleton_vertex_ = nullptr;
		std::unique_ptr<SpheresOptimizerType> spheres_optimizer_;
		std::unique_ptr<SkeletonTopologyType> skeleton_topology_;

		// Skeleton
		NONMANIFOLD* skeleton_ = nullptr;
		bool skeleton_invalidated_ = false;
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
		if (p.samples_mesh_ && p.samples_kdtree_)
			geometry::compute_samples_area(
				*p.samples_mesh_, *p.samples_position_, *p.samples_normal_, *p.samples_area_, *p.samples_knn_,
				*p.samples_kdtree_, p.samples_kdtree_vertices_, p.knn_k_, Scalar(p.alpha_));
		if (p.samples_mesh_)
			geometry::compute_quadrics(*p.samples_mesh_, *p.samples_position_, *p.samples_normal_, *p.samples_area_,
										*p.samples_knn_, *p.samples_quadric_, *p.samples_line_quadric_, p.knn_k_);
	}

	void headless_compute_fitting_primitives_prepared(POINTS& points)
	{
		headless_compute_fitting_primitives_prepared(points_parameters_[&points]);
	}

	void headless_init_spheres_prepared(POINTS& points)
	{
		auto& p = points_parameters_[&points];
		if (optimizer_for(p).initialize_from_samples(Scalar(p.init_dilation_constant_)))
		{
			p.skeleton_invalidated_ = true;
			if (p.skeleton_)
			{
				clear(*p.skeleton_);
				p.skeleton_topology_.reset();
				if (p.spheres_skeleton_vertex_)
					p.spheres_skeleton_vertex_->fill(NMVertex());
			}
		}
		sync_optimizer_metrics(p);
	}

	HeadlessOptimizationStats headless_optimize_spheres_prepared(PointsParameters& p, bool verbose = false)
	{
		const OutputVerbosity output_verbosity = verbose ? OUTPUT_VERBOSE : p.output_verbosity_;
		HeadlessOptimizationStats stats;
		p.running_ = true;
		p.stopping_ = false;
		auto& optimizer = optimizer_for(p);
		optimizer.reset_optimization();
		sync_optimizer_metrics(p);
		p.pending_full_refresh_after_stop_ = false;
		auto optimization_start = std::chrono::high_resolution_clock::now();
		SpheresOptimizerStatus status = SpheresOptimizerStatus::running;
		while (status == SpheresOptimizerStatus::running)
		{
			status = optimizer.update_once(Scalar(p.sqem_update_lambda_line_plane_));
			sync_optimizer_metrics(p);
			if (optimizer.metrics().sphere_topology_changed)
				p.skeleton_invalidated_ = true;
			log_basic(output_verbosity, "[HeadlessOptimize] iteration=", p.iteration_count_, " spheres=",
					  p.nb_spheres_, " error=", p.total_error_, " diff=", p.total_error_diff_, '\n');
		}

		auto optimization_end = std::chrono::high_resolution_clock::now();
		stats.optimization_total_ms_ =
			std::chrono::duration<float64, std::milli>(optimization_end - optimization_start).count();
		stats.cluster_total_ms_ = optimizer.metrics().cluster_total_ms;
		stats.sphere_update_total_ms_ = optimizer.metrics().sphere_update_total_ms;
		stats.error_total_ms_ = optimizer.metrics().error_total_ms;
		stats.optimization_iterations_ = optimizer.metrics().iteration;
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
		const SkeletonTopologyStatus status = skeleton_topology_for(p).build_from_spheres();
		if (status == SkeletonTopologyStatus::success)
		{
			p.skeleton_invalidated_ = false;
			log_basic(p, "[SkeletonBuild] tets=", skeleton_topology_for(p).metrics().tet_count, '\n');
			refresh_skeleton_topology_colors(p);
		}
	}

	void headless_build_skeleton_prepared(POINTS& points)
	{
		headless_build_skeleton_prepared(points_parameters_[&points]);
	}

	void headless_run_topology_fix_prepared(PointsParameters& p)
	{
		if (p.skeleton_invalidated_)
		{
			log_error(p, "Topology fix requires a skeleton built from the current sphere set.", '\n');
			return;
		}
		auto evaluator = [this, &p](const std::vector<Vec3>& points, std::vector<Scalar>& values) {
			return eval_topology_score_values(p, points, values);
		};
		auto& topology = skeleton_topology_for(p);
		const SkeletonTopologyStatus status = topology.fix_topology(evaluator);
		if (status == SkeletonTopologyStatus::success)
		{
			log_skeleton_topology_summary(p, "[TopologyFull]");
			apply_skeleton_topology_results(p);
		}
		else
			log_error(p, "[TopologyFull] score evaluation failed or topology data is invalid.", '\n');
	}

	void headless_run_topology_fix_prepared(POINTS& points)
	{
		headless_run_topology_fix_prepared(points_parameters_[&points]);
	}
	void headless_run_completion_residual_prune_prepared(PointsParameters& p)
	{
		if (skeleton_topology_for(p).prune_residual_sheets() == SkeletonTopologyStatus::success)
		{
			log_skeleton_topology_summary(p, "[ResidualSheetPrune]");
			apply_skeleton_topology_results(p);
			if (colorize_skeleton_face_components(p, "[FaceComponentsColor-Residual-Final]") &&
				non_manifold_provider_ && p.skeleton_ && p.skeleton_face_component_color_)
				non_manifold_provider_->emit_attribute_changed(*p.skeleton_, p.skeleton_face_component_color_.get());
		}
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
		counts.skeleton_vertices_ = p.skeleton_ && !p.skeleton_invalidated_ ? nb_cells<NMVertex>(*p.skeleton_) : 0;
		counts.skeleton_edges_ = p.skeleton_ && !p.skeleton_invalidated_ ? nb_cells<NMEdge>(*p.skeleton_) : 0;
		counts.skeleton_faces_ = p.skeleton_ && !p.skeleton_invalidated_ ? nb_cells<NMFace>(*p.skeleton_) : 0;
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
		if (p.skeleton_invalidated_ || !p.skeleton_ || !p.skeleton_position_ || !p.skeleton_radius_ || !p.spheres_ || !p.spheres_radius_)
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
			skeleton_topology_for(p).prune_fully_non_manifold_triangles();
			if (skeleton_topology_for(p).compute_skeleton_face_components_union_find() == SkeletonTopologyStatus::success)
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
			p.input_mode_ = INPUT_NEURAL_UDF;
			log_basic(p, "Loaded neural UDF model from: ", model_path, '\n');
		}
		catch (const c10::Error& e)
		{
			log_error(p, "Error loading Neural UDF model: ", e.what(), '\n');
			p.neural_udf_loaded_ = false;
		}
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
					if (p.skeleton_invalidated_ && p.skeleton_)
					{
						clear(*p.skeleton_);
						p.skeleton_topology_.reset();
						if (p.spheres_skeleton_vertex_)
							p.spheres_skeleton_vertex_->fill(NMVertex());
					}
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
			geometry::compute_input_normals(*p.points_, *p.position_, *p.normal_, *p.knn_, *p.input_kdtree_,
										p.input_kdtree_vertices_, p.knn_k_);
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
		if (p.samples_mesh_ && p.samples_kdtree_)
			geometry::compute_samples_area(
				*p.samples_mesh_, *p.samples_position_, *p.samples_normal_, *p.samples_area_, *p.samples_knn_,
				*p.samples_kdtree_, p.samples_kdtree_vertices_, p.knn_k_, Scalar(p.alpha_));
		log_basic(p, "Computing Quadrics...", '\n');
		if (p.samples_mesh_)
			geometry::compute_quadrics(*p.samples_mesh_, *p.samples_position_, *p.samples_normal_, *p.samples_area_,
										*p.samples_knn_, *p.samples_quadric_, *p.samples_line_quadric_, p.knn_k_);

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
		SpheresOptimizerType& optimizer = optimizer_for(p);
		if (optimizer.initialize_from_samples(Scalar(p.init_dilation_constant_)))
		{
			p.skeleton_invalidated_ = true;
			if (p.skeleton_)
			{
				clear(*p.skeleton_);
				p.skeleton_topology_.reset();
				if (p.spheres_skeleton_vertex_)
					p.spheres_skeleton_vertex_->fill(NMVertex());
			}
		}
		sync_optimizer_metrics(p);

		if (!p.running_)
		{
			set_post_init_sphere_render_state(p);
		}
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

	SpheresOptimizerData make_spheres_optimizer_data(PointsParameters& p)
	{
		SpheresOptimizerData data;
		data.samples_mesh = p.samples_mesh_;
		data.spheres = p.spheres_;
		data.sample_position = p.samples_position_;
		data.sample_area = p.samples_area_;
		data.sample_knn = p.samples_knn_;
		data.sample_quadric = p.samples_quadric_;
		data.sample_line_quadric = p.samples_line_quadric_;
		data.sample_ma_position = p.samples_ma_position_;
		data.sample_ma_radius = p.samples_ma_radius_;
		data.sample_ma_secondary_vertex = p.samples_ma_secondary_vertex_;
		data.sample_sphere = p.samples_sphere_;
		data.sample_error = p.samples_error_;
		data.sphere_position = p.spheres_position_;
		data.sphere_radius = p.spheres_radius_;
		data.sphere_cluster = p.spheres_cluster_;
		data.sphere_cluster_area = p.spheres_cluster_area_;
		data.sphere_cluster_color = p.spheres_cluster_color_;
		data.sphere_neighbors = p.spheres_neighbor_clusters_;
		data.sphere_error = p.spheres_error_;
		data.sphere_error_not_normalized = p.spheres_error_not_normalized_;
		data.sample_kdtree = p.samples_kdtree_;
		data.sample_kdtree_vertices = &p.samples_kdtree_vertices_;
		data.sqem_update_lambda_line_plane = Scalar(p.sqem_update_lambda_line_plane_);
		data.sphere_count = &p.nb_spheres_;
		return data;
	}

	SpheresOptimizerType& optimizer_for(PointsParameters& p)
	{
		SpheresOptimizerData data = make_spheres_optimizer_data(p);
		if (!p.spheres_optimizer_)
			p.spheres_optimizer_ = std::make_unique<SpheresOptimizerType>(std::move(data));
		else
			p.spheres_optimizer_->data() = std::move(data);
		return *p.spheres_optimizer_;
	}

	SkeletonTopologyType& skeleton_topology_for(PointsParameters& p)
	{
		SkeletonTopologyData data;
		data.skeleton = p.skeleton_;
		data.spheres_optimizer = &optimizer_for(p);
		data.sphere_skeleton_vertex = p.spheres_skeleton_vertex_;
		data.skeleton_position = p.skeleton_position_;
		data.skeleton_radius = p.skeleton_radius_;
		data.skeleton_source_sphere = p.skeleton_source_sphere_;
		data.face_incident_tets = p.incident_tets_;
		data.face_component_id = p.skeleton_face_component_id_;
		data.edge_degree = p.edge_degree_;
		if (!p.skeleton_topology_)
			p.skeleton_topology_ = std::make_unique<SkeletonTopologyType>(std::move(data));
		else
			p.skeleton_topology_->data() = std::move(data);
		return *p.skeleton_topology_;
	}

	void log_skeleton_topology_summary(const PointsParameters& p, const char* prefix) const
	{
		if (!p.skeleton_topology_)
			return;
		const SkeletonTopologyMetrics& metrics = p.skeleton_topology_->metrics();
		log_basic(p, prefix, " rounds=", metrics.topology_rounds,
				  " boundary_removed_tets=", metrics.boundary_removed_tets,
				  " simple_removed_tets=", metrics.simple_removed_tets,
				  " nonsimple_removed_tets=", metrics.nonsimple_removed_tets,
				  " removed_faces=", metrics.removed_faces,
				  " removed_edges=", metrics.removed_edges,
				  " removed_tets=", metrics.removed_tets,
				  " remaining_tets=", metrics.tet_count,
				  " repaired_regions=", metrics.repaired_dense_regions,
				  " added_spheres=", metrics.added_spheres.size(),
				  " added_faces=", metrics.added_faces,
				  " residual_removed_sheets=", metrics.residual_removed_sheets,
				  " residual_removed_faces=", metrics.residual_removed_faces, '\n');
	}

	void apply_skeleton_topology_results(PointsParameters& p)
	{
		if (!p.skeleton_topology_)
			return;
		sync_optimizer_metrics(p);
		for (PVertex sphere : p.skeleton_topology_->metrics().added_spheres)
		{
			const uint32 id = p.spheres_ ? index_of(*p.spheres_, sphere) : INVALID_INDEX;
			if (id == INVALID_INDEX)
				continue;
			if (p.spheres_color_)
				(*p.spheres_color_)[id] = Vec4(0.95, 0.25, 0.15, 1.0);
			if (p.spheres_cluster_color_)
				(*p.spheres_cluster_color_)[id] = Vec4(0.95, 0.25, 0.15, 1.0);
		}
		refresh_skeleton_topology_colors(p);
		if (non_manifold_provider_ && p.skeleton_)
		{
			non_manifold_provider_->emit_connectivity_changed(*p.skeleton_);
			if (p.skeleton_position_) non_manifold_provider_->emit_attribute_changed(*p.skeleton_, p.skeleton_position_.get());
			if (p.skeleton_radius_) non_manifold_provider_->emit_attribute_changed(*p.skeleton_, p.skeleton_radius_.get());
		}
		if (points_provider_ && p.spheres_)
		{
			points_provider_->emit_connectivity_changed(*p.spheres_);
			if (p.spheres_position_) points_provider_->emit_attribute_changed(*p.spheres_, p.spheres_position_.get());
			if (p.spheres_radius_) points_provider_->emit_attribute_changed(*p.spheres_, p.spheres_radius_.get());
			if (p.spheres_color_) points_provider_->emit_attribute_changed(*p.spheres_, p.spheres_color_.get());
			if (p.spheres_cluster_color_)
				points_provider_->emit_attribute_changed(*p.spheres_, p.spheres_cluster_color_.get());
		}
	}

	void sync_optimizer_metrics(PointsParameters& p)
	{
		if (!p.spheres_optimizer_)
			return;
		const SpheresOptimizerMetrics& metrics = p.spheres_optimizer_->metrics();
		p.nb_spheres_ = metrics.sphere_count;
		p.iteration_count_ = metrics.iteration;
		p.total_error_ = metrics.total_error;
		p.total_error_not_normalized_ = metrics.total_error_not_normalized;
		p.last_total_error_ = metrics.last_total_error;
		p.total_error_diff_ = metrics.error_difference;
		p.min_error_ = metrics.minimum_error;
		p.max_error_ = metrics.maximum_error;
	}
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
	void refresh_skeleton_topology_colors(PointsParameters& p)
	{
		skeleton_topology_for(p).compute_edge_degree();
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
		{
			const SkeletonTopologyStatus status = skeleton_topology_for(p).build_from_spheres();
			p.skeleton_invalidated_ = status != SkeletonTopologyStatus::success;
			if (status == SkeletonTopologyStatus::success)
				refresh_skeleton_topology_colors(p);
		}
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
		SpheresOptimizerType& optimizer = optimizer_for(p);
		optimizer.reset_optimization();
		sync_optimizer_metrics(p);
		p.pending_full_refresh_after_stop_ = false;

		launch_thread([this, &p, &optimizer]() {
			auto start = std::chrono::high_resolution_clock::now();
			SpheresOptimizerStatus status = SpheresOptimizerStatus::running;
			while (status == SpheresOptimizerStatus::running)
			{
				{
					std::lock_guard<std::mutex> lock(p.mutex_);
					if (!p.stopping_)
					{
						log_basic(p, "Start Sphere update", '\n');
						status = optimizer.update_once(Scalar(p.sqem_update_lambda_line_plane_));
						sync_optimizer_metrics(p);
						if (optimizer.metrics().sphere_topology_changed)
							p.skeleton_invalidated_ = true;
					}
				}
				if (p.stopping_)
					break;
				if (p.slow_down_)
					std::this_thread::sleep_for(std::chrono::microseconds(1000000 / p.update_rate_));
				else
					std::this_thread::yield();

				log_basic(p, "Iteration: ", p.iteration_count_, " | Spheres: ", p.nb_spheres_, " | Error: ",
						  p.total_error_, " | Diff: ", p.total_error_diff_, '\n');
				if (status == SpheresOptimizerStatus::converged)
				{
					log_basic(p, "Auto stop: error converged after post-convergence iterations.", '\n');
					p.stopping_ = true;
				}
				else if (status == SpheresOptimizerStatus::max_iterations)
				{
					log_basic(p, "Stop: reached max iterations (150).", '\n');
					p.stopping_ = true;
				}
				else if (status == SpheresOptimizerStatus::failed)
				{
					log_error(p, "Sphere optimizer is missing required data.", '\n');
					p.stopping_ = true;
				}
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
							SpheresOptimizerType& optimizer = optimizer_for(p);
							optimizer.update_once(Scalar(p.sqem_update_lambda_line_plane_));
							if (optimizer.metrics().sphere_topology_changed)
							{
								p.skeleton_invalidated_ = true;
								if (p.skeleton_)
								{
									clear(*p.skeleton_);
									p.skeleton_topology_.reset();
									if (p.spheres_skeleton_vertex_)
										p.spheres_skeleton_vertex_->fill(NMVertex());
								}
							}
							sync_optimizer_metrics(p);
							update_render_data(p, false, true, true);
						}
					}


					if (ImGui::Button("Build Skeleton"))
					{
						if (!p.running_)
						{
							std::lock_guard<std::mutex> lock(p.mutex_);
							const SkeletonTopologyStatus status = skeleton_topology_for(p).build_from_spheres();
							p.skeleton_invalidated_ = status != SkeletonTopologyStatus::success;
							if (status == SkeletonTopologyStatus::success)
							{
								log_basic(p, "[SkeletonBuild] tets=", skeleton_topology_for(p).metrics().tet_count, '\n');
								refresh_skeleton_topology_colors(p);
							}
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
							skeleton_topology_for(p).prune_fully_non_manifold_triangles();
							if (skeleton_topology_for(p).compute_skeleton_face_components_union_find() == SkeletonTopologyStatus::success &&
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
							if (p.skeleton_invalidated_)
							{
								log_error(p, "Topology fix requires a skeleton built from the current sphere set.", '\n');
								return;
							}
							auto evaluator = [this, &p](const std::vector<Vec3>& points, std::vector<Scalar>& values) {
								return eval_topology_score_values(p, points, values);
							};
							auto& topology = skeleton_topology_for(p);
							const SkeletonTopologyStatus status = topology.fix_topology(evaluator);
							if (status == SkeletonTopologyStatus::success)
							{
								log_skeleton_topology_summary(p, "[TopologyFull]");
								apply_skeleton_topology_results(p);
							}
							else
								log_error(p, "[TopologyFull] score evaluation failed or topology data is invalid.", '\n');
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

					ImGui::Separator();

					ImGui::Text("Total error: %f", p.total_error_ / p.nb_spheres_);
					ImGui::Text("Min error: %f", p.min_error_);
					ImGui::Text("Max error: %f", p.max_error_);

					ImGui::Separator();

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
	std::shared_ptr<boost::synapse::connection> timer_connection_;

	torch::Device device_ = torch::kCPU;
};

} // namespace ui

} // namespace cgogn

#endif // CGOGN_MODULE_UDF_TRAINING_H_
