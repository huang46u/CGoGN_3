#ifndef CGOGN_MODULE_UDF_TRAINING_H_
#define CGOGN_MODULE_UDF_TRAINING_H_

#include <cgogn/ui/app.h>
#include <cgogn/ui/module.h>

#include <cgogn/core/ui_modules/mesh_provider.h>
#include <cgogn/core/utils/thread.h>
#if !defined(CGOGN_UDF_WITH_TORCH) || CGOGN_UDF_WITH_TORCH
#include <cgogn/geometry/algos/udf/neural/neural_distance_field.h>
#endif
#include <cgogn/geometry/algos/udf/neural_model_type.h>
#include <cgogn/geometry/algos/udf/reconstruction.h>
#include <cgogn/geometry/algos/udf/reconstruction_options.h>
#include <cgogn/geometry/algos/udf/reconstruction_profiles.h>
#include <cgogn/geometry/algos/udf/skeleton_face_components.h>
#include <cgogn/geometry/algos/udf/skeleton_ply_export.h>
#include <cgogn/geometry/types/vector_traits.h>
#include <cgogn/ui/portable-file-dialogs.h>

#include <cgogn/rendering/ui_modules/point_cloud_render.h>
#include <cgogn/rendering/ui_modules/surface_render.h>

#include <GLFW/glfw3.h>
#include <algorithm>
#include <atomic>
#include <chrono>
#include <cmath>
#include <filesystem>
#include <iostream>
#include <limits>
#include <map>
#include <memory>
#include <mutex>
#include <set>
#include <stdexcept>
#include <string>
#include <thread>
#include <unordered_set>
#if !defined(CGOGN_UDF_WITH_TORCH) || CGOGN_UDF_WITH_TORCH
#include <torch/torch.h>
#endif
#include <utility>

namespace cgogn
{
namespace ui
{
using geometry::Scalar;
using geometry::Vec3;
using geometry::Vec4;

template <typename SURFACE, typename POINTS, typename NONMANIFOLD, typename RaySamplerTag>
class UDFTraining : public ViewModule
{
public:
	using PVertex = typename mesh_traits<POINTS>::Vertex;
	using NMVertex = typename mesh_traits<NONMANIFOLD>::Vertex;
	using NMFace = typename mesh_traits<NONMANIFOLD>::Face;
	using NMEdge = typename mesh_traits<NONMANIFOLD>::Edge;

	template <typename T>
	using PAttribute = typename mesh_traits<POINTS>::template Attribute<T>;
	template <typename T>
	using NMAttribute = typename mesh_traits<NONMANIFOLD>::template Attribute<T>;
	using SpheresOptimizerType = geometry::SpheresOptimizer<POINTS>;
	using SpheresOptimizerMetrics = typename SpheresOptimizerType::Metrics;
	using SpheresOptimizerStatus = typename SpheresOptimizerType::Status;
	using ReconstructionType = geometry::UDFReconstruction<SURFACE, POINTS, NONMANIFOLD, RaySamplerTag>;
	using ReconstructionStatus = typename ReconstructionType::Status;
	using ReconstructionStatusInfo = typename ReconstructionType::StatusInfo;

	enum OutputVerbosity : uint32
	{
		OUTPUT_MUTE,
		OUTPUT_NORMAL
	};

private:
	// "<action>: <message>" for the log and the panel; the colon is left out when the message is empty
	static std::string with_message(const std::string& action, const std::string& message)
	{
		return message.empty() ? action : action + ": " + message;
	}

	struct PointsParameters
	{
		bool initialized_ = false;
		bool fitting_data_computed_ = false; // UI workflow gate for the sphere-fitting controls.
		bool spheres_initialized_ = false; // UI workflow gate for the sphere-optimization controls; reset with the fitting data.
		std::unique_ptr<ReconstructionType> reconstruction_;
		// Input Points
		POINTS* points_ = nullptr;
		std::shared_ptr<PAttribute<Vec3>> position_ = nullptr;

		// Neural model information shown in the UI; the model itself is owned by ReconstructionType.
		bool neural_udf_loaded_ = false;
		std::string neural_udf_model_path_ = "";

		// Reason of the last failed stage, shown in red in the Neural UDF panel; cleared when the stage succeeds.
		std::string last_error_;

		// Dataset profile currently applied to the options (see reconstruction_profile_names()).
		std::string profile_ = "default";

		// Sampling & Fitting Data
		POINTS* samples_mesh_ = nullptr;
		std::shared_ptr<PAttribute<Vec3>> samples_position_ = nullptr;
		std::shared_ptr<PAttribute<Vec3>> samples_normal_ = nullptr;

		// Medial Axis (on samples)

		// Clustering (on samples)
		std::shared_ptr<PAttribute<PVertex>> samples_sphere_ = nullptr; // Cluster ID for each sample
		std::shared_ptr<PAttribute<Vec4>> samples_color_ = nullptr;
		std::shared_ptr<PAttribute<Vec4>> samples_normal_color_ = nullptr;

		// Spheres
		POINTS* spheres_ = nullptr;
		std::shared_ptr<PAttribute<Vec3>> spheres_position_ = nullptr;
		std::shared_ptr<PAttribute<Scalar>> spheres_radius_ = nullptr;
		std::shared_ptr<PAttribute<Vec4>> spheres_color_ = nullptr;
		std::shared_ptr<PAttribute<Vec4>> spheres_cluster_color_ = nullptr;
		std::shared_ptr<PAttribute<Scalar>> spheres_error_ = nullptr;

		// Skeleton
		NONMANIFOLD* skeleton_ = nullptr;
		std::atomic<bool> skeleton_invalidated_{false};
		std::shared_ptr<NMAttribute<Vec3>> skeleton_position_ = nullptr;
		std::shared_ptr<NMAttribute<Scalar>> skeleton_radius_ = nullptr;
		std::shared_ptr<NMAttribute<std::set<std::size_t>>> incident_tets_ = nullptr;

		std::shared_ptr<NMAttribute<Vec3>> skeleton_face_color_ = nullptr;
		std::shared_ptr<NMAttribute<Vec3>> skeleton_edge_color_ = nullptr;

		bool error_as_spheres_color_ = false;
		float32 spheres_transparency_ = 0.5f;

		// UI state
		std::atomic<OutputVerbosity> output_verbosity_{OUTPUT_NORMAL};

		// Threading
		std::mutex mutex_;
		std::atomic<bool> running_{false};
		std::atomic<bool> stop_requested_{false};
		std::atomic<bool> render_refresh_requested_{false};
		std::atomic<bool> preview_render_during_sphere_update_{true};
		std::atomic<bool> slow_down_{true};
		std::atomic<uint32> update_rate_{20};
		std::atomic<float32> sphere_update_lambda_{0.2f};
		SpheresOptimizerMetrics metrics_snapshot_{};
	};

public:
	UDFTraining(const App& app) : ViewModule(app, "UDFTraining")
	{
	}

	~UDFTraining() override
	{
		stop_and_join_spheres_update();
	}

	static int output_verbosity_index(OutputVerbosity value)
	{
		switch (value)
		{
		case OUTPUT_MUTE:
			return 0;
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
		case 1:
		default:
			return OUTPUT_NORMAL;
		}
	}

	bool is_basic_logging_enabled(const PointsParameters& p) const
	{
		return p.output_verbosity_.load() != OUTPUT_MUTE;
	}

	bool is_basic_logging_enabled(OutputVerbosity output_verbosity) const
	{
		return output_verbosity != OUTPUT_MUTE;
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
	void log_error(Args&&... args) const
	{
		(std::cerr << ... << std::forward<Args>(args));
	}

	void set_selected_surface(SURFACE& s)
	{
		if (!surface_provider_)
			throw std::runtime_error("UDFTraining surface provider is not initialized before set_selected_surface().");
		selected_surface_ = &s;
		if (selected_points_ && points_parameters_[selected_points_].reconstruction_)
			points_parameters_[selected_points_].reconstruction_->data().surface = &s;
	}

	void set_selected_points(POINTS& p)
	{
		if (!points_provider_)
			throw std::runtime_error("UDFTraining points provider is not initialized before set_selected_points().");
		if (!non_manifold_provider_)
			throw std::runtime_error(
				"UDFTraining non-manifold provider is not initialized before set_selected_points().");
		selected_points_ = &p;
		init_points_data(p);
		PointsParameters& params = points_parameters_[selected_points_];
		if (params.reconstruction_)
			params.reconstruction_->data().surface = selected_surface_;
	}

	// Reset every option to the table defaults, then the model-type defaults, then the dataset profile.
	// The neural model type of the current reconstruction is kept.
	void reset_options_to_profile(PointsParameters& p, const std::string& profile)
	{
		if (!p.reconstruction_)
			return;
		auto& options = p.reconstruction_->options();
		const cgogn::geometry::NeuralModelType model_type = options.neural_model_type;
		options = cgogn::geometry::ReconstructionOptions{};
		options.neural_model_type = model_type;
		options.apply_filtering = cgogn::geometry::default_apply_filtering(model_type);
		options.recompute_normals_after_sampling =
			cgogn::geometry::default_recompute_normals_after_sampling(model_type);
		if (!cgogn::geometry::apply_reconstruction_profile(profile, options))
			log_error("Unknown reconstruction profile: ", profile, '\n');
		p.profile_ = profile;
		p.sphere_update_lambda_.store(options.sqem_update_lambda_line_plane);
	}

	ReconstructionStatusInfo load_neural_udf_model(POINTS& points, const std::string& model_path,
											   cgogn::geometry::NeuralModelType model_type)
	{
		PointsParameters& p = points_parameters_[&points];
		init_points_data(points);
		p.neural_udf_loaded_ = false;
		p.neural_udf_model_path_.clear();
		if (!p.initialized_ || !p.reconstruction_)
		{
			p.last_error_ = "Failed to initialize UDF reconstruction data.";
			log_error(p.last_error_, '\n');
			return ReconstructionStatus::invalid_data;
		}
		try
		{
			if (!std::filesystem::exists(model_path))
			{
				p.last_error_ = "Neural UDF model file does not exist: " + model_path;
				log_error(p.last_error_, '\n');
				return ReconstructionStatus::model_load_failed;
			}
			log_basic(p, "Loading Neural UDF model from: ", model_path, '\n');
			const auto status = p.reconstruction_->load_neural_model(model_path, model_type);
			p.neural_udf_loaded_ = status.code == ReconstructionStatus::success;
			if (!p.neural_udf_loaded_)
			{
				p.last_error_ = with_message("Error loading Neural UDF model: loader rejected the model", status.message);
				log_error(p.last_error_, '\n');
				return status;
			}
			p.neural_udf_model_path_ = model_path;
			p.last_error_.clear();
			// dataset profile guessed from the model type and path; the Sampling panel can change it
			reset_options_to_profile(p, cgogn::geometry::guess_reconstruction_profile(true, model_type, model_path));
			log_basic(p, "Reconstruction profile: ", p.profile_, '\n');
			log_basic(p, "Loaded neural UDF model from: ", model_path, '\n');
			return ReconstructionStatus::success;
		}
		catch (const std::exception& e)
		{
			p.last_error_ = std::string("Error loading Neural UDF model: ") + e.what();
			log_error(p.last_error_, '\n');
			p.neural_udf_loaded_ = false;
			return ReconstructionStatus::model_load_failed;
		}
		catch (...)
		{
			p.last_error_ = "Error loading Neural UDF model: unknown exception.";
			log_error(p.last_error_, '\n');
			p.neural_udf_loaded_ = false;
			return ReconstructionStatus::model_load_failed;
		}
	}

	void load_alpha_samples_to_mesh(PointsParameters& p)
	{
		if (!p.reconstruction_)
			return;
		const auto status = p.reconstruction_->sample();
		if (status.code != ReconstructionStatus::success)
		{
			p.last_error_ = with_message("Alpha level set sampling failed", status.message);
			log_error(p.last_error_, '\n');
			return;
		}
		p.last_error_.clear();
		if (p.samples_color_)
			p.samples_color_->fill(Vec4(0.0, 0.0, 0.0, 1.0));
		refresh_sample_normals_color(p);
		points_provider_->emit_connectivity_changed(*p.samples_mesh_);
		points_provider_->emit_attribute_changed(*p.samples_mesh_, p.samples_position_.get());
		points_provider_->emit_attribute_changed(*p.samples_mesh_, p.samples_normal_.get());
		if (p.samples_color_)
			points_provider_->emit_attribute_changed(*p.samples_mesh_, p.samples_color_.get());
		p.fitting_data_computed_ = false;
		log_basic(p, "Alpha level set sampling complete. Ready for fitting.", '\n');
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

		pcr_ = static_cast<PointCloudRender<POINTS>*>(
			app_.module("PointCloudRender (" + std::string{mesh_traits<POINTS>::name} + ")"));
		skeleton_render_ = static_cast<SurfaceRender<NONMANIFOLD>*>(
			app_.module("SurfaceRender (" + std::string{mesh_traits<NONMANIFOLD>::name} + ")"));

		timer_connection_ = boost::synapse::connect<App::timer_tick>(&app_, [this]() { handle_spheres_update_tick(); });
		const OutputVerbosity startup_output_verbosity =
			selected_points_ ? points_parameters_[selected_points_].output_verbosity_.load() : OUTPUT_NORMAL;
#if !defined(CGOGN_UDF_WITH_TORCH) || CGOGN_UDF_WITH_TORCH
		if (torch::cuda::is_available())
		{
			log_basic(startup_output_verbosity, "CUDA is available! Using GPU device 0.", '\n');
		}
		else
		{
			log_basic(startup_output_verbosity, "CUDA is not available! Using the CPU.", '\n');
		}
#else
		log_basic(startup_output_verbosity, "Built without LibTorch: neural UDF input is not available.", '\n');
#endif
	}

private:
	void init_points_data(POINTS& m)
	{
		PointsParameters& p = points_parameters_[&m];
		if (p.initialized_)
			return;

		p.points_ = &m;
		p.position_ = get_attribute<Vec3, PVertex>(m, "position");

		const std::string sample_name = points_provider_->mesh_name(m) + "_samples";
		p.samples_mesh_ = points_provider_->has_mesh(sample_name) ? points_provider_->mesh(sample_name)
																  : points_provider_->add_mesh(sample_name);
		points_provider_->clear_mesh(*p.samples_mesh_);
		p.samples_position_ = get_or_add_attribute<Vec3, PVertex>(*p.samples_mesh_, "position");
		p.samples_normal_ = get_or_add_attribute<Vec3, PVertex>(*p.samples_mesh_, "normal");
		p.samples_sphere_ = get_or_add_attribute<PVertex, PVertex>(*p.samples_mesh_, "sphere");
		p.samples_color_ = get_or_add_attribute<Vec4, PVertex>(*p.samples_mesh_, "color");
		p.samples_normal_color_ = get_or_add_attribute<Vec4, PVertex>(*p.samples_mesh_, "normal_color");

		const std::string sphere_name = points_provider_->mesh_name(m) + "_spheres";
		p.spheres_ = points_provider_->has_mesh(sphere_name) ? points_provider_->mesh(sphere_name)
															 : points_provider_->add_mesh(sphere_name);
		p.spheres_position_ = get_or_add_attribute<Vec3, PVertex>(*p.spheres_, "position");
		p.spheres_radius_ = get_or_add_attribute<Scalar, PVertex>(*p.spheres_, "radius");
		p.spheres_color_ = get_or_add_attribute<Vec4, PVertex>(*p.spheres_, "color");
		p.spheres_cluster_color_ = get_or_add_attribute<Vec4, PVertex>(*p.spheres_, "cluster_color");
		p.spheres_error_ = get_or_add_attribute<Scalar, PVertex>(*p.spheres_, "error");

		const std::string skeleton_name = points_provider_->mesh_name(m) + "_skeleton";
		p.skeleton_ = non_manifold_provider_->has_mesh(skeleton_name) ? non_manifold_provider_->mesh(skeleton_name)
																	  : non_manifold_provider_->add_mesh(skeleton_name);
		p.skeleton_position_ = get_or_add_attribute<Vec3, NMVertex>(*p.skeleton_, "position");
		p.skeleton_radius_ = get_or_add_attribute<Scalar, NMVertex>(*p.skeleton_, "radius");
		p.incident_tets_ = get_or_add_attribute<std::set<std::size_t>, NMFace>(*p.skeleton_, "incident_tets");
		p.skeleton_face_color_ = get_or_add_attribute<Vec3, NMFace>(*p.skeleton_, "color");
		p.skeleton_edge_color_ = get_or_add_attribute<Vec3, NMEdge>(*p.skeleton_, "color");

		typename ReconstructionType::Data data{selected_surface_, p.points_, p.samples_mesh_, p.spheres_, p.skeleton_};
		p.reconstruction_ = std::make_unique<ReconstructionType>(data);
		p.sphere_update_lambda_.store(p.reconstruction_->options().sqem_update_lambda_line_plane);
		const auto status = p.reconstruction_->initialize_data();
		p.initialized_ = status.code == ReconstructionStatus::success;
		if (p.initialized_)
			p.last_error_.clear();
		else
			p.last_error_ = with_message("Failed to initialize reconstruction data", status.message);
	}

	void compute_fitting_data(PointsParameters& p)
	{
		if (!p.reconstruction_ || p.fitting_data_computed_)
			return;
		const auto status = p.reconstruction_->compute_fitting_data();
		if (status.code != ReconstructionStatus::success)
		{
			p.last_error_ = with_message("Failed to compute fitting data", status.message);
			log_error(p.last_error_, '\n');
			return;
		}
		p.last_error_.clear();
		p.fitting_data_computed_ = true;
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

	void update_spheres_color(PointsParameters& p)
	{
		parallel_foreach_cell(*p.spheres_, [&](PVertex v) -> bool {
			uint32 v_index = index_of(*p.spheres_, v);
			if (p.error_as_spheres_color_)
				(*p.spheres_color_)[v_index] =
					color_map((*p.spheres_error_)[v_index], p.metrics_snapshot_.minimum_error,
							  p.metrics_snapshot_.maximum_error, p.spheres_transparency_);
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
		if (!(max > min))
			return Vec4(0.0, 0.0, 1.0, transparency);
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

	void apply_skeleton_topology_results(PointsParameters& p)
	{
		if (!p.reconstruction_)
			return;
		const auto& metrics = p.reconstruction_->skeleton_topology().metrics();
		for (PVertex sphere : metrics.added_spheres)
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
			if (p.skeleton_position_)
				non_manifold_provider_->emit_attribute_changed(*p.skeleton_, p.skeleton_position_.get());
			if (p.skeleton_radius_)
				non_manifold_provider_->emit_attribute_changed(*p.skeleton_, p.skeleton_radius_.get());
			if (p.skeleton_edge_color_)
				non_manifold_provider_->emit_attribute_changed(*p.skeleton_, p.skeleton_edge_color_.get());
			if (p.skeleton_face_color_)
				non_manifold_provider_->emit_attribute_changed(*p.skeleton_, p.skeleton_face_color_.get());
		}
		if (points_provider_ && p.spheres_)
		{
			points_provider_->emit_connectivity_changed(*p.spheres_);
			if (p.spheres_position_)
				points_provider_->emit_attribute_changed(*p.spheres_, p.spheres_position_.get());
			if (p.spheres_radius_)
				points_provider_->emit_attribute_changed(*p.spheres_, p.spheres_radius_.get());
			if (p.spheres_color_)
				points_provider_->emit_attribute_changed(*p.spheres_, p.spheres_color_.get());
			if (p.spheres_cluster_color_)
				points_provider_->emit_attribute_changed(*p.spheres_, p.spheres_cluster_color_.get());
		}
	}
	void refresh_skeleton_topology_colors(PointsParameters& p)
	{
		p.reconstruction_->skeleton_topology().compute_edge_degree();
		// Edge colors as in V0: degree 1 blue, degree 2 green, degree > 2 red, all others dim gray.
		foreach_cell(*p.skeleton_, [&](NMEdge e) {
			const std::size_t degree = incident_faces(*p.skeleton_, e).size();
			Vec3 color(0.15, 0.15, 0.15);
			if (degree == 1)
				color = Vec3(0.0, 0.0, 1.0);
			else if (degree == 2)
				color = Vec3(0.0, 1.0, 0.0);
			else if (degree > 2)
				color = Vec3(1.0, 0.2, 0.1);
			value<Vec3>(*p.skeleton_, p.skeleton_edge_color_, e) = color;
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
	}

	// The skeleton coloring helpers below only change colors; they never change the skeleton topology.
	bool skeleton_colors_ready(const PointsParameters& p) const
	{
		return p.skeleton_ && p.skeleton_edge_color_ && p.skeleton_face_color_;
	}

	void emit_skeleton_colors_changed(PointsParameters& p)
	{
		if (!non_manifold_provider_ || !p.skeleton_)
			return;
		non_manifold_provider_->emit_attribute_changed(*p.skeleton_, p.skeleton_edge_color_.get());
		non_manifold_provider_->emit_attribute_changed(*p.skeleton_, p.skeleton_face_color_.get());
	}

	// Shows the skeleton with the colors that a color button just wrote: per-face color vectors and per-edge colors.
	void show_skeleton_colors(PointsParameters& p)
	{
		using SkeletonRender = SurfaceRender<NONMANIFOLD>;
		if (!skeleton_render_ || !p.skeleton_ || !p.skeleton_face_color_ || !p.skeleton_edge_color_)
			return;
		for (View* v : linked_views_)
		{
			skeleton_render_->set_face_color(*v, *p.skeleton_, p.skeleton_face_color_);
			skeleton_render_->set_edge_color(*v, *p.skeleton_, p.skeleton_edge_color_);
			skeleton_render_->set_face_color_mode(*v, *p.skeleton_, SkeletonRender::PER_FACE, SkeletonRender::VECTOR);
			skeleton_render_->set_edge_color_mode(*v, *p.skeleton_, SkeletonRender::PER_EDGE);
		}
	}

	// V0 "Boundary tet mask": faces of boundary tets yellow, other faces dim gray. The longest degree-2 edge of
	// each boundary tet is orange, other edges dim gray. The V0 K5 mask is not included.
	void mark_boundary_tets_color(PointsParameters& p)
	{
		const auto mask = p.reconstruction_->skeleton_topology().boundary_tet_mask();
		const std::unordered_set<uint32> boundary_faces(mask.face_ids.begin(), mask.face_ids.end());
		const std::unordered_set<uint32> boundary_edges(mask.edge_ids.begin(), mask.edge_ids.end());
		foreach_cell(*p.skeleton_, [&](NMFace f) -> bool {
			const bool hit = boundary_faces.count(index_of(*p.skeleton_, f)) > 0;
			value<Vec3>(*p.skeleton_, p.skeleton_face_color_, f) = hit ? Vec3(1.0, 0.85, 0.1) : Vec3(0.12, 0.12, 0.12);
			return true;
		});
		foreach_cell(*p.skeleton_, [&](NMEdge e) -> bool {
			const bool hit = boundary_edges.count(index_of(*p.skeleton_, e)) > 0;
			value<Vec3>(*p.skeleton_, p.skeleton_edge_color_, e) = hit ? Vec3(1.0, 0.45, 0.05) : Vec3(0.1, 0.1, 0.1);
			return true;
		});
	}

	// V0 palette helpers, copied from V0 (golden-ratio hue spacing).
	static Vec3 hsv_to_rgb(Scalar h, Scalar s, Scalar v)
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

	static Vec3 component_palette_color(uint32 component_id)
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

	// V0 "Face components (UF)": color each skeleton face by its sheet. Faces are only recolored, not removed.
	void mark_skeleton_face_components_color(PointsParameters& p)
	{
		const auto sheet_of_face = geometry::compute_skeleton_face_components(*p.skeleton_);
		foreach_cell(*p.skeleton_, [&](NMFace f) -> bool {
			const auto it = sheet_of_face.find(index_of(*p.skeleton_, f));
			value<Vec3>(*p.skeleton_, p.skeleton_face_color_, f) =
				(it == sheet_of_face.end()) ? Vec3(0.0, 0.0, 0.0) : component_palette_color(it->second);
			return true;
		});
	}

	// V0 colormaps: cells with a score are colored by color_map over the range of those scores (0.5 if the range is
	// empty, as in V0); NaN (cell not adjacent to a tet) gets the dim fallback color.
	template <typename CELL>
	void color_cells_by_scores(PointsParameters& p, const std::vector<Scalar>& scores,
							   const std::shared_ptr<NMAttribute<Vec3>>& attribute, const Vec3& fallback)
	{
		Scalar score_min = std::numeric_limits<Scalar>::max();
		Scalar score_max = std::numeric_limits<Scalar>::lowest();
		for (const Scalar v : scores)
		{
			if (std::isnan(v))
				continue;
			score_min = std::min(score_min, v);
			score_max = std::max(score_max, v);
		}
		const bool has_range = score_max > score_min;
		foreach_cell(*p.skeleton_, [&](CELL c) -> bool {
			const uint32 id = index_of(*p.skeleton_, c);
			Vec3 color = fallback;
			if (id < scores.size() && !std::isnan(scores[id]))
			{
				const Vec4 rgba = has_range ? color_map(std::clamp(scores[id], score_min, score_max), score_min,
														score_max, 1.0f)
											: color_map(Scalar(0.5), Scalar(0), Scalar(1), 1.0f);
				color = Vec3(rgba[0], rgba[1], rgba[2]);
			}
			value<Vec3>(*p.skeleton_, attribute, c) = color;
			return true;
		});
	}

	// V0 "Tet face UDF colormap (Integral / Normalized)": read-only, only the face colors change.
	bool mark_tet_face_udf_color(PointsParameters& p, bool normalize_by_area)
	{
		if (!skeleton_colors_ready(p))
		{
			p.last_error_ = "No skeleton: build the skeleton first.";
			log_error(p.last_error_, '\n');
			return false;
		}
		const auto scores = p.reconstruction_->skeleton_face_scores(normalize_by_area);
		if (!scores)
		{
			p.last_error_ = "Failed to evaluate the face UDF scores (is the distance field loaded?).";
			log_error(p.last_error_, '\n');
			return false;
		}
		color_cells_by_scores<NMFace>(p, *scores, p.skeleton_face_color_, Vec3(0.2, 0.2, 0.2));
		p.last_error_.clear();
		return true;
	}

	// V0 "Edge UDF score colormap (Gauss3)": read-only, only the edge colors change.
	bool mark_tet_edge_gauss3_color(PointsParameters& p)
	{
		if (!skeleton_colors_ready(p))
		{
			p.last_error_ = "No skeleton: build the skeleton first.";
			log_error(p.last_error_, '\n');
			return false;
		}
		const auto scores = p.reconstruction_->skeleton_edge_scores_gauss3();
		if (!scores)
		{
			p.last_error_ = "Failed to evaluate the edge UDF scores (is the distance field loaded?).";
			log_error(p.last_error_, '\n');
			return false;
		}
		color_cells_by_scores<NMEdge>(p, *scores, p.skeleton_edge_color_, Vec3(0.08, 0.08, 0.08));
		p.last_error_.clear();
		return true;
	}

protected:
	enum class RenderUpdateMode : uint8
	{
		spheres_only,
		all_current_data,
		rebuild_skeleton_and_update_all
	};

	void update_render_data(PointsParameters& p, RenderUpdateMode mode)
	{
		std::unique_lock<std::mutex> lock(p.mutex_, std::defer_lock);
		if (p.running_)
			lock.lock();

		points_provider_->emit_connectivity_changed(*p.spheres_);
		points_provider_->emit_attribute_changed(*p.spheres_, p.spheres_position_.get());
		points_provider_->emit_attribute_changed(*p.spheres_, p.spheres_radius_.get());

		update_spheres_color(p);
		points_provider_->emit_attribute_changed(*p.spheres_, p.spheres_color_.get());

		if (mode == RenderUpdateMode::spheres_only)
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

		if (mode == RenderUpdateMode::rebuild_skeleton_and_update_all)
		{
			const auto status = p.reconstruction_->build_skeleton();
			p.skeleton_invalidated_ = status.code != ReconstructionStatus::success;
			if (status.code == ReconstructionStatus::success)
			{
				p.last_error_.clear();
				refresh_skeleton_topology_colors(p);
			}
			else
				p.last_error_ = with_message("Failed to build skeleton", status.message);
		}
		if (non_manifold_provider_ && p.skeleton_)
		{
			non_manifold_provider_->emit_connectivity_changed(*p.skeleton_);
			if (p.skeleton_position_)
				non_manifold_provider_->emit_attribute_changed(*p.skeleton_, p.skeleton_position_.get());
			if (p.skeleton_radius_)
				non_manifold_provider_->emit_attribute_changed(*p.skeleton_, p.skeleton_radius_.get());
			if (p.skeleton_edge_color_)
				non_manifold_provider_->emit_attribute_changed(*p.skeleton_, p.skeleton_edge_color_.get());
			if (p.skeleton_face_color_)
				non_manifold_provider_->emit_attribute_changed(*p.skeleton_, p.skeleton_face_color_.get());
		}
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
		if (spheres_update_thread_.joinable())
			spheres_update_thread_.join();

		active_spheres_update_ = &p;
		p.stop_requested_.store(false);
		p.render_refresh_requested_.store(false);
		{
			std::lock_guard<std::mutex> lock(p.mutex_);
			p.reconstruction_->spheres_optimizer().reset_optimization();
			p.metrics_snapshot_ = p.reconstruction_->spheres_optimizer().metrics();
			p.sphere_update_lambda_.store(p.reconstruction_->options().sqem_update_lambda_line_plane);
		}
		p.running_.store(true);

		spheres_update_thread_ = std::thread([this, active = &p]() {
			cgogn::thread_start(cgogn::max_nb_threads() - 1);
			auto start = std::chrono::high_resolution_clock::now();
			auto last_preview_refresh = start;
			SpheresOptimizerType& optimizer = active->reconstruction_->spheres_optimizer();
			SpheresOptimizerStatus status = SpheresOptimizerStatus::running;
			SpheresOptimizerMetrics final_metrics;
			try
			{
				while (status == SpheresOptimizerStatus::running && !active->stop_requested_.load())
				{
					SpheresOptimizerMetrics metrics;
					bool emit_preview_refresh = false;
					{
						std::lock_guard<std::mutex> lock(active->mutex_);
						if (active->stop_requested_.load())
							break;
						const Scalar lambda = Scalar(active->sphere_update_lambda_.load());
						status = active->reconstruction_->update_spheres_once(lambda);
						metrics = optimizer.metrics();
						if (metrics.sphere_topology_changed)
							active->skeleton_invalidated_ = true;
						const auto now = std::chrono::high_resolution_clock::now();
						if (active->preview_render_during_sphere_update_.load() &&
							std::chrono::duration_cast<std::chrono::milliseconds>(now - last_preview_refresh).count() >=
								100)
						{
							emit_preview_refresh = !active->render_refresh_requested_.exchange(true);
							if (emit_preview_refresh)
								last_preview_refresh = now;
						}
					}

					log_basic(*active, "Iteration: ", metrics.iteration, " | Spheres: ", metrics.sphere_count,
							  " | Error: ", metrics.total_error, " | Diff: ", metrics.error_difference, '\n');
					if (emit_preview_refresh)
						boost::synapse::emit<App::timer_tick>(&app_);

					if (status == SpheresOptimizerStatus::converged)
						log_basic(*active, "Auto stop: error converged after post-convergence iterations.", '\n');
					else if (status == SpheresOptimizerStatus::max_iterations)
						log_basic(*active, "Stop: reached max iterations (150).", '\n');
					else if (status == SpheresOptimizerStatus::failed)
						log_error("Sphere optimizer is missing required data.", '\n');

					if (status != SpheresOptimizerStatus::running)
						break;
					if (active->stop_requested_.load())
						break;
					if (active->slow_down_.load())
						std::this_thread::sleep_for(
							std::chrono::microseconds(1000000 / std::max<uint32>(1, active->update_rate_.load())));
					else
						std::this_thread::yield();
				}
			}
			catch (const std::exception& e)
			{
				log_error("Sphere update failed: ", e.what(), '\n');
			}
			catch (...)
			{
				log_error("Sphere update failed with an unknown exception.", '\n');
			}

			{
				std::lock_guard<std::mutex> lock(active->mutex_);
				final_metrics = optimizer.metrics();
				active->render_refresh_requested_.store(true);
				active->stop_requested_.store(false);
				active->running_.store(false);
			}
			const auto end = std::chrono::high_resolution_clock::now();
			log_basic(*active, "Sphere optimizations time: ", std::chrono::duration<Scalar>(end - start).count(), "s",
					  '\n');
			log_basic(*active, "Nb iterations: ", final_metrics.iteration, '\n');
			cgogn::thread_stop();
			boost::synapse::emit<App::timer_tick>(&app_);
		});
	}

	void stop_spheres_update(PointsParameters& p)
	{
		p.stop_requested_.store(true);
	}

	void handle_spheres_update_tick()
	{
		PointsParameters* active = active_spheres_update_;
		if (!active)
			return;

		if (!active->render_refresh_requested_.load())
			return;

		bool running;
		bool should_render;
		{
			std::lock_guard<std::mutex> lock(active->mutex_);
			running = active->running_.load();
			should_render = (running && active->preview_render_during_sphere_update_.load()) ||
							(!running && active->render_refresh_requested_.load());
			active->metrics_snapshot_ = active->reconstruction_->spheres_optimizer().metrics();
			active->render_refresh_requested_.store(false);
		}
		if (!should_render)
			return;

		if (!running && spheres_update_thread_.joinable())
			spheres_update_thread_.join();

		// The preview rebuilds the skeleton from the current spheres, like the final refresh;
		// update_render_data holds the worker mutex while it runs.
		if (active->skeleton_invalidated_ && active->skeleton_)
			clear(*active->skeleton_);
		update_render_data(*active, RenderUpdateMode::rebuild_skeleton_and_update_all);
		for (View* v : linked_views_)
			v->request_update();
		if (!running)
			active_spheres_update_ = nullptr;
	}

	void stop_and_join_spheres_update()
	{
		PointsParameters* active = active_spheres_update_;
		if (active)
			active->stop_requested_.store(true);
		if (spheres_update_thread_.joinable())
			spheres_update_thread_.join();
		if (active)
		{
			active->running_.store(false);
			active->render_refresh_requested_.store(false);
			active_spheres_update_ = nullptr;
		}
	}

	void close_event() override
	{
		stop_and_join_spheres_update();
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
				update_render_data(p, RenderUpdateMode::all_current_data);
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
		PointsParameters* active_update = active_spheres_update_;
		const bool update_running = active_update && active_update->running_.load();
		ImGui::BeginDisabled(update_running);
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
		ImGui::EndDisabled();

		if (!selected_points_)
			return;
		PointsParameters& p = points_parameters_[selected_points_];
		auto& options = p.reconstruction_->options();
		int output_verbosity = output_verbosity_index(p.output_verbosity_.load());
		if (ImGui::Combo("Output", &output_verbosity, "Mute\0Normal\0"))
			p.output_verbosity_.store(output_verbosity_from_index(output_verbosity));

		ImGui::Separator();
		if (ImGui::CollapsingHeader("Neural UDF", ImGuiTreeNodeFlags_DefaultOpen))
		{
			if (p.neural_udf_loaded_)
			{
				ImGui::TextColored(ImVec4(0, 1, 0, 1), "Model loaded: %s", p.neural_udf_model_path_.c_str());

				ImGui::Separator();
				ImGui::Text("Alpha Level Set Sampling Settings");
			}
			if (!p.last_error_.empty())
				ImGui::TextColored(ImVec4(1, 0.3f, 0.3f, 1), "%s", p.last_error_.c_str());
		}
		// Sampling
		const bool running = p.running_.load();
		const bool has_samples = p.samples_mesh_ && nb_cells<PVertex>(*p.samples_mesh_) > 0;
		if (ImGui::CollapsingHeader("Sampling", ImGuiTreeNodeFlags_DefaultOpen))
		{
			ImGui::BeginDisabled(running);
			const auto& profile_names = cgogn::geometry::reconstruction_profile_names();
			if (ImGui::BeginCombo("Profile", p.profile_.c_str()))
			{
				for (const std::string& name : profile_names)
					if (ImGui::Selectable(name.c_str(), name == p.profile_))
						reset_options_to_profile(p, name);
				ImGui::EndCombo();
			}
			if (ImGui::InputFloat("Alpha", &options.alpha, 0.001f, 0.1f, "%.4f"))
				options.alpha = std::clamp(options.alpha, float(cgogn::geometry::option_info::alpha.min),
										   float(cgogn::geometry::option_info::alpha.max));
			if (ImGui::InputFloat("Sample Radius", &options.sample_radius, 0.001f, 0.01f, "%.4f"))
				options.sample_radius =
					std::clamp(options.sample_radius, float(cgogn::geometry::option_info::sample_radius.min),
							   float(cgogn::geometry::option_info::sample_radius.max));
			if (ImGui::InputInt("KNN K", &options.knn_k, 1, 5))
				options.knn_k = std::clamp(options.knn_k, int(cgogn::geometry::option_info::knn_k.min),
										   int(cgogn::geometry::option_info::knn_k.max));
			ImGui::Checkbox("Filter samples by distance to input", &options.apply_filtering);
			ImGui::Checkbox("Recompute normals after sampling", &options.recompute_normals_after_sampling);
			if (ImGui::InputInt("Bridson candidates per parent", &options.alpha_sampling_candidates_per_parent, 1, 10))
				options.alpha_sampling_candidates_per_parent =
					std::clamp(options.alpha_sampling_candidates_per_parent,
							   int(cgogn::geometry::option_info::alpha_sampling_candidates_per_parent.min),
							   int(cgogn::geometry::option_info::alpha_sampling_candidates_per_parent.max));
			if (ImGui::InputInt("Bridson parent batch size", &options.alpha_sampling_parent_batch_size, 1, 64))
				options.alpha_sampling_parent_batch_size =
					std::clamp(options.alpha_sampling_parent_batch_size,
							   int(cgogn::geometry::option_info::alpha_sampling_parent_batch_size.min),
							   int(cgogn::geometry::option_info::alpha_sampling_parent_batch_size.max));
			if (ImGui::InputInt("Seed", &options.seed, 1, 10))
				options.seed = std::clamp(options.seed, int(cgogn::geometry::option_info::seed.min),
										  int(cgogn::geometry::option_info::seed.max));
			ImGui::EndDisabled();
			if (!running)
			{
				if (ImGui::Button("Sample UDF"))
				{
					load_alpha_samples_to_mesh(p);
					p.fitting_data_computed_ = false;
					p.spheres_initialized_ = false;
				}
				if (has_samples)
				{
					ImGui::SameLine();
					if (ImGui::Button("Clear Samples"))
					{
						if (p.samples_mesh_)
							points_provider_->clear_mesh(*p.samples_mesh_);
						p.fitting_data_computed_ = false;
						p.spheres_initialized_ = false;
					}
				}
			}
			if (p.samples_mesh_)
				ImGui::Text("Number of samples: %zu", nb_cells<PVertex>(*p.samples_mesh_));
		}

		if (!has_samples)
		{
			ImGui::TextColored(ImVec4(1, 1, 0, 1), "Please sample points first.");
			return;
		}

		// Fitting Data
		if (ImGui::CollapsingHeader("Fitting Data", ImGuiTreeNodeFlags_DefaultOpen))
		{
			ImGui::BeginDisabled(running);
			ImGui::Checkbox("Enable MAFlipPrune", &options.ma_flip_prune);
			float ma_flip_prune_alpha_factor = static_cast<float>(options.ma_flip_prune_alpha_factor);
			if (ImGui::SliderFloat("MAFlipPrune MF Alpha Factor", &ma_flip_prune_alpha_factor,
								   float(cgogn::geometry::option_info::ma_flip_prune_alpha_factor.ui_min),
								   float(cgogn::geometry::option_info::ma_flip_prune_alpha_factor.ui_max), "%.2f"))
				options.ma_flip_prune_alpha_factor = Scalar(ma_flip_prune_alpha_factor);
			ImGui::EndDisabled();
			if (!running && ImGui::Button("Compute Fitting Data"))
			{
				p.fitting_data_computed_ = false;
				p.spheres_initialized_ = false;
				compute_fitting_data(p);
				update_render_data(p, RenderUpdateMode::all_current_data);
			}
		}
		if (!p.fitting_data_computed_)
		{
			ImGui::TextColored(ImVec4(1, 1, 0, 1), "Compute fitting data to enable sphere fitting.");
			return;
		}

		// Sphere Fitting
		if (ImGui::CollapsingHeader("Sphere Fitting", ImGuiTreeNodeFlags_DefaultOpen))
		{
			ImGui::BeginDisabled(running);
			ImGui::SliderFloat("Init dilation constant", &options.init_dilation_constant,
							   float(cgogn::geometry::option_info::init_dilation_constant.ui_min),
							   float(cgogn::geometry::option_info::init_dilation_constant.ui_max), "%.4f");
			ImGui::EndDisabled();
			if (!running && ImGui::Button("Init spheres"))
			{
				std::lock_guard<std::mutex> lock(p.mutex_);
				ReconstructionStatusInfo status{ReconstructionStatus::invalid_data};
				if (p.reconstruction_)
					status = p.reconstruction_->initialize_spheres();
				p.spheres_initialized_ = status.code == ReconstructionStatus::success;
				if (status.code != ReconstructionStatus::success)
				{
					p.last_error_ = with_message("Failed to initialize spheres", status.message);
					log_error(p.last_error_, '\n');
				}
				else
				{
					p.last_error_.clear();
					p.skeleton_invalidated_ = true;
					if (p.skeleton_)
						clear(*p.skeleton_);
					set_post_init_sphere_render_state(p);
					p.metrics_snapshot_ = p.reconstruction_->spheres_optimizer().metrics();
					update_render_data(p, RenderUpdateMode::rebuild_skeleton_and_update_all);
				}
			}
			if (!p.spheres_initialized_)
				ImGui::TextColored(ImVec4(1, 1, 0, 1), "Init spheres to continue.");
		}

		// Sphere Optimization
		const bool spheres_ready = p.spheres_initialized_;
		if (spheres_ready && ImGui::CollapsingHeader("Sphere Optimization", ImGuiTreeNodeFlags_DefaultOpen))
		{
			float update_lambda = p.sphere_update_lambda_.load();
			if (ImGui::SliderFloat("Update lambda", &update_lambda,
								   float(cgogn::geometry::option_info::sqem_update_lambda_line_plane.ui_min),
								   float(cgogn::geometry::option_info::sqem_update_lambda_line_plane.ui_max), "%.6f"))
			{
				p.sphere_update_lambda_.store(update_lambda);
				options.sqem_update_lambda_line_plane = update_lambda;
			}
			if (!running && ImGui::Button("Update spheres"))
			{
				std::lock_guard<std::mutex> lock(p.mutex_);
				SpheresOptimizerType& optimizer = p.reconstruction_->spheres_optimizer();
				p.reconstruction_->update_spheres_once(Scalar(options.sqem_update_lambda_line_plane));
				p.metrics_snapshot_ = optimizer.metrics();
				if (optimizer.metrics().sphere_topology_changed)
				{
					p.skeleton_invalidated_ = true;
					if (p.skeleton_)
					{
						clear(*p.skeleton_);
					}
				}
				update_render_data(p, RenderUpdateMode::rebuild_skeleton_and_update_all);
			}
			if (!running)
			{
				if (ImGui::Button("Start spheres update"))
					start_spheres_update(p);
			}
			else if (ImGui::Button("Stop spheres update"))
				stop_spheres_update(p);

			bool slow_down = p.slow_down_.load();
			if (ImGui::Checkbox("Slow down", &slow_down))
				p.slow_down_.store(slow_down);
			if (slow_down)
			{
				int update_rate = static_cast<int>(p.update_rate_.load());
				if (ImGui::SliderInt("Update rate", &update_rate, 1, 100))
					p.update_rate_.store(static_cast<uint32>(update_rate));
			}
			bool preview = p.preview_render_during_sphere_update_.load();
			if (ImGui::Checkbox("Preview during update", &preview))
				p.preview_render_during_sphere_update_.store(preview);
			if (ImGui::Checkbox("Error as color", &p.error_as_spheres_color_))
			{
				if (!running)
				{
					std::lock_guard<std::mutex> lock(p.mutex_);
					update_render_data(p, RenderUpdateMode::spheres_only);
				}
			}
			if (ImGui::SliderFloat("Transparency", &p.spheres_transparency_, 0.0f, 1.0f))
			{
				if (!running)
				{
					std::lock_guard<std::mutex> lock(p.mutex_);
					update_render_data(p, RenderUpdateMode::all_current_data);
				}
			}

			ImGui::Separator();

			const SpheresOptimizerMetrics& metrics = p.metrics_snapshot_;
			ImGui::Text("Total error: %f", metrics.total_error / std::max<uint32>(1, metrics.sphere_count));
			ImGui::Text("Min error: %f", metrics.minimum_error);
			ImGui::Text("Max error: %f", metrics.maximum_error);
		}

		// Skeleton
		const bool skeleton_built = p.skeleton_ && !p.skeleton_invalidated_ && nb_cells<NMVertex>(*p.skeleton_) > 0;
		if (spheres_ready && ImGui::CollapsingHeader("Skeleton", ImGuiTreeNodeFlags_DefaultOpen))
		{
			if (!running && ImGui::Button("Build Skeleton"))
			{
				std::lock_guard<std::mutex> lock(p.mutex_);
				const auto status = p.reconstruction_->build_skeleton();
				p.skeleton_invalidated_ = status.code != ReconstructionStatus::success;
				if (status.code == ReconstructionStatus::success)
				{
					p.last_error_.clear();
					refresh_skeleton_topology_colors(p);
				}
				else
					p.last_error_ = with_message("Failed to build skeleton", status.message);
				update_render_data(p, RenderUpdateMode::all_current_data);
			}
		}

		// Topology Fix. While the sphere update runs, the topology is not read (the update thread writes it), so
		// the step stays visible with its parameters disabled, as before.
		// Read the existing topology object: skeleton_topology() would rebuild its data (or create it) every frame.
		const auto& topology = p.reconstruction_->state().topology;
		const bool has_tets = !running && spheres_ready && skeleton_built && topology && topology->has_tets();
		const bool topology_fix_step = spheres_ready && skeleton_built && (running || has_tets);
		if (topology_fix_step && ImGui::CollapsingHeader("Topology Fix", ImGuiTreeNodeFlags_DefaultOpen))
		{
			ImGui::BeginDisabled(running);
			ImGui::Checkbox("Residual prune", &options.residual_prune);
			ImGui::SameLine();
			int residual_mode = static_cast<int>(options.residual_prune_mode);
			if (ImGui::Combo("##ResidualPruneMode", &residual_mode,
							 "current\0keep_fully_non_manifold\0collapsible_only\0"))
				options.residual_prune_mode = static_cast<cgogn::geometry::ResidualPruneOption>(residual_mode);
			ImGui::EndDisabled();
			if (!running && ImGui::Button("Topology fix full pipeline"))
			{
				std::lock_guard<std::mutex> lock(p.mutex_);
				if (p.skeleton_invalidated_)
				{
					log_error("Topology fix requires a skeleton built from the current sphere set.\n");
				}
				else
				{
					// Same order as run(): residual sheet prune after fix_topology, if enabled.
					const auto status = p.reconstruction_->fix_topology_and_prune();
					if (status.code == ReconstructionStatus::success)
					{
						p.last_error_.clear();
						p.metrics_snapshot_ = p.reconstruction_->spheres_optimizer().metrics();
						apply_skeleton_topology_results(p);
					}
					else
					{
						p.last_error_ = with_message("[TopologyFull] score evaluation failed or topology data is invalid", status.message);
						log_error(p.last_error_, '\n');
					}
				}
			}
		}

		// Skeleton Colors. The tet-based colormaps need tets, so they are hidden once topology fix removed them all.
		if (spheres_ready && skeleton_built && !running && ImGui::CollapsingHeader("Skeleton Colors", ImGuiTreeNodeFlags_DefaultOpen))
		{
			bool first_in_row = true;
			auto next_in_row = [&first_in_row]() {
				if (!first_in_row)
					ImGui::SameLine();
				first_in_row = false;
			};
			next_in_row();
			if (ImGui::Button("Topology colors") && skeleton_colors_ready(p))
			{
				std::lock_guard<std::mutex> lock(p.mutex_);
				refresh_skeleton_topology_colors(p);
				emit_skeleton_colors_changed(p);
				show_skeleton_colors(p);
			}
			if (has_tets)
			{
				next_in_row();
				if (ImGui::Button("Boundary tet mask") && skeleton_colors_ready(p))
				{
					std::lock_guard<std::mutex> lock(p.mutex_);
					mark_boundary_tets_color(p);
					emit_skeleton_colors_changed(p);
					show_skeleton_colors(p);
				}
			}
			next_in_row();
			if (ImGui::Button("Face components (UF)") && skeleton_colors_ready(p))
			{
				std::lock_guard<std::mutex> lock(p.mutex_);
				mark_skeleton_face_components_color(p);
				emit_skeleton_colors_changed(p);
				show_skeleton_colors(p);
			}
			if (has_tets)
			{
				next_in_row();
				if (ImGui::Button("Tet face UDF colormap (Integral)"))
				{
					std::lock_guard<std::mutex> lock(p.mutex_);
					if (mark_tet_face_udf_color(p, false))
					{
						emit_skeleton_colors_changed(p);
						show_skeleton_colors(p);
					}
				}
				next_in_row();
				if (ImGui::Button("Tet face UDF colormap (Normalized)"))
				{
					std::lock_guard<std::mutex> lock(p.mutex_);
					if (mark_tet_face_udf_color(p, true))
					{
						emit_skeleton_colors_changed(p);
						show_skeleton_colors(p);
					}
				}
				next_in_row();
				if (ImGui::Button("Edge UDF score colormap (Gauss3)"))
				{
					std::lock_guard<std::mutex> lock(p.mutex_);
					if (mark_tet_edge_gauss3_color(p))
					{
						emit_skeleton_colors_changed(p);
						show_skeleton_colors(p);
					}
				}
			}
		}

		// Export
		const bool skeleton_export_available = !running && !p.skeleton_invalidated_ && p.skeleton_ &&
											   p.skeleton_position_ && p.skeleton_radius_ &&
											   nb_cells<NMVertex>(*p.skeleton_) > 0;
		if (skeleton_export_available && ImGui::CollapsingHeader("Export", ImGuiTreeNodeFlags_DefaultOpen))
		{
			if (ImGui::Button("Export skeleton PLY"))
			{
				const std::string filename =
					pfd::save_file("Export skeleton PLY", "skeleton.ply", {"PLY", "*.ply"}).result();
				if (!filename.empty())
				{
					std::lock_guard<std::mutex> lock(p.mutex_);
					if (!(p.reconstruction_->skeleton() &&
						  geometry::export_skeleton_ply(
							  *p.reconstruction_->skeleton(), p.reconstruction_->skeleton_position_attribute(),
							  p.reconstruction_->skeleton_radius_attribute(), filename)))
						log_error("Failed to export skeleton PLY: ", filename, '\n');
				}
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

	POINTS* selected_points_ = nullptr;
	std::map<POINTS*, PointsParameters> points_parameters_;
	std::shared_ptr<boost::synapse::connection> timer_connection_;
	std::thread spheres_update_thread_;
	PointsParameters* active_spheres_update_ = nullptr;
};

} // namespace ui

} // namespace cgogn

#endif // CGOGN_MODULE_UDF_TRAINING_H_
