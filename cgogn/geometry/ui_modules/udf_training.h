#ifndef CGOGN_MODULE_UDF_TRAINING_H_
#define CGOGN_MODULE_UDF_TRAINING_H_

#include <cgogn/ui/app.h>
#include <cgogn/ui/module.h>

#include <cgogn/core/ui_modules/mesh_provider.h>
#include <cgogn/geometry/algos/udf/reconstruction.h>
#include <cgogn/geometry/types/vector_traits.h>
#include <cgogn/ui/portable-file-dialogs.h>

#include <cgogn/rendering/ui_modules/point_cloud_render.h>
#include <cgogn/rendering/ui_modules/surface_render.h>

#include <GLFW/glfw3.h>
#include <algorithm>
#include <chrono>
#include <cmath>
#include <filesystem>
#include <iostream>
#include <map>
#include <memory>
#include <mutex>
#include <set>
#include <stdexcept>
#include <string>
#include <thread>
#include <torch/torch.h>
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
	using SkeletonTopologyType = geometry::SkeletonTopology<POINTS, NONMANIFOLD>;
	using ReconstructionType = geometry::UDFReconstruction<SURFACE, POINTS, NONMANIFOLD, RaySamplerTag>;
	using ReconstructionStatus = typename ReconstructionType::Status;

	enum NeuralModelType : uint32
	{
		NEURAL_MODEL_UDF,
		NEURAL_MODEL_MF
	};
	enum OutputVerbosity : uint32
	{
		OUTPUT_MUTE,
		OUTPUT_NORMAL
	};

private:
	struct PointsParameters;

	struct PointsParameters
	{
		bool initialized_ = false;
		bool fitting_data_computed_ = false; // UI workflow gate for the sphere-fitting controls.
		std::unique_ptr<ReconstructionType> reconstruction_;
		// Input Points
		POINTS* points_ = nullptr;
		std::shared_ptr<PAttribute<Vec3>> position_ = nullptr;

		// Neural model information shown in the UI; the model itself is owned by ReconstructionType.
		bool neural_udf_loaded_ = false;
		std::string neural_udf_model_path_ = "";

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
		bool skeleton_invalidated_ = false;
		std::shared_ptr<NMAttribute<Vec3>> skeleton_position_ = nullptr;
		std::shared_ptr<NMAttribute<Scalar>> skeleton_radius_ = nullptr;
		std::shared_ptr<NMAttribute<std::set<std::size_t>>> incident_tets_ = nullptr;

		std::shared_ptr<NMAttribute<Vec3>> skeleton_face_color_ = nullptr;
		std::shared_ptr<NMAttribute<Vec3>> skeleton_face_component_color_ = nullptr;
		std::shared_ptr<NMAttribute<Vec3>> skeleton_edge_color_ = nullptr;

		bool error_as_spheres_color_ = false;
		float32 spheres_transparency_ = 0.5f;

		// UI state
		OutputVerbosity output_verbosity_ = OUTPUT_NORMAL;

		// Threading
		std::mutex mutex_;
		bool running_ = false;
		bool stopping_ = false;
		bool preview_render_during_sphere_update_ = true;
		bool pending_full_refresh_after_stop_ = false;
		bool slow_down_ = true;
		uint32 update_rate_ = 20;

	};

public:
	UDFTraining(const App& app) : ViewModule(app, "UDFTraining")
	{
	}

	~UDFTraining()
	{
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
		return p.output_verbosity_ != OUTPUT_MUTE;
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
	void log_error(const PointsParameters& p, Args&&... args) const
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
			throw std::runtime_error("UDFTraining non-manifold provider is not initialized before set_selected_points().");
		selected_points_ = &p;
		init_points_data(p);
		PointsParameters& params = points_parameters_[selected_points_];
		if (params.reconstruction_)
			params.reconstruction_->data().surface = selected_surface_;
	}

	void load_neural_udf_model(POINTS& points, const std::string& model_path, NeuralModelType model_type)
	{
		PointsParameters& p = points_parameters_[&points];
		init_points_data(points);
		if (!std::filesystem::exists(model_path))
		{
			log_error(p, "Neural UDF model file does not exist: ", model_path, '\n');
			return;
		}
		try
		{
			log_basic(p, "Loading Neural UDF model from: ", model_path, '\n');
			const auto type = model_type == NEURAL_MODEL_MF ? ReconstructionType::NeuralModelType::mf
													: ReconstructionType::NeuralModelType::udf;
			p.neural_udf_loaded_ = p.reconstruction_ &&
				p.reconstruction_->load_neural_model(model_path, type) == ReconstructionStatus::success;
			if (!p.neural_udf_loaded_)
				throw std::runtime_error("The reconstruction model loader rejected the model.");
			p.neural_udf_model_path_ = model_path;
			log_basic(p, "Loaded neural UDF model from: ", model_path, '\n');
		}
		catch (const std::exception& e)
		{
			log_error(p, "Error loading Neural UDF model: ", e.what(), '\n');
			p.neural_udf_loaded_ = false;
		}
	}


	void load_alpha_samples_to_mesh(PointsParameters& p)
	{
		if (!p.reconstruction_)
			return;
		const ReconstructionStatus status = p.reconstruction_->sample_alpha_level_set();
		if (status != ReconstructionStatus::success)
		{
			log_error(p, "Alpha level set sampling failed.", '\n');
			return;
		}
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

	void apply_sampling_preprocess_filtering(PointsParameters& p)
	{
		if (!p.reconstruction_ || p.running_)
			return;
		const uint32 before = nb_cells<PVertex>(*p.samples_mesh_);
		const ReconstructionStatus status = p.reconstruction_->apply_sampling_filtering();
		if (status != ReconstructionStatus::success)
		{
			log_error(p, "Sampling filtering failed.", '\n');
			return;
		}
		const uint32 after = nb_cells<PVertex>(*p.samples_mesh_);
		if (after == before)
			return;
		p.fitting_data_computed_ = false;
		if (p.samples_color_)
			p.samples_color_->fill(Vec4(0.0, 0.0, 0.0, 1.0));
		points_provider_->emit_connectivity_changed(*p.samples_mesh_);
		points_provider_->emit_attribute_changed(*p.samples_mesh_, p.samples_position_.get());
		points_provider_->emit_attribute_changed(*p.samples_mesh_, p.samples_normal_.get());
		if (p.samples_color_)
			points_provider_->emit_attribute_changed(*p.samples_mesh_, p.samples_color_.get());
		log_basic(p, "Sampling filtering applied: ", before, " -> ", after, " points.", '\n');
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
					}
					update_render_data(p, false, true, true);
					request_linked_views_update();
					p.pending_full_refresh_after_stop_ = false;
				}
			}
		});
		const OutputVerbosity startup_output_verbosity =
			selected_points_ ? points_parameters_[selected_points_].output_verbosity_ : OUTPUT_NORMAL;
		if (torch::cuda::is_available())
		{
			log_basic(startup_output_verbosity, "CUDA is available! Using GPU device 0.", '\n');
		}
		else
		{
			log_basic(startup_output_verbosity, "CUDA is not available! Using the CPU.", '\n');
		}
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
		p.skeleton_face_component_color_ = get_or_add_attribute<Vec3, NMFace>(*p.skeleton_, "face_component_color");
		p.skeleton_edge_color_ = get_or_add_attribute<Vec3, NMEdge>(*p.skeleton_, "color");

		typename ReconstructionType::Data data{selected_surface_, p.points_, p.samples_mesh_, p.spheres_, p.skeleton_};
		p.reconstruction_ = std::make_unique<ReconstructionType>(data);
		p.reconstruction_->options().ma_flip_prune = false;
		p.reconstruction_->options().ray_sampler_batch_size = 8192;
		const ReconstructionStatus status = p.reconstruction_->initialize_data();
		p.initialized_ = status == ReconstructionStatus::success;
	}

	void compute_fitting_data(PointsParameters& p)
	{
		if (!p.reconstruction_ || p.fitting_data_computed_)
			return;
		ReconstructionType& reconstruction = *p.reconstruction_;
		if (reconstruction.build_kdtree_and_normals() != ReconstructionStatus::success ||
			reconstruction.compute_fitting_primitives() != ReconstructionStatus::success ||
			reconstruction.compute_initial_medial_axis() != ReconstructionStatus::success)
		{
			log_error(p, "Failed to compute fitting data.", '\n');
			return;
		}
		p.fitting_data_computed_ = true;
	}

	void init_spheres(PointsParameters& p)
	{
		if (!p.reconstruction_ ||
			p.reconstruction_->initialize_spheres() != ReconstructionStatus::success)
		{
			log_error(p, "Failed to initialize spheres.", '\n');
			return;
		}
		p.skeleton_invalidated_ = true;
		if (p.skeleton_)
			clear(*p.skeleton_);
		if (!p.running_)
			set_post_init_sphere_render_state(p);
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
					color_map((*p.spheres_error_)[v_index], p.reconstruction_->spheres_optimizer().metrics().minimum_error,
											  p.reconstruction_->spheres_optimizer().metrics().maximum_error, p.spheres_transparency_);
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

	SpheresOptimizerType& optimizer_for(PointsParameters& p) { return p.reconstruction_->spheres_optimizer(); }
	SkeletonTopologyType& skeleton_topology_for(PointsParameters& p) { return p.reconstruction_->skeleton_topology(); }
	void log_skeleton_topology_summary(const PointsParameters& p, const char* prefix) const
	{
		if (!p.reconstruction_) return;
		const auto& metrics = p.reconstruction_->skeleton_topology().metrics();
		log_basic(p, prefix, " rounds=", metrics.topology_rounds,
				  " boundary_removed_tets=", metrics.boundary_removed_tets,
				  " simple_removed_tets=", metrics.simple_removed_tets,
				  " nonsimple_removed_tets=", metrics.nonsimple_removed_tets,
				  " removed_faces=", metrics.removed_faces, " removed_edges=", metrics.removed_edges,
				  " removed_tets=", metrics.removed_tets, " remaining_tets=", metrics.tet_count,
				  " repaired_regions=", metrics.repaired_dense_regions, " added_spheres=", metrics.added_spheres.size(),
				  " added_faces=", metrics.added_faces, " residual_removed_sheets=", metrics.residual_removed_sheets,
				  " residual_removed_faces=", metrics.residual_removed_faces, '\n');
	}
	void apply_skeleton_topology_results(PointsParameters& p)
	{
		if (!p.reconstruction_) return;
		const auto& metrics = p.reconstruction_->skeleton_topology().metrics();
		for (PVertex sphere : metrics.added_spheres)
		{
			const uint32 id = p.spheres_ ? index_of(*p.spheres_, sphere) : INVALID_INDEX;
			if (id == INVALID_INDEX) continue;
			if (p.spheres_color_) (*p.spheres_color_)[id] = Vec4(0.95, 0.25, 0.15, 1.0);
			if (p.spheres_cluster_color_) (*p.spheres_cluster_color_)[id] = Vec4(0.95, 0.25, 0.15, 1.0);
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
			if (p.spheres_cluster_color_) points_provider_->emit_attribute_changed(*p.spheres_, p.spheres_cluster_color_.get());
		}
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
			const ReconstructionStatus status = p.reconstruction_->build_skeleton();
			p.skeleton_invalidated_ = status != ReconstructionStatus::success;
			if (status == ReconstructionStatus::success)
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
						status = optimizer.update_once(Scalar(p.reconstruction_->options().sqem_update_lambda_line_plane));
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

				log_basic(p, "Iteration: ", optimizer.metrics().iteration, " | Spheres: ", optimizer.metrics().sphere_count,
						  " | Error: ", optimizer.metrics().total_error, " | Diff: ", optimizer.metrics().error_difference, '\n');
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
			log_basic(p, "Nb iterations: ", optimizer.metrics().iteration, '\n');
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
		auto& options = p.reconstruction_->options();
		int output_verbosity = output_verbosity_index(p.output_verbosity_);
		if (ImGui::Combo("Output", &output_verbosity, "Mute\0Normal\0"))
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
			ImGui::InputFloat("Alpha", &options.alpha, 0.001f, 0.1f, "%.4f");
			ImGui::InputFloat("Sample Radius", &options.sample_radius, 0.001f, 0.01f, "%.4f");
			ImGui::InputInt("Eval Batch Size", &options.batch_size, 256, 1024);
			ImGui::InputInt("Ray Batch Size", &options.ray_sampler_batch_size, 256, 2048);
			ImGui::InputInt("Max Iterations", &options.udf_max_iterations, 1000, 8000);
			ImGui::InputFloat("Tolerance", &options.tolerance, 0.0f, 0.0f, "%.6f");
			ImGui::InputInt("KNN K", &options.knn_k, 1, 5);
			ImGui::Checkbox("Recompute Normals In Fitting Data", &options.recompute_normals_after_sampling);
			if (ImGui::Button("Apply Sampling Filtering"))
			{
				std::lock_guard<std::mutex> lock(p.mutex_);
				options.apply_filtering = true;
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
			ImGui::Checkbox("Enable MAFlipPrune", &options.ma_flip_prune);
			float ma_flip_prune_alpha_factor = static_cast<float>(options.ma_flip_prune_alpha_factor);
			if (ImGui::SliderFloat("MAFlipPrune MF Alpha Factor", &ma_flip_prune_alpha_factor, 1.0f, 5.0f, "%.2f"))
				options.ma_flip_prune_alpha_factor = Scalar(ma_flip_prune_alpha_factor);
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
					ImGui::SliderFloat("Init dilation constant", &options.init_dilation_constant, 0.001f, 0.01f, "%.4f");
					if (ImGui::Button("Init spheres"))
					{
						std::lock_guard<std::mutex> lock(p.mutex_);
						init_spheres(p);
						update_render_data(p, false, true, true);
					}
					ImGui::SliderFloat("Update lambda", &options.sqem_update_lambda_line_plane, 0.0f, 4.0f, "%.6f");
					if (ImGui::Button("Update spheres"))
					{
						if (!p.running_)
						{
							std::lock_guard<std::mutex> lock(p.mutex_);
							SpheresOptimizerType& optimizer = optimizer_for(p);
							optimizer.update_once(Scalar(options.sqem_update_lambda_line_plane));
							if (optimizer.metrics().sphere_topology_changed)
							{
								p.skeleton_invalidated_ = true;
								if (p.skeleton_)
								{
									clear(*p.skeleton_);
								}
							}
												update_render_data(p, false, true, true);
						}
					}


					if (ImGui::Button("Build Skeleton"))
					{
						if (!p.running_)
						{
							std::lock_guard<std::mutex> lock(p.mutex_);
							const ReconstructionStatus status = p.reconstruction_->build_skeleton();
							p.skeleton_invalidated_ = status != ReconstructionStatus::success;
							if (status == ReconstructionStatus::success)
							{
								log_basic(p, "[SkeletonBuild] tets=", skeleton_topology_for(p).metrics().tet_count, '\n');
								refresh_skeleton_topology_colors(p);
							}
							update_render_data(p, false, true, false);
						}
					}
					const bool skeleton_export_available =
						!p.skeleton_invalidated_ && p.skeleton_ && p.skeleton_position_ && p.skeleton_radius_ &&
						nb_cells<NMVertex>(*p.skeleton_) > 0;
					const bool can_export_skeleton = !p.running_ && skeleton_export_available;
					ImGui::TextDisabled(skeleton_export_available ? "ready" : "missing");
					if (!can_export_skeleton)
						ImGui::BeginDisabled();
					if (ImGui::Button("Export skeleton PLY"))
					{
						const std::string filename = pfd::save_file("Export skeleton PLY", "skeleton.ply", {"PLY", "*.ply"}).result();
						if (!filename.empty())
						{
							std::lock_guard<std::mutex> lock(p.mutex_);
							if (!p.reconstruction_->export_skeleton_ply(filename, false))
								log_error(p, "Failed to export skeleton PLY: ", filename, '\n');
						}
					}
					if (!can_export_skeleton)
						ImGui::EndDisabled();
					if (ImGui::Button("Face components (UF)"))
					{
						if (!p.running_)
						{
							std::lock_guard<std::mutex> lock(p.mutex_);
							if (p.reconstruction_->prepare_face_components() == ReconstructionStatus::success &&
								non_manifold_provider_ && p.skeleton_ &&
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
							const ReconstructionStatus status = p.reconstruction_->fix_topology();
							if (status == ReconstructionStatus::success)
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

					const SpheresOptimizerMetrics& metrics = optimizer_for(p).metrics();
					ImGui::Text("Total error: %f", metrics.total_error / std::max<uint32>(1, metrics.sphere_count));
					ImGui::Text("Min error: %f", metrics.minimum_error);
					ImGui::Text("Max error: %f", metrics.maximum_error);

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

	POINTS* selected_points_ = nullptr;
	std::map<POINTS*, PointsParameters> points_parameters_;
	std::shared_ptr<boost::synapse::connection> timer_connection_;

};

} // namespace ui

} // namespace cgogn

#endif // CGOGN_MODULE_UDF_TRAINING_H_
