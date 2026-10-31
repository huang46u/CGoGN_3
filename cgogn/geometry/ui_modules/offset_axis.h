/*******************************************************************************
 * CGoGN: Combinatorial and Geometric modeling with Generic N-dimensional Maps  *
 * Copyright (C), IGG Group, ICube, University of Strasbourg, France            *
 *                                                                              *
 * This library is free software; you can redistribute it and/or modify it      *
 * under the terms of the GNU Lesser General Public License as published by the *
 * Free Software Foundation; either version 2.1 of the License, or (at your     *
 * option) any later version.                                                   *
 *                                                                              *
 * This library is distributed in the hope that it will be useful, but WITHOUT *
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or        *
 * FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License *
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
#ifndef CGOGN_MODULE_OFFSET_AXIS_H_
#define CGOGN_MODULE_OFFSET_AXIS_H_

#include <cgogn/ui/app.h>
#include <cgogn/ui/module.h>

#include <cgogn/core/ui_modules/mesh_provider.h>
#include <cgogn/core/utils/thread.h>
#if !defined(CGOGN_OFFSET_AXIS_WITH_TORCH) || CGOGN_OFFSET_AXIS_WITH_TORCH
#include <cgogn/geometry/algos/offset_axis/neural/neural_distance_field.h>
#endif
#include <cgogn/geometry/algos/offset_axis/neural_model_type.h>
#include <cgogn/geometry/algos/offset_axis/alpha_medial.h>
#include <cgogn/geometry/algos/offset_axis/alpha_sampling.h>
#include <cgogn/geometry/algos/offset_axis/distance_field.h>
#include <cgogn/geometry/algos/offset_axis/input_geometry.h>
#include <cgogn/geometry/algos/offset_axis/topology_scores.h>
#include <cgogn/geometry/algos/offset_axis/reconstruction_options.h>
#include <cgogn/geometry/algos/offset_axis/sample_processing.h>
#include <cgogn/geometry/algos/offset_axis/sampling_filter.h>
#include <cgogn/geometry/algos/offset_axis/skeleton_face_components.h>
#include <cgogn/geometry/algos/offset_axis/skeleton_topology.h>
#include <cgogn/geometry/algos/offset_axis/spheres_optimizer.h>
#include <cgogn/io/surface/skeleton_ply.h>
#include <cgogn/geometry/types/line_quadric.h>
#include <cgogn/geometry/types/spherical_quadric.h>
#include <cgogn/geometry/types/vector_traits.h>
#include <cgogn/ui/portable-file-dialogs.h>

#include <cgogn/rendering/ui_modules/point_cloud_render.h>
#include <cgogn/rendering/ui_modules/surface_render.h>

#include <libacc/kd_tree.h>

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
#include <utility>
#include <vector>

namespace cgogn
{

namespace ui
{

using geometry::Scalar;
using geometry::Vec3;
using geometry::Vec4;

template <typename SURFACE, typename POINTS, typename NONMANIFOLD, template <typename, typename> class FIELD>
class OffsetAxis : public ViewModule
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
	using SpheresOptimizerType = cgogn::geometry::offset_axis::SpheresOptimizer<POINTS>;
	using SpheresOptimizerMetrics = typename SpheresOptimizerType::Metrics;
	using SpheresOptimizerStatus = typename SpheresOptimizerType::Status;
	using TopologyType = cgogn::geometry::offset_axis::SkeletonTopology<POINTS, NONMANIFOLD>;
	using Field = FIELD<SURFACE, POINTS>;

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
		// Reconstruction
		cgogn::geometry::offset_axis::ReconstructionOptions options_;
		cgogn::geometry::offset_axis::InputGeometry<SURFACE, POINTS> input_;
		Field field_{input_, options_};
		std::unique_ptr<acc::KDTree<3, uint32>> samples_kdtree_, samples_ma_kdtree_;
		std::vector<PVertex> samples_kdtree_vertices_, samples_ma_kdtree_vertices_;
		std::unique_ptr<SpheresOptimizerType> optimizer_;
		std::unique_ptr<TopologyType> topology_;

		// Input
		SURFACE* surface_ = nullptr;
		POINTS* points_ = nullptr;
		std::shared_ptr<PAttribute<Vec3>> position_ = nullptr;
		std::shared_ptr<PAttribute<Vec3>> normal_ = nullptr;
		std::shared_ptr<PAttribute<std::vector<PVertex>>> knn_ = nullptr;

		// Neural model information shown in the UI; the model itself is owned by field_.
		bool neural_udf_loaded_ = false;
		std::string neural_udf_model_path_ = "";

		// Reason of the last failed stage, shown in red in the Neural UDF panel; cleared when the stage succeeds.
		std::string last_error_;

		// Sampling & Fitting Data
		POINTS* samples_mesh_ = nullptr;
		std::shared_ptr<PAttribute<Vec3>> samples_position_ = nullptr;
		std::shared_ptr<PAttribute<Vec3>> samples_normal_ = nullptr;
		std::shared_ptr<PAttribute<Scalar>> samples_area_ = nullptr;
		std::shared_ptr<PAttribute<std::vector<PVertex>>> samples_knn_ = nullptr;
		std::shared_ptr<PAttribute<geometry::Spherical_Quadric>> samples_quadric_ = nullptr;
		std::shared_ptr<PAttribute<geometry::Line_Quadric>> samples_line_quadric_ = nullptr;
		std::shared_ptr<PAttribute<Scalar>> samples_error_ = nullptr;

		// Medial Axis (on samples)
		std::shared_ptr<PAttribute<Vec3>> samples_ma_position_ = nullptr;
		std::shared_ptr<PAttribute<Scalar>> samples_ma_radius_ = nullptr;
		std::shared_ptr<PAttribute<PVertex>> samples_ma_secondary_ = nullptr;

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
		std::shared_ptr<PAttribute<Scalar>> spheres_error_raw_ = nullptr;
		std::shared_ptr<PAttribute<std::vector<PVertex>>> spheres_cluster_ = nullptr;
		std::shared_ptr<PAttribute<Scalar>> spheres_cluster_area_ = nullptr;
		std::shared_ptr<PAttribute<std::vector<PVertex>>> spheres_neighbors_ = nullptr;
		std::shared_ptr<PAttribute<NMVertex>> spheres_skeleton_vertex_ = nullptr;

		// Skeleton
		NONMANIFOLD* skeleton_ = nullptr;
		std::atomic<bool> skeleton_invalidated_{false};
		std::shared_ptr<NMAttribute<Vec3>> skeleton_position_ = nullptr;
		std::shared_ptr<NMAttribute<Scalar>> skeleton_radius_ = nullptr;
		std::shared_ptr<NMAttribute<std::set<std::size_t>>> incident_tets_ = nullptr;
		std::shared_ptr<NMAttribute<PVertex>> skeleton_source_sphere_ = nullptr;
		std::shared_ptr<NMAttribute<uint32>> skeleton_edge_degree_ = nullptr;
		std::shared_ptr<NMAttribute<Scalar>> skeleton_face_score_ = nullptr;
		std::shared_ptr<NMAttribute<Scalar>> skeleton_edge_score_ = nullptr;

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
	OffsetAxis(const App& app) : ViewModule(app, "OffsetAxis")
	{
	}

	~OffsetAxis() override
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
			throw std::runtime_error("OffsetAxis surface provider is not initialized before set_selected_surface().");
		selected_surface_ = &s;
		if (selected_points_ && points_parameters_[selected_points_].initialized_)
			bind_input(points_parameters_[selected_points_], &s);
	}

	void set_selected_points(POINTS& p)
	{
		if (!points_provider_)
			throw std::runtime_error("OffsetAxis points provider is not initialized before set_selected_points().");
		if (!non_manifold_provider_)
			throw std::runtime_error(
				"OffsetAxis non-manifold provider is not initialized before set_selected_points().");
		selected_points_ = &p;
		init_points_data(p);
		bind_input(points_parameters_[selected_points_], selected_surface_);
	}

	bool load_neural_udf_model(POINTS& points, const std::string& model_path,
							   cgogn::geometry::offset_axis::NeuralModelType model_type)
	{
		PointsParameters& p = points_parameters_[&points];
		init_points_data(points);
		p.neural_udf_loaded_ = false;
		p.neural_udf_model_path_.clear();
		try
		{
			if (!std::filesystem::exists(model_path))
			{
				p.last_error_ = "Neural UDF model file does not exist: " + model_path;
				log_error(p.last_error_, '\n');
				return false;
			}
			log_basic(p, "Loading Neural UDF model from: ", model_path, '\n');
			p.options_.neural_model_type = model_type;
			p.neural_udf_loaded_ = p.field_.load_model(model_path, model_type);
			if (!p.neural_udf_loaded_)
			{
				p.last_error_ =
					with_message("Error loading Neural UDF model: loader rejected the model", p.field_.load_error());
				log_error(p.last_error_, '\n');
				return false;
			}
			p.neural_udf_model_path_ = model_path;
			p.last_error_.clear();
			p.options_ = cgogn::geometry::offset_axis::default_options(model_type);
			p.sphere_update_lambda_.store(p.options_.sqem_update_lambda_line_plane);
			log_basic(p, "Loaded neural UDF model from: ", model_path, '\n');
			return true;
		}
		catch (const std::exception& e)
		{
			p.last_error_ = std::string("Error loading Neural UDF model: ") + e.what();
			log_error(p.last_error_, '\n');
			p.neural_udf_loaded_ = false;
			return false;
		}
		catch (...)
		{
			p.last_error_ = "Error loading Neural UDF model: unknown exception.";
			log_error(p.last_error_, '\n');
			p.neural_udf_loaded_ = false;
			return false;
		}
	}

	void load_alpha_samples_to_mesh(PointsParameters& p)
	{
		std::string reason;
		if (!sample(p, reason))
		{
			p.last_error_ = with_message("Alpha level set sampling failed", reason);
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
#if !defined(CGOGN_OFFSET_AXIS_WITH_TORCH) || CGOGN_OFFSET_AXIS_WITH_TORCH
		if (geometry::offset_axis::neural_cuda_available())
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
		p.position_ = get_or_add_attribute<Vec3, PVertex>(m, "position");
		p.normal_ = get_or_add_attribute<Vec3, PVertex>(m, "normal");
		p.knn_ = get_or_add_attribute<std::vector<PVertex>, PVertex>(m, "knn");

		const std::string sample_name = points_provider_->mesh_name(m) + "_samples";
		p.samples_mesh_ = points_provider_->has_mesh(sample_name) ? points_provider_->mesh(sample_name)
																  : points_provider_->add_mesh(sample_name);
		points_provider_->clear_mesh(*p.samples_mesh_);
		p.samples_position_ = get_or_add_attribute<Vec3, PVertex>(*p.samples_mesh_, "position");
		p.samples_normal_ = get_or_add_attribute<Vec3, PVertex>(*p.samples_mesh_, "normal");
		p.samples_sphere_ = get_or_add_attribute<PVertex, PVertex>(*p.samples_mesh_, "sphere");
		p.samples_color_ = get_or_add_attribute<Vec4, PVertex>(*p.samples_mesh_, "color");
		p.samples_normal_color_ = get_or_add_attribute<Vec4, PVertex>(*p.samples_mesh_, "normal_color");
		p.samples_area_ = get_or_add_attribute<Scalar, PVertex>(*p.samples_mesh_, "area");
		p.samples_knn_ = get_or_add_attribute<std::vector<PVertex>, PVertex>(*p.samples_mesh_, "knn");
		p.samples_quadric_ = get_or_add_attribute<geometry::Spherical_Quadric, PVertex>(*p.samples_mesh_, "quadric");
		p.samples_line_quadric_ =
			get_or_add_attribute<geometry::Line_Quadric, PVertex>(*p.samples_mesh_, "line_quadric");
		p.samples_ma_position_ = get_or_add_attribute<Vec3, PVertex>(*p.samples_mesh_, "ma_position");
		p.samples_ma_radius_ = get_or_add_attribute<Scalar, PVertex>(*p.samples_mesh_, "ma_radius");
		p.samples_ma_secondary_ = get_or_add_attribute<PVertex, PVertex>(*p.samples_mesh_, "ma_secondary_vertex");
		p.samples_error_ = get_or_add_attribute<Scalar, PVertex>(*p.samples_mesh_, "error");

		const std::string sphere_name = points_provider_->mesh_name(m) + "_spheres";
		p.spheres_ = points_provider_->has_mesh(sphere_name) ? points_provider_->mesh(sphere_name)
															 : points_provider_->add_mesh(sphere_name);
		p.spheres_position_ = get_or_add_attribute<Vec3, PVertex>(*p.spheres_, "position");
		p.spheres_radius_ = get_or_add_attribute<Scalar, PVertex>(*p.spheres_, "radius");
		p.spheres_color_ = get_or_add_attribute<Vec4, PVertex>(*p.spheres_, "color");
		p.spheres_cluster_color_ = get_or_add_attribute<Vec4, PVertex>(*p.spheres_, "cluster_color");
		p.spheres_error_ = get_or_add_attribute<Scalar, PVertex>(*p.spheres_, "error");
		p.spheres_error_raw_ = get_or_add_attribute<Scalar, PVertex>(*p.spheres_, "error_not_normalized");
		p.spheres_cluster_ = get_or_add_attribute<std::vector<PVertex>, PVertex>(*p.spheres_, "cluster");
		p.spheres_cluster_area_ = get_or_add_attribute<Scalar, PVertex>(*p.spheres_, "cluster_area");
		p.spheres_neighbors_ = get_or_add_attribute<std::vector<PVertex>, PVertex>(*p.spheres_, "neighbor_clusters");
		p.spheres_skeleton_vertex_ = get_or_add_attribute<NMVertex, PVertex>(*p.spheres_, "skeleton_vertex");

		const std::string skeleton_name = points_provider_->mesh_name(m) + "_skeleton";
		p.skeleton_ = non_manifold_provider_->has_mesh(skeleton_name) ? non_manifold_provider_->mesh(skeleton_name)
																	  : non_manifold_provider_->add_mesh(skeleton_name);
		p.skeleton_position_ = get_or_add_attribute<Vec3, NMVertex>(*p.skeleton_, "position");
		p.skeleton_radius_ = get_or_add_attribute<Scalar, NMVertex>(*p.skeleton_, "radius");
		p.incident_tets_ = get_or_add_attribute<std::set<std::size_t>, NMFace>(*p.skeleton_, "incident_tets");
		p.skeleton_face_color_ = get_or_add_attribute<Vec3, NMFace>(*p.skeleton_, "color");
		p.skeleton_edge_color_ = get_or_add_attribute<Vec3, NMEdge>(*p.skeleton_, "color");
		p.skeleton_source_sphere_ = get_or_add_attribute<PVertex, NMVertex>(*p.skeleton_, "source_sphere");
		p.skeleton_edge_degree_ = get_or_add_attribute<uint32, NMEdge>(*p.skeleton_, "degree");
		p.skeleton_face_score_ = get_or_add_attribute<Scalar, NMFace>(*p.skeleton_, "udf_score");
		p.skeleton_edge_score_ = get_or_add_attribute<Scalar, NMEdge>(*p.skeleton_, "udf_score");

		typename SpheresOptimizerType::Data optimizer_data{*p.samples_mesh_, *p.spheres_, p.samples_kdtree_,
														   p.samples_kdtree_vertices_};
		optimizer_data.sample_position = p.samples_position_.get();
		optimizer_data.sample_area = p.samples_area_.get();
		optimizer_data.sample_knn = p.samples_knn_.get();
		optimizer_data.sample_quadric = p.samples_quadric_.get();
		optimizer_data.sample_line_quadric = p.samples_line_quadric_.get();
		optimizer_data.sample_ma_position = p.samples_ma_position_.get();
		optimizer_data.sample_ma_radius = p.samples_ma_radius_.get();
		optimizer_data.sample_ma_secondary_vertex = p.samples_ma_secondary_.get();
		optimizer_data.sample_sphere = p.samples_sphere_.get();
		optimizer_data.sample_error = p.samples_error_.get();
		optimizer_data.sphere_position = p.spheres_position_.get();
		optimizer_data.sphere_radius = p.spheres_radius_.get();
		optimizer_data.sphere_cluster = p.spheres_cluster_.get();
		optimizer_data.sphere_cluster_area = p.spheres_cluster_area_.get();
		optimizer_data.sphere_cluster_color = p.spheres_cluster_color_.get();
		optimizer_data.sphere_neighbors = p.spheres_neighbors_.get();
		optimizer_data.sphere_error = p.spheres_error_.get();
		optimizer_data.sphere_error_not_normalized = p.spheres_error_raw_.get();
		p.optimizer_ = std::make_unique<SpheresOptimizerType>(optimizer_data);

		typename TopologyType::Data topology_data{*p.skeleton_, *p.optimizer_};
		topology_data.sphere_skeleton_vertex = p.spheres_skeleton_vertex_.get();
		topology_data.skeleton_position = p.skeleton_position_.get();
		topology_data.skeleton_radius = p.skeleton_radius_.get();
		topology_data.skeleton_source_sphere = p.skeleton_source_sphere_.get();
		topology_data.face_incident_tets = p.incident_tets_.get();
		topology_data.edge_degree = p.skeleton_edge_degree_.get();
		topology_data.face_score = p.skeleton_face_score_.get();
		topology_data.edge_score = p.skeleton_edge_score_.get();
		p.topology_ = std::make_unique<TopologyType>(topology_data);

		p.sphere_update_lambda_.store(p.options_.sqem_update_lambda_line_plane);
		bind_input(p, selected_surface_);
		p.initialized_ = true;
		p.last_error_.clear();
	}

	// Points the input caches at the surface (the selected surface may change) and clears them.
	static void bind_input(PointsParameters& p, SURFACE* surface)
	{
		p.surface_ = surface;
		p.input_.bind(surface, p.points_, p.position_, p.normal_);
		p.input_.invalidate();
	}

	static void build_samples_kdtrees(PointsParameters& p)
	{
		geometry::offset_axis::build_point_kdtree(*p.samples_mesh_, *p.samples_position_, p.samples_kdtree_,
												  p.samples_kdtree_vertices_);
		geometry::offset_axis::build_medial_kdtree(*p.samples_mesh_, *p.samples_ma_position_, *p.samples_ma_radius_,
												   p.samples_ma_kdtree_, p.samples_ma_kdtree_vertices_);
	}

	// Alpha level set sampling, then (apply_filtering) the removal of the samples far from the input. reason is set
	// when the filtering removes every sample.
	static bool sample(PointsParameters& p, std::string& reason)
	{
		p.input_.invalidate();
		p.input_.build_input_kdtree();
		const geometry::offset_axis::PrepareInputs<SURFACE, POINTS> inputs{p.surface_, p.points_, p.position_,
																			p.normal_, p.knn_};
		if (p.field_.prepare(inputs) != geometry::offset_axis::FieldPrepareResult::ok)
			return false;
		if (geometry::offset_axis::sample_alpha_level_set(p.field_, inputs,
														   geometry::offset_axis::alpha_sampling_params(p.options_),
														   *p.samples_mesh_, p.samples_position_.get(),
														   p.samples_normal_.get()) != geometry::offset_axis::SamplingStatus::ok)
			return false;
		build_samples_kdtrees(p);
		if (!p.options_.apply_filtering ||
			geometry::offset_axis::filter_samples<Field>(p.input_, Scalar(p.options_.alpha), *p.samples_mesh_,
														 p.samples_position_.get(), p.samples_normal_.get()) == 0)
			return true;
		if (nb_cells<PVertex>(*p.samples_mesh_) == 0)
		{
			reason = "every sample is farther than alpha + filter_margin from the input (the field does not match the "
					 "input)";
			return false;
		}
		build_samples_kdtrees(p);
		return true;
	}

	void compute_fitting_data(PointsParameters& p)
	{
		if (p.fitting_data_computed_)
			return;
		if (nb_cells<PVertex>(*p.samples_mesh_) == 0)
		{
			p.last_error_ = "Failed to compute fitting data: there is no sample.";
			log_error(p.last_error_, '\n');
			return;
		}
		const cgogn::geometry::offset_axis::ReconstructionOptions& options = p.options_;
		// PCA for every input: neural samples already carry the network gradient, and for MF models that gradient
		// is not orthogonal to the medial axis.
		if (options.recompute_normals_after_sampling)
			geometry::offset_axis::recompute_normals_pca(*p.samples_mesh_, *p.samples_position_, *p.samples_normal_,
														 *p.samples_kdtree_, p.samples_kdtree_vertices_, options.knn_k);
		// the medial flip prune changes normals and removes samples, so it runs before the primitives that use them
		geometry::offset_axis::compute_initial_medial_axis(
			p.field_, *p.samples_mesh_, *p.samples_position_, *p.samples_normal_, p.samples_kdtree_,
			p.samples_kdtree_vertices_, *p.samples_ma_position_, *p.samples_ma_radius_, *p.samples_ma_secondary_,
			Scalar(options.alpha), options.ma_flip_prune, Scalar(options.ma_flip_prune_alpha_factor), options.knn_k);
		build_samples_kdtrees(p);
		if (nb_cells<PVertex>(*p.samples_mesh_) == 0)
		{
			p.last_error_ = "Failed to compute fitting data: the medial flip prune removed every sample.";
			log_error(p.last_error_, '\n');
			return;
		}
		geometry::offset_axis::compute_samples_area(*p.samples_mesh_, *p.samples_position_, *p.samples_normal_,
													*p.samples_area_, *p.samples_knn_, *p.samples_kdtree_,
													p.samples_kdtree_vertices_, options.knn_k, Scalar(options.alpha));
		geometry::offset_axis::compute_quadrics(*p.samples_mesh_, *p.samples_position_, *p.samples_normal_,
												*p.samples_area_, *p.samples_knn_, *p.samples_quadric_,
												*p.samples_line_quadric_, options.knn_k);
		p.last_error_.clear();
		p.fitting_data_computed_ = true;
	}

	static bool initialize_spheres(PointsParameters& p)
	{
		if (!p.fitting_data_computed_)
			return false;
		p.optimizer_->initialize_from_samples(Scalar(p.options_.init_dilation_constant),
											  p.options_.init_owner_euclidean,
											  Scalar(p.options_.sqem_update_lambda_line_plane));
		return true;
	}

	// The caller passes lambda so that the update thread can feed a live value.
	static SpheresOptimizerStatus update_spheres_once(PointsParameters& p, Scalar sqem_update_lambda_line_plane)
	{
		if (nb_cells<PVertex>(*p.spheres_) == 0)
			return SpheresOptimizerStatus::failed;
		return p.optimizer_->update_once(sqem_update_lambda_line_plane);
	}

	static bool build_skeleton(PointsParameters& p)
	{
		if (nb_cells<PVertex>(*p.spheres_) == 0)
			return false;
		p.optimizer_->compute_sphere_neighbors();
		p.topology_->build_from_spheres();
		return true;
	}

	static auto topology_scores_evaluator(PointsParameters& p)
	{
		return [&p](const std::vector<Vec3>& points, std::vector<Scalar>& values) {
			return geometry::offset_axis::evaluate_topology_scores(p.field_, *p.samples_mesh_, *p.samples_ma_position_,
																   *p.samples_ma_kdtree_, p.samples_ma_kdtree_vertices_,
																   points, values);
		};
	}

	// Residual sheet prune after the topology fix, if enabled.
	static bool fix_topology_and_prune(PointsParameters& p)
	{
		if (nb_cells<NMVertex>(*p.skeleton_) == 0 || !p.topology_->fix_topology(topology_scores_evaluator(p)))
			return false;
		if (p.options_.residual_prune)
			p.topology_->prune_residual_sheets();
		return true;
	}

	// The scores indexed by index_of(cell); the indices without a cell get NaN.
	template <typename CELL>
	static std::vector<Scalar> dense_scores(PointsParameters& p, const NMAttribute<Scalar>& scores)
	{
		uint32 size = 0;
		foreach_cell(*p.skeleton_, [&](CELL c) -> bool {
			size = std::max(size, index_of(*p.skeleton_, c) + 1);
			return true;
		});
		std::vector<Scalar> values(size, std::numeric_limits<Scalar>::quiet_NaN());
		foreach_cell(*p.skeleton_, [&](CELL c) -> bool {
			const uint32 id = index_of(*p.skeleton_, c);
			values[id] = scores[id];
			return true;
		});
		return values;
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
		for (PVertex sphere : p.topology_->added_spheres())
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
		p.topology_->compute_edge_degree();
		// Edge colors: degree 1 blue, degree 2 green, degree > 2 red, all others dim gray.
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

	// "Boundary tet mask": faces of boundary tets yellow, other faces dim gray. The longest degree-2 edge of each
	// boundary tet is orange, other edges dim gray.
	void mark_boundary_tets_color(PointsParameters& p)
	{
		const auto mask = p.topology_->boundary_tet_mask();
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

	// Palette helpers (golden-ratio hue spacing).
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

	// "Face components (UF)": color each skeleton face by its sheet. Faces are only recolored, not removed.
	void mark_skeleton_face_components_color(PointsParameters& p)
	{
		const auto sheet_of_face = cgogn::geometry::offset_axis::compute_skeleton_face_components(*p.skeleton_);
		foreach_cell(*p.skeleton_, [&](NMFace f) -> bool {
			const auto it = sheet_of_face.find(index_of(*p.skeleton_, f));
			value<Vec3>(*p.skeleton_, p.skeleton_face_color_, f) =
				(it == sheet_of_face.end()) ? Vec3(0.0, 0.0, 0.0) : component_palette_color(it->second);
			return true;
		});
	}

	// Colormaps: cells with a score are colored by color_map over the range of those scores (0.5 if the range is
	// empty); NaN (cell not adjacent to a tet) gets the dim fallback color.
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

	// "Tet face UDF colormap (Integral / Normalized)": read-only, only the face colors change.
	bool mark_tet_face_udf_color(PointsParameters& p, bool normalize_by_area)
	{
		if (!skeleton_colors_ready(p))
		{
			p.last_error_ = "No skeleton: build the skeleton first.";
			log_error(p.last_error_, '\n');
			return false;
		}
		if (!p.topology_->compute_face_scores(topology_scores_evaluator(p), normalize_by_area))
		{
			p.last_error_ = "Failed to evaluate the face UDF scores (is the distance field loaded?).";
			log_error(p.last_error_, '\n');
			return false;
		}
		color_cells_by_scores<NMFace>(p, dense_scores<NMFace>(p, *p.skeleton_face_score_), p.skeleton_face_color_,
									  Vec3(0.2, 0.2, 0.2));
		p.last_error_.clear();
		return true;
	}

	// "Edge UDF score colormap (Gauss3)": read-only, only the edge colors change.
	bool mark_tet_edge_gauss3_color(PointsParameters& p)
	{
		if (!skeleton_colors_ready(p))
		{
			p.last_error_ = "No skeleton: build the skeleton first.";
			log_error(p.last_error_, '\n');
			return false;
		}
		if (!p.topology_->compute_edge_scores(topology_scores_evaluator(p)))
		{
			p.last_error_ = "Failed to evaluate the edge UDF scores (is the distance field loaded?).";
			log_error(p.last_error_, '\n');
			return false;
		}
		color_cells_by_scores<NMEdge>(p, dense_scores<NMEdge>(p, *p.skeleton_edge_score_), p.skeleton_edge_color_,
									  Vec3(0.08, 0.08, 0.08));
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
			const bool built = build_skeleton(p);
			p.skeleton_invalidated_ = !built;
			if (built)
			{
				p.last_error_.clear();
				refresh_skeleton_topology_colors(p);
			}
			else
				p.last_error_ = "Failed to build skeleton: there is no sphere.";
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
			p.optimizer_->reset_optimization();
			p.metrics_snapshot_ = p.optimizer_->metrics();
			p.sphere_update_lambda_.store(p.options_.sqem_update_lambda_line_plane);
		}
		p.running_.store(true);

		spheres_update_thread_ = std::thread([this, active = &p]() {
			cgogn::thread_start(cgogn::max_nb_threads() - 1);
			auto start = std::chrono::high_resolution_clock::now();
			auto last_preview_refresh = start;
			SpheresOptimizerType& optimizer = *active->optimizer_;
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
						status = update_spheres_once(*active, lambda);
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
			active->metrics_snapshot_ = active->optimizer_->metrics();
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
		auto& options = p.options_;
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
			if (ImGui::InputFloat("Alpha", &options.alpha, 0.001f, 0.1f, "%.4f"))
				options.alpha = std::clamp(options.alpha, float(cgogn::geometry::offset_axis::option_info::alpha.min),
										   float(cgogn::geometry::offset_axis::option_info::alpha.max));
			if (ImGui::InputFloat("Sample Radius", &options.sample_radius, 0.001f, 0.01f, "%.4f"))
				options.sample_radius =
					std::clamp(options.sample_radius, float(cgogn::geometry::offset_axis::option_info::sample_radius.min),
							   float(cgogn::geometry::offset_axis::option_info::sample_radius.max));
			if (ImGui::InputInt("KNN K", &options.knn_k, 1, 5))
				options.knn_k = std::clamp(options.knn_k, int(cgogn::geometry::offset_axis::option_info::knn_k.min),
										   int(cgogn::geometry::offset_axis::option_info::knn_k.max));
			ImGui::Checkbox("Filter samples by distance to input", &options.apply_filtering);
			ImGui::Checkbox("Recompute normals after sampling", &options.recompute_normals_after_sampling);
			if (ImGui::InputInt("Bridson candidates per parent", &options.alpha_sampling_candidates_per_parent, 1, 10))
				options.alpha_sampling_candidates_per_parent =
					std::clamp(options.alpha_sampling_candidates_per_parent,
							   int(cgogn::geometry::offset_axis::option_info::alpha_sampling_candidates_per_parent.min),
							   int(cgogn::geometry::offset_axis::option_info::alpha_sampling_candidates_per_parent.max));
			if (ImGui::InputInt("Bridson parent batch size", &options.alpha_sampling_parent_batch_size, 1, 64))
				options.alpha_sampling_parent_batch_size =
					std::clamp(options.alpha_sampling_parent_batch_size,
							   int(cgogn::geometry::offset_axis::option_info::alpha_sampling_parent_batch_size.min),
							   int(cgogn::geometry::offset_axis::option_info::alpha_sampling_parent_batch_size.max));
			if (ImGui::InputInt("Reprojection iterations (MF, 0 = off)", &options.reprojection_iterations, 1, 5))
				options.reprojection_iterations =
					std::clamp(options.reprojection_iterations,
							   int(cgogn::geometry::offset_axis::option_info::reprojection_iterations.min),
							   int(cgogn::geometry::offset_axis::option_info::reprojection_iterations.max));
			if (ImGui::InputInt("Seed", &options.seed, 1, 10))
				options.seed = std::clamp(options.seed, int(cgogn::geometry::offset_axis::option_info::seed.min),
										  int(cgogn::geometry::offset_axis::option_info::seed.max));
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
								   float(cgogn::geometry::offset_axis::option_info::ma_flip_prune_alpha_factor.min),
								   float(cgogn::geometry::offset_axis::option_info::ma_flip_prune_alpha_factor.max), "%.2f"))
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
							   float(cgogn::geometry::offset_axis::option_info::init_dilation_constant.min),
							   float(cgogn::geometry::offset_axis::option_info::init_dilation_constant.max), "%.4f");
			ImGui::Checkbox("Euclidean initial owners", &options.init_owner_euclidean);
			ImGui::EndDisabled();
			if (!running && ImGui::Button("Init spheres"))
			{
				std::lock_guard<std::mutex> lock(p.mutex_);
				p.spheres_initialized_ = initialize_spheres(p);
				if (!p.spheres_initialized_)
				{
					p.last_error_ = "Failed to initialize spheres: compute the fitting data first.";
					log_error(p.last_error_, '\n');
				}
				else
				{
					p.last_error_.clear();
					p.skeleton_invalidated_ = true;
					if (p.skeleton_)
						clear(*p.skeleton_);
					set_post_init_sphere_render_state(p);
					p.metrics_snapshot_ = p.optimizer_->metrics();
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
								   float(cgogn::geometry::offset_axis::option_info::sqem_update_lambda_line_plane.min),
								   float(cgogn::geometry::offset_axis::option_info::sqem_update_lambda_line_plane.max), "%.6f"))
			{
				p.sphere_update_lambda_.store(update_lambda);
				options.sqem_update_lambda_line_plane = update_lambda;
			}
			if (!running && ImGui::Button("Update spheres"))
			{
				std::lock_guard<std::mutex> lock(p.mutex_);
				SpheresOptimizerType& optimizer = *p.optimizer_;
				update_spheres_once(p, Scalar(options.sqem_update_lambda_line_plane));
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
				const bool built = build_skeleton(p);
				p.skeleton_invalidated_ = !built;
				if (built)
				{
					p.last_error_.clear();
					refresh_skeleton_topology_colors(p);
				}
				else
					p.last_error_ = "Failed to build skeleton: there is no sphere.";
				update_render_data(p, RenderUpdateMode::all_current_data);
			}
		}

		// Topology Fix. While the sphere update runs, the topology is not read (the update thread writes it), so the
		// step stays visible with its parameters disabled.
		const bool has_tets = !running && spheres_ready && skeleton_built && p.topology_->has_tets();
		const bool topology_fix_step = spheres_ready && skeleton_built && (running || has_tets);
		if (topology_fix_step && ImGui::CollapsingHeader("Topology Fix", ImGuiTreeNodeFlags_DefaultOpen))
		{
			ImGui::BeginDisabled(running);
			ImGui::Checkbox("Residual prune", &options.residual_prune);
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
					if (fix_topology_and_prune(p))
					{
						p.last_error_.clear();
						p.metrics_snapshot_ = p.optimizer_->metrics();
						apply_skeleton_topology_results(p);
					}
					else
					{
						p.last_error_ = "[TopologyFull] score evaluation failed or topology data is invalid";
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
					if (!cgogn::io::export_skeleton_ply(*p.skeleton_, p.skeleton_position_, p.skeleton_radius_, filename))
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

#endif // CGOGN_MODULE_OFFSET_AXIS_H_
