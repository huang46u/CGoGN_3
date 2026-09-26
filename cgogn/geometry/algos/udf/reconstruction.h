/*******************************************************************************
 * CGoGN: Combinatorial and Geometric modeling with Generic N-dimensional Maps  *
 * Copyright (C), IGG Group, ICube, University of Strasbourg, France            *
 *******************************************************************************/
#ifndef CGOGN_GEOMETRY_ALGOS_UDF_RECONSTRUCTION_H_
#define CGOGN_GEOMETRY_ALGOS_UDF_RECONSTRUCTION_H_

#include <cgogn/core/functions/attributes.h>
#include <cgogn/core/functions/mesh_info.h>
#include <cgogn/core/functions/traversals/edge.h>
#include <cgogn/core/types/cell_marker.h>
#include <cgogn/geometry/algos/area.h>
#include <cgogn/geometry/algos/fitting.h>
#include <cgogn/geometry/algos/medial_axis.h>
#include <cgogn/geometry/algos/normal.h>
#include <cgogn/geometry/algos/udf/alpha_projection.h>
#include <cgogn/geometry/algos/udf/alpha_sampling.h>
#include <cgogn/geometry/algos/udf/neural_alpha_projection.h>
#include <cgogn/geometry/algos/udf/neural_field_query.h>
#include <cgogn/geometry/algos/udf/sample_processing.h>
#include <cgogn/geometry/algos/udf/skeleton_topology.h>
#include <cgogn/geometry/algos/udf/spatial_query.h>
#include <cgogn/geometry/algos/udf/spheres_optimizer.h>
#include <cgogn/geometry/functions/bounding_box.h>
#include <cgogn/geometry/functions/normal.h>
#include <cgogn/geometry/types/line_quadric.h>
#include <cgogn/geometry/types/neural_field_forward.h>
#include <cgogn/geometry/types/ray_level_set_sampler.h>
#include <cgogn/geometry/types/ray_level_set_sampler_traits.h>
#include <cgogn/geometry/types/spatial_grid.h>
#include <cgogn/geometry/types/spherical_quadric.h>
#include <cgogn/geometry/types/vector_traits.h>
#include <cgogn/io/surface/export_options.h>
#include <cgogn/io/surface/ply.h>

#include <libacc/bvh_tree.h>
#include <libacc/kd_tree.h>

#include <Eigen/Dense>

#include <algorithm>
#include <chrono>
#include <cmath>
#include <filesystem>
#include <functional>
#include <memory>
#include <set>
#include <stdexcept>
#include <string>
#include <tuple>
#include <type_traits>
#include <utility>
#include <vector>

namespace cgogn::geometry
{

template <typename SURFACE, typename POINTS, typename NONMANIFOLD, typename RaySamplerTag>
class UDFReconstruction
{
public:
	using PVertex = typename mesh_traits<POINTS>::Vertex;
	using SVertex = typename mesh_traits<SURFACE>::Vertex;
	using SFace = typename mesh_traits<SURFACE>::Face;
	using NMVertex = typename mesh_traits<NONMANIFOLD>::Vertex;
	using NMEdge = typename mesh_traits<NONMANIFOLD>::Edge;
	using NMFace = typename mesh_traits<NONMANIFOLD>::Face;
	using SpheresOptimizerType = SpheresOptimizer<POINTS>;
	using SpheresOptimizerMetrics = typename SpheresOptimizerType::Metrics;
	using TopologyType = SkeletonTopology<POINTS, NONMANIFOLD>;
	using RayConfig = RaySamplerTraits<RaySamplerTag, SURFACE, POINTS>;
	using RayTraits = typename RayConfig::SamplerTraits;
	using RaySampler = RayLevelSetSampler<RayTraits>;
	using RayParams = typename RaySampler::Params;
	using PAttribute = typename mesh_traits<POINTS>::template Attribute<Vec3>;
	template <typename T>
	using PointAttribute = typename mesh_traits<POINTS>::template Attribute<T>;
	template <typename T>
	using SurfaceAttribute = typename mesh_traits<SURFACE>::template Attribute<T>;
	template <typename T>
	using SkeletonAttribute = typename mesh_traits<NONMANIFOLD>::template Attribute<T>;

	enum class NeuralModelType
	{
		udf,
		mf
	};
	enum class Status
	{
		success,
		invalid_data,
		model_load_failed,
		sampling_failed,
		fitting_failed,
		medial_axis_failed,
		sphere_initialization_failed,
		optimization_failed,
		skeleton_build_failed,
		topology_fix_failed
	};

	struct Data
	{
		SURFACE* surface = nullptr;
		POINTS* input_points = nullptr;
		POINTS* samples = nullptr;
		POINTS* spheres = nullptr;
		NONMANIFOLD* skeleton = nullptr;
	};

	struct Options
	{
		NeuralModelType neural_model_type = NeuralModelType::udf;
		bool ma_flip_prune = true;
		float ma_flip_prune_alpha_factor = 1.0f;
		float alpha = 0.005f;
		int knn_k = 10;
		bool apply_filtering = false;
		bool recompute_normals_after_sampling = false;
		int ray_sampler_batch_size = 4096;
		int batch_size = 1310640;
		float tolerance = 1e-6f;
		int udf_max_iterations = 3000;
		float sample_radius = 0.0025f;
		float sqem_update_lambda_line_plane = 0.2f;
		int seed = 42;
		float init_dilation_constant = 0.001f;
		bool residual_prune = false;
	};

	struct Timing
	{
		float64 preprocessing_ms = 0.0;
		float64 sampling_ms = 0.0;
		float64 sample_filtering_ms = 0.0;
		float64 kdtree_and_normals_ms = 0.0;
		float64 fitting_primitives_ms = 0.0;
		float64 initial_medial_axis_ms = 0.0;
		float64 sphere_initialization_ms = 0.0;
		float64 optimization_ms = 0.0;
		float64 cluster_ms = 0.0;
		float64 sphere_update_ms = 0.0;
		float64 error_ms = 0.0;
		float64 skeleton_construction_ms = 0.0;
		float64 topology_processing_ms = 0.0;
	};

	struct Counts
	{
		uint32 input_vertices = 0;
		uint32 input_points = 0;
		uint32 samples_before_filtering = 0;
		uint32 samples = 0;
		uint32 spheres = 0;
		uint32 skeleton_vertices = 0;
		uint32 skeleton_edges = 0;
		uint32 skeleton_faces = 0;
		uint32 optimization_iterations = 0;
	};

	struct Result
	{
		Status status = Status::invalid_data;
		Timing timing;
		Counts counts;
	};

	using IterationCallback = std::function<void(const SpheresOptimizerMetrics&)>;

	explicit UDFReconstruction(Data data, Options options = {}) : data_(data), options_(options)
	{
	}
	UDFReconstruction(const UDFReconstruction&) = delete;
	UDFReconstruction& operator=(const UDFReconstruction&) = delete;

	Data& data()
	{
		return data_;
	}
	const Data& data() const
	{
		return data_;
	}
	Options& options()
	{
		return options_;
	}
	const Options& options() const
	{
		return options_;
	}

	void reset()
	{
		input_kdtree_.reset();
		sample_kdtree_.reset();
		sample_ma_kdtree_.reset();
		sample_vertices_.clear();
		sample_ma_vertices_.clear();
		input_vertices_.clear();
		surface_bvh_.reset();
		surface_faces_.clear();
		surface_vertices_.clear();
		surface_positions_.clear();
		ray_sampler_.reset();
		optimizer_.reset();
		topology_.reset();
		neural_model_loaded_ = false;
		fitting_data_computed_ = false;
		face_components_prepared_ = false;
		samples_before_filtering_ = 0;
		initialized_ = false;
	}

	Status initialize_data()
	{
		if (initialized_)
			return Status::success;
		if (!data_.input_points || !data_.samples || !data_.spheres || !data_.skeleton)
			return Status::invalid_data;
		input_position_ = get_or_add_attribute<Vec3, PVertex>(*data_.input_points, "position");
		input_normal_ = get_or_add_attribute<Vec3, PVertex>(*data_.input_points, "normal");
		input_knn_ = get_or_add_attribute<std::vector<PVertex>, PVertex>(*data_.input_points, "knn");
		sample_position_ = get_or_add_attribute<Vec3, PVertex>(*data_.samples, "position");
		sample_normal_ = get_or_add_attribute<Vec3, PVertex>(*data_.samples, "normal");
		sample_area_ = get_or_add_attribute<Scalar, PVertex>(*data_.samples, "area");
		sample_knn_ = get_or_add_attribute<std::vector<PVertex>, PVertex>(*data_.samples, "knn");
		sample_quadric_ = get_or_add_attribute<Spherical_Quadric, PVertex>(*data_.samples, "quadric");
		sample_line_quadric_ = get_or_add_attribute<Line_Quadric, PVertex>(*data_.samples, "line_quadric");
		sample_ma_position_ = get_or_add_attribute<Vec3, PVertex>(*data_.samples, "ma_position");
		sample_ma_radius_ = get_or_add_attribute<Scalar, PVertex>(*data_.samples, "ma_radius");
		sample_ma_secondary_ = get_or_add_attribute<PVertex, PVertex>(*data_.samples, "ma_secondary_vertex");
		sample_sphere_ = get_or_add_attribute<PVertex, PVertex>(*data_.samples, "sphere");
		sample_error_ = get_or_add_attribute<Scalar, PVertex>(*data_.samples, "error");
		sphere_position_ = get_or_add_attribute<Vec3, PVertex>(*data_.spheres, "position");
		sphere_radius_ = get_or_add_attribute<Scalar, PVertex>(*data_.spheres, "radius");
		sphere_cluster_ = get_or_add_attribute<std::vector<PVertex>, PVertex>(*data_.spheres, "cluster");
		sphere_cluster_area_ = get_or_add_attribute<Scalar, PVertex>(*data_.spheres, "cluster_area");
		sphere_cluster_color_ = get_or_add_attribute<Vec4, PVertex>(*data_.spheres, "cluster_color");
		sphere_neighbors_ = get_or_add_attribute<std::set<PVertex>, PVertex>(*data_.spheres, "neighbor_clusters");
		sphere_error_ = get_or_add_attribute<Scalar, PVertex>(*data_.spheres, "error");
		sphere_error_raw_ = get_or_add_attribute<Scalar, PVertex>(*data_.spheres, "error_not_normalized");
		sphere_skeleton_vertex_ = get_or_add_attribute<NMVertex, PVertex>(*data_.spheres, "skeleton_vertex");
		skeleton_position_ = get_or_add_attribute<Vec3, NMVertex>(*data_.skeleton, "position");
		skeleton_radius_ = get_or_add_attribute<Scalar, NMVertex>(*data_.skeleton, "radius");
		skeleton_source_sphere_ = get_or_add_attribute<PVertex, NMVertex>(*data_.skeleton, "source_sphere");
		skeleton_incident_tets_ = get_or_add_attribute<std::set<std::size_t>, NMFace>(*data_.skeleton, "incident_tets");
		skeleton_component_id_ = get_or_add_attribute<uint32, NMFace>(*data_.skeleton, "face_component_id");
		skeleton_edge_degree_ = get_or_add_attribute<uint32, NMEdge>(*data_.skeleton, "degree");
		skeleton_face_component_color_ = get_or_add_attribute<Vec3, NMFace>(*data_.skeleton, "face_component_color");
		initialized_ = input_position_ && sample_position_ && sample_normal_ && sphere_position_ && skeleton_position_;
		return initialized_ ? Status::success : Status::invalid_data;
	}

	Status load_neural_model(const std::filesystem::path& path, NeuralModelType type)
	{
		if constexpr (!std::is_same_v<RaySamplerTag, RaySamplerNeural>)
		{
			(void)path;
			(void)type;
			return Status::model_load_failed;
		}
		else
		{
			try
			{
				neural_model_ = torch::jit::load(path.string(), device_);
				neural_model_.eval();
				neural_model_loaded_ = true;
				options_.neural_model_type = type;
				return Status::success;
			}
			catch (const c10::Error&)
			{
				neural_model_loaded_ = false;
				return Status::model_load_failed;
			}
		}
	}

	Status sample_alpha_level_set();
	Status apply_sampling_filtering();
	Status build_kdtree_and_normals();
	Status compute_fitting_primitives();
	Status compute_initial_medial_axis();
	Status initialize_spheres();
	Status optimize_spheres(const IterationCallback& iteration_callback = {});
	Status build_skeleton();
	Status fix_topology();
	Status prune_residual_sheets();
	Status prepare_face_components();
	Result run(const IterationCallback& iteration_callback = {});

	Counts counts() const
	{
		Counts c;
		if (data_.input_points)
			c.input_points = nb_cells<PVertex>(*data_.input_points);
		c.input_vertices = c.input_points;
		if (data_.surface)
			c.input_vertices = nb_cells<SVertex>(*data_.surface);
		if (data_.samples)
			c.samples = nb_cells<PVertex>(*data_.samples);
		c.samples_before_filtering = samples_before_filtering_;
		if (data_.spheres)
			c.spheres = nb_cells<PVertex>(*data_.spheres);
		if (data_.skeleton)
		{
			c.skeleton_vertices = nb_cells<NMVertex>(*data_.skeleton);
			c.skeleton_edges = nb_cells<NMEdge>(*data_.skeleton);
			c.skeleton_faces = nb_cells<NMFace>(*data_.skeleton);
		}
		if (optimizer_)
			c.optimization_iterations = optimizer_->metrics().iteration;
		return c;
	}

	bool export_skeleton_ply(const std::filesystem::path& path, bool include_face_components) const
	{
		if (!data_.skeleton || !skeleton_position_ || !skeleton_radius_)
			return false;
		if (include_face_components && !face_components_prepared_)
			return false;
		std::error_code error;
		if (!path.parent_path().empty())
		{
			std::filesystem::create_directories(path.parent_path(), error);
			if (error)
				return false;
		}
		if (std::filesystem::exists(path, error))
		{
			std::filesystem::remove(path, error);
			if (error)
				return false;
		}
		io::SurfaceExportAttributeSelection<NONMANIFOLD> attributes;
		attributes.vertex_attributes.push_back(skeleton_radius_);
		if (include_face_components && skeleton_face_component_color_)
			attributes.face_attributes.push_back(skeleton_face_component_color_);
		io::export_PLY(*data_.skeleton, skeleton_position_.get(), path.string(), &attributes);
		const bool exists = std::filesystem::exists(path, error);
		if (error || !exists)
			return false;
		const auto size = std::filesystem::file_size(path, error);
		return !error && size > 0;
	}

	SpheresOptimizerType& spheres_optimizer()
	{
		return optimizer_for();
	}
	TopologyType& skeleton_topology()
	{
		return topology_for();
	}

private:
	Data data_;
	Options options_;
	torch::Device device_ = torch::cuda::is_available() ? torch::Device(torch::kCUDA, 0) : torch::Device(torch::kCPU);
	torch::jit::Module neural_model_;
	bool neural_model_loaded_ = false;
	NeuralProjectionWorkspace projection_workspace_;
	std::unique_ptr<RaySampler> ray_sampler_;
	std::unique_ptr<acc::KDTree<3, uint32>> input_kdtree_;
	std::unique_ptr<acc::KDTree<3, uint32>> sample_kdtree_;
	std::unique_ptr<acc::KDTree<3, uint32>> sample_ma_kdtree_;
	std::vector<PVertex> input_vertices_;
	std::vector<PVertex> sample_vertices_;
	std::vector<PVertex> sample_ma_vertices_;
	std::unique_ptr<acc::BVHTree<uint32, Vec3>> surface_bvh_;
	std::vector<SFace> surface_faces_;
	std::vector<SVertex> surface_vertices_;
	std::vector<Vec3> surface_positions_;
	std::unique_ptr<SpheresOptimizerType> optimizer_;
	std::unique_ptr<TopologyType> topology_;
	bool fitting_data_computed_ = false;
	bool initialized_ = false;
	bool face_components_prepared_ = false;
	uint32 sphere_count_ = 0;
	uint32 samples_before_filtering_ = 0;

	std::shared_ptr<PointAttribute<Vec3>> input_position_, input_normal_;
	std::shared_ptr<PointAttribute<std::vector<PVertex>>> input_knn_;
	std::shared_ptr<PointAttribute<Vec3>> sample_position_, sample_normal_;
	std::shared_ptr<PointAttribute<Scalar>> sample_area_;
	std::shared_ptr<PointAttribute<std::vector<PVertex>>> sample_knn_;
	std::shared_ptr<PointAttribute<Spherical_Quadric>> sample_quadric_;
	std::shared_ptr<PointAttribute<Line_Quadric>> sample_line_quadric_;
	std::shared_ptr<PointAttribute<Vec3>> sample_ma_position_;
	std::shared_ptr<PointAttribute<Scalar>> sample_ma_radius_;
	std::shared_ptr<PointAttribute<PVertex>> sample_ma_secondary_, sample_sphere_;
	std::shared_ptr<PointAttribute<Scalar>> sample_error_;
	std::shared_ptr<PointAttribute<Vec3>> sphere_position_;
	std::shared_ptr<PointAttribute<Scalar>> sphere_radius_;
	std::shared_ptr<PointAttribute<std::vector<PVertex>>> sphere_cluster_;
	std::shared_ptr<PointAttribute<Scalar>> sphere_cluster_area_;
	std::shared_ptr<PointAttribute<Vec4>> sphere_cluster_color_;
	std::shared_ptr<PointAttribute<std::set<PVertex>>> sphere_neighbors_;
	std::shared_ptr<PointAttribute<Scalar>> sphere_error_, sphere_error_raw_;
	std::shared_ptr<PointAttribute<NMVertex>> sphere_skeleton_vertex_;
	std::shared_ptr<SkeletonAttribute<Vec3>> skeleton_position_;
	std::shared_ptr<SkeletonAttribute<Scalar>> skeleton_radius_;
	std::shared_ptr<SkeletonAttribute<PVertex>> skeleton_source_sphere_;
	std::shared_ptr<SkeletonAttribute<std::set<std::size_t>>> skeleton_incident_tets_;
	std::shared_ptr<SkeletonAttribute<uint32>> skeleton_component_id_, skeleton_edge_degree_;
	std::shared_ptr<SkeletonAttribute<Vec3>> skeleton_face_component_color_;
	std::shared_ptr<SurfaceAttribute<Vec3>> surface_vertex_normal_;

	void build_surface_bvh();
	void build_input_kdtree();
	void build_sample_kdtree();
	bool query_surface_projection_info(const Vec3& query, AlphaProjectionInfo& info);
	bool query_point_cloud_projection_info(const Vec3& query, AlphaProjectionInfo& info);
	void recompute_samples_normals_from_input();
	bool recompute_sample_normal_pca(PVertex vertex);
	NeuralFieldForward neural_field()
	{
		return NeuralFieldForward(&neural_model_, neural_model_loaded_, device_);
	}
	SpheresOptimizerType& optimizer_for();
	TopologyType& topology_for();
	bool evaluate_topology_scores(const std::vector<Vec3>& points, std::vector<Scalar>& values);
	Vec3 hsv_to_rgb(Scalar h, Scalar s, Scalar v) const;
	Vec3 component_palette_color(uint32 component_id) const;
};

template <typename SURFACE, typename POINTS, typename NONMANIFOLD, typename RaySamplerTag>
void UDFReconstruction<SURFACE, POINTS, NONMANIFOLD, RaySamplerTag>::build_input_kdtree()
{
	input_kdtree_.reset();
	input_vertices_.clear();
	if (!data_.input_points || !input_position_)
		return;
	std::vector<Vec3> positions;
	foreach_cell(*data_.input_points, [&](PVertex v) {
		const uint32 idx = index_of(*data_.input_points, v);
		positions.push_back((*input_position_)[idx]);
		input_vertices_.push_back(v);
		return true;
	});
	if (!positions.empty())
		input_kdtree_ = std::make_unique<acc::KDTree<3, uint32>>(positions);
}

template <typename SURFACE, typename POINTS, typename NONMANIFOLD, typename RaySamplerTag>
void UDFReconstruction<SURFACE, POINTS, NONMANIFOLD, RaySamplerTag>::build_sample_kdtree()
{
	sample_kdtree_.reset();
	sample_ma_kdtree_.reset();
	sample_vertices_.clear();
	sample_ma_vertices_.clear();
	if (!data_.samples || !sample_position_)
		return;
	std::vector<Vec3> positions;
	std::vector<Vec3> ma_positions;
	foreach_cell(*data_.samples, [&](PVertex v) {
		const uint32 idx = index_of(*data_.samples, v);
		positions.push_back((*sample_position_)[idx]);
		sample_vertices_.push_back(v);
		const Vec3 ma = (*sample_ma_position_)[idx];
		if (ma.allFinite() && (*sample_ma_radius_)[idx] > Scalar(0))
		{
			ma_positions.push_back(ma);
			sample_ma_vertices_.push_back(v);
		}
		return true;
	});
	if (!positions.empty())
		sample_kdtree_ = std::make_unique<acc::KDTree<3, uint32>>(positions);
	if (!ma_positions.empty())
		sample_ma_kdtree_ = std::make_unique<acc::KDTree<3, uint32>>(ma_positions);
}

template <typename SURFACE, typename POINTS, typename NONMANIFOLD, typename RaySamplerTag>
void UDFReconstruction<SURFACE, POINTS, NONMANIFOLD, RaySamplerTag>::build_surface_bvh()
{
	surface_bvh_.reset();
	surface_faces_.clear();
	surface_vertices_.clear();
	surface_positions_.clear();
	if (!data_.surface)
		return;
	auto position = get_attribute<Vec3, SVertex>(*data_.surface, "position");
	if (!position)
		return;
	surface_vertex_normal_ = get_or_add_attribute<Vec3, SVertex>(*data_.surface, "normal");
	geometry::compute_normal<SVertex>(*data_.surface, position.get(), surface_vertex_normal_.get());
	auto vertex_index = get_or_add_attribute<uint32, SVertex>(*data_.surface, "__udf_bvh_vertex_index");
	uint32 next = 0;
	foreach_cell(*data_.surface, [&](SVertex v) {
		surface_vertices_.push_back(v);
		surface_positions_.push_back((*position)[index_of(*data_.surface, v)]);
		value<uint32>(*data_.surface, vertex_index, v) = next++;
		return true;
	});
	std::vector<uint32> face_indices;
	foreach_cell(*data_.surface, [&](typename mesh_traits<SURFACE>::Face f) {
		surface_faces_.push_back(f);
		foreach_incident_vertex(*data_.surface, f, [&](SVertex v) {
			face_indices.push_back(value<uint32>(*data_.surface, vertex_index, v));
			return true;
		});
		return true;
	});
	remove_attribute<SVertex>(*data_.surface, vertex_index);
	if (!face_indices.empty() && !surface_positions_.empty())
		surface_bvh_ = std::make_unique<acc::BVHTree<uint32, Vec3>>(face_indices, surface_positions_);
}

template <typename SURFACE, typename POINTS, typename NONMANIFOLD, typename RaySamplerTag>
bool UDFReconstruction<SURFACE, POINTS, NONMANIFOLD, RaySamplerTag>::query_surface_projection_info(
	const Vec3& query, AlphaProjectionInfo& info)
{
	info = AlphaProjectionInfo{};
	if (!data_.surface || !surface_bvh_ || !surface_vertex_normal_)
		return false;
	uint32 primitive;
	if (!query_surface_closest_point(*surface_bvh_, query, primitive, info.closest_point, info.distance))
		return false;
	info.normal = query - info.closest_point;
	if (info.normal.squaredNorm() < Scalar(1e-12) && primitive < surface_faces_.size())
		info.normal = geometry::normal(*data_.surface, surface_faces_[primitive], surface_vertex_normal_.get());
	if (info.normal.squaredNorm() < Scalar(1e-12))
		info.normal = Vec3(0, 0, 1);
	else
		info.normal.normalize();
	info.valid =
		info.closest_point.allFinite() && info.normal.allFinite() && std::isfinite(static_cast<double>(info.distance));
	return info.valid;
}

template <typename SURFACE, typename POINTS, typename NONMANIFOLD, typename RaySamplerTag>
bool UDFReconstruction<SURFACE, POINTS, NONMANIFOLD, RaySamplerTag>::query_point_cloud_projection_info(
	const Vec3& query, AlphaProjectionInfo& info)
{
	info = AlphaProjectionInfo{};
	if (!input_kdtree_ || !data_.input_points || input_vertices_.empty())
		return false;
	uint32 point_index;
	if (!query_point_cloud_nearest_point(*input_kdtree_, query, point_index, info.distance) ||
		point_index >= input_vertices_.size())
		return false;
	const uint32 vertex_index = index_of(*data_.input_points, input_vertices_[point_index]);
	if (vertex_index == INVALID_INDEX)
		return false;
	info.closest_point = (*input_position_)[vertex_index];
	info.normal = query - info.closest_point;
	if (info.normal.squaredNorm() < Scalar(1e-12))
		info.normal = (*input_normal_)[vertex_index];
	if (info.normal.squaredNorm() < Scalar(1e-12))
		info.normal = Vec3(0, 0, 1);
	else
		info.normal.normalize();
	info.valid =
		info.closest_point.allFinite() && info.normal.allFinite() && std::isfinite(static_cast<double>(info.distance));
	return info.valid;
}

template <typename SURFACE, typename POINTS, typename NONMANIFOLD, typename RaySamplerTag>
void UDFReconstruction<SURFACE, POINTS, NONMANIFOLD, RaySamplerTag>::recompute_samples_normals_from_input()
{
	if constexpr (std::is_same_v<RaySamplerTag, RaySamplerNeural>)
	{
		std::vector<Vec3> positions, normals;
		positions.reserve(nb_cells<PVertex>(*data_.samples));
		foreach_cell(*data_.samples, [&](PVertex v) {
			positions.push_back((*sample_position_)[index_of(*data_.samples, v)]);
			return true;
		});
		if (evaluate_udf_normals(neural_field(), positions, static_cast<size_t>(std::max(1, options_.batch_size)),
								 projection_workspace_, normals) &&
			normals.size() == positions.size())
		{
			uint32 i = 0;
			foreach_cell(*data_.samples, [&](PVertex v) {
				(*sample_normal_)[index_of(*data_.samples, v)] = normals[i++];
				return true;
			});
			return;
		}
	}
	if (!input_kdtree_ && data_.surface)
	{
		build_surface_bvh();
		if (surface_bvh_)
		{
			foreach_cell(*data_.samples, [&](PVertex v) {
				const uint32 idx = index_of(*data_.samples, v);
				const Vec3& position = (*sample_position_)[idx];
				std::pair<uint32, Vec3> closest;
				Vec3 normal(0, 0, 1);
				if (surface_bvh_->closest_point(position, &closest))
				{
					normal = position - closest.second;
					if (normal.squaredNorm() < Scalar(1e-12))
						normal = Vec3(0, 0, 1);
					else
						normal.normalize();
				}
				(*sample_normal_)[idx] = normal;
				return true;
			});
			return;
		}
	}
	if (input_kdtree_ && data_.input_points)
	{
		foreach_cell(*data_.samples, [&](PVertex v) {
			const uint32 idx = index_of(*data_.samples, v);
			const Vec3& position = (*sample_position_)[idx];
			std::pair<uint32, Scalar> nearest;
			Vec3 normal(0, 0, 1);
			if (input_kdtree_->find_nn(position, &nearest) && nearest.first < input_vertices_.size())
			{
				const uint32 input_idx = index_of(*data_.input_points, input_vertices_[nearest.first]);
				if (input_idx != INVALID_INDEX)
				{
					normal = position - (*input_position_)[input_idx];
					if (normal.squaredNorm() < Scalar(1e-12))
						normal = (*input_normal_)[input_idx];
					if (normal.squaredNorm() < Scalar(1e-12))
						normal = Vec3(0, 0, 1);
					else
						normal.normalize();
				}
			}
			(*sample_normal_)[idx] = normal;
			return true;
		});
		return;
	}
	foreach_cell(*data_.samples, [&](PVertex v) {
		(*sample_normal_)[index_of(*data_.samples, v)] = Vec3(0, 0, 1);
		return true;
	});
}

template <typename SURFACE, typename POINTS, typename NONMANIFOLD, typename RaySamplerTag>
bool UDFReconstruction<SURFACE, POINTS, NONMANIFOLD, RaySamplerTag>::recompute_sample_normal_pca(PVertex vertex)
{
	if (!data_.samples || !sample_kdtree_ || !vertex.is_valid())
		return false;
	const int k = std::max(3, options_.knn_k);
	const uint32 vertex_index = index_of(*data_.samples, vertex);
	if (vertex_index == INVALID_INDEX)
		return false;
	const Vec3& center = (*sample_position_)[vertex_index];
	std::vector<std::pair<uint32, Scalar>> neighbors_index;
	sample_kdtree_->find_nns(center, k + 1, &neighbors_index);
	std::vector<Vec3> neighbors;
	neighbors.reserve(k);
	for (const auto& item : neighbors_index)
	{
		const PVertex neighbor = sample_vertices_[item.first];
		if (neighbor == vertex)
			continue;
		const uint32 neighbor_index = index_of(*data_.samples, neighbor);
		neighbors.push_back((*sample_position_)[neighbor_index]);
		if (static_cast<int>(neighbors.size()) >= k)
			break;
	}
	if (neighbors.size() < 3)
		return false;
	Vec3 mean(0, 0, 0);
	for (const Vec3& p : neighbors)
		mean += p;
	mean /= Scalar(neighbors.size());
	Eigen::Matrix<Scalar, 3, 3> covariance = Eigen::Matrix<Scalar, 3, 3>::Zero();
	for (const Vec3& p : neighbors)
	{
		const Vec3 d = p - mean;
		covariance(0, 0) += d.x() * d.x();
		covariance(0, 1) += d.x() * d.y();
		covariance(0, 2) += d.x() * d.z();
		covariance(1, 1) += d.y() * d.y();
		covariance(1, 2) += d.y() * d.z();
		covariance(2, 2) += d.z() * d.z();
	}
	covariance(1, 0) = covariance(0, 1);
	covariance(2, 0) = covariance(0, 2);
	covariance(2, 1) = covariance(1, 2);
	Eigen::SelfAdjointEigenSolver<Eigen::Matrix<Scalar, 3, 3>> solver(covariance);
	if (solver.info() != Eigen::Success)
		return false;
	const auto eigenvector = solver.eigenvectors().col(0);
	Vec3 normal(eigenvector(0), eigenvector(1), eigenvector(2));
	if (normal.squaredNorm() < Scalar(1e-12))
		return false;
	const Vec3 current = (*sample_normal_)[vertex_index];
	if (current.squaredNorm() > Scalar(1e-12) && normal.dot(current) < Scalar(0))
		normal = -normal;
	normal.normalize();
	(*sample_normal_)[vertex_index] = normal;
	return true;
}

template <typename SURFACE, typename POINTS, typename NONMANIFOLD, typename RaySamplerTag>
typename UDFReconstruction<SURFACE, POINTS, NONMANIFOLD, RaySamplerTag>::SpheresOptimizerType& UDFReconstruction<
	SURFACE, POINTS, NONMANIFOLD, RaySamplerTag>::optimizer_for()
{
	typename SpheresOptimizerType::Data d;
	d.samples_mesh = data_.samples;
	d.spheres = data_.spheres;
	d.sample_position = sample_position_;
	d.sample_area = sample_area_;
	d.sample_knn = sample_knn_;
	d.sample_quadric = sample_quadric_;
	d.sample_line_quadric = sample_line_quadric_;
	d.sample_ma_position = sample_ma_position_;
	d.sample_ma_radius = sample_ma_radius_;
	d.sample_ma_secondary_vertex = sample_ma_secondary_;
	d.sample_sphere = sample_sphere_;
	d.sample_error = sample_error_;
	d.sphere_position = sphere_position_;
	d.sphere_radius = sphere_radius_;
	d.sphere_cluster = sphere_cluster_;
	d.sphere_cluster_area = sphere_cluster_area_;
	d.sphere_cluster_color = sphere_cluster_color_;
	d.sphere_neighbors = sphere_neighbors_;
	d.sphere_error = sphere_error_;
	d.sphere_error_not_normalized = sphere_error_raw_;
	d.sample_kdtree = sample_kdtree_.get();
	d.sample_kdtree_vertices = &sample_vertices_;
	d.sqem_update_lambda_line_plane = Scalar(options_.sqem_update_lambda_line_plane);
	d.sphere_count = &sphere_count_;
	if (!optimizer_)
		optimizer_ = std::make_unique<SpheresOptimizerType>(std::move(d));
	else
		optimizer_->data() = std::move(d);
	return *optimizer_;
}

template <typename SURFACE, typename POINTS, typename NONMANIFOLD, typename RaySamplerTag>
typename UDFReconstruction<SURFACE, POINTS, NONMANIFOLD, RaySamplerTag>::TopologyType& UDFReconstruction<
	SURFACE, POINTS, NONMANIFOLD, RaySamplerTag>::topology_for()
{
	typename TopologyType::Data d;
	d.skeleton = data_.skeleton;
	d.spheres_optimizer = &optimizer_for();
	d.sphere_skeleton_vertex = sphere_skeleton_vertex_;
	d.skeleton_position = skeleton_position_;
	d.skeleton_radius = skeleton_radius_;
	d.skeleton_source_sphere = skeleton_source_sphere_;
	d.face_incident_tets = skeleton_incident_tets_;
	d.face_component_id = skeleton_component_id_;
	d.edge_degree = skeleton_edge_degree_;
	if (!topology_)
		topology_ = std::make_unique<TopologyType>(std::move(d));
	else
		topology_->data() = std::move(d);
	return *topology_;
}

template <typename SURFACE, typename POINTS, typename NONMANIFOLD, typename RaySamplerTag>
typename UDFReconstruction<SURFACE, POINTS, NONMANIFOLD, RaySamplerTag>::Status UDFReconstruction<
	SURFACE, POINTS, NONMANIFOLD, RaySamplerTag>::sample_alpha_level_set()
{
	if (initialize_data() != Status::success)
		return Status::invalid_data;
	build_input_kdtree();
	if constexpr (std::is_same_v<RaySamplerTag, RaySamplerPointCloud>)
	{
		if (!input_kdtree_ || !input_normal_ || !input_knn_)
			return Status::invalid_data;
		compute_input_normals(*data_.input_points, *input_position_, *input_normal_, *input_knn_, *input_kdtree_,
							  input_vertices_, options_.knn_k);
	}
	else if constexpr (std::is_same_v<RaySamplerTag, RaySamplerSurface>)
		build_surface_bvh();

	Vec3 bbox_min(0, 0, 0), bbox_max(1, 1, 1);
	if constexpr (std::is_same_v<RaySamplerTag, RaySamplerNeural>)
	{
		if (!neural_model_loaded_)
			return Status::model_load_failed;
		if (options_.neural_model_type == NeuralModelType::udf)
		{
			bbox_min = Vec3(-0.5, -0.5, -0.5);
			bbox_max = Vec3(0.5, 0.5, 0.5);
		}
		if (data_.surface && nb_cells<SVertex>(*data_.surface) > 0)
		{
			auto position = get_attribute<Vec3, SVertex>(*data_.surface, "position");
			if (position)
				std::tie(bbox_min, bbox_max) = geometry::bounding_box(*position);
		}
		else if (data_.input_points && nb_cells<PVertex>(*data_.input_points) > 0)
			std::tie(bbox_min, bbox_max) = geometry::bounding_box(*input_position_);
	}
	else if constexpr (std::is_same_v<RaySamplerTag, RaySamplerSurface>)
	{
		if (!data_.surface || !surface_bvh_)
			return Status::invalid_data;
		auto s_pos = get_attribute<Vec3, SVertex>(*data_.surface, "position");
		std::tie(bbox_min, bbox_max) = geometry::bounding_box(*s_pos);
	}
	else if (data_.input_points && nb_cells<PVertex>(*data_.input_points) > 0)
		std::tie(bbox_min, bbox_max) = geometry::bounding_box(*input_position_);
	const Scalar bbox_expand = Scalar(0.1);
	bbox_min -= Vec3(bbox_expand, bbox_expand, bbox_expand);
	bbox_max += Vec3(bbox_expand, bbox_expand, bbox_expand);

	RayTraits traits;
	// NeuralFieldRayTraits stores a pointer, so keep this object alive through
	// the complete sample_alpha_level_set call below.
	NeuralFieldForward sampling_field = neural_field();
	if constexpr (std::is_same_v<RaySamplerTag, RaySamplerNeural>)
		traits = RayConfig::make(sampling_field);
	else if constexpr (std::is_same_v<RaySamplerTag, RaySamplerSurface>)
	{
		auto s_pos = get_attribute<Vec3, SVertex>(*data_.surface, "position");
		traits = RayConfig::make(*data_.surface, s_pos.get(), surface_bvh_.get(), &surface_faces_);
	}
	else
		traits = RayConfig::make(*data_.input_points, input_position_.get(), input_normal_.get(), input_knn_.get(),
								 input_kdtree_.get(), &input_vertices_);
	if (!traits.is_ready())
		return Status::invalid_data;
	RayParams ray_params;
	ray_params.bbox_expand = Scalar(0.05);
	ray_params.alpha = Scalar(options_.alpha);
	ray_params.tol = Scalar(options_.tolerance);
	ray_params.step_bound = Scalar(2.0);
	ray_params.batch_size = std::max(1, options_.ray_sampler_batch_size);
	ray_params.max_iterations = options_.udf_max_iterations;
	ray_params.max_outer_iterations = 10;
	ray_params.seed = options_.seed;
	torch::Device sampling_device = torch::kCPU;
	if constexpr (std::is_same_v<RaySamplerTag, RaySamplerNeural>)
		sampling_device = device_;
	if (!ray_sampler_)
		ray_sampler_ = std::make_unique<RaySampler>(ray_params, sampling_device);
	else
	{
		ray_sampler_->set_params(ray_params);
		ray_sampler_->set_device(sampling_device);
	}
	AlphaSamplingParameters params;
	params.sample_radius = Scalar(options_.sample_radius);
	params.seed = static_cast<uint32>(options_.seed);
	params.ray_sampler_batch_size = static_cast<size_t>(std::max(1, options_.ray_sampler_batch_size));
	params.bbox_min = bbox_min;
	params.bbox_max = bbox_max;
	auto compute_normals = [&](const std::vector<Vec3>& points, std::vector<Vec3>& normals,
							   std::vector<uint8_t>& keep) {
		keep.assign(points.size(), uint8_t(0));
		normals.assign(points.size(), Vec3(0, 0, 1));
		if constexpr (std::is_same_v<RaySamplerTag, RaySamplerNeural>)
		{
			if (!evaluate_udf_normals(neural_field(), points, static_cast<size_t>(options_.batch_size),
									  projection_workspace_, normals))
				return false;
			std::fill(keep.begin(), keep.end(), uint8_t(1));
		}
		else
		{
			for (size_t i = 0; i < points.size(); ++i)
			{
				AlphaProjectionInfo info;
				const bool query_ok = [&]() {
					if constexpr (std::is_same_v<RaySamplerTag, RaySamplerSurface>)
						return query_surface_projection_info(points[i], info);
					else
						return query_point_cloud_projection_info(points[i], info);
				}();
				if (!query_ok || !info.valid || !points[i].allFinite() || !info.normal.allFinite())
					continue;
				Vec3 n = info.normal;
				if (n.squaredNorm() < Scalar(1e-12))
					n = Vec3(0, 0, 1);
				else
					n.normalize();
				normals[i] = n;
				keep[i] = uint8_t(1);
			}
		}
		return true;
	};
	auto project_points = [&](const std::vector<Vec3>& points, std::vector<Vec3>& projected, std::vector<Vec3>& normals,
							  std::vector<uint8_t>& keep) {
		if constexpr (std::is_same_v<RaySamplerTag, RaySamplerNeural>)
		{
			AlphaProjectionResult projection;
			const AlphaProjectionParameters pp{Scalar(options_.alpha), Scalar(options_.tolerance), bbox_min, bbox_max,
											   1};
			const bool filter_positive = options_.neural_model_type == NeuralModelType::mf;
			if (!geometry::project_points_to_alpha(neural_field(), points, pp, filter_positive, projection_workspace_,
												   projection))
				return false;
			projected = std::move(projection.projected_points);
			normals = std::move(projection.normals);
			keep = std::move(projection.keep_mask);
			return true;
		}
		AlphaProjectionResult projection;
		const AlphaProjectionParameters pp{Scalar(options_.alpha), Scalar(options_.tolerance), bbox_min, bbox_max, 1};
		const bool ok = geometry::project_points_to_alpha(
			points, pp,
			[&](const Vec3& p, AlphaProjectionInfo& info) {
				if constexpr (std::is_same_v<RaySamplerTag, RaySamplerSurface>)
					return query_surface_projection_info(p, info);
				else
					return query_point_cloud_projection_info(p, info);
			},
			projection);
		projected = std::move(projection.projected_points);
		normals = std::move(projection.normals);
		keep = std::move(projection.keep_mask);
		return ok;
	};
	AlphaSamplingResult result =
		geometry::sample_alpha_level_set(traits, *ray_sampler_, params, compute_normals, project_points);
	if (!result.success || result.positions.empty())
		return Status::sampling_failed;
	clear(*data_.samples);
	samples_before_filtering_ = static_cast<uint32>(result.positions.size());
	for (size_t i = 0; i < result.positions.size(); ++i)
	{
		PVertex v = add_vertex(*data_.samples);
		const uint32 idx = index_of(*data_.samples, v);
		(*sample_position_)[idx] = result.positions[i];
		(*sample_normal_)[idx] = result.normals[i];
	}
	build_sample_kdtree();
	return Status::success;
}

template <typename SURFACE, typename POINTS, typename NONMANIFOLD, typename RaySamplerTag>
typename UDFReconstruction<SURFACE, POINTS, NONMANIFOLD, RaySamplerTag>::Status UDFReconstruction<
	SURFACE, POINTS, NONMANIFOLD, RaySamplerTag>::apply_sampling_filtering()
{
	if (!options_.apply_filtering || !data_.samples || nb_cells<PVertex>(*data_.samples) == 0)
		return Status::success;
	build_input_kdtree();
	if (!input_kdtree_)
		build_surface_bvh();
	std::vector<Vec3> kept;
	kept.reserve(nb_cells<PVertex>(*data_.samples));
	foreach_cell(*data_.samples, [&](PVertex v) {
		const Vec3 p = (*sample_position_)[index_of(*data_.samples, v)];
		bool within_distance = false;
		const Scalar max_distance = Scalar(options_.alpha) + Scalar(1e-2);
		if (input_kdtree_)
		{
			std::pair<uint32, Scalar> nearest;
			within_distance = input_kdtree_->find_nn(p, &nearest, max_distance);
		}
		else if (surface_bvh_)
		{
			std::pair<uint32, Vec3> closest;
			within_distance = surface_bvh_->closest_point(p, &closest, max_distance);
		}
		if (within_distance)
			kept.push_back(p);
		return true;
	});
	const size_t before = nb_cells<PVertex>(*data_.samples);
	if (kept.size() == before)
		return Status::success;
	clear(*data_.samples);
	for (const Vec3& p : kept)
	{
		PVertex v = add_vertex(*data_.samples);
		(*sample_position_)[index_of(*data_.samples, v)] = p;
	}
	fitting_data_computed_ = false;
	if (kept.empty())
	{
		sample_kdtree_.reset();
		sample_vertices_.clear();
		return Status::success;
	}
	recompute_samples_normals_from_input();
	build_sample_kdtree();
	return Status::success;
}

template <typename SURFACE, typename POINTS, typename NONMANIFOLD, typename RaySamplerTag>
typename UDFReconstruction<SURFACE, POINTS, NONMANIFOLD, RaySamplerTag>::Status UDFReconstruction<
	SURFACE, POINTS, NONMANIFOLD, RaySamplerTag>::build_kdtree_and_normals()
{
	if (!data_.samples || nb_cells<PVertex>(*data_.samples) == 0)
		return Status::invalid_data;
	build_sample_kdtree();
	if (!sample_kdtree_)
		return Status::invalid_data;
	if (options_.recompute_normals_after_sampling)
	{
		if constexpr (std::is_same_v<RaySamplerTag, RaySamplerNeural>)
		{
			std::vector<Vec3> positions, normals;
			foreach_cell(*data_.samples, [&](PVertex v) {
				positions.push_back((*sample_position_)[index_of(*data_.samples, v)]);
				return true;
			});
			if (!evaluate_udf_normals(neural_field(), positions, static_cast<size_t>(options_.batch_size),
									  projection_workspace_, normals))
				return Status::fitting_failed;
			for (uint32 i = 0; i < normals.size(); ++i)
				(*sample_normal_)[i] = normals[i];
		}
		else
			parallel_foreach_cell(*data_.samples, [&](PVertex v) {
				recompute_sample_normal_pca(v);
				return true;
			});
	}
	return Status::success;
}

template <typename SURFACE, typename POINTS, typename NONMANIFOLD, typename RaySamplerTag>
typename UDFReconstruction<SURFACE, POINTS, NONMANIFOLD, RaySamplerTag>::Status UDFReconstruction<
	SURFACE, POINTS, NONMANIFOLD, RaySamplerTag>::compute_fitting_primitives()
{
	if (!data_.samples || !sample_kdtree_)
		return Status::invalid_data;
	compute_samples_area(*data_.samples, *sample_position_, *sample_normal_, *sample_area_, *sample_knn_,
						 *sample_kdtree_, sample_vertices_, options_.knn_k, Scalar(options_.alpha));
	compute_quadrics(*data_.samples, *sample_position_, *sample_normal_, *sample_area_, *sample_knn_, *sample_quadric_,
					 *sample_line_quadric_, options_.knn_k);
	fitting_data_computed_ = true;
	return Status::success;
}

template <typename SURFACE, typename POINTS, typename NONMANIFOLD, typename RaySamplerTag>
typename UDFReconstruction<SURFACE, POINTS, NONMANIFOLD, RaySamplerTag>::Status UDFReconstruction<
	SURFACE, POINTS, NONMANIFOLD, RaySamplerTag>::compute_initial_medial_axis()
{
	if (!data_.samples || !sample_kdtree_)
		return Status::invalid_data;
	const Scalar fallback_radius = Scalar(options_.alpha);
	const Scalar min_norm = Scalar(1e-12);
	const bool neural_input = std::is_same_v<RaySamplerTag, RaySamplerNeural> && neural_model_loaded_;
	const bool mf_model = neural_input && options_.neural_model_type == NeuralModelType::mf;
	const bool udf_model = neural_input && options_.neural_model_type == NeuralModelType::udf;
	const Scalar initial_radius = std::max<Scalar>(fallback_radius * Scalar(10), Scalar(0));
	auto run_shrinking_ball_for_vertex = [&](PVertex v) -> bool {
		if (!v.is_valid())
			return false;
		const uint32 idx = index_of(*data_.samples, v);
		if (idx == INVALID_INDEX)
			return false;
		const Vec3& pt = (*sample_position_)[idx];
		Vec3 n = (*sample_normal_)[idx];
		const bool normal_finite = n.allFinite();
		const bool can_shrink = normal_finite && n.squaredNorm() > min_norm;
		if (can_shrink)
			n.normalize();
		else
			n = Vec3(0, 0, 1);
		Vec3 center = pt - n * fallback_radius;
		Scalar radius = fallback_radius;
		PVertex secondary;
		if (can_shrink)
		{
			auto [c1, r1, q1] =
				geometry::shrinking_ball_center<PVertex>(pt, n, sample_kdtree_.get(), sample_vertices_, initial_radius);
			if (c1.allFinite() && std::isfinite(static_cast<double>(r1)) && r1 > Scalar(0))
			{
				center = c1;
				radius = r1;
				secondary = q1;
			}
		}
		if (!center.allFinite())
			center = pt - n * fallback_radius;
		if (!std::isfinite(static_cast<double>(radius)) || radius <= Scalar(0))
			radius = fallback_radius;
		if (!secondary.is_valid() && sample_kdtree_ && !sample_vertices_.empty())
		{
			std::pair<uint32, Scalar> nn;
			if (sample_kdtree_->find_nn(center - n * radius, &nn) && nn.first < sample_vertices_.size())
				secondary = sample_vertices_[nn.first];
		}
		(*sample_ma_position_)[idx] = center;
		(*sample_ma_radius_)[idx] = std::max(radius, fallback_radius);
		(*sample_ma_secondary_)[idx] = secondary;
		return true;
	};
	auto run_all = [&] {
		parallel_foreach_cell(*data_.samples, [&](PVertex v) {
			run_shrinking_ball_for_vertex(v);
			return true;
		});
	};
	run_all();
	if (options_.ma_flip_prune)
	{
		uint32 flipped = 0, deleted = 0;
		std::vector<uint32> flipped_indices;
		const Scalar udf_threshold = udf_model ? Scalar(1.2) * fallback_radius : fallback_radius;
		const Scalar mf_threshold =
			fallback_radius * std::max(Scalar(1), std::min(Scalar(5), Scalar(options_.ma_flip_prune_alpha_factor)));
		const bool supported = mf_model || std::is_same_v<RaySamplerTag, RaySamplerSurface> ||
							   std::is_same_v<RaySamplerTag, RaySamplerPointCloud> || neural_input;
		auto scores = [&](const std::vector<Vec3>& points, std::vector<Scalar>& values,
						  std::vector<Scalar>& sdf) -> bool {
			sdf.clear();
			if (!mf_model)
				return evaluate_topology_scores(points, values);
			NeuralFieldForward field = neural_field();
			if (!field.is_loaded())
				return false;
			values.assign(points.size(), Scalar(0));
			sdf.assign(points.size(), Scalar(0));
			const size_t batch = std::max<size_t>(1, static_cast<size_t>(options_.batch_size));
			for (size_t offset = 0; offset < points.size(); offset += batch)
			{
				const size_t count = std::min(batch, points.size() - offset);
				auto cpu = torch::TensorOptions().dtype(torch::kFloat32).device(torch::kCPU);
				if (device_.is_cuda())
					cpu = cpu.pinned_memory(true);
				torch::Tensor input = torch::empty({static_cast<long>(count), 3}, cpu);
				auto acc = input.accessor<float, 2>();
				for (size_t i = 0; i < count; ++i)
					for (int j = 0; j < 3; ++j)
						acc[static_cast<long>(i)][j] = static_cast<float>(points[offset + i][j]);
				auto [v, s] = field.forward_values_sdf_gpu(input);
				if (!v.defined() || !s.defined() || v.numel() != static_cast<long>(count) ||
					s.numel() != static_cast<long>(count))
					return false;
				if (v.dim() == 2 && v.size(1) == 1)
					v = v.squeeze(1);
				if (s.dim() == 2 && s.size(1) == 1)
					s = s.squeeze(1);
				v = v.to(torch::kCPU).contiguous();
				s = s.to(torch::kCPU).contiguous();
				auto va = v.accessor<float, 1>();
				auto sa = s.accessor<float, 1>();
				for (size_t i = 0; i < count; ++i)
				{
					values[offset + i] = Scalar(va[static_cast<long>(i)]);
					sdf[offset + i] = Scalar(sa[static_cast<long>(i)]);
				}
			}
			return true;
		};
		if (supported)
		{
			std::vector<PVertex> vertices;
			std::vector<Vec3> centers;
			vertices.reserve(nb_cells<PVertex>(*data_.samples));
			centers.reserve(vertices.capacity());
			foreach_cell(*data_.samples, [&](PVertex v) {
				const uint32 i = index_of(*data_.samples, v);
				vertices.push_back(v);
				centers.push_back((*sample_ma_position_)[i]);
				return true;
			});
			std::vector<Scalar> values, sdf;
			if (scores(centers, values, sdf) && values.size() == centers.size())
			{
				std::vector<PVertex> retry;
				for (size_t i = 0; i < vertices.size(); ++i)
					if ((mf_model && i < sdf.size() && sdf[i] > Scalar(0)) ||
						values[i] > (mf_model ? mf_threshold : udf_threshold))
						retry.push_back(vertices[i]);
				flipped_indices.reserve(retry.size());
				for (PVertex v : retry)
				{
					const uint32 i = index_of(*data_.samples, v);
					if (i == INVALID_INDEX)
						continue;
					Vec3 n = (*sample_normal_)[i];
					if (n.squaredNorm() > min_norm)
					{
						n.normalize();
						(*sample_normal_)[i] = -n;
						flipped_indices.push_back(i);
						++flipped;
					}
					run_shrinking_ball_for_vertex(v);
				}
				auto prune = [&](const std::vector<PVertex>& candidates) {
					std::vector<Vec3> p;
					std::vector<PVertex> active;
					for (PVertex v : candidates)
					{
						uint32 i = index_of(*data_.samples, v);
						if (i != INVALID_INDEX)
						{
							active.push_back(v);
							p.push_back((*sample_ma_position_)[i]);
						}
					}
					std::vector<Scalar> val, sd;
					if (!scores(p, val, sd) || val.size() != p.size())
						return uint32(0);
					uint32 removed = 0;
					for (size_t i = 0; i < active.size(); ++i)
						if ((mf_model && i < sd.size() && sd[i] > Scalar(0)) ||
							val[i] > (mf_model ? mf_threshold : udf_threshold))
						{
							remove_vertex(*data_.samples, active[i]);
							++removed;
						}
					return removed;
				};
				deleted += prune(retry);
				if (deleted > 0 && !flipped_indices.empty() && nb_cells<PVertex>(*data_.samples) > 0)
				{
					build_sample_kdtree();
					std::vector<PVertex> survivors;
					survivors.reserve(flipped_indices.size());
					for (uint32 i : flipped_indices)
					{
						PVertex v = of_index<PVertex>(*data_.samples, i);
						if (v.is_valid())
							survivors.push_back(v);
					}
					for (PVertex v : survivors)
						recompute_sample_normal_pca(v);
					for (PVertex v : survivors)
						run_shrinking_ball_for_vertex(v);
					deleted += prune(survivors);
				}
			}
		}
	}
	build_sample_kdtree();
	return Status::success;
}

template <typename SURFACE, typename POINTS, typename NONMANIFOLD, typename RaySamplerTag>
Vec3 UDFReconstruction<SURFACE, POINTS, NONMANIFOLD, RaySamplerTag>::hsv_to_rgb(Scalar h, Scalar s, Scalar v) const
{
	const Scalar hh = h - std::floor(h);
	if (!(s > Scalar(0)))
		return Vec3(v, v, v);
	const Scalar scaled = hh * Scalar(6);
	const int sector = static_cast<int>(std::floor(scaled));
	const Scalar f = scaled - Scalar(sector), p = v * (Scalar(1) - s), q = v * (Scalar(1) - s * f),
				 t = v * (Scalar(1) - s * (Scalar(1) - f));
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

template <typename SURFACE, typename POINTS, typename NONMANIFOLD, typename RaySamplerTag>
Vec3 UDFReconstruction<SURFACE, POINTS, NONMANIFOLD, RaySamplerTag>::component_palette_color(uint32 id) const
{
	static constexpr double golden_ratio_conjugate = 0.6180339887498949;
	const Scalar hue = Scalar(std::fmod(0.17 + golden_ratio_conjugate * static_cast<double>(id), 1.0));
	const Scalar saturation = (id % 3 == 0) ? Scalar(0.82) : (id % 3 == 1) ? Scalar(0.72) : Scalar(0.90);
	const Scalar value = (id % 4 == 0)	 ? Scalar(0.95)
						 : (id % 4 == 1) ? Scalar(0.88)
						 : (id % 4 == 2) ? Scalar(0.80)
										 : Scalar(0.92);
	return hsv_to_rgb(hue, saturation, value);
}

template <typename SURFACE, typename POINTS, typename NONMANIFOLD, typename RaySamplerTag>
typename UDFReconstruction<SURFACE, POINTS, NONMANIFOLD, RaySamplerTag>::Status UDFReconstruction<
	SURFACE, POINTS, NONMANIFOLD, RaySamplerTag>::initialize_spheres()
{
	if (!fitting_data_computed_ || !data_.samples || !data_.spheres)
		return Status::invalid_data;
	if (!optimizer_for().initialize_from_samples(Scalar(options_.init_dilation_constant)))
		return Status::sphere_initialization_failed;
	return Status::success;
}

template <typename SURFACE, typename POINTS, typename NONMANIFOLD, typename RaySamplerTag>
typename UDFReconstruction<SURFACE, POINTS, NONMANIFOLD, RaySamplerTag>::Status UDFReconstruction<
	SURFACE, POINTS, NONMANIFOLD, RaySamplerTag>::optimize_spheres(const IterationCallback& callback)
{
	auto& optimizer = optimizer_for();
	optimizer.reset_optimization();
	typename SpheresOptimizerType::Status status = SpheresOptimizerType::Status::running;
	while (status == SpheresOptimizerType::Status::running)
	{
		status = optimizer.update_once(Scalar(options_.sqem_update_lambda_line_plane));
		if (callback)
			callback(optimizer.metrics());
	}
	return status == SpheresOptimizerType::Status::failed ? Status::optimization_failed : Status::success;
}

template <typename SURFACE, typename POINTS, typename NONMANIFOLD, typename RaySamplerTag>
typename UDFReconstruction<SURFACE, POINTS, NONMANIFOLD, RaySamplerTag>::Status UDFReconstruction<
	SURFACE, POINTS, NONMANIFOLD, RaySamplerTag>::build_skeleton()
{
	face_components_prepared_ = false;
	return topology_for().build_from_spheres() == TopologyType::Status::success ? Status::success
																				: Status::skeleton_build_failed;
}

template <typename SURFACE, typename POINTS, typename NONMANIFOLD, typename RaySamplerTag>
typename UDFReconstruction<SURFACE, POINTS, NONMANIFOLD, RaySamplerTag>::Status UDFReconstruction<
	SURFACE, POINTS, NONMANIFOLD, RaySamplerTag>::fix_topology()
{
	face_components_prepared_ = false;
	return topology_for().fix_topology([this](const std::vector<Vec3>& p, std::vector<Scalar>& v) {
		return evaluate_topology_scores(p, v);
	}) == TopologyType::Status::success
			   ? Status::success
			   : Status::topology_fix_failed;
}

template <typename SURFACE, typename POINTS, typename NONMANIFOLD, typename RaySamplerTag>
typename UDFReconstruction<SURFACE, POINTS, NONMANIFOLD, RaySamplerTag>::Status UDFReconstruction<
	SURFACE, POINTS, NONMANIFOLD, RaySamplerTag>::prune_residual_sheets()
{
	face_components_prepared_ = false;
	return topology_for().prune_residual_sheets() == TopologyType::Status::success ? Status::success
																				   : Status::topology_fix_failed;
}

template <typename SURFACE, typename POINTS, typename NONMANIFOLD, typename RaySamplerTag>
typename UDFReconstruction<SURFACE, POINTS, NONMANIFOLD, RaySamplerTag>::Status UDFReconstruction<
	SURFACE, POINTS, NONMANIFOLD, RaySamplerTag>::prepare_face_components()
{
	face_components_prepared_ = false;
	auto& topology = topology_for();
	// Keep the legacy face-component export contract explicit: fully non-manifold
	// triangles are removed before component ids and colors are computed.
	if (topology.prune_fully_non_manifold_triangles() != TopologyType::Status::success ||
		topology.compute_skeleton_face_components_union_find() != TopologyType::Status::success)
		return Status::topology_fix_failed;
	skeleton_face_component_color_ = get_or_add_attribute<Vec3, NMFace>(*data_.skeleton, "face_component_color");
	const uint32 count = nb_cells<NMFace>(*data_.skeleton);
	if (count == 0)
	{
		face_components_prepared_ = true;
		return Status::success;
	}
	foreach_cell(*data_.skeleton, [&](NMFace f) {
		const uint32 idx = index_of(*data_.skeleton, f);
		const uint32 id = (*skeleton_component_id_)[idx];
		if (id == INVALID_INDEX)
		{
			(*skeleton_face_component_color_)[idx] = Vec3(0, 0, 0);
			return true;
		}
		(*skeleton_face_component_color_)[idx] = component_palette_color(id);
		return true;
	});
	face_components_prepared_ = true;
	return Status::success;
}

template <typename SURFACE, typename POINTS, typename NONMANIFOLD, typename RaySamplerTag>
bool UDFReconstruction<SURFACE, POINTS, NONMANIFOLD, RaySamplerTag>::evaluate_topology_scores(
	const std::vector<Vec3>& points, std::vector<Scalar>& values)
{
	if constexpr (std::is_same_v<RaySamplerTag, RaySamplerNeural>)
	{
		if (options_.neural_model_type == NeuralModelType::mf)
		{
			if (!sample_ma_kdtree_)
				return false;
			values.resize(points.size(), Scalar(0));
			for (size_t i = 0; i < points.size(); ++i)
			{
				std::pair<uint32, Scalar> nearest;
				if (!sample_ma_kdtree_->find_nn(points[i], &nearest) || nearest.first >= sample_ma_vertices_.size())
					continue;
				const uint32 sid = index_of(*data_.samples, sample_ma_vertices_[nearest.first]);
				if (sid != INVALID_INDEX)
					values[i] = (points[i] - (*sample_ma_position_)[sid]).norm();
			}
			return true;
		}
		return evaluate_udf_values(neural_field(), points, static_cast<size_t>(options_.batch_size), values);
	}
	else if constexpr (std::is_same_v<RaySamplerTag, RaySamplerSurface>)
	{
		build_surface_bvh();
		if (!surface_bvh_)
			return false;
		values.resize(points.size());
		for (size_t i = 0; i < points.size(); ++i)
		{
			uint32 primitive;
			Vec3 closest;
			Scalar distance;
			if (!query_surface_closest_point(*surface_bvh_, points[i], primitive, closest, distance))
				return false;
			values[i] = distance;
		}
		return true;
	}
	else
	{
		build_input_kdtree();
		if (!input_kdtree_)
			return false;
		values.resize(points.size());
		for (size_t i = 0; i < points.size(); ++i)
		{
			uint32 point;
			Scalar distance;
			if (!query_point_cloud_nearest_point(*input_kdtree_, points[i], point, distance))
				return false;
			values[i] = distance;
		}
		return true;
	}
}

template <typename SURFACE, typename POINTS, typename NONMANIFOLD, typename RaySamplerTag>
typename UDFReconstruction<SURFACE, POINTS, NONMANIFOLD, RaySamplerTag>::Result UDFReconstruction<
	SURFACE, POINTS, NONMANIFOLD, RaySamplerTag>::run(const IterationCallback& callback)
{
	Result r;
	if (data_.samples)
		clear(*data_.samples);
	if (data_.spheres)
		clear(*data_.spheres);
	if (data_.skeleton)
		clear(*data_.skeleton);
	sample_kdtree_.reset();
	sample_ma_kdtree_.reset();
	sample_vertices_.clear();
	sample_ma_vertices_.clear();
	optimizer_.reset();
	topology_.reset();
	fitting_data_computed_ = false;
	face_components_prepared_ = false;
	samples_before_filtering_ = 0;
	sphere_count_ = 0;
	auto run_stage = [](float64& time, auto&& stage) {
		const auto start = std::chrono::high_resolution_clock::now();
		const Status status = stage();
		time = std::chrono::duration<float64, std::milli>(std::chrono::high_resolution_clock::now() - start).count();
		return status;
	};
	if ((r.status = run_stage(r.timing.preprocessing_ms, [&] { return initialize_data(); })) != Status::success)
		return r;
	if ((r.status = run_stage(r.timing.sampling_ms, [&] { return sample_alpha_level_set(); })) != Status::success)
		return r;
	r.counts.samples_before_filtering = samples_before_filtering_;
	if ((r.status = run_stage(r.timing.sample_filtering_ms, [&] { return apply_sampling_filtering(); })) !=
		Status::success)
		return r;
	if ((r.status = run_stage(r.timing.kdtree_and_normals_ms, [&] { return build_kdtree_and_normals(); })) !=
		Status::success)
		return r;
	if ((r.status = run_stage(r.timing.fitting_primitives_ms, [&] { return compute_fitting_primitives(); })) !=
		Status::success)
		return r;
	if ((r.status = run_stage(r.timing.initial_medial_axis_ms, [&] { return compute_initial_medial_axis(); })) !=
		Status::success)
		return r;
	if ((r.status = run_stage(r.timing.sphere_initialization_ms, [&] { return initialize_spheres(); })) !=
		Status::success)
		return r;
	if ((r.status = run_stage(r.timing.optimization_ms, [&] { return optimize_spheres(callback); })) != Status::success)
		return r;
	if (optimizer_)
	{
		r.timing.cluster_ms = optimizer_->metrics().cluster_total_ms;
		r.timing.sphere_update_ms = optimizer_->metrics().sphere_update_total_ms;
		r.timing.error_ms = optimizer_->metrics().error_total_ms;
	}
	if ((r.status = run_stage(r.timing.skeleton_construction_ms, [&] { return build_skeleton(); })) != Status::success)
		return r;
	if ((r.status = run_stage(r.timing.topology_processing_ms, [&] { return fix_topology(); })) != Status::success)
		return r;
	if (options_.residual_prune)
	{
		float64 residual_ms = 0.0;
		if ((r.status = run_stage(residual_ms, [&] { return prune_residual_sheets(); })) != Status::success)
			return r;
		r.timing.topology_processing_ms += residual_ms;
	}
	r.status = Status::success;
	r.counts = counts();
	r.counts.samples_before_filtering = samples_before_filtering_;
	return r;
}

} // namespace cgogn::geometry

#endif // CGOGN_GEOMETRY_ALGOS_UDF_RECONSTRUCTION_H_
