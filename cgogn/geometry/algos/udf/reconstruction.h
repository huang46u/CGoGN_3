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
#include <cgogn/core/utils/thread_pool.h>
#include <cgogn/geometry/algos/area.h>
#include <cgogn/geometry/algos/fitting.h>
#include <cgogn/geometry/algos/medial_axis.h>
#include <cgogn/geometry/algos/normal.h>
#include <cgogn/geometry/algos/udf/alpha_projection.h>
#include <cgogn/geometry/algos/udf/alpha_sampling.h>
#include <cgogn/geometry/algos/udf/neural_alpha_projection.h>
#include <cgogn/geometry/algos/udf/neural_field_query.h>
#include <cgogn/geometry/algos/udf/reconstruction_stats.h>
#include <cgogn/geometry/algos/udf/sample_processing.h>
#include <cgogn/geometry/algos/udf/skeleton_topology.h>
#include <cgogn/geometry/algos/udf/spatial_query.h>
#include <cgogn/geometry/algos/udf/neural_model_type.h>
#include <cgogn/geometry/algos/udf/normal_utils.h>
#include <cgogn/geometry/algos/udf/parallel_for.h>
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
#include <atomic>
#include <chrono>
#include <cmath>
#include <filesystem>
#include <functional>
#include <future>
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

// Fixed constants of the reconstruction (not user options). Lengths are absolute, in model units (the inputs are
// rescaled to [0,1] before the reconstruction).
namespace udf_constants
{
// Padding added on each side of the input bounding box. The padded box is the sampling domain and the box the
// projected points are clamped to.
inline constexpr Scalar sampling_padding = Scalar(0.1);
// Relative enlargement of the padded sampling box from which the seed rays are cast (its half size is multiplied by
// 1 + ray_box_relative_expand). Overrides the ray sampler default (0.1).
inline constexpr Scalar ray_box_relative_expand = Scalar(0.05);
// Samples farther than alpha + filter_margin from the input are removed by the sampling filter.
inline constexpr Scalar filter_margin = Scalar(1e-2);
// MAFlipPrune, neural UDF: a medial centre is retried when its UDF value exceeds this factor times alpha.
inline constexpr Scalar udf_flip_threshold_factor = Scalar(1.2);
// MAFlipPrune, MF: range the ma_flip_prune_alpha_factor option is clamped to.
inline constexpr Scalar mf_flip_alpha_factor_min = Scalar(1);
inline constexpr Scalar mf_flip_alpha_factor_max = Scalar(5);
} // namespace udf_constants

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

	using NeuralModelType = cgogn::geometry::NeuralModelType;
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
		bool apply_filtering = true;
		bool recompute_normals_after_sampling = false;
		int ray_sampler_batch_size = 8192;
		int batch_size = 1048576;
		float tolerance = 1e-6f;
		int udf_max_iterations = 3000;
		// neural seed rays: end a ray batch once fewer than this many rays are still marching (0 = off)
		int udf_min_active_rays = 100;
		float sample_radius = 0.0025f;
		float sqem_update_lambda_line_plane = 0.2f;
		int seed = 42;
		float init_dilation_constant = 0.001f;
		bool residual_prune = false;
		int residual_prune_mode = 2; // 0 current, 1 keep_fully_non_manifold, 2 collapsible_only
		int alpha_sampling_candidates_per_parent = 12;
		int alpha_sampling_parent_batch_size = 1024;
	};

	using Timing = ReconstructionTiming;
	using Counts = ReconstructionCounts;

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
		invalidate_input_caches();
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
		neural_field_ = NeuralFieldForward(&neural_model_, false, device_);
		fitting_data_computed_ = false;
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
		sphere_neighbors_ = get_or_add_attribute<std::vector<PVertex>, PVertex>(*data_.spheres, "neighbor_clusters");
		sphere_error_ = get_or_add_attribute<Scalar, PVertex>(*data_.spheres, "error");
		sphere_error_raw_ = get_or_add_attribute<Scalar, PVertex>(*data_.spheres, "error_not_normalized");
		sphere_skeleton_vertex_ = get_or_add_attribute<NMVertex, PVertex>(*data_.spheres, "skeleton_vertex");
		skeleton_position_ = get_or_add_attribute<Vec3, NMVertex>(*data_.skeleton, "position");
		skeleton_radius_ = get_or_add_attribute<Scalar, NMVertex>(*data_.skeleton, "radius");
		skeleton_source_sphere_ = get_or_add_attribute<PVertex, NMVertex>(*data_.skeleton, "source_sphere");
		skeleton_incident_tets_ = get_or_add_attribute<std::set<std::size_t>, NMFace>(*data_.skeleton, "incident_tets");
		skeleton_edge_degree_ = get_or_add_attribute<uint32, NMEdge>(*data_.skeleton, "degree");
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
				neural_field_ = NeuralFieldForward(&neural_model_, true, device_);
				options_.neural_model_type = type;
				return Status::success;
			}
			catch (const c10::Error&)
			{
				neural_model_loaded_ = false;
				neural_field_ = NeuralFieldForward(&neural_model_, false, device_);
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
	Result run(const IterationCallback& iteration_callback = {});

	Counts counts() const
	{
		Counts c;
		if (data_.input_points)
			c.input_points = nb_cells<PVertex>(*data_.input_points);
		if (data_.surface)
			c.input_vertices = nb_cells<SVertex>(*data_.surface);
		else
			c.input_vertices = c.input_points;
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

	bool export_skeleton_ply(const std::filesystem::path& path) const
	{
		if (!data_.skeleton || !skeleton_position_ || !skeleton_radius_)
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
	// one wrapper for the loaded model: it remembers the output layout found by its first forward, so MF models
	// do not go through autograd again on every call
	NeuralFieldForward neural_field_{&neural_model_, false, device_};
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
	std::shared_ptr<PointAttribute<std::vector<PVertex>>> sphere_neighbors_;
	std::shared_ptr<PointAttribute<Scalar>> sphere_error_, sphere_error_raw_;
	std::shared_ptr<PointAttribute<NMVertex>> sphere_skeleton_vertex_;
	std::shared_ptr<SkeletonAttribute<Vec3>> skeleton_position_;
	std::shared_ptr<SkeletonAttribute<Scalar>> skeleton_radius_;
	std::shared_ptr<SkeletonAttribute<PVertex>> skeleton_source_sphere_;
	std::shared_ptr<SkeletonAttribute<std::set<std::size_t>>> skeleton_incident_tets_;
	std::shared_ptr<SkeletonAttribute<uint32>> skeleton_edge_degree_;
	std::shared_ptr<SurfaceAttribute<Vec3>> surface_vertex_normal_;

	void invalidate_input_caches();
	void build_surface_bvh();
	void build_input_kdtree();
	void build_sample_kdtree();
	bool query_surface_projection_info(const Vec3& query, AlphaProjectionInfo& info);
	bool query_point_cloud_projection_info(const Vec3& query, AlphaProjectionInfo& info);
	void recompute_samples_normals_from_input();
	bool recompute_sample_normal_pca(PVertex vertex);
	NeuralFieldForward& neural_field()
	{
		return neural_field_;
	}
	SpheresOptimizerType& optimizer_for();
	TopologyType& topology_for();
	bool evaluate_topology_scores(const std::vector<Vec3>& points, std::vector<Scalar>& values);
};

template <typename SURFACE, typename POINTS, typename NONMANIFOLD, typename RaySamplerTag>
void UDFReconstruction<SURFACE, POINTS, NONMANIFOLD, RaySamplerTag>::invalidate_input_caches()
{
	input_kdtree_.reset();
	input_vertices_.clear();
	surface_bvh_.reset();
	surface_faces_.clear();
	surface_vertices_.clear();
	surface_positions_.clear();
	surface_vertex_normal_.reset();
}

template <typename SURFACE, typename POINTS, typename NONMANIFOLD, typename RaySamplerTag>
void UDFReconstruction<SURFACE, POINTS, NONMANIFOLD, RaySamplerTag>::build_input_kdtree()
{
	if (input_kdtree_)
		return;
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
	if (surface_bvh_)
		return;
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
	info.normal = unit_normal_or_z(info.normal);
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
	info.normal = unit_normal_or_z(info.normal);
	info.valid =
		info.closest_point.allFinite() && info.normal.allFinite() && std::isfinite(static_cast<double>(info.distance));
	return info.valid;
}

template <typename SURFACE, typename POINTS, typename NONMANIFOLD, typename RaySamplerTag>
void UDFReconstruction<SURFACE, POINTS, NONMANIFOLD, RaySamplerTag>::recompute_samples_normals_from_input()
{
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
					normal = unit_normal_or_z(position - closest.second);
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
					normal = unit_normal_or_z(normal);
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
	invalidate_input_caches();
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
	const Scalar padding = udf_constants::sampling_padding;
	bbox_min -= Vec3(padding, padding, padding);
	bbox_max += Vec3(padding, padding, padding);

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
	ray_params.bbox_expand = udf_constants::ray_box_relative_expand;
	ray_params.alpha = Scalar(options_.alpha);
	ray_params.tol = Scalar(options_.tolerance);
	ray_params.step_bound = Scalar(2.0);
	ray_params.batch_size = std::max(1, options_.ray_sampler_batch_size);
	ray_params.max_iterations = options_.udf_max_iterations;
	ray_params.min_active_rays = options_.udf_min_active_rays;
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
	params.candidates_per_parent = static_cast<size_t>(std::max(1, options_.alpha_sampling_candidates_per_parent));
	params.parent_batch_size = static_cast<size_t>(std::max(1, options_.alpha_sampling_parent_batch_size));
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
			detail::parallel_for_index(points.size(), [&](std::size_t i) {
				AlphaProjectionInfo info;
				const bool query_ok = [&]() {
					if constexpr (std::is_same_v<RaySamplerTag, RaySamplerSurface>)
						return query_surface_projection_info(points[i], info);
					else
						return query_point_cloud_projection_info(points[i], info);
				}();
				if (!query_ok || !info.valid || !points[i].allFinite() || !info.normal.allFinite())
					return;
				normals[i] = unit_normal_or_z(info.normal);
				keep[i] = uint8_t(1);
			});
		}
		return true;
	};
	auto project_points = [&](const std::vector<Vec3>& points, AlphaProjectionResult& projection) {
		if constexpr (std::is_same_v<RaySamplerTag, RaySamplerNeural>)
		{
			const AlphaProjectionParameters pp{Scalar(options_.alpha), Scalar(options_.tolerance), bbox_min, bbox_max, 1};
			const bool filter_positive = options_.neural_model_type == NeuralModelType::mf;
			if (!geometry::project_points_to_alpha(neural_field(), points, pp, filter_positive,
												  projection_workspace_, projection))
				return false;
			return true;
		}
		const AlphaProjectionParameters pp{Scalar(options_.alpha), Scalar(options_.tolerance), bbox_min, bbox_max, 1};
		return geometry::project_points_to_alpha(
			points, pp,
			[&](const Vec3& p, AlphaProjectionInfo& info) {
				if constexpr (std::is_same_v<RaySamplerTag, RaySamplerSurface>)
					return query_surface_projection_info(p, info);
				else
					return query_point_cloud_projection_info(p, info);
			},
			projection);
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
	std::vector<Vec3> kept, kept_normals;
	kept.reserve(nb_cells<PVertex>(*data_.samples));
	foreach_cell(*data_.samples, [&](PVertex v) {
		const Vec3 p = (*sample_position_)[index_of(*data_.samples, v)];
		bool within_distance = false;
		const Scalar max_distance = Scalar(options_.alpha) + udf_constants::filter_margin;
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
		{
			kept.push_back(p);
			kept_normals.push_back((*sample_normal_)[index_of(*data_.samples, v)]);
		}
		return true;
	});
	const size_t before = nb_cells<PVertex>(*data_.samples);
	if (kept.size() == before)
		return Status::success;
	clear(*data_.samples);
	for (size_t i = 0; i < kept.size(); ++i)
	{
		PVertex v = add_vertex(*data_.samples);
		(*sample_position_)[index_of(*data_.samples, v)] = kept[i];
		(*sample_normal_)[index_of(*data_.samples, v)] = kept_normals[i];
	}
	fitting_data_computed_ = false;
	if (kept.empty())
	{
		sample_kdtree_.reset();
		sample_vertices_.clear();
		return Status::success;
	}
	// neural samples keep their projection normals (the network gradient); for point cloud / surface input the
	// normal is recomputed at the projected position from the input
	if constexpr (!std::is_same_v<RaySamplerTag, RaySamplerNeural>)
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
	if (!sample_kdtree_)
		return Status::invalid_data;
	// PCA for every input: neural samples already carry the network gradient, and for MF models that gradient
	// is not orthogonal to the medial axis.
	if (options_.recompute_normals_after_sampling)
		parallel_foreach_cell(*data_.samples, [&](PVertex v) {
			recompute_sample_normal_pca(v);
			return true;
		});
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
		const Scalar udf_threshold =
			udf_model ? udf_constants::udf_flip_threshold_factor * fallback_radius : fallback_radius;
		const Scalar mf_threshold =
			fallback_radius * std::max(udf_constants::mf_flip_alpha_factor_min,
									   std::min(udf_constants::mf_flip_alpha_factor_max,
												Scalar(options_.ma_flip_prune_alpha_factor)));
		const bool supported = mf_model || std::is_same_v<RaySamplerTag, RaySamplerSurface> ||
							   std::is_same_v<RaySamplerTag, RaySamplerPointCloud> || neural_input;
		auto scores = [&](const std::vector<Vec3>& points, std::vector<Scalar>& values,
						  std::vector<Scalar>& sdf) -> bool {
			sdf.clear();
			if (!mf_model)
				return evaluate_topology_scores(points, values);
			return evaluate_udf_values_sdf(neural_field(), points, static_cast<size_t>(options_.batch_size), values,
										   sdf);
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
	return topology_for().build_from_spheres() == TopologyType::Status::success ? Status::success
																				: Status::skeleton_build_failed;
}

template <typename SURFACE, typename POINTS, typename NONMANIFOLD, typename RaySamplerTag>
typename UDFReconstruction<SURFACE, POINTS, NONMANIFOLD, RaySamplerTag>::Status UDFReconstruction<
	SURFACE, POINTS, NONMANIFOLD, RaySamplerTag>::fix_topology()
{
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
	topology_for().set_residual_prune_mode(static_cast<typename TopologyType::ResidualPruneMode>(options_.residual_prune_mode));
	return topology_for().prune_residual_sheets() == TopologyType::Status::success ? Status::success
																	   : Status::topology_fix_failed;
}

template <typename SURFACE, typename POINTS, typename NONMANIFOLD, typename RaySamplerTag>
bool UDFReconstruction<SURFACE, POINTS, NONMANIFOLD, RaySamplerTag>::evaluate_topology_scores(
	const std::vector<Vec3>& points, std::vector<Scalar>& values)
{
	if constexpr (std::is_same_v<RaySamplerTag, RaySamplerNeural>)
	{
		if (options_.neural_model_type == NeuralModelType::mf)
		{
			// MF: mf - |sdf| is not a stable medial axis (its minimum is above zero and varies a lot over a shape),
			// so the score is the distance to the medial point of the nearest sample instead of the network value
			if (!sample_ma_kdtree_)
				return false;
			values.resize(points.size(), Scalar(0));
			detail::parallel_for_index(points.size(), [&](std::size_t i) {
				std::pair<uint32, Scalar> nearest;
				if (!sample_ma_kdtree_->find_nn(points[i], &nearest) || nearest.first >= sample_ma_vertices_.size())
					return;
				const uint32 sid = index_of(*data_.samples, sample_ma_vertices_[nearest.first]);
				if (sid != INVALID_INDEX)
					values[i] = (points[i] - (*sample_ma_position_)[sid]).norm();
			});
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
		std::atomic<bool> all_found{true};
		detail::parallel_for_index(points.size(), [&](std::size_t i) {
			uint32 primitive;
			Vec3 closest;
			Scalar distance;
			if (!query_surface_closest_point(*surface_bvh_, points[i], primitive, closest, distance))
			{
				all_found = false;
				return;
			}
			values[i] = distance;
		});
		if (!all_found)
			return false;
		return true;
	}
	else
	{
		build_input_kdtree();
		if (!input_kdtree_)
			return false;
		values.resize(points.size());
		std::atomic<bool> all_found{true};
		detail::parallel_for_index(points.size(), [&](std::size_t i) {
			uint32 point;
			Scalar distance;
			if (!query_point_cloud_nearest_point(*input_kdtree_, points[i], point, distance))
			{
				all_found = false;
				return;
			}
			values[i] = distance;
		});
		if (!all_found)
			return false;
		return true;
	}
}

template <typename SURFACE, typename POINTS, typename NONMANIFOLD, typename RaySamplerTag>
typename UDFReconstruction<SURFACE, POINTS, NONMANIFOLD, RaySamplerTag>::Result UDFReconstruction<
	SURFACE, POINTS, NONMANIFOLD, RaySamplerTag>::run(const IterationCallback& callback)
{
	Result r;
	invalidate_input_caches();
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
	samples_before_filtering_ = 0;
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
	return r;
}

} // namespace cgogn::geometry

#endif // CGOGN_GEOMETRY_ALGOS_UDF_RECONSTRUCTION_H_
