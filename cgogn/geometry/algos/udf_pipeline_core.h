#ifndef CGOGN_GEOMETRY_ALGOS_UDF_PIPELINE_CORE_H_
#define CGOGN_GEOMETRY_ALGOS_UDF_PIPELINE_CORE_H_

#include <cgogn/core/functions/attributes.h>
#include <cgogn/geometry/algos/udf_pipeline_context.h>
#include <cgogn/geometry/functions/bounding_box.h>

#include <stdexcept>

namespace cgogn
{

namespace geometry
{

template <typename RaySamplerTag>
class UDFPipelineCore
{
public:
	using Context = UDFPipelineContext<RaySamplerTag>;
	using Surface = typename Context::Surface;
	using Points = typename Context::Points;
	using NonManifold = typename Context::NonManifold;
	using Training = typename Context::Training;
	using PVertex = typename mesh_traits<Points>::Vertex;
	using SVertex = typename mesh_traits<Surface>::Vertex;

	explicit UDFPipelineCore(Context& context) : context_(context)
	{
	}

	Points& load_point_cloud_input(const std::string& input_path, bool normalize_for_udf)
	{
		std::cerr << "[BenchmarkInput] load_point_cloud_input begin: " << input_path << std::endl;
		Points* points = context_.points_provider.load_points_from_file(input_path);
		std::cerr << "[BenchmarkInput] load_points_from_file returned" << std::endl;
		if (!points)
			throw std::runtime_error("Failed to load point cloud: " + input_path);

		auto position = get_attribute<Vec3, PVertex>(*points, "position");
		std::cerr << "[BenchmarkInput] get_attribute(position) returned" << std::endl;
		if (!position)
			throw std::runtime_error("Point cloud has no `position` attribute: " + input_path);

		if (normalize_for_udf)
		{
			geometry::normalize_centered(*position.get());
			context_.points_provider.emit_attribute_changed(*points, position.get());
		}

		std::cerr << "[BenchmarkInput] set_mesh_bb_vertex_position begin" << std::endl;
		context_.points_provider.set_mesh_bb_vertex_position(*points, position);
		std::cerr << "[BenchmarkInput] set_mesh_bb_vertex_position end" << std::endl;
		context_.selected_points = points;
		std::cerr << "[BenchmarkInput] set_selected_points begin" << std::endl;
		context_.udf_training.set_selected_points(*points);
		std::cerr << "[BenchmarkInput] set_selected_points end" << std::endl;
		return *points;
	}

	Surface& load_surface_input(const std::string& surface_path, bool normalize_for_udf)
	{
		Surface* surface = context_.surface_provider.load_surface_from_file(surface_path);
		if (!surface)
			throw std::runtime_error("Failed to load surface mesh: " + surface_path);
		auto position = get_attribute<Vec3, SVertex>(*surface, "position");
		if (!position)
			throw std::runtime_error("Surface mesh has no `position` attribute: " + surface_path);
		if (normalize_for_udf)
		{
			geometry::normalize_centered(*position.get());
			context_.surface_provider.emit_attribute_changed(*surface, position.get());
		}
		context_.surface_provider.set_mesh_bb_vertex_position(*surface, position);
		context_.selected_surface = surface;
		context_.udf_training.set_selected_surface(*surface);
		return *surface;
	}

	Points& create_empty_points_input(const std::string& mesh_name)
	{
		Points* points = context_.points_provider.add_mesh(mesh_name);
		if (!points)
			throw std::runtime_error("Failed to create point container: " + mesh_name);
		auto position = get_or_add_attribute<Vec3, PVertex>(*points, "position");
		context_.points_provider.set_mesh_bb_vertex_position(*points, position);
		context_.selected_points = points;
		context_.udf_training.set_selected_points(*points);
		return *points;
	}

	void load_neural_model(Points& points, const std::string& model_path, typename Training::NeuralModelType model_type)
	{
		context_.udf_training.load_neural_udf_model(points, model_path, model_type);
	}

	void apply_options(Points& points, const typename Training::HeadlessBenchmarkOptions& options)
	{
		context_.udf_training.apply_headless_benchmark_options(points, options);
	}

	void sample_alpha_level_set(Points& points, const typename Training::HeadlessBenchmarkOptions& options)
	{
		context_.udf_training.headless_sample_alpha_level_set(points, options);
	}

	void apply_poisson_eliminate(Points& points, size_t target_num)
	{
		context_.udf_training.headless_apply_poisson_eliminate(points, target_num);
	}

	void apply_sampling_filtering(Points& points)
	{
		context_.udf_training.headless_apply_sampling_filtering(points);
	}

	void build_kdtree_and_normals(Points& points)
	{
		context_.udf_training.headless_build_sample_kdtree_and_normals(points);
	}

	void compute_samples_area(Points& points)
	{
		context_.udf_training.headless_compute_samples_area(points);
	}

	void compute_winding_numbers(Points& points)
	{
		context_.udf_training.headless_compute_winding_numbers(points);
	}

	void compute_quadrics(Points& points)
	{
		context_.udf_training.headless_compute_quadrics(points);
	}

	void compute_initial_medial_axis(Points& points)
	{
		context_.udf_training.headless_compute_initial_medial_axis(points);
	}

	void init_spheres(Points& points, unsigned int max_nb_spheres)
	{
		context_.udf_training.headless_init_spheres(points, max_nb_spheres);
	}

	typename Training::HeadlessOptimizationStats optimize(Points& points, bool verbose)
	{
		return context_.udf_training.headless_optimize_spheres(points, verbose);
	}

	void build_skeleton(Points& points)
	{
		context_.udf_training.headless_build_skeleton(points);
	}

	void run_topology_fix(Points& points, bool run_deg_face_deletion)
	{
		context_.udf_training.headless_run_topology_fix(points, run_deg_face_deletion);
	}

	void run_deg_face_deletion(Points& points)
	{
		context_.udf_training.headless_run_deg_face_deletion(points);
	}

	void run_face_post_processing(Points& points)
	{
		context_.udf_training.headless_run_face_post_processing(points);
	}

	typename Training::HeadlessCounts collect_counts(const Points& points) const
	{
		return context_.udf_training.headless_collect_counts(points);
	}

	void export_skeleton_ply(Points& points, const std::string& filename)
	{
		context_.udf_training.headless_export_skeleton_ply(points, filename);
	}

private:
	Context& context_;
};

} // namespace geometry

} // namespace cgogn

#endif // CGOGN_GEOMETRY_ALGOS_UDF_PIPELINE_CORE_H_
