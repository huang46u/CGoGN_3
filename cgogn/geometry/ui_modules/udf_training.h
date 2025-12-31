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
#include <cgogn/geometry/types/spherical_quadric.h>
#include <cgogn/geometry/types/vector_traits.h>

#include <cgogn/rendering/ui_modules/point_cloud_render.h>
#include <cgogn/rendering/ui_modules/surface_render.h>

#include <Eigen/Dense>
#include <Eigen/Eigenvalues>
#include <Eigen/Sparse>
#include <libacc/bvh_tree.h>
#include <libacc/bvh_tree_spheres.h>
#include <libacc/kd_tree.h>

//import CGAL
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Polygon_mesh_processing/distance.h>
#include <CGAL/Surface_mesh.h>

#include <GLFW/glfw3.h>
#include <numeric>
#include <algorithm>
#include <random>
#include <thread>
#include <mutex>
#include <array>
#include <set>
#include <unordered_map>

namespace cgogn
{

namespace ui
{

using geometry::Vec3;
using geometry::Vec4;
using geometry::Scalar;
using geometry::Spherical_Quadric;

template <typename SURFACE, typename POINTS, typename NONMANIFOLD>
class UDFTraining : public ViewModule
{
	using PVertex = typename mesh_traits<POINTS>::Vertex;
	using NMVertex = typename mesh_traits<NONMANIFOLD>::Vertex;
	using SVertex = typename mesh_traits<SURFACE>::Vertex;

	using SFace = typename mesh_traits<SURFACE>::Face;
	using NMEdge = typename mesh_traits<NONMANIFOLD>::Edge;
	
	template <typename T>
	using PAttribute = typename mesh_traits<POINTS>::template Attribute<T>;
	template <typename T>
	using NMAttribute = typename mesh_traits<NONMANIFOLD>::template Attribute<T>;
	using K = CGAL::Exact_predicates_inexact_constructions_kernel;
	using CGAL_Mesh = CGAL::Surface_mesh<K::Point_3>;
	using Vector_3 = typename K::Vector_3;
	using Point_3 = typename K::Point_3;
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

private:
	struct SpatialGrid {
		struct GridKey {
			int x, y, z;
			bool operator==(const GridKey& o) const { return x == o.x && y == o.y && z == o.z; }
		};
		struct GridKeyHash {
			std::size_t operator()(const GridKey& k) const {
				return std::hash<int>()(k.x) ^ (std::hash<int>()(k.y) << 1) ^ (std::hash<int>()(k.z) << 2);
			}
		};

		Scalar cell_size_;
		std::unordered_map<GridKey, std::vector<uint32>, GridKeyHash> grid_;

		SpatialGrid(Scalar cell_size) : cell_size_(cell_size) {}

		GridKey get_key(const Vec3& v) const {
			return {
				(int)std::floor(v[0] / cell_size_),
				(int)std::floor(v[1] / cell_size_),
				(int)std::floor(v[2] / cell_size_)
			};
		}

		void insert(const Vec3& pos, uint32 idx) {
			grid_[get_key(pos)].push_back(idx);
		}

		template<typename PositionAttribute>
		bool is_valid_sample(const Vec3& pos, Scalar radius, const PositionAttribute& positions) const {
			GridKey k = get_key(pos);
			Scalar r2 = radius * radius;
			// Check 3x3x3 neighborhood
			for (int dx = -1; dx <= 1; ++dx) {
				for (int dy = -1; dy <= 1; ++dy) {
					for (int dz = -1; dz <= 1; ++dz) {
						GridKey neighbor_key = { k.x + dx, k.y + dy, k.z + dz };
						auto it = grid_.find(neighbor_key);
						if (it != grid_.end()) {
							for (uint32 idx : it->second) {
								if ((positions[idx] - pos).squaredNorm() < r2)
									return false;
							}
						}
					}
				}
			}
			return true;
		}
	};

	struct PointsParameters
	{
		bool initialized_ = false;
		bool fitting_data_computed_ = false;
		// Input Points
		POINTS* points_ = nullptr;
		std::shared_ptr<PAttribute<Vec3>> position_ = nullptr;
		std::shared_ptr<PAttribute<Vec3>> normal_ = nullptr; // Computed on input for sampling
		std::shared_ptr<PAttribute<std::vector<PVertex>>> knn_ = nullptr; // For input normals

		// Sampling & Fitting Data
		POINTS* samples_mesh_ = nullptr;
		std::shared_ptr<PAttribute<Vec3>> samples_position_ = nullptr;
		std::shared_ptr<PAttribute<Vec3>> samples_normal_ = nullptr;
		std::shared_ptr<PAttribute<Scalar>> samples_area_ = nullptr;
		std::shared_ptr<PAttribute<std::vector<PVertex>>> samples_knn_ = nullptr;
		std::shared_ptr<PAttribute<Spherical_Quadric>> samples_quadric_ = nullptr;

		// Medial Axis (on samples)
		std::shared_ptr<PAttribute<Vec3>> samples_ma_position_ = nullptr;
		std::shared_ptr<PAttribute<Scalar>> samples_ma_radius_ = nullptr;
		std::shared_ptr<PAttribute<PVertex>> samples_ma_secondary_vertex_ = nullptr;

		// Clustering (on samples)
		std::shared_ptr<PAttribute<PVertex>> samples_sphere_ = nullptr; // Cluster ID for each sample
		std::shared_ptr<PAttribute<Scalar>> samples_error_ = nullptr;
		std::shared_ptr<PAttribute<Vec4>> samples_color_ = nullptr;
		std::shared_ptr<PAttribute<Vec4>> samples_normal_color_ = nullptr;

		acc::KDTree<3, uint32>* samples_kdtree_ = nullptr; // KDTree of samples
		std::vector<PVertex> samples_kdtree_vertices_; // Vertices of samples in KDTree order

		acc::KDTree<3, uint32>* input_kdtree_ = nullptr; // KDTree of input points
		std::vector<PVertex> input_kdtree_vertices_; // Vertices of input points in KDTree order

		// Spheres (Fitted)
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

		// Skeleton
		NONMANIFOLD* skeleton_ = nullptr;
		std::shared_ptr<NMAttribute<Vec3>> skeleton_position_ = nullptr;
		std::shared_ptr<NMAttribute<Scalar>> skeleton_radius_ = nullptr;

		float32 filter_radius_threshold_ = 0.0f;
		float32 init_dilation_factor_ = 3.0f;
		bool use_sqem_term_ = true;
		bool sphere_correction_ = true;
		CorrectionMode sphere_correction_mode_ = CORRECT_ALWAYS;
		bool auto_stop_ = false;
		bool auto_split_ = false;
		AutoSplitMode auto_split_mode_ = ERROR_THRESHOLD;
		float32 auto_split_error_threshold_ = 0.00025f;
		uint32 auto_split_max_nb_spheres_ = 50;
		bool error_as_spheres_color_ = false;
		float32 spheres_transparency_ = 0.5f;
		float32 sqem_update_lambda_ = 0.02f;
		float32 sqem_clustering_lambda_ = 0.02f;
		
		
		// Filtering
		float32 target_radius_ = 0.1f;
		float32 radius_tolerance_ = 0.01f;

		// Sampling Parameters
		float epsilon_ = 0.01f;
		float sample_radius_ = 0.0025f;
		int sample_iterations_ = 30;  // Max attempts per point
		int knn_k_ = 10;
		int seed_ = 42;
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
			if (samples_kdtree_) delete samples_kdtree_;
			if (input_kdtree_) delete input_kdtree_;
		}
	};

public:
	UDFTraining(const App& app)
		: ViewModule(app, "UDFTraining")
	{
	}

	~UDFTraining()
	{
	}

	void set_selected_surface(SURFACE& s) {} // Compatibility

	void set_selected_points(POINTS& p)
	{
		selected_points_ = &p;
		init_points_data(p);
	}

	void set_point_mesh_provider(MeshProvider<POINTS>* mp) { points_provider_ = mp; }
	void set_point_cloud_render(PointCloudRender<POINTS>* pcr) { pcr_ = pcr; }
	void set_non_manifold_mesh_provider(MeshProvider<NONMANIFOLD>* mp) { non_manifold_provider_ = mp; }
	void set_non_manifold_render(SurfaceRender<NONMANIFOLD>* sr) { sr_nm_ = sr; }

public:
	// --- Surface Sampling (CGAL) ---
	void sample_surface_to_points(SURFACE& surface, POINTS& points, int num_samples)
	{
		points_provider_->clear_mesh(points);

		// 1. Convert CGoGN SURFACE to CGAL::Surface_mesh
		CGAL_Mesh cgal_mesh;
		
		auto pos = cgogn::get_attribute<Vec3, SVertex>(surface, "position");

		std::unordered_map<uint32, CGAL_Mesh::Vertex_index> v_map;

		// Add vertices
		foreach_cell(surface, [&](SVertex v) {
			uint32 v_idx = index_of(surface, v);
			const Vec3& p = (*pos)[v_idx];
			v_map[v_idx] = cgal_mesh.add_vertex(Point_3(p[0], p[1], p[2]));
			return true;
		});

		// Add faces
		foreach_cell(surface, [&](SFace f) {
			std::vector<CGAL_Mesh::Vertex_index> face_v;
			foreach_incident_vertex(surface, f, [&](SVertex v) {
				face_v.push_back(v_map[index_of(surface, v)]);
				return true;
			});
				
			cgal_mesh.add_face(face_v);
			return true;
		});
	

		// 2. Sample mesh
		std::vector<Point_3> sampled_points;

		// Using simple random sampling on mesh
		CGAL::Polygon_mesh_processing::sample_triangle_mesh(cgal_mesh, std::back_inserter(sampled_points),
															CGAL::parameters::use_grid_sampling(true).grid_spacing(0.01));

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
		timer_connection_ = boost::synapse::connect<App::timer_tick>(&app_, [this]() {
			if (selected_points_)
			{
				PointsParameters& p = points_parameters_[selected_points_];
				update_render_data(p);
			}
		});
	}


private:
	// --- Initialization ---

	void init_points_data(POINTS& m)
	{
		PointsParameters& p = points_parameters_[&m];
		if (p.initialized_) return;

		// Init Input Points 
		p.points_ = &m;
		p.position_ = get_attribute<Vec3, PVertex>(m, "position");
		p.normal_ = get_or_add_attribute<Vec3, PVertex>(m, "normal");
		p.knn_ = get_or_add_attribute<std::vector<PVertex>, PVertex>(m, "knn");

		// Init Samples Mesh
		std::string sample_name = points_provider_->mesh_name(*p.points_) + "_samples";
		if (!p.samples_mesh_)
			p.samples_mesh_ =
				points_provider_->has_mesh(sample_name) ? points_provider_->mesh(sample_name) : points_provider_->add_mesh(sample_name);
		else points_provider_->clear_mesh(*p.samples_mesh_);

		// Initialize attributes for samples
		p.samples_position_ = get_or_add_attribute<Vec3, PVertex>(*p.samples_mesh_, "position");
		p.samples_normal_ = get_or_add_attribute<Vec3, PVertex>(*p.samples_mesh_, "normal");
		p.samples_area_ = get_or_add_attribute<Scalar, PVertex>(*p.samples_mesh_, "area");
		p.samples_knn_ = get_or_add_attribute<std::vector<PVertex>, PVertex>(*p.samples_mesh_, "knn");
		p.samples_quadric_ = get_or_add_attribute<Spherical_Quadric, PVertex>(*p.samples_mesh_, "quadric");
		p.samples_ma_position_ = get_or_add_attribute<Vec3, PVertex>(*p.samples_mesh_, "ma_position");
		p.samples_ma_radius_ = get_or_add_attribute<Scalar, PVertex>(*p.samples_mesh_, "ma_radius");
		p.samples_ma_secondary_vertex_ = get_or_add_attribute<PVertex, PVertex>(*p.samples_mesh_, "ma_secondary_vertex");
		p.samples_sphere_ = get_or_add_attribute<PVertex, PVertex>(*p.samples_mesh_, "sphere");
		p.samples_error_ = get_or_add_attribute<Scalar, PVertex>(*p.samples_mesh_, "error");
		p.samples_color_ = get_or_add_attribute<Vec4, PVertex>(*p.samples_mesh_, "color");
		p.samples_normal_color_ = get_or_add_attribute<Vec4, PVertex>(*p.samples_mesh_, "normal_color");

		// Init Spheres Mesh
		std::string sphere_name = points_provider_->mesh_name(m) + "_spheres";
		if (!p.spheres_)
			p.spheres_ =
				points_provider_->has_mesh(sphere_name) ? points_provider_->mesh(sphere_name) : points_provider_->add_mesh(sphere_name);
		
		p.spheres_position_ = get_or_add_attribute<Vec3, PVertex>(*p.spheres_, "position");
		p.spheres_radius_ = get_or_add_attribute<Scalar, PVertex>(*p.spheres_, "radius");
		p.spheres_color_ = get_or_add_attribute<Vec4, PVertex>(*p.spheres_, "color");
		p.spheres_cluster_ = get_or_add_attribute<std::vector<PVertex>, PVertex>(*p.spheres_, "cluster");
		p.spheres_cluster_color_ = get_or_add_attribute<Vec4, PVertex>(*p.spheres_, "cluster_color");
		p.spheres_cluster_area_ = get_or_add_attribute<Scalar, PVertex>(*p.spheres_, "cluster_area");
		p.spheres_neighbor_clusters_ = get_or_add_attribute<std::set<PVertex>, PVertex>(*p.spheres_, "neighbor_clusters");
		p.spheres_parent_ = get_or_add_attribute<PVertex, PVertex>(*p.spheres_, "parent");
		p.spheres_do_not_split_ = get_or_add_attribute<bool, PVertex>(*p.spheres_, "do_not_split");
		p.spheres_error_ = get_or_add_attribute<Scalar, PVertex>(*p.spheres_, "error");
		p.spheres_error_not_normalized_ = get_or_add_attribute<Scalar, PVertex>(*p.spheres_, "error_not_normalized");

		// Init Skeleton Mesh
		std::string skel_name = points_provider_->mesh_name(m) + "_skeleton";
		if (!p.skeleton_)
			p.skeleton_ = non_manifold_provider_->has_mesh(skel_name) ? non_manifold_provider_->mesh(skel_name) : non_manifold_provider_->add_mesh(skel_name);
		p.skeleton_position_ = get_or_add_attribute<Vec3, NMVertex>(*p.skeleton_, "position");
		p.skeleton_radius_ = get_or_add_attribute<Scalar, NMVertex>(*p.skeleton_, "radius");

		p.initialized_ = true;
	}

	void compute_fitting_data(PointsParameters& p)
	{
		if (p.fitting_data_computed_) return;
		if (!p.samples_mesh_) {
			std::cerr << "Error: No sampled points found. Please sample points first." << std::endl;
			return;
		}

		std::cout << "Building KDTree..." << std::endl;
		build_kdtree(p);
		std::cout << "Computing KNN and Area..." << std::endl;
		compute_samples_area(p); // Compute KNN and Area for samples
		std::cout << "Computing Quadrics..." << std::endl;
		compute_quadrics(p);
		std::cout << "Computing Initial Medial Axis..." << std::endl;
		compute_initial_medial_axis(p);

		std::cout << "Fitting Data Computed." << std::endl;

		p.fitting_data_computed_ = true;
	}

	void build_kdtree(PointsParameters& p)
	{
		if (p.samples_kdtree_) delete p.samples_kdtree_;

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
	template<typename MESH, typename POS_ATTR>
	Vec3 compute_pca_normal(const MESH& mesh, const POS_ATTR& position, 
							const std::vector<uint32>& indices, 
							const std::vector<typename mesh_traits<MESH>::Vertex>& kdtree_vertices)
	{
		using Vertex = typename mesh_traits<MESH>::Vertex;
		if (indices.size() < 3) return Vec3(0, 0, 1);
		
		Vec3 centroid(0, 0, 0);
		for (uint32 idx : indices) {
			Vertex v = kdtree_vertices[idx];
			uint32 v_idx = index_of(mesh, v);
			centroid += position[v_idx];
		}
		centroid /= Scalar(indices.size());
		
		Eigen::Matrix3d covariance = Eigen::Matrix3d::Zero();
		for (uint32 idx : indices) {
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
			for(auto& res : knn_res) indices.push_back(res.first);
			
			(*p.normal_)[v_idx] = compute_pca_normal(*p.points_, *p.position_, indices, p.input_kdtree_vertices_);
			return true;
		});
	}

	Vec3 compute_avg_normal(PointsParameters& p, PVertex v, acc::KDTree<3, uint32>& kdtree, const std::vector<PVertex>& vertices)
	{
		uint32 v_idx = index_of(*p.points_, v);
		const Vec3& pt = (*p.position_)[v_idx];
		const Vec3& n = (*p.normal_)[v_idx];
		
		std::vector<std::pair<uint32, Scalar>> knn_res;
		kdtree.find_nns(pt, p.knn_k_, &knn_res);

		Vec3 avg_n(0,0,0);
		for(auto& res : knn_res) {
			PVertex neighbor = vertices[res.first];
			uint32 n_idx = index_of(*p.points_, neighbor);
			Vec3 nn = (*p.normal_)[n_idx];
			if (nn.dot(n) < 0) nn = -nn;
			avg_n += nn;
		}
		return avg_n.normalized();
	}
	void compute_samples_area(PointsParameters& p)
	{
		// Compute KNN and Area on samples_mesh_
		if (!p.samples_mesh_ || !p.samples_kdtree_) return;

		parallel_foreach_cell(*p.samples_mesh_, [&](PVertex v) {
			uint32 v_idx = index_of(*p.samples_mesh_, v);
			const Vec3& pt = (*p.samples_position_)[v_idx];
			std::vector<std::pair<uint32, Scalar>> knn_res;
			p.samples_kdtree_->find_nns(pt, p.knn_k_, &knn_res);
			
			(*p.samples_knn_)[v_idx].clear();
			Scalar sum_dist = 0.0;
			for(auto& res : knn_res) {
				if (p.samples_kdtree_vertices_[res.first] != v) {
					(*p.samples_knn_)[v_idx].push_back(p.samples_kdtree_vertices_[res.first]);
					sum_dist += res.second;
				}
			}
			// Normals are already computed/oriented in sample_points
			(*p.samples_area_)[v_idx] = (sum_dist * sum_dist) / (2.0 * p.knn_k_); // Rough area estimate
			return true;
		});
	}

	void compute_quadrics(PointsParameters& p)
	{	
		parallel_foreach_cell(*p.samples_mesh_, [&](PVertex v) {
			uint32 v_idx = index_of(*p.samples_mesh_, v);
			Spherical_Quadric& q = (*p.samples_quadric_)[v_idx];
			q.clear();
			const Vec3& pos = (*p.samples_position_)[v_idx];
			const Vec3& n = (*p.samples_normal_)[v_idx];
			Scalar a = (*p.samples_area_)[v_idx] / (p.knn_k_ + 1.0);
			q += Spherical_Quadric(Vec4(pos.x(), pos.y(), pos.z(), 0), Vec4(n.x(), n.y(), n.z(), 1)) * a;
			for (PVertex vn : (*p.samples_knn_)[v_idx]) {
				uint32 vn_idx = index_of(*p.samples_mesh_, vn);
				const Vec3& pn = (*p.samples_position_)[vn_idx];
				const Vec3& nn = (*p.samples_normal_)[vn_idx];
				Scalar an = (*p.samples_area_)[vn_idx] / (p.knn_k_ + 1.0);
				q += Spherical_Quadric(Vec4(pn.x(), pn.y(), pn.z(), 0), Vec4(nn.x(), nn.y(), nn.z(), 1)) * an;
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

	std::pair<Vec3, Vec3> project_and_normal(PointsParameters& p, const Vec3& sample_pos)
	{
		std::pair<uint32, Scalar> knn_res;
		p.input_kdtree_->find_nn(sample_pos, &knn_res);
		PVertex nn = p.input_kdtree_vertices_[knn_res.first];
		Vec3 nn_pos = (*p.position_)[index_of(*p.points_, nn)];
		Vec3 n = (sample_pos - nn_pos).normalized();
		Vec3 query = nn_pos + n * p.epsilon_;
		PVertex last_nn = PVertex();
		do {
			p.input_kdtree_->find_nn(query, &knn_res);
			nn = p.input_kdtree_vertices_[knn_res.first];
			Vec3 nn_pos = (*p.position_)[index_of(*p.points_, nn)];
			Vec3 n = (query - nn_pos).normalized();
			query = nn_pos + n * p.epsilon_;
			last_nn = nn;
		}while(index_of(*p.points_, nn) != index_of(*p.points_, last_nn));
		
		return {query, n};
	}

	void sample_points(PointsParameters& p)
	{
		// Build KDTree for input points
		if (p.input_kdtree_) delete p.input_kdtree_;
		std::vector<Vec3> points;
		p.input_kdtree_vertices_.clear();
		points.reserve(nb_cells<PVertex>(*p.points_));
		p.input_kdtree_vertices_.reserve(points.size());
		foreach_cell(*p.points_, [&](PVertex v) {
			uint32 v_idx = index_of(*p.points_, v);
			points.push_back((*p.position_)[v_idx]);
			p.input_kdtree_vertices_.push_back(v);
			return true;
		});
		p.input_kdtree_ = new acc::KDTree<3, uint32>(points);

		std::uniform_real_distribution<Scalar> uniform(0.0, 1.0);
		std::mt19937 gen(p.seed_);
		
		// Pick a random point from input as generator seed
		std::uniform_int_distribution<uint32> uniform_idx(0, uint32(p.input_kdtree_vertices_.size() - 1));
		uint32 rand_start_idx = uniform_idx(gen);
		PVertex start_seed_vertex = p.input_kdtree_vertices_[rand_start_idx];
		Vec3 generator = (*p.position_)[index_of(*p.points_, start_seed_vertex)];

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
					if (input_nn_res.second < p.epsilon_)
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
						std::cout << "Sampled point " << count  << "\r" << std::flush;
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
			(*p.samples_normal_color_)[v_idx] = Vec4((n.x() + 1.0) * 0.5, (n.y() + 1.0) * 0.5, (n.z() + 1.0) * 0.5, 1.0);
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
			
			auto [c1, r1, q1] = cgogn::geometry::shrinking_ball_center<PVertex>(
				pt, n,
				p.samples_kdtree_, p.samples_kdtree_vertices_,
				true
			);
			
			(*p.samples_ma_position_)[v_idx] = c1;
			(*p.samples_ma_radius_)[v_idx] = r1;
			(*p.samples_ma_secondary_vertex_)[v_idx] = *reinterpret_cast<PVertex*>(&q1);
	
			return true;
		});

		// Filter points with ball radius >= 1.1 * epsilon (ball may be outside the surface)
		// First try to fix by re-estimating normal from input cloud
		std::vector<PVertex> to_remove;
		foreach_cell(*p.samples_mesh_, [&](PVertex v) {
			uint32 v_idx = index_of(*p.samples_mesh_, v);
			Scalar radius = (*p.samples_ma_radius_)[v_idx];
			
			if (radius >= 1.1 * p.epsilon_)
			{
				// Try to fix: re-estimate normal using local PCA on samples
				Vec3 pt = (*p.samples_position_)[v_idx];
				
				// Find K nearest neighbors in samples mesh
				std::vector<std::pair<uint32, Scalar>> knn_res;
				p.samples_kdtree_->find_nns(pt, p.knn_k_, &knn_res);
				
				// Extract indices for PCA
				std::vector<uint32> indices;
				for (auto& res : knn_res)
					indices.push_back(res.first);
				
				// Compute PCA normal from local neighborhood
				if (indices.size() >= 3)
				{
					Vec3 new_normal = compute_pca_normal(*p.samples_mesh_, *p.samples_position_, indices, p.samples_kdtree_vertices_);
					
					// Orient normal: should point away from surface (same side as current position relative to input cloud)
					std::pair<uint32, Scalar> input_nn_res;
					p.input_kdtree_->find_nn(pt, &input_nn_res);
					PVertex input_nn = p.input_kdtree_vertices_[input_nn_res.first];
					Vec3 input_nn_pos = (*p.position_)[index_of(*p.points_, input_nn)];
					Vec3 to_pt = (pt - input_nn_pos).normalized();
					if (new_normal.dot(to_pt) < 0)
						new_normal = -new_normal;
					
					(*p.samples_normal_)[v_idx] = new_normal;
					
					// Re-compute shrinking ball with new normal
					auto [c2, r2, q2] = cgogn::geometry::shrinking_ball_center<PVertex>(
						pt, new_normal,
						p.samples_kdtree_, p.samples_kdtree_vertices_,
						true
					);
					
					(*p.samples_ma_position_)[v_idx] = c2;
					(*p.samples_ma_radius_)[v_idx] = r2;
					(*p.samples_ma_secondary_vertex_)[v_idx] = *reinterpret_cast<PVertex*>(&q2);
					
					// Still bad? Mark for removal
					if (r2 >= 2 * p.epsilon_)
						to_remove.push_back(v);
				}
				else
				{
					// Not enough neighbors for PCA, mark for removal
					to_remove.push_back(v);
				}
			}
			return true;
		});
		if (!to_remove.empty())
		{
			std::cout << "Removed " << to_remove.size() << " points with ball outside surface (after retry)."
					  << std::endl;
			
		}
		for (PVertex v : to_remove)
			remove_vertex(*p.samples_mesh_, v);
		build_kdtree(p); // Rebuild KDTree after removing vertices
		points_provider_->emit_connectivity_changed(*p.samples_mesh_);
		
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
				for(PVertex u : (*p.samples_knn_)[w_idx]) {
					uint32 u_idx = index_of(*p.samples_mesh_, u);
					if (!(*covered)[u_idx] &&
						((*p.samples_position_)[u_idx] - vp).norm() < p.init_dilation_factor_ * vr)
						stack.push_back(u);
				}
			}
			
			// Also check secondary vertex logic if needed, but for now stick to KNN propagation
			PVertex secondary = (*p.samples_ma_secondary_vertex_)[v_index];
			if (secondary.is_valid()) {
				stack.push_back(secondary);
				while (!stack.empty())
				{
					PVertex w = stack.back();
					stack.pop_back();
					uint32 w_idx = index_of(*p.samples_mesh_, w);
					(*covered)[w_idx] = true;
					
					for(PVertex u : (*p.samples_knn_)[w_idx]) {
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
				Scalar dist_eucl = ((vp - center).norm() - radius);
				dist_eucl *= dist_eucl;
				dist_eucl *= a;
				Scalar dist;
				if (p.use_sqem_term_)
				{
					Scalar dist_sqem =
						(*p.samples_quadric_)[v_index].eval(Vec4(center.x(), center.y(), center.z(), radius));
					dist = dist_sqem + p.sqem_clustering_lambda_ * dist_eucl;
				}
				else
					dist = dist_eucl;
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
		// remove small clusters
		foreach_cell(*p.spheres_, [&](PVertex v) -> bool {
			uint32 v_idx = index_of(*p.spheres_, v);
			std::vector<PVertex>& cluster = (*p.spheres_cluster_)[v_idx];
			if (cluster.size() < 4)
			{
				for (PVertex sv : cluster) {
					uint32 sv_idx = index_of(*p.samples_mesh_, sv);
					(*p.samples_sphere_)[sv_idx] = PVertex();
				}
				remove_vertex(*p.spheres_, v);
				p.nb_spheres_--;
			}
			return true;
		});
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
		
		// remove small clusters
		foreach_cell(*p.spheres_, [&](PVertex v) -> bool {
			uint32 v_idx = index_of(*p.spheres_, v);
			std::vector<PVertex>& cluster = (*p.spheres_cluster_)[v_idx];
			if (cluster.size() < 4)
			{
				for (PVertex sv : cluster) {
					uint32 sv_idx = index_of(*p.samples_mesh_, sv);
					(*p.samples_sphere_)[sv_idx] = PVertex();
				}
				remove_vertex(*p.spheres_, v);
				p.nb_spheres_--;
			}
			return true;
		});
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

				Scalar a = (*p.samples_area_)[sv_index];

				Scalar dist_eucl = ((*p.samples_position_)[sv_index] - center).norm() - radius;
				dist_eucl *= dist_eucl;
				dist_eucl *= a;
				Scalar dist;
				if (p.use_sqem_term_)
				{
					Scalar dist_sqem =
						(*p.samples_quadric_)[sv_index].eval(Vec4(center.x(), center.y(), center.z(), radius));
					dist = dist_sqem + p.sqem_clustering_lambda_ * dist_eucl;
				}
				else
					dist = dist_eucl;
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
			
			if (error < p.min_error_) p.min_error_ = error;
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
		if (cluster.empty()) return;

		Vec3 c = (*p.spheres_position_)[sphere_index];
		Scalar r = (*p.spheres_radius_)[sphere_index];

		Eigen::MatrixXd J(cluster.size(), 4);
		J.setZero();
		Eigen::VectorXd b(cluster.size());
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

				Vec3 d = pos - Vec3(s(0), s(1), s(2));
				Scalar l = d.norm();
				
				Scalar a = sqrt((*p.samples_area_)[v_index]);
				J.row(idx) = Eigen::Vector4d(-(d[0] / l), -(d[1] / l), -(d[2] / l), -1.0) * a;
				b(idx) = -(l - s(3)) * a;
				
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

		(*p.spheres_position_)[sphere_index] = c;
		(*p.spheres_radius_)[sphere_index] = r;
	}

	void update_sphere_sqem(PointsParameters& p, PVertex sphere)
	{
		uint32 sphere_index = index_of(*p.spheres_, sphere);

		const std::vector<PVertex>& cluster = (*p.spheres_cluster_)[sphere_index];
		if (cluster.empty()) return;

		Vec3 c = (*p.spheres_position_)[sphere_index];
		Scalar r = (*p.spheres_radius_)[sphere_index];

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

		(*p.spheres_position_)[sphere_index] = c;
		(*p.spheres_radius_)[sphere_index] = r;
	}

	void correct_sphere(PointsParameters& p, PVertex v)
	{
		uint32 v_index = index_of(*p.spheres_, v);

		Vec3& c = (*p.spheres_position_)[v_index];
		Scalar& r = (*p.spheres_radius_)[v_index];

		std::pair<uint32, Scalar> k_res;
		p.samples_kdtree_->find_nn(c, &k_res);
		uint32 k_idx = index_of(*p.samples_mesh_, p.samples_kdtree_vertices_[k_res.first]);
		Vec3 closest_pos = (*p.samples_position_)[k_idx];
		Vec3 dir = (closest_pos - c).normalized();
		
		auto [nc, nr, nq] = cgogn::geometry::shrinking_ball_center<PVertex>(
			closest_pos, dir,
			p.samples_kdtree_, p.samples_kdtree_vertices_,
			true
		);
		c = nc; r = nr;
	}

	void update_spheres(PointsParameters& p)
	{
		compute_clusters(p);

		parallel_foreach_cell(*p.spheres_, [&](PVertex v) -> bool {
			if (p.use_sqem_term_)
				update_sphere_sqem(p, v);
			else
				update_sphere_euclidean(p, v);
			
			if (p.sphere_correction_ && p.sphere_correction_mode_ == CORRECT_ALWAYS)
				correct_sphere(p, v);
			
			(*p.spheres_do_not_split_)[index_of(*p.spheres_, v)] = false;
			return true;
		});

		compute_spheres_error(p);

		if (p.auto_split_ &&
			(p.total_error_diff_ < 1e-5 || p.iteration_count_ % 10 == 0))
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
					
					uint32 to_split_max = std::min(uint32(std::ceil(p.nb_spheres_ * 0.2)), 10u);
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

		auto spheres_skeleton_vertex_map =
			add_attribute<NMVertex, PVertex>(*p.spheres_, "__spheres_skeleton_vertex_map");

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
			uint32 pv_index = index_of(*p.spheres_, pv);
			NMVertex nmv1 = (*spheres_skeleton_vertex_map)[pv_index];
			const std::set<PVertex>& n_pv = (*p.spheres_neighbor_clusters_)[pv_index];
			for (const PVertex& ne1 : n_pv)
			{
				uint32 ne1_index = index_of(*p.spheres_, ne1);
				NMVertex nmv2 = (*spheres_skeleton_vertex_map)[ne1_index];
				const std::set<PVertex>& ne_ne1 = (*p.spheres_neighbor_clusters_)[ne1_index];
				for (const PVertex& ne2 : ne_ne1)
				{
					if (std::find(n_pv.begin(), n_pv.end(), ne2) != n_pv.end())
					{
						uint32 ne2_index = index_of(*p.spheres_, ne2);
						NMVertex nmv3 = (*spheres_skeleton_vertex_map)[ne2_index];
						if (index_of(*p.skeleton_, nmv1) < index_of(*p.skeleton_, nmv2) &&
							index_of(*p.skeleton_, nmv2) < index_of(*p.skeleton_, nmv3))
						{
							std::vector<NMEdge> edges;
							edges.reserve(3);
							edges.push_back(edge_indices[{index_of(*p.skeleton_, nmv1), index_of(*p.skeleton_, nmv2)}]);
							edges.push_back(edge_indices[{index_of(*p.skeleton_, nmv2), index_of(*p.skeleton_, nmv3)}]);
							edges.push_back(edge_indices[{index_of(*p.skeleton_, nmv1), index_of(*p.skeleton_, nmv3)}]);
							add_face(*p.skeleton_, edges);
						}
					}
				}
			}
			return true;
		});

		remove_attribute<PVertex>(*p.spheres_, spheres_skeleton_vertex_map);
	}

protected:
	void split_sphere(PointsParameters& p, PVertex sphere)
	{
		if (!sphere.is_valid()) return;
		uint32 s_index = index_of(*p.spheres_, sphere);
		Vec3 c = (*p.spheres_position_)[s_index];
		Scalar r = (*p.spheres_radius_)[s_index];

		// find the point in the cluster with the max error
		const std::vector<PVertex>& cluster = (*p.spheres_cluster_)[s_index];
		if (cluster.empty()) return;

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

		if (!max_err_v.is_valid()) return;

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
		if (!selected_points_) return;
		PointsParameters& p = points_parameters_[selected_points_];

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
		if (ImGui::BeginCombo("Input Point Cloud", selected_points_ ? points_provider_->mesh_name(*selected_points_).c_str() : "None"))
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

		if (!selected_points_) return;
		PointsParameters& p = points_parameters_[selected_points_];

		ImGui::Separator();

		// Sampling
		if (ImGui::CollapsingHeader("Sampling", ImGuiTreeNodeFlags_DefaultOpen))
		{
			ImGui::InputFloat("Epsilon", &p.epsilon_, 0.001f, 0.1f, "%.4f");
			ImGui::InputFloat("Sample Radius", &p.sample_radius_, 0.0001f, 0.01f, "%.4f");
			ImGui::SliderInt("Sample Iterations", &p.sample_iterations_, 10, 100);
			ImGui::SliderInt("KNN for Normal", &p.knn_k_, 3, 50);

			if (ImGui::Button("Sample Points"))
				sample_points(p);
			ImGui::SameLine();
			if (ImGui::Button("Clear Samples"))
			{
				if (p.samples_mesh_) points_provider_->clear_mesh(*p.samples_mesh_);
				p.fitting_data_computed_ = false;
			}
		}

		bool has_samples = p.samples_mesh_ && nb_cells<PVertex>(*p.samples_mesh_) > 0;

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

					static bool sync_lambda = false;
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

					ImGui::Checkbox("Use SQEM term", &p.use_sqem_term_);

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

					ImGui::Text("Total error: %f", p.total_error_);
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

	MeshProvider<POINTS>* points_provider_ = nullptr;
	PointCloudRender<POINTS>* pcr_ = nullptr;
	MeshProvider<NONMANIFOLD>* non_manifold_provider_ = nullptr;
	SurfaceRender<NONMANIFOLD>* sr_nm_ = nullptr;
	
	POINTS* selected_points_ = nullptr;
	std::map<POINTS*, PointsParameters> points_parameters_;
	PVertex picked_sphere_;
	std::shared_ptr<boost::synapse::connection> timer_connection_;
	std::array<std::mutex, 43> spheres_mutex_;
};

} // namespace ui

} // namespace cgogn

#endif // CGOGN_MODULE_UDF_TRAINING_H_
