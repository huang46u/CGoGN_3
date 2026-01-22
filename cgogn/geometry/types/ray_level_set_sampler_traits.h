#ifndef CGOGN_GEOMETRY_TYPES_RAY_LEVEL_SET_SAMPLER_TRAITS_H_
#define CGOGN_GEOMETRY_TYPES_RAY_LEVEL_SET_SAMPLER_TRAITS_H_

#include <cgogn/core/functions/mesh_info.h>
#include <cgogn/geometry/types/neural_field_forward.h>
#include <cgogn/geometry/types/vector_traits.h>

#include <libacc/bvh_tree.h>
#include <libacc/kd_tree.h>

#include <cmath>
#include <limits>
#include <utility>
#include <vector>

namespace cgogn
{

namespace geometry
{

struct RaySamplerNeural
{
};

struct RaySamplerSurface
{
};

struct RaySamplerPointCloud
{
};

struct NeuralFieldRayTraits
{
	static constexpr bool kUsesTorch = true;
	NeuralFieldForward* field = nullptr;

	bool is_ready() const
	{
		return field && field->is_loaded();
	}

	const torch::Device& device() const
	{
		return field->device();
	}

	std::pair<torch::Tensor, torch::Tensor> eval_values_sdf(const torch::Tensor& points) const
	{
		if (!field)
			return {torch::Tensor(), torch::Tensor()};
		return field->forward_values_sdf_gpu(points);
	}
};

template <typename RaySamplerTag, typename SURFACE, typename POINTS>
struct RaySamplerTraits
{
};

template <typename SURFACE, typename POINTS>
struct RaySamplerTraits<RaySamplerNeural, SURFACE, POINTS>
{
	using SamplerTraits = NeuralFieldRayTraits;

	static SamplerTraits make(NeuralFieldForward& field)
	{
		SamplerTraits traits;
		traits.field = &field;
		return traits;
	}
};

template <typename SURFACE>
struct SurfaceMeshRayTraits
{
	static constexpr bool kUsesTorch = false;
	using Face = typename mesh_traits<SURFACE>::Face;
	template <typename T>
	using Attribute = typename mesh_traits<SURFACE>::template Attribute<T>;

	SURFACE* mesh = nullptr;
	Attribute<Vec3>* position = nullptr;
	acc::BVHTree<uint32, Vec3>* bvh = nullptr;
	std::vector<Face>* bvh_faces = nullptr;

	bool is_ready() const
	{
		return bvh != nullptr;
	}

	Scalar eval_distance(const Vec3& pos) const
	{
		if (!bvh)
			return std::numeric_limits<Scalar>::max();
		std::pair<uint32, Vec3> cp;
		if (!bvh->closest_point(pos, &cp))
			return std::numeric_limits<Scalar>::max();
		return (pos - cp.second).norm();
	}
};

template <typename SURFACE, typename POINTS>
struct RaySamplerTraits<RaySamplerSurface, SURFACE, POINTS>
{
	using SamplerTraits = SurfaceMeshRayTraits<SURFACE>;
	template <typename T>
	using Attribute = typename mesh_traits<SURFACE>::template Attribute<T>;
	using Face = typename mesh_traits<SURFACE>::Face;

	static SamplerTraits make(SURFACE& mesh, Attribute<Vec3>* position, acc::BVHTree<uint32, Vec3>* bvh,
							  std::vector<Face>* bvh_faces)
	{
		SamplerTraits traits;
		traits.mesh = &mesh;
		traits.position = position;
		traits.bvh = bvh;
		traits.bvh_faces = bvh_faces;
		return traits;
	}
};

template <typename POINTS>
struct PointCloudRayTraits
{
	static constexpr bool kUsesTorch = false;
	using Vertex = typename mesh_traits<POINTS>::Vertex;
	template <typename T>
	using Attribute = typename mesh_traits<POINTS>::template Attribute<T>;

	POINTS* mesh = nullptr;
	Attribute<Vec3>* position = nullptr;
	Attribute<Vec3>* normal = nullptr;
	Attribute<std::vector<Vertex>>* knn = nullptr;
	acc::KDTree<3, uint32>* kdtree = nullptr;
	const std::vector<Vertex>* kdtree_vertices = nullptr;

	bool is_ready() const
	{
		return mesh && position && normal && knn && kdtree && kdtree_vertices;
	}

	Scalar eval_distance(const Vec3& pos) const
	{
		std::pair<uint32, Scalar> knn_res;
		kdtree->find_nn(pos, &knn_res);
		return knn_res.second;

	}
};

template <typename SURFACE, typename POINTS>
struct RaySamplerTraits<RaySamplerPointCloud, SURFACE, POINTS>
{
	using SamplerTraits = PointCloudRayTraits<POINTS>;
	using Vertex = typename mesh_traits<POINTS>::Vertex;
	template <typename T>
	using Attribute = typename mesh_traits<POINTS>::template Attribute<T>;

	static SamplerTraits make(POINTS& mesh, Attribute<Vec3>* position, Attribute<Vec3>* normal,
							  Attribute<std::vector<Vertex>>* knn, acc::KDTree<3, uint32>* kdtree,
							  const std::vector<Vertex>* kdtree_vertices)
	{
		SamplerTraits traits;
		traits.mesh = &mesh;
		traits.position = position;
		traits.normal = normal;
		traits.knn = knn;
		traits.kdtree = kdtree;
		traits.kdtree_vertices = kdtree_vertices;
		return traits;
	}
};

} // namespace geometry

} // namespace cgogn

#endif // CGOGN_GEOMETRY_TYPES_RAY_LEVEL_SET_SAMPLER_TRAITS_H_
