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
#ifndef CGOGN_GEOMETRY_ALGOS_OFFSET_AXIS_DISTANCE_FIELD_H_
#define CGOGN_GEOMETRY_ALGOS_OFFSET_AXIS_DISTANCE_FIELD_H_

#include <cgogn/core/functions/attributes.h>
#include <cgogn/core/functions/mesh_info.h>
#include <cgogn/geometry/algos/offset_axis/ray_level_set_sampler_traits.h>
#include <cgogn/geometry/types/vector_traits.h>

#include <memory>
#include <vector>

namespace cgogn::geometry
{

// Distance field: the function the UDF reconstruction samples and scores. The reconstruction is written once and
// the input kind selects a field class through DistanceFieldFor<RaySamplerTag, SURFACE, POINTS>::type:
//   - SurfaceDistanceField: the unsigned distance to a triangle surface (exact closest point on the BVH);
//   - PointCloudDistanceField: the unsigned distance to an oriented point cloud (k-d tree queries);
//   - NeuralDistanceField: a TorchScript network that predicts the UDF or the MF value of the input.
// Every field is constructed from the InputGeometry caches and the ReconstructionOptions of the reconstruction, and
// provides the following interface (the CPU fields hold no model, see the notes on each method):
//   RayConfig, RayTraits, RaySampler, is_neural_type  types of the ray sampler and compile-time kind of the field
//   prepare(inputs)                              build the caches, or report that the field is not usable yet
//   sampling_box(inputs, min, max)               bounding box of the input, before padding
//   ray_traits(inputs)                           traits of the ray sampler (valid until the sampling ends)
//   make_ray_sampler(params)                     new ray sampler for the field (on the device of the field)
//   configure_ray_sampler(sampler, params)       update a ray sampler created by make_ray_sampler
//   normals(points, normals, keep)               normals of the candidate points, keep flags the usable ones
//   project(points, parameters, result)          projection of the points onto the alpha level set
//   load_model(path, type) / reset() / loaded()  neural model life cycle (CPU fields: no model, loaded() is false)
//   distances(points, values)                    topology scores of the points (unsigned distance to the input)
//   values_sdf(points, values, sdf)              MF value and sdf of the points (CPU fields: returns false)
//   recompute_normals_from_input_after_filtering constant: whether the normals are recomputed from the input after
//                                                the sampling filter
//   forward() / workspace() / device()           neural state (neural fields only)

// Result of the preparation and sampling-box steps of a field.
enum class FieldPrepareResult
{
	ok,
	invalid_data,
	model_not_loaded
};

// Reconstruction data a field may need. The pointers are non-owning: the owner is ReconstructionState::data and the
// attributes are the ones created by UDFReconstruction::initialize_data().
template <typename SURFACE, typename POINTS>
struct PrepareInputs
{
	using PVertex = typename mesh_traits<POINTS>::Vertex;
	template <typename T>
	using PointAttribute = typename mesh_traits<POINTS>::template Attribute<T>;

	SURFACE* surface = nullptr;
	POINTS* input_points = nullptr;
	std::shared_ptr<PointAttribute<Vec3>> input_position;
	std::shared_ptr<PointAttribute<Vec3>> input_normal;
	std::shared_ptr<PointAttribute<std::vector<PVertex>>> input_knn;
};

template <typename SURFACE, typename POINTS>
class SurfaceDistanceField;

template <typename SURFACE, typename POINTS>
class PointCloudDistanceField;

// Maps a ray sampler tag to the distance field class that implements it.
template <typename RaySamplerTag, typename SURFACE, typename POINTS>
struct DistanceFieldFor;

template <typename SURFACE, typename POINTS>
struct DistanceFieldFor<RaySamplerSurface, SURFACE, POINTS>
{
	using type = SurfaceDistanceField<SURFACE, POINTS>;
};

template <typename SURFACE, typename POINTS>
struct DistanceFieldFor<RaySamplerPointCloud, SURFACE, POINTS>
{
	using type = PointCloudDistanceField<SURFACE, POINTS>;
};

// the neural specialization is in neural/neural_distance_field.h; include it to reconstruct from a neural field

} // namespace cgogn::geometry

#endif // CGOGN_GEOMETRY_ALGOS_OFFSET_AXIS_DISTANCE_FIELD_H_
