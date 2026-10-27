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
#include <cgogn/geometry/algos/offset_axis/neural_model_type.h>
#include <cgogn/geometry/types/vector_traits.h>

#include <memory>
#include <vector>

namespace cgogn
{

namespace geometry
{

namespace offset_axis
{

// Distance field: the function the UDF reconstruction samples and scores, the FIELD template parameter of the
// reconstruction (SurfaceDistanceField, PointCloudDistanceField, NeuralDistanceField). Every field is constructed
// from the InputGeometry and the ReconstructionOptions and provides:
//   prepare(inputs)                        build the caches, or report that the field is not usable yet
//   sampling_box(inputs, min, max)         bounding box of the input, before padding
//   ray_traits(inputs)                     traits of the ray sampler (valid until the sampling ends)
//   normals(points, normals, keep)         normals of the candidate points, keep flags the usable ones
//   project(points, parameters, result)    projection of the points onto the alpha level set
//   distances(points, values)              topology scores of the points
//   values_sdf(points, values, sdf)        MF value and sdf of the points (CPU fields: returns false)
//   load_model / reset / loaded            neural model life cycle (CPU fields: no model)

enum class FieldPrepareResult
{
	ok,
	invalid_data,
	model_not_loaded
};

// Non-owning pointers to the reconstruction data a field may need.
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

template <typename FIELD>
bool is_neural_model(const FIELD& field, NeuralModelType type)
{
	if constexpr (FIELD::is_neural_type)
		return field.model_type() == type;
	else
		return false;
}

} // namespace offset_axis

} // namespace geometry

} // namespace cgogn

#endif // CGOGN_GEOMETRY_ALGOS_OFFSET_AXIS_DISTANCE_FIELD_H_
