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
#ifndef CGOGN_GEOMETRY_ALGOS_OFFSET_AXIS_NEURAL_MODEL_TYPE_H_
#define CGOGN_GEOMETRY_ALGOS_OFFSET_AXIS_NEURAL_MODEL_TYPE_H_

namespace cgogn
{

namespace geometry
{

namespace offset_axis
{

// output layout of a neural field: unsigned distance only, or medial field (udf + sdf)
enum class NeuralModelType
{
	udf,
	mf
};

// MF fields vanish on the medial axis, not on the surface: their samples lie away from the input and must not be
// filtered by their distance to it, and their normals are re-estimated by PCA over neighbouring samples (the network
// gradient is not orthogonal to the medial axis)
inline bool default_apply_filtering(NeuralModelType type)
{
	return type != NeuralModelType::mf;
}

inline bool default_recompute_normals_after_sampling(NeuralModelType type)
{
	return type == NeuralModelType::mf;
}

} // namespace offset_axis

} // namespace geometry

} // namespace cgogn

#endif // CGOGN_GEOMETRY_ALGOS_OFFSET_AXIS_NEURAL_MODEL_TYPE_H_
