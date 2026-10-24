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
#ifndef CGOGN_GEOMETRY_ALGOS_OFFSET_AXIS_NEURAL_NEURAL_RAY_TRAITS_H_
#define CGOGN_GEOMETRY_ALGOS_OFFSET_AXIS_NEURAL_NEURAL_RAY_TRAITS_H_

#include <cgogn/geometry/algos/offset_axis/neural/neural_field_forward.h>
#include <cgogn/geometry/algos/offset_axis/ray_level_set_sampler_traits.h>

#include <torch/torch.h>

#include <utility>

namespace cgogn
{

namespace geometry
{

namespace offset_axis
{

struct NeuralFieldRayTraits
{
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

} // namespace offset_axis

} // namespace geometry

} // namespace cgogn

#endif // CGOGN_GEOMETRY_ALGOS_OFFSET_AXIS_NEURAL_NEURAL_RAY_TRAITS_H_
