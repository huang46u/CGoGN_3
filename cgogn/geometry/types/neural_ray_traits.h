#ifndef CGOGN_GEOMETRY_TYPES_NEURAL_RAY_TRAITS_H_
#define CGOGN_GEOMETRY_TYPES_NEURAL_RAY_TRAITS_H_

#include <cgogn/geometry/types/neural_field_forward.h>
#include <cgogn/geometry/types/ray_level_set_sampler_traits.h>

#include <torch/torch.h>

#include <utility>

namespace cgogn
{

namespace geometry
{

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

} // namespace geometry

} // namespace cgogn

#endif // CGOGN_GEOMETRY_TYPES_NEURAL_RAY_TRAITS_H_
