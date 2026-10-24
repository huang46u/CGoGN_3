#ifndef CGOGN_GEOMETRY_ALGOS_OFFSET_AXIS_NEURAL_MODEL_TYPE_H_
#define CGOGN_GEOMETRY_ALGOS_OFFSET_AXIS_NEURAL_MODEL_TYPE_H_

namespace cgogn::geometry
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

} // namespace cgogn::geometry

#endif // CGOGN_GEOMETRY_ALGOS_OFFSET_AXIS_NEURAL_MODEL_TYPE_H_
