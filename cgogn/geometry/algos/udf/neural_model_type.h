#ifndef CGOGN_GEOMETRY_ALGOS_UDF_NEURAL_MODEL_TYPE_H_
#define CGOGN_GEOMETRY_ALGOS_UDF_NEURAL_MODEL_TYPE_H_

namespace cgogn::geometry
{

// output layout of a neural field: unsigned distance only, or medial field (udf + sdf)
enum class NeuralModelType
{
	udf,
	mf
};

} // namespace cgogn::geometry

#endif // CGOGN_GEOMETRY_ALGOS_UDF_NEURAL_MODEL_TYPE_H_
