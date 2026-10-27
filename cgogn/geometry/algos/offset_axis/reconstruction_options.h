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
#ifndef CGOGN_GEOMETRY_ALGOS_OFFSET_AXIS_RECONSTRUCTION_OPTIONS_H_
#define CGOGN_GEOMETRY_ALGOS_OFFSET_AXIS_RECONSTRUCTION_OPTIONS_H_

#include <cgogn/geometry/algos/offset_axis/neural_model_type.h>

#include <string>

namespace cgogn
{

namespace geometry
{

namespace offset_axis
{

inline std::string to_string(NeuralModelType value)
{
	return value == NeuralModelType::mf ? "mf" : "udf";
}

inline bool parse(const std::string& value, NeuralModelType& out)
{
	if (value == "udf")
	{
		out = NeuralModelType::udf;
		return true;
	}
	if (value == "mf")
	{
		out = NeuralModelType::mf;
		return true;
	}
	return false;
}

// X(type, name, default, json key, min, max, ui_min, ui_max)
//   json key : path in the benchmark config / timing JSON (dotted)
//   min, max : accepted range (checked by the benchmark config loader), unused for bool/enum
//   ui_min, ui_max : slider range in the UI
// neural_model_type: type of the neural model to load; the loaded NeuralDistanceField keeps its own copy
// udf_min_active_rays: neural seed rays, a ray batch ends once fewer than this many rays are still marching (0 = off)
// init_owner_euclidean: initial sphere owners from a sphere BVH, the closest sphere surface (default); false scans all spheres
// clang-format off
#define CGOGN_OFFSET_AXIS_OPTIONS(X) \
  X(NeuralModelType,     neural_model_type,              NeuralModelType::udf,                  "input.neural_model_type",                   0, 0, 0, 0) \
  X(bool,                ma_flip_prune,                  true,                                  "input.ma_flip_prune",                       0, 1, 0, 1) \
  X(float,               ma_flip_prune_alpha_factor,     1.5f,                                "input.ma_flip_prune_alpha_factor",          1.0, 5.0, 1.0, 5.0) \
  X(float,               alpha,                          0.002f,                                "sampling.alpha",                            1e-6, 1.0, 0.0005, 0.05) \
  X(int,                 knn_k,                          10,                                    "sampling.knn_k",                            3, 1000, 3, 50) \
  X(bool,                apply_filtering,                true,                                  "sampling.apply_filtering",                  0, 1, 0, 1) \
  X(bool,                recompute_normals_after_sampling, false,                               "sampling.recompute_normals_after_sampling", 0, 1, 0, 1) \
  X(float,               sample_radius,                  0.002f,                                "sampling.bridson.sample_radius",            1e-6, 1.0, 0.0005, 0.02) \
  X(int,                 alpha_sampling_candidates_per_parent, 12,                              "sampling.bridson.candidates_per_parent",    1, 100000, 1, 200) \
  X(int,                 alpha_sampling_parent_batch_size, 1024,                                "sampling.bridson.parent_batch_size",        1, 1 << 24, 1, 4096) \
  X(int,                 ray_sampler_batch_size,         8192,                                  "sampling.ray.batch_size",                   1, 1 << 24, 1, 65536) \
  X(int,                 batch_size,                     1048576,                               "sampling.ray.query_batch_size",             1, 1 << 30, 1, 1 << 22) \
  X(float,               tolerance,                      1e-6f,                                 "sampling.ray.tolerance",                    1e-12, 1.0, 1e-9, 1e-3) \
  X(int,                 udf_max_iterations,             3000,                                  "sampling.ray.max_iterations",               1, 1 << 24, 1, 10000) \
  X(int,                 udf_min_active_rays,            100,                                   "sampling.ray.min_active_rays",              0, 1 << 24, 0, 10000) \
  X(float,               sqem_update_lambda_line_plane,  0.2f,                                  "optimization.sqem_update_lambda_line_plane", 0.0, 100.0, 0.0, 4.0) \
  X(int,                 seed,                           42,                                    "initialization.seed",                       0, 2147483647, 0, 1000) \
  X(float,               init_dilation_constant,         0.01f,                                 "initialization.init_dilation_constant",     1e-6, 1.0, 0.001, 0.01) \
  X(bool,                init_owner_euclidean,           true,                                  "initialization.init_owner_euclidean",       0, 1, 0, 1) \
  X(bool,                residual_prune,                 true,                                  "postprocess.residual_prune",                0, 1, 0, 1)
// clang-format on

struct ReconstructionOptions
{
#define CGOGN_OFFSET_AXIS_OPTION_MEMBER(T, name, def, key, lo, hi, ui_lo, ui_hi) T name = def;
	CGOGN_OFFSET_AXIS_OPTIONS(CGOGN_OFFSET_AXIS_OPTION_MEMBER)
#undef CGOGN_OFFSET_AXIS_OPTION_MEMBER
};

struct OptionInfo
{
	const char* name;
	const char* json_key;
	double min, max, ui_min, ui_max;
};

namespace option_info
{

#define CGOGN_OFFSET_AXIS_OPTION_INFO(T, name, def, key, lo, hi, ui_lo, ui_hi)                                                 \
	inline constexpr OptionInfo name{#name, key, lo, hi, ui_lo, ui_hi};
CGOGN_OFFSET_AXIS_OPTIONS(CGOGN_OFFSET_AXIS_OPTION_INFO)
#undef CGOGN_OFFSET_AXIS_OPTION_INFO
} // namespace option_info

// Options of a reconstruction with a field of the given type: the table defaults, which are the UDF ones, and for
// MF models the following changes. MF fields vanish on the medial axis, not on the surface: their samples lie away
// from the input and must not be filtered by their distance to it, and their normals are re-estimated by PCA over
// neighbouring samples (the network gradient is not orthogonal to the medial axis).
inline ReconstructionOptions default_options(NeuralModelType type)
{
	ReconstructionOptions options;
	options.neural_model_type = type;
	if (type == NeuralModelType::mf)
	{
		options.alpha = 0.005f;
		options.sample_radius = 0.0025f;
		options.alpha_sampling_candidates_per_parent = 60;
		options.alpha_sampling_parent_batch_size = 102;
		options.apply_filtering = false;
		options.recompute_normals_after_sampling = true;
	}
	return options;
}

// f(const OptionInfo&, T& value) for every option, in table order.
template <typename Options, typename F>
void for_each_option(Options& options, F&& f)
{
#define CGOGN_OFFSET_AXIS_OPTION_VISIT(T, name, def, key, lo, hi, ui_lo, ui_hi) f(option_info::name, options.name);
	CGOGN_OFFSET_AXIS_OPTIONS(CGOGN_OFFSET_AXIS_OPTION_VISIT)
#undef CGOGN_OFFSET_AXIS_OPTION_VISIT
}

} // namespace offset_axis

} // namespace geometry

} // namespace cgogn

#endif // CGOGN_GEOMETRY_ALGOS_OFFSET_AXIS_RECONSTRUCTION_OPTIONS_H_
