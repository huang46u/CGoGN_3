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
#ifndef CGOGN_GEOMETRY_ALGOS_OFFSET_AXIS_RECONSTRUCTION_PROFILES_H_
#define CGOGN_GEOMETRY_ALGOS_OFFSET_AXIS_RECONSTRUCTION_PROFILES_H_

#include <cgogn/geometry/algos/offset_axis/neural_model_type.h>
#include <cgogn/geometry/algos/offset_axis/reconstruction_options.h>

#include <algorithm>
#include <string>
#include <vector>

namespace cgogn
{

namespace geometry
{

namespace offset_axis
{

// Dataset presets: named sets of sampling / postprocess defaults for one dataset family.
// A profile only overrides the fields listed for it. Profiles never change neural_model_type or seed.

inline const std::vector<std::string>& reconstruction_profile_names()
{
	static const std::vector<std::string> names = {"default", "mf_organic", "mf_abc", "udf"};
	return names;
}

// Applies the named profile on top of the current options.
// Returns false (and leaves options untouched) if the name is unknown.
inline bool apply_reconstruction_profile(const std::string& name, ReconstructionOptions& options)
{
	if (name == "default")
		return true;

	// Source: apps/configs/offset_axis_benchmark_batch_organic.json
	if (name == "mf_organic")
	{
		options.apply_filtering = false;
		options.recompute_normals_after_sampling = true;
		options.sample_radius = 0.002f;
		options.alpha_sampling_candidates_per_parent = 120;
		options.alpha_sampling_parent_batch_size = 102;
		options.init_dilation_constant = 0.01f;
		options.residual_prune = true;
		return true;
	}

	// Source: apps/configs/offset_axis_benchmark_batch_ABC.json
	if (name == "mf_abc")
	{
		options.apply_filtering = false;
		options.recompute_normals_after_sampling = true;
		options.sample_radius = 0.002f;
		options.alpha_sampling_candidates_per_parent = 60;
		options.alpha_sampling_parent_batch_size = 102;
		options.init_dilation_constant = 0.01f;
		options.residual_prune = true;
		return true;
	}

	// Source: apps/configs/offset_axis_benchmark_batch_deepfashion.json
	if (name == "udf")
	{
		options.alpha = 0.002f;
		options.sample_radius = 0.002f;
		options.init_dilation_constant = 0.01f;
		options.residual_prune = true;
		return true;
	}

	return false;
}

// Options of a reconstruction: the table defaults, then the defaults of the model type, then the dataset profile.
// Returns false if the profile is unknown; options then hold the defaults of the model type.
inline bool make_reconstruction_options(NeuralModelType type, const std::string& profile, ReconstructionOptions& options)
{
	options = ReconstructionOptions{};
	options.neural_model_type = type;
	options.apply_filtering = default_apply_filtering(type);
	options.recompute_normals_after_sampling = default_recompute_normals_after_sampling(type);
	return apply_reconstruction_profile(profile, options);
}

// Profile suggested for a neural model (used by the UI to pick a preset when a model is loaded).
// Neural MF model with "/ABC/" in its path -> mf_abc, other MF -> mf_organic, neural UDF -> udf, no model -> default.
inline std::string guess_reconstruction_profile(bool neural, NeuralModelType type, const std::string& model_path)
{
	if (!neural)
		return "default";

	if (type == NeuralModelType::mf)
	{
		std::string normalized = model_path;
		std::replace(normalized.begin(), normalized.end(), '\\', '/');
		if (normalized.find("/ABC/") != std::string::npos)
			return "mf_abc";
		return "mf_organic";
	}

	return "udf";
}

} // namespace offset_axis

} // namespace geometry

} // namespace cgogn

#endif // CGOGN_GEOMETRY_ALGOS_OFFSET_AXIS_RECONSTRUCTION_PROFILES_H_
