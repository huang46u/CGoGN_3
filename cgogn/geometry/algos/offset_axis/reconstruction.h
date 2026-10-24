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
#ifndef CGOGN_GEOMETRY_ALGOS_OFFSET_AXIS_RECONSTRUCTION_H_
#define CGOGN_GEOMETRY_ALGOS_OFFSET_AXIS_RECONSTRUCTION_H_

#include <cgogn/core/functions/attributes.h>
#include <cgogn/core/functions/mesh_info.h>
#include <cgogn/geometry/algos/offset_axis/reconstruction_stages.h>
#include <cgogn/geometry/algos/offset_axis/reconstruction_state.h>
#include <cgogn/geometry/algos/offset_axis/reconstruction_stats.h>
#include <cgogn/geometry/types/vector_traits.h>

#include <algorithm>
#include <filesystem>
#include <limits>
#include <memory>
#include <optional>
#include <unordered_map>
#include <vector>

namespace cgogn
{

namespace geometry
{

namespace offset_axis
{

template <typename SURFACE, typename POINTS, typename NONMANIFOLD, template <typename, typename> class FIELD>
class OffsetAxisReconstruction
{
public:
	using State = ReconstructionState<SURFACE, POINTS, NONMANIFOLD, FIELD>;

	using PVertex = typename State::PVertex;
	using SVertex = typename State::SVertex;
	using SFace = typename State::SFace;
	using NMVertex = typename State::NMVertex;
	using NMEdge = typename State::NMEdge;
	using NMFace = typename State::NMFace;
	using SpheresOptimizerType = typename State::SpheresOptimizerType;
	using SpheresOptimizerMetrics = typename State::SpheresOptimizerMetrics;
	using TopologyType = typename State::TopologyType;
	using Field = typename State::Field;
	using RayTraits = typename State::RayTraits;
	using RaySampler = typename State::RaySampler;
	using RayParams = typename State::RayParams;
	using PAttribute = typename State::PAttribute;
	template <typename T>
	using PointAttribute = typename State::template PointAttribute<T>;
	template <typename T>
	using SurfaceAttribute = typename State::template SurfaceAttribute<T>;
	template <typename T>
	using SkeletonAttribute = typename State::template SkeletonAttribute<T>;

	using NeuralModelType = typename State::NeuralModelType;
	using Status = typename State::Status;
	using StatusInfo = typename State::StatusInfo;
	using Data = typename State::Data;
	using Options = typename State::Options;
	using Timing = typename State::Timing;
	using Counts = typename State::Counts;
	using Result = typename State::Result;
	using IterationCallback = typename State::IterationCallback;

	explicit OffsetAxisReconstruction(Data data, Options options = {})
	{
		s_.data = data;
		s_.options = options;
	}
	OffsetAxisReconstruction(const OffsetAxisReconstruction&) = delete;
	OffsetAxisReconstruction& operator=(const OffsetAxisReconstruction&) = delete;

	State& state()
	{
		return s_;
	}
	const State& state() const
	{
		return s_;
	}

	Data& data()
	{
		s_.input.invalidate();
		return s_.data;
	}
	const Data& data() const
	{
		return s_.data;
	}
	Options& options()
	{
		return s_.options;
	}
	const Options& options() const
	{
		return s_.options;
	}

	void reset()
	{
		stages::reset(s_);
	}
	StatusInfo initialize_data()
	{
		return stages::initialize_data(s_);
	}
	StatusInfo load_neural_model(const std::filesystem::path& path, NeuralModelType type)
	{
		return stages::load_neural_model(s_, path, type);
	}
	StatusInfo sample_alpha_level_set()
	{
		return stages::sample_alpha_level_set(s_);
	}
	StatusInfo apply_sampling_filtering()
	{
		return stages::apply_sampling_filtering(s_);
	}
	StatusInfo build_kdtree_and_normals()
	{
		return stages::build_kdtree_and_normals(s_);
	}
	StatusInfo compute_fitting_primitives()
	{
		return stages::compute_fitting_primitives(s_);
	}
	StatusInfo compute_initial_medial_axis()
	{
		return stages::compute_initial_medial_axis(s_);
	}
	StatusInfo initialize_spheres()
	{
		return stages::initialize_spheres(s_);
	}
	StatusInfo optimize_spheres(const IterationCallback& iteration_callback = {})
	{
		return stages::optimize_spheres(s_, iteration_callback);
	}
	StatusInfo build_skeleton()
	{
		return stages::build_skeleton(s_);
	}
	StatusInfo fix_topology()
	{
		return stages::fix_topology(s_);
	}
	StatusInfo prune_residual_sheets()
	{
		return stages::prune_residual_sheets(s_);
	}
	StatusInfo sample()
	{
		return stages::sample(s_);
	}
	StatusInfo compute_fitting_data()
	{
		return stages::compute_fitting_data(s_);
	}
	typename SpheresOptimizerType::Status update_spheres_once(Scalar sqem_update_lambda_line_plane)
	{
		return stages::update_spheres_once(s_, sqem_update_lambda_line_plane);
	}
	StatusInfo fix_topology_and_prune()
	{
		return stages::fix_topology_and_prune(s_);
	}
	Result run(const IterationCallback& iteration_callback = {})
	{
		return stages::reconstruct(s_, iteration_callback);
	}

	Counts counts() const
	{
		return stages::counts(s_);
	}
	ReconstructionRepairs repairs() const
	{
		return stages::repairs(s_);
	}

	NONMANIFOLD* skeleton() const
	{
		return s_.data.skeleton;
	}
	const std::shared_ptr<SkeletonAttribute<Vec3>>& skeleton_position_attribute() const
	{
		return s_.skeleton_position;
	}
	const std::shared_ptr<SkeletonAttribute<Scalar>>& skeleton_radius_attribute() const
	{
		return s_.skeleton_radius;
	}
	SpheresOptimizerType& spheres_optimizer()
	{
		return stages::optimizer_for(s_);
	}
	TopologyType& skeleton_topology()
	{
		return stages::topology_for(s_);
	}

	// Read-only display scores for the training UI, indexed by index_of(face) / index_of(edge). Cells that are not
	// adjacent to a current tet get NaN. Returns nullopt if the skeleton is missing or the UDF cannot be evaluated.
	std::optional<std::vector<Scalar>> skeleton_face_scores(bool normalize_by_area)
	{
		if (!s_.data.skeleton)
			return std::nullopt;
		std::unordered_map<uint32, Scalar> scores;
		auto evaluator = [this](const std::vector<Vec3>& points, std::vector<Scalar>& values) {
			return stages::evaluate_topology_scores(s_, points, values);
		};
		if (!skeleton_topology().tet_face_udf_scores(evaluator, normalize_by_area, scores))
			return std::nullopt;
		return dense_scores<NMFace>(scores);
	}

	std::optional<std::vector<Scalar>> skeleton_edge_scores_gauss3()
	{
		if (!s_.data.skeleton)
			return std::nullopt;
		std::unordered_map<uint32, Scalar> scores;
		auto evaluator = [this](const std::vector<Vec3>& points, std::vector<Scalar>& values) {
			return stages::evaluate_topology_scores(s_, points, values);
		};
		if (!skeleton_topology().tet_edge_gauss3_scores(evaluator, scores))
			return std::nullopt;
		return dense_scores<NMEdge>(scores);
	}

private:
	template <typename CELL>
	std::vector<Scalar> dense_scores(const std::unordered_map<uint32, Scalar>& scores)
	{
		NONMANIFOLD& skeleton = *s_.data.skeleton;
		uint32 size = 0;
		foreach_cell(skeleton, [&](CELL c) -> bool {
			const uint32 id = index_of(skeleton, c);
			if (id != INVALID_INDEX)
				size = std::max(size, id + 1);
			return true;
		});
		std::vector<Scalar> values(size, std::numeric_limits<Scalar>::quiet_NaN());
		for (const auto& kv : scores)
		{
			if (kv.first < size)
				values[kv.first] = kv.second;
		}
		return values;
	}

	State s_;
};

} // namespace offset_axis

} // namespace geometry

} // namespace cgogn

#endif // CGOGN_GEOMETRY_ALGOS_OFFSET_AXIS_RECONSTRUCTION_H_
