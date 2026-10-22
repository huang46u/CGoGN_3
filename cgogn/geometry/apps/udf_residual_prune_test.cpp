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

#include <cgogn/core/functions/attributes.h>
#include <cgogn/core/functions/mesh_info.h>
#include <cgogn/core/types/incidence_graph/incidence_graph.h>
#include <cgogn/core/types/maps/cmap/cmap0.h>
#include <cgogn/core/utils/thread.h>
// skeleton_topology.h is not self-contained; reconstruction.h pulls in the traversals it needs
#include <cgogn/geometry/algos/udf/reconstruction.h>

#include <Eigen/Dense>

#include <cstddef>
#include <cstdint>
#include <iostream>
#include <map>
#include <set>
#include <string>
#include <utility>
#include <vector>

namespace
{

using Graph = cgogn::IncidenceGraph;
using Points = cgogn::CMap0;
using GVertex = cgogn::mesh_traits<Graph>::Vertex;
using GEdge = cgogn::mesh_traits<Graph>::Edge;
using GFace = cgogn::mesh_traits<Graph>::Face;
using PVertex = cgogn::mesh_traits<Points>::Vertex;
using Topology = cgogn::geometry::SkeletonTopology<Points, Graph>;
using Optimizer = cgogn::geometry::SpheresOptimizer<Points>;
using Mode = Topology::ResidualPruneMode;

struct Counts
{
	long long vertices = 0;
	long long edges = 0;
	long long faces = 0;
	long long chi() const { return vertices - edges + faces; }
};

Counts counts(Graph& g)
{
	Counts c;
	c.vertices = static_cast<long long>(cgogn::nb_cells<GVertex>(g));
	c.edges = static_cast<long long>(cgogn::nb_cells<GEdge>(g));
	c.faces = static_cast<long long>(cgogn::nb_cells<GFace>(g));
	return c;
}

// builds triangles on a graph, reusing edges by their vertex pair
class Builder
{
public:
	explicit Builder(Graph& g) : g_(g) {}

	GVertex vertex() { return cgogn::add_vertex(g_); }

	GEdge edge(GVertex a, GVertex b)
	{
		const std::uint32_t ia = cgogn::index_of(g_, a);
		const std::uint32_t ib = cgogn::index_of(g_, b);
		const std::pair<std::uint32_t, std::uint32_t> key = ia < ib ? std::make_pair(ia, ib) : std::make_pair(ib, ia);
		auto it = edges_.find(key);
		if (it != edges_.end())
			return it->second;
		GEdge e = cgogn::add_edge(g_, a, b);
		edges_[key] = e;
		return e;
	}

	GFace triangle(GVertex a, GVertex b, GVertex c)
	{
		std::vector<GEdge> edges{edge(a, b), edge(b, c), edge(c, a)};
		return cgogn::add_face(g_, edges);
	}

	// strip (u,v,x), (v,x,y), (x,y,z) glued on edge uv, 3 faces, internal edges shared by 2 faces
	void strip_on(GVertex u, GVertex v)
	{
		GVertex x = vertex();
		GVertex y = vertex();
		GVertex z = vertex();
		triangle(u, v, x);
		triangle(v, x, y);
		triangle(x, y, z);
	}

	// two independent support strips on edge uv, so that uv is shared by 3 faces in total
	void support_edge(GVertex u, GVertex v)
	{
		strip_on(u, v);
		strip_on(u, v);
	}

private:
	Graph& g_;
	std::map<std::pair<std::uint32_t, std::uint32_t>, GEdge> edges_;
};

// scenario A: isolated fully non-manifold triangle T = (a,b,c)
void build_scenario_a(Graph& g, Builder& builder)
{
	(void)g;
	GVertex va = builder.vertex();
	GVertex vb = builder.vertex();
	GVertex vc = builder.vertex();
	builder.triangle(va, vb, vc);
	builder.support_edge(va, vb);
	builder.support_edge(vb, vc);
	builder.support_edge(vc, va);
}

// scenario B: two-face sheet T1 = (a,b,c), T2 = (a,c,d) sharing only edge ac, no boundary edge
void build_scenario_b(Graph& g, Builder& builder)
{
	(void)g;
	GVertex va = builder.vertex();
	GVertex vb = builder.vertex();
	GVertex vc = builder.vertex();
	GVertex vd = builder.vertex();
	builder.triangle(va, vb, vc);
	builder.triangle(va, vc, vd);
	builder.support_edge(va, vb);
	builder.support_edge(vb, vc);
	builder.support_edge(vc, vd);
	builder.support_edge(vd, va);
}

const char* mode_name(Mode mode)
{
	switch (mode)
	{
	case Mode::current:
		return "current";
	case Mode::keep_fully_non_manifold:
		return "keep_fully_non_manifold";
	case Mode::collapsible_only:
		return "collapsible_only";
	}
	return "unknown";
}

// returns true when every check of this (scenario, mode) pair passes
bool run_case(bool scenario_b, Mode mode)
{
	const char* scenario = scenario_b ? "B" : "A";

	Graph skeleton;
	Points spheres;

	// attributes are created before the cells, as in reconstruction.h
	auto sphere_skeleton_vertex = cgogn::add_attribute<GVertex, PVertex>(spheres, "skeleton_vertex");
	auto skeleton_position = cgogn::add_attribute<Eigen::Vector3d, GVertex>(skeleton, "position");
	auto skeleton_radius = cgogn::add_attribute<double, GVertex>(skeleton, "radius");
	auto skeleton_source_sphere = cgogn::add_attribute<PVertex, GVertex>(skeleton, "source_sphere");
	auto face_incident_tets = cgogn::add_attribute<std::set<std::size_t>, GFace>(skeleton, "incident_tets");
	auto edge_degree = cgogn::add_attribute<std::uint32_t, GEdge>(skeleton, "degree");

	Builder builder(skeleton);
	if (scenario_b)
		build_scenario_b(skeleton, builder);
	else
		build_scenario_a(skeleton, builder);

	Optimizer optimizer{Optimizer::Data{}};

	Topology::Data data;
	data.skeleton = &skeleton;
	data.spheres_optimizer = &optimizer;
	data.sphere_skeleton_vertex = sphere_skeleton_vertex;
	data.skeleton_position = skeleton_position;
	data.skeleton_radius = skeleton_radius;
	data.skeleton_source_sphere = skeleton_source_sphere;
	data.face_incident_tets = face_incident_tets;
	data.edge_degree = edge_degree;

	Topology topology(data);
	topology.set_residual_prune_mode(mode);

	const Counts before = counts(skeleton);
	const Topology::Status status = topology.prune_residual_sheets();
	const Counts after = counts(skeleton);
	const Topology::Metrics& m = topology.metrics();

	const long long chi_before = before.chi();
	const long long chi_after = after.chi();
	const long long deleted_faces = before.faces - after.faces;
	const long long removed_edges = before.edges - after.edges;
	const bool current = mode == Mode::current;
	const bool collapsible = mode == Mode::collapsible_only;

	// expectations
	long long expected_deleted_faces = 0;
	long long expected_removed_edges = 0;
	long long expected_fully_nm = 0;
	long long expected_residual_removed = 0;
	if (!scenario_b)
	{
		// T is removed by the fully non-manifold step only in current mode
		expected_deleted_faces = current ? 1 : 0;
		expected_fully_nm = current ? 1 : 0;
	}
	else
	{
		// T1, T2 form a 2-face sheet without boundary edge: deleted unless collapsible_only
		const bool deleted = !collapsible;
		expected_deleted_faces = deleted ? 2 : 0;
		expected_residual_removed = deleted ? 2 : 0;
		// edge ac loses both faces and is removed as an orphan edge, so chi drops by 1 not 2
		expected_removed_edges = deleted ? 1 : 0;
	}

	std::vector<std::string> errors;
	if (status != Topology::Status::success)
		errors.push_back("prune_residual_sheets did not return success");
	if (before.vertices != after.vertices)
		errors.push_back("vertex count changed (vertex removal path was not expected)");
	if (deleted_faces != expected_deleted_faces)
		errors.push_back("deleted faces: expected " + std::to_string(expected_deleted_faces) + ", got " +
						 std::to_string(deleted_faces));
	if (removed_edges != expected_removed_edges)
		errors.push_back("removed edges: expected " + std::to_string(expected_removed_edges) + ", got " +
						 std::to_string(removed_edges));
	if (chi_after != chi_before - expected_deleted_faces + expected_removed_edges)
		errors.push_back("chi_after: expected " + std::to_string(chi_before - expected_deleted_faces + expected_removed_edges) +
						 ", got " + std::to_string(chi_after));
	if (m.residual_chi_before != chi_before || m.residual_chi_after != chi_after)
		errors.push_back("metrics chi values do not match the skeleton");
	if (static_cast<long long>(m.residual_fully_nm_removed_faces) != expected_fully_nm)
		errors.push_back("residual_fully_nm_removed_faces: expected " + std::to_string(expected_fully_nm) + ", got " +
						 std::to_string(m.residual_fully_nm_removed_faces));
	if (scenario_b && static_cast<long long>(m.residual_removed_faces) != expected_residual_removed)
		errors.push_back("residual_removed_faces: expected " + std::to_string(expected_residual_removed) + ", got " +
						 std::to_string(m.residual_removed_faces));
	if (scenario_b && m.residual_small_sheets_without_boundary_edge < 1)
		errors.push_back("residual_small_sheets_without_boundary_edge is 0, expected >= 1");

	const bool ok = errors.empty();
	std::cout << "[scenario " << scenario << "][" << mode_name(mode) << "] chi_before=" << chi_before
			  << " chi_after=" << chi_after << " deleted_faces=" << deleted_faces
			  << " removed_edges=" << removed_edges << (ok ? " PASS" : " FAIL") << '\n';
	for (const std::string& e : errors)
		std::cout << "  FAIL: " << e << '\n';
	return ok;
}

} // namespace

int main()
{
	cgogn::thread_start(0);
	int failed_cases = 0;
	const Mode modes[] = {Mode::current, Mode::keep_fully_non_manifold, Mode::collapsible_only};
	for (bool scenario_b : {false, true})
		for (Mode mode : modes)
			if (!run_case(scenario_b, mode))
				++failed_cases;
	cgogn::thread_stop();

	if (failed_cases == 0)
	{
		std::cout << "PASS\n";
		return 0;
	}
	std::cout << "FAIL: " << failed_cases << " case(s) failed\n";
	return 1;
}
