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
#include <cgogn/geometry/algos/offset_axis/skeleton_topology.h>
#include <cgogn/geometry/algos/offset_axis/spheres_optimizer.h>

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
using Topology = cgogn::geometry::offset_axis::SkeletonTopology<Points, Graph>;
using Optimizer = cgogn::geometry::offset_axis::SpheresOptimizer<Points>;
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

// scenario C: complete 2-complex on 5 vertices (all 10 triangles). Every edge is shared by 3 triangles, so every
// triangle is fully non-manifold and all of them are removed in the same batch; every edge then becomes an orphan
// edge and every vertex an isolated vertex (covers the edge and vertex removal of the orphan cleanup).
void build_scenario_c(Graph& g, Builder& builder)
{
	(void)g;
	std::vector<GVertex> v;
	for (int i = 0; i < 5; ++i)
		v.push_back(builder.vertex());
	for (int i = 0; i < 5; ++i)
		for (int j = i + 1; j < 5; ++j)
			for (int k = j + 1; k < 5; ++k)
				builder.triangle(v[i], v[j], v[k]);
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
bool run_case(char scenario, Mode mode)
{
	const bool scenario_b = scenario == 'B';
	const bool scenario_c = scenario == 'C';

	Graph skeleton;
	Points spheres;

	// attributes are created before the cells
	auto sphere_skeleton_vertex = cgogn::add_attribute<GVertex, PVertex>(spheres, "skeleton_vertex");
	auto skeleton_position = cgogn::add_attribute<Eigen::Vector3d, GVertex>(skeleton, "position");
	auto skeleton_radius = cgogn::add_attribute<double, GVertex>(skeleton, "radius");
	auto skeleton_source_sphere = cgogn::add_attribute<PVertex, GVertex>(skeleton, "source_sphere");
	auto face_incident_tets = cgogn::add_attribute<std::set<std::size_t>, GFace>(skeleton, "incident_tets");
	auto edge_degree = cgogn::add_attribute<std::uint32_t, GEdge>(skeleton, "degree");

	Builder builder(skeleton);
	if (scenario_b)
		build_scenario_b(skeleton, builder);
	else if (scenario_c)
		build_scenario_c(skeleton, builder);
	else
		build_scenario_a(skeleton, builder);

	// the skeleton has no source spheres, so the optimizer attributes are never read
	Points samples;
	std::unique_ptr<acc::KDTree<3, std::uint32_t>> sample_kdtree;
	std::vector<PVertex> sample_kdtree_vertices;
	Optimizer optimizer{Optimizer::Data{samples, spheres, sample_kdtree, sample_kdtree_vertices}};

	Topology::Data data{skeleton, optimizer};
	data.sphere_skeleton_vertex = sphere_skeleton_vertex.get();
	data.skeleton_position = skeleton_position.get();
	data.skeleton_radius = skeleton_radius.get();
	data.skeleton_source_sphere = skeleton_source_sphere.get();
	data.face_incident_tets = face_incident_tets.get();
	data.edge_degree = edge_degree.get();

	Topology topology(data);

	const Counts before = counts(skeleton);
	topology.prune_residual_sheets(mode);
	const Counts after = counts(skeleton);

	const long long chi_before = before.chi();
	const long long chi_after = after.chi();
	const long long deleted_faces = before.faces - after.faces;
	const long long removed_edges = before.edges - after.edges;
	const long long removed_vertices = before.vertices - after.vertices;
	const bool current = mode == Mode::current;
	const bool collapsible = mode == Mode::collapsible_only;

	// expectations
	long long expected_deleted_faces = 0;
	long long expected_removed_edges = 0;
	long long expected_removed_vertices = 0;
	if (scenario_c)
	{
		// the 10 triangles are removed by the fully non-manifold step only in current mode, together with the 10
		// orphan edges and the 5 isolated vertices
		expected_deleted_faces = current ? 10 : 0;
		expected_removed_edges = current ? 10 : 0;
		expected_removed_vertices = current ? 5 : 0;
	}
	else if (!scenario_b)
	{
		// T is removed by the fully non-manifold step only in current mode
		expected_deleted_faces = current ? 1 : 0;
	}
	else
	{
		// T1, T2 form a 2-face sheet without boundary edge: deleted unless collapsible_only
		const bool deleted = !collapsible;
		expected_deleted_faces = deleted ? 2 : 0;
		// edge ac loses both faces and is removed as an orphan edge, so chi drops by 1 not 2
		expected_removed_edges = deleted ? 1 : 0;
	}

	std::vector<std::string> errors;
	if (removed_vertices != expected_removed_vertices)
		errors.push_back("removed vertices: expected " + std::to_string(expected_removed_vertices) + ", got " +
						 std::to_string(removed_vertices));
	if (deleted_faces != expected_deleted_faces)
		errors.push_back("deleted faces: expected " + std::to_string(expected_deleted_faces) + ", got " +
						 std::to_string(deleted_faces));
	if (removed_edges != expected_removed_edges)
		errors.push_back("removed edges: expected " + std::to_string(expected_removed_edges) + ", got " +
						 std::to_string(removed_edges));
	const long long expected_chi_after =
		chi_before - expected_deleted_faces + expected_removed_edges - expected_removed_vertices;
	if (chi_after != expected_chi_after)
		errors.push_back("chi_after: expected " + std::to_string(expected_chi_after) + ", got " +
						 std::to_string(chi_after));

	const bool ok = errors.empty();
	std::cout << "[scenario " << scenario << "][" << mode_name(mode) << "] chi_before=" << chi_before
			  << " chi_after=" << chi_after << " deleted_faces=" << deleted_faces
			  << " removed_edges=" << removed_edges << " removed_vertices=" << removed_vertices
			  << (ok ? " PASS" : " FAIL") << '\n';
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
	for (char scenario : {'A', 'B', 'C'})
		for (Mode mode : modes)
			if (!run_case(scenario, mode))
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
