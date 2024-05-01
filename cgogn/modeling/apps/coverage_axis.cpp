/*******************************************************************************
 * CGoGN: Combinatorial and Geometric modeling with Generic N-dimensional Maps  *
 * Copyright (C), IGG Group, ICube, University of Strasbourg, France            *
 *                                                                              *
 * This library is free software; you can redistribute it and/or modify it      *
 * under the terms of the GNU Lesser General Public License as published by the *
 * Free Software Foundation; either version 2.1 of the License, or (at your     *
 * option) any later version.                                                   *
 *                                                                              *
 * This library is distributed in the hope that it will be useful, but WITHOUT  *
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or        *
 * FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License  *
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

#include <cgogn/core/types/incidence_graph/incidence_graph.h>
#include <cgogn/core/types/maps/cmap/cmap0.h>
#include <cgogn/core/types/maps/cmap/cmap2.h>

#include <cgogn/geometry/types/vector_traits.h>

#include <cgogn/ui/app.h>
#include <cgogn/ui/view.h>

#include <cgogn/core/ui_modules/mesh_provider.h>
#include <cgogn/geometry/ui_modules/surface_differential_properties.h>
#include <cgogn/geometry/ui_modules/surface_selection.h>
#include <cgogn/modeling/ui_modules/coverage_axis.h>
#include <cgogn/modeling/ui_modules/surface_modeling.h>
#include <cgogn/rendering/ui_modules/point_cloud_render.h>
#include <cgogn/rendering/ui_modules/surface_render.h>

#define DEFAULT_MESH_PATH CGOGN_STR(CGOGN_DATA_PATH) "/meshes/"

using Point = cgogn::CMap0;
using NonManifold = cgogn::IncidenceGraph;
using Surface = cgogn::CMap2;

template <typename T>
using Attribute = typename cgogn::mesh_traits<Surface>::Attribute<T>;
using Vertex = typename cgogn::mesh_traits<Surface>::Vertex;
using Edge = typename cgogn::mesh_traits<Surface>::Edge;
using Face = typename cgogn::mesh_traits<Surface>::Face;

using Vec3 = cgogn::geometry::Vec3;
using Scalar = cgogn::geometry::Scalar;

int main(int argc, char** argv)
{
	std::string filename;
	if (argc > 1)
		filename = std::string(argv[1]);

	cgogn::thread_start();

	cgogn::ui::App app;
	app.set_window_title("Coverage axis viewer");
	app.set_window_size(1000, 800);

	cgogn::ui::MeshProvider<Point> mp(app);
	cgogn::ui::MeshProvider<Surface> ms(app);
	cgogn::ui::MeshProvider<NonManifold> mpnm(app);
	cgogn::ui::PointCloudRender<Point> pcr(app);
	cgogn::ui::SurfaceModeling<Surface> sms(app);
	cgogn::ui::SurfaceSelection<Surface> sl(app);
	cgogn::ui::SurfaceDifferentialProperties<Surface> sdp(app);
	cgogn::ui::SurfaceRender<Surface> sr(app);
	cgogn::ui::SurfaceRender<NonManifold> srnm(app);
	cgogn::ui::CoverageAxis<Point, Surface, NonManifold> ca(app);

	app.init_modules();

	cgogn::ui::View* v1 = app.current_view();

	v1->link_module(&mp);
	v1->link_module(&ms);
	v1->link_module(&mpnm);
	v1->link_module(&sl);
	v1->link_module(&sr);
	v1->link_module(&srnm);
	v1->link_module(&pcr);
	v1->link_module(&ca);

	if (filename.length() > 0)
	{
		Surface* m = ms.load_surface_from_file(filename);
		if (!m)
		{
			std::cout << "File could not be loaded" << std::endl;
			return 1;
		}

		auto surface_vertex_position = cgogn::get_attribute<Vec3, Vertex>(*m, "position");
		auto surface_vertex_normal = cgogn::get_or_add_attribute<Vec3, Vertex>(*m, "normal");

		sdp.compute_normal(*m, surface_vertex_position.get(), surface_vertex_normal.get());
		sr.set_vertex_position(*v1, *m, surface_vertex_position);
		sr.set_vertex_normal(*v1, *m, surface_vertex_normal);
		sr.set_render_edges(*v1, *m, false);
		sr.set_render_vertices(*v1, *m, false);

		sr.set_ghost_mode(*v1, *m, true);
		pcr.set_render_vertices(*v1, *m, true);
	}

	return app.launch();
}