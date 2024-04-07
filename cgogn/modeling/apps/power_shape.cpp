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
#include <cgogn/modeling/ui_modules/power_shape.h>
#include <cgogn/rendering/ui_modules/point_cloud_render.h>
#include <cgogn/rendering/ui_modules/surface_render.h>
#include <cgogn/geometry/ui_modules/surface_selection.h>
#include <cgogn/geometry/ui_modules/surface_differential_properties.h>
#include <cgogn/modeling/ui_modules/surface_modeling.h>

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
	app.set_window_title("Power shape viewer");
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
	cgogn::ui::PowerShape<Point, Surface, NonManifold> pw(app);
	

	app.init_modules();

	cgogn::ui::View* v1 = app.current_view();

	v1->link_module(&mp);
	v1->link_module(&ms);
	v1->link_module(&mpnm);
	v1->link_module(&sl);
	v1->link_module(&sr);
	v1->link_module(&srnm);
	v1->link_module(&pcr);
	v1->link_module(&pw);

	
	if (filename.length() > 0)
	{
		Surface* m = ms.load_surface_from_file(filename);
		if (!m)
		{
			std::cout << "File could not be loaded" << std::endl;
			return 1;
		}

	
		auto surface_vertex_position = cgogn::get_attribute<Vec3, Vertex>(*m, "position");
		auto surface_vertex_normal =
			cgogn::get_or_add_attribute<Vec3, Vertex>(*m, "normal");
		auto surface_vertex_kmax = cgogn::get_or_add_attribute<Scalar, Vertex>(*m, "kmax");
		auto surface_vertex_kmin = cgogn::get_or_add_attribute<Scalar, Vertex>(*m, "kmin");
		auto surface_vertex_kgaussian = cgogn::get_or_add_attribute<Scalar, Vertex>(*m, "kgaussian");
		auto surface_vertex_Kmax = cgogn::get_or_add_attribute<Vec3, Vertex>(*m, "Kmax");
		auto surface_vertex_Kmin = cgogn::get_or_add_attribute<Vec3, Vertex>(*m, "Kmin");
		auto surface_vertex_Knormal = cgogn::get_or_add_attribute<Vec3, Vertex>(*m, "Knormal");
		
		sdp.compute_normal(*m, surface_vertex_position.get(), surface_vertex_normal.get());
		sr.set_vertex_position(*v1, *m, surface_vertex_position);
		sr.set_vertex_normal(*v1, *m, surface_vertex_normal);
		sr.set_render_edges(*v1, *m, false);
		sr.set_render_vertices(*v1, *m, false);

		sr.set_ghost_mode(*v1, *m, true);
		pcr.set_render_vertices(*v1, *m, true);

		/*std::shared_ptr<Attribute<Scalar>> edge_angle = cgogn::get_or_add_attribute<Scalar, Edge>(*m,
		   "__edge_angle");
		cgogn::geometry::compute_angle(*m, surface_vertex_position.get(), edge_angle.get());

		Scalar mean_edge_length = cgogn::geometry::mean_edge_length(*m, surface_vertex_position.get());
		
		sdp.compute_curvature(*m, mean_edge_length * 2.5, surface_vertex_position.get(), surface_vertex_normal.get(),
							  edge_angle.get(), surface_vertex_kmax.get(), surface_vertex_kmin.get(),
							  surface_vertex_kgaussian.get(), surface_vertex_Kmax.get(), surface_vertex_Kmin.get(),
							  surface_vertex_Knormal.get());*/

		
	}

	return app.launch();
}
