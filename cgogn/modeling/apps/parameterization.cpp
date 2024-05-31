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
#include <cgogn/modeling/ui_modules/parameterization.h>
#include <cgogn/geometry/types/vector_traits.h>


#include <cgogn/ui/app.h>
#include <cgogn/ui/view.h>

#include <cgogn/core/ui_modules/mesh_provider.h>
#include <cgogn/rendering/ui_modules/surface_render.h>
#include <cgogn/rendering/ui_modules/surface_tex_render.h>
#define DEFAULT_MESH_PATH CGOGN_STR(CGOGN_DATA_PATH) "/meshes/"

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
	app.set_window_title("Parameterization viewer");
	app.set_window_size(1000, 800);

	cgogn::ui::MeshProvider<Surface> ms(app);
	cgogn::ui::SurfaceRender<Surface> sr(app);
	cgogn::ui::SurfaceTexRender<Surface> str(app);
	cgogn::ui::Parameterization<Surface> para(app);

	app.init_modules();

	cgogn::ui::View* v1 = app.current_view();

	v1->link_module(&ms);
	v1->link_module(&sr);
	v1->link_module(&str);
	v1->link_module(&para);

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

		sr.set_vertex_position(*v1, *m, surface_vertex_position);
		sr.set_vertex_normal(*v1, *m, surface_vertex_normal);
		str.chekered_texture();
		cgogn::index_cells<Edge>(*m);
		cgogn::index_cells<Face>(*m);
	}

	return app.launch();
}
