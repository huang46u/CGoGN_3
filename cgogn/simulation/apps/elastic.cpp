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

#include <cgogn/core/types/maps/gmap/gmap3.h>
#include <cgogn/core/types/maps/gmap/gmap2.h>
// #include <cgogn/core/types/maps/cmap/cmap3.h>

#include <cgogn/geometry/types/vector_traits.h>

#include <cgogn/ui/app.h>
#include <cgogn/ui/view.h>

#include <cgogn/core/ui_modules/mesh_provider.h>
#include <cgogn/rendering/ui_modules/volume_render.h>
#include <cgogn/rendering/ui_modules/surface_render.h>
#include <cgogn/simulation/ui_modules/elastic.h>


#define DEFAULT_MESH_PATH CGOGN_STR(CGOGN_DATA_PATH) "/meshes/"

using Volume = cgogn::GMap3;
using Surface = cgogn::GMap2;
// using Mesh = cgogn::CMap3;

template <typename T>
using VolumeAttribute = typename cgogn::mesh_traits<Volume>::Attribute<T>;
template <typename T>
using SurfaceAttribute = typename cgogn::mesh_traits<Surface>::Attribute<T>;

int main(int argc, char** argv)
{

	using Vec3 = cgogn::geometry::Vec3;
	using Scalar = cgogn::geometry::Scalar;

	std::string objfilename;
	if (argc < 2)
	{
		objfilename = std::string(DEFAULT_MESH_PATH) + std::string("tet/hand.tet");
	}
	else
		objfilename = std::string(argv[1]);

	std::string plan = std::string(DEFAULT_MESH_PATH) + std::string("off/plan.off");

	cgogn::thread_start();
	cgogn::ui::App app;
	app.set_window_title("Simple volume viewer");
	app.set_window_size(1000, 800);

	cgogn::ui::MeshProvider<Volume> mp(app);
	cgogn::ui::MeshProvider<Surface> mp_surf(app);
	cgogn::ui::VolumeRender<Volume> vr(app);
	cgogn::ui::SurfaceRender<Surface> sr(app);
	//cgogn::ui::Elastic<Volume> el(app);

	app.init_modules();

	cgogn::ui::View* v1 = app.current_view();
	v1->link_module(&mp);
	v1->link_module(&vr);

	v1->link_module(&mp_surf);
	v1->link_module(&sr);

	Volume* m = mp.load_volume_from_file(objfilename);
	Surface* s = mp_surf.load_surface_from_file(plan);
	if (!m)
	{
		std::cout << "Gmap3 could not be loaded" << std::endl;
		return 1;
	}
	if (!s)
	{
		std::cout << "Gmap2 could not be loaded" << std::endl;
		return 1;
	}

	return app.launch();
}
