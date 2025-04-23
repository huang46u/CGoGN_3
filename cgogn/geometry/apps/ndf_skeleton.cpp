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

#include <algorithm>

#include <cgogn/core/types/maps/cmap/cmap2.h>
#include <cgogn/geometry/functions/bounding_box.h>
#include <cgogn/geometry/types/vector_traits.h>

#include <cgogn/ui/app.h>
#include <cgogn/ui/view.h>

#include <cgogn/core/ui_modules/mesh_provider.h>
#include <cgogn/geometry/ui_modules/ndf_skeleton.h>
#include <cgogn/rendering/ui_modules/point_cloud_render.h>
#include <cgogn/rendering/ui_modules/surface_render.h>

#define DEFAULT_MESH_PATH CGOGN_STR(CGOGN_DATA_PATH) "/meshes/"

using namespace cgogn::numerics;

using Surface = cgogn::CMap2;
using Points = cgogn::CMap0;

template <typename T>
using SurfaceAttribute = typename cgogn::mesh_traits<Surface>::Attribute<T>;
using SurfaceVertex = typename cgogn::mesh_traits<Surface>::Vertex;

template <typename T>
using PointsAttribute = typename cgogn::mesh_traits<Points>::Attribute<T>;
using PointsVertex = typename cgogn::mesh_traits<Points>::Vertex;

using cgogn::geometry::Scalar;
using cgogn::geometry::Vec3;
using cgogn::geometry::Vec4;

template <template <typename VEC> typename CONTAINER, typename VEC>
void center_at(CONTAINER<VEC>& container, const VEC& center)
{
	using Scalar = typename cgogn::geometry::vector_traits<VEC>::Scalar;
	const std::size_t dimension = cgogn::geometry::vector_traits<VEC>::SIZE;

	VEC centroid;
	for (std::size_t i = 0; i < dimension; ++i)
		centroid[i] = 0;
	uint32 count = 0;
	for (const VEC& v : container)
	{
		for (std::size_t i = 0; i < dimension; ++i)
			centroid[i] += v[i];
		count++;
	}
	centroid /= Scalar(count);
	for (VEC& v : container)
	{
		for (std::size_t i = 0; i < dimension; ++i)
			v[i] -= centroid[i];
	}
}

int main(int argc, char** argv)
{
	std::string surface_filename;
	std::string model_filename;
	if (argc < 3)
	{
		std::cout << "Usage: " << argv[0] << " <surface_file> <model_file>" << std::endl;
		std::exit(1);
	}
	else
	{
		surface_filename = std::string(argv[1]);
		model_filename = std::string(argv[2]);
	}

	cgogn::thread_start();

	cgogn::ui::App app;
	app.set_window_title("NDF Skeleton");
	app.set_window_size(1000, 800);

	cgogn::ui::MeshProvider<Surface> mps(app);
	cgogn::ui::MeshProvider<Points> mpp(app);

	cgogn::ui::SurfaceRender<Surface> sr(app);
	cgogn::ui::PointCloudRender<Points> pcr(app);

	cgogn::ui::LipNDF<Surface, Points> lip_ndf(app);

	app.init_modules();

	cgogn::ui::View* v1 = app.current_view();
	v1->link_module(&mps);
	v1->link_module(&mpp);
	v1->link_module(&sr);
	v1->link_module(&pcr);

	Surface* s = mps.load_surface_from_file(surface_filename);
	if (!s)
	{
		std::cout << "File could not be loaded" << std::endl;
		return 1;
	}

	auto surface_vertex_position = cgogn::get_attribute<Vec3, SurfaceVertex>(*s, "position");

	center_at(*surface_vertex_position.get(), {0, 0, 0});
	mps.emit_attribute_changed(*s, surface_vertex_position.get());
	v1->pivot_around_scene_center();

	lip_ndf.set_selected_surface(*s);
	lip_ndf.set_surface_vertex_position(*s, surface_vertex_position);
	lip_ndf.init_surface_data(*s, model_filename);

	sr.set_vertex_position(*v1, *s, surface_vertex_position);
	sr.set_render_vertices(*v1, *s, false);
	sr.set_render_edges(*v1, *s, false);

	Points* pc = mpp.mesh(mps.mesh_name(*s) + "_spheres");
	auto points_vertex_position = cgogn::get_attribute<Vec3, PointsVertex>(*pc, "position");
	auto points_vertex_radius = cgogn::get_attribute<Scalar, PointsVertex>(*pc, "radius");
	auto points_vertex_color = cgogn::get_attribute<Vec4, PointsVertex>(*pc, "color");

	mpp.set_mesh_bb_vertex_position(*pc, points_vertex_position);

	pcr.set_vertex_position(*v1, *pc, points_vertex_position);
	pcr.set_vertex_radius(*v1, *pc, points_vertex_radius);
	pcr.set_vertex_color(*v1, *pc, points_vertex_color);
	pcr.set_vertex_color_per_cell(*v1, *pc, cgogn::ui::PointCloudRender<Points>::AttributePerCell::PER_VERTEX);

	return app.launch();
}
