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

#include <cgogn/core/types/maps/cmap/cmap2.h>
#include <cgogn/core/types/maps/cmap/cmap3.h>
#include <cgogn/geometry/types/vector_traits.h>

#include <cgogn/ui/app.h>
#include <cgogn/ui/view.h>

#include <cgogn/core/ui_modules/mesh_provider.h>
#include <cgogn/rendering/ui_modules/surface_render.h>
#include <cgogn/rendering/ui_modules/volume_render.h>

#include <cgogn/io/scene/scene_loader.h> // Include the scene loader

using namespace cgogn::numerics;

using Surface = cgogn::CMap2;
using Volume = cgogn::CMap3;

using SurfaceVertex = typename cgogn::mesh_traits<Surface>::Vertex;
using VolumeVertex = typename cgogn::mesh_traits<Volume>::Vertex;
using Vec3 = cgogn::geometry::Vec3;
using Scalar = cgogn::geometry::Scalar;

int main(int argc, char** argv)
{
    std::string scene_file;
    if (argc < 2)
        scene_file = "D:/Code/CGoGN_3/data/scenes/simulation/test_scene.json";
    else
        scene_file = std::string(argv[1]);

    cgogn::thread_start();

    // Create application and view
    cgogn::ui::App app;
    app.set_window_title("Local-Global Solver");
    app.set_window_size(1200, 800);

    // Create necessary rendering modules
    cgogn::ui::SurfaceRender<Surface> surface_render(app);
    cgogn::ui::VolumeRender<Volume> volume_render(app);
    
    // Create scene loader
    cgogn::io::SceneLoader scene_loader(app);

    // Set up view
    cgogn::ui::View* view = app.current_view();
    view->link_module(&surface_render);
    view->link_module(&volume_render);
	view->link_module(&scene_loader.surface_provider());
	view->link_module(&scene_loader.volume_provider());
    
    app.init_modules();

    // Load scene
    std::cout << "Starting to load scene: " << scene_file << std::endl;
    cgogn::io::SceneInfo scene = scene_loader.load_scene(scene_file);
	for (const auto& model : scene.models)
	{
		if (model.type == "surface")
		{
			auto surface_mesh = scene_loader.get_surface_mesh(model.name, scene);
			if (surface_mesh)
			{
				auto vertex_position = cgogn::get_attribute<Vec3, SurfaceVertex>(*surface_mesh, "position");
				if (vertex_position)
				{
					surface_render.set_vertex_position(*view, *surface_mesh, vertex_position);
				}
			}
		}
		if (model.type == " volume")
		{
			auto volume_mesh = scene_loader.get_volume_mesh(model.name, scene);
			if (volume_mesh)
			{
				auto vertex_position = cgogn::get_attribute<Vec3, VolumeVertex>(*volume_mesh, "position");
				if (vertex_position)
				{
					volume_render.set_vertex_position(*view, *volume_mesh, vertex_position);
				}
			}
        }
		
    }
    // Output scene loading results
    if (scene.models.empty())
    {
        std::cout << "Scene loading failed or scene contains no models" << std::endl;
        return 1;
    }
    
    std::cout << "Successfully loaded scene: " << scene.scene_name << std::endl;
    std::cout << "Scene contains " << scene.models.size() << " models" << std::endl;
    std::cout << "Scene gravity setting: [" << scene.gravity[0] << ", " 
                                          << scene.gravity[1] << ", " 
                                          << scene.gravity[2] << "]" << std::endl;
    std::cout << "Scene time step: " << scene.time_step << std::endl;
    
    // Output information about each model
    for (const auto& model : scene.models)
    {
        std::cout << "Model: " << model.name << " (" << model.type << ")" << std::endl;
        std::cout << "  - Mass: " << model.mass << std::endl;
        std::cout << "  - Position: [" << model.position[0] << ", " 
                                      << model.position[1] << ", " 
                                      << model.position[2] << "]" << std::endl;
        
        if (model.is_static)
            std::cout << "  - Static object" << std::endl;
            
        if (!model.fixed_vertices.empty())
        {
            std::cout << "  - Fixed vertices count: " << model.fixed_vertices.size() << std::endl;
        }
    }

    // Launch application
    return app.launch();
}
