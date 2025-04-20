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

#include <cgogn/io/scene/scene_loader.h> 

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

    cgogn::ui::MeshProvider<Surface> sp(app);
	cgogn::ui::MeshProvider<Volume> vp(app);
    cgogn::ui::SurfaceRender<Surface> sr(app);
    cgogn::ui::VolumeRender<Volume> vr(app);
  
    
    // Create scene loader
    cgogn::io::SceneLoader scene_loader(app, sp, vp);

    // Set up view
    cgogn::ui::View* view = app.current_view();
	view->link_module(&sp);
	view->link_module(&vp);
    view->link_module(&sr);
    view->link_module(&vr);
	
    
    app.init_modules();

    // Load scene
    std::cout << "Starting to load scene: " << scene_file << std::endl;
    cgogn::io::SceneInfo scene = scene_loader.load_scene(scene_file);
	
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
