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

#ifndef CGOGN_SIMULATION_SCENE_LOADER_H_
#define CGOGN_SIMULATION_SCENE_LOADER_H_

#include <cgogn/core/types/maps/cmap/cmap2.h>
#include <cgogn/core/types/maps/cmap/cmap3.h>
#include <cgogn/geometry/types/vector_traits.h>
#include <cgogn/ui/app.h>
#include <cgogn/core/ui_modules/mesh_provider.h>

#include <nlohmann/json.hpp>
#include <fstream>
#include <string>
#include <vector>
#include <unordered_map>
#include <memory>
#include <cmath>

// Define DEFAULT_MESH_PATH if not already defined
#ifndef DEFAULT_MESH_PATH
#define DEFAULT_MESH_PATH CGOGN_STR(CGOGN_DATA_PATH)"meshes/"
#endif

namespace cgogn
{

namespace io
{

using Surface = cgogn::CMap2;
using Volume = cgogn::CMap3;
using SurfaceVertex = typename mesh_traits<Surface>::Vertex;
using VolumeVertex = typename mesh_traits<Volume>::Vertex;	
using Vec3 = cgogn::geometry::Vec3;
using Scalar = cgogn::geometry::Scalar;

struct ModelInfo
{
    std::string name;
    std::string type;  // "surface" or "volume"
    std::string file_path;
    Vec3 position{0.0f, 0.0f, 0.0f};
    Vec3 rotation{0.0f, 0.0f, 0.0f};
    Vec3 scale{1.0f, 1.0f, 1.0f};
    
    // Physics properties
    Scalar mass{1.0f};
    bool is_static{false};
    std::vector<uint32_t> fixed_vertices;
    Scalar stiffness{1.0f};

    // Pointers to loaded meshes
    Surface* surface_mesh{nullptr};
    Volume* volume_mesh{nullptr};
};

struct SceneInfo
{
    std::string scene_name;
    std::vector<ModelInfo> models;
    Vec3 gravity{0.0f, -9.81f, 0.0f};
    Scalar time_step{0.016f};
};

class SceneLoader
{
public:
	SceneLoader(cgogn::ui::App& app) : 
        app_(app), 
        surface_provider_(app), 
        volume_provider_(app)
    {
    }

    SceneInfo load_scene(const std::string& json_path)
    {
        SceneInfo scene;
        
        try
        {
            // Read and parse JSON file
            std::ifstream file(json_path);
            if (!file.is_open())
            {
                std::cerr << "Unable to open scene file: " << json_path << std::endl;
                return scene;
            }

            nlohmann::json json_data;
            file >> json_data;

            // Read scene basic information
            scene.scene_name = json_data.value("scene_name", "Unnamed Scene");
            
            // Read simulation parameters
            if (json_data.contains("simulation"))
            {
                auto& sim = json_data["simulation"];
                if (sim.contains("gravity") && sim["gravity"].is_array() && sim["gravity"].size() == 3)
                {
                    scene.gravity[0] = sim["gravity"][0];
                    scene.gravity[1] = sim["gravity"][1];
                    scene.gravity[2] = sim["gravity"][2];
                }
                scene.time_step = sim.value("time_step", 0.016f);
            }

            // Load models
            if (json_data.contains("models") && json_data["models"].is_array())
            {
                for (auto& model_json : json_data["models"])
                {
                    ModelInfo model;
                    model.name = model_json.value("name", "unnamed_model");
                    model.type = model_json.value("type", "surface");
                    model.file_path = model_json.value("file_path", "");

                    // Process DEFAULT_MESH_PATH macro
                    if (model.file_path.find("DEFAULT_MESH_PATH") == 0)
                    {
                        model.file_path.replace(0, 17, std::string(DEFAULT_MESH_PATH));
                    }

                    // Read position, rotation, scale
                    if (model_json.contains("position") && model_json["position"].is_array() && model_json["position"].size() == 3)
                    {
                        model.position[0] = model_json["position"][0];
                        model.position[1] = model_json["position"][1];
                        model.position[2] = model_json["position"][2];
                    }
                    
                    if (model_json.contains("rotation") && model_json["rotation"].is_array() && model_json["rotation"].size() == 3)
                    {
                        model.rotation[0] = model_json["rotation"][0];
                        model.rotation[1] = model_json["rotation"][1];
                        model.rotation[2] = model_json["rotation"][2];
                    }
                    
                    if (model_json.contains("scale") && model_json["scale"].is_array() && model_json["scale"].size() == 3)
                    {
                        model.scale[0] = model_json["scale"][0];
                        model.scale[1] = model_json["scale"][1];
                        model.scale[2] = model_json["scale"][2];
                    }

                    // Read physics properties
                    if (model_json.contains("physics"))
                    {
                        auto& physics = model_json["physics"];
                        model.mass = physics.value("mass", 1.0f);
                        model.is_static = physics.value("is_static", false);
                        model.stiffness = physics.value("stiffness", 1.0f);
                        
                        if (physics.contains("fixed_vertices") && physics["fixed_vertices"].is_array())
                        {
                            for (auto& v : physics["fixed_vertices"])
                            {
                                model.fixed_vertices.push_back(v);
                            }
                        }
                    }

                    // Load model mesh
                    if (model.type == "surface")
                    {
                        model.surface_mesh = surface_provider_.load_surface_from_file(model.file_path);
                        if (!model.surface_mesh)
                        {
                            std::cerr << "Failed to load surface mesh: " << model.file_path << std::endl;
                        }
						
                    }
                    else if (model.type == "volume")
                    {
                        model.volume_mesh = volume_provider_.load_volume_from_file(model.file_path);
                        if (!model.volume_mesh)
                        {
                            std::cerr << "Failed to load volume mesh: " << model.file_path << std::endl;
                        }
                    }

                    // Apply transformation (position, rotation, scale)
                    apply_transformation(model);

                    // Add to scene
                    scene.models.push_back(model);
                }
            }
        }
        catch (const std::exception& e)
        {
            std::cerr << "Error loading scene: " << e.what() << std::endl;
        }

        return scene;
    }

    Surface* get_surface_mesh(const std::string& model_name, const SceneInfo& scene)
    {
        for (const auto& model : scene.models)
        {
            if (model.name == model_name && model.type == "surface")
            {
                return model.surface_mesh;
            }
        }
        return nullptr;
    }

    Volume* get_volume_mesh(const std::string& model_name, const SceneInfo& scene)
    {
        for (const auto& model : scene.models)
        {
            if (model.name == model_name && model.type == "volume")
            {
                return model.volume_mesh;
            }
        }
        return nullptr;
    }

    cgogn::ui::MeshProvider<Surface>& surface_provider()
	{
		return surface_provider_;
	}
	cgogn::ui::MeshProvider<Volume>& volume_provider()
	{
		return volume_provider_;
    }
private:
    void apply_transformation(ModelInfo& model)
    {
        // Apply transformation to mesh vertices
        if (model.type == "surface" && model.surface_mesh)
        {
            auto position_attribute = get_attribute<Vec3, SurfaceVertex>(*model.surface_mesh, "position");
            if (position_attribute)
            {
                // Use CGOGN's foreach_cell method to iterate through all vertices and apply transformation
                parallel_foreach_cell(*model.surface_mesh, [&](SurfaceVertex v) -> bool {
					uint32 index_v = index_of(*model.surface_mesh, v);
					Vec3& pos = (*position_attribute)[index_v];
                    
                    // Apply scaling
                    pos[0] *= model.scale[0];
                    pos[1] *= model.scale[1];
                    pos[2] *= model.scale[2];

                    // Implement rotation (based on XYZ Euler angles)
                    // Convert angles to radians
                    Scalar rx = model.rotation[0] * M_PI / 180.0f;
                    Scalar ry = model.rotation[1] * M_PI / 180.0f;
                    Scalar rz = model.rotation[2] * M_PI / 180.0f;

                    // Save original position
                    Scalar x = pos[0];
                    Scalar y = pos[1];
                    Scalar z = pos[2];

                    // Rotate around X-axis
                    Scalar y1 = y * std::cos(rx) - z * std::sin(rx);
                    Scalar z1 = y * std::sin(rx) + z * std::cos(rx);
                    
                    // Rotate around Y-axis
                    Scalar x2 = x * std::cos(ry) + z1 * std::sin(ry);
                    Scalar z2 = -x * std::sin(ry) + z1 * std::cos(ry);
                    
                    // Rotate around Z-axis
                    Scalar x3 = x2 * std::cos(rz) - y1 * std::sin(rz);
                    Scalar y3 = x2 * std::sin(rz) + y1 * std::cos(rz);
                    
                    // Update position
                    pos[0] = x3;
                    pos[1] = y3;
                    pos[2] = z2;
                    
                    // Apply translation
                    pos[0] += model.position[0];
                    pos[1] += model.position[1];
                    pos[2] += model.position[2];
                    
                    return true; // Continue iteration
                });
            }
        }
        else if (model.type == "volume" && model.volume_mesh)
        {
            auto position_attribute = get_attribute<Vec3, VolumeVertex>(*model.volume_mesh, "position");
            if (position_attribute)
            {
                // Use CGOGN's foreach_cell method to iterate through all vertices and apply transformation
                parallel_foreach_cell(*model.volume_mesh, [&](VolumeVertex v) -> bool {
					uint32 index_v = index_of(*model.volume_mesh, v);
					Vec3& pos = (*position_attribute)[index_v];
                    
                    // Apply scaling
                    pos[0] *= model.scale[0];
                    pos[1] *= model.scale[1];
                    pos[2] *= model.scale[2];

                    // Implement rotation (based on XYZ Euler angles)
                    // Convert angles to radians
                    Scalar rx = model.rotation[0] * M_PI / 180.0f;
                    Scalar ry = model.rotation[1] * M_PI / 180.0f;
                    Scalar rz = model.rotation[2] * M_PI / 180.0f;

                    // Save original position
                    Scalar x = pos[0];
                    Scalar y = pos[1];
                    Scalar z = pos[2];

                    // Rotate around X-axis
                    Scalar y1 = y * std::cos(rx) - z * std::sin(rx);
                    Scalar z1 = y * std::sin(rx) + z * std::cos(rx);
                    
                    // Rotate around Y-axis
                    Scalar x2 = x * std::cos(ry) + z1 * std::sin(ry);
                    Scalar z2 = -x * std::sin(ry) + z1 * std::cos(ry);
                    
                    // Rotate around Z-axis
                    Scalar x3 = x2 * std::cos(rz) - y1 * std::sin(rz);
                    Scalar y3 = x2 * std::sin(rz) + y1 * std::cos(rz);
                    
                    // Update position
                    pos[0] = x3;
                    pos[1] = y3;
                    pos[2] = z2;
                    
                    // Apply translation
                    pos[0] += model.position[0];
                    pos[1] += model.position[1];
                    pos[2] += model.position[2];
                    
                    return true; // Continue iteration
                });
            }
        }
    }

private:
    cgogn::ui::App& app_;
    cgogn::ui::MeshProvider<Surface> surface_provider_;
    cgogn::ui::MeshProvider<Volume> volume_provider_;
};

} // namespace io
} // namespace cgogn

#endif // CGOGN_SIMULATION_SCENE_LOADER_H_
