/*******************************************************************************
 * CGoGN                                                                        *
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

#ifndef CGOGN_SIMULATION_LOCAL_GLOBAL_H_
#define CGOGN_SIMULATION_LOCAL_GLOBAL_H_

#include <cgogn/ui/app.h>
#include <cgogn/ui/module.h>
#include <cgogn/ui/view.h>

#include <cgogn/core/ui_modules/mesh_provider.h>

#include <cgogn/geometry/types/vector_traits.h>
#include <cgogn/geometry/algos/area.h>
#include <cgogn/geometry/algos/length.h>
#include <cgogn/geometry/functions/distance.h>

#include <cgogn/io/scene/scene_loader.h>

#include <cgogn/simulation/type/constraint/base_constraint.h>

#include <cgogn/simulation/type/constraint/triangle_spring_constraint.h>
#include <cgogn/simulation/solver/projective_dynamics_solver.h>

#include <GLFW/glfw3.h>

#include <boost/synapse/connect.hpp>
#include <memory>
#include <mutex>
#include <vector>

#include <cgogn/simulation/type/constraint/pin_constraint.h>

namespace cgogn
{

namespace ui
{

using geometry::Scalar;
using geometry::Vec3;
using geometry::Vec4;

template <typename MESH>
class LocalGlobalSimulator : public ui::ViewModule
{
    template <typename T>
    using MeshAttribute = typename mesh_traits<MESH>::template Attribute<T>;
    
    using Vertex = typename mesh_traits<MESH>::Vertex;
    using Edge = typename mesh_traits<MESH>::Edge;
    using Face = typename mesh_traits<MESH>::Face;

public:
    struct SimulationParameters
    {
        bool initialized_ = false;
        bool running_ = false;
        bool stopping_ = false;
        bool slow_down_ = true;
        uint32 update_rate_ = 20;
        
        // Simulation parameters
        Vec3 gravity_ = Vec3(0.0, -9.81, 0.0);
        Scalar time_step_ = 0.016f;
        uint32 iteration_count_ = 0;
        
        // Mesh and attributes
        MESH* mesh_ = nullptr;
        std::shared_ptr<MeshAttribute<Vec3>> vertex_position_ = nullptr;
        std::shared_ptr<MeshAttribute<Vec3>> vertex_position_prev_ = nullptr;
        std::shared_ptr<MeshAttribute<Vec3>> vertex_velocity_ = nullptr;
        std::shared_ptr<MeshAttribute<Scalar>> vertex_mass_ = nullptr;
        std::shared_ptr<MeshAttribute<bool>> vertex_fixed_ = nullptr;
        
        // Constraints
        std::vector<std::unique_ptr<Constraint<MESH>>> constraints_;
        
        // Solver
        ProjectiveDynamicsSolver<MESH> solver_;
        
        std::mutex mutex_;
    };

    LocalGlobalSimulator(const ui::App& app)
        : ViewModule(app, "LocalGlobalSimulator (" + std::string{mesh_traits<MESH>::name} + ")")
    {
    }
    
    ~LocalGlobalSimulator()
    {
    }
    
    void set_selected_mesh(MESH& m)
    {
        selected_mesh_ = &m;
        
        // Set signal connections to update the data when the mesh connectivity or position changes
        if (mesh_connections_.find(&m) == mesh_connections_.end())
        {
            mesh_connections_[&m].push_back(
                boost::synapse::connect<typename ui::MeshProvider<MESH>::connectivity_changed>(
                    &m, [this, &m = m]() {
                        SimulationParameters& p = simulation_parameters_[&m];
                        init_simulation_data(m);
                    }));
                    
            mesh_connections_[&m].push_back(
                boost::synapse::connect<typename ui::MeshProvider<MESH>::template attribute_changed_t<Vec3>>(
                    &m, [this, &m = m](MeshAttribute<Vec3>* attribute) {
                        SimulationParameters& p = simulation_parameters_[&m];
                        if (attribute == p.vertex_position_.get())
                            init_simulation_data(m);
                    }));
        }
    }
    
    void set_mesh_vertex_position(MESH& m, const std::shared_ptr<MeshAttribute<Vec3>>& vertex_position)
    {
        SimulationParameters& p = simulation_parameters_[&m];
        p.vertex_position_ = vertex_position;
    }

    void init_simulation_data(MESH& m)
    {
        SimulationParameters& p = simulation_parameters_[&m];
        p.mesh_ = &m;
        
        if (!p.vertex_position_)
        {
            std::cout << "No vertex position attribute set" << std::endl;
            return;
        }
        
        // Initialize attributes
        p.vertex_position_prev_ = get_or_add_attribute<Vec3, Vertex>(m, "position_prev");
        p.vertex_position_prev_->copy(p.vertex_position_.get());
        
        p.vertex_velocity_ = get_or_add_attribute<Vec3, Vertex>(m, "velocity");
        p.vertex_velocity_->fill(Vec3(0.0, 0.0, 0.0));
        
        // Initialize mass attribute (default to 1.0)
        p.vertex_mass_ = get_or_add_attribute<Scalar, Vertex>(m, "mass");
        p.vertex_mass_->fill(1.0);
        
        // Initialize pin constraint for fixed vertices
        p.pin_constraint_ = std::make_shared<PinConstraint<MESH>>(1e10);
        p.pin_constraint_->init(m);
        
        // Reset solver
        p.solver_.reset();
        
        p.initialized_ = true;
        std::cout << "Simulation data initialized" << std::endl;
    }
    
    void load_fixed_vertices_from_scene(MESH& m, const io::SceneInfo& scene)
    {
        SimulationParameters& p = simulation_parameters_[&m];
        if (!p.initialized_)
        {
            std::cout << "Mesh not initialized yet" << std::endl;
            return;
        }
        
        // Find the model in the scene that corresponds to this mesh
        for (const auto& model : scene.models)
        {
            if (model.surface_mesh == &m)
            {
                // Reset all fixed flags
                p.vertex_fixed_->fill(false);
                
                // Set fixed vertices from the model
                for (uint32_t v_idx : model.fixed_vertices)
                {
                    if (v_idx < nb_vertices(m))
                    {
                        // Find vertex by index (actual implementation may differ)
                        foreach_cell(m, [&](Vertex v) -> bool {
                            if (index_of(m, v) == v_idx)
                            {
                                (*p.vertex_fixed_)[v_idx] = true;
                                return false; // Stop iteration
                            }
                            return true;
                        });
                    }
                }
                
                std::cout << "Fixed " << model.fixed_vertices.size() << " vertices from scene" << std::endl;
                break;
            }
        }
    }
    
    void add_constraint(MESH& m, std::unique_ptr<Constraint<MESH>> constraint)
    {
        SimulationParameters& p = simulation_parameters_[&m];
        if (!p.initialized_)
        {
            std::cout << "Mesh not initialized yet" << std::endl;
            return;
        }
        
        constraint->init(m);
        p.constraints_.push_back(std::move(constraint));
        
        // Reset system matrix since constraints changed
        p.solver_.reset();
    }
    
    void add_triangle_spring_constraint(MESH& m, Scalar weight = 1.0)
    {
        add_constraint(m, std::make_unique<TriangleSpringConstraint<MESH>>(weight));
    }
    
    void add_triangle_arap_constraint(MESH& m, Scalar weight = 1.0)
    {
        add_constraint(m, std::make_unique<TriangleARAPConstraint<MESH>>(weight));
    }
    
    void update_fixed_vertices_from_selection(MESH& m, CellsSet<MESH, Vertex>* selected_vertices)
    {
        SimulationParameters& p = simulation_parameters_[&m];
        if (!p.initialized_ || !selected_vertices)
            return;
        
        // Clear existing fixed vertices
        p.pin_constraint_->clear();
        
        // Add selected vertices as fixed
        p.pin_constraint_->add_from_cells_set(m, selected_vertices, p.vertex_position_.get());
        
        // Reset solver
        p.solver_.reset();
        
        std::cout << "Updated fixed vertices from selection: " << p.pin_constraint_->size() << " vertices fixed" << std::endl;
    }
    
    void clear_constraints(MESH& m)
    {
        SimulationParameters& p = simulation_parameters_[&m];
        p.constraints_.clear();
        p.solver_.reset();
    }
    
    void simulation_step(SimulationParameters& p)
    {
        if (!p.initialized_)
            return;
        
        // Store current positions as previous
        p.vertex_position_prev_->copy(p.vertex_position_.get());
        
        // Prepare all constraints for solver
        std::vector<std::reference_wrapper<std::unique_ptr<Constraint<MESH>>>> all_constraints;
        
        // First, add the pin constraint if it has fixed vertices
        if (p.pin_constraint_->size() > 0)
        {
            // Create a temporary unique_ptr that doesn't own the object
            auto pin_constraint_ptr = std::unique_ptr<Constraint<MESH>>(
                p.pin_constraint_.get(), [](Constraint<MESH>*) {}
            );
            all_constraints.push_back(std::ref(pin_constraint_ptr));
        }
        
        // Add the regular constraints
        for (auto& constraint : p.constraints_)
        {
            all_constraints.push_back(std::ref(constraint));
        }
        
        // Prepare a vector of raw unique_ptrs for the solver
        std::vector<std::unique_ptr<Constraint<MESH>>> solver_constraints;
        for (auto& constraint_ref : all_constraints)
        {
            solver_constraints.push_back(std::move(constraint_ref.get()));
        }
        
        // Use the solver to perform a simulation step
        p.solver_.step(
            *(p.mesh_),
            p.vertex_position_,
            p.vertex_position_prev_,
            p.vertex_velocity_,
            p.vertex_mass_,
            solver_constraints,
            p.time_step_,
            p.gravity_
        );
        
        // Restore the unique_ptrs
        for (size_t i = 0; i < all_constraints.size(); ++i)
        {
            all_constraints[i].get() = std::move(solver_constraints[i]);
        }
        
        // Update iteration count
        p.iteration_count_++;
    }
    
    void start_simulation(SimulationParameters& p)
    {
        if (!p.initialized_)
        {
            std::cout << "Mesh not initialized yet" << std::endl;
            return;
        }
        
        if (p.running_)
            return;
        
        p.running_ = true;
        p.stopping_ = false;
        p.iteration_count_ = 0;
        
        launch_thread([&]() {
            auto start = std::chrono::high_resolution_clock::now();
            
            while (true)
            {
                {
                    std::lock_guard<std::mutex> lock(p.mutex_);
                    simulation_step(p);
                }
                
                // Update rendering
                mesh_provider_->emit_attribute_changed(*(p.mesh_), p.vertex_position_.get());
                
                // Control simulation speed
                if (p.slow_down_)
                    std::this_thread::sleep_for(std::chrono::microseconds(1000000 / p.update_rate_));
                else
                    std::this_thread::yield();
                
                if (p.stopping_)
                {
                    p.stopping_ = false;
                    p.running_ = false;
                    break;
                }
            }
            
            auto end = std::chrono::high_resolution_clock::now();
            std::cout << "Simulation time: " << std::chrono::duration<Scalar>(end - start).count() << "s" << std::endl;
            std::cout << "Number of iterations: " << p.iteration_count_ << std::endl;
        });
        
        app_.start_timer(100, [&]() -> bool { return !p.running_; });
    }
    
    void stop_simulation(SimulationParameters& p)
    {
        p.stopping_ = true;
    }
    
    void reset_simulation(SimulationParameters& p)
    {
        if (p.running_)
            stop_simulation(p);
        
        // Wait for simulation to stop
        while (p.running_)
            std::this_thread::sleep_for(std::chrono::milliseconds(10));
            
        // Reset velocities
        p.vertex_velocity_->fill(Vec3(0.0, 0.0, 0.0));
        
        // Reset iteration count
        p.iteration_count_ = 0;
        
        // Update rendering
        mesh_provider_->emit_attribute_changed(*(p.mesh_), p.vertex_position_.get());
    }

protected:
    void init() override
    {
        mesh_provider_ = static_cast<ui::MeshProvider<MESH>*>(
            app_.module("MeshProvider (" + std::string{mesh_traits<MESH>::name} + ")"));
            
        // Try to get the SurfaceSelection module if available
        surface_selection_ = static_cast<ui::SurfaceSelection<MESH>*>(
            app_.module("SurfaceSelection (" + std::string{mesh_traits<MESH>::name} + ")"));
            
        timer_connection_ = boost::synapse::connect<ui::App::timer_tick>(&app_, [this]() {
            SimulationParameters& p = simulation_parameters_[selected_mesh_];
            
            if (!p.running_ && p.initialized_)
                mesh_provider_->emit_attribute_changed(*(p.mesh_), p.vertex_position_.get());
        });
    }
    
    void key_press_event(ui::View* view, int32 key_code) override
    {
        SimulationParameters& p = simulation_parameters_[selected_mesh_];
        
        if (!p.initialized_ || !selected_mesh_)
            return;
            
        if (key_code == GLFW_KEY_SPACE)
        {
            if (!p.running_)
                start_simulation(p);
            else
                stop_simulation(p);
        }
        else if (key_code == GLFW_KEY_R)
        {
            reset_simulation(p);
        }
        else if (key_code == GLFW_KEY_S)
        {
            if (!p.running_)
            {
                std::lock_guard<std::mutex> lock(p.mutex_);
                simulation_step(p);
                mesh_provider_->emit_attribute_changed(*(p.mesh_), p.vertex_position_.get());
            }
        }
    }
    
    void left_panel() override
    {
        ImGui::Text("Local-Global Projective Dynamics Solver");
        ImGui::Separator();
        
        imgui_mesh_selector(mesh_provider_, selected_mesh_, "Mesh", 
                           [&](MESH& m) { set_selected_mesh(m); });
                           
        if (selected_mesh_)
        {
            SimulationParameters& p = simulation_parameters_[selected_mesh_];
            
            imgui_combo_attribute<Vertex, Vec3>(
                *selected_mesh_, p.vertex_position_, "Position",
                [&](const std::shared_ptr<MeshAttribute<Vec3>>& attribute) { 
                    p.vertex_position_ = attribute;
                });
                
            if (p.vertex_position_ && !p.initialized_)
            {
                if (ImGui::Button("Initialize Simulation"))
                    init_simulation_data(*selected_mesh_);
            }
            
            if (p.initialized_)
            {
                ImGui::Separator();
                
                ImGui::Text("Physics Parameters");
                ImGui::InputFloat3("Gravity", &p.gravity_[0]);
                ImGui::InputFloat("Time Step", &p.time_step_, 0.001f, 0.01f);
                
                ImGui::Separator();
                
                ImGui::Text("Fixed Vertices");
                if (surface_selection_)
                {
                    MeshData<MESH>& md = mesh_provider_->mesh_data(*selected_mesh_);
                    CellsSet<MESH, Vertex>* selected_vertices_set = nullptr;
                    
                    imgui_combo_cells_set(md, selected_vertices_set, "Selection Set", 
                        [&](CellsSet<MESH, Vertex>* cs) {
                            selected_vertices_set = cs;
                        });
                    
                    if (selected_vertices_set)
                    {
                        if (ImGui::Button("Fix Selected Vertices"))
                        {
                            std::lock_guard<std::mutex> lock(p.mutex_);
                            update_fixed_vertices_from_selection(*selected_mesh_, selected_vertices_set);
                        }
                        
                        ImGui::Text("Fixed Vertices: %zu", p.pin_constraint_->size());
                        
                        if (ImGui::Button("Clear Fixed Vertices"))
                        {
                            std::lock_guard<std::mutex> lock(p.mutex_);
                            p.pin_constraint_->clear();
                            p.solver_.reset();
                        }
                    }
                    else
                    {
                        ImGui::TextColored(ImVec4(1,0.5,0,1), "Please create a vertex selection set");
                        ImGui::TextColored(ImVec4(1,0.5,0,1), "using the SurfaceSelection module");
                    }
                }
                else
                {
                    ImGui::TextColored(ImVec4(1,0,0,1), "SurfaceSelection module not available");
                }
                
                ImGui::Separator();
                
                ImGui::Text("Constraints");
                if (ImGui::Button("Add Spring Constraint"))
                {
                    add_triangle_spring_constraint(*selected_mesh_, 1.0);
                }
                ImGui::SameLine();
                if (ImGui::Button("Add ARAP Constraint"))
                {
                    add_triangle_arap_constraint(*selected_mesh_, 1.0);
                }
                
                if (ImGui::Button("Clear Constraints"))
                {
                    clear_constraints(*selected_mesh_);
                }
                
                ImGui::Text("Active Constraints: %zu", p.constraints_.size());
                
                ImGui::Separator();
                
                ImGui::Text("Simulation Controls");
                ImGui::Checkbox("Slow Down", &p.slow_down_);
                if (p.slow_down_)
                    ImGui::SliderInt("Update Rate", (int*)&p.update_rate_, 1, 100);
                    
                if (!p.running_)
                {
                    if (ImGui::Button("Start Simulation"))
                        start_simulation(p);
                        
                    ImGui::SameLine();
                    if (ImGui::Button("Single Step"))
                    {
                        std::lock_guard<std::mutex> lock(p.mutex_);
                        simulation_step(p);
                        mesh_provider_->emit_attribute_changed(*(p.mesh_), p.vertex_position_.get());
                    }
                }
                else
                {
                    if (ImGui::Button("Stop Simulation"))
                        stop_simulation(p);
                }
                
                if (ImGui::Button("Reset Simulation"))
                    reset_simulation(p);
                    
                ImGui::Separator();
                
                ImGui::Text("Statistics");
                ImGui::Text("Iteration Count: %u", p.iteration_count_);
                ImGui::Text("Vertices: %u", nb_vertices(*selected_mesh_));
                ImGui::Text("Edges: %u", nb_edges(*selected_mesh_));
                ImGui::Text("Faces: %u", nb_faces(*selected_mesh_));
                
                ImGui::Separator();
                
                ImGui::Text("Keyboard Controls:");
                ImGui::Text("SPACE - Start/Stop Simulation");
                ImGui::Text("S - Single Step");
                ImGui::Text("R - Reset Simulation");
            }
        }
    }
    
private:
    ui::MeshProvider<MESH>* mesh_provider_ = nullptr;
    ui::SurfaceSelection<MESH>* surface_selection_ = nullptr;
    MESH* selected_mesh_ = nullptr;
    
    std::unordered_map<MESH*, SimulationParameters> simulation_parameters_;
    std::unordered_map<MESH*, std::vector<std::shared_ptr<boost::synapse::connection>>> mesh_connections_;
    std::shared_ptr<boost::synapse::connection> timer_connection_;
};

} // namespace ui

} // namespace cgogn

#endif // CGOGN_SIMULATION_LOCAL_GLOBAL_H_

