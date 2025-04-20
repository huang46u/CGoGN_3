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

#ifndef CGOGN_SIMULATION_SOLVER_PROJECTIVE_DYNAMICS_SOLVER_H_
#define CGOGN_SIMULATION_SOLVER_PROJECTIVE_DYNAMICS_SOLVER_H_

#include <cgogn/core/types/mesh_traits.h>
#include <cgogn/geometry/types/vector_traits.h>
#include <cgogn/simulation/type/constraint/base_constraint.h>

#include <Eigen/Sparse>
#include <Eigen/SparseCholesky>
#include <vector>
#include <memory>

namespace cgogn
{

namespace simulation
{

using geometry::Vec3;
using geometry::Scalar;

/**
 * @brief Projective Dynamics solver for physical simulations
 * 
 * This solver implements the Projective Dynamics algorithm for physical simulations.
 * It handles both the local step (projection) and global step (linear system solving).
 * 
 * @tparam MESH The mesh type
 */
template <typename MESH>
class ProjectiveDynamicsSolver
{
public:
    using Vertex = typename mesh_traits<MESH>::Vertex;
    using VertexAttribute = typename mesh_traits<MESH>::template Attribute<Vec3>;
    
    /**
     * @brief Constructor
     */
    ProjectiveDynamicsSolver() : 
        system_matrix_initialized_(false)
    {}

    /**
     * @brief Destructor
     */
    ~ProjectiveDynamicsSolver() = default;

    /**
     * @brief Initialize the system matrix for efficient solving
     * 
     * @param mesh The mesh reference
     * @param vertex_position Current vertex positions attribute
     * @param vertex_mass Vertex mass attribute
     * @param vertex_fixed Fixed vertices flag attribute
     * @param constraints List of constraints
     * @param time_step Simulation time step (for inertia terms)
     */
    void init_system_matrix(
        MESH& mesh,
        const std::shared_ptr<VertexAttribute>& vertex_position,
        const std::shared_ptr<typename mesh_traits<MESH>::template Attribute<Scalar>>& vertex_mass,
        const std::shared_ptr<typename mesh_traits<MESH>::template Attribute<bool>>& vertex_fixed,
        const std::vector<std::unique_ptr<Constraint<MESH>>>& constraints,
        Scalar time_step)
    {
        uint32_t nb_vertices = nb_cells<Vertex>(m);
        uint32_t system_size = 3 * nb_vertices;
        
        // Prepare triplets for system matrix
        std::vector<Eigen::Triplet<Scalar>> system_triplets;
        
        // Add mass matrix contribution (M/h²)
        Scalar h2_inv = 1.0 / (time_step * time_step);
        
        foreach_cell(mesh, [&](Vertex v) -> bool {
            uint32_t v_idx = index_of(mesh, v);
            Scalar mass = (*vertex_mass)[v_idx];
            
            // Add diagonal entries for each coordinate (x, y, z)
            for (int d = 0; d < 3; ++d)
            {
                if ((*vertex_fixed)[v_idx])
                {
                    // For fixed vertices, use a large value to constrain position
                    system_triplets.emplace_back(3 * v_idx + d, 3 * v_idx + d, 1e10);
                }
                else
                {
                    // For normal vertices, add mass/h² term
                    system_triplets.emplace_back(3 * v_idx + d, 3 * v_idx + d, mass * h2_inv);
                }
            }
            return true;
        });
        
        // Add contributions from all constraints
        for (const auto& constraint : constraints)
        {
            constraint->accumulation_matrix(mesh, system_triplets);
        }
        
        // Build system matrix
        system_matrix_.resize(system_size, system_size);
        system_matrix_.setFromTriplets(system_triplets.begin(), system_triplets.end());
        
        // Initialize solver with system matrix
        solver_.analyzePattern(system_matrix_);
        solver_.factorize(system_matrix_);
        
        // Initialize RHS and solution vectors
        rhs_.resize(system_size);
        solution_.resize(system_size);
        
        system_matrix_initialized_ = true;
    }

    /**
     * @brief Perform one simulation step using Projective Dynamics
     * 
     * @param mesh The mesh reference
     * @param vertex_position Current vertex positions (will be updated)
     * @param vertex_position_prev Previous vertex positions
     * @param vertex_velocity Vertex velocities (will be updated)
     * @param vertex_mass Vertex mass attribute
     * @param vertex_fixed Fixed vertices flag attribute
     * @param constraints List of constraints
     * @param time_step Simulation time step
     * @param gravity Gravity vector
     */
    void step(
        MESH& mesh,
        const std::shared_ptr<VertexAttribute>& vertex_position,
        const std::shared_ptr<VertexAttribute>& vertex_position_prev,
        const std::shared_ptr<VertexAttribute>& vertex_velocity,
        const std::shared_ptr<typename mesh_traits<MESH>::template Attribute<Scalar>>& vertex_mass,
        const std::shared_ptr<typename mesh_traits<MESH>::template Attribute<bool>>& vertex_fixed,
        const std::vector<std::unique_ptr<Constraint<MESH>>>& constraints,
        Scalar time_step,
        const Vec3& gravity)
    {
        // Initialize system matrix if not already done
        if (!system_matrix_initialized_)
        {
            init_system_matrix(mesh, vertex_position, vertex_mass, vertex_fixed, constraints, time_step);
        }
        
        // Local step: compute projections
        local_step(mesh, vertex_position, vertex_position_prev, vertex_velocity, 
                  vertex_mass, vertex_fixed, constraints, time_step, gravity);
        
        // Global step: solve linear system
        global_step(mesh, vertex_position, vertex_position_prev, 
                   vertex_velocity, vertex_fixed, time_step);
    }

    /**
     * @brief Compute projections in the local step
     */
    void local_step(
        MESH& mesh,
        const std::shared_ptr<VertexAttribute>& vertex_position,
        const std::shared_ptr<VertexAttribute>& vertex_position_prev,
        const std::shared_ptr<VertexAttribute>& vertex_velocity,
        const std::shared_ptr<typename mesh_traits<MESH>::template Attribute<Scalar>>& vertex_mass,
        const std::vector<std::unique_ptr<Constraint<MESH>>>& constraints,
        Scalar time_step,
        const Vec3& gravity)
    {
        uint32_t nb_vertices = cgogn::nb_vertices(mesh);
        uint32_t system_size = 3 * nb_vertices;
        
        // 清空RHS向量
        rhs_.setZero(system_size);
        
        // 收集约束投影（本地步骤）- 直接更新RHS
        for (const auto& constraint : constraints)
        {
            constraint->project(mesh, vertex_position.get(), rhs_);
        }
        
        // 添加惯性项: M/h² * q̃
        Scalar h2_inv = 1.0 / (time_step * time_step);
        foreach_cell(mesh, [&](Vertex v) -> bool {
            uint32_t v_idx = index_of(mesh, v);
            
            Scalar mass = (*vertex_mass)[v_idx];
            const Vec3& pos_prev = (*vertex_position_prev)[v_idx];
            const Vec3& velocity = (*vertex_velocity)[v_idx];
            
            // 计算惯性位置: q̃ = qᵗ + h * vᵗ + h² * f_ext / m
            Vec3 sn = pos_prev + time_step * velocity + 
                               time_step * time_step * gravity / mass;
            
            // 直接使用segment更新RHS向量: M/h² * q̃
            rhs_.segment<3>(3 * v_idx) += mass * h2_inv * sn;
            
            return true;
        });
    }

    /**
     * @brief Solve the linear system in the global step
     * 
     * @param mesh The mesh reference
     * @param vertex_position Current vertex positions (will be updated)
     * @param vertex_position_prev Previous vertex positions
     * @param vertex_velocity Vertex velocities (will be updated)
     * @param vertex_fixed Fixed vertices flag attribute
     * @param time_step Simulation time step
     */
    void global_step(
        MESH& mesh,
        const std::shared_ptr<VertexAttribute>& vertex_position,
        const std::shared_ptr<VertexAttribute>& vertex_position_prev,
        const std::shared_ptr<VertexAttribute>& vertex_velocity,
        const std::shared_ptr<typename mesh_traits<MESH>::template Attribute<bool>>& vertex_fixed,
        Scalar time_step)
    {
        // Solve the linear system
        solution_ = solver_.solve(rhs_);
        
        // Update positions and velocities
        foreach_cell(mesh, [&](Vertex v) -> bool {
            uint32_t v_idx = index_of(mesh, v);
            
            if (!(*vertex_fixed)[v_idx])
            {
                // Update position
                Vec3& pos = (*vertex_position)[v_idx];
                for (int d = 0; d < 3; ++d) {
                    pos[d] = solution_(3 * v_idx + d);
                }
                
                // Update velocity
                Vec3& vel = (*vertex_velocity)[v_idx];
                const Vec3& pos_prev = (*vertex_position_prev)[v_idx];
                vel = (pos - pos_prev) / time_step;
            }
            
            return true;
        });
    }

    /**
     * @brief Check if the system matrix is initialized
     * 
     * @return true if initialized, false otherwise
     */
    bool is_initialized() const { 
        return system_matrix_initialized_; 
    }

    /**
     * @brief Reset the solver, clearing any previous initialization
     */
    void reset() {
        system_matrix_initialized_ = false;
        system_matrix_.resize(0, 0);
        rhs_.resize(0);
        solution_.resize(0);
    }

private:
    bool system_matrix_initialized_;
    Eigen::SparseMatrix<Scalar> system_matrix_;
    Eigen::VectorXd rhs_;
    Eigen::VectorXd solution_;
    Eigen::SimplicialLDLT<Eigen::SparseMatrix<Scalar>> solver_;　
};

} // namespace simulation

} // namespace cgogn

#endif // CGOGN_SIMULATION_SOLVER_PROJECTIVE_DYNAMICS_SOLVER_H_
