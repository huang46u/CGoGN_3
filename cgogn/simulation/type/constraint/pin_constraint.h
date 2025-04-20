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

#ifndef CGOGN_SIMULATION_CONSTRAINT_PIN_CONSTRAINT_H_
#define CGOGN_SIMULATION_CONSTRAINT_PIN_CONSTRAINT_H_

#include <cgogn/simulation/type/constraint/base_constraint.h>
#include <Eigen/Core>
#include <Eigen/Sparse>
#include <vector>
#include <unordered_set>

namespace cgogn
{

namespace simulation
{

/**
 * @brief Constraint for fixing vertices in place
 * 
 * This constraint fixes selected vertices to their initial positions.
 */
template <typename MESH>
class PinConstraint : public Constraint<MESH>
{
public:
    using Vertex = typename mesh_traits<MESH>::Vertex;
    using VertexAttribute = typename mesh_traits<MESH>::template Attribute<Vec3>;
    
    /**
     * @brief Constructor with weight parameter
     * 
     * @param weight Weight of the constraint in the global system (typically very large)
     */
    PinConstraint(Scalar weight = 1e10) : 
        Constraint<MESH>(weight), 
        initialized_(false)
    {}
    
    /**
     * @brief Initialize the constraint with a set of fixed vertices
     * 
     * @param mesh The mesh reference
     */
    void init(MESH& mesh) override
    {
        initialized_ = true;
    }
    
    /**
     * @brief Add a vertex to the fixed set
     * 
     * @param v The vertex to fix
     * @param position The initial position of the vertex
     */
    void add_fixed_vertex(const Vertex& v, const Vec3& position)
    {
        fixed_vertices_.insert(v);
        fixed_positions_[v] = position;
    }
    
    /**
     * @brief Add vertices from a cells set
     * 
     * @param mesh The mesh reference
     * @param cells_set The set of vertices to fix
     * @param position_attribute The attribute containing vertex positions
     */
    void add_from_cells_set(MESH& mesh, 
                           CellsSet<MESH, Vertex>* cells_set, 
                           const VertexAttribute* position_attribute)
    {
        if (!cells_set) return;
        
        cells_set->foreach_cell([&](Vertex v) {
            uint32 v_idx = index_of(mesh, v);
            add_fixed_vertex(v, (*position_attribute)[v_idx]);
        });
    }
    
    /**
     * @brief Clear all fixed vertices
     */
    void clear()
    {
        fixed_vertices_.clear();
        fixed_positions_.clear();
    }
    
    /**
     * @brief Get number of fixed vertices
     */
    size_t size() const
    {
        return fixed_vertices_.size();
    }
    
    /**
     * @brief Compute RHS contributions and add directly to the RHS vector
     * 
     * @param mesh The mesh reference
     * @param positions Current positions attribute
     * @param rhs The RHS vector to update directly
     */
    void project(MESH& mesh, 
                const VertexAttribute* positions,
                Eigen::VectorXd& rhs) override
    {
        if (!initialized_)
            return;
        
        // For each fixed vertex, add a constraint to maintain its initial position
        for (const Vertex& v : fixed_vertices_)
        {
            uint32 v_idx = index_of(mesh, v);
            const Vec3& target_pos = fixed_positions_[v];
            
            // Add the target position to the RHS with high weight - directly update using segment
            rhs.segment<3>(3 * v_idx) += this->weight_ * target_pos;
        }
    }
    
    /**
     * @brief Add the constraint contribution to the global system matrix
     * 
     * @param mesh The mesh reference
     * @param system_triplets The vector of triplets to accumulate the system matrix
     */
    void accumulation_matrix(MESH& mesh, 
                            std::vector<Eigen::Triplet<Scalar>>& system_triplets) override
    {
        if (!initialized_)
            return;
        
        // For each fixed vertex, add a high weight to the diagonal to fix its position
        for (const Vertex& v : fixed_vertices_)
        {
            uint32 v_idx = index_of(mesh, v);
            
            for (int d = 0; d < 3; ++d)
            {
                system_triplets.emplace_back(3 * v_idx + d, 3 * v_idx + d, this->weight_);
            }
        }
    }
    
private:
    bool initialized_;
    std::unordered_set<Vertex> fixed_vertices_;
    std::unordered_map<Vertex, Vec3> fixed_positions_;
};

} // namespace simulation

} // namespace cgogn

#endif // CGOGN_SIMULATION_CONSTRAINT_PIN_CONSTRAINT_H_
