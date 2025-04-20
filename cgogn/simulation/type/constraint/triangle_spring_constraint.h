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

#ifndef CGOGN_SIMULATION_CONSTRAINT_TRIANGLE_SPRING_CONSTRAINT_H_
#define CGOGN_SIMULATION_CONSTRAINT_TRIANGLE_SPRING_CONSTRAINT_H_

#include <cgogn/simulation/type/constraint/base_constraint.h>
#include <Eigen/Core>
#include <Eigen/Sparse>
#include <cgogn/geometry/types/vector_traits.h>

namespace cgogn
{

namespace simulation
{

/**
 * @brief Spring constraint for triangular meshes
 * 
 * This constraint implements a spring-based energy for triangle mesh edges,
 * which tries to maintain the rest length of each edge.
 */
template <typename MESH>
class TriangleSpringConstraint : public Constraint<MESH>
{
public:
    using Vertex = typename mesh_traits<MESH>::Vertex;
    using Edge = typename mesh_traits<MESH>::Edge;
    using Face = typename mesh_traits<MESH>::Face;
    using VertexAttribute = typename mesh_traits<MESH>::template Attribute<Vec3>;
    using ScalarAttribute = typename mesh_traits<MESH>::template Attribute<Scalar>;

    /**
     * @brief Constructor with weight parameter
     * 
     * @param weight Weight of the constraint in the global system
     */
    TriangleSpringConstraint(Scalar weight = 1.0) : 
        Constraint<MESH>(weight), 
        initialized_(false)
    {}

    /**
     * @brief Initialize the constraint data from the m
     * 
     * @param m The mesh reference
     */
    void init(MESH& m) override
    {
        auto position = get_or_add_attribute<Vec3, Vertex>(m, "position");
        auto rest_length = get_or_add_attribute<Scalar, Edge>(m, "rest_length");
        

        parallel_foreach_cell(m, [&](Edge e) -> bool {
            uint32 e_idx = index_of(m, e);
            
            auto& vertices = incident_vertices(m, e);
            uint32 v1 = vertices[0];
            uint32 v2 = vertices[1];
            auto edge_pair = std::make_pair(std::min(v1, v2), std::max(v1, v2));
            
            Vec3 p1 = (*position)[index_of(m, v1)];
            Vec3 p2 = (*position)[index_of(m, v2)];
            
            value<Scalar>(m, rest_length, e) = (p2 - p1).norm();
            
            return true;
        });

        initialized_ = true;
        
    }

    /**
     * @brief Compute RHS matrix elements to satisfy the constraint
     * 
     * @param m The mesh reference
     * @param positions Current positions attribute
     * @param RHS RHS matrix
     */
    void project(MESH& m, 
                const VertexAttribute* positions,
                Eigen::Matrix3d RHS) override
    {
        
        auto rest_length = get_attribute<Scalar, Edge>(m, "rest_length");
        
        foreach_cell(m, [&](Edge e) -> bool {

            uint32 e_idx = index_of(m, e);
            auto& vertices = incident_vertices(m, e);
            uint32 v1_idx = index_of(m, vertices[0]);
            uint32 v2_idx = index_of(m, vertices[1]);
            Vec3 p1 = (*positions)[v1_idx];
            Vec3 p2 = (*positions)[v2_idx];
            
            Scalar rest_len = (*rest_length)[e_idx];
        
            Vec3 direction = (p2 - p1).normalized();
            
            rhs.segment<3>(3 *v1_idx) += - weight_ * 0.5f * direction * rest_len ;
            rhs.segment<3>(3 *v2_idx) += weight_ * 0.5f * direction * rest_len ;
            
            return true;
        });
    }

    /**
     * @brief Add the constraint contribution to the global system matrix
     * 
     * @param m The mesh reference
     * @param system_triplets The vector of triplets to accumulate the system matrix
     */
    void accumulation_matrix(MESH& m, 
                            std::vector<Eigen::Triplet<Scalar>>& system_triplets) override
    {

        foreach_cell(m, [&](Edge e) -> bool {
          
            auto& vertices = incident_vertices(m, e);
            uint32 v1_ind = index_of(m,vertices[0]);
            uint32 v2_ind = index_of(m,vertices[1]);
           
            for (int i = 0; i < 3; ++i) {
                system_triplets.emplace_back(3 * v1_ind + i, 3 * v1_ind + i, weight_ * 0.5);
                system_triplets.emplace_back(3 * v1_ind + i, 3 * v2_ind + i, -weight_ * 0.5);
                system_triplets.emplace_back(3 * v2_ind + i, 3 * v1_ind + i, -weight_ * 0.5);
                system_triplets.emplace_back(3 * v2_ind + i, 3 * v2_ind + i, weight_ * 0.5);
            }
            
            return true;
        });
    }

};

} // namespace simulation

} // namespace cgogn

#endif // CGOGN_SIMULATION_CONSTRAINT_TRIANGLE_SPRING_CONSTRAINT_H_
