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

#ifndef CGOGN_SIMULATION_CONSTRAINTS_CONSTRAINT_H_
#define CGOGN_SIMULATION_CONSTRAINTS_CONSTRAINT_H_

#include <cgogn/core/types/mesh_traits.h>
#include <cgogn/geometry/types/vector_traits.h>
#include <cgogn/geometry/types/eigen.h>

namespace cgogn
{

namespace simulation
{

using geometry::Vec3;
using geometry::Scalar;

/**
 * @brief Base class for all constraints in the projective dynamics framework
 * 
 * @tparam MESH The mesh type
 */
template <typename MESH>
class Constraint
{
public:
   
    
    /**
     * @brief Constructor with weight parameter
     * 
     * @param weight The weight of this constraint in the global system
     */
    Constraint(Scalar weight = 1.0) : weight_(weight) {}
    
    /**
     * @brief Virtual destructor
     */
    virtual ~Constraint() {}

    /**
     * @brief Get the weight of this constraint
     * 
     * @return Scalar The weight value
     */
    Scalar weight() const { return weight_; }
    
    /**
     * @brief Set the weight of this constraint
     * 
     * @param weight The new weight value
     */
    void set_weight(Scalar weight) { weight_ = weight; }

    /**
     * @brief Initialize the constraint data
     * 
     * @param m The mesh reference
     */
    virtual void init(MESH& m) = 0;
    
    /**
     * @brief Compute RHS matrix elements to satisfy the constraint
     * 
     * @param m The mesh reference
     * @param positions Current positions attribute
     * @param rhs_triplets Output vector of triplets for the RHS matrix
     */
    virtual void project(MESH& m, 
                        const typename mesh_traits<MESH>::template Attribute<Vec3>* positions, 
                        std::vector<geometry::Triplet>& rhs_triplets) = 0;
    
 
    /**
     * @brief Add the constraint contribution to the global system matrix using triplets
     * 
     * @param m The mesh reference
     * @param system_triplets The vector of triplets to accumulate the system matrix
     */
    virtual void accumulation_matrix(MESH& m, 
                                std::vector<geometry::Triplet>& system_triplets) = 0;
    

protected:
    Scalar weight_;
};

} // namespace simulation

} // namespace cgogn

#endif // CGOGN_SIMULATION_CONSTRAINTS_CONSTRAINT_H_