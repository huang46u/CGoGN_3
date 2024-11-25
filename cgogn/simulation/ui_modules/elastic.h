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

#ifndef CGOGN_MODULE_ELASTIC_H_
#define CGOGN_MODULE_ELASTIC_H_

#include <cgogn/ui/app.h>
#include <cgogn/ui/module.h>

#include <boost/synapse/connect.hpp>

namespace cgogn
{

namespace ui
{

using geometry::Scalar;
using geometry::Vec3;

template <typename MESH>
class Elastic : public Module
{
	static_assert(mesh_traits<MESH>::dimension == 3, "Elastic can only be used with meshes of dimension 3");

	template <typename T>
	using Attribute = typename mesh_traits<MESH>::template Attribute<T>;

	using Vertex = typename mesh_traits<MESH>::Vertex;
	using Edge = typename mesh_traits<MESH>::Edge;
	using Face = typename mesh_traits<MESH>::Face;

public:
	Elastic(const App& app) : Module(app, "Elastic (" + std::string{mesh_traits<MESH>::name} + ")")
	{
	}
	~Elastic()
	{
	}

protected:
	void init() override
	{
		
	}

	void left_panel() override
	{
		
	}

private:
	MeshProvider<MESH>* mesh_provider_ = nullptr;
	MESH* domain_ = nullptr;

	bool domain_initialized_ = false;
	bool running_ = false;

	std::shared_ptr<Attribute<Vec3>> vertex_water_position_;
	std::shared_ptr<Attribute<Vec3>> vertex_water_flux_;


	CellsSet<MESH, Face>* selected_faces_set_ = nullptr;
	CellsSet<MESH, Edge>* selected_edges_set_ = nullptr;

	std::vector<std::shared_ptr<boost::synapse::connection>> domain_connections_;
	std::shared_ptr<boost::synapse::connection> timer_connection_;
};

} // namespace ui

} // namespace cgogn

#endif // CGOGN_MODULE_ELASTIC_H_
