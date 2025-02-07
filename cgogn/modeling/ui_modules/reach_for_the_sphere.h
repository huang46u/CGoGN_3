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

#ifndef CGOGN_MODULE_SKELETON_EXTRACTOR_H_
#define CGOGN_MODULE_SKELETON_EXTRACTOR_H_

#include <cgogn/core/ui_modules/mesh_provider.h>
#include <cgogn/ui/app.h>
#include <cgogn/ui/module.h>

#include <cgogn/core/functions/convert.h>

#include <cgogn/core/functions/traversals/face.h>
#include <cgogn/core/functions/traversals/vertex.h>

#include <cgogn/geometry/types/vector_traits.h>

#include <cgogn/geometry/algos/length.h>
#include <cgogn/geometry/algos/medial_axis.h>

#include <cgogn/modeling/algos/remeshing/pliant_remeshing.h>
#include <cgogn/modeling/algos/skeleton.h>

#include <boost/synapse/connect.hpp>

#include <unordered_map>

namespace cgogn
{

namespace ui
{

using geometry::Mat3;
using geometry::Scalar;
using geometry::Vec3;

template <typename POINT, typename SURFACE>
class ReachForTheSphere : public Module
{
	template <typename T>
	using PointAttribute = typename mesh_traits<POINT>::template Attribute<T>;
	template <typename T>
	using SurfaceAttribute = typename mesh_traits<SURFACE>::template Attribute<T>;

	using PointVertex = typename mesh_traits<POINT>::Vertex;

	using SurfaceVertex = typename mesh_traits<SURFACE>::Vertex;
	using SurfaceEdge = typename mesh_traits<SURFACE>::Edge;
	using SurfaceFace = typename mesh_traits<SURFACE>::Face;


	struct SurfaceParameters
	{
		SURFACE* mesh_;
		std::shared_ptr<SurfaceAttribute<Vec3>> surface_vertex_positions = nullptr;
		
		POINT* sdf_;
		std::shared_ptr<PointAttribute<Vec3>> point_positions = nullptr;
		std::shared_ptr<PointAttribute<Scalar>> point_sdf = nullptr;

	};

public:
	ReachForTheSphere(const App& app)
		: Module(app, "ReachForTheSphere"),
		  selected_surface_(nullptr), selected_surface_vertex_position_(nullptr)
	{
	}

	~ReachForTheSphere()
	{
	}

private:
	void init_surface_mesh(SURFACE* s)
	{
		SurfaceParameters& p = surface_parameters_[s];
		p.mesh_ = s;
		MeshData<SURFACE>& md = surface_provider_->mesh_data(*s);
	}

public:
	
protected:
	void init() override
	{
		
		surface_provider_ = static_cast<ui::MeshProvider<SURFACE>*>(
			app_.module("MeshProvider (" + std::string{mesh_traits<SURFACE>::name} + ")"));
		point_provider_ = static_cast<ui::MeshProvider<POINT>*>(
			app_.module("PointProvider (" + std::string{mesh_traits<SURFACE>::name} + ")"));

		surface_provider_->foreach_mesh([this](SURFACE& m, const std::string&) { init_surface_mesh(&m); });


		connections_.push_back(boost::synapse::connect<typename MeshProvider<SURFACE>::mesh_added>(
			surface_provider_, this, &ReachForTheSphere<POINT, SURFACE>::init_surface_mesh));
	}

	void left_panel() override
	{
		imgui_mesh_selector(surface_provider_, selected_surface_, "Surface", [&](SURFACE& m) {
			selected_surface_ = &m;
			surface_provider_->mesh_data(m).outlined_until_ = App::frame_time_ + 1.0;
		});

		if (selected_surface_)
		{
			imgui_combo_attribute<SurfaceVertex, Vec3>(*selected_surface_, selected_surface_vertex_position_,
													   "Position",
													   [&](const std::shared_ptr<SurfaceAttribute<Vec3>>& attribute) {
														   selected_surface_vertex_position_ = attribute;
													   });

			if (selected_surface_vertex_position_)
			{
				static float wL = 0.6f, wH = 0.1f, wM = 0.1f, resampling_ratio = 0.9f;
				ImGui::SliderFloat("Smoothness", &wL, 0.01f, 1.0f);
				ImGui::SliderFloat("Velocity", &wH, 0.01f, 1.0f);
				ImGui::SliderFloat("Medial attraction", &wM, 0.01f, 1.0f);
				ImGui::SliderFloat("Resampling ratio", &resampling_ratio, 0.01f, 2.0f);

			
			}
			
		}
	}

private:
	MeshProvider<POINT>* point_provider_ = nullptr;
	MeshProvider<SURFACE>* surface_provider_ = nullptr;

	SURFACE* selected_surface_ = nullptr;
	std::shared_ptr<SurfaceAttribute<Vec3>> selected_surface_vertex_position_ = nullptr;

	std::unordered_map<const SURFACE*, SurfaceParameters> surface_parameters_;

	std::vector<std::shared_ptr<boost::synapse::connection>> connections_;
};

} // namespace ui

} // namespace cgogn

#endif // CGOGN_MODULE_SKELETON_EXTRACTOR_H_
