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

#ifndef CGOGN_MODULE_PARAMETERIZATION_H_
#define CGOGN_MODULE_PARAMETERIZATION_H_
#include <cgogn/core/types/maps/gmap/gmap_base.h>
#include <cgogn/core/ui_modules/mesh_provider.h>
#include <cgogn/ui/app.h>
#include <cgogn/ui/module.h>
#include <cgogn/ui/imgui_helpers.h>
#include <cmath>
#include <filesystem>
#include <libacc/bvh_trees_spheres.h>
#include <libacc/bvh_tree.h>
#include <libacc/kd_tree.h>
#include <cgogn/rendering/frame_manipulator.h>

#include <GLFW/glfw3.h>
#include <GL/gl3w.h>

#include <iomanip>
#include <limits>

namespace cgogn
{

namespace ui
{
template <typename SURFACE>
class Parameterization : public ViewModule
{
	static_assert(mesh_traits<SURFACE>::dimension >= 2, "Coverage axis can only be used with meshes of dimension >= 2");
	template <typename T>
	using SAttribute = typename mesh_traits<SURFACE>::template Attribute<T>;

	using SVertex = typename mesh_traits<SURFACE>::Vertex;
	using SEdge = typename mesh_traits<SURFACE>::Edge;
	using SFace = typename mesh_traits<SURFACE>::Face;

	using Vec4 = geometry::Vec4;
	using Vec3 = geometry::Vec3;
	using Scalar = geometry::Scalar;

public:
	Parameterization(const App& app) : ViewModule(app, "Parameterization")
	{
	}
	~Parameterization()
	{
	}


public:
	struct SurfaceParameter
	{
		SURFACE* surface_;
		bool initialized_ = false;
		bool first_point_saved_ = false;
		bool activate_plane_selection_ = false;
		bool point_selected_ = false;

		rendering::GLVec3d first_point_, second_point_, selected_point_;
		rendering::GLVec3d plane_normal_;
		std::shared_ptr<SAttribute<Vec3>> surface_vertex_position_ = nullptr;
		std::shared_ptr<SAttribute<Vec3>> surface_edge_color_ = nullptr;
		std::shared_ptr<SAttribute<Vec3>> surface_face_color_ = nullptr;
		std::shared_ptr<SAttribute<Vec3>> surface_vertex_color_ = nullptr;
		std::shared_ptr<SAttribute<Vec3>> surface_tc_position_ = nullptr;
		std::vector<SEdge> boundary_edge;


		rendering::FrameManipulator frame_manipulator_;
		bool show_frame_manipulator_;
		bool manipulating_frame_;
		bool clipping_plane_;
	};

	void parameterization_init(SURFACE& surface)
	{
		MeshData<SURFACE>& md = surface_provider_->mesh_data(surface);
		SurfaceParameter& p = surface_parameters_[&surface];
		p.surface_= &surface;
		p.surface_vertex_position_ = get_or_add_attribute<Vec3, SVertex>(*p.surface_, "position");
		p.surface_edge_color_ = get_or_add_attribute<Vec3, SEdge>(*p.surface_, "surface_edge_color");
		p.surface_face_color_ = get_or_add_attribute<Vec3, SFace>(*p.surface_, "surface_face_color");
		p.surface_vertex_color_ = get_or_add_attribute<Vec3, SVertex>(*p.surface_, "surface_vert_color");
		p.surface_tc_position_ = get_or_add_attribute<Vec3, SVertex>(*p.surface_, "tc_position");
		
		foreach_cell(*p.surface_, [&](SEdge v) {
			value<Vec3>(*p.surface_, p.surface_edge_color_, v) = Vec3(0.0, 0.0, 0.0);
			return true;
		});
		p.initialized_ = true;

		Scalar size = (md.bb_max_ - md.bb_min_).norm() / 25;
		Vec3 position = 0.2 * md.bb_min_ + 0.8 * md.bb_max_;
		
		p.frame_manipulator_.set_size(size);
		p.frame_manipulator_.set_position(position);
		
	}


	//Only applicable for disk-like surface
	void tutte_embedding(SurfaceParameter& p)
	{
		//Identify bord vertices
		std::vector<Edge> boundary_edges;
		foreach_cell(*p.surface_, [&](SEdge e) {
			if (incident_faces(*p.surface_, e).size()==1)
			{
				
				boundary_edges.push_back(e);
				
			}
		
			return true;
		});
		for (auto e : boundary_edges)
		{
			auto& vs = incident_vertices(*p.surface_, e);
			for (auto v : vs)
			{
				value<Vec3>(*p.surface_, p.surface_vertex_color_, v) = Vec3(1.0, 0.0, 0.0);
			}
			value<Vec3>(*p.surface_, p.surface_edge_color_, e) = Vec3(1.0, 0.0, 0.0);
		}
		surface_provider_->emit_attribute_changed(*p.surface_, p.surface_vertex_color_.get());
		surface_provider_->emit_attribute_changed(*p.surface_, p.surface_edge_color_.get());
		surface_provider_->emit_attribute_changed(*p.surface_, p.surface_face_color_.get());
		//Fix the boundary vertices on a square of size 1

	}

	void cut_mesh(SurfaceParameter& p)
	{
		
		foreach_cell(*p.surface_, [&](SEdge se) {
			value<Vec3>(*p.surface_, p.surface_edge_color_, se) = Vec3(0.0, 0.0, 0.0);
			return true;
		});
		foreach_cell(*p.surface_, [&](SVertex se) {
			value<Vec3>(*p.surface_, p.surface_vertex_color_, se) = Vec3(0.0, 0.0, 0.0);
			return true;
		});
		foreach_cell(*p.surface_, [&](SFace se) {
			value<Vec3>(*p.surface_, p.surface_face_color_, se) = Vec3(0.0, 0.0, 0.0);
			return true;
		});
		MeshData<SURFACE>& md = surface_provider_->mesh_data(*p.surface_);
		Vec3 position;
		p.frame_manipulator_.get_position(position);
		Vec3 axis_z;
		p.frame_manipulator_.get_axis(rendering::FrameManipulator::Zt, axis_z);
		Scalar d = -(position.dot(axis_z));
		Vec4 plane = Vec4(axis_z.x(), axis_z.y(), axis_z.z(), d);
		//Find the edges that intersected with the plane
		foreach_cell(*p.surface_, [&](SEdge se) {
			auto& vs = incident_vertices(*p.surface_, se);
			auto& fs = incident_faces(*p.surface_, se);
			Vec3 pos1 = value<Vec3>(*p.surface_, p.surface_vertex_position_, vs[0]);
			Vec3 pos2 = value<Vec3>(*p.surface_, p.surface_vertex_position_, vs[1]);
			Scalar d1 = Vec4(pos1.x(), pos1.y(), pos1.z(), 1).dot(plane);
			Scalar d2 = Vec4(pos2.x(), pos2.y(), pos2.z(), 1).dot(plane);

			if (std::signbit(d1) != std::signbit(d2))
			{
				for (SVertex& v : vs)
				{
					value<Vec3>(*p.surface_, p.surface_vertex_color_, v) = Vec3(1.0, 0.0, 0.0);
				}
				for (SFace& f : fs)
				{
					value<Vec3>(*p.surface_, p.surface_face_color_, f) = Vec3(1.0, 0.0, 0.0);
				}
				value<Vec3>(*p.surface_, p.surface_edge_color_, se) = Vec3(1.0, 0.0, 0.0);
			}

			return true;
		}); 
		/*foreach_cell(*p.surface_, [&](SVertex se) {
			Vec3 pos = value<Vec3>(*p.surface_, p.surface_vertex_position_, se);
			Scalar d1 = Vec4(pos.x(), pos.y(), pos.z(), 1).dot(plane);
			if (std::signbit(d1))
			{
				auto& es = incident_edges(*p.surface_, se);
				for (SEdge& e : es)
				{
					value<Vec3>(*p.surface_, p.surface_edge_color_, e) = Vec3(1.0, 0.0, 0.0);
				}
			
			}

			return true;
		}); */
		surface_provider_->emit_attribute_changed(*p.surface_, p.surface_vertex_color_.get());
		surface_provider_->emit_attribute_changed(*p.surface_, p.surface_edge_color_.get());
		surface_provider_->emit_attribute_changed(*p.surface_, p.surface_face_color_.get());
	}

	void compute_plane_from_points(SurfaceParameter& p)
	{
		rendering::GLVec3d a = p.second_point_ - p.first_point_;
		a.normalize();
		rendering::GLVec3d b = p.plane_normal_.cross(a);
		b.normalize();
		rendering::GLVec3 X = a.cast<float>(); // rendering::GLVec3(a.x(), a.y(), a.z());
		rendering::GLVec3 Y = b.cast<float>(); // rendering::GLVec3(b.x(), b.y(), b.z());
		p.frame_manipulator_.set_orientation(X,Y);
		p.frame_manipulator_.set_position(p.first_point_);
	} 

protected:
	void init() override
	{
		surface_provider_ = static_cast<ui::MeshProvider<SURFACE>*>(
			app_.module("MeshProvider (" + std::string{mesh_traits<SURFACE>::name} + ")"));
	}
	void draw(View* view) override
	{
		SurfaceParameter& p = surface_parameters_[selected_surface_mesh_];
		auto& proj_matrix = view->projection_matrix();
		auto& view_matrix = view->modelview_matrix();
		if (p.show_frame_manipulator_)
			p.frame_manipulator_.draw(true, true, proj_matrix, view_matrix);
	}

	void mouse_press_event(View* view, int32 button, int32 x, int32 y) override
	{
		SurfaceParameter& p = surface_parameters_[selected_surface_mesh_];
		rendering::GLVec3d first_pos_screen, second_pos_screen;
		if (button == GLFW_MOUSE_BUTTON_LEFT && p.activate_plane_selection_ )
		{
			rendering::GLVec3d near_point = view->unproject(x, y, 0.25);
			rendering::GLVec3d far_point = view->unproject(x, y, 1.0);
			
			if (!p.first_point_saved_)
			{
				first_pos_screen = rendering::GLVec3d(x, y, 0);
				p.first_point_ = near_point;
				p.first_point_saved_ = true;
				p.selected_point_ = near_point;
				p.point_selected_ = true;
				std::cout << "pick first point: " << p.first_point_.x() << ", " << p.first_point_.y() << ", "
						  << p.first_point_.z() << std::endl;

			}
			else
			{
				second_pos_screen = rendering::GLVec3d(x, y, 0);
				p.second_point_ = near_point;
				p.first_point_saved_ = false;
				p.point_selected_ = false;
				std::cout << "pick second point: " << p.second_point_.x() << ", " << p.second_point_.y() << ", "
						  << p.second_point_.z() << std::endl;
				rendering::GLVec3d v =  second_pos_screen-first_pos_screen;
				p.plane_normal_ = view->unproject(-v.y(), v.x(), 0.0) - view->unproject(0.0,0.0,0.0);
				p.plane_normal_.normalize();
				
				compute_plane_from_points(p);
			}
		}
		if (p.manipulating_frame_)
		{
			auto [P, Q] = view->pixel_ray(x, y);
			p.frame_manipulator_.pick(x, y, P, Q);
			view->request_update();
		}
	}
	void mouse_release_event(View* view, int32, int32, int32) override
	{
		SurfaceParameter& p = surface_parameters_[selected_surface_mesh_];
		p.frame_manipulator_.release();
		view->request_update();
		
	}
	void mouse_move_event(View* view, int32 x, int32 y) override
	{
		SurfaceParameter& p = surface_parameters_[selected_surface_mesh_];
		bool leftpress = view->mouse_button_pressed(GLFW_MOUSE_BUTTON_LEFT);
		bool rightpress = view->mouse_button_pressed(GLFW_MOUSE_BUTTON_RIGHT);
		if (p.manipulating_frame_ && (rightpress || leftpress))
		{
			p.frame_manipulator_.drag(leftpress, x, y);
			if (p.clipping_plane_)
			{
				cut_mesh(p);
			}
			view->stop_event();
			view->request_update();
		}
	}


	void key_press_event(View* view, int32 keycode) override
	{
		SurfaceParameter& p = surface_parameters_[selected_surface_mesh_];
		if (keycode == GLFW_KEY_C)
		{
			p.activate_plane_selection_ = true;
			
		}
		if (keycode == GLFW_KEY_D)
		{

			if (p.show_frame_manipulator_)
				p.manipulating_frame_ = true;
		}
	}
	void key_release_event(View* view, int32 keycode) override
	{
		SurfaceParameter& p = surface_parameters_[selected_surface_mesh_];
		if (keycode == GLFW_KEY_C)
		{
			p.activate_plane_selection_ = false;
			
		}
		if (keycode == GLFW_KEY_D)
		{
			p.manipulating_frame_ = false;
		}
	}

	void left_panel() override
	{
		imgui_mesh_selector(surface_provider_, selected_surface_mesh_, "Surface", [&](SURFACE& m) {
			selected_surface_mesh_ = &m;
			surface_provider_->mesh_data(m).outlined_until_ = App::frame_time_ + 1.0;
		});
		
		ImGui::Separator();
		if (selected_surface_mesh_){
		
			SurfaceParameter& p = surface_parameters_[selected_surface_mesh_];
			imgui_combo_attribute<SVertex, Vec3>(*selected_surface_mesh_, p.surface_vertex_position_, "Position",
													   [&](const std::shared_ptr<SAttribute<Vec3>>& attribute) {
														   p.surface_vertex_position_ = attribute;
													   });
			
			if (ImGui::Button("Parameterization init"))
			{
				parameterization_init(*selected_surface_mesh_);
			}
			
			if (p.initialized_)
			{
				if (ImGui::Button("Cut mesh")){
					cut_mesh(p);
				}
				if (ImGui::Button("Tutte embedding"))
				{
					tutte_embedding(p);
				}
				if (ImGui::Checkbox("Clipping plane ", &p.clipping_plane_))
				{
					for (View* v : linked_views_)
						v->request_update();
				}
				if (ImGui::Checkbox("Show clipping plane", &p.show_frame_manipulator_))
				{
					for (View* v : linked_views_)
						v->request_update();
				}
				if (p.show_frame_manipulator_)
					ImGui::TextUnformatted("Press D to manipulate the plane");
			}
		}
	}

private:
	SURFACE* selected_surface_mesh_ = nullptr;

	std::unordered_map<const SURFACE*, SurfaceParameter> surface_parameters_;

	MeshProvider<SURFACE>* surface_provider_;
};

} // namespace ui

} // namespace cgogn

#endif // CGOGN_MODULE_COVERAGE_AXIS_H_