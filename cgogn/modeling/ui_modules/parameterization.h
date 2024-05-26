#pragma once

#ifndef CGOGN_MODULE_COVERAGE_AXIS_H_
#define CGOGN_MODULE_COVERAGE_AXIS_H_
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
// import CGAL

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
	// Kernel for construct Delaunay
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

		Vec3 first_point_, second_point_, selected_point_;
		Vec4 cut_plane_;
		Vec3 plane_normal_;
		std::shared_ptr<SAttribute<Vec3>> surface_vertex_position_ = nullptr;
		std::shared_ptr<SAttribute<Vec3>> surface_edge_color_ = nullptr;
		std::vector<SEdge> boundary_edge;
	};

	void parameterization_init(SURFACE& surface)
	{
		SurfaceParameter& p = surface_parameters_[&surface];
		p.surface_= &surface;
		p.surface_vertex_position_ = get_or_add_attribute<Vec3, SVertex>(*p.surface_, "position");
		p.surface_edge_color_ = get_or_add_attribute<Vec3, SEdge>(*p.surface_, "surface_edge_color");
		parallel_foreach_cell(*p.surface_, [&](SEdge se) {
			value<Vec3>(*p.surface_, p.surface_edge_color_, se) = Vec3(0.0, 0.0, 0.0);
			return true;
		});
		p.initialized_ = true;
		


		
	}

	void cut_mesh(SurfaceParameter& p)
	{
		//Find the edges that intersected with the plane
		parallel_foreach_cell(*p.surface_, [&](SEdge se) {
			auto& vs = incident_vertices(*p.surface_, se);
			Vec3 pos1 = value<Vec3>(*p.surface_, p.surface_vertex_position_, vs[0]);
			Vec3 pos2 = value<Vec3>(*p.surface_, p.surface_vertex_position_, vs[1]);
			Scalar d1 = Vec4(pos1.x(), pos1.y(), pos1.z(), 1).dot(p.cut_plane_);
			Scalar d2 = Vec4(pos2.x(), pos2.y(), pos2.z(), 1).dot(p.cut_plane_);
			if((d1 < 0 && d2 > 0) || (d1 > 0 && d2 < 0))
			{
				std::cout << "find intersection" << std::endl;
				value<Vec3>(*p.surface_, p.surface_edge_color_, se) = Vec3(1.0, 0.0, 0.0);
				p.boundary_edge.push_back(se);
			}
			return true;
			});
	}

	void compute_plane_from_points(SurfaceParameter& p)
	{
		Scalar d = -p.plane_normal_.dot(p.first_point_);
		p.cut_plane_ = Vec4(p.plane_normal_.x(), p.plane_normal_.y(), p.plane_normal_.z(), d);
		std::cout << "plane defined: " << p.cut_plane_ <<std::endl;
	} 

protected:
	void init() override
	{
		surface_provider_ = static_cast<ui::MeshProvider<SURFACE>*>(
			app_.module("MeshProvider (" + std::string{mesh_traits<SURFACE>::name} + ")"));
	}
	void draw(View* view) override
	{
		
	}

	void mouse_press_event(View* view, int32 button, int32 x, int32 y) override
	{
		SurfaceParameter& p = surface_parameters_[selected_surface_mesh_];
		Vec3 first_pos_screen, second_pos_screen;
		if (button == GLFW_MOUSE_BUTTON_LEFT && p.activate_plane_selection_ )
		{
			
			rendering::GLVec3d near_point = view->unproject(x, y, 0.0);
			rendering::GLVec3d far_point = view->unproject(x, y, 1.0);
			
			if (!p.first_point_saved_)
			{
				first_pos_screen = Vec3(x, y, 0);
				p.first_point_ = near_point;
				p.first_point_saved_ = true;
				p.selected_point_ = near_point;
				p.point_selected_ = true;
				std::cout << "pick first point: " << p.first_point_.x() << ", " << p.first_point_.y() << ", "
						  << p.first_point_.z() << std::endl;

			}
			else
			{
				second_pos_screen = Vec3(x, y, 0);
				p.second_point_ = near_point;
				p.first_point_saved_ = false;
				p.point_selected_ = false;
				std::cout << "pick second point: " << p.second_point_.x() << ", " << p.second_point_.y() << ", "
						  << p.second_point_.z() << std::endl;
				/*Vec3 v = first_pos_screen - second_pos_screen;
				
				p.plane_normal_ = view->unproject(-v.y(), v.x(), 0);
				p.plane_normal_.normalize();*/
				p.plane_normal_= Vec3(0.5,0.5,1);
				compute_plane_from_points(p);
			}
		}
		

	}
	void key_press_event(View* view, int32 keycode) override
	{
		SurfaceParameter& p = surface_parameters_[selected_surface_mesh_];
		if (keycode == GLFW_KEY_C)
		{
			p.activate_plane_selection_ = true;
		}
		
	}
	void key_release_event(View* view, int32 keycode) override
	{
		SurfaceParameter& p = surface_parameters_[selected_surface_mesh_];
		if (keycode == GLFW_KEY_C)
		{
			p.activate_plane_selection_ = false;
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
			
			
			if (ImGui::Button("Parameterization init"))
			{
				parameterization_init(*selected_surface_mesh_);
			}
			
			if (p.initialized_)
			{
				if (ImGui::Button("Cut mesh")){
					cut_mesh(p);
				}
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