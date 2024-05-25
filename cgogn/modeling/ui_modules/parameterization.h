#pragma once

#ifndef CGOGN_MODULE_COVERAGE_AXIS_H_
#define CGOGN_MODULE_COVERAGE_AXIS_H_
#include <cgogn/core/types/maps/gmap/gmap_base.h>
#include <cgogn/core/ui_modules/mesh_provider.h>
#include <cgogn/ui/app.h>
#include <cgogn/ui/module.h>
#include <cmath>
#include <filesystem>
#include <libacc/bvh_trees_spheres.h>
#include <libacc/bvh_tree.h>
#include <libacc/kd_tree.h>
// import CGAL

#include <GLFW/glfw3.h>


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
	template <typename T>
	using SVertex = typename mesh_traits<SURFACE>::Vertex;
	using SEdge = typename mesh_traits<SURFACE>::Edge;
	using SFace = typename mesh_traits<SURFACE>::Face;

	using Vec4 = geometry::Vec4;
	using Vec3 = geometry::Vec3;
	using Scalar = geometry::Scalar;

public:
	Parameterization(const App& app) : ViewModule(app, "CoverageAxis")
	{
	}
	~Parameterization()
	{
	}


public:
	struct SufaceParameter
	{
		
	};


protected:
	void init() override
	{
		surface_provider_ = static_cast<ui::MeshProvider<SURFACE>*>(
			app_.module("MeshProvider (" + std::string{mesh_traits<SURFACE>::name} + ")"));
	}
	void draw(View* view) override
	{
	}
	void left_panel() override
	{
		imgui_mesh_selector(surface_provider_, selected_surface_mesh_, "Surface", [&](SURFACE& m) {
			selected_surface_mesh_ = &m;
			surface_provider_->mesh_data(m).outlined_until_ = App::frame_time_ + 1.0;
		});

		ImGui::Separator();
		if (selected_surface_mesh_)
		{
		}
	}

private:
	SURFACE* selected_surface_mesh_ = nullptr;

	std::unordered_map<const SURFACE*, SufaceParameter> surface_parameters_;

	MeshProvider<SURFACE>* surface_provider_;
};

} // namespace ui

} // namespace cgogn

#endif // CGOGN_MODULE_COVERAGE_AXIS_H_