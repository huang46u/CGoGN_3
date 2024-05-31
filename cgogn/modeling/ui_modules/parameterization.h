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

#include <libacc/primitives.h>
#include <cgogn/rendering/frame_manipulator.h>
#include <Eigen/Sparse>
#include <cgogn/geometry/algos/angle.h>


#include <GLFW/glfw3.h>
#include <GL/gl3w.h>

#include <iomanip>
#include <limits>

namespace cgogn
{

namespace ui
{

using Vec4 = geometry::Vec4;
using Vec3 = geometry::Vec3;
using Vec2 = geometry::Vec2;
using Scalar = geometry::Scalar;

enum boundary_form
{
	BOUNDARY_FORM_CIRCLE,
	BOUNDARY_FORM_SQUARE
};

enum parameterization_weight
{
	WEIGHT_UNIFORM,
	WEIGHT_FLOATER
};

void fix_boundary_vertices(CMap2& m, boundary_form form,
						   typename mesh_traits<CMap2>::template Attribute<bool>* boundary,
						   typename mesh_traits<CMap2>::template Attribute<Vec2>* tc_position)
{
	CMap2::Edge start_edge;
	foreach_cell(m, [&](CMap2::Edge e) {
		if (is_incident_to_boundary(m, e))
		{
			start_edge = e;
			return false;
		}
		return true;
		});
	Dart b = is_boundary(m, start_edge.dart_) ? start_edge.dart_ : phi2(m, start_edge.dart_);
	CMap2::Face bf(b);
	std::vector<CMap2::Vertex> boundary_vertices;
	foreach_incident_vertex(m, bf, [&](CMap2::Vertex bv) {
		value<bool>(m, boundary, bv) = true;
		boundary_vertices.push_back(bv);
		return true;
	});
	switch (form)
	{
	case cgogn::ui::BOUNDARY_FORM_CIRCLE: {
		Scalar step = 2 * M_PI / boundary_vertices.size();
		for (uint32 idx = 0; idx < boundary_vertices.size(); idx++)
		{
			Scalar x = std::cos(idx * step)+1;
			Scalar y = std::sin(idx * step)+1;
			value<Vec2>(m, tc_position, boundary_vertices[idx]) = Vec2(x, y);
		}
	}
		break;
	case cgogn::ui::BOUNDARY_FORM_SQUARE: {
		Scalar x, y;
		uint32 boundarySize = boundary_vertices.size();
		uint32 sideLength = boundarySize / 4;
		for (uint32 i = 0; i < boundarySize; ++i)
		{
			uint32 idx = index_of(m, boundary_vertices[i]);
			if (i < sideLength)
			{
				// Bottom edge
				x = static_cast<Scalar>(i) / (sideLength - 1);
				y = 0.0;
			}
			else if (i < 2 * sideLength)
			{
				// Right edge
				x = 1.0;
				y = static_cast<Scalar>(i - sideLength) / (sideLength - 1);
			}
			else if (i < 3 * sideLength)
			{
				// Top edge
				x = 1.0 - static_cast<Scalar>(i - 2 * sideLength) / (sideLength - 1);
				y = 1.0;
			}
			else
			{
				// Left edge
				x = 0.0;
				y = 1.0 - static_cast<Scalar>(i - 3 * sideLength) / (sideLength - 1);
			}
			value<Vec2>(m, tc_position, boundary_vertices[i]) = Vec2(x, y);
		}	
	}
		break;
	default:
		break;
	}
}

Vec3 barycentric_coordinates(const Vec2& a, const Vec2& b, const Vec2& c, const Vec2& p)
{
	Vec2 v0 = b - a, v1 = c - a, v2 = p - a;
	Scalar d00 = v0.dot(v0);
	Scalar d01 = v0.dot(v1);
	Scalar d11 = v1.dot(v1);
	Scalar d20 = v2.dot(v0);
	Scalar d21 = v2.dot(v1);
	Scalar denom = d00 * d11 - d01 * d01;
	Scalar v = (d11 * d20 - d01 * d21) / denom;
	Scalar w = (d00 * d21 - d01 * d20) / denom;
	Scalar u = 1.0 - v - w;
	return Vec3(u, v, 1.0 - u - v);
}
template <typename SURFACE>
class Parameterization : public ViewModule
{
	static_assert(mesh_traits<SURFACE>::dimension >= 2, "Coverage axis can only be used with meshes of dimension >= 2");
	template <typename T>
	using SAttribute = typename mesh_traits<SURFACE>::template Attribute<T>;

	using SVertex = typename mesh_traits<SURFACE>::Vertex;
	using SEdge = typename mesh_traits<SURFACE>::Edge;
	using SFace = typename mesh_traits<SURFACE>::Face;



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
		std::shared_ptr<SAttribute<Vec3>> surface_tc3_position_ = nullptr;
		std::shared_ptr<SAttribute<Vec2>> surface_tc2_position_ = nullptr;
		std::shared_ptr<SAttribute<bool>> is_boundary_vertex_= nullptr;
		std::vector<SEdge> boundary_edge;
		boundary_form boundary_form_;
		parameterization_weight weight_;
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
		p.surface_tc2_position_ = get_or_add_attribute<Vec2, SVertex>(*p.surface_, "tc2_position");
		p.surface_tc3_position_ = get_or_add_attribute<Vec3, SVertex>(*p.surface_, "tc3_position");
		p.is_boundary_vertex_ = get_or_add_attribute<bool, SVertex>(*p.surface_, "is_boundary_vertex");

		p.boundary_form_ = BOUNDARY_FORM_CIRCLE;
		foreach_cell(*p.surface_, [&](SVertex v) {
			value<bool>(*p.surface_, p.is_boundary_vertex_, v) = false;
			return true;
		});
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
		fix_boundary_vertices(*p.surface_, p.boundary_form_, p.is_boundary_vertex_.get(), p.surface_tc2_position_.get());
		typedef Eigen::SparseMatrix<Scalar> SparseMatrix;
		typedef Eigen::Triplet<Scalar> Triplet;
		std::vector<Triplet> triplets;
		uint32 nb_vertices = nb_cells<SVertex>(*p.surface_);
		SparseMatrix A(nb_vertices, nb_vertices);
		Eigen::MatrixXd B(nb_vertices, 2);
		B.setZero();
		auto geodesic_polar_map_coord = add_attribute<Vec2, SVertex>(*p.surface_, "geodesic_polar_map_coord");
		foreach_cell(*p.surface_, [&](SVertex v) {
			uint32 v_index = index_of(*p.surface_, v);
			if ((*p.is_boundary_vertex_)[v_index])
			{
				Vec2 pos = value<Vec2>(*p.surface_, p.surface_tc2_position_, v);
				B(v_index, 0) = pos.x();
				B(v_index, 1) = pos.y();
				triplets.push_back(Triplet(v_index, v_index, 1.0));
			}
			else
			{
				switch (p.weight_)
				{
				case WEIGHT_UNIFORM: {
					foreach_adjacent_vertex_through_edge(*p.surface_, v, [&](SVertex iv) {
						uint32 iv_index = index_of(*p.surface_, iv);
						triplets.push_back(Triplet(v_index, iv_index, -1.0));
						return true;
					});
					triplets.push_back(Triplet(v_index, v_index, degree(*p.surface_, v)));
				}
				break;
				// Derived from "Parametrization and smooth approximation of surface triangulations", Section 6
				// Shape prevsering parameterization
				case WEIGHT_FLOATER: {
					Scalar sum = 0.0;
					(*geodesic_polar_map_coord)[v_index] = Vec2(0.0, 0.0);
					std::vector<SVertex> one_ring = adjacent_vertices_through_edge(*p.surface_, v);
					Scalar angle_sum = geometry::vertex_angle_sum(*p.surface_, p.surface_vertex_position_.get(), v);

					uint32 degree = one_ring.size();
					// geodesic polar mapping
					// ||p_{k} - p|| = ||x_{jk}-x_{i}||, angle(p_{k}, p, p_{k+}) = 2*Pi * angle(x_{jk}, x_{i}, x_{jk+1})/angle_sum
					Scalar sum_angle = 0;
					for (uint32 i = 0; i < degree; i++)
					{
						Vec3 pos = (*p.surface_vertex_position_)[v_index];

						SVertex iv = one_ring[i];
						uint32 iv_index = index_of(*p.surface_, iv);
						Vec3 pos_v = (*p.surface_vertex_position_)[iv_index];
						
						if (i == 0)
						{
							(*geodesic_polar_map_coord)[iv_index] = Vec2((pos_v - pos).norm(), 0);
						}
						else
						{
							SVertex iv_prev = one_ring[(i - 1)];
							Vec3 pos_prev = value<Vec3>(*p.surface_, p.surface_vertex_position_, iv_prev);
							sum_angle += geometry::angle(pos_v - pos, pos_prev - pos);
							Scalar percent = (sum_angle / angle_sum) * 2 * M_PI;
							(*geodesic_polar_map_coord)[iv_index] =
								Vec2(std::cos(percent), std::sin(percent)) * (pos_v - pos).norm();
						}
					}
					if (degree == 3)
					{
						SVertex va = one_ring[0];
						SVertex vb = one_ring[1];
						SVertex vc = one_ring[2];
						uint32 va_index = index_of(*p.surface_, va);
						uint32 vb_index = index_of(*p.surface_, vb);
						uint32 vc_index = index_of(*p.surface_, vc);
						Vec2 pos = (*geodesic_polar_map_coord)[v_index];
						Vec2 pos_a = (*geodesic_polar_map_coord)[va_index];
						Vec2 pos_b = (*geodesic_polar_map_coord)[vb_index];
						Vec2 pos_c = (*geodesic_polar_map_coord)[vc_index];
						Vec3 bary_coord = barycentric_coordinates(pos_a, pos_b, pos_c, pos);
						triplets.push_back(Triplet(v_index, va_index, bary_coord[0]));
						triplets.push_back(Triplet(v_index, vb_index, bary_coord[1]));
						triplets.push_back(Triplet(v_index, vc_index, bary_coord[2]));
						triplets.push_back(Triplet(v_index, v_index, -bary_coord[0] - bary_coord[1] - bary_coord[2]));
					}
					else
					{
						Eigen::MatrixXd mu(degree, degree);
						mu.setZero();
						for (uint32 l = 0; l < degree; l++)
						{
							Vec2 pos = (*geodesic_polar_map_coord)[v_index];
							SVertex lv = one_ring[l];
							uint32 lv_index = index_of(*p.surface_, lv);
							Vec2 pos_lv = (*geodesic_polar_map_coord)[lv_index];
							for (uint32 k = 0; k < degree; k++)
							{
								if (k == l || (k+1)%degree == l)
									continue;
								SVertex kv = one_ring[k];
								SVertex kv_next = one_ring[(k + 1) % degree];

								uint32 kv_index = index_of(*p.surface_, kv);
								uint32 kv_next_index = index_of(*p.surface_, kv_next);
								Vec2 pos_kv = (*geodesic_polar_map_coord)[kv_index];
								Vec2 pos_kv_next = (*geodesic_polar_map_coord)[kv_next_index];
								Vec3 bary_coord = barycentric_coordinates(pos_lv, pos_kv, pos_kv_next, pos);
								if (bary_coord[0] > 0 && bary_coord[1] > 0 && bary_coord[2] > 0)
								{
									mu(l, l) = bary_coord[0];
									mu(k, l) = bary_coord[1];
									mu((k + 1) % degree, l) = bary_coord[2];
									break;
								}
							}
						}
						Eigen::VectorXd sum_mu = mu.rowwise().sum();
						double sum_weight = 0;
						for (uint32 idx = 0; idx < degree; idx++)
						{
							triplets.push_back(Triplet(v_index, index_of(*p.surface_,one_ring[idx]), sum_mu(idx) / degree));
							sum_weight += sum_mu(idx) / degree;
						}
						triplets.push_back(Triplet(v_index, v_index, -sum_weight));
					}
				}
				break;
				default:
					break;
				}
			}
			return true;
		});
		A.setFromTriplets(triplets.begin(), triplets.end());
		Eigen::SparseLU<SparseMatrix> solver;
		solver.compute(A);
		Eigen::MatrixXd X = solver.solve(B);
		foreach_cell(*p.surface_, [&](SVertex v) {
			uint32 v_index = index_of(*p.surface_, v);
			value<Vec2>(*p.surface_, p.surface_tc2_position_, v) = Vec2(X(v_index, 0), X(v_index, 1));
			value<Vec3>(*p.surface_, p.surface_tc3_position_,v) = Vec3(X(v_index, 0), X(v_index, 1), 0);
			return true;
		});
		remove_attribute<SVertex>(*p.surface_, geodesic_polar_map_coord.get());
		surface_provider_->emit_attribute_changed(*p.surface_, p.surface_tc2_position_.get());
		surface_provider_->emit_attribute_changed(*p.surface_, p.surface_tc3_position_.get());

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
			
			if (ImGui::Button("Parameterization init"))
			{
				parameterization_init(*selected_surface_mesh_);
			}
			
			if (p.initialized_)
			{
				if (ImGui::Button("Cut mesh"))
				{
					cut_mesh(p);
				}
				ImGui::RadioButton("Circle", (int*)&p.boundary_form_, BOUNDARY_FORM_CIRCLE);
				ImGui::SameLine();
				ImGui::RadioButton("Square", (int*)&p.boundary_form_, BOUNDARY_FORM_SQUARE);

				ImGui::RadioButton("Uniform", (int*)&p.weight_, WEIGHT_UNIFORM);
				ImGui::SameLine();
				ImGui::RadioButton("Floater", (int*)&p.weight_, WEIGHT_FLOATER);
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