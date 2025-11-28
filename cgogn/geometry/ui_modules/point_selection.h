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

#ifndef CGOGN_MODULE_POINT_SELECTION_H_
#define CGOGN_MODULE_POINT_SELECTION_H_

#include <cgogn/ui/app.h>
#include <cgogn/ui/imgui_helpers.h>
#include <cgogn/ui/module.h>
#include <cgogn/ui/view.h>

#include <cgogn/core/functions/attributes.h>
#include <cgogn/core/functions/traversals/vertex.h>
#include <cgogn/core/ui_modules/mesh_provider.h>

#include <cgogn/geometry/functions/distance.h>
#include <cgogn/geometry/types/vector_traits.h>

#include <cgogn/rendering/shaders/shader_point_sprite.h>
#include <cgogn/rendering/vbo_update.h>

#include <boost/synapse/connect.hpp>

#include <GLFW/glfw3.h>

#include <algorithm>
#include <limits>
#include <optional>
#include <unordered_map>
#include <memory>

#undef near
#undef far

namespace cgogn
{

namespace ui
{

using geometry::Scalar;
using geometry::Vec3;
using geometry::Vec4;

/**
 * Simple vertex selection module for point-cloud meshes (CMap0).
 * Selection is toggled with the S key (left click to select, right click to unselect).
 */
template <typename MESH>
class PointSelection : public ViewModule
{
	static_assert(mesh_traits<MESH>::dimension == 0, "PointSelection can only be used with meshes of dimension 0");

	template <typename T>
	using Attribute = typename mesh_traits<MESH>::template Attribute<T>;

	using Vertex = typename mesh_traits<MESH>::Vertex;

	enum SelectionMethod
	{
		SingleCell = 0,
		WithinSphere
	};

	enum SelectionColor
	{
		Uniform = 0,
		PerVertex
	};

	struct Parameters
	{
		Parameters()
			: mesh_(nullptr), vertex_position_(nullptr), vertex_scale_factor_(1.0f), sphere_scale_factor_(10.0f),
			  selected_vertices_set_(nullptr), selection_method_(SingleCell), selection_color_(Uniform)
		{
			param_point_sprite_ = rendering::ShaderPointSprite::generate_param();
			param_point_sprite_->color_ = rendering::GLColor(1, 0, 0, 0.65f);
			param_point_sprite_->set_vbos({&selected_vertices_vbo_});

			param_point_sprite_color_ = rendering::ShaderPointSpriteColor::generate_param();
			param_point_sprite_color_->set_vbos({&selected_vertices_vbo_, &selected_vertices_color_vbo_});
		}

		CGOGN_NOT_COPYABLE_NOR_MOVABLE(Parameters);

		void update_selected_vertices_vbo()
		{
			if (selected_vertices_set_)
			{
				std::vector<Vec3> selected_vertices_position;
				selected_vertices_position.reserve(selected_vertices_set_->size());
				selected_vertices_set_->foreach_cell(
					[&](Vertex v) { selected_vertices_position.push_back(value<Vec3>(*mesh_, vertex_position_, v)); });
				rendering::update_vbo(selected_vertices_position, &selected_vertices_vbo_);

				if (vertex_color_)
				{
					std::vector<Vec4> selected_vertices_color;
					selected_vertices_color.reserve(selected_vertices_set_->size());
					selected_vertices_set_->foreach_cell(
						[&](Vertex v) { selected_vertices_color.push_back(value<Vec4>(*mesh_, vertex_color_, v)); });
					rendering::update_vbo(selected_vertices_color, &selected_vertices_color_vbo_);
				}
				else
				{
					// clear color buffer so shader reads a defined value (unused if Uniform)
					std::vector<Vec4> dummy;
					rendering::update_vbo(dummy, &selected_vertices_color_vbo_);
				}
			}
		}

		MESH* mesh_;
		std::shared_ptr<Attribute<Vec3>> vertex_position_;
		std::shared_ptr<Attribute<Vec4>> vertex_color_;

		std::unique_ptr<rendering::ShaderPointSprite::Param> param_point_sprite_;
		std::unique_ptr<rendering::ShaderPointSpriteColor::Param> param_point_sprite_color_;

		float32 vertex_scale_factor_;
		float32 vertex_base_size_;
		float32 sphere_scale_factor_;

		rendering::VBO selected_vertices_vbo_;
		rendering::VBO selected_vertices_color_vbo_;

		CellsSet<MESH, Vertex>* selected_vertices_set_;

		SelectionMethod selection_method_;
		SelectionColor selection_color_;

		struct DisplayedSet
		{
			CellsSet<MESH, Vertex>* set_;
			std::unique_ptr<rendering::VBO> pos_vbo_;
			std::unique_ptr<rendering::VBO> color_vbo_;
		};
		std::vector<DisplayedSet> displayed_sets_;
	};

public:
	explicit PointSelection(const App& app)
		: ViewModule(app, "PointSelection (" + std::string{mesh_traits<MESH>::name} + ")"), selected_mesh_(nullptr),
		  selecting_(false), pick_max_distance_factor_(3.0f)
	{
		param_point_sprite_selecting_ = rendering::ShaderPointSprite::generate_param();
		param_point_sprite_selecting_->color_ = rendering::GLColor(1, 0.5, 0, 0.65f);
		param_point_sprite_selecting_->set_vbos({&selecting_vertices_vbo_});
	}

	~PointSelection()
	{
	}

	void set_vertex_position(const MESH& m, const std::shared_ptr<Attribute<Vec3>>& vertex_position)
	{
		Parameters& p = parameters_[&m];

		p.vertex_position_ = vertex_position;
		if (p.vertex_position_)
		{
			update_vertex_base_size(p);
			p.update_selected_vertices_vbo();
		}
		else
			selecting_ = false;

		for (View* v : linked_views_)
			v->request_update();

		// refresh VAO bindings (matches PointCloudRender behavior)
		p.param_point_sprite_->set_vbos({&p.selected_vertices_vbo_});
		p.param_point_sprite_color_->set_vbos({&p.selected_vertices_vbo_, &p.selected_vertices_color_vbo_});
		param_point_sprite_selecting_->set_vbos({&selecting_vertices_vbo_});

		for (auto& ds : p.displayed_sets_)
		{
			update_displayed_set_vbo(p, ds);
		}
	}

	void set_vertex_color(const MESH& m, const std::shared_ptr<Attribute<Vec4>>& vertex_color)
	{
		Parameters& p = parameters_[&m];

		p.vertex_color_ = vertex_color;
		p.selection_color_ = vertex_color ? PerVertex : Uniform;
		p.update_selected_vertices_vbo();

		// refresh VAO bindings to use the (possibly new) color buffer
		p.param_point_sprite_color_->set_vbos({&p.selected_vertices_vbo_, &p.selected_vertices_color_vbo_});

		for (View* v : linked_views_)
			v->request_update();

		for (auto& ds : p.displayed_sets_)
		{
			update_displayed_set_vbo(p, ds);
		}
	}

protected:
	void init() override
	{
		mesh_provider_ = static_cast<ui::MeshProvider<MESH>*>(
			app_.module("MeshProvider (" + std::string{mesh_traits<MESH>::name} + ")"));
		mesh_provider_->foreach_mesh([this](MESH& m, const std::string&) { init_mesh(&m); });
		connections_.push_back(boost::synapse::connect<typename MeshProvider<MESH>::mesh_added>(
			mesh_provider_, this, &PointSelection<MESH>::init_mesh));
	}

	void key_press_event(View*, int32 key_code) override
	{
		if (selected_mesh_ && key_code == GLFW_KEY_S)
		{
			Parameters& p = parameters_[selected_mesh_];

			if (p.vertex_position_ && p.selected_vertices_set_)
			{
				selecting_ = true;
				selecting_vertices_.clear();
			}
		}
	}

	void key_release_event(View*, int32 key_code) override
	{
		if (key_code == GLFW_KEY_S)
		{
			selecting_ = false;
			selecting_vertices_.clear();

			for (View* v : linked_views_)
				v->request_update();
		}
	}

	void mouse_wheel_event(View* view, int32, int32 dy) override
	{
		if (selecting_)
		{
			Parameters& p = parameters_[selected_mesh_];
			if (p.selection_method_ == WithinSphere)
			{
				p.sphere_scale_factor_ += dy > 0 ? 1.0f : -1.0f;
				p.sphere_scale_factor_ = std::clamp(p.sphere_scale_factor_, 1.0f, 100.0f);
				view->stop_event();
				view->request_update();
			}
		}
	}

	void mouse_move_event(View* view, int32 x, int32 y) override
	{
		update_selecting_vertices(view, x, y);
	}

	void mouse_press_event(View*, int32 button, int32, int32) override
	{
		if (!selecting_ || selecting_vertices_.empty())
			return;

		Parameters& p = parameters_[selected_mesh_];
		switch (p.selection_method_)
		{
		case SingleCell:
			switch (button)
			{
			case 0:
				p.selected_vertices_set_->select(selecting_vertices_.front());
				break;
			case 1:
				p.selected_vertices_set_->unselect(selecting_vertices_.front());
				break;
			}
			mesh_provider_->emit_cells_set_changed(*selected_mesh_, p.selected_vertices_set_);
			break;
		case WithinSphere: {
			Scalar radius = p.vertex_base_size_ * p.sphere_scale_factor_;
			const Vec3& center = value<Vec3>(*p.mesh_, p.vertex_position_, selecting_vertices_.front());
			Scalar r2 = radius * radius;
			switch (button)
			{
			case 0:
				foreach_cell(*p.mesh_, [&](Vertex v) -> bool {
					if ((value<Vec3>(*p.mesh_, p.vertex_position_, v) - center).squaredNorm() <= r2)
						p.selected_vertices_set_->select(v);
					return true;
				});
				break;
			case 1:
				foreach_cell(*p.mesh_, [&](Vertex v) -> bool {
					if ((value<Vec3>(*p.mesh_, p.vertex_position_, v) - center).squaredNorm() <= r2)
						p.selected_vertices_set_->unselect(v);
					return true;
				});
				break;
			}
			mesh_provider_->emit_cells_set_changed(*selected_mesh_, p.selected_vertices_set_);
		}
		break;
		}
	}

	void draw(View* view) override
	{
		const rendering::GLMat4& proj_matrix = view->projection_matrix();
		const rendering::GLMat4& view_matrix = view->modelview_matrix();

		for (auto& [m, p] : parameters_)
		{
			// draw extra displayed sets
			for (auto& ds : p.displayed_sets_)
			{
				if (!ds.set_ || ds.set_->size() == 0)
					continue;

				if (p.vertex_color_ && p.selection_color_ == PerVertex &&
					p.param_point_sprite_color_->attributes_initialized() && ds.pos_vbo_ && ds.color_vbo_)
				{
					p.param_point_sprite_color_->set_vbos({ds.pos_vbo_.get(), ds.color_vbo_.get()});
					p.param_point_sprite_color_->point_size_ = p.vertex_base_size_ * p.vertex_scale_factor_;
					p.param_point_sprite_color_->bind(proj_matrix, view_matrix);
					glEnable(GL_BLEND);
					glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
					glDrawArrays(GL_POINTS, 0, ds.set_->size());
					glDisable(GL_BLEND);
					p.param_point_sprite_color_->release();
				}
				else if (p.param_point_sprite_->attributes_initialized() && ds.pos_vbo_)
				{
					p.param_point_sprite_->set_vbos({ds.pos_vbo_.get()});
					p.param_point_sprite_->point_size_ = p.vertex_base_size_ * p.vertex_scale_factor_;
					p.param_point_sprite_->bind(proj_matrix, view_matrix);
					glEnable(GL_BLEND);
					glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
					glDrawArrays(GL_POINTS, 0, ds.set_->size());
					glDisable(GL_BLEND);
					p.param_point_sprite_->release();
				}
			}

			if (p.selected_vertices_set_ && p.selected_vertices_set_->size() > 0)
			{
				if (p.vertex_color_ && p.selection_color_ == PerVertex &&
					p.param_point_sprite_color_->attributes_initialized())
				{
					p.param_point_sprite_color_->point_size_ = p.vertex_base_size_ * p.vertex_scale_factor_;
					p.param_point_sprite_color_->bind(proj_matrix, view_matrix);
					glEnable(GL_BLEND);
					glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
					glDrawArrays(GL_POINTS, 0, p.selected_vertices_set_->size());
					glDisable(GL_BLEND);
					p.param_point_sprite_color_->release();
				}
				else if (p.param_point_sprite_->attributes_initialized())
				{
					p.param_point_sprite_->point_size_ = p.vertex_base_size_ * p.vertex_scale_factor_;
					p.param_point_sprite_->bind(proj_matrix, view_matrix);
					glEnable(GL_BLEND);
					glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
					glDrawArrays(GL_POINTS, 0, p.selected_vertices_set_->size());
					glDisable(GL_BLEND);
					p.param_point_sprite_->release();
				}
			}
		}

		if (selecting_ && !selecting_vertices_.empty() && param_point_sprite_selecting_->attributes_initialized())
		{
			param_point_sprite_selecting_->point_size_ =
				parameters_[selected_mesh_].vertex_base_size_ * parameters_[selected_mesh_].vertex_scale_factor_;
			param_point_sprite_selecting_->bind(proj_matrix, view_matrix);
			glEnable(GL_BLEND);
			glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
			glDrawArrays(GL_POINTS, 0, selecting_vertices_.size());
			glDisable(GL_BLEND);
			param_point_sprite_selecting_->release();
		}
	}

	void left_panel() override
	{
		bool need_update = false;

		imgui_mesh_selector(mesh_provider_, selected_mesh_, "Points", [&](MESH& m) {
			selected_mesh_ = &m;
			mesh_provider_->mesh_data(m).outlined_until_ = App::frame_time_ + 1.0;
		});

		if (selected_mesh_)
		{
			Parameters& p = parameters_[selected_mesh_];

			imgui_combo_attribute<Vertex, Vec3>(*selected_mesh_, p.vertex_position_, "Position",
												[&](const std::shared_ptr<Attribute<Vec3>>& attribute) {
													set_vertex_position(*selected_mesh_, attribute);
												});

			if (p.vertex_position_)
			{
				ImGui::Separator();
				ImGui::RadioButton("Single", reinterpret_cast<int*>(&p.selection_method_), SingleCell);
				ImGui::SameLine();
				ImGui::RadioButton("Sphere", reinterpret_cast<int*>(&p.selection_method_), WithinSphere);

				if (p.selection_method_ == WithinSphere)
					ImGui::SliderFloat("Sphere radius", &(p.sphere_scale_factor_), 1.0f, 100.0f);

				MeshData<MESH>& md = mesh_provider_->mesh_data(*selected_mesh_);

				if (ImGui::Button("Create set##vertices_set"))
					md.template add_cells_set<Vertex>();
				imgui_combo_cells_set(md, p.selected_vertices_set_, "Sets", [&](CellsSet<MESH, Vertex>* cs) {
					p.selected_vertices_set_ = cs;
					p.update_selected_vertices_vbo();
					need_update = true;
				});

				ImGui::Separator();
				ImGui::TextUnformatted("Displayed sets");
				if (p.mesh_ && p.vertex_position_ && ImGui::Button("Add new set##add_display_set"))
				{
					p.displayed_sets_.push_back(
						{nullptr, std::make_unique<rendering::VBO>(3), std::make_unique<rendering::VBO>(4)});
				}
				if (!p.displayed_sets_.empty())
				{
					uint32 idx = 0;
					for (auto it = p.displayed_sets_.begin(); it != p.displayed_sets_.end();)
					{
						std::string label = "Set##disp" + std::to_string(idx);
						imgui_combo_cells_set(md, it->set_, label, [&](CellsSet<MESH, Vertex>* cs) {
							it->set_ = cs;
							update_displayed_set_vbo(p, *it);
							need_update = true;
						});
						ImGui::SameLine();
						ImGui::Text("(nb: %u)", it->set_ ? it->set_->size() : 0u);
						ImGui::SameLine();
						if (ImGui::SmallButton((std::string("Remove##disp") + std::to_string(idx)).c_str()))
							it = p.displayed_sets_.erase(it);
						else
						{
							++it;
						}
						++idx;
					}
				}

				if (p.selected_vertices_set_)
				{
					ImGui::TextUnformatted("Press S to select vertices");
					ImGui::Text("(nb elements: %d)", p.selected_vertices_set_->size());
					if (ImGui::Button("Clear##vertices_set"))
					{
						p.selected_vertices_set_->clear();
						mesh_provider_->emit_cells_set_changed(*selected_mesh_, p.selected_vertices_set_);
					}
					if (ImGui::Button("Select all##vertices_set"))
					{
						p.selected_vertices_set_->select_if([&](Vertex) { return true; });
						mesh_provider_->emit_cells_set_changed(*selected_mesh_, p.selected_vertices_set_);
					}
				}
				ImGui::TextUnformatted("Colors");
				ImGui::BeginGroup();
				if (ImGui::RadioButton("Global##color", p.selection_color_ == Uniform))
				{
					p.selection_color_ = Uniform;
					need_update = true;
				}
				ImGui::SameLine();
				if (ImGui::RadioButton("Per vertex##color", p.selection_color_ == PerVertex))
				{
					p.selection_color_ = PerVertex;
					need_update = true;
				}
				ImGui::EndGroup();

				if (p.selection_color_ == Uniform)
				{
					need_update |= ImGui::ColorEdit4("Color##vertices", p.param_point_sprite_->color_.data(),
													 ImGuiColorEditFlags_NoInputs);
				}
				else if (p.selection_color_ == PerVertex)
				{
					imgui_combo_attribute<Vertex, Vec4>(*selected_mesh_, p.vertex_color_, "Attribute##vertexcolor",
														[&](const std::shared_ptr<Attribute<Vec4>>& attribute) {
															set_vertex_color(*selected_mesh_, attribute);
															need_update = true;
														});
				}

				ImGui::TextUnformatted("Drawing parameters");
				need_update |= ImGui::SliderFloat("size##vertices", &(p.vertex_scale_factor_), 0.1f, 2.0f);
			}
		}

		if (need_update)
			for (View* v : linked_views_)
				v->request_update();
	}

private:
	void init_mesh(MESH* m)
	{
		Parameters& p = parameters_[m];
		p.mesh_ = m;
		mesh_connections_[m].push_back(
			boost::synapse::connect<typename MeshProvider<MESH>::template attribute_changed_t<Vec3>>(
				m, [this, m](Attribute<Vec3>* attribute) {
					Parameters& p = parameters_[m];
					if (p.vertex_position_.get() == attribute)
					{
						update_vertex_base_size(p);
						p.update_selected_vertices_vbo();
						for (auto& ds : p.displayed_sets_)
							update_displayed_set_vbo(p, ds);
					}

					for (View* v : linked_views_)
						v->request_update();
				}));
		mesh_connections_[m].push_back(
			boost::synapse::connect<typename MeshProvider<MESH>::template attribute_changed_t<Vec4>>(
				m, [this, m](Attribute<Vec4>* attribute) {
					Parameters& p = parameters_[m];
					if (p.vertex_color_.get() == attribute)
					{
						p.update_selected_vertices_vbo();
						for (auto& ds : p.displayed_sets_)
							update_displayed_set_vbo(p, ds);
						for (View* v : linked_views_)
							v->request_update();
					}
				}));
		mesh_connections_[m].push_back(
			boost::synapse::connect<typename MeshProvider<MESH>::template cells_set_changed<Vertex>>(
				m, [this, m](CellsSet<MESH, Vertex>* set) {
					Parameters& p = parameters_[m];
					if (p.selected_vertices_set_ == set && p.vertex_position_)
					{
						p.update_selected_vertices_vbo();
						for (View* v : linked_views_)
							v->request_update();
					}
					for (auto& ds : p.displayed_sets_)
					{
						if (ds.set_ == set && p.vertex_position_)
							update_displayed_set_vbo(p, ds);
					}
				}));

		// Default to a "color" attribute if present to mirror point cloud renderer behavior
		if (!p.vertex_color_)
		{
			auto default_color = get_attribute<Vec4, Vertex>(*m, "color");
			if (default_color)
				set_vertex_color(*m, default_color);
		}
	}

	void update_vertex_base_size(Parameters& p)
	{
		MeshData<MESH>& md = mesh_provider_->mesh_data(*p.mesh_);
		Vec3 diag = md.bb_max_ - md.bb_min_;
		Scalar diag_norm = diag.norm();
		p.vertex_base_size_ = float32(diag_norm > 0.0 ? diag_norm / 200.0 : 0.01);
	}

	void update_displayed_set_vbo(Parameters& p, typename Parameters::DisplayedSet& ds)
	{
		if (!p.mesh_ || !ds.set_ || !p.vertex_position_)
			return;

		std::vector<Vec3> pos;
		pos.reserve(ds.set_->size());
		ds.set_->foreach_cell([&](Vertex v) -> bool {
			pos.push_back(value<Vec3>(*p.mesh_, p.vertex_position_, v));
			return true;
		});
		if (ds.pos_vbo_)
			rendering::update_vbo(pos, ds.pos_vbo_.get());

		if (p.vertex_color_ && p.selection_color_ == PerVertex)
		{
			std::vector<Vec4> cols;
			cols.reserve(ds.set_->size());
			ds.set_->foreach_cell([&](Vertex v) -> bool {
				cols.push_back(value<Vec4>(*p.mesh_, p.vertex_color_, v));
				return true;
			});
			if (ds.color_vbo_)
				rendering::update_vbo(cols, ds.color_vbo_.get());
		}
		else
		{
			std::vector<Vec4> dummy;
			if (ds.color_vbo_)
				rendering::update_vbo(dummy, ds.color_vbo_.get());
		}
	}

	void update_selecting_vertices(View* view, int32 x, int32 y)
	{
		if (!selecting_ || !selected_mesh_)
			return;

		Parameters& p = parameters_[selected_mesh_];
		if (!p.vertex_position_ || !p.selected_vertices_set_)
			return;

		rendering::GLVec3d near = view->unproject(x, y, 0.0);
		rendering::GLVec3d far_d = view->unproject(x, y, 1.0);
		Vec3 A{near.x(), near.y(), near.z()};
		Vec3 B{far_d.x(), far_d.y(), far_d.z()};

		selecting_vertices_.clear();
		if (auto picked = pick_closest_vertex(A, B, p))
		{
			selecting_vertices_.push_back(*picked);
			std::vector<Vec3> selecting_points;
			selecting_points.reserve(selecting_vertices_.size());
			for (Vertex v : selecting_vertices_)
				selecting_points.push_back(value<Vec3>(*p.mesh_, p.vertex_position_, v));
			rendering::update_vbo(selecting_points, &selecting_vertices_vbo_);
		}

		for (View* v : linked_views_)
			v->request_update();
	}

	std::optional<Vertex> pick_closest_vertex(const Vec3& A, const Vec3& B, Parameters& p) const
	{
		Vec3 AB = B - A;
		Scalar ab2 = AB.squaredNorm();
		if (ab2 <= std::numeric_limits<Scalar>::epsilon())
			return std::nullopt;
		AB.normalize();

		Scalar best_d2 = std::numeric_limits<Scalar>::max();
		Vertex best_vertex;
		Scalar max_d2 =
			(p.vertex_base_size_ * pick_max_distance_factor_) * (p.vertex_base_size_ * pick_max_distance_factor_);

		foreach_cell(*p.mesh_, [&](Vertex v) -> bool {
			const Vec3& pos = value<Vec3>(*p.mesh_, p.vertex_position_, v);
			Scalar along_ray = AB.dot(pos - A);
			if (along_ray < 0.0)
				return true;
			Scalar d2 = geometry::squared_distance_normalized_line_point(A, AB, pos);
			if (d2 < best_d2)
			{
				best_d2 = d2;
				best_vertex = v;
			}
			return true;
		});

		if (best_d2 < std::min(max_d2, std::numeric_limits<Scalar>::max()))
			return best_vertex;
		else
			return std::nullopt;
	}

private:
	const MESH* selected_mesh_;
	std::unordered_map<const MESH*, Parameters> parameters_;
	std::vector<std::shared_ptr<boost::synapse::connection>> connections_;
	std::unordered_map<const MESH*, std::vector<std::shared_ptr<boost::synapse::connection>>> mesh_connections_;
	MeshProvider<MESH>* mesh_provider_;

	bool selecting_;
	std::vector<Vertex> selecting_vertices_;
	rendering::VBO selecting_vertices_vbo_;
	std::unique_ptr<rendering::ShaderPointSprite::Param> param_point_sprite_selecting_;
	float32 pick_max_distance_factor_;
};

} // namespace ui

} // namespace cgogn

#endif // CGOGN_MODULE_POINT_SELECTION_H_
