#ifndef CGOGN_MODULE_SHRINKING_BALL_DEBUGGER_H_
#define CGOGN_MODULE_SHRINKING_BALL_DEBUGGER_H_

#include <cgogn/ui/app.h>
#include <cgogn/ui/imgui_helpers.h>
#include <cgogn/ui/module.h>
#include <cgogn/ui/view.h>

#include <cgogn/core/functions/attributes.h>
#include <cgogn/core/functions/traversals/vertex.h>
#include <cgogn/core/ui_modules/mesh_provider.h>
#include <cgogn/core/utils/numerics.h>

#include <cgogn/geometry/algos/medial_axis.h>
#include <cgogn/geometry/functions/angle.h>
#include <cgogn/geometry/functions/distance.h>
#include <cgogn/geometry/types/vector_traits.h>

#include <cgogn/rendering/ui_modules/point_cloud_render.h>
#include <cgogn/rendering/shaders/shader_bold_line.h>
#include <cgogn/rendering/vbo_update.h>

#include <libacc/kd_tree.h>

#include <GLFW/glfw3.h>

#include <algorithm>
#include <cmath>
#include <functional>
#include <limits>
#include <string>
#include <unordered_map>
#include <vector>

namespace cgogn
{

namespace ui
{

using geometry::Scalar;
using geometry::Vec3;
using geometry::Vec4;

template <typename POINTS>
class ShrinkingBallDebugger : public ViewModule
{
public:
	template <typename T>
	using PAttribute = typename mesh_traits<POINTS>::template Attribute<T>;
	using PVertex = typename mesh_traits<POINTS>::Vertex;
	using CenterUDFQuery = std::function<bool(POINTS&, const Vec3&, Scalar&)>;

	ShrinkingBallDebugger(const App& app)
		: ViewModule(app, "ShrinkingBallDebugger (" + std::string{mesh_traits<POINTS>::name} + ")"),
		  selected_view_(app.current_view())
	{
	}

	~ShrinkingBallDebugger()
	{
		for (auto& [mesh, p] : parameters_)
		{
			if (p.kdtree_)
				delete p.kdtree_;
		}
	}

	void set_center_udf_query(CenterUDFQuery query)
	{
		center_udf_query_ = std::move(query);
	}

private:
	struct Parameters
	{
		Parameters()
			: points_(nullptr), position_(nullptr), normal_(nullptr), color_(nullptr), kdtree_(nullptr),
			  ball_center_(0, 0, 0), ball_radius_(0), ball_prev_contact_pos_(0, 0, 0), ball_iter_(0),
			  ball_initialized_(false), ball_finished_(false), normal_length_(1.0f), initial_radius_(0.5f),
			  sphere_alpha_(0.5f), sphere_mesh_(nullptr), sphere_position_(nullptr), sphere_radius_(nullptr),
			  sphere_color_(nullptr), sphere_ready_(false), normals_count_(0)
		{
			param_normals_ = rendering::ShaderBoldLine::generate_param();
			param_normals_->color_ = rendering::GLColor(1.0f, 1.0f, 0.0f, 1.0f);
			param_normals_->width_ = 2.0f;
			param_normals_->set_vbos({&normals_vbo_});
		}

		POINTS* points_ = nullptr;
		std::shared_ptr<PAttribute<Vec3>> position_ = nullptr;
		std::shared_ptr<PAttribute<Vec3>> normal_ = nullptr;
		std::shared_ptr<PAttribute<Vec4>> color_ = nullptr;

		acc::KDTree<3, uint32>* kdtree_ = nullptr;
		std::vector<PVertex> kdtree_vertices_;

		PVertex picked_vertex_;
		PVertex contact_vertex_;
		PVertex last_picked_;
		PVertex last_contact_;

		Vec3 ball_center_ = Vec3(0, 0, 0);
		Scalar ball_radius_ = Scalar(0);
		Vec3 ball_prev_contact_pos_ = Vec3(0, 0, 0);
		PVertex ball_prev_contact_vertex_;
		uint32 ball_iter_ = 0;
		bool ball_initialized_ = false;
		bool ball_finished_ = false;

		float32 normal_length_ = 1.0f;

		float32 initial_radius_ = 0.5f;
		float32 sphere_alpha_ = 0.5f;

		POINTS* sphere_mesh_ = nullptr;
		std::shared_ptr<PAttribute<Vec3>> sphere_position_ = nullptr;
		std::shared_ptr<PAttribute<Scalar>> sphere_radius_ = nullptr;
		std::shared_ptr<PAttribute<Vec4>> sphere_color_ = nullptr;
		PVertex sphere_vertex_;
		bool sphere_ready_ = false;

		rendering::VBO normals_vbo_;
		std::unique_ptr<rendering::ShaderBoldLine::Param> param_normals_;
		uint32 normals_count_ = 0;
	};

protected:
	void init() override
	{
		points_provider_ = static_cast<ui::MeshProvider<POINTS>*>(
			app_.module("MeshProvider (" + std::string{mesh_traits<POINTS>::name} + ")"));
		pcr_ = static_cast<PointCloudRender<POINTS>*>(
			app_.module("PointCloudRender (" + std::string{mesh_traits<POINTS>::name} + ")"));
		if (points_provider_)
			points_provider_->foreach_mesh([this](POINTS& m, const std::string&) { parameters_[&m]; });
	}

	void key_press_event(View* view, int32 key_code) override
	{
		if (!selected_points_)
			return;
		if (key_code == GLFW_KEY_I)
		{
			Parameters& p = parameters_[selected_points_];
			if (pick_vertex_under_cursor(view, p))
			{
				init_ball_from_pick(p);
				update_sphere_visual(p);
				update_normals_vbo(p);
			}
		}
	}

	void left_panel() override
	{
		if (!points_provider_)
			return;

		if (app_.nb_views() > 1)
			imgui_view_selector(this, selected_view_, [&](View* v) { selected_view_ = v; });
		else
			selected_view_ = app_.current_view();

		imgui_mesh_selector(points_provider_, selected_points_, "Cmap0", [&](POINTS& m) { set_selected_points(m); });

		if (!selected_points_)
			return;

		Parameters& p = parameters_[selected_points_];

		imgui_combo_attribute<PVertex, Vec3>(
			*selected_points_, p.position_, "Position",
			[&](const std::shared_ptr<PAttribute<Vec3>>& attribute) { set_position_attribute(p, attribute); });

		imgui_combo_attribute<PVertex, Vec3>(
			*selected_points_, p.normal_, "Normal",
			[&](const std::shared_ptr<PAttribute<Vec3>>& attribute) {
				p.normal_ = attribute;
				update_normals_vbo(p);
			});
		if (!p.normal_)
			ImGui::TextUnformatted("Normal attribute required for shrinking ball");

		if (ImGui::SliderFloat("Normal length", &p.normal_length_, 0.1f, 5.0f))
			update_normals_vbo(p);

		ImGui::Separator();
		ImGui::InputFloat("Initial radius", &p.initial_radius_, 0.0f, 0.0f, "%.6f");
		p.initial_radius_ = std::max(0.0f, p.initial_radius_);
		if (ImGui::SliderFloat("Sphere alpha", &p.sphere_alpha_, 0.0f, 1.0f))
			update_sphere_visual(p);

		if (ImGui::Button("Build KDTree"))
		{
			build_kdtree(p);
		}

		ImGui::Separator();
		ImGui::TextUnformatted("Pick point under cursor with I");
		if (p.picked_vertex_.is_valid())
		{
			uint32 idx = index_of(*selected_points_, p.picked_vertex_);
			ImGui::Text("Picked id: %u", idx);
			Scalar picked_udf = Scalar(0);
			if (center_udf_query_ && selected_points_ && p.position_ &&
				center_udf_query_(*selected_points_, (*p.position_)[idx], picked_udf))
				ImGui::Text("Picked udf: %f", static_cast<double>(picked_udf));
			else
				ImGui::TextUnformatted("Picked udf: N/A");
		}
		if (p.contact_vertex_.is_valid())
		{
			uint32 idx = index_of(*selected_points_, p.contact_vertex_);
			ImGui::Text("Contact id: %u", idx);
		}
		if (p.ball_initialized_)
		{
			ImGui::Text("Ball: r=%f", static_cast<double>(p.ball_radius_));
			Scalar center_udf = Scalar(0);
			if (center_udf_query_ && selected_points_ && center_udf_query_(*selected_points_, p.ball_center_, center_udf))
				ImGui::Text("Ball center udf: %f", static_cast<double>(center_udf));
			else
				ImGui::TextUnformatted("Ball center udf: N/A");
		}

		if (ImGui::Button("Shrink step"))
		{
			if (p.ball_finished_)
			{
				ImGui::OpenPopup("Shrinking ball finished");
			}
			else
			{
				step_shrinking_ball(p);
			}
		}

		if (ImGui::BeginPopupModal("Shrinking ball finished", nullptr, ImGuiWindowFlags_AlwaysAutoResize))
		{
			ImGui::TextUnformatted("shrinking ball finished");
			if (ImGui::Button("OK"))
				ImGui::CloseCurrentPopup();
			ImGui::EndPopup();
		}
	}

	void draw(View* view) override
	{
		if (!selected_points_)
			return;
		Parameters& p = parameters_[selected_points_];
		if (p.normals_count_ == 0 || !p.param_normals_ || !p.param_normals_->attributes_initialized())
			return;

		const rendering::GLMat4& proj_matrix = view->projection_matrix();
		const rendering::GLMat4& view_matrix = view->modelview_matrix();

		p.param_normals_->bind(proj_matrix, view_matrix);
		glDrawArrays(GL_LINES, 0, p.normals_count_);
		p.param_normals_->release();
	}

private:
	void set_selected_points(POINTS& m)
	{
		selected_points_ = &m;
		Parameters& p = parameters_[selected_points_];
		p.points_ = selected_points_;
		reset_state(p);
		if (p.position_ && selected_view_)
			apply_point_render_settings(p);
		update_normals_vbo(p);
	}

	void set_position_attribute(Parameters& p, const std::shared_ptr<PAttribute<Vec3>>& attribute)
	{
		p.position_ = attribute;
		if (!p.position_)
			return;

		ensure_color_attribute(p);
		apply_point_render_settings(p);
		build_kdtree(p);
		reset_ball(p);
		update_normals_vbo(p);
	}

	void apply_point_render_settings(Parameters& p)
	{
		if (!pcr_ || !selected_view_ || !p.points_ || !p.position_)
			return;

		pcr_->set_vertex_position(*selected_view_, *p.points_, p.position_);
		if (p.color_)
		{
			pcr_->set_vertex_color(*selected_view_, *p.points_, p.color_);
			pcr_->set_vertex_color_per_cell(*selected_view_, *p.points_, PointCloudRender<POINTS>::PER_VERTEX);
		}
	}

	void ensure_color_attribute(Parameters& p)
	{
		if (!p.points_)
			return;
		if (!p.color_)
			p.color_ = get_or_add_attribute<Vec4, PVertex>(*p.points_, "sbd_color");
		p.color_->fill(Vec4(1, 1, 1, 1));
		if (points_provider_ && p.points_)
			points_provider_->emit_attribute_changed(*p.points_, p.color_.get());
	}

	void reset_state(Parameters& p)
	{
		p.picked_vertex_ = PVertex();
		p.contact_vertex_ = PVertex();
		p.last_picked_ = PVertex();
		p.last_contact_ = PVertex();
		reset_ball(p);
		update_highlights(p);
		update_normals_vbo(p);
	}

	void update_normals_vbo(Parameters& p)
	{
		if (!p.points_ || !p.position_ || !p.normal_ || !p.picked_vertex_.is_valid())
		{
			p.normals_count_ = 0;
			rendering::update_vbo(std::vector<Vec3>{}, &p.normals_vbo_);
			for (View* v : linked_views_)
				v->request_update();
			return;
		}

		const Scalar min_norm = Scalar(1e-12);
		std::vector<Vec3> lines;
		lines.reserve(2);

		uint32 idx = index_of(*p.points_, p.picked_vertex_);
		const Vec3& pos = (*p.position_)[idx];
		Vec3 n = (*p.normal_)[idx];
		if (n.squaredNorm() > min_norm)
		{
			n.normalize();
			lines.push_back(pos);
			lines.push_back(pos + n * Scalar(p.normal_length_));
		}

		rendering::update_vbo(lines, &p.normals_vbo_);
		p.normals_count_ = static_cast<uint32>(lines.size());
		for (View* v : linked_views_)
			v->request_update();
	}

	void reset_ball(Parameters& p)
	{
		p.ball_initialized_ = false;
		p.ball_finished_ = false;
		p.ball_iter_ = 0;
		p.ball_center_ = Vec3(0, 0, 0);
		p.ball_radius_ = Scalar(0);
		p.ball_prev_contact_pos_ = Vec3(0, 0, 0);
		p.ball_prev_contact_vertex_ = PVertex();
		p.contact_vertex_ = PVertex();
		update_sphere_visual(p);
	}

	void build_kdtree(Parameters& p)
	{
		if (!p.points_ || !p.position_)
			return;
		if (p.kdtree_)
		{
			delete p.kdtree_;
			p.kdtree_ = nullptr;
		}
		std::vector<Vec3> points;
		points.reserve(nb_cells<PVertex>(*p.points_));
		p.kdtree_vertices_.clear();
		p.kdtree_vertices_.reserve(points.size());

		foreach_cell(*p.points_, [&](PVertex v) {
			uint32 idx = index_of(*p.points_, v);
			points.push_back((*p.position_)[idx]);
			p.kdtree_vertices_.push_back(v);
			return true;
		});

		p.kdtree_ = new acc::KDTree<3, uint32>(points);
	}

	bool pick_vertex_under_cursor(View* view, Parameters& p)
	{
		if (!view || !p.points_ || !p.position_)
			return false;

		int32 x = view->mouse_x();
		int32 y = view->mouse_y();

		rendering::GLVec3d near_ = view->unproject(x, y, 0.0);
		rendering::GLVec3d far_d = view->unproject(x, y, 1.0);
		Vec3 A{near_.x(), near_.y(), near_.z()};
		Vec3 B{far_d.x(), far_d.y(), far_d.z()};

		Scalar best_d2 = std::numeric_limits<Scalar>::max();
		PVertex best;
		foreach_cell(*p.points_, [&](PVertex v) -> bool {
			uint32 idx = index_of(*p.points_, v);
			const Vec3& pos = (*p.position_)[idx];
			Scalar d2 = geometry::squared_distance_line_point(A, B, pos);
			if (d2 < best_d2)
			{
				best_d2 = d2;
				best = v;
			}
			return true;
		});

		if (!best.is_valid())
			return false;

		p.picked_vertex_ = best;
		p.ball_finished_ = false;
		p.ball_initialized_ = false;
		update_highlights(p);
		return true;
	}

	void update_highlights(Parameters& p)
	{
		if (!p.color_ || !p.points_ || !points_provider_)
			return;

		const Vec4 white(1, 1, 1, 1);
		const Vec4 picked_color(1, 0.2f, 0.2f, 1);
		const Vec4 contact_color(0.2f, 1, 0.2f, 1);

		auto reset_vertex = [&](PVertex v) {
			if (!v.is_valid())
				return;
			uint32 idx = index_of(*p.points_, v);
			(*p.color_)[idx] = white;
		};

		if (p.last_picked_.is_valid() && p.last_picked_ != p.picked_vertex_ && p.last_picked_ != p.contact_vertex_)
			reset_vertex(p.last_picked_);
		if (p.last_contact_.is_valid() && p.last_contact_ != p.picked_vertex_ && p.last_contact_ != p.contact_vertex_)
			reset_vertex(p.last_contact_);

		if (p.picked_vertex_.is_valid())
			(*p.color_)[index_of(*p.points_, p.picked_vertex_)] = picked_color;
		if (p.contact_vertex_.is_valid() && p.contact_vertex_ != p.picked_vertex_)
			(*p.color_)[index_of(*p.points_, p.contact_vertex_)] = contact_color;

		p.last_picked_ = p.picked_vertex_;
		p.last_contact_ = p.contact_vertex_;

		points_provider_->emit_attribute_changed(*p.points_, p.color_.get());
	}

	bool init_ball_from_pick(Parameters& p)
	{
		if (!p.picked_vertex_.is_valid() || !p.position_ || !p.normal_)
			return false;
		if (!p.kdtree_)
			build_kdtree(p);
		if (!p.kdtree_)
			return false;

		const uint32 v_idx = index_of(*p.points_, p.picked_vertex_);
		const Vec3& pt = (*p.position_)[v_idx];
		Vec3 n = (*p.normal_)[v_idx];
		if (n.squaredNorm() <= Scalar(1e-12))
			return false;
		n.normalize();

		Scalar r = std::max<Scalar>(Scalar(0), Scalar(p.initial_radius_));
		p.ball_center_ = pt - (r * n);
		p.ball_radius_ = r;
		p.ball_prev_contact_pos_ = pt - (Scalar(2) * r * n);
		p.ball_prev_contact_vertex_ = PVertex();
		p.ball_iter_ = 0;
		p.ball_initialized_ = true;
		p.ball_finished_ = false;
		p.contact_vertex_ = PVertex();
		update_highlights(p);
		return true;
	}

	bool step_shrinking_ball(Parameters& p)
	{
		if (!p.picked_vertex_.is_valid() || !p.position_ || !p.normal_)
			return false;
		if (!p.ball_initialized_)
			if (!init_ball_from_pick(p))
				return false;
		if (p.ball_finished_)
			return false;
		if (!p.kdtree_)
			build_kdtree(p);
		if (!p.kdtree_)
			return false;

		const uint32 v_idx = index_of(*p.points_, p.picked_vertex_);
		const Vec3& pt = (*p.position_)[v_idx];
		Vec3 n = (*p.normal_)[v_idx];
		if (n.squaredNorm() <= Scalar(1e-12))
			return false;
		n.normalize();

		std::pair<uint32, Scalar> k_res;
		p.kdtree_->find_nn(p.ball_center_, &k_res);
		const Vec3 q_next = p.kdtree_->vertex(k_res.first);
		const Scalar d = k_res.second;
		const PVertex q_next_v = p.kdtree_vertices_[k_res.first];

		if (std::abs(d - p.ball_radius_) <= geometry::delta_convergence ||
			(q_next - p.ball_prev_contact_pos_).norm() < geometry::delta_convergence ||
			(q_next - pt).norm() < geometry::delta_convergence)
		{
			p.ball_finished_ = true;
			p.contact_vertex_ = q_next_v;
			update_highlights(p);
			update_sphere_visual(p);
			return false;
		}

		const Scalar r_next = geometry::compute_radius(pt, n, q_next);
		if (!std::isfinite(static_cast<double>(r_next)) || r_next <= Scalar(0))
		{
			p.ball_finished_ = true;
			update_sphere_visual(p);
			return false;
		}
		const Vec3 c_next = pt - (r_next * n);

		const Scalar separation_angle = geometry::angle(pt - c_next, q_next - c_next);
		if (p.ball_iter_ > 0 && separation_angle < geometry::denoise_preserve)
		{
			p.ball_finished_ = true;
			update_sphere_visual(p);
			return false;
		}

		p.ball_center_ = c_next;
		p.ball_radius_ = r_next;
		p.ball_prev_contact_pos_ = q_next;
		p.ball_prev_contact_vertex_ = q_next_v;
		p.contact_vertex_ = q_next_v;
		p.ball_iter_++;
		if (p.ball_iter_ > geometry::iteration_limit)
			p.ball_finished_ = true;

		update_highlights(p);
		update_sphere_visual(p);
		return true;
	}

	void ensure_sphere_mesh(Parameters& p)
	{
		if (!points_provider_ || !p.points_)
			return;

		if (!p.sphere_mesh_)
		{
			const std::string name = points_provider_->mesh_name(*p.points_) + "_shrinking_ball";
			p.sphere_mesh_ = points_provider_->has_mesh(name) ? points_provider_->mesh(name) : points_provider_->add_mesh(name);
			p.sphere_position_ = get_or_add_attribute<Vec3, PVertex>(*p.sphere_mesh_, "position");
			p.sphere_radius_ = get_or_add_attribute<Scalar, PVertex>(*p.sphere_mesh_, "radius");
			p.sphere_color_ = get_or_add_attribute<Vec4, PVertex>(*p.sphere_mesh_, "color");
		}

		if (nb_cells<PVertex>(*p.sphere_mesh_) == 0)
		{
			p.sphere_vertex_ = add_vertex(*p.sphere_mesh_);
			points_provider_->emit_connectivity_changed(*p.sphere_mesh_);
		}
		else
		{
			p.sphere_vertex_ = of_index<PVertex>(*p.sphere_mesh_, 0);
		}

		p.sphere_ready_ = p.sphere_vertex_.is_valid();

		if (pcr_ && selected_view_ && p.sphere_mesh_ && p.sphere_position_ && p.sphere_radius_ && p.sphere_color_)
		{
			pcr_->set_vertex_position(*selected_view_, *p.sphere_mesh_, p.sphere_position_);
			pcr_->set_vertex_radius(*selected_view_, *p.sphere_mesh_, p.sphere_radius_);
			pcr_->set_vertex_color(*selected_view_, *p.sphere_mesh_, p.sphere_color_);
			pcr_->set_vertex_color_per_cell(*selected_view_, *p.sphere_mesh_, PointCloudRender<POINTS>::PER_VERTEX);
		}
	}

	void update_sphere_visual(Parameters& p)
	{
		ensure_sphere_mesh(p);
		if (!p.sphere_ready_ || !p.sphere_position_ || !p.sphere_radius_ || !p.sphere_color_)
			return;

		const Vec4 sphere_color(0.2f, 0.6f, 1.0f, p.sphere_alpha_);
		const uint32 idx = index_of(*p.sphere_mesh_, p.sphere_vertex_);

		if (p.ball_initialized_)
		{
			(*p.sphere_position_)[idx] = p.ball_center_;
			(*p.sphere_radius_)[idx] = p.ball_radius_;
			(*p.sphere_color_)[idx] = sphere_color;
		}
		else
		{
			(*p.sphere_radius_)[idx] = Scalar(0);
			(*p.sphere_color_)[idx] = sphere_color;
		}

		if (points_provider_)
		{
			points_provider_->emit_attribute_changed(*p.sphere_mesh_, p.sphere_position_.get());
			points_provider_->emit_attribute_changed(*p.sphere_mesh_, p.sphere_radius_.get());
			points_provider_->emit_attribute_changed(*p.sphere_mesh_, p.sphere_color_.get());
		}
	}

private:
	ui::MeshProvider<POINTS>* points_provider_ = nullptr;
	PointCloudRender<POINTS>* pcr_ = nullptr;
	View* selected_view_ = nullptr;
	POINTS* selected_points_ = nullptr;
	CenterUDFQuery center_udf_query_;
	std::unordered_map<POINTS*, Parameters> parameters_;
};

} // namespace ui

} // namespace cgogn

#endif // CGOGN_MODULE_SHRINKING_BALL_DEBUGGER_H_
