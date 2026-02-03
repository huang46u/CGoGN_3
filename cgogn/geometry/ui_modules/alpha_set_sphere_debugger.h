#ifndef CGOGN_MODULE_ALPHA_SET_SPHERE_DEBUGGER_H_
#define CGOGN_MODULE_ALPHA_SET_SPHERE_DEBUGGER_H_

#include <cgogn/ui/app.h>
#include <cgogn/ui/imgui_helpers.h>
#include <cgogn/ui/module.h>
#include <cgogn/ui/view.h>

#include <cgogn/core/functions/attributes.h>
#include <cgogn/core/functions/traversals/vertex.h>
#include <cgogn/core/ui_modules/mesh_provider.h>

#include <cgogn/geometry/functions/distance.h>
#include <cgogn/geometry/types/line_quadric.h>
#include <cgogn/geometry/types/spherical_quadric.h>
#include <cgogn/geometry/types/vector_traits.h>

#include <cgogn/rendering/ui_modules/point_cloud_render.h>

#include <GLFW/glfw3.h>

#include <Eigen/Dense>

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <mutex>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>

namespace cgogn
{

namespace ui
{

using geometry::Line_Quadric;
using geometry::Mat3;
using geometry::Mat4;
using geometry::Scalar;
using geometry::Spherical_Quadric;
using geometry::Vec3;
using geometry::Vec4;

template <typename POINTS>
class AlphaSetSphereDebugger : public ViewModule
{
public:
	template <typename T>
	using PAttribute = typename mesh_traits<POINTS>::template Attribute<T>;
	using PVertex = typename mesh_traits<POINTS>::Vertex;

	AlphaSetSphereDebugger(const App& app)
		: ViewModule(app, "AlphaSetSphereDebugger (" + std::string{mesh_traits<POINTS>::name} + ")"),
		  selected_view_(app.current_view())
	{
	}

private:
	struct Parameters
	{
		POINTS* points_ = nullptr;
		POINTS* samples_mesh_ = nullptr;
		POINTS* alpha_inside_mesh_ = nullptr;
		POINTS* spheres_mesh_ = nullptr;
		POINTS* debug_cluster_mesh_ = nullptr;
		POINTS* debug_sphere_mesh_ = nullptr;

		std::shared_ptr<PAttribute<Vec3>> samples_position_ = nullptr;
		std::shared_ptr<PAttribute<Vec3>> samples_normal_ = nullptr;
		std::shared_ptr<PAttribute<Scalar>> samples_area_ = nullptr;
		std::shared_ptr<PAttribute<std::vector<PVertex>>> samples_knn_ = nullptr;
		std::shared_ptr<PAttribute<Spherical_Quadric>> samples_quadric_ = nullptr;
		std::shared_ptr<PAttribute<Line_Quadric>> samples_line_quadric_ = nullptr;

		std::shared_ptr<PAttribute<Vec3>> alpha_inside_position_ = nullptr;
		std::shared_ptr<PAttribute<Vec4>> alpha_inside_color_ = nullptr;
		std::shared_ptr<PAttribute<Vec3>> alpha_inside_projected_position_ = nullptr;
		std::shared_ptr<PAttribute<Vec3>> alpha_inside_projected_normal_ = nullptr;
		std::shared_ptr<PAttribute<Scalar>> alpha_inside_projected_area_ = nullptr;
		std::shared_ptr<PAttribute<std::vector<PVertex>>> alpha_inside_projected_knn_ = nullptr;
		std::shared_ptr<PAttribute<Spherical_Quadric>> alpha_inside_projected_quadric_ = nullptr;
		std::shared_ptr<PAttribute<Line_Quadric>> alpha_inside_projected_line_quadric_ = nullptr;

		std::shared_ptr<PAttribute<Vec3>> spheres_position_ = nullptr;
		std::shared_ptr<PAttribute<Scalar>> spheres_radius_ = nullptr;
		std::shared_ptr<PAttribute<Vec4>> spheres_cluster_color_ = nullptr;

		std::shared_ptr<PAttribute<Vec3>> debug_cluster_position_ = nullptr;
		std::shared_ptr<PAttribute<Vec3>> debug_cluster_projected_position_ = nullptr;
		std::shared_ptr<PAttribute<Vec4>> debug_cluster_color_ = nullptr;
		std::shared_ptr<PAttribute<Vec3>> debug_sphere_position_ = nullptr;
		std::shared_ptr<PAttribute<Scalar>> debug_sphere_radius_ = nullptr;
		std::shared_ptr<PAttribute<Vec4>> debug_sphere_color_ = nullptr;

		enum SourceMode
		{
			SOURCE_SAMPLES,
			SOURCE_ALPHA_INSIDE
		};
		SourceMode source_mode_ = SOURCE_ALPHA_INSIDE;
		std::vector<std::vector<PVertex>> clusters_;
		std::vector<uint32> cluster_component_counts_;
		std::vector<uint32> cluster_component_totals_;
		std::vector<std::vector<uint32>> cluster_component_sizes_;
		uint32 cluster_points_count_ = 0;
		uint32 cluster_spheres_count_ = 0;
		uint32 cluster_spheres_capacity_ = 0;

		PVertex picked_sphere_;
		bool show_normals_ = false;
		bool solo_view_ = false;
		bool solo_active_ = false;
		std::unordered_map<const POINTS*, bool> render_backup_;
		float32 sqem_update_lambda_ = 0.20f;
		float32 alpha_ = 0.005f;
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
		if (key_code != GLFW_KEY_O)
			return;

		Parameters& p = parameters_[selected_points_];
		if (!p.spheres_mesh_ || !p.spheres_position_)
			return;

		int32 x = view->mouse_x();
		int32 y = view->mouse_y();

		rendering::GLVec3d near_ = view->unproject(x, y, 0.0);
		rendering::GLVec3d far_d = view->unproject(x, y, 1.0);
		Vec3 A{near_.x(), near_.y(), near_.z()};
		Vec3 B{far_d.x(), far_d.y(), far_d.z()};

		Vec3 picked_sphere_center;
		bool first = true;
		foreach_cell(*p.spheres_mesh_, [&](PVertex v) -> bool {
			const Vec3& sp = (*p.spheres_position_)[index_of(*p.spheres_mesh_, v)];
			if (first)
			{
				first = false;
				p.picked_sphere_ = v;
				picked_sphere_center = sp;
				return true;
			}
			if (geometry::squared_distance_line_point(A, B, sp) <
				geometry::squared_distance_line_point(A, B, picked_sphere_center))
			{
				p.picked_sphere_ = v;
				picked_sphere_center = sp;
			}
			return true;
		});

		if (!p.picked_sphere_.is_valid())
			return;

		const uint32 s_idx = index_of(*p.spheres_mesh_, p.picked_sphere_);
		SourceBasic src;
		bool needs_cluster = true;
		if (get_source_basic(p, src))
		{
			const uint32 nb_spheres = nb_cells<PVertex>(*p.spheres_mesh_);
			const uint32 cluster_size =
				p.spheres_position_ ? static_cast<uint32>(p.spheres_position_->size()) : nb_spheres;
			const uint32 nb_points = nb_cells<PVertex>(*src.mesh);
			needs_cluster = (p.cluster_spheres_count_ != nb_spheres) ||
							(p.cluster_spheres_capacity_ != cluster_size) ||
							(p.cluster_points_count_ != nb_points);
			if (!needs_cluster && s_idx < p.clusters_.size())
				needs_cluster = p.clusters_[s_idx].empty();
		}
		if (needs_cluster)
			compute_power_cluster_source(p);
		update_debug_meshes(p, p.picked_sphere_);
		p.solo_view_ = true;
		set_solo_view(p, true);
		view->request_update();
	}

	void left_panel() override
	{
		if (!points_provider_)
			return;

		if (app_.nb_views() > 1)
			imgui_view_selector(this, selected_view_, [&](View* v) { selected_view_ = v; });
		else
			selected_view_ = app_.current_view();

		imgui_mesh_selector(points_provider_, selected_points_, "Cmap0",
							[&](POINTS& m) { set_selected_points(m); });

		if (!selected_points_)
			return;

		auto it = parameters_.find(selected_points_);
		if (it == parameters_.end())
			return;
		Parameters& p = parameters_[selected_points_];
		const bool has_samples = p.samples_mesh_ && p.samples_position_;
		const bool has_alpha_inside = p.alpha_inside_mesh_ && p.alpha_inside_projected_position_;
		const bool has_spheres = p.spheres_mesh_ && p.spheres_position_ && p.spheres_radius_;

		if (!has_spheres)
		{
			ImGui::TextColored(ImVec4(1, 1, 0, 1), "Spheres mesh not found.");
			return;
		}
		if (p.source_mode_ == Parameters::SOURCE_SAMPLES && !has_samples)
		{
			ImGui::TextColored(ImVec4(1, 1, 0, 1), "Samples mesh not found.");
			return;
		}
		if (p.source_mode_ == Parameters::SOURCE_ALPHA_INSIDE && !has_alpha_inside)
		{
			ImGui::TextColored(ImVec4(1, 1, 0, 1), "Alpha inside projected data not found.");
			return;
		}

		ImGui::Checkbox("Solo Debug View", &p.solo_view_);
		set_solo_view(p, p.solo_view_);

		ImGui::Separator();
		ImGui::Text("Source: %s", p.source_mode_ == Parameters::SOURCE_SAMPLES ? "Samples" : "Alpha Inside (Projected)");
		if (ImGui::Button("Power Cluster"))
		{
			compute_power_cluster_source(p);
			update_debug_meshes(p, p.picked_sphere_);
		}
		if (p.source_mode_ == Parameters::SOURCE_SAMPLES)
			ImGui::Text("Sample points: %zu", nb_cells<PVertex>(*p.samples_mesh_));
		else
			ImGui::Text("Alpha inside points: %zu", nb_cells<PVertex>(*p.alpha_inside_mesh_));

		ImGui::Separator();
		ImGui::Text("Press O to pick a sphere and show its cluster.");
		if (p.picked_sphere_.is_valid())
		{
			const uint32 s_idx = index_of(*p.spheres_mesh_, p.picked_sphere_);
			const Vec3& sp = (*p.spheres_position_)[s_idx];
			ImGui::Text("Picked sphere center: (%f, %f, %f)", sp[0], sp[1], sp[2]);
			ImGui::Text("Picked sphere radius: %f", (*p.spheres_radius_)[s_idx]);
			ImGui::Text("Picked sphere index: %u", s_idx);
			if (s_idx < p.clusters_.size())
			{
				const auto& cluster = p.clusters_[s_idx];
				ImGui::Text("Cluster size: %zu", cluster.size());
			}
			if (s_idx < p.cluster_component_counts_.size())
			{
				ImGui::Text("Components kept (>3): %u", p.cluster_component_counts_[s_idx]);
				if (s_idx < p.cluster_component_totals_.size())
					ImGui::Text("Components total: %u", p.cluster_component_totals_[s_idx]);
				if (s_idx < p.cluster_component_sizes_.size())
				{
					const auto& sizes = p.cluster_component_sizes_[s_idx];
					if (!sizes.empty())
					{
						std::string size_str;
						for (size_t i = 0; i < sizes.size(); ++i)
						{
							if (i > 0)
								size_str += ", ";
							size_str += std::to_string(sizes[i]);
						}
						ImGui::Text("Component sizes: %s", size_str.c_str());
					}
				}
			}
			SourceData src;
			const bool src_ready = get_source_data(p, src);
			ImGui::Text("Source data: %s", src_ready ? "ready" : "not ready");
		}

		if (ImGui::Button("Show Cluster Colors"))
		{
			p.show_normals_ = false;
			update_debug_meshes(p, p.picked_sphere_);
		}
		ImGui::SameLine();
		if (ImGui::Button("Show Cluster Normals"))
		{
			p.show_normals_ = true;
			update_debug_meshes(p, p.picked_sphere_);
		}

		ImGui::Separator();
		ImGui::InputFloat("Alpha (fixed radius)", &p.alpha_, 0.0f, 0.0f, "%.6f");
		if (p.alpha_ < 0.0f)
			p.alpha_ = 0.0f;
		ImGui::SliderFloat("Update lambda", &p.sqem_update_lambda_, 0.0f, 2.0f, "%.6f");
		SourceData src;
		const bool src_ready = get_source_data(p, src);
		if (!src_ready)
			ImGui::TextColored(ImVec4(1, 1, 0, 1), "Source data not ready.");

		if (ImGui::Button("Update SQEM (one-step)"))
		{
			if (p.picked_sphere_.is_valid() && src_ready)
			{
				update_sphere_sqem_one_step(p, p.picked_sphere_);
				update_debug_meshes(p, p.picked_sphere_);
			}
		}
		if (ImGui::Button("Update SQEM + Euclidean (iter)"))
		{
			if (p.picked_sphere_.is_valid() && src_ready)
			{
				update_sphere_sqem_euclidean(p, p.picked_sphere_);
				update_debug_meshes(p, p.picked_sphere_);
			}
		}
		if (ImGui::Button("Update SQEM + Line Quadric (fixed alpha)"))
		{
			if (p.picked_sphere_.is_valid() && src_ready)
			{
				update_sphere_sqem_line_quadric_fix_alpha(p, p.picked_sphere_);
				update_debug_meshes(p, p.picked_sphere_);
			}
		}
		if (ImGui::Button("Update SQEM + Line Quadric (free radius)"))
		{
			if (p.picked_sphere_.is_valid() && src_ready)
			{
				update_sphere_sqem_line_quadric_free_radius(p, p.picked_sphere_);
				update_debug_meshes(p, p.picked_sphere_);
			}
		}
	}

private:
	void set_selected_points(POINTS& m)
	{
		if (selected_points_)
		{
			Parameters& prev = parameters_[selected_points_];
			if (prev.solo_active_)
				set_solo_view(prev, false);
		}
		selected_points_ = &m;
		Parameters& p = parameters_[selected_points_];
		init_parameters(m, p);
		p.solo_view_ = false;
		set_solo_view(p, false);
	}

	void init_parameters(POINTS& m, Parameters& p)
	{
		p.points_ = &m;
		const std::string selected_name = points_provider_->mesh_name(m);
		auto ends_with = [](const std::string& s, const std::string& suffix) {
			return s.size() >= suffix.size() && s.compare(s.size() - suffix.size(), suffix.size(), suffix) == 0;
		};

		std::string base_name = selected_name;
		if (ends_with(selected_name, "_samples"))
		{
			base_name = selected_name.substr(0, selected_name.size() - std::string("_samples").size());
			p.source_mode_ = Parameters::SOURCE_SAMPLES;
		}
		else if (ends_with(selected_name, "_alpha_inside"))
		{
			base_name = selected_name.substr(0, selected_name.size() - std::string("_alpha_inside").size());
			p.source_mode_ = Parameters::SOURCE_ALPHA_INSIDE;
		}
		else
		{
			const std::string guess_samples = base_name + "_samples";
			if (points_provider_->has_mesh(guess_samples))
				p.source_mode_ = Parameters::SOURCE_SAMPLES;
			else
				p.source_mode_ = Parameters::SOURCE_ALPHA_INSIDE;
		}

		const std::string samples_name = base_name + "_samples";
		p.samples_mesh_ = points_provider_->has_mesh(samples_name) ? points_provider_->mesh(samples_name) : nullptr;

		const std::string inside_name = base_name + "_alpha_inside";
		p.alpha_inside_mesh_ =
			points_provider_->has_mesh(inside_name) ? points_provider_->mesh(inside_name) : nullptr;

		const std::string spheres_name = base_name + "_spheres";
		p.spheres_mesh_ =
			points_provider_->has_mesh(spheres_name) ? points_provider_->mesh(spheres_name) : nullptr;

		const std::string debug_cluster_name = base_name + "_alpha_inside_debug";
		if (!p.debug_cluster_mesh_)
			p.debug_cluster_mesh_ = points_provider_->has_mesh(debug_cluster_name)
										? points_provider_->mesh(debug_cluster_name)
										: points_provider_->add_mesh(debug_cluster_name);

		const std::string debug_sphere_name = base_name + "_alpha_inside_debug_sphere";
		if (!p.debug_sphere_mesh_)
			p.debug_sphere_mesh_ = points_provider_->has_mesh(debug_sphere_name)
									   ? points_provider_->mesh(debug_sphere_name)
									   : points_provider_->add_mesh(debug_sphere_name);

		if (p.samples_mesh_)
		{
			p.samples_position_ = get_attribute<Vec3, PVertex>(*p.samples_mesh_, "position");
			p.samples_normal_ = get_attribute<Vec3, PVertex>(*p.samples_mesh_, "normal");
			p.samples_area_ = get_attribute<Scalar, PVertex>(*p.samples_mesh_, "area");
			p.samples_knn_ = get_attribute<std::vector<PVertex>, PVertex>(*p.samples_mesh_, "knn");
			p.samples_quadric_ = get_attribute<Spherical_Quadric, PVertex>(*p.samples_mesh_, "quadric");
			p.samples_line_quadric_ = get_attribute<Line_Quadric, PVertex>(*p.samples_mesh_, "line_quadric");
		}

		if (p.alpha_inside_mesh_)
		{
			p.alpha_inside_position_ = get_attribute<Vec3, PVertex>(*p.alpha_inside_mesh_, "position");
			p.alpha_inside_color_ = get_or_add_attribute<Vec4, PVertex>(*p.alpha_inside_mesh_, "color");
			p.alpha_inside_projected_position_ =
				get_attribute<Vec3, PVertex>(*p.alpha_inside_mesh_, "projected_alpha_inside_position");
			p.alpha_inside_projected_normal_ =
				get_attribute<Vec3, PVertex>(*p.alpha_inside_mesh_, "projected_alpha_inside_normal");
			p.alpha_inside_projected_area_ =
				get_attribute<Scalar, PVertex>(*p.alpha_inside_mesh_, "projected_alpha_inside_area");
			p.alpha_inside_projected_knn_ =
				get_attribute<std::vector<PVertex>, PVertex>(*p.alpha_inside_mesh_, "projected_alpha_inside_knn");
			p.alpha_inside_projected_quadric_ =
				get_attribute<Spherical_Quadric, PVertex>(*p.alpha_inside_mesh_, "projected_alpha_inside_quadric");
			p.alpha_inside_projected_line_quadric_ =
				get_attribute<Line_Quadric, PVertex>(*p.alpha_inside_mesh_, "projected_alpha_inside_line_quadric");
		}

		if (p.spheres_mesh_)
		{
			p.spheres_position_ = get_attribute<Vec3, PVertex>(*p.spheres_mesh_, "position");
			p.spheres_radius_ = get_attribute<Scalar, PVertex>(*p.spheres_mesh_, "radius");
			p.spheres_cluster_color_ = get_or_add_attribute<Vec4, PVertex>(*p.spheres_mesh_, "cluster_color");
		}

		if (p.debug_cluster_mesh_)
		{
			p.debug_cluster_position_ = get_or_add_attribute<Vec3, PVertex>(*p.debug_cluster_mesh_, "position");
			p.debug_cluster_projected_position_ =
				get_or_add_attribute<Vec3, PVertex>(*p.debug_cluster_mesh_, "projected_position");
			p.debug_cluster_color_ = get_or_add_attribute<Vec4, PVertex>(*p.debug_cluster_mesh_, "color");
		}

		if (p.debug_sphere_mesh_)
		{
			p.debug_sphere_position_ = get_or_add_attribute<Vec3, PVertex>(*p.debug_sphere_mesh_, "position");
			p.debug_sphere_radius_ = get_or_add_attribute<Scalar, PVertex>(*p.debug_sphere_mesh_, "radius");
			p.debug_sphere_color_ = get_or_add_attribute<Vec4, PVertex>(*p.debug_sphere_mesh_, "color");
		}

		setup_debug_render(p);
	}

	void setup_debug_render(Parameters& p)
	{
		if (!pcr_ || !selected_view_)
			return;

		if (p.samples_mesh_ && p.samples_position_)
		{
			pcr_->set_vertex_position(*selected_view_, *p.samples_mesh_, p.samples_position_);
			if(p.source_mode_ == Parameters::SOURCE_SAMPLES)
				pcr_->set_render_vertices(*selected_view_, *p.samples_mesh_, true);
		}

		if (p.alpha_inside_mesh_ && p.alpha_inside_position_ && p.alpha_inside_color_)
		{
			pcr_->set_vertex_position(*selected_view_, *p.alpha_inside_mesh_, p.alpha_inside_position_);
			pcr_->set_vertex_color(*selected_view_, *p.alpha_inside_mesh_, p.alpha_inside_color_);
			pcr_->set_vertex_color_per_cell(*selected_view_, *p.alpha_inside_mesh_,
											PointCloudRender<POINTS>::PER_VERTEX);
			if (p.source_mode_ == Parameters::SOURCE_ALPHA_INSIDE)
				pcr_->set_render_vertices(*selected_view_, *p.alpha_inside_mesh_, true);
		}

		if (p.spheres_mesh_ && p.spheres_position_ && p.spheres_radius_ && p.spheres_cluster_color_)
		{
			pcr_->set_vertex_position(*selected_view_, *p.spheres_mesh_, p.spheres_position_);
			pcr_->set_vertex_radius(*selected_view_, *p.spheres_mesh_, p.spheres_radius_);
			pcr_->set_vertex_color(*selected_view_, *p.spheres_mesh_, p.spheres_cluster_color_);
			pcr_->set_vertex_color_per_cell(*selected_view_, *p.spheres_mesh_, PointCloudRender<POINTS>::PER_VERTEX);
			pcr_->set_render_vertices(*selected_view_, *p.spheres_mesh_, true);
		}

		if (p.debug_cluster_mesh_ && p.debug_cluster_position_ && p.debug_cluster_color_)
		{
			auto pos_attr = p.debug_cluster_projected_position_ ? p.debug_cluster_projected_position_ : p.debug_cluster_position_;
			pcr_->set_vertex_position(*selected_view_, *p.debug_cluster_mesh_, pos_attr);
			pcr_->set_vertex_color(*selected_view_, *p.debug_cluster_mesh_, p.debug_cluster_color_);
			pcr_->set_vertex_color_per_cell(*selected_view_, *p.debug_cluster_mesh_,
											PointCloudRender<POINTS>::PER_VERTEX);
			pcr_->set_render_vertices(*selected_view_, *p.debug_cluster_mesh_, true);
		}

		if (p.debug_sphere_mesh_ && p.debug_sphere_position_ && p.debug_sphere_radius_ && p.debug_sphere_color_)
		{
			pcr_->set_vertex_position(*selected_view_, *p.debug_sphere_mesh_, p.debug_sphere_position_);
			pcr_->set_vertex_radius(*selected_view_, *p.debug_sphere_mesh_, p.debug_sphere_radius_);
			pcr_->set_vertex_color(*selected_view_, *p.debug_sphere_mesh_, p.debug_sphere_color_);
			pcr_->set_vertex_color_per_cell(*selected_view_, *p.debug_sphere_mesh_,
											PointCloudRender<POINTS>::PER_VERTEX);
			pcr_->set_render_vertices(*selected_view_, *p.debug_sphere_mesh_, true);
		}
	}

	void set_solo_view(Parameters& p, bool enable)
	{
		if (!pcr_ || !selected_view_)
			return;
		if (enable == p.solo_active_)
			return;

		if (enable)
		{
			p.render_backup_.clear();
			points_provider_->foreach_mesh([&](POINTS& m, const std::string&) {
				const bool current = pcr_->render_vertices(*selected_view_, m);
				p.render_backup_[&m] = current;
				const bool keep = (&m == p.debug_cluster_mesh_ || &m == p.debug_sphere_mesh_);
					pcr_->set_render_vertices(*selected_view_, m, keep);
			});
			p.solo_active_ = true;
		}
		else
		{
			for (const auto& it : p.render_backup_)
			{
				if (it.first)
					pcr_->set_render_vertices(*selected_view_, *it.first, it.second);
			}
			p.render_backup_.clear();
			p.solo_active_ = false;
		}
	}

	struct SourceData
	{
		POINTS* mesh = nullptr;
		std::shared_ptr<PAttribute<Vec3>> position = nullptr;
		std::shared_ptr<PAttribute<Vec3>> projected_position = nullptr;
		std::shared_ptr<PAttribute<Vec3>> cluster_position = nullptr;
		std::shared_ptr<PAttribute<Vec3>> normal = nullptr;
		std::shared_ptr<PAttribute<Scalar>> area = nullptr;
		std::shared_ptr<PAttribute<std::vector<PVertex>>> knn = nullptr;
		std::shared_ptr<PAttribute<Spherical_Quadric>> quadric = nullptr;
		std::shared_ptr<PAttribute<Line_Quadric>> line_quadric = nullptr;
	};

	struct SourceBasic
	{
		POINTS* mesh = nullptr;
		std::shared_ptr<PAttribute<Vec3>> position = nullptr;
		std::shared_ptr<PAttribute<Vec3>> projected_position = nullptr;
		std::shared_ptr<PAttribute<Vec3>> cluster_position = nullptr;
		std::shared_ptr<PAttribute<Vec3>> normal = nullptr;
	};

	bool get_source_basic(Parameters& p, SourceBasic& out)
	{
		if (p.source_mode_ == Parameters::SOURCE_SAMPLES)
		{
			out.mesh = p.samples_mesh_;
			out.position = p.samples_position_;
			out.projected_position = p.samples_position_;
			out.cluster_position = p.samples_position_;
			out.normal = p.samples_normal_;
		}
		else
		{
			out.mesh = p.alpha_inside_mesh_;
			out.position = p.alpha_inside_position_;
			out.projected_position = p.alpha_inside_projected_position_;
			out.cluster_position = p.alpha_inside_position_;
			out.normal = p.alpha_inside_projected_normal_;
		}
		if (!out.mesh || !out.cluster_position || !out.projected_position)
			return false;
		return true;
	}

	bool get_source_data(Parameters& p, SourceData& out)
	{
		if (p.source_mode_ == Parameters::SOURCE_SAMPLES)
		{
			out.mesh = p.samples_mesh_;
			out.position = p.samples_position_;
			out.projected_position = p.samples_position_;
			out.cluster_position = p.samples_position_;
			out.normal = p.samples_normal_;
			out.area = p.samples_area_;
			out.knn = p.samples_knn_;
			out.quadric = p.samples_quadric_;
			out.line_quadric = p.samples_line_quadric_;
		}
		else
		{
			out.mesh = p.alpha_inside_mesh_;
			out.position = p.alpha_inside_projected_position_;
			out.projected_position = p.alpha_inside_projected_position_;
			out.cluster_position = p.alpha_inside_position_;
			out.normal = p.alpha_inside_projected_normal_;
			out.area = p.alpha_inside_projected_area_;
			out.knn = p.alpha_inside_projected_knn_;
			out.quadric = p.alpha_inside_projected_quadric_;
			out.line_quadric = p.alpha_inside_projected_line_quadric_;
		}

		if (!out.mesh || !out.position || !out.projected_position || !out.cluster_position || !out.normal || !out.area ||
			!out.knn || !out.quadric || !out.line_quadric)
			return false;
		return true;
	}

	void filter_cluster_components(Parameters& p, const SourceData& src, std::vector<PVertex>& cluster,
								   uint32& kept_components, uint32& total_components,
								   std::vector<uint32>& kept_sizes)
	{
		kept_components = 0;
		total_components = 0;
		kept_sizes.clear();
		if (cluster.empty())
			return;

		std::unordered_set<uint32> in_cluster;
		in_cluster.reserve(cluster.size() * 2);
		for (PVertex v : cluster)
			in_cluster.insert(index_of(*src.mesh, v));

		std::unordered_set<uint32> visited;
		visited.reserve(cluster.size() * 2);

		std::vector<PVertex> kept;
		kept.reserve(cluster.size());

		std::vector<PVertex> stack;
		std::vector<PVertex> component;
		for (PVertex v : cluster)
		{
			uint32 v_idx = index_of(*src.mesh, v);
			if (visited.find(v_idx) != visited.end())
				continue;

			total_components++;
			stack.clear();
			component.clear();
			stack.push_back(v);
			visited.insert(v_idx);

			while (!stack.empty())
			{
				PVertex cur = stack.back();
				stack.pop_back();
				component.push_back(cur);

				uint32 cur_idx = index_of(*src.mesh, cur);
				for (PVertex nb : (*src.knn)[cur_idx])
				{
					uint32 nb_idx = index_of(*src.mesh, nb);
					if (in_cluster.find(nb_idx) == in_cluster.end())
						continue;
					if (visited.insert(nb_idx).second)
						stack.push_back(nb);
				}
			}

			if (component.size() > 3)
			{
				kept_components++;
				kept_sizes.push_back(static_cast<uint32>(component.size()));
				kept.insert(kept.end(), component.begin(), component.end());
			}
		}

		cluster.swap(kept);
	}

	void compute_power_cluster_source(Parameters& p)
	{
		SourceBasic src;
		if (!get_source_basic(p, src))
			return;
		if (!p.spheres_mesh_ || !p.spheres_position_ || !p.spheres_radius_)
			return;

		const uint32 nb_spheres = nb_cells<PVertex>(*p.spheres_mesh_);
		const uint32 cluster_size =
			p.spheres_position_ ? static_cast<uint32>(p.spheres_position_->size()) : nb_spheres;
		if (p.clusters_.size() != cluster_size)
			p.clusters_.assign(cluster_size, {});
		else
		{
			for (auto& c : p.clusters_)
				c.clear();
		}

		if (nb_spheres == 0)
			return;

		parallel_foreach_cell(*src.mesh, [&](PVertex v) -> bool {
			uint32 v_index = index_of(*src.mesh, v);
			const Vec3& vp = (*src.cluster_position)[v_index];

			Scalar min_power_distance = std::numeric_limits<Scalar>::max();
			PVertex closest_sphere;
			uint32 closest_sphere_index = 0;

			foreach_cell(*p.spheres_mesh_, [&](PVertex pv) {
				uint32 pv_index = index_of(*p.spheres_mesh_, pv);
				const Vec3& center = (*p.spheres_position_)[pv_index];
				Scalar radius = (*p.spheres_radius_)[pv_index];
				Scalar dist_sq = (vp - center).squaredNorm();
				Scalar power_dist = dist_sq - radius * radius;

				if (power_dist < min_power_distance)
				{
					min_power_distance = power_dist;
					closest_sphere = pv;
					closest_sphere_index = pv_index;
				}
				return true;
			});

			std::lock_guard<std::mutex> lock(spheres_mutex_[closest_sphere_index % spheres_mutex_.size()]);
			p.clusters_[closest_sphere_index].push_back(v);
			return true;
		});

		SourceData src_full;
		if (!get_source_data(p, src_full))
			return;

		p.cluster_component_counts_.assign(cluster_size, 0);
		p.cluster_component_totals_.assign(cluster_size, 0);
		p.cluster_component_sizes_.assign(cluster_size, {});
		for (uint32 i = 0; i < cluster_size; ++i)
		{
			uint32 kept = 0;
			uint32 total = 0;
			filter_cluster_components(p, src_full, p.clusters_[i], kept, total, p.cluster_component_sizes_[i]);
			p.cluster_component_counts_[i] = kept;
			p.cluster_component_totals_[i] = total;
		}
		p.cluster_spheres_count_ = nb_spheres;
		p.cluster_spheres_capacity_ = cluster_size;
		p.cluster_points_count_ = nb_cells<PVertex>(*src.mesh);
	}


	void update_debug_meshes(Parameters& p, PVertex sphere)
	{
		if (!p.debug_cluster_mesh_ || !p.debug_cluster_position_ || !p.debug_cluster_color_)
			return;
		if (!p.debug_sphere_mesh_ || !p.debug_sphere_position_ || !p.debug_sphere_radius_ || !p.debug_sphere_color_)
			return;
		if (!p.spheres_position_ || !p.spheres_radius_)
			return;
		SourceBasic src;
		if (!get_source_basic(p, src))
			return;

		points_provider_->clear_mesh(*p.debug_cluster_mesh_);
		points_provider_->clear_mesh(*p.debug_sphere_mesh_);

		if (!sphere.is_valid())
			return;

		const Vec4 base_color = p.spheres_cluster_color_
									? value<Vec4>(*p.spheres_mesh_, p.spheres_cluster_color_, sphere)
									: Vec4(0.2f, 0.8f, 0.9f, 1.0f);

		const uint32 s_idx = index_of(*p.spheres_mesh_, sphere);
		if (s_idx >= p.clusters_.size())
			return;
		const std::vector<PVertex>& cluster = p.clusters_[s_idx];
		for (PVertex v : cluster)
		{
			uint32 v_idx = index_of(*src.mesh, v);
			PVertex dv = add_vertex(*p.debug_cluster_mesh_);
			uint32 dv_idx = index_of(*p.debug_cluster_mesh_, dv);
			(*p.debug_cluster_position_)[dv_idx] = (*src.cluster_position)[v_idx];
			if (p.debug_cluster_projected_position_)
				(*p.debug_cluster_projected_position_)[dv_idx] = (*src.projected_position)[v_idx];

			Vec4 c = base_color;
			if (p.show_normals_ && src.normal)
			{
				const Vec3& n = (*src.normal)[v_idx];
				c = Vec4((n.x() + 1.0f) * 0.5f, (n.y() + 1.0f) * 0.5f, (n.z() + 1.0f) * 0.5f, 1.0f);
			}
			(*p.debug_cluster_color_)[dv_idx] = c;
		}

		PVertex sv = add_vertex(*p.debug_sphere_mesh_);
		uint32 sv_idx = index_of(*p.debug_sphere_mesh_, sv);
		const Vec3& sp = (*p.spheres_position_)[s_idx];
		(*p.debug_sphere_position_)[sv_idx] = sp;
		(*p.debug_sphere_radius_)[sv_idx] = (*p.spheres_radius_)[s_idx];
		(*p.debug_sphere_color_)[sv_idx] = base_color;

		points_provider_->emit_connectivity_changed(*p.debug_cluster_mesh_);
		points_provider_->emit_connectivity_changed(*p.debug_sphere_mesh_);
		points_provider_->emit_attribute_changed(*p.debug_cluster_mesh_, p.debug_cluster_position_.get());
		if (p.debug_cluster_projected_position_)
			points_provider_->emit_attribute_changed(*p.debug_cluster_mesh_, p.debug_cluster_projected_position_.get());
		points_provider_->emit_attribute_changed(*p.debug_cluster_mesh_, p.debug_cluster_color_.get());
		points_provider_->emit_attribute_changed(*p.debug_sphere_mesh_, p.debug_sphere_position_.get());
		points_provider_->emit_attribute_changed(*p.debug_sphere_mesh_, p.debug_sphere_radius_.get());
		points_provider_->emit_attribute_changed(*p.debug_sphere_mesh_, p.debug_sphere_color_.get());
	}

	void update_sphere_sqem_one_step(Parameters& p, PVertex sphere)
	{
		SourceData src;
		if (!get_source_data(p, src))
			return;
		if (!p.spheres_position_ || !p.spheres_radius_)
			return;
		const uint32 s_idx = index_of(*p.spheres_mesh_, sphere);
		if (s_idx >= p.clusters_.size() || p.clusters_[s_idx].empty())
			compute_power_cluster_source(p);
		if (s_idx >= p.clusters_.size())
			return;
		const std::vector<PVertex>& cluster = p.clusters_[s_idx];
		if (cluster.empty())
			return;

		Spherical_Quadric q;
		Scalar weight_sum = Scalar(0);
		for (PVertex v : cluster)
		{
			uint32 v_index = index_of(*src.mesh, v);
			Scalar weight = value<Scalar>(*src.mesh, src.area, v);
			if (weight <= 0.0)
				continue;
			q += (*src.quadric)[v_index] * weight;
			weight_sum += weight;
		}

		if (weight_sum <= Scalar(0))
			return;

		const Vec4 s = q._A.completeOrthogonalDecomposition().solve(q._b);
		if (!s.allFinite())
			return;

		(*p.spheres_position_)[s_idx] = s.head<3>();
		(*p.spheres_radius_)[s_idx] = s[3];
		if (points_provider_)
		{
			points_provider_->emit_attribute_changed(*p.spheres_mesh_, p.spheres_position_.get());
			points_provider_->emit_attribute_changed(*p.spheres_mesh_, p.spheres_radius_.get());
		}
	}

	void update_sphere_sqem_euclidean(Parameters& p, PVertex sphere)
	{
		SourceData src;
		if (!get_source_data(p, src))
			return;
		if (!p.spheres_position_ || !p.spheres_radius_)
			return;
		const uint32 s_idx = index_of(*p.spheres_mesh_, sphere);
		if (s_idx >= p.clusters_.size() || p.clusters_[s_idx].empty())
			compute_power_cluster_source(p);
		if (s_idx >= p.clusters_.size())
			return;
		const std::vector<PVertex>& cluster = p.clusters_[s_idx];
		if (cluster.empty())
			return;

		Vec3 c = (*p.spheres_position_)[s_idx];
		Scalar r = (*p.spheres_radius_)[s_idx];

		Eigen::MatrixXd J(2 * cluster.size(), 4);
		J.setZero();
		Eigen::VectorXd b(2 * cluster.size());
		b.setZero();
		uint32 idx = 0;
		Eigen::VectorXd s(4);
		s << c[0], c[1], c[2], r;

		for (uint32 i = 0; i < 10; ++i)
		{
			idx = 0;
			for (PVertex v : cluster)
			{
				uint32 v_index = index_of(*src.mesh, v);
				const Vec3& pos = (*src.position)[v_index];

				Eigen::Vector4d lhs = Eigen::Vector4d::Zero();
				Scalar rhs = 0.0;
				const Vec3& n = (*src.normal)[v_index];
				Vec4 n4 = Vec4(n.x(), n.y(), n.z(), 1.0);
				const int k = std::max<int>(1, static_cast<int>((*src.knn)[v_index].size()));
				Scalar a = sqrt((*src.area)[v_index] / (k + 1.0));
				lhs += -n4 * a;
				rhs += -1.0 * ((pos - Vec3(s(0), s(1), s(2))).dot(n) - s(3)) * a;

				for (PVertex vn : (*src.knn)[v_index])
				{
					uint32 vn_index = index_of(*src.mesh, vn);
					const Vec3& pn = (*src.position)[vn_index];
					const Vec3& nn = (*src.normal)[vn_index];
					Vec4 nn4 = Vec4(nn.x(), nn.y(), nn.z(), 1.0);
					const int kn = std::max<int>(1, static_cast<int>((*src.knn)[vn_index].size()));
					Scalar an = sqrt((*src.area)[vn_index] / (kn + 1.0));
					lhs += -nn4 * an;
					rhs += -1.0 * ((pn - Vec3(s(0), s(1), s(2))).dot(nn) - s(3)) * an;
				}
				J.row(idx) = lhs;
				b(idx) = rhs;
				++idx;

				Vec3 d = pos - Vec3(s(0), s(1), s(2));
				Scalar l = d.norm();
				if (l > Scalar(1e-12))
				{
					Scalar a_dist = sqrt((*src.area)[v_index]);
					J.row(idx) =
						Eigen::Vector4d(-(d[0] / l), -(d[1] / l), -(d[2] / l), -1.0) * a_dist * p.sqem_update_lambda_;
					b(idx) = -(l - s(3)) * a_dist * p.sqem_update_lambda_;
				}
				else
				{
					J.row(idx).setZero();
					b(idx) = 0.0;
				}
				++idx;
			};

			Eigen::LDLT<Eigen::MatrixXd> solver(J.transpose() * J);
			Eigen::VectorXd delta_s = solver.solve(J.transpose() * b);
			s += delta_s;
			if (delta_s.norm() < 1e-6)
				break;
		}
		c = s.head<3>();
		r = s[3];

		(*p.spheres_position_)[s_idx] = c;
		(*p.spheres_radius_)[s_idx] = r;
		if (points_provider_)
		{
			points_provider_->emit_attribute_changed(*p.spheres_mesh_, p.spheres_position_.get());
			points_provider_->emit_attribute_changed(*p.spheres_mesh_, p.spheres_radius_.get());
		}
	}

	void update_sphere_sqem_line_quadric_fix_alpha(Parameters& p, PVertex sphere)
	{
		SourceData src;
		if (!get_source_data(p, src))
			return;
		if (!p.spheres_position_ || !p.spheres_radius_)
			return;
		const uint32 s_idx = index_of(*p.spheres_mesh_, sphere);
		if (s_idx >= p.clusters_.size() || p.clusters_[s_idx].empty())
			compute_power_cluster_source(p);
		if (s_idx >= p.clusters_.size())
			return;
		const std::vector<PVertex>& cluster = p.clusters_[s_idx];
		if (cluster.empty())
			return;

		Spherical_Quadric q;
		Line_Quadric lq;
		for (PVertex v : cluster)
		{
			uint32 v_index = index_of(*src.mesh, v);
			Scalar weight = value<Scalar>(*src.mesh, src.area, v);
			if (weight <= 0.0)
				continue;
			q += (*src.quadric)[v_index] * weight;
			lq += (*src.line_quadric)[v_index] * weight;
		}

		Mat4 Ql = lq.get_quadric().matrix();
		Mat3 Al = Ql.block<3, 3>(0, 0);
		Vec3 bl = -Ql.block<3, 1>(0, 3);

		Mat3 As = q._A.block<3, 3>(0, 0);
		Vec3 bs = q._b.head<3>();
		Vec3 Asr = q._A.block<3, 1>(0, 3);

		Mat3 A = As + p.sqem_update_lambda_ * Al;
		Vec3 b = (bs + p.sqem_update_lambda_ * bl) - Asr * static_cast<Scalar>(p.alpha_);

		Vec3 c = A.ldlt().solve(b);
		(*p.spheres_position_)[s_idx] = c;
		(*p.spheres_radius_)[s_idx] = static_cast<Scalar>(p.alpha_);
		if (points_provider_)
		{
			points_provider_->emit_attribute_changed(*p.spheres_mesh_, p.spheres_position_.get());
			points_provider_->emit_attribute_changed(*p.spheres_mesh_, p.spheres_radius_.get());
		}
	}

	void update_sphere_sqem_line_quadric_free_radius(Parameters& p, PVertex sphere)
	{
		SourceData src;
		if (!get_source_data(p, src))
			return;
		if (!p.spheres_position_ || !p.spheres_radius_)
			return;
		const uint32 s_idx = index_of(*p.spheres_mesh_, sphere);
		if (s_idx >= p.clusters_.size() || p.clusters_[s_idx].empty())
			compute_power_cluster_source(p);
		if (s_idx >= p.clusters_.size())
			return;
		const std::vector<PVertex>& cluster = p.clusters_[s_idx];
		if (cluster.empty())
			return;

		Spherical_Quadric q;
		Line_Quadric lq;
		Scalar weight_sum = Scalar(0);
		for (PVertex v : cluster)
		{
			uint32 v_index = index_of(*src.mesh, v);
			Scalar weight = value<Scalar>(*src.mesh, src.area, v);
			if (weight <= 0.0)
				continue;
			q += (*src.quadric)[v_index] * weight;
			lq += (*src.line_quadric)[v_index] * weight;
			weight_sum += weight;
		}
		if (weight_sum <= Scalar(0))
			return;

		Mat4 A = q._A;
		Vec4 b = q._b;

		Mat4 Ql = lq.get_quadric().matrix();
		A.block<3, 3>(0, 0) += p.sqem_update_lambda_ * Ql.block<3, 3>(0, 0);
		b.head<3>() += p.sqem_update_lambda_ * (-Ql.block<3, 1>(0, 3));

		const Vec4 s = A.completeOrthogonalDecomposition().solve(b);
		if (!s.allFinite())
			return;

		(*p.spheres_position_)[s_idx] = s.head<3>();
		(*p.spheres_radius_)[s_idx] = s[3];
		if (points_provider_)
		{
			points_provider_->emit_attribute_changed(*p.spheres_mesh_, p.spheres_position_.get());
			points_provider_->emit_attribute_changed(*p.spheres_mesh_, p.spheres_radius_.get());
		}
	}

private:
	MeshProvider<POINTS>* points_provider_ = nullptr;
	PointCloudRender<POINTS>* pcr_ = nullptr;
	View* selected_view_ = nullptr;
	POINTS* selected_points_ = nullptr;
	std::map<POINTS*, Parameters> parameters_;
	std::array<std::mutex, 43> spheres_mutex_;
};

} // namespace ui

} // namespace cgogn

#endif // CGOGN_MODULE_ALPHA_SET_SPHERE_DEBUGGER_H_
