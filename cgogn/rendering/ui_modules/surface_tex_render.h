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

#ifndef CGOGN_MODULE_SURFACE_TEX_RENDER_H_
#define CGOGN_MODULE_SURFACE_TEX_RENDER_H_

#include <cgogn/core/ui_modules/mesh_provider.h>
#include <cgogn/ui/app.h>
#include <cgogn/ui/imgui_helpers.h>
#include <cgogn/ui/module.h>
#include <cgogn/ui/view.h>
#include <cgogn/geometry/types/vector_traits.h>


#include <cgogn/rendering/shaders/shader_flat_texture.h>
#include <cgogn/rendering/shaders/shader_flat_color_per_vertex.h>

#include <cgogn/rendering/texture.h>

#include <boost/synapse/connect.hpp>
#include <unordered_map>

namespace cgogn
{

namespace ui
{

using geometry::Scalar;
using geometry::Vec3;
using geometry::Vec2;

template <typename MESH>
class SurfaceTexRender : public ViewModule
{
	static_assert(mesh_traits<MESH>::dimension >= 2, "SurfaceRender can only be used with meshes of dimension >= 2");


	template <typename T>
	using Attribute = typename mesh_traits<MESH>::template Attribute<T>;

	using Vertex = typename mesh_traits<MESH>::Vertex;
	using Edge = typename mesh_traits<MESH>::Edge;
	using Face = typename mesh_traits<MESH>::Face;

	struct Parameters
	{
		Parameters()
			: vertex_position_(nullptr), vertex_position_vbo_(nullptr), vertex_tc_(nullptr), vertex_tc_vbo_(nullptr), 
            render_texture_(true)
		{
			param_textured_ = rendering::ShaderFlatTexture::generate_param();

            //param_textured_->draw_param_ = true;
		}

		CGOGN_NOT_COPYABLE_NOR_MOVABLE(Parameters);

		std::unique_ptr<rendering::ShaderFlatTexture::Param> param_textured_;
        
		std::shared_ptr<Attribute<Vec3>> vertex_position_;
		rendering::VBO* vertex_position_vbo_;

		std::shared_ptr<Attribute<Vec2>> vertex_tc_;
		rendering::VBO* vertex_tc_vbo_;
        bool render_texture_;

	};

public:
	SurfaceTexRender(const App& app)
		: ViewModule(app, "SurfaceTextureRender (" + std::string{mesh_traits<MESH>::name} + ")"),
		  selected_view_(app.current_view()), selected_mesh_(nullptr)
	{
	}

	~SurfaceTexRender()
	{
	}

private:


	void init_mesh(MESH* m)
	{
        for (View* v : linked_views_)
		{
			parameters_[v][m];
			std::shared_ptr<Attribute<Vec3>> vertex_position = cgogn::get_attribute<Vec3, Vertex>(*m, "position");
			if (vertex_position)
				set_vertex_position(*v, *m, vertex_position);
            std::shared_ptr<Attribute<Vec2>> vertex_tc = cgogn::get_attribute<Vec2, Vertex>(*m, "tex_coord");
            if (vertex_tc)
				set_vertex_texcoord(*v, *m, vertex_tc);
		}
	}



public:
	void load_texture(const std::string& img_name)
	{
		rendering::GLImage img(img_name);
		tex_->load(img);
	}

	void load_texture(const rendering::GLImage& img)
	{
		tex_->load(img);
	}

	void chekered_texture()
	{
		cgogn::rendering::GLImage img(16, 16, 3);
		std::vector<std::array<cgogn::uint8, 3>> pix;
		pix.reserve(16*16);
		for (int i = 0; i < 16; ++i)
			for (int j = 0; j < 16; ++j)
				if ((i + j) % 2 == 0)
					pix.push_back({0u, 0u, 0u});
				else
					pix.push_back({255u, 255u, 255u});
		img.copy_pixels_data(pix.data()->data());
		load_texture(img);
	}

	void set_vertex_position(View& v, const MESH& m, const std::shared_ptr<Attribute<Vec3>>& vertex_position)
	{
		Parameters& p = parameters_[&v][&m];
		if (p.vertex_position_ == vertex_position) 
			return;

		MeshData<MESH>& md = mesh_provider_->mesh_data(m);

		p.vertex_position_ = vertex_position;
		if (p.vertex_position_)
			p.vertex_position_vbo_ = md.update_vbo(p.vertex_position_.get(), true);
		else
			p.vertex_position_vbo_ = nullptr;

		p.param_textured_->set_vbos({p.vertex_position_vbo_, p.vertex_tc_vbo_});
		v.request_update();	
	}


	void set_vertex_texcoord(View& v, const MESH& m, const std::shared_ptr<Attribute<Vec2>>& vertex_tc)
	{
		Parameters& p = parameters_[&v][&m];
		if (p.vertex_tc_ == vertex_tc)
			return;

		MeshData<MESH>& md = mesh_provider_->mesh_data(m);

		p.vertex_tc_ = vertex_tc;
		if (p.vertex_tc_)
			p.vertex_tc_vbo_ = md.update_vbo(p.vertex_tc_.get(), true);
		else
			p.vertex_tc_vbo_ = nullptr;

		p.param_textured_->set_vbos({p.vertex_position_vbo_, p.vertex_tc_vbo_});
		v.request_update();	
	}

	void init() override
	{
		mesh_provider_ = static_cast<ui::MeshProvider<MESH>*>(
			app_.module("MeshProvider (" + std::string{mesh_traits<MESH>::name} + ")"));
		mesh_provider_->foreach_mesh([this](MESH& m, const std::string&) { init_mesh(&m); });
		connections_.push_back(boost::synapse::connect<typename MeshProvider<MESH>::mesh_added>(
			mesh_provider_, this, &SurfaceTexRender<MESH>::init_mesh));
		tex_ = std::make_shared<rendering::Texture2D>(
			std::vector<std::pair<GLenum, GLint>>{
			{GL_TEXTURE_MIN_FILTER, GL_LINEAR}, {GL_TEXTURE_MAG_FILTER, GL_NEAREST}, {GL_TEXTURE_WRAP_S, GL_REPEAT},
													   {GL_TEXTURE_WRAP_T, GL_REPEAT}});
	}


	void draw(View* view) override
	{
		for (auto& [m, pp] : parameters_[view])
		{
			Parameters& p = pp;
            if (p.render_texture_)
            {
                MeshData<MESH>& md = mesh_provider_->mesh_data(*m);

                const rendering::GLMat4& proj_matrix = view->projection_matrix();
                const rendering::GLMat4& view_matrix = view->modelview_matrix();
                    if (p.param_textured_->attributes_initialized())
                {
                        p.param_textured_->texture_ = tex_;
                        p.param_textured_->bind(proj_matrix, view_matrix);
                        md.draw(rendering::TRIANGLES, p.vertex_position_);
                        p.param_textured_->release();

                }
            }
		}
	}



	void left_panel() override
	{
	    bool need_update = false;

		if (app_.nb_views() > 1)
			imgui_view_selector(this, selected_view_, [&](View* v) { selected_view_ = v; });

		imgui_mesh_selector(mesh_provider_, selected_mesh_, "Surface", [&](MESH& m) {
			selected_mesh_ = &m;
//			mesh_provider_->mesh_data(m).outlined_until_ = App::frame_time_ + 1.0;
		});

		if (selected_view_ && selected_mesh_)
		{
			Parameters& p = parameters_[selected_view_][selected_mesh_];

			imgui_combo_attribute<Vertex, Vec3>(*selected_mesh_, p.vertex_position_, "Position",
												[&](const std::shared_ptr<Attribute<Vec3>>& attribute) {
													set_vertex_position(*selected_view_, *selected_mesh_, attribute);
												});
			ImGui::Separator();
			imgui_combo_attribute<Vertex, Vec2>(*selected_mesh_, p.vertex_tc_, "TexCoord",
												[&](const std::shared_ptr<Attribute<Vec2>>& attribute) {
													set_vertex_texcoord(*selected_view_, *selected_mesh_, attribute);
												});

            need_update |= ImGui::Checkbox("Render texture", &p.render_texture_);
        }

        if (need_update)
            for (View* v : linked_views_)
                v->request_update();
	}

private:
	View* selected_view_;
	const MESH* selected_mesh_;
	std::unordered_map<View*, std::unordered_map<const MESH*, Parameters>> parameters_;
	std::vector<std::shared_ptr<boost::synapse::connection>> connections_;
	std::unordered_map<const MESH*, std::vector<std::shared_ptr<boost::synapse::connection>>> mesh_connections_;
	MeshProvider<MESH>* mesh_provider_;
	std::shared_ptr<rendering::Texture2D> tex_;
};

} // namespace ui

} // namespace cgogn

#endif // CGOGN_MODULE_SURFACE_RENDER_H_
