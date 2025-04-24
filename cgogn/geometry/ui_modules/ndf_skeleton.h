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

#ifndef CGOGN_MODULE_NDF_SKELETON_H_
#define CGOGN_MODULE_NDF_SKELETON_H_

#include <cgogn/ui/app.h>
#include <cgogn/ui/module.h>

#include <cgogn/core/ui_modules/mesh_provider.h>
#include <cgogn/geometry/types/vector_traits.h>
#include <filesystem>
#include <torch/script.h>
#include <torch/torch.h>

namespace cgogn
{

namespace ui
{

using geometry::Scalar;
using geometry::Vec3;
using geometry::Vec4;

Vec4 color_map(Scalar x, Scalar min, Scalar max, float32 transparency = 1.0)
{
	x = (x - min) / (max - min);
	x = std::clamp(x, 0.0, 1.0);

	Scalar x2 = 2.0 * x;
	switch (int(std::floor(std::max(0.0, x2 + 1.0))))
	{
	case 0:
		return Vec4(0.0, 0.0, 1.0, transparency);
	case 1:
		return Vec4(x2, x2, 1.0, transparency);
	case 2:
		return Vec4(1.0, 2.0 - x2, 2.0 - x2, transparency);
	}
	return Vec4(1.0, 0.0, 0.0, transparency);
}

template <typename SURFACE, typename POINTS>
class LipNDF : public Module
{
	template <typename T>
	using PointsAttribute = typename mesh_traits<POINTS>::template Attribute<T>;
	template <typename T>
	using SurfaceAttribute = typename mesh_traits<SURFACE>::template Attribute<T>;

	using PointsVertex = typename mesh_traits<POINTS>::Vertex;
	using SurfaceVertex = typename mesh_traits<SURFACE>::Vertex;

	struct Parameters
	{
		bool initialized_ = false;

		SURFACE* surface_ = nullptr;
		std::shared_ptr<SurfaceAttribute<Vec3>> surface_vertex_position_ = nullptr;
		POINTS* spheres_ = nullptr;
		std::shared_ptr<PointsAttribute<Vec3>> spheres_position_ = nullptr;
		std::shared_ptr<PointsAttribute<Scalar>> spheres_radius_ = nullptr;
		std::shared_ptr<PointsAttribute<Vec4>> spheres_color_ = nullptr;
		torch::jit::Module model_;
	};

public:
	LipNDF(const App& app)
		: Module(app, "LipNDF (" + std::string{mesh_traits<SURFACE>::name} + "," +
						  std::string{mesh_traits<POINTS>::name} + ")")
	{
	}
	~LipNDF()
	{
	}

	void set_selected_surface(SURFACE& s)
	{
		selected_surface_ = &s;
	}

	void set_surface_vertex_position(const SURFACE& s, const std::shared_ptr<SurfaceAttribute<Vec3>>& vertex_position)
	{
		Parameters& p = parameters_[&s];
		p.surface_vertex_position_ = vertex_position;
	}

	void init_surface_data(SURFACE& s, const std::string& model_path)
	{
		Parameters& p = parameters_[&s];
		p.surface_ = &s;

		if (!std::filesystem::exists(model_path))
		{
			std::cout << "Model file does not exist: " << model_path << std::endl;
			return;
		}

		if (!p.surface_vertex_position_)
		{
			std::cout << "No surface vertex position attribute set" << std::endl;
			return;
		}

		p.model_ = torch::jit::load(model_path, device_);

		p.spheres_ = points_provider_->add_mesh(surface_provider_->mesh_name(s) + "_spheres");
		p.spheres_position_ = get_or_add_attribute<Vec3, PointsVertex>(*p.spheres_, "position");
		p.spheres_radius_ = get_or_add_attribute<Scalar, PointsVertex>(*p.spheres_, "radius");
		p.spheres_color_ = get_or_add_attribute<Vec4, PointsVertex>(*p.spheres_, "color");

		p.initialized_ = true;
	}

	void add_random_sample(SURFACE& s)
	{
		Parameters& p = parameters_[&s];

		Vec3 rp = Vec3::Random(); // random point in [-1, 1]^3
		rp /= Scalar(2);		  // contract to [-0.5, 0.5]^3

		at::Tensor point = torch::tensor({rp[0], rp[1], rp[2]}, torch::kFloat32);
		at::Tensor output = p.model_.forward({point}).toTensor();
		Scalar radius = output.item<Scalar>();

		PointsVertex v = cgogn::add_vertex(*p.spheres_);
		cgogn::value<Vec3>(*p.spheres_, p.spheres_position_, v) = rp;
		cgogn::value<Scalar>(*p.spheres_, p.spheres_radius_, v) = radius;

		points_provider_->emit_connectivity_changed(*p.spheres_);
		points_provider_->emit_attribute_changed(*p.spheres_, p.spheres_position_.get());
		points_provider_->emit_attribute_changed(*p.spheres_, p.spheres_radius_.get());

		update_spheres_color(s);
	}

	void sample_sdf_on_grid(SURFACE& s, uint32 N)
	{
		Parameters& p = parameters_[&s];

		clear(*p.spheres_);

		Scalar step = 1.0f / (N - 1);
		for (size_t i = 0; i < N; i++)
		{
			for (size_t j = 0; j < N; j++)
			{
				for (size_t k = 0; k < N; k++)
				{
					Scalar x = -0.5f + i * step;
					Scalar y = -0.5f + j * step;
					Scalar z = -0.5f + k * step;

					at::Tensor point = torch::tensor({x, y, z}, torch::kFloat32);
					at::Tensor output = p.model_.forward({point}).toTensor();
					Scalar radius = output.item<Scalar>();

					PointsVertex v = cgogn::add_vertex(*p.spheres_);
					cgogn::value<Vec3>(*p.spheres_, p.spheres_position_, v) = {x, y, z};
					cgogn::value<Scalar>(*p.spheres_, p.spheres_radius_, v) = radius;
				}
			}
		}

		points_provider_->emit_connectivity_changed(*p.spheres_);
		points_provider_->emit_attribute_changed(*p.spheres_, p.spheres_position_.get());
		points_provider_->emit_attribute_changed(*p.spheres_, p.spheres_radius_.get());

		update_spheres_color(s);
	}

	void clear_samples(SURFACE& s)
	{
		Parameters& p = parameters_[&s];
		clear(*p.spheres_);
		points_provider_->emit_connectivity_changed(*p.spheres_);
	}

	void update_spheres_color(SURFACE& s)
	{
		Parameters& p = parameters_[&s];

		Scalar min_radius = std::numeric_limits<Scalar>::max();
		Scalar max_radius = std::numeric_limits<Scalar>::min();
		for (Scalar r : *p.spheres_radius_)
		{
			min_radius = std::min(min_radius, r);
			max_radius = std::max(max_radius, r);
		}
		cgogn::parallel_foreach_cell(*p.spheres_, [&](PointsVertex v) -> bool {
			cgogn::value<Vec4>(*p.spheres_, p.spheres_color_, v) =
				color_map(cgogn::value<Scalar>(*p.spheres_, p.spheres_radius_, v), min_radius, max_radius, 0.2);
			return true;
		});

		points_provider_->emit_attribute_changed(*p.spheres_, p.spheres_color_.get());
	}

protected:
	void init() override
	{
		points_provider_ = static_cast<ui::MeshProvider<POINTS>*>(
			app_.module("MeshProvider (" + std::string{mesh_traits<POINTS>::name} + ")"));
		surface_provider_ = static_cast<ui::MeshProvider<SURFACE>*>(
			app_.module("MeshProvider (" + std::string{mesh_traits<SURFACE>::name} + ")"));

		if (torch::cuda::is_available())
		{
			std::cout << "CUDA is available! Using the GPU." << std::endl;
			device_ = torch::kCUDA;
		}
		std::cout << "CUDA is not available! Using the CPU." << std::endl;
		device_ = torch::kCPU;
	}

	void left_panel() override
	{
		imgui_mesh_selector(surface_provider_, selected_surface_, "Surface",
							[&](SURFACE& s) { set_selected_surface(s); });

		if (selected_surface_)
		{
			Parameters& p = parameters_[selected_surface_];

			imgui_combo_attribute<SurfaceVertex, Vec3>(*selected_surface_, p.surface_vertex_position_, "Position",
													   [&](const std::shared_ptr<SurfaceAttribute<Vec3>>& attribute) {
														   set_surface_vertex_position(*selected_surface_, attribute);
													   });

			if (p.initialized_)
			{
				if (ImGui::Button("Add random sample"))
					add_random_sample(*selected_surface_);

				static uint32 N = 5;
				ImGui::InputScalar("Grid size", ImGuiDataType_U32, &N);
				if (ImGui::Button("Sample grid"))
					sample_sdf_on_grid(*selected_surface_, N);

				if (ImGui::Button("Clear"))
					clear_samples(*selected_surface_);
			}
		}
	}

private:
	SURFACE* selected_surface_ = nullptr;
	std::unordered_map<const SURFACE*, Parameters> parameters_;
	std::vector<std::shared_ptr<boost::synapse::connection>> connections_;

	MeshProvider<POINTS>* points_provider_;
	MeshProvider<SURFACE>* surface_provider_;

	torch::Device device_ = torch::kCPU;
};

} // namespace ui

} // namespace cgogn

#endif // CGOGN_MODULE_NDF_SKELETON_H_
