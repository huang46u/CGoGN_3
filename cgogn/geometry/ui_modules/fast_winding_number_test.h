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

#ifndef CGOGN_MODULE_FAST_WINDING_NUMBER_TEST_H_
#define CGOGN_MODULE_FAST_WINDING_NUMBER_TEST_H_

#include <cgogn/ui/app.h>
#include <cgogn/ui/module.h>

#include <cgogn/core/ui_modules/mesh_provider.h>
#include <cgogn/core/utils/numerics.h>

#include <cgogn/core/functions/attributes.h>
//#include <cgogn/core/functions/traversals/cell.h>
#include <cgogn/core/functions/traversals/vertex.h>

#include <cgogn/geometry/algos/area.h>
#include <cgogn/geometry/algos/centroid.h>
#include <cgogn/geometry/algos/normal.h>
#include <cgogn/geometry/types/fast_winding_number.h>
#include <cgogn/geometry/types/fast_winding_number_traits.h>
#include <cgogn/geometry/types/vector_traits.h>

#include <libacc/bvh_tree.h>
#include <libacc/bvh_tree_spheres.h>
#include <libacc/kd_tree.h>

#include <cmath>
#include <memory>
#include <random>
#include <unordered_map>

namespace cgogn
{

namespace ui
{

using cgogn::numerics::uint32;
using geometry::Scalar;
using geometry::Vec3;
using geometry::Vec4;

template <typename SURFACE, typename POINTS>
class FastWindingNumberTest : public ViewModule
{
	template <typename T>
	using SAttribute = typename mesh_traits<SURFACE>::template Attribute<T>;
	using SVertex = typename mesh_traits<SURFACE>::Vertex;
	using SFace = typename mesh_traits<SURFACE>::Face;

	template <typename T>
	using PAttribute = typename mesh_traits<POINTS>::template Attribute<T>;
	using PVertex = typename mesh_traits<POINTS>::Vertex;

	using TriangleTraits = geometry::FWN_Triangle_Traits<SURFACE>;
	using PointTraits = geometry::FWN_Point_Traits<SURFACE>;

	template <int ORDER>
	using TriangleFWN = geometry::Fast_Winding_Number<TriangleTraits, ORDER>;
	template <int ORDER>
	using PointFWN = geometry::Fast_Winding_Number<PointTraits, ORDER>;

	enum class EvalMode
	{
		Triangles,
		Points
	};

	struct SurfaceParameters
	{
		bool initialized_ = false;

		SURFACE* surface_ = nullptr;
		std::shared_ptr<SAttribute<Vec3>> surface_vertex_position_ = nullptr;
		std::shared_ptr<SAttribute<Vec3>> surface_vertex_normal_ = nullptr;
		std::shared_ptr<SAttribute<Vec3>> surface_face_normal_ = nullptr;
		std::shared_ptr<SAttribute<Scalar>> surface_face_area_ = nullptr;
		std::shared_ptr<SAttribute<Vec3>> surface_face_centroid_ = nullptr;
		std::shared_ptr<SAttribute<Scalar>> surface_vertex_area_pc_ = nullptr;
		std::shared_ptr<SAttribute<Scalar>> surface_vertex_area_surf_ = nullptr;

		std::unique_ptr<acc::BVHTree<uint32, Vec3>> surface_bvh_;
		std::vector<SFace> surface_bvh_faces_;
		std::vector<SVertex> surface_bvh_vertices_;
		std::vector<Vec3> surface_bvh_vertex_positions_;

		std::unique_ptr<acc::BVHTreeSpheres<uint32, Vec3>> point_bvh_;
		std::vector<Vec3> point_bvh_centers_;
		std::vector<double> point_bvh_radii_;

		std::unique_ptr<TriangleFWN<1>> tri_fwn_1_;
		std::unique_ptr<TriangleFWN<2>> tri_fwn_2_;
		std::unique_ptr<TriangleFWN<3>> tri_fwn_3_;

		std::unique_ptr<PointFWN<1>> point_fwn_1_;
		std::unique_ptr<PointFWN<2>> point_fwn_2_;
		std::unique_ptr<PointFWN<3>> point_fwn_3_;

		POINTS* samples_ = nullptr;
		std::shared_ptr<PAttribute<Vec3>> samples_position_ = nullptr;
		std::shared_ptr<PAttribute<Vec4>> samples_color_ = nullptr;
		std::shared_ptr<PAttribute<Scalar>> samples_radius_ = nullptr;

		uint32 sample_count_ = 200;
		Scalar sample_radius_ = Scalar(0.01);
		Scalar beta_ = Scalar(2.0);
	};

public:
	FastWindingNumberTest(const App& app)
		: ViewModule(app, "FastWindingNumberTest (" + std::string{mesh_traits<SURFACE>::name} + ")")
	{
	}
	~FastWindingNumberTest()
	{
	}

	void set_selected_surface(SURFACE& s)
	{
		selected_surface_ = &s;
	}

	void set_surface_vertex_position(SURFACE& s, const std::shared_ptr<SAttribute<Vec3>>& surface_vertex_position)
	{
		SurfaceParameters& p = surface_parameters_[&s];
		p.surface_vertex_position_ = surface_vertex_position;
	}

	void set_ghost_mode(View& v, const SURFACE& m, bool b)
	{
		Parameters& p = surface_parameters_[&v][&m];
		p.ghost_mode_ = b;
		p.param_flat_->ghost_mode_ = p.ghost_mode_;
		p.param_phong_->ghost_mode_ = p.ghost_mode_;
		v.request_update();
	}

	void init_surface_data(SURFACE& s)
	{
		SurfaceParameters& p = surface_parameters_[&s];
		p.surface_ = &s;

		if (!p.surface_vertex_position_)
		{
			std::cout << "No surface vertex position attribute set" << std::endl;
			return;
		}

		compute_surface_attributes(p);
		build_bvhs(p);
		build_fwn_instances(p);
		init_sample_mesh(p);

		p.initialized_ = true;
	}

protected:
	void init() override
	{
		surface_provider_ = static_cast<ui::MeshProvider<SURFACE>*>(
			app_.module("MeshProvider (" + std::string{mesh_traits<SURFACE>::name} + ")"));

		points_provider_ = static_cast<ui::MeshProvider<POINTS>*>(
			app_.module("MeshProvider (" + std::string{mesh_traits<POINTS>::name} + ")"));
	}

	void left_panel() override
	{
		imgui_mesh_selector(surface_provider_, selected_surface_, "Surface",
							[&](SURFACE& s) { set_selected_surface(s); });

		if (!selected_surface_)
			return;

		SurfaceParameters& p = surface_parameters_[selected_surface_];

		imgui_combo_attribute<SVertex, Vec3>(
			*selected_surface_, p.surface_vertex_position_, "Position",
			[&](const std::shared_ptr<SAttribute<Vec3>>& attribute) { p.surface_vertex_position_ = attribute; });

		if (p.surface_vertex_position_ && !p.initialized_)
		{
			if (ImGui::Button("Init surface data"))
				init_surface_data(*selected_surface_);
		}

		if (!p.initialized_)
			return;

		ImGui::InputScalar("Sample count", ImGuiDataType_U32, &p.sample_count_);

		if (ImGui::Button("Sample + Eval (Triangles)"))
			sample_and_evaluate(p, EvalMode::Triangles);
		if (ImGui::Button("Sample + Eval (Points)"))
			sample_and_evaluate(p, EvalMode::Points);
	}

private:
	void compute_surface_attributes(SurfaceParameters& p)
	{
		p.surface_face_normal_ = get_or_add_attribute<Vec3, SFace>(*p.surface_, "normal");
		geometry::compute_normal<SFace>(*p.surface_, p.surface_vertex_position_.get(), p.surface_face_normal_.get());

		p.surface_vertex_normal_ = get_or_add_attribute<Vec3, SVertex>(*p.surface_, "normal");
		geometry::compute_normal<SVertex>(*p.surface_, p.surface_vertex_position_.get(), p.surface_vertex_normal_.get());

		p.surface_face_area_ = get_or_add_attribute<Scalar, SFace>(*p.surface_, "area");
		geometry::compute_area<SFace>(*p.surface_, p.surface_vertex_position_.get(), p.surface_face_area_.get());

		p.surface_face_centroid_ = get_or_add_attribute<Vec3, SFace>(*p.surface_, "centroid");
		geometry::compute_centroid<Vec3,SFace>(*p.surface_, p.surface_vertex_position_.get(), p.surface_face_centroid_.get());

		p.surface_vertex_area_surf_ = get_or_add_attribute<Scalar, SVertex>(*p.surface_, "area_surf");
		geometry::compute_area<SVertex>(*p.surface_, p.surface_vertex_position_.get(), p.surface_vertex_area_surf_.get(),
									geometry::VertexAreaPolicy::THIRD);
	}

	void build_bvhs(SurfaceParameters& p)
	{
		MeshData<SURFACE>& md = surface_provider_->mesh_data(*p.surface_);
		uint32 nb_vertices = md.template nb_cells<SVertex>();
		uint32 nb_faces = md.template nb_cells<SFace>();

		auto bvh_vertex_index = get_or_add_attribute<uint32, SVertex>(*p.surface_, "__bvh_vertex_index");

		p.surface_bvh_vertices_.clear();
		p.surface_bvh_vertices_.reserve(nb_vertices);
		p.surface_bvh_vertex_positions_.clear();
		p.surface_bvh_vertex_positions_.reserve(nb_vertices);

		uint32 idx = 0;
		foreach_cell(*p.surface_, [&](SVertex v) -> bool {
			p.surface_bvh_vertices_.push_back(v);
			value<uint32>(*p.surface_, bvh_vertex_index, v) = idx++;
			p.surface_bvh_vertex_positions_.push_back(value<Vec3>(*p.surface_, p.surface_vertex_position_, v));
			return true;
		});

		p.surface_bvh_faces_.clear();
		p.surface_bvh_faces_.reserve(nb_faces);
		std::vector<uint32> face_vertex_indices;
		face_vertex_indices.reserve(nb_faces * 3);
		foreach_cell(*p.surface_, [&](SFace f) -> bool {
			p.surface_bvh_faces_.push_back(f);
			foreach_incident_vertex(*p.surface_, f, [&](SVertex v) -> bool {
				face_vertex_indices.push_back(value<uint32>(*p.surface_, bvh_vertex_index, v));
				return true;
			});
			return true;
		});

		p.surface_bvh_ = std::make_unique<acc::BVHTree<uint32, Vec3>>(face_vertex_indices, p.surface_bvh_vertex_positions_);

		compute_point_cloud_areas(p);

		p.point_bvh_centers_ = p.surface_bvh_vertex_positions_;
		p.point_bvh_radii_.clear();
		p.point_bvh_radii_.reserve(p.surface_bvh_vertices_.size());
		for (SVertex v : p.surface_bvh_vertices_)
		{
			Scalar area = value<Scalar>(*p.surface_, p.surface_vertex_area_pc_, v);
			Scalar radius = std::sqrt(area / M_PI);
			p.point_bvh_radii_.push_back(static_cast<double>(radius));
		}

		p.point_bvh_ = std::make_unique<acc::BVHTreeSpheres<uint32, Vec3>>(p.point_bvh_centers_, p.point_bvh_radii_);

		remove_attribute<SVertex>(*p.surface_, bvh_vertex_index);
	}

	void compute_point_cloud_areas(SurfaceParameters& p)
	{
		static constexpr uint32 k = 6;
		p.surface_vertex_area_pc_ = get_or_add_attribute<Scalar, SVertex>(*p.surface_, "area_pc");

		acc::KDTree<3, uint32> kdt(p.surface_bvh_vertex_positions_);
		for (uint32 i = 0; i < p.surface_bvh_vertices_.size(); ++i)
		{
			const Vec3& pos = p.surface_bvh_vertex_positions_[i];
			std::vector<std::pair<uint32, double>> k_res;
			kdt.find_nns(pos, k, &k_res);
			Scalar sum = 0;
			for (const auto& [idx, dist] : k_res)
			{
				if (idx == i)
					continue;
				sum += Scalar(dist);
			}
			value<Scalar>(*p.surface_, p.surface_vertex_area_pc_, p.surface_bvh_vertices_[i]) =
				(sum * sum) / (Scalar(2.0) * Scalar(k));
		}
	}

	void build_fwn_instances(SurfaceParameters& p)
	{
		TriangleTraits tri_traits(*p.surface_, p.surface_bvh_.get(), p.surface_bvh_faces_, p.surface_vertex_position_.get(),
							 p.surface_face_normal_.get(), p.surface_face_area_.get(), p.surface_face_centroid_.get());

		p.tri_fwn_1_ = std::make_unique<TriangleFWN<1>>(*p.surface_bvh_, tri_traits, p.beta_);
		p.tri_fwn_2_ = std::make_unique<TriangleFWN<2>>(*p.surface_bvh_, tri_traits, p.beta_);
		p.tri_fwn_3_ = std::make_unique<TriangleFWN<3>>(*p.surface_bvh_, tri_traits, p.beta_);

		PointTraits point_traits(*p.surface_, p.point_bvh_.get(), p.surface_bvh_vertices_, p.surface_vertex_position_.get(),
							 p.surface_vertex_normal_.get(), p.surface_vertex_area_pc_.get());

		p.point_fwn_1_ = std::make_unique<PointFWN<1>>(*p.point_bvh_, point_traits, p.beta_);
		p.point_fwn_2_ = std::make_unique<PointFWN<2>>(*p.point_bvh_, point_traits, p.beta_);
		p.point_fwn_3_ = std::make_unique<PointFWN<3>>(*p.point_bvh_, point_traits, p.beta_);
	}

	void init_sample_mesh(SurfaceParameters& p)
	{
		if (!p.samples_)
			p.samples_ = points_provider_->add_mesh(surface_provider_->mesh_name(*p.surface_) + "_fwn_samples");

		p.samples_position_ = get_or_add_attribute<Vec3, PVertex>(*p.samples_, "position");
		p.samples_color_ = get_or_add_attribute<Vec4, PVertex>(*p.samples_, "color");
		p.samples_radius_ = get_or_add_attribute<Scalar, PVertex>(*p.samples_, "radius");

		points_provider_->set_mesh_bb_vertex_position(*p.samples_, p.samples_position_);
	}

	void sample_and_evaluate(SurfaceParameters& p, EvalMode mode)
	{
		if (!p.initialized_)
			return;

		clear(*p.samples_);

		std::mt19937 rng(std::random_device{}());
		std::uniform_real_distribution<Scalar> dist(Scalar(0.0), Scalar(1.0));

		for (uint32 i = 0; i < p.sample_count_; ++i)
		{
			PVertex v = add_vertex(*p.samples_);
			uint32 v_index = index_of(*p.samples_, v);
			Vec3 pos(dist(rng), dist(rng), dist(rng));
			(*p.samples_position_)[v_index] = pos;
			(*p.samples_radius_)[v_index] = p.sample_radius_;

			Scalar wn1 = 0;
			Scalar wn2 = 0;
			Scalar wn3 = 0;
			Scalar exact = 0;

			if (mode == EvalMode::Triangles)
			{
				wn1 = p.tri_fwn_1_->evaluate_fast_winding_number(pos);
				wn2 = p.tri_fwn_2_->evaluate_fast_winding_number(pos);
				wn3 = p.tri_fwn_3_->evaluate_fast_winding_number(pos);
				exact = p.tri_fwn_3_->exact_winding_number(pos);
			}
			else
			{
				wn1 = p.point_fwn_1_->evaluate_fast_winding_number(pos);
				wn2 = p.point_fwn_2_->evaluate_fast_winding_number(pos);
				wn3 = p.point_fwn_3_->evaluate_fast_winding_number(pos);
				exact = p.point_fwn_3_->exact_winding_number(pos);
			}

			const Scalar diff = wn3 - exact;

			std::cout << "Sample " << i << " pos (" << pos[0] << ", " << pos[1] << ", " << pos[2] << ")" << std::endl;
			std::cout << " order1=" << wn1 << " order2=" << wn2 << " order3=" << wn3
					  << " exact=" << exact << std::endl;
			std::cout << "  order3-exact=" << diff << std::endl;

			if (wn3 > Scalar(0.5))
				(*p.samples_color_)[v_index] = Vec4(Scalar(0.1), Scalar(0.8), Scalar(0.1), Scalar(1.0));
			else
				(*p.samples_color_)[v_index] = Vec4(Scalar(0.9), Scalar(0.2), Scalar(0.2), Scalar(1.0));
		}

		points_provider_->emit_connectivity_changed(*p.samples_);
		points_provider_->emit_attribute_changed(*p.samples_, p.samples_position_.get());
		points_provider_->emit_attribute_changed(*p.samples_, p.samples_radius_.get());
		points_provider_->emit_attribute_changed(*p.samples_, p.samples_color_.get());
	}

private:
	MeshProvider<SURFACE>* surface_provider_ = nullptr;
	MeshProvider<POINTS>* points_provider_ = nullptr;

	std::unordered_map<const SURFACE*, SurfaceParameters> surface_parameters_;
	SURFACE* selected_surface_ = nullptr;
};

} // namespace ui

} // namespace cgogn

#endif // CGOGN_MODULE_FAST_WINDING_NUMBER_TEST_H_
