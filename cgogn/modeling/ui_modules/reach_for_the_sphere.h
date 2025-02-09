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

#ifndef CGOGN_MODULE_REACH_FOR_THE_SPHERE_
#define CGOGN_MODULE_REACH_FOR_THE_SPHERE_

#include <random>

#include <cgogn/core/ui_modules/mesh_provider.h>
#include <cgogn/ui/app.h>
#include <cgogn/ui/module.h>

#include <cgogn/core/functions/convert.h>

#include <cgogn/core/functions/traversals/face.h>
#include <cgogn/core/functions/traversals/vertex.h>

#include <cgogn/geometry/types/vector_traits.h>

#include <cgogn/geometry/algos/length.h>
#include <cgogn/geometry/algos/medial_axis.h>

#include <cgogn/modeling/algos/remeshing/pliant_remeshing.h>
#include <cgogn/modeling/algos/skeleton.h>

#include <boost/synapse/connect.hpp>

#include <unordered_map>

namespace cgogn
{

namespace ui
{

using geometry::Mat3;
using geometry::Scalar;
using geometry::Vec3;
using geometry::Vec4;

template <typename POINT, typename SURFACE>
class ReachForTheSphere : public Module
{
	template <typename T>
	using PointAttribute = typename mesh_traits<POINT>::template Attribute<T>;
	template <typename T>
	using SurfaceAttribute = typename mesh_traits<SURFACE>::template Attribute<T>;

	using PointVertex = typename mesh_traits<POINT>::Vertex;

	using SurfaceVertex = typename mesh_traits<SURFACE>::Vertex;
	using SurfaceEdge = typename mesh_traits<SURFACE>::Edge;
	using SurfaceFace = typename mesh_traits<SURFACE>::Face;

	enum SamplingMode : uint32
	{
		GRILLE,
		RANDOM
	};

	struct SurfaceParameters
	{
		bool initialized_;

		SURFACE* surface_;
		std::shared_ptr<SurfaceAttribute<Vec3>> surface_vertex_position_ = nullptr;

		acc::BVHTree<uint32, Vec3>* surface_bvh_ = nullptr;
		std::vector<SurfaceFace> surface_bvh_faces_;

		SURFACE* flow_mesh_;

		std::shared_ptr<SurfaceAttribute<Vec3>> flow_vertex_position_ = nullptr;
		std::shared_ptr<SurfaceAttribute<Vec3>> flow_vertex_normal_ = nullptr;
		std::shared_ptr<SurfaceAttribute<uint32>> flow_vertex_id_ = nullptr;

		acc::BVHTree<uint32, Vec3>* flow_mesh_bvh_ = nullptr;
		std::vector<SurfaceFace> flow_mesh_bvh_faces_;


		POINT* sdf_;
		std::shared_ptr<PointAttribute<Vec3>> sdf_sample_position_ = nullptr;
		std::shared_ptr<PointAttribute<Scalar>> sdf_sample_value_ = nullptr;
		std::shared_ptr<PointAttribute<Scalar>> sdf_sample_radius_ = nullptr;
		std::shared_ptr<PointAttribute<Vec4>> sdf_sample_color_ = nullptr;

		Scalar t_min = 1e-6;
		Scalar t_max = 50;

		SamplingMode sampling_mode_ = GRILLE;

		uint32 grille_sample_resolution_ = 10;
		uint32 random_sample_number_ = 10000; 
	};

public:
	ReachForTheSphere(const App& app)
		: Module(app, "ReachForTheSphere"),
		  selected_surface_(nullptr), selected_surface_vertex_position_(nullptr)
	{
	}

	~ReachForTheSphere()
	{
	}

private:
	void init_surface_data(SURFACE& s)
	{
		SurfaceParameters& p = surface_parameters_[&s];
		p.surface_ = &s;

		if (!p.surface_vertex_position_)
		{
			std::cout << "No surface vertex position attribute set" << std::endl;
			return;
		}

		// set signal connections to update the data when the surface connectivity or position changes

		if (surface_connections_.find(&s) == surface_connections_.end())
		{
			surface_connections_[&s].push_back(
				boost::synapse::connect<typename MeshProvider<SURFACE>::connectivity_changed>(&s, [this, &s = s]() {
					SurfaceParameters& p = surface_parameters_[&s];
					init_surface_data(s);
				}));
			surface_connections_[&s].push_back(
				boost::synapse::connect<typename MeshProvider<SURFACE>::template attribute_changed_t<Vec3>>(
					&s, [this, &s = s](SurfaceAttribute<Vec3>* attribute) {
						SurfaceParameters& p = surface_parameters_[&s];
						if (attribute == p.surface_vertex_position_.get())
							init_surface_data(s);
					}));
		}
		// Build the BVH for the surface mesh
		
		build_bvh(s, p.surface_vertex_position_, p.surface_bvh_, p.surface_bvh_faces_);

		// Create sdf mesh
		if (!p.sdf_)
			p.sdf_ = points_provider_->add_mesh(surface_provider_->mesh_name(s) + "_sdf");
		
		clear(*p.sdf_);

		p.sdf_sample_position_ = get_or_add_attribute<Vec3, PointVertex>(*p.sdf_, "position");
		p.sdf_sample_value_ = get_or_add_attribute<Scalar, PointVertex>(*p.sdf_, "sdf");
		p.sdf_sample_color_ = get_or_add_attribute<Vec4, PointVertex>(*p.sdf_, "color");
		p.sdf_sample_radius_ = get_or_add_attribute<Scalar, PointVertex>(*p.sdf_, "radius");

		// Sample the surface to create the SDF mesh
		if (p.sampling_mode_ == GRILLE)
			grille_samping_mesh(s, p.grille_sample_resolution_);
		if (p.sampling_mode_ == RANDOM)
			random_sampling_mesh(s, p.random_sample_number_);

		// Compute the convex hull of the sampled points
		compute_convex_hull(s);

		// Ensure the flow mesh is created before building its BVH
		if (!p.flow_mesh_)
		{
			std::cout << "Flow mesh not created" << std::endl;
			return;
		}

		// Build the BVH for the flow mesh
		build_bvh(*p.flow_mesh_, p.flow_vertex_position_, p.flow_mesh_bvh_, p.flow_mesh_bvh_faces_);

		p.initialized_ = true;
		
	}

	bool is_inside(SURFACE& s, Vec3& pos, std::vector<SurfaceFace>& bvh_faces_,
				   std::shared_ptr<SurfaceAttribute<Vec3>>& vertex_position_, std::pair<uint32, Vec3>& cp)
	{

		Vec3 dir = (cp.second - pos).normalized();
		Vec3 n = geometry::normal(s, bvh_faces_[cp.first], vertex_position_.get());
		// If point is inside of the mesh
		return dir.dot(n) >= 0.0; // TODO: not reliable, better use general winding number.
	}

public:

	void build_bvh(SURFACE& s, std::shared_ptr<SurfaceAttribute<Vec3>>& surface_vertex_position,
					   acc::BVHTree<uint32, Vec3>*& surface_bvh, std::vector<SurfaceFace>& surface_bvh_faces)
		{
			MeshData<SURFACE>& md = surface_provider_->mesh_data(s);
			uint32 nb_vertices = md.template nb_cells<SurfaceVertex>();
			uint32 nb_faces = md.template nb_cells<SurfaceFace>();

			auto bvh_vertex_index = get_or_add_attribute<uint32, SurfaceVertex>(s, "__bvh_vertex_index");

			std::vector<Vec3> vertex_position_vector;
			vertex_position_vector.reserve(nb_vertices);
			uint32 idx = 0;
			foreach_cell(s, [&](SurfaceVertex v) -> bool {
				value<uint32>(s, bvh_vertex_index, v) = idx++;
				vertex_position_vector.push_back(value<Vec3>(s, surface_vertex_position, v));
				return true;
			});

			surface_bvh_faces.clear();
			surface_bvh_faces.reserve(nb_faces);
			std::vector<uint32> face_vertex_indices;
			face_vertex_indices.reserve(nb_faces * 3);
			foreach_cell(s, [&](SurfaceFace f) -> bool {
				surface_bvh_faces.push_back(f);
				foreach_incident_vertex(s, f, [&](SurfaceVertex v) -> bool {
					face_vertex_indices.push_back(value<uint32>(s, bvh_vertex_index, v));
					return true;
				});
				return true;
			});

			if (surface_bvh)
				delete surface_bvh;
			surface_bvh = new acc::BVHTree<uint32, Vec3>(face_vertex_indices, vertex_position_vector);

			remove_attribute<SurfaceVertex>(s, bvh_vertex_index);
		}

	void random_sampling_mesh(SURFACE& surface, uint32 numbers)
	{
		
		SurfaceParameters& p = surface_parameters_[&surface];

		// Compute bounding box of the surface mesh
		Vec3 bb_min, bb_max;
		std::tie(bb_min, bb_max) = geometry::bounding_box(*p.surface_vertex_position_);

		std::random_device rd;
		std::mt19937 gen(rd());
		std::uniform_real_distribution<> dis_x(bb_min[0], bb_max[0]);
		std::uniform_real_distribution<> dis_y(bb_min[1], bb_max[1]);
		std::uniform_real_distribution<> dis_z(bb_min[2], bb_max[2]);

		for (uint32 i = 0; i < numbers; ++i)
		{
			Vec3 sample_point(dis_x(gen), dis_y(gen), dis_z(gen));

			std::pair<uint32, Vec3> cp;
			p.surface_bvh_->closest_point(sample_point, &cp);
			float distance = (cp.second - sample_point).norm();
			float signed_distance = distance;
			if (is_inside(surface, sample_point, p.surface_bvh_faces_, p.surface_vertex_position_, cp))
			{
				signed_distance = -signed_distance;
			}

			// Add sample point to the point mesh
			PointVertex pv = add_vertex(*p.sdf_);
			value<Scalar>(*p.sdf_, p.sdf_sample_radius_, pv) = distance;
			value<Scalar>(*p.sdf_, p.sdf_sample_value_, pv) = signed_distance;
			value<Vec3>(*p.sdf_, p.sdf_sample_position_, pv) = sample_point;

			// Assign color based on signed distance
			Vec4 color;
			if (signed_distance < 0)
			{
				color = Vec4(1.0, 0.0, 0.0, 0.1); // Red for inside
			}
			else
			{
				color = Vec4(0.0, 0.0, 1.0, 0.1); // Blue for outside
			}
			value<Vec4>(*p.sdf_, p.sdf_sample_color_, pv) = color;
		}

		std::cout << "Finished random sampling" << std::endl;
		points_provider_->emit_connectivity_changed(*p.sdf_);
		points_provider_->emit_attribute_changed(*p.sdf_, p.sdf_sample_position_.get());

		points_provider_->set_mesh_bb_vertex_position(*p.sdf_, p.sdf_sample_position_);
		

	}

	//marching cube like sampling
	void grille_samping_mesh(SURFACE& surface, uint32 step = 50)
	{
		SurfaceParameters& p = surface_parameters_[&surface];
		// Compute bounding box of the surface mesh
		Vec3 bb_min, bb_max;
		std::tie(bb_min, bb_max) = geometry::bounding_box(*p.surface_vertex_position_);

		// Compute the size of each division
		Vec3 size = (bb_max - bb_min) / float(step);

		// Sample points within the bounding box using cubic grid approach
		for (int i = 0; i < step; ++i)
		{
			for (int j = 0; j < step; ++j)
			{
				for (int k = 0; k < step; ++k)
				{
					Vec3 sample_point = bb_min + Vec3((i + 0.5f) * size(0), (j + 0.5f) * size(1), (k + 0.5f) * size(2));

					std::pair<uint32, Vec3> cp;
					p.surface_bvh_->closest_point(sample_point, &cp);
					float distance = (cp.second - sample_point).norm();
					float signed_distance = distance;
					if (is_inside(surface, sample_point, p.surface_bvh_faces_, p.surface_vertex_position_, cp))
					{
						signed_distance = -signed_distance;
					}
					// Add sample point to the point mesh
					PointVertex pv = add_vertex(*p.sdf_);
					value<Scalar>(*p.sdf_, p.sdf_sample_radius_, pv) = distance;
					value<Scalar>(*p.sdf_, p.sdf_sample_value_, pv) = signed_distance;
					value<Vec3>(*p.sdf_, p.sdf_sample_position_, pv) = sample_point;
					// Assign color based on signed distance
					Vec4 color;
					if (signed_distance < 0)
					{
						color = Vec4(1.0, 0.0, 0.0, 0.1); // Red for inside
					}
					else
					{
						color = Vec4(0.0, 0.0, 1.0, 0.1); // Blue for outside
					}
					value<Vec4>(*p.sdf_, p.sdf_sample_color_, pv) = color;
				}
			}
		}
		std::cout << "Finished grille sampling" << std::endl;
		points_provider_->emit_connectivity_changed(*p.sdf_);
		points_provider_->emit_attribute_changed(*p.sdf_, p.sdf_sample_position_.get());

		points_provider_->set_mesh_bb_vertex_position(*p.sdf_, p.sdf_sample_position_);
	}
	
	void compute_convex_hull(SURFACE& surface)	
	{
		SurfaceParameters& p = surface_parameters_[&surface];

		if (!p.sdf_)
		{
			std::cout << "No SDF mesh available" << std::endl;
			return;
		}

		std::vector<Vec3> points;
		points.reserve(nb_cells<SurfaceVertex>(*p.surface_));
		foreach_cell(*p.surface_, [&](SurfaceVertex v) -> bool {
			points.push_back(value<Vec3>(*p.surface_, p.surface_vertex_position_, v));
			return true;
		});

		if (!p.flow_mesh_)
		{
			p.flow_mesh_ = surface_provider_->add_mesh(surface_provider_->mesh_name(*selected_surface_) + "_convex_hull");
		}
		clear(*p.flow_mesh_);
		p.flow_vertex_position_ = get_or_add_attribute<Vec3, SurfaceVertex>(*p.flow_mesh_, "position");
		p.flow_vertex_id_ = get_or_add_attribute<uint32, SurfaceVertex>(*p.flow_mesh_, "id");

		cgogn::modeling::convex_hull(points, *p.flow_mesh_, p.flow_vertex_position_.get());

		cgogn::modeling::pliant_remeshing(*p.flow_mesh_, p.flow_vertex_position_, 3, false, true);
		cgogn::modeling::pliant_remeshing(*p.flow_mesh_, p.flow_vertex_position_, 3, false, true);
		uint32 vertex_id = 0;
		foreach_cell(*p.flow_mesh_, [&](SurfaceVertex sv) -> bool {
			value<uint32>(*p.flow_mesh_, p.flow_vertex_id_, sv) = vertex_id++;
			return true;
		});

		surface_provider_->emit_connectivity_changed(*p.flow_mesh_);
		surface_provider_->emit_attribute_changed(*p.flow_mesh_, p.flow_vertex_position_.get());

		surface_provider_->set_mesh_bb_vertex_position(*p.flow_mesh_, p.flow_vertex_position_);
	}

	inline Vec3 barycentric_coordinates(Vec3 const& a, Vec3 const& b, Vec3 const& c, Vec3 const& p)
	{
		/* Derived from the book "Real-Time Collision Detection"
		 * by Christer Ericson published by Morgan Kaufmann in 2005 */
		Vec3 ab = b - a;
		Vec3 ac = c - a;
		Vec3 ap = p - a;

		double d00 = ab.dot(ab);
		double d01 = ab.dot(ac);
		double d11 = ac.dot(ac);
		double d20 = ap.dot(ab);
		double d21 = ap.dot(ac);
		double denom = d00 * d11 - d01 * d01;

		Vec3 bcoords;
		bcoords[1] = (d11 * d20 - d01 * d21) / denom;
		bcoords[2] = (d00 * d21 - d01 * d20) / denom;
		bcoords[0] = 1.0f - bcoords[1] - bcoords[2];

		return bcoords;
	}

	void reach_for_the_sphere_iteration(SURFACE& s)
	{
		SurfaceParameters& p = surface_parameters_[&s];

		uint32 nb_vertices = nb_cells<SurfaceVertex>(*p.flow_mesh_);
		uint32 nb_samples = nb_cells<PointVertex>(*p.sdf_);

		Eigen::SparseMatrix<Scalar, Eigen::ColMajor> A(nb_samples, nb_vertices);
		A.setZero();
		std::cout << "(" << A.rows() << ", " << A.cols() << ") " << std::endl;
		Eigen::SparseMatrix<Scalar, Eigen::ColMajor> M(nb_vertices, nb_vertices);
		Eigen::MatrixXd V(nb_vertices, 3);
		Eigen::MatrixXd S(nb_samples, 3);

		M.setIdentity();

		std::vector<Eigen::Triplet<Scalar>> Acoeffs;
		Acoeffs.reserve(nb_samples * 3);

		foreach_cell(*p.sdf_, [&](PointVertex pv) -> bool {
			// get sdf sample info
			uint32 pv_index = index_of(*p.sdf_, pv);
			Vec3 pv_pos = value<Vec3>(*p.sdf_, p.sdf_sample_position_, pv);
			Scalar sdf = value<Scalar>(*p.sdf_, p.sdf_sample_value_, pv);
			Scalar radius = value<Scalar>(*p.sdf_, p.sdf_sample_radius_, pv);

			// find the closest triangle for each sphere
			std::pair<uint32, Vec3> cp;
			p.flow_mesh_bvh_->closest_point(pv_pos, &cp);
			Scalar distance = std::abs((cp.second - pv_pos).norm() - radius);
			// compute the barycentric coordinate
			auto vertices = incident_vertices(*p.flow_mesh_, p.flow_mesh_bvh_faces_[cp.first]);
			Vec3 v0 = value<Vec3>(*p.flow_mesh_, p.flow_vertex_position_, vertices[0]);
			Vec3 v1 = value<Vec3>(*p.flow_mesh_, p.flow_vertex_position_, vertices[1]);
			Vec3 v2 = value<Vec3>(*p.flow_mesh_, p.flow_vertex_position_, vertices[2]);

			Vec3 bary_coord = barycentric_coordinates(v0, v1, v2, cp.second);

			for (uint32 i = 0; i < 3; i++)
			{
				uint32 sv_index = value<uint32>(*p.flow_mesh_, p.flow_vertex_id_, vertices[i]);
				Acoeffs.push_back(Eigen::Triplet<Scalar>((int)pv_index, (int)sv_index, bary_coord[i]));
			}

			// compute the projection points on the spheres

			Vec3 projection;

			if (distance < 1e-5)
				projection = cp.second;
			else
			{
				bool inside = is_inside(*p.flow_mesh_, pv_pos, p.flow_mesh_bvh_faces_, p.flow_vertex_position_, cp);
				int sigma = (inside == (sdf < 0)) ? 1 : -1;
				projection = pv_pos + sigma * (cp.second - pv_pos).normalized() * radius;
			}
			// Build S matrix
			S.row(pv_index) = projection.transpose();
			return true;
		});

		A.setFromTriplets(Acoeffs.begin(), Acoeffs.end());
		cgogn::index_cells<SurfaceVertex>(*p.flow_mesh_);

		std::cout << "Finished build A and S" << std::endl;
		// Build V matrix
		foreach_cell(*p.flow_mesh_, [&](SurfaceVertex sv) -> bool {
			uint32 sv_index = value<uint32>(*p.flow_mesh_, p.flow_vertex_id_, sv);
			Vec3 sv_pos = value<Vec3>(*p.flow_mesh_, p.flow_vertex_position_, sv);
			V.row(sv_index) = sv_pos.transpose();
			return true;
		});

		Scalar rho = 1.0 / nb_samples;

		std::vector<Eigen::Triplet<Scalar>> triplets;
		triplets.reserve(nb_samples);

		for (int i = 0; i < nb_samples; ++i)
		{
			triplets.emplace_back(i, i, rho);
		}

		Eigen::SparseMatrix<Scalar> R(nb_samples, nb_samples);
		R.setFromTriplets(triplets.begin(), triplets.end());

		Scalar n_c = 0.01;

		// Compute n_p = -A^T * R * (A * V - S)
		Eigen::MatrixXd P = -A.transpose() * R * (A * V - S);

		Scalar numerator = (A * V - S).cwiseProduct(R * (A * P)).sum() + n_c * P.squaredNorm();
		Scalar denominator = (A * P).cwiseProduct(R * (A * P)).sum();
		std::cout << "numerator :" << numerator << ", denominator: " << denominator << std::endl;
		Scalar t = -numerator / denominator; // Avoid division by zero

		// Handle NaN and Inf values
		if (std::isnan(t) || std::isinf(t))
		{
			std::cout << "t not correct" << std::endl;
			t = 0.0;
		}
		std::cout << "computed time step: " << t << std::endl;
		// Clamp t within [min_t, max_t]
		t = std::min(p.t_max, std::max(t, p.t_min));
		std::cout << "selected time step: " << t << std::endl;

		Eigen::SparseMatrix<Scalar> lhs = M + t * A.transpose() * A;
		Eigen::SimplicialCholesky<Eigen::SparseMatrix<Scalar>> chol(lhs);

		if (chol.info() != Eigen::Success)
		{
			std::cerr << "Cholesky decomposition failed!" << std::endl;
		}
		else
		{
			Eigen::MatrixXd V_t = chol.solve(M * V + t * A.transpose() * S);
			std::cout << "Solved the system" << std::endl;
			// Assign the new postition of each vertex
			foreach_cell(*p.flow_mesh_, [&](SurfaceVertex sv) -> bool {
				uint32 sv_index = value<uint32>(*p.flow_mesh_, p.flow_vertex_id_, sv);
				value<Vec3>(*p.flow_mesh_, p.flow_vertex_position_, sv) = V_t.row(sv_index);
				return true;
			});

			uint32 vertex_id = 0;
			foreach_cell(*p.flow_mesh_, [&](SurfaceVertex sv) -> bool {
				value<uint32>(*p.flow_mesh_, p.flow_vertex_id_, sv) = vertex_id++;
				return true;
			});
			build_bvh(*p.flow_mesh_, p.flow_vertex_position_, p.flow_mesh_bvh_, p.flow_mesh_bvh_faces_);
			surface_provider_->emit_attribute_changed(*p.flow_mesh_, p.flow_vertex_position_.get());
		}
	}
	
	void remeshing(SURFACE& s)
	{
		SurfaceParameters& p = surface_parameters_[&s];

		cgogn::modeling::pliant_remeshing(*p.flow_mesh_, p.flow_vertex_position_, 1, false, true, true);
	
		uint32 vertex_id = 0;
		foreach_cell(*p.flow_mesh_, [&](SurfaceVertex sv) -> bool {
			value<uint32>(*p.flow_mesh_, p.flow_vertex_id_, sv) = vertex_id++;
			return true;
		});
		build_bvh(*p.flow_mesh_, p.flow_vertex_position_, p.flow_mesh_bvh_, p.flow_mesh_bvh_faces_);
		surface_provider_->emit_connectivity_changed(*p.flow_mesh_);
		surface_provider_->emit_attribute_changed(*p.flow_mesh_, p.flow_vertex_position_.get());

		surface_provider_->set_mesh_bb_vertex_position(*p.flow_mesh_, p.flow_vertex_position_);
	}

protected:
	void init() override
	{
		surface_provider_ = static_cast<ui::MeshProvider<SURFACE>*>(
			app_.module("MeshProvider (" + std::string{mesh_traits<SURFACE>::name} + ")"));

		points_provider_ = static_cast<ui::MeshProvider<POINT>*>(
			app_.module("MeshProvider (" + std::string{mesh_traits<POINT>::name} + ")"));
	}

	void left_panel() override
	{
		imgui_mesh_selector(surface_provider_, selected_surface_, "Surface",
							[&](SURFACE& s) {selected_surface_ = &s; });

		if (selected_surface_)
		{
			SurfaceParameters& p = surface_parameters_[selected_surface_];

			imgui_combo_attribute<SurfaceVertex, Vec3>(
				*selected_surface_, p.surface_vertex_position_, "Position",
				[&](const std::shared_ptr<SurfaceAttribute<Vec3>>& attribute) { p.surface_vertex_position_ = attribute; });

			if (p.surface_vertex_position_)
			{
				ImGui::RadioButton("Grille", (int*)&p.sampling_mode_, GRILLE);
				ImGui::SameLine();
				ImGui::RadioButton("Random", (int*)&p.sampling_mode_, RANDOM);
				if (p.sampling_mode_ == GRILLE)
				{
					ImGui::InputScalar("Grille resolution", ImGuiDataType_U32, &p.grille_sample_resolution_);
				}
				if (p.sampling_mode_ == RANDOM)
				{
					ImGui::InputScalar("Samples number", ImGuiDataType_U32, &p.random_sample_number_);
				}
				if (ImGui::Button("Init surface data"))
					init_surface_data(*selected_surface_);
				if (p.initialized_)
				{

					if (ImGui::Button("Reach for the sphere"))
						reach_for_the_sphere_iteration(*selected_surface_);
					if (ImGui::Button("Remshing"))
						remeshing(*selected_surface_);
				}
			}
			
			
			
		}
	}

private:
	MeshProvider<POINT>* points_provider_ = nullptr;
	MeshProvider<SURFACE>* surface_provider_ = nullptr;

	SURFACE* selected_surface_ = nullptr;
	std::shared_ptr<SurfaceAttribute<Vec3>> selected_surface_vertex_position_ = nullptr;

	std::unordered_map<const SURFACE*, SurfaceParameters> surface_parameters_;

	std::unordered_map<const SURFACE*, std::vector<std::shared_ptr<boost::synapse::connection>>> surface_connections_;
};

} // namespace ui

} // namespace cgogn

#endif // CGOGN_MODULE_REACH_FOR_THE_SPHERE_
