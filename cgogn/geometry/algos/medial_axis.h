/*******************************************************************************
 * CGoGN: Combinatorial and Geometric modeling with Generic N-dimensional Maps  *
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

#ifndef CGOGN_GEOMETRY_ALGOS_MEDIAL_AXIS_H_
#define CGOGN_GEOMETRY_ALGOS_MEDIAL_AXIS_H_

#include <cgogn/core/functions/traversals/vertex.h>

#include <cgogn/geometry/functions/angle.h>
#include <cgogn/geometry/types/vector_traits.h>

#include <array>
#include <vector>
#include <nlopt.hpp>
#include <libacc/bvh_tree.h>
#include <libacc/kd_tree.h>

namespace cgogn
{

namespace geometry
{

using geometry::Scalar;
using geometry::Vec3;

inline Scalar compute_radius(const Vec3& p, const Vec3& n, const Vec3& q)
{
	// Compute radius of the ball that touches points p and q and whose center falls on the normal n from p
	Vec3 qp = p - q;
	Scalar d = qp.norm();
	// Scalar cos_theta = n.dot(p - q) / d;
	Scalar cos_theta = geometry::cos_angle(n, qp);
	return Scalar(d / (2 * cos_theta));
}

inline Vec3 projection_point_on_triangle(const Vec3& p, const Vec3 p1, const Vec3 p2, const Vec3 p3)
{
	Vec3 e1 = p1 - p2;
	Vec3 e2 = p1 - p3;
	Vec3 norm = e1.cross(e2);
	norm.normalize();
	return  p - (p - p2).dot(norm) * norm;
}
bool SameSide(const Vec3& v1, const Vec3& v2, const Vec3& a, const Vec3& b)
{
	Vec3 cv1 = (b - a).cross(v1 - a);
	Vec3 cv2 = (b - a).cross(v2 - a);
	if (cv1.dot(cv2) >= 1e-12)
		return true;
	else
		return false;
}

bool InsideTriangle(const Vec3& p, const Vec3& v0, const Vec3& v1, const Vec3& v2)
{
	if (SameSide(p, v0, v1, v2) && SameSide(p, v1, v0, v2) && SameSide(p, v2, v0, v1))
		return true;
	else
		return false;
}
const Scalar denoise_preserve = 30.0 * M_PI / 180.0;
const Scalar denoise_planar = 32.0 * M_PI / 180.0;
const Scalar delta_convergence = 1e-5;
const uint32 iteration_limit = 30;

template <typename MESH>
std::tuple<Vec3, Scalar, typename mesh_traits<MESH>::Vertex> shrinking_ball_center(
	const MESH& m, const Vec3& p, const Vec3& n,
	const typename mesh_traits<MESH>::template Attribute<Vec3>* vertex_position,
	const acc::BVHTree<uint32, Vec3>* surface_bvh, const std::vector<typename mesh_traits<MESH>::Face>& bvh_faces,
	const acc::KDTree<3, uint32>* surface_kdt, const std::vector<typename mesh_traits<MESH>::Vertex>& kdt_vertices)
{
	using Vertex = typename mesh_traits<MESH>::Vertex;
	using Face = typename mesh_traits<MESH>::Face;

	uint32 j = 0;
	Scalar r = 0.;

	acc::Ray<Vec3> ray{p, -n, 1e-10, acc::inf};
	acc::BVHTree<uint32, Vec3>::Hit h1;
	acc::BVHTree<uint32, Vec3>::Hit h2;
	if (surface_bvh->intersect(ray, &h1))
	{
		Face f = bvh_faces[h1.idx];
		std::vector<Vertex> vertices = incident_vertices(m, f);
		Vec3 ip = h1.bcoords[0] * value<Vec3>(m, vertex_position, vertices[0]) +
				  h1.bcoords[1] * value<Vec3>(m, vertex_position, vertices[1]) +
				  h1.bcoords[2] * value<Vec3>(m, vertex_position, vertices[2]);
		r = (p - ip).norm() * 0.75;
	}
	else
	{
		std::cout << "intersection point not found !!!";
		r = 0;
		return
		{
			Vec3(0, 0, 0), 0, Vertex()
		};
	}
	Vec3 c = p - (r * n);
	Vec3 q = p - (2 * r * n);
	Vertex q_v;

	while (true)
	{
		// Find closest point to c
		Vertex q_next_v;
		Vec3 q_next;
		Scalar d = std::numeric_limits<Scalar>::max();
		std::pair<uint32, Scalar> k_res;
		if (!surface_kdt->find_nn(c, &k_res))
		{
			std::cout << "closest point not found !!!";
			return
			{
				Vec3(0, 0, 0), 0, Vertex()
			};
		}
		else
		{
			q_next_v = kdt_vertices[k_res.first];
			foreach_incident_face(m, q_next_v, [&](Face f) {
				std::vector<Vertex> vertices = incident_vertices(m, f);
				Vec3 p1 = value<Vec3>(m, vertex_position, vertices[0]);
				Vec3 p2 = value<Vec3>(m, vertex_position, vertices[1]);
				Vec3 p3 = value<Vec3>(m, vertex_position, vertices[2]);
				Vec3 proj = projection_point_on_triangle(c, p1, p2, p3);
				Scalar dis = (proj - c).norm();
				if (dis < d)
				{
					d = dis;
					q_next = proj;
				}
				
				return true;
			});
		}
		/*const Vec3& q_next = surface_kdt->vertex(k_res.first);
		Scalar d = k_res.second;
		Vertex q_next_v = kdt_vertices[k_res.first];*/

		 /* std::pair<uint32, Vec3> cp_res;
		 surface_bvh->closest_point(c, &cp_res);
		 Vec3 q_next = cp_res.second;
		 Scalar d = (q_next - c).norm();
		 Vertex q_next_v;

		ray = acc::Ray<Vec3>{c, q_next - c, 1e-10, acc::inf};
		 if(surface_bvh->intersect(ray, &h2)){
			Face f = bvh_faces[h2.idx];
			std::vector<Vertex> vertices = incident_vertices(m, f);
			int max_idx = 0;
			if (h2.bcoords[1] > h2.bcoords[max_idx])
				max_idx = 1;
			if (h2.bcoords[2] > h2.bcoords[max_idx])
				max_idx = 2;
			q_next_v = vertices[max_idx];
		 }
		 else
		 {
			std::cout << "intersection point of the center not found !!!";
			 return
			 {
				 Vec3(0, 0, 0), 0, Vertex()
			 };
		 }*/
		// If the closest point is (almost) the same as the previous one, or if the ball no longer shrinks, we stop
		if ((d >= r - delta_convergence) || (q_next - q).norm() < delta_convergence)
		{
			if (j == 0)
			{
				std::cout << "the closest point is (almost) the same as the previous one, or if the ball no longer shrinks, we stop!!!"<<std::endl;
				q_v = q_next_v;
			}
			break;
		}
		// Compute next ball center
		Scalar r_next = compute_radius(p, n, q_next);
		Vec3 c_next = p - (r_next * n);

		// Denoising
		Scalar separation_angle = geometry::angle(p - c_next, q_next - c_next);
		if (j > 0 && separation_angle < denoise_preserve)
		{ // && r_next > // (q_next - p).norm())
			// std::cout << "Denoising preserve" << std::endl;
			break;
		}

		c = c_next;
		r = r_next;
		q = q_next;
		q_v = q_next_v;

		j++;
		
		if (j > iteration_limit)
			break;
	}

	return {c, r, q_v};
}

// adapted from https://github.com/tudelft3d/masbcpp

template <typename MESH>
void shrinking_ball_centers(
	MESH& m, const typename mesh_traits<MESH>::template Attribute<Vec3>* vertex_position,
	const typename mesh_traits<MESH>::template Attribute<Vec3>* vertex_normal,
	typename mesh_traits<MESH>::template Attribute<Vec3>* vertex_shrinking_ball_center,
	typename mesh_traits<MESH>::template Attribute<Scalar>* vertex_shrinking_ball_radius,
typename mesh_traits<MESH>::template Attribute<std::pair<typename mesh_traits<MESH>::Vertex, 
	typename mesh_traits<MESH>::Vertex>>* vertex_shrinking_ball_closest_points)
{
	using Vertex = typename mesh_traits<MESH>::Vertex;
	using Face = typename mesh_traits<MESH>::Face;

	uint32 nb_vertices = nb_cells<Vertex>(m);
	uint32 nb_faces = nb_cells<Face>(m);

	auto bvh_vertex_index = get_or_add_attribute<uint32, Vertex>(m, "__bvh_vertex_index");

	std::vector<Vec3> vertex_position_vector;
	vertex_position_vector.reserve(nb_vertices);
	std::vector<Vertex> kdt_vertices;

	kdt_vertices.reserve(nb_vertices);

	uint32 idx = 0;
	foreach_cell(m, [&](Vertex v) -> bool {
		value<uint32>(m, bvh_vertex_index, v) = idx++;
		vertex_position_vector.push_back(value<Vec3>(m, vertex_position, v));
		kdt_vertices.push_back(v);
		return true;
	});

	std::vector<Face> bvh_faces;
	bvh_faces.reserve(nb_faces);
	std::vector<uint32> face_vertex_indices;
	face_vertex_indices.reserve(nb_faces * 3);
	foreach_cell(m, [&](Face f) -> bool {
		bvh_faces.push_back(f);
		foreach_incident_vertex(m, f, [&](Vertex v) -> bool {
			face_vertex_indices.push_back(value<uint32>(m, bvh_vertex_index, v));
			return true;
		});
		return true;
	});

	acc::BVHTree<uint32, Vec3>* surface_bvh =
		new acc::BVHTree<uint32, Vec3>(face_vertex_indices, vertex_position_vector);

	acc::KDTree<3, uint32>* surface_kdt = new acc::KDTree<3, uint32>(vertex_position_vector);

	parallel_foreach_cell(m, [&](Vertex v) -> bool {
		auto [c, r, q] = shrinking_ball_center(m, value<Vec3>(m, vertex_position, v), value<Vec3>(m, vertex_normal, v),
												vertex_position, surface_bvh, bvh_faces, surface_kdt, kdt_vertices);
		value<Vec3>(m, vertex_shrinking_ball_center, v) = c;
		value<Scalar>(m, vertex_shrinking_ball_radius, v) = r;
		value<std::pair<Vertex, Vertex>>(m, vertex_shrinking_ball_closest_points, v) = {v, q};
		return true;
	});

	delete surface_kdt;

	remove_attribute<Vertex>(m, bvh_vertex_index);
	delete surface_bvh;
}


template <typename MESH>
std::tuple<Scalar, typename mesh_traits<MESH>::Vertex, typename mesh_traits<MESH>::Vertex> non_linear_solver(
	MESH& mesh, const typename mesh_traits<MESH>::template Attribute<Vec3>* vertex_position,
	const typename mesh_traits<MESH>::template Attribute<Vec3>* vertex_normal,
	const std::vector<typename mesh_traits<MESH>::Vertex>& vertices, Vec3& pos, Scalar radius)
{
	struct OptimizationData
	{
		const MESH* mesh;
		const typename mesh_traits<MESH>::template Attribute<Vec3>* vertex_position;
		std::vector<typename mesh_traits<MESH>::Vertex> vertices;
		Vec3 p;
		Scalar rad;
	};
	/*struct ConstraintData
	{
		const MESH* mesh;
		const typename mesh_traits<MESH>::template Attribute<Vec3>* vertex_position;
		const typename mesh_traits<MESH>::template Attribute<Vec3>* vertex_position;
		typename mesh_traits<MESH>::Vertex v;
		
	};*/
	using Vertex = typename mesh_traits<MESH>::Vertex;

	auto objective = [](const std::vector<double>& x, std::vector<double>& grad, void* data) -> double { 
		auto *optData = reinterpret_cast<OptimizationData*>(data);
		double energy = 0;
		Vec3 f;
		f.setZero();
		Scalar dEdr = 0;
		// compute the mean distance to the sphere center
		for (const auto& v : optData->vertices)
		{
			Vec3 position = value<Vec3>(*optData->mesh, optData->vertex_position, v);
			double distance = (position - Vec3(x[0], x[1], x[2])).norm() - x[3];	
			energy += distance * distance;
			//compute gradient
			Vec3 x01 = value<Vec3>(*optData->mesh, optData->vertex_position, v) - Vec3(x[0], x[1], x[2]);
			Scalar norm = x01.norm();
			Vec3 x01_normalised = x01.normalized();
			f -= 2 * (norm - x[3]) * x01_normalised;
			dEdr -= 2 * (norm - x[3]);
		}
		if (!grad.empty())
		{
			grad[0] = f.x();
			grad[1] = f.y();
			grad[2] = f.z();
			grad[3] = dEdr;
		}
		optData->p = Vec3(x[0], x[1], x[2]);
		optData->rad = x[3];
		std::cout << "current position: " << x[0] << ", " << x[1] << ", " << x[2] << ", "
				  << "current radius: " << x[3] << std::endl;
		std::cout << "current energy: " << energy << std::endl;
		return energy;
	}; 
	
	nlopt::opt opt(nlopt::LD_SLSQP, 4);
	OptimizationData data;
	data.mesh = &mesh;
	data.vertex_position = vertex_position;
	data.vertices = vertices;
	opt.set_min_objective(objective, &data);

	auto constraint = [](const std::vector<double>& x, std::vector<double>& grad, void* data) -> double {
		OptimizationData* constraintData = reinterpret_cast<OptimizationData*>(data);

		Vec3 center(x[0], x[1], x[2]); // 球心位置
		double r = x[3];			   // 球半径

		double minDistance = std::numeric_limits<double>::max();

		Vec3 cloest_vertex;
		for (const auto& v : constraintData->vertices)
		{
			Vec3 position = value<Vec3>(*constraintData->mesh, constraintData->vertex_position, v);
			double distance = (position - center).norm();
			if (distance < minDistance)
			{
				minDistance = distance;
				cloest_vertex = position;
			}
		}
		if (!grad.empty())
		{
			Vec3 gradient = (cloest_vertex - center).normalized();
			grad[0] = -gradient.x();
			grad[1] = -gradient.y();
			grad[2] = -gradient.z();
			grad[3] = 1;
		}
		return r - minDistance;
	};
	
	opt.add_inequality_constraint(constraint, &data, 1e-8);
	opt.set_ftol_rel(1e-5);
	opt.set_upper_bounds({1, 1, 1, 1});
	opt.set_lower_bounds({-1, -1, -1, 0});
	/* opt.set_maxtime(0.1);
	*/
		
	std::vector<double> x = {pos.x(), pos.y(), pos.z(), radius};  
	double min_f;				   

	try
	{
		nlopt::result result = opt.optimize(x, min_f);
		// 输出结果
		std::cout << "after optimization: " << x[0] << ", " << x[1] << ", " << x[2] << ", " << x[3]
				  << 
				  ", with min energy:" << min_f<< std::endl;
		pos = data.p;
		std::cout << "after optimization: " << pos.x() << ", " << pos.y() << ", " << pos.z() << ", " << data.rad
				  << ", with min energy:" << min_f << std::endl;
		return std::make_tuple(data.rad, vertices[0], vertices[1]);

	}
	catch (std::exception& e)
	{
		std::cerr << "NLopt failed: " << e.what() << std::endl;
		return std::make_tuple(0.0, vertices[0], vertices[1]); 
	}
}





	} // namespace geometry

} // namespace cgogn

#endif // CGOGN_GEOMETRY_ALGOS_MEDIAL_AXIS_H_
