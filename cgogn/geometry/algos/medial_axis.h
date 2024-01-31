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

const Scalar denoise_preserve = 20.0 * M_PI / 180.0;
const Scalar denoise_planar = 32.0 * M_PI / 180.0;
const Scalar delta_convergence = 1e-5;
const uint32 iteration_limit = 30;

template <typename MESH>
std::tuple<Vec3, Scalar, Vec3> shrinking_ball_center(
	MESH& m, typename mesh_traits<MESH>::Vertex v,
	const typename mesh_traits<MESH>::template Attribute<Vec3>* vertex_position,
	const typename mesh_traits<MESH>::template Attribute<Vec3>* vertex_normal,
	const acc::BVHTree<uint32, Vec3>* surface_bvh, const std::vector<typename mesh_traits<MESH>::Face>& bvh_faces,
	const acc::KDTree<3, uint32>* surface_kdt)
{
	using Vertex = typename mesh_traits<MESH>::Vertex;
	using Face = typename mesh_traits<MESH>::Face;

	const Vec3& p = value<Vec3>(m, vertex_position, v);
	const Vec3& n = value<Vec3>(m, vertex_normal, v);

	uint32 j = 0;
	Scalar r = 0.;

	acc::Ray<Vec3> ray{p, -n, 1e-5, acc::inf};
	acc::BVHTree<uint32, Vec3>::Hit h;
	if (surface_bvh->intersect(ray, &h))
	{
		Face f = bvh_faces[h.idx];
		std::vector<Vertex> vertices = incident_vertices(m, f);
		Vec3 ip = h.bcoords[0] * value<Vec3>(m, vertex_position, vertices[0]) +
				  h.bcoords[1] * value<Vec3>(m, vertex_position, vertices[1]) +
				  h.bcoords[2] * value<Vec3>(m, vertex_position, vertices[2]);
		r = (p - ip).norm() * 0.75;
	}
	// else
	// 	std::cout << "intersection point not found !!!";

	Vec3 c = p - (r * n);
	Vec3 q = p - (2 * r * n);

	while (true)
	{
		// find closest point to c
		std::pair<uint32, Scalar> k_res;
		bool found = surface_kdt->find_nn(c, &k_res);
		// std::pair<uint32, Vec3> cp_res;
		// bool found = surface_bvh->closest_point(c, &cp_res);
		if (!found)
			std::cout << "closest point not found !!!";

		const Vec3& q_next = surface_kdt->vertex(k_res.first);
		Scalar d = k_res.second;
		// Vec3 q_next = cp_res.second;
		// Scalar d = (q_next - c).norm();

		// This should handle all (special) cases where we want to break the loop
		// - normal case when ball no longer shrinks
		// - the case where q == p
		// - any duplicate point cases
		if ((d >= r - delta_convergence) || (p == q_next))
			break;

		// Compute next ball center
		r = compute_radius(p, n, q_next);
		Vec3 c_next = p - (r * n);

		// Denoising
		if (denoise_preserve > 0 || denoise_planar > 0)
		{
			Scalar separation_angle = geometry::angle(p - c_next, q_next - c_next);

			// if (j == 0 && denoise_planar > 0 && separation_angle < denoise_planar)
			// 	break;
			if (j > 0 && denoise_preserve > 0 && (separation_angle < denoise_preserve && r > (q_next - p).norm()))
				break;
		}

		// Stop iteration if this looks like an infinite loop:
		if (j > iteration_limit)
			break;

		c = c_next;
		q = q_next;
		j++;
	}

	return {c, r, q};
}

// adapted from https://github.com/tudelft3d/masbcpp

template <typename MESH, bool ReturnCloestVertexPosition = true>
void shrinking_ball_centers(
	MESH& m, const typename mesh_traits<MESH>::template Attribute<Vec3>* vertex_position,
	const typename mesh_traits<MESH>::template Attribute<Vec3>* vertex_normal,
	typename mesh_traits<MESH>::template Attribute<Vec3>* vertex_shrinking_ball_center,
	typename mesh_traits<MESH>::template Attribute<Scalar>* vertex_shrinking_ball_radius,
	typename std::conditional<
		ReturnCloestVertexPosition, 
		typename mesh_traits<MESH>::template Attribute<std::pair<Vec3, Vec3>>,
		typename mesh_traits<MESH>::template Attribute<
			std::pair<typename mesh_traits<MESH>::Vertex, typename mesh_traits<MESH>::Vertex>>>::type*
		vertex_shrinking_ball_closest_points)
{
	using Vertex = typename mesh_traits<MESH>::Vertex;
	using Face = typename mesh_traits<MESH>::Face;

	uint32 nb_vertices = nb_cells<Vertex>(m);
	uint32 nb_faces = nb_cells<Face>(m);

	auto bvh_vertex_index = add_attribute<uint32, Vertex>(m, "__bvh_vertex_index");

	std::vector<Vec3> vertex_position_vector;
	vertex_position_vector.reserve(nb_vertices);
	std::vector<Vertex> kdt_vertices;
	if constexpr (!ReturnCloestVertexPosition)
	{
		kdt_vertices.reserve(nb_vertices);
	}
	uint32 idx = 0;
	foreach_cell(m, [&](Vertex v) -> bool {
		value<uint32>(m, bvh_vertex_index, v) = idx++;
		vertex_position_vector.push_back(value<Vec3>(m, vertex_position, v));
		if constexpr (!ReturnCloestVertexPosition)
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
		auto [c, r, q] =
			shrinking_ball_center(m, v, vertex_position, vertex_normal, surface_bvh, bvh_faces, surface_kdt);
		value<Vec3>(m, vertex_shrinking_ball_center, v) = c;
		value<Scalar>(m, vertex_shrinking_ball_radius, v) = r;
		if constexpr (ReturnCloestVertexPosition)
		{
			value<std::pair<Vec3, Vec3>>(m, vertex_shrinking_ball_closest_points,
										 v) = {value<Vec3>(m, vertex_position, v), q};
		}
		else
		{
			std::pair<uint32, Scalar> k_res;
			bool found = surface_kdt->find_nn(q, &k_res);
			if (found)
			{
				Vertex closest_vertex = kdt_vertices[k_res.first];
				value<std::pair<Vertex, Vertex>>(m, vertex_shrinking_ball_closest_points, v) = {v, closest_vertex};
			}
			else
			{
				std::cout << "closest point not found !!!";
			}
		}
		return true;
	});

	delete surface_kdt;

	remove_attribute<Vertex>(m, bvh_vertex_index);
	delete surface_bvh;
}


template <typename MESH, bool UseDisMatrix = true>
typename std::enable_if<UseDisMatrix, Scalar>::type surface_medial_distance_variance(
	MESH& m, const typename mesh_traits<MESH>::Vertex medial_axis_sample,
	const typename mesh_traits<MESH>::template Attribute<Vec3>* vertex_position,
	const typename mesh_traits<MESH>::template Attribute<Vec3>* vertex_normal,
	const typename mesh_traits<MESH>::template Attribute<Vec3>* vertex_shrinking_ball_center,
	std::vector<typename mesh_traits<MESH>::Vertex>& clusters, Eigen::MatrixXd& dis_matrix)
{
	using Vertex = typename mesh_traits<MESH>::Vertex;
	std::vector<Scalar> cosines;
	Scalar sum_dist = 0;
	for (Vertex& v : clusters)
	{
		Vec3 dir =
			(value<Vec3>(m, vertex_position, v) - value<Vec3>(m, vertex_shrinking_ball_center, medial_axis_sample))
				.normalized();
		Scalar cosine = value<Vec3>(m, vertex_normal, v).dot(dir);
		cosines.push_back(cosine);
		sum_dist += -cosine * dis_matrix(index_of(m, v), index_of(m, medial_axis_sample));
	}
	Scalar average_dist = sum_dist / clusters.size();
	Scalar variance = 0;
	size_t idx = 0;
	for (Vertex& v : clusters)
	{
		Scalar cosine = cosines[idx++];
		variance += (-cosine * dis_matrix(index_of(m, v), index_of(m, medial_axis_sample)) - average_dist) *
					(-cosine * dis_matrix(index_of(m, v), index_of(m, medial_axis_sample)) - average_dist);
	}
	return variance / clusters.size();
}

template <typename MESH1, typename MESH2> 
Scalar surface_medial_distance_variance(
	MESH1& m1, MESH2& m2, const typename mesh_traits<MESH2>::Vertex medial_axis_sample,
	const typename mesh_traits<MESH1>::template Attribute<Vec3>* vertex_position,
	const typename mesh_traits<MESH1>::template Attribute<Vec3>* vertex_normal,
	const typename mesh_traits<MESH2>::template Attribute<Vec3>* vertex_cluster_center,
	const typename mesh_traits<MESH1>::template Attribute<Scalar>* weight,
	std::vector<typename mesh_traits<MESH1>::Vertex>& clusters_points, Scalar median_sphere_radius)
{
	using PointVertex = typename mesh_traits<MESH2>::Vertex;
	using SurfaceVertex = typename mesh_traits<MESH1>::Vertex;
	std::vector<Scalar> dist;
	Scalar sum_dist = 0;
	for (SurfaceVertex& v : clusters_points)
	{
		Vec3 vec = (value<Vec3>(m1, vertex_position, v) - value<Vec3>(m2, vertex_cluster_center, medial_axis_sample));
		Vec3 dir = vec.normalized();
		/*Scalar cosine = value<Vec3>(m1, vertex_normal, v).dot(dir);*/
		Scalar length = vec.norm() - median_sphere_radius;
		Scalar distance = length * value<Scalar>(m1, weight, v);
		dist.push_back(distance);
		sum_dist += distance;
		
	}
	Scalar average_dist = sum_dist / clusters_points.size();
	Scalar variance = 0;
	size_t idx = 0;
	for (SurfaceVertex& v : clusters_points)
	{
		Scalar distance = dist[idx++];
		variance += (distance - average_dist) * (distance - average_dist);
	}
	return variance / clusters_points.size();
}

template <typename MESH>
	std::tuple < Scalar, 
		typename mesh_traits<MESH>::Vertex,
		typename mesh_traits<MESH>::Vertex> move_point_to_medial_axis(
	MESH& mesh, 
	const typename mesh_traits<MESH>::template Attribute<Vec3>* vertex_position,
	const typename mesh_traits<MESH>::template Attribute<Vec3>* vertex_normal,
	const std::vector<typename mesh_traits<MESH>::Vertex>& vertices, 
	Vec3& pos, const acc::KDTree<3, uint32>* surface_kdt,
	const acc::KDTree<3, uint32>* medial_kdt, 
	acc::BVHTree<uint32, Vec3>* surface_bvh)
{
	using Vertex = typename mesh_traits<MESH>::Vertex;
	Scalar distance_to_nearest = 0;
	Scalar distance_to_new_nearest = 0;
	Vertex nearest_medial_vertex, other_nearest_medial_vertex;
	Vec3 nearset_point_pos, new_nearest_point_pos;
	Vec3 old_dir, new_dir;
	Scalar cosine;
	Scalar step = 0.02;
	bool found = false;
	int iteration = 0;
	do{
		iteration++;
		if (iteration > 200)
		{
			std::cout << "no solution" << std::endl;
			break;
		}
		nearset_point_pos = surface_bvh->closest_point(pos);
		distance_to_nearest = (nearset_point_pos - pos).norm();
		old_dir = (pos - nearset_point_pos)/*.normalized()*/;
		std::pair<uint32, Scalar> k_res;
		found = surface_kdt->find_nn(nearset_point_pos, &k_res);
		if (found)
		{
			nearest_medial_vertex = vertices[k_res.first];
		}
		else
		{
			std::cout << "closest point not found !!!";
		}
		pos += step * old_dir;

		new_nearest_point_pos = surface_bvh->closest_point(pos);
		std::pair<uint32, Scalar> new_k_res;
		found = surface_kdt->find_nn(new_nearest_point_pos, &new_k_res);
		if (found)
		{
			other_nearest_medial_vertex = vertices[new_k_res.first];
		}
		else
		{
			std::cout << "closest point not found !!!";
		}
		distance_to_new_nearest = (new_nearest_point_pos - pos).norm();
		new_dir = (pos - new_nearest_point_pos)/*.normalized()*/;
		cosine = old_dir.dot(new_dir) / new_dir.norm()/old_dir.norm();
		/*if (cosine<0)
		{
			step *= 0.1;
		}*/
		if (new_k_res.first!= k_res.first)
		{
			step *= 0.1;
		}
		/*std::cout << "distance_to_nearest: " << distance_to_nearest << ", "
				  << "distance_to_new_nearest: " << distance_to_new_nearest << std::endl;
		
		std::cout << "cosine: " << cosine << ", iteration: "<<iteration
				  << std::endl;
		std::cout << "cosine2: " << old_dir.normalized().dot(value<Vec3>(mesh, vertex_normal, nearest_medial_vertex))
				  << std::endl;*/
	} while (/*cosine > 0||*/ 
			 std::fabs(distance_to_nearest - distance_to_new_nearest) > 1e-5);
	std::pair<uint32, Scalar> k_res;
	/*bool found = medial_kdt->find_nn(pos, &k_res);
	if (!found)
	{
		std::cout << "Closest point not found !!!" << std::endl;
		return {-1, nearest_medial_vertex};
	}
	nearest_medial_vertex = vertices[k_res.first];*/
	//std::cout << "---------------------------------------------------" << std::endl;
	return
	{
		distance_to_nearest, nearest_medial_vertex, other_nearest_medial_vertex};
}
} // namespace geometry

} // namespace cgogn

#endif // CGOGN_GEOMETRY_ALGOS_MEDIAL_AXIS_H_
