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

#ifndef CGOGN_MODELING_DECIMATION_EDGE_QUEUE_SQEM_H_
#define CGOGN_MODELING_DECIMATION_EDGE_QUEUE_SQEM_H_

#include <cgogn/core/functions/attributes.h>
#include <cgogn/core/functions/traversals/global.h>
#include <cgogn/core/functions/traversals/vertex.h>
#include <cgogn/core/functions/traversals/edge.h>
#include <cgogn/core/types/incidence_graph/incidence_graph_ops.h>

#include <cgogn/geometry/types/slab_quadric.h>
#include <cgogn/geometry/types/vector_traits.h>
#include <cgogn/geometry/algos/invertion.h>

namespace cgogn
{

namespace modeling
{

template <typename NONMANIFOLD>
struct DecimationSQEM_Helper
{
	template <typename T>
	using Attribute = typename mesh_traits<NONMANIFOLD>::template Attribute<T>;

	using Vertex = typename NONMANIFOLD::Vertex;
	using Edge = typename NONMANIFOLD::Edge;
	using Face = typename NONMANIFOLD::Face;

	using Vec3 = typename geometry::Vec3;
	using Vec4 = typename geometry::Vec4;
	using Scalar = typename geometry::Scalar;
	using Slab_Quadric = typename geometry::Slab_Quadric;
	using EdgeQueue = std::multimap<Scalar, Edge>;
	using EdgeQueueIt = typename EdgeQueue::const_iterator;
	using EdgeInfo = std::pair<bool, EdgeQueueIt>; // {valid, iterator}
	
	friend class Slab_Quadric;

	DecimationSQEM_Helper(float k, NONMANIFOLD& m, std::shared_ptr<Attribute<Vec3>>& position,
						   std::shared_ptr<Attribute<Vec4>> & sphere_info,
						  std::shared_ptr<Attribute<Vec3>>& stability_color,
		std::shared_ptr<Attribute<double>>& stability_ratio)
		:k_(k), m_(m), position_(position), sphere_info_(sphere_info), stability_color_(stability_color), stability_ratio_(stability_ratio)
	{
		edge_queue_it_ = add_attribute<EdgeInfo, Edge>(m_, "non_manifold_edge_queue_it");
		sphere_opt_ = add_attribute<Vec4, Edge>(m_, "sphere_opt");
		vertex_slab_quadric_ = add_attribute<Slab_Quadric, Vertex>(m_, "vertex_slab_quadric");
		edge_slab_quadric_ = add_attribute<Slab_Quadric, Edge>(m_, "edge_slab_quadric");
		slab_normals_ = add_attribute<std::pair<Vec4, Vec4>, Face>(m_, "slab_normals");
		queue_.clear();
	}

	~DecimationSQEM_Helper()
	{
		remove_attribute<Edge>(m_, edge_queue_it_);
		remove_attribute<Edge>(m_, sphere_opt_);
		remove_attribute<Vertex>(m_, vertex_slab_quadric_);
		remove_attribute<Edge>(m_, edge_slab_quadric_);
		remove_attribute<Face>(m_, slab_normals_);

	}

	void initial_slab_mesh()
	{
		parallel_foreach_cell(m_, [&](Vertex v) -> bool {
			value<Slab_Quadric>(m_, vertex_slab_quadric_, v).clear();
			value<Slab_Quadric>(m_, edge_slab_quadric_, v).clear();
			return true;
		});
		foreach_cell(m_, [&](Face f) -> bool {
			Vec4 n1 = Vec4(0, 0, 0, 1), n2 = Vec4(0, 0, 0, 1);
			auto ifv = incident_vertices(m_, f);
			Vec4 p1 = value<Vec4>(m_, sphere_info_, ifv[0]);
			Vec4 p2 = value<Vec4>(m_, sphere_info_, ifv[1]);
			Vec4 p3 = value<Vec4>(m_, sphere_info_, ifv[2]);
			if (slab_normal(p1, p2, p3, n1, n2) && n1 != Vec4(0, 0, 0, 1) && n2 != Vec4(0, 0, 0, 1))
			{
				value<std::pair<Vec4, Vec4>>(m_, slab_normals_, f) = {n1, n2};
				Slab_Quadric q1(p1, p2, p3, n1, n2);
				Slab_Quadric q2(p2, p1, p3, n1, n2);
				Slab_Quadric q3(p3, p1, p2, n1, n2);
				value<Slab_Quadric>(m_, vertex_slab_quadric_, ifv[0]) += q1;
				value<Slab_Quadric>(m_, vertex_slab_quadric_, ifv[1]) += q2;
				value<Slab_Quadric>(m_, vertex_slab_quadric_, ifv[2]) += q3;
			}
			return true;
		});
	}

	void initial_boundary_mesh()
	{
		foreach_cell(m_, [&](Edge e) {
			if (is_boundary(m_, e))
			{
				auto ief = incident_faces(m_, e);
				auto iev = incident_vertices(m_, e);
				auto iv1_e = incident_edges(m_, iev[0]);
				auto iv2_e = incident_edges(m_, iev[1]);
				double sr = value<double>(m_, stability_ratio_, e);
				Vec4 p1 = value<Vec4>(m_, sphere_info_, iev[0]);
				Vec4 p2 = value<Vec4>(m_, sphere_info_, iev[1]);
				Vec3 v1 = p1.head<3>();
				Vec3 v2 = p2.head<3>();
				Vec3 v1v2 = v2 - v1;
				double r1 = p1.w();
				double r2 = p2.w();
				Vec3 v3;
				if (ief.size() == 1)
				{
					auto iefv = incident_vertices(m_, ief[0]);

					for (auto vert : iefv)
					{
						auto vert_pos = value<Vec4>(m_, sphere_info_, vert).head<3>();
						if (vert_pos != v1 && vert_pos != v2)
						{
							v3 = vert_pos;
							break;
						}
					}
					Vec3 v1v3 = v3 - v1;
					Vec3 iefv_nor = v1v2.cross(v1v3);
					Vec3 normal_boundary = iefv_nor.cross(v1v2);
					double temp_cos = acos(normal_boundary.dot(v1v3) / v1v3.norm() / normal_boundary.norm());
					bool dir = temp_cos > M_PI / 2.0 ? true : false;
					if (!dir)
						normal_boundary *= -1;
					Vec4 n1 = Vec4(normal_boundary.x(), normal_boundary.y(), normal_boundary.z(), 1.0);
					Slab_Quadric boundary_v1(p1, n1, sr);
					Slab_Quadric boundary_v2(p2, n1, sr);
					value<Slab_Quadric>(m_, vertex_slab_quadric_, iev[0]) += boundary_v1;
					value<Slab_Quadric>(m_, vertex_slab_quadric_, iev[1]) += boundary_v2;
				}
				else if (ief.size() == 0)
				{
					Vec4 p3 = Vec4(v1.x() + v2.x(), v1.y() + v2.y(), v1.z() + v2.z(), (r1 + r2) * 0.5);
					Vec4 n1 = Vec4(0, 0, 0, 1), n2 = Vec4(0, 0, 0, 1);
					if (slab_normal(p1, p2, p3, n1, n2))
					{
						if (n1 != Vec4(0, 0, 0, 1) && n2 != Vec4(0, 0, 0, 1))
						{
							// Add slab quadric of current slab
							Slab_Quadric q1(p1, p2, p3, n1, n2, true);
							Slab_Quadric q2(p2, p1, p3, n1, n2, true);
							value<Slab_Quadric>(m_, vertex_slab_quadric_, iev[0]) += q1;
							value<Slab_Quadric>(m_, vertex_slab_quadric_, iev[1]) += q2;
							Vec3 t1 = v1v2.cross(n1.head<3>());
							Vec3 t2 = v1v2.cross(n2.head<3>());
							// Add slab quadric of two orthogonal plane of the slab
							Vec4 tnormal1 = Vec4(t1.x(), t1.y(), t1.z(), 1.0);
							Vec4 tnormal2 = Vec4(t2.x(), t2.y(), t2.z(), 1.0);
							Slab_Quadric q3(p1, p2, p3, tnormal1, tnormal2, true);
							Slab_Quadric q4(p2, p1, p3, tnormal1, tnormal2, true);
							value<Slab_Quadric>(m_, vertex_slab_quadric_, iev[0]) += q3;
							value<Slab_Quadric>(m_, vertex_slab_quadric_, iev[1]) += q4;
						}
					}
					if (iv1_e.size() == 1)
					{
						Vec4 boundary_vertex_norm = Vec4(v1v2.x(), v1v2.y(), v1v2.z(), 1.0);
						Slab_Quadric boundary_v1(p1, boundary_vertex_norm, sr);
						value<Slab_Quadric>(m_, vertex_slab_quadric_, iev[0]) += boundary_v1;
					}
					if (iv2_e.size() == 1)
					{
						Vec4 boundary_vertex_norm = Vec4(-v1v2.x(), -v1v2.y(), -v1v2.z(), 1.0);
						Slab_Quadric boundary_v2(p2, boundary_vertex_norm, sr);
						value<Slab_Quadric>(m_, vertex_slab_quadric_, iev[1]) += boundary_v2;
					}
				}
			}
			return true;
			});
	}

	void initial_collapse_queue()
	{
		// Initialize the queue with all the edges and their cost
		foreach_cell(m_, [&](Edge e) -> bool {
			Vec4 opt;
			Scalar cost_opt = edge_optimal(e,opt);
			value<Vec4>(m_, sphere_opt_, e) = opt;
			
			/* double cost_opt = value<double>(nm, stability_ratio, e);
			std::cout << std::fixed;
			std::cout  << "cost: " << std::setprecision(9)<< cost;
			std::cout << ", cost optimal: " << std::setprecision(9) <<cost_opt;
			std::cout << ", stability ratio: " << std::setprecision(9)<< value<double>(nm, stability_ratio, e) <<
			std::endl;*/

			value<EdgeInfo>(m_, edge_queue_it_, e) = {true, queue_.emplace(cost_opt, e)};
			return true;
		});
	}

	void simplify(int threshold)
	{
		int count = 0;
		while (!queue_.empty() && count < threshold)
		{
			auto it = queue_.begin();
			Edge e = (*it).second;
			queue_.erase(it);
			value<EdgeInfo>(m_, edge_queue_it_, e).first = false;

			Vec4 opt = value<Vec4>(m_, sphere_opt_, e);
			Slab_Quadric eq = value<Slab_Quadric>(m_, edge_slab_quadric_, e);
			auto [v, removed_edges] = collapse_edge_qmat(m_, e);
			value<Vec4>(m_, sphere_info_, v) = opt;
			value<Slab_Quadric>(m_, vertex_slab_quadric_, v) = eq;
			for (Edge re : removed_edges)
			{
				EdgeInfo einfo = value<EdgeInfo>(m_, edge_queue_it_, re);
				if (einfo.first)
					queue_.erase(einfo.second);
			}

			// recompute the cost of the edges incident to v and update the queue
			foreach_incident_edge(m_, v, [&](Edge ie) -> bool {
				auto iv = incident_vertices(m_, ie);
				Vertex v1 = iv[0];
				Vertex v2 = iv[1];
				const Vec3& v1_p = value<Vec4>(m_, sphere_info_, v1).head<3>();
				const Vec3& v2_p = value<Vec4>(m_, sphere_info_, v2).head<3>();
				const double& r1 = value<Vec4>(m_, sphere_info_, v1).w();
				const double& r2 = value<Vec4>(m_, sphere_info_, v2).w();
				const double center_dist = (v1_p - v2_p).norm();
				double dis = std::max(0.0, (center_dist - std::abs(r1 - r2)));
				if (center_dist == 0.0)
				{
					(*stability_ratio_)[ie.index_] = 0.0;
					(*stability_color_)[ie.index_] = Vec3(0, 0, 0.5);
				}
				else
				{
					double stability = dis / center_dist;
					value<double>(m_, stability_ratio_, ie) = stability;
					value<Vec3>(m_, stability_color_, ie) = (stability <= 0.5)
																? Vec3(0, stability, (0.5 - stability))
																: Vec3(stability - 0.5, (1 - stability), 0);
				}
				// recompute the cost of the edge
				Vec4 opt;
				Scalar cost_opt = edge_optimal(ie, opt);
				value<EdgeInfo>(m_, edge_queue_it_, ie) = {true, queue_.emplace(cost_opt, ie)};

				return true;
			});
			++count;
		}
	}

	bool slab_normal(Vec4& p1, Vec4& p2, Vec4& p3, Vec4& n1, Vec4& n2)
	{

		Vec3 c1 = p1.head<3>();
		Vec3 c2 = p2.head<3>();
		Vec3 c3 = p3.head<3>();
		double r1 = p1[3];
		double r2 = p2[3];
		double r3 = p3[3];

		Vec3 c12 = c2 - c1;
		Vec3 c13 = c3 - c1;
		Vec3 c23 = c3 - c2;
		double dc12 = c12.norm();
		double dc13 = c13.norm();
		double dc23 = c23.norm();

		double dr12 = fabs(r2 - r1);
		double dr13 = fabs(r3 - r1);
		double dr23 = fabs(r3 - r2);

		if ((dc12 < 1e-8) || (dc13 < 1e-8) || (dc23 < 1e-8))
			return false;

		Vec3 n = c12.cross(c13).normalized();

		if ((dr12 < 1e-8) && (dr13 < 1e-8) && (dr23 < 1e-8))
		{
			n1 = Vec4(n[0], n[1], n[2], 1);
			n2 = Vec4(-n[0], -n[1], -n[2], 1);
			return true;
		}
		else
		{
			// two points on the tangent plane
			Vec3 contact_point1, contact_point2;

			if (dr12 < 1e-8)
			{
				contact_point1 = (r3 * c1 - r1 * c3) / (r3 - r1);
				contact_point2 = (r3 * c2 - r2 * c3) / (r3 - r2);
			}
			else if (dr13 < 1e-8)
			{
				contact_point1 = (r2 * c1 - r1 * c2) / (r2 - r1);
				contact_point2 = (r3 * c2 - r2 * c3) / (r3 - r2);
			}
			else if (dr23 < 1e-8)
			{
				contact_point1 = (r3 * c1 - r1 * c3) / (r3 - r1);
				contact_point2 = (r1 * c2 - r2 * c1) / (r1 - r2);
			}
			else
			{
				contact_point1 = (r3 * c1 - r1 * c3) / (r3 - r1);
				contact_point2 = (r3 * c2 - r2 * c3) / (r3 - r2);
			}

			Vec3 inter_point;
			double dist;
			DistanceToLine(c1, contact_point1, contact_point2, dist, inter_point);
			double sangle = r1 / dist;
			if (fabs(sangle) > 1.)
				return false;
			double cangle = sqrt(1. - sangle * sangle);
			Vec3 c1inter = (c1 - inter_point ).normalized();
			Vec3 n1_vec3 = n * cangle - c1inter * sangle;
			Vec3 n2_vec3 = -n * cangle - c1inter * sangle;
			n1_vec3.normalize();
			n2_vec3.normalize();
			Vec3 v11 = c1 + r1 * n1_vec3;
			Vec3 v12 = c2 + r2 * n1_vec3;
			Vec3 v13 = c3 + r3 * n1_vec3;

			Vec3 nv1 = (v12- v11).cross(v13-v11).normalized();

			Vec3 v21 = c1 + r1 * n2_vec3;
			Vec3 v22 = c2 + r2 * n2_vec3;
			Vec3 v23 = c3 + r3 * n2_vec3;

			Vec3 nv2 = (v22 - v21).cross(v23 - v21).normalized() * -1;

			n1 = Vec4(nv1[0], nv1[1], nv1[2], 1);
			n2 = Vec4(nv2[0], nv2[1], nv2[2], 1);
		}
		return true;
	}

	bool DistanceToLine(const Vec3& p, const Vec3& v0, const Vec3& v1, double& dist, Vec3& fp)
	{
		Vec3 v0v1(v1 - v0), pv0(v0 - p), pv1(v1 - p);
		double area = fabs(v0v1.cross(pv0).norm());
		if (v0v1.norm() > 1e-12)
		{
			dist = area / v0v1.norm();
			double t = (pv0.dot(pv0) - pv0.dot(pv1)) / (pv0.dot(pv0) + pv1.dot(pv1) - 2 * pv0.dot(pv1));
			fp = (1 - t) * v0 + t * v1;
			return true;
		}
		else
			return false;
	}

	Scalar edge_optimal(Edge e, Vec4& opt)
	{
		Scalar cost = 0;
		std::vector<Vertex> iv = incident_vertices(m_, e);
		auto p1 = value<Vec4>(m_, sphere_info_, iv[0]);
		auto p2 = value<Vec4>(m_, sphere_info_, iv[1]);
		auto eq = value<Slab_Quadric>(m_, edge_slab_quadric_, e);
		auto q1 = value<Slab_Quadric>(m_, vertex_slab_quadric_, iv[0]);
		auto q2 = value<Slab_Quadric>(m_, vertex_slab_quadric_, iv[1]);
		double st = value<double>(m_, stability_ratio_, e);
		eq._A = q1._A + q2._A; 
		eq._b = q1._b + q2._b;
		eq._c = q1._c + q2._c;
		Vec4 p;
		if (eq.optimized(opt) || is_boundary(m_, e))
		{
			eq._A = q1._A + q2._A + q1._add_A + q2._add_A;
			eq._b = q1._b + q2._b + q1._add_b + q2._add_b;
			eq._c = q1._c + q2._c + q1._add_c + q2._add_c;
			if (eq.optimized(opt))
			{
				if (opt[3] < 0)
					opt = Scalar(0.5) * (p1 + p2);
				}
			else {
				opt = Scalar(0.5) * (p1 + p2);
			}
		}
		else
		{
			Vec4 spheres[3];
			Scalar cost_collapse[3];
			Vec4 p3 = (p1 + p2) * Scalar(0.5);
			int min_index = 0;
			spheres[0] = p1;
			spheres[1] = p2;
			spheres[2] = p3;
			
			cost_collapse[0] = eq.eval(p1);
			cost_collapse[1] = eq.eval(p2);
			cost_collapse[2] = eq.eval(p3);
			if (cost_collapse[0] >= cost_collapse[1])
				min_index = 1;
			min_index = cost_collapse[min_index] > cost_collapse[2] ? 2 : min_index;
			cost = cost_collapse[min_index];
			opt = spheres[min_index];
			return (cost + k_) * st * st;
		}
		cost = eq.eval(opt);
		
		return (cost + k_) * st * st;
	}

	EdgeQueue queue_;
	float k_;
	NONMANIFOLD& m_;
	std::shared_ptr<Attribute<Vec4>> sphere_info_;
	std::shared_ptr<Attribute<Vec3>> stability_color_;
	std::shared_ptr<Attribute<double>> stability_ratio_;
	std::shared_ptr<Attribute<Vec3>> position_;

	std::shared_ptr<Attribute<EdgeInfo>> edge_queue_it_;
	std::shared_ptr<Attribute<Vec4>> sphere_opt_;
	std::shared_ptr<Attribute<Slab_Quadric>> edge_slab_quadric_;
	std::shared_ptr<Attribute<Slab_Quadric>> vertex_slab_quadric_;
	std::shared_ptr<Attribute<std::pair<Vec4, Vec4>>> slab_normals_;


};

} // namespace modeling

} // namespace cgogn

#endif // CGOGN_MODELING_DECIMATION_EDGE_QUEUE_SQEM_H_
