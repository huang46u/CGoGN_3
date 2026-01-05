/*******************************************************************************
 * CGoGN: Combinatorial and Geometric modeling with Generic N-dimensional Maps  *
 * Copyright (C), IGG Group, ICube, University of Strasbourg, France            *
 *******************************************************************************/
#pragma once

#include <cgogn/core/types/mesh_traits.h>
#include <cgogn/core/functions/attributes.h>
#include <cgogn/core/functions/mesh_ops.h>

#include <cgogn/geometry/algos/area.h>
#include <cgogn/geometry/algos/centroid.h>
#include <cgogn/geometry/algos/normal.h>
#include <cgogn/geometry/types/fast_winding_number.h>

#include <cgogn/io/surface/surface_import.h>
#include <cgogn/io/surface/obj.h>

#include <cgogn/ui/app.h>
#include <cgogn/ui/view.h>
#include <cgogn/ui/modules/mesh_provider/mesh_provider.h>
#include <cgogn/ui/modules/surface_render/surface_render.h>

#include <cgogn/rendering/shaders/shader_point_sprite.h>

#include <libacc/bvh_tree.h>

#include <random>
#include <iostream>
#include <iomanip>

using namespace cgogn::numerics;

using Mesh = cgogn::CMap2;
using PointCloud = cgogn::CMap0;

using Vertex = typename cgogn::mesh_traits<Mesh>::Vertex;
using Face = typename cgogn::mesh_traits<Mesh>::Face;

template <typename MESH>
using Attribute = typename cgogn::mesh_traits<MESH>::template Attribute<cgogn::geometry::Vec3>;

// Build BVH for triangle mesh
acc::BVHTree<std::size_t, cgogn::geometry::Vec3> build_mesh_bvh(
	const Mesh& mesh,
	const cgogn::mesh_traits<Mesh>::template Attribute<cgogn::geometry::Vec3>* vertex_position,
	std::vector<Face>& out_faces)
{
	out_faces.clear();
	std::vector<acc::BVHTree<std::size_t, cgogn::geometry::Vec3>::Primitive> primitives;

	cgogn::foreach_cell(mesh, [&](Face f) -> bool {
		auto vertices = cgogn::incident_vertices(mesh, f);
		if (vertices.size() >= 3)
		{
			cgogn::geometry::Vec3 p0 = cgogn::value<cgogn::geometry::Vec3>(mesh, vertex_position, vertices[0]);
			cgogn::geometry::Vec3 p1 = cgogn::value<cgogn::geometry::Vec3>(mesh, vertex_position, vertices[1]);
			cgogn::geometry::Vec3 p2 = cgogn::value<cgogn::geometry::Vec3>(mesh, vertex_position, vertices[2]);

			acc::BVHTree<std::size_t, cgogn::geometry::Vec3>::AABB box;
			box.min = p0.cwiseMin(p1).cwiseMin(p2);
			box.max = p0.cwiseMax(p1).cwiseMax(p2);

			primitives.push_back({out_faces.size(), box});
			out_faces.push_back(f);
		}
		return true;
	});

	return acc::BVHTree<std::size_t, cgogn::geometry::Vec3>(primitives);
}

// Build BVH for point cloud
acc::BVHTree<std::size_t, cgogn::geometry::Vec3> build_pointcloud_bvh(
	const PointCloud& cloud,
	const cgogn::mesh_traits<PointCloud>::template Attribute<cgogn::geometry::Vec3>* vertex_position,
	std::vector<Vertex>& out_vertices)
{
	out_vertices.clear();
	std::vector<acc::BVHTree<std::size_t, cgogn::geometry::Vec3>::Primitive> primitives;

	cgogn::foreach_cell(cloud, [&](Vertex v) -> bool {
		cgogn::geometry::Vec3 p = cgogn::value<cgogn::geometry::Vec3>(cloud, vertex_position, v);

		acc::BVHTree<std::size_t, cgogn::geometry::Vec3>::AABB box;
		box.min = p - cgogn::geometry::Vec3(1e-6, 1e-6, 1e-6);
		box.max = p + cgogn::geometry::Vec3(1e-6, 1e-6, 1e-6);

		primitives.push_back({out_vertices.size(), box});
		out_vertices.push_back(v);
		return true;
	});

	return acc::BVHTree<std::size_t, cgogn::geometry::Vec3>(primitives);
}

// Compute bounding box
struct BBox
{
	cgogn::geometry::Vec3 min, max;

	BBox() : min(std::numeric_limits<cgogn::Scalar>::max(), 
				 std::numeric_limits<cgogn::Scalar>::max(), 
				 std::numeric_limits<cgogn::Scalar>::max()),
			 max(std::numeric_limits<cgogn::Scalar>::lowest(),
				 std::numeric_limits<cgogn::Scalar>::lowest(),
				 std::numeric_limits<cgogn::Scalar>::lowest()) {}

	void expand(const cgogn::geometry::Vec3& p)
	{
		min = min.cwiseMin(p);
		max = max.cwiseMax(p);
	}

	cgogn::geometry::Vec3 center() const { return 0.5 * (min + max); }
	cgogn::geometry::Vec3 extent() const { return max - min; }
};

template <typename MESH>
BBox compute_bbox(const MESH& mesh,
				  const typename cgogn::mesh_traits<MESH>::template Attribute<cgogn::geometry::Vec3>* vertex_position)
{
	BBox bbox;
	cgogn::foreach_cell(mesh, [&](typename cgogn::mesh_traits<MESH>::Vertex v) -> bool {
		bbox.expand(cgogn::value<cgogn::geometry::Vec3>(mesh, vertex_position, v));
		return true;
	});
	return bbox;
}

// Test result structure
struct TestResult
{
	std::string name;
	double mean_error;
	double max_error;
	double rmse;
	int num_samples;
};

// Main application
int main(int argc, char** argv)
{
	if (argc < 3)
	{
		std::cerr << "Usage: " << argv[0] << " <mesh.obj> <pointcloud.obj>" << std::endl;
		return 1;
	}

	std::string mesh_file = argv[1];
	std::string cloud_file = argv[2];

	// Create UI application
	auto app = cgogn::ui::App::init_app("Fast Winding Number Test");
	auto& view = app->add_view();

	// Load mesh
	Mesh mesh;
	cgogn::io::import_surface<Mesh>(mesh_file, mesh);

	auto mesh_vertex_position = cgogn::get_attribute<cgogn::geometry::Vec3, Vertex>(mesh, "position");
	if (!mesh_vertex_position)
	{
		std::cerr << "Mesh has no position attribute!" << std::endl;
		return 1;
	}

	// Compute mesh attributes
	auto mesh_face_normal = cgogn::add_attribute<cgogn::geometry::Vec3, Face>(mesh, "normal");
	auto mesh_face_area = cgogn::add_attribute<cgogn::Scalar, Face>(mesh, "area");
	auto mesh_face_centroid = cgogn::add_attribute<cgogn::geometry::Vec3, Face>(mesh, "centroid");

	cgogn::geometry::compute_normal<Face>(mesh, mesh_vertex_position.get(), mesh_face_normal.get());
	cgogn::geometry::compute_area<Face>(mesh, mesh_vertex_position.get(), mesh_face_area.get());
	cgogn::geometry::compute_centroid<cgogn::geometry::Vec3, Face>(mesh, mesh_vertex_position.get(), mesh_face_centroid.get());

	// Build mesh BVH
	std::vector<Face> mesh_faces;
	auto mesh_bvh = build_mesh_bvh(mesh, mesh_vertex_position.get(), mesh_faces);
	std::cout << "Mesh BVH built with " << mesh_faces.size() << " faces" << std::endl;

	// Load point cloud
	PointCloud cloud;
	cgogn::io::import_surface<PointCloud>(cloud_file, cloud);

	auto cloud_vertex_position = cgogn::get_attribute<cgogn::geometry::Vec3, Vertex>(cloud, "position");
	if (!cloud_vertex_position)
	{
		std::cerr << "Point cloud has no position attribute!" << std::endl;
		return 1;
	}

	// Compute point cloud attributes (estimate normals and areas)
	auto cloud_vertex_normal = cgogn::add_attribute<cgogn::geometry::Vec3, Vertex>(cloud, "normal");
	auto cloud_vertex_area = cgogn::add_attribute<cgogn::Scalar, Vertex>(cloud, "area");

	// Simple normal estimation: average of neighboring normals (placeholder)
	cgogn::Scalar avg_area = 0.01; // Placeholder: estimate from point density
	cgogn::foreach_cell(cloud, [&](Vertex v) -> bool {
		// For now, use simple estimation
		cgogn::value<cgogn::geometry::Vec3>(cloud, cloud_vertex_normal, v) = cgogn::geometry::Vec3(0, 0, 1);
		cgogn::value<cgogn::Scalar>(cloud, cloud_vertex_area, v) = avg_area;
		return true;
	});

	// Build point cloud BVH
	std::vector<Vertex> cloud_vertices;
	auto cloud_bvh = build_pointcloud_bvh(cloud, cloud_vertex_position.get(), cloud_vertices);
	std::cout << "Point cloud BVH built with " << cloud_vertices.size() << " points" << std::endl;

	// Create FWN instances
	cgogn::Scalar beta = 2.0;

	// Mesh FWN
	cgogn::geometry::Fast_Winding_Number_Triangle_Traits<Mesh> mesh_traits(
		mesh, mesh_vertex_position.get(), mesh_face_normal.get(), 
		mesh_face_area.get(), mesh_face_centroid.get(), mesh_faces);

	cgogn::geometry::Fast_Winding_Number<cgogn::geometry::Fast_Winding_Number_Triangle_Traits<Mesh>, 1> mesh_fwn1(mesh_traits, mesh_bvh, beta);
	cgogn::geometry::Fast_Winding_Number<cgogn::geometry::Fast_Winding_Number_Triangle_Traits<Mesh>, 2> mesh_fwn2(mesh_traits, mesh_bvh, beta);
	cgogn::geometry::Fast_Winding_Number<cgogn::geometry::Fast_Winding_Number_Triangle_Traits<Mesh>, 3> mesh_fwn3(mesh_traits, mesh_bvh, beta);

	// Point cloud FWN
	cgogn::geometry::Fast_Winding_Number__PointCloud_Traits<PointCloud> cloud_traits(
		cloud, cloud_vertex_position.get(), cloud_vertex_normal.get(), 
		cloud_vertex_area.get(), cloud_vertices);

	cgogn::geometry::Fast_Winding_Number<cgogn::geometry::Fast_Winding_Number__PointCloud_Traits<PointCloud>, 1> cloud_fwn1(cloud_traits, cloud_bvh, beta);
	cgogn::geometry::Fast_Winding_Number<cgogn::geometry::Fast_Winding_Number__PointCloud_Traits<PointCloud>, 2> cloud_fwn2(cloud_traits, cloud_bvh, beta);

	// Compute bounding box
	BBox mesh_bbox = compute_bbox(mesh, mesh_vertex_position.get());
	BBox cloud_bbox = compute_bbox(cloud, cloud_vertex_position.get());

	// Expand bbox slightly
	cgogn::geometry::Vec3 padding = 0.1 * mesh_bbox.extent();
	mesh_bbox.min -= padding;
	mesh_bbox.max += padding;

	// Sample points in bbox
	const int num_samples = 10000;
	std::vector<cgogn::geometry::Vec3> sample_points;
	std::vector<cgogn::geometry::Vec3> sample_colors;

	std::random_device rd;
	std::mt19937 gen(rd());
	std::uniform_real_distribution<cgogn::Scalar> dist_x(mesh_bbox.min.x(), mesh_bbox.max.x());
	std::uniform_real_distribution<cgogn::Scalar> dist_y(mesh_bbox.min.y(), mesh_bbox.max.y());
	std::uniform_real_distribution<cgogn::Scalar> dist_z(mesh_bbox.min.z(), mesh_bbox.max.z());

	for (int i = 0; i < num_samples; ++i)
	{
		sample_points.push_back(cgogn::geometry::Vec3(dist_x(gen), dist_y(gen), dist_z(gen)));
	}

	std::cout << "\n=== Testing Mesh FWN ===" << std::endl;
	std::cout << "Sampling " << num_samples << " points..." << std::endl;

	// Test mesh FWN
	std::vector<TestResult> mesh_results;
	
	auto test_mesh_order = [&](int order, auto& fwn, const std::string& name) {
		TestResult result;
		result.name = name;
		result.num_samples = num_samples;
		
		double sum_error = 0, sum_sq_error = 0, max_err = 0;
		
		for (const auto& q : sample_points)
		{
			cgogn::Scalar exact = mesh_fwn3.exact_winding_number(q);
			cgogn::Scalar approx = fwn.evaluate_fast_winding_number(q);
			
			double error = std::abs(exact - approx);
			sum_error += error;
			sum_sq_error += error * error;
			max_err = std::max(max_err, error);
		}
		
		result.mean_error = sum_error / num_samples;
		result.max_error = max_err;
		result.rmse = std::sqrt(sum_sq_error / num_samples);
		
		mesh_results.push_back(result);
	};

	test_mesh_order(1, mesh_fwn1, "Mesh ORDER=1");
	test_mesh_order(2, mesh_fwn2, "Mesh ORDER=2");
	test_mesh_order(3, mesh_fwn3, "Mesh ORDER=3 (Fast)");

	// Print mesh results
	std::cout << "\n--- Mesh Results ---" << std::endl;
	std::cout << std::setw(20) << "Method" 
			  << std::setw(15) << "Mean Error" 
			  << std::setw(15) << "Max Error"
			  << std::setw(15) << "RMSE" << std::endl;
	std::cout << std::string(65, '-') << std::endl;
	
	for (const auto& r : mesh_results)
	{
		std::cout << std::setw(20) << r.name
				  << std::setw(15) << std::scientific << std::setprecision(4) << r.mean_error
				  << std::setw(15) << r.max_error
				  << std::setw(15) << r.rmse << std::endl;
	}

	// Test point cloud FWN
	std::cout << "\n=== Testing Point Cloud FWN ===" << std::endl;
	std::vector<TestResult> cloud_results;
	
	auto test_cloud_order = [&](int order, auto& fwn, const std::string& name) {
		TestResult result;
		result.name = name;
		result.num_samples = num_samples;
		
		double sum_error = 0, sum_sq_error = 0, max_err = 0;
		
		for (const auto& q : sample_points)
		{
			cgogn::Scalar exact = cloud_fwn2.exact_winding_number(q);
			cgogn::Scalar approx = fwn.evaluate_fast_winding_number(q);
			
			double error = std::abs(exact - approx);
			sum_error += error;
			sum_sq_error += error * error;
			max_err = std::max(max_err, error);
		}
		
		result.mean_error = sum_error / num_samples;
		result.max_error = max_err;
		result.rmse = std::sqrt(sum_sq_error / num_samples);
		
		cloud_results.push_back(result);
	};

	test_cloud_order(1, cloud_fwn1, "Cloud ORDER=1");
	test_cloud_order(2, cloud_fwn2, "Cloud ORDER=2");

	// Print cloud results
	std::cout << "\n--- Point Cloud Results ---" << std::endl;
	std::cout << std::setw(20) << "Method" 
			  << std::setw(15) << "Mean Error" 
			  << std::setw(15) << "Max Error"
			  << std::setw(15) << "RMSE" << std::endl;
	std::cout << std::string(65, '-') << std::endl;
	
	for (const auto& r : cloud_results)
	{
		std::cout << std::setw(20) << r.name
				  << std::setw(15) << std::scientific << std::setprecision(4) << r.mean_error
				  << std::setw(15) << r.max_error
				  << std::setw(15) << r.rmse << std::endl;
	}

	// Colorize sample points based on winding number (mesh ORDER=3)
	const cgogn::Scalar threshold = 0.5; // Inside/outside threshold
	sample_colors.resize(num_samples);
	
	for (int i = 0; i < num_samples; ++i)
	{
		cgogn::Scalar w = mesh_fwn3.evaluate_fast_winding_number(sample_points[i]);
		
		if (w > threshold)
		{
			// Inside: red
			sample_colors[i] = cgogn::geometry::Vec3(1.0, 0.0, 0.0);
		}
		else
		{
			// Outside: green
			sample_colors[i] = cgogn::geometry::Vec3(0.0, 1.0, 0.0);
		}
	}

	// Visualization setup
	std::cout << "\n=== Rendering ===" << std::endl;
	std::cout << "Red points: inside (w > 0.5)" << std::endl;
	std::cout << "Green points: outside (w <= 0.5)" << std::endl;

	// Create point cloud for visualization
	PointCloud vis_cloud;
	auto vis_position = cgogn::add_attribute<cgogn::geometry::Vec3, Vertex>(vis_cloud, "position");
	auto vis_color = cgogn::add_attribute<cgogn::geometry::Vec3, Vertex>(vis_cloud, "color");

	for (int i = 0; i < num_samples; ++i)
	{
		Vertex v = cgogn::add_vertex(vis_cloud);
		cgogn::value<cgogn::geometry::Vec3>(vis_cloud, vis_position, v) = sample_points[i];
		cgogn::value<cgogn::geometry::Vec3>(vis_cloud, vis_color, v) = sample_colors[i];
	}

	// Setup rendering (simplified - you may need to adapt to your UI framework)
	auto* mesh_provider = static_cast<cgogn::ui::MeshProvider<Mesh>*>(
		view.add_module(new cgogn::ui::MeshProvider<Mesh>(mesh, "input_mesh")));
	mesh_provider->set_vertex_position(mesh_vertex_position);

	auto* vis_provider = static_cast<cgogn::ui::MeshProvider<PointCloud>*>(
		view.add_module(new cgogn::ui::MeshProvider<PointCloud>(vis_cloud, "sample_points")));
	vis_provider->set_vertex_position(vis_position);

	// Run app
	return app->launch();
}