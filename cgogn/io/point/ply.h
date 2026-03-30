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

#ifndef CGOGN_IO_POINT_PLY_H_
#define CGOGN_IO_POINT_PLY_H_

#include <cgogn/io/point/export_options.h>
#include <cgogn/io/point/point_import.h>
#include <cgogn/io/utils.h>

#include <cgogn/core/functions/attributes.h>
#include <cgogn/core/functions/mesh_info.h>

#include <thirdparty/happly/happly.h>

#include <algorithm>
#include <cctype>
#include <cmath>

namespace cgogn
{

namespace io
{

namespace point_internal
{

inline std::string sanitize_ply_property_name(const std::string& name)
{
	std::string result = name;
	for (char& c : result)
	{
		if (!(std::isalnum(static_cast<unsigned char>(c)) || c == '_'))
			c = '_';
	}
	return result;
}

inline unsigned char ply_color_channel(float64 value)
{
	return static_cast<unsigned char>(std::lround(std::clamp(value, 0.0, 1.0) * 255.0));
}

template <typename T>
void add_scalar_ply_property(happly::Element& element, const std::string& name, const std::vector<T>& values)
{
	element.addProperty<T>(name, values);
}

template <typename MESH, typename CELL>
bool add_selected_ply_color_attribute(
	happly::Element& element, MESH& m, const std::shared_ptr<typename mesh_traits<MESH>::AttributeGen>& attribute_gen,
	std::integral_constant<uint32, 3>)
{
	using Vec3 = geometry::Vec3;
	using Attribute = typename mesh_traits<MESH>::template Attribute<Vec3>;
	auto attribute = std::dynamic_pointer_cast<Attribute>(attribute_gen);
	if (!attribute)
		return false;

	std::vector<unsigned char> red;
	std::vector<unsigned char> green;
	std::vector<unsigned char> blue;
	red.reserve(nb_cells<CELL>(m));
	green.reserve(nb_cells<CELL>(m));
	blue.reserve(nb_cells<CELL>(m));
	foreach_cell(m, [&](CELL c) -> bool {
		const Vec3& color = value<Vec3>(m, attribute.get(), c);
		red.push_back(ply_color_channel(color[0]));
		green.push_back(ply_color_channel(color[1]));
		blue.push_back(ply_color_channel(color[2]));
		return true;
	});

	element.addProperty<unsigned char>("red", red);
	element.addProperty<unsigned char>("green", green);
	element.addProperty<unsigned char>("blue", blue);
	return true;
}

template <typename MESH, typename CELL>
bool add_selected_ply_color_attribute(
	happly::Element& element, MESH& m, const std::shared_ptr<typename mesh_traits<MESH>::AttributeGen>& attribute_gen,
	std::integral_constant<uint32, 4>)
{
	using Vec4 = geometry::Vec4;
	using Attribute = typename mesh_traits<MESH>::template Attribute<Vec4>;
	auto attribute = std::dynamic_pointer_cast<Attribute>(attribute_gen);
	if (!attribute)
		return false;

	std::vector<unsigned char> red;
	std::vector<unsigned char> green;
	std::vector<unsigned char> blue;
	red.reserve(nb_cells<CELL>(m));
	green.reserve(nb_cells<CELL>(m));
	blue.reserve(nb_cells<CELL>(m));
	foreach_cell(m, [&](CELL c) -> bool {
		const Vec4& color = value<Vec4>(m, attribute.get(), c);
		red.push_back(ply_color_channel(color[0]));
		green.push_back(ply_color_channel(color[1]));
		blue.push_back(ply_color_channel(color[2]));
		return true;
	});

	element.addProperty<unsigned char>("red", red);
	element.addProperty<unsigned char>("green", green);
	element.addProperty<unsigned char>("blue", blue);
	return true;
}

template <typename MESH, typename CELL>
bool add_selected_ply_color_attribute(
	happly::Element&, MESH&, const std::shared_ptr<typename mesh_traits<MESH>::AttributeGen>&, std::integral_constant<uint32, 2>)
{
	return false;
}

template <typename MESH, typename CELL>
bool add_selected_ply_color_attribute(happly::Element&, MESH&, const std::shared_ptr<typename mesh_traits<MESH>::AttributeGen>&)
{
	return false;
}

template <typename VEC>
void add_vector_ply_property_components(happly::Element& element, const std::string& name,
										 const std::vector<VEC>& values, const char* const* suffixes, uint32 dim)
{
	std::vector<double> components[4];
	for (uint32 i = 0; i < dim; ++i)
		components[i].reserve(values.size());
	for (const VEC& v : values)
	{
		for (uint32 i = 0; i < dim; ++i)
			components[i].push_back(static_cast<double>(v[i]));
	}
	for (uint32 i = 0; i < dim; ++i)
		element.addProperty<double>(name + suffixes[i], components[i]);
}

template <typename MESH, typename CELL, typename T>
bool try_add_selected_ply_attribute(happly::Element& element, MESH& m,
									 const std::shared_ptr<typename mesh_traits<MESH>::AttributeGen>& attribute_gen)
{
	using Attribute = typename mesh_traits<MESH>::template Attribute<T>;
	auto attribute = std::dynamic_pointer_cast<Attribute>(attribute_gen);
	if (!attribute)
		return false;

	std::vector<T> values;
	values.reserve(nb_cells<CELL>(m));
	foreach_cell(m, [&](CELL c) -> bool {
		values.push_back(value<T>(m, attribute.get(), c));
		return true;
	});

	add_scalar_ply_property(element, sanitize_ply_property_name(attribute->name()), values);
	return true;
}

template <typename MESH, typename CELL>
bool try_add_selected_ply_attribute(happly::Element& element, MESH& m,
									 const std::shared_ptr<typename mesh_traits<MESH>::AttributeGen>& attribute_gen,
									 std::integral_constant<uint32, 2>)
{
	using Vec2 = geometry::Vec2;
	using Attribute = typename mesh_traits<MESH>::template Attribute<Vec2>;
	auto attribute = std::dynamic_pointer_cast<Attribute>(attribute_gen);
	if (!attribute)
		return false;

	std::vector<Vec2> values;
	values.reserve(nb_cells<CELL>(m));
	foreach_cell(m, [&](CELL c) -> bool {
		values.push_back(value<Vec2>(m, attribute.get(), c));
		return true;
	});

	static const char* suffixes[2] = {"_x", "_y"};
	add_vector_ply_property_components(element, sanitize_ply_property_name(attribute->name()), values, suffixes, 2);
	return true;
}

template <typename MESH, typename CELL>
bool try_add_selected_ply_attribute(happly::Element& element, MESH& m,
									 const std::shared_ptr<typename mesh_traits<MESH>::AttributeGen>& attribute_gen,
									 std::integral_constant<uint32, 3>)
{
	using Vec3 = geometry::Vec3;
	using Attribute = typename mesh_traits<MESH>::template Attribute<Vec3>;
	auto attribute = std::dynamic_pointer_cast<Attribute>(attribute_gen);
	if (!attribute)
		return false;

	std::vector<Vec3> values;
	values.reserve(nb_cells<CELL>(m));
	foreach_cell(m, [&](CELL c) -> bool {
		values.push_back(value<Vec3>(m, attribute.get(), c));
		return true;
	});

	static const char* suffixes[3] = {"_x", "_y", "_z"};
	add_vector_ply_property_components(element, sanitize_ply_property_name(attribute->name()), values, suffixes, 3);
	return true;
}

template <typename MESH, typename CELL>
bool try_add_selected_ply_attribute(happly::Element& element, MESH& m,
									 const std::shared_ptr<typename mesh_traits<MESH>::AttributeGen>& attribute_gen,
									 std::integral_constant<uint32, 4>)
{
	using Vec4 = geometry::Vec4;
	using Attribute = typename mesh_traits<MESH>::template Attribute<Vec4>;
	auto attribute = std::dynamic_pointer_cast<Attribute>(attribute_gen);
	if (!attribute)
		return false;

	std::vector<Vec4> values;
	values.reserve(nb_cells<CELL>(m));
	foreach_cell(m, [&](CELL c) -> bool {
		values.push_back(value<Vec4>(m, attribute.get(), c));
		return true;
	});

	static const char* suffixes[4] = {"_x", "_y", "_z", "_w"};
	add_vector_ply_property_components(element, sanitize_ply_property_name(attribute->name()), values, suffixes, 4);
	return true;
}

template <typename MESH, typename CELL>
void add_selected_ply_attributes(happly::Element& element, MESH& m,
								 const std::vector<std::shared_ptr<typename mesh_traits<MESH>::AttributeGen>>& attributes,
								 const std::shared_ptr<typename mesh_traits<MESH>::AttributeGen>& color_attribute = nullptr)
{
	for (const auto& attribute_gen : attributes)
	{
		if (!attribute_gen)
			continue;
		if (color_attribute && attribute_gen.get() == color_attribute.get())
			continue;

		bool exported = false;
		exported = try_add_selected_ply_attribute<MESH, CELL, int32>(element, m, attribute_gen);
		if (!exported)
			exported = try_add_selected_ply_attribute<MESH, CELL, uint32>(element, m, attribute_gen);
		if (!exported)
			exported = try_add_selected_ply_attribute<MESH, CELL, float32>(element, m, attribute_gen);
		if (!exported)
			exported = try_add_selected_ply_attribute<MESH, CELL, float64>(element, m, attribute_gen);
		if (!exported)
			exported =
				try_add_selected_ply_attribute<MESH, CELL>(element, m, attribute_gen, std::integral_constant<uint32, 2>{});
		if (!exported)
			exported =
				try_add_selected_ply_attribute<MESH, CELL>(element, m, attribute_gen, std::integral_constant<uint32, 3>{});
		if (!exported)
			exported =
				try_add_selected_ply_attribute<MESH, CELL>(element, m, attribute_gen, std::integral_constant<uint32, 4>{});
	}
}

} // namespace point_internal

template <typename MESH>
typename std::enable_if<mesh_traits<MESH>::dimension == 0, bool>::type import_PLY(MESH& m, const std::string& filename)
{
	static_assert(mesh_traits<MESH>::dimension == 0, "MESH dimension should be 0");

	using Vertex = typename MESH::Vertex;

	Scoped_C_Locale loc;

	PointImportData point_data;

	happly::PLYData plyData(filename);
	std::vector<std::array<double, 3>> position = plyData.getVertexPositions();

	const uint32 nb_vertices = position.size();

	point_data.reserve(nb_vertices);

	for (uint32 i = 0u; i < nb_vertices; ++i)
	{
		const std::array<double, 3>& p = position[i];
		point_data.vertex_position_.push_back({p[0], p[1], p[2]});
	}

	import_point_data(m, point_data);

	return true;
}

template <typename MESH>
typename std::enable_if<mesh_traits<MESH>::dimension == 0, void>::type export_PLY(
	MESH& m, const typename mesh_traits<MESH>::template Attribute<geometry::Vec3>* vertex_position,
	const std::string& filename, const PointExportAttributeSelection<MESH>* export_attributes = nullptr)
{
	static_assert(mesh_traits<MESH>::dimension == 0, "MESH dimension should be 0");

	using Vertex = typename mesh_traits<MESH>::Vertex;
	using Vec3 = geometry::Vec3;

	Scoped_C_Locale loc;

	std::vector<std::array<double, 3>> position;
	uint32 nb_vertices = nb_cells<Vertex>(m);
	position.reserve(nb_vertices);

	foreach_cell(m, [&](Vertex v) -> bool {
		const Vec3& p = value<geometry::Vec3>(m, vertex_position, v);
		position.push_back({p.x(), p.y(), p.z()});
		return true;
	});

	happly::PLYData plyOut;
	plyOut.addVertexPositions(position);
	if (export_attributes)
	{
		std::shared_ptr<typename mesh_traits<MESH>::AttributeGen> exported_color_attribute = nullptr;
		if (export_attributes->vertex_color_attribute)
		{
			bool exported_color = false;
			exported_color = point_internal::add_selected_ply_color_attribute<MESH, Vertex>(
				plyOut.getElement("vertex"), m, export_attributes->vertex_color_attribute, std::integral_constant<uint32, 3>{});
			if (!exported_color)
				exported_color = point_internal::add_selected_ply_color_attribute<MESH, Vertex>(
					plyOut.getElement("vertex"), m, export_attributes->vertex_color_attribute, std::integral_constant<uint32, 4>{});
			if (exported_color)
				exported_color_attribute = export_attributes->vertex_color_attribute;
		}
		point_internal::add_selected_ply_attributes<MESH, Vertex>(
			plyOut.getElement("vertex"), m, export_attributes->vertex_attributes, exported_color_attribute);
	}
	plyOut.write(filename);
}

} // namespace io

} // namespace cgogn

#endif // CGOGN_IO_POINT_PLY_H_
