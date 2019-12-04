/*******************************************************************************
* CGoGN: Combinatorial and Geometric modeling with Generic N-dimensional Maps  *
* Copyright (C) 2015, IGG Group, ICube, University of Strasbourg, France       *
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

#ifndef CGOGN_CORE_TYPES_CMAP_CMAP1_H_
#define CGOGN_CORE_TYPES_CMAP_CMAP1_H_

#include <cgogn/core/cgogn_core_export.h>

#include <cgogn/core/types/cmap/cmap0.h>

namespace cgogn
{

struct CGOGN_CORE_EXPORT CMap1 : public CMap0
{
	
	using Self = CMap1;
	using Inherit = CMap0;
	
	using AttributeContainer = Inherit::AttributeContainer;

	template <typename T>
	using Attribute = Inherit::Attribute<T>;
	using AttributeGen = Inherit::AttributeGen;
	using MarkAttribute = Inherit::MarkAttribute;
	
	std::shared_ptr<Attribute<Dart>> phi1_;
	std::shared_ptr<Attribute<Dart>> phi_1_;

	using Vertex = Cell<DART>;
	using Edge = Cell<DART>;
	using Face = Cell<PHI1>;
	using CC = Face;

	using Cells = std::tuple<Vertex, Edge, Face>;

	CMap1() : CMap0()
	{
		phi1_ = base_map_->add_or_get_relation("phi1");
		phi_1_ = base_map_->add_relation("phi_1");
	}

	CMap1(std::shared_ptr<CMapBase> m) : CMap0(m)
	{
		phi1_ = base_map_->add_or_get_relation("phi1");
		phi_1_ = base_map_->add_relation("phi_1");
	}
};

} // namespace cgogn

#endif // CGOGN_CORE_TYPES_CMAP_CMAP1_H_
