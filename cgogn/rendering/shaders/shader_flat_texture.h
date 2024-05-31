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

#ifndef CGOGN_RENDERING_SHADERS_FLAT_TEXTURE_H_
#define CGOGN_RENDERING_SHADERS_FLAT_TEXTURE_H_

#include <cgogn/rendering/cgogn_rendering_export.h>
#include <cgogn/rendering/shader_program.h>
#include <cgogn/rendering/texture.h>
namespace cgogn
{

namespace rendering
{
DECLARE_SHADER_CLASS(FlatTexture, false, CGOGN_STR(FlatTexture))

class CGOGN_RENDERING_EXPORT ShaderParamFlatTexture : public ShaderParam
{
	void set_uniforms() override;

	enum VBOName : uint32
	{
		VERTEX_POSITION = 0,
		VERTEX_TC
	};

public:
	GLVec3 light_position_;
	std::shared_ptr<Texture2D> texture_;
	bool draw_param_;

	using ShaderType = ShaderFlatTexture;

	ShaderParamFlatTexture(ShaderType* sh)
		: ShaderParam(sh), light_position_(1000, 10000, 100000), texture_(nullptr), draw_param_(false)
	{
	}
};

} // namespace rendering

} // namespace cgogn

#endif // CGOGN_RENDERING_SHADERS_FLAT_TEXTURE_H__
