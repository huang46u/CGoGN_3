#ifndef CGOGN_GEOMETRY_TYPES_SPATIAL_GRID_H_
#define CGOGN_GEOMETRY_TYPES_SPATIAL_GRID_H_

#include <cgogn/core/utils/numerics.h>
#include <cgogn/geometry/types/vector_traits.h>

#include <cmath>
#include <functional>
#include <unordered_map>
#include <vector>

namespace cgogn
{

namespace geometry
{

struct SpatialGrid
{
	struct GridKey
	{
		int x, y, z;
		bool operator==(const GridKey& o) const
		{
			return x == o.x && y == o.y && z == o.z;
		}
	};
	struct GridKeyHash
	{
		std::size_t operator()(const GridKey& k) const
		{
			return std::hash<int>()(k.x) ^ (std::hash<int>()(k.y) << 1) ^ (std::hash<int>()(k.z) << 2);
		}
	};

	Scalar cell_size_;
	std::unordered_map<GridKey, std::vector<uint32>, GridKeyHash> grid_;

	explicit SpatialGrid(Scalar cell_size) : cell_size_(cell_size)
	{
	}

	GridKey get_key(const Vec3& v) const
	{
		return {(int)std::floor(v[0] / cell_size_), (int)std::floor(v[1] / cell_size_),
				(int)std::floor(v[2] / cell_size_)};
	}

	void insert(const Vec3& pos, uint32 idx)
	{
		grid_[get_key(pos)].push_back(idx);
	}

	template <typename PositionAttribute>
	bool is_valid_sample(const Vec3& pos, Scalar radius, const PositionAttribute& positions) const
	{
		GridKey k = get_key(pos);
		Scalar r2 = radius * radius;
		// Check 3x3x3 neighborhood
		for (int dx = -1; dx <= 1; ++dx)
		{
			for (int dy = -1; dy <= 1; ++dy)
			{
				for (int dz = -1; dz <= 1; ++dz)
				{
					GridKey neighbor_key = {k.x + dx, k.y + dy, k.z + dz};
					auto it = grid_.find(neighbor_key);
					if (it != grid_.end())
					{
						for (uint32 idx : it->second)
						{
							if ((positions[idx] - pos).squaredNorm() < r2)
								return false;
						}
					}
				}
			}
		}
		return true;
	}
};

} // namespace geometry

} // namespace cgogn

#endif // CGOGN_GEOMETRY_TYPES_SPATIAL_GRID_H_
