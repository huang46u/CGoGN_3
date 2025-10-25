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

#ifndef CGOGN_GEOMETRY_ALGOS_VOLUME_POISSON_DISK_SAMPLER_H_
#define CGOGN_GEOMETRY_ALGOS_VOLUME_POISSON_DISK_SAMPLER_H_
#include <algorithm>
#include <cmath>
#include <random>
#include <unordered_map>
namespace cgogn
{

namespace geometry
{
using Vec3 = cgogn::geometry::Vec3;
using Scalar = cgogn::geometry::Scalar;

// Spread 10 bits of n so that there are 3 zeros between each bit
inline uint32 spreadbits(uint32 n)
{
	n &= 0x000003ff; // we only look at the first 10 bits
	n = (n ^ (n << 16)) & 0x030000ff;
	n = (n ^ (n << 8)) & 0x0300f00f;
	n = (n ^ (n << 4)) & 0x030c30c3;
	n = (n ^ (n << 2)) & 0x09249249;
	return n;
}

inline uint32 morton3D(uint32 x, uint32 y, uint32 z)
{
	return (spreadbits(x) | (spreadbits(y) << 1) | (spreadbits(z) << 2));
}

class OctreeLocationCode
{
public:
	static OctreeLocationCode root()
	{
		return OctreeLocationCode{};
	}

	void clear()
	{
		bits_ = 0;
		depth_ = 0;
	}

	void push_octant(uint8 oct)
	{
		oct &= 7u;
		bits_ |= (uint64(oct) << (3u * depth_));
		++depth_;
	}

	OctreeLocationCode child(uint8 oct) const
	{
		OctreeLocationCode c = *this;
		c.push_octant(oct);
		return c;
	}

	void parent()
	{
		if (depth_ == 0)
			return;
		--depth_;
		bits_ &= mask_for_depth(depth_);
	}

	uint8 octant_at(uint8 level) const
	{
		return uint8((bits_ >> (3u * level)) & 7u);
	}

	OctreeLocationCode truncated(uint8 d) const
	{
		if (d >= depth_)
			return *this;
		OctreeLocationCode c;
		c.bits_ = bits_ & mask_for_depth(d);
		c.depth_ = d;
		return c;
	}

	static OctreeLocationCode from_ijk(uint32 ix, uint32 iy, uint32 iz, uint8 depth)
	{
		OctreeLocationCode c;
		c.depth_ = depth;
		uint64 b = 0;
		for (uint8 k = 0; k < depth; ++k)
		{
			uint64 oct = ((ix >> k) & 1u) | (((iy >> k) & 1u) << 1u) | (((iz >> k) & 1u) << 2u);
			b |= (oct << (3ull * k));
		}
		c.bits_ = b;
		return c;
	}

	std::tuple<uint32, uint32, uint32> to_ijk() const
	{
		uint32 ix = 0;
		uint32 iy = 0;
		uint32 iz = 0;
		for (uint8 k = 0; k < depth_; ++k)
		{
			uint8 oct = octant_at(k);
			ix |= ((oct & 1u) << k);
			iy |= (((oct >> 1u) & 1u) << k);
			iz |= (((oct >> 2u) & 1u) << k);
		}
		return {ix, iy, iz};
	}

	static constexpr uint32 reverse_bits_depth(uint32 v, uint8 depth) noexcept
	{
		uint32 r = 0;
		for (uint8 k = 0; k < depth; ++k)
		{
			r |= ((v >> k) & 1u) << (depth - 1 - k);
		}
		return r;
	}

	inline std::tuple<uint32, uint32, uint32> to_grid_ijk() const
	{
		auto [ix, iy, iz] = to_ijk();
		const uint8 d = depth_;
		return {reverse_bits_depth(ix, d), 
				reverse_bits_depth(iy, d), reverse_bits_depth(iz, d)};
	}

	static inline OctreeLocationCode from_grid_ijk(uint32 gx, uint32 gy, uint32 gz, uint8 depth)
	{
		return from_ijk(reverse_bits_depth(gx, depth), reverse_bits_depth(gy, depth), reverse_bits_depth(gz, depth),
						depth);
	}
	uint8 depth() const
	{
		return depth_;
	}
	static constexpr uint64 mask_for_depth(uint8 d)
	{
		return (d == 0) ? 0ull : ((1ull << (3ull * d)) - 1ull);
	}

	std::string to_string() const
	{
		std::string s = "d=" + std::to_string(depth_) + " [";
		for (uint8 k = 0; k < depth_; ++k)
		{
			s += std::to_string(octant_at(k));
			if (k + 1 < depth_)
				s += ",";
		}
		s += "]";
		return s;
	}

	bool operator==(const OctreeLocationCode& o) const noexcept
	{
		return bits_ == o.bits_ && depth_ == o.depth_;
	}
	bool operator!=(const OctreeLocationCode& o) const noexcept
	{
		return !(*this == o);
	}

	struct Hasher
	{
		size_t operator()(const OctreeLocationCode& c) const noexcept
		{
			// mix bits and depth (put depth in top bits to avoid collisions across depths)
			uint64_t h = c.bits_ ^ (uint64_t(c.depth_) << 61);
			// final avalanching
			h ^= h >> 33;
			h *= 0xff51afd7ed558ccdULL;
			h ^= h >> 33;
			h *= 0xc4ceb9fe1a85ec53ULL;
			h ^= h >> 33;
			return size_t(h);
		}
	};

	uint64 bits() const { return bits_; }
	uint8 depth_value() const { return depth_; }

private:
	uint64 bits_ = 0;
	uint8 depth_ = 0;
};

class SparseGrid
{
public:
	using Key = uint32;
	using Grid_Index = std::tuple<int, int, int>;

	SparseGrid(uint32 res) : res_(res)
	{
	}

	// New method: check neighbors INCLUDING the center cell
	template <class Func>
	bool adjacent_hit_local(const Vec3& local01, Scalar r_world, Scalar cell_size, Func&& func) const
	{
		Grid_Index q = cell_of_local(local01);
		int w = int(std::ceil(r_world / cell_size));
		
		const int ix =std::get<0>(q);
		const int iy =std::get<1>(q);
		const int iz =std::get<2>(q);
		const int loX = std::clamp(ix - w, 0, int(res_) - 1);
		const int hiX = std::clamp(ix + w, 0, int(res_) - 1);
		const int loY = std::clamp(iy - w, 0, int(res_) - 1);
		const int hiY = std::clamp(iy + w, 0, int(res_) - 1);
		const int loZ = std::clamp(iz - w, 0, int(res_) - 1);
		const int hiZ = std::clamp(iz + w, 0, int(res_) - 1);

		for (int z = loZ; z <= hiZ; ++z)
		{
			for (int y = loY; y <= hiY; ++y)
			{
				for (int x = loX; x <= hiX; ++x)
				{
					if (!is_legal(Grid_Index{x, y, z}))
						continue;
					Key key = morton3D(uint32(x), uint32(y), uint32(z));
					auto it = cells_.find(key);
					if (it == cells_.end())
						continue;
					if (func(it->second))
						return true;
				}
			}
		}
		return false;
	}
	
	bool occupied(const Grid_Index& idx) const
	{
		return occupied(morton_of(idx));
	}
	bool occupied(const uint32 key) const
	{
		return cells_.find(key) != cells_.end();
	}
	int get(const Grid_Index& idx) const
	{
		return get(morton_of(idx));
	}
	int get(const uint32 key) const
	{
		auto it = cells_.find(key);
		if (it != cells_.end())
			return it->second;
		return -1;
	}
	void set(const Grid_Index& idx, uint32 value)
	{
		set(morton_of(idx), value);
	}
	void set(const uint32 key, uint32 value)
	{
		cells_[key] = value;
	}

	Grid_Index cell_of_local(const Vec3& local01) const
	{
		Scalar fx = local01.x() * Scalar(res_);
		Scalar fy = local01.y() * Scalar(res_);
		Scalar fz = local01.z() * Scalar(res_);
		int x = int(std::floor(fx));
		int y = int(std::floor(fy));
		int z = int(std::floor(fz));
		return {x, y, z};
	}

	Key morton_of(const Grid_Index& q) const
	{
		return morton3D(uint32(std::get<0>(q)), std::get<1>(q), std::get<2>(q));
	}

	uint32 resolution() const
	{
		return res_;
	}

	bool empty() const
	{
		return cells_.empty();
	}

	std::size_t cell_count() const
	{
		return cells_.size();
	}

	void clear()
	{
		cells_.clear();
	}

	bool is_legal(Grid_Index idx) const
	{
		return (std::get<0>(idx) >= 0 && std::get<0>(idx) < int(res_) && std::get<1>(idx) >= 0 &&
				std::get<1>(idx) < int(res_) && std::get<2>(idx) >= 0 && std::get<2>(idx) < int(res_));
	}

private:
	

	std::unordered_map<uint32, uint32> cells_; // morton code -> sample ID
	uint32 res_;
};

class OctreeNode
{
public:
	OctreeNode(const OctreeLocationCode& loc, const Vec3& center, Scalar half)
		: location_(loc), center_(center), half_size_(half), grid_(16)
	{
		cell_size_ = (2.0 * half_size_) / Scalar(grid_.resolution());
		local_radius_ = cell_size_ * std::sqrt(3.0) / 2.0;
	}

	template <typename ConflictPred>
	bool can_insert(const Vec3& p, ConflictPred&& is_conflict)
	{
		Vec3 local = to_local01(p);
		auto idx = grid_.cell_of_local(local);
		if (!grid_.is_legal(idx))
			return false;
		if (grid_.occupied(idx))
			return false;
		if (grid_.adjacent_hit_local(local, local_radius_, cell_size_,
									 [&](uint32 sampleId) { return is_conflict(sampleId, p, local_radius_); }))
			return false;
		
		return true;
	}

	Scalar cell_size() const
	{
		return cell_size_;
	}

	Vec3 to_local01(const Vec3& p) const
	{
		const Scalar inv = 1.0 / (2.0 * half_size_);
		Vec3 local = (p - (center_ - Vec3(half_size_, half_size_, half_size_))) * inv;
		return local;
	}

	SparseGrid& grid()
	{
		return grid_;
	}

	Scalar radius() const
	{
		return local_radius_;
	}

	Vec3 random_sample_around(const Vec3& p, std::mt19937& rng, std::uniform_real_distribution<Scalar>& Uni)
	{
		Scalar u = Uni(rng);
		Scalar v = Uni(rng);
		Scalar w = Uni(rng);

		Scalar R3 = local_radius_ * local_radius_ * local_radius_;
		Scalar r = std::cbrt(R3 + u * (8 * R3 - R3)); // r = pow((R^3 + u(8R^3 - R^3)), 1/3)

		Scalar phi = v * 2.0 * M_PI;
		Scalar theta = std::acos(1.0 - 2.0 * w);

		Scalar x = r * std::sin(theta) * std::cos(phi);
		Scalar y = r * std::sin(theta) * std::sin(phi);
		Scalar z = r * std::cos(theta);

		return p + Vec3(x, y, z);
	}
	inline uint8 octant_of(const Vec3& p) const
	{
		
		return (p.x() >= center_.x()) | ((p.y() >= center_.y()) << 1) | ((p.z() >= center_.z()) << 2);
		
	}

	bool sphere_overlap_node(const Vec3& q, Scalar r)
	{
		const Scalar ax = std::abs(q.x() - center_.x());
		const Scalar ay = std::abs(q.y() - center_.y());
		const Scalar az = std::abs(q.z() - center_.z());
		const Scalar dx = std::max(Scalar(0), ax - half_size_);
		const Scalar dy = std::max(Scalar(0), ay - half_size_);
		const Scalar dz = std::max(Scalar(0), az - half_size_);
		return (dx * dx + dy * dy + dz * dz) <= r * r;
	}

private:
	OctreeLocationCode location_;
	Vec3 center_;
	Scalar half_size_;
	Scalar cell_size_;
	Scalar local_radius_;
	SparseGrid grid_;
};
class NestedOctree
{

public:
	using NodeVecIdx = std::size_t;
	using Node = OctreeNode;
	using Grid_Index = typename SparseGrid::Grid_Index;

	NestedOctree() = default;

	Node& get_or_add_node(const OctreeLocationCode& loc)
	{
		auto it = node_indices_.find(loc);
		if (it != node_indices_.end())
		{
			return nodes_[it->second];
		}
		// create new node
		Vec3 center(0.5, 0.5, 0.5);
		Scalar half = 0.5;
		for (uint8 d = 0; d < loc.depth(); ++d)
		{
			uint8 oct = loc.octant_at(d);
			half *= 0.5;
			if (oct & 1u)
				center.x() += half;
			else
				center.x() -= half;
			if (oct & 2u)
				center.y() += half;
			else
				center.y() -= half;
			if (oct & 4u)
				center.z() += half;
			else
				center.z() -= half;
		}
		OctreeNode node(loc, center, half);
		nodes_.push_back(node);
		NodeVecIdx idx = nodes_.size() - 1;
		node_indices_[loc] = idx;
		return nodes_[idx];
	}

	Vec3 next_pos(const Vec3& pos, uint8 depth, std::mt19937& rng, std::uniform_real_distribution<Scalar>& Uni)
	{
		OctreeLocationCode loc = OctreeLocationCode::root();

		for (uint8 d = 0; d < depth; ++d)
		{
			Node& node = get_or_add_node(loc);
			uint8 oct = node.octant_of(pos);
			loc.push_octant(oct);
		}
		Node& node = get_or_add_node(loc);
		return node.random_sample_around(pos, rng, Uni);
	}
	template <typename ConflictPred, typename AcceptCallBack>
	bool insert(const Vec3& p, uint8 depth, ConflictPred&& is_conflict, AcceptCallBack&& on_accept)
	{
		OctreeLocationCode loc = OctreeLocationCode::root();
		for (uint8 d = 0; d < depth; ++d)
		{
			Node& node = get_or_add_node(loc);
			uint8 oct = node.octant_of(p);
			loc.push_octant(oct);
		}
		Node& node = get_or_add_node(loc);
		
		if (!node.can_insert(p, is_conflict))
			return false;
		
		Vec3 local = node.to_local01(p);
		auto idx = node.grid().cell_of_local(local);
		const Scalar R = node.radius();
		const Scalar cs = node.cell_size();
		
		if (conflict_ancestors(loc, p, R, cs, is_conflict))
			return false;
		
		if (conflict_same_depth_neighbors(loc, p, R, cs, is_conflict))
			return false;
		
		uint32 sampleId = on_accept(p, depth);
		node.grid().set(idx, sampleId);
		return true;
	}
	void clear()
	{
		nodes_.clear();
		node_indices_.clear();
	}

private:
	template <class ConflictPred>
	bool conflict_ancestors(const OctreeLocationCode& loc, const Vec3& p, Scalar R, Scalar cell_size, ConflictPred&& is_conflict)
	{
		const uint8 depth = loc.depth();
		if (depth == 0)
			return false;

		for (int ancestor_depth = int(depth); ancestor_depth > 0; --ancestor_depth)
		{
			OctreeLocationCode ac = loc.truncated(uint8(ancestor_depth - 1));

			if (auto it = node_indices_.find(ac); it != node_indices_.end())
			{
				Node& an = nodes_[it->second];
				
				const Vec3 local = an.to_local01(p);
				
				const Scalar ancestor_cell_size = an.cell_size();
				
				if (an.grid().adjacent_hit_local(local, R, ancestor_cell_size,
												 [&](int sid) { 
													return is_conflict(uint32(sid), p, R); 
												}))
				{
					return true;
				}
			}

			const uint8 k = ac.depth();
			const int N = 1 << k;
			auto t = ac.to_grid_ijk();
			const int ix = static_cast<int>(std::get<0>(t));
			const int iy = static_cast<int>(std::get<1>(t));
			const int iz = static_cast<int>(std::get<2>(t));

			for (int dz = -1; dz <= 1; ++dz)
				for (int dy = -1; dy <= 1; ++dy)
					for (int dx = -1; dx <= 1; ++dx)
					{
						if (dx == 0 && dy == 0 && dz == 0)
							continue;
						int nx = ix + dx, ny = iy + dy, nz = iz + dz;
						if (nx < 0 || ny < 0 || nz < 0 || nx >= N || ny >= N || nz >= N)
							continue;

						OctreeLocationCode neighbor_code = OctreeLocationCode::from_grid_ijk(nx, ny, nz, k);
						auto neighbor_it = node_indices_.find(neighbor_code);
						if (neighbor_it == node_indices_.end())
							continue;

						Node& neighbor = nodes_[neighbor_it->second];
						Vec3 nlocal = neighbor.to_local01(p);
						
						if (!neighbor.sphere_overlap_node(p, R))
							continue;
						const Scalar neighbor_cell_size = neighbor.cell_size();
						
						if (neighbor.grid().adjacent_hit_local(nlocal, R, neighbor_cell_size,
															   [&](int sid) {
								return is_conflict(uint32(sid), p, R);
							}))
							return true;
					}
		}
		return false;
	}

	template <typename ConflictPred>
	bool conflict_same_depth_neighbors(const OctreeLocationCode& loc, const Vec3& p, Scalar R, Scalar cell_size,
									   ConflictPred&& is_conflict)
	{
		auto t = loc.to_grid_ijk();
		const int ix = static_cast<int>(std::get<0>(t));
		const int iy = static_cast<int>(std::get<1>(t));
		const int iz = static_cast<int>(std::get<2>(t));
		const int N = 1 << loc.depth();

		for (int dz = -1; dz <= 1; ++dz)
			for (int dy = -1; dy <= 1; ++dy)
				for (int dx = -1; dx <= 1; ++dx)
				{
					if (dx == 0 && dy == 0 && dz == 0)
						continue;
					int nx = ix + dx, ny = iy + dy, nz = iz + dz;
					if (nx < 0 || ny < 0 || nz < 0 || nx >= N || ny >= N || nz >= N)
						continue;

					OctreeLocationCode ncode =
						OctreeLocationCode::from_grid_ijk(uint32(nx), uint32(ny), uint32(nz), loc.depth());
					auto it = node_indices_.find(ncode);
					if (it == node_indices_.end())
						continue;

					Node& nb = nodes_[it->second];
					if (!nb.sphere_overlap_node(p, R))
						continue;
					
					Vec3 nlocal = nb.to_local01(p);
					
				const Scalar neighbor_cell_size = nb.cell_size();
					if (nb.grid().adjacent_hit_local(nlocal, R, neighbor_cell_size,
															 [&](int sid) { return is_conflict(sid, p, R); }))
					{
						return true;
					}
				}
		return false;
	}

private:
	uint8 max_depth_;
	std::vector<OctreeNode> nodes_;
	std::unordered_map<OctreeLocationCode, NodeVecIdx, OctreeLocationCode::Hasher> node_indices_;
};

template <typename MESH>
class VolumePoissonDiskSampler
{

	using Vertex = typename mesh_traits<MESH>::Vertex;
	template <typename T>
	using Attribute = typename mesh_traits<MESH>::template Attribute<T>;
	static constexpr uint8 MAX_DEPTH = 21;

public:
	VolumePoissonDiskSampler(MESH& m, uint32 max_trials = 30, uint32 seed = 42)
		: mesh_(m), max_trials_(max_trials), rng_(seed), uni_(0.0, 1.0)
	{
		sample_position_ = get_or_add_attribute<Vec3, Vertex>(mesh_, "position");
		poisson_sample_depth_ = get_or_add_attribute<uint8, Vertex>(mesh_, "poisson_sample_depth");
		poisson_sample_radius_ = get_or_add_attribute<Scalar, Vertex>(mesh_, "poisson_sample_radius");
	}

	template <class Domain, class OnAccept>
	uint32 sample_fill_at_depth(uint8 depth, OnAccept&& post, Domain&& domain,
								uint32 target_count = (std::numeric_limits<uint32>::max)())
	{

		uint32 added = 0;
		auto is_conflict = [&](uint32 sid, const Vec3& q, Scalar R_candidate) -> bool {
			Vec3 s_pos = (*sample_position_)[sid];
			Scalar R_existing = (*poisson_sample_radius_)[sid];
			// Use the smaller radius for conflict detection (cross-layer constraint)
			Scalar R_threshold = std::min(R_candidate, R_existing);
			Scalar dist_sq = (s_pos - q).dot(s_pos - q);
			bool conflict = dist_sq < (R_threshold * R_threshold);
			
			return conflict;
		};
		auto on_accept = [&](const Vec3& pos, const uint8 depth) -> uint32 {
			Vertex v = add_vertex(mesh_);
			uint32 vid = index_of(mesh_, v);
			active_list_.push_back(vid);
			
			
			Scalar node_size = 1.0 / Scalar(1u << depth);
			Scalar cell_size = node_size / 16.0; 
			Scalar radius = cell_size * std::sqrt(3.0) / 2.0;
			(*poisson_sample_radius_)[vid] = radius;
			
			post(pos, v, depth);
			return vid;
		};
		const uint32 max_seed_trials = std::max<uint32>(256, 10 * max_trials_);
		if (active_list_.empty())
		{
			Vec3 seed_pos;
			bool accetped = false;
			for (int i = 0; i < max_trials_; ++i)
			{
				if (pick_seed(domain, seed_pos, max_seed_trials))
				{

					if (octree_.insert(seed_pos, depth, is_conflict, on_accept))
					{
						accetped = true;
						++added;
						break;
					}
				}
			}
			if (!accetped)
			{
				return 0;
			}
			if (added >= target_count)
				return added;
		}
		while (!active_list_.empty() && added < target_count)
		{
			size_t idx = uni_(rng_) * active_list_.size();
			const uint32 seed_vid = active_list_[idx];
			const Vec3 current_pos = (*sample_position_)[seed_vid];
			bool accepted = false;
			for (uint32 i = 0; i < max_trials_; i++)
			{
				Vec3 candidate = octree_.next_pos(current_pos, depth, rng_, uni_);
				if (!domain(candidate))
					continue;
				if (octree_.insert(candidate, depth, is_conflict, on_accept))
				{
					++added;
					accepted = true;
					break;
				}
			}
			if (!accepted)
			{
				active_list_[idx] = active_list_.back();
				active_list_.pop_back();
				if (active_list_.empty() && added < target_count)
				{
					Vec3 seed_pos;
					bool accetped = false;
					for (int i = 0; i < max_trials_; ++i)
					{
						if (pick_seed(domain, seed_pos, max_seed_trials))
						{

							if (octree_.insert(seed_pos, depth, is_conflict, on_accept))
							{
								accetped = true;
								++added;
								break;
							}
						}
					}
					if (!accetped)
					{
						break;
					}
				}
			}
		}
		return added;
	}

	template <class Domain, class ClusterDomain, class OnAccept>
	uint32 sample_cluster_at_depth(uint8 depth, OnAccept&& post, Domain&& domain, ClusterDomain&& cluster_domain,
								   Vec3 bb_min, Vec3 bb_max, uint32 target_count = (std::numeric_limits<uint32>::max)())
	{

		uint32 added = 0;
		auto is_conflict = [&](uint32 sid, const Vec3& q, Scalar R_candidate) -> bool {
			Vec3 s_pos = (*sample_position_)[sid];
			Scalar R_existing = (*poisson_sample_radius_)[sid];
			// Use the smaller radius for conflict detection (cross-layer constraint)
			Scalar R_threshold = std::min(R_candidate, R_existing);
			return (s_pos - q).dot(s_pos - q) < (R_threshold * R_threshold);
		};
		auto on_accept = [&](const Vec3& pos, const uint8 depth) -> uint32 {
			Vertex v = add_vertex(mesh_);
			uint32 vid = index_of(mesh_, v);
			active_list_.push_back(vid);
			
			Scalar node_size = 1.0 / Scalar(1u << depth);
			Scalar cell_size = node_size / 16.0; 
			Scalar radius = cell_size * std::sqrt(3.0) / 2.0;
			(*poisson_sample_radius_)[vid] = radius;
			
			post(pos, v, depth);
			return vid;
		};
		const uint32 max_seed_trials = std::max<uint32>(256, 10 * max_trials_);
		if (active_list_.empty())
		{
			Vec3 seed_pos;
			bool accetped = false;
			for (int i = 0; i < max_trials_; ++i)
			{
				if (pick_seed_in_bbox(bb_min, bb_max, domain, cluster_domain, seed_pos, max_seed_trials))
				{

					if (octree_.insert(seed_pos, depth, is_conflict, on_accept))
					{
						accetped = true;
						++added;
						break;
					}
				}
			}
			if (!accetped)
			{
				return 0;
			}
			if (added >= target_count)
				return added;
		}
		while (!active_list_.empty() && added < target_count)
		{
			size_t idx = uni_(rng_) * active_list_.size();
			const uint32 seed_vid = active_list_[idx];
			const Vec3 current_pos = (*sample_position_)[seed_vid];
			bool accepted = false;
			for (uint32 i = 0; i < max_trials_; i++)
			{
				Vec3 candidate = octree_.next_pos(current_pos, depth, rng_, uni_);
				if (!domain(candidate) || !cluster_domain(candidate))
					continue;
				if (octree_.insert(candidate, depth, is_conflict, on_accept))
				{
					++added;
					accepted = true;
					break;
				}
			}
			if (!accepted)
			{
				active_list_[idx] = active_list_.back();
				active_list_.pop_back();
				if (active_list_.empty() && added < target_count)
				{
					Vec3 seed_pos;
					bool accetped = false;
					for (int i = 0; i < max_trials_; ++i)
					{
						if (pick_seed_in_bbox(bb_min, bb_max, domain, cluster_domain, seed_pos, max_seed_trials))
						{

							if (octree_.insert(seed_pos, depth, is_conflict, on_accept))
							{
								accetped = true;
								++added;
								break;
							}
						}
					}
					if (!accetped)
					{
						break;
					}
				}
			}
		}
		return added;
	}

	template <typename OnAccept, typename Domain, typename ClusterDomain>
	uint32 sample_cluster(Vec3& center, Scalar r, std::vector<Vertex>& cluster, OnAccept&& post, Domain&& domain,
						  ClusterDomain&& cluster_domain, uint32 target_count = (std::numeric_limits<uint32>::max)(),
						  uint8 default_depth = 1)
	{
		active_list_.clear();
		Vec3 bb_min, bb_max;
		if (cluster.empty())
		{
			bb_min = center - Vec3(r, r, r);
			bb_max = center + Vec3(r, r, r);
		}
		else
		{
			auto bbox = compute_cluster_bbox(cluster);
			bb_min = bbox.first;
			bb_max = bbox.second;
		}
		uint8 work_depth = default_depth;
		
		uint32 remaining = target_count;
		while (remaining > 0)
		{
			active_list_.clear();
			for (Vertex v : cluster)
			{
				const uint32 v_index = index_of(mesh_, v);
				if ((*poisson_sample_depth_)[v_index] != work_depth)
					continue;
				const Vec3& p = (*sample_position_)[v_index];
				active_list_.push_back(v_index);
			}
			const uint32 added =
				sample_cluster_at_depth(work_depth, post, domain, cluster_domain, bb_min, bb_max, remaining);
			if (added > remaining)
				return target_count;
			remaining -= added;
			if (work_depth >= MAX_DEPTH)
			{
				return target_count - remaining;
			}
			++work_depth;
		}
		return target_count;
	}

	void clear()
	{
		octree_.clear();
		active_list_.clear();
	}

private:
	std::pair<Vec3, Vec3> compute_cluster_bbox(std::vector<Vertex>& cluster)
	{
		Vec3 bbox_min(std::numeric_limits<Scalar>::max(), std::numeric_limits<Scalar>::max(),
					  std::numeric_limits<Scalar>::max());
		Vec3 bbox_max(std::numeric_limits<Scalar>::lowest(), std::numeric_limits<Scalar>::lowest(),
					  std::numeric_limits<Scalar>::lowest());
		for (Vertex v : cluster)
		{
			const Vec3& p = (*sample_position_)[index_of(mesh_, v)];
			bbox_min.x() = std::min(bbox_min.x(), p.x());
			bbox_min.y() = std::min(bbox_min.y(), p.y());
			bbox_min.z() = std::min(bbox_min.z(), p.z());
			bbox_max.x() = std::max(bbox_max.x(), p.x());
			bbox_max.y() = std::max(bbox_max.y(), p.y());
			bbox_max.z() = std::max(bbox_max.z(), p.z());
		}
		return {bbox_min, bbox_max};
	}

	template <typename Domain>
	bool pick_seed(Domain&& domain, Vec3& out, uint32 max_trials = 1000)
	{
		while ((max_trials--) > 0)
		{
			Scalar x = uni_(rng_);
			Scalar y = uni_(rng_);
			Scalar z = uni_(rng_);
			Vec3 p(x, y, z);
			if (domain(p))
			{
				out = p;
				return true;
			}
		}
		return false;
	}

	template <typename Domain, typename LocalDomain>
	bool pick_seed_in_bbox(const Vec3 bbox_min, const Vec3 bbox_max, Domain&& domain, LocalDomain&& local_domain,
						   Vec3& out, uint32 max_trials = 1000)
	{
		const Vec3 extent = bbox_max - bbox_min;

		for (uint32 i = 0; i < max_trials; ++i)
		{
			Scalar x = uni_(rng_) * extent.x() + bbox_min.x();
			Scalar y = uni_(rng_) * extent.y() + bbox_min.y();
			Scalar z = uni_(rng_) * extent.z() + bbox_min.z();
			Vec3 p(x, y, z);
			if (domain(p) && local_domain(p))
			{
				out = p;
				return true;
			}
		}
		return false;
	}

private:
	MESH& mesh_;
	NestedOctree octree_;
	std::shared_ptr<Attribute<uint8>> poisson_sample_depth_;
	std::shared_ptr<Attribute<Vec3>> sample_position_;
	std::shared_ptr<Attribute<Scalar>> poisson_sample_radius_;
	uint32 max_trials_;
	std::mt19937 rng_;
	std::uniform_real_distribution<Scalar> uni_;

	std::vector<uint32> active_list_;
};
} // namespace geometry
} // namespace cgogn

#endif // CGOGN_GEOMETRY_ALGOS_VOLUME_POISSON_DISK_SAMPLER_H_
