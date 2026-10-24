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
#ifndef CGOGN_GEOMETRY_ALGOS_OFFSET_AXIS_PARALLEL_FOR_H_
#define CGOGN_GEOMETRY_ALGOS_OFFSET_AXIS_PARALLEL_FOR_H_

#include <cgogn/core/utils/thread_pool.h>

#include <algorithm>
#include <cstddef>
#include <future>
#include <vector>

namespace cgogn
{

namespace geometry
{

namespace detail
{

// runs f(i) for i in [0, n), in chunks on the thread pool; must not be called from a pool worker (nested wait)
template <typename FUNC>
void parallel_for_index(std::size_t n, const FUNC& f)
{
	constexpr std::size_t min_parallel_size = 1024;
	ThreadPool* pool = thread_pool();
	const std::size_t nb_workers = pool->nb_workers();
	if (n < min_parallel_size || nb_workers == 0)
	{
		for (std::size_t i = 0; i < n; ++i)
			f(i);
		return;
	}

	const std::size_t nb_chunks = nb_workers * 4;
	const std::size_t chunk_size = (n + nb_chunks - 1) / nb_chunks;
	std::vector<std::future<void>> futures;
	futures.reserve(nb_chunks);
	for (std::size_t begin = 0; begin < n; begin += chunk_size)
	{
		const std::size_t end = std::min(n, begin + chunk_size);
		futures.push_back(pool->enqueue([&f, begin, end]() {
			for (std::size_t i = begin; i < end; ++i)
				f(i);
		}));
	}
	for (auto& future : futures)
		future.wait();
}

} // namespace detail

} // namespace geometry

} // namespace cgogn

#endif // CGOGN_GEOMETRY_ALGOS_OFFSET_AXIS_PARALLEL_FOR_H_
