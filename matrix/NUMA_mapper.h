#ifndef __NUMA_MAPPER_H__
#define __NUMA_MAPPER_H__

/*
 * Copyright 2014 Open Connectome Project (http://openconnecto.me)
 * Written by Da Zheng (zhengda1936@gmail.com)
 *
 * This file is part of FlashMatrix.
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

#include <stdlib.h>
#include <math.h>
#include <assert.h>

#include <vector>

namespace fm
{

namespace detail
{

/*
 * This class maps elements in a container to multiple physical containers.
 * The number of physical containers has to be 2^n.
 * It's mainly used for distributing data to multiple NUMA nodes.
 */
class NUMA_mapper
{
	// TODO We need to try different range size to get better performance.
	static const size_t range_size_log = 20;
	static const size_t range_size = 1 << range_size_log;
	static const size_t range_mask = range_size - 1;

	const size_t numa_log;
	const size_t numa_mask;
public:
	NUMA_mapper(size_t num_nodes): numa_log(log2(num_nodes)), numa_mask(
				(1 << numa_log) - 1) {
		assert(num_nodes == 1UL << numa_log);
	}

	size_t get_num_nodes() const {
		return 1 << numa_log;
	}

	size_t get_range_size() const {
		return range_size;
	}

	size_t get_logical_range_id(size_t idx) const {
		return idx >> range_size_log;
	}

	/*
	 * This function maps the logical location in the vector to the physical
	 * location of the data. This function is a little expensive if it's
	 * invoked for every element.
	 */
	std::pair<int, size_t> map2physical(size_t idx) const {
		size_t tmp_idx = idx >> range_size_log;
		size_t range_idx = idx & range_mask;
		int node_id = tmp_idx & numa_mask;
		tmp_idx = tmp_idx >> numa_log;
		return std::pair<int, size_t>(node_id,
				(tmp_idx << range_size_log) + range_idx);
	}

	/*
	 * This method maps the offset in a raw array to the logical location
	 * in the vector.
	 */
	size_t map2logical(int node_id, size_t local_off) const {
		size_t off_in_range = local_off & range_mask;
		size_t range_id = local_off >> range_size_log;
		// The number of elements in the previous ranges on all NUMA nodes.
		return range_id * range_size * get_num_nodes()
			// The number of elements in the same range but on the NUMA nodes
			// in front of the current node.
			+ range_size * node_id
			// The number of elements in this range on the same NUMA node.
			+ off_in_range;
	}

	/*
	 * Given the number of elements, this method calculates the number of
	 * elements assigned to each NUMA node.
	 */
	std::vector<size_t> cal_local_lengths(size_t len) const;
};

}

}

#endif
