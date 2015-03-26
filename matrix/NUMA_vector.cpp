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

#include <boost/format.hpp>

#include "common.h"

#include "NUMA_vector.h"
#include "data_frame.h"

namespace fm
{

NUMA_vector::NUMA_vector(size_t length, size_t num_nodes,
		const scalar_type &_type): vector(length, _type.get_size(),
			true), numa_log(log2(num_nodes)), numa_mask(
				(1 << numa_log) - 1), type(_type)
{
	assert(num_nodes == 1UL << numa_log);
	data.resize(num_nodes);
	size_t num_eles_per_node = ceil(((double) length) / num_nodes);
	num_eles_per_node = ROUNDUP(num_eles_per_node, range_size);
	size_t size_per_node = num_eles_per_node * type.get_size();
	for (size_t node_id = 0; node_id < num_nodes; node_id++)
		data[node_id] = detail::raw_data_array(size_per_node, node_id);
}

void NUMA_vector::reset_data()
{
	for (size_t i = 0; i < data.size(); i++)
		data[i].reset_data();
}

const char *NUMA_vector::get_sub_arr(off_t start, off_t end) const
{
	off_t rid1 = start >> range_size_log;
	off_t rid2 = (end - 1) >> range_size_log;
	// The start and end needs to fall into the same range.
	if (rid1 != rid2) {
		BOOST_LOG_TRIVIAL(error) << boost::format(
				"[%1%, %2%) isn't in the same range") % start % end;
		return NULL;
	}

	std::pair<int, size_t> loc = map2data(start);
	size_t off = loc.second * get_entry_size();
	return data[loc.first].get_raw() + off;
}

char *NUMA_vector::get_sub_arr(off_t start, off_t end)
{
	const NUMA_vector *const_this = this;
	return (char *) const_this->get_sub_arr(start, end);
}

vector::const_ptr NUMA_vector::get_sub_vec(off_t start, size_t length) const
{
	// TODO
	assert(0);
}

bool NUMA_vector::expose_sub_vec(off_t start, size_t length)
{
	// TODO
	assert(0);
}

bool NUMA_vector::append(std::vector<vector::ptr>::const_iterator vec_it,
		std::vector<vector::ptr>::const_iterator vec_end)
{
	// TODO
	assert(0);
}

bool NUMA_vector::append(const vector &vec)
{
	// TODO
	assert(0);
}

void NUMA_vector::sort()
{
	// TODO
	assert(0);
}

vector::ptr NUMA_vector::sort_with_index()
{
	// TODO
	assert(0);
}

bool NUMA_vector::is_sorted() const
{
	// TODO
	assert(0);
}

// It should return data frame instead of vector.
data_frame::ptr NUMA_vector::groupby(const gr_apply_operate<mem_vector> &op,
		bool with_val) const
{
	// TODO
	assert(0);
}

vector::ptr NUMA_vector::deep_copy() const
{
	// TODO
	assert(0);
}

vector::ptr NUMA_vector::shallow_copy()
{
	// TODO
	assert(0);
}

vector::const_ptr NUMA_vector::shallow_copy() const
{
	// TODO
	assert(0);
}

}
