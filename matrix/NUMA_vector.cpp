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
#include "thread.h"

#include "generic_type.h"
#include "NUMA_vector.h"
#include "data_frame.h"
#include "mem_worker_thread.h"

namespace fm
{

template<class T>
T ceil_divide(T v1, T v2)
{
	return ceil(((double) v1) / v2);
}

NUMA_vector::NUMA_vector(size_t length, size_t num_nodes,
		const scalar_type &_type): vector(length, _type.get_size(),
			true), mapper(num_nodes), type(_type)
{
	data.resize(num_nodes);
	size_t num_eles_per_node = ceil_divide(length, num_nodes);
	num_eles_per_node = ROUNDUP(num_eles_per_node, mapper.get_range_size());
	size_t size_per_node = num_eles_per_node * type.get_size();
	for (size_t node_id = 0; node_id < num_nodes; node_id++)
		data[node_id] = detail::raw_data_array(size_per_node, node_id);
}

void NUMA_vector::reset_data()
{
	reset_arrays(data);
}

const char *NUMA_vector::get_sub_arr(off_t start, off_t end) const
{
	// The start and end needs to fall into the same range.
	if (mapper.get_logical_range_id(start)
			!= mapper.get_logical_range_id(end - 1)) {
		BOOST_LOG_TRIVIAL(error) << boost::format(
				"[%1%, %2%) isn't in the same range") % start % end;
		return NULL;
	}

	std::pair<int, size_t> loc = mapper.map2physical(start);
	size_t off = loc.second * get_entry_size();
	return data[loc.first].get_raw() + off;
}

char *NUMA_vector::get_sub_arr(off_t start, off_t end)
{
	const NUMA_vector *const_this = this;
	return (char *) const_this->get_sub_arr(start, end);
}

bool NUMA_vector::is_sub_vec() const
{
	for (size_t i = 0; i < data.size(); i++)
		if (!data[i].has_entire_array())
			return true;
	return false;
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

namespace
{

class sort_task: public thread_task
{
	const scalar_type &type;
	char *buf;
	size_t num_eles;
public:
	sort_task(char *buf, size_t num_eles,
			const scalar_type &_type): type(_type) {
		this->buf = buf;
		this->num_eles = num_eles;
	}

	void run() {
		type.get_sorter().serial_sort(buf, num_eles, false);
	}
};

/*
 * This task is to copy some data in the from buffer to the right ranges.
 */
class copy_operate: public detail::set_range_operate
{
	const char *from_buf;
	size_t entry_size;
	const detail::NUMA_mapper &mapper;
public:
	copy_operate(const detail::NUMA_mapper &_mapper, size_t entry_size,
			const char *from_buf): mapper(_mapper) {
		this->entry_size = entry_size;
		this->from_buf = from_buf;
	}

	virtual void set(char *buf, size_t size, off_t local_off,
			int node_id) const {
		assert(size % entry_size == 0);
		assert(local_off % entry_size == 0);
		size_t from_off = mapper.map2logical(node_id, local_off / entry_size);
		memcpy(buf, from_buf + from_off * entry_size, size);
	}
};

}

void NUMA_vector::sort()
{
	assert(!is_sub_vec());
	// The number of threads per NUMA node.
	size_t nthreads_per_node
		= matrix_conf.get_num_threads() / matrix_conf.get_num_nodes();
	typedef std::pair<char *, char *> arr_pair;
	std::vector<arr_pair> arrs;
	detail::mem_thread_pool::ptr mem_threads
		= detail::mem_thread_pool::get_global_mem_threads();
	std::vector<size_t> local_lens = mapper.cal_local_lengths(get_length());
	// We first split data into multiple partitions and sort each partition
	// independently.
	for (size_t i = 0; i < data.size(); i++) {
		// The average number of elements processed in a thread.
		size_t len_per_thread = ceil_divide(local_lens[i], nthreads_per_node);
		for (size_t j = 0; j < nthreads_per_node; j++) {
			char *start = data[i].get_raw() + len_per_thread * j * get_entry_size();
			assert(local_lens[i] >= len_per_thread * j);
			// The number of elements processed in this thread.
			size_t local_len = std::min(len_per_thread,
					local_lens[i] - len_per_thread * j);
			arrs.push_back(std::pair<char *, char *>(start,
						start + local_len * get_entry_size()));
			mem_threads->process_task(data[i].get_node_id(),
					new sort_task(start, local_len, get_type()));
		}
	}
	mem_threads->wait4complete();
	assert(arrs.size() <= (size_t) matrix_conf.get_num_threads());

	for (size_t i = 0; i < arrs.size(); i++) {
		size_t len = ((long) (arrs[i].second - arrs[i].first)) / get_entry_size();
		assert(get_type().get_sorter().is_sorted(arrs[i].first, len, false));
	}

	// Now we should merge the array parts together.
	size_t tot_num_bytes = get_length() * get_entry_size();
	std::unique_ptr<char[]> tmp(new char[tot_num_bytes]);
	get_type().get_sorter().merge(arrs, tmp.get(), get_length());
	assert(get_type().get_sorter().is_sorted(tmp.get(), get_length(), false));

	// We need to copy the result back.
	copy_from(tmp.get(), tot_num_bytes);
}

void NUMA_vector::copy_from(const char *buf, size_t num_bytes)
{
	assert(num_bytes % get_entry_size() == 0);
	assert(num_bytes / get_entry_size() == get_length());

	copy_operate cp(mapper, get_entry_size(), buf);
	set_array_ranges(mapper, get_length(), get_entry_size(), cp, data);
}

vector::ptr NUMA_vector::sort_with_index()
{
	assert(!is_sub_vec());
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
