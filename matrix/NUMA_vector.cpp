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
#include "local_vec_store.h"
#include "matrix_store.h"
#include "NUMA_dense_matrix.h"

namespace fm
{

namespace detail
{

template<class T>
T ceil_divide(T v1, T v2)
{
	return ceil(((double) v1) / v2);
}

NUMA_vec_store::ptr NUMA_vec_store::cast(vec_store::ptr vec)
{
	if (!vec->is_in_mem()
			|| std::static_pointer_cast<mem_vec_store>(vec)->get_num_nodes() < 0) {
		BOOST_LOG_TRIVIAL(error) << "The vector isn't a NUMA vector";
		return NUMA_vec_store::ptr();
	}
	return std::static_pointer_cast<NUMA_vec_store>(vec);
}

NUMA_vec_store::NUMA_vec_store(size_t length, size_t num_nodes,
		const scalar_type &type): mem_vec_store(length, type), mapper(num_nodes)
{
	data.resize(num_nodes);
	size_t num_eles_per_node = ceil_divide(length, num_nodes);
	num_eles_per_node = ROUNDUP(num_eles_per_node, mapper.get_range_size());
	size_t size_per_node = num_eles_per_node * type.get_size();
	for (size_t node_id = 0; node_id < num_nodes; node_id++)
		data[node_id] = detail::raw_data_array(size_per_node, node_id);
}

void NUMA_vec_store::reset_data()
{
	reset_arrays(data);
}

namespace
{

class set_ele_operate: public detail::set_range_operate
{
	const set_vec_operate &op;
	const detail::NUMA_mapper &mapper;
	size_t entry_size;
public:
	set_ele_operate(const set_vec_operate &_op, const detail::NUMA_mapper &_mapper,
			size_t entry_size): op(_op), mapper(_mapper) {
		this->entry_size = entry_size;
	}

	void set(char *buf, size_t size, off_t local_off, int node_id) const {
		assert(local_off % entry_size == 0);
		assert(size % entry_size == 0);
		size_t num_eles = size / entry_size;
		size_t local_ele_idx = local_off / entry_size;
		size_t global_ele_idx = mapper.map2logical(node_id, local_ele_idx);
		op.set(buf, num_eles, global_ele_idx);
	}
};

}

void NUMA_vec_store::set_data(const set_vec_operate &op)
{
	set_ele_operate set_ele(op, mapper, get_entry_size());
	set_array_ranges(mapper, get_length(), get_entry_size(), set_ele, data);
}

const char *NUMA_vec_store::get_sub_arr(off_t start, off_t end) const
{
	// The start and end needs to fall into the same range.
	if (mapper.get_logical_range_id(start)
			!= mapper.get_logical_range_id(end - 1)) {
		BOOST_LOG_TRIVIAL(error) << boost::format(
				"[%1%, %2%) isn't in the same range") % start % end;
		return NULL;
	}
	if ((size_t) start >= get_length() || (size_t) end > get_length()) {
		BOOST_LOG_TRIVIAL(error)
			<< boost::format("[%1%, %2%) is out of boundary") % start % end;
		return NULL;
	}

	std::pair<int, size_t> loc = mapper.map2physical(start);
	size_t off = loc.second * get_entry_size();
	return data[loc.first].get_raw() + off;
}

char *NUMA_vec_store::get_sub_arr(off_t start, off_t end)
{
	const NUMA_vec_store *const_this = this;
	return (char *) const_this->get_sub_arr(start, end);
}

bool NUMA_vec_store::is_sub_vec() const
{
	return false;
}

bool NUMA_vec_store::append(std::vector<vec_store::const_ptr>::const_iterator vec_it,
		std::vector<vec_store::const_ptr>::const_iterator vec_end)
{
	// TODO
	assert(0);
}

bool NUMA_vec_store::append(const vec_store &vec)
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

void NUMA_vec_store::sort()
{
	assert(!is_sub_vec());
	// The number of threads per NUMA node.
	typedef std::pair<const char *, const char *> arr_pair;
	std::vector<arr_pair> arrs;
	detail::mem_thread_pool::ptr mem_threads
		= detail::mem_thread_pool::get_global_mem_threads();
	size_t nthreads_per_node
		= mem_threads->get_num_threads() / matrix_conf.get_num_nodes();
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
			arrs.push_back(arr_pair(start, start + local_len * get_entry_size()));
			mem_threads->process_task(data[i].get_node_id(),
					new sort_task(start, local_len, get_type()));
		}
	}
	mem_threads->wait4complete();
	assert(arrs.size() <= (size_t) mem_threads->get_num_threads());

	// Now we should merge the array parts together.
	size_t tot_num_bytes = get_length() * get_entry_size();
	std::unique_ptr<char[]> tmp(new char[tot_num_bytes]);
	get_type().get_sorter().merge(arrs, tmp.get(), get_length());

	// We need to copy the result back.
	copy_from(tmp.get(), tot_num_bytes);
}

bool NUMA_vec_store::copy_from(const char *buf, size_t num_bytes)
{
	if (num_bytes % get_entry_size() != 0
			|| num_bytes / get_entry_size() != get_length())
		return false;

	copy_operate cp(mapper, get_entry_size(), buf);
	set_array_ranges(mapper, get_length(), get_entry_size(), cp, data);
	return true;
}

bool NUMA_vec_store::copy_from(const NUMA_vec_store &vec)
{
	if (get_type() != vec.get_type()) {
		BOOST_LOG_TRIVIAL(error)
			<< "can't copy from a vector of different type";
		return false;
	}

	if (mapper != vec.mapper) {
		BOOST_LOG_TRIVIAL(error) << "The vector has different NUMA setup";
		return false;
	}

	// TODO we should do it in parallel.
	for (size_t i = 0; i < data.size(); i++)
		if (!data[i].copy_from(vec.data[i]))
			return false;
	return true;
}

vec_store::ptr NUMA_vec_store::sort_with_index()
{
	assert(!is_sub_vec());
	// TODO
	assert(0);
}

bool NUMA_vec_store::is_sorted() const
{
	// TODO
	assert(0);
}

vec_store::ptr NUMA_vec_store::deep_copy() const
{
	NUMA_vec_store *ret = new NUMA_vec_store(*this);
	for (size_t i = 0; i < data.size(); i++)
		ret->data[i] = data[i].deep_copy();
	return vec_store::ptr(ret);
}

local_vec_store::const_ptr NUMA_vec_store::get_portion(off_t start,
		size_t length) const
{
	const char *raw_data = get_sub_arr(start, start + length);
	if (raw_data == NULL)
		return local_vec_store::const_ptr();

	std::pair<int, size_t> loc = mapper.map2physical(start);
	return local_vec_store::const_ptr(new local_cref_vec_store(raw_data,
				start, length, get_type(), data[loc.first].get_node_id()));
}

local_vec_store::ptr NUMA_vec_store::get_portion(off_t start, size_t length)
{
	char *raw_data = get_sub_arr(start, start + length);
	if (raw_data == NULL)
		return local_vec_store::ptr();

	std::pair<int, size_t> loc = mapper.map2physical(start);
	return local_vec_store::ptr(new local_ref_vec_store(raw_data,
				start, length, get_type(), data[loc.first].get_node_id()));
}

matrix_store::const_ptr NUMA_vec_store::conv2mat(size_t nrow, size_t ncol,
			bool byrow) const
{
	std::vector<NUMA_vec_store::ptr> cols(1);
	cols[0] = NUMA_vec_store::ptr(new NUMA_vec_store(*this));
	return NUMA_col_tall_matrix_store::create(cols);
}

}

}
