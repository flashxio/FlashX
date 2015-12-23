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

/*
 * This interface is to set the data in the vector of arrays.
 */
class set_range_operate
{
public:
	/*
	 * @buf: the buffer where to initialize data.
	 * @size: the size of the buffer (the number of bytes).
	 * @local_off: the offset of the buffer in an array (the number of bytes).
	 * @arr_id: the id of the array where the buffer is.
	 */
	virtual void set(char *buf, size_t size, off_t local_off,
			int arr_id) const = 0;
};

class reset_data_task: public thread_task
{
	chunked_raw_array arr;
public:
	reset_data_task(chunked_raw_array &_arr): arr(_arr) {
	}

	void run() {
		arr.reset_data();
	}
};

class set_data_task: public thread_task
{
	chunked_raw_array &to_buf;
	off_t to_off;
	size_t to_size;
	int node_id;
	const set_range_operate &set_range;
	size_t range_size;
public:
	// `to_off', `to_size' and `range_size' are in the number of bytes.
	set_data_task(chunked_raw_array &_arr, off_t to_off, size_t to_size, int node_id,
			const set_range_operate &_set_range, size_t range_size): to_buf(
				_arr), set_range(_set_range) {
		this->to_off = to_off;
		this->to_size = to_size;
		this->node_id = node_id;
		this->range_size = range_size;
	}

	void run() {
		for (size_t rel_off = 0; rel_off < to_size; rel_off += range_size) {
			size_t size = std::min(to_size - rel_off, range_size);
			off_t off = to_off + rel_off;
			set_range.set(to_buf.get_raw(off), size, off, node_id);
		}
	}
};

/*
 * Reset all elements in the arrays.
 */
void reset_arrays(std::vector<chunked_raw_array> &arrs);
/*
 * Set the elements in the arrays.
 * @mapper: defines how the elements are distributed in the arrays.
 * @length: the total number of elements in the arrays.
 * @entry_size: the size of an element.
 * @set_range: the function to set elements.
 * @arrs: the arrays.
 */
void set_array_ranges(const NUMA_mapper &mapper, size_t length,
		size_t entry_size, const set_range_operate &set_range,
		std::vector<chunked_raw_array> &arrs);

void reset_arrays(std::vector<chunked_raw_array> &arrs)
{
	// TODO This only runs a thread for each NUMA node.
	mem_thread_pool::ptr mem_threads
		= mem_thread_pool::get_global_mem_threads();
	for (size_t i = 0; i < arrs.size(); i++) {
		assert(arrs[i].get_node_id() >= 0);
		mem_threads->process_task(arrs[i].get_node_id(),
				new reset_data_task(arrs[i]));
	}
	mem_threads->wait4complete();
}

void set_array_ranges(const NUMA_mapper &mapper, size_t length,
		size_t entry_size, const set_range_operate &set_range,
		std::vector<chunked_raw_array> &arrs)
{
	// The number of threads per NUMA node.
	detail::mem_thread_pool::ptr mem_threads
		= detail::mem_thread_pool::get_global_mem_threads();
	size_t nthreads_per_node
		= mem_threads->get_num_threads() / matrix_conf.get_num_nodes();
	std::vector<size_t> local_lens = mapper.cal_local_lengths(length);
	for (size_t i = 0; i < arrs.size(); i++) {
		size_t num_local_bytes = local_lens[i] * entry_size;
		// The number of ranges in array of the node.
		size_t nranges = ceil_divide(local_lens[i], mapper.get_range_size());
		// The number of ranges a thread should get.
		size_t nranges_per_thread = ceil_divide(nranges, nthreads_per_node);
		// The number of bytes a thread should get.
		size_t nbytes_per_thread
			= nranges_per_thread * mapper.get_range_size() * entry_size;
		for (size_t j = 0; j < nthreads_per_node; j++) {
			if (num_local_bytes <= nbytes_per_thread * j)
				continue;
			// The number of bytes a thread actually gets.
			size_t local_nbytes = std::min(nbytes_per_thread,
					num_local_bytes - nbytes_per_thread * j);
			assert(arrs[i].get_node_id() >= 0);
			mem_threads->process_task(arrs[i].get_node_id(),
					new set_data_task(arrs[i], nbytes_per_thread * j,
						local_nbytes, arrs[i].get_node_id(), set_range,
						mapper.get_range_size() * entry_size));
		}
	}
	mem_threads->wait4complete();
}

NUMA_vec_store::ptr NUMA_vec_store::create(size_t length, const scalar_type &type,
		const std::vector<chunked_raw_array> &data, const NUMA_mapper &mapper)
{
	if (mapper.get_num_nodes() != data.size()) {
		BOOST_LOG_TRIVIAL(error)
			<< "can't create a NUMA vector. Each NUMA node should have a data array";
		return ptr();
	}
	size_t off = length - 1;
	// We should check if the last ranges can be mapped to the data arrays.
	for (size_t i = 0; i < data.size(); i++) {
		auto phy_off = mapper.map2physical(off);
		if (phy_off.second * type.get_size() >= data[phy_off.first].get_num_bytes()) {
			BOOST_LOG_TRIVIAL(error) << boost::format(
					"Data array %ld doesn't have enough bytes") % phy_off.first;
			return ptr();
		}

		if (off < mapper.get_range_size())
			break;
		off -= mapper.get_range_size();
	}
	return ptr(new NUMA_vec_store(length, type, data, mapper));
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

NUMA_vec_store::NUMA_vec_store(const NUMA_vec_store &vec): mem_vec_store(
		vec.get_length(), vec.get_type()), mapper(vec.mapper)
{
	data = vec.data;
}

NUMA_vec_store::NUMA_vec_store(size_t length, size_t num_nodes,
		const scalar_type &type): mem_vec_store(length, type), mapper(num_nodes,
			NUMA_range_size_log)
{
	data.resize(num_nodes);
	size_t num_eles_per_node = ceil_divide(length, num_nodes);
	num_eles_per_node = ROUNDUP(num_eles_per_node, mapper.get_range_size());
	size_t size_per_node = num_eles_per_node * type.get_size();
	size_t block_bytes = mapper.get_range_size() * type.get_size();
	for (size_t node_id = 0; node_id < num_nodes; node_id++)
		data[node_id] = detail::chunked_raw_array(size_per_node, block_bytes,
				node_id);
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
	const NUMA_mapper &mapper;
	size_t entry_size;
public:
	set_ele_operate(const set_vec_operate &_op, const NUMA_mapper &_mapper,
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
	return data[loc.first].get_raw(off);
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
	const NUMA_mapper &mapper;
public:
	copy_operate(const NUMA_mapper &_mapper, size_t entry_size,
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
	std::vector<size_t> local_lens = mapper.cal_local_lengths(get_length());
	// We first split data into multiple partitions and sort each partition
	// independently.
	size_t entry_size = get_type().get_size();
	for (size_t i = 0; i < data.size(); i++) {
		for (size_t j = 0; j < data[i].get_num_chunks(); j++) {
			auto chunk = data[i].get_chunk(j);
			assert(chunk.second % entry_size == 0);
			size_t size = std::min(local_lens[i] * entry_size, chunk.second);
			local_lens[i] -= size / entry_size;
			arrs.push_back(arr_pair(chunk.first, chunk.first + size));
			mem_threads->process_task(data[i].get_node_id(),
					new sort_task(chunk.first, size / entry_size, get_type()));
		}
		assert(local_lens[i] == 0);
	}
	mem_threads->wait4complete();

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
	// TODO I need to make it in parallel.
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

namespace
{

class vec2mat_set_operate: public set_operate
{
	const NUMA_vec_store &vec;
	bool byrow;
	// #rows and cols of the output matrix.
	size_t num_rows;
	size_t num_cols;
public:
	vec2mat_set_operate(const NUMA_vec_store &_vec, bool byrow,
			size_t num_rows, size_t num_cols): vec(_vec) {
		this->byrow = byrow;
		this->num_rows = num_rows;
		this->num_cols = num_cols;
	}
	virtual void set(void *arr, size_t num_eles, off_t row_idx,
			off_t col_idx) const;
	virtual const scalar_type &get_type() const {
		return vec.get_type();
	}
};

void vec2mat_set_operate::set(void *arr, size_t num_eles,
		off_t row_idx, off_t col_idx) const
{
	// The start offset in the vector.
	off_t vec_off;
	if (byrow)
		vec_off = num_cols * row_idx + col_idx;
	else
		vec_off = num_rows * col_idx + row_idx;

	while (num_eles > 0) {
		off_t next_range_start = ROUNDUP(vec_off + 1,
				vec.get_mapper().get_range_size());
		off_t end = std::min(next_range_start, (off_t) vec.get_length());
		size_t num_avails = end - vec_off;
		num_avails = std::min(num_eles, num_avails);
		const void *sub_arr = vec.get_sub_arr(vec_off, end);
		assert(sub_arr);
		memcpy(arr, sub_arr, num_avails * vec.get_entry_size());

		arr = ((char *) arr) +  num_avails * vec.get_entry_size();
		vec_off += num_avails;
		num_eles -= num_avails;
	}
}

}

matrix_store::ptr NUMA_vec_store::conv2mat(size_t nrow, size_t ncol, bool byrow)
{
	if (nrow > 1 && ncol > 1) {
		matrix_store::ptr res;
		if (byrow)
			res = NUMA_matrix_store::create(nrow, ncol, get_num_nodes(),
					matrix_layout_t::L_ROW, get_type());
		else
			res = NUMA_matrix_store::create(nrow, ncol, get_num_nodes(),
					matrix_layout_t::L_COL, get_type());
		res->set_data(vec2mat_set_operate(*this, byrow, nrow, ncol));
		return res;
	}

	matrix_store::ptr mat;
	if (nrow == 1 && byrow) {
		std::vector<NUMA_vec_store::ptr> cols(1);
		cols[0] = NUMA_vec_store::ptr(new NUMA_vec_store(*this));
		NUMA_col_tall_matrix_store::ptr tmp
			= NUMA_col_tall_matrix_store::create(cols);
		mat = NUMA_row_wide_matrix_store::create_transpose(*tmp);
	}
	else if (ncol == 1 && !byrow) {
		std::vector<NUMA_vec_store::ptr> cols(1);
		cols[0] = NUMA_vec_store::ptr(new NUMA_vec_store(*this));
		mat = NUMA_col_tall_matrix_store::create(cols);
	}
	else if (ncol == 1 && byrow)
		mat = NUMA_row_tall_matrix_store::create(data, nrow, ncol, mapper,
				get_type());
	else {
		assert(nrow == 1 && !byrow);
		NUMA_row_tall_matrix_store::ptr tmp = NUMA_row_tall_matrix_store::create(
				data, ncol, nrow, mapper, get_type());
		mat = NUMA_col_wide_matrix_store::create_transpose(*tmp);
	}
	return mat;
}

}

}
