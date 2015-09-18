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

#include <numa.h>
#include <malloc.h>

#include "thread.h"

#include "mem_worker_thread.h"
#include "raw_data_array.h"
#include "NUMA_mapper.h"
#include "matrix_config.h"
#include "bulk_operate.h"
#include "local_mem_buffer.h"

namespace fm
{

namespace detail
{

namespace
{

template<class T>
T ceil_divide(T v1, T v2)
{
	return ceil(((double) v1) / v2);
}

class NUMA_deleter
{
	size_t size;
public:
	NUMA_deleter(size_t size) {
		this->size = size;
	}

	void operator()(char *addr) {
		numa_free(addr, size);
	}
};

class aligned_deleter
{
public:
	void operator()(char *addr) {
		free(addr);
	}
};

std::shared_ptr<char> memalloc_node(int node_id, size_t num_bytes)
{
	if (num_bytes == 0)
		return std::shared_ptr<char>();

	std::shared_ptr<char> ret;
	if (node_id >= 0) {
		void *addr = numa_alloc_onnode(num_bytes, node_id);
		ret = std::shared_ptr<char>((char *) addr, NUMA_deleter(num_bytes));
	}
	else {
		ret = local_mem_buffer::alloc(num_bytes);
		if (ret == NULL)
			ret = std::shared_ptr<char>((char *) memalign(PAGE_SIZE, num_bytes),
					aligned_deleter());
	}
	assert(ret);
	return ret;
}

class reset_data_task: public thread_task
{
	char *raw_arr;
	size_t num_bytes;
public:
	reset_data_task(char *raw_arr, size_t num_bytes) {
		this->raw_arr = raw_arr;
		this->num_bytes = num_bytes;
	}

	void run() {
		memset(raw_arr, 0, num_bytes);
	}
};

class set_data_task: public thread_task
{
	char *to_buf;
	off_t to_off;
	size_t to_size;
	int node_id;
	const set_range_operate &set_range;
	size_t range_size;
public:
	// `to_off', `to_size' and `range_size' are in the number of bytes.
	set_data_task(char *to_buf, off_t to_off, size_t to_size, int node_id,
			const set_range_operate &_set_range, size_t range_size): set_range(
				_set_range) {
		this->to_buf = to_buf;
		this->to_off = to_off;
		this->to_size = to_size;
		this->node_id = node_id;
		this->range_size = range_size;
	}

	void run() {
		for (size_t rel_off = 0; rel_off < to_size; rel_off += range_size) {
			size_t size = std::min(to_size - rel_off, range_size);
			off_t off = to_off + rel_off;
			set_range.set(to_buf + off, size, off, node_id);
		}
	}
};

}

raw_data_array::raw_data_array(size_t num_bytes, int node_id)
{
	this->node_id = node_id;
	this->num_bytes = num_bytes;
	data = memalloc_node(node_id, num_bytes);
}

raw_data_array raw_data_array::deep_copy() const
{
	raw_data_array ret(*this);
	ret.data = memalloc_node(get_node_id(), num_bytes);
	memcpy(ret.data.get(), this->data.get(), num_bytes);
	return ret;
}

bool raw_data_array::copy_from(const raw_data_array &arr)
{
	if (num_bytes != arr.num_bytes) {
		BOOST_LOG_TRIVIAL(error) << "copy from an array of different length";
		return false;
	}
	memcpy(data.get(), arr.data.get(), num_bytes);
	return true;
}

void raw_data_array::expand(size_t min)
{
	size_t new_num_bytes = num_bytes;
	for (; new_num_bytes < min; new_num_bytes *= 2);
	std::shared_ptr<char> new_data;
	new_data = memalloc_node(get_node_id(), new_num_bytes);
	memcpy(new_data.get(), data.get(), num_bytes);
	num_bytes = new_num_bytes;
	data = new_data;
}

bool raw_data_array::set_sub_arr(off_t start, const char *buf, size_t size)
{
	if (start + size > get_num_bytes()) {
		BOOST_LOG_TRIVIAL(error) << "set_sub_arr: out of range";
		return false;
	}
	memcpy(data.get() + start, buf, size);
	return true;
}

void reset_arrays(std::vector<raw_data_array> &arrs)
{
	mem_thread_pool::ptr mem_threads
		= mem_thread_pool::get_global_mem_threads();
	for (size_t i = 0; i < arrs.size(); i++) {
		assert(arrs[i].get_node_id() >= 0);
		mem_threads->process_task(arrs[i].get_node_id(),
				new reset_data_task(arrs[i].get_raw(), arrs[i].get_num_bytes()));
	}
	mem_threads->wait4complete();
}

void set_array_ranges(const NUMA_mapper &mapper, size_t length,
		size_t entry_size, const set_range_operate &set_range,
		std::vector<raw_data_array> &arrs)
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
					new set_data_task(arrs[i].get_raw(), nbytes_per_thread * j,
						local_nbytes, arrs[i].get_node_id(), set_range,
						mapper.get_range_size() * entry_size));
		}
	}
	mem_threads->wait4complete();
}

}

}
