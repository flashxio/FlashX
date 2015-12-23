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

std::shared_ptr<char> memalloc_node(int node_id, bool is_local, size_t num_bytes)
{
	if (num_bytes == 0)
		return std::shared_ptr<char>();

	std::shared_ptr<char> ret;
	if (node_id >= 0) {
		void *addr = numa_alloc_onnode(num_bytes, node_id);
		ret = std::shared_ptr<char>((char *) addr, NUMA_deleter(num_bytes));
	}
	else {
		// We allocate memory for the local buffer differently to reduce
		// overhead.
		if (is_local)
			ret = local_mem_buffer::alloc(num_bytes);
		if (ret == NULL)
			ret = std::shared_ptr<char>((char *) memalign(PAGE_SIZE, num_bytes),
					aligned_deleter());
	}
	assert(ret);
	return ret;
}

}

local_raw_array::local_raw_array(size_t num_bytes): raw_array(num_bytes, -1)
{
	data = memalloc_node(-1, true, num_bytes);
}

void local_raw_array::expand(size_t min)
{
	BOOST_LOG_TRIVIAL(error) << "local raw array doesn't support expanding";
	assert(0);
}

simple_raw_array::simple_raw_array(size_t num_bytes,
		int node_id): raw_array(num_bytes, node_id)
{
	data = memalloc_node(node_id, false, num_bytes);
}

void simple_raw_array::reset_data()
{
	assert(0);
}

simple_raw_array simple_raw_array::deep_copy() const
{
	simple_raw_array ret(*this);
	ret.data = memalloc_node(get_node_id(), false, get_num_bytes());
	memcpy(ret.data.get(), this->data.get(), get_num_bytes());
	return ret;
}

bool simple_raw_array::copy_from(const simple_raw_array &arr)
{
	if (get_num_bytes() != arr.get_num_bytes()) {
		BOOST_LOG_TRIVIAL(error) << "copy from an array of different length";
		return false;
	}
	memcpy(data.get(), arr.data.get(), get_num_bytes());
	return true;
}

void simple_raw_array::expand(size_t min)
{
	size_t new_num_bytes = get_num_bytes();
	for (; new_num_bytes < min; new_num_bytes *= 2);
	std::shared_ptr<char> new_data;
	new_data = memalloc_node(get_node_id(), false, new_num_bytes);
	memcpy(new_data.get(), data.get(), get_num_bytes());
	resize(new_num_bytes);
	data = new_data;
}

bool simple_raw_array::set_sub_arr(off_t start, const char *buf, size_t size)
{
	if (start + size > get_num_bytes()) {
		BOOST_LOG_TRIVIAL(error) << "set_sub_arr: out of range";
		return false;
	}
	memcpy(data.get() + start, buf, size);
	return true;
}

/*
 * We need the memory chunk size to be at least a certain size.
 * We partition a matrix in one dimension (either horizontally or vertically).
 * We keep many rows or columns in the same range and map ranges to NUMA nodes.
 * As such, we need to store all data in a range to contiguous memory so that
 * we can access them easily. If a chunk size is too small, we can only allow
 * a tall dense matrix with very few columns.
 */
static const size_t ARR_CHUNK_SIZE = 64 * 1024 * 1024;

chunked_raw_array::chunked_raw_array(size_t num_bytes, size_t block_size,
		int node_id): raw_array(num_bytes, node_id)
{
	assert(block_size <= ARR_CHUNK_SIZE);
	// Total number of blocks to store data. We need to round it up.
	size_t num_blocks = ceil(((double) num_bytes) / block_size);
	// The number of blocks in a memory chunk. We need to round it down.
	size_t num_blocks_chunk = ARR_CHUNK_SIZE / block_size;
	// The number of memory chunks required to store blocks.
	size_t num_chunks = ceil(((double) num_blocks) / num_blocks_chunk);
	chunk_data_size = num_blocks_chunk * block_size;
	data = std::shared_ptr<std::vector<mem_chunk_t> >(
			new std::vector<mem_chunk_t>(num_chunks));
	size_t remain_size = num_bytes;
	for (size_t i = 0; i < num_chunks; i++) {
		data->at(i).first = memalloc_node(node_id, false, ARR_CHUNK_SIZE);
		assert(remain_size > 0);
		data->at(i).second = std::min(chunk_data_size, remain_size);
		assert(data->at(i).second <= ARR_CHUNK_SIZE);
		remain_size -= data->at(i).second;
	}
	printf("actual bytes: %ld, alloc mem: %ld\n", num_bytes,
			ARR_CHUNK_SIZE * num_chunks);
}

chunked_raw_array chunked_raw_array::deep_copy() const
{
	chunked_raw_array ret(*this);
	ret.data = std::shared_ptr<std::vector<mem_chunk_t> >(
			new std::vector<mem_chunk_t>(get_num_chunks()));
	for (size_t i = 0; i < ret.data->size(); i++) {
		ret.data->at(i).first = memalloc_node(get_node_id(), false,
				ARR_CHUNK_SIZE);
		ret.data->at(i).second = data->at(i).second;
		memcpy(ret.data->at(i).first.get(), data->at(i).first.get(),
				data->at(i).second);
	}
	return ret;
}

void chunked_raw_array::reset_data()
{
	for (size_t i = 0; i < data->size(); i++)
		memset(data->at(i).first.get(), 0, data->at(i).second);
}

}

}
