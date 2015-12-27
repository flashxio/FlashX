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

#include <boost/format.hpp>

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
	size_t new_num_bytes = get_num_bytes();
	for (; new_num_bytes < min; new_num_bytes *= 2);
	std::shared_ptr<char> new_data;
	// TODO should we allocate memory in the local buffer when we expand
	// memory.
	new_data = memalloc_node(-1, true, new_num_bytes);
	memcpy(new_data.get(), data.get(), get_num_bytes());
	resize(new_num_bytes);
	data = new_data;
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

namespace
{

/*
 * We need the memory chunk size to be at least a certain size.
 * We partition a matrix in one dimension (either horizontally or vertically).
 * We keep many rows or columns in the same range and map ranges to NUMA nodes.
 * As such, we need to store all data in a range to contiguous memory so that
 * we can access them easily. If a chunk size is too small, we can only allow
 * a tall dense matrix with very few columns.
 */
static const size_t ARR_CHUNK_SIZE = 64 * 1024 * 1024;

typedef std::vector<char *> chunk_set_t;

/*
 * To avoid the overhead of memory allocation, we keep some memory chunks that
 * have been populated with pages. These memory chunks were allocated previously.
 * Instead of deallocating them after they are used, we keep them to serve
 * next memory allocation requests. All memory chunks have the same size.
 *
 * smp_reserved_chunks maintains memory chunks without associating to any NUMA
 * node.
 * reserved_chunks maintains memory chunks for NUMA nodes. Each entry maintains
 * memory chunks in the corresponding NUMA node.
 * I would expect smp_reserved_chunks has very few memory chunks because
 * this scheme is mainly used for NUMA data containers.
 */
static std::vector<chunk_set_t> reserved_chunks;
static chunk_set_t smp_reserved_chunks;

static std::atomic<size_t> reserved_bytes;
static std::atomic<size_t> smp_reserved_bytes;

class smp_reserved_deleter
{
public:
	void operator()(char *addr) {
		smp_reserved_chunks.push_back(addr);
	}
};

class NUMA_reserved_deleter
{
	int node_id;
public:
	NUMA_reserved_deleter(int node_id) {
		this->node_id = node_id;
	}
	void operator()(char *addr) {
		reserved_chunks[node_id].push_back(addr);
	}
};

static std::shared_ptr<char> memchunk_alloc(int node_id, size_t num_bytes)
{
	if (num_bytes != ARR_CHUNK_SIZE)
		return memalloc_node(node_id, false, num_bytes);
	if (node_id < 0 && !smp_reserved_chunks.empty()) {
		auto ret = smp_reserved_chunks.back();
		smp_reserved_chunks.pop_back();
		return std::shared_ptr<char>(ret, smp_reserved_deleter());
	}
	else if (node_id < 0) {
		smp_reserved_bytes += num_bytes;
		return std::shared_ptr<char>(
				(char *) memalign(PAGE_SIZE, num_bytes), smp_reserved_deleter());
	}

	assert((size_t) node_id < reserved_chunks.size());
	if (!reserved_chunks[node_id].empty()) {
		auto ret = reserved_chunks[node_id].back();
		reserved_chunks[node_id].pop_back();
		return std::shared_ptr<char>(ret, NUMA_reserved_deleter(node_id));
	}
	else {
		reserved_bytes += num_bytes;
		void *addr = numa_alloc_onnode(num_bytes, node_id);
		return std::shared_ptr<char>((char *) addr, NUMA_reserved_deleter(node_id));
	}
}

}

void destroy_memchunk_reserve()
{
	BOOST_LOG_TRIVIAL(error) << boost::format(
			"reserved %ld bytes are associated with NUMA nodes and %ld bytes aren't\n")
		% reserved_bytes.load() % smp_reserved_bytes.load();
	for (size_t i = 0; i < smp_reserved_chunks.size(); i++)
		free(smp_reserved_chunks[i]);
	for (size_t i = 0; i < reserved_chunks.size(); i++)
		for (size_t j = 0; j < reserved_chunks[i].size(); j++)
			numa_free(reserved_chunks[i][j], ARR_CHUNK_SIZE);
	reserved_bytes = 0;
	smp_reserved_bytes = 0;
}

void init_memchunk_reserve(int num_nodes)
{
	destroy_memchunk_reserve();
	if (num_nodes > 0)
		reserved_chunks.resize(num_nodes);
}

chunked_raw_array::chunked_raw_array(size_t num_bytes, size_t block_size,
		int node_id): raw_array(num_bytes, node_id)
{
	this->contig_block_size = block_size;
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
		data->at(i).first = memchunk_alloc(node_id, ARR_CHUNK_SIZE);
		assert(remain_size > 0);
		data->at(i).second = std::min(chunk_data_size, remain_size);
		assert(data->at(i).second <= ARR_CHUNK_SIZE);
		remain_size -= data->at(i).second;
	}
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
