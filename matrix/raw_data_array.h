#ifndef __RAW_DATA_ARRAY_H__
#define __RAW_DATA_ARRAY_H__

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

#include <string.h>
#include <assert.h>

#include <memory>
#include <vector>

#include "NUMA_mapper.h"

namespace fm
{

class set_operate;
class bulk_operate;
class scalar_variable;

namespace detail
{

/*
 * This class holds the allocated memory and helps other classes access
 * the data in the memory. We shouldn't pass the object of this class
 * around, which is relatively expensive.
 */
class raw_array
{
	int node_id;
	// The total number of bytes in the allocated memory.
	size_t num_bytes;
protected:
	void resize(size_t num_bytes) {
		this->num_bytes = num_bytes;
	}
public:
	raw_array() {
		node_id = -1;
		num_bytes = 0;
	}

	raw_array(size_t num_bytes, int node_id) {
		this->num_bytes = num_bytes;
		this->node_id = node_id;
	}

	bool is_empty() const {
		return num_bytes == 0;
	}

	/*
	 * Get the number of bytes allocated.
	 */
	size_t get_num_bytes() const {
		return num_bytes;
	}

	int get_node_id() const {
		return node_id;
	}
};

/*
 * This raw array is used as a local buffer for computation on a vector
 * or matrix.
 */
class local_raw_array: public raw_array
{
	// This points to the beginning of the allocated memory.
	std::shared_ptr<char> data;
public:
	local_raw_array() {
	}

	local_raw_array(size_t num_bytes);

	char *get_raw() {
		return data.get();
	}

	const char *get_raw() const {
		return data.get();
	}

	void reset_data() {
		memset(data.get(), 0, get_num_bytes());
	}

	void expand(size_t min);
};


/*
 * This raw array is used to hold the memory for a simple version of vector
 * and matrix implementations.
 */
class simple_raw_array: public raw_array
{
	// This points to the beginning of the allocated memory.
	std::shared_ptr<char> data;
public:
	simple_raw_array() {
	}

	simple_raw_array(size_t num_bytes, int node_id);

	char *get_raw() {
		return data.get();
	}

	const char *get_raw() const {
		return data.get();
	}

	std::shared_ptr<char> get_raw_data() {
		return data;
	}

	std::shared_ptr<const char> get_raw_data() const {
		return data;
	}

	void reset_data();

	simple_raw_array deep_copy() const;

	bool copy_from(const simple_raw_array &arr);
	/*
	 * Copy the data in the buffer to the start location.
	 * @start and @size are in the number of bytes.
	 */
	bool set_sub_arr(off_t start, const char *arr, size_t size);

	void expand(size_t min);
};

/*
 * This raw array holds data for large vectors and matrices in NUMA machines.
 * Data in this raw array isn't stored in a piece of contiguous memory.
 * Instead, the raw array is composed of many pieces of memory (memory chunks).
 * We can reuse memory chunks when they are free'd. In this way, we can reduce memory
 * allocation overhead because it's expensive to allocate a large piece memory
 * and populate the memory with pages.
 *
 * Even though data in this raw array isn't stored in contiguous memory, users
 * still want to access contiguous memory to some extent. As such, users need
 * to tell us their data access pattern. We allow users to access data aligned
 * with the pre-specified block size.
 *
 * In summary, the raw array is composed of many physical memory chunks and
 * each memory chunk is divided into multiple logical blocks.
 */
class chunked_raw_array: public raw_array
{
	// The pointer to the memory chunk, the actual number of bytes
	// in the memory chunk.
	typedef std::pair<std::shared_ptr<char>, size_t> mem_chunk_t;
	std::shared_ptr<std::vector<mem_chunk_t> > data;
	// The actual number of bytes stored in a chunk.
	size_t chunk_data_size;
	// The block size that data is guaranteed to be contiguous.
	size_t contig_block_size;
public:
	chunked_raw_array() {
		chunk_data_size = 0;
		contig_block_size = 0;
	}

	chunked_raw_array(size_t num_bytes, size_t block_size, int node_id);

	void reset_data();

	char *get_raw(off_t start) {
		off_t first_chunk = start / chunk_data_size;
		return data->at(first_chunk).first.get() + (start % chunk_data_size);
	}

	const char *get_raw(off_t start) const {
		off_t first_chunk = start / chunk_data_size;
		return data->at(first_chunk).first.get() + (start % chunk_data_size);
	}

	char *get_raw(off_t start, off_t end) {
		off_t first_chunk = start / chunk_data_size;
		off_t last_chunk = (end - 1) / chunk_data_size;
		// The data in the range isn't in the same memory chunk.
		if (first_chunk != last_chunk)
			return NULL;
		return data->at(first_chunk).first.get() + (start % chunk_data_size);
	}

	const char *get_raw(off_t start, off_t end) const {
		off_t first_chunk = start / chunk_data_size;
		off_t last_chunk = (end - 1) / chunk_data_size;
		// The data in the range isn't in the same memory chunk.
		if (first_chunk != last_chunk)
			return NULL;
		return data->at(first_chunk).first.get() + (start % chunk_data_size);
	}

	std::pair<const char *, size_t> get_chunk(size_t off) const {
		return std::pair<const char *, size_t>(data->at(off).first.get(),
				data->at(off).second);
	}
	std::pair<char *, size_t> get_chunk(size_t off) {
		return std::pair<char *, size_t>(data->at(off).first.get(),
				data->at(off).second);
	}
	size_t get_num_chunks() const {
		return data->size();
	}

	size_t get_contig_block_size() const {
		return contig_block_size;
	}

	chunked_raw_array deep_copy() const;
};

void init_memchunk_reserve(int num_nodes);
void destroy_memchunk_reserve();

}

}

#endif
