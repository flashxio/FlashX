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

namespace fm
{

namespace detail
{

/*
 * This class holds the allocated memory and helps other classes access
 * the data in the memory. We shouldn't pass the object of this class
 * around, which is relatively expensive.
 */
class raw_data_array
{
	int node_id;
	// The total number of bytes in the allocated memory.
	size_t num_bytes;
	// This points to the beginning of the allocated memory.
	std::shared_ptr<char> data;
public:
	raw_data_array() {
		node_id = -1;
		num_bytes = 0;
	}

	raw_data_array(size_t num_bytes, int node_id);

	/*
	 * Get the number of bytes allocated.
	 */
	size_t get_num_bytes() const {
		return num_bytes;
	}

	void reset_data() {
		memset(data.get(), 0, num_bytes);
	}

	char *get_raw() {
		return data.get();
	}

	const char *get_raw() const {
		return data.get();
	}

	int get_node_id() const {
		return node_id;
	}

	raw_data_array deep_copy() const;

	bool copy_from(const raw_data_array &arr);
};

class NUMA_mapper;

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

/*
 * Reset all elements in the arrays.
 */
void reset_arrays(std::vector<raw_data_array> &arrs);
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
		std::vector<raw_data_array> &arrs);

}

}

#endif
