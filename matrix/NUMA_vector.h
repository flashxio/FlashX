#ifndef __NUMA_VECTOR_H__
#define __NUMA_VECTOR_H__

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
#include <assert.h>

#include <vector>

#include "bulk_operate.h"
#include "vector.h"
#include "matrix_config.h"
#include "raw_data_array.h"
#include "NUMA_mapper.h"

namespace fm
{

class scalar_type;

/*
 * This is another implementation of in-memory vector. In this implementation,
 * the data in the vector is split and stored on multiple NUMA nodes. All
 * operations on the vector are optimized accordingly.
 */
class NUMA_vector: public vector
{
	detail::NUMA_mapper mapper;

	std::vector<detail::raw_data_array> data;

	// The copy constructor performs shallow copy.
	NUMA_vector(const NUMA_vector &vec): vector(vec.get_length(),
			vec.get_type(), true), mapper(vec.mapper) {
		data = vec.data;
	}

	NUMA_vector(size_t length, size_t num_nodes, const scalar_type &type);
public:
	typedef std::shared_ptr<NUMA_vector> ptr;

	static ptr create(size_t length, const scalar_type &type) {
		return ptr(new NUMA_vector(length, matrix_conf.get_num_nodes(), type));
	}

	static ptr create(size_t length, size_t num_nodes, const scalar_type &type) {
		return ptr(new NUMA_vector(length, num_nodes, type));
	}

	static ptr cast(vector::ptr vec);

	virtual vector::ptr get_sub_vec(off_t start, size_t length) const;
	virtual bool expose_sub_vec(off_t start, size_t length);

	virtual bool append(std::vector<vector::ptr>::const_iterator vec_it,
			std::vector<vector::ptr>::const_iterator vec_end);
	virtual bool append(const vector &vec);

	virtual void sort();
	virtual vector::ptr sort_with_index();
	virtual bool is_sorted() const;

	// It should return data frame instead of vector.
	virtual std::shared_ptr<data_frame> groupby(
			const gr_apply_operate<mem_vector> &op, bool with_val) const;
	virtual bool dot_prod(const NUMA_vector &vec, scalar_variable &res) const;

	virtual vector::ptr deep_copy() const;
	virtual vector::ptr shallow_copy() const;

	virtual void reset_data();
	void set_data(const set_operate &op);

	bool is_sub_vec() const;

	/*
	 * Get the node id where the data specified by `off' is located.
	 */
	int get_node_id(off_t off) const {
		auto phy_loc = mapper.map2physical(off);
		return phy_loc.first;
	}

	/*
	 * Get a subarray in [start, end), which must be in the same range.
	 */
	const char *get_sub_arr(off_t start, off_t end) const;
	char *get_sub_arr(off_t start, off_t end);

	/*
	 * This copies a piece of contiguous memory to the NUMA vector.
	 */
	void copy_from(const char *buf, size_t num_bytes);
	bool copy_from(const NUMA_vector &vec);

	size_t get_num_nodes() const {
		return data.size();
	}

	char *get(off_t idx) {
		std::pair<int, size_t> loc = mapper.map2physical(idx);
		size_t off = loc.second * get_entry_size();
		return data[loc.first].get_raw() + off;
	}

	const char *get(off_t idx) const {
		std::pair<int, size_t> loc = mapper.map2physical(idx);
		size_t off = loc.second * get_entry_size();
		return data[loc.first].get_raw() + off;
	}

	template<class T>
	T get(off_t idx) const {
		return *(const T*) get(idx);
	}

	template<class T>
	void set(off_t idx, T v) {
		*(T*) get(idx) = v;
	}

	const detail::NUMA_mapper &get_mapper() const {
		return mapper;
	}

	size_t get_portion_size() const {
		return mapper.get_range_size();
	}
};

}

#endif
