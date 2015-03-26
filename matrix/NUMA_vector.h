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
	// TODO We need to try different range size to get better performance.
	static const size_t range_size_log = 20;
	static const size_t range_size = 1 << range_size_log;
	static const size_t range_mask = range_size - 1;

	const size_t numa_log;
	const size_t numa_mask;

	std::vector<detail::raw_data_array> data;
	const scalar_type &type;

	/*
	 * This function is to map the linear index in the vector to the physical
	 * location of the data. This function is a little expensive if it's
	 * invoked for every element.
	 */
	std::pair<int, size_t> map2data(size_t idx) const {
		size_t tmp_idx = idx >> range_size_log;
		size_t range_idx = idx & range_mask;
		int node_id = tmp_idx & numa_mask;
		tmp_idx = tmp_idx >> numa_log;
		return std::pair<int, size_t>(node_id,
				(tmp_idx << range_size_log) + range_idx);
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

	virtual const scalar_type &get_type() const {
		return type;
	}

	virtual bool set_sub_vec(off_t start, const vector &vec) {
		// TODO
		assert(0);
	}

	virtual vector::const_ptr get_sub_vec(off_t start, size_t length) const {
		// TODO
		assert(0);
	}

	virtual bool expose_sub_vec(off_t start, size_t length) {
		// TODO
		assert(0);
	}

	virtual bool append(std::vector<vector::ptr>::const_iterator vec_it,
			std::vector<vector::ptr>::const_iterator vec_end) {
		// TODO
		assert(0);
	}

	virtual bool append(const vector &vec) {
		// TODO
		assert(0);
	}

	virtual void sort() {
		// TODO
		assert(0);
	}

	virtual vector::ptr sort_with_index() {
		// TODO
		assert(0);
	}

	virtual bool is_sorted() const {
		// TODO
		assert(0);
	}

	// It should return data frame instead of vector.
	virtual std::shared_ptr<data_frame> groupby(
			const gr_apply_operate<mem_vector> &op, bool with_val) const {
		// TODO
		assert(0);
	}

	virtual vector::ptr deep_copy() const {
		// TODO
		assert(0);
	}

	virtual vector::ptr shallow_copy() {
		// TODO
		assert(0);
	}

	virtual vector::const_ptr shallow_copy() const {
		// TODO
		assert(0);
	}

	virtual void reset_data();

	/*
	 * Get a subarray in [start, end), which must be in the same range.
	 */
	const char *get_sub_arr(off_t start, off_t end) const;
	char *get_sub_arr(off_t start, off_t end);

	size_t get_num_nodes() const {
		return data.size();
	}

	char *get(off_t idx) {
		std::pair<int, size_t> loc = map2data(idx);
		size_t off = loc.second * get_entry_size();
		return data[loc.first].get_raw() + off;
	}

	const char *get(off_t idx) const {
		std::pair<int, size_t> loc = map2data(idx);
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
};

}

#endif
