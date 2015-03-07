#ifndef __VECTOR_VECTOR_H__
#define __VECTOR_VECTOR_H__

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

#include <memory>

#include "comm_exception.h"
#include "log.h"

#include "vector.h"

namespace fm
{

class mem_vector;
class scalar_type;
class sub_vector_vector;
class factor_vector;

/*
 * This stores a vector of vectors. It's similar to the row-wise matrix,
 * but this allows each vector to have different lengths.
 */
class vector_vector: public vector
{
public:
	vector_vector(size_t length, bool in_mem): vector(length, 0, in_mem) {
	}

	typedef std::shared_ptr<vector_vector> ptr;

	static ptr cast(vector::ptr vec) {
		if (vec->get_entry_size() != 0) {
			BOOST_LOG_TRIVIAL(error) << "This isn't a vector of vectors";
			return ptr();
		}
		return std::static_pointer_cast<vector_vector>(vec);
	}

	size_t get_num_vecs() const {
		return vector::get_length();
	}

	virtual size_t get_tot_num_entries() const = 0;
	virtual size_t get_length(off_t idx) const = 0;
	virtual const char*get_raw_arr(off_t idx) const = 0;

	/*
	 * Catenate all vectors into a single vector.
	 */
	virtual std::shared_ptr<vector> cat() const = 0;

	virtual const scalar_type &get_type() const = 0;

	virtual bool resize(size_t new_length) {
		throw unsupported_exception("resize");
	}

	virtual bool set_sub_vec(off_t start, const vector &vec) {
		throw unsupported_exception("set_sub_vec");
	}

	virtual vector::const_ptr get_sub_vec(off_t start, size_t length) const {
		throw unsupported_exception("get_sub_vec");
	}
	virtual bool expose_sub_vec(off_t start, size_t length) {
		throw unsupported_exception("expose_sub_vec");
	}

	virtual bool append(std::vector<vector::ptr>::const_iterator vec_it,
			std::vector<vector::ptr>::const_iterator vec_end) = 0;
	// We can assume each vector can be stored in memory.
	virtual bool append(const vector &vec) = 0;

	virtual void sort() {
		throw unsupported_exception("sort");
	}
	virtual vector::ptr sort_with_index() {
		throw unsupported_exception("sort_with_index");
	}
	virtual bool is_sorted() const {
		throw unsupported_exception("is_sorted");
	}

	virtual std::shared_ptr<data_frame> groupby(
			const gr_apply_operate<mem_vector> &op, bool with_val) const {
		throw unsupported_exception("groupby");
	}
	virtual vector_vector::ptr groupby(const factor_vector &labels,
			const gr_apply_operate<sub_vector_vector> &op) const = 0;

	virtual vector::ptr deep_copy() const {
		// TODO
		throw unsupported_exception("deep_copy");
	}
	virtual vector::ptr shallow_copy() {
		// TODO
		throw unsupported_exception("shallow_copy");
	}
	virtual vector::const_ptr shallow_copy() const {
		// TODO
		throw unsupported_exception("shallow_copy");
	}
};

/*
 * This is used to access some vectors in vector_vector.
 * It is always in memory.
 */
class sub_vector_vector
{
	const vector_vector &vv;
	std::vector<off_t> vec_idxs;
public:
	sub_vector_vector(const vector_vector &_vv,
			const std::vector<off_t> &vec_idxs): vv(_vv) {
		this->vec_idxs = vec_idxs;
	}

	const scalar_type &get_type() const {
		return vv.get_type();
	}

	size_t get_num_vecs() const {
		return vec_idxs.size();
	}

	size_t get_length(off_t idx) const {
		return vv.get_length(vec_idxs[idx]);
	}

	const char *get_raw_arr(off_t idx) const {
		return vv.get_raw_arr(vec_idxs[idx]);
	}
};

}

#endif
