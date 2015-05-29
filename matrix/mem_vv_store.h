#ifndef __MEM_VV_STORE_H__
#define __MEM_VV_STORE_H__

/*
 * Copyright 2015 Open Connectome Project (http://openconnecto.me)
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

#include "mem_vec_store.h"

namespace fm
{

namespace detail
{

class matrix_store;

class mem_vv_store: public mem_vec_store
{
	// The offsets (in #bytes) of the vectors in the data array.
	// The last offset is the end of the last vector.
	std::vector<off_t> vec_offs;
	mem_vec_store::ptr store;

	mem_vv_store(const scalar_type &type);
	mem_vv_store(const detail::raw_data_array &data,
			const std::vector<off_t> &offs, const scalar_type &type);
public:
	typedef std::shared_ptr<mem_vv_store> ptr;
	typedef std::shared_ptr<const mem_vv_store> const_ptr;

	static ptr cast(vec_store::ptr store);

	static ptr create(const scalar_type &type) {
		return ptr(new mem_vv_store(type));
	}

	static ptr create(const detail::raw_data_array &data,
			const std::vector<off_t> &offs, const scalar_type &type) {
		return ptr(new mem_vv_store(data, offs, type));
	}

	static bool is_vector_vector(const vec_store &vec) {
		return vec.get_entry_size() == 0;
	}

	size_t get_num_bytes() const {
		return vec_offs.back() - vec_offs.front();
	}

	size_t get_num_bytes(off_t idx) const {
		return vec_offs[idx + 1] - vec_offs[idx];
	}

	size_t get_length(off_t idx) const {
		return (vec_offs[idx + 1] - vec_offs[idx]) / get_type().get_size();
	}

	size_t get_num_vecs() const {
		return vec_offs.size() - 1;
	}

	virtual char *get_raw_arr() {
		return store->get_raw_arr();
	}
	virtual const char *get_raw_arr() const {
		return store->get_raw_arr();
	}

	const char*get_raw_arr(off_t idx) const {
		return store->get_raw_arr() + vec_offs[idx];
	}

	virtual vec_store::const_ptr cat() const;

	virtual mem_vv_store::const_ptr get_sub_vec_vec(off_t start,
			size_t len) const;

	virtual vec_store::ptr deep_copy() const {
		mem_vv_store::ptr ret = mem_vv_store::ptr(new mem_vv_store(*this));
		ret->store = mem_vec_store::cast(this->store->deep_copy());
		return ret;
	}
	virtual vec_store::ptr shallow_copy() {
		// TODO the vector offsets may also be very large.
		return mem_vv_store::ptr(new mem_vv_store(*this));
	}
	virtual vec_store::const_ptr shallow_copy() const {
		// TODO the vector offsets may also be very large.
		return mem_vv_store::ptr(new mem_vv_store(*this));
	}

	bool append(const vec_store &vec);
	virtual bool append(std::vector<vec_store::const_ptr>::const_iterator vec_it,
			std::vector<vec_store::const_ptr>::const_iterator vec_end);

	/*
	 * The following methods aren't supported in the vv_store.
	 */

	virtual const char *get_sub_arr(off_t start, off_t end) const {
		assert(0);
		return NULL;
	}
	virtual char *get_sub_arr(off_t start, off_t end) {
		assert(0);
		return NULL;
	}
	virtual bool resize(size_t new_length) {
		assert(0);
		return false;
	}
	virtual bool expose_sub_vec(off_t start, size_t length) {
		assert(0);
		return false;
	}
	virtual size_t get_portion_size() const {
		assert(0);
		return -1;
	}
	virtual std::shared_ptr<local_vec_store> get_portion(off_t loc,
			size_t size) {
		assert(0);
		return std::shared_ptr<local_vec_store>();
	}
	virtual std::shared_ptr<const local_vec_store> get_portion(off_t loc,
			size_t size) const {
		return std::shared_ptr<local_vec_store>();
	}
	virtual vec_store::ptr sort_with_index() {
		assert(0);
		return vec_store::ptr();
	}
	virtual void sort() {
		assert(0);
	}
	virtual bool is_sorted() const {
		return false;
	}

	virtual void reset_data() {
		assert(0);
	}
	virtual void set_data(const set_operate &op) {
		assert(0);
	}

	virtual std::shared_ptr<const matrix_store> conv2mat(size_t nrow,
			size_t ncol, bool byrow) const {
		return std::shared_ptr<const matrix_store>();
	}
};

}

}

#endif
