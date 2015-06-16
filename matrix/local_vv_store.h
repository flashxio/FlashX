#ifndef __LOCAL_VV_STORE_H__
#define __LOCAL_VV_STORE_H__

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

#include "vv_store.h"
#include "local_vec_store.h"

namespace fm
{

class local_vv_store: public local_vec_store
{
	std::vector<off_t> offs;
	local_vec_store::ptr vec;

	const char *get_orig_const_data() const {
		return static_cast<const local_vec_store &>(*vec).get_raw_arr();
	}

	char *get_orig_data() {
		return vec->get_raw_arr();
	}
public:
	typedef std::shared_ptr<local_vv_store> ptr;
	typedef std::shared_ptr<const local_vv_store> const_ptr;

	static ptr cast(local_vec_store::ptr vec) {
		assert(vec->get_entry_size() == 0);
		return std::static_pointer_cast<local_vv_store>(vec);
	}

	static const_ptr cast(local_vec_store::const_ptr vec) {
		assert(vec->get_entry_size() == 0);
		return std::static_pointer_cast<const local_vv_store>(vec);
	}

	local_vv_store(off_t global_start, const std::vector<off_t> &offs,
			local_vec_store::ptr vec): local_vec_store(
				static_cast<const local_vec_store &>(*vec).get_raw_arr(),
				vec->get_raw_arr(), global_start, offs.size() - 1, 0,
				vec->get_type(), -1) {
		assert(!detail::vv_store::is_vector_vector(*vec));
		this->offs = offs;
		this->vec = vec;
	}

	bool is_whole() const {
		return offs.front() == 0
			&& offs.size() - 1 == local_vec_store::get_length();
	}

	virtual bool expose_sub_vec(off_t start, size_t length) {
		assert(start + length < offs.size());
		char *new_data = this->get_orig_data() + offs[start];
		const char *const_new_data = this->get_orig_const_data() + offs[start];
		local_vec_store::set_data(const_new_data, new_data);
		return local_vec_store::expose_sub_vec(start, length);
	}
	virtual void reset_expose() {
		local_vec_store::set_data(this->get_orig_const_data(),
				this->get_orig_data());
		local_vec_store::reset_expose();
	}

	virtual const char *get_sub_arr(off_t start, off_t end) const {
		start += get_local_start();
		assert((size_t) start < offs.size() && (size_t) end < offs.size());
		return this->get_orig_const_data() + offs[start];
	}
	virtual char *get_sub_arr(off_t start, off_t end) {
		start += get_local_start();
		assert((size_t) start < offs.size() && (size_t) end < offs.size());
		return this->get_orig_data() + offs[start];
	}

	virtual void reset_data() {
		assert(is_whole());
		memset(this->get_orig_data(), 0, get_num_bytes());
	}

	size_t get_num_bytes(off_t idx) const {
		idx += get_local_start();
		return offs[idx + 1] - offs[idx];
	}

	using local_vec_store::get_length;
	size_t get_length(off_t idx) const {
		return get_num_bytes(idx) / get_type().get_size();
	}

	size_t get_num_vecs() const {
		return local_vec_store::get_length();
	}

	size_t get_num_bytes() const {
		off_t start = get_local_start();
		off_t end = start + local_vec_store::get_length();
		return offs[end] - offs[start];
	}

	using local_vec_store::get_raw_arr;
	const char*get_raw_arr(off_t idx) const {
		idx += get_local_start();
		return this->get_orig_const_data() + offs[idx];
	}

	virtual bool resize(size_t new_length) {
		assert(0);
		return false;
	}

	virtual detail::vec_store::ptr shallow_copy() {
		return ptr(new local_vv_store(*this));
	}
	virtual detail::vec_store::const_ptr shallow_copy() const {
		return ptr(new local_vv_store(*this));
	}

	virtual local_vec_store::ptr get_portion(off_t loc, size_t size);
	virtual local_vec_store::const_ptr get_portion(off_t loc,
			size_t size) const;

	/*
	 * The methods below aren't used.
	 */

	virtual void set_data(const set_vec_operate &op) {
		assert(0);
	}

	virtual detail::vec_store::ptr sort_with_index() {
		assert(0);
	}
	virtual void sort() {
		assert(0);
	}
	virtual bool is_sorted() const {
		return false;
	}
	virtual std::shared_ptr<const detail::matrix_store> conv2mat(
			size_t nrow, size_t ncol, bool byrow) const {
		assert(0);
	}

	local_vec_store::ptr get(std::vector<off_t> &idxs) const;
};

namespace detail
{
class mem_vv_store;

std::shared_ptr<detail::mem_vv_store> apply(const local_vv_store &store,
		const arr_apply_operate &op);
}

}

#endif
