#ifndef __LOCAL_VEC_STORE_H__
#define __LOCAL_VEC_STORE_H__

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

#include <string.h>

#include "mem_vec_store.h"
#include "bulk_operate.h"
#include "raw_data_array.h"

namespace fm
{

class data_frame;

namespace detail
{

class matrix_store;

}

class local_vec_store: public detail::mem_vec_store
{
	const off_t orig_global_start;
	const size_t orig_length;
	const char *const_data;
	char *data;
	off_t global_start;
	int node_id;
protected:
	void set_data(const char *const_data, char *data) {
		this->data = data;
		this->const_data = const_data;
	}

public:
	local_vec_store(const char *const_data, char *data, off_t _global_start,
			size_t length, const scalar_type &type,
			int node_id): detail::mem_vec_store(length, type), orig_global_start(
				_global_start), orig_length(length) {
		this->data = data;
		this->const_data = const_data;
		this->global_start = _global_start;
		this->node_id = node_id;
	}
	local_vec_store(const char *const_data, char *data, off_t _global_start,
			size_t length, size_t entry_size, const scalar_type &type,
			int node_id): detail::mem_vec_store(length, entry_size,
				type), orig_global_start(_global_start), orig_length(length) {
		this->data = data;
		this->const_data = const_data;
		this->global_start = _global_start;
		this->node_id = node_id;
	}

	virtual bool expose_sub_vec(off_t start, size_t length) {
		assert(start + length <= orig_length);
		global_start = orig_global_start + start;
		return mem_vec_store::resize(length);
	}
	virtual void reset_expose() {
		global_start = orig_global_start;
		mem_vec_store::resize(orig_length);
	}

	int get_node_id() const {
		return node_id;
	}

	off_t get_global_start() const {
		return global_start;
	}

	off_t get_local_start() const {
		return global_start - orig_global_start;
	}

	typedef std::shared_ptr<local_vec_store> ptr;
	typedef std::shared_ptr<const local_vec_store> const_ptr;

	virtual const char *get_raw_arr() const {
		return const_data;
	}
	virtual char *get_raw_arr() {
		return data;
	}
	virtual const char *get_sub_arr(off_t start, off_t end) const {
		assert(start < end && (size_t) end <= get_length());
		return const_data + start * get_type().get_size();
	}
	virtual char *get_sub_arr(off_t start, off_t end) {
		assert(start < end && (size_t) end <= get_length());
		return data + start * get_type().get_size();
	}

	std::shared_ptr<data_frame> groupby(
			const gr_apply_operate<local_vec_store> &op, bool with_val) const;

	virtual void reset_data() {
		memset(get_raw_arr(), 0, get_length() * get_entry_size());
	}

	virtual void set_data(const set_vec_operate &op) {
		// We always view a vector as a single-column matrix.
		op.set(get_raw_arr(), get_length(), global_start);
	}

	virtual detail::vec_store::ptr sort_with_index();
	virtual void sort();
	virtual bool is_sorted() const;
	virtual std::shared_ptr<detail::matrix_store> conv2mat(
			size_t nrow, size_t ncol, bool byrow);

	local_vec_store::ptr get(std::vector<off_t> &idxs) const;

	const char *get(off_t idx) const {
		return const_data + idx * get_entry_size();
	}

	char *get(off_t idx) {
		return data + idx * get_entry_size();
	}

	virtual local_vec_store::ptr get_portion(off_t loc, size_t size);
	virtual local_vec_store::const_ptr get_portion(off_t loc,
			size_t size) const;
	virtual bool set_portion(std::shared_ptr<const local_vec_store> store,
			off_t loc);

	void set_raw(off_t idx, const char *val) {
		assert(const_data == data);
		char *dst = data + idx * get_entry_size();
		memcpy(dst, val, get_entry_size());
	}

	template<class T>
	T get(off_t idx) const {
		return *(const T *) (const_data + idx * sizeof(T));
	}
	template<class T>
	void set(off_t idx, T val) {
		assert(const_data == data);
		*(T *) (data + idx * sizeof(T)) = val;
	}

	/*
	 * The following methods aren't needed by this class and its child classes.
	 */

	virtual size_t get_portion_size() const {
		assert(0);
		return 0;
	}
	virtual bool append(std::vector<detail::vec_store::const_ptr>::const_iterator vec_it,
			std::vector<detail::vec_store::const_ptr>::const_iterator vec_end) {
		assert(0);
		return false;
	}
	virtual bool append(const detail::vec_store &vec) {
		assert(0);
		return false;
	}
	virtual detail::vec_store::ptr deep_copy() const {
		assert(0);
		return detail::vec_store::ptr();
	}
};

class local_ref_vec_store: public local_vec_store
{
	char *orig_data;
public:
	local_ref_vec_store(char *data, off_t global_start, size_t length,
			const scalar_type &type, int node_id): local_vec_store(data,
				data, global_start, length, type, node_id) {
		this->orig_data = data;
	}
	virtual bool resize(size_t new_length);

	virtual detail::vec_store::ptr shallow_copy() {
		return detail::vec_store::ptr(new local_ref_vec_store(*this));
	}
	virtual detail::vec_store::const_ptr shallow_copy() const {
		return detail::vec_store::const_ptr(new local_ref_vec_store(*this));
	}
	virtual bool expose_sub_vec(off_t start, size_t length) {
		char *new_data = orig_data + start * get_type().get_size();
		set_data(new_data, new_data);
		return local_vec_store::expose_sub_vec(start, length);
	}
	virtual void reset_expose() {
		set_data(orig_data, orig_data);
		local_vec_store::reset_expose();
	}
};

class local_cref_vec_store: public local_vec_store
{
	const char *orig_data;
public:
	local_cref_vec_store(const char *data, off_t global_start, size_t length,
			const scalar_type &type, int node_id): local_vec_store(data,
				NULL, global_start, length, type, node_id) {
		this->orig_data = data;
	}
	virtual bool resize(size_t new_length);

	virtual detail::vec_store::ptr shallow_copy() {
		return detail::vec_store::ptr(new local_cref_vec_store(*this));
	}
	virtual detail::vec_store::const_ptr shallow_copy() const {
		return detail::vec_store::const_ptr(new local_cref_vec_store(*this));
	}
	virtual bool expose_sub_vec(off_t start, size_t length) {
		const char *new_data = orig_data + start * get_type().get_size();
		set_data(new_data, NULL);
		return local_vec_store::expose_sub_vec(start, length);
	}
	virtual void reset_expose() {
		set_data(orig_data, NULL);
		local_vec_store::reset_expose();
	}
};

class local_buf_vec_store: public local_vec_store
{
	detail::local_raw_array arr;
public:
	local_buf_vec_store(off_t global_start, size_t length, const scalar_type &type,
			// Let's not cache the memory used for local vectors first.
			// It's very frequent that we need to allocate memory for vectors
			// of different lengths. We need to deallocate the memory for these
			// vectors once they aren't used. Otherwise, we get memory leak.
			// TODO are there cases we actually need to keep allocated memory.
			int node_id, bool cached = false): local_vec_store(NULL, NULL,
				global_start, length, type, node_id),
			arr(length * type.get_size(), cached) {
		set_data(arr.get_raw(), arr.get_raw());
	}

	virtual bool resize(size_t new_length);

	virtual detail::vec_store::ptr shallow_copy() {
		return detail::vec_store::ptr(new local_buf_vec_store(*this));
	}
	virtual detail::vec_store::const_ptr shallow_copy() const {
		return detail::vec_store::const_ptr(new local_buf_vec_store(*this));
	}
	virtual bool expose_sub_vec(off_t start, size_t length) {
		char *new_data = arr.get_raw() + start * get_type().get_size();
		set_data(new_data, new_data);
		return local_vec_store::expose_sub_vec(start, length);
	}
	virtual void reset_expose() {
		set_data(arr.get_raw(), arr.get_raw());
		local_vec_store::reset_expose();
	}
};

}

#endif
