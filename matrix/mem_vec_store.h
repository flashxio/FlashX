#ifndef __MEM_VEC_STORE_H__
#define __MEM_VEC_STORE_H__

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

#include "log.h"

#include "vec_store.h"
#include "raw_data_array.h"
#include "bulk_operate.h"

namespace fm
{

class local_vec_store;

namespace detail
{

class matrix_store;

class mem_vec_store: public vec_store
{
public:
	mem_vec_store(size_t length, const scalar_type &type): vec_store(
			length, type, true) {
	}
	mem_vec_store(size_t length, size_t entry_size,
			const scalar_type &type): vec_store(length, entry_size, type, true) {
	}

	typedef std::shared_ptr<mem_vec_store> ptr;
	typedef std::shared_ptr<const mem_vec_store> const_ptr;

	static ptr cast(vec_store::ptr store) {
		assert(store->is_in_mem());
		return std::static_pointer_cast<mem_vec_store>(store);
	}

	static const_ptr cast(vec_store::const_ptr store) {
		assert(store->is_in_mem());
		return std::static_pointer_cast<const mem_vec_store>(store);
	}

	static ptr create(size_t length, int num_nodes, const scalar_type &type);

	virtual int get_num_nodes() const {
		return -1;
	}

	virtual char *get_raw_arr() = 0;
	virtual const char *get_raw_arr() const = 0;
	virtual const char *get_sub_arr(off_t start, off_t end) const = 0;
	virtual char *get_sub_arr(off_t start, off_t end) = 0;
	virtual char *get(off_t idx) = 0;
	virtual const char *get(off_t idx) const = 0;

	virtual std::shared_ptr<local_vec_store> get_portion(off_t loc,
			size_t size) = 0;
	virtual std::shared_ptr<const local_vec_store> get_portion(off_t loc,
			size_t size) const = 0;
	virtual bool copy_from(const char *buf, size_t num_bytes);

	template<class T>
	T get(off_t idx) const {
		return *(const T *) get(idx);
	}

	template<class T>
	void set(off_t idx, T val) {
		*(T *) get(idx) = val;
	}
};

class smp_vec_store: public mem_vec_store
{
	char *arr;
	detail::simple_raw_array data;

	smp_vec_store(const detail::simple_raw_array &data, size_t length,
			const scalar_type &type);
	smp_vec_store(size_t length, const scalar_type &type);
public:
	typedef std::shared_ptr<smp_vec_store> ptr;
	typedef std::shared_ptr<const smp_vec_store> const_ptr;

	static ptr create(size_t length, const scalar_type &type) {
		return ptr(new smp_vec_store(length, type));
	}
	static ptr create(const detail::simple_raw_array &data,
			size_t length, const scalar_type &type);
	static ptr create(const detail::simple_raw_array &data,
			const scalar_type &type);

	static ptr cast(vec_store::ptr store) {
		assert(store->is_in_mem() && std::static_pointer_cast<mem_vec_store>(
					store)->get_num_nodes() < 0);
		return std::static_pointer_cast<smp_vec_store>(store);
	}

	static const_ptr cast(vec_store::const_ptr store) {
		assert(store->is_in_mem());
		return std::static_pointer_cast<const smp_vec_store>(store);
	}

	std::shared_ptr<char> get_raw_data() {
		return data.get_raw_data();
	}

	std::shared_ptr<const char> get_raw_data() const {
		return data.get_raw_data();
	}

	virtual char *get_raw_arr() {
		return arr;
	}
	virtual const char *get_raw_arr() const {
		return arr;
	}
	virtual const char *get_sub_arr(off_t start, off_t end) const {
		assert(start < end && (size_t) end <= get_length());
		return arr + start * get_type().get_size();
	}
	virtual char *get_sub_arr(off_t start, off_t end) {
		assert(start < end && (size_t) end <= get_length());
		return arr + start * get_type().get_size();
	}

	size_t get_sub_start() const {
		return (arr - data.get_raw()) / get_entry_size();
	}

	bool expose_sub_vec(off_t start, size_t length);

	virtual smp_vec_store::ptr get(const smp_vec_store &idxs) const;

	virtual bool reserve(size_t length);
	virtual size_t get_reserved_size() const;

	virtual bool append(std::vector<vec_store::const_ptr>::const_iterator vec_it,
			std::vector<vec_store::const_ptr>::const_iterator vec_end);
	virtual bool append(const vec_store &vec);
	virtual bool resize(size_t new_length);
	virtual vec_store::ptr deep_copy() const;
	virtual vec_store::ptr shallow_copy() {
		return vec_store::ptr(new smp_vec_store(*this));
	}
	virtual vec_store::const_ptr shallow_copy() const {
		return vec_store::ptr(new smp_vec_store(*this));
	}

	virtual bool set_portion(std::shared_ptr<const local_vec_store> store,
			off_t loc);
	virtual std::shared_ptr<local_vec_store> get_portion(off_t loc, size_t size);
	virtual std::shared_ptr<const local_vec_store> get_portion(off_t loc,
			size_t size) const;
	virtual size_t get_portion_size() const;

	virtual void reset_data() {
		data.reset_data();
	}

	void set_data(const set_vec_operate &op);

	void set(const std::vector<const char *> &locs);

	virtual vec_store::ptr sort_with_index();
	virtual void sort() {
		get_type().get_sorter().sort(arr, get_length(), false);
	}
	virtual bool is_sorted() const {
		return get_type().get_sorter().is_sorted(get_raw_arr(),
				get_length(), false);
	}

	virtual std::shared_ptr<matrix_store> conv2mat(size_t nrow,
			size_t ncol, bool byrow);

	char *get(off_t idx) {
		return arr + idx * get_entry_size();
	}

	const char *get(off_t idx) const {
		return arr + idx * get_entry_size();
	}

	template<class T>
	T get(off_t idx) const {
		return *(const T *) get(idx);
	}

	template<class T>
	void set(off_t idx, T val) {
		*(T *) get(idx) = val;
	}
};

}

}

#endif
