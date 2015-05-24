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

	virtual int get_num_nodes() const {
		return -1;
	}

	virtual char *get_raw_arr() = 0;
	virtual const char *get_raw_arr() const = 0;
	virtual const char *get_sub_arr(off_t start, off_t end) const = 0;
	virtual char *get_sub_arr(off_t start, off_t end) = 0;

	virtual std::shared_ptr<local_vec_store> get_portion(off_t loc,
			size_t size) = 0;
	virtual std::shared_ptr<const local_vec_store> get_portion(off_t loc,
			size_t size) const = 0;
};

class smp_vec_store: public mem_vec_store
{
	char *arr;
	detail::raw_data_array data;

	smp_vec_store(const detail::raw_data_array &data, const scalar_type &type);
	smp_vec_store(size_t length, const scalar_type &type);
public:
	typedef std::shared_ptr<smp_vec_store> ptr;
	typedef std::shared_ptr<const smp_vec_store> const_ptr;

	static ptr create(size_t length, const scalar_type &type) {
		return ptr(new smp_vec_store(length, type));
	}
	static ptr create(const detail::raw_data_array &data,
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

	virtual std::shared_ptr<local_vec_store> get_portion(off_t loc, size_t size);
	virtual std::shared_ptr<const local_vec_store> get_portion(off_t loc,
			size_t size) const;
	virtual size_t get_portion_size() const;

	virtual void reset_data() {
		data.reset_data();
	}

	void set_data(const set_operate &op);

	void set(const std::vector<const char *> &locs);

	virtual vec_store::ptr sort_with_index();
	virtual void sort() {
		get_type().get_sorter().sort(arr, get_length(), false);
	}
	virtual bool is_sorted() const {
		return get_type().get_sorter().is_sorted(get_raw_arr(),
				get_length(), false);
	}

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

template<class T>
class seq_set_operate: public type_set_operate<T>
{
	long n;
	T from;
	T by;
public:
	seq_set_operate(long n, T from, T by) {
		this->n = n;
		this->from = from;
		this->by = by;
	}

	virtual void set(T *arr, size_t num_eles, off_t row_idx,
			off_t col_idx) const {
		// We are initializing a single-column matrix.
		T v = from + row_idx * by;
		for (size_t i = 0; i < num_eles; i++) {
			arr[i] = v;
			v += by;
		}
	}
};

template<class EntryType>
class set_const_operate: public type_set_operate<EntryType>
{
	EntryType v;
public:
	set_const_operate(EntryType v) {
		this->v = v;
	}

	virtual void set(EntryType *arr, size_t num_eles, off_t row_idx,
			off_t col_idx) const {
		for (size_t i = 0; i < num_eles; i++)
			arr[i] = v;
	}
};

/*
 * Create a sequence of values in [start, end]. `end' is inclusive.
 */
template<class EntryType>
vec_store::ptr create_vec_store(EntryType start, EntryType end, EntryType stride)
{
	if ((end < start && stride > 0) || (end > stride && stride < 0)) {
		BOOST_LOG_TRIVIAL(error)
			<< "There are a negative number of elements in the sequence";
		return vec_store::ptr();
	}
	// TODO let's just use in-memory dense matrix first.
	long n = (end - start) / stride;
	// We need to count the start element.
	n++;
	detail::smp_vec_store::ptr v = detail::smp_vec_store::create(n,
			get_scalar_type<EntryType>());
	v->set_data(seq_set_operate<EntryType>(n, start, stride));
	return v;
}

/*
 * Create a vector filled with a constant value.
 */
template<class EntryType>
vec_store::ptr create_vec_store(size_t length, EntryType initv)
{
	// TODO let's just use in-memory dense matrix first.
	detail::smp_vec_store::ptr v = detail::smp_vec_store::create(length,
			get_scalar_type<EntryType>());
	v->set_data(set_const_operate<EntryType>(initv));
	return v;
}

template<>
vec_store::ptr create_vec_store<double>(double start, double end, double stride);

}

}

#endif
