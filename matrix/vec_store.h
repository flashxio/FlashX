#ifndef __VEC_STORE_H__
#define __VEC_STORE_H__

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

#include "generic_type.h"
#include "bulk_operate.h"

namespace safs
{
class file_io_factory;
};

namespace fm
{

class local_vec_store;
class set_vec_operate;

namespace detail
{

class matrix_store;

class vec_store
{
	size_t length;
	const scalar_type &type;
	bool in_mem;
	// This is used to avoid virtual function call.
	int entry_size;
public:
	typedef std::shared_ptr<vec_store> ptr;
	typedef std::shared_ptr<const vec_store> const_ptr;

	static ptr create(size_t length, const scalar_type &_type, int num_nodes,
			bool in_mem);
	/*
	 * This creates a vector from an existing file or in-memory buffer.
	 * The I/O factory may point to a memory buffer or a file.
	 */
	static ptr create(std::shared_ptr<safs::file_io_factory> factory,
			const scalar_type &type);

	vec_store(size_t length, const scalar_type &_type, bool in_mem): type(_type) {
		this->length = length;
		this->in_mem = in_mem;
		this->entry_size = type.get_size();
	}
	/*
	 * This is used for vector vector because its entry size isn't fixed.
	 */
	vec_store(size_t length, size_t entry_size, const scalar_type &_type,
			bool in_mem): type(_type) {
		this->length = length;
		this->entry_size = entry_size;
		this->in_mem = in_mem;
	}

	virtual ~vec_store() {
	}

	bool is_in_mem() const {
		return in_mem;
	}

	size_t get_length() const {
		return length;
	}

	const scalar_type &get_type() const {
		return type;
	}

	template<class T>
	bool is_type() const {
		return get_type().get_type() == fm::get_type<T>();
	}

	size_t get_entry_size() const {
		return entry_size;
	}

	virtual int get_num_nodes() const {
		return -1;
	}

	virtual size_t get_num_bytes() const {
		return get_length() * get_type().get_size();
	}

	virtual bool resize(size_t new_length) {
		this->length = new_length;
		return true;
	}

	virtual void clear() {
		vec_store::resize(0);
	}

	// Copy #eles to the data array and return #eles copied.
	virtual size_t copy_to(char *data, size_t num_eles) const;

	/*
	 * This is almost the same as `resize' except that it doesn't change
	 * the length of the vector.
	 */
	virtual bool reserve(size_t num_eles) = 0;

	virtual size_t get_reserved_size() const = 0;

	/*
	 * These two functions need to be thread-safe.
	 */
	virtual bool append(std::vector<vec_store::const_ptr>::const_iterator vec_it,
			std::vector<vec_store::const_ptr>::const_iterator vec_end) = 0;
	virtual bool append(const vec_store &vec) = 0;

	virtual vec_store::ptr deep_copy() const = 0;
	virtual vec_store::ptr shallow_copy() = 0;
	virtual vec_store::const_ptr shallow_copy() const = 0;

	virtual bool set_portion(std::shared_ptr<const local_vec_store> store,
			off_t loc) = 0;
	virtual std::shared_ptr<local_vec_store> get_portion(off_t loc,
			size_t size) = 0;
	virtual std::shared_ptr<const local_vec_store> get_portion(off_t loc,
			size_t size) const = 0;
	virtual size_t get_portion_size() const = 0;
	virtual size_t get_num_portions() const {
		double len = get_length();
		return ceil(len / get_portion_size());
	}

	virtual void reset_data() = 0;
	virtual void set_data(const set_vec_operate &op) = 0;

	virtual vec_store::ptr sort_with_index() = 0;
	virtual void sort() = 0;
	virtual bool is_sorted() const = 0;
	virtual std::shared_ptr<const matrix_store> conv2mat(size_t nrow,
			size_t ncol, bool byrow) const;
	virtual std::shared_ptr<matrix_store> conv2mat(size_t nrow,
			size_t ncol, bool byrow) = 0;
};

template<class T>
class seq_set_vec_operate: public type_set_vec_operate<T>
{
	long n;
	T from;
	T by;
public:
	seq_set_vec_operate(long n, T from, T by) {
		this->n = n;
		this->from = from;
		this->by = by;
	}

	virtual void set(T *arr, size_t num_eles, off_t start_idx) const {
		// We are initializing a single-column matrix.
		T v = from + start_idx * by;
		for (size_t i = 0; i < num_eles; i++) {
			arr[i] = v;
			v += by;
		}
	}
};

/*
 * Create a sequence of values in [start, end]. `end' is inclusive.
 */
template<class EntryType>
vec_store::ptr create_seq_vec_store(EntryType start, EntryType end, EntryType stride,
		int num_nodes = -1, bool in_mem = true)
{
	if ((end < start && stride > 0) || (end > stride && stride < 0)) {
		BOOST_LOG_TRIVIAL(error)
			<< "There are a negative number of elements in the sequence";
		return vec_store::ptr();
	}
	long n = (end - start) / stride;
	// We need to count the start element.
	n++;
	detail::vec_store::ptr v = detail::vec_store::create(n,
			get_scalar_type<EntryType>(), num_nodes, in_mem);
	v->set_data(seq_set_vec_operate<EntryType>(n, start, stride));
	return v;
}

/*
 * Create a vector filled with a constant value.
 */
template<class EntryType>
vec_store::ptr create_rep_vec_store(size_t length, EntryType initv,
		int num_nodes = -1, bool in_mem = true)
{
	detail::vec_store::ptr v = detail::vec_store::create(length,
			get_scalar_type<EntryType>(), num_nodes, in_mem);
	v->set_data(const_set_vec_operate<EntryType>(initv));
	return v;
}

template<>
vec_store::ptr create_seq_vec_store<double>(double start, double end,
		double stride, int num_nodes, bool in_mem);

}

}

#endif
