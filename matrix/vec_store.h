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

#include "generic_type.h"

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

	static ptr create(size_t length, const scalar_type &_type, bool in_mem);

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

	virtual bool resize(size_t new_length) {
		this->length = new_length;
		return true;
	}
	virtual bool append(std::vector<vec_store::const_ptr>::const_iterator vec_it,
			std::vector<vec_store::const_ptr>::const_iterator vec_end) = 0;
	virtual bool append(const vec_store &vec) = 0;
	virtual vec_store::ptr deep_copy() const = 0;
	virtual vec_store::ptr shallow_copy() = 0;
	virtual vec_store::const_ptr shallow_copy() const = 0;

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
			size_t ncol, bool byrow) const = 0;
};

}

}

#endif
