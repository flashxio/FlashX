#ifndef __MEM_VECTOR_H__
#define __MEM_VECTOR_H__

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

#include <boost/format.hpp>

#include "log.h"
#include "common.h"

#include "vector.h"
#include "mem_vec_store.h"

namespace fm
{

class data_frame;
class scalar_type;

class mem_vector: public vector
{
protected:
	mem_vector(detail::mem_vec_store::const_ptr store): vector(store) {
	}
	mem_vector(const mem_vector &vec): vector(vec) {
	}
public:
	typedef std::shared_ptr<mem_vector> ptr;
	typedef std::shared_ptr<const mem_vector> const_ptr;

	static ptr cast(vector::ptr vec);

	static ptr create(detail::mem_vec_store::const_ptr store) {
		return ptr(new mem_vector(store));
	}

	static ptr create(size_t length, const scalar_type &type,
			const set_operate &op);

	const char *get_raw_arr() const {
		return static_cast<const detail::mem_vec_store &>(
				get_data()).get_raw_arr();
	}

	const char *get(off_t idx) const {
		return get_raw_arr() + idx * get_data().get_entry_size();
	}

	template<class T>
	T get(off_t idx) const {
		return *(const T *) get(idx);
	}

	virtual bool equals(const mem_vector &vec) const;

	bool verify_groupby(const gr_apply_operate<local_vec_store> &op) const;
	using vector::groupby;
	virtual std::shared_ptr<data_frame> groupby(
			const gr_apply_operate<local_vec_store> &op,
			bool with_val) const;
	virtual scalar_variable::ptr aggregate(const bulk_operate &op) const;
	virtual scalar_variable::ptr dot_prod(const vector &vec) const;

	virtual vector::ptr sort() const;
	virtual std::shared_ptr<data_frame> sort_with_index() const;

	bool export2(FILE *f) const;
};

/*
 * Create a sequence of values in [start, end]. `end' is inclusive.
 */
template<class EntryType>
vector::ptr create_vector(EntryType start, EntryType end, EntryType stride)
{
	detail::vec_store::ptr store = detail::create_vec_store(start, end, stride);
	if (store == NULL)
		return vector::ptr();
	return mem_vector::create(detail::mem_vec_store::cast(store));
}

/*
 * Create a vector filled with a constant value.
 */
template<class EntryType>
vector::ptr create_vector(size_t length, EntryType initv)
{
	detail::vec_store::ptr store = detail::create_vec_store(length, initv);
	if (store == NULL)
		return vector::ptr();
	return mem_vector::create(detail::mem_vec_store::cast(store));
}

}

#endif
