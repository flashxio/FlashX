#ifndef __MEM_VECTOR_VECTOR_H__
#define __MEM_VECTOR_VECTOR_H__

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

#include <vector>

#include "generic_type.h"
#include "vector_vector.h"
#include "raw_data_array.h"
#include "mem_vv_store.h"

namespace fm
{

class mem_vector_vector: public vector_vector
{
	detail::mem_vv_store::const_ptr store;

protected:
	mem_vector_vector(detail::mem_vv_store::const_ptr store): vector_vector(store) {
		this->store = store;
	}
public:
	typedef std::shared_ptr<mem_vector_vector> ptr;
	typedef std::shared_ptr<const mem_vector_vector> const_ptr;

	static ptr cast(vector_vector::ptr vec) {
		if (!vec->is_in_mem()) {
			BOOST_LOG_TRIVIAL(error) << "This vector vector isn't in memory";
			return ptr();
		}
		return std::static_pointer_cast<mem_vector_vector>(vec);
	}

	static ptr create(detail::mem_vv_store::const_ptr store) {
		return ptr(new mem_vector_vector(store));
	}

	static ptr create(const scalar_type &type) {
		return ptr(new mem_vector_vector(detail::mem_vv_store::create(type)));
	}

	static ptr create(const detail::raw_data_array &data,
			const std::vector<off_t> &offs, const scalar_type &type) {
		detail::mem_vv_store::ptr vec = detail::mem_vv_store::create(data,
				offs, type);
		return ptr(new mem_vector_vector(vec));
	}

	const char *get_raw_data() const {
		return store->get_raw_arr();
	}

	virtual size_t get_length() const {
		return store->get_num_vecs();
	}

	virtual size_t get_length(off_t idx) const {
		return store->get_length(idx);
	}

	virtual size_t get_tot_num_entries() const {
		return store->get_num_bytes() / get_type().get_size();
	}

	virtual const char*get_raw_arr(off_t idx) const {
		return store->get_raw_arr(idx);
	}

	virtual std::shared_ptr<vector> cat() const;
	virtual vector_vector::ptr groupby(const factor_vector &labels,
			const gr_apply_operate<sub_vector_vector> &op) const;
	virtual vector_vector::ptr apply(const arr_apply_operate &op) const;
	virtual vector_vector::ptr serial_apply(const arr_apply_operate &op) const;
};

}

#endif
