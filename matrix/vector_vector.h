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

#include "vec_store.h"
#include "vector.h"
#include "mem_vv_store.h"

namespace fm
{

class mem_vector;
class scalar_type;
class factor_vector;
class local_vv_store;
class data_frame;

/*
 * This stores a vector of vectors. It's similar to the row-wise matrix,
 * but this allows each vector to have different lengths.
 */
class vector_vector: public vector
{
	vector_vector(detail::mem_vv_store::const_ptr store): vector(store) {
	}
public:
	typedef std::shared_ptr<vector_vector> ptr;

	static bool is_vector_vector(const vector &vec) {
		return vec.get_entry_size() == 0;
	}

	static ptr cast(vector::ptr vec) {
		if (!is_vector_vector(*vec)) {
			BOOST_LOG_TRIVIAL(error) << "This isn't a vector of vectors";
			return ptr();
		}
		return std::static_pointer_cast<vector_vector>(vec);
	}

	static ptr create(detail::mem_vv_store::ptr store) {
		return ptr(new vector_vector(store));
	}

	static ptr create(const detail::raw_data_array &data,
			const std::vector<off_t> &offs, const scalar_type &type) {
		detail::mem_vv_store::ptr vec = detail::mem_vv_store::create(data,
				offs, type);
		return ptr(new vector_vector(vec));
	}

	virtual size_t get_entry_size() const {
		return 0;
	}

	size_t get_num_vecs() const {
		return get_length();
	}

	/*
	 * This return the number of vectors in the vector vector.
	 */
	virtual size_t get_length() const {
		return static_cast<const detail::mem_vv_store &>(
				get_data()).get_num_vecs();
	}

	virtual size_t get_length(off_t idx) const {
		return static_cast<const detail::mem_vv_store &>(
				get_data()).get_length(idx);
	}

	virtual size_t get_tot_num_entries() const {
		return static_cast<const detail::mem_vv_store &>(
				get_data()).get_num_bytes() / get_type().get_size();
	}

	/*
	 * Catenate all vectors into a single vector.
	 */
	virtual std::shared_ptr<vector> cat() const;

	virtual vector_vector::ptr groupby(const factor_vector &labels,
			const gr_apply_operate<local_vv_store> &op) const;
	virtual vector_vector::ptr apply(const arr_apply_operate &op) const;

	virtual vector::ptr sort() const {
		return vector::ptr();
	}
	virtual std::shared_ptr<data_frame> sort_with_index() const {
		return std::shared_ptr<data_frame>();
	}
	virtual std::shared_ptr<data_frame> groupby(
			const gr_apply_operate<local_vec_store> &op,
			bool with_val) const {
		return std::shared_ptr<data_frame>();
	}
	virtual scalar_variable::ptr aggregate(const bulk_operate &op) const {
		return scalar_variable::ptr();
	}
	virtual scalar_variable::ptr dot_prod(const vector &vec) const {
		return scalar_variable::ptr();
	}
	virtual std::shared_ptr<dense_matrix> conv2mat(size_t nrow, size_t ncol,
			bool byrow) const {
		return std::shared_ptr<dense_matrix>();
	}
};

}

#endif
