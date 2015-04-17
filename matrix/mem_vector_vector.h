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

namespace fm
{

class mem_vector_vector: public vector_vector
{
	// The offsets (in #bytes) of the vectors in the data array.
	// The last offset is the end of the last vector.
	std::vector<off_t> vec_offs;

	detail::raw_data_array data;

	size_t get_num_bytes() const {
		return vec_offs.back();
	}

	char *get_end() {
		return data.get_raw() + get_num_bytes();
	}

protected:
	mem_vector_vector(const scalar_type &type): vector_vector(0, type,
			true), data(1024) {
		vec_offs.push_back(0);
	}

	mem_vector_vector(const detail::raw_data_array &data,
			const std::vector<off_t> &offs, const scalar_type &type): vector_vector(
				offs.size() - 1, type, true) {
		this->vec_offs = offs;
		this->data = data;
	}

	const char *get_raw_data() const {
		return data.get_raw();
	}

	size_t get_num_bytes(off_t idx) const {
		return vec_offs[idx + 1] - vec_offs[idx];
	}

	void append_mem_vector(const mem_vector &vec);
	void append_mem_vectors(std::vector<vector::ptr>::const_iterator vec_it,
			std::vector<vector::ptr>::const_iterator vec_end);
	void append_mem_vv(const mem_vector_vector &vv);
	void append_mem_vvs(std::vector<vector_vector::ptr>::const_iterator vec_it,
			std::vector<vector_vector::ptr>::const_iterator vec_end);
public:
	typedef std::shared_ptr<mem_vector_vector> ptr;

	static ptr cast(vector_vector::ptr vec) {
		if (!vec->is_in_mem()) {
			BOOST_LOG_TRIVIAL(error) << "This vector vector isn't in memory";
			return ptr();
		}
		return std::static_pointer_cast<mem_vector_vector>(vec);
	}

	static ptr create(const scalar_type &type) {
		return ptr(new mem_vector_vector(type));
	}

	static ptr create(const detail::raw_data_array &data,
			const std::vector<off_t> &offs, const scalar_type &type) {
		return ptr(new mem_vector_vector(data, offs, type));
	}

	virtual size_t get_tot_num_entries() const {
		return get_num_bytes() / get_type().get_size();
	}

	virtual size_t get_length(off_t idx) const {
		return (vec_offs[idx + 1] - vec_offs[idx]) / get_type().get_size();
	}
	virtual const char*get_raw_arr(off_t idx) const {
		return data.get_raw() + vec_offs[idx];
	}
	virtual mem_vector_vector::ptr get_sub_vec_vec(off_t start,
			size_t len) const;

	bool append(const vector &vec);
	virtual bool append(std::vector<vector::ptr>::const_iterator vec_it,
			std::vector<vector::ptr>::const_iterator vec_end);

	virtual std::shared_ptr<vector> cat() const;
	virtual vector_vector::ptr groupby(const factor_vector &labels,
			const gr_apply_operate<sub_vector_vector> &op) const;
	virtual vector_vector::ptr apply(const arr_apply_operate &op) const;
	virtual vector_vector::ptr serial_apply(const arr_apply_operate &op) const;
	virtual void reset_data() {
		data.reset_data();
	}
	virtual vector::ptr flatten() const;
};

}

#endif
