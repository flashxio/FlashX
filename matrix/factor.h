#ifndef __FACTOR_H__
#define __FACTOR_H__

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

#include "mem_vector.h"

/*
 * This factor is the same as the one in R.
 */

namespace fm
{

typedef int factor_value_t;

class factor
{
	size_t num_levels;
public:
	factor(size_t num_levels) {
		this->num_levels = num_levels;
	}

	size_t get_num_levels() const {
		return num_levels;
	}

	factor_value_t get_max() const {
		return num_levels - 1;
	}

	factor_value_t get_min() const {
		return 0;
	}

	bool is_valid_level(factor_value_t v) const {
		return v >= 0 && (size_t) v < num_levels;
	}
};

class factor_vector: public mem_vector
{
	factor f;

	factor_vector(const factor &_f, size_t len, const set_operate &op): mem_vector(
			detail::mem_vec_store::create(len,
				get_scalar_type<factor_value_t>())), f(_f) {
		const_cast<detail::vec_store &>(get_data()).set_data(op);
	}
public:
	typedef std::shared_ptr<factor_vector> ptr;

	static ptr create(const factor &f, size_t length, const set_operate &op) {
		return ptr(new factor_vector(f, length, op));
	}

	const factor &get_factor() const {
		return f;
	}
};

}

#endif
