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
#include "data_frame.h"

namespace fm
{

class agg_gr_apply: public gr_apply_operate<mem_vector>
{
	const agg_operate &agg_op;
public:
	agg_gr_apply(const agg_operate &_agg_op): agg_op(_agg_op) {
	}

	virtual void run(const void *key, const mem_vector &val,
			mem_vector &vec) const {
		agg_op.run(val.get_length(), val.get_raw_arr(), vec.get_raw_arr());
	}
	virtual const scalar_type &get_key_type() const {
		return agg_op.get_input_type();
	}
	virtual const scalar_type &get_output_type() const {
		return agg_op.get_output_type();
	}
	virtual size_t get_num_out_eles() const {
		return 1;
	}
};

data_frame::ptr vector::groupby(const agg_operate &op, bool with_val) const
{
	// TODO this is a temporary solution for groupby with aggregation.
	// We can certainly implement it more efficiently.
	return groupby(agg_gr_apply(op), with_val);
}

}
