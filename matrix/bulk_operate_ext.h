#ifndef __BULK_OPERATE_EXT_H__
#define __BULK_OPERATE_EXT_H__

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

#include <boost/format.hpp>

#include "comm_exception.h"
#include "bulk_operate.h"

namespace fm
{

/*
 * This is the interface for aggregating elements into a single element.
 */
class agg_operate
{
	bulk_operate::const_ptr agg;
	bulk_operate::const_ptr combine;

	agg_operate(bulk_operate::const_ptr agg, bulk_operate::const_ptr combine) {
		this->agg = agg;
		this->combine = combine;
	}
public:
	typedef std::shared_ptr<const agg_operate> const_ptr;

	static const_ptr create(bulk_operate::const_ptr agg);

	static const_ptr create(bulk_operate::const_ptr agg,
			bulk_operate::const_ptr combine);

	bool is_same() const {
		return agg == combine;
	}

	bool has_combine() const {
		return combine != NULL;
	}

	bulk_operate::const_ptr get_agg_ptr() const {
		return agg;
	}

	bulk_operate::const_ptr get_combine_ptr() const {
		return combine;
	}

	const bulk_operate &get_agg() const {
		return *agg;
	}

	const bulk_operate &get_combine() const {
		return *combine;
	}

	void runAgg(size_t num_eles, const void *left_arr, void *output) const {
		agg->runAgg(num_eles, left_arr, output);
	}

	/*
	 * This combines the partially aggregated result.
	 * The input and output types of this method are both the output type of
	 * the aggregation.
	 */
	void runCombine(size_t num_eles, const void *in, void *out) const {
		assert(combine);
		combine->runAgg(num_eles, in, out);
	}

	const scalar_type &get_input_type() const {
		return agg->get_left_type();
	}
	const scalar_type &get_output_type() const {
		return agg->get_output_type();
	}
};


class agg_ops
{
protected:
	std::vector<agg_operate::const_ptr> ops;
public:
	enum op_idx {
		COUNT,
		FIND_NEXT,
		FIND_PREV,
		ARGMIN,
		ARGMAX,
		MIN,
		MAX,
		SUM,
		PROD,
		AND,
		OR,
		NUM_OPS,
	};

	typedef std::shared_ptr<agg_ops> ptr;

	agg_ops() {
		ops.resize(NUM_OPS);
	}

	agg_operate::const_ptr get_count() const {
		return ops[COUNT];
	}
	agg_operate::const_ptr get_find_next() const {
		return ops[FIND_NEXT];
	}
	agg_operate::const_ptr get_find_prev() const {
		return ops[FIND_PREV];
	}

	agg_operate::const_ptr get_op(op_idx idx) const {
		if (idx >= NUM_OPS)
			return agg_operate::const_ptr();
		return ops[(int) idx];
	}
};

}

#endif
