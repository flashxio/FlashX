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

#include "comm_exception.h"
#include "bulk_operate.h"

namespace fm
{

/*
 * This is the interface for aggregating elements into a single element.
 */
class agg_operate: public bulk_operate
{
	struct empty_deleter {
		void operator()(const agg_operate *addr) {
		}
	};

	virtual void runAA(size_t num_eles, const void *left_arr,
			const void *right_arr, void *output_arr) const {
		throw unsupported_exception();
	}
	virtual void runAE(size_t num_eles, const void *left_arr,
			const void *right, void *output_arr) const {
		throw unsupported_exception();
	}
	virtual void runEA(size_t num_eles, const void *left,
			const void *right_arr, void *output_arr) const {
		throw unsupported_exception();
	}
public:
	typedef std::shared_ptr<const agg_operate> const_ptr;

	static const_ptr conv2ptr(const agg_operate &op) {
		return const_ptr(&op, empty_deleter());
	}

	/*
	 * This combines the partially aggregated result.
	 * The input and output types of this method are both the output type of
	 * the aggregation.
	 */
	virtual void runCombine(size_t num_eles, const void *in, const void *orig,
			void *out) const = 0;
	const scalar_type &get_input_type() const {
		return get_left_type();
	}
};

/*
 * Find the first location where its element is different from the first
 * element in the array.
 */
template<class T>
class find_next_impl: public agg_operate
{
public:
	/*
	 * The first element in the array is pointed by `left_arr' and there are
	 * `num_eles' in the array.
	 * If all elements are the same, it returns the number of elements
	 * in the array.
	 */
	virtual void runAgg(size_t num_eles, const void *left_arr, const void *orig,
			void *output) const {
		const T *curr = (const T *) left_arr;
		T val = *curr;
		size_t loc = 1;
		for (; loc < num_eles && curr[loc] == val; loc++);
		*(size_t *) output = loc;
	}
	virtual void runCombine(size_t num_eles, const void *in, const void *orig,
			void *out) const {
		throw unsupported_exception();
	}

	virtual const scalar_type &get_left_type() const {
		return get_scalar_type<T>();
	}
	virtual const scalar_type &get_right_type() const {
		return get_scalar_type<T>();
	}
	virtual const scalar_type &get_output_type() const {
		return get_scalar_type<size_t>();
	}
};

/*
 * Search backwards and find the first location where its element is different
 * from the last element in the array.
 */
template<class T>
class find_prev_impl: public agg_operate
{
public:
	/*
	 * The end of the array is indicated by `arr_end'. The last element
	 * is right before `arr_end'. There are `num_eles' elements in the array.
	 * If all elements are the same, it returns the number of elements
	 * in the array.
	 */
	virtual void runAgg(size_t num_eles, const void *arr_end, const void *orig,
			void *output) const {
		const T *curr = ((const T *) arr_end) - 1;
		T val = *curr;
		const T *first = ((const T *) arr_end) - num_eles;
		for (; curr > first && *curr == val; curr--);
		if (*curr != val)
			curr++;
		assert(*curr == val);
		*(size_t *) output = ((const T *) arr_end) - curr;
	}
	virtual void runCombine(size_t num_eles, const void *in, const void *orig,
			void *out) const {
		throw unsupported_exception();
	}

	virtual const scalar_type &get_left_type() const {
		return get_scalar_type<T>();
	}
	virtual const scalar_type &get_right_type() const {
		return get_scalar_type<T>();
	}
	virtual const scalar_type &get_output_type() const {
		return get_scalar_type<size_t>();
	}
};

template<class T>
class count_operate: public agg_operate
{
public:
	virtual void runAgg(size_t num_eles, const void *in, const void *orig,
			void *output) const {
		size_t *t_out = (size_t *) output;
		if (orig == NULL)
			t_out[0] = num_eles;
		else
			t_out[0] = (*(const size_t *) orig) + num_eles;
	}
	virtual void runCombine(size_t num_eles, const void *in, const void *orig,
			void *out) const {
		assert(num_eles > 0);
		const size_t *t_in = (const size_t *) in;
		size_t res;
		size_t i;
		if (orig == NULL) {
			i = 1;
			res = t_in[0];
		}
		else {
			i = 0;
			res = *(const size_t *) orig;
		}
		for (; i < num_eles; i++)
			res += t_in[i];
		*(size_t *) out = res;
	}

	virtual const scalar_type &get_left_type() const {
		return get_scalar_type<T>();
	}
	virtual const scalar_type &get_right_type() const {
		return get_scalar_type<T>();
	}
	virtual const scalar_type &get_output_type() const {
		return get_scalar_type<size_t>();
	}
};

class agg_ops
{
public:
	typedef std::shared_ptr<agg_ops> ptr;

	virtual const agg_operate &get_count() const = 0;
	virtual const agg_operate &get_find_next() const = 0;
	virtual const agg_operate &get_find_prev() const = 0;
};

template<class InType, class OutType>
class agg_ops_impl: public agg_ops
{
	count_operate<InType> count;
	find_next_impl<InType> find_next;
	find_prev_impl<InType> find_prev;
public:
	virtual const agg_operate &get_count() const {
		return count;
	}
	virtual const agg_operate &get_find_next() const {
		return find_next;
	}
	virtual const agg_operate &get_find_prev() const {
		return find_prev;
	}
};

}

#endif
