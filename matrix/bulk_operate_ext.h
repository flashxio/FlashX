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

	void runAgg(size_t num_eles, const void *left_arr, const void *orig,
			void *output) const {
		agg->runAgg(num_eles, left_arr, orig, output);
	}

	/*
	 * This combines the partially aggregated result.
	 * The input and output types of this method are both the output type of
	 * the aggregation.
	 */
	void runCombine(size_t num_eles, const void *in, const void *orig,
			void *out) const {
		assert(combine);
		combine->runAgg(num_eles, in, orig, out);
	}

	const scalar_type &get_input_type() const {
		return agg->get_left_type();
	}
	const scalar_type &get_output_type() const {
		return agg->get_output_type();
	}
};

/*
 * Find the first location where its element is different from the first
 * element in the array.
 */
template<class T>
class find_next_impl: public bulk_operate
{
public:
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
class find_prev_impl: public bulk_operate
{
public:
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
class count_operate: public bulk_operate
{
public:
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

	virtual void runAgg(size_t num_eles, const void *in, const void *orig,
			void *output) const {
		size_t *t_out = (size_t *) output;
		if (orig == NULL)
			t_out[0] = num_eles;
		else
			t_out[0] = (*(const size_t *) orig) + num_eles;
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
protected:
	agg_operate::const_ptr count;
	agg_operate::const_ptr find_next;
	agg_operate::const_ptr find_prev;
public:
	typedef std::shared_ptr<agg_ops> ptr;

	agg_operate::const_ptr get_count() const {
		return count;
	}
	agg_operate::const_ptr get_find_next() const {
		return find_next;
	}
	agg_operate::const_ptr get_find_prev() const {
		return find_prev;
	}
};

template<class InType, class OutType>
class agg_ops_impl: public agg_ops
{
public:
	agg_ops_impl() {
		bulk_operate::const_ptr count_agg
			= bulk_operate::const_ptr(new count_operate<InType>());
		count = agg_operate::create(count_agg,
				bulk_operate::conv2ptr(
					count_agg->get_output_type().get_basic_ops().get_add()));
		find_next = agg_operate::create(
				bulk_operate::const_ptr(new find_next_impl<InType>()),
				bulk_operate::const_ptr());
		find_prev = agg_operate::create(
				bulk_operate::const_ptr(new find_prev_impl<InType>()),
				bulk_operate::const_ptr());
	}
};

template<class T1, class T2>
class type_cast: public bulk_uoperate
{
public:
	virtual void runA(size_t num, const void *in, void *out) const {
		const T1 *t_in = (const T1 *) in;
		T2 *t_out = (T2 *) out;
		for (size_t i = 0; i < num; i++)
			t_out[i] = t_in[i];
	}
	virtual const scalar_type &get_input_type() const {
		return get_scalar_type<T1>();
	}
	virtual const scalar_type &get_output_type() const {
		return get_scalar_type<T2>();
	}
};

}

#endif
