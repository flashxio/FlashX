#ifndef __FM_BULK_OPERATE_IMPL_H__
#define __FM_BULK_OPERATE_IMPL_H__

/*
 * Copyright 2017 Open Connectome Project (http://openconnecto.me)
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

#include "bulk_operate.h"

namespace fm
{

/*
 * This template implements the interface of a bulk unary operator.
 */
template<class OpType, class InType, class OutType>
class bulk_uoperate_impl: public bulk_uoperate
{
	OpType op;
public:
	bulk_uoperate_impl() {
	}

	bulk_uoperate_impl(const OpType &_op): op(_op) {
	}

	virtual void runA(size_t num_eles, const void *in_arr1,
			void *output_arr1) const {
		const InType *in_arr = (const InType *) in_arr1;
		OutType *output_arr = (OutType *) output_arr1;
		for (size_t i = 0; i < num_eles; i++)
			output_arr[i] = op(in_arr[i]);
	}

	virtual const scalar_type &get_input_type() const {
		return get_scalar_type<InType>();
	}
	virtual const scalar_type &get_output_type() const {
		return get_scalar_type<OutType>();
	}
	virtual std::string get_name() const {
		return OpType::get_name();
	}
};

/*
 * This template implements the interface of the bulk binary operator.
 */
template<class OpType, class LeftType, class RightType, class ResType>
class bulk_operate_impl: public bulk_operate
{
	static const size_t BULK_LEN = 128;
	OpType op;

	void runAA_short(size_t num_eles, const LeftType *left_arr,
			const RightType *right_arr, ResType *output_arr) const {
		for (size_t i = 0; i < num_eles; i++)
			output_arr[i] = op(left_arr[i], right_arr[i]);
	}
	void runAA_fixed(const LeftType *left_arr,
			const RightType *right_arr, ResType *output_arr) const {
		for (size_t i = 0; i < BULK_LEN; i++)
			output_arr[i] = op(left_arr[i], right_arr[i]);
	}

	void runAE_short(size_t num_eles, const LeftType *left_arr,
			const RightType &entry, ResType *output_arr) const {
		for (size_t i = 0; i < num_eles; i++)
			output_arr[i] = op(left_arr[i], entry);
	}
	void runAE_fixed(const LeftType *left_arr, const RightType &entry,
			ResType *output_arr) const {
		for (size_t i = 0; i < BULK_LEN; i++)
			output_arr[i] = op(left_arr[i], entry);
	}

	void runEA_short(size_t num_eles, const LeftType &entry,
			const RightType *right_arr, ResType *output_arr) const {
		for (size_t i = 0; i < num_eles; i++)
			output_arr[i] = op(entry, right_arr[i]);
	}
	void runEA_fixed(const LeftType &entry, const RightType *right_arr,
			ResType *output_arr) const {
		for (size_t i = 0; i < BULK_LEN; i++)
			output_arr[i] = op(entry, right_arr[i]);
	}

	ResType runAgg_short(size_t num_eles, const LeftType *left_arr,
			const ResType &orig) const {
		ResType res = orig;
		for (size_t i = 0; i < num_eles; i++)
			res = op(left_arr[i], res);
		return res;
	}
	ResType runAgg_fixed(const LeftType *left_arr, const ResType &orig) const {
		// GCC cannot optimize reduction well for many operations.
		// Because some optimizations may change the result of the reduction.
		// For example, float-point operations aren't strictly commutative or
		// associative; even though min/max is commutative and associative
		// semantically, GCC cannot discover this property by default.
		// Arranging the operation as below helps GCC apply more optimizations
		// and vectorize operations.
		// NOTE: such arranging will change the result of float-point operations
		// slightly. We may not care because parallelization changes the result
		// anyway.
		const size_t tmp_len = 8;
		ResType tmp[tmp_len];
		for (size_t i = 0; i < tmp_len; i++)
			tmp[i] = OpType::get_agg_init();
		for (size_t i = 0; i < BULK_LEN; i += tmp_len) {
			for (size_t j = 0; j < tmp_len; j++)
				tmp[j] = op(left_arr[i + j], tmp[j]);
		}
		ResType res = orig;
		for (size_t i = 0; i < tmp_len; i++)
			res = op(tmp[i], res);
		return res;
	}
public:
	virtual void runAA(size_t num_eles, const void *left_arr1,
			const void *right_arr1, void *output_arr1) const {
		const LeftType *left_arr = (const LeftType *) left_arr1;
		const RightType *right_arr = (const RightType *) right_arr1;
		ResType *output_arr = (ResType *) output_arr1;
		if (num_eles < BULK_LEN)
			runAA_short(num_eles, left_arr, right_arr, output_arr);
		else {
			size_t num_eles_align = ROUND(num_eles, BULK_LEN);
			for (size_t i = 0; i < num_eles_align; i += BULK_LEN)
				runAA_fixed(left_arr + i, right_arr + i, output_arr + i);
			runAA_short(num_eles - num_eles_align, left_arr + num_eles_align,
					right_arr + num_eles_align, output_arr + num_eles_align);
		}
	}

	virtual void runAE(size_t num_eles, const void *left_arr1,
			const void *right, void *output_arr1) const {
		const LeftType *left_arr = (const LeftType *) left_arr1;
		ResType *output_arr = (ResType *) output_arr1;
		RightType entry = *(const RightType *) right;
		if (num_eles < BULK_LEN)
			runAE_short(num_eles, left_arr, entry, output_arr);
		else {
			size_t num_eles_align = ROUND(num_eles, BULK_LEN);
			for (size_t i = 0; i < num_eles_align; i += BULK_LEN)
				runAE_fixed(left_arr + i, entry, output_arr + i);
			runAE_short(num_eles - num_eles_align, left_arr + num_eles_align,
					entry, output_arr + num_eles_align);
		}
	}

	virtual void runEA(size_t num_eles, const void *left,
			const void *right_arr1, void *output_arr1) const {
		LeftType entry = *(const LeftType *) left;
		const RightType *right_arr = (const RightType *) right_arr1;
		ResType *output_arr = (ResType *) output_arr1;
		if (num_eles < BULK_LEN)
			runEA_short(num_eles, entry, right_arr, output_arr);
		else {
			size_t num_eles_align = ROUND(num_eles, BULK_LEN);
			for (size_t i = 0; i < num_eles_align; i += BULK_LEN)
				runEA_fixed(entry, right_arr + i, output_arr + i);
			runEA_short(num_eles - num_eles_align, entry,
					right_arr + num_eles_align, output_arr + num_eles_align);
		}
	}

	virtual void runAgg(size_t num_eles, const void *left_arr1,
			void *output) const {
		const LeftType *left_arr = (const LeftType *) left_arr1;
		ResType res = OpType::get_agg_init();
		if (num_eles < BULK_LEN)
			*(ResType *) output = runAgg_short(num_eles, left_arr, res);
		else {
			size_t num_eles_align = ROUND(num_eles, BULK_LEN);
			for (size_t i = 0; i < num_eles_align; i += BULK_LEN)
				res = runAgg_fixed(left_arr + i, res);
			*(ResType *) output = runAgg_short(num_eles - num_eles_align,
					left_arr + num_eles_align, res);
		}
	}

	virtual void runCum(size_t num_eles, const void *left_arr1,
			const void *prev1, void *output_arr1) const {
		const LeftType *left_arr = (const LeftType *) left_arr1;
		const RightType *prev = (const RightType *) prev1;
		ResType *output_arr = (ResType *) output_arr1;
		if (prev1)
			output_arr[0] = op(left_arr[0], *prev);
		else
			output_arr[0] = left_arr[0];
		for (size_t i = 1; i < num_eles; i++)
			output_arr[i] = op(left_arr[i], output_arr[i - 1]);
	}

	virtual const scalar_type &get_left_type() const {
		return get_scalar_type<LeftType>();
	}
	virtual const scalar_type &get_right_type() const {
		return get_scalar_type<RightType>();
	}
	virtual const scalar_type &get_output_type() const {
		return get_scalar_type<ResType>();
	}
	virtual std::string get_name() const {
		return OpType::get_name();
	}
};

}

#endif
