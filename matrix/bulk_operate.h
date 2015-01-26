#ifndef __BULK_OPERATE_H__
#define __BULK_OPERATE_H__

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
#include <vector>
#include <cmath>

namespace fm
{

class bulk_operate
{
public:
	/*
	 * This performs element-wise operation on two input arrays, and stores
	 * the result on the output array.
	 */
	virtual void runAA(size_t num_eles, const void *left_arr,
			const void *right_arr, void *output_arr) const = 0;
	/*
	 * This performs operations on the left input array and the right element,
	 * and stores the result on the output array.
	 */
	virtual void runAE(size_t num_eles, const void *left_arr,
			const void *right, void *output_arr) const = 0;

	/*
	 * This performs aggregation on the input array and stores the result
	 * on output.
	 */
	virtual void runA(size_t num_eles, const void *left_arr,
			void *output) const = 0;

	virtual size_t left_entry_size() const = 0;
	virtual size_t right_entry_size() const = 0;
	virtual size_t output_entry_size() const = 0;
};

template<class OpType, class LeftType, class RightType, class ResType>
class bulk_operate_impl: public bulk_operate
{
public:
	virtual void runAA(size_t num_eles, const void *left_arr1,
			const void *right_arr1, void *output_arr1) const {
		const LeftType *left_arr = (const LeftType *) left_arr1;
		const RightType *right_arr = (const RightType *) right_arr1;
		ResType *output_arr = (ResType *) output_arr1;
		OpType op;
		for (size_t i = 0; i < num_eles; i++)
			output_arr[i] = op(left_arr[i], right_arr[i]);
	}

	virtual void runAE(size_t num_eles, const void *left_arr1,
			const void *right, void *output_arr1) const {
		const LeftType *left_arr = (const LeftType *) left_arr1;
		ResType *output_arr = (ResType *) output_arr1;
		RightType entry = *(const RightType *) right;
		OpType op;
		for (size_t i = 0; i < num_eles; i++)
			output_arr[i] = op(left_arr[i], entry);
	}

	virtual void runA(size_t num_eles, const void *left_arr1,
			void *output) const {
		const LeftType *left_arr = (const LeftType *) left_arr1;
		ResType res = 0;
		OpType op;
		for (size_t i = 0; i < num_eles; i++)
			res = op(left_arr[i], res);
		*(ResType *) output = res;
	}

	virtual size_t left_entry_size() const {
		return sizeof(LeftType);
	}

	virtual size_t right_entry_size() const {
		return sizeof(RightType);
	}

	virtual size_t output_entry_size() const {
		return sizeof(ResType);
	}
};

class basic_ops
{
public:
	typedef std::shared_ptr<basic_ops> ptr;

	enum op_idx {
		ADD,
		SUB,
		MUL,
		DIV,
		MIN,
		MAX,
		POW,
		NUM_OPS,
	};

	virtual const bulk_operate *get_op(op_idx idx) const = 0;
	virtual const bulk_operate &get_add() const = 0;
	virtual const bulk_operate &get_sub() const = 0;
	virtual const bulk_operate &get_multiply() const = 0;
	virtual const bulk_operate &get_divide() const = 0;
};

template<class LeftType, class RightType, class ResType>
class basic_ops_impl: public basic_ops
{
	struct multiply {
		ResType operator()(const LeftType &e1, const RightType &e2) {
			return e1 * e2;
		}
	};

	struct add {
		ResType operator()(const LeftType &e1, const RightType &e2) {
			return e1 + e2;
		}
	};

	struct sub {
		ResType operator()(const LeftType &e1, const RightType &e2) {
			return e1 - e2;
		}
	};

	struct divide {
		ResType operator()(const LeftType &e1, const RightType &e2) {
			return e1 / e2;
		}
	};

	struct min {
		ResType operator()(const LeftType &e1, const RightType &e2) {
			if (e1 < e2)
				return e1;
			else
				return e2;
		}
	};

	struct max {
		ResType operator()(const LeftType &e1, const RightType &e2) {
			if (e1 > e2)
				return e1;
			else
				return e2;
		}
	};

	struct pow {
		ResType operator()(const LeftType &e1, const RightType &e2) {
			return std::pow(e1, e2);
		}
	};

	bulk_operate_impl<add, LeftType, RightType, ResType> add_op;
	bulk_operate_impl<sub, LeftType, RightType, ResType> sub_op;
	bulk_operate_impl<multiply, LeftType, RightType, ResType> mul_op;
	bulk_operate_impl<divide, LeftType, RightType, ResType> div_op;
	bulk_operate_impl<min, LeftType, RightType, ResType> min_op;
	bulk_operate_impl<max, LeftType, RightType, ResType> max_op;
	bulk_operate_impl<pow, LeftType, RightType, ResType> pow_op;

	std::vector<bulk_operate *> ops;
public:
	basic_ops_impl() {
		ops.resize(NUM_OPS);
		ops[ADD] = &add_op;
		ops[SUB] = &sub_op;
		ops[MUL] = &mul_op;
		ops[DIV] = &div_op;
		ops[MIN] = &min_op;
		ops[MAX] = &max_op;
		ops[POW] = &pow_op;
	}

	virtual const bulk_operate *get_op(op_idx idx) const {
		return ops[idx];
	}

	virtual const bulk_operate &get_add() const {
		return add_op;
	}

	virtual const bulk_operate &get_sub() const {
		return sub_op;
	}

	virtual const bulk_operate &get_multiply() const {
		return mul_op;
	}

	virtual const bulk_operate &get_divide() const {
		return div_op;
	}
};

/*
 * This operate is used to set values on a matrix.
 */
class set_operate
{
public:
	virtual void set(void *arr, size_t num_eles, off_t row_idx,
			off_t col_idx) const = 0;
	virtual size_t entry_size() const = 0;
};

template<class T>
class type_set_operate: public set_operate
{
	virtual void set(void *arr, size_t num_eles, off_t row_idx,
			off_t col_idx) const {
		set((T *) arr, num_eles, row_idx, col_idx);
	}
public:
	virtual void set(T *arr, size_t num_eles, off_t row_idx,
			off_t col_idx) const = 0;
	virtual size_t entry_size() const {
		return sizeof(T);
	}
};

}

#endif
