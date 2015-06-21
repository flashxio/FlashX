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

#include <assert.h>

#include <memory>
#include <vector>
#include <cmath>

#include "generic_type.h"

namespace fm
{

/**
 * This is a bulk version of a unary operator that takes one input and
 * generates an output. The bulk version is to amortize the overhead of
 * invoking a virtual method.
 */
class bulk_uoperate
{
	struct empty_deleter {
		void operator()(const bulk_uoperate *addr) {
		}
	};
public:
	typedef std::shared_ptr<const bulk_uoperate> const_ptr;

	static const_ptr conv2ptr(const bulk_uoperate &op) {
		return const_ptr(&op, empty_deleter());
	}

	virtual void runA(size_t num_eles, const void *in_arr,
			void *out_arr) const = 0;
	virtual const scalar_type &get_input_type() const = 0;
	virtual const scalar_type &get_output_type() const = 0;

	size_t input_entry_size() const {
		return get_input_type().get_size();
	}

	size_t output_entry_size() const {
		return get_output_type().get_size();
	}
};

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

	virtual const scalar_type &get_input_type() const;
	virtual const scalar_type &get_output_type() const;
};

/**
 * This is a bulk version of a binary operator that takes two inputs and
 * generates an output. The bulk version is to amortize the overhead of
 * invoking a virtual method.
 * This class defines the interface of a binary operator. A subclass has to
 * implement three forms of the operator.
 */
class bulk_operate
{
	struct empty_deleter {
		void operator()(const bulk_operate *addr) {
		}
	};
public:
	typedef std::shared_ptr<const bulk_operate> const_ptr;

	static const_ptr conv2ptr(const bulk_operate &op) {
		return const_ptr(&op, empty_deleter());
	}

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
	 * This performs operations on the left element array and the right array,
	 * and stores the result on the output array.
	 */
	virtual void runEA(size_t num_eles, const void *left,
			const void *right_arr, void *output_arr) const = 0;

	/*
	 * This performs aggregation on the input array and stores the result
	 * on output.
	 */
	virtual void runA(size_t num_eles, const void *left_arr,
			void *output) const = 0;

	virtual const scalar_type &get_left_type() const = 0;
	virtual const scalar_type &get_right_type() const = 0;
	virtual const scalar_type &get_output_type() const = 0;

	size_t left_entry_size() const {
		return get_left_type().get_size();
	}

	size_t right_entry_size() const {
		return get_right_type().get_size();
	}

	size_t output_entry_size() const {
		return get_output_type().get_size();
	}
};

/*
 * This template implements the interface of the bulk binary operator.
 */
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

	virtual void runEA(size_t num_eles, const void *left,
			const void *right_arr1, void *output_arr1) const {
		LeftType entry = *(const LeftType *) left;
		const RightType *right_arr = (const RightType *) right_arr1;
		ResType *output_arr = (ResType *) output_arr1;
		OpType op;
		for (size_t i = 0; i < num_eles; i++)
			output_arr[i] = op(entry, right_arr[i]);
	}

	virtual void runA(size_t num_eles, const void *left_arr1,
			void *output) const {
		const LeftType *left_arr = (const LeftType *) left_arr1;
		if (num_eles == 0)
			return;

		ResType res = left_arr[0];
		OpType op;
		for (size_t i = 1; i < num_eles; i++)
			res = op(left_arr[i], res);
		*(ResType *) output = res;
	}

	virtual const scalar_type &get_left_type() const;
	virtual const scalar_type &get_right_type() const;
	virtual const scalar_type &get_output_type() const;
};

/*
 * This is the interface for aggregating elements into a single element.
 */
class agg_operate
{
public:
	virtual void run(size_t num_eles, const void *in, void *output) const = 0;
	virtual const scalar_type &get_input_type() const = 0;
	virtual const scalar_type &get_output_type() const = 0;

	size_t input_entry_size() const {
		return get_input_type().get_size();
	}

	size_t output_entry_size() const {
		return get_output_type().get_size();
	}
};

/*
 * This interface defines a collection of basic unary operators.
 */
class basic_uops
{
public:
	typedef std::shared_ptr<basic_uops> ptr;

	enum op_idx {
		NEG,
		SQRT,
		ABS,
		NOT,
		SQ,
		NUM_OPS,
	};

	virtual const bulk_uoperate *get_op(op_idx idx) const = 0;
};

/*
 * This template implements all basic binary operators for different types.
 */
template<class InType, class OutType>
class basic_uops_impl: public basic_uops
{
	struct uop_neg {
		OutType operator()(const InType &e) const {
			return -e;
		}
	};

	struct uop_sqrt {
		double operator()(const InType &e) const {
			return std::sqrt(e);
		}
	};

	struct uop_abs {
		OutType operator()(const InType &e) const {
			return std::abs(e);
		}
	};

	struct uop_not {
		bool operator()(const bool &e) const {
			return !e;
		}
	};

	struct sq {
		OutType operator()(const InType &e) const {
			return e * e;
		}
	};

	bulk_uoperate_impl<uop_neg, InType, OutType> neg_op;
	bulk_uoperate_impl<uop_sqrt, InType, double> sqrt_op;
	bulk_uoperate_impl<uop_abs, InType, OutType> abs_op;
	bulk_uoperate_impl<uop_not, bool, bool> not_op;
	bulk_uoperate_impl<sq, InType, OutType> sq_op;

	std::vector<bulk_uoperate *> ops;
public:
	basic_uops_impl() {
		ops.push_back(&neg_op);
		ops.push_back(&sqrt_op);
		ops.push_back(&abs_op);
		ops.push_back(&not_op);
		ops.push_back(&sq_op);
	}

	virtual const bulk_uoperate *get_op(op_idx idx) const {
		if (idx >= op_idx::NUM_OPS)
			return NULL;
		return ops[idx];
	}
};

/*
 * This interface defines a collection of basic binary operators.
 */
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
		EQ,
		GT,
		GE,
		NUM_OPS,
	};

	virtual const bulk_operate *get_op(op_idx idx) const = 0;
	virtual const bulk_operate &get_add() const {
		return *get_op(ADD);
	}

	virtual const bulk_operate &get_sub() const {
		return *get_op(SUB);
	}

	virtual const bulk_operate &get_multiply() const {
		return *get_op(MUL);
	}

	virtual const bulk_operate &get_divide() const {
		return *get_op(DIV);
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
	virtual void run(size_t num_eles, const void *left_arr,
			void *output) const {
		const T *curr = (const T *) left_arr;
		T val = *curr;
		size_t loc = 1;
		for (; loc < num_eles && curr[loc] == val; loc++);
		*(size_t *) output = loc;
	}

	virtual const scalar_type &get_input_type() const;
	virtual const scalar_type &get_output_type() const;
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
	virtual void run(size_t num_eles, const void *arr_end,
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

	virtual const scalar_type &get_input_type() const;
	virtual const scalar_type &get_output_type() const;
};

class agg_ops
{
public:
	typedef std::shared_ptr<agg_ops> ptr;

	virtual const agg_operate &get_find_next() const = 0;
	virtual const agg_operate &get_find_prev() const = 0;
};

template<class InType, class OutType>
class agg_ops_impl: public agg_ops
{
	find_next_impl<InType> find_next;
	find_prev_impl<InType> find_prev;
public:
	virtual const agg_operate &get_find_next() const {
		return find_next;
	}
	virtual const agg_operate &get_find_prev() const {
		return find_prev;
	}
};

/*
 * This template implements all basic binary operators for different types.
 */
template<class LeftType, class RightType, class ResType>
class basic_ops_impl: public basic_ops
{
	struct multiply {
		ResType operator()(const LeftType &e1, const RightType &e2) const {
			return e1 * e2;
		}
	};

	struct add {
		ResType operator()(const LeftType &e1, const RightType &e2) const {
			return e1 + e2;
		}
	};

	struct sub {
		ResType operator()(const LeftType &e1, const RightType &e2) const {
			return e1 - e2;
		}
	};

	// Division is special. Its output should be float point.
	// Therefore, we convert both input values to float point.
	struct divide {
		double operator()(const LeftType &e1, const RightType &e2) const {
			double d1 = e1;
			double d2 = e2;
			return d1 / d2;
		}
	};

	struct min {
		ResType operator()(const LeftType &e1, const RightType &e2) const {
			if (e1 < e2)
				return e1;
			else
				return e2;
		}
	};

	struct max {
		ResType operator()(const LeftType &e1, const RightType &e2) const {
			if (e1 > e2)
				return e1;
			else
				return e2;
		}
	};

	struct pow {
		ResType operator()(const LeftType &e1, const RightType &e2) const {
			return std::pow(e1, e2);
		}
	};

	struct eq {
		bool operator()(const LeftType &e1, const RightType &e2) const {
			return e1 == e2;
		}
	};

	struct gt {
		bool operator()(const LeftType &e1, const RightType &e2) const {
			return e1 > e2;
		}
	};

	struct ge {
		bool operator()(const LeftType &e1, const RightType &e2) const {
			return e1 >= e2;
		}
	};

	bulk_operate_impl<add, LeftType, RightType, ResType> add_op;
	bulk_operate_impl<sub, LeftType, RightType, ResType> sub_op;
	bulk_operate_impl<multiply, LeftType, RightType, ResType> mul_op;
	bulk_operate_impl<divide, LeftType, RightType, double> div_op;
	bulk_operate_impl<min, LeftType, RightType, ResType> min_op;
	bulk_operate_impl<max, LeftType, RightType, ResType> max_op;
	bulk_operate_impl<pow, LeftType, RightType, ResType> pow_op;
	bulk_operate_impl<eq, LeftType, RightType, bool> eq_op;
	bulk_operate_impl<gt, LeftType, RightType, bool> gt_op;
	bulk_operate_impl<ge, LeftType, RightType, bool> ge_op;

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
		ops[EQ] = &eq_op;
		ops[GT] = &gt_op;
		ops[GE] = &ge_op;
	}

	virtual const bulk_operate *get_op(op_idx idx) const {
		if (idx >= op_idx::NUM_OPS)
			return NULL;
		return ops[idx];
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
	virtual const scalar_type &get_type() const = 0;
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

	virtual const scalar_type &get_type() const {
		return get_scalar_type<T>();
	}
};

template<class T>
class const_set_operate: public type_set_operate<T>
{
	T val;
public:
	const_set_operate(T val) {
		this->val = val;
	}

	virtual void set(T *arr, size_t num_eles, off_t row_idx,
			off_t col_idx) const {
		for (size_t i = 0; i < num_eles; i++)
			arr[i] = val;
	}
};

/*
 * This operate set values in a vector.
 */

class set_vec_operate
{
public:
	typedef std::shared_ptr<const set_vec_operate> const_ptr;

	virtual void set(void *arr, size_t num_eles, off_t start_idx) const = 0;
	virtual const scalar_type &get_type() const = 0;
};

template<class T>
class type_set_vec_operate: public set_vec_operate
{
	virtual void set(void *arr, size_t num_eles, off_t start_idx) const {
		set((T *) arr, num_eles, start_idx);
	}
public:
	virtual void set(T *arr, size_t num_eles, off_t start_idx) const = 0;
	virtual const scalar_type &get_type() const {
		return get_scalar_type<T>();
	}
};

template<class T>
class const_set_vec_operate: public type_set_vec_operate<T>
{
	T val;
public:
	const_set_vec_operate(T val) {
		this->val = val;
	}
	virtual void set(T *arr, size_t num_eles, off_t start_idx) const {
		for (size_t i = 0; i < num_eles; i++)
			arr[i] = val;
	}
};

template<class OpType, class InType, class OutType>
const scalar_type &bulk_uoperate_impl<OpType, InType, OutType>::get_input_type() const
{
	return get_scalar_type<InType>();
}

template<class OpType, class InType, class OutType>
const scalar_type &bulk_uoperate_impl<OpType, InType, OutType>::get_output_type() const
{
	return get_scalar_type<OutType>();
}

template<class OpType, class LeftType, class RightType, class ResType>
const scalar_type &bulk_operate_impl<OpType, LeftType, RightType, ResType>::get_left_type() const
{
	return get_scalar_type<LeftType>();
}

template<class OpType, class LeftType, class RightType, class ResType>
const scalar_type &bulk_operate_impl<OpType, LeftType, RightType, ResType>::get_right_type() const
{
	return get_scalar_type<RightType>();
}

template<class OpType, class LeftType, class RightType, class ResType>
const scalar_type &bulk_operate_impl<OpType, LeftType, RightType, ResType>::get_output_type() const
{
	return get_scalar_type<ResType>();
}

template<class T>
const scalar_type &find_next_impl<T>::get_input_type() const
{
	return get_scalar_type<T>();
}

template<class T>
const scalar_type &find_next_impl<T>::get_output_type() const
{
	return get_scalar_type<size_t>();
}

template<class T>
const scalar_type &find_prev_impl<T>::get_input_type() const
{
	return get_scalar_type<T>();
}

template<class T>
const scalar_type &find_prev_impl<T>::get_output_type() const
{
	return get_scalar_type<size_t>();
}

class local_vec_store;

/*
 * This operator is different from bulk_uoperate. It treats an array
 * as a single input and outputs an array of potentially different length.
 */
class arr_apply_operate
{
	size_t num_out_eles;
public:
	typedef std::shared_ptr<const arr_apply_operate> const_ptr;

	arr_apply_operate(size_t num_out_eles) {
		this->num_out_eles = num_out_eles;
	}
	/*
	 * This virtual method accepts an input array and stores the result
	 * in an output array.
	 */
	virtual void run(const local_vec_store &in,
			local_vec_store &out) const = 0;

	virtual const scalar_type &get_input_type() const = 0;
	virtual const scalar_type &get_output_type() const = 0;

	size_t input_entry_size() const {
		return get_input_type().get_size();
	}

	size_t output_entry_size() const {
		return get_output_type().get_size();
	}

	size_t get_num_out_eles() const {
		return num_out_eles;
	}
};

/*
 * This operator applies to an entry generated by groupby.
 * It takes a key and a value as input and the type of the value varies.
 * TODO let's not worry about the output yet.
 */
template<class T>
class gr_apply_operate
{
public:
	virtual void run(const void *key, const T &val,
			local_vec_store &vec) const = 0;

	virtual const scalar_type &get_key_type() const = 0;
	virtual const scalar_type &get_output_type() const = 0;
	// This method tells how many elements output for a key. If it's unknown
	// or isn't fixed, it should return 0.
	virtual size_t get_num_out_eles() const = 0;
};

class scatter_gather
{
public:
	virtual void scatter(const char *arr, std::vector<char *> &arrs) const = 0;
	virtual void gather(const std::vector<const char *> &arrs,
			char *arr) const = 0;
};

template<class T>
class type_scatter_gather: public scatter_gather
{
public:
	virtual void scatter(const char *arr, std::vector<char *> &arrs) const {
		const T *t_arr = (const T *) arr;
		for (size_t i = 0; i < arrs.size(); i++)
			*(T *) arrs[i] = t_arr[i];
	}

	virtual void gather(const std::vector<const char *> &arrs,
			char *arr) const {
		T *t_arr = (T *) arr;
		for (size_t i = 0; i < arrs.size(); i++)
			t_arr[i] = *(const T *) arrs[i];
	}
};

}

#endif
