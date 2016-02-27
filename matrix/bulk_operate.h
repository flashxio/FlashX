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
	virtual std::string get_name() const = 0;

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
	virtual std::string get_name() const {
		return OpType::get_name();
	}
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
	 * This performs aggregation on the input array, combines the agg result
	 * with the original agg result and stores the result on output.
	 */
	virtual void runAgg(size_t num_eles, const void *left_arr,
			void *output) const = 0;

	virtual const scalar_type &get_left_type() const = 0;
	virtual const scalar_type &get_right_type() const = 0;
	virtual const scalar_type &get_output_type() const = 0;
	virtual std::string get_name() const = 0;

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
	OpType op;
public:
	virtual void runAA(size_t num_eles, const void *left_arr1,
			const void *right_arr1, void *output_arr1) const {
		const LeftType *left_arr = (const LeftType *) left_arr1;
		const RightType *right_arr = (const RightType *) right_arr1;
		ResType *output_arr = (ResType *) output_arr1;
		for (size_t i = 0; i < num_eles; i++)
			output_arr[i] = op(left_arr[i], right_arr[i]);
	}

	virtual void runAE(size_t num_eles, const void *left_arr1,
			const void *right, void *output_arr1) const {
		const LeftType *left_arr = (const LeftType *) left_arr1;
		ResType *output_arr = (ResType *) output_arr1;
		RightType entry = *(const RightType *) right;
		for (size_t i = 0; i < num_eles; i++)
			output_arr[i] = op(left_arr[i], entry);
	}

	virtual void runEA(size_t num_eles, const void *left,
			const void *right_arr1, void *output_arr1) const {
		LeftType entry = *(const LeftType *) left;
		const RightType *right_arr = (const RightType *) right_arr1;
		ResType *output_arr = (ResType *) output_arr1;
		for (size_t i = 0; i < num_eles; i++)
			output_arr[i] = op(entry, right_arr[i]);
	}

	virtual void runAgg(size_t num_eles, const void *left_arr1,
			void *output) const {
		const LeftType *left_arr = (const LeftType *) left_arr1;
		ResType res = left_arr[0];
		for (size_t i = 1; i < num_eles; i++)
			res = op(left_arr[i], res);
		*(ResType *) output = res;
	}

	virtual const scalar_type &get_left_type() const;
	virtual const scalar_type &get_right_type() const;
	virtual const scalar_type &get_output_type() const;
	virtual std::string get_name() const {
		return OpType::get_name();
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
		CEIL,
		FLOOR,
		ROUND,
		LOG,
		LOG2,
		LOG10,
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
		static std::string get_name() {
			return "-";
		}
		OutType operator()(const InType &e) const {
			return -e;
		}
	};

	struct uop_sqrt {
		static std::string get_name() {
			return "sqrt";
		}
		double operator()(const InType &e) const {
			return std::sqrt(e);
		}
	};

	struct uop_abs {
		static std::string get_name() {
			return "abs";
		}
		OutType operator()(const InType &e) const {
			return (OutType) std::abs(e);
		}
	};

	struct uop_not {
		static std::string get_name() {
			return "not";
		}
		bool operator()(const InType &e) const {
			return !e;
		}
	};

	struct sq {
		static std::string get_name() {
			return "sq";
		}
		OutType operator()(const InType &e) const {
			return e * e;
		}
	};

	struct ceil {
		static std::string get_name() {
			return "ceil";
		}
		OutType operator()(const InType &e) const {
			return std::ceil(e);
		}
	};

	struct floor {
		static std::string get_name() {
			return "floor";
		}
		OutType operator()(const InType &e) const {
			return std::floor(e);
		}
	};

	struct round {
		static std::string get_name() {
			return "round";
		}
		OutType operator()(const InType &e) const {
			return std::round(e);
		}
	};

	struct log {
		static std::string get_name() {
			return "log";
		}
		OutType operator()(const InType &e) const {
			return std::log(e);
		}
	};

	struct log2 {
		static std::string get_name() {
			return "log2";
		}
		OutType operator()(const InType &e) const {
			return std::log2(e);
		}
	};

	struct log10 {
		static std::string get_name() {
			return "log10";
		}
		OutType operator()(const InType &e) const {
			return std::log10(e);
		}
	};

	bulk_uoperate_impl<uop_neg, InType, OutType> neg_op;
	bulk_uoperate_impl<uop_sqrt, InType, double> sqrt_op;
	bulk_uoperate_impl<uop_abs, InType, OutType> abs_op;
	bulk_uoperate_impl<uop_not, InType, bool> not_op;
	bulk_uoperate_impl<sq, InType, OutType> sq_op;
	bulk_uoperate_impl<ceil, InType, OutType> ceil_op;
	bulk_uoperate_impl<floor, InType, OutType> floor_op;
	bulk_uoperate_impl<round, InType, OutType> round_op;
	bulk_uoperate_impl<log, InType, OutType> log_op;
	bulk_uoperate_impl<log2, InType, OutType> log2_op;
	bulk_uoperate_impl<log10, InType, OutType> log10_op;

	std::vector<bulk_uoperate *> ops;
public:
	basic_uops_impl() {
		ops.push_back(&neg_op);
		ops.push_back(&sqrt_op);
		ops.push_back(&abs_op);
		ops.push_back(&not_op);
		ops.push_back(&sq_op);
		ops.push_back(&ceil_op);
		ops.push_back(&floor_op);
		ops.push_back(&round_op);
		ops.push_back(&log_op);
		ops.push_back(&log2_op);
		ops.push_back(&log10_op);
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
		NEQ,
		GT,
		GE,
		LT,
		LE,
		OR,
		AND,
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

template<class LeftType1, class RightType1, class ResType1>
struct multiply
{
	static std::string get_name() {
		return "*";
	}
	ResType1 operator()(const LeftType1 &e1, const RightType1 &e2) const {
		return e1 * e2;
	}
};

template<>
struct multiply<double, double, double>
{
	static std::string get_name() {
		return "*";
	}
	double operator()(const double &e1, const double &e2) const {
		long double first = e1;
		long double second = e2;
		return first * second;
	}
};

template<class LeftType, class RightType, class ResType>
struct min
{
	static std::string get_name() {
		return "min";
	}
	ResType operator()(const LeftType &e1, const RightType &e2) const {
		return std::min(e1, e2);
	}
};

template<>
struct min<double, double, double>
{
	static std::string get_name() {
		return "min";
	}
	double operator()(const double &e1, const double &e2) const {
		if (std::isnan(e1))
			return e1;
		else if (std::isnan(e2))
			return e2;
		else
			return std::min(e1, e2);
	}
};

template<class LeftType, class RightType, class ResType>
struct max
{
	static std::string get_name() {
		return "max";
	}
	ResType operator()(const LeftType &e1, const RightType &e2) const {
		return std::max(e1, e2);
	}
};

template<>
struct max<double, double, double>
{
	static std::string get_name() {
		return "max";
	}
	double operator()(const double &e1, const double &e2) const {
		if (std::isnan(e1))
			return e1;
		else if (std::isnan(e2))
			return e2;
		else
			return std::max(e1, e2);
	}
};

/*
 * This template implements all basic binary operators for different types.
 */
template<class LeftType, class RightType, class ResType>
class basic_ops_impl: public basic_ops
{
	struct add {
		static std::string get_name() {
			return "+";
		}
		ResType operator()(const LeftType &e1, const RightType &e2) const {
			return e1 + e2;
		}
	};

	struct sub {
		static std::string get_name() {
			return "-";
		}
		ResType operator()(const LeftType &e1, const RightType &e2) const {
			return e1 - e2;
		}
	};

	// Division is special. Its output should be float point.
	// Therefore, we convert both input values to float point.
	struct divide {
		static std::string get_name() {
			return "/";
		}
		double operator()(const LeftType &e1, const RightType &e2) const {
			double d1 = e1;
			double d2 = e2;
			return d1 / d2;
		}
	};
	struct divide_float {
		static std::string get_name() {
			return "/";
		}
		float operator()(const float &e1, const float &e2) const {
			return e1 / e2;
		}
	};

	struct pow {
		static std::string get_name() {
			return "pow";
		}
		ResType operator()(const LeftType &e1, const RightType &e2) const {
			return (ResType) std::pow(e1, e2);
		}
	};

	struct eq {
		static std::string get_name() {
			return "==";
		}
		bool operator()(const LeftType &e1, const RightType &e2) const {
			return e1 == e2;
		}
	};

	struct neq {
		static std::string get_name() {
			return "!=";
		}
		bool operator()(const LeftType &e1, const RightType &e2) const {
			return e1 != e2;
		}
	};

	struct gt {
		static std::string get_name() {
			return ">";
		}
		bool operator()(const LeftType &e1, const RightType &e2) const {
			return e1 > e2;
		}
	};

	struct ge {
		static std::string get_name() {
			return ">=";
		}
		bool operator()(const LeftType &e1, const RightType &e2) const {
			return e1 >= e2;
		}
	};

	struct lt {
		static std::string get_name() {
			return "<";
		}
		bool operator()(const LeftType &e1, const RightType &e2) const {
			return e1 < e2;
		}
	};

	struct le {
		static std::string get_name() {
			return "<=";
		}
		bool operator()(const LeftType &e1, const RightType &e2) const {
			return e1 <= e2;
		}
	};

	struct logic_or {
		static std::string get_name() {
			return "||";
		}
		bool operator()(const LeftType &e1, const RightType &e2) const {
			return e1 || e2;
		}
	};

	struct logic_and {
		static std::string get_name() {
			return "&&";
		}
		bool operator()(const LeftType &e1, const RightType &e2) const {
			return e1 && e2;
		}
	};

	bulk_operate_impl<add, LeftType, RightType, ResType> add_op;
	bulk_operate_impl<sub, LeftType, RightType, ResType> sub_op;
	bulk_operate_impl<multiply<LeftType, RightType, ResType>,
		LeftType, RightType, ResType> mul_op;
	bulk_operate_impl<divide, LeftType, RightType, double> div_op;
	bulk_operate_impl<divide_float, float, float, float> div_float_op;
	bulk_operate_impl<min<LeftType, RightType, ResType>,
		LeftType, RightType, ResType> min_op;
	bulk_operate_impl<max<LeftType, RightType, ResType>,
		LeftType, RightType, ResType> max_op;
	bulk_operate_impl<pow, LeftType, RightType, ResType> pow_op;
	bulk_operate_impl<eq, LeftType, RightType, bool> eq_op;
	bulk_operate_impl<neq, LeftType, RightType, bool> neq_op;
	bulk_operate_impl<gt, LeftType, RightType, bool> gt_op;
	bulk_operate_impl<ge, LeftType, RightType, bool> ge_op;
	bulk_operate_impl<lt, LeftType, RightType, bool> lt_op;
	bulk_operate_impl<le, LeftType, RightType, bool> le_op;
	bulk_operate_impl<logic_or, LeftType, RightType, bool> or_op;
	bulk_operate_impl<logic_and, LeftType, RightType, bool> and_op;

	std::vector<bulk_operate *> ops;
public:
	basic_ops_impl() {
		ops.resize(NUM_OPS);
		ops[ADD] = &add_op;
		ops[SUB] = &sub_op;
		ops[MUL] = &mul_op;
		// We should treat float differently, because we want to output float.
		if (get_type<LeftType>() == prim_type::P_FLOAT
				&& get_type<RightType>() == prim_type::P_FLOAT)
			ops[DIV] = &div_float_op;
		else
			ops[DIV] = &div_op;
		ops[MIN] = &min_op;
		ops[MAX] = &max_op;
		ops[POW] = &pow_op;
		ops[EQ] = &eq_op;
		ops[NEQ] = &neq_op;
		ops[GT] = &gt_op;
		ops[GE] = &ge_op;
		ops[LT] = &lt_op;
		ops[LE] = &le_op;
		ops[OR] = &or_op;
		ops[AND] = &and_op;
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

class set_vv_operate
{
public:
	typedef std::shared_ptr<const set_vv_operate> const_ptr;

	virtual void set(off_t arr_idx, void *arr, size_t num_eles) const = 0;
	virtual const scalar_type &get_type() const = 0;
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

class conv_layout
{
public:
	/*
	 * Move data in a matrix stored in `arrs' to a single piece of memory
	 * pointed by `contig_arr'.
	 * The number of elements in each array of the vector is `arr_len'.
	 */
	virtual void conv(const std::vector<const char *> &arrs, size_t arr_len,
			char *contig_arr) const = 0;
	/*
	 * Move data in a single piece of memory to a matrix stored in `arrs'.
	 * The number of elements in the contiguous memory is `arr_len'.
	 */
	virtual void conv(const char *contig_arr, size_t arr_len,
			const std::vector<char *> &arrs) const = 0;
};

template<class T>
class type_conv_layout: public conv_layout
{
public:
	virtual void conv(const std::vector<const char *> &arrs, size_t arr_len,
			char *contig_arr) const {
		std::vector<const T *> t_arrs(arrs.size());
		for (size_t i = 0; i < arrs.size(); i++)
			t_arrs[i] = reinterpret_cast<const T *>(arrs[i]);
		T *t_res = reinterpret_cast<T *>(contig_arr);
		size_t res_idx = 0;
		for (size_t i = 0; i < arr_len; i++)
			for (size_t j = 0; j < t_arrs.size(); j++)
				t_res[res_idx++] = t_arrs[j][i];

	}
	virtual void conv(const char *contig_arr, size_t arr_len,
			const std::vector<char *> &arrs) const {
		const T *t_arr = reinterpret_cast<const T *>(contig_arr);
		std::vector<T *> t_res(arrs.size());
		for (size_t i = 0; i < t_res.size(); i++)
			t_res[i] = reinterpret_cast<T *>(arrs[i]);
		size_t each_len = arr_len / arrs.size();
		size_t src_idx = 0;
		for (size_t i = 0; i < each_len; i++)
			for (size_t j = 0; j < t_res.size(); j++)
				t_res[j][i] = t_arr[src_idx++];
	}
};

}

#endif
