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

#include "common.h"

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
