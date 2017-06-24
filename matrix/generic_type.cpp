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

#include <string.h>

#include "log.h"
#include "comm_exception.h"

#include "generic_type.h"
#include "rand_gen.h"
#include "bulk_operate.h"
#include "bulk_operate_ext.h"
#include "generic_hashtable.h"

namespace fm
{

std::vector<scalar_type::ptr> scalar_type::types;
std::vector<basic_ops::ptr> scalar_type::basic_ops_impls;
std::vector<basic_uops::ptr> scalar_type::basic_uops_impls;
std::vector<agg_ops::ptr> scalar_type::agg_ops_impls;

/*
 * Here we implement the scalar type.
 */
template<class T>
class scalar_type_impl: public scalar_type
{
public:
	virtual std::string get_name() const {
		return get_type_str<T>();
	}
	virtual bool is_floating_point() const {
		return std::is_floating_point<T>::value;
	}
	virtual std::string conv2str(const char *arr, size_t num_eles,
			const std::string &sep) const;
	virtual std::shared_ptr<generic_hashtable> create_hashtable(
			const scalar_type &val_type) const;
	virtual const basic_uops &get_basic_uops() const {
		int type = (int) fm::get_type<T>();
		return *basic_uops_impls[type];
	}
	virtual const basic_ops &get_basic_ops() const {
		int type = (int) fm::get_type<T>();
		return *basic_ops_impls[type];
	}
	virtual const agg_ops &get_agg_ops() const {
		int type = (int) fm::get_type<T>();
		return *agg_ops_impls[type];
	}

	virtual std::shared_ptr<scalar_variable> create_scalar() const;
	virtual std::shared_ptr<rand_gen> create_randu_gen(const scalar_variable &min,
			const scalar_variable &max) const;
	virtual std::shared_ptr<rand_gen> create_randu_gen(const scalar_variable &min,
			const scalar_variable &max, const scalar_variable &seed) const;
	virtual std::shared_ptr<rand_gen> create_randn_gen(const scalar_variable &mean,
			const scalar_variable &var) const;
	virtual std::shared_ptr<rand_gen> create_randn_gen(const scalar_variable &mean,
			const scalar_variable &var, const scalar_variable &seed) const;

	virtual const sorter &get_sorter() const {
		static type_sorter<T> sort;
		return sort;
	}

	virtual const scatter_gather &get_sg() const;
	virtual const ifelse &get_ifelse() const;
	virtual const conv_layout &get_conv() const;
	virtual const stl_algs &get_stl_algs() const {
		static stl_algs_impl<T> algs;
		return algs;
	}
	virtual std::shared_ptr<const set_operate> get_set_const(
			const scalar_variable &val) const;
	virtual std::shared_ptr<const set_vec_operate> get_set_vec_const(
			const scalar_variable &val) const;
	virtual std::shared_ptr<const set_operate> get_set_seq(
			const scalar_variable &start, const scalar_variable &stride,
			size_t num_rows, size_t num_cols, bool byrow,
			matrix_layout_t layout) const;
	virtual const bulk_uoperate &get_type_cast(const scalar_type &type) const;

	virtual prim_type get_type() const {
		return fm::get_type<T>();
	}

	virtual size_t get_size() const {
		return sizeof(T);
	}
};

void scalar_type::init()
{
	/*
	 * We need to initialize all the types first before we can initialize
	 * the operations on the types.
	 */
	types.resize((int) prim_type::NUM_TYPES);
	types[fm::get_type<char>()] = scalar_type::ptr(new scalar_type_impl<char>());
	types[fm::get_type<short>()] = scalar_type::ptr(new scalar_type_impl<short>());
	types[fm::get_type<int>()] = scalar_type::ptr(new scalar_type_impl<int>());
	types[fm::get_type<long>()] = scalar_type::ptr(new scalar_type_impl<long>());
	types[fm::get_type<float>()] = scalar_type::ptr(new scalar_type_impl<float>());
	types[fm::get_type<double>()] = scalar_type::ptr(new scalar_type_impl<double>());
	types[fm::get_type<long double>()] = scalar_type::ptr(new scalar_type_impl<long double>());
	types[fm::get_type<bool>()] = scalar_type::ptr(new scalar_type_impl<bool>());
	types[fm::get_type<unsigned short>()] = scalar_type::ptr(new scalar_type_impl<unsigned short>());
	types[fm::get_type<unsigned int>()] = scalar_type::ptr(new scalar_type_impl<unsigned int>());
	types[fm::get_type<unsigned long>()] = scalar_type::ptr(new scalar_type_impl<unsigned long>());

	for (size_t i = 0; i < types.size(); i++)
		if (types[i] == NULL)
			throw unsupported_exception("find an unsupported type");

	init_ops();
}

// This initializes all of the scalar types supported by FlashMatrix.
class scalar_type_initializer
{
public:
	scalar_type_initializer() {
		scalar_type::init();
	}
};
static scalar_type_initializer initializer;

template<class T>
generic_hashtable::ptr scalar_type_impl<T>::create_hashtable(
		const scalar_type &val_type) const
{
	switch(val_type.get_size()) {
		case 4: return generic_hashtable::ptr(
						new generic_hashtable_impl<T, 4>(val_type));
		case 8: return generic_hashtable::ptr(
						new generic_hashtable_impl<T, 8>(val_type));
		default: BOOST_LOG_TRIVIAL(error) << "Can't create a generic hashtable";
				 return generic_hashtable::ptr();
	}
}

scalar_variable::ptr scalar_variable::cast_type(const scalar_type &type) const
{
	const scalar_type &vtype = get_type();
	const bulk_uoperate &cast_op = vtype.get_type_cast(type);
	scalar_variable::ptr ret = type.create_scalar();
	cast_op.runA(1, get_raw(), ret->get_raw());
	return ret;
}

template<class T>
std::string conv2str(T val)
{
	return std::to_string(val);
}

template<>
std::string conv2str<float>(float val)
{
	char str[20];
	snprintf(str, 20, "%g", val);
	return str;
}

template<>
std::string conv2str<double>(double val)
{
	char str[20];
	snprintf(str, 20, "%g", val);
	return str;
}

template<class T>
std::string scalar_type_impl<T>::conv2str(const char *arr, size_t num_eles,
		const std::string &sep) const
{
	const T *tarr = reinterpret_cast<const T *>(arr);
	assert(num_eles > 0);
	std::string ret = fm::conv2str<T>(tarr[0]);
	for (size_t i = 1; i < num_eles; i++)
		ret += sep + fm::conv2str<T>(tarr[i]);
	return ret;
}

template<class T>
scalar_variable::ptr scalar_type_impl<T>::create_scalar() const
{
	return scalar_variable::ptr(new scalar_variable_impl<T>());
}

template<class T>
rand_gen::ptr scalar_type_impl<T>::create_randu_gen(const scalar_variable &min,
		const scalar_variable &max) const
{
	scalar_variable_impl<T> &t_min = (scalar_variable_impl<T> &) min;
	scalar_variable_impl<T> &t_max = (scalar_variable_impl<T> &) max;
	return rand_gen::create_randu<T>(t_min.get(), t_max.get());
}

template<class T>
rand_gen::ptr scalar_type_impl<T>::create_randu_gen(const scalar_variable &min,
		const scalar_variable &max, const scalar_variable &seed) const
{
	scalar_variable_impl<T> &t_min = (scalar_variable_impl<T> &) min;
	scalar_variable_impl<T> &t_max = (scalar_variable_impl<T> &) max;
	scalar_variable_impl<T> &t_seed = (scalar_variable_impl<T> &) seed;
	return rand_gen::create_randu<T>(t_min.get(), t_max.get(), t_seed.get());
}

template<>
rand_gen::ptr scalar_type_impl<char>::create_randn_gen(const scalar_variable &mean,
		const scalar_variable &var) const
{
	assert(0);
	return rand_gen::ptr();
}

template<>
rand_gen::ptr scalar_type_impl<char>::create_randn_gen(const scalar_variable &mean,
		const scalar_variable &var, const scalar_variable &seed) const
{
	assert(0);
	return rand_gen::ptr();
}

template<>
rand_gen::ptr scalar_type_impl<int>::create_randn_gen(const scalar_variable &mean,
		const scalar_variable &var) const
{
	assert(0);
	return rand_gen::ptr();
}

template<>
rand_gen::ptr scalar_type_impl<int>::create_randn_gen(const scalar_variable &mean,
		const scalar_variable &var, const scalar_variable &seed) const
{
	assert(0);
	return rand_gen::ptr();
}

template<>
rand_gen::ptr scalar_type_impl<long>::create_randn_gen(const scalar_variable &mean,
		const scalar_variable &var) const
{
	assert(0);
	return rand_gen::ptr();
}

template<>
rand_gen::ptr scalar_type_impl<long>::create_randn_gen(const scalar_variable &mean,
		const scalar_variable &var, const scalar_variable &seed) const
{
	assert(0);
	return rand_gen::ptr();
}

template<>
rand_gen::ptr scalar_type_impl<short>::create_randn_gen(const scalar_variable &mean,
		const scalar_variable &var) const
{
	assert(0);
	return rand_gen::ptr();
}

template<>
rand_gen::ptr scalar_type_impl<short>::create_randn_gen(const scalar_variable &mean,
		const scalar_variable &var, const scalar_variable &seed) const
{
	assert(0);
	return rand_gen::ptr();
}

template<>
rand_gen::ptr scalar_type_impl<unsigned short>::create_randn_gen(const scalar_variable &mean,
		const scalar_variable &var) const
{
	assert(0);
	return rand_gen::ptr();
}

template<>
rand_gen::ptr scalar_type_impl<unsigned short>::create_randn_gen(const scalar_variable &mean,
		const scalar_variable &var, const scalar_variable &seed) const
{
	assert(0);
	return rand_gen::ptr();
}

template<>
rand_gen::ptr scalar_type_impl<unsigned int>::create_randn_gen(const scalar_variable &mean,
		const scalar_variable &var) const
{
	assert(0);
	return rand_gen::ptr();
}

template<>
rand_gen::ptr scalar_type_impl<unsigned int>::create_randn_gen(const scalar_variable &mean,
		const scalar_variable &var, const scalar_variable &seed) const
{
	assert(0);
	return rand_gen::ptr();
}

template<>
rand_gen::ptr scalar_type_impl<unsigned long>::create_randn_gen(const scalar_variable &mean,
		const scalar_variable &var) const
{
	assert(0);
	return rand_gen::ptr();
}

template<>
rand_gen::ptr scalar_type_impl<unsigned long>::create_randn_gen(const scalar_variable &mean,
		const scalar_variable &var, const scalar_variable &seed) const
{
	assert(0);
	return rand_gen::ptr();
}

template<>
rand_gen::ptr scalar_type_impl<bool>::create_randn_gen(const scalar_variable &mean,
		const scalar_variable &var) const
{
	assert(0);
	return rand_gen::ptr();
}

template<>
rand_gen::ptr scalar_type_impl<bool>::create_randn_gen(const scalar_variable &mean,
		const scalar_variable &var, const scalar_variable &seed) const
{
	assert(0);
	return rand_gen::ptr();
}

template<class T>
rand_gen::ptr scalar_type_impl<T>::create_randn_gen(const scalar_variable &mean,
		const scalar_variable &var) const
{
	scalar_variable_impl<T> &t_mean = (scalar_variable_impl<T> &) mean;
	scalar_variable_impl<T> &t_var = (scalar_variable_impl<T> &) var;
	return rand_gen::create_randn<T>(t_mean.get(), t_var.get());
}

template<class T>
rand_gen::ptr scalar_type_impl<T>::create_randn_gen(const scalar_variable &mean,
		const scalar_variable &var, const scalar_variable &seed) const
{
	scalar_variable_impl<T> &t_mean = (scalar_variable_impl<T> &) mean;
	scalar_variable_impl<T> &t_var = (scalar_variable_impl<T> &) var;
	scalar_variable_impl<T> &t_seed = (scalar_variable_impl<T> &) seed;
	return rand_gen::create_randn<T>(t_mean.get(), t_var.get(), t_seed.get());
}

namespace
{

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

template<class T>
class type_ifelse: public ifelse
{
public:
	virtual void run(const bool *cond, size_t len, const char *arr1,
			const char *arr2, char *arr3) const {
		const T *A = reinterpret_cast<const T *>(arr1);
		const T *B = reinterpret_cast<const T *>(arr2);
		T *C = reinterpret_cast<T *>(arr3);
		for (size_t i = 0; i < len; i++) {
			if (cond[i])
				C[i] = A[i];
			else
				C[i] = B[i];
		}
	}
};

template<class T>
class type_conv_layout: public conv_layout
{
public:
	virtual void conv1(const std::vector<const char *> &arrs, size_t arr_len,
			char *contig_arr) const {
		std::vector<const T *> t_arrs(arrs.size());
		for (size_t i = 0; i < arrs.size(); i++)
			t_arrs[i] = reinterpret_cast<const T *>(arrs[i]);
		T *t_res = reinterpret_cast<T *>(contig_arr);
		for (size_t j = 0; j < t_arrs.size(); j++)
			for (size_t i = 0; i < arr_len; i++)
				t_res[i * arrs.size() + j] = t_arrs[j][i];
	}
	virtual void conv2(const char *contig_arr, size_t arr_len,
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

/*
 * Set the data of a matrix with sequence numbers.
 */
template<class T>
class set_seq: public type_set_operate<T>
{
	T start;
	// The stride between two adjacent elements in the sequence.
	T stride;
	// The stride between two elements stored contiguously.
	// If the sequence number is placed in matrix by rows,
	// * For row-major matrices, it's the same as `stride'.
	// * For col-major matrices, it's `stride * round_len'.
	// If the sequence number is placed by cols,
	// * For row-major matrices, it's `stride * round_len'.
	// * For col-major matrices, it's `stride'.
	T seq_ele_stride;
	size_t round_len;
	bool byrow;
	set_seq(T start, T stride, T seq_ele_stride, size_t round_len, bool byrow) {
		this->start = start;
		this->stride = stride;
		this->seq_ele_stride = seq_ele_stride;
		this->round_len = round_len;
		this->byrow = byrow;
	}
public:
	set_seq(T start, T stride, size_t round_len, bool byrow,
			matrix_layout_t layout) {
		this->start = start;
		this->stride = stride;
		this->round_len = round_len;
		this->byrow = byrow;

		if (layout == matrix_layout_t::L_ROW && byrow)
			this->seq_ele_stride = stride;
		else if (layout == matrix_layout_t::L_COL && byrow)
			this->seq_ele_stride = stride * round_len;
		else if (layout == matrix_layout_t::L_ROW)
			this->seq_ele_stride = stride * round_len;
		else
			this->seq_ele_stride = stride;
	}

	void set(T *arr, size_t num_eles, off_t row_idx, off_t col_idx) const {
		T curr_start;
		if (byrow)
			curr_start = start + (row_idx * round_len + col_idx) * stride;
		else
			curr_start = start + (col_idx * round_len + row_idx) * stride;

		for (size_t i = 0; i < num_eles; i++)
			arr[i] = curr_start + i * seq_ele_stride;
	}

	virtual set_operate::const_ptr transpose() const {
		return set_operate::const_ptr(new set_seq<T>(start, stride,
					seq_ele_stride, round_len, !byrow));
	}
};

}

template<class T>
const scatter_gather &scalar_type_impl<T>::get_sg() const
{
	static type_scatter_gather<T> sg;
	return sg;
}

template<class T>
const ifelse &scalar_type_impl<T>::get_ifelse() const
{
	static type_ifelse<T> ie;
	return ie;
}

template<class T>
const conv_layout &scalar_type_impl<T>::get_conv() const
{
	static type_conv_layout<T> sg;
	return sg;
}

template<class T>
set_operate::const_ptr scalar_type_impl<T>::get_set_const(
		const scalar_variable &val) const
{
	assert(val.get_type() == get_scalar_type<T>());
	const scalar_variable_impl<T> &t_val
		= static_cast<const scalar_variable_impl<T> &>(val);
	return set_operate::const_ptr(new const_set_operate<T>(t_val.get()));
}

template<class T>
set_vec_operate::const_ptr scalar_type_impl<T>::get_set_vec_const(
		const scalar_variable &val) const
{
	assert(val.get_type() == get_scalar_type<T>());
	const scalar_variable_impl<T> &t_val
		= static_cast<const scalar_variable_impl<T> &>(val);
	return set_vec_operate::const_ptr(new const_set_vec_operate<T>(t_val.get()));
}

template<class T>
set_operate::const_ptr scalar_type_impl<T>::get_set_seq(
		const scalar_variable &start, const scalar_variable &stride,
		size_t num_rows, size_t num_cols, bool byrow,
		matrix_layout_t layout) const
{
	assert(start.get_type() == get_scalar_type<T>());
	assert(stride.get_type() == get_scalar_type<T>());
	const scalar_variable_impl<T> &t_start
		= static_cast<const scalar_variable_impl<T> &>(start);
	const scalar_variable_impl<T> &t_stride
		= static_cast<const scalar_variable_impl<T> &>(stride);
	if (byrow)
		return set_operate::const_ptr(new set_seq<T>(t_start.get(),
					t_stride.get(), num_cols, byrow, layout));
	else
		return set_operate::const_ptr(new set_seq<T>(t_start.get(),
					t_stride.get(), num_rows, byrow, layout));
}

namespace
{

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
	virtual std::string get_name() const {
		return boost::str(boost::format("cast_%1%2%2%") % get_type_str<T1>()
				% get_type_str<T2>());
	}
};

template<class T1, class T2>
const bulk_uoperate &get_type_cast()
{
	static type_cast<T1, T2> cast;
	return cast;
}

}

template<class T>
const bulk_uoperate &scalar_type_impl<T>::get_type_cast(const scalar_type &type) const
{
	switch(type.get_type()) {
		case P_CHAR:
			return fm::get_type_cast<T, char>();
		case P_SHORT:
			return fm::get_type_cast<T, short>();
		case P_INTEGER:
			return fm::get_type_cast<T, int>();
		case P_LONG:
			return fm::get_type_cast<T, long>();
		case P_FLOAT:
			return fm::get_type_cast<T, float>();
		case P_DOUBLE:
			return fm::get_type_cast<T, double>();
		case P_LDOUBLE:
			return fm::get_type_cast<T, long double>();
		case P_BOOL:
			return fm::get_type_cast<T, bool>();
		case P_USHORT:
			return fm::get_type_cast<T, unsigned short>();
		case P_UINT:
			return fm::get_type_cast<T, unsigned int>();
		case P_ULONG:
			return fm::get_type_cast<T, unsigned long>();
		default:
			throw invalid_arg_exception("invalid prim type");
	}
}

bool require_cast(const scalar_type &t1, const scalar_type &t2)
{
	if (t1 == t2)
		return false;
	// If the two types require different memory storage size, we definitely
	// need to cast them.
	if (t1.get_size() != t2.get_size())
		return true;

	if ((t1.get_type() == prim_type::P_SHORT
				&& t2.get_type() == prim_type::P_USHORT)
			|| (t1.get_type() == prim_type::P_INTEGER
				&& t2.get_type() == prim_type::P_UINT)
			|| (t1.get_type() == prim_type::P_LONG
				&& t2.get_type() == prim_type::P_ULONG))
		return false;

	if ((t2.get_type() == prim_type::P_SHORT
				&& t1.get_type() == prim_type::P_USHORT)
			|| (t2.get_type() == prim_type::P_INTEGER
				&& t1.get_type() == prim_type::P_UINT)
			|| (t2.get_type() == prim_type::P_LONG
				&& t1.get_type() == prim_type::P_ULONG))
		return false;

	// We need to cast the rest of the type pairs.
	return true;
}

const scalar_type &get_larger_type(const scalar_type &t1, const scalar_type &t2)
{
	if (t1.get_type() > t2.get_type())
		return t1;
	else
		return t2;
}

}
