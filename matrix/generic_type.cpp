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

namespace fm
{

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

template<class T>
const scatter_gather &scalar_type_impl<T>::get_sg() const
{
	static type_scatter_gather<T> sg;
	return sg;
}

template<class T>
const set_operate &scalar_type_impl<T>::get_set_const(const scalar_variable &val) const
{
	assert(val.get_type() == get_scalar_type<T>());
	const scalar_variable_impl<T> &t_val
		= static_cast<const scalar_variable_impl<T> &>(val);
	static const_set_operate<T> op(t_val.get());
	return op;
}

template<class T>
const basic_uops &scalar_type_impl<T>::get_basic_uops() const
{
	static basic_uops_impl<T, T> uops;
	return uops;
}

template<class T>
const basic_ops &scalar_type_impl<T>::get_basic_ops() const
{
	static basic_ops_impl<T, T, T> ops;
	return ops;
}

template<class T>
const agg_ops &scalar_type_impl<T>::get_agg_ops() const
{
	static agg_ops_impl<T, T> aops;
	return aops;
}

template<class T>
const type_cast &scalar_type_impl<T>::get_type_cast(const scalar_type &type) const
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

const scalar_type &get_scalar_type(prim_type type)
{
	switch(type) {
		case P_CHAR:
			return get_scalar_type<char>();
		case P_SHORT:
			return get_scalar_type<short>();
		case P_INTEGER:
			return get_scalar_type<int>();
		case P_LONG:
			return get_scalar_type<long>();
		case P_FLOAT:
			return get_scalar_type<float>();
		case P_DOUBLE:
			return get_scalar_type<double>();
		case P_LDOUBLE:
			return get_scalar_type<long double>();
		case P_BOOL:
			return get_scalar_type<bool>();
		case P_USHORT:
			return get_scalar_type<unsigned short>();
		case P_UINT:
			return get_scalar_type<unsigned int>();
		case P_ULONG:
			return get_scalar_type<unsigned long>();
		default:
			throw invalid_arg_exception("invalid prim type");
	}
}

}
