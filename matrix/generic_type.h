#ifndef __GENERAL_TYPE_H__
#define __GENERAL_TYPE_H__

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
#include <stdlib.h>

#include <sstream>
#include <string>
#include <memory>

#include "sorter.h"
#include "stl_algs.h"

namespace fm
{

enum matrix_layout_t
{
	L_COL,
	L_ROW,
	L_ROW_2D,
	// It indicates that the layout isn't defined.
	L_NONE,
};

enum matrix_margin
{
	MAR_ROW = 1,
	MAR_COL = 2,
	BOTH,
};

/**
 * Here defines the primitive types.
 * The types in the list are also ordered to help arithmetic conversion.
 * The types in the front has lower size than the ones in the end.
 */
enum prim_type
{
	P_BOOL,
	P_CHAR,
	P_SHORT,
	P_USHORT,
	P_INTEGER,
	P_UINT,
	P_LONG,
	P_ULONG,
	P_FLOAT,
	P_DOUBLE,
	P_LDOUBLE,
	NUM_TYPES,
};

template<class T>
prim_type get_type()
{
	return prim_type::NUM_TYPES;
}

template<>
inline prim_type get_type<char>()
{
	return prim_type::P_CHAR;
}

template<>
inline prim_type get_type<short>()
{
	return prim_type::P_SHORT;
}

template<>
inline prim_type get_type<int>()
{
	return prim_type::P_INTEGER;
}

template<>
inline prim_type get_type<long>()
{
	return prim_type::P_LONG;
}

template<>
inline prim_type get_type<float>()
{
	return prim_type::P_FLOAT;
}

template<>
inline prim_type get_type<double>()
{
	return prim_type::P_DOUBLE;
}

template<>
inline prim_type get_type<long double>()
{
	return prim_type::P_LDOUBLE;
}

template<>
inline prim_type get_type<bool>()
{
	return prim_type::P_BOOL;
}

template<>
inline prim_type get_type<unsigned short>()
{
	return prim_type::P_USHORT;
}

template<>
inline prim_type get_type<unsigned int>()
{
	return prim_type::P_UINT;
}

template<>
inline prim_type get_type<unsigned long>()
{
	return prim_type::P_ULONG;
}

template<class T>
std::string get_type_str()
{
	return "unknown";
}

template<>
inline std::string get_type_str<char>()
{
	return "char";
}

template<>
inline std::string get_type_str<short>()
{
	return "short";
}

template<>
inline std::string get_type_str<int>()
{
	return "int";
}

template<>
inline std::string get_type_str<long>()
{
	return "long";
}

template<>
inline std::string get_type_str<float>()
{
	return "float";
}

template<>
inline std::string get_type_str<double>()
{
	return "double";
}

template<>
inline std::string get_type_str<long double>()
{
	return "ldouble";
}

template<>
inline std::string get_type_str<bool>()
{
	return "bool";
}

template<>
inline std::string get_type_str<unsigned short>()
{
	return "ushort";
}

template<>
inline std::string get_type_str<unsigned int>()
{
	return "uint";
}

template<>
inline std::string get_type_str<unsigned long>()
{
	return "ulong";
}

class basic_uops;
class basic_ops;
class agg_ops;
class scatter_gather;
class conv_layout;
class scalar_variable;
class rand_gen;
class set_operate;
class set_vec_operate;
class generic_hashtable;
class bulk_uoperate;
class ifelse;

/**
 * This interface defines a scalar type and the operations related to the type.
 */
class scalar_type
{
	static std::vector<std::shared_ptr<scalar_type> > types;

	static std::vector<std::shared_ptr<basic_ops> > basic_ops_impls;
	static std::vector<std::shared_ptr<basic_uops> > basic_uops_impls;
	static std::vector<std::shared_ptr<agg_ops> > agg_ops_impls;
public:
	typedef std::shared_ptr<scalar_type> ptr;

	// This registers all scalar types supported by FlashMatrix.
	static void init();
	static void init_ops();

	static const scalar_type &get_type(prim_type type) {
		return *types[(int) type];
	}

	virtual std::string get_name() const = 0;
	virtual bool is_floating_point() const = 0;
	virtual std::string conv2str(const char *arr, size_t num_eles,
			const std::string &sep) const = 0;

	/*
	 * The operators that work on this type.
	 */
	virtual std::shared_ptr<generic_hashtable> create_hashtable(
			const scalar_type &val_type) const = 0;
	virtual const basic_uops &get_basic_uops() const = 0;
	virtual const basic_ops &get_basic_ops() const = 0;
	virtual const agg_ops &get_agg_ops() const = 0;
	virtual prim_type get_type() const = 0;
	virtual size_t get_size() const = 0;
	virtual const sorter &get_sorter() const = 0;
	virtual const scatter_gather &get_sg() const = 0;
	virtual const ifelse &get_ifelse() const = 0;
	virtual const conv_layout &get_conv() const = 0;
	virtual const stl_algs &get_stl_algs() const = 0;
	virtual std::shared_ptr<const set_operate> get_set_const(
			const scalar_variable &val) const = 0;
	virtual std::shared_ptr<const set_vec_operate> get_set_vec_const(
			const scalar_variable &val) const = 0;
	virtual std::shared_ptr<const set_operate> get_set_seq(
			const scalar_variable &start, const scalar_variable &stride,
			size_t num_rows, size_t num_cols, bool byrow,
			matrix_layout_t layout) const = 0;
	virtual std::shared_ptr<scalar_variable> create_scalar() const = 0;
	// Create Random generator with the uniform distribution.
	virtual std::shared_ptr<rand_gen> create_randu_gen(const scalar_variable &min,
			const scalar_variable &max) const = 0;
	virtual std::shared_ptr<rand_gen> create_randu_gen(const scalar_variable &min,
			const scalar_variable &max, const scalar_variable &seed) const = 0;
	// Create Random generator with the normal distribution.
	virtual std::shared_ptr<rand_gen> create_randn_gen(const scalar_variable &mean,
			const scalar_variable &var) const = 0;
	virtual std::shared_ptr<rand_gen> create_randn_gen(const scalar_variable &mean,
			const scalar_variable &var, const scalar_variable &seed) const = 0;
	virtual const bulk_uoperate &get_type_cast(const scalar_type &type) const = 0;

	virtual bool operator==(const scalar_type &type) const {
		return get_type() == type.get_type();
	}

	virtual bool operator!=(const scalar_type &type) const {
		return get_type() != type.get_type();
	}
};

template<class T>
const scalar_type &get_scalar_type()
{
	return scalar_type::get_type(get_type<T>());
}

static inline const scalar_type &get_scalar_type(prim_type type)
{
	return scalar_type::get_type(type);
}

/**
 * This class defines a generic type for a scalar variable.
 * It shouldn't be used in an array because it has a lot of overhead.
 */
class scalar_variable
{
public:
	typedef std::shared_ptr<scalar_variable> ptr;
	typedef std::shared_ptr<const scalar_variable> const_ptr;

	template<class T>
	static T get_val(const scalar_variable &var) {
		assert(var.get_type() == get_scalar_type<T>());
		const T *val = reinterpret_cast<const T *>(var.get_raw());
		return *val;
	}
	/**
	 * Get the raw representation of the type.
	 */
	virtual const char *get_raw() const = 0;
	virtual char *get_raw() = 0;
	/**
	 * The type.
	 */
	virtual const scalar_type &get_type() const = 0;
	/**
	 * Set the value of the scalar variable in the raw representation.
	 */
	virtual bool set_raw(const char *v, int size) = 0;
	/*
	 * Test if the value is equal to the one stored in the given address.
	 */
	virtual bool equals(const char *addr) const = 0;

	/*
	 * Get the text representation of the variable.
	 */
	virtual std::string get_name() const = 0;

	virtual size_t get_size() const {
		return get_type().get_size();
	}

	scalar_variable::ptr cast_type(const scalar_type &type) const;
};

template<class T>
class scalar_variable_impl: public scalar_variable
{
	T v;
public:
	scalar_variable_impl() {
		v = 0;
	}

	scalar_variable_impl(T v) {
		this->v = v;
	}

	virtual const char *get_raw() const {
		return (const char *) &v;
	}
	virtual char *get_raw() {
		return (char *) &v;
	}

	virtual const scalar_type &get_type() const {
		return get_scalar_type<T>();
	}

	virtual bool set_raw(const char *v, int size) {
		if (sizeof(T) != size)
			return false;

		memcpy(&this->v, v, size);
		return true;
	}

	virtual bool equals(const char *addr) const {
		return v == *(const T *) addr;
	}

	T get() const {
		return v;
	}

	void set(T v) {
		this->v = v;
	}

	virtual std::string get_name() const {
		std::ostringstream conv;
		conv << v;
		return conv.str();
	}
};

bool require_cast(const scalar_type &t1, const scalar_type &t2);
/*
 * This function returns the type that can contain more values.
 */
const scalar_type &get_larger_type(const scalar_type &t1, const scalar_type &t2);

class set_operate;

template<class T>
std::shared_ptr<const set_operate> create_urand_init(T _min, T _max)
{
	extern std::shared_ptr<const set_operate> create_urand_init(
			scalar_variable::const_ptr min, scalar_variable::const_ptr max);
	scalar_variable::const_ptr min(new scalar_variable_impl<T>(_min));
	scalar_variable::const_ptr max(new scalar_variable_impl<T>(_max));
	return create_urand_init(min, max);
}

template<class T>
std::shared_ptr<const set_operate> create_nrand_init(T _mean, T _var)
{
	extern std::shared_ptr<const set_operate> create_nrand_init(
			scalar_variable::const_ptr mean, scalar_variable::const_ptr var);
	scalar_variable::const_ptr mean(new scalar_variable_impl<T>(_mean));
	scalar_variable::const_ptr var(new scalar_variable_impl<T>(_var));
	return create_nrand_init(mean, var);
}

}

#endif
