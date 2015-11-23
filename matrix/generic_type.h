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

#include <memory>

#include "sorter.h"
#include "stl_algs.h"

namespace fm
{

/**
 * Here defines the primitive types.
 */
enum prim_type
{
	P_CHAR,
	P_SHORT,
	P_INTEGER,
	P_LONG,
	P_FLOAT,
	P_DOUBLE,
	P_LDOUBLE,
	P_BOOL,
	P_USHORT,
	P_UINT,
	P_ULONG,
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

class basic_uops;
class basic_ops;
class agg_ops;
class scatter_gather;
class scalar_variable;
class rand_gen;
class set_operate;
class generic_hashtable;
class bulk_uoperate;

/**
 * This interface defines a scalar type and the operations related to the type.
 */
class scalar_type
{
public:
	/**
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
	virtual const stl_algs &get_stl_algs() const = 0;
	virtual const set_operate &get_set_const(const scalar_variable &val) const = 0;
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

/*
 * Here we implement the scalar type.
 */
template<class T>
class scalar_type_impl: public scalar_type
{
public:
	virtual std::shared_ptr<generic_hashtable> create_hashtable(
			const scalar_type &val_type) const;
	virtual const basic_uops &get_basic_uops() const;
	virtual const basic_ops &get_basic_ops() const;
	virtual const agg_ops &get_agg_ops() const;

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
	virtual const stl_algs &get_stl_algs() const {
		static stl_algs_impl<T> algs;
		return algs;
	}
	virtual const set_operate &get_set_const(const scalar_variable &val) const;
	virtual const bulk_uoperate &get_type_cast(const scalar_type &type) const;

	virtual prim_type get_type() const {
		return fm::get_type<T>();
	}

	virtual size_t get_size() const {
		return sizeof(T);
	}
};

template<class T>
const scalar_type &get_scalar_type()
{
	static scalar_type_impl<T> t;
	return t;
}

const scalar_type &get_scalar_type(prim_type type);

/**
 * This class defines a generic type for a scalar variable.
 * It shouldn't be used in an array because it has a lot of overhead.
 */
class scalar_variable
{
public:
	typedef std::shared_ptr<scalar_variable> ptr;
	typedef std::shared_ptr<const scalar_variable> const_ptr;
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

	virtual size_t get_size() const {
		return get_type().get_size();
	}
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
};

bool require_cast(const scalar_type &t1, const scalar_type &t2);

}

#endif
