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

#include "log.h"

#include "generic_type.h"
#include "mem_vector.h"
#include "mem_vector_vector.h"

namespace fm
{

template<class T>
mem_vector::ptr scalar_type_impl<T>::create_mem_vec(
		size_t length) const
{
	return mem_vector::create(length, get_scalar_type<T>());
}

template<class T>
mem_vector::ptr scalar_type_impl<T>::create_mem_vec(std::shared_ptr<char> data,
			size_t num_bytes) const
{
	return mem_vector::create(data, num_bytes, get_scalar_type<T>());
}

template<class T>
scalar_variable::ptr scalar_type_impl<T>::create_scalar() const
{
	return scalar_variable::ptr(new scalar_variable_impl<T>());
}

template<class T>
mem_vector_vector::ptr scalar_type_impl<T>::create_mem_vec_vec() const
{
	return type_mem_vector_vector<T>::create();
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
