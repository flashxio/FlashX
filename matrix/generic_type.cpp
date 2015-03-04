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
	return type_mem_vector<T>::create(length);
}

template<class T>
mem_vector::ptr scalar_type_impl<T>::create_mem_vec(std::shared_ptr<char> data,
			size_t num_bytes) const
{
	return type_mem_vector<T>::create(data, num_bytes);
}

template<class T>
mem_vector_vector::ptr scalar_type_impl<T>::create_mem_vec_vec() const
{
	return type_mem_vector_vector<T>::create();
}

template class scalar_type_impl<bool>;
template class scalar_type_impl<char>;
template class scalar_type_impl<short>;
template class scalar_type_impl<int>;
template class scalar_type_impl<long>;
template class scalar_type_impl<float>;
template class scalar_type_impl<double>;
template class scalar_type_impl<unsigned short>;
template class scalar_type_impl<unsigned int>;
template class scalar_type_impl<unsigned long>;

}
