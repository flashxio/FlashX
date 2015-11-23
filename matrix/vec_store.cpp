/*
 * Copyright 2015 Open Connectome Project (http://openconnecto.me)
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

#include "vec_store.h"
#include "mem_vec_store.h"
#include "NUMA_vector.h"
#include "EM_vector.h"

namespace fm
{

namespace detail
{

vec_store::ptr vec_store::create(size_t length, const scalar_type &type,
		int num_nodes, bool in_mem)
{
	if (in_mem && num_nodes < 0)
		return smp_vec_store::create(length, type);
	else if (in_mem)
		return NUMA_vec_store::create(length, num_nodes, type);
	else
		return EM_vec_store::create(length, type);
}

size_t vec_store::copy_to(char *data, size_t num_eles) const
{
	size_t num_copy_eles = std::min(get_length(), num_eles);
	size_t portion_size = get_portion_size();
	size_t entry_size = get_type().get_size();
	for (size_t idx = 0; idx < num_copy_eles; idx += portion_size) {
		size_t len = std::min(portion_size, num_copy_eles - idx);
		local_vec_store::const_ptr store = get_portion(idx, len);
		assert(store);
		memcpy(data + idx * entry_size, store->get_raw_arr(), len * entry_size);
	}
	return num_copy_eles;
}

template<>
vec_store::ptr create_seq_vec_store<double>(double start, double end,
		double stride, int num_nodes, bool in_mem)
{
	// The result of division may generate a real number slightly smaller than
	// what we want because of the representation precision in the machine.
	// When we convert the real number to an integer, we may find the number
	// is smaller than exepcted. We need to add a very small number to
	// the real number to correct the problem.
	// TODO is it the right way to correct the problem?
	long n = (end - start) / stride + 1e-9;
	if (n < 0) {
		BOOST_LOG_TRIVIAL(error) <<"wrong sign in 'by' argument";
		return vec_store::ptr();
	}
	// We need to count the start element.
	n++;

	detail::vec_store::ptr v = detail::vec_store::create(n,
			get_scalar_type<double>(), num_nodes, in_mem);
	v->set_data(seq_set_vec_operate<double>(n, start, stride));
	return v;
}

}
}
