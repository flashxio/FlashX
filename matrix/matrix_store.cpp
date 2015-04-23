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

#include "matrix_store.h"
#include "local_matrix_store.h"

namespace fm
{

namespace detail
{

size_t matrix_store::get_num_portions() const
{
	std::pair<size_t, size_t> chunk_size = get_portion_size();
	if (is_wide()) {
		assert(chunk_size.first == get_num_rows());
		return ceil(((double) get_num_cols()) / chunk_size.second);
	}
	else {
		assert(chunk_size.second == get_num_cols());
		return ceil(((double) get_num_rows()) / chunk_size.first);
	}
}

void matrix_store::reset_data()
{
	size_t num_chunks = get_num_portions();
	if (is_wide()) {
#pragma omp parallel for
		for (size_t i = 0; i < num_chunks; i++) {
			local_matrix_store::ptr local_store = get_portion(i);
			local_store->reset_data();
		}
	}
	else {
#pragma omp parallel for
		for (size_t i = 0; i < num_chunks; i++) {
			local_matrix_store::ptr local_store = get_portion(i);
			local_store->reset_data();
		}
	}
}

void matrix_store::set_data(const set_operate &op)
{
	size_t num_chunks = get_num_portions();
	if (is_wide()) {
#pragma omp parallel for
		for (size_t i = 0; i < num_chunks; i++) {
			local_matrix_store::ptr local_store = get_portion(i);
			local_store->set_data(op);
		}
	}
	else {
#pragma omp parallel for
		for (size_t i = 0; i < num_chunks; i++) {
			local_matrix_store::ptr local_store = get_portion(i);
			local_store->set_data(op);
		}
	}
}

}

}
