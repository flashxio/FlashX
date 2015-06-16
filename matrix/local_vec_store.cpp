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

#include "log.h"
#include "local_vec_store.h"
#include "data_frame.h"
#include "mem_vv_store.h"
#include "mem_vec_store.h"
#include "matrix_store.h"

namespace fm
{

local_vec_store::ptr local_vec_store::get_portion(off_t loc, size_t size)
{
	assert(loc + size <= get_length());
	size_t entry_size = get_type().get_size();
	if (get_raw_arr())
		return local_vec_store::ptr(new local_ref_vec_store(
					get_raw_arr() + loc * entry_size, get_global_start() + loc,
					size, get_type(), get_node_id()));
	else {
		const local_vec_store *const_this = this;
		return local_vec_store::ptr(new local_cref_vec_store(
					const_this->get_raw_arr() + loc * entry_size,
					get_global_start() + loc, size, get_type(), get_node_id()));
	}
}

local_vec_store::const_ptr local_vec_store::get_portion(off_t loc,
			size_t size) const
{
	assert(get_raw_arr());
	assert(loc + size <= get_length());
	size_t entry_size = get_type().get_size();
	return local_vec_store::ptr(new local_cref_vec_store(
				get_raw_arr() + loc * entry_size, get_global_start() + loc,
				size, get_type(), get_node_id()));
}

detail::vec_store::ptr local_vec_store::sort_with_index()
{
	char *data_ptr = get_raw_arr();
	if (data_ptr == NULL) {
		BOOST_LOG_TRIVIAL(error) << "This local vec store is read-only";
		return detail::vec_store::ptr();
	}

	local_vec_store::ptr buf(new local_buf_vec_store(get_global_start(),
				get_length(), get_scalar_type<off_t>(), get_node_id()));
	// TODO we need to make is serial.
	get_type().get_sorter().sort_with_index(data_ptr,
			(off_t *) buf->get_raw_arr(), get_length(), false);
	return buf;
}

void local_vec_store::sort()
{
	char *data_ptr = get_raw_arr();
	assert(data_ptr);
	get_type().get_sorter().serial_sort(data_ptr, get_length(), false);
}

bool local_vec_store::is_sorted() const
{
	// Test if the array is sorted in the ascending order.
	return get_type().get_sorter().is_sorted(get_raw_arr(), get_length(), false);
}

data_frame::ptr local_vec_store::groupby(
		const gr_apply_operate<local_vec_store> &op, bool with_val) const
{
	const scalar_type &output_type = op.get_output_type();
	const agg_operate &find_next = get_type().get_agg_ops().get_find_next();

	assert(is_sorted());
	detail::vec_store::ptr agg;
	// TODO it might not be a good idea to create a mem_vector
	if (op.get_num_out_eles() == 1)
		agg = detail::smp_vec_store::create(0, output_type);
	else
		agg = detail::mem_vv_store::create(output_type);
	detail::smp_vec_store::ptr val;
	if (with_val)
		val = detail::smp_vec_store::create(0, get_type());

	size_t out_size;
	// If the user can predict the number of output elements, we can create
	// a buffer of the expected size.
	if (op.get_num_out_eles() > 0)
		out_size = op.get_num_out_eles();
	else
		// If the user can't, we create a small buffer.
		out_size = 16;
	local_buf_vec_store one_agg(0, out_size, output_type, -1);
	std::vector<const char *> val_locs;
	off_t init_global_start = get_global_start();
	size_t loc = 0;
	while (loc < get_length()) {
		size_t curr_length = get_length() - loc;
		const char *curr_ptr = get_raw_arr() + get_entry_size() * loc;
		size_t rel_end;
		find_next.run(curr_length, curr_ptr, &rel_end);
		local_cref_vec_store lcopy(curr_ptr, init_global_start + loc, rel_end,
				get_type(), get_node_id());
		op.run(curr_ptr, lcopy, one_agg);
		agg->append(one_agg);
		if (with_val)
			val_locs.push_back(curr_ptr);
		loc += rel_end;
	}
	if (with_val) {
		val->resize(val_locs.size());
		val->set(val_locs);
	}
	data_frame::ptr ret = data_frame::create();
	if (with_val)
		ret->add_vec("val", val);
	ret->add_vec("agg", agg);
	return ret;
}

local_vec_store::ptr local_vec_store::get(std::vector<off_t> &idxs) const
{
	local_vec_store::ptr ret(new local_buf_vec_store(-1, idxs.size(), get_type(), -1));
	std::vector<const char *> ptrs(idxs.size());
	for (size_t i = 0; i < idxs.size(); i++)
		ptrs[i] = this->get(idxs[i]);
	get_type().get_sg().gather(ptrs, ret->get_raw_arr());
	return ret;
}

detail::matrix_store::const_ptr local_vec_store::conv2mat(size_t nrow,
		size_t ncol, bool byrow) const
{
	BOOST_LOG_TRIVIAL(error)
		<< "can't convert a local vector to a matrix";
	return detail::matrix_store::ptr();
}

bool local_ref_vec_store::resize(size_t new_length)
{
	BOOST_LOG_TRIVIAL(error)
		<< "can't resize a local reference vector store";
	assert(0);
	return false;
}

bool local_cref_vec_store::resize(size_t new_length)
{
	BOOST_LOG_TRIVIAL(error)
		<< "can't resize a local reference vector store";
	assert(0);
	return false;
}

bool local_buf_vec_store::resize(size_t new_length)
{
	if (get_length() < new_length) {
		arr.expand(new_length * get_type().get_size());
		set_data(arr.get_raw(), arr.get_raw());
	}
	detail::vec_store::resize(new_length);
	return true;
}

}
