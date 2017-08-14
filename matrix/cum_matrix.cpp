/*
 * Copyright 2017
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

#include "cum_matrix.h"
#include "local_matrix_store.h"

namespace fm
{

namespace detail
{

cum_long_dim_op::cum_long_dim_op(matrix_margin margin, agg_operate::const_ptr op,
		size_t portion_size, size_t num_rows, size_t num_cols): portion_mapply_op(
			num_rows, num_cols, op->get_output_type())
{
	this->margin = margin;
	this->op = op;
	this->portion_size = portion_size;
	pthread_mutex_init(&lock, NULL);
	pthread_cond_init(&cond, NULL);
}

void cum_long_dim_op::run(const std::vector<local_matrix_store::const_ptr> &ins,
		local_matrix_store &out) const
{
	assert(ins[0]->get_num_rows() == out.get_num_rows());
	assert(ins[0]->get_num_cols() == out.get_num_cols());
	assert(ins[0]->get_type() == out.get_type());

	part_dim_t dim = margin
		== matrix_margin::MAR_COL ? part_dim_t::PART_DIM1 : part_dim_t::PART_DIM2;
	cum(*ins[0], NULL, *op, margin, dim, out);
	local_vec_store::ptr local_cum;
	if (margin == matrix_margin::MAR_COL) {
		local_cum = local_vec_store::ptr(new local_buf_vec_store(0,
					out.get_num_cols(), get_output_type(), -1));
		copy_last_row(out, *local_cum);
	}
	else {
		local_cum = local_vec_store::ptr(new local_buf_vec_store(0,
					out.get_num_rows(), get_output_type(), -1));
		copy_last_col(out, *local_cum);
	}

	// Insert the local accumulated result to the map.
	// start_row or start_col is 0.
	off_t curr_off = ins[0]->get_global_start_row() * get_out_num_cols()
		+ ins[0]->get_global_start_col() * get_out_num_rows();

	cum_long_dim_op *mutable_this = const_cast<cum_long_dim_op *>(this);
	pthread_mutex_lock(&mutable_this->lock);
	// If this is the first portion, we only need to insert the portion to
	// the map.
	mutable_this->curr_cums.insert(std::pair<off_t, cum_res_t>(curr_off,
				cum_res_t(local_cum, curr_off == 0)));

	// Compute global accumulated results.
	std::map<off_t, cum_res_t>::iterator it = mutable_this->curr_cums.begin();
	off_t prev_off = it->first;
	auto prev_cum = it->second.data;
	bool prev_used = it->second.is_used;
	bool prev_global = it->second.is_global;
	for (it++; it != curr_cums.end()
			&& prev_off + portion_size == (size_t) it->first; it++) {
		if (prev_global && !it->second.is_global) {
			auto curr = it->second.data;
			op->get_agg().runAA(curr->get_length(), curr->get_raw_arr(),
					prev_cum->get_raw_arr(), curr->get_raw_arr());
			it->second.is_global = true;
		}
		prev_off = it->first;
		prev_cum = it->second.data;
		prev_used = it->second.is_used;
		prev_global = it->second.is_global;
	}

	// We should delete the accumulated results that won't be used again.
	// The accumulated results are the ones that have been used and
	// whose next has global accumulated results.
	std::map<off_t, cum_res_t>::iterator del_it = mutable_this->curr_cums.begin();
	it = del_it;
	prev_off = it->first;
	prev_used = it->second.is_used;
	for (it++; it != curr_cums.end() && prev_used
			&& prev_off + portion_size == (size_t) it->first; it++) {
		// We count the number of accumulated results we should delete.
		if (prev_used && it->second.is_global)
			del_it = it;
		prev_off = it->first;
		prev_used = it->second.is_used;
	}
	if (del_it != curr_cums.begin())
		mutable_this->curr_cums.erase(curr_cums.begin(), del_it);

	// Find the global accumulated result for the current portion.
	local_vec_store::const_ptr global_cum;
	// TODO we should do it asynchronously.
	while (global_cum == NULL && (size_t) curr_off >= portion_size) {
		// We need to find the accumulated result from the previous portion.
		it = mutable_this->curr_cums.find(curr_off - portion_size);
		if (it != curr_cums.end() && it->second.is_global) {
			global_cum = it->second.data;
			it->second.is_used = true;
		}
		if (global_cum == NULL)
			pthread_cond_wait(&mutable_this->cond, &mutable_this->lock);
	}
	pthread_mutex_unlock(&mutable_this->lock);
	pthread_cond_broadcast(&mutable_this->cond);

	if (global_cum) {
		if (margin == matrix_margin::MAR_COL)
			mapply_rows(out, *global_cum, op->get_agg(), out);
		else
			mapply_cols(out, *global_cum, op->get_agg(), out);
	}
}

}

}
