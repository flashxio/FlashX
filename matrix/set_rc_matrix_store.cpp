/*
 * Copyright 2017 Open Connectome Project (http://openconnecto.me)
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

#include "set_rc_matrix_store.h"
#include "local_matrix_store.h"

namespace fm
{

namespace detail
{

portion_mapply_op::const_ptr set_row_mapply_op::transpose() const
{
	auto t = rows->transpose();
	auto cols = std::dynamic_pointer_cast<const mem_col_matrix_store>(t);
	return portion_mapply_op::const_ptr(new set_col_mapply_op(idx, cols,
				get_out_num_cols(), get_out_num_rows(), get_output_type()));
}

void set_row_mapply_op::run(
		const std::vector<local_matrix_store::const_ptr> &ins,
		local_matrix_store &out) const
{
	assert(ins[0]->store_layout() == out.store_layout());
	off_t row_start = ins[0]->get_global_start_row();
	off_t row_end = ins[0]->get_global_start_row() + ins[0]->get_num_rows();
	auto start = std::lower_bound(idx->begin(), idx->end(), row_start);
	if (start == idx->end() || *start >= row_end)
		return;
	auto end = std::lower_bound(start, idx->end(), row_end);
	out.copy_from(*ins[0]);
	if (ins[0]->store_layout() == matrix_layout_t::L_ROW) {
		auto row_in
			= std::dynamic_pointer_cast<const local_row_matrix_store>(ins[0]);
		local_row_matrix_store &row_out
			= dynamic_cast<local_row_matrix_store &>(out);
		for (auto it = start; it != end; it++) {
			const char *src_row = rows->get_row(it - idx->begin());
			char *dst_row = row_out.get_row(*it - row_start);
			memcpy(dst_row, src_row, out.get_num_cols() * out.get_entry_size());
		}
	}
	else {
		std::vector<char *> first_dst_row(out.get_num_cols());
		std::vector<char *> dst_row(out.get_num_cols());
		local_col_matrix_store &col_out
			= dynamic_cast<local_col_matrix_store &>(out);
		for (size_t i = 0; i < out.get_num_cols(); i++)
			first_dst_row[i] = col_out.get_col(i);

		const scatter_gather &sg = out.get_type().get_sg();
		for (auto it = start; it != end; it++) {
			const char *src_row = rows->get_row(it - idx->begin());
			for (size_t i = 0; i < out.get_num_cols(); i++)
				dst_row[i] = first_dst_row[i]
					+ (*it - row_start) * out.get_entry_size();
			sg.scatter(src_row, dst_row);
		}
	}
}

portion_mapply_op::const_ptr set_col_mapply_op::transpose() const
{
	auto t = cols->transpose();
	auto rows = std::dynamic_pointer_cast<const mem_row_matrix_store>(t);
	return portion_mapply_op::const_ptr(new set_row_mapply_op(idx, rows,
				get_out_num_cols(), get_out_num_rows(), get_output_type()));
}

void set_col_mapply_op::run(
		const std::vector<local_matrix_store::const_ptr> &ins,
		local_matrix_store &out) const
{
	assert(ins[0]->store_layout() == out.store_layout());
	off_t col_start = ins[0]->get_global_start_col();
	off_t col_end = ins[0]->get_global_start_col() + ins[0]->get_num_cols();
	auto start = std::lower_bound(idx->begin(), idx->end(), col_start);
	if (start == idx->end() || *start >= col_end)
		return;
	auto end = std::lower_bound(start, idx->end(), col_end);
	out.copy_from(*ins[0]);
	if (ins[0]->store_layout() == matrix_layout_t::L_COL) {
		auto col_in
			= std::dynamic_pointer_cast<const local_col_matrix_store>(ins[0]);
		local_col_matrix_store &col_out
			= dynamic_cast<local_col_matrix_store &>(out);
		for (auto it = start; it != end; it++) {
			const char *src_col = cols->get_col(it - idx->begin());
			char *dst_col = col_out.get_col(*it - col_start);
			memcpy(dst_col, src_col, out.get_num_rows() * out.get_entry_size());
		}
	}
	else {
		std::vector<char *> first_dst_col(out.get_num_rows());
		std::vector<char *> dst_col(out.get_num_rows());
		local_row_matrix_store &row_out
			= dynamic_cast<local_row_matrix_store &>(out);
		for (size_t i = 0; i < out.get_num_rows(); i++)
			first_dst_col[i] = row_out.get_row(i);

		const scatter_gather &sg = out.get_type().get_sg();
		for (auto it = start; it != end; it++) {
			const char *src_col = cols->get_col(it - idx->begin());
			for (size_t i = 0; i < out.get_num_rows(); i++)
				dst_col[i] = first_dst_col[i]
					+ (*it - col_start) * out.get_entry_size();
			sg.scatter(src_col, dst_col);
		}
	}
}

}

}
