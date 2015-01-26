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

#include "dense_matrix.h"
#include "bulk_operate.h"
#include "mem_dense_matrix.h"
#include "EM_dense_matrix.h"
#include "generic_type.h"

namespace fm
{

bool dense_matrix::verify_inner_prod(const dense_matrix &m,
		const bulk_operate &left_op, const bulk_operate &right_op) const
{
	if (this->get_entry_size() != left_op.left_entry_size()
			|| m.get_entry_size() != left_op.right_entry_size()) {
		BOOST_LOG_TRIVIAL(error)
			<< "The left operator isn't compatible with input matrices";
		return false;
	}

	if (left_op.output_entry_size() != right_op.left_entry_size()) {
		BOOST_LOG_TRIVIAL(error)
			<< "The type of the left operator doesn't match the right operator";
		return false;
	}

	if (right_op.left_entry_size() != right_op.right_entry_size()
			|| right_op.left_entry_size() != right_op.output_entry_size()) {
		BOOST_LOG_TRIVIAL(error)
			<< "The input and output of the right operator has different types";
		return false;
	}

	if (get_num_cols() != m.get_num_rows()) {
		BOOST_LOG_TRIVIAL(error) << "The matrix size doesn't match";
		return false;
	}
	return true;
}

bool dense_matrix::verify_aggregate(const bulk_operate &op,
		scalar_type &res) const
{
	if (op.left_entry_size() != op.right_entry_size()
			|| op.left_entry_size() != op.output_entry_size()) {
		BOOST_LOG_TRIVIAL(error)
			<< "The input and output type of the operator is different";
		return false;
	}

	if (this->get_entry_size() != op.left_entry_size()) {
		BOOST_LOG_TRIVIAL(error)
			<< "The matrix entry size is different from the operator";
		return false;
	}

	if (res.get_size() != op.left_entry_size()) {
		BOOST_LOG_TRIVIAL(error)
			<< "The result entry size is different from the operator";
		return false;
	}
	return true;
}

bool dense_matrix::verify_mapply2(const dense_matrix &m,
			const bulk_operate &op) const
{
	if (this->get_num_rows() != m.get_num_rows()
			|| this->get_num_cols() != m.get_num_cols()) {
		BOOST_LOG_TRIVIAL(error)
			<< "two matrices in mapply2 don't have the same shape";
		return false;
	}

	if (this->store_layout() != m.store_layout()) {
		BOOST_LOG_TRIVIAL(error)
			<< "two matrices in mapply2 don't have the same data layout";
		return false;
	}

	if (get_entry_size() != op.left_entry_size()
			|| m.get_entry_size() != op.right_entry_size()) {
		BOOST_LOG_TRIVIAL(error)
			<< "the element type in the matrices isn't compatible with the operator";
		return false;
	}

	return true;
}

dense_matrix::ptr dense_matrix::create(size_t nrow, size_t ncol,
		size_t entry_size, matrix_layout_t layout, bool in_mem)
{
	if (in_mem) {
		if (layout == matrix_layout_t::L_ROW)
			return mem_row_dense_matrix::create(nrow, ncol, entry_size);
		else
			return mem_col_dense_matrix::create(nrow, ncol, entry_size);
	}
	else {
		if (layout == matrix_layout_t::L_ROW) {
			fprintf(stderr, "EM row dense matrix isn't supported\n");
			return dense_matrix::ptr();
		}
		else
			return EM_col_dense_matrix::create(nrow, ncol, entry_size);
	}
}

}
