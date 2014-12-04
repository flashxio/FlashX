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

#include "mem_dense_matrix.h"

namespace fm
{

mem_dense_matrix::ptr mem_col_dense_matrix::inner_prod(const mem_dense_matrix &m,
		const bulk_operate &left_op, const bulk_operate &right_op) const
{
	if (this->get_entry_size() != left_op.left_entry_size()
			|| m.get_entry_size() != left_op.right_entry_size()) {
		BOOST_LOG_TRIVIAL(error)
			<< "The left operator isn't compatible with input matrices";
		return mem_dense_matrix::ptr();
	}

	size_t ncol = this->get_num_cols();
	size_t nrow = this->get_num_rows();
	if (ncol != m.get_num_rows()) {
		BOOST_LOG_TRIVIAL(error) << "The matrix size doesn't match";
		return mem_dense_matrix::ptr();
	}

	mem_col_dense_matrix::ptr res = mem_col_dense_matrix::create(nrow,
			m.get_num_cols(), right_op.output_entry_size());
	res->reset_data();

#pragma omp parallel
	{
		char *tmp_res = (char *) malloc(SUB_CHUNK_SIZE * res->get_entry_size());
#pragma omp for
		for (size_t k = 0; k < nrow; k += SUB_CHUNK_SIZE) {
			sub_matrix subm(k, std::min(SUB_CHUNK_SIZE, nrow - k), 0, ncol, *this);
			for (size_t i = 0; i < ncol; i++) {
				for (size_t j = 0; j < m.get_num_cols(); j++) {
					left_op.runAE(subm.get_num_rows(), subm.get_col(i),
							m.get(i, j), tmp_res);
					char *store_col = res->get_col(j) + k * res->get_entry_size();
					right_op.runAA(subm.get_num_rows(), tmp_res, store_col,
							store_col);
				}
			}
		}
		free(tmp_res);
	}
	return res;
}

mem_dense_matrix::ptr mem_row_dense_matrix::inner_prod(const mem_dense_matrix &m,
		const bulk_operate &left_op, const bulk_operate &right_op) const
{
	if (m.store_layout() != matrix_layout_t::L_COL) {
		BOOST_LOG_TRIVIAL(error)
			<< "The layout of the right matrix has to be column matrix";
		return mem_dense_matrix::ptr();
	}
	const mem_col_dense_matrix &col_m = (const mem_col_dense_matrix &) m;

	if (this->get_entry_size() != left_op.left_entry_size()
			|| m.get_entry_size() != left_op.right_entry_size()) {
		BOOST_LOG_TRIVIAL(error)
			<< "The left operator isn't compatible with input matrices";
		return mem_dense_matrix::ptr();
	}

	size_t ncol = this->get_num_cols();
	size_t nrow = this->get_num_rows();
	if (ncol != m.get_num_rows()) {
		BOOST_LOG_TRIVIAL(error) << "The matrix size doesn't match";
		return mem_dense_matrix::ptr();
	}

	mem_row_dense_matrix::ptr res = mem_row_dense_matrix::create(nrow,
			m.get_num_cols(), right_op.output_entry_size());
	res->reset_data();

#pragma omp parallel
	{
		char *tmp_res = (char *) malloc(ncol * res->get_entry_size());
#pragma omp for
		for (size_t i = 0; i < nrow; i++) {
			for (size_t j = 0; j < m.get_num_cols(); j++) {
				left_op.runAA(ncol, get_row(i), col_m.get_col(j), tmp_res);
				right_op.runA(ncol, tmp_res, res->get(i, j));
			}
		}
		free(tmp_res);
	}
	return res;
}

}
