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

#include "mem_dense_matrix.h"
#include "EM_dense_matrix.h"
#include "bulk_operate.h"

namespace fm
{

void EM_col_dense_matrix::resize(size_t nrow, size_t ncol)
{
	cols.resize(ncol);
	for (size_t i = 0; i < cols.size(); i++) {
		if (cols[i] == NULL)
			cols[i] = EM_vector::ptr(new EM_vector(nrow));
		else
			cols[i]->resize(nrow);
	}
}

struct subcol_struct
{
	mem_col_dense_matrix::ptr subm;
	size_t count;

	subcol_struct() {
		count = 0;
	}

	typedef std::shared_ptr<subcol_struct> ptr;
};

/*
 * This subvec compute is to help read data from a column-wise matrix
 * on disks and store the result on a column-wise matrix in memory.
 * Once the in-memory matrix gets complete data, it invokes the computation
 * on the in-memory matrix.
 */
class fetch_subcol_compute: public subvec_compute
{
	subcol_struct::ptr subcol;
	size_t col_idx;
	submatrix_compute::ptr compute;
public:
	fetch_subcol_compute(subcol_struct::ptr subcol, size_t col_idx,
			submatrix_compute::ptr compute) {
		this->subcol = subcol;
		this->col_idx = col_idx;
		this->compute = compute;
	}

	virtual void run(char *buf, size_t size) {
		subcol->subm->set_col(buf, size, col_idx);
		subcol->count++;
		if (subcol->count == subcol->subm->get_num_cols())
			compute->run(*subcol->subm);
	}
};

bool EM_col_dense_matrix::fetch_submatrix(size_t start_row, size_t sub_nrow,
		size_t start_col, size_t sub_ncol, submatrix_compute::ptr compute) const
{
	if (start_row + sub_nrow > this->get_num_rows()
			|| start_col + sub_ncol > this->get_num_cols()) {
		BOOST_LOG_TRIVIAL(error) << "fetch submatrix is out of range";
		return false;
	}

	subcol_struct::ptr subcol = subcol_struct::ptr(new subcol_struct());
	subcol->subm = mem_col_dense_matrix::create(sub_nrow, sub_ncol,
			get_entry_size());
	for (size_t i = 0; i < sub_ncol; i++) {
		cols[start_col + i]->fetch_subvec(start_row, sub_nrow, subvec_compute::ptr(
					new fetch_subcol_compute(subcol, i, compute)));
	}
	return true;
}

/*
 * The subvec compute doesn't do anything. The only function of the class
 * is to keep a reference to the sub matrix. As long as the compute object
 * exists, the sub matrix won't be deallocated.
 */
class store_subcol_compute: public subvec_compute
{
	mem_col_dense_matrix::ptr subm;
public:
	store_subcol_compute(mem_col_dense_matrix::ptr subm) {
		this->subm = subm;
	}

	virtual void run(char *buf, size_t size) {
	}
};

bool EM_col_dense_matrix::set_submatrix(size_t start_row, size_t start_col,
		mem_dense_matrix::ptr subm)
{
	if (start_row + subm->get_num_rows() > this->get_num_rows()
			|| start_col + subm->get_num_cols() > this->get_num_cols()) {
		BOOST_LOG_TRIVIAL(error) << "store submatrix is out of range";
		return false;
	}

	assert(subm->store_layout() == matrix_layout_t::COL);
	mem_col_dense_matrix::ptr sub_colm
		= std::static_pointer_cast<mem_col_dense_matrix>(subm);
	size_t sub_nrow = subm->get_num_rows();
	size_t sub_ncol = subm->get_num_cols();
	for (size_t i = 0; i < sub_ncol; i++) {
		cols[start_col + i]->set_subvec(sub_colm->get_col(i), start_row, sub_nrow,
				subvec_compute::ptr(new store_subcol_compute(sub_colm)));
	}
	return true;
}

class submatrix_inner_prod_compute: public submatrix_compute
{
	const mem_dense_matrix &m;
	EM_dense_matrix &res_m;
	const bulk_operate &left_op;
	const bulk_operate &right_op;
	size_t start_out_row;
	size_t start_out_col;
public:
	submatrix_inner_prod_compute(size_t start_row, size_t start_col,
			size_t start_out_row, size_t start_out_col,
			const bulk_operate &_left_op, const bulk_operate &_right_op,
			const mem_dense_matrix &_m, EM_dense_matrix &_res_m): submatrix_compute(
				start_row, start_col), m(_m), res_m(_res_m), left_op(_left_op),
			right_op(_right_op) {
		this->start_out_row = start_out_row;
		this->start_out_col = start_out_col;
	}

	void run(const mem_dense_matrix &subm) {
		mem_dense_matrix::ptr sub_res = subm.inner_prod(m, left_op,
				right_op);
		res_m.set_submatrix(start_out_row, start_out_col, sub_res);
	}
};

EM_dense_matrix::ptr EM_col_dense_matrix::inner_prod(const mem_dense_matrix &m,
		const bulk_operate &left_op, const bulk_operate &right_op)
{
	size_t nrow = get_num_rows();
	size_t ncol = get_num_cols();
	if (ncol != m.get_num_rows()) {
		BOOST_LOG_TRIVIAL(error) << boost::format(
				"#row (%1%) of right matrix doesn't match #columns (%1%) of left matrix")
			% m.get_num_rows() % ncol;
		return EM_dense_matrix::ptr();
	}

	EM_col_dense_matrix::ptr res = EM_col_dense_matrix::create(
			nrow, m.get_num_cols(), right_op.output_entry_size());
	for (size_t i = 0; i < nrow; i += COL_CHUNK_SIZE) {
		size_t sub_nrow = std::min(COL_CHUNK_SIZE, nrow - i);
		fetch_submatrix(i, sub_nrow, 0, ncol, submatrix_compute::ptr(
					new submatrix_inner_prod_compute(i, 0, i, 0, left_op,
						right_op, m, *res)));
		// TODO wait for some computation complete.
	}
	wait4all();
	return res;
}

}
