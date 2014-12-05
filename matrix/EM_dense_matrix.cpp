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

#include <boost/foreach.hpp>

#include "io_interface.h"
#include "safs_file.h"

#include "matrix_config.h"
#include "mem_dense_matrix.h"
#include "EM_dense_matrix.h"
#include "bulk_operate.h"

namespace fm
{

bool EM_dense_matrix::exist(const std::string &name)
{
	safs::safs_file f(safs::get_sys_RAID_conf(), name);
	return f.exist();
}

void EM_col_matrix_accessor::flush()
{
	if (!fetch_reqs.empty()) {
		accessor->fetch_subvecs(fetch_reqs.data(), fetch_reqs.size());
		fetch_reqs.clear();
	}
	if (!set_reqs.empty()) {
		accessor->set_subvecs(set_reqs.data(), set_reqs.size());
		set_reqs.clear();
	}
}

class EM_subvec_accessor
{
	std::vector<fetch_vec_request> &fetch_reqs;
	std::vector<set_vec_request> &set_reqs;
	size_t start;
	size_t end;
public:
	typedef std::shared_ptr<EM_subvec_accessor> ptr;

	EM_subvec_accessor(std::vector<fetch_vec_request> &_fetch_reqs,
			std::vector<set_vec_request> &_set_reqs, size_t start,
			size_t end): fetch_reqs(_fetch_reqs), set_reqs(_set_reqs) {
		this->start = start;
		this->end = end;
	}

	void fetch_subvec(char *buf, size_t start, size_t length,
			subvec_compute::ptr compute) {
		assert(this->start + start + length <= end);
		fetch_vec_request req;
		req.buf = buf;
		req.start = this->start + start;
		req.length = length;
		req.compute = compute;
		fetch_reqs.push_back(req);
	}

	void set_subvec(const char *buf, size_t start, size_t length,
			subvec_compute::ptr compute) {
		assert(this->start + start + length <= end);
		set_vec_request req;
		req.buf = buf;
		req.start = this->start + start;
		req.length = length;
		req.compute = compute;
		set_reqs.push_back(req);
	}
};

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
	submatrix_compute::ptr compute;
	EM_col_matrix_accessor &accessor;
public:
	fetch_subcol_compute(EM_col_matrix_accessor &_accessor,
			subcol_struct::ptr subcol,
			submatrix_compute::ptr compute): accessor(_accessor) {
		this->subcol = subcol;
		this->compute = compute;
	}

	virtual void run(char *buf, size_t size) {
		subcol->count++;
		accessor.complete_req();
		if (subcol->count == subcol->subm->get_num_cols())
			compute->run(*subcol->subm);
	}
};

bool EM_col_matrix_accessor::fetch_submatrix(size_t start_row, size_t sub_nrow,
		size_t start_col, size_t sub_ncol, submatrix_compute::ptr compute)
{
	if (start_row + sub_nrow > m.get_num_rows()
			|| start_col + sub_ncol > m.get_num_cols()) {
		BOOST_LOG_TRIVIAL(error) << "fetch submatrix is out of range";
		return false;
	}

	subcol_struct::ptr subcol = subcol_struct::ptr(new subcol_struct());
	subcol->subm = mem_col_dense_matrix::create(sub_nrow, sub_ncol,
			m.get_entry_size());
	for (size_t i = 0; i < sub_ncol; i++) {
		sub_accessors[start_col + i]->fetch_subvec(subcol->subm->get_col(i),
				start_row, sub_nrow,
				subvec_compute::ptr(new fetch_subcol_compute(*this,
						subcol, compute)));
	}
	flush();
	pending_reqs += sub_ncol;
	return true;
}

/*
 * The subvec compute doesn't do anything. The only function of the class
 * is to keep a reference to the sub matrix. As long as the compute object
 * exists, the sub matrix won't be deallocated.
 */
class store_subcol_compute: public subvec_compute
{
	subcol_struct::ptr subcol;
	EM_col_matrix_accessor &accessor;
public:
	store_subcol_compute(EM_col_matrix_accessor &_accessor,
			subcol_struct::ptr subcol): accessor(_accessor) {
		this->subcol = subcol;
	}

	virtual void run(char *buf, size_t size) {
		subcol->count++;
		accessor.complete_req();
	}
};

bool EM_col_matrix_accessor::set_submatrix(size_t start_row, size_t start_col,
		mem_dense_matrix::ptr subm)
{
	if (start_row + subm->get_num_rows() > m.get_num_rows()
			|| start_col + subm->get_num_cols() > m.get_num_cols()) {
		BOOST_LOG_TRIVIAL(error) << "store submatrix is out of range";
		return false;
	}

	assert(subm->store_layout() == matrix_layout_t::L_COL);
	mem_col_dense_matrix::ptr sub_colm
		= std::static_pointer_cast<mem_col_dense_matrix>(subm);
	size_t sub_nrow = subm->get_num_rows();
	size_t sub_ncol = subm->get_num_cols();
	subcol_struct::ptr subcol = subcol_struct::ptr(new subcol_struct());
	subcol->subm = sub_colm;
	for (size_t i = 0; i < sub_ncol; i++) {
		sub_accessors[start_col + i]->set_subvec(sub_colm->get_col(i), start_row,
				sub_nrow, subvec_compute::ptr(new store_subcol_compute(*this,
						subcol)));
	}
	flush();
	pending_reqs += sub_ncol;
	return true;
}

void EM_col_matrix_accessor::wait4complete(int num)
{
	int num2wait = std::min(num, pending_reqs);
	if (num2wait == 0)
		return;
	accessor->wait4complete(num2wait);
}

void EM_col_matrix_accessor::wait4all()
{
	accessor->wait4all();
}

EM_col_matrix_accessor::EM_col_matrix_accessor(EM_col_dense_matrix &_m,
		EM_vector::ptr data): m(_m)
{
	pending_reqs = 0;
	accessor = data->create_accessor();
	sub_accessors.resize(m.get_num_cols());
	for (size_t i = 0; i < m.get_num_cols(); i++)
		sub_accessors[i] = EM_subvec_accessor::ptr(new EM_subvec_accessor(
					fetch_reqs, set_reqs, i * m.get_num_rows(),
					(i + 1) * m.get_num_rows()));
}

EM_col_matrix_accessor::~EM_col_matrix_accessor()
{
	wait4all();
	assert(pending_reqs == 0);
}

/*
 * The task to perform inner product once the submatrix is read from disks.
 */
class submatrix_inner_prod_compute: public submatrix_compute
{
	const mem_dense_matrix &m;
	EM_dense_matrix_accessor &res_m;
	const bulk_operate &left_op;
	const bulk_operate &right_op;
	size_t start_out_row;
	size_t start_out_col;
public:
	submatrix_inner_prod_compute(size_t start_row, size_t start_col,
			size_t start_out_row, size_t start_out_col,
			const bulk_operate &_left_op, const bulk_operate &_right_op,
			const mem_dense_matrix &_m, EM_dense_matrix_accessor &_res_m): submatrix_compute(
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
	EM_dense_matrix_accessor::ptr accessor = create_accessor();
	EM_dense_matrix_accessor::ptr res_accessor = res->create_accessor();
	for (size_t i = 0; i < nrow; i += COL_CHUNK_SIZE) {
		size_t sub_nrow = std::min(COL_CHUNK_SIZE, nrow - i);
		accessor->fetch_submatrix(i, sub_nrow, 0, ncol, submatrix_compute::ptr(
					new submatrix_inner_prod_compute(i, 0, i, 0, left_op,
						right_op, m, *res_accessor)));
		while (accessor->num_pending_reqs() > 8 * ncol)
			accessor->wait4complete(ncol);
		while (res_accessor->num_pending_reqs() > 2 * ncol)
			res_accessor->wait4complete(1);
	}
	accessor->wait4all();
	res_accessor->wait4all();
	return res;
}

EM_dense_matrix_accessor::ptr EM_col_dense_matrix::create_accessor()
{
	return EM_dense_matrix_accessor::ptr(new EM_col_matrix_accessor(*this,
				data));
}

void EM_col_dense_matrix::set_data(const set_operate &op)
{
	EM_dense_matrix_accessor::ptr accessor = create_accessor();
	size_t nrow = get_num_rows();
	size_t ncol = get_num_cols();
	for (size_t i = 0; i < nrow; i += COL_CHUNK_SIZE) {
		size_t chunk_size = std::min(COL_CHUNK_SIZE, nrow);
		mem_dense_matrix::ptr mem_m = mem_col_dense_matrix::create(
				chunk_size, ncol, get_entry_size());
		mem_m->set_data(op);
		accessor->set_submatrix(i, 0, mem_m);
		while ((size_t) accessor->num_pending_reqs() > 2 * ncol)
			accessor->wait4complete(1);
	}
	accessor->wait4all();
}

}
