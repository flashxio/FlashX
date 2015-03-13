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

#if defined(_OPENMP)
#include <omp.h>
#endif

#include <atomic>

#include <boost/format.hpp>

#include "log.h"

#include "mem_dense_matrix.h"
#include "generic_type.h"

namespace fm
{

const size_t SUB_CHUNK_SIZE = 1024;

/*
 * This class contains the information of a submatrix in the original matrix.
 * It is mainly used in inner product.
 */
class sub_matrix_info
{
	size_t start_row;
	size_t start_col;
	size_t nrow;
	size_t ncol;
public:
	sub_matrix_info(size_t start_row, size_t nrow, size_t start_col,
			size_t ncol) {
		this->start_row = start_row;
		this->start_col = start_col;
		this->nrow = nrow;
		this->ncol = ncol;
	}

	size_t get_num_rows() const {
		return nrow;
	}

	size_t get_num_cols() const {
		return ncol;
	}

	size_t get_start_row() const {
		return start_row;
	}

	size_t get_start_col() const {
		return start_col;
	}
};

/*
 * This class contains the information of a submatrix in the a column-wise matrix.
 * It is mainly used in inner product.
 */
class sub_col_matrix_info: public sub_matrix_info
{
	const mem_col_dense_matrix &m;
public:
	sub_col_matrix_info(size_t start_row, size_t nrow, size_t start_col,
			size_t ncol, const mem_col_dense_matrix &_m): sub_matrix_info(
				start_row, nrow, start_col, ncol), m(_m) {
		assert(start_row + nrow <= m.get_num_rows());
		assert(start_col + ncol <= m.get_num_cols());
	}

	const char *get_col(size_t col) const {
		return m.get_col(get_start_col() + col)
			+ get_start_row() * m.get_entry_size();
	}
};

/*
 * This class contains the information of a submatrix in the a row-wise matrix.
 * It is mainly used in inner product.
 */
class sub_row_matrix_info: public sub_matrix_info
{
	const mem_row_dense_matrix &m;
public:
	sub_row_matrix_info(size_t start_row, size_t nrow, size_t start_col,
			size_t ncol, const mem_row_dense_matrix &_m): sub_matrix_info(
				start_row, nrow, start_col, ncol), m(_m) {
		assert(start_row + nrow <= m.get_num_rows());
		assert(start_col + ncol <= m.get_num_cols());
	}

	const char *get_row(size_t row) const {
		return m.get_row(get_start_row() + row) + get_start_col() * m.get_entry_size();
	}
};

/*
 * This class defines a submatrix from a column-wise matrix.
 * Users can use the submatrix as if it is a normal matrix.
 */
class mem_sub_col_dense_matrix: public mem_col_dense_matrix
{
	// The data buffer is referenced in the parent class.
	// but this class also needs to access the data buffer.
	std::shared_ptr<char> data;
	std::vector<off_t> orig_col_idxs;

	mem_sub_col_dense_matrix(size_t nrow, size_t ncol, const scalar_type &type,
			std::shared_ptr<char> data,
			const std::vector<off_t> &col_idxs): mem_col_dense_matrix(nrow, ncol,
				type, data) {
		this->orig_col_idxs = col_idxs;
		this->data = data;
	}
public:
	typedef std::shared_ptr<mem_sub_col_dense_matrix> ptr;

	static ptr create(size_t nrow, size_t ncol, const scalar_type &type,
			std::shared_ptr<char> data, const std::vector<off_t> &col_idxs) {
		return ptr(new mem_sub_col_dense_matrix(nrow, ncol, type, data,
					col_idxs));
	}

	~mem_sub_col_dense_matrix() {
	}

	virtual dense_matrix::ptr shallow_copy() const {
		// TODO
		BOOST_LOG_TRIVIAL(error)
			<< "shallow_copy() isn't supported in a sub_col_matrix";
		return dense_matrix::ptr();
	}

	virtual dense_matrix::ptr deep_copy() const {
		// TODO
		BOOST_LOG_TRIVIAL(error)
			<< "deep_copy() isn't supported in a sub_col_matrix";
		return dense_matrix::ptr();
	}

	virtual dense_matrix::ptr conv2(size_t nrow, size_t ncol, bool byrow) const {
		// TODO
		BOOST_LOG_TRIVIAL(error) << "conv2() isn't supported in a sub_col_matrix";
		return dense_matrix::ptr();
	}

	virtual void reset_data() {
		// TODO
		BOOST_LOG_TRIVIAL(error)
			<< "reset_data() isn't supported in a sub_col_matrix";
	}
	virtual void set_data(const set_operate &op) {
		// TODO
		BOOST_LOG_TRIVIAL(error)
			<< "set_data() isn't supported in a sub_col_matrix";
	}
	virtual void serial_reset_data() {
		// TODO
		BOOST_LOG_TRIVIAL(error)
			<< "serial_reset_data() isn't supported in a sub_col_matrix";
	}
	virtual void serial_set_data(const set_operate &op) {
		// TODO
		BOOST_LOG_TRIVIAL(error)
			<< "serial_set_data() isn't supported in a sub_col_matrix";
	}

	virtual bool set_cols(const mem_col_dense_matrix &m,
			const std::vector<off_t> &idxs) {
		// TODO
		BOOST_LOG_TRIVIAL(error)
			<< "set_cols() isn't supported in a sub_col_matrix";
		return false;
	}
	virtual bool set_col(const char *buf, size_t size, size_t col) {
		// TODO
		BOOST_LOG_TRIVIAL(error)
			<< "set_col() isn't supported in a sub_col_matrix";
		return false;
	}

	virtual dense_matrix::ptr get_cols(const std::vector<off_t> &idxs) const {
		// TODO
		BOOST_LOG_TRIVIAL(error)
			<< "get_cols() isn't supported in a sub_col_matrix";
		return false;
	}

	virtual dense_matrix::ptr apply(apply_margin margin,
			const arr_apply_operate &op) const {
		// TODO
		BOOST_LOG_TRIVIAL(error) << "apply() isn't supported in a sub_col_matrix";
		return dense_matrix::ptr();
	}

	virtual bool aggregate(const bulk_operate &op, scalar_variable &res) const;
	virtual dense_matrix::ptr mapply2(const dense_matrix &m,
			const bulk_operate &op) const;
	virtual dense_matrix::ptr sapply(const bulk_uoperate &op) const;

	virtual dense_matrix::ptr transpose() const;

	virtual char *get_col(size_t col) {
		off_t orig_col = orig_col_idxs[col];
		return data.get() + orig_col * get_num_rows() * get_entry_size();
	}

	virtual const char *get_col(size_t col) const {
		off_t orig_col = orig_col_idxs[col];
		return data.get() + orig_col * get_num_rows() * get_entry_size();
	}
};

/*
 * This class defines a submatrix from a row-wise matrix.
 * Users can use the submatrix as if it is a normal matrix.
 */
class mem_sub_row_dense_matrix: public mem_row_dense_matrix
{
	// The data buffer is referenced in the parent class.
	// but this class also needs to access the data buffer.
	std::shared_ptr<char> data;
	std::vector<off_t> orig_row_idxs;

	mem_sub_row_dense_matrix(size_t nrow, size_t ncol, const scalar_type &type,
			std::shared_ptr<char> data,
			const std::vector<off_t> &row_idxs): mem_row_dense_matrix(nrow, ncol,
				type, data) {
		this->orig_row_idxs = row_idxs;
		this->data = data;
	}
public:
	typedef std::shared_ptr<mem_sub_row_dense_matrix> ptr;

	static ptr create(size_t nrow, size_t ncol, const scalar_type &type,
			std::shared_ptr<char> data, const std::vector<off_t> &row_idxs) {
		return ptr(new mem_sub_row_dense_matrix(nrow, ncol, type, data,
					row_idxs));
	}

	~mem_sub_row_dense_matrix() {
	}

	virtual dense_matrix::ptr shallow_copy() const {
		// TODO
		BOOST_LOG_TRIVIAL(error)
			<< "shallow_copy() isn't supported in a sub_col_matrix";
		return dense_matrix::ptr();
	}

	virtual dense_matrix::ptr deep_copy() const {
		// TODO
		BOOST_LOG_TRIVIAL(error)
			<< "deep_copy() isn't supported in a sub_col_matrix";
		return dense_matrix::ptr();
	}

	virtual dense_matrix::ptr conv2(size_t nrow, size_t ncol, bool byrow) const {
		// TODO
		BOOST_LOG_TRIVIAL(error)
			<< "conv2() isn't supported in a sub_col_matrix";
		return dense_matrix::ptr();
	}

	virtual void reset_data() {
		// TODO
		BOOST_LOG_TRIVIAL(error)
			<< "reset_data() isn't supported in a sub_col_matrix";
	}
	virtual void set_data(const set_operate &op) {
		// TODO
		BOOST_LOG_TRIVIAL(error)
			<< "set_data() isn't supported in a sub_col_matrix";
	}
	virtual void serial_reset_data() {
		// TODO
		BOOST_LOG_TRIVIAL(error)
			<< "serial_reset_data() isn't supported in a sub_col_matrix";
	}
	virtual void serial_set_data(const set_operate &op) {
		// TODO
		BOOST_LOG_TRIVIAL(error)
			<< "serial_set_data() isn't supported in a sub_col_matrix";
	}

	virtual dense_matrix::ptr apply(apply_margin margin,
			const arr_apply_operate &op) const {
		// TODO
		BOOST_LOG_TRIVIAL(error) << "apply() isn't supported in a sub_row_matrix";
		return dense_matrix::ptr();
	}

	virtual bool aggregate(const bulk_operate &op, scalar_variable &res) const;
	virtual dense_matrix::ptr mapply2(const dense_matrix &m,
			const bulk_operate &op) const;
	virtual dense_matrix::ptr sapply(const bulk_uoperate &op) const;

	virtual dense_matrix::ptr transpose() const;

	char *get_row(size_t row) {
		size_t orig_row = orig_row_idxs[row];
		return data.get() + orig_row * get_num_cols() * get_entry_size();
	}

	const char *get_row(size_t row) const {
		size_t orig_row = orig_row_idxs[row];
		return data.get() + orig_row * get_num_cols() * get_entry_size();
	}
};

bool mem_sub_col_dense_matrix::aggregate(const bulk_operate &op,
		scalar_variable &res) const
{
	if (!verify_aggregate(op, res))
		return false;

	size_t ncol = this->get_num_cols();
	size_t nrow = this->get_num_rows();
	char *raw_arr = (char *) malloc(res.get_size() * ncol);
	// TODO parallel
	for (size_t i = 0; i < ncol; i++)
		op.runA(nrow, get_col(i), raw_arr + res.get_size() * i);
	char raw_res[res.get_size()];
	op.runA(ncol, raw_arr, raw_res);
	free(raw_arr);
	res.set_raw(raw_res, res.get_size());
	return true;
}

dense_matrix::ptr mem_sub_col_dense_matrix::mapply2(const dense_matrix &m,
		const bulk_operate &op) const
{
	// The same shape and the same data layout.
	if (!verify_mapply2(m, op))
		return dense_matrix::ptr();

	const mem_col_dense_matrix &col_m = (const mem_col_dense_matrix &) m;
	size_t ncol = this->get_num_cols();
	size_t nrow = this->get_num_rows();
	mem_col_dense_matrix::ptr res = mem_col_dense_matrix::create(nrow, ncol,
			op.get_output_type());
	// TODO parallel
	for (size_t i = 0; i < ncol; i++)
		op.runAA(nrow, get_col(i), col_m.get_col(i), res->get_col(i));
	return res;
}


dense_matrix::ptr mem_sub_col_dense_matrix::sapply(const bulk_uoperate &op) const
{
	// TODO parallel
	size_t ncol = get_num_cols();
	size_t nrow = get_num_rows();
	mem_col_dense_matrix::ptr res = mem_col_dense_matrix::create(nrow, ncol,
			op.get_output_type());
	for (size_t i = 0; i < ncol; i++)
		op.runA(nrow, get_col(i), res->get_col(i));
	return res;
}

dense_matrix::ptr mem_sub_col_dense_matrix::transpose() const
{
	return mem_sub_row_dense_matrix::create(get_num_cols(), get_num_rows(),
			get_type(), data, orig_col_idxs);
}

bool mem_sub_row_dense_matrix::aggregate(const bulk_operate &op,
		scalar_variable &res) const
{
	if (!verify_aggregate(op, res))
		return false;

	size_t ncol = this->get_num_cols();
	size_t nrow = this->get_num_rows();
	char *raw_arr = (char *) malloc(res.get_size() * nrow);
	// TODO parallel
	for (size_t i = 0; i < nrow; i++)
		op.runA(ncol, get_row(i), raw_arr + res.get_size() * i);
	char raw_res[res.get_size()];
	op.runA(nrow, raw_arr, raw_res);
	free(raw_arr);
	res.set_raw(raw_res, res.get_size());
	return true;
}

dense_matrix::ptr mem_sub_row_dense_matrix::mapply2(const dense_matrix &m,
		const bulk_operate &op) const
{
	// The same shape and the same data layout.
	if (!verify_mapply2(m, op))
		return dense_matrix::ptr();

	const mem_row_dense_matrix &row_m = (const mem_row_dense_matrix &) m;
	size_t ncol = this->get_num_cols();
	size_t nrow = this->get_num_rows();
	mem_row_dense_matrix::ptr res = mem_row_dense_matrix::create(nrow, ncol,
			op.get_output_type());
	// TODO parallel
	for (size_t i = 0; i < nrow; i++)
		op.runAA(ncol, get_row(i), row_m.get_row(i), res->get_row(i));
	return res;
}

dense_matrix::ptr mem_sub_row_dense_matrix::sapply(const bulk_uoperate &op) const
{
	// TODO parallel
	size_t ncol = get_num_cols();
	size_t nrow = get_num_rows();
	mem_row_dense_matrix::ptr res = mem_row_dense_matrix::create(nrow, ncol,
			op.get_output_type());
	for (size_t i = 0; i < nrow; i++)
		op.runA(ncol, get_row(i), res->get_row(i));
	return res;
}

dense_matrix::ptr mem_sub_row_dense_matrix::transpose() const
{
	return mem_sub_col_dense_matrix::create(get_num_cols(),
			get_num_rows(), get_type(), data, orig_row_idxs);
}

bool mem_dense_matrix::verify_inner_prod(const dense_matrix &m,
		const bulk_operate &left_op, const bulk_operate &right_op) const
{
	if (!m.is_in_mem()) {
		BOOST_LOG_TRIVIAL(error) << "The right matrix isn't in memory";
		return false;
	}
	return dense_matrix::verify_inner_prod(m, left_op, right_op);
}

mem_dense_matrix::ptr mem_dense_matrix::cast(dense_matrix::ptr m)
{
	if (!m->is_in_mem()) {
		BOOST_LOG_TRIVIAL(error)
			<< "Can't cast an EM matrix to mem_dense_matrix";
		return mem_dense_matrix::ptr();
	}
	return std::static_pointer_cast<mem_dense_matrix>(m);
}

mem_row_dense_matrix::ptr mem_row_dense_matrix::cast(dense_matrix::ptr m)
{
	if (!m->is_in_mem()) {
		BOOST_LOG_TRIVIAL(error)
			<< "can't cast an EM matrix to mem_row_dense_matrix";
		return mem_row_dense_matrix::ptr();
	}
	if (m->store_layout() != matrix_layout_t::L_ROW) {
		BOOST_LOG_TRIVIAL(error)
			<< "the matrix to be cast isn't row-wise";
		return mem_row_dense_matrix::ptr();
	}

	return std::static_pointer_cast<mem_row_dense_matrix>(m);
}

mem_row_dense_matrix::ptr mem_row_dense_matrix::cast(mem_dense_matrix::ptr m)
{
	if (m->store_layout() != matrix_layout_t::L_ROW) {
		BOOST_LOG_TRIVIAL(error)
			<< "the matrix to be cast isn't row-wise";
		return mem_row_dense_matrix::ptr();
	}

	return std::static_pointer_cast<mem_row_dense_matrix>(m);
}

mem_col_dense_matrix::ptr mem_col_dense_matrix::cast(dense_matrix::ptr m)
{
	if (!m->is_in_mem()) {
		BOOST_LOG_TRIVIAL(error)
			<< "can't cast an EM matrix to mem_col_dense_matrix";
		return mem_col_dense_matrix::ptr();
	}
	if (m->store_layout() != matrix_layout_t::L_COL) {
		BOOST_LOG_TRIVIAL(error)
			<< "the matrix to be cast isn't col-wise";
		return mem_col_dense_matrix::ptr();
	}

	return std::static_pointer_cast<mem_col_dense_matrix>(m);
}

mem_col_dense_matrix::ptr mem_col_dense_matrix::cast(mem_dense_matrix::ptr m)
{
	if (m->store_layout() != matrix_layout_t::L_COL) {
		BOOST_LOG_TRIVIAL(error)
			<< "the matrix to be cast isn't col-wise";
		return mem_col_dense_matrix::ptr();
	}

	return std::static_pointer_cast<mem_col_dense_matrix>(m);
}

dense_matrix::ptr mem_col_dense_matrix::shallow_copy() const
{
	// The data array is read-only. It's safe to have two matrices reference
	// the same data array.
	return dense_matrix::ptr(new mem_col_dense_matrix(get_num_rows(),
				get_num_cols(), get_type(), data));
}

dense_matrix::ptr mem_col_dense_matrix::deep_copy() const
{
	size_t num_bytes = get_num_rows() * get_num_cols() * get_entry_size();
	std::shared_ptr<char> new_data = std::shared_ptr<char>(
			(char *) memalign(PAGE_SIZE, num_bytes), deleter());
	memcpy(new_data.get(), data.get(), num_bytes);
	return dense_matrix::ptr(new mem_col_dense_matrix(get_num_rows(),
				get_num_cols(), get_type(), new_data));
}

dense_matrix::ptr mem_col_dense_matrix::conv2(size_t nrow, size_t ncol,
		bool byrow) const
{
	if (nrow * ncol > get_num_rows() * get_num_cols()) {
		BOOST_LOG_TRIVIAL(error)
			<< "can't convert to a matrix larger than the original matrix";
		return dense_matrix::ptr();
	}

	if (!byrow) {
		// The new matrix has the same data layout as the original matrix.
		// so we can simply clone the original matrix.
		dense_matrix::ptr new_mat = shallow_copy();
		new_mat->resize(nrow, ncol);
		return new_mat;
	}
	else
		return dense_matrix::ptr(new mem_row_dense_matrix(nrow, ncol,
					get_type(), data));
}

dense_matrix::ptr mem_col_dense_matrix::transpose() const
{
	// When we transpose a column-wise matrix, the matrix becomes row-wise,
	// and no data needs to be copied.
	return dense_matrix::ptr(new mem_row_dense_matrix(get_num_cols(),
				get_num_rows(), get_type(), data));
}

void mem_col_dense_matrix::serial_reset_data()
{
	size_t tot_bytes = get_num_rows() * get_num_cols() * get_entry_size();
	memset(data.get(), 0, tot_bytes);
}

void mem_col_dense_matrix::serial_set_data(const set_operate &op)
{
	size_t ncol = get_num_cols();
	size_t nrow = get_num_rows();
	for (size_t i = 0; i < ncol; i++)
		op.set(get_col(i), nrow, 0, i);
}

void mem_col_dense_matrix::reset_data()
{
	size_t tot_bytes = get_num_rows() * get_num_cols() * get_entry_size();
#pragma omp parallel for
	for (size_t i = 0; i < tot_bytes; i += PAGE_SIZE)
		memset(data.get() + i, 0, std::min(tot_bytes - i, (size_t) PAGE_SIZE));
}

void mem_col_dense_matrix::set_data(const set_operate &op)
{
	size_t ncol = get_num_cols();
	size_t nrow = get_num_rows();
#pragma omp parallel for
	for (size_t i = 0; i < ncol; i++)
		op.set(get_col(i), nrow, 0, i);
}

dense_matrix::ptr mem_col_dense_matrix::serial_inner_prod(const dense_matrix &m,
		const bulk_operate &left_op, const bulk_operate &right_op) const
{
	if (!verify_inner_prod(m, left_op, right_op))
		return dense_matrix::ptr();

	size_t ncol = this->get_num_cols();
	size_t nrow = this->get_num_rows();
	assert(nrow > ncol);
	// TODO we need to determine the layout of the output matrix smartly.
	mem_col_dense_matrix::ptr res = mem_col_dense_matrix::create(nrow,
			m.get_num_cols(), right_op.get_output_type());
	res->serial_reset_data();
	const mem_dense_matrix &mem_m = (const mem_dense_matrix &) m;

	char *tmp_res = (char *) malloc(SUB_CHUNK_SIZE * res->get_entry_size());
	for (size_t k = 0; k < nrow; k += SUB_CHUNK_SIZE) {
		sub_col_matrix_info subm(k, std::min(SUB_CHUNK_SIZE, nrow - k), 0, ncol, *this);
		for (size_t i = 0; i < ncol; i++) {
			for (size_t j = 0; j < m.get_num_cols(); j++) {
				left_op.runAE(subm.get_num_rows(), subm.get_col(i),
						mem_m.get(i, j), tmp_res);
				char *store_col = res->get_col(j) + k * res->get_entry_size();
				right_op.runAA(subm.get_num_rows(), tmp_res, store_col,
						store_col);
			}
		}
	}
	free(tmp_res);
	return std::static_pointer_cast<dense_matrix>(res);
}

dense_matrix::ptr mem_col_dense_matrix::inner_prod(const dense_matrix &m,
		const bulk_operate &left_op, const bulk_operate &right_op) const
{
	if (!verify_inner_prod(m, left_op, right_op))
		return dense_matrix::ptr();

	size_t ncol = this->get_num_cols();
	size_t nrow = this->get_num_rows();
	assert(nrow > ncol);
	mem_col_dense_matrix::ptr res = mem_col_dense_matrix::create(nrow,
			m.get_num_cols(), right_op.get_output_type());
	res->reset_data();
	const mem_dense_matrix &mem_m = (const mem_dense_matrix &) m;

#pragma omp parallel
	{
		char *tmp_res = (char *) malloc(SUB_CHUNK_SIZE * res->get_entry_size());
#pragma omp for
		for (size_t k = 0; k < nrow; k += SUB_CHUNK_SIZE) {
			sub_col_matrix_info subm(k, std::min(SUB_CHUNK_SIZE, nrow - k), 0, ncol, *this);
			for (size_t i = 0; i < ncol; i++) {
				for (size_t j = 0; j < m.get_num_cols(); j++) {
					left_op.runAE(subm.get_num_rows(), subm.get_col(i),
							mem_m.get(i, j), tmp_res);
					char *store_col = res->get_col(j) + k * res->get_entry_size();
					right_op.runAA(subm.get_num_rows(), tmp_res, store_col,
							store_col);
				}
			}
		}
		free(tmp_res);
	}
	return std::static_pointer_cast<dense_matrix>(res);
}

bool mem_col_dense_matrix::aggregate(const bulk_operate &op,
		scalar_variable &res) const
{
	if (!verify_aggregate(op, res))
		return false;

	size_t ncol = this->get_num_cols();
	size_t nrow = this->get_num_rows();
	char raw_arr[res.get_size()];
	// TODO parallel
	op.runA(nrow * ncol, data.get(), raw_arr);
	res.set_raw(raw_arr, res.get_size());
	return true;
}

dense_matrix::ptr mem_col_dense_matrix::mapply2(const dense_matrix &m,
		const bulk_operate &op) const
{
	// The same shape and the same data layout.
	if (!verify_mapply2(m, op))
		return dense_matrix::ptr();

	// TODO parallel
	const mem_col_dense_matrix &col_m = (const mem_col_dense_matrix &) m;
	size_t ncol = get_num_cols();
	size_t nrow = get_num_rows();
	mem_col_dense_matrix::ptr res = mem_col_dense_matrix::create(nrow, ncol,
			op.get_output_type());
	op.runAA(ncol * nrow, data.get(), col_m.data.get(), res->data.get());
	return res;
}

dense_matrix::ptr mem_col_dense_matrix::sapply(const bulk_uoperate &op) const
{
	// TODO parallel
	size_t ncol = get_num_cols();
	size_t nrow = get_num_rows();
	mem_col_dense_matrix::ptr res = mem_col_dense_matrix::create(nrow, ncol,
			op.get_output_type());
	op.runA(ncol * nrow, data.get(), res->data.get());
	return res;
}

dense_matrix::ptr mem_col_dense_matrix::apply(apply_margin margin,
		const arr_apply_operate &op) const
{
	if (!verify_apply(margin, op))
		return dense_matrix::ptr();

	// TODO parallel
	size_t ncol = get_num_cols();
	size_t nrow = get_num_rows();
	// Each operation runs on a row
	if (margin == apply_margin::MAR_ROW) {
		size_t out_nrow = nrow;
		size_t out_ncol = op.get_num_out_eles();
		mem_row_dense_matrix::ptr res = mem_row_dense_matrix::create(out_nrow,
				out_ncol, op.get_output_type());
		// TODO this might be a very large array.
		char *tmp_arr = (char *) malloc(ncol * get_entry_size());
		for (size_t i = 0; i < nrow; i++) {
			get_row(i, tmp_arr);
			op.run(tmp_arr, ncol, res->get_row(i));
		}
		free(tmp_arr);
		return res;
	}
	// Each operation runs on a column
	else {
		size_t out_nrow = op.get_num_out_eles();
		size_t out_ncol = ncol;
		mem_col_dense_matrix::ptr res = mem_col_dense_matrix::create(out_nrow,
				out_ncol, op.get_output_type());
		for (size_t i = 0; i < ncol; i++)
			op.run(get_col(i), nrow, res->get_col(i));
		return res;
	}
}

bool mem_col_dense_matrix::set_cols(const mem_col_dense_matrix &m,
		const std::vector<off_t> &idxs)
{
	if (m.get_num_cols() != idxs.size()) {
		BOOST_LOG_TRIVIAL(error)
			<< "# vectors is different from # column indices";
		return false;
	}

	for (size_t i = 0; i < idxs.size(); i++) {
		// If idxs[i] is negative, the unsigned value will be very large.
		if ((size_t) idxs[i] >= get_num_cols()) {
			BOOST_LOG_TRIVIAL(error)
				<< "a column index is out of bounds\n";
			return false;
		}
	}

	if (get_num_rows() != m.get_num_rows()) {
		BOOST_LOG_TRIVIAL(error)
			<< "The length of the vectors is different from #rows of the matrix";
		return false;
	}

	if (get_type() != m.get_type()) {
		BOOST_LOG_TRIVIAL(error)
			<< "The matrix and the vectors have different types";
		return false;
	}

	size_t num_cols = idxs.size();
	for (size_t i = 0; i < num_cols; i++)
		set_col(m.get_col(i), m.get_num_rows() * m.get_entry_size(), idxs[i]);
	return true;
}

bool mem_col_dense_matrix::set_col(const char *buf, size_t size, size_t col)
{
	if (size != get_entry_size() * get_num_rows()) {
		BOOST_LOG_TRIVIAL(error)
			<< "set_col: has a different column length";
		return false;
	}

	memcpy(get_col(col), buf, size);
	return true;
}

dense_matrix::ptr mem_col_dense_matrix::get_cols(const std::vector<off_t> &idxs) const
{
	for (size_t i = 0; i < idxs.size(); i++) {
		if ((size_t) idxs[i] >= get_num_cols()) {
			BOOST_LOG_TRIVIAL(error)
				<< "a column index is out of bounds\n";
			return dense_matrix::ptr();
		}
	}

	return mem_sub_col_dense_matrix::create(get_num_rows(),
			idxs.size(), get_type(), data, idxs);
}

bool mem_col_dense_matrix::write2file(const std::string &file_name) const
{
	FILE *f = fopen(file_name.c_str(), "w");
	if (f == NULL) {
		BOOST_LOG_TRIVIAL(error)
			<< boost::format("can't open %1%: %2%") % file_name % strerror(errno);
		return false;
	}
	if (!write_header(f))
		return false;

	size_t ncol = get_num_cols();
	size_t col_size = get_num_rows() * get_entry_size();
	for (size_t i = 0; i < ncol; i++) {
		const char *col = get_col(i);
		size_t ret = fwrite(col, col_size, 1, f);
		if (ret == 0) {
			BOOST_LOG_TRIVIAL(error)
				<< boost::format("can't write to %1%: %2%")
				% file_name % strerror(errno);
			return false;
		}
	}
	fclose(f);
	return true;
}

mem_col_dense_matrix::ptr mem_col_dense_matrix::create(size_t nrow, size_t ncol,
		const scalar_type &type, FILE *f)
{
	size_t mat_size = nrow * ncol * type.get_size();
	std::shared_ptr<char> data = std::shared_ptr<char>((char *) memalign(
				PAGE_SIZE, mat_size), deleter());
	if (data == NULL) {
		BOOST_LOG_TRIVIAL(error) << "can't allocate memory for the matrix";
		return mem_col_dense_matrix::ptr();
	}
	size_t ret = fread(data.get(), mat_size, 1, f);
	if (ret == 0) {
		BOOST_LOG_TRIVIAL(error)
			<< boost::format("can't read %1% bytes from the file") % mat_size;
		return mem_col_dense_matrix::ptr();
	}

	return mem_col_dense_matrix::ptr(new mem_col_dense_matrix(nrow, ncol,
				type, data));
}

void mem_col_dense_matrix::get_row(size_t row_idx, char *arr) const
{
	size_t ncol = get_num_cols();
	size_t entry_size = get_entry_size();
	for (size_t i = 0; i < ncol; i++)
		memcpy(arr + i * entry_size, get_col(i) + row_idx * entry_size,
				entry_size);
}

dense_matrix::ptr mem_row_dense_matrix::shallow_copy() const
{
	// The data array is read-only. It's safe to have two matrices reference
	// the same data array.
	return dense_matrix::ptr(new mem_row_dense_matrix(get_num_rows(),
				get_num_cols(), get_type(), data));
}

dense_matrix::ptr mem_row_dense_matrix::deep_copy() const
{
	size_t num_bytes = get_num_rows() * get_num_cols() * get_entry_size();
	std::shared_ptr<char> new_data = std::shared_ptr<char>(
			(char *) memalign(PAGE_SIZE, num_bytes), deleter());
	memcpy(new_data.get(), data.get(), num_bytes);
	return dense_matrix::ptr(new mem_row_dense_matrix(get_num_rows(),
				get_num_cols(), get_type(), new_data));
}

dense_matrix::ptr mem_row_dense_matrix::conv2(size_t nrow, size_t ncol,
		bool byrow) const
{
	if (nrow * ncol > get_num_rows() * get_num_cols()) {
		BOOST_LOG_TRIVIAL(error)
			<< "can't convert to a matrix larger than the original matrix";
		return dense_matrix::ptr();
	}

	if (byrow) {
		// The new matrix has the same data layout as the original matrix.
		// so we can simply clone the original matrix.
		dense_matrix::ptr new_mat = shallow_copy();
		new_mat->resize(nrow, ncol);
		return new_mat;
	}
	else
		return dense_matrix::ptr(new mem_col_dense_matrix(nrow, ncol,
					get_type(), data));
}

dense_matrix::ptr mem_row_dense_matrix::transpose() const
{
	// When we transpose a row-wise matrix, the matrix becomes column-wise,
	// and no data needs to be copied.
	return dense_matrix::ptr(new mem_col_dense_matrix(get_num_cols(),
				get_num_rows(), get_type(), data));
}

void mem_row_dense_matrix::serial_reset_data()
{
	size_t tot_bytes = get_num_rows() * get_num_cols() * get_entry_size();
	memset(data.get(), 0, tot_bytes);
}

void mem_row_dense_matrix::serial_set_data(const set_operate &op)
{
	size_t ncol = get_num_cols();
	size_t nrow = get_num_rows();
	for (size_t i = 0; i < nrow; i++)
		op.set(get_row(i), ncol, i, 0);
}

void mem_row_dense_matrix::reset_data()
{
	size_t tot_bytes = get_num_rows() * get_num_cols() * get_entry_size();
#pragma omp parallel for
	for (size_t i = 0; i < tot_bytes; i += PAGE_SIZE)
		memset(data.get() + i, 0, std::min(tot_bytes - i, (size_t) PAGE_SIZE));
}

void mem_row_dense_matrix::set_data(const set_operate &op)
{
	size_t ncol = get_num_cols();
	size_t nrow = get_num_rows();
#pragma omp parallel for
	for (size_t i = 0; i < nrow; i++)
		op.set(get_row(i), ncol, i, 0);
}

bool mem_row_dense_matrix::verify_inner_prod(const dense_matrix &m,
		const bulk_operate &left_op, const bulk_operate &right_op) const
{
	if (!mem_dense_matrix::verify_inner_prod(m, left_op, right_op))
		return false;

	const mem_dense_matrix &mem_m = (const mem_dense_matrix &) m;
	if (mem_m.store_layout() != matrix_layout_t::L_COL) {
		BOOST_LOG_TRIVIAL(error)
			<< "The layout of the right matrix has to be column matrix";
		return false;
	}
	else
		return true;
}

/*
 * In this case, this matrix is wide. The right matrix is tall and its data
 * is stored column-wise.
 */
void mem_row_dense_matrix::serial_inner_prod_wide(const dense_matrix &m,
		const bulk_operate &left_op, const bulk_operate &right_op,
		mem_row_dense_matrix &res) const
{
	size_t ncol = this->get_num_cols();
	size_t nrow = this->get_num_rows();
	const mem_col_dense_matrix &col_m = (const mem_col_dense_matrix &) m;
	char *tmp_res = (char *) malloc(SUB_CHUNK_SIZE * left_op.output_entry_size());
	char *tmp_res2 = (char *) malloc(res.get_num_cols() * res.get_entry_size());
	for (size_t k = 0; k < ncol; k += SUB_CHUNK_SIZE) {
		size_t sub_ncol = std::min(SUB_CHUNK_SIZE, ncol - k);
		sub_row_matrix_info sub_left(0, nrow, k, sub_ncol, *this);
		sub_col_matrix_info sub_right(k, sub_ncol, 0, m.get_num_cols(), col_m);
		for (size_t i = 0; i < sub_left.get_num_rows(); i++) {
			for (size_t j = 0; j < sub_right.get_num_cols(); j++) {
				left_op.runAA(sub_ncol, sub_left.get_row(i),
						sub_right.get_col(j), tmp_res);
				right_op.runA(sub_ncol, tmp_res,
						tmp_res2 + res.get_entry_size() * j);
			}
			// This is fine because we assume the input type of the right operator
			// should be the same as the type of the output matrix.
			right_op.runAA(sub_right.get_num_cols(), tmp_res2, res.get_row(i),
					res.get_row(i));
		}
	}
}

/*
 * In this case, this matrix is tall, and we assume the right matrix is
 * small and its shape is close to a square. We don't need to consider
 * the case that the right matrix is wide because the product would
 * be too large to be stored in any storage media.
 */
void mem_row_dense_matrix::serial_inner_prod_tall(const dense_matrix &m,
		const bulk_operate &left_op, const bulk_operate &right_op,
		mem_row_dense_matrix &res) const
{
	size_t ncol = this->get_num_cols();
	size_t nrow = this->get_num_rows();
	const mem_col_dense_matrix &col_m = (const mem_col_dense_matrix &) m;
	char *tmp_res = (char *) malloc(ncol * res.get_entry_size());
	for (size_t i = 0; i < nrow; i++) {
		for (size_t j = 0; j < m.get_num_cols(); j++) {
			left_op.runAA(ncol, get_row(i), col_m.get_col(j), tmp_res);
			right_op.runA(ncol, tmp_res, res.get(i, j));
		}
	}
	free(tmp_res);
}

dense_matrix::ptr mem_row_dense_matrix::serial_inner_prod(const dense_matrix &m,
		const bulk_operate &left_op, const bulk_operate &right_op) const
{
	if (!verify_inner_prod(m, left_op, right_op))
		return dense_matrix::ptr();

	// TODO we need to determine the layout of the output matrix smartly.
	mem_row_dense_matrix::ptr res = mem_row_dense_matrix::create(
			get_num_rows(), m.get_num_cols(), right_op.get_output_type());
	res->serial_reset_data();

	if (is_wide())
		serial_inner_prod_wide(m, left_op, right_op, *res);
	else
		serial_inner_prod_tall(m, left_op, right_op, *res);

	return std::static_pointer_cast<dense_matrix>(res);
}

void mem_row_dense_matrix::inner_prod_wide(const dense_matrix &m,
		const bulk_operate &left_op, const bulk_operate &right_op,
		mem_row_dense_matrix &res) const
{
	const mem_col_dense_matrix &col_m = (const mem_col_dense_matrix &) m;
	size_t ncol = this->get_num_cols();
	size_t nrow = this->get_num_rows();
	int nthreads = get_num_omp_threads();
	std::vector<mem_row_dense_matrix::ptr> local_ms(nthreads);

#pragma omp parallel
	{
		char *tmp_res = (char *) malloc(
				SUB_CHUNK_SIZE * left_op.output_entry_size());
		char *tmp_res2 = (char *) malloc(
				res.get_num_cols() * res.get_entry_size());
		mem_row_dense_matrix::ptr local_m = mem_row_dense_matrix::create(nrow,
				m.get_num_cols(), right_op.get_output_type());
		local_m->serial_reset_data();
#pragma omp for
		for (size_t k = 0; k < ncol; k += SUB_CHUNK_SIZE) {
			size_t sub_ncol = std::min(SUB_CHUNK_SIZE, ncol - k);
			sub_row_matrix_info sub_left(0, nrow, k, sub_ncol, *this);
			sub_col_matrix_info sub_right(k, sub_ncol, 0, m.get_num_cols(), col_m);
			for (size_t i = 0; i < sub_left.get_num_rows(); i++) {
				for (size_t j = 0; j < sub_right.get_num_cols(); j++) {
					left_op.runAA(sub_ncol, sub_left.get_row(i),
							sub_right.get_col(j), tmp_res);
					right_op.runA(sub_ncol, tmp_res,
							tmp_res2 + res.get_entry_size() * j);
				}
				// This is fine because we assume the input type and output type of
				// the right operator should be the same as the type of the output
				// matrix.
				right_op.runAA(sub_right.get_num_cols(), tmp_res2,
						local_m->get_row(i),
						local_m->get_row(i));
			}
		}
		local_ms[get_omp_thread_num()] = local_m;
	}

	// Aggregate the results from omp threads.
	for (size_t i = 0; i < res.get_num_rows(); i++) {
		for (int j = 0; j < nthreads; j++) {
			right_op.runAA(res.get_num_cols(), local_ms[j]->get_row(i),
					res.get_row(i), res.get_row(i));
		}
	}
}

void mem_row_dense_matrix::inner_prod_tall(const dense_matrix &m,
		const bulk_operate &left_op, const bulk_operate &right_op,
		mem_row_dense_matrix &res) const
{
	const mem_col_dense_matrix &col_m = (const mem_col_dense_matrix &) m;
	size_t ncol = this->get_num_cols();
	size_t nrow = this->get_num_rows();

#pragma omp parallel
	{
		char *tmp_res = (char *) malloc(ncol * res.get_entry_size());
#pragma omp for
		for (size_t i = 0; i < nrow; i++) {
			for (size_t j = 0; j < m.get_num_cols(); j++) {
				left_op.runAA(ncol, get_row(i), col_m.get_col(j), tmp_res);
				right_op.runA(ncol, tmp_res, res.get(i, j));
			}
		}
		free(tmp_res);
	}
}

dense_matrix::ptr mem_row_dense_matrix::inner_prod(const dense_matrix &m,
		const bulk_operate &left_op, const bulk_operate &right_op) const
{
	if (!verify_inner_prod(m, left_op, right_op))
		return dense_matrix::ptr();

	// TODO we need to determine the layout of the output matrix smartly.
	mem_row_dense_matrix::ptr res = mem_row_dense_matrix::create(get_num_rows(),
			m.get_num_cols(), right_op.get_output_type());
	res->reset_data();

	if (is_wide())
		inner_prod_wide(m, left_op, right_op, *res);
	else
		inner_prod_tall(m, left_op, right_op, *res);

	return std::static_pointer_cast<dense_matrix>(res);
}

bool mem_row_dense_matrix::aggregate(const bulk_operate &op,
		scalar_variable &res) const
{
	return get_t_mat().aggregate(op, res);
}

dense_matrix::ptr mem_row_dense_matrix::mapply2(const dense_matrix &m,
		const bulk_operate &op) const
{
	// The same shape and the same data layout.
	if (!verify_mapply2(m, op))
		return dense_matrix::ptr();

	const mem_row_dense_matrix &row_m = (const mem_row_dense_matrix &) m;
	return get_t_mat().mapply2(row_m.get_t_mat(), op)->transpose();
}

dense_matrix::ptr mem_row_dense_matrix::sapply(const bulk_uoperate &op) const
{
	return get_t_mat().sapply(op)->transpose();
}

dense_matrix::ptr mem_row_dense_matrix::apply(apply_margin margin,
		const arr_apply_operate &op) const
{
	if (margin == apply_margin::MAR_ROW)
		margin = apply_margin::MAR_COL;
	else
		margin = apply_margin::MAR_ROW;
	return get_t_mat().apply(margin, op)->transpose();
}

bool mem_row_dense_matrix::write2file(const std::string &file_name) const
{
	FILE *f = fopen(file_name.c_str(), "w");
	if (f == NULL) {
		BOOST_LOG_TRIVIAL(error)
			<< boost::format("can't open %1%: %2%") % file_name % strerror(errno);
		return false;
	}
	if (!write_header(f))
		return false;

	size_t nrow = get_num_rows();
	size_t row_size = get_num_cols() * get_entry_size();
	for (size_t i = 0; i < nrow; i++) {
		const char *row = get_row(i);
		size_t ret = fwrite(row, row_size, 1, f);
		if (ret == 0) {
			BOOST_LOG_TRIVIAL(error)
				<< boost::format("can't write to %1%: %2%")
				% file_name % strerror(errno);
			return false;
		}
	}
	fclose(f);
	return true;
}

mem_row_dense_matrix::ptr mem_row_dense_matrix::create(size_t nrow, size_t ncol,
		const scalar_type &type, FILE *f)
{
	size_t mat_size = nrow * ncol * type.get_size();
	std::shared_ptr<char> data = std::shared_ptr<char>((char *) memalign(
				PAGE_SIZE, mat_size), deleter());
	if (data == NULL) {
		BOOST_LOG_TRIVIAL(error) << "can't allocate memory for the matrix";
		return mem_row_dense_matrix::ptr();
	}
	size_t ret = fread(data.get(), mat_size, 1, f);
	if (ret == 0) {
		BOOST_LOG_TRIVIAL(error)
			<< boost::format("can't read %1% bytes from the file") % mat_size;
		return mem_row_dense_matrix::ptr();
	}

	return mem_row_dense_matrix::ptr(new mem_row_dense_matrix(nrow, ncol,
				type, data));
}

mem_col_dense_matrix &mem_row_dense_matrix::get_t_mat()
{
	if (t_mat == NULL)
		t_mat = mem_col_dense_matrix::cast(transpose());
	return *t_mat;
}

const mem_col_dense_matrix &mem_row_dense_matrix::get_t_mat() const
{
	mem_row_dense_matrix *mutable_this = (mem_row_dense_matrix *) this;
	if (mutable_this->t_mat == NULL)
		mutable_this->t_mat = mem_col_dense_matrix::cast(transpose());
	return *t_mat;
}

}
