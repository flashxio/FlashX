#ifndef __EM_DENSE_MATRIX_H__
#define __EM_DENSE_MATRIX_H__

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

#include <memory>
#include <boost/format.hpp>

#include "log.h"

#include "bulk_operate.h"
#include "EM_vector.h"
#include "dense_matrix.h"

namespace fm
{

class mem_dense_matrix;
class mem_col_dense_matrix;
class bulk_operate;

struct submatrix_loc
{
	size_t start_row;
	size_t start_col;
	size_t nrow;
	size_t ncol;
};

class submatrix_compute
{
	size_t start_row;
	size_t start_col;
public:
	typedef std::shared_ptr<submatrix_compute> ptr;
	submatrix_compute(size_t start_row, size_t start_col) {
		this->start_row = start_row;
		this->start_col = start_col;
	}

	virtual void run(const mem_dense_matrix &subm) = 0;
};

class EM_dense_matrix_accessor
{
public:
	typedef std::shared_ptr<EM_dense_matrix_accessor> ptr;

	virtual ~EM_dense_matrix_accessor() {
	}

	virtual bool fetch_submatrix(size_t start_row, size_t nrow,
			size_t start_col, size_t ncol,
			submatrix_compute::ptr compute) = 0;
	virtual bool set_submatrix(size_t start_row, size_t start_col,
			std::shared_ptr<mem_dense_matrix> subm) = 0;
	virtual int num_pending_reqs() const = 0;
	virtual void wait4complete(int num) = 0;
	virtual void wait4all() = 0;
};

class EM_dense_matrix: public dense_matrix
{
protected:
	EM_dense_matrix(size_t nrow, size_t ncol,
			const scalar_type &type): dense_matrix(nrow, ncol, type, false) {
	}
public:
	typedef std::shared_ptr<EM_dense_matrix> ptr;

	static ptr cast(dense_matrix::ptr m);
	static bool exist(const std::string &name);

	virtual ~EM_dense_matrix() {
	}

	virtual EM_dense_matrix_accessor::ptr create_accessor() = 0;
};

/*
 * In this matrix class, data is stored in columns.
 * This has to be a very narrow matrix, i.e., the column length must be much
 * larger than the row length.
 */
class EM_col_dense_matrix: public EM_dense_matrix
{
	// The number of elements.
	static const size_t COL_CHUNK_SIZE;

	EM_vector::ptr data;

	EM_col_dense_matrix(size_t nrow, size_t ncol,
			const scalar_type &type): EM_dense_matrix(nrow, ncol, type) {
		data = EM_vector::create(nrow * ncol, type.get_size());
	}

	EM_col_dense_matrix(size_t nrow, size_t ncol, const scalar_type &type,
			const std::string &name): EM_dense_matrix(nrow, ncol, type) {
		data = EM_vector::create(nrow * ncol, type.get_size(), name);
	}

	void split_matrix(std::vector<submatrix_loc> &locs) const;
public:
	static ptr create(size_t nrow, size_t ncol, const scalar_type &type) {
		return ptr(new EM_col_dense_matrix(nrow, ncol, type));
	}

	static ptr create(size_t nrow, size_t ncol, const scalar_type &type,
			const std::string &name) {
		return ptr(new EM_col_dense_matrix(nrow, ncol, type, name));
	}

	virtual bool write2file(const std::string&) const {
		BOOST_LOG_TRIVIAL(error) << "can't write an EM matrix to a file";
		return false;
	}

	virtual matrix_layout_t store_layout() const {
		return matrix_layout_t::L_COL;
	}

	virtual dense_matrix::ptr inner_prod(const dense_matrix &m,
			const bulk_operate &left_op, const bulk_operate &right_op) const;
	virtual bool aggregate(const bulk_operate &op, scalar_variable &res) const;
	virtual dense_matrix::ptr mapply2(const dense_matrix &m,
			const bulk_operate &op) const;
	virtual dense_matrix::ptr sapply(const bulk_uoperate &op) const;
	virtual dense_matrix::ptr apply(apply_margin margin,
			const arr_apply_operate &op) const;

	virtual void set_data(const set_operate &op);
	virtual void reset_data();
	virtual dense_matrix::ptr shallow_copy() const;
	virtual dense_matrix::ptr deep_copy() const;
	virtual dense_matrix::ptr conv2(size_t nrow, size_t ncol, bool byrow) const;
	virtual dense_matrix::ptr transpose() const;

	EM_dense_matrix_accessor::ptr create_accessor();
};

class EM_subvec_accessor;

class EM_col_matrix_accessor: public EM_dense_matrix_accessor
{
	EM_col_dense_matrix &m;
	EM_vector_accessor::ptr accessor;

	// The buffers for requests issued by sub accessors.
	std::vector<fetch_vec_request> fetch_reqs;
	std::vector<set_vec_request> set_reqs;
	std::vector<std::shared_ptr<EM_subvec_accessor> > sub_accessors;

	int pending_reqs;

	void flush();
public:
	typedef std::shared_ptr<EM_col_matrix_accessor> ptr;

	EM_col_matrix_accessor(EM_col_dense_matrix &_m, EM_vector::ptr cols);
	~EM_col_matrix_accessor();

	bool fetch_submatrix(size_t start_row, size_t nrow,
			size_t start_col, size_t ncol, submatrix_compute::ptr compute);
	bool set_submatrix(size_t start_row, size_t start_col,
			std::shared_ptr<mem_dense_matrix> subm);
	/*
	 * The number of pending I/O requests.
	 */
	int num_pending_reqs() const {
		return pending_reqs;
	}
	/*
	 * Notify the completion of an I/O request.
	 */
	void complete_req() {
		pending_reqs--;
	}
	/*
	 * Wait for the specified number of I/O requests to be completed.
	 */
	virtual void wait4complete(int num);
	virtual void wait4all();
};

template<class LeftType, class RightType, class ResType>
EM_dense_matrix::ptr multiply(EM_dense_matrix &m1, dense_matrix &m2)
{
	basic_ops_impl<LeftType, RightType, ResType> ops;
	return EM_dense_matrix::cast(m1.inner_prod(m2, ops.get_multiply(), ops.get_add()));
}

}

#endif
