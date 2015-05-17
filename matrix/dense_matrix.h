#ifndef __DENSE_MATRIX_H__
#define __DENSE_MATRIX_H__

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

#include <stdlib.h>

#include <vector>
#include <memory>

#include "generic_type.h"
#include "matrix_header.h"
#include "bulk_operate.h"
#include "matrix_store.h"

namespace fm
{

class bulk_operate;
class bulk_uoperate;
class set_operate;
class arr_apply_operate;

enum apply_margin
{
	MAR_ROW = 1,
	MAR_COL = 2,
};

class dense_matrix;

namespace detail
{

class local_matrix_store;

class portion_mapply_op
{
	size_t out_num_rows;
	size_t out_num_cols;
	const scalar_type &type;
public:
	typedef std::shared_ptr<const portion_mapply_op> const_ptr;

	portion_mapply_op(size_t out_num_rows, size_t out_num_cols,
			const scalar_type &_type): type(_type) {
		this->out_num_rows = out_num_rows;
		this->out_num_cols = out_num_cols;
	}

	virtual portion_mapply_op::const_ptr transpose() const = 0;

	virtual void run(
			const std::vector<std::shared_ptr<const local_matrix_store> > &ins,
			local_matrix_store &out) const = 0;

	size_t get_out_num_rows() const {
		return out_num_rows;
	}
	size_t get_out_num_cols() const {
		return out_num_cols;
	}
	const scalar_type &get_output_type() const {
		return type;
	}
};

std::shared_ptr<dense_matrix> mapply_portion(
		const std::vector<std::shared_ptr<const dense_matrix> > &mats,
		// A user can specify the layout of the output dense matrix.
		portion_mapply_op::const_ptr op, matrix_layout_t out_layout);

matrix_store::ptr __mapply_portion(
		const std::vector<matrix_store::const_ptr> &mats,
		portion_mapply_op::const_ptr op, matrix_layout_t out_layout);

matrix_store::ptr __mapply_portion_virtual(
		const std::vector<matrix_store::const_ptr> &store,
		portion_mapply_op::const_ptr op, matrix_layout_t out_layout);

}

/*
 * This class represents a dense matrix and is able to perform computation
 * on the matrix. However, this class can't modify the matrix data. The only
 * way to modify the matrix is to have the pointer point to another matrix.
 */
class dense_matrix
{
	detail::matrix_store::const_ptr store;

	class multiply_scalar_op: public detail::portion_mapply_op {
		scalar_variable::ptr var;
		const bulk_operate &op;
	public:
		multiply_scalar_op(scalar_variable::ptr var, size_t out_num_rows,
				size_t out_num_cols): detail::portion_mapply_op(out_num_rows,
					out_num_cols, var->get_type()),
				op(var->get_type().get_basic_ops().get_multiply()) {
			this->var = var;
		}
		void run(const std::vector<std::shared_ptr<const detail::local_matrix_store> > &ins,
				detail::local_matrix_store &out) const;
		detail::portion_mapply_op::const_ptr transpose() const {
			return detail::portion_mapply_op::const_ptr(new multiply_scalar_op(
						var, get_out_num_cols(), get_out_num_rows()));
		}
	};

protected:
	dense_matrix(detail::matrix_store::const_ptr store) {
		this->store = store;
	}
	virtual bool verify_aggregate(const bulk_operate &op) const;
	virtual bool verify_inner_prod(const dense_matrix &m,
		const bulk_operate &left_op, const bulk_operate &right_op) const;
	virtual bool verify_mapply2(const dense_matrix &m,
			const bulk_operate &op) const;
	virtual bool verify_apply(apply_margin margin, const arr_apply_operate &op) const;
public:
	typedef std::shared_ptr<dense_matrix> ptr;
	typedef std::shared_ptr<const dense_matrix> const_ptr;

	static ptr create(size_t nrow, size_t ncol, const scalar_type &type,
			matrix_layout_t layout, bool in_mem);
	static ptr create(detail::matrix_store::const_ptr store);

	virtual ~dense_matrix() {
	}

	const detail::matrix_store &get_data() const {
		return *store;
	}

	detail::matrix_store::const_ptr get_raw_store() const {
		return store;
	}

	size_t get_entry_size() const {
		return store->get_entry_size();
	}

	size_t get_num_rows() const {
		return store->get_num_rows();
	}

	size_t get_num_cols() const {
		return store->get_num_cols();
	}

	const scalar_type &get_type() const {
		return store->get_type();
	}

	matrix_layout_t store_layout() const {
		return store->store_layout();
	}

	bool is_in_mem() const {
		return store->is_in_mem();
	}

	bool is_wide() const {
		return store->is_wide();
	}

	bool is_virtual() const {
		return store->is_virtual();
	}

	void materialize_self() const;

	template<class T>
	bool is_type() const {
		return get_type() == get_scalar_type<T>();
	}

	/*
	 * We can't change the matrix data that it points to, but we can change
	 * the pointer in the class so that it can point to another matrix data.
	 */
	void assign(const dense_matrix &mat) {
		store = mat.store;
	}

	virtual dense_matrix::ptr get_cols(
			const std::vector<off_t> &idxs) const = 0;
	virtual dense_matrix::ptr get_rows(
			const std::vector<off_t> &idxs) const = 0;
	/*
	 * Clone the matrix.
	 * The class can't modify the matrix data that it points to, but it
	 * can modify the pointer. If someone changes in the pointer in the cloned
	 * matrix, it doesn't affect the current matrix.
	 */
	virtual dense_matrix::ptr clone() const = 0;

	virtual dense_matrix::ptr transpose() const = 0;

	virtual dense_matrix::ptr inner_prod(const dense_matrix &m,
			const bulk_operate &left_op, const bulk_operate &right_op,
			matrix_layout_t out_layout) const = 0;
	virtual std::shared_ptr<scalar_variable> aggregate(
			const bulk_operate &op) const = 0;
	/*
	 * A subclass should define this method for element-wise operations.
	 */
	virtual dense_matrix::ptr mapply2(const dense_matrix &m,
			const bulk_operate &op) const = 0;
	virtual dense_matrix::ptr sapply(const bulk_uoperate &op) const = 0;
	virtual dense_matrix::ptr apply(apply_margin margin,
			arr_apply_operate::const_ptr op) const = 0;
	virtual dense_matrix::ptr scale_cols(
			std::shared_ptr<const mem_vector> vals) const = 0;
	virtual dense_matrix::ptr scale_rows(
			std::shared_ptr<const mem_vector> vals) const = 0;

	dense_matrix::ptr multiply(const dense_matrix &mat,
			matrix_layout_t out_layout) const {
		return inner_prod(mat, get_type().get_basic_ops().get_multiply(),
				get_type().get_basic_ops().get_add(), out_layout);
	}

	dense_matrix::ptr add(const dense_matrix &mat) const {
		return this->mapply2(mat, get_type().get_basic_ops().get_add());
	}

	template<class T>
	dense_matrix::ptr multiply_scalar(T val) const {
		assert(get_type() == get_scalar_type<T>());
		scalar_variable::ptr scal_var(new scalar_variable_impl<T>(val));
		std::vector<detail::matrix_store::const_ptr> stores(1);
		stores[0] = store;
		detail::portion_mapply_op::const_ptr op(new multiply_scalar_op(
					scal_var, get_num_rows(), get_num_cols()));
		detail::matrix_store::ptr ret = __mapply_portion_virtual(stores, op,
				store_layout());
		return dense_matrix::create(ret);
	}

	double norm2() const;
};

}

#endif
