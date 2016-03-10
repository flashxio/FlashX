#ifndef __FM_GROUP_MATRIX_H__
#define __FM_GROUP_MATRIX_H__

/*
 * Copyright 2016 Open Connectome Project (http://openconnecto.me)
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

#include "dense_matrix.h"
#include "combined_matrix_store.h"

namespace fm
{

/*
 * This class specializes the computation on a group of matrices of the same
 * size (except the last one). Such a group of matrices are used to store
 * a tall matrix with many columns or a wide matrix with many rows. We need
 * to explicitly optimize the computation on it to improve performance.
 */
class block_matrix: public dense_matrix
{
	size_t block_size;
	detail::combined_matrix_store::const_ptr store;

	block_matrix(detail::combined_matrix_store::const_ptr store): dense_matrix(
			store) {
		if (store->get_mat_ref(0).is_wide())
			this->block_size = store->get_mat_ref(0).get_num_rows();
		else
			this->block_size = store->get_mat_ref(0).get_num_cols();
		this->store = store;
	}

	static dense_matrix::ptr create(
			detail::combined_matrix_store::const_ptr store);
protected:
	virtual dense_matrix::ptr inner_prod_tall(const dense_matrix &m,
			bulk_operate::const_ptr left_op, bulk_operate::const_ptr right_op,
			matrix_layout_t out_layout) const;
	virtual dense_matrix::ptr inner_prod_wide(const dense_matrix &m,
			bulk_operate::const_ptr left_op, bulk_operate::const_ptr right_op,
			matrix_layout_t out_layout) const;
public:
	typedef std::shared_ptr<block_matrix> ptr;

	static dense_matrix::ptr create(size_t num_rows, size_t num_cols,
			size_t block_size, const scalar_type &type, const set_operate &op,
			int num_nodes = -1, bool in_mem = true,
			safs::safs_file_group::ptr group = NULL);

	size_t get_block_size() const {
		return block_size;
	}

	virtual matrix_layout_t store_layout() const;

	virtual bool is_virtual() const;

	virtual void materialize_self() const;

	virtual void set_materialize_level(materialize_level level,
			detail::matrix_store::ptr materialize_buf);

	virtual void assign(const dense_matrix &mat);

	virtual std::shared_ptr<vector> get_col(off_t idx) const;
	virtual std::shared_ptr<vector> get_row(off_t idx) const;
	virtual dense_matrix::ptr get_cols(const std::vector<off_t> &idxs) const;
	virtual dense_matrix::ptr get_rows(const std::vector<off_t> &idxs) const;

#if 0
	virtual dense_matrix::ptr aggregate(matrix_margin margin,
			agg_operate::const_ptr op) const;
	virtual dense_matrix::ptr groupby_row(std::shared_ptr<const factor_col_vector> labels,
			agg_operate::const_ptr) const;

	virtual dense_matrix::ptr deep_copy() const;
	virtual dense_matrix::ptr conv2(matrix_layout_t layout) const;
	virtual dense_matrix::ptr conv_store(bool in_mem, int num_nodes) const;
	virtual bool move_store(bool in_mem, int num_nodes) const;
	virtual bool drop_cache();
	virtual size_t get_num_cached() const;
#endif

	virtual dense_matrix::ptr clone() const;

	virtual dense_matrix::ptr transpose() const;

	virtual dense_matrix::ptr multiply(const dense_matrix &mat,
			matrix_layout_t out_layout, bool use_blas) const;


	virtual dense_matrix::ptr mapply_cols(std::shared_ptr<const vector> vals,
			bulk_operate::const_ptr op) const;
	virtual dense_matrix::ptr mapply_rows(std::shared_ptr<const vector> vals,
			bulk_operate::const_ptr op) const;
	virtual dense_matrix::ptr mapply2(const dense_matrix &m,
			bulk_operate::const_ptr op) const;
	virtual dense_matrix::ptr sapply(bulk_uoperate::const_ptr op) const;
	virtual dense_matrix::ptr apply(matrix_margin margin,
			arr_apply_operate::const_ptr op) const;
};

}

#endif
