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

#include "matrix_config.h"
#include "block_matrix.h"
#include "vector.h"
#include "local_matrix_store.h"
#include "one_val_matrix_store.h"
#include "col_vec.h"
#include "sink_matrix.h"
#include "agg_matrix_store.h"
#include "project_matrix_store.h"
#include "IPW_matrix_store.h"
#include "materialize.h"
#include "mem_matrix_store.h"

namespace fm
{

dense_matrix::ptr block_matrix::create(
		detail::combined_matrix_store::const_ptr store)
{
	if (store->get_mat_ref(0).is_wide()) {
		for (size_t i = 1; i < store->get_num_mats() - 1; i++)
			if (store->get_mat_ref(i).get_num_rows()
					!= store->get_mat_ref(i - 1).get_num_rows()) {
				BOOST_LOG_TRIVIAL(error)
					<< "The matrices have different block sizes";
				return dense_matrix::ptr();
			}
	}
	else {
		for (size_t i = 1; i < store->get_num_mats() - 1; i++)
			if (store->get_mat_ref(i).get_num_cols()
					!= store->get_mat_ref(i - 1).get_num_cols()) {
				BOOST_LOG_TRIVIAL(error)
					<< "The matrices have different block sizes";
				return dense_matrix::ptr();
			}
	}
	return dense_matrix::ptr(new block_matrix(store));
}

dense_matrix::ptr block_matrix::create_layout(scalar_variable::ptr val,
		size_t num_rows, size_t num_cols, matrix_layout_t layout,
		size_t block_size, int num_nodes, bool in_mem,
		safs::safs_file_group::ptr group)
{
	// For tall matrices
	if (num_rows > num_cols) {
		// If there is only one block.
		if (num_cols < block_size)
			return dense_matrix::create_const(val, num_rows, num_cols,
					layout, num_nodes, in_mem, group);

		std::vector<detail::matrix_store::const_ptr> stores(div_ceil<size_t>(
					num_cols, block_size));
		for (size_t i = 0; i < stores.size(); i++) {
			size_t local_num_cols = std::min(num_cols - i * block_size,
					block_size);
			stores[i] = detail::matrix_store::ptr(new detail::one_val_matrix_store(
						val, num_rows, local_num_cols, layout, num_nodes));
		}
		return block_matrix::create(detail::combined_matrix_store::create(
					stores, layout));
	}
	else {
		// If there is only one block.
		if (num_rows < block_size)
			return dense_matrix::create_const(val, num_rows, num_cols,
					layout, num_nodes, in_mem, group);

		std::vector<detail::matrix_store::const_ptr> stores(div_ceil<size_t>(
					num_rows, block_size));
		for (size_t i = 0; i < stores.size(); i++) {
			size_t local_num_rows = std::min(num_rows - i * block_size,
					block_size);
			stores[i] = detail::matrix_store::ptr(new detail::one_val_matrix_store(
						val, local_num_rows, num_cols, layout, num_nodes));
		}
		return block_matrix::create(detail::combined_matrix_store::create(
					stores, layout));
	}
}

dense_matrix::ptr block_matrix::create_layout(size_t num_rows, size_t num_cols,
		matrix_layout_t layout, size_t block_size, const scalar_type &type,
		const set_operate &op, int num_nodes, bool in_mem,
		safs::safs_file_group::ptr group)
{
	detail::combined_matrix_store::ptr combined;
	// For tall matrices
	if (num_rows > num_cols) {
		// If there is only one block.
		if (num_cols < block_size)
			return dense_matrix::create(num_rows, num_cols, layout, type, op,
					num_nodes, in_mem, group);

		std::vector<detail::matrix_store::ptr> stores(div_ceil<size_t>(
					num_cols, block_size));
		for (size_t i = 0; i < stores.size(); i++) {
			size_t local_num_cols = std::min(num_cols - i * block_size,
					block_size);
			stores[i] = detail::matrix_store::create(num_rows, local_num_cols,
					layout, type, num_nodes, in_mem, group);
			// TODO we lose some information if we initialize them individually.
			stores[i]->set_data(op);
		}
		std::vector<detail::matrix_store::const_ptr> const_stores(
				stores.begin(), stores.end());
		combined = detail::combined_matrix_store::create(const_stores, layout);
	}
	else {
		// If there is only one block.
		if (num_rows < block_size)
			return dense_matrix::create(num_rows, num_cols, layout, type,
					op, num_nodes, in_mem, group);

		std::vector<detail::matrix_store::ptr> stores(div_ceil<size_t>(
					num_rows, block_size));
		for (size_t i = 0; i < stores.size(); i++) {
			size_t local_num_rows = std::min(num_rows - i * block_size,
					block_size);
			stores[i] = detail::matrix_store::create(local_num_rows, num_cols,
					layout, type, num_nodes, in_mem, group);
			// TODO we lose some information if we initialize them individually.
			stores[i]->set_data(op);
		}
		std::vector<detail::matrix_store::const_ptr> const_stores(
				stores.begin(), stores.end());
		combined = detail::combined_matrix_store::create(const_stores, layout);
	}
	if (combined)
		return block_matrix::create(combined);
	else
		return dense_matrix::ptr();
}

dense_matrix::ptr block_matrix::create_seq_layout(scalar_variable::ptr start,
		scalar_variable::ptr stride, scalar_variable::ptr seq_ele_stride,
		size_t num_rows, size_t num_cols, matrix_layout_t layout,
		size_t block_size, bool byrow, int num_nodes, bool in_mem,
		safs::safs_file_group::ptr group)
{
	if (num_rows > num_cols && num_cols < block_size)
		return dense_matrix::create_seq(start, stride, seq_ele_stride,
				num_rows, num_cols, layout, byrow, num_nodes, in_mem, group);
	else if (num_rows <= num_cols && num_rows < block_size)
		return dense_matrix::create_seq(start, stride, seq_ele_stride,
				num_rows, num_cols, layout, byrow, num_nodes, in_mem, group);

	const scalar_type &type = start->get_type();
	assert(type == stride->get_type());
	assert(type == seq_ele_stride->get_type());

	size_t num_blocks;
	if (num_rows > num_cols)
		num_blocks = div_ceil<size_t>(num_cols, block_size);
	else
		num_blocks = div_ceil<size_t>(num_rows, block_size);

	const bulk_operate &mul = type.get_basic_ops().get_multiply();
	const bulk_operate &add = type.get_basic_ops().get_add();
	std::vector<detail::matrix_store::const_ptr> stores(num_blocks);
	for (size_t i = 0; i < stores.size(); i++) {
		// Calculate the value of the first element of a block.
		size_t len;
		if ((num_rows > num_cols && byrow) || (num_rows <= num_cols && !byrow))
			len = block_size * i;
		else if (num_rows > num_cols)
			len = block_size * num_rows * i;
		else
			len = block_size * num_cols * i;
		scalar_variable::ptr len_var(new scalar_variable_impl<size_t>(len));
		len_var = len_var->cast_type(type);
		scalar_variable::ptr block_stride = type.create_scalar();
		mul.runAA(1, len_var->get_raw(), stride->get_raw(), block_stride->get_raw());
		scalar_variable::ptr lstart = type.create_scalar();
		add.runAA(1, start->get_raw(), block_stride->get_raw(),
				lstart->get_raw());

		size_t local_num_cols, local_num_rows;
		if (num_rows > num_cols) {
			local_num_rows = num_rows;
			local_num_cols = std::min(num_cols - i * block_size, block_size);
		}
		else {
			local_num_rows = std::min(num_rows - i * block_size, block_size);
			local_num_cols = num_cols;
		}


		detail::matrix_store::ptr store = detail::matrix_store::create(
				local_num_rows, local_num_cols, layout, type, num_nodes,
				in_mem, group);
		auto op = type.get_set_seq(*lstart, *stride, *seq_ele_stride,
				num_rows, num_cols, byrow);
		store->set_data(*op);
		stores[i] = store;
	}
	return block_matrix::create(detail::combined_matrix_store::create(
				stores, layout));

}

matrix_layout_t block_matrix::store_layout() const
{
	// All matrices in the group should have the same data layout.
	return store->get_mat_ref(0).store_layout();
}

bool block_matrix::is_virtual() const
{
	// If one matrix is virtual, all other matrices should also be virtual.
	return store->get_mat_ref(0).is_virtual();
}

bool block_matrix::materialize_self() const
{
	if (!is_virtual())
		return true;

	std::vector<detail::matrix_store::const_ptr> res_stores(store->get_num_mats());
	bool ret = true;
	// TODO materializing individual matrices in serial may hurt performance.
	for (size_t i = 0; i < store->get_num_mats(); i++) {
		dense_matrix::ptr mat = dense_matrix::create(store->get_mat(i));
		ret = ret && mat->materialize_self();
		res_stores[i] = mat->get_raw_store();
	}
	if (!ret)
		return false;

	block_matrix *mutable_this = const_cast<block_matrix *>(this);
	mutable_this->store = detail::combined_matrix_store::create(res_stores,
			res_stores[0]->store_layout());

	// This is only way to change the store pointer in dense_matrix.
	dense_matrix::ptr tmp = dense_matrix::create(this->store);
	mutable_this->dense_matrix::assign(*tmp);
	return true;
}

void block_matrix::set_materialize_level(materialize_level level,
		detail::matrix_store::ptr materialize_buf)
{
	for (size_t i = 0; i < store->get_num_mats(); i++) {
		detail::matrix_store::const_ptr mat = store->get_mat(i);
		if (mat->is_virtual()) {
			const detail::virtual_matrix_store *tmp
				= dynamic_cast<const detail::virtual_matrix_store *>(mat.get());
			detail::virtual_matrix_store *tmp1
				= const_cast<detail::virtual_matrix_store *>(tmp);
			// TODO we ignore the customization of materialized matrix store.
			tmp1->set_materialize_level(level, NULL);
		}
	}
}

std::vector<detail::virtual_matrix_store::const_ptr> block_matrix::get_compute_matrices() const
{
	if (!is_virtual())
		return std::vector<detail::virtual_matrix_store::const_ptr>();

	// Double check that a block matrix shouldn't be a sink matrix.
	detail::sink_store::const_ptr sink
		= std::dynamic_pointer_cast<const detail::sink_store>(get_raw_store());
	assert(sink == NULL);

	std::vector<detail::virtual_matrix_store::const_ptr> vmats;
	for (size_t i = 0; i < store->get_num_mats(); i++) {
		detail::virtual_matrix_store::const_ptr vmat
			= std::dynamic_pointer_cast<const detail::virtual_matrix_store>(
					store->get_mat(i));
		if (vmat)
			vmats.push_back(vmat);
	}
	return vmats;
}

void block_matrix::assign(const dense_matrix &mat)
{
	// The input matrix must be a block matrix. Otherwise, dynamic
	// casting will throw an exception.
	const block_matrix &gmat = dynamic_cast<const block_matrix &>(mat);
	this->store = gmat.store;
	dense_matrix::assign(mat);
}

dense_matrix::ptr block_matrix::get_cols(const std::vector<off_t> &idxs) const
{
	auto cols = store->get_cols(idxs);
	if (cols == NULL)
		return dense_matrix::ptr();
	else
		return dense_matrix::create(cols);
}

dense_matrix::ptr block_matrix::get_rows(const std::vector<off_t> &idxs) const
{
	auto rows = store->get_rows(idxs);
	if (rows == NULL)
		return dense_matrix::ptr();
	else
		return dense_matrix::create(rows);
}

dense_matrix::ptr block_matrix::clone() const
{
	return block_matrix::create(store);
}

dense_matrix::ptr block_matrix::transpose() const
{
	detail::matrix_store::const_ptr tmp = store->transpose();
	return block_matrix::create(detail::combined_matrix_store::cast(tmp));
}

static detail::mem_matrix_store::const_ptr get_sub_mat(
		detail::mem_matrix_store::const_ptr mat, size_t start_row,
		size_t start_col, size_t num_rows, size_t num_cols)
{
	detail::local_matrix_store::const_ptr portion = mat->get_portion(
			start_row, start_col, num_rows, num_cols);
	assert(portion);

	detail::mem_matrix_store::ptr ret = detail::mem_matrix_store::create(
			num_rows, num_cols, mat->store_layout(), portion->get_type(), -1);
	ret->write_portion_async(portion, 0, 0);
	return ret;
}

namespace
{

/*
 * This is a summation with generalized operator.
 */
class gsum_op: public fm::detail::portion_mapply_op
{
	bulk_operate::const_ptr op;
public:
	gsum_op(bulk_operate::const_ptr op, size_t out_num_rows,
			size_t out_num_cols): fm::detail::portion_mapply_op(
				out_num_rows, out_num_cols, op->get_output_type()) {
		this->op = op;
	}

	virtual void run(
			const std::vector<fm::detail::local_matrix_store::const_ptr> &ins,
			fm::detail::local_matrix_store &out) const {
		assert(!ins.empty());
		if (ins.size() == 1)
			out.copy_from(*ins[0]);
		else {
			detail::part_dim_t dim = get_out_num_rows() > get_out_num_cols()
				? detail::part_dim_t::PART_DIM1 : detail::part_dim_t::PART_DIM2;
			mapply2(*ins[0], *ins[1], *op, dim, out);
			for (size_t i = 2; i < ins.size(); i++)
				mapply2(*ins[i], out, *op, dim, out);
		}
	}

	virtual fm::detail::portion_mapply_op::const_ptr transpose() const;

	virtual std::string to_string(
			const std::vector<detail::matrix_store::const_ptr> &mats) const {
		std::string str;
		assert(mats.size() > 1);
		str = "(";
		for (size_t i = 0; i < mats.size() - 1; i++)
			str += mats[i]->get_name() + op->get_name();
		str += mats[mats.size() - 1]->get_name() + ")";
		return str;
	}

	virtual bool is_agg() const {
		return true;
	}

	bulk_operate::const_ptr get_op() const {
		return op;
	}
};

class t_gsum_op: public fm::detail::portion_mapply_op
{
	bulk_operate::const_ptr op;
	gsum_op portion_op;
public:
	t_gsum_op(const gsum_op &_op): fm::detail::portion_mapply_op(
			_op.get_out_num_cols(), _op.get_out_num_rows(),
			_op.get_output_type()), portion_op(_op) {
		this->op = _op.get_op();
	}

	virtual void run(
			const std::vector<fm::detail::local_matrix_store::const_ptr> &ins,
			fm::detail::local_matrix_store &out) const {
		assert(!ins.empty());
		if (ins.size() == 1)
			out.copy_from(*ins[0]);
		else {
			detail::part_dim_t dim = get_out_num_rows() > get_out_num_cols()
				? detail::part_dim_t::PART_DIM1 : detail::part_dim_t::PART_DIM2;
			mapply2(*ins[0], *ins[1], *op, dim, out);
			for (size_t i = 2; i < ins.size(); i++)
				mapply2(*ins[i], out, *op, dim, out);
		}
	}

	virtual fm::detail::portion_mapply_op::const_ptr transpose() const {
		return fm::detail::portion_mapply_op::const_ptr(new gsum_op(portion_op));
	}
	virtual std::string to_string(
			const std::vector<detail::matrix_store::const_ptr> &mats) const {
		return portion_op.to_string(mats);
	}

	virtual bool is_agg() const {
		return true;
	}
};

fm::detail::portion_mapply_op::const_ptr gsum_op::transpose() const
{
	return fm::detail::portion_mapply_op::const_ptr(new t_gsum_op(*this));
}

}

dense_matrix::ptr block_matrix::inner_prod_tall(const dense_matrix &m,
			bulk_operate::const_ptr left_op, bulk_operate::const_ptr right_op,
			matrix_layout_t out_layout) const
{
	// Get the right matrix in memory.
	// We don't want the right matrix to be a block matrix either.
	dense_matrix::ptr m2 = dense_matrix::create(m.get_raw_store());
	m2 = m2->conv_store(true, -1);
	detail::mem_matrix_store::const_ptr mem_m2
		= std::dynamic_pointer_cast<const detail::mem_matrix_store>(
				m2->get_raw_store());
	assert(mem_m2);

	// Here is to reuse the code for matrix multiplication with BLAS.
	bool use_blas = left_op == NULL;
	if (use_blas) {
		assert(get_type() == get_scalar_type<double>()
				|| get_type() == get_scalar_type<float>());
		assert(m.get_type() == get_scalar_type<double>()
				|| m.get_type() == get_scalar_type<float>());
		right_op = bulk_operate::conv2ptr(get_type().get_basic_ops().get_add());
	}

	// This contains the blocks for the final output.
	std::vector<detail::matrix_store::const_ptr> res_blocks(
			div_ceil<size_t>(mem_m2->get_num_cols(), get_block_size()));
	for (size_t m2_col = 0; m2_col < mem_m2->get_num_cols();
			m2_col += get_block_size()) {
		// We multiply with individual matrices and output a vector of
		// temporary matrices. Later on, we need to sum the temporary matrices.
		std::vector<dense_matrix::const_ptr> tmp_mats(store->get_num_mats());
		for (size_t m2_row = 0; m2_row < mem_m2->get_num_rows();
				m2_row += get_block_size()) {
			size_t i = m2_row / get_block_size();
			dense_matrix::ptr left = dense_matrix::create(store->get_mat(i));
			// Get the submatrix in the right matrix
			size_t part_num_rows = std::min(get_block_size(),
					m2->get_num_rows() - m2_row);
			size_t part_num_cols = std::min(get_block_size(),
					m2->get_num_cols() - m2_col);
			detail::mem_matrix_store::const_ptr part = get_sub_mat(mem_m2,
					m2_row, m2_col, part_num_rows, part_num_cols);
			dense_matrix::ptr right = dense_matrix::create(part);

			// Compute the temporary matrix.
			// TODO maybe we can perform multiply with multiple block matrices
			// together in the cost of more memory consumption.
			if (use_blas)
				tmp_mats[i] = left->multiply(*right, out_layout);
			else
				tmp_mats[i] = left->inner_prod(*right, left_op, right_op, out_layout);
			// We really don't need to cache the portion in this intermediate
			// matrix and the EM matrix beneath it in the hierarchy.
			const_cast<detail::matrix_store &>(tmp_mats[i]->get_data()).set_cache_portion(false);
		}

		// We then sum all of the temp matrices.
		detail::portion_mapply_op::const_ptr op(new gsum_op(right_op,
					tmp_mats[0]->get_num_rows(), tmp_mats[0]->get_num_cols()));
		// We will materialize the mapply in a hierarchical way,
		// so the intermediate matrices should be read from SSDs in serial.
		size_t i = m2_col / get_block_size();
		dense_matrix::ptr res = mapply_portion(tmp_mats, op,
				tmp_mats[0]->store_layout(), false);
		res->materialize_self();
		res_blocks[i] = res->get_raw_store();
	}

	// TODO we need to restore the original caching policy in the EM matrix.

	if (res_blocks.size() == 1)
		return dense_matrix::create(res_blocks[0]);
	else
		return block_matrix::create(detail::combined_matrix_store::create(
					res_blocks, res_blocks[0]->store_layout()));
}

dense_matrix::ptr block_matrix::inner_prod_wide(const dense_matrix &m,
			bulk_operate::const_ptr left_op, bulk_operate::const_ptr right_op,
			matrix_layout_t out_layout) const
{
	std::vector<detail::matrix_store::const_ptr> right_mats;
	const block_matrix *block_m = dynamic_cast<const block_matrix *>(&m);
	if (block_m == NULL)
		right_mats.push_back(m.get_raw_store());
	else {
		for (size_t i = 0; i < block_m->store->get_num_mats(); i++)
			right_mats.push_back(block_m->store->get_mat(i));
	}

	detail::matrix_store::ptr res;
	// If the left operator isn't defined, we assume it uses BLAS for matrix
	// multiplication.
	bool use_blas = left_op == NULL;
	if (use_blas) {
		assert(get_type() == m.get_type());
		assert(get_type() == get_scalar_type<double>()
				|| get_type() == get_scalar_type<float>());
		res = detail::matrix_store::create(get_num_rows(), m.get_num_cols(),
				out_layout, get_type(), -1, true);
	}
	else
		res = detail::matrix_store::create(get_num_rows(), m.get_num_cols(),
				out_layout, right_op->get_output_type(), -1, true);
	// Each time we take one matrix in the right group and perform inner product
	// with all matrices in the left group.
	size_t right_block_size = right_mats[0]->get_num_cols();
	for (size_t i = 0; i < right_mats.size(); i++) {
		dense_matrix::ptr right = dense_matrix::create(right_mats[i]);
		std::vector<dense_matrix::ptr> tmp_mats(store->get_num_mats());
		for (size_t j = 0; j < store->get_num_mats(); j++) {
			dense_matrix::ptr left = dense_matrix::create(store->get_mat(j));
			if (use_blas)
				tmp_mats[j] = left->multiply(*right, out_layout);
			else
				tmp_mats[j] = left->inner_prod(*right, left_op, right_op,
						out_layout);
			const_cast<detail::matrix_store &>(left->get_data()).set_cache_portion(false);
		}
		materialize(tmp_mats, false);

		// We now copy the inner product result to the final matrix.
		size_t col_idx = i * right_block_size;
		for (size_t j = 0; j < tmp_mats.size(); j++) {
			size_t row_idx = j * block_size;
			size_t num_rows = std::min(block_size, get_num_rows() - row_idx);
			size_t num_cols = std::min(right_block_size,
					m.get_num_cols() - col_idx);
			assert(num_rows == tmp_mats[j]->get_num_rows());
			assert(num_cols == tmp_mats[j]->get_num_cols());
			detail::local_matrix_store::ptr res_part = res->get_portion(row_idx,
					col_idx, num_rows, num_cols);
			detail::local_matrix_store::const_ptr src_part
				= tmp_mats[j]->get_data().get_portion(0);
			res_part->copy_from(*src_part);
		}
	}
	return dense_matrix::create(res);
}

dense_matrix::ptr block_matrix::multiply_tall(const dense_matrix &m,
		matrix_layout_t out_layout) const
{
	return block_matrix::inner_prod_tall(m, NULL, NULL, out_layout);
}

static void get_wider_matrices(detail::combined_matrix_store::const_ptr in,
		std::vector<detail::matrix_store::const_ptr> &mats, size_t block_size)
{
	size_t short_dim = std::min(in->get_mat_ref(0).get_num_rows(),
			in->get_mat_ref(0).get_num_cols());
	block_size = std::min(block_size, matrix_conf.get_max_multiply_block_size());
	// We prefer to round it up.
	size_t num_block_mats = div_ceil(block_size, short_dim);
	if (num_block_mats == 0)
		num_block_mats = 1;
	if (num_block_mats <= 1) {
		for (size_t i = 0; i < in->get_num_mats(); i++)
			mats.push_back(in->get_mat(i));
	}
	else {
		for (size_t i = 0; i < in->get_num_mats(); i += num_block_mats) {
			std::vector<detail::matrix_store::const_ptr> tmp(std::min(
						in->get_num_mats() - i, num_block_mats));
			for (size_t j = 0; j < tmp.size(); j++)
				tmp[j] = in->get_mat(i + j);
			mats.push_back(detail::combined_matrix_store::create(tmp,
						in->store_layout()));
		}
	}
}

dense_matrix::ptr block_matrix::multiply_wide(const dense_matrix &m,
		matrix_layout_t out_layout) const
{
	std::vector<detail::matrix_store::const_ptr> right_mats;
	const block_matrix *block_m = dynamic_cast<const block_matrix *>(&m);
	if (block_m == NULL)
		right_mats.push_back(m.get_raw_store());
	else
		get_wider_matrices(block_m->store, right_mats,
				std::min(get_num_rows(), m.get_num_cols()));

	std::vector<detail::matrix_store::const_ptr> left_mats;
	get_wider_matrices(store, left_mats,
			std::min(get_num_rows(), m.get_num_cols()));

	assert(get_type() == m.get_type());
	assert(get_type() == get_scalar_type<double>()
			|| get_type() == get_scalar_type<float>());
	std::vector<detail::matrix_store::const_ptr> blocks(
			left_mats.size() * right_mats.size());
	// Each time we take one matrix in the right group and perform inner product
	// with all matrices in the left group.
	for (size_t i = 0; i < right_mats.size(); i++) {
		dense_matrix::ptr right = dense_matrix::create(right_mats[i]);
		for (size_t j = 0; j < left_mats.size(); j++) {
			dense_matrix::ptr left = dense_matrix::create(left_mats[j]);
			dense_matrix::ptr res = left->multiply(*right, out_layout);
			assert(res);
			blocks[j * right_mats.size() + i] = res->get_raw_store();
			// This is only necessary for EM matrices.
			if (!left->is_in_mem()) {
				detail::matrix_store &tmp
					= const_cast<detail::matrix_store &>(left->get_data());
				tmp.set_cache_portion(false);
			}
		}
		if (!right->is_in_mem()) {
			detail::matrix_store &tmp
				= const_cast<detail::matrix_store &>(right->get_data());
			tmp.set_cache_portion(true);
		}
	}
	return dense_matrix::create(detail::block_sink_store::create(blocks,
				left_mats.size(), right_mats.size()));
}

/*
 * Multiply a wide block matrix with a sparse matrix.
 * Here we assume the sparse matrix is small and can be stored in memory.
 */
dense_matrix::ptr block_matrix::multiply_sparse_wide(const dense_matrix &m,
		matrix_layout_t out_layout) const
{
	detail::sparse_project_matrix_store::const_ptr right
		= std::dynamic_pointer_cast<const detail::sparse_project_matrix_store>(
				m.get_raw_store());
	assert(right);

	if (out_layout == matrix_layout_t::L_NONE)
		out_layout = matrix_layout_t::L_COL;

	assert(get_type() == m.get_type());
	assert(get_type() == get_scalar_type<double>()
			|| get_type() == get_scalar_type<float>());

	std::vector<dense_matrix::ptr> res_mats;
	std::vector<dense_matrix::ptr> tmps;
	for (size_t i = 0; i < this->store->get_num_mats(); i++) {
		// We have to make sure the matrix is column major.
		detail::matrix_store::const_ptr left = this->store->get_mat(i);
		if (left->store_layout() == matrix_layout_t::L_ROW) {
			dense_matrix::ptr tmp = dense_matrix::create(left);
			tmp = tmp->conv2(matrix_layout_t::L_COL);
			left = tmp->get_raw_store();
		}
		tmps.push_back(dense_matrix::create(detail::matrix_store::ptr(
						new detail::IPW_matrix_store(left, right, NULL, NULL,
							out_layout))));
		// TODO It might be better if we can perform computation on more
		// matrices together. However, computation on more matrices requires
		// more memory allocation.
		if (tmps.size() >= 8) {
			materialize(tmps, false);
			res_mats.insert(res_mats.end(), tmps.begin(), tmps.end());
			tmps.clear();
		}
	}
	if (tmps.size() > 0) {
		materialize(tmps, false);
		res_mats.insert(res_mats.end(), tmps.begin(), tmps.end());
		tmps.clear();
	}

	return dense_matrix::rbind(res_mats);
}

dense_matrix::ptr block_matrix::multiply(const dense_matrix &mat,
		matrix_layout_t out_layout) const
{
	if (mat.get_data().is_sparse()) {
		assert(is_wide());

		// TODO we need to deal with tall matrix.
		detail::combined_matrix_store::const_ptr combined
			= std::dynamic_pointer_cast<const detail::combined_matrix_store>(
					mat.get_raw_store());
		if (combined == NULL)
			return multiply_sparse_wide(mat, out_layout);

		std::vector<dense_matrix::ptr> res_mats(combined->get_num_mats());
		for (size_t i = 0; i < combined->get_num_mats(); i++) {
			dense_matrix::ptr right = dense_matrix::create(combined->get_mat(i));
			res_mats[i] = multiply(*right, out_layout);
		}
		materialize(res_mats, false);
		return dense_matrix::cbind(res_mats);
	}
	else if ((get_type() == get_scalar_type<double>()
				|| get_type() == get_scalar_type<float>())) {
		assert(get_type() == mat.get_type());
		size_t long_dim1 = std::max(get_num_rows(), get_num_cols());
		size_t long_dim2 = std::max(mat.get_num_rows(), mat.get_num_cols());
		// We prefer to perform computation on the larger matrix.
		// If the matrix in the right operand is larger, we transpose
		// the entire computation.
		if (long_dim2 > long_dim1) {
			dense_matrix::ptr t_mat1 = this->transpose();
			dense_matrix::ptr t_mat2 = mat.transpose();
			matrix_layout_t t_layout = out_layout;
			if (t_layout == matrix_layout_t::L_ROW)
				t_layout = matrix_layout_t::L_COL;
			else if (t_layout == matrix_layout_t::L_COL)
				t_layout = matrix_layout_t::L_ROW;
			dense_matrix::ptr t_res = t_mat2->multiply(*t_mat1, t_layout);
			return t_res->transpose();
		}

		if (is_wide())
			return multiply_wide(mat, out_layout);
		else
			return multiply_tall(mat, out_layout);
	}
	else {
		// This relies on inner product to compute matrix multiplication.
		bulk_operate::const_ptr multiply = bulk_operate::conv2ptr(
				get_type().get_basic_ops().get_multiply());
		bulk_operate::const_ptr add = bulk_operate::conv2ptr(
				get_type().get_basic_ops().get_add());
		return inner_prod(mat, multiply, add, out_layout);
	}
}

dense_matrix::ptr block_matrix::mapply_cols(col_vec::const_ptr vals,
		bulk_operate::const_ptr op) const
{
	dense_matrix::ptr tmat = transpose();
	return tmat->mapply_rows(vals, op)->transpose();
}

dense_matrix::ptr block_matrix::mapply_rows(col_vec::const_ptr vals,
		bulk_operate::const_ptr op) const
{
	if (!vals->is_in_mem()) {
		BOOST_LOG_TRIVIAL(error) << "Can't scale rows with an EM vector";
		return dense_matrix::ptr();
	}
	if (get_num_cols() != vals->get_length()) {
		BOOST_LOG_TRIVIAL(error)
			<< "The vector's length needs to equal to #columns";
		return dense_matrix::ptr();
	}
	
	std::vector<detail::matrix_store::const_ptr> res_stores(
			store->get_num_mats());
	if (is_wide()) {
		for (size_t i = 0; i < res_stores.size(); i++) {
			dense_matrix::ptr mat = dense_matrix::create(store->get_mat(i));
			dense_matrix::ptr res = mat->mapply_rows(vals, op);
			if (res == NULL)
				return dense_matrix::ptr();
			res_stores[i] = res->get_raw_store();
		}
		return block_matrix::create(detail::combined_matrix_store::create(
					res_stores, store->store_layout()));
	}
	else {
		size_t val_start = 0;
		vals->move_store(true, -1);
		detail::mem_matrix_store::const_ptr mem_vals
			= detail::mem_matrix_store::cast(vals->get_raw_store());
		for (size_t i = 0; i < res_stores.size(); i++) {
			// Get part of the vector.
			size_t llen = store->get_mat_ref(i).get_num_cols();
			detail::mem_col_matrix_store::ptr vals_store
				= detail::mem_col_matrix_store::create(llen, 1, vals->get_type());
			// mem_vals is stored in a contiguous memory, so we can get
			// the starting address of the memory we want by getting the address
			// of the element.
			memcpy(vals_store->get_col(0), mem_vals->get(val_start, 0),
					llen * vals->get_entry_size());
			col_vec::ptr vals_part = col_vec::create(vals_store);

			// Perform computation.
			dense_matrix::ptr mat = dense_matrix::create(store->get_mat(i));
			dense_matrix::ptr res = mat->mapply_rows(vals_part, op);
			if (res == NULL)
				return dense_matrix::ptr();
			res_stores[i] = res->get_raw_store();

			val_start += llen;
		}
		assert(val_start == vals->get_length());
		return block_matrix::create(detail::combined_matrix_store::create(
					res_stores, store->store_layout()));
	}
}

dense_matrix::ptr block_matrix::mapply2(const dense_matrix &m,
		bulk_operate::const_ptr op) const
{
	if (get_num_rows() != m.get_num_rows()
			|| get_num_cols() != m.get_num_cols()) {
		BOOST_LOG_TRIVIAL(error) << "The matrix size isn't compatible";
		return dense_matrix::ptr();
	}
	const block_matrix *block_m = dynamic_cast<const block_matrix *>(&m);
	if (block_m == NULL) {
		BOOST_LOG_TRIVIAL(error) << "The input matrix isn't a block matrix";
		return dense_matrix::ptr();
	}
	if (block_m->get_block_size() != get_block_size()) {
		BOOST_LOG_TRIVIAL(error)
			<< "The input matrix has a different block size";
		return dense_matrix::ptr();
	}

	std::vector<detail::matrix_store::const_ptr> res_stores(
			store->get_num_mats());
	for (size_t i = 0; i < res_stores.size(); i++) {
		dense_matrix::ptr mat1 = dense_matrix::create(store->get_mat(i));
		dense_matrix::ptr mat2 = dense_matrix::create(block_m->store->get_mat(i));
		dense_matrix::ptr lres = mat1->mapply2(*mat2, op);
		if (lres == NULL)
			return dense_matrix::ptr();

		res_stores[i] = lres->get_raw_store();
	}
	return block_matrix::create(detail::combined_matrix_store::create(
				res_stores, store->store_layout()));
}

dense_matrix::ptr block_matrix::sapply(bulk_uoperate::const_ptr op) const
{
	std::vector<detail::matrix_store::const_ptr> res_stores(
			store->get_num_mats());
	for (size_t i = 0; i < res_stores.size(); i++) {
		dense_matrix::ptr mat = dense_matrix::create(store->get_mat(i));
		mat = mat->sapply(op);
		if (mat == NULL)
			return dense_matrix::ptr();

		res_stores[i] = mat->get_raw_store();
	}
	return block_matrix::create(detail::combined_matrix_store::create(
				res_stores, store->store_layout()));
}

dense_matrix::ptr block_matrix::apply(matrix_margin margin,
		arr_apply_operate::const_ptr op) const
{
	// When applying the function on the shorter dimension, we have to combine all matrices
	// first. We can't optimize it.
	// TODO we currently don't support applying on the longer dimension.
	return dense_matrix::apply(margin, op);
}

dense_matrix::ptr block_matrix::conv_store(bool in_mem, int num_nodes) const
{
	std::vector<detail::matrix_store::const_ptr> new_stores(store->get_num_mats());
	// TODO it's better to convert them all together.
	for (size_t i = 0; i < store->get_num_mats(); i++) {
		dense_matrix::ptr mat = dense_matrix::create(store->get_mat(i));
		mat = mat->conv_store(in_mem, num_nodes);
		if (mat == NULL)
			return dense_matrix::ptr();
		new_stores[i] = mat->get_raw_store();
	}
	return block_matrix::create(detail::combined_matrix_store::create(
				new_stores, store->store_layout()));
}

bool block_matrix::move_store(bool in_mem, int num_nodes) const
{
	std::vector<detail::matrix_store::const_ptr> new_stores(store->get_num_mats());
	// TODO it's better to convert them all together.
	for (size_t i = 0; i < store->get_num_mats(); i++) {
		dense_matrix::ptr mat = dense_matrix::create(store->get_mat(i));
		mat = mat->conv_store(in_mem, num_nodes);
		if (mat == NULL)
			return false;
		new_stores[i] = mat->get_raw_store();
	}
	const_cast<block_matrix *>(this)->store
		= detail::combined_matrix_store::create(new_stores, store->store_layout());
	return true;
}

namespace
{

class agg_block_sink_store: public detail::sink_store
{
	detail::block_sink_store::const_ptr bsink;
	agg_operate::const_ptr op;

	static bool is_all_in_mem(
			const std::vector<detail::matrix_store::const_ptr> &stores) {
		for (size_t i = 0; i < stores.size(); i++)
			if (!stores[i]->is_in_mem())
				return false;
		return true;
	}
public:
	agg_block_sink_store(
			const std::vector<detail::matrix_store::const_ptr> &stores,
			agg_operate::const_ptr op): detail::sink_store(1, 1,
				is_all_in_mem(stores), op->get_output_type()) {
		bsink = detail::block_sink_store::create(stores, stores.size(), 1);
		this->op = op;
	}

	virtual bool has_materialized() const {
		return bsink->has_materialized();
	}

	virtual matrix_store::const_ptr get_result() const {
		detail::matrix_store::const_ptr part_res = bsink->get_result();
		detail::mem_matrix_store::const_ptr mem_part
			= std::dynamic_pointer_cast<const detail::mem_matrix_store>(part_res);
		assert(mem_part);
		const bulk_operate &combine = op->get_combine();
		detail::mem_matrix_store::ptr ret = detail::mem_matrix_store::create(
				1, 1, matrix_layout_t::L_ROW, combine.get_output_type(), -1);
		combine.runAgg(mem_part->get_num_rows() * mem_part->get_num_cols(),
				mem_part->get_raw_arr(), ret->get_raw_arr());
		return ret;
	}

	virtual std::vector<detail::virtual_matrix_store::const_ptr> get_compute_matrices() const {
		return bsink->get_compute_matrices();
	}
	virtual void materialize_self() const {
		bsink->materialize_self();
	}
	virtual matrix_store::const_ptr materialize(bool in_mem, int num_nodes) const {
		return get_result();
	}

	virtual detail::matrix_store::const_ptr transpose() const {
		return bsink->transpose();
	}

	virtual matrix_layout_t store_layout() const {
		return bsink->store_layout();
	}

	virtual std::string get_name() const {
		return bsink->get_name();
	}
	virtual std::unordered_map<size_t, size_t> get_underlying_mats() const {
		return bsink->get_underlying_mats();
	}
};

}

dense_matrix::ptr block_matrix::aggregate(matrix_margin margin,
			agg_operate::const_ptr op) const
{
	bool is_wide = store->get_mat_ref(0).is_wide();
	// TODO we might want to optimize the first case specially because
	// the current implementation still cause a lot of CPU cache misses.
	if ((margin == matrix_margin::MAR_ROW && !is_wide)
			|| (margin == matrix_margin::MAR_COL && is_wide)) {
		// If the agg operation doesn't have combine, we have to concatenate
		// rows first before running agg on it.
		if (!op->has_combine())
			return dense_matrix::aggregate(margin, op);

		std::vector<dense_matrix::ptr> vecs(store->get_num_mats());
		for (size_t i = 0; i < store->get_num_mats(); i++) {
			dense_matrix::ptr mat = dense_matrix::create(store->get_mat(i));
			// If the matrix has only one row or one column, we don't need
			// to run agg on it.
			if (mat->is_wide() && mat->get_num_rows() == 1) {
				// TODO What if agg and combine are different?
				assert(op->is_same());
				// we want a single col matrix.
				vecs[i] = mat->transpose();
			}
			else if (!mat->is_wide() && mat->get_num_cols() == 1) {
				// TODO What if agg and combine are different?
				assert(op->is_same());
				vecs[i] = mat;
			}
			else
				vecs[i] = mat->aggregate(margin, op);
		}
		// aggregate always returns a column matrix.
		dense_matrix::ptr combined = dense_matrix::cbind(vecs);
		return combined->aggregate(matrix_margin::MAR_ROW, op);
	}
	else {
		std::vector<detail::matrix_store::const_ptr> sinks(store->get_num_mats());
		for (size_t i = 0; i < store->get_num_mats(); i++) {
			dense_matrix::ptr mat = dense_matrix::create(store->get_mat(i));
			dense_matrix::ptr res = mat->aggregate(margin, op);
			sinks[i] = res->get_raw_store();
		}
		detail::matrix_store::ptr ret;
		if (margin == matrix_margin::BOTH) {
			ret = detail::matrix_store::ptr(new agg_block_sink_store(sinks, op));
		}
		else if (margin == matrix_margin::MAR_ROW)
			ret = detail::block_sink_store::create(sinks, sinks.size(), 1);
		else
			ret = detail::block_sink_store::create(sinks, 1, sinks.size());
		return dense_matrix::create(ret);
	}
}

}
