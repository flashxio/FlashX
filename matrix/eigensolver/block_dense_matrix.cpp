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

#include <cblas.h>

#include "block_dense_matrix.h"
#include "dotp_matrix_store.h"
#include "matrix_stats.h"
#include "collected_col_matrix_store.h"

namespace fm
{

namespace eigen
{

bool cache_recent = true;
dense_matrix::ptr cached_mat;

size_t num_col_writes = 0;
size_t num_col_writes_concept = 0;
size_t num_col_reads_concept = 0;
size_t num_multiply_concept = 0;

class mirror_block_multi_vector: public block_multi_vector
{
	// When mirroring the original matrix, we need to make sure that
	// the pointers point to the original matrix.
	mirror_block_multi_vector(const std::vector<fm::dense_matrix::ptr> &mats,
			bool in_mem): block_multi_vector(mats, in_mem) {
	}
public:
	virtual void set_block(off_t block_idx, fm::dense_matrix::const_ptr mat) {
		assert(mat->get_num_cols() == get_block_size());
		mats[block_idx]->assign(*mat);
	}

	static ptr create(const std::vector<fm::dense_matrix::ptr> &mats,
			bool in_mem) {
		return ptr(new mirror_block_multi_vector(mats, in_mem));
	}
};

/*
 * This is a simply way of replacing some columns in a matrix with columns
 * from another matrix.
 */
class set_cols_op: public detail::portion_mapply_op
{
	std::vector<off_t> offs_in_mirror;
public:
	set_cols_op(const std::vector<off_t> &offs_in_mirror, size_t num_rows,
			size_t num_cols, const scalar_type &type): detail::portion_mapply_op(
				num_rows, num_cols, type) {
		this->offs_in_mirror = offs_in_mirror;
	}

	virtual portion_mapply_op::const_ptr transpose() const;
	virtual void run(const std::vector<detail::local_matrix_store::const_ptr> &ins,
			detail::local_matrix_store &out) const;
	virtual std::string to_string(
			const std::vector<detail::matrix_store::const_ptr> &mats) const {
		assert(mats.size() == 2);
		return std::string("set_cols(") + mats[0]->get_name()
			+ ", " + mats[1]->get_name() + ")";
	}
};

class set_rows_op: public detail::portion_mapply_op
{
	std::vector<off_t> offs_in_mirror;
public:
	set_rows_op(const std::vector<off_t> &offs_in_mirror, size_t num_rows,
			size_t num_cols, const scalar_type &type): detail::portion_mapply_op(
				num_rows, num_cols, type) {
		this->offs_in_mirror = offs_in_mirror;
	}

	virtual portion_mapply_op::const_ptr transpose() const;
	virtual void run(const std::vector<detail::local_matrix_store::const_ptr> &ins,
			detail::local_matrix_store &out) const;
	virtual std::string to_string(
			const std::vector<detail::matrix_store::const_ptr> &mats) const {
		assert(mats.size() == 2);
		return std::string("set_rows(") + mats[0]->get_name()
			+ ", " + mats[1]->get_name() + ")";
	}
};

void set_cols_op::run(
		const std::vector<detail::local_matrix_store::const_ptr> &ins,
		detail::local_matrix_store &out) const
{
	assert(ins.size() == 2);
	assert(ins[0]->store_layout() == matrix_layout_t::L_COL);
	assert(ins[1]->store_layout() == matrix_layout_t::L_COL);
	assert(out.store_layout() == matrix_layout_t::L_COL);
	const detail::local_col_matrix_store &col_in1
		= static_cast<const detail::local_col_matrix_store &>(*ins[1]);
	detail::local_col_matrix_store &col_out
		= static_cast<detail::local_col_matrix_store &>(out);
	out.copy_from(*ins[0]);
	for (size_t i = 0; i < offs_in_mirror.size(); i++)
		memcpy(col_out.get_col(offs_in_mirror[i]), col_in1.get_col(i),
				col_out.get_num_rows() * col_out.get_entry_size());
}

void set_rows_op::run(
		const std::vector<detail::local_matrix_store::const_ptr> &ins,
		detail::local_matrix_store &out) const
{
	assert(ins.size() == 2);
	assert(ins[0]->store_layout() == matrix_layout_t::L_ROW);
	assert(ins[1]->store_layout() == matrix_layout_t::L_ROW);
	assert(out.store_layout() == matrix_layout_t::L_ROW);
	const detail::local_row_matrix_store &row_in1
		= static_cast<const detail::local_row_matrix_store &>(*ins[1]);
	detail::local_row_matrix_store &row_out
		= static_cast<detail::local_row_matrix_store &>(out);
	out.copy_from(*ins[0]);
	for (size_t i = 0; i < offs_in_mirror.size(); i++)
		memcpy(row_out.get_row(offs_in_mirror[i]), row_in1.get_row(i),
				row_out.get_num_cols() * row_out.get_entry_size());
}

detail::portion_mapply_op::const_ptr set_cols_op::transpose() const
{
	return detail::portion_mapply_op::const_ptr(new set_rows_op(offs_in_mirror,
				get_out_num_cols(), get_out_num_rows(), get_output_type()));
}

detail::portion_mapply_op::const_ptr set_rows_op::transpose() const
{
	return detail::portion_mapply_op::const_ptr(new set_cols_op(offs_in_mirror,
				get_out_num_cols(), get_out_num_rows(), get_output_type()));
}

/*
 * This mirrors some columns in the block of a block_multi_vector.
 */
class mirror_cols_block_multi_vector: public block_multi_vector
{
	fm::dense_matrix::ptr mirrored_mat;
	std::vector<off_t> offs_in_mirror;

	mirror_cols_block_multi_vector(
			const std::vector<fm::dense_matrix::ptr> &mats,
			const std::vector<off_t> &offs_in_mirror,
			fm::dense_matrix::ptr mirrored_mat, bool in_mem): block_multi_vector(
				mats, in_mem) {
		this->mirrored_mat = mirrored_mat;
		this->offs_in_mirror = offs_in_mirror;
	}
public:
	static ptr create(fm::dense_matrix::ptr mat, const std::vector<off_t> &offs,
			fm::dense_matrix::ptr mirrored_mat, bool in_mem) {
		std::vector<fm::dense_matrix::ptr> mats(1);
		mats[0] = mat;
		return ptr(new mirror_cols_block_multi_vector(mats, offs, mirrored_mat,
					in_mem));
	}

	virtual void set_block(off_t block_idx, fm::dense_matrix::const_ptr mat);
};

void mirror_cols_block_multi_vector::set_block(off_t block_idx,
		fm::dense_matrix::const_ptr mat)
{
	assert(block_idx == 0);
	assert(mat->get_num_cols() == get_block_size());
	mats[block_idx]->assign(*mat);
	if (mirrored_mat->is_virtual()) {
		num_col_writes += mirrored_mat->get_num_cols();
		printf("set_block1: materialize %s\n",
				mirrored_mat->get_data().get_name().c_str());
		detail::matrix_stats_t orig_stats = detail::matrix_stats;
		mirrored_mat->materialize_self();
		detail::matrix_stats.print_diff(orig_stats);
	}
	if (mat->is_virtual()) {
		num_col_writes += mat->get_num_cols();
		printf("set_block2: materialize %s\n",
				mat->get_data().get_name().c_str());
		detail::matrix_stats_t orig_stats = detail::matrix_stats;
		mat->materialize_self();
		detail::matrix_stats.print_diff(orig_stats);
	}

	std::vector<fm::dense_matrix::const_ptr> in_mats(2);
	in_mats[0] = mirrored_mat;
	in_mats[1] = mat;
	dense_matrix::ptr new_mat = fm::detail::mapply_portion(in_mats,
			detail::portion_mapply_op::const_ptr(new set_cols_op(
					offs_in_mirror, mirrored_mat->get_num_rows(),
					mirrored_mat->get_num_cols(), mirrored_mat->get_type())),
			mirrored_mat->store_layout());
	mirrored_mat->assign(*new_mat);
}

block_multi_vector::block_multi_vector(
		const std::vector<fm::dense_matrix::ptr> &mats,
		bool in_mem): type(mats[0]->get_type())
{
	this->MAX_MUL_BLOCKS = 8;
	this->in_mem = in_mem;
	this->num_rows = mats[0]->get_num_rows();
	this->num_cols = mats[0]->get_num_cols() * mats.size();
	this->block_size = mats[0]->get_num_cols();
	this->mats = mats;
	for (size_t i = 1; i < mats.size(); i++) {
		assert(this->block_size == mats[i]->get_num_cols());
		assert(this->num_rows == mats[i]->get_num_rows());
	}
}

block_multi_vector::block_multi_vector(size_t nrow, size_t ncol,
		size_t block_size, const fm::scalar_type &_type, bool in_mem): type(_type)
{
	this->MAX_MUL_BLOCKS = 8;
	this->in_mem = in_mem;
	this->num_rows = nrow;
	this->num_cols = ncol;
	this->block_size = block_size;
	mats.resize(ncol / block_size);
	for (size_t i = 0; i < mats.size(); i++)
		mats[i] = dense_matrix::create(nrow, block_size,
				matrix_layout_t::L_COL, type, get_num_nodes(), in_mem);
}

dense_matrix::const_ptr block_multi_vector::get_col(off_t col_idx) const
{
	off_t block_idx = col_idx / block_size;
	off_t local_col_idx = col_idx % block_size;
	std::vector<off_t> offs(1);
	offs[0] = local_col_idx;
	fm::dense_matrix::const_ptr block = get_block(block_idx);
	fm::dense_matrix::const_ptr ret;
	if (block->is_virtual()) {
		// We need to handle the special case explicitly.
		const dotp_matrix_store *dotp
			= dynamic_cast<const dotp_matrix_store *>(
					block->get_raw_store().get());
		if (dotp)
			ret = dense_matrix::create(dotp->get_cols(offs));
		else {
			num_col_writes += block->get_num_cols();
			printf("get_col: materialize %s\n",
					block->get_data().get_name().c_str());
			detail::matrix_stats_t orig_stats = detail::matrix_stats;
			block->materialize_self();
			detail::matrix_stats.print_diff(orig_stats);
			ret = block->get_cols(offs);
		}
	}
	else
		ret = block->get_cols(offs);
	assert(ret);
	return ret;
}

bool is_same_block(const std::vector<int> &index, size_t block_size)
{
	if (index.size() > block_size)
		return false;
	size_t block_start = index[0] / block_size;
	for (size_t i = 1; i < index.size(); i++) {
		if ((index[i] / block_size) != block_start)
			return false;
	}
	return true;
}

block_multi_vector::ptr block_multi_vector::get_cols(const std::vector<int> &index)
{
	// get entire blocks.
	if (index.size() % get_block_size() == 0 && index[0] % get_block_size() == 0) {
		block_multi_vector::ptr ret = block_multi_vector::create(get_num_rows(),
				index.size(), get_block_size(), get_type(), in_mem);
		size_t num_blocks = index.size() / get_block_size();
		size_t block_start = index[0] / get_block_size();
		for (size_t i = 0; i < num_blocks; i++)
			ret->set_block(i, get_block(i + block_start));
		return ret;
	}
	else if (is_same_block(index, block_size)) {
		// We can only fetch columns from the same block.
		size_t block_start = index[0] / get_block_size();
		std::vector<off_t> local_offs(index.size());
		for (size_t i = 0; i < local_offs.size(); i++)
			local_offs[i] = index[i] - block_start * get_block_size();

		block_multi_vector::ptr ret = block_multi_vector::create(get_num_rows(),
				index.size(), index.size(), get_type(), in_mem);
		dense_matrix::ptr block = get_block(block_start);
		dense_matrix::ptr ret1;
		if (block->is_virtual()) {
			// We need to handle the special case explicitly.
			const dotp_matrix_store *dotp
				= dynamic_cast<const dotp_matrix_store *>(
					block->get_raw_store().get());
			if (dotp)
				ret1 = dense_matrix::create(dotp->get_cols(local_offs));
			else {
				num_col_writes += block->get_num_cols();
				printf("get_cols: materialize %s\n",
						block->get_data().get_name().c_str());
				detail::matrix_stats_t orig_stats = detail::matrix_stats;
				block->materialize_self();
				detail::matrix_stats.print_diff(orig_stats);
				ret1 = block->get_cols(local_offs);
			}
		}
		else
			ret1 = block->get_cols(local_offs);
		ret->set_block(0, ret1);
		return ret;
	}
	else {
		block_multi_vector::ptr ret = block_multi_vector::create(get_num_rows(),
				index.size(), 1, get_type(), in_mem);
		for (size_t i = 0; i < index.size(); i++)
			ret->set_block(i, get_col(index[i]));
		return ret;
	}
}

block_multi_vector::ptr block_multi_vector::get_cols_mirror(
		const std::vector<int> &index)
{
	// get entire blocks.
	if (index.size() % get_block_size() == 0 && index[0] % get_block_size() == 0) {
		size_t num_blocks = index.size() / get_block_size();
		size_t block_start = index[0] / get_block_size();
		std::vector<fm::dense_matrix::ptr> mats(num_blocks);
		for (size_t i = 0; i < num_blocks; i++)
			mats[i] = get_block(i + block_start);
		mirror_block_multi_vector::ptr ret
			= mirror_block_multi_vector::create(mats, in_mem);
		// We want the dense matrices in the basis materialized.
		// This can reduce significant computation.
		if (num_blocks == 1) {
			if (!in_mem && cached_mat && cached_mat->is_virtual()) {
				printf("materialize the old cached mat %s to disks\n",
						cached_mat->get_data().get_name().c_str());
				num_col_writes += cached_mat->get_num_cols();
				detail::matrix_stats_t orig_stats = detail::matrix_stats;
				bool ret = cached_mat->move_store(false, -1);
				assert(ret);
				detail::matrix_stats.print_diff(orig_stats);
			}
			else if (cached_mat && cached_mat->is_virtual()) {
				printf("materialize the old cached mat %s\n",
						cached_mat->get_data().get_name().c_str());
				detail::matrix_stats_t orig_stats = detail::matrix_stats;
				cached_mat->materialize_self();
				detail::matrix_stats.print_diff(orig_stats);
			}
			cached_mat = mats[0];
		}
		return ret;
	}
	else if (is_same_block(index, get_block_size())) {
		size_t block_start = index[0] / get_block_size();

		// Get the relative offsets in a block.
		std::vector<off_t> idxs_in_block(index.size());
		idxs_in_block[0] = index[0] % get_block_size();
		for (size_t i = 1; i < index.size(); i++) {
			// All columns has to be in the same block.
			assert(block_start == index[i] / get_block_size());
			idxs_in_block[i] = index[i] % get_block_size();
		}

		// Get the columns in the block.
		fm::dense_matrix::ptr block = get_block(block_start);
		if (block->is_virtual()) {
			printf("get_cols_mirror: materialize %s\n",
					block->get_data().get_name().c_str());
			num_col_writes += block->get_num_cols();
			block->materialize_self();
		}
		mirror_cols_block_multi_vector::ptr ret
			= mirror_cols_block_multi_vector::create(
					block->get_cols(idxs_in_block), idxs_in_block, block,
					in_mem);
		return ret;
	}
	else {
		assert(0);
		return block_multi_vector::ptr();
	}
#if 0
	else {
		// otherwise, we can only fetch columns from the same block.
		assert(is_same_block(index, block_size));
		size_t block_start = index[0] / get_block_size();
		std::vector<off_t> local_offs(index.size());
		for (size_t i = 0; i < local_offs.size(); i++)
			local_offs[i] = index[i] - block_start * get_block_size();

		std::vector<fm::dense_matrix::ptr> mats(1);
		mats[0] = get_block(block_start)->get_cols(local_offs);
		mirror_block_multi_vector::ptr ret = mirror_block_multi_vector::create(mats);
		return ret;
	}
#endif
}

std::atomic<size_t> num_spmm;

void block_multi_vector::sparse_matrix_multiply(const spm_function &multiply,
		const block_multi_vector &X, block_multi_vector &Y)
{
	assert(multiply.get_num_cols() == X.get_num_rows());
	assert(multiply.get_num_rows() == Y.get_num_rows());
	assert(Y.get_num_cols() == X.get_num_cols());
	assert(Y.get_num_blocks() == X.get_num_blocks());
	assert(X.get_block_size() == Y.get_block_size());
	assert(X.get_type() == Y.get_type());
	size_t num_blocks = X.get_num_blocks();
	int num_nodes = matrix_conf.get_num_nodes();
	bool in_mem = X.in_mem;
	printf("SpMM %ld: input matirx is in-mem: %d\n", num_spmm.load(), in_mem);
	num_spmm++;
	for (size_t i = 0; i < num_blocks; i++) {
		dense_matrix::ptr in = X.get_block(i);
		dense_matrix::ptr res;
		num_col_reads_concept += in->get_num_cols();
		// If we run in the EM mode, we should materialize the cached matrix
		// and write out the most recently cached matrix.
		if (!in_mem && in->is_virtual() && ((cached_mat == in)
					|| (cached_mat
						&& cached_mat->get_raw_store() == in->get_raw_store()))) {
			printf("spmm: materialize in mat %s to disks\n",
					in->get_data().get_name().c_str());
			num_col_writes += cached_mat->get_num_cols();
			detail::matrix_stats_t orig_stats = detail::matrix_stats;
			bool ret = cached_mat->move_store(false, -1);
			assert(ret);
			detail::matrix_stats.print_diff(orig_stats);
		}
		// Otherwise, we still want to materialize the matrix.
		else if (in->is_virtual() && cached_mat
				&& cached_mat->get_raw_store() == in->get_raw_store()) {
			printf("spmm: materialize in mat %s\n",
					in->get_data().get_name().c_str());
			detail::matrix_stats_t orig_stats = detail::matrix_stats;
			in->materialize_self();
			detail::matrix_stats.print_diff(orig_stats);
		}

		assert(in->store_layout() == matrix_layout_t::L_COL);
		dense_matrix::ptr row_in = in->conv2(matrix_layout_t::L_ROW);
		// If the input matrix isn't in memory, we should load it first.
		// When loading, if the input matrix is a virtual matrix, materialize
		// it in memory.
		if (!row_in->is_in_mem()) {
			printf("load the input matrix for SpMM to memory\n");
			detail::matrix_stats_t orig_stats = detail::matrix_stats;
			row_in = row_in->conv_store(true, num_nodes);
			detail::matrix_stats.print_diff(orig_stats);
		}
		// If the input matrix is in memory and is virtual, we should
		// materialize it.
		else if (row_in->is_virtual()) {
			printf("spmm: materialize %s\n",
					row_in->get_data().get_name().c_str());
			detail::matrix_stats_t orig_stats = detail::matrix_stats;
			row_in->materialize_self();
			detail::matrix_stats.print_diff(orig_stats);
		}
		res = multiply.run(row_in);
		if (res->store_layout() == matrix_layout_t::L_ROW)
			res = res->conv2(matrix_layout_t::L_COL);
		// If the input matrix isn't in memory, we should convert it to
		// EM matrix.
		if (!in_mem)
			num_col_writes_concept += res->get_num_cols();
		if (!in_mem && !cache_recent) {
			printf("write the output matrix of SpMM to disks\n");
			num_col_writes += res->get_num_cols();
			detail::matrix_stats_t orig_stats = detail::matrix_stats;
			res = res->conv_store(false, -1);
			detail::matrix_stats.print_diff(orig_stats);
		}
		assert(res->store_layout() == matrix_layout_t::L_COL);
		Y.set_block(i, res);
	}
}

block_multi_vector::ptr block_multi_vector::clone() const
{
	block_multi_vector::ptr vecs= block_multi_vector::create(
			get_num_rows(), get_num_cols(), block_size, type, in_mem);
	size_t num_blocks = get_num_blocks();
	for (size_t i = 0; i < num_blocks; i++)
		vecs->mats[i] = this->get_block(i)->clone();
	return vecs;
}

namespace
{

template<class T>
class gemm_op: public fm::detail::portion_mapply_op
{
	T alpha;
	T beta;
	size_t A_num_blocks;
	size_t C_num_blocks;
	detail::mem_col_matrix_store::const_ptr Bstore;
public:
	gemm_op(detail::mem_col_matrix_store::const_ptr B, size_t A_num_blocks,
			size_t C_num_blocks, T alpha, T beta, size_t out_num_rows,
			size_t out_num_cols): fm::detail::portion_mapply_op(
				out_num_rows, out_num_cols, get_scalar_type<T>()) {
		this->alpha = alpha;
		this->beta = beta;
		this->A_num_blocks = A_num_blocks;
		this->C_num_blocks = C_num_blocks;
		this->Bstore = B;
	}

	virtual void run(
			const std::vector<fm::detail::local_matrix_store::const_ptr> &ins,
			fm::detail::local_matrix_store &out) const;
	virtual fm::detail::portion_mapply_op::const_ptr transpose() const;

	virtual std::string to_string(
			const std::vector<detail::matrix_store::const_ptr> &mats) const {
		std::string str;
		if (A_num_blocks == 1)
			str = mats[0]->get_name();
		else {
			str = "cat(";
			for (size_t i = 0; i < A_num_blocks - 1; i++)
				str += mats[i]->get_name() + ",";
			str += mats[A_num_blocks - 1]->get_name() + ")";
		}
		if (mats.size() == A_num_blocks)
			return (boost::format("(%1% * %2% * %3%)") % alpha % str
					% Bstore->get_name()).str();
		else
			return (boost::format("(%1% * %2% * %3% + %4% * %5%)") % alpha % str
					% Bstore->get_name() % beta % mats[A_num_blocks]->get_name()).str();
	}
};

template<class T>
class t_gemm_op: public fm::detail::portion_mapply_op
{
	gemm_op<T> op;
public:
	t_gemm_op(const gemm_op<T> &_op): fm::detail::portion_mapply_op(
			_op.get_out_num_cols(), _op.get_out_num_rows(),
			_op.get_output_type()), op(_op) {
	}

	virtual void run(
			const std::vector<fm::detail::local_matrix_store::const_ptr> &ins,
			fm::detail::local_matrix_store &out) const {
		std::vector<fm::detail::local_matrix_store::const_ptr> t_ins(ins.size());
		for (size_t i = 0; i < ins.size(); i++)
			t_ins[i] = std::static_pointer_cast<const fm::detail::local_matrix_store>(
					ins[i]->transpose());
		fm::detail::local_matrix_store::ptr t_out
			= std::static_pointer_cast<fm::detail::local_matrix_store>(
					out.transpose());
		op.run(t_ins, *t_out);
	}
	virtual fm::detail::portion_mapply_op::const_ptr transpose() const {
		return fm::detail::portion_mapply_op::const_ptr(new gemm_op<T>(op));
	}
	virtual std::string to_string(
			const std::vector<detail::matrix_store::const_ptr> &mats) const {
		return op.to_string(mats);
	}
};

template<class T>
class sum_op: public fm::detail::portion_mapply_op
{
public:
	sum_op(size_t out_num_rows, size_t out_num_cols): fm::detail::portion_mapply_op(
				out_num_rows, out_num_cols, get_scalar_type<T>()) {
	}

	virtual void run(
			const std::vector<fm::detail::local_matrix_store::const_ptr> &ins,
			fm::detail::local_matrix_store &out) const {
		assert(!ins.empty());
		if (ins.size() == 1)
			out.copy_from(*ins[0]);
		else {
			mapply2(*ins[0], *ins[1],
					out.get_type().get_basic_ops().get_add(), out);
			for (size_t i = 2; i < ins.size(); i++)
				mapply2(*ins[i], out,
						out.get_type().get_basic_ops().get_add(), out);
		}
	}

	virtual fm::detail::portion_mapply_op::const_ptr transpose() const;

	virtual std::string to_string(
			const std::vector<detail::matrix_store::const_ptr> &mats) const {
		std::string str;
		assert(mats.size() > 1);
		str = "(";
		for (size_t i = 0; i < mats.size() - 1; i++)
			str += mats[i]->get_name() + "+";
		str += mats[mats.size() - 1]->get_name() + ")";
		return str;
	}

	virtual bool is_agg() const {
		return true;
	}
};

template<class T>
class t_sum_op: public fm::detail::portion_mapply_op
{
	sum_op<T> op;
public:
	t_sum_op(const sum_op<T> &_op): fm::detail::portion_mapply_op(
			_op.get_out_num_cols(), _op.get_out_num_rows(),
			_op.get_output_type()), op(_op) {
	}

	virtual void run(
			const std::vector<fm::detail::local_matrix_store::const_ptr> &ins,
			fm::detail::local_matrix_store &out) const {
		std::vector<fm::detail::local_matrix_store::const_ptr> t_ins(ins.size());
		for (size_t i = 0; i < ins.size(); i++)
			t_ins[i] = std::static_pointer_cast<const fm::detail::local_matrix_store>(
					ins[i]->transpose());
		fm::detail::local_matrix_store::ptr t_out
			= std::static_pointer_cast<fm::detail::local_matrix_store>(
					out.transpose());
		op.run(t_ins, *t_out);
	}
	virtual fm::detail::portion_mapply_op::const_ptr transpose() const {
		return fm::detail::portion_mapply_op::const_ptr(new sum_op<T>(op));
	}
	virtual std::string to_string(
			const std::vector<detail::matrix_store::const_ptr> &mats) const {
		return op.to_string(mats);
	}

	virtual bool is_agg() const {
		return true;
	}
};

template<class T>
fm::detail::portion_mapply_op::const_ptr gemm_op<T>::transpose() const
{
	return fm::detail::portion_mapply_op::const_ptr(new t_gemm_op<T>(*this));
}

template<class T>
fm::detail::portion_mapply_op::const_ptr sum_op<T>::transpose() const
{
	return fm::detail::portion_mapply_op::const_ptr(new t_sum_op<T>(*this));
}

typedef std::vector<fm::detail::local_matrix_store::const_ptr>::const_iterator block_iterator;

static void copy_from_blocks(block_iterator begin, block_iterator end,
		fm::detail::local_col_matrix_store &res_store)
{
	std::vector<const fm::detail::local_col_matrix_store *> col_ins(
			end - begin);
	size_t i = 0;
	for (auto it = begin; it != end; it++, i++)
		col_ins[i] = static_cast<const fm::detail::local_col_matrix_store *>(
				it->get());
	size_t type_size = res_store.get_type().get_size();
	for (size_t i = 0; i < col_ins.size(); i++) {
		assert(res_store.get_global_start_row()
				== col_ins[i]->get_global_start_row());
		assert(res_store.get_global_start_col()
				== col_ins[i]->get_global_start_col());
		assert(res_store.get_type() == col_ins[i]->get_type());
		assert(res_store.get_num_rows() == col_ins[i]->get_num_rows());

		size_t block_size = col_ins[i]->get_num_cols();
		for (size_t j = 0; j < block_size; j++) {
			size_t col_idx = i * block_size + j;
			memcpy(res_store.get_col(col_idx), col_ins[i]->get_col(j),
					res_store.get_num_rows() * type_size);
		}
	}
}

template<class T>
void gemm_op<T>::run(
		const std::vector<fm::detail::local_matrix_store::const_ptr> &ins,
		fm::detail::local_matrix_store &out) const
{
	detail::matrix_stats.inc_multiplies(
			ins[0]->get_num_rows() * Bstore->get_num_cols() * Bstore->get_num_cols());

	assert(A_num_blocks + C_num_blocks == ins.size());
	off_t global_start_row = ins.front()->get_global_start_row();
	off_t global_start_col = ins.front()->get_global_start_col();

	T *res_mat;
	res_mat = (T *) out.get_raw_arr();
	fm::detail::local_matrix_store::ptr tmp_res;
	if (res_mat == NULL) {
		fm::detail::local_col_matrix_store *raw_tmp_res
			= new detail::local_buf_col_matrix_store(
					global_start_row, global_start_col,
					out.get_num_rows(), out.get_num_cols(),
					get_scalar_type<T>(), -1);
		// TODO we don't need to allocate this every time.
		tmp_res = fm::detail::local_col_matrix_store::ptr(raw_tmp_res);
		if (beta && C_num_blocks == 1)
			tmp_res->copy_from(*ins.back());
		else if (beta)
			copy_from_blocks(ins.begin() + A_num_blocks, ins.end(),
					*raw_tmp_res);
		res_mat = (T *) tmp_res->get_raw_arr();
	}
	else {
		if (beta && C_num_blocks == 1)
			out.copy_from(*ins.back());
		else if (beta)
			copy_from_blocks(ins.begin() + A_num_blocks, ins.end(),
					static_cast<fm::detail::local_col_matrix_store &>(out));
	}
	assert(res_mat);

	fm::detail::local_col_matrix_store::const_ptr Astore;
	// We need to merge multiple matrix portions.
	if (A_num_blocks >= 2) {
		size_t block_size = ins.front()->get_num_cols();
		size_t num_rows = ins.front()->get_num_rows();
		size_t num_cols = block_size * A_num_blocks;
		assert(num_cols == Bstore->get_num_rows());
		// TODO we don't need to allocate this every time.
		fm::detail::local_col_matrix_store::ptr in_buf(
				new fm::detail::local_buf_col_matrix_store(
					global_start_row, global_start_col,
					num_rows, num_cols, get_scalar_type<T>(), -1));
		copy_from_blocks(ins.begin(), ins.begin() + A_num_blocks, *in_buf);
		Astore = in_buf;
	}
	// All data in this portion isn't contiguous.
	else if (ins[0]->get_raw_arr() == NULL) {
		size_t num_rows = ins.front()->get_num_rows();
		size_t num_cols = ins.front()->get_num_cols();
		// TODO we don't need to allocate this every time.
		fm::detail::local_col_matrix_store::ptr in_buf(
				new fm::detail::local_buf_col_matrix_store(
					global_start_row, global_start_col,
					num_rows, num_cols, get_scalar_type<T>(), -1));
		in_buf->copy_from(*ins.front());
		Astore = in_buf;
	}
	else
		Astore = fm::detail::local_col_matrix_store::cast(ins[0]);

	const T *Amat = (const T *) Astore->get_raw_arr();
	assert(Bstore->get_data().get_num_bytes()
			== Bstore->get_num_rows() * Bstore->get_num_cols() *
			Bstore->get_entry_size());
	const T *Bmat = (const T *) Bstore->get_data().get_raw();
	assert(Amat);
	assert(Bmat);
	cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans,
			Astore->get_num_rows(), Bstore->get_num_cols(),
			Astore->get_num_cols(), alpha, Amat,
			Astore->get_num_rows(), Bmat, Bstore->get_num_rows(),
			beta, res_mat, out.get_num_rows());
	if (tmp_res)
		out.copy_from(*tmp_res);
}

/*
 * This split a col-major matrix into multiple col-major matrices.
 */
class split_op: public detail::portion_mapply_op
{
public:
	split_op(): portion_mapply_op(0, 0, get_scalar_type<int>()) {
	}
	void run(const std::vector<detail::local_matrix_store::const_ptr> &ins,
			const std::vector<detail::local_matrix_store::ptr> &outs) const;

	virtual portion_mapply_op::const_ptr transpose() const {
		return portion_mapply_op::const_ptr();
	}
	virtual std::string to_string(
			const std::vector<detail::matrix_store::const_ptr> &mats) const {
		return std::string();
	}
};

void split_op::run(
		const std::vector<detail::local_matrix_store::const_ptr> &ins,
		const std::vector<detail::local_matrix_store::ptr> &outs) const
{
	assert(ins.size() == 1);
	assert(ins[0]->store_layout() == matrix_layout_t::L_COL);
	for (size_t i = 0; i < outs.size(); i++)
		assert(outs[i]->store_layout() == matrix_layout_t::L_COL);
	const detail::local_col_matrix_store &col_in
		= dynamic_cast<const detail::local_col_matrix_store &>(*ins[0]);
	size_t col_idx = 0;
	for (size_t j = 0; j < outs.size(); j++) {
		detail::local_col_matrix_store &col_out
			= dynamic_cast<detail::local_col_matrix_store &>(*outs[j]);
		for (size_t i = 0; i < col_out.get_num_cols(); i++) {
			assert(col_idx < col_in.get_num_cols());
			memcpy(col_out.get_col(i), col_in.get_col(col_idx),
					col_in.get_num_rows() * col_in.get_entry_size());
			col_idx++;
		}
	}
	assert(col_idx == ins[0]->get_num_cols());
}

detail::mem_col_matrix_store::const_ptr get_sub_mat(
		detail::mem_col_matrix_store::const_ptr mat, size_t start_row,
		size_t start_col, size_t num_rows, size_t num_cols)
{
	detail::local_matrix_store::const_ptr portion = mat->get_portion(
			start_row, start_col, num_rows, num_cols);
	assert(portion);
	detail::local_col_matrix_store::const_ptr col_portion
		= detail::local_col_matrix_store::cast(portion);
	assert(col_portion);

	detail::mem_col_matrix_store::ptr ret = detail::mem_col_matrix_store::create(
			num_rows, num_cols, col_portion->get_type());
	for (size_t i = 0; i < num_cols; i++)
		memcpy(ret->get_col(i), col_portion->get_col(i),
				num_rows * ret->get_entry_size());
	return ret;
}

}

static void set_caching(
		const std::vector<detail::matrix_store::const_ptr> &mats, bool caching)
{
	for (size_t i = 0; i < mats.size(); i++)
		const_cast<detail::matrix_store &>(*mats[i]).set_cache_portion(caching);
}

static std::vector<detail::matrix_store::const_ptr> get_input_matrices(
		const block_multi_vector &mv)
{
	std::vector<detail::matrix_store::const_ptr> mats(mv.get_num_blocks());
	for (size_t i = 0; i < mv.get_num_blocks(); i++)
		mats[i] = mv.get_block(i)->get_raw_store();
	return mats;
}

block_multi_vector::ptr block_multi_vector::gemm(const block_multi_vector &A,
		detail::mem_col_matrix_store::const_ptr B, const scalar_variable &alpha,
		const scalar_variable &beta) const
{
	if (A.get_num_cols() != B->get_num_rows()
			|| A.get_num_rows() != get_num_rows()
			|| B->get_num_cols() != get_num_cols()) {
		BOOST_LOG_TRIVIAL(error) << "Wrong dimension for GEMM";
		return block_multi_vector::ptr();
	}

	assert(A.get_num_rows() == this->get_num_rows());

	double d_alpha
		= dynamic_cast<const scalar_variable_impl<double> &>(alpha).get();
	double d_beta
		= dynamic_cast<const scalar_variable_impl<double> &>(beta).get();

	// This works for double float points.
	assert(type == fm::get_scalar_type<double>());
	size_t A_num_blocks = A.get_num_blocks();
	size_t C_num_blocks = d_beta != 0 ? this->get_num_blocks() : 0;
	dense_matrix::ptr block;
	bool use_hierarchy = false;
	if (A_num_blocks <= MAX_MUL_BLOCKS) {
		std::vector<dense_matrix::const_ptr> mats;
		if (d_beta) {
			mats.resize(A.get_num_blocks() + this->get_num_blocks());
			for (size_t i = A.get_num_blocks(); i < mats.size(); i++)
				mats[i] = this->get_block(i - A.get_num_blocks());
		}
		else
			mats.resize(A.get_num_blocks());
		for (size_t i = 0; i < A.get_num_blocks(); i++)
			mats[i] = A.get_block(i);
		assert(A_num_blocks + C_num_blocks == mats.size());

		// I assume B is small enough.
		detail::portion_mapply_op::const_ptr op(new gemm_op<double>(B,
					A_num_blocks, C_num_blocks, d_alpha, d_beta,
					this->get_num_rows(), this->get_num_cols()));
		block = mapply_portion(mats, op, matrix_layout_t::L_COL);
	}
	else {
		printf("compute gemm hierarchically\n");
		use_hierarchy = true;
		// If there are too many blocks in the subspace, we need to perform
		// gemm in a hierarchical fashion.
		std::vector<dense_matrix::const_ptr> tmp_mats;
		// We first perform multiplication on the subset of the matrices
		// and generate some temp matrices.
		for (size_t i = 0; i < A_num_blocks; i += MAX_MUL_BLOCKS) {
			size_t num_sub_mats = std::min(MAX_MUL_BLOCKS, A_num_blocks - i);
			std::vector<detail::matrix_store::const_ptr> sub_mats(num_sub_mats);
			for (size_t j = 0; j < num_sub_mats; j++)
				sub_mats[j] = A.get_block(i + j)->get_raw_store();
			detail::portion_mapply_op::const_ptr op;
			detail::mem_col_matrix_store::const_ptr sub_B = get_sub_mat(B,
					i * A.get_block_size(), 0, num_sub_mats * A.get_block_size(),
					B->get_num_cols());
			if (i + MAX_MUL_BLOCKS < A_num_blocks)
				op = detail::portion_mapply_op::const_ptr(new gemm_op<double>(
							sub_B, sub_mats.size(), 0, d_alpha, 0,
							this->get_num_rows(), this->get_num_cols()));
			else {
				// This is the last group.
				size_t num_sub_mats = sub_mats.size();
				if (d_beta) {
					sub_mats.resize(num_sub_mats + C_num_blocks);
					for (size_t j = num_sub_mats; j < sub_mats.size(); j++)
						sub_mats[j] = this->get_block(
								j - num_sub_mats)->get_raw_store();
				}
				op = detail::portion_mapply_op::const_ptr(new gemm_op<double>(
							sub_B, num_sub_mats, C_num_blocks, d_alpha, d_beta,
							this->get_num_rows(), this->get_num_cols()));
			}
			detail::matrix_store::ptr tmp_store = detail::__mapply_portion_virtual(
					sub_mats, op, matrix_layout_t::L_COL);
			// We really don't need to cache the portion in this intermediate
			// matrix in the hierarchy.
			tmp_store->set_cache_portion(false);
			tmp_mats.push_back(dense_matrix::create(tmp_store));
		}

		// We then sum all of the temp matrices with the original alpha
		// and beta.
		detail::portion_mapply_op::const_ptr op(new sum_op<double>(
					this->get_num_rows(), this->get_num_cols()));
		// We will materialize the mapply in a hierarchical way,
		// so the intermediate matrices should be read from SSDs in serial.
		block = mapply_portion(tmp_mats, op, matrix_layout_t::L_COL, false);
	}

	std::unordered_map<size_t, size_t> bytes
		= block->get_data().get_underlying_mats();
	block_multi_vector::ptr vecs;
	// We are going to assign the result to the current MV. The result
	// of gemm is a single matrix. If the current MV stores data in multiple
	// matrices, its block size will be different from the number of columns
	// of the gemm result. In this case, we need to split the gemm result.
	if (get_block_size() != block->get_num_cols()) {
		size_t num_out_cols = block->get_num_cols();
		assert(num_out_cols % get_block_size() == 0);

		printf("Split matrix\n");
		printf("There are %ld underlying matrices\n", bytes.size());
		num_col_writes += block->get_num_cols();
		printf("materialize %s\n", block->get_data().get_name().c_str());
		detail::matrix_stats_t orig_stats = detail::matrix_stats;
		std::vector<detail::matrix_store::const_ptr> orig_input_mats
			= get_input_matrices(A);
		if (use_hierarchy)
			set_caching(orig_input_mats, false);
		block->materialize_self();
		if (use_hierarchy)
			set_caching(orig_input_mats, true);
		detail::matrix_stats.print_diff(orig_stats);

		vecs = block_multi_vector::create(get_num_rows(), B->get_num_cols(),
				get_block_size(), type, in_mem);
		off_t col_off = 0;
		for (size_t i = 0; i < vecs->mats.size(); i++) {
			std::vector<off_t> idxs(get_block_size());
			for (size_t j = 0; j < get_block_size(); j++)
				idxs[j] = col_off++;
			vecs->mats[i] = block->get_cols(idxs);
		}
	}
	else if (bytes.size() > 2) {
		printf("There are %ld underlying matrices\n", bytes.size());
		printf("materialize %s\n", block->get_data().get_name().c_str());
		detail::matrix_stats_t orig_stats = detail::matrix_stats;
		std::vector<detail::matrix_store::const_ptr> orig_input_mats
			= get_input_matrices(A);
		if (use_hierarchy)
			set_caching(orig_input_mats, false);
		// If we want to cache the most recently materialized matrix.
		if (!in_mem && cache_recent)
			block->move_store(true, matrix_conf.get_num_nodes());
		else {
			num_col_writes += block->get_num_cols();
			block->materialize_self();
		}
		if (use_hierarchy)
			set_caching(orig_input_mats, true);
		detail::matrix_stats.print_diff(orig_stats);

		vecs = block_multi_vector::create(get_num_rows(), B->get_num_cols(),
				B->get_num_cols(), type, in_mem);
		vecs->mats[0] = block;
	}
	else {
		vecs = block_multi_vector::create(get_num_rows(), B->get_num_cols(),
				B->get_num_cols(), type, in_mem);
		vecs->mats[0] = block;
	}
	return vecs;
}

void block_multi_vector::assign(const block_multi_vector &vecs)
{
	assert(num_rows == vecs.get_num_rows());
	assert(num_cols == vecs.get_num_cols());
	assert(type == vecs.get_type());
	if (block_size != vecs.get_block_size()) {
		for (size_t i = 0; i < mats.size(); i++)
			assert(mats[i].use_count() == 1);
		block_size = vecs.get_block_size();
		mats.resize(vecs.get_num_blocks());
	}
	size_t num_blocks = vecs.get_num_blocks();
	for (size_t i = 0; i < num_blocks; i++)
		set_block(i, vecs.get_block(i));
}

block_multi_vector::ptr block_multi_vector::add(
		const block_multi_vector &vecs) const
{
	block_multi_vector::ptr ret= block_multi_vector::create(
			get_num_rows(), get_num_cols(), block_size, type, in_mem);
	size_t num_blocks = get_num_blocks();
	for (size_t i = 0; i < num_blocks; i++)
		ret->mats[i] = this->get_block(i)->add(*vecs.get_block(i));
	return ret;
}

template<class T>
class multiply_wide_op: public detail::portion_mapply_op
{
	std::vector<detail::local_matrix_store::ptr> res_bufs;
	size_t out_num_rows;
	size_t out_num_cols;
	matrix_layout_t Alayout;
	matrix_layout_t Blayout;
public:
	multiply_wide_op(size_t num_threads, size_t out_num_rows,
			size_t out_num_cols, size_t fake_num_rows,
			matrix_layout_t required_layout): detail::portion_mapply_op(
				// We have to output something to make it work for the upper
				// level of the mapply hierarchy.
				fake_num_rows, 1, get_scalar_type<T>()) {
		res_bufs.resize(num_threads);
		this->out_num_rows = out_num_rows;
		this->out_num_cols = out_num_cols;
		// We need to transpose the A matrix, so we want the data in the A matrix
		// to be organized in the opposite layout to the required.
		if (required_layout == matrix_layout_t::L_COL)
			Alayout = matrix_layout_t::L_ROW;
		else
			Alayout = matrix_layout_t::L_COL;
		Blayout = required_layout;
	}

	const std::vector<detail::local_matrix_store::ptr> &get_partial_results(
			) const {
		return res_bufs;
	}

	virtual void run(
			const std::vector<detail::local_matrix_store::const_ptr> &ins,
			detail::local_matrix_store &out) const;

	virtual detail::portion_mapply_op::const_ptr transpose() const {
		assert(0);
		return detail::portion_mapply_op::const_ptr();
	}

	virtual std::string to_string(
			const std::vector<detail::matrix_store::const_ptr> &mats) const {
		assert(mats.size() == 2);
		return std::string("(") + (mats[0]->get_name()
					+ "*") + mats[1]->get_name() + std::string(")");
	}
};

class empty_op: public detail::portion_mapply_op
{
public:
	empty_op(): detail::portion_mapply_op(0, 0,
			get_scalar_type<double>()) {
	}

	virtual void run(
			const std::vector<detail::local_matrix_store::const_ptr> &ins) const {
		for (size_t i = 0; i < ins.size(); i++)
			ins[i]->get_raw_arr();
	}

	virtual detail::portion_mapply_op::const_ptr transpose() const {
		assert(0);
		return detail::portion_mapply_op::const_ptr();
	}

	virtual std::string to_string(
			const std::vector<detail::matrix_store::const_ptr> &mats) const {
		return std::string();
	}

	virtual bool is_agg() const {
		return true;
	}
};

template<class T>
void multiply_wide_op<T>::run(
		const std::vector<detail::local_matrix_store::const_ptr> &ins,
		detail::local_matrix_store &out) const
{
	assert(ins.size() >= 2);
	detail::pool_task_thread *thread = dynamic_cast<detail::pool_task_thread *>(
			thread::get_curr_thread());
	int thread_id = thread->get_pool_thread_id();

	const T *Amat = NULL;
	size_t num_Arows;
	size_t num_Acols;
	if (ins.size() == 2) {
		assert(ins[0]->store_layout() == matrix_layout_t::L_COL);
		Amat = (const T *) ins[0]->get_raw_arr();
		num_Arows = ins[0]->get_num_rows();
		num_Acols = ins[0]->get_num_cols();
	}
	detail::local_col_matrix_store::ptr Abuf;
	if (Amat == NULL) {
		size_t A_num_blocks = ins.size() - 1;
		size_t block_size = ins.front()->get_num_cols();
		size_t num_rows = ins.front()->get_num_rows();
		size_t num_cols = block_size * A_num_blocks;
		off_t global_start_row = ins.front()->get_global_start_row();
		off_t global_start_col = ins.front()->get_global_start_col();
		Abuf = fm::detail::local_col_matrix_store::ptr(
				new fm::detail::local_buf_col_matrix_store(
					global_start_row, global_start_col,
					num_rows, num_cols, get_scalar_type<T>(), -1));
		copy_from_blocks(ins.begin(), ins.begin() + ins.size() - 1, *Abuf);
		num_Arows = Abuf->get_num_rows();
		num_Acols = Abuf->get_num_cols();
		Amat = (const T *) Abuf->get_raw_arr();
	}
	assert(Amat);

	detail::local_matrix_store::const_ptr Bstore = ins.back();
	const T *Bmat = (const T *) Bstore->get_raw_arr();
	detail::local_matrix_store::ptr Bbuf;
	assert(Bstore->store_layout() == Blayout);
	if (Bmat == NULL) {
		if (Blayout == matrix_layout_t::L_COL)
			Bbuf = detail::local_matrix_store::ptr(
					new fm::detail::local_buf_col_matrix_store(0, 0,
						Bstore->get_num_rows(), Bstore->get_num_cols(),
						Bstore->get_type(), -1));
		else
			Bbuf = detail::local_matrix_store::ptr(
					new fm::detail::local_buf_row_matrix_store(0, 0,
						Bstore->get_num_rows(), Bstore->get_num_cols(),
						Bstore->get_type(), -1));
		Bbuf->copy_from(*Bstore);
		Bmat = (const T *) Bbuf->get_raw_arr();
	}
	assert(Bmat);

	if (res_bufs[thread_id] == NULL) {
		if (Blayout == matrix_layout_t::L_COL)
			const_cast<multiply_wide_op<T> *>(this)->res_bufs[thread_id]
				= detail::local_matrix_store::ptr(
						new fm::detail::local_buf_col_matrix_store(0, 0,
							out_num_rows, out_num_cols, get_scalar_type<T>(), -1));
		else
			const_cast<multiply_wide_op<T> *>(this)->res_bufs[thread_id]
				= detail::local_matrix_store::ptr(
						new fm::detail::local_buf_row_matrix_store(0, 0,
							out_num_rows, out_num_cols, get_scalar_type<T>(), -1));
		res_bufs[thread_id]->reset_data();
	}
	assert(res_bufs[thread_id]->store_layout() == Blayout);
	T *res_mat = (T *) res_bufs[thread_id]->get_raw_arr();
	// The A matrix is the transpose of the matrix we need. Since the A matrix
	// is stored in contiguous memory and is organized in row major, we can
	// easily interpret it as its transpose by switching its #rows and #cols.
	if (Blayout == matrix_layout_t::L_COL)
		cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, num_Acols,
				Bstore->get_num_cols(), num_Arows, 1, Amat, num_Acols, Bmat,
				Bstore->get_num_rows(), 1, res_mat, out_num_rows);
	else
		cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, num_Acols,
				Bstore->get_num_cols(), num_Arows, 1, Amat, num_Arows, Bmat,
				Bstore->get_num_cols(), 1, res_mat, out_num_cols);
}

dense_matrix::ptr MvTransMv_wide(
		const std::vector<detail::matrix_store::const_ptr> &blocks1,
		detail::matrix_store::const_ptr in2, const size_t MAX_MUL_BLOCKS)
{
	detail::mem_thread_pool::ptr threads
		= detail::mem_thread_pool::get_global_mem_threads();
	size_t num_threads = threads->get_num_threads();
	std::vector<detail::matrix_store::const_ptr> tmp_res;
	std::vector<detail::portion_mapply_op::const_ptr> mul_ops;
	for (size_t i = 0; i < blocks1.size(); i += MAX_MUL_BLOCKS) {
		size_t num = std::min(blocks1.size() - i, MAX_MUL_BLOCKS);
		std::vector<detail::matrix_store::const_ptr> tmp_ins(
				blocks1.begin() + i, blocks1.begin() + i + num);
		tmp_ins.push_back(in2);
		size_t out_num_rows = num * blocks1.front()->get_num_cols();
		size_t out_num_cols = in2->get_num_cols();
		detail::portion_mapply_op::const_ptr op(new multiply_wide_op<double>(
					num_threads, out_num_rows, out_num_cols,
					// The left matrix is larger. After it is transposed,
					// it is stored in row-major order. Let's multiply
					// the matrices in row-major order.
					in2->get_num_rows(), matrix_layout_t::L_ROW));
		tmp_res.push_back(__mapply_portion_virtual(tmp_ins, op,
					// This is a layout of the fake output. It doesn't
					// really matter.
					matrix_layout_t::L_COL));
		mul_ops.push_back(op);
	}

	set_caching(blocks1, false);
	set_caching(tmp_res, false);
	const_cast<detail::matrix_store &>(*in2).set_cache_portion(true);
	__mapply_portion(tmp_res, detail::portion_mapply_op::const_ptr(new empty_op()),
			matrix_layout_t::L_COL, false);
	set_caching(blocks1, true);

	// There we put all result matrices together.
	detail::matrix_store::ptr res = detail::matrix_store::create(
			blocks1.size() * blocks1.front()->get_num_cols(),
			in2->get_num_cols(), matrix_layout_t::L_ROW,
			tmp_res.front()->get_type(), -1, true);
	res->reset_data();
	size_t start_row = 0;
	for (size_t i = 0; i < mul_ops.size(); i++) {
		std::shared_ptr<const multiply_wide_op<double> > op
			= std::static_pointer_cast<const multiply_wide_op<double> >(mul_ops[i]);
		std::vector<detail::local_matrix_store::ptr> local_ms
			= op->get_partial_results();
		assert(local_ms.size() == num_threads);

		// Get a non-empty result generated in a worker thread.
		detail::local_matrix_store::ptr local_m;
		for (size_t i = 0; i < local_ms.size(); i++)
			if (local_ms[i])
				local_m = local_ms[i];
		assert(local_m);

		// Aggregate the results from omp threads.
		// This matrix is small. We can always keep it in memory.
		detail::local_matrix_store::ptr local_res = res->get_portion(
				start_row, 0, local_m->get_num_rows(), local_m->get_num_cols());
		const bulk_operate &add = local_res->get_type().get_basic_ops().get_add();
		for (size_t j = 0; j < local_ms.size(); j++) {
			// It's possible that the local matrix store doesn't exist
			// because the input matrix is very small.
			if (local_ms[j])
				detail::mapply2(*local_res, *local_ms[j], add, *local_res);
		}
		start_row += local_res->get_num_rows();
	}
	return dense_matrix::create(res);
}

dense_matrix::ptr block_multi_vector::MvTransMv(
		const block_multi_vector &mv) const
{
	detail::mem_col_matrix_store::ptr res = detail::mem_col_matrix_store::create(
			mv.get_num_cols(), this->get_num_cols(), type);
	std::vector<detail::matrix_store::const_ptr> blocks1(mv.get_num_blocks());
	for (size_t i = 0; i < blocks1.size(); i++) {
		blocks1[i] = mv.get_block(i)->get_raw_store();
		if (blocks1[i]->is_virtual())
			printf("materialize %s on the fly\n", blocks1[i]->get_name().c_str());
	}
	std::vector<detail::matrix_store::const_ptr> blocks2(this->get_num_blocks());
	for (size_t i = 0; i < blocks2.size(); i++) {
		blocks2[i] = this->get_block(i)->get_raw_store();
		if (blocks2[i]->is_virtual())
			printf("materialize %s on the fly\n", blocks2[i]->get_name().c_str());
	}

	if (blocks1.size() <= MAX_MUL_BLOCKS && blocks2.size() <= MAX_MUL_BLOCKS) {
		dense_matrix::ptr in1;
		if (blocks1.size() == 1)
			in1 = mv.get_block(0);
		else
			in1 = dense_matrix::create(collected_matrix_store::create(blocks1,
						mv.get_num_cols()));
		dense_matrix::ptr in2;
		if (blocks2.size() == 1)
			in2 = this->get_block(0);
		else
			in2 = dense_matrix::create(collected_matrix_store::create(blocks2,
						this->get_num_cols()));
		detail::matrix_stats_t orig_stats = detail::matrix_stats;
		dense_matrix::ptr ret = in1->transpose()->multiply(*in2,
				matrix_layout_t::L_NONE, true);
		detail::matrix_stats.print_diff(orig_stats);
		return ret;
	}
	else if (blocks1.size() > MAX_MUL_BLOCKS && blocks2.size() <= MAX_MUL_BLOCKS) {
		detail::matrix_store::const_ptr in2;
		assert(blocks1.front()->store_layout() == matrix_layout_t::L_COL);
		// We should convert the data layout in advance so that multiple groups
		// can share the data with the converted layout.
		if (blocks2.size() == 1) {
			dense_matrix::ptr tmp
				= this->get_block(0)->conv2(matrix_layout_t::L_ROW);
			in2 = tmp->get_raw_store();
		}
		else {
			dense_matrix::ptr tmp = dense_matrix::create(
					collected_matrix_store::create(blocks2, this->get_num_cols()));
			tmp = tmp->conv2(matrix_layout_t::L_ROW);
			in2 = tmp->get_raw_store();
		}
		detail::matrix_stats_t orig_stats = detail::matrix_stats;
		dense_matrix::ptr ret = MvTransMv_wide(blocks1, in2, MAX_MUL_BLOCKS);
		detail::matrix_stats.print_diff(orig_stats);
		return ret;
	}
	else if (blocks1.size() <= MAX_MUL_BLOCKS && blocks2.size() > MAX_MUL_BLOCKS) {
		detail::matrix_store::const_ptr in1;
		if (blocks1.size() == 1) {
			dense_matrix::ptr tmp
				= mv.get_block(0)->conv2(matrix_layout_t::L_ROW);
			in1 = tmp->get_raw_store();
		}
		else {
			dense_matrix::ptr tmp = dense_matrix::create(
					collected_matrix_store::create(blocks1, mv.get_num_cols()));
			tmp = tmp->conv2(matrix_layout_t::L_ROW);
			in1 = tmp->get_raw_store();
		}
		detail::matrix_stats_t orig_stats = detail::matrix_stats;
		dense_matrix::ptr ret = MvTransMv_wide(blocks2, in1, MAX_MUL_BLOCKS);
		detail::matrix_stats.print_diff(orig_stats);
		return ret->transpose();
	}
	else {
		// I don't need to handle this case right now.
		assert(0);
		return dense_matrix::ptr();
	}
}

std::vector<double> block_multi_vector::MvDot(const block_multi_vector &mv) const
{
	std::vector<double> ret;

	assert(mv.get_num_cols() == this->get_num_cols());
	assert(mv.get_block_size() == this->get_block_size());
	for (size_t i = 0; i < mv.get_num_blocks(); i++) {
		dense_matrix::ptr mat1 = get_block(i);
		dense_matrix::ptr mat2 = mv.get_block(i);
		dense_matrix::ptr res = mat1->multiply_ele(*mat2);
		vector::ptr sum = res->col_sum();
		std::vector<double> tmp = sum->conv2std<double>();
		ret.insert(ret.end(), tmp.begin(), tmp.end());
	}
	return ret;
}

typedef std::vector<int> block_col_set_t;

/*
 * This function splits the columns into blocks.
 */
std::vector<block_col_set_t> get_col_index_blocks(const std::vector<int>& index,
		size_t block_size)
{
	std::vector<block_col_set_t> col_blocks;
	col_blocks.push_back(block_col_set_t());
	for (size_t i = 0; i < index.size(); i++) {
		if (col_blocks.back().empty())
			col_blocks.back().push_back(index[i]);
		else if (col_blocks.back().back() / block_size == index[i] / block_size)
			col_blocks.back().push_back(index[i]);
		else {
			col_blocks.push_back(block_col_set_t());
			col_blocks.back().push_back(index[i]);
		}
	}
	return col_blocks;
}

fm::dense_matrix::ptr block_multi_vector::get_col_mat(
		const std::vector<off_t> &_index) const
{
	std::vector<int> index(_index.begin(), _index.end());
	std::vector<block_col_set_t> col_index_blocks
		= get_col_index_blocks(index, get_block_size());
	std::vector<dense_matrix::ptr> mats(col_index_blocks.size());
	for (size_t i = 0; i < mats.size(); i++) {
		block_col_set_t block_cols = col_index_blocks[i];
		int block_idx = block_cols[0] / get_block_size();
		std::vector<off_t> offs_in_block(block_cols.size());
		for (size_t i = 0; i < block_cols.size(); i++) {
			assert((size_t) block_idx == block_cols[i] / get_block_size());
			offs_in_block[i] = block_cols[i] % get_block_size();
		}
		mats[i] = get_block(block_idx)->get_cols(offs_in_block);
	}

	if (mats.size() == 1)
		return mats[0];
	else {
		std::vector<detail::matrix_store::const_ptr> const_mats(mats.size());
		for (size_t i = 0; i < mats.size(); i++)
			const_mats[i] = mats[i]->get_raw_store();
		return dense_matrix::create(collected_matrix_store::create(
					const_mats, index.size()));
	}
}

dense_matrix::ptr materialize_block(dense_matrix::ptr mat)
{
	assert(mat->is_virtual());
	num_col_writes += mat->get_num_cols();
	printf("set_block: materialize %s\n", mat->get_data().get_name().c_str());
	detail::matrix_stats_t orig_stats = detail::matrix_stats;
	mat->materialize_self();
	detail::matrix_stats.print_diff(orig_stats);
	return mat;
}

void block_multi_vector::set_block(const block_multi_vector &mv,
		const std::vector<int>& index)
{
	// We have to set the entire block.
	if (index[0] % get_block_size() == 0 && index.size() % get_block_size() == 0) {
		size_t num_blocks = index.size() / get_block_size();
		size_t block_start = index[0] / get_block_size();
		for (size_t i = 0; i < num_blocks; i++)
			this->set_block(block_start + i, mv.get_block(i));
	}
	else {
		assert(std::is_sorted(index.begin(), index.end()));
		std::vector<block_col_set_t> col_index_blocks
			= get_col_index_blocks(index, get_block_size());
		size_t mv_col_idx = 0;
		for (size_t k = 0; k < col_index_blocks.size(); k++) {
			// The column index on the current MV.
			block_col_set_t block_cols = col_index_blocks[k];
			// The block index on the current MV.
			int this_block_idx = block_cols[0] / get_block_size();
			// The column index inside a block of the current MV.
			std::vector<off_t> offs_in_block(block_cols.size());
			for (size_t i = 0; i < block_cols.size(); i++) {
				assert((size_t) this_block_idx == block_cols[i] / get_block_size());
				offs_in_block[i] = block_cols[i] % get_block_size();
			}
			std::vector<off_t> mv_block_cols(block_cols.size());
			for (size_t i = 0; i < block_cols.size(); i++)
				mv_block_cols[i] = mv_col_idx++;

			// Get all columns in the block.
			fm::dense_matrix::ptr this_block = get_block(this_block_idx);
			fm::dense_matrix::ptr new_mat;
			// We can use the entire block in the input MV to replace the block
			// in the current MV.
			if (offs_in_block.front() == 0
					&& offs_in_block.size() == get_block_size())
				new_mat = mv.get_col_mat(mv_block_cols);
			// The front part of the block in the current MV is replaced.
			else if (offs_in_block.front() == 0) {
				std::vector<detail::matrix_store::const_ptr> mats(2);
				mats[0] = mv.get_col_mat(mv_block_cols)->get_raw_store();
				std::vector<off_t> remain_cols(
						get_block_size() - mv_block_cols.size());
				for (size_t i = 0; i < remain_cols.size(); i++)
					remain_cols[i] = offs_in_block.back() + i + 1;
				assert((size_t) remain_cols.back() == get_block_size() - 1);
				mats[1] = this_block->get_cols(remain_cols)->get_raw_store();
				if (mats[1] == NULL)
					mats[1] = materialize_block(this_block)->get_cols(
							remain_cols)->get_raw_store();
				new_mat = fm::dense_matrix::create(
						collected_matrix_store::create(mats, get_block_size()));
			}
			else if ((size_t) offs_in_block.back() == get_block_size() - 1) {
				std::vector<detail::matrix_store::const_ptr> mats(2);
				std::vector<off_t> remain_cols(
						get_block_size() - mv_block_cols.size());
				for (size_t i = 0; i < remain_cols.size(); i++)
					remain_cols[i] = i;
				assert(remain_cols.back() == offs_in_block.front() - 1);
				mats[0] = this_block->get_cols(remain_cols)->get_raw_store();
				if (mats[0] == NULL)
					mats[0] = materialize_block(this_block)->get_cols(
							remain_cols)->get_raw_store();
				mats[1] = mv.get_col_mat(mv_block_cols)->get_raw_store();
				new_mat = fm::dense_matrix::create(
						collected_matrix_store::create(mats, get_block_size()));
			}
			else {
				std::vector<detail::matrix_store::const_ptr> mats(3);
				std::vector<off_t> front_cols(offs_in_block.front());
				for (size_t i = 0; i < front_cols.size(); i++)
					front_cols[i] = i;
				assert(front_cols.back() == offs_in_block.front() - 1);
				mats[0] = this_block->get_cols(front_cols)->get_raw_store();
				if (mats[0] == NULL)
					mats[0] = materialize_block(this_block)->get_cols(
							front_cols)->get_raw_store();
				mats[1] = mv.get_col_mat(mv_block_cols)->get_raw_store();
				std::vector<off_t> back_cols(
						get_block_size() - offs_in_block.back() - 1);
				for (size_t i = 0; i < back_cols.size(); i++)
					back_cols[i] = offs_in_block.back() + i + 1;
				assert((size_t) back_cols.back() == get_block_size() - 1);
				mats[2] = this_block->get_cols(back_cols)->get_raw_store();
				assert(mats[2]);
				assert(front_cols.size() + mv_block_cols.size() + back_cols.size()
						== get_block_size());
				new_mat = fm::dense_matrix::create(
						collected_matrix_store::create(mats, get_block_size()));
			}

			set_block(this_block_idx, new_mat);
		}
	}
}

fm::dense_matrix::ptr block_multi_vector::conv2matrix() const
{
	if (mats.size() == 1)
		return mats[0];
	else {
		std::vector<detail::matrix_store::const_ptr> const_mats(mats.size());
		for (size_t i = 0; i < mats.size(); i++)
			const_mats[i] = mats[i]->get_raw_store();
		size_t tot_num_cols = 0;
		for (size_t i = 0; i < mats.size(); i++)
			tot_num_cols += mats[i]->get_num_cols();
		dense_matrix::ptr new_mat
			= dense_matrix::create(collected_matrix_store::create(
						const_mats, tot_num_cols));
		printf("conv2matrix: materialize %s\n",
				new_mat->get_data().get_name().c_str());
		num_col_writes += new_mat->get_num_cols();
		new_mat->materialize_self();
		return new_mat;
	}
}

}

}
