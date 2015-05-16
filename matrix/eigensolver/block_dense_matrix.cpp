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

using namespace fm;
size_t num_col_writes = 0;

class mirror_block_multi_vector: public block_multi_vector
{
	// When mirroring the original matrix, we need to make sure that
	// the pointers point to the original matrix.
	mirror_block_multi_vector(
			const std::vector<fm::dense_matrix::ptr> &mats): block_multi_vector(mats) {
	}
public:
	virtual void set_block(off_t block_idx, fm::dense_matrix::const_ptr mat) {
		assert(mat->get_num_cols() == get_block_size());
		mats[block_idx]->assign(*mat);
	}

	static ptr create(const std::vector<fm::dense_matrix::ptr> &mats) {
		return ptr(new mirror_block_multi_vector(mats));
	}
};

/*
 * This mirrors some columns in the block of a block_multi_vector.
 * It only works if each block use NUMA_matrix_store.
 */
class mirror_cols_block_multi_vector: public block_multi_vector
{
	fm::dense_matrix::ptr mirrored_mat;
	std::vector<off_t> offs_in_mirror;

	mirror_cols_block_multi_vector(
			const std::vector<fm::dense_matrix::ptr> &mats,
			const std::vector<off_t> &offs_in_mirror,
			fm::dense_matrix::ptr mirrored_mat): block_multi_vector(mats) {
		this->mirrored_mat = mirrored_mat;
		this->offs_in_mirror = offs_in_mirror;
	}
public:
	static ptr create(const std::vector<NUMA_vector::ptr> &vecs,
			const std::vector<off_t> &offs_in_mirror,
			fm::dense_matrix::ptr mirrored_mat) {
		std::vector<fm::dense_matrix::ptr> mats(1);
		mats[0] = fm::mem_dense_matrix::create(
				detail::NUMA_col_tall_matrix_store::create(vecs));
		return ptr(new mirror_cols_block_multi_vector(mats, offs_in_mirror,
					mirrored_mat));
	}

	virtual void set_block(off_t block_idx, fm::dense_matrix::const_ptr mat) {
		assert(block_idx == 0);
		assert(mat->get_num_cols() == get_block_size());
		mats[block_idx]->assign(*mat);
		if (mirrored_mat->is_virtual()) {
			printf("materialize %ld cols\n", mirrored_mat->get_num_cols());
			num_col_writes += mirrored_mat->get_num_cols();
		}
		if (mat->is_virtual()) {
			printf("materialize %ld cols\n", mat->get_num_cols());
			num_col_writes += mat->get_num_cols();
		}
		mirrored_mat->materialize_self();
		mat->materialize_self();

		detail::NUMA_col_tall_matrix_store &mirrored_numa_mat
			= const_cast<detail::NUMA_col_tall_matrix_store &>(
					dynamic_cast<const detail::NUMA_col_tall_matrix_store &>(
						mirrored_mat->get_data()));
		std::vector<NUMA_vector::ptr> cols(mirrored_mat->get_num_cols());
		for (size_t i = 0; i < cols.size(); i++)
			cols[i] = mirrored_numa_mat.get_col_vec(i);

		detail::NUMA_col_tall_matrix_store &numa_in_mat
			= const_cast<detail::NUMA_col_tall_matrix_store &>(
					dynamic_cast<const detail::NUMA_col_tall_matrix_store &>(
						mat->get_data()));
		for (size_t i = 0; i < mat->get_num_cols(); i++)
			cols[offs_in_mirror[i]] = numa_in_mat.get_col_vec(i);

		fm::dense_matrix::ptr new_mat = fm::mem_dense_matrix::create(
				fm::detail::NUMA_col_tall_matrix_store::create(cols));
		mirrored_mat->assign(*new_mat);
	}
};

block_multi_vector::block_multi_vector(
		const std::vector<fm::dense_matrix::ptr> &mats): type(mats[0]->get_type())
{
	this->num_rows = mats[0]->get_num_rows();
	this->num_cols = mats[0]->get_num_cols() * mats.size();
	this->block_size = mats[0]->get_num_cols();
	this->mats = mats;
	for (size_t i = 1; i < mats.size(); i++)
		assert(this->block_size == mats[i]->get_num_cols());
}

block_multi_vector::block_multi_vector(size_t nrow, size_t ncol,
		size_t block_size, const fm::scalar_type &_type): type(_type)
{
	this->num_rows = nrow;
	this->num_cols = ncol;
	this->block_size = block_size;
	mats.resize(ncol / block_size);
	for (size_t i = 0; i < mats.size(); i++)
		mats[i] = mem_dense_matrix::create(nrow, block_size,
				matrix_layout_t::L_COL, type, get_num_nodes());
}

dense_matrix::const_ptr block_multi_vector::get_col(off_t col_idx) const
{
	off_t block_idx = col_idx / block_size;
	off_t local_col_idx = col_idx % block_size;
	std::vector<off_t> offs(1);
	offs[0] = local_col_idx;
	fm::dense_matrix::const_ptr block = get_block(block_idx);
	if (block->is_virtual()) {
		printf("materialize %ld cols\n", block->get_num_cols());
		num_col_writes += block->get_num_cols();
	}
	block->materialize_self();
	return block->get_cols(offs);
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
				index.size(), get_block_size(), get_type());
		size_t num_blocks = index.size() / get_block_size();
		size_t block_start = index[0] / get_block_size();
		for (size_t i = 0; i < num_blocks; i++)
			ret->set_block(i, get_block(i + block_start));
		printf("get block [%ld-%ld]\n", block_start, block_start + num_blocks - 1);
		return ret;
	}
	else if (is_same_block(index, block_size)) {
		// We can only fetch columns from the same block.
		size_t block_start = index[0] / get_block_size();
		std::vector<off_t> local_offs(index.size());
		for (size_t i = 0; i < local_offs.size(); i++)
			local_offs[i] = index[i] - block_start * get_block_size();

		block_multi_vector::ptr ret = block_multi_vector::create(get_num_rows(),
				index.size(), index.size(), get_type());
		dense_matrix::ptr block = get_block(block_start);
		if (block->is_virtual()) {
			printf("materialize %ld cols\n", block->get_num_cols());
			num_col_writes += block->get_num_cols();
		}
		block->materialize_self();
		ret->set_block(0, block->get_cols(local_offs));
		printf("get block [%ld]\n", block_start);
		return ret;
	}
	else {
		block_multi_vector::ptr ret = block_multi_vector::create(get_num_rows(),
				index.size(), 1, get_type());
		for (size_t i = 0; i < index.size(); i++)
			ret->set_block(i, get_col(index[i]));
		printf("get column [%d-%d]\n", index.front(), index.back());
		return ret;
	}
}

bool block_multi_vector::resize_block(size_t new_block_size)
{
	assert(0);
	for (size_t i = 0; i < mats.size(); i++)
		assert(mats[i].use_count() == 1);
	assert(get_block_size() % new_block_size == 0);
	size_t num_subblocks = get_block_size() / new_block_size;
	std::vector<fm::dense_matrix::ptr> new_mats(
			num_subblocks * get_num_blocks());
	for (size_t i = 0; i < new_mats.size(); i++) {
		size_t subblock_idx = i % num_subblocks;
		size_t old_block_idx = i / num_subblocks;
		std::vector<off_t> subblock_offs(new_block_size);
		for (size_t j = 0; j < new_block_size; j++)
			subblock_offs[j] = subblock_idx * new_block_size + j;
		const fm::detail::mem_matrix_store &mem_store
			= dynamic_cast<const fm::detail::mem_matrix_store &>(
					mats[old_block_idx]->get_data());
		fm::detail::mem_matrix_store::const_ptr sub_store
			= fm::detail::mem_matrix_store::cast(mem_store.get_cols(
						subblock_offs));
		new_mats[i] = fm::mem_dense_matrix::create(sub_store);
	}

	mats = new_mats;
	block_size = new_block_size;
	return true;
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
			= mirror_block_multi_vector::create(mats);
		printf("mirror block [%ld-%ld]\n", block_start, block_start + num_blocks - 1);
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
		assert(!block->is_virtual());
		detail::NUMA_col_tall_matrix_store &numa_block
			= const_cast<detail::NUMA_col_tall_matrix_store &>(
					dynamic_cast<const detail::NUMA_col_tall_matrix_store &>(
						block->get_data()));
		std::vector<NUMA_vector::ptr> cols(index.size());
		for (size_t i = 0; i < cols.size(); i++)
			cols[i] = numa_block.get_col_vec(idxs_in_block[i]);

		mirror_cols_block_multi_vector::ptr ret
			= mirror_cols_block_multi_vector::create(cols, idxs_in_block, block);
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

void block_multi_vector::sparse_matrix_multiply(const sp_multiply &multiply,
		const sparse_matrix &A, const block_multi_vector &X,
		block_multi_vector &Y)
{
	assert(A.get_num_cols() == X.get_num_rows());
	assert(A.get_num_rows() == Y.get_num_rows());
	assert(Y.get_num_cols() == X.get_num_cols());
	assert(Y.get_num_blocks() == X.get_num_blocks());
	assert(X.get_block_size() == Y.get_block_size());
	assert(X.get_type() == Y.get_type());
	size_t num_blocks = X.get_num_blocks();
	for (size_t i = 0; i < num_blocks; i++) {
		dense_matrix::const_ptr block = X.get_block(i);
		if (block->is_virtual()) {
			printf("materialize %ld cols\n", block->get_num_cols());
			num_col_writes += block->get_num_cols();
		}
		block->materialize_self();
		const detail::mem_matrix_store &in
			= dynamic_cast<const detail::mem_matrix_store &>(X.get_block(i)->get_data());
		detail::mem_matrix_store::ptr res = detail::mem_matrix_store::create(
				Y.get_num_rows(), X.get_block_size(), matrix_layout_t::L_COL,
				X.get_type(), in.get_num_nodes());
		if (in.store_layout() == matrix_layout_t::L_COL) {
			detail::matrix_store::const_ptr tmp
				= in.conv2(matrix_layout_t::L_ROW);
			multiply.run(A, *tmp, *res);
		}
		else
			multiply.run(A, in, *res);
		Y.set_block(i, mem_dense_matrix::create(res));
	}
}

block_multi_vector::ptr block_multi_vector::clone() const
{
	block_multi_vector::ptr vecs= block_multi_vector::create(
			get_num_rows(), get_num_cols(), block_size, type);
	size_t num_blocks = get_num_blocks();
	for (size_t i = 0; i < num_blocks; i++)
		vecs->mats[i] = this->get_block(i)->clone();
	return vecs;
}

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
};

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
		assert(res_store.get_node_id() == col_ins[i]->get_node_id());
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
	assert(A_num_blocks + C_num_blocks == ins.size());
	int node_id = ins.front()->get_node_id();
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
					get_scalar_type<T>(), node_id);
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
					num_rows, num_cols, get_scalar_type<T>(), node_id));
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
					num_rows, num_cols, get_scalar_type<T>(), node_id));
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

block_multi_vector::ptr block_multi_vector::gemm(const block_multi_vector &A,
		detail::mem_col_matrix_store::const_ptr B, const scalar_variable &alpha,
		const scalar_variable &beta) const
{
	// The product can only have one block.
	block_multi_vector::ptr vecs= block_multi_vector::create(
			// mapply_portion can only output one matrix, so `vecs' can
			// only have one block.
			get_num_rows(), B->get_num_cols(), B->get_num_cols(), type);
	assert(A.get_num_rows() == this->get_num_rows());

	double d_alpha
		= dynamic_cast<const scalar_variable_impl<double> &>(alpha).get();
	double d_beta
		= dynamic_cast<const scalar_variable_impl<double> &>(beta).get();

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

	// This works for double float points.
	assert(type == fm::get_scalar_type<double>());
	size_t A_num_blocks = A.get_num_blocks();
	size_t C_num_blocks = d_beta != 0 ? this->get_num_blocks() : 0;
	assert(A_num_blocks + C_num_blocks == mats.size());
	// I assume B is small enough.
	detail::portion_mapply_op::const_ptr op(new gemm_op<double>(B,
				A_num_blocks, C_num_blocks, d_alpha, d_beta,
				this->get_num_rows(), this->get_num_cols()));
	vecs->mats[0] = mapply_portion(mats, op, matrix_layout_t::L_COL);
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
			get_num_rows(), get_num_cols(), block_size, type);
	size_t num_blocks = get_num_blocks();
	for (size_t i = 0; i < num_blocks; i++)
		ret->mats[i] = this->get_block(i)->add(*vecs.get_block(i));
	return ret;
}

mem_dense_matrix::ptr block_multi_vector::MvTransMv(
		const block_multi_vector &mv) const
{
	// TODO this isn't an efficient implementation.
	detail::mem_row_matrix_store::ptr res = detail::mem_row_matrix_store::create(
			mv.get_num_cols(), this->get_num_cols(), type);
	size_t block_num_rows = mv.get_block_size();
	size_t block_num_cols = get_block_size();
	for (size_t i = 0; i < get_num_blocks(); i++) {
		for (size_t j = 0; j < mv.get_num_blocks(); j++) {
			dense_matrix::const_ptr mv_block = mv.get_block(j);
			if (mv_block->is_virtual()) {
				printf("materialize %ld cols\n", mv_block->get_num_cols());
				num_col_writes += mv_block->get_num_cols();
			}
			mv_block->materialize_self();
			dense_matrix::const_ptr block = get_block(i);
			if (block->is_virtual()) {
				printf("materialize %ld cols\n", block->get_num_cols());
				num_col_writes += block->get_num_cols();
			}
			block->materialize_self();
			fm::mem_dense_matrix::ptr tA = fm::mem_dense_matrix::cast(
					mv.get_block(j)->transpose());
			fm::mem_dense_matrix::ptr res1 = fm::mem_dense_matrix::cast(
					tA->multiply(*get_block(i), matrix_layout_t::L_ROW));
			assert(res->store_layout() == res1->store_layout());
			detail::local_matrix_store::ptr part = res->get_portion(
					block_num_rows * j, block_num_cols * i,
					block_num_rows, block_num_cols);
			detail::local_matrix_store::const_ptr local_store
				= static_cast<const detail::mem_matrix_store &>(
						res1->get_data()).get_portion(0);
			part->copy_from(*local_store);
		}
	}
	return mem_dense_matrix::create(res);
}

typedef std::vector<std::pair<int, NUMA_vector::ptr> > block_col_set_t;

std::vector<block_col_set_t> get_col_index_blocks(const block_multi_vector &mv,
		const std::vector<int>& index, size_t block_size)
{
	std::vector<block_col_set_t> col_blocks;

	std::set<int> block_set;
	for (size_t i = 0; i < index.size(); i++) {
		int col_idx = index[i];
		size_t block_idx = col_idx / block_size;
		block_set.insert(block_idx);

		NUMA_vector::ptr col = const_cast<detail::NUMA_col_tall_matrix_store &>(
					dynamic_cast<const detail::NUMA_col_tall_matrix_store &>(
						mv.get_col(i)->get_data())).get_col_vec(0);
		// Not in the same block
		if (col_blocks.empty() || col_blocks.back()[0].first / block_size
				!= block_idx) {
			block_col_set_t set(1);
			set[0] = std::pair<int, NUMA_vector::ptr>(col_idx, col);
			col_blocks.push_back(set);
		}
		else
			// In the same block.
			col_blocks.back().push_back(std::pair<int, NUMA_vector::ptr>(
						col_idx, col));
	}
	assert(block_set.size() == col_blocks.size());

	return col_blocks;
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
		// They should all be in the same block.
		std::vector<block_col_set_t> col_index_blocks
			= get_col_index_blocks(mv, index, get_block_size());
		for (size_t k = 0; k < col_index_blocks.size(); k++) {
			block_col_set_t block_cols = col_index_blocks[k];
			int block_idx = block_cols[0].first / get_block_size();
			std::vector<off_t> offs_in_block(block_cols.size());
			for (size_t i = 0; i < block_cols.size(); i++) {
				assert((size_t) block_idx == block_cols[i].first / get_block_size());
				offs_in_block[i] = block_cols[i].first % get_block_size();
			}

			// Get all columns in the block.
			fm::dense_matrix::ptr block = get_block(block_idx);
			detail::NUMA_col_tall_matrix_store &numa_mat
				= const_cast<detail::NUMA_col_tall_matrix_store &>(
						dynamic_cast<const detail::NUMA_col_tall_matrix_store &>(
							block->get_data()));
			std::vector<NUMA_vector::ptr> cols(numa_mat.get_num_cols());
			for (size_t i = 0; i < cols.size(); i++)
				cols[i] = numa_mat.get_col_vec(i);

			// changes some of the columns accordingly.
			for (size_t i = 0; i < block_cols.size(); i++)
				cols[offs_in_block[i]] = block_cols[i].second;

			// Create a dense matrix for the block and assign it to the original
			// block.
			fm::dense_matrix::ptr new_mat = fm::mem_dense_matrix::create(
					fm::detail::NUMA_col_tall_matrix_store::create(cols));
			set_block(block_idx, new_mat);
		}
	}
}
