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

using namespace fm;
size_t num_col_writes = 0;
size_t num_col_writes_concept = 0;
size_t num_col_reads_concept = 0;

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
			fm::dense_matrix::ptr mirrored_mat): block_multi_vector(mats) {
		this->mirrored_mat = mirrored_mat;
		this->offs_in_mirror = offs_in_mirror;
	}
public:
	static ptr create(fm::dense_matrix::ptr mat, const std::vector<off_t> &offs,
			fm::dense_matrix::ptr mirrored_mat) {
		std::vector<fm::dense_matrix::ptr> mats(1);
		mats[0] = mat;
		return ptr(new mirror_cols_block_multi_vector(mats, offs, mirrored_mat));
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
		mats[i] = dense_matrix::create(nrow, block_size,
				matrix_layout_t::L_COL, type, get_num_nodes(), true);
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
		const detail::dotp_matrix_store *dotp
			= dynamic_cast<const detail::dotp_matrix_store *>(
					block->get_raw_store().get());
		if (dotp)
			ret = dense_matrix::create(dotp->get_cols(offs));
		else {
			num_col_writes += block->get_num_cols();
			printf("materialize %s\n", block->get_data().get_name().c_str());
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
				index.size(), get_block_size(), get_type());
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
				index.size(), index.size(), get_type());
		dense_matrix::ptr block = get_block(block_start);
		dense_matrix::ptr ret1;
		if (block->is_virtual()) {
			// We need to handle the special case explicitly.
			const detail::dotp_matrix_store *dotp
				= dynamic_cast<const detail::dotp_matrix_store *>(
					block->get_raw_store().get());
			if (dotp)
				ret1 = dense_matrix::create(dotp->get_cols(local_offs));
			else {
				num_col_writes += block->get_num_cols();
				printf("materialize %s\n", block->get_data().get_name().c_str());
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
				index.size(), 1, get_type());
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
			= mirror_block_multi_vector::create(mats);
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
					block->get_cols(idxs_in_block), idxs_in_block, block);
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
	for (size_t i = 0; i < num_blocks; i++) {
		dense_matrix::ptr in = X.get_block(i);
		dense_matrix::ptr res;
		if (in->store_layout() == matrix_layout_t::L_COL)
			in = in->conv2(matrix_layout_t::L_ROW);
		if (in->is_virtual()) {
			num_col_writes += in->get_num_cols();
			printf("materialize %s\n", in->get_data().get_name().c_str());
			detail::matrix_stats_t orig_stats = detail::matrix_stats;
			in->materialize_self();
			detail::matrix_stats.print_diff(orig_stats);
		}
		res = multiply.run(in);
		if (res->store_layout() == matrix_layout_t::L_ROW)
			res = res->conv2(matrix_layout_t::L_COL);
		assert(res->store_layout() == matrix_layout_t::L_COL);
		Y.set_block(i, res);
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
fm::detail::portion_mapply_op::const_ptr gemm_op<T>::transpose() const
{
	return fm::detail::portion_mapply_op::const_ptr(new t_gemm_op<T>(*this));
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
	if (A.get_num_blocks() > 2) {
		num_col_writes += vecs->mats[0]->get_num_cols();
		printf("materialize %s\n", vecs->mats[0]->get_data().get_name().c_str());
		detail::matrix_stats_t orig_stats = detail::matrix_stats;
		vecs->mats[0]->materialize_self();
		detail::matrix_stats.print_diff(orig_stats);
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
			get_num_rows(), get_num_cols(), block_size, type);
	size_t num_blocks = get_num_blocks();
	for (size_t i = 0; i < num_blocks; i++)
		ret->mats[i] = this->get_block(i)->add(*vecs.get_block(i));
	return ret;
}

dense_matrix::ptr block_multi_vector::MvTransMv(
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
			if (mv_block->is_virtual())
				printf("materialize %s on the fly\n",
						mv_block->get_data().get_name().c_str());
			dense_matrix::const_ptr block = get_block(i);
			if (block->is_virtual())
				printf("materialize %s on the fly\n",
						block->get_data().get_name().c_str());
			detail::matrix_stats_t orig_stats = detail::matrix_stats;
			fm::dense_matrix::ptr tA = mv.get_block(j)->transpose();
			fm::dense_matrix::ptr res1 = tA->multiply(*get_block(i),
					// I should use BLAS for multiplication here.
					matrix_layout_t::L_ROW, true);
			detail::matrix_stats.print_diff(orig_stats);
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
	return dense_matrix::create(res);
}

typedef std::vector<std::pair<int, detail::NUMA_vec_store::ptr> > block_col_set_t;

std::vector<block_col_set_t> get_col_index_blocks(const block_multi_vector &mv,
		const std::vector<int>& index, size_t block_size)
{
	std::vector<block_col_set_t> col_blocks;

	std::set<int> block_set;
	for (size_t i = 0; i < index.size(); i++) {
		int col_idx = index[i];
		size_t block_idx = col_idx / block_size;
		block_set.insert(block_idx);

		dense_matrix::const_ptr tmp_mat = mv.get_col(i);
		detail::vec_store::ptr col = const_cast<detail::NUMA_col_tall_matrix_store &>(
					dynamic_cast<const detail::NUMA_col_tall_matrix_store &>(
						tmp_mat->get_data())).get_col_vec(0);
		// Not in the same block
		if (col_blocks.empty() || col_blocks.back()[0].first / block_size
				!= block_idx) {
			block_col_set_t set(1);
			set[0] = std::pair<int, detail::NUMA_vec_store::ptr>(col_idx,
					detail::NUMA_vec_store::cast(col));
			col_blocks.push_back(set);
		}
		else
			// In the same block.
			col_blocks.back().push_back(std::pair<int, detail::NUMA_vec_store::ptr>(
						col_idx, detail::NUMA_vec_store::cast(col)));
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
			// The block might be a one-value matrix, so we have to materialize
			// the matrix to change its columns. This rarely happens.
			// TODO We might want the uninitialized columns to be virtualized.
			if (block->is_virtual()) {
				num_col_writes += block->get_num_cols();
				printf("materialize %s\n", block->get_data().get_name().c_str());
				detail::matrix_stats_t orig_stats = detail::matrix_stats;
				block->materialize_self();
				detail::matrix_stats.print_diff(orig_stats);
			}
			detail::NUMA_col_tall_matrix_store &numa_mat
				= const_cast<detail::NUMA_col_tall_matrix_store &>(
						dynamic_cast<const detail::NUMA_col_tall_matrix_store &>(
							block->get_data()));
			std::vector<detail::NUMA_vec_store::ptr> cols(numa_mat.get_num_cols());
			for (size_t i = 0; i < cols.size(); i++)
				cols[i] = detail::NUMA_vec_store::cast(numa_mat.get_col_vec(i));

			// changes some of the columns accordingly.
			for (size_t i = 0; i < block_cols.size(); i++)
				cols[offs_in_block[i]] = block_cols[i].second;

			// Create a dense matrix for the block and assign it to the original
			// block.
			fm::dense_matrix::ptr new_mat = fm::dense_matrix::create(
					fm::detail::NUMA_col_tall_matrix_store::create(cols));
			set_block(block_idx, new_mat);
		}
	}
}

class merge_cols_op: public detail::portion_mapply_op
{
public:
	merge_cols_op(size_t num_rows, size_t num_cols,
			const scalar_type &type): detail::portion_mapply_op(
				num_rows, num_cols, type) {
	}

	virtual void run(const std::vector<detail::local_matrix_store::const_ptr> &ins,
			detail::local_matrix_store &out) const;

	// This is only used once, so it's not necessary to implement these
	// two methods.
	virtual portion_mapply_op::const_ptr transpose() const {
		assert(0);
		return portion_mapply_op::const_ptr();
	}
	virtual std::string to_string(
			const std::vector<detail::matrix_store::const_ptr> &mats) const {
		return std::string();
	}
};

void merge_cols_op::run(
		const std::vector<detail::local_matrix_store::const_ptr> &ins,
		detail::local_matrix_store &out) const
{
	assert(out.store_layout() == matrix_layout_t::L_COL);
	detail::local_col_matrix_store &col_out
		= static_cast<detail::local_col_matrix_store &>(out);

	size_t out_col_idx = 0;
	for (size_t i = 0; i < ins.size(); i++) {
		assert(ins[i]->store_layout() == matrix_layout_t::L_COL);
		const detail::local_col_matrix_store &col_in
			= static_cast<const detail::local_col_matrix_store &>(*ins[i]);
		for (size_t col_idx = 0; col_idx < col_in.get_num_cols(); col_idx++) {
			assert(out_col_idx < out.get_num_cols());
			memcpy(col_out.get_col(out_col_idx), col_in.get_col(col_idx),
					col_in.get_num_rows() * col_in.get_entry_size());
			out_col_idx++;
		}
	}
	assert(out_col_idx == out.get_num_cols());
}

fm::dense_matrix::ptr block_multi_vector::conv2matrix() const
{
	if (mats.size() == 1)
		return mats[0];
	else {
		std::vector<dense_matrix::const_ptr> const_mats(mats.begin(),
				mats.end());
		size_t tot_num_cols = 0;
		for (size_t i = 0; i < mats.size(); i++)
			tot_num_cols += mats[i]->get_num_cols();
		dense_matrix::ptr new_mat = fm::detail::mapply_portion(const_mats,
				detail::portion_mapply_op::const_ptr(new merge_cols_op(
						mats[0]->get_num_rows(), tot_num_cols,
						mats[0]->get_type())),
				mats[0]->store_layout());
		printf("conv2matrix: materialize %s\n",
				new_mat->get_data().get_name().c_str());
		num_col_writes += new_mat->get_num_cols();
		new_mat->materialize_self();
		return new_mat;
	}
}
