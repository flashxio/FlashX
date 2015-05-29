/*
 * Copyright 2015 Open Connectome Project (http://openconnecto.me)
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

#include "mapply_matrix_store.h"
#include "local_matrix_store.h"
#include "vec_store.h"

namespace fm
{

namespace detail
{

namespace
{

size_t SUB_CHUNK_SIZE = 4 * 1024;

class mapply_store
{
	std::vector<local_matrix_store::const_ptr> ins;
	const portion_mapply_op &op;
	raw_data_array buf;
	bool materialized;
	local_matrix_store *lstore;
public:
	mapply_store(const std::vector<local_matrix_store::const_ptr> &ins,
			const portion_mapply_op &_op, local_matrix_store *lstore): op(_op) {
		this->lstore = lstore;
		this->ins = ins;
		materialized = false;
	}

	const char *get_raw_arr() const {
		return buf.get_raw();
	}

	bool is_materialized() const {
		return materialized;
	}
	void materialize() const;
	void materialize_whole();

	void resize(off_t local_start_row, off_t local_start_col,
			size_t local_num_rows, size_t local_num_cols);
	void reset_size();
};

void mapply_store::reset_size()
{
	// When the matrix store is resized, we should drop the materialized data.
	materialized = false;
	buf = raw_data_array();
	for (size_t i = 0; i < ins.size(); i++)
		const_cast<local_matrix_store *>(ins[i].get())->reset_size();
}

void mapply_store::resize(off_t local_start_row, off_t local_start_col,
		size_t local_num_rows, size_t local_num_cols)
{
	// When the matrix store is resized, we should drop the materialized data.
	materialized = false;
	// To avoid unnecessary memory allocation, if the new exposed part isn't
	// larger than the buffer size, we can keep the buffer.
	if (local_num_rows * local_num_cols * lstore->get_entry_size()
			> buf.get_num_bytes())
		buf = raw_data_array();

	if (lstore->is_wide()) {
		assert(local_start_row == 0 && local_num_rows == lstore->get_num_rows());
		for (size_t i = 0; i < ins.size(); i++)
			const_cast<local_matrix_store *>(ins[i].get())->resize(0,
					local_start_col, ins[i]->get_num_rows(), local_num_cols);
	}
	else {
		assert(local_start_col == 0 && local_num_cols == lstore->get_num_cols());
		for (size_t i = 0; i < ins.size(); i++)
			const_cast<local_matrix_store *>(ins[i].get())->resize(
					local_start_row, 0, local_num_rows, ins[i]->get_num_cols());
	}
}

void mapply_store::materialize_whole()
{
	local_matrix_store::ptr res;
	if (lstore->store_layout() == matrix_layout_t::L_COL) {
		local_buf_col_matrix_store *tmp = new local_buf_col_matrix_store(
				lstore->get_global_start_row(),
				lstore->get_global_start_col(), lstore->get_num_rows(),
				lstore->get_num_cols(), lstore->get_type(), -1);
		res = local_matrix_store::ptr(tmp);
		this->buf = tmp->get_data();
	}
	else {
		local_buf_row_matrix_store *tmp = new local_buf_row_matrix_store(
				lstore->get_global_start_row(),
				lstore->get_global_start_col(), lstore->get_num_rows(),
				lstore->get_num_cols(), lstore->get_type(), -1);
		res = local_matrix_store::ptr(tmp);
		this->buf = tmp->get_data();
	}

	std::vector<local_matrix_store *> mutable_ins(ins.size());
	for (size_t i = 0; i < ins.size(); i++)
		mutable_ins[i] = const_cast<local_matrix_store *>(ins[i].get());

	if (lstore->is_wide()) {
		for (size_t local_start_col = 0; local_start_col < lstore->get_num_cols();
				local_start_col += SUB_CHUNK_SIZE) {
			size_t local_num_cols = std::min(SUB_CHUNK_SIZE,
					lstore->get_num_cols() - local_start_col);
			for (size_t i = 0; i < mutable_ins.size(); i++) {
				mutable_ins[i]->resize(0, local_start_col,
						mutable_ins[i]->get_num_rows(), local_num_cols);
			}
			res->resize(0, local_start_col, res->get_num_rows(), local_num_cols);
			op.run(ins, *res);
		}
	}
	else {
		// If this is a tall matrix.
		for (size_t local_start_row = 0; local_start_row < lstore->get_num_rows();
				local_start_row += SUB_CHUNK_SIZE) {
			size_t local_num_rows = std::min(SUB_CHUNK_SIZE,
					lstore->get_num_rows() - local_start_row);
			for (size_t i = 0; i < mutable_ins.size(); i++) {
				mutable_ins[i]->resize(local_start_row, 0, local_num_rows,
						mutable_ins[i]->get_num_cols());
			}
			res->resize(local_start_row, 0, local_num_rows, res->get_num_cols());
			op.run(ins, *res);
		}
	}
	for (size_t i = 0; i < mutable_ins.size(); i++)
		mutable_ins[i]->reset_size();
}

void mapply_store::materialize() const
{
	mapply_store *mutable_this = const_cast<mapply_store *>(this);
	if (lstore->is_whole())
		mutable_this->materialize_whole();
	else if (lstore->store_layout() == matrix_layout_t::L_COL
			&& this->buf.is_empty()) {
		local_buf_col_matrix_store res(lstore->get_global_start_row(),
				lstore->get_global_start_col(), lstore->get_num_rows(),
				lstore->get_num_cols(), lstore->get_type(), -1);
		op.run(ins, res);
		mutable_this->buf = res.get_data();
	}
	else if (lstore->store_layout() == matrix_layout_t::L_COL) {
		local_ref_contig_col_matrix_store res(mutable_this->buf.get_raw(),
				lstore->get_global_start_row(), lstore->get_global_start_col(),
				lstore->get_num_rows(), lstore->get_num_cols(),
				lstore->get_type(), -1);
		op.run(ins, res);
	}
	else if (this->buf.is_empty()) {
		local_buf_row_matrix_store res(lstore->get_global_start_row(),
				lstore->get_global_start_col(), lstore->get_num_rows(),
				lstore->get_num_cols(), lstore->get_type(), -1);
		op.run(ins, res);
		mutable_this->buf = res.get_data();
	}
	else {
		local_ref_contig_row_matrix_store res(mutable_this->buf.get_raw(),
				lstore->get_global_start_row(), lstore->get_global_start_col(),
				lstore->get_num_rows(), lstore->get_num_cols(),
				lstore->get_type(), -1);
		op.run(ins, res);
	}
	mutable_this->materialized = true;
}

class lmapply_col_matrix_store: public lvirtual_col_matrix_store
{
	mapply_store store;
public:
	lmapply_col_matrix_store(
			const std::vector<local_matrix_store::const_ptr> &ins,
			const portion_mapply_op &op,
			off_t global_start_row, off_t global_start_col,
			size_t nrow, size_t ncol, const scalar_type &type,
			int node_id): lvirtual_col_matrix_store(global_start_row,
				global_start_col, nrow, ncol, type, node_id),
			store(ins, op, this) {
	}

	virtual bool resize(off_t local_start_row, off_t local_start_col,
			size_t local_num_rows, size_t local_num_cols) {
		store.resize(local_start_row, local_start_col, local_num_rows,
				local_num_cols);
		return local_matrix_store::resize(local_start_row, local_start_col,
				local_num_rows, local_num_cols);
	}
	virtual void reset_size() {
		store.reset_size();
		local_matrix_store::reset_size();
	}

	virtual const char *get_raw_arr() const {
		if (!store.is_materialized())
			store.materialize();
		return store.get_raw_arr();
	}

	virtual matrix_store::const_ptr transpose() const {
		if (!store.is_materialized())
			store.materialize();
		// TODO we somehow need to make sure the current matrix store isn't
		// destroy when the reference matrix store exists.
		return matrix_store::const_ptr(
				new local_cref_contig_row_matrix_store(store.get_raw_arr(),
					get_global_start_col(), get_global_start_row(),
					get_num_cols(), get_num_rows(), get_type(), get_node_id()));
	}

	virtual const char *get_col(size_t col) const {
		if (!store.is_materialized())
			store.materialize();
		return store.get_raw_arr() + get_num_rows() * col * get_type().get_size();
	}
};

class lmapply_row_matrix_store: public lvirtual_row_matrix_store
{
	mapply_store store;
public:
	lmapply_row_matrix_store(
			const std::vector<local_matrix_store::const_ptr> &ins,
			const portion_mapply_op &op,
			off_t global_start_row, off_t global_start_col,
			size_t nrow, size_t ncol, const scalar_type &type,
			int node_id): lvirtual_row_matrix_store(global_start_row,
				global_start_col, nrow, ncol, type, node_id),
			store(ins, op, this) {
	}

	virtual bool resize(off_t local_start_row, off_t local_start_col,
			size_t local_num_rows, size_t local_num_cols) {
		store.resize(local_start_row, local_start_col, local_num_rows,
				local_num_cols);
		return local_matrix_store::resize(local_start_row, local_start_col,
				local_num_rows, local_num_cols);
	}
	virtual void reset_size() {
		store.reset_size();
		local_matrix_store::reset_size();
	}

	virtual const char *get_raw_arr() const {
		if (!store.is_materialized())
			store.materialize();
		return store.get_raw_arr();
	}

	virtual matrix_store::const_ptr transpose() const {
		if (!store.is_materialized())
			store.materialize();
		// TODO we somehow need to make sure the current matrix store isn't
		// destroy when the reference matrix store exists.
		return matrix_store::const_ptr(
				new local_cref_contig_col_matrix_store(store.get_raw_arr(),
					get_global_start_col(), get_global_start_row(),
					get_num_cols(), get_num_rows(), get_type(), get_node_id()));
	}

	virtual const char *get_row(size_t row) const {
		if (!store.is_materialized())
			store.materialize();
		return store.get_raw_arr() + get_num_cols() * row * get_type().get_size();
	}

	virtual const char *get_rows(size_t row_start, size_t row_end) const {
		// TODO
		assert(0);
		return NULL;
	}
};

}

mapply_matrix_store::mapply_matrix_store(
		const std::vector<mem_matrix_store::const_ptr> &in_mats,
		portion_mapply_op::const_ptr op, matrix_layout_t layout,
		size_t nrow, size_t ncol): virtual_matrix_store(nrow, ncol,
			op->get_output_type())
{
	this->layout = layout;
	this->in_mats = in_mats;
	this->op = op;
}

void mapply_matrix_store::materialize_self() const
{
	// Materialize the matrix store if it hasn't.
	if (res)
		return;

	mapply_matrix_store *mutable_this = const_cast<mapply_matrix_store *>(this);
	mutable_this->res = mem_matrix_store::cast(materialize());
}

matrix_store::ptr mapply_matrix_store::materialize() const
{
	std::vector<matrix_store::const_ptr> tmp(in_mats.begin(), in_mats.end());
	return __mapply_portion(tmp, op, layout);
}

vec_store::const_ptr mapply_matrix_store::get_col_vec(off_t idx) const
{
	if (res == NULL)
		materialize_self();
	return res->get_col_vec(idx);
}

vec_store::const_ptr mapply_matrix_store::get_row_vec(off_t idx) const
{
	if (res == NULL)
		materialize_self();
	return res->get_row_vec(idx);
}

matrix_store::const_ptr mapply_matrix_store::get_cols(
		const std::vector<off_t> &idxs) const
{
	if (res == NULL)
		materialize_self();
	return res->get_cols(idxs);
}

matrix_store::const_ptr mapply_matrix_store::get_rows(
		const std::vector<off_t> &idxs) const
{
	if (res == NULL)
		materialize_self();
	return res->get_rows(idxs);
}

local_matrix_store::const_ptr mapply_matrix_store::get_portion(
			size_t start_row, size_t start_col, size_t num_rows,
			size_t num_cols) const
{
	// If the virtual matrix store has been materialized, we should return
	// the portion from the materialized store directly.
	if (res)
		return res->get_portion(start_row, start_col, num_rows, num_cols);

	std::vector<local_matrix_store::const_ptr> parts(in_mats.size());
	if (is_wide()) {
		assert(start_row == 0);
		assert(num_rows == get_num_rows());
		for (size_t i = 0; i < in_mats.size(); i++)
			parts[i] = in_mats[i]->get_portion(start_row, start_col,
					in_mats[i]->get_num_rows(), num_cols);
	}
	else {
		assert(start_col == 0);
		assert(num_cols == get_num_cols());
		for (size_t i = 0; i < in_mats.size(); i++)
			parts[i] = in_mats[i]->get_portion(start_row, start_col,
					num_rows, in_mats[i]->get_num_cols());
	}

	if (store_layout() == matrix_layout_t::L_ROW)
		return local_matrix_store::const_ptr(new lmapply_row_matrix_store(
					parts, *op, start_row, start_col, num_rows, num_cols,
					get_type(), parts.front()->get_node_id()));
	else
		return local_matrix_store::const_ptr(new lmapply_col_matrix_store(
					parts, *op, start_row, start_col, num_rows, num_cols,
					get_type(), parts.front()->get_node_id()));
}

local_matrix_store::const_ptr mapply_matrix_store::get_portion(
			size_t id) const
{
	// If the virtual matrix store has been materialized, we should return
	// the portion from the materialized store directly.
	if (res)
		return res->get_portion(id);

	size_t start_row;
	size_t start_col;
	size_t num_rows;
	size_t num_cols;
	std::pair<size_t, size_t> chunk_size = get_portion_size();
	if (is_wide()) {
		start_row = 0;
		start_col = chunk_size.second * id;
		num_rows = get_num_rows();
		num_cols = std::min(chunk_size.second, get_num_cols() - start_col);
	}
	else {
		start_row = chunk_size.first * id;
		start_col = 0;
		num_rows = std::min(chunk_size.first, get_num_rows() - start_row);
		num_cols = get_num_cols();
	}
	return get_portion(start_row, start_col, num_rows, num_cols);
}

matrix_store::const_ptr mapply_matrix_store::transpose() const
{
	std::vector<mem_matrix_store::const_ptr> t_in_mats(in_mats.size());
	for (size_t i = 0; i < in_mats.size(); i++)
		t_in_mats[i] = mem_matrix_store::cast(in_mats[i]->transpose());
	matrix_layout_t t_layout;
	if (layout == matrix_layout_t::L_COL)
		t_layout = matrix_layout_t::L_ROW;
	else
		t_layout = matrix_layout_t::L_COL;
	mapply_matrix_store *ret = new mapply_matrix_store(t_in_mats,
			op->transpose(), t_layout, get_num_cols(), get_num_rows());
	if (this->res)
		ret->res = mem_matrix_store::cast(this->res->transpose());
	return matrix_store::const_ptr(ret);
}

}

}
