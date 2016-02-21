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

#include <numeric>

#include "agg_matrix_store.h"
#include "dense_matrix.h"
#include "local_vec_store.h"
#include "local_matrix_store.h"
#include "mem_worker_thread.h"
#include "EM_object.h"

namespace fm
{
namespace detail
{

/*
 * This aggregates on the longer dimension.
 * It outputs a very short vector, so the result is materialized immediately.
 */
class matrix_long_agg_op: public portion_mapply_op
{
	matrix_margin margin;
	agg_operate::const_ptr op;
	// Each row stores the local aggregation results on a thread.
	mem_row_matrix_store::ptr partial_res;
	std::vector<local_vec_store::ptr> local_bufs;
	// This indicates the number of elements that have been aggregated to
	// the partial results.
	std::vector<size_t> num_aggs;
public:
	matrix_long_agg_op(mem_row_matrix_store::ptr partial_res,
			matrix_margin margin, agg_operate::const_ptr &op): portion_mapply_op(
				0, 0, partial_res->get_type()) {
		this->partial_res = partial_res;
		this->margin = margin;
		this->op = op;
		local_bufs.resize(partial_res->get_num_rows());
		num_aggs.resize(partial_res->get_num_rows());
	}

	size_t get_num_aggs() const {
		return std::accumulate(num_aggs.begin(), num_aggs.end(), 0);
	}

	const agg_operate &get_agg_op() const {
		return *op;
	}

	bool valid_row(size_t off) const {
		return num_aggs[off] > 0;
	}

	size_t get_num_valid_rows() const {
		size_t ret = 0;
		for (size_t i = 0; i < num_aggs.size(); i++)
			if (num_aggs[i] > 0)
				ret++;
		return ret;
	}

	virtual void run(const std::vector<local_matrix_store::const_ptr> &ins) const;

	virtual portion_mapply_op::const_ptr transpose() const {
		assert(0);
		return portion_mapply_op::const_ptr();
	}

	virtual std::string to_string(
			const std::vector<matrix_store::const_ptr> &mats) const {
		return std::string();
	}
};

void matrix_long_agg_op::run(
		const std::vector<local_matrix_store::const_ptr> &ins) const
{
	assert(ins.size() == 1);
	int thread_id = mem_thread_pool::get_curr_thread_id();
	if (local_bufs[thread_id] == NULL)
		const_cast<matrix_long_agg_op *>(this)->local_bufs[thread_id]
			= local_vec_store::ptr(new local_buf_vec_store(0,
						partial_res->get_num_cols(), partial_res->get_type(),
						ins[0]->get_node_id()));
	aggregate(*ins[0], *op, margin, *local_bufs[thread_id]);

	// If this is the first time, we should copy the local results to
	// the corresponding row.
	if (num_aggs[thread_id] == 0)
		memcpy(partial_res->get_row(thread_id),
				local_bufs[thread_id]->get_raw_arr(),
				partial_res->get_num_cols() * partial_res->get_entry_size());
	else
		op->get_combine().runAA(partial_res->get_num_cols(),
				partial_res->get_row(thread_id),
				local_bufs[thread_id]->get_raw_arr(),
				partial_res->get_row(thread_id));
	const_cast<matrix_long_agg_op *>(this)->num_aggs[thread_id]
		+= ins[0]->get_num_rows() * ins[0]->get_num_cols();
}

agg_matrix_store::agg_matrix_store(matrix_store::const_ptr data,
		matrix_margin margin, agg_operate::const_ptr op): virtual_matrix_store(
			data->get_num_rows(), data->get_num_cols(), data->is_in_mem(),
			data->get_type())
{
	this->data = data;

	size_t num_threads = detail::mem_thread_pool::get_global_num_threads();
	detail::mem_row_matrix_store::ptr partial_res;
	if (margin == matrix_margin::BOTH)
		partial_res = detail::mem_row_matrix_store::create(num_threads,
				1, op->get_output_type());
	// For the next two cases, I assume the partial result is small enough
	// to be kept in memory.
	else if (margin == matrix_margin::MAR_ROW)
		partial_res = detail::mem_row_matrix_store::create(num_threads,
				data->get_num_rows(), op->get_output_type());
	else if (margin == matrix_margin::MAR_COL)
		partial_res = detail::mem_row_matrix_store::create(num_threads,
				data->get_num_cols(), op->get_output_type());
	else
		// This shouldn't happen.
		assert(0);
	partial_res->reset_data();

	portion_op = std::shared_ptr<matrix_long_agg_op>(new matrix_long_agg_op(
				partial_res, margin, op));
	this->partial_res = partial_res;
}

matrix_store::ptr agg_matrix_store::get_agg_res() const
{
	std::shared_ptr<matrix_long_agg_op> agg_op
		= std::static_pointer_cast<matrix_long_agg_op>(portion_op);
	// The last step is to aggregate the partial results from all portions.
	// It runs in serial. I hope it's not a bottleneck.
	detail::local_matrix_store::const_ptr local_res;
	size_t num_valid_rows = agg_op->get_num_valid_rows();
	// If there is only one row that is valid, we have the final agg results.
	if (num_valid_rows == 1) {
		detail::mem_matrix_store::ptr res = detail::mem_matrix_store::create(
				partial_res->get_num_cols(), 1, matrix_layout_t::L_COL,
				partial_res->get_type(), -1);
		memcpy(res->get_raw_arr(), partial_res->get_row(0),
				partial_res->get_num_cols() * partial_res->get_entry_size());
		return res;
	}
	// If not, we need to combine the partial aggregation results.
	// If all rows are valid.
	else if (num_valid_rows == partial_res->get_num_rows())
		local_res = partial_res->get_portion(
				0, 0, partial_res->get_num_rows(), partial_res->get_num_cols());
	else {
		// Otherwise, we have to pick all the valid rows out.
		detail::local_row_matrix_store::ptr tmp(
				new detail::local_buf_row_matrix_store(0, 0, num_valid_rows,
					partial_res->get_num_cols(), partial_res->get_type(), -1));
		size_t entry_size = partial_res->get_entry_size();
		size_t copy_row = 0;
		for (size_t i = 0; i < partial_res->get_num_rows(); i++)
			if (agg_op->valid_row(i)) {
				memcpy(tmp->get_row(copy_row), partial_res->get_row(i),
						partial_res->get_num_cols() * entry_size);
				copy_row++;
			}
		assert(copy_row == num_valid_rows);
		local_res = tmp;
	}

	detail::mem_matrix_store::ptr res = detail::mem_matrix_store::create(
			partial_res->get_num_cols(), 1, matrix_layout_t::L_COL,
			partial_res->get_type(), -1);
	local_ref_vec_store local_vec(res->get_raw_arr(), 0, res->get_num_rows(),
			res->get_type(), -1);
	// I need to create new aggregation with the combine operation
	// to run aggregation on the columns of the matrix.
	agg_operate::const_ptr combine_agg = agg_operate::create(
			agg_op->get_agg_op().get_combine_ptr(), bulk_operate::const_ptr());
	detail::aggregate(*local_res, *combine_agg, matrix_margin::MAR_COL,
			local_vec);
	return res;
}

vec_store::const_ptr agg_matrix_store::get_col_vec(off_t idx) const
{
	// An agg matrix is always a one-column matrix.
	// We should support this method.
	matrix_store::const_ptr ret = materialize(true, -1);
	return ret->get_col_vec(idx);
}

vec_store::const_ptr agg_matrix_store::get_row_vec(off_t idx) const
{
	assert(0);
	return vec_store::const_ptr();
}

matrix_store::const_ptr agg_matrix_store::get_cols(
		const std::vector<off_t> &idxs) const
{
	assert(0);
	return matrix_store::const_ptr();
}

matrix_store::const_ptr agg_matrix_store::get_rows(
		const std::vector<off_t> &idxs) const
{
	assert(0);
	return matrix_store::const_ptr();
}

bool agg_matrix_store::has_materialized() const
{
	std::shared_ptr<matrix_long_agg_op> agg_op
		= std::static_pointer_cast<matrix_long_agg_op>(portion_op);
	return agg_op->get_num_aggs() == get_num_rows() * get_num_cols();
}

void agg_matrix_store::materialize_self() const
{
	if (!has_materialized()) {
		// This computes the partial aggregation result.
		std::vector<detail::matrix_store::const_ptr> ins(1);
		ins[0] = data;
		__mapply_portion(ins, portion_op, matrix_layout_t::L_ROW);
	}
}

matrix_store::const_ptr agg_matrix_store::materialize(bool in_mem,
		int num_nodes) const
{
	materialize_self();
	return get_agg_res();
}


namespace
{

/*
 * The role of these two matrices is to materialize the underlying local matrix
 * piece by piece so that we can keep data in the CPU cache when computing
 * aggregation.
 */

class lmaterialize_col_matrix_store: public lvirtual_col_matrix_store
{
	std::vector<local_matrix_store::const_ptr> parts;
	local_matrix_store &mutable_part;
	portion_mapply_op::const_ptr portion_op;
public:
	lmaterialize_col_matrix_store(local_matrix_store::const_ptr part,
			portion_mapply_op::const_ptr portion_op): lvirtual_col_matrix_store(
				part->get_global_start_row(), part->get_global_start_col(),
				part->get_num_rows(), part->get_num_cols(), part->get_type(),
				part->get_node_id()), parts(1, part),
			mutable_part(const_cast<local_matrix_store &>(*part)) {
		this->portion_op = portion_op;
	}

	virtual bool resize(off_t local_start_row, off_t local_start_col,
			size_t local_num_rows, size_t local_num_cols) {
		mutable_part.resize(local_start_row, local_start_col, local_num_rows,
				local_num_cols);
		return local_matrix_store::resize(local_start_row, local_start_col,
				local_num_rows, local_num_cols);
	}
	virtual void reset_size() {
		mutable_part.reset_size();
		local_matrix_store::reset_size();
	}

	using lvirtual_col_matrix_store::get_raw_arr;
	virtual const char *get_raw_arr() const {
		assert(0);
		return NULL;
	}

	using lvirtual_col_matrix_store::transpose;
	virtual matrix_store::const_ptr transpose() const {
		assert(0);
		return matrix_store::const_ptr();
	}

	using lvirtual_col_matrix_store::get_col;
	virtual const char *get_col(size_t col) const {
		assert(0);
		return NULL;
	}

	virtual local_matrix_store::const_ptr get_portion(
			size_t local_start_row, size_t local_start_col, size_t num_rows,
			size_t num_cols) const {
		assert(0);
		return local_matrix_store::const_ptr();
	}

	virtual void materialize_self() const {
		portion_op->run(parts);
	}
};

class lmaterialize_row_matrix_store: public lvirtual_row_matrix_store
{
	std::vector<local_matrix_store::const_ptr> parts;
	local_matrix_store &mutable_part;
	portion_mapply_op::const_ptr portion_op;
public:
	lmaterialize_row_matrix_store(local_matrix_store::const_ptr part,
			portion_mapply_op::const_ptr portion_op): lvirtual_row_matrix_store(
				part->get_global_start_row(), part->get_global_start_col(),
				part->get_num_rows(), part->get_num_cols(), part->get_type(),
				part->get_node_id()), parts(1, part),
			mutable_part(const_cast<local_matrix_store &>(*part)) {
		this->portion_op = portion_op;
	}

	virtual bool resize(off_t local_start_row, off_t local_start_col,
			size_t local_num_rows, size_t local_num_cols) {
		mutable_part.resize(local_start_row, local_start_col, local_num_rows,
				local_num_cols);
		return local_matrix_store::resize(local_start_row, local_start_col,
				local_num_rows, local_num_cols);
	}
	virtual void reset_size() {
		mutable_part.reset_size();
		local_matrix_store::reset_size();
	}

	using lvirtual_row_matrix_store::get_raw_arr;
	virtual const char *get_raw_arr() const {
		assert(0);
		return NULL;
	}

	using lvirtual_row_matrix_store::transpose;
	virtual matrix_store::const_ptr transpose() const {
		assert(0);
		return matrix_store::const_ptr();
	}

	using lvirtual_row_matrix_store::get_row;
	virtual const char *get_row(size_t row) const {
		assert(0);
		return NULL;
	}

	virtual local_matrix_store::const_ptr get_portion(
			size_t local_start_row, size_t local_start_col, size_t num_rows,
			size_t num_cols) const {
		assert(0);
		return local_matrix_store::const_ptr();
	}

	virtual void materialize_self() const {
		portion_op->run(parts);
	}
};

}

static local_matrix_store::const_ptr create_lmaterialize_matrix(
		local_matrix_store::const_ptr part,
		portion_mapply_op::const_ptr portion_op)
{
	if (part->store_layout() == matrix_layout_t::L_ROW)
		return local_matrix_store::const_ptr(new lmaterialize_row_matrix_store(
					part, portion_op));
	else
		return local_matrix_store::const_ptr(new lmaterialize_col_matrix_store(
					part, portion_op));
}

local_matrix_store::const_ptr agg_matrix_store::get_portion(
		size_t start_row, size_t start_col, size_t num_rows,
		size_t num_cols) const
{
	local_matrix_store::const_ptr part = data->get_portion(start_row,
			start_col, num_rows, num_cols);
	return create_lmaterialize_matrix(part, portion_op);
}

local_matrix_store::const_ptr agg_matrix_store::get_portion(size_t id) const
{
	local_matrix_store::const_ptr part = data->get_portion(id);
	return create_lmaterialize_matrix(part, portion_op);
}

async_cres_t agg_matrix_store::get_portion_async(
		size_t start_row, size_t start_col, size_t num_rows,
		size_t num_cols, std::shared_ptr<portion_compute> compute) const
{
	async_cres_t ret = data->get_portion_async(start_row, start_col,
			num_rows, num_cols, compute);
	ret.second = create_lmaterialize_matrix(ret.second, portion_op);
	return ret;
}

matrix_store::const_ptr agg_matrix_store::transpose() const
{
	// TODO This method should also be implemented.
	assert(0);
	return matrix_store::const_ptr();
}

std::vector<safs::io_interface::ptr> agg_matrix_store::create_ios() const
{
	const EM_object *obj = dynamic_cast<const EM_object *>(data.get());
	if (obj)
		return obj->create_ios();
	else
		return std::vector<safs::io_interface::ptr>();
}

}

}
