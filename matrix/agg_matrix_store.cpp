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

local_matrix_store::const_ptr agg_matrix_store::get_portion(
		size_t start_row, size_t start_col, size_t num_rows,
		size_t num_cols) const
{
	local_matrix_store::const_ptr part = data->get_portion(start_row,
			start_col, num_rows, num_cols);
	std::vector<local_matrix_store::const_ptr> ins(1, part);
	portion_op->run(ins);
	return part;
}

local_matrix_store::const_ptr agg_matrix_store::get_portion(size_t id) const
{
	local_matrix_store::const_ptr part = data->get_portion(id);
	std::vector<local_matrix_store::const_ptr> ins(1, part);
	portion_op->run(ins);
	return part;
}

class agg_portion_compute: public portion_compute
{
	std::vector<local_matrix_store::const_ptr> portions;
	std::shared_ptr<portion_mapply_op> op;
	portion_compute::ptr orig_compute;
public:
	agg_portion_compute(std::shared_ptr<portion_mapply_op> op,
			portion_compute::ptr orig_compute) {
		this->orig_compute = orig_compute;
		this->op = op;
		portions.resize(1);
	}

	void set_buf(local_matrix_store::const_ptr portion) {
		portions[0] = portion;
	}

	virtual void run(char *buf, size_t size) {
		assert(!portions.empty());
		op->run(portions);
		// We should also invoke the original compute.
		orig_compute->run(buf, size);
		// This is to make sure the portion compute is only invoked once.
		portions.clear();
	}
};

async_cres_t agg_matrix_store::get_portion_async(
		size_t start_row, size_t start_col, size_t num_rows,
		size_t num_cols, std::shared_ptr<portion_compute> compute) const
{
	agg_portion_compute *_agg_compute = new agg_portion_compute(portion_op,
			compute);
	portion_compute::ptr agg_compute(_agg_compute);
	async_cres_t ret = data->get_portion_async(start_row, start_col,
			num_rows, num_cols, agg_compute);
	// If the data in the portion available, we should run on the data directly.
	if (ret.first) {
		std::vector<local_matrix_store::const_ptr> ins(1, ret.second);
		portion_op->run(ins);
	}
	// Otherwise, agg_compute will run on the portion when its data becomes
	// available.
	else
		_agg_compute->set_buf(ret.second);
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
