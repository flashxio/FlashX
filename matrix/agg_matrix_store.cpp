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
#include "materialize.h"
#include "local_matrix_store.h"
#include "mem_worker_thread.h"

namespace fm
{
namespace detail
{

/*
 * This contains the partial result of aggregation.
 */
class partial_matrix
{
	size_t vec_len;
	// The number of effective elements in a vector.
	size_t num_effect;
	// Each vector contains the local aggregation result in a thread.
	std::vector<smp_vec_store::ptr> res;
	const scalar_type &type;
public:
	typedef std::shared_ptr<partial_matrix> ptr;
	typedef std::shared_ptr<const partial_matrix> const_ptr;

	partial_matrix(size_t num_rows, size_t num_cols,
			const scalar_type &_type): type(_type) {
		res.resize(num_rows);
		num_effect = num_cols;
		// We want each vector has at least 32 elements to avoid false sharing
		// in a CPU cache line.
		vec_len = std::max(num_cols, 32UL);
	}

	bool has_agg(size_t idx) const {
		return res[idx] != NULL;
	}

	bool has_materialized() const {
		bool materialized = false;
		for (size_t i = 0; i < res.size(); i++)
			if (res[i])
				materialized = true;
		return materialized;
	}

	char *get_row(size_t idx) {
		if (res[idx] == NULL) {
			res[idx] = smp_vec_store::create(vec_len, type);
			res[idx]->reset_data();
		}
		return res[idx]->get_raw_arr();
	}

	const char *get_row(size_t idx) const {
		assert(res[idx]);
		return res[idx]->get_raw_arr();
	}

	const scalar_type &get_type() const {
		return type;
	}

	const size_t get_entry_size() const {
		return get_type().get_size();
	}

	size_t get_num_rows() const {
		return res.size();
	}

	size_t get_num_valid_rows() const {
		size_t num_valid_rows = 0;
		for (size_t i = 0; i < res.size(); i++)
			if (res[i])
				num_valid_rows++;
		return num_valid_rows;
	}

	size_t get_num_cols() const {
		return num_effect;
	}

	local_matrix_store::const_ptr get_local_matrix() const;
};

local_matrix_store::const_ptr partial_matrix::get_local_matrix() const
{
	size_t num_valid_rows = get_num_valid_rows();
	size_t entry_size = get_entry_size();
	local_row_matrix_store::ptr tmp(new local_buf_row_matrix_store(0, 0,
				num_valid_rows, num_effect, type, -1));
	size_t copy_row = 0;
	for (size_t i = 0; i < res.size(); i++)
		if (res[i]) {
			memcpy(tmp->get_row(copy_row), res[i]->get_raw_arr(),
					num_effect * entry_size);
			copy_row++;
		}

	return tmp;
}

/*
 * This aggregates on the longer dimension.
 * It outputs a very short vector, so the result is materialized immediately.
 */
class matrix_long_agg_op: public portion_mapply_op
{
	matrix_margin margin;
	agg_operate::const_ptr op;
	// Each row stores the local aggregation results on a thread.
	partial_matrix::ptr partial_res;
	std::vector<local_matrix_store::ptr> local_bufs;
	const size_t data_id;
public:
	typedef std::shared_ptr<const matrix_long_agg_op> const_ptr;

	matrix_long_agg_op(partial_matrix::ptr partial_res,
			matrix_margin margin, agg_operate::const_ptr &op): portion_mapply_op(
				0, 0, partial_res->get_type()), data_id(matrix_store::mat_counter++) {
		this->partial_res = partial_res;
		this->margin = margin;
		this->op = op;
		local_bufs.resize(partial_res->get_num_rows());
	}

	size_t get_data_id() const {
		return data_id;
	}

	partial_matrix::const_ptr get_partial_res() const {
		return partial_res;
	}

	matrix_margin get_margin() const {
		return margin;
	}

	agg_operate::const_ptr get_agg_op() const {
		return op;
	}

	virtual void run(const std::vector<local_matrix_store::const_ptr> &ins) const;

	virtual portion_mapply_op::const_ptr transpose() const {
		matrix_long_agg_op *op = new matrix_long_agg_op(*this);
		op->local_bufs.clear();
		op->local_bufs.resize(partial_res->get_num_rows());
		if (get_margin() == matrix_margin::MAR_ROW)
			op->margin = matrix_margin::MAR_COL;
		else if (get_margin() == matrix_margin::MAR_COL)
			op->margin = matrix_margin::MAR_ROW;
		else
			op->margin = get_margin();
		return portion_mapply_op::const_ptr(op);
	}

	virtual std::string to_string(
			const std::vector<matrix_store::const_ptr> &mats) const {
		assert(mats.size() > 0);
		std::string name;
		if (margin == matrix_margin::MAR_ROW)
			name = "agg_row";
		else
			name = "agg_col";
		name += op->get_agg().get_name() + "(" + mats[0]->get_name();
		for (size_t i = 1; i < mats.size(); i++)
			name += ", " + mats[i]->get_name();
		name += ")";
		return name;
	}
};

void matrix_long_agg_op::run(
		const std::vector<local_matrix_store::const_ptr> &ins) const
{
	assert(ins.size() == 1);
	int thread_id = mem_thread_pool::get_curr_thread_id();
	if (local_bufs[thread_id] == NULL)
		const_cast<matrix_long_agg_op *>(this)->local_bufs[thread_id]
			= local_matrix_store::ptr(new local_buf_col_matrix_store(0, 0,
						partial_res->get_num_cols(), 1, partial_res->get_type(),
						ins[0]->get_node_id()));
	part_dim_t dim = margin
		== matrix_margin::MAR_ROW ? part_dim_t::PART_DIM2 : part_dim_t::PART_DIM1;
	aggregate(*ins[0], *op, margin, dim, *local_bufs[thread_id]);

	// If this is the first time, we should copy the local results to
	// the corresponding row.
	if (!partial_res->has_agg(thread_id))
		memcpy(partial_res->get_row(thread_id),
				local_bufs[thread_id]->get_raw_arr(),
				partial_res->get_num_cols() * partial_res->get_entry_size());
	else
		op->get_combine().runAA(partial_res->get_num_cols(),
				partial_res->get_row(thread_id),
				local_bufs[thread_id]->get_raw_arr(),
				partial_res->get_row(thread_id));
}

static size_t get_num_rows_agg(const matrix_store &data, matrix_margin margin)
{
	if (margin == matrix_margin::MAR_ROW)
		return data.get_num_rows();
	else
		return 1;
}

static size_t get_num_cols_agg(const matrix_store &data, matrix_margin margin)
{
	if (margin == matrix_margin::MAR_COL)
		return data.get_num_cols();
	else
		return 1;
}

agg_matrix_store::agg_matrix_store(matrix_store::const_ptr data,
		matrix_margin margin, agg_operate::const_ptr op): sink_store(
			get_num_rows_agg(*data, margin), get_num_cols_agg(*data, margin),
			data->is_in_mem(), op->get_output_type())
{
	this->data = data;

	size_t num_threads = detail::mem_thread_pool::get_global_num_threads();
	partial_matrix::ptr partial_res;
	if (margin == matrix_margin::BOTH)
		partial_res = partial_matrix::ptr(new partial_matrix(num_threads,
				1, op->get_output_type()));
	// For the next two cases, I assume the partial result is small enough
	// to be kept in memory.
	else if (margin == matrix_margin::MAR_ROW)
		partial_res = partial_matrix::ptr(new partial_matrix(num_threads,
				data->get_num_rows(), op->get_output_type()));
	else if (margin == matrix_margin::MAR_COL)
		partial_res = partial_matrix::ptr(new partial_matrix(num_threads,
				data->get_num_cols(), op->get_output_type()));
	else
		// This shouldn't happen.
		assert(0);
	agg_op = matrix_long_agg_op::const_ptr(new matrix_long_agg_op(
				partial_res, margin, op));
	this->underlying = get_underlying_mats();
}

agg_matrix_store::agg_matrix_store(matrix_store::const_ptr data,
		matrix_long_agg_op::const_ptr portion_op): sink_store(
			get_num_rows_agg(*data, portion_op->get_margin()),
			get_num_cols_agg(*data, portion_op->get_margin()),
			data->is_in_mem(), portion_op->get_agg_op()->get_output_type())
{
	this->data = data;
	this->agg_op = portion_op;
	this->underlying = get_underlying_mats();
}

matrix_store::const_ptr agg_matrix_store::get_agg_res() const
{
	partial_matrix::const_ptr partial_res = agg_op->get_partial_res();
	// The last step is to aggregate the partial results from all portions.
	// It runs in serial. I hope it's not a bottleneck.
	detail::local_matrix_store::const_ptr local_res;
	size_t num_valid_rows = partial_res->get_num_valid_rows();
	// If there is only one row that is valid, we have the final agg results.
	if (num_valid_rows == 1) {
		detail::mem_matrix_store::ptr res = detail::mem_matrix_store::create(
				partial_res->get_num_cols(), 1, matrix_layout_t::L_COL,
				partial_res->get_type(), -1);
		off_t valid_row = -1;
		for (size_t i = 0; i < partial_res->get_num_rows(); i++) {
			if (partial_res->has_agg(i)) {
				valid_row = i;
				break;
			}
		}
		assert(valid_row >= 0);
		memcpy(res->get_raw_arr(), partial_res->get_row(valid_row),
				partial_res->get_num_cols() * partial_res->get_entry_size());
		if ((get_num_rows() == 1 && get_num_cols() == 1)
				|| (get_num_rows() > 1 && get_num_cols() == 1)) {
			assert(res->get_num_rows() == get_num_rows());
			assert(res->get_num_cols() == get_num_cols());
			return res;
		}
		else {
			assert(res->get_num_rows() == get_num_cols());
			assert(res->get_num_cols() == get_num_rows());
			return res->transpose();
		}
	}
	else
		local_res = partial_res->get_local_matrix();

	detail::mem_matrix_store::ptr res = detail::mem_matrix_store::create(
			partial_res->get_num_cols(), 1, matrix_layout_t::L_COL,
			partial_res->get_type(), -1);
	local_matrix_store::ptr portion = res->get_portion(0);
	assert(res->get_num_rows() == portion->get_num_rows());
	assert(res->get_num_cols() == portion->get_num_cols());
	// I need to create new aggregation with the combine operation
	// to run aggregation on the columns of the matrix.
	agg_operate::const_ptr combine_agg = agg_operate::create(
			agg_op->get_agg_op()->get_combine_ptr(), bulk_operate::const_ptr());
	detail::aggregate(*local_res, *combine_agg, matrix_margin::MAR_COL,
			part_dim_t::PART_NONE, *portion);
	if ((get_num_rows() == 1 && get_num_cols() == 1)
			|| (get_num_rows() > 1 && get_num_cols() == 1)) {
		assert(res->get_num_rows() == get_num_rows());
		assert(res->get_num_cols() == get_num_cols());
		return res;
	}
	else {
		assert(res->get_num_rows() == get_num_cols());
		assert(res->get_num_cols() == get_num_rows());
		return res->transpose();
	}
}

bool agg_matrix_store::has_materialized() const
{
	return agg_op->get_partial_res()->has_materialized();
}

void agg_matrix_store::materialize_self() const
{
	if (!has_materialized()) {
		// This computes the partial aggregation result.
		std::vector<detail::matrix_store::const_ptr> ins(1);
		ins[0] = data;
		__mapply_portion(ins, agg_op, matrix_layout_t::L_ROW);
	}
}

size_t agg_matrix_store::get_data_id() const
{
	return agg_op->get_data_id();
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
			const scalar_type &type,
			portion_mapply_op::const_ptr portion_op): lvirtual_col_matrix_store(
				part->get_global_start_row(), part->get_global_start_col(),
				part->get_num_rows(), part->get_num_cols(), type,
				part->get_node_id()), parts(1, part),
			mutable_part(const_cast<local_matrix_store &>(*part)) {
		this->portion_op = portion_op;
	}

	virtual bool resize(off_t local_start_row, off_t local_start_col,
			size_t local_num_rows, size_t local_num_cols) {
		bool success = mutable_part.resize(local_start_row, local_start_col,
				local_num_rows, local_num_cols);
		if (!success)
			return false;
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
			const scalar_type &type,
			portion_mapply_op::const_ptr portion_op): lvirtual_row_matrix_store(
				part->get_global_start_row(), part->get_global_start_col(),
				part->get_num_rows(), part->get_num_cols(), type,
				part->get_node_id()), parts(1, part),
			mutable_part(const_cast<local_matrix_store &>(*part)) {
		this->portion_op = portion_op;
	}

	virtual bool resize(off_t local_start_row, off_t local_start_col,
			size_t local_num_rows, size_t local_num_cols) {
		bool success = mutable_part.resize(local_start_row, local_start_col,
				local_num_rows, local_num_cols);
		if (!success)
			return false;
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

matrix_store::const_ptr agg_matrix_store::transpose() const
{
	if (has_materialized()) {
		matrix_store::const_ptr res = get_agg_res();
		return res->transpose();
	}
	matrix_store::const_ptr tdata = data->transpose();
	portion_mapply_op::const_ptr top = agg_op->transpose();
	matrix_long_agg_op::const_ptr tagg_op
		= std::dynamic_pointer_cast<const matrix_long_agg_op>(top);
	return matrix_store::const_ptr(new agg_matrix_store(tdata, tagg_op));
}

class agg_compute_store: public sink_compute_store, public EM_object
{
	matrix_long_agg_op::const_ptr portion_op;
	matrix_store::const_ptr data;
public:
	agg_compute_store(matrix_store::const_ptr data,
			portion_mapply_op::const_ptr portion_op): sink_compute_store(
				data->get_num_rows(), data->get_num_cols(), data->is_in_mem(),
				portion_op->get_output_type()) {
		this->data = data;
		this->portion_op
			= std::dynamic_pointer_cast<const matrix_long_agg_op>(portion_op);
		assert(this->portion_op);
	}
	virtual size_t get_data_id() const {
		return portion_op->get_data_id();
	}
	using virtual_matrix_store::get_portion;
	virtual std::shared_ptr<const local_matrix_store> get_portion(
			size_t start_row, size_t start_col, size_t num_rows,
			size_t num_cols) const;
	virtual std::shared_ptr<const local_matrix_store> get_portion(
			size_t id) const;
	using virtual_matrix_store::get_portion_async;
	virtual async_cres_t get_portion_async(
			size_t start_row, size_t start_col, size_t num_rows,
			size_t num_cols, std::shared_ptr<portion_compute> compute) const;

	virtual int get_portion_node_id(size_t id) const {
		return data->get_portion_node_id(id);
	}

	virtual std::pair<size_t, size_t> get_portion_size() const {
		return data->get_portion_size();
	}

	virtual int get_num_nodes() const {
		return data->get_num_nodes();
	}

	virtual matrix_layout_t store_layout() const {
		return data->store_layout();
	}
	virtual std::string get_name() const {
		std::vector<matrix_store::const_ptr> mats(1);
		mats[0] = data;
		return portion_op->to_string(mats);
	}

	virtual std::vector<safs::io_interface::ptr> create_ios() const;
	virtual std::unordered_map<size_t, size_t> get_underlying_mats() const {
		return data->get_underlying_mats();
	}
};

static local_matrix_store::const_ptr create_lmaterialize_matrix(
		local_matrix_store::const_ptr part, const scalar_type &type,
		portion_mapply_op::const_ptr portion_op)
{
	if (part->store_layout() == matrix_layout_t::L_ROW)
		return local_matrix_store::const_ptr(new lmaterialize_row_matrix_store(
					part, type, portion_op));
	else
		return local_matrix_store::const_ptr(new lmaterialize_col_matrix_store(
					part, type, portion_op));
}

local_matrix_store::const_ptr agg_compute_store::get_portion(
		size_t start_row, size_t start_col, size_t num_rows,
		size_t num_cols) const
{
	local_matrix_store::const_ptr part = data->get_portion(start_row,
			start_col, num_rows, num_cols);
	return create_lmaterialize_matrix(part, get_type(), portion_op);
}

local_matrix_store::const_ptr agg_compute_store::get_portion(size_t id) const
{
	local_matrix_store::const_ptr part = data->get_portion(id);
	return create_lmaterialize_matrix(part, get_type(), portion_op);
}

async_cres_t agg_compute_store::get_portion_async(
		size_t start_row, size_t start_col, size_t num_rows,
		size_t num_cols, std::shared_ptr<portion_compute> compute) const
{
	async_cres_t ret = data->get_portion_async(start_row, start_col,
			num_rows, num_cols, compute);
	ret.second = create_lmaterialize_matrix(ret.second, get_type(), portion_op);
	return ret;
}

std::vector<safs::io_interface::ptr> agg_compute_store::create_ios() const
{
	const EM_object *obj = dynamic_cast<const EM_object *>(data.get());
	if (obj)
		return obj->create_ios();
	else
		return std::vector<safs::io_interface::ptr>();
}

std::vector<virtual_matrix_store::const_ptr> agg_matrix_store::get_compute_matrices() const
{
	// If the agg matrix has been materialized, we don't need to do
	// anything.
	if (has_materialized())
		return std::vector<virtual_matrix_store::const_ptr>();
	else if (data->is_wide()) {
		portion_mapply_op::const_ptr top = agg_op->transpose();
		return std::vector<virtual_matrix_store::const_ptr>(1,
				virtual_matrix_store::const_ptr(new agg_compute_store(
						data->transpose(), top)));
	}
	else
		return std::vector<virtual_matrix_store::const_ptr>(1,
				virtual_matrix_store::const_ptr(new agg_compute_store(data,
						agg_op)));
}

std::string agg_matrix_store::get_name() const
{
	std::vector<matrix_store::const_ptr> mats(1);
	mats[0] = data;
	return agg_op->to_string(mats);
}

}

}
