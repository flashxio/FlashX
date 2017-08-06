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

#include "agg_matrix_store.h"
#include "groupby_matrix_store.h"
#include "local_matrix_store.h"
#include "factor.h"

namespace fm
{

namespace detail
{

namespace
{

class groupby_op: public detail::portion_mapply_op
{
	// This contains a bool vector for each thread.
	// The bool vector indicates whether a label gets partially aggregated data.
	std::vector<std::vector<bool> > part_agg;
	// This contains a local matrix for each thread.
	// Each row/column of a local matrix contains partially aggregated data
	// for a label.
	std::vector<detail::local_matrix_store::ptr> part_results;
	std::vector<bool> part_status;
	size_t num_levels;
	matrix_margin margin;
	agg_operate::const_ptr op;
	const size_t data_id;
public:
	groupby_op(agg_operate::const_ptr op, size_t num_levels,
			matrix_margin margin): detail::portion_mapply_op(0, 0,
				op->get_output_type()), data_id(matrix_store::mat_counter++) {
		size_t num_threads = detail::mem_thread_pool::get_global_num_threads();
		part_results.resize(num_threads);
		part_agg.resize(num_threads);
		part_status.resize(num_threads, true);
		this->num_levels = num_levels;
		this->margin = margin;
		this->op = op;
	}
	size_t get_data_id() const {
		return data_id;
	}

	virtual detail::portion_mapply_op::const_ptr transpose() const {
		throw unsupported_exception("Don't support transpose of groupby_op");
	}

	virtual void run(
			const std::vector<detail::local_matrix_store::const_ptr> &ins) const;

	virtual std::string to_string(
			const std::vector<detail::matrix_store::const_ptr> &mats) const {
		std::string str = margin == matrix_margin::MAR_ROW ? "row" : "col";
		return std::string("groupby_") + str + "(" + mats[0]->get_name() + ")";
	}

	detail::matrix_store::ptr get_agg() const;
	bool has_materialized() const;
};

detail::matrix_store::ptr groupby_op::get_agg() const
{
	for (size_t i = 0; i < part_status.size(); i++) {
		if (!part_status[i]) {
			BOOST_LOG_TRIVIAL(error) << "groupby fails on a partition";
			return detail::matrix_store::ptr();
		}
	}
	size_t first_idx;
	for (first_idx = 0; first_idx < part_results.size(); first_idx++)
		if (part_results[first_idx] != NULL)
			break;
	assert(first_idx < part_results.size());

	size_t nrow = part_results[first_idx]->get_num_rows();
	size_t ncol = part_results[first_idx]->get_num_cols();
	const scalar_type &type = part_results[first_idx]->get_type();
	if (part_results[first_idx]->store_layout() == matrix_layout_t::L_ROW) {
		detail::mem_row_matrix_store::ptr res
			= detail::mem_row_matrix_store::create(nrow, ncol, type);
		for (size_t i = 0; i < res->get_num_rows(); i++) {
			const local_row_matrix_store &row_res
				= static_cast<const local_row_matrix_store &>(
						*part_results[first_idx]);
			memcpy(res->get_row(i), row_res.get_row(i),
					res->get_num_cols() * res->get_entry_size());
			for (size_t j = first_idx + 1; j < part_results.size(); j++) {
				if (part_results[j] != NULL) {
					const local_row_matrix_store &row_res
						= static_cast<const local_row_matrix_store &>(
								*part_results[j]);
					op->get_combine().runAA(res->get_num_cols(),
							row_res.get_row(i), res->get_row(i), res->get_row(i));
				}
			}
		}
		return res;
	}
	else {
		detail::mem_col_matrix_store::ptr res
			= detail::mem_col_matrix_store::create(nrow, ncol, type);
		for (size_t i = 0; i < res->get_num_cols(); i++) {
			const local_col_matrix_store &col_res
				= static_cast<const local_col_matrix_store &>(
						*part_results[first_idx]);
			memcpy(res->get_col(i), col_res.get_col(i),
					res->get_num_rows() * res->get_entry_size());
			for (size_t j = first_idx + 1; j < part_results.size(); j++) {
				if (part_results[j] != NULL) {
					const local_col_matrix_store &col_res
						= static_cast<const local_col_matrix_store &>(
								*part_results[j]);
					op->get_combine().runAA(res->get_num_rows(),
							col_res.get_col(i), res->get_col(i), res->get_col(i));
				}
			}
		}
		return res;
	}
}

bool groupby_op::has_materialized() const
{
	// If we have materialized the groupby operation, we should have
	// partial results from at least one thread.
	for (size_t i = 0; i < part_results.size(); i++)
		if (part_results[i])
			return true;
	return false;
}

void groupby_op::run(
		const std::vector<detail::local_matrix_store::const_ptr> &ins) const
{
	assert(ins.size() == 2);
	detail::local_matrix_store::const_ptr labels = ins[0];
	detail::local_matrix_store::const_ptr in = ins[1];

	groupby_op *mutable_this = const_cast<groupby_op *>(this);
	// Prepare for the output result.
	int thread_id = detail::mem_thread_pool::get_curr_thread_id();
	if (part_results[thread_id] == NULL) {
		assert(part_agg[thread_id].empty());
		if (margin == matrix_margin::MAR_ROW) {
			mutable_this->part_results[thread_id]
				= detail::local_matrix_store::ptr(
					new detail::local_buf_row_matrix_store(0, 0, num_levels,
						in->get_num_cols(), op->get_output_type(), -1));
			assert(in->store_layout() == matrix_layout_t::L_ROW);
		}
		else {
			mutable_this->part_results[thread_id]
				= detail::local_matrix_store::ptr(
					new detail::local_buf_col_matrix_store(0, 0,
						in->get_num_rows(), num_levels,
						op->get_output_type(), -1));
			assert(in->store_layout() == matrix_layout_t::L_COL);
		}
		mutable_this->part_results[thread_id]->reset_data();
		mutable_this->part_agg[thread_id].resize(num_levels);
	}
	// If there was a failure in this thread, we don't need to perform more
	// computation.
	if (!part_status[thread_id])
		return;

	detail::part_dim_t dim = margin == matrix_margin::MAR_ROW
		? detail::part_dim_t::PART_DIM1 : detail::part_dim_t::PART_DIM2;
	bool ret = detail::groupby(*labels, *in, *op, margin, dim,
			*part_results[thread_id], mutable_this->part_agg[thread_id]);
	if (!ret)
		mutable_this->part_status[thread_id] = false;
}

}

static size_t get_num_rows_groupby(const matrix_store &mat,
		size_t num_levels, matrix_margin margin)
{
	return margin == matrix_margin::MAR_ROW ? num_levels : mat.get_num_rows();
}

static size_t get_num_cols_groupby(const matrix_store &mat,
		size_t num_levels, matrix_margin margin)
{
	return margin == matrix_margin::MAR_COL ? num_levels : mat.get_num_cols();
}

static matrix_store::const_ptr conv_layout(matrix_store::const_ptr data,
		matrix_layout_t layout)
{
	if ((layout == matrix_layout_t::L_ROW
				&& data->store_layout() == matrix_layout_t::L_ROW)
			|| (layout == matrix_layout_t::L_COL
				&& data->store_layout() == matrix_layout_t::L_COL))
		return data;
	else {
		dense_matrix::ptr tmp = dense_matrix::create(data);
		tmp = tmp->conv2(layout);
		return tmp->get_raw_store();
	}
}

groupby_matrix_store::groupby_matrix_store(matrix_store::const_ptr data,
		matrix_store::const_ptr label_store, const factor &_f,
		matrix_margin margin, agg_operate::const_ptr op): sink_store(
			get_num_rows_groupby(*data, _f.get_num_levels(), margin),
			get_num_cols_groupby(*data, _f.get_num_levels(), margin),
			data->is_in_mem(), op->get_output_type()), f(_f)
{
	this->label_store = label_store;
	if (margin == matrix_margin::MAR_ROW) {
		assert(this->label_store->get_num_rows() == data->get_num_rows());
		this->data = conv_layout(data, matrix_layout_t::L_ROW);
		assert(this->data->get_num_rows() >= this->data->get_num_cols());
	}
	else {
		assert(this->label_store->get_num_cols() == data->get_num_cols());
		this->data = conv_layout(data, matrix_layout_t::L_COL);
		assert(this->data->get_num_cols() >= this->data->get_num_rows());
	}
	this->f = f;
	this->margin = margin;
	portion_op = std::shared_ptr<groupby_op>(new groupby_op(op,
				f.get_num_levels(), margin));
	agg_op = op;
	this->underlying = get_underlying_mats();
}

groupby_matrix_store::groupby_matrix_store(matrix_store::const_ptr data,
		factor_col_vector::const_ptr labels, matrix_margin margin,
		agg_operate::const_ptr op): sink_store(
			get_num_rows_groupby(*data, labels->get_num_levels(), margin),
			get_num_cols_groupby(*data, labels->get_num_levels(), margin),
			data->is_in_mem(), op->get_output_type()), f(labels->get_factor())
{
	this->data = data;
	// labels is a col matrix. If we group by columns, we need to transpose
	// the col matrix.
	if (margin == matrix_margin::MAR_ROW) {
		this->data = conv_layout(data, matrix_layout_t::L_ROW);
		this->label_store = labels->get_raw_store();
		assert(this->label_store->get_num_rows() == data->get_num_rows());
		assert(this->data->get_num_rows() >= this->data->get_num_cols());
	}
	else {
		this->data = conv_layout(data, matrix_layout_t::L_COL);
		this->label_store = labels->get_raw_store()->transpose();
		assert(data->store_layout() == matrix_layout_t::L_COL);
		assert(this->label_store->get_num_cols() == data->get_num_cols());
		assert(this->data->get_num_cols() >= this->data->get_num_rows());
	}
	this->margin = margin;
	portion_op = std::shared_ptr<groupby_op>(new groupby_op(op,
				labels->get_factor().get_num_levels(), margin));
	agg_op = op;
	this->underlying = get_underlying_mats();
}

matrix_store::ptr groupby_matrix_store::get_agg_res() const
{
	return std::static_pointer_cast<groupby_op>(portion_op)->get_agg();
}

size_t groupby_matrix_store::get_data_id() const
{
	return std::static_pointer_cast<groupby_op>(portion_op)->get_data_id();
}

bool groupby_matrix_store::has_materialized() const
{
	return std::static_pointer_cast<groupby_op>(portion_op)->has_materialized();
}

void groupby_matrix_store::materialize_self() const
{
	if (!has_materialized()) {
		// This computes the partial aggregation result.
		std::vector<detail::matrix_store::const_ptr> ins(2);
		ins[0] = label_store;
		ins[1] = data;
		__mapply_portion(ins, portion_op, matrix_layout_t::L_ROW);
	}
}

matrix_store::const_ptr groupby_matrix_store::materialize(bool in_mem,
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
	local_matrix_store &mutable_data;
	local_matrix_store &mutable_labels;
	portion_mapply_op::const_ptr portion_op;
public:
	lmaterialize_col_matrix_store(local_matrix_store::const_ptr data,
			local_matrix_store::const_ptr labels, const scalar_type &type,
			portion_mapply_op::const_ptr portion_op): lvirtual_col_matrix_store(
				data->get_global_start_row(), data->get_global_start_col(),
				data->get_num_rows(), data->get_num_cols(), type,
				data->get_node_id()), mutable_data(
				const_cast<local_matrix_store &>(*data)), mutable_labels(
				const_cast<local_matrix_store &>(*labels)) {
		this->portion_op = portion_op;
		parts.resize(2);
		parts[0] = labels;
		parts[1] = data;
	}

	virtual bool resize(off_t local_start_row, off_t local_start_col,
			size_t local_num_rows, size_t local_num_cols) {
		// We should resize the number of columns.
		assert(local_start_row == 0
				&& local_num_rows == mutable_data.get_num_rows());
		// labels is a column vector, so we can only resize the number of rows.
		local_matrix_store::exposed_area orig_labels
			= mutable_labels.get_exposed_area();
		bool success = mutable_labels.resize(local_start_col, 0, local_num_cols, 1);
		if (!success)
			return false;

		success = mutable_data.resize(local_start_row, local_start_col, local_num_rows,
				local_num_cols);
		if (!success) {
			mutable_labels.restore_size(orig_labels);
			return false;
		}
		return local_matrix_store::resize(local_start_row, local_start_col,
				local_num_rows, local_num_cols);
	}
	virtual void reset_size() {
		mutable_data.reset_size();
		mutable_labels.reset_size();
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
	local_matrix_store &mutable_data;
	local_matrix_store &mutable_labels;
	portion_mapply_op::const_ptr portion_op;
public:
	lmaterialize_row_matrix_store(local_matrix_store::const_ptr data,
			local_matrix_store::const_ptr labels, const scalar_type &type,
			portion_mapply_op::const_ptr portion_op): lvirtual_row_matrix_store(
				data->get_global_start_row(), data->get_global_start_col(),
				data->get_num_rows(), data->get_num_cols(), type,
				data->get_node_id()), mutable_data(
				const_cast<local_matrix_store &>(*data)), mutable_labels(
				const_cast<local_matrix_store &>(*labels)) {
		this->portion_op = portion_op;
		parts.resize(2);
		parts[0] = labels;
		parts[1] = data;
	}

	virtual bool resize(off_t local_start_row, off_t local_start_col,
			size_t local_num_rows, size_t local_num_cols) {
		// We should resize the number of rows.
		assert(local_start_col == 0
				&& local_num_cols == mutable_data.get_num_cols());
		local_matrix_store::exposed_area orig_labels
			= mutable_labels.get_exposed_area();
		bool success = mutable_labels.resize(local_start_row, 0, local_num_rows, 1);
		if (!success)
			return false;

		success = mutable_data.resize(local_start_row, local_start_col, local_num_rows,
				local_num_cols);
		if (!success) {
			mutable_labels.restore_size(orig_labels);
			return false;
		}
		return local_matrix_store::resize(local_start_row, local_start_col,
				local_num_rows, local_num_cols);
	}
	virtual void reset_size() {
		mutable_data.reset_size();
		mutable_labels.reset_size();
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

matrix_store::const_ptr groupby_matrix_store::transpose() const
{
	matrix_margin new_margin;
	if (margin == matrix_margin::MAR_ROW)
		new_margin = matrix_margin::MAR_COL;
	else
		new_margin = matrix_margin::MAR_ROW;
	return matrix_store::const_ptr(new groupby_matrix_store(data->transpose(),
				label_store->transpose(), f, new_margin, agg_op));
}

class groupby_compute_store: public sink_compute_store, public EM_object
{
	std::shared_ptr<groupby_op> portion_op;
	matrix_store::const_ptr data;
	matrix_store::const_ptr label_store;
	matrix_margin margin;
public:
	groupby_compute_store(matrix_store::const_ptr data,
			matrix_store::const_ptr label_store,
			std::shared_ptr<portion_mapply_op> portion_op,
			matrix_margin margin): sink_compute_store(data->get_num_rows(),
				data->get_num_cols(),
				data->is_in_mem() && label_store->is_in_mem(),
				portion_op->get_output_type()) {
		this->data = data;
		this->label_store = label_store;
		this->portion_op = std::dynamic_pointer_cast<groupby_op>(portion_op);
		assert(this->portion_op);
		this->margin = margin;
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
		local_matrix_store::const_ptr data, local_matrix_store::const_ptr labels,
		const scalar_type &type, portion_mapply_op::const_ptr portion_op)
{
	if (data->store_layout() == matrix_layout_t::L_ROW)
		return local_matrix_store::const_ptr(new lmaterialize_row_matrix_store(
					data, labels, type, portion_op));
	else
		return local_matrix_store::const_ptr(new lmaterialize_col_matrix_store(
					data, labels, type, portion_op));
}

local_matrix_store::const_ptr groupby_compute_store::get_portion(
		size_t start_row, size_t start_col, size_t num_rows,
		size_t num_cols) const
{
	local_matrix_store::const_ptr data_part = data->get_portion(start_row,
			start_col, num_rows, num_cols);
	local_matrix_store::const_ptr label_part;
	// `label_store' is a col_vec.
	if (margin == matrix_margin::MAR_ROW) {
		assert(label_store->get_num_cols() == 1);
		label_part = label_store->get_portion(start_row, 0, num_rows, 1);
	}
	else {
		assert(label_store->get_num_rows() == 1);
		label_part = label_store->get_portion(0, start_col, 1, num_cols);
	}
	return create_lmaterialize_matrix(data_part, label_part, get_type(),
			portion_op);
}

local_matrix_store::const_ptr groupby_compute_store::get_portion(size_t id) const
{
	local_matrix_store::const_ptr data_part = data->get_portion(id);
	local_matrix_store::const_ptr label_part = label_store->get_portion(id);
	return create_lmaterialize_matrix(data_part, label_part, get_type(),
			portion_op);
}

namespace
{

class collect_portion_compute: public portion_compute
{
	size_t num_EM_parts;
	size_t num_reads;
	portion_compute::ptr orig_compute;
public:
	typedef std::shared_ptr<collect_portion_compute> ptr;

	collect_portion_compute(portion_compute::ptr orig_compute) {
		this->num_EM_parts = 0;
		this->num_reads = 0;
		this->orig_compute = orig_compute;
	}

	void set_EM_count(size_t num_EM_parts) {
		this->num_EM_parts = num_EM_parts;
	}

	virtual void run(char *buf, size_t size) {
		num_reads++;
		if (num_reads == num_EM_parts) {
			orig_compute->run(NULL, 0);
			// This only runs once.
			// Let's remove all user's portion compute to indicate that it has
			// been invoked.
			orig_compute = NULL;
		}
	}
};

}

async_cres_t groupby_compute_store::get_portion_async(
		size_t start_row, size_t start_col, size_t num_rows,
		size_t num_cols, std::shared_ptr<portion_compute> compute) const
{
	collect_portion_compute *_compute = new collect_portion_compute(compute);
	portion_compute::ptr collect_compute(_compute);

	async_cres_t ret = data->get_portion_async(start_row, start_col,
			num_rows, num_cols, collect_compute);
	local_matrix_store::const_ptr label_part;
	async_cres_t label_ret;
	// `label_store' is a col_vec.
	if (margin == matrix_margin::MAR_ROW) {
		assert(label_store->get_num_cols() == 1);
		label_ret = label_store->get_portion_async(start_row, 0, num_rows, 1,
				collect_compute);
	}
	else {
		assert(label_store->get_num_rows() == 1);
		label_ret = label_store->get_portion_async(0, start_col, 1, num_cols,
				collect_compute);
	}
	// It's two if both are unavailable.
	_compute->set_EM_count(!ret.first + !label_ret.first);
	return async_cres_t(ret.first && label_ret.first,
			create_lmaterialize_matrix(ret.second, label_ret.second,
			get_type(), portion_op));
}

std::vector<safs::io_interface::ptr> groupby_compute_store::create_ios() const
{
	std::vector<safs::io_interface::ptr> ios;
	const EM_object *obj = dynamic_cast<const EM_object *>(data.get());
	if (obj) {
		auto tmp = obj->create_ios();
		ios.insert(ios.end(), tmp.begin(), tmp.end());
	}

	obj = dynamic_cast<const EM_object *>(label_store.get());
	if (obj) {
		auto tmp = obj->create_ios();
		ios.insert(ios.end(), tmp.begin(), tmp.end());
	}
	return ios;
}

std::vector<virtual_matrix_store::const_ptr> groupby_matrix_store::get_compute_matrices() const
{
	// If the groupby matrix has been materialized, we don't need to do
	// anything.
	if (has_materialized())
		return std::vector<virtual_matrix_store::const_ptr>();
	else
		return std::vector<virtual_matrix_store::const_ptr>(1,
				virtual_matrix_store::const_ptr(new groupby_compute_store(data,
						label_store, portion_op, margin)));
}

}

}
