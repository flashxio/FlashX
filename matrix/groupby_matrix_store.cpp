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
	// Each row of a local matrix contains partially aggregated data for a label.
	std::vector<detail::local_row_matrix_store::ptr> part_results;
	std::vector<bool> part_status;
	size_t num_levels;
	matrix_margin margin;
	agg_operate::const_ptr op;
public:
	groupby_op(agg_operate::const_ptr op, size_t num_levels,
			matrix_margin margin): detail::portion_mapply_op(0, 0,
				op->get_output_type()) {
		size_t num_threads = detail::mem_thread_pool::get_global_num_threads();
		part_results.resize(num_threads);
		part_agg.resize(num_threads);
		part_status.resize(num_threads, true);
		this->num_levels = num_levels;
		this->margin = margin;
		this->op = op;
	}

	virtual detail::portion_mapply_op::const_ptr transpose() const {
		throw unsupported_exception("Don't support transpose of groupby_op");
	}

	virtual void run(
			const std::vector<detail::local_matrix_store::const_ptr> &ins) const;

	virtual std::string to_string(
			const std::vector<detail::matrix_store::const_ptr> &mats) const {
		throw unsupported_exception("Don't support to_string of groupby_op");
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
	detail::mem_matrix_store::ptr res = detail::mem_matrix_store::create(nrow,
			ncol, matrix_layout_t::L_ROW, type, -1);
	for (size_t i = 0; i < res->get_num_rows(); i++) {
		memcpy(res->get_row(i), part_results[first_idx]->get_row(i),
				res->get_num_cols() * res->get_entry_size());
		for (size_t j = first_idx + 1; j < part_results.size(); j++) {
			if (part_results[j] != NULL)
				op->get_combine().runAA(res->get_num_cols(),
						part_results[j]->get_row(i), res->get_row(i),
						res->get_row(i));
		}
	}
	return res;
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
		mutable_this->part_results[thread_id] = detail::local_row_matrix_store::ptr(
				new detail::local_buf_row_matrix_store(0, 0, num_levels,
					in->get_num_cols(), op->get_output_type(), -1));
		mutable_this->part_agg[thread_id].resize(num_levels);
	}
	// If there was a failure in this thread, we don't need to perform more
	// computation.
	if (!part_status[thread_id])
		return;

	assert(in->store_layout() == matrix_layout_t::L_ROW);
	bool ret = detail::groupby_row(*labels,
			static_cast<const detail::local_row_matrix_store &>(*in),
			*op, detail::part_dim_t::PART_DIM1, *part_results[thread_id],
			mutable_this->part_agg[thread_id]);
	if (!ret)
		mutable_this->part_status[thread_id] = false;
}

}

static size_t get_num_rows_groupby(const factor_col_vector &labels,
		matrix_margin margin)
{
	return margin == matrix_margin::MAR_ROW ? labels.get_length() : 1;
}

static size_t get_num_cols_groupby(const factor_col_vector &labels,
		matrix_margin margin)
{
	return margin == matrix_margin::MAR_COL ? labels.get_length() : 1;
}

groupby_matrix_store::groupby_matrix_store(matrix_store::const_ptr data,
		factor_col_vector::const_ptr labels, matrix_margin margin,
		agg_operate::const_ptr op): sink_store(get_num_rows_groupby(*labels, margin),
			get_num_cols_groupby(*labels, margin), data->is_in_mem(),
			op->get_output_type())
{
	this->data = data;
	this->label_store = labels->get_raw_store();
	this->margin = margin;
	portion_op = std::shared_ptr<groupby_op>(new groupby_op(op,
				labels->get_factor().get_num_levels(), margin));
}

matrix_store::ptr groupby_matrix_store::get_agg_res() const
{
	return std::static_pointer_cast<groupby_op>(portion_op)->get_agg();
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
		mutable_labels.resize(local_start_col, 0, local_num_cols, 1);
		mutable_data.resize(local_start_row, local_start_col, local_num_rows,
				local_num_cols);
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
		mutable_labels.resize(local_start_row, 0, local_num_rows, 1);
		mutable_data.resize(local_start_row, local_start_col, local_num_rows,
				local_num_cols);
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
	// TODO This method should also be implemented.
	assert(0);
	return matrix_store::const_ptr();
}

class groupby_compute_store: public sink_compute_store, public EM_object
{
	std::shared_ptr<portion_mapply_op> portion_op;
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
		this->portion_op = portion_op;
		this->margin = margin;
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
	assert(label_store->get_num_cols() == 1);
	// `label_store' is a col_vec.
	if (margin == matrix_margin::MAR_ROW)
		label_part = label_store->get_portion(start_row, 0, num_rows, 1);
	else
		label_part = label_store->get_portion(start_col, 0, num_cols, 1);
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

async_cres_t groupby_compute_store::get_portion_async(
		size_t start_row, size_t start_col, size_t num_rows,
		size_t num_cols, std::shared_ptr<portion_compute> compute) const
{
	async_cres_t ret = data->get_portion_async(start_row, start_col,
			num_rows, num_cols, compute);
	local_matrix_store::const_ptr label_part;
	assert(label_store->get_num_cols() == 1);
	assert(label_store->is_in_mem());
	// `label_store' is a col_vec.
	if (margin == matrix_margin::MAR_ROW)
		label_part = label_store->get_portion(start_row, 0, num_rows, 1);
	else
		label_part = label_store->get_portion(start_col, 0, num_cols, 1);
	ret.second = create_lmaterialize_matrix(ret.second, label_part,
			get_type(), portion_op);
	return ret;
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

virtual_matrix_store::const_ptr groupby_matrix_store::get_compute_matrix() const
{
	return virtual_matrix_store::const_ptr(new groupby_compute_store(data,
				label_store, portion_op, margin));
}

}

}
