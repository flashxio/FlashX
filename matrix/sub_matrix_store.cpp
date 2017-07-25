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

#include "EM_dense_matrix.h"
#include "sub_matrix_store.h"
#include "local_mem_buffer.h"

namespace fm
{

namespace detail
{

namespace
{

class collect_rc_compute: public portion_compute
{
	std::vector<portion_compute::ptr> computes;
	local_matrix_store::ptr collected;
	std::vector<local_matrix_store::const_ptr> orig_portions;
	size_t num_collects;
public:
	typedef std::shared_ptr<collect_rc_compute> ptr;

	collect_rc_compute(portion_compute::ptr compute,
			local_matrix_store::ptr collected) {
		computes.push_back(compute);
		this->collected = collected;
		num_collects = 0;
	}

	void add_orig_compute(portion_compute::ptr compute) {
		computes.push_back(compute);
	}

	void set_bufs(const std::vector<local_matrix_store::const_ptr> &orig) {
		this->orig_portions = orig;
	}
	virtual void run(char *buf, size_t size);
	void run_complete();
};

void collect_rc_compute::run_complete()
{
	if (collected->store_layout() == matrix_layout_t::L_COL) {
		size_t num_cols = 0;
		for (size_t i = 0; i < orig_portions.size(); i++) {
			collected->resize(0, num_cols, collected->get_num_rows(),
					orig_portions[i]->get_num_cols());
			collected->copy_from(*orig_portions[i]);
			num_cols += orig_portions[i]->get_num_cols();
		}
	}
	else {
		size_t num_rows = 0;
		for (size_t i = 0; i < orig_portions.size(); i++) {
			collected->resize(num_rows, 0, orig_portions[i]->get_num_rows(),
					collected->get_num_cols());
			collected->copy_from(*orig_portions[i]);
			num_rows += orig_portions[i]->get_num_rows();
		}
	}

	collected->reset_size();
	// I can't pass the raw array to the compute run now.
	// TODO Maybe I should pass the portion object to the run.
	for (size_t i = 0; i < computes.size(); i++)
		computes[i]->run(NULL, 0);
}

void collect_rc_compute::run(char *buf, size_t size)
{
	num_collects++;
	// After we collect all columns.
	if (num_collects == orig_portions.size())
		run_complete();
}

class local_collected_buf_col_matrix_store: public local_buf_col_matrix_store
{
	std::weak_ptr<collect_rc_compute> compute;
public:
	typedef std::shared_ptr<local_collected_buf_col_matrix_store> ptr;

	local_collected_buf_col_matrix_store(off_t global_start_row,
			off_t global_start_col, size_t nrow, size_t ncol,
			const scalar_type &type, int node_id): local_buf_col_matrix_store(
				global_start_row, global_start_col, nrow, ncol, type, node_id) {
	}

	void set_compute(collect_rc_compute::ptr compute) {
		this->compute = compute;
	}

	collect_rc_compute::ptr get_compute() const {
		return compute.lock();
	}
};

class local_collected_buf_row_matrix_store: public local_buf_row_matrix_store
{
	std::weak_ptr<collect_rc_compute> compute;
public:
	typedef std::shared_ptr<local_collected_buf_row_matrix_store> ptr;

	local_collected_buf_row_matrix_store(off_t global_start_row,
			off_t global_start_col, size_t nrow, size_t ncol,
			const scalar_type &type, int node_id): local_buf_row_matrix_store(
				global_start_row, global_start_col, nrow, ncol, type, node_id) {
	}

	void set_compute(collect_rc_compute::ptr compute) {
		this->compute = compute;
	}

	collect_rc_compute::ptr get_compute() const {
		return compute.lock();
	}
};

}

static std::vector<std::pair<off_t, off_t> > split_idxs(
		std::vector<off_t>::const_iterator it,
		std::vector<off_t>::const_iterator end)
{
	assert(it != end);
	std::vector<std::pair<off_t, off_t> > vecs(1);
	vecs[0].first = *it;
	vecs[0].second = *it + 1;
	for (it++; it != end; it++) {
		// It's contiguous.
		if (vecs.back().second == *it)
			vecs.back().second++;
		else {
			std::pair<off_t, off_t> new_range(*it, *it + 1);
			vecs.push_back(new_range);
		}
	}
	return vecs;
}

matrix_store::const_ptr sub_col_matrix_store::transpose() const
{
	matrix_store::const_ptr t_mat = get_orig().transpose();
	return matrix_store::const_ptr(new sub_row_matrix_store(rc_idxs,
				t_mat, data_id));
}

async_cres_t sub_col_matrix_store::get_mem_portion_async(
		size_t start_row, size_t start_col, size_t num_rows,
		size_t num_cols, portion_compute::ptr compute) const
{
	// Get the entire portion from the original matrix.
	local_matrix_store::const_ptr portion
		= get_orig().get_portion(start_row, 0, num_rows,
				get_orig().get_num_cols());
	assert(portion->store_layout() == matrix_layout_t::L_COL);
	local_col_matrix_store::const_ptr col_portion
		= std::static_pointer_cast<const local_col_matrix_store>(portion);
	std::vector<const char *> cols(num_cols);
	assert(cols.size() <= rc_idxs.size());
	for (size_t i = 0; i < num_cols; i++)
		cols[i] = col_portion->get_col(rc_idxs[i]);
	// TODO we might be able to return a local matrix with contiguous memory.
	return async_cres_t(true, local_matrix_store::const_ptr(
				new local_cref_col_matrix_store(
					cols, start_row, start_col, num_rows, num_cols,
					get_type(), portion->get_node_id())));
}

async_cres_t sub_col_matrix_store::get_EM_portion_async(
		size_t start_row, size_t start_col, size_t num_rows,
		size_t num_cols, portion_compute::ptr compute) const
{
	auto wanted_it = rc_idxs.begin() + start_col;
	auto wanted_end = wanted_it + num_cols;
	std::vector<std::pair<off_t, off_t> > ranges = split_idxs(wanted_it,
			wanted_end);

	size_t start_fetched_row;
	size_t num_fetched_rows;
	if (is_wide()) {
		start_fetched_row = start_row;
		num_fetched_rows = num_rows;
	}
	else {
		auto chunk_size
			= is_wide() ? get_portion_size().second : get_portion_size().first;
		start_fetched_row = start_row - (start_row % chunk_size);
		num_fetched_rows = std::min(chunk_size,
				get_num_rows() - start_fetched_row);
	}

	// Fetch the portion from the cache.
	local_matrix_store::const_ptr ret1 = local_mem_buffer::get_mat_portion(
			data_id->get_id());

	bool match = false;
	bool match_trans = false;
	if (ret1) {
		// If it's in the same portion.
		match = (size_t) ret1->get_global_start_row() == start_fetched_row
			&& (size_t) ret1->get_global_start_col() == start_col
			&& ret1->get_num_rows() == num_fetched_rows
			&& ret1->get_num_cols() == num_cols;
		// If it's in the corresponding portion in the transposed matrix.
		match_trans = (size_t) ret1->get_global_start_row() == start_col
			&& (size_t) ret1->get_global_start_col() == start_fetched_row
			&& ret1->get_num_rows() == num_cols
			&& ret1->get_num_cols() == num_fetched_rows;
	}

	if (match || match_trans) {
		assert(ret1->get_local_start_row() == 0);
		assert(ret1->get_local_start_col() == 0);
		collect_rc_compute::ptr collect_compute;
		// In the asynchronous version, data in the portion isn't ready when
		// the method is called. We should add the user's portion computation
		// to the queue. When the data is ready, all user's portion computations
		// will be invoked.
		if (match) {
			local_matrix_store *tmp = const_cast<local_matrix_store *>(ret1.get());
			local_collected_buf_col_matrix_store *store
				= dynamic_cast<local_collected_buf_col_matrix_store *>(tmp);
			assert(store);
			compute = store->get_compute();
		}
		else {
			local_matrix_store *tmp = const_cast<local_matrix_store *>(ret1.get());
			local_collected_buf_row_matrix_store *store
				= dynamic_cast<local_collected_buf_row_matrix_store *>(tmp);
			assert(store);
			compute = store->get_compute();
		}
		// If collect_rc_compute doesn't exist, it mean the data has been read
		// from disks.
		bool valid_data = collect_compute == NULL;
		if (!valid_data)
			collect_compute->add_orig_compute(compute);

		// We need its transpose.
		if (match_trans)
			ret1 = std::static_pointer_cast<const local_matrix_store>(
					ret1->transpose());

		local_matrix_store::const_ptr ret;
		if (start_fetched_row < start_row || num_rows < num_fetched_rows)
			ret = ret1->get_portion(start_row - start_fetched_row, 0,
					num_rows, num_cols);
		else
			ret = ret1;
		return async_cres_t(valid_data, ret);
	}

	local_collected_buf_col_matrix_store::ptr ret(
			new local_collected_buf_col_matrix_store(start_fetched_row,
				start_col, num_fetched_rows, num_cols, get_type(), -1));
	collect_rc_compute::ptr collect_compute(new collect_rc_compute(compute, ret));
	ret->set_compute(collect_compute);

	std::vector<local_matrix_store::const_ptr> bufs(ranges.size());
	size_t collected_cols = 0;
	size_t num_ready = 0;
	for (size_t i = 0; i < ranges.size(); i++) {
		size_t local_num_cols = ranges[i].second - ranges[i].first;
		collected_cols += local_num_cols;
		async_cres_t res = get_orig().get_portion_async(
				start_fetched_row, ranges[i].first, num_fetched_rows,
				local_num_cols, collect_compute);
		if (res.first)
			num_ready++;
		bufs[i] = res.second;
	}
	assert(collected_cols == num_cols);
	collect_compute->set_bufs(bufs);
	// Here we tell the collected compute how many of its buffers are ready.
	for (size_t i = 0; i < num_ready; i++)
		collect_compute->run(NULL, 0);
	if (get_orig().is_cache_portion())
		local_mem_buffer::cache_portion(data_id->get_id(), ret);

	bool ready = num_ready == bufs.size();
	if (num_fetched_rows == num_rows)
		return async_cres_t(ready, ret);
	else
		return async_cres_t(ready, ret->get_portion(
					start_row - start_fetched_row, 0, num_rows, num_cols));
}

matrix_store::const_ptr sub_row_matrix_store::transpose() const
{
	matrix_store::const_ptr t_mat = get_orig().transpose();
	return matrix_store::const_ptr(new sub_col_matrix_store(rc_idxs,
				t_mat, data_id));
}

async_cres_t sub_row_matrix_store::get_mem_portion_async(
		size_t start_row, size_t start_col, size_t num_rows,
		size_t num_cols, portion_compute::ptr compute) const
{
	// Get the entire portion from the original matrix.
	local_matrix_store::const_ptr portion = get_orig().get_portion(0,
			start_col, get_orig().get_num_rows(), num_cols);
	assert(portion->store_layout() == matrix_layout_t::L_ROW);
	local_row_matrix_store::const_ptr row_portion
		= std::static_pointer_cast<const local_row_matrix_store>(portion);
	std::vector<const char *> rows(num_rows);
	assert(rows.size() <= rc_idxs.size());
	for (size_t i = 0; i < num_rows; i++)
		rows[i] = row_portion->get_row(rc_idxs[i]);
	// TODO we might be able to return a local matrix with contiguous memory.
	return async_cres_t(true, local_matrix_store::const_ptr(
				new local_cref_row_matrix_store(
					rows, start_row, start_col, num_rows, num_cols,
					get_type(), portion->get_node_id())));
}

async_cres_t sub_row_matrix_store::get_EM_portion_async(
		size_t start_row, size_t start_col, size_t num_rows,
		size_t num_cols, portion_compute::ptr compute) const
{
	auto wanted_it = rc_idxs.begin() + start_row;
	auto wanted_end = wanted_it + num_rows;
	std::vector<std::pair<off_t, off_t> > ranges = split_idxs(wanted_it,
			wanted_end);

	size_t start_fetched_col;
	size_t num_fetched_cols;
	if (is_wide()) {
		auto chunk_size
			= is_wide() ? get_portion_size().second : get_portion_size().first;
		start_fetched_col = start_col - (start_col % chunk_size);
		num_fetched_cols = std::min(chunk_size,
				get_num_cols() - start_fetched_col);
	}
	else {
		start_fetched_col = start_col;
		num_fetched_cols = num_cols;
	}

	// Fetch the portion from the cache.
	local_matrix_store::const_ptr ret1 = local_mem_buffer::get_mat_portion(
			data_id->get_id());
	bool match = false;
	bool match_trans = false;
	if (ret1) {
		// If it's in the same portion.
		match = (size_t) ret1->get_global_start_row() == start_row
			&& (size_t) ret1->get_global_start_col() == start_fetched_col
			&& ret1->get_num_rows() == num_rows
			&& ret1->get_num_cols() == num_fetched_cols;
		// If it's in the corresponding portion in the transposed matrix.
		match_trans = (size_t) ret1->get_global_start_row() == start_fetched_col
			&& (size_t) ret1->get_global_start_col() == start_row
			&& ret1->get_num_rows() == num_fetched_cols
			&& ret1->get_num_cols() == num_rows;
	}
	if (match || match_trans) {
		assert(ret1->get_local_start_row() == 0);
		assert(ret1->get_local_start_col() == 0);
		// In the asynchronous version, data in the portion isn't ready when
		// the method is called. We should add the user's portion computation
		// to the queue. When the data is ready, all user's portion computations
		// will be invoked.
		collect_rc_compute::ptr collect_compute;
		if (match) {
			local_matrix_store *tmp = const_cast<local_matrix_store *>(ret1.get());
			local_collected_buf_row_matrix_store *store
				= dynamic_cast<local_collected_buf_row_matrix_store *>(tmp);
			assert(store);
			collect_compute = store->get_compute();
		}
		else {
			local_matrix_store *tmp = const_cast<local_matrix_store *>(ret1.get());
			local_collected_buf_col_matrix_store *store
				= dynamic_cast<local_collected_buf_col_matrix_store *>(tmp);
			assert(store);
			collect_compute = store->get_compute();
		}
		bool valid_data = collect_compute == NULL;
		if (!valid_data)
			collect_compute->add_orig_compute(compute);

		// We need its transpose.
		if (match_trans)
			ret1 = std::static_pointer_cast<const local_matrix_store>(
					ret1->transpose());

		local_matrix_store::const_ptr ret;
		if (start_fetched_col < start_col || num_cols < num_fetched_cols)
			ret = ret1->get_portion(0, start_col - start_fetched_col,
					num_rows, num_cols);
		else
			ret = ret1;
		return async_cres_t(valid_data, ret);
	}

	local_collected_buf_row_matrix_store::ptr ret(
			new local_collected_buf_row_matrix_store(start_row,
				start_fetched_col, num_rows, num_fetched_cols, get_type(), -1));
	collect_rc_compute::ptr collect_compute(new collect_rc_compute(compute, ret));
	ret->set_compute(collect_compute);

	std::vector<local_matrix_store::const_ptr> bufs(ranges.size());
	size_t collected_rows = 0;
	size_t num_ready = 0;
	for (size_t i = 0; i < ranges.size(); i++) {
		size_t local_num_rows = ranges[i].second - ranges[i].first;
		collected_rows += local_num_rows;
		async_cres_t res = get_orig().get_portion_async(
				ranges[i].first, start_fetched_col, local_num_rows,
				num_fetched_cols, collect_compute);
		if (res.first)
			num_ready++;
		bufs[i] = res.second;
	}
	assert(collected_rows == num_rows);
	collect_compute->set_bufs(bufs);
	// Here we tell the collected compute how many of its buffers are ready.
	for (size_t i = 0; i < num_ready; i++)
		collect_compute->run(NULL, 0);
	if (get_orig().is_cache_portion())
		local_mem_buffer::cache_portion(data_id->get_id(), ret);

	bool ready = num_ready == bufs.size();
	if (num_fetched_cols == num_cols)
		return async_cres_t(ready, ret);
	else
		return async_cres_t(ready, ret->get_portion(
					0, start_col - start_fetched_col, num_rows, num_cols));
}

local_matrix_store::const_ptr sub_matrix_store::get_portion(
		size_t start_row, size_t start_col, size_t num_rows, size_t num_cols) const
{
	// If it's in memory, we can get the portion directly.
	if (get_orig().is_in_mem()) {
		auto ret = get_portion_async(start_row, start_col,
				num_rows, num_cols, portion_compute::ptr());
//		printf("get portion %ld,%ld,%ld,%ld, node %d\n",
//				start_row, start_col, num_rows, num_cols,
//				ret.second->get_node_id());
		assert(ret.first);
		return ret.second;
	}

	bool ready = false;
	portion_compute::ptr compute(new sync_read_compute(ready));
	async_cres_t ret = get_portion_async(start_row, start_col,
			num_rows, num_cols, compute);
	// If we can't get the specified portion or the portion already has
	// the valid data.
	if (ret.second == NULL || ret.first)
		return ret.second;

	// TODO this is an ugly implementation.
	const EM_matrix_store *em_orig
		= dynamic_cast<const EM_matrix_store *>(orig.get());
	while (!ready)
		em_orig->wait4complete(1);
	return ret.second;
}

std::vector<safs::io_interface::ptr> sub_matrix_store::create_ios() const
{
	const EM_object *obj = dynamic_cast<const EM_object *>(orig.get());
	if (obj)
		return obj->create_ios();
	else
		return std::vector<safs::io_interface::ptr>();
}

async_cres_t sub_matrix_store::get_portion_async(
		size_t start_row, size_t start_col, size_t num_rows,
		size_t num_cols, portion_compute::ptr compute) const
{
	if (start_row + num_rows > get_num_rows()
			|| start_col + num_cols > get_num_cols()) {
		BOOST_LOG_TRIVIAL(error) << "get portion async: out of boundary";
		return async_cres_t();
	}

	if (!is_in_mem())
		return get_EM_portion_async(start_row, start_col,
				num_rows, num_cols, compute);
	else
		return get_mem_portion_async(start_row, start_col,
				num_rows, num_cols, compute);
}

}

}
