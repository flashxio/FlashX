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

#include <boost/foreach.hpp>

#include "io_interface.h"
#include "safs_file.h"

#include "matrix_config.h"
#include "EM_dense_matrix.h"
#include "local_matrix_store.h"
#include "raw_data_array.h"
#include "mem_matrix_store.h"
#include "matrix_stats.h"
#include "local_mem_buffer.h"

namespace fm
{

namespace detail
{

static const size_t CHUNK_SIZE = 1024 * 1024;

/*
 * These two functions define the length and portion size for 1D partitioning
 * on a matrix.
 */

static inline size_t get_tot_len(const matrix_store &mat)
{
	return mat.is_wide() ? mat.get_num_cols() : mat.get_num_rows();
}

static inline size_t get_portion_size(const matrix_store &mat)
{
	return mat.is_wide() ? mat.get_portion_size().second : mat.get_portion_size().first;
}

namespace
{

/*
 * When we write data to disks, we need to have something to hold the buffer.
 * This holds the local buffer until the write completes.
 */
class portion_write_complete: public portion_compute
{
	local_matrix_store::const_ptr store;
public:
	portion_write_complete(local_matrix_store::const_ptr store) {
		this->store = store;
	}

	virtual void run(char *buf, size_t size) {
	}
};

}

EM_matrix_store::EM_matrix_store(size_t nrow, size_t ncol, matrix_layout_t layout,
		const scalar_type &type): matrix_store(nrow, ncol, false,
			type), mat_id(mat_counter++), data_id(mat_id)
{
	this->layout = layout;
	holder = file_holder::create_temp("mat", nrow * ncol * type.get_size());
	safs::file_io_factory::shared_ptr factory = safs::create_io_factory(
			holder->get_name(), safs::REMOTE_ACCESS);
	ios = io_set::ptr(new io_set(factory));
}

EM_matrix_store::EM_matrix_store(file_holder::ptr holder, io_set::ptr ios,
		size_t nrow, size_t ncol, matrix_layout_t layout,
		const scalar_type &type, size_t _data_id): matrix_store(nrow, ncol,
			false, type), mat_id(mat_counter++), data_id(_data_id)
{
	this->layout = layout;
	this->holder = holder;
	this->ios = ios;
}

void EM_matrix_store::reset_data()
{
	assert(0);
}

namespace
{

class EM_mat_setdata_dispatcher: public EM_portion_dispatcher
{
	const set_operate &op;
	EM_matrix_store &to_mat;
public:
	EM_mat_setdata_dispatcher(EM_matrix_store &store, const set_operate &_op);

	virtual void create_task(off_t global_start, size_t length);
};

EM_mat_setdata_dispatcher::EM_mat_setdata_dispatcher(EM_matrix_store &store,
		const set_operate &_op): EM_portion_dispatcher(get_tot_len(store),
			fm::detail::get_portion_size(store)), op(_op), to_mat(store)
{
}

void EM_mat_setdata_dispatcher::create_task(off_t global_start,
		size_t length)
{
	size_t global_start_row, global_start_col;
	size_t num_rows, num_cols;
	if (to_mat.is_wide()) {
		global_start_row = 0;
		global_start_col = global_start;
		num_rows = to_mat.get_num_rows();
		num_cols = length;
	}
	else {
		global_start_row = global_start;
		global_start_col = 0;
		num_rows = length;
		num_cols = to_mat.get_num_cols();
	}
	local_matrix_store::ptr buf;
	if (to_mat.store_layout() == matrix_layout_t::L_COL)
		buf = local_matrix_store::ptr(new local_buf_col_matrix_store(
					global_start_row, global_start_col, num_rows, num_cols,
					to_mat.get_type(), -1));
	else
		buf = local_matrix_store::ptr(new local_buf_row_matrix_store(
					global_start_row, global_start_col, num_rows, num_cols,
					to_mat.get_type(), -1));
	buf->set_data(op);
	to_mat.write_portion_async(buf, global_start_row, global_start_col);
}

}

void EM_matrix_store::set_data(const set_operate &op)
{
	mem_thread_pool::ptr threads = mem_thread_pool::get_global_mem_threads();
	EM_mat_setdata_dispatcher::ptr dispatcher(
			new EM_mat_setdata_dispatcher(*this, op));
	for (size_t i = 0; i < threads->get_num_threads(); i++) {
		io_worker_task *task = new io_worker_task(dispatcher);
		task->register_EM_obj(this);
		threads->process_task(i % threads->get_num_nodes(), task);
	}
	threads->wait4complete();
}

std::pair<size_t, size_t> EM_matrix_store::get_portion_size() const
{
	return std::pair<size_t, size_t>(std::min(get_num_rows(), CHUNK_SIZE),
			std::min(get_num_cols(), CHUNK_SIZE));
}

local_matrix_store::ptr EM_matrix_store::get_portion(
		size_t start_row, size_t start_col, size_t num_rows,
		size_t num_cols)
{
	// This doesn't need to be used. Changing the data in the local portion
	// doesn't affect the data in the disks.
	assert(0);
	return local_matrix_store::ptr();
}

local_matrix_store::const_ptr EM_matrix_store::get_portion(
		size_t start_row, size_t start_col, size_t num_rows,
		size_t num_cols) const
{
	safs::io_interface &io = ios->get_curr_io();
	bool ready = false;
	portion_compute::ptr compute(new sync_read_compute(ready));
	local_matrix_store::const_ptr ret = get_portion_async(start_row, start_col,
			num_rows, num_cols, compute);
	if (ret == NULL)
		return ret;

	while (!ready)
		io.wait4complete(1);
	return ret;
}

local_matrix_store::ptr EM_matrix_store::get_portion_async(
		size_t start_row, size_t start_col, size_t num_rows,
		size_t num_cols, portion_compute::ptr compute)
{
	// This doesn't need to be used. Changing the data in the local portion
	// doesn't affect the data in the disks.
	assert(0);
	return local_matrix_store::ptr();
}

local_matrix_store::const_ptr EM_matrix_store::get_portion_async(
		size_t start_row, size_t start_col, size_t num_rows,
		size_t num_cols, portion_compute::ptr compute) const
{
	size_t local_start_row = start_row % CHUNK_SIZE;
	size_t local_start_col = start_col % CHUNK_SIZE;
	// For now, we only allow this method to fetch data from a portion.
	if (local_start_row + num_rows > CHUNK_SIZE
			|| local_start_col + num_cols > CHUNK_SIZE
			|| start_row + num_rows > get_num_rows()
			|| start_col + num_cols > get_num_cols()) {
		BOOST_LOG_TRIVIAL(error) << "Out of boundary of a portion";
		return local_matrix_store::ptr();
	}

	safs::io_interface &io = ios->get_curr_io();
	size_t entry_size = get_type().get_size();
	size_t portion_start_row = start_row - local_start_row;
	size_t portion_start_col = start_col - local_start_col;
	size_t portion_num_rows;
	size_t portion_num_cols;
	if (portion_start_row + CHUNK_SIZE > get_num_rows())
		portion_num_rows = get_num_rows() - portion_start_row;
	else
		portion_num_rows = CHUNK_SIZE;
	if (portion_start_col + CHUNK_SIZE > get_num_cols())
		portion_num_cols = get_num_cols() - portion_start_col;
	else
		portion_num_cols = CHUNK_SIZE;

	// Location of the portion on the disks.
	off_t off = (get_num_cols() * portion_start_row
		+ portion_num_rows * portion_start_col) * entry_size;
	// If this is the very last portion (the bottom right portion), the data
	// size may not be aligned with the page size.
	size_t num_bytes
		= ROUNDUP(portion_num_rows * portion_num_cols * entry_size, PAGE_SIZE);

	// We should try to get the portion from the local thread memory buffer
	// first.
	local_matrix_store::const_ptr ret1 = local_mem_buffer::get_mat_portion(
			data_id);
	// If it's in the same portion.
	if (ret1 && (((size_t) ret1->get_global_start_row() == start_row
					&& (size_t) ret1->get_global_start_col() == start_col
					&& ret1->get_num_rows() == num_rows
					&& ret1->get_num_cols() == num_cols)
				// If it's in the corresponding portion in the transposed matrix.
				|| ((size_t) ret1->get_global_start_row() == start_col
					&& (size_t) ret1->get_global_start_col() == start_row
					&& ret1->get_num_rows() == num_cols
					&& ret1->get_num_cols() == num_rows))) {
		assert(ret1->get_local_start_row() == 0);
		assert(ret1->get_local_start_col() == 0);
		// In the asynchronous version, data in the portion isn't ready when
		// the method is called. We should add the user's portion computation
		// to the queue. When the data is ready, all user's portion computations
		// will be invoked.
		safs::data_loc_t loc(io.get_file_id(), off);
		safs::io_request req(const_cast<char *>(ret1->get_raw_arr()), loc,
				num_bytes, READ);
		portion_callback &cb = static_cast<portion_callback &>(io.get_callback());
		// If there isn't a portion compute related to the I/O request,
		// it means the data in the portion is ready. TODO should we invoke
		// user's portion compute directly?
		assert(cb.has_callback(req));
		cb.add(req, compute);
		return ret1;
	}

	raw_data_array data_arr(num_bytes, -1);
	// Read the portion in a single I/O request.
	local_matrix_store::ptr buf;
	if (store_layout() == matrix_layout_t::L_ROW)
		buf = local_matrix_store::ptr(new local_buf_row_matrix_store(data_arr,
					portion_start_row, portion_start_col, portion_num_rows,
					portion_num_cols, get_type(), data_arr.get_node_id()));
	else
		buf = local_matrix_store::ptr(new local_buf_col_matrix_store(data_arr,
					portion_start_row, portion_start_col, portion_num_rows,
					portion_num_cols, get_type(), data_arr.get_node_id()));

	safs::data_loc_t loc(io.get_file_id(), off);
	safs::io_request req(buf->get_raw_arr(), loc, num_bytes, READ);
	static_cast<portion_callback &>(io.get_callback()).add(req, compute);
	io.access(&req, 1);
	io.flush_requests();

	// Count the number of bytes really read from disks.
	detail::matrix_stats.inc_read_bytes(
			num_rows * num_cols * get_entry_size(), false);

	if (local_start_row > 0 || local_start_col > 0
			|| num_rows < portion_num_rows || num_cols < portion_num_cols)
		buf->resize(local_start_row, local_start_col, num_rows, num_cols);
	local_mem_buffer::cache_portion(data_id, buf);
	return buf;
}

void EM_matrix_store::write_portion_async(
		local_matrix_store::const_ptr portion, off_t start_row,
		off_t start_col)
{
	assert(store_layout() == portion->store_layout());
	assert(start_row % CHUNK_SIZE == 0 && start_col % CHUNK_SIZE == 0);
	assert(start_row + portion->get_num_rows() <= get_num_rows()
			&& start_col + portion->get_num_cols() <= get_num_cols());
	// Make sure the portion is stored contiguously on disks.
	// Otherwise, we need to read data from the disk first, modify it and
	// write it back.
	if (portion->get_global_start_row() + CHUNK_SIZE > get_num_rows())
		assert(portion->get_num_rows()
				== get_num_rows() - portion->get_global_start_row());
	else
		assert(portion->get_num_rows() == CHUNK_SIZE);
	if (portion->get_global_start_col() + CHUNK_SIZE > get_num_cols())
		assert(portion->get_num_cols()
				== get_num_cols() - portion->get_global_start_col());
	else
		assert(portion->get_num_cols() == CHUNK_SIZE);

	// And data in memory is also stored contiguously.
	// This constraint can be relaxed in the future.
	assert(portion->get_raw_arr());

	safs::io_interface &io = ios->get_curr_io();
	size_t entry_size = get_type().get_size();

	// Location of the portion on the disks.
	// I need to compute the location here because I might need to change
	// the size of a portion later.
	off_t off = (get_num_cols() * portion->get_global_start_row()
		+ portion->get_num_rows() * portion->get_global_start_col()) * entry_size;

	size_t num_bytes
		= portion->get_num_rows() * portion->get_num_cols() * entry_size;
	// If this is the very last portion (the bottom right portion), the data
	// size may not be aligned with the page size.
	if (num_bytes % PAGE_SIZE != 0) {
		raw_data_array data_arr(ROUNDUP(num_bytes, PAGE_SIZE),
				portion->get_node_id());
		// if the data layout is row wise, we should align the number of
		// rows, so the rows are still stored contiguously.
		if (store_layout() == matrix_layout_t::L_ROW) {
			local_buf_row_matrix_store::ptr tmp_buf(new local_buf_row_matrix_store(
						data_arr, portion->get_global_start_row(),
						portion->get_global_start_col(),
						portion->get_num_rows(), portion->get_num_cols(),
						portion->get_type(), portion->get_node_id()));
			memcpy(tmp_buf->get_raw_arr(), portion->get_raw_arr(), num_bytes);
			portion = tmp_buf;
		}
		// if the data layout is column wise, we should align the number of
		// columns, so the columns are still stored contiguously.
		else {
			local_buf_col_matrix_store::ptr tmp_buf(new local_buf_col_matrix_store(
						portion->get_global_start_row(),
						portion->get_global_start_col(),
						portion->get_num_rows(), portion->get_num_cols(),
						portion->get_type(), portion->get_node_id()));
			memcpy(tmp_buf->get_raw_arr(), portion->get_raw_arr(), num_bytes);
			portion = tmp_buf;
		}
		num_bytes = ROUNDUP(num_bytes, PAGE_SIZE);
	}

	safs::data_loc_t loc(io.get_file_id(), off);
	safs::io_request req(const_cast<char *>(portion->get_raw_arr()),
			loc, num_bytes, WRITE);
	portion_compute::ptr compute(new portion_write_complete(portion));
	static_cast<portion_callback &>(io.get_callback()).add(req, compute);
	io.access(&req, 1);
	io.flush_requests();
}

matrix_store::const_ptr EM_matrix_store::append_cols(
		const std::vector<matrix_store::const_ptr> &mats) const
{
	// TODO
	assert(0);
}

std::vector<safs::io_interface::ptr> EM_matrix_store::create_ios() const
{
	std::vector<safs::io_interface::ptr> ret(1);
	ret[0] = ios->create_io();
	return ret;
}

vec_store::const_ptr EM_matrix_store::get_col_vec(off_t idx) const
{
	if ((size_t) idx >= get_num_cols()) {
		BOOST_LOG_TRIVIAL(error) << "Out of boundary";
		return vec_store::const_ptr();
	}

	size_t entry_size = get_type().get_size();
	safs::io_interface &io = ios->get_curr_io();
	if (get_num_cols() == 1) {
		size_t len = roundup_ele(get_num_rows(), PAGE_SIZE, entry_size);
		smp_vec_store::ptr vec = smp_vec_store::create(len, get_type());
		vec->expose_sub_vec(0, get_num_rows());

		safs::data_loc_t loc(io.get_file_id(), 0);
		safs::io_request req(vec->get_raw_arr(), loc, len * entry_size, READ);
		io.access(&req, 1);
		io.wait4complete(1);
		return vec;
	}
	else if (store_layout() == matrix_layout_t::L_COL) {
		smp_vec_store::ptr vec = smp_vec_store::create(get_num_rows(),
				get_type());

		// Read all portions of data with the full length.
		std::vector<safs::io_request> reqs;
		off_t row_idx = 0;
		for (; row_idx + CHUNK_SIZE < get_num_rows(); row_idx += CHUNK_SIZE) {
			off_t off = (get_num_cols() * row_idx
					+ CHUNK_SIZE * idx) * entry_size;

			safs::data_loc_t loc(io.get_file_id(), off);
			safs::io_request req(vec->get_raw_arr() + row_idx * entry_size, loc,
					CHUNK_SIZE * entry_size, READ);
			reqs.push_back(req);
		}

		// Read the data in the last portion. It's possible that the offset
		// and the length of the data isn't aligned with the page size.
		size_t last_nrows = get_num_rows() - row_idx;
		assert(last_nrows <= CHUNK_SIZE);
		off_t ele_start = (get_num_cols() * row_idx + last_nrows * idx);
		off_t read_start = round_ele(ele_start, PAGE_SIZE, entry_size);
		size_t num_read_eles = roundup_ele(last_nrows, PAGE_SIZE, entry_size);
		assert(read_start <= ele_start
				&& ele_start + last_nrows <= read_start + num_read_eles);
		smp_vec_store::ptr tmp = smp_vec_store::create(num_read_eles, get_type());
		safs::data_loc_t loc(io.get_file_id(), read_start * entry_size);
		safs::io_request req(tmp->get_raw_arr(), loc,
				num_read_eles * entry_size, READ);
		reqs.push_back(req);

		io.access(reqs.data(), reqs.size());
		io.wait4complete(reqs.size());
		memcpy(vec->get_raw_arr() + row_idx * entry_size,
				tmp->get_raw_arr() + (ele_start - read_start) * entry_size,
				last_nrows * entry_size);
		return vec;
	}
	else {
		return vec_store::const_ptr();
	}
}

vec_store::const_ptr EM_matrix_store::get_row_vec(off_t idx) const
{
	if ((size_t) idx >= get_num_rows()) {
		BOOST_LOG_TRIVIAL(error) << "Out of boundary";
		return vec_store::const_ptr();
	}
	return transpose()->get_col_vec(idx);
}

namespace
{

class EM_mat_load_dispatcher: public detail::EM_portion_dispatcher
{
	const detail::EM_matrix_store &from_mat;
	detail::mem_matrix_store::ptr to_mat;
public:
	EM_mat_load_dispatcher(const detail::EM_matrix_store &_from_mat,
			detail::mem_matrix_store::ptr to_mat): detail::EM_portion_dispatcher(
				get_tot_len(_from_mat), fm::detail::get_portion_size(_from_mat)),
			from_mat(_from_mat) {
		this->to_mat = to_mat;
	}

	virtual void create_task(off_t global_start, size_t length);
};

class load_portion_compute: public portion_compute
{
	local_matrix_store::const_ptr from_store;
	local_matrix_store::ptr to_store;
public:
	void set_buf(local_matrix_store::const_ptr from_store,
			local_matrix_store::ptr to_store) {
		this->from_store = from_store;
		this->to_store = to_store;
	}

	virtual void run(char *buf, size_t size) {
		assert(from_store);
		assert(to_store);
		to_store->copy_from(*from_store);
	}
};

void EM_mat_load_dispatcher::create_task(off_t global_start, size_t length)
{
	load_portion_compute *load_compute = new load_portion_compute();
	load_portion_compute::ptr compute(load_compute);
	size_t global_start_row;
	size_t global_start_col;
	size_t num_rows;
	size_t num_cols;
	if (from_mat.is_wide()) {
		global_start_row = 0;
		global_start_col = global_start;
		num_rows = from_mat.get_num_rows();
		num_cols = length;
	}
	else {
		global_start_row = global_start;
		global_start_col = 0;
		num_rows = length;
		num_cols = from_mat.get_num_cols();
	}
	detail::local_matrix_store::const_ptr store1 = from_mat.get_portion_async(
			global_start_row, global_start_col, num_rows, num_cols, compute);
	detail::local_matrix_store::ptr store2 = to_mat->get_portion(
			global_start_row, global_start_col, num_rows, num_cols);
	load_compute->set_buf(store1, store2);
}

}

mem_matrix_store::ptr EM_matrix_store::load() const
{
	mem_matrix_store::ptr ret = mem_matrix_store::create(get_num_rows(),
			get_num_cols(), store_layout(), get_type(), -1);

	EM_mat_load_dispatcher::ptr dispatcher(
			new EM_mat_load_dispatcher(*this, ret));
	detail::mem_thread_pool::ptr threads
		= detail::mem_thread_pool::get_global_mem_threads();
	for (size_t i = 0; i < threads->get_num_threads(); i++) {
		detail::io_worker_task *task = new detail::io_worker_task(dispatcher);
		const detail::EM_object *obj = this;
		task->register_EM_obj(const_cast<detail::EM_object *>(obj));
		threads->process_task(i % threads->get_num_nodes(), task);
	}
	threads->wait4complete();
	return ret;
}

namespace
{

class EM_mat_copy_dispatcher: public EM_portion_dispatcher
{
	const detail::matrix_store &from_mat;
	EM_matrix_store &to_mat;
public:
	EM_mat_copy_dispatcher(const matrix_store &from, EM_matrix_store &to);

	virtual void create_task(off_t global_start, size_t length);
};

EM_mat_copy_dispatcher::EM_mat_copy_dispatcher(const matrix_store &from,
		EM_matrix_store &to): EM_portion_dispatcher(get_tot_len(to),
			fm::detail::get_portion_size(to)), from_mat(from), to_mat(to)
{
}

void EM_mat_copy_dispatcher::create_task(off_t global_start, size_t length)
{
	size_t global_start_row, global_start_col;
	size_t num_rows, num_cols;
	if (to_mat.is_wide()) {
		global_start_row = 0;
		global_start_col = global_start;
		num_rows = to_mat.get_num_rows();
		num_cols = length;
	}
	else {
		global_start_row = global_start;
		global_start_col = 0;
		num_rows = length;
		num_cols = to_mat.get_num_cols();
	}
	local_matrix_store::ptr buf;
	if (to_mat.store_layout() == matrix_layout_t::L_COL)
		buf = local_matrix_store::ptr(new local_buf_col_matrix_store(
					global_start_row, global_start_col, num_rows, num_cols,
					to_mat.get_type(), -1));
	else
		buf = local_matrix_store::ptr(new local_buf_row_matrix_store(
					global_start_row, global_start_col, num_rows, num_cols,
					to_mat.get_type(), -1));
	auto from_portion_size = from_mat.get_portion_size();
	auto to_portion_size = to_mat.get_portion_size();
	if (from_portion_size.first == to_portion_size.first
			&& from_portion_size.second == to_portion_size.second) {
		local_matrix_store::const_ptr lstore = from_mat.get_portion(
				global_start_row, global_start_col, num_rows, num_cols);
		assert(lstore);
		buf->copy_from(*lstore);
	}
	else if (from_portion_size.first < to_portion_size.first
			&& from_portion_size.second == to_portion_size.second) {
		size_t to_len = std::min(to_portion_size.first,
				to_mat.get_num_rows() - global_start_row);
		size_t from_len = from_portion_size.first;
		for (size_t lstart = 0; lstart < to_len; lstart += from_len) {
			size_t llen = std::min(from_len, to_len - lstart);
			local_matrix_store::const_ptr lstore = from_mat.get_portion(
					global_start_row + lstart, global_start_col, llen, num_cols);
			assert(lstore);
			buf->resize(lstart, 0, llen, num_cols);
			buf->copy_from(*lstore);
		}
		buf->reset_size();
	}
	else if (from_portion_size.first == to_portion_size.first
			&& from_portion_size.second < to_portion_size.second) {
		size_t to_len = std::min(to_portion_size.second,
				to_mat.get_num_cols() - global_start_col);
		size_t from_len = from_portion_size.second;
		for (size_t lstart = 0; lstart < to_len; lstart += from_len) {
			size_t llen = std::min(from_len, to_len - lstart);
			llen = std::min(llen,
					from_mat.get_num_cols() - global_start_col - lstart);
			local_matrix_store::const_ptr lstore = from_mat.get_portion(
					global_start_row, global_start_col + lstart, num_rows, llen);
			assert(lstore);
			buf->resize(0, lstart, num_rows, llen);
			buf->copy_from(*lstore);
		}
		buf->reset_size();
	}
	else {
		// We shouldn't reach here.
		abort();
	}
	to_mat.write_portion_async(buf, global_start_row, global_start_col);
}

}

bool EM_matrix_store::copy_from(matrix_store::const_ptr mat)
{
	assert(mat->is_in_mem());
	if (mat->get_num_rows() != get_num_rows()
			|| mat->get_num_cols() != get_num_cols()
			|| mat->get_type() != get_type()) {
		BOOST_LOG_TRIVIAL(error) << "copy from an incompatible matrix";
		return false;
	}

	mem_thread_pool::ptr threads = mem_thread_pool::get_global_mem_threads();
	EM_mat_copy_dispatcher::ptr dispatcher(
			new EM_mat_copy_dispatcher(*mat, *this));
	for (size_t i = 0; i < threads->get_num_threads(); i++) {
		io_worker_task *task = new io_worker_task(dispatcher, 1);
		task->register_EM_obj(this);
		threads->process_task(i % threads->get_num_nodes(), task);
	}
	threads->wait4complete();
	return true;
}

matrix_store::const_ptr EM_matrix_store::transpose() const
{
	matrix_layout_t new_layout;
	if (layout == matrix_layout_t::L_ROW)
		new_layout = matrix_layout_t::L_COL;
	else
		new_layout = matrix_layout_t::L_ROW;
	return matrix_store::const_ptr(new EM_matrix_store(holder, ios,
				get_num_cols(), get_num_rows(), new_layout, get_type(),
				data_id));
}

}

}
