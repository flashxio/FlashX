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

namespace fm
{

namespace detail
{

static const size_t CHUNK_SIZE = 1024 * 1024;

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
		const scalar_type &type): matrix_store(nrow, ncol, false, type)
{
	this->layout = layout;
	holder = file_holder::create_temp("mat", nrow * ncol * type.get_size());
	safs::file_io_factory::shared_ptr factory = safs::create_io_factory(
			holder->get_name(), safs::REMOTE_ACCESS);
	ios = io_set::ptr(new io_set(factory));
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
		const set_operate &_op): EM_portion_dispatcher(
			store.is_wide() ? store.get_num_cols() : store.get_num_rows(),
			store.is_wide() ? store.get_portion_size(
				).second : store.get_portion_size().first), op(
				_op), to_mat(store)
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

local_matrix_store::const_ptr EM_matrix_store::get_portion(
		size_t start_row, size_t start_col, size_t num_rows,
		size_t num_cols) const
{
	return const_cast<EM_matrix_store *>(this)->get_portion(start_row,
			start_col, num_rows, num_cols);
}

local_matrix_store::ptr EM_matrix_store::get_portion(
		size_t start_row, size_t start_col, size_t num_rows,
		size_t num_cols)
{
	safs::io_interface &io = ios->get_curr_io();
	bool ready = false;
	portion_compute::ptr compute(new sync_read_compute(ready));
	local_matrix_store::ptr ret = get_portion_async(start_row, start_col,
			num_rows, num_cols, compute);
	if (ret == NULL)
		return ret;

	while (!ready)
		io.wait4complete(1);
	return ret;
}

local_matrix_store::const_ptr EM_matrix_store::get_portion_async(
		size_t start_row, size_t start_col, size_t num_rows,
		size_t num_cols, portion_compute::ptr compute) const
{
	return const_cast<EM_matrix_store *>(this)->get_portion_async(
			start_row, start_col, num_rows, num_cols, compute);
}

local_matrix_store::ptr EM_matrix_store::get_portion_async(
		size_t start_row, size_t start_col, size_t num_rows,
		size_t num_cols, portion_compute::ptr compute)
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

	if (local_start_row > 0 || local_start_col > 0
			|| num_rows < portion_num_rows || num_cols < portion_num_cols)
		buf->resize(local_start_row, local_start_col, num_rows, num_cols);
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

safs::io_interface::ptr EM_matrix_store::create_io()
{
	return ios->create_io();
}

}

}
