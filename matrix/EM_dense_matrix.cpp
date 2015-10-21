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

#include <unordered_map>
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

const size_t EM_matrix_store::CHUNK_SIZE = 256 * 1024;

static std::unordered_map<std::string, EM_object::file_holder::ptr> file_holders;

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
		const scalar_type &type, safs::safs_file_group::ptr group): matrix_store(
			nrow, ncol, false, type), mat_id(mat_counter++), data_id(mat_id)
{
	this->cache_portion = true;
	this->orig_num_rows = nrow;
	this->orig_num_cols = ncol;
	this->layout = layout;
	holder = file_holder::create_temp("mat", nrow * ncol * type.get_size(),
			group);
	safs::file_io_factory::shared_ptr factory = safs::create_io_factory(
			holder->get_name(), safs::REMOTE_ACCESS);
	ios = io_set::ptr(new io_set(factory));

	// Store the header as the metadata.
	std::vector<char> header_buf(matrix_header::get_header_size());
	new (header_buf.data()) matrix_header(matrix_type::DENSE, type.get_size(),
			nrow, ncol, layout, type.get_type());
	safs::safs_file f(safs::get_sys_RAID_conf(), holder->get_name());
	bool ret = f.set_user_metadata(header_buf);
	assert(ret);
}

EM_matrix_store::EM_matrix_store(file_holder::ptr holder, io_set::ptr ios,
		size_t nrow, size_t ncol, size_t orig_nrow, size_t orig_ncol,
		matrix_layout_t layout, const scalar_type &type,
		size_t _data_id): matrix_store(nrow, ncol, false, type), mat_id(
			mat_counter++), data_id(_data_id)
{
	this->cache_portion = true;
	this->orig_num_rows = orig_nrow;
	this->orig_num_cols = orig_ncol;
	this->layout = layout;
	this->holder = holder;
	this->ios = ios;
}

EM_matrix_store::ptr EM_matrix_store::create(const std::string &mat_file)
{
	// The file holder might already exist in the hashtable.
	// We should create only one holder for a matrix file in the system.
	file_holder::ptr holder;
	auto it = file_holders.find(mat_file);
	if (it != file_holders.end())
		holder = it->second;
	else {
		holder = file_holder::create(mat_file);
		if (holder) {
			auto ret = file_holders.insert(
					std::pair<std::string, file_holder::ptr>(mat_file, holder));
			assert(ret.second);
		}
	}
	if (holder == NULL) {
		BOOST_LOG_TRIVIAL(error) << mat_file + " doesn't exist";
		return EM_matrix_store::ptr();
	}

	safs::file_io_factory::shared_ptr factory = safs::create_io_factory(
			holder->get_name(), safs::REMOTE_ACCESS);
	io_set::ptr ios(new io_set(factory));

	// Read the matrix header.
	safs::safs_file f(safs::get_sys_RAID_conf(), holder->get_name());
	std::vector<char> header_buf = f.get_user_metadata();
	if (header_buf.size() != matrix_header::get_header_size()) {
		BOOST_LOG_TRIVIAL(error) << "Cannot get the matrix header";
		return EM_matrix_store::ptr();
	}

	matrix_header *header = (matrix_header *) header_buf.data();
	EM_matrix_store::ptr ret_mat(new EM_matrix_store(holder, ios,
				header->get_num_rows(), header->get_num_cols(),
				header->get_num_rows(), header->get_num_cols(),
				// TODO we should save the data Id in the matrix header.
				header->get_layout(), header->get_data_type(), mat_counter++));
	return ret_mat;
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
	if (is_wide())
		return std::pair<size_t, size_t>(get_num_rows(), CHUNK_SIZE);
	else
		return std::pair<size_t, size_t>(CHUNK_SIZE, get_num_cols());
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
	async_cres_t ret = get_portion_async(start_row, start_col,
			num_rows, num_cols, compute);
	// If we can't get the specified portion or the portion already has
	// the valid data.
	if (ret.second == NULL || ret.first)
		return ret.second;

	while (!ready)
		io.wait4complete(1);
	return ret.second;
}

async_res_t EM_matrix_store::get_portion_async(
		size_t start_row, size_t start_col, size_t num_rows,
		size_t num_cols, portion_compute::ptr compute)
{
	// This doesn't need to be used. Changing the data in the local portion
	// doesn't affect the data in the disks.
	assert(0);
	return async_res_t();
}

async_cres_t EM_matrix_store::get_portion_async(
		size_t start_row, size_t start_col, size_t num_rows,
		size_t num_cols, portion_compute::ptr compute) const
{
	size_t local_start_row = start_row % CHUNK_SIZE;
	size_t local_start_col = start_col % CHUNK_SIZE;
	// For now, we only allow this method to fetch data from a portion.
	if (local_start_row + num_rows > CHUNK_SIZE
			|| local_start_col + num_cols > CHUNK_SIZE
			|| start_row + num_rows > get_orig_num_rows()
			|| start_col + num_cols > get_orig_num_cols()) {
		BOOST_LOG_TRIVIAL(error) << "Out of boundary of a portion";
		return async_cres_t();
	}

	safs::io_interface &io = ios->get_curr_io();
	size_t entry_size = get_type().get_size();
	// The information of a portion.
	size_t portion_start_row = start_row - local_start_row;
	size_t portion_start_col = start_col - local_start_col;
	size_t portion_num_rows;
	size_t portion_num_cols;
	if (store_layout() == matrix_layout_t::L_COL && num_rows == CHUNK_SIZE)
		portion_num_rows = num_rows;
	else if (portion_start_row + CHUNK_SIZE > get_orig_num_rows())
		portion_num_rows = get_orig_num_rows() - portion_start_row;
	else
		portion_num_rows = CHUNK_SIZE;

	if (store_layout() == matrix_layout_t::L_ROW && num_cols == CHUNK_SIZE)
		portion_num_cols = num_cols;
	else if (portion_start_col + CHUNK_SIZE > get_orig_num_cols())
		portion_num_cols = get_orig_num_cols() - portion_start_col;
	else
		portion_num_cols = CHUNK_SIZE;

	// The information of the part of data fetched in a portion.
	size_t fetch_start_row = portion_start_row;
	size_t fetch_start_col = portion_start_col;
	size_t fetch_num_rows = portion_num_rows;
	size_t fetch_num_cols = portion_num_cols;
	// Location of the portion on the disks.
	// The number of elements above the portion row
	off_t off = (get_orig_num_cols() * portion_start_row
			// The number of elements in front of the wanted portion
			// in the same portion row.
		+ portion_num_rows * portion_start_col) * entry_size;
	if (portion_num_rows < CHUNK_SIZE && portion_num_cols < CHUNK_SIZE) {
		// This is the very last portion, we have to fetch the entire portion.
	}
	// If we fetch data from a col-major matrix and fetch the entire cols,
	// we only need to fetch the wanted columns
	else if (store_layout() == matrix_layout_t::L_COL) {
		fetch_start_col += local_start_col;
		off += local_start_col * portion_num_rows * entry_size;
		local_start_col = 0;
	}
	// For the same reason, we only need to fetch the wanted rows.
	else {
		fetch_start_row += local_start_row;
		off += local_start_row * portion_num_cols * entry_size;
		local_start_row = 0;
	}

	// If this is the very last portion (the bottom right portion), the data
	// size may not be aligned with the page size.
	size_t num_bytes;
	if (portion_num_rows < CHUNK_SIZE && portion_num_cols < CHUNK_SIZE) {
		// We have to fetch the entire portion for the very last portion.
		num_bytes = ROUNDUP(portion_num_rows * portion_num_cols * entry_size,
				PAGE_SIZE);
	}
	else if (store_layout() == matrix_layout_t::L_COL) {
		num_bytes = ROUNDUP(portion_num_rows * num_cols * entry_size, PAGE_SIZE);
		fetch_num_cols = num_cols;
	}
	else {
		num_bytes = ROUNDUP(num_rows * portion_num_cols * entry_size, PAGE_SIZE);
		fetch_num_rows = num_rows;
	}

	// We should try to get the portion from the local thread memory buffer
	// first.
	local_matrix_store::const_ptr ret1 = local_mem_buffer::get_mat_portion(
			data_id);
	// If it's in the same portion.
	if (ret1 && (((size_t) ret1->get_global_start_row() == fetch_start_row
					&& (size_t) ret1->get_global_start_col() == fetch_start_col
					&& ret1->get_num_rows() == fetch_num_rows
					&& ret1->get_num_cols() == fetch_num_cols)
				// If it's in the corresponding portion in the transposed matrix.
				|| ((size_t) ret1->get_global_start_row() == fetch_start_col
					&& (size_t) ret1->get_global_start_col() == fetch_start_row
					&& ret1->get_num_rows() == fetch_num_cols
					&& ret1->get_num_cols() == fetch_num_rows))) {
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
		// it means the data in the portion is ready.
		bool valid_data = !cb.has_callback(req);
		if (!valid_data)
			cb.add(req, compute);

		local_matrix_store::const_ptr ret;
		if (local_start_row > 0 || local_start_col > 0
				|| num_rows < fetch_num_rows || num_cols < fetch_num_cols)
			ret = ret1->get_portion(local_start_row, local_start_col,
					num_rows, num_cols);
		else
			ret = ret1;
		return async_cres_t(valid_data, ret);
	}

	raw_data_array data_arr(num_bytes, -1);
	// Read the portion in a single I/O request.
	local_matrix_store::ptr buf;
	if (store_layout() == matrix_layout_t::L_ROW)
		buf = local_matrix_store::ptr(new local_buf_row_matrix_store(data_arr,
					fetch_start_row, fetch_start_col, fetch_num_rows,
					fetch_num_cols, get_type(), data_arr.get_node_id()));
	else
		buf = local_matrix_store::ptr(new local_buf_col_matrix_store(data_arr,
					fetch_start_row, fetch_start_col, fetch_num_rows,
					fetch_num_cols, get_type(), data_arr.get_node_id()));

	safs::data_loc_t loc(io.get_file_id(), off);
	safs::io_request req(buf->get_raw_arr(), loc, num_bytes, READ);
	static_cast<portion_callback &>(io.get_callback()).add(req, compute);
	io.access(&req, 1);
	io.flush_requests();

	// Count the number of bytes really read from disks.
	detail::matrix_stats.inc_read_bytes(
			buf->get_num_rows() * buf->get_num_cols() * get_entry_size(), false);

	if (cache_portion)
		local_mem_buffer::cache_portion(data_id, buf);
	local_matrix_store::const_ptr ret;
	if (local_start_row > 0 || local_start_col > 0
			|| num_rows < fetch_num_rows || num_cols < fetch_num_cols)
		ret = buf->get_portion(local_start_row, local_start_col,
				num_rows, num_cols);
	else
		ret = buf;
	return async_cres_t(false, ret);
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

matrix_store::const_ptr EM_matrix_store::transpose() const
{
	matrix_layout_t new_layout;
	if (layout == matrix_layout_t::L_ROW)
		new_layout = matrix_layout_t::L_COL;
	else
		new_layout = matrix_layout_t::L_ROW;
	return matrix_store::const_ptr(new EM_matrix_store(holder, ios,
				get_num_cols(), get_num_rows(), get_orig_num_cols(),
				get_orig_num_rows(), new_layout, get_type(), data_id));
}

/*
 * This matrix store accesses a set of rows of a wide matrix or a set of columns
 * of a tall matrix.
 */
class sub_EM_matrix_store: public EM_matrix_store
{
	const size_t mat_id;
	const size_t data_id;
	std::vector<off_t> rc_idxs;

	sub_EM_matrix_store(const std::vector<off_t> idxs,
			const EM_matrix_store &store, size_t _data_id): EM_matrix_store(
				store), mat_id(mat_counter++), data_id(_data_id) {
		this->rc_idxs = idxs;
		assert(idxs.size() > 0);
		if (store_layout() == matrix_layout_t::L_COL)
			matrix_store::resize(get_num_rows(), idxs.size());
		else
			matrix_store::resize(idxs.size(), get_num_cols());
	}
public:
	sub_EM_matrix_store(const std::vector<off_t> idxs,
			const EM_matrix_store &store): EM_matrix_store(store), mat_id(
				// The sub matrix has a different matrix ID and data ID from
				// its parent matrix.
				mat_counter++), data_id(mat_id) {
		this->rc_idxs = idxs;
		assert(idxs.size() > 0);
		if (store_layout() == matrix_layout_t::L_COL)
			matrix_store::resize(get_num_rows(), idxs.size());
		else
			matrix_store::resize(idxs.size(), get_num_cols());
	}

	virtual std::string get_name() const {
		return (boost::format("sub_EM_mat-%1%(%2%,%3%)") % mat_id
				% get_num_rows() % get_num_cols()).str();
	}

	virtual matrix_store::const_ptr transpose() const {
		matrix_store::const_ptr t_mat = EM_matrix_store::transpose();
		return matrix_store::const_ptr(new sub_EM_matrix_store(rc_idxs,
					dynamic_cast<const EM_matrix_store &>(*t_mat), data_id));
	}

	virtual async_cres_t get_col_portion_async(
			size_t start_row, size_t start_col, size_t num_rows,
			size_t num_cols, portion_compute::ptr compute) const;
	virtual async_cres_t get_row_portion_async(
			size_t start_row, size_t start_col, size_t num_rows,
			size_t num_cols, portion_compute::ptr compute) const;
	virtual async_cres_t get_portion_async(
			size_t start_row, size_t start_col, size_t num_rows,
			size_t num_cols, portion_compute::ptr compute) const {
		if (start_row + num_rows > get_num_rows()
				|| start_col + num_cols > get_num_cols()) {
			BOOST_LOG_TRIVIAL(error) << "get portion async: out of boundary";
			return async_cres_t();
		}

		if (store_layout() == matrix_layout_t::L_COL)
			return get_col_portion_async(start_row, start_col,
					num_rows, num_cols, compute);
		else
			return get_row_portion_async(start_row, start_col,
					num_rows, num_cols, compute);
	}

	virtual matrix_store::const_ptr get_cols(
			const std::vector<off_t> &idxs) const {
		std::vector<off_t> orig_idxs(idxs.size());
		for (size_t i = 0; i < idxs.size(); i++)
			orig_idxs[i] = rc_idxs[idxs[i]];
		return EM_matrix_store::get_cols(orig_idxs);
	}
	virtual matrix_store::const_ptr get_rows(
			const std::vector<off_t> &idxs) const {
		std::vector<off_t> orig_idxs(idxs.size());
		for (size_t i = 0; i < idxs.size(); i++)
			orig_idxs[i] = rc_idxs[idxs[i]];
		return EM_matrix_store::get_rows(orig_idxs);
	}

	virtual vec_store::const_ptr get_col_vec(off_t idx) const {
		if (store_layout() == matrix_layout_t::L_COL)
			return EM_matrix_store::get_col_vec(rc_idxs[idx]);
		else {
			BOOST_LOG_TRIVIAL(error)
				<< "Can't get a col from a sub EM row matrix";
			return vec_store::const_ptr();
		}
	}
	virtual vec_store::const_ptr get_row_vec(off_t idx) const {
		if (store_layout() == matrix_layout_t::L_COL) {
			BOOST_LOG_TRIVIAL(error)
				<< "Can't get a row from a sub EM col matrix";
			return vec_store::const_ptr();
		}
		else
			return EM_matrix_store::get_row_vec(rc_idxs[idx]);
	}

	virtual async_res_t get_portion_async(
			size_t start_row, size_t start_col, size_t num_rows,
			size_t num_cols, portion_compute::ptr compute) {
		BOOST_LOG_TRIVIAL(error)
			<< "Don't support non-const get_portion_async in a sub EM matrix";
		return async_res_t();
	}

	virtual void write_portion_async(local_matrix_store::const_ptr portion,
			off_t start_row, off_t start_col) {
		BOOST_LOG_TRIVIAL(error)
			<< "Don't support write_portion_async in a sub EM matrix";
	}
	virtual void reset_data() {
		BOOST_LOG_TRIVIAL(error)
			<< "Don't support reset_data in a sub EM matrix";
	}
	virtual void set_data(const set_operate &op) {
		BOOST_LOG_TRIVIAL(error) << "Don't support set_data in a sub EM matrix";
	}
};

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

async_cres_t sub_EM_matrix_store::get_col_portion_async(
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
		start_fetched_row = start_row - (start_row % CHUNK_SIZE);
		num_fetched_rows = std::min(CHUNK_SIZE,
				get_num_rows() - start_fetched_row);
	}

	// Fetch the portion from the cache.
	local_matrix_store::const_ptr ret1 = local_mem_buffer::get_mat_portion(
			data_id);
	// If it's in the same portion.
	if (ret1 && (((size_t) ret1->get_global_start_row() == start_fetched_row
					&& (size_t) ret1->get_global_start_col() == start_col
					&& ret1->get_num_rows() == num_fetched_rows
					&& ret1->get_num_cols() == num_cols)
				// If it's in the corresponding portion in the transposed matrix.
				|| ((size_t) ret1->get_global_start_row() == start_col
					&& (size_t) ret1->get_global_start_col() == start_fetched_row
					&& ret1->get_num_rows() == num_cols
					&& ret1->get_num_cols() == num_fetched_rows))) {
		assert(ret1->get_local_start_row() == 0);
		assert(ret1->get_local_start_col() == 0);
		// In the asynchronous version, data in the portion isn't ready when
		// the method is called. We should add the user's portion computation
		// to the queue. When the data is ready, all user's portion computations
		// will be invoked.
		local_matrix_store *tmp = const_cast<local_matrix_store *>(ret1.get());
		local_collected_buf_col_matrix_store *store
			= dynamic_cast<local_collected_buf_col_matrix_store *>(tmp);
		collect_rc_compute::ptr collect_compute = store->get_compute();
		// If collect_rc_compute doesn't exist, it mean the data has been read
		// from disks.
		bool valid_data = collect_compute == NULL;
		if (!valid_data)
			collect_compute->add_orig_compute(compute);

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
		async_cres_t res = EM_matrix_store::get_portion_async(
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
	if (is_cache_portion())
		local_mem_buffer::cache_portion(data_id, ret);

	bool ready = num_ready == bufs.size();
	if (num_fetched_rows == num_rows)
		return async_cres_t(ready, ret);
	else
		return async_cres_t(ready, ret->get_portion(
					start_row - start_fetched_row, 0, num_rows, num_cols));
}

async_cres_t sub_EM_matrix_store::get_row_portion_async(
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
		start_fetched_col = start_col - (start_col % CHUNK_SIZE);
		num_fetched_cols = std::min(CHUNK_SIZE,
				get_num_cols() - start_fetched_col);
	}
	else {
		start_fetched_col = start_col;
		num_fetched_cols = num_cols;
	}

	// Fetch the portion from the cache.
	local_matrix_store::const_ptr ret1 = local_mem_buffer::get_mat_portion(
			data_id);
	// If it's in the same portion.
	if (ret1 && (((size_t) ret1->get_global_start_row() == start_row
					&& (size_t) ret1->get_global_start_col() == start_fetched_col
					&& ret1->get_num_rows() == num_rows
					&& ret1->get_num_cols() == num_fetched_cols)
				// If it's in the corresponding portion in the transposed matrix.
				|| ((size_t) ret1->get_global_start_row() == start_fetched_col
					&& (size_t) ret1->get_global_start_col() == start_row
					&& ret1->get_num_rows() == num_fetched_cols
					&& ret1->get_num_cols() == num_rows))) {
		assert(ret1->get_local_start_row() == 0);
		assert(ret1->get_local_start_col() == 0);
		// In the asynchronous version, data in the portion isn't ready when
		// the method is called. We should add the user's portion computation
		// to the queue. When the data is ready, all user's portion computations
		// will be invoked.
		local_matrix_store *tmp = const_cast<local_matrix_store *>(ret1.get());
		local_collected_buf_row_matrix_store *store
			= dynamic_cast<local_collected_buf_row_matrix_store *>(tmp);
		collect_rc_compute::ptr collect_compute = store->get_compute();
		bool valid_data = collect_compute == NULL;
		if (!valid_data)
			collect_compute->add_orig_compute(compute);

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
	for (size_t i = 0; i < ranges.size(); i++) {
		size_t local_num_rows = ranges[i].second - ranges[i].first;
		collected_rows += local_num_rows;
		async_cres_t res = EM_matrix_store::get_portion_async(
				ranges[i].first, start_fetched_col, local_num_rows,
				num_fetched_cols, collect_compute);
		// We assume the requested portion doesn't have data, so the callback
		// function is invoked when the data is ready.
		assert(!res.first);
		bufs[i] = res.second;
	}
	assert(collected_rows == num_rows);
	collect_compute->set_bufs(bufs);
	if (is_cache_portion())
		local_mem_buffer::cache_portion(data_id, ret);
	if (num_fetched_cols == num_cols)
		return async_cres_t(false, ret);
	else
		return async_cres_t(false, ret->get_portion(
					0, start_col - start_fetched_col, num_rows, num_cols));
}

matrix_store::const_ptr EM_matrix_store::get_cols(
			const std::vector<off_t> &idxs) const
{
	if (store_layout() == matrix_layout_t::L_ROW || is_wide()) {
		BOOST_LOG_TRIVIAL(error)
			<< "can't support get cols from a row-major or wide matrix";
		return matrix_store::const_ptr();
	}

	return matrix_store::const_ptr(new sub_EM_matrix_store(idxs, *this));
}

matrix_store::const_ptr EM_matrix_store::get_rows(
			const std::vector<off_t> &idxs) const
{
	if (store_layout() == matrix_layout_t::L_COL || !is_wide()) {
		BOOST_LOG_TRIVIAL(error)
			<< "can't support get cols from a col-major or tall matrix";
		return matrix_store::const_ptr();
	}

	return matrix_store::const_ptr(new sub_EM_matrix_store(idxs, *this));
}

bool EM_matrix_store::set_persistent(const std::string &name) const
{
	// We need to keep holder in a global hashtable, so later on
	// when someone else creates a matrix to access the file, he can
	// use the holder directly from the hashtable.
	auto ret = file_holders.insert(
			std::pair<std::string, file_holder::ptr>(name, holder));
	if (!ret.second) {
		BOOST_LOG_TRIVIAL(error) << "The matrix name already exists";
		return false;
	}
	return holder->set_persistent(name);
}

void EM_matrix_store::unset_persistent() const
{
	size_t ret = file_holders.erase(holder->get_name());
	if (ret == 0) {
		assert(!holder->is_persistent());
		return;
	}
	holder->unset_persistent();
}

}

}
