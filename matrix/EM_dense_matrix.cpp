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

const size_t EM_matrix_store::CHUNK_SIZE = 16 * 1024;

static std::unordered_map<std::string, EM_object::file_holder::ptr> file_holders;

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

EM_matrix_store::ptr EM_matrix_store::load(const std::string &ext_file,
		size_t nrow, size_t ncol, matrix_layout_t layout,
		const scalar_type &type, safs::safs_file_group::ptr group)
{
	EM_matrix_store::ptr store = EM_matrix_store::create(nrow, ncol, layout,
			type, group);
	if (store == NULL)
		return store;

	safs::safs_file f(safs::get_sys_RAID_conf(), store->holder->get_name());
	bool ret = f.load_data(ext_file);
	if (!ret)
		return EM_matrix_store::ptr();

	return store;
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

		// We need its transpose.
		if (ret1->store_layout() != store_layout())
			ret1 = std::static_pointer_cast<const local_matrix_store>(
					ret1->transpose());

		local_matrix_store::const_ptr ret;
		if (local_start_row > 0 || local_start_col > 0
				|| num_rows < fetch_num_rows || num_cols < fetch_num_cols)
			ret = ret1->get_portion(local_start_row, local_start_col,
					num_rows, num_cols);
		else
			ret = ret1;
		assert((size_t) ret->get_global_start_row() == start_row);
		assert((size_t) ret->get_global_start_col() == start_col);
		assert(ret->get_num_rows() == num_rows);
		assert(ret->get_num_cols() == num_cols);
		return async_cres_t(valid_data, ret);
	}

	local_raw_array data_arr(num_bytes);
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

void EM_matrix_store::_write_portion_async(
		local_matrix_store::const_ptr portion, off_t start_row,
		off_t start_col)
{
	assert(store_layout() == portion->store_layout());
	assert(start_row + portion->get_num_rows() <= get_num_rows()
			&& start_col + portion->get_num_cols() <= get_num_cols());
	// If this is a tall column-major matrix with more than one column,
	// we need to make sure that the portion to be written is either aligned
	// with the portion size or is at the end of the matrix.
	if (!is_wide() && store_layout() == matrix_layout_t::L_COL
			&& get_num_cols() > 1) {
		assert(start_row % CHUNK_SIZE == 0 && start_col % CHUNK_SIZE == 0);
		if (portion->get_num_rows() % CHUNK_SIZE > 0)
			assert(portion->get_num_rows() == get_num_rows() - start_row);
	}
	// If this is a wide row-major matrix with more than one row,
	// we need to make sure that the portion to be written is either aligned
	// with the portion size or is at the end of the matrix.
	if (is_wide() && store_layout() == matrix_layout_t::L_ROW
			&& get_num_rows() > 1) {
		assert(start_row % CHUNK_SIZE == 0 && start_col % CHUNK_SIZE == 0);
		if (portion->get_num_cols() % CHUNK_SIZE > 0)
			assert(portion->get_num_cols() == get_num_cols() - start_col);
	}

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
	assert(off % PAGE_SIZE == 0);

	size_t num_bytes
		= portion->get_num_rows() * portion->get_num_cols() * entry_size;
	// If this is the very last portion (the bottom right portion), the data
	// size may not be aligned with the page size.
	if (num_bytes % PAGE_SIZE != 0) {
		local_raw_array data_arr(ROUNDUP(num_bytes, PAGE_SIZE));
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
	detail::matrix_stats.inc_write_bytes(num_bytes, false);
	portion_compute::ptr compute(new portion_write_complete(portion));
	static_cast<portion_callback &>(io.get_callback()).add(req, compute);
	io.access(&req, 1);
	io.flush_requests();
}

void EM_matrix_store::write_portion_async(
		local_matrix_store::const_ptr portion, off_t start_row,
		off_t start_col)
{
	if (stream)
		stream->write_async(portion, start_row, start_col);
	else
		_write_portion_async(portion, start_row, start_col);
}

namespace
{

struct empty_deleter {
	void operator()(EM_matrix_store *addr) {
	}
};

}

void EM_matrix_store::start_stream()
{
	stream = EM_matrix_stream::create(std::shared_ptr<EM_matrix_store>(this,
				empty_deleter()));
}

void EM_matrix_store::end_stream()
{
	wait4complete();
	assert(stream->is_complete());
	stream = NULL;
}

void EM_matrix_store::wait4complete()
{
	safs::io_interface &io = ios->get_curr_io();
	io.wait4complete(io.num_pending_ios());
}

std::vector<safs::io_interface::ptr> EM_matrix_store::create_ios() const
{
	std::vector<safs::io_interface::ptr> ret(1);
	ret[0] = ios->create_io();
	return ret;
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
class sub_EM_matrix_store: public matrix_store, public EM_object
{
	const size_t mat_id;
	const size_t data_id;
	EM_matrix_store::const_ptr orig;
	std::vector<off_t> rc_idxs;

	/*
	 * We can only get sub matrix from a col-major tall matrix and a row-major
	 * wide matrix.
	 */
	static size_t cal_num_rows(const std::vector<off_t> idxs,
			const EM_matrix_store &store) {
		if (store.is_wide())
			return idxs.size();
		else
			return store.get_num_rows();
	}
	static size_t cal_num_cols(const std::vector<off_t> idxs,
			const EM_matrix_store &store) {
		if (store.is_wide())
			return store.get_num_cols();
		else
			return idxs.size();
	}

	sub_EM_matrix_store(const std::vector<off_t> idxs,
			EM_matrix_store::const_ptr store, size_t _data_id): matrix_store(
				cal_num_rows(idxs, *store), cal_num_cols(idxs, *store), false,
				store->get_type()), mat_id(mat_counter++), data_id(_data_id) {
		this->rc_idxs = idxs;
		this->orig = store;
		assert(idxs.size() > 0);
	}
public:
	sub_EM_matrix_store(const std::vector<off_t> idxs,
			EM_matrix_store::const_ptr store): matrix_store(cal_num_rows(idxs, *store),
				cal_num_cols(idxs, *store), false, store->get_type()), mat_id(
				// The sub matrix has a different matrix ID and data ID from
				// its parent matrix.
				mat_counter++), data_id(mat_id) {
		this->rc_idxs = idxs;
		this->orig = store;
		assert(idxs.size() > 0);
	}

	virtual std::unordered_map<size_t, size_t> get_underlying_mats() const {
		return orig->get_underlying_mats();
	}
	virtual matrix_layout_t store_layout() const {
		return orig->store_layout();
	}
	virtual std::pair<size_t, size_t> get_portion_size() const {
		return orig->get_portion_size();
	}
	virtual int get_portion_node_id(size_t id) const {
		return -1;
	}

	virtual std::string get_name() const {
		return (boost::format("sub_EM_mat-%1%(%2%,%3%)") % mat_id
				% get_num_rows() % get_num_cols()).str();
	}

	virtual matrix_store::const_ptr transpose() const {
		EM_matrix_store::const_ptr t_mat
			= std::static_pointer_cast<const EM_matrix_store>(orig->transpose());
		return matrix_store::const_ptr(new sub_EM_matrix_store(rc_idxs,
					t_mat, data_id));
	}

	virtual std::shared_ptr<const local_matrix_store> get_portion(
			size_t start_row, size_t start_col, size_t num_rows,
			size_t num_cols) const {
		// TODO let's not support this first.
		assert(0);
		return local_matrix_store::ptr();
	}
	virtual std::shared_ptr<local_matrix_store> get_portion(
			size_t start_row, size_t start_col, size_t num_rows,
			size_t num_cols) {
		// This doesn't need to be used. Changing the data in the local portion
		// doesn't affect the data in the disks.
		assert(0);
		return local_matrix_store::ptr();
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
		if (is_wide()) {
			matrix_store::const_ptr cols = orig->get_cols(idxs);
			return cols->get_rows(rc_idxs);
		}
		else {
			std::vector<off_t> orig_idxs(idxs.size());
			for (size_t i = 0; i < idxs.size(); i++)
				orig_idxs[i] = rc_idxs[idxs[i]];
			return orig->get_cols(orig_idxs);
		}
	}
	virtual matrix_store::const_ptr get_rows(
			const std::vector<off_t> &idxs) const {
		if (!is_wide()) {
			matrix_store::const_ptr rows = orig->get_rows(idxs);
			return rows->get_cols(rc_idxs);
		}
		else {
			std::vector<off_t> orig_idxs(idxs.size());
			for (size_t i = 0; i < idxs.size(); i++)
				orig_idxs[i] = rc_idxs[idxs[i]];
			return orig->get_rows(orig_idxs);
		}
	}

	virtual void write_portion_async(local_matrix_store::const_ptr portion,
			off_t start_row, off_t start_col) {
		BOOST_LOG_TRIVIAL(error)
			<< "Don't support write_portion_async in a sub EM matrix";
	}

	virtual std::vector<safs::io_interface::ptr> create_ios() const {
		return orig->create_ios();
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
		start_fetched_row = start_row - (start_row % EM_matrix_store::CHUNK_SIZE);
		num_fetched_rows = std::min(EM_matrix_store::CHUNK_SIZE,
				get_num_rows() - start_fetched_row);
	}

	// Fetch the portion from the cache.
	local_matrix_store::const_ptr ret1 = local_mem_buffer::get_mat_portion(
			data_id);

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
		async_cres_t res = orig->get_portion_async(
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
	if (orig->is_cache_portion())
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
		start_fetched_col = start_col - (start_col % EM_matrix_store::CHUNK_SIZE);
		num_fetched_cols = std::min(EM_matrix_store::CHUNK_SIZE,
				get_num_cols() - start_fetched_col);
	}
	else {
		start_fetched_col = start_col;
		num_fetched_cols = num_cols;
	}

	// Fetch the portion from the cache.
	local_matrix_store::const_ptr ret1 = local_mem_buffer::get_mat_portion(
			data_id);
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
	for (size_t i = 0; i < ranges.size(); i++) {
		size_t local_num_rows = ranges[i].second - ranges[i].first;
		collected_rows += local_num_rows;
		async_cres_t res = orig->get_portion_async(
				ranges[i].first, start_fetched_col, local_num_rows,
				num_fetched_cols, collect_compute);
		// We assume the requested portion doesn't have data, so the callback
		// function is invoked when the data is ready.
		assert(!res.first);
		bufs[i] = res.second;
	}
	assert(collected_rows == num_rows);
	collect_compute->set_bufs(bufs);
	if (orig->is_cache_portion())
		local_mem_buffer::cache_portion(data_id, ret);
	if (num_fetched_cols == num_cols)
		return async_cres_t(false, ret);
	else
		return async_cres_t(false, ret->get_portion(
					0, start_col - start_fetched_col, num_rows, num_cols));
}

struct comp_row_idx
{
	bool operator()(const std::pair<off_t, off_t> &a,
			const std::pair<off_t, off_t> &b) const {
		return a.first < b.first;
	}
};

matrix_store::const_ptr EM_matrix_store::get_rows(
			const std::vector<off_t> &idxs) const
{
	// We can't get rows from a wide col-major matrix.
	// We'll rely on dense_matrix to get the rows.
	if (store_layout() == matrix_layout_t::L_COL && is_wide())
		return matrix_store::const_ptr();

	// If this is a tall matrix, we read portions that contains the specified
	// rows. This operation is slow, should be used to read a very small number
	// of rows.
	if (!is_wide()) {
		std::pair<size_t, size_t> chunk_size = get_portion_size();
		mem_matrix_store::ptr ret = mem_matrix_store::create(idxs.size(),
				get_num_cols(), matrix_layout_t::L_ROW, get_type(), -1);
		local_row_matrix_store::const_ptr prev;
		// We sort the row idxs, so we only need to keep one portion in memory.
		// But we need to know the offset in the original row idx array before
		// sorted.
		std::vector<std::pair<off_t, off_t> > sorted_idxs(idxs.size());
		for (size_t i = 0; i < idxs.size(); i++) {
			sorted_idxs[i].first = idxs[i];
			sorted_idxs[i].second = i;
		}
		std::sort(sorted_idxs.begin(), sorted_idxs.end(), comp_row_idx());

		for (size_t i = 0; i < sorted_idxs.size(); i++) {
			off_t portion_id = sorted_idxs[i].first / chunk_size.first;
			local_row_matrix_store::const_ptr row_portion;
			// If we don't have a portion previously or the previous portion
			// doesn't match with the required portion, we read the portion
			// from the matrix.
			if (prev == NULL
					|| (size_t) prev->get_global_start_row()
					!= portion_id * chunk_size.first) {
				size_t lnum_rows = std::min(chunk_size.first,
						get_num_rows() - portion_id * chunk_size.first);
				local_matrix_store::const_ptr portion = get_portion(
						portion_id * chunk_size.first, 0, lnum_rows,
						get_num_cols());
				if (portion->store_layout() != matrix_layout_t::L_ROW)
					portion = portion->conv2(matrix_layout_t::L_ROW);
				row_portion = std::static_pointer_cast<const local_row_matrix_store>(
							portion);
				prev = row_portion;
			}
			else
				row_portion = prev;
			off_t lrow_idx = sorted_idxs[i].first - row_portion->get_global_start_row();
			assert(lrow_idx >= 0 && (size_t) lrow_idx < row_portion->get_num_rows());
			memcpy(ret->get_row(sorted_idxs[i].second), row_portion->get_row(lrow_idx),
					row_portion->get_num_cols() * row_portion->get_entry_size());
		}
		return ret;
	}

	return matrix_store::const_ptr(new sub_EM_matrix_store(idxs, shallow_copy()));
}

bool EM_matrix_store::set_persistent(const std::string &name) const
{
	bool persist_ret = holder->set_persistent(name);
	if (!persist_ret)
		return persist_ret;

	// We need to keep holder in a global hashtable, so later on
	// when someone else creates a matrix to access the file, he can
	// use the holder directly from the hashtable.
	auto ret = file_holders.insert(
			std::pair<std::string, file_holder::ptr>(name, holder));
	if (!ret.second) {
		BOOST_LOG_TRIVIAL(error) << "The matrix name already exists";
		holder->unset_persistent();
		return false;
	}
	return true;
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

EM_matrix_stream::filled_portion::filled_portion(local_raw_array &arr,
		off_t start_row, off_t start_col, size_t num_rows, size_t num_cols,
		const scalar_type &type, matrix_layout_t layout, size_t portion_size,
		bool is_wide)
{
	size_t num_bytes_portion;
	if (is_wide)
		num_bytes_portion = portion_size * num_rows * type.get_size();
	else
		num_bytes_portion = portion_size * num_cols * type.get_size();
	assert(arr.get_num_bytes() % num_bytes_portion == 0);
	size_t num_portions = arr.get_num_bytes() / num_bytes_portion;

	if (layout == matrix_layout_t::L_ROW)
		this->data = local_matrix_store::ptr(new local_buf_row_matrix_store(
					arr, start_row, start_col, num_rows, num_cols, type, -1));
	else
		this->data = local_matrix_store::ptr(new local_buf_col_matrix_store(
					arr, start_row, start_col, num_rows, num_cols, type, -1));

	num_filled = 0;
	this->tot_num_eles = num_rows * num_cols;
	this->portion_size = portion_size;
	this->is_wide = is_wide;
	portions.resize(num_portions);

	// create the portions where we can fill data.
	char *raw_arr = arr.get_raw();
	if (is_wide) {
		for (size_t i = 0; i < portions.size(); i++) {
			size_t portion_start_row = start_row;
			size_t portion_start_col = start_col + portion_size * i;
			size_t portion_num_rows = num_rows;
			size_t portion_num_cols = std::min(num_cols - portion_size * i,
					portion_size);
			portions[i] = create_portion(raw_arr + num_bytes_portion * i,
					portion_start_row, portion_start_col, portion_num_rows,
					portion_num_cols, layout, type);
		}
	}
	else {
		for (size_t i = 0; i < portions.size(); i++) {
			size_t portion_start_row = start_row + portion_size * i;
			size_t portion_start_col = start_col;
			size_t portion_num_rows = std::min(num_rows - portion_size * i,
					portion_size);
			size_t portion_num_cols = num_cols;
			portions[i] = create_portion(raw_arr + num_bytes_portion * i,
					portion_start_row, portion_start_col, portion_num_rows,
					portion_num_cols, layout, type);
		}
	}
}

local_matrix_store::ptr EM_matrix_stream::filled_portion::create_portion(
		char *data, off_t start_row, off_t start_col, size_t num_rows,
		size_t num_cols, matrix_layout_t layout, const scalar_type &type)
{
	if (layout == matrix_layout_t::L_COL)
		return local_matrix_store::ptr(new local_ref_contig_col_matrix_store(data,
					start_row, start_col, num_rows, num_cols, type, -1));
	else
		return local_matrix_store::ptr(new local_ref_contig_row_matrix_store(data,
					start_row, start_col, num_rows, num_cols, type, -1));
}

local_matrix_store::ptr EM_matrix_stream::filled_portion::get_portion(
		off_t start_row, off_t start_col)
{
	size_t id;
	if (is_wide)
		id = (start_col - data->get_global_start_col()) / portion_size;
	else
		id = (start_row - data->get_global_start_row()) / portion_size;
	assert(id < portions.size());
	return portions[id];
}

bool EM_matrix_stream::filled_portion::write(
		local_matrix_store::const_ptr in,
		off_t global_start_row, off_t global_start_col)
{
	local_matrix_store::ptr portion = get_portion(global_start_row,
			global_start_col);
	assert(portion);
	local_matrix_store::ptr part = portion->get_portion(
			global_start_row - portion->get_global_start_row(),
			global_start_col - portion->get_global_start_col(),
			in->get_num_rows(), in->get_num_cols());
	part->copy_from(*in);
	size_t ret = num_filled.fetch_add(
			in->get_num_rows() * in->get_num_cols());
	// I assume that no region is filled multiple times.
	return ret + in->get_num_rows() * in->get_num_cols() == tot_num_eles;
}

EM_matrix_stream::EM_matrix_stream(EM_matrix_store::ptr mat)
{
	pthread_spin_init(&lock, PTHREAD_PROCESS_PRIVATE);
	this->mat = mat;
	auto psize = mat->get_portion_size();

	min_io_portions = div_ceil(matrix_conf.get_write_io_buf_size(),
			psize.first * psize.second * mat->get_entry_size());
	// We want the min number of I/O portions to be at least 1 and the smallest
	// number 2^i that is larger than min_io_portions.
	if (min_io_portions == 0)
		min_io_portions = 1;
	else
		min_io_portions = 1 << (size_t) ceil(log2(min_io_portions));
}

void EM_matrix_stream::write_async(local_matrix_store::const_ptr portion,
		off_t start_row, off_t start_col)
{
	if (mat->is_wide())
		assert(start_row == 0 && portion->get_num_rows() == mat->get_num_rows());
	else
		assert(start_col == 0 && portion->get_num_cols() == mat->get_num_cols());

	const size_t CHUNK_SIZE = min_io_portions * EM_matrix_store::CHUNK_SIZE;
	if (mat->is_wide()) {
		// If the portion is aligned with the default EM matrix portion size.
		if (start_col % CHUNK_SIZE == 0
				&& (portion->get_num_cols() % CHUNK_SIZE == 0
				// If this is the last portion.
				|| start_col + portion->get_num_cols() == mat->get_num_cols())) {
			mat->_write_portion_async(portion, start_row, start_col);
			return;
		}
	}
	else {
		// If the portion is aligned with the default EM matrix portion size.
		if (start_row % CHUNK_SIZE == 0
				&& (portion->get_num_rows() % CHUNK_SIZE == 0
				// If this is the last portion.
				|| start_row + portion->get_num_rows() == mat->get_num_rows())) {
			mat->_write_portion_async(portion, start_row, start_col);
			return;
		}
	}

	off_t EM_portion_start;
	if (mat->is_wide())
		EM_portion_start = (start_col / CHUNK_SIZE) * CHUNK_SIZE;
	else
		EM_portion_start = (start_row / CHUNK_SIZE) * CHUNK_SIZE;
	// TODO we might want to a thread-safe hashtable.
	pthread_spin_lock(&lock);
	auto it = portion_bufs.find(EM_portion_start);
	filled_portion::ptr buf;
	if (it == portion_bufs.end()) {
		size_t portion_num_rows, portion_num_cols;
		off_t portion_start_row, portion_start_col;
		size_t num_bytes;
		if (mat->is_wide()) {
			portion_start_row = start_row;
			portion_start_col = EM_portion_start;
			portion_num_rows = mat->get_num_rows();
			portion_num_cols = std::min(CHUNK_SIZE,
					mat->get_num_cols() - EM_portion_start);
			num_bytes = portion_num_rows * CHUNK_SIZE * mat->get_entry_size();
		}
		else  {
			portion_start_row = EM_portion_start;
			portion_start_col = start_col;
			portion_num_rows = std::min(CHUNK_SIZE,
					mat->get_num_rows() - EM_portion_start);
			portion_num_cols = mat->get_num_cols();
			num_bytes = CHUNK_SIZE * portion_num_cols * mat->get_entry_size();
		}
		// We don't want to allocate memory from the local memory buffers
		// because it's not clear which thread will destroy the raw array.
		local_raw_array arr(num_bytes, false);
		buf = filled_portion::ptr(new filled_portion(arr, portion_start_row,
					portion_start_col, portion_num_rows, portion_num_cols,
					mat->get_type(), mat->store_layout(),
					EM_matrix_store::CHUNK_SIZE, mat->is_wide()));
		auto ret = portion_bufs.insert(std::pair<off_t, filled_portion::ptr>(
					EM_portion_start, buf));
		assert(ret.second);
	}
	else
		buf = it->second;
	pthread_spin_unlock(&lock);

	bool ret = buf->write(portion, start_row, start_col);
	// If we fill the buffer, we should flush it to disks.
	if (ret) {
		local_matrix_store::const_ptr data = buf->get_whole_portion();
		if (mat->is_wide())
			assert(data->get_global_start_col() == EM_portion_start);
		else
			assert(data->get_global_start_row() == EM_portion_start);

		mat->_write_portion_async(data, data->get_global_start_row(),
				data->get_global_start_col());
		pthread_spin_lock(&lock);
		portion_bufs.erase(EM_portion_start);
		pthread_spin_unlock(&lock);
	}
}

EM_matrix_stream::~EM_matrix_stream()
{
	pthread_spin_lock(&lock);
	assert(portion_bufs.empty());
	pthread_spin_unlock(&lock);
}

bool EM_matrix_stream::is_complete() const
{
	EM_matrix_stream *mutable_this = const_cast<EM_matrix_stream *>(this);
	pthread_spin_lock(&mutable_this->lock);
	bool ret = portion_bufs.empty();
	pthread_spin_unlock(&mutable_this->lock);
	return ret;
}

}

}
