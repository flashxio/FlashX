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

#include <boost/format.hpp>

#include "NUMA_dense_matrix.h"
#include "mem_worker_thread.h"
#include "local_matrix_store.h"
#include "matrix_stats.h"
#include "sub_matrix_store.h"

namespace fm
{

namespace detail
{

NUMA_matrix_store::ptr NUMA_matrix_store::cast(matrix_store::ptr store)
{
	mem_matrix_store::ptr mem_store = mem_matrix_store::cast(store);
	if (mem_store == NULL)
		return NUMA_matrix_store::ptr();
	if (mem_store->get_num_nodes() < 0) {
		BOOST_LOG_TRIVIAL(error)
			<< "cast to NUMA matrix: the matrix isn't stored in NUMA memory";
		return NUMA_matrix_store::ptr();
	}
	return std::static_pointer_cast<NUMA_matrix_store>(store);
}

NUMA_matrix_store::const_ptr NUMA_matrix_store::cast(matrix_store::const_ptr store)
{
	mem_matrix_store::const_ptr mem_store = mem_matrix_store::cast(store);
	if (mem_store == NULL)
		return NUMA_matrix_store::const_ptr();
	if (mem_store->get_num_nodes() < 0) {
		BOOST_LOG_TRIVIAL(error)
			<< "cast to NUMA matrix: the matrix isn't stored in NUMA memory";
		return NUMA_matrix_store::const_ptr();
	}
	return std::static_pointer_cast<const NUMA_matrix_store>(store);
}

NUMA_matrix_store::ptr NUMA_matrix_store::create(size_t nrow, size_t ncol,
		int num_nodes, matrix_layout_t layout, const scalar_type &type)
{
	if (layout == matrix_layout_t::L_ROW) {
		if (nrow > ncol)
			return NUMA_row_tall_matrix_store::create(nrow, ncol, num_nodes,
					type);
		else
			return NUMA_row_wide_matrix_store::create(nrow, ncol, num_nodes,
					type);
	}
	else {
		if (nrow > ncol)
			return NUMA_col_tall_matrix_store::create(nrow, ncol, num_nodes,
					type);
		else
			return NUMA_col_wide_matrix_store::create(nrow, ncol, num_nodes,
					type);
	}
}

NUMA_row_tall_matrix_store::NUMA_row_tall_matrix_store(
		const NUMA_row_tall_matrix_store &mat): NUMA_matrix_store(
			mat.get_num_rows(), mat.get_num_cols(),
			mat.get_type(), mat.get_data_id()), mapper(mat.get_num_nodes(),
			NUMA_range_size_log)
{
	this->data = mat.data;
}

NUMA_row_tall_matrix_store::NUMA_row_tall_matrix_store(
		const std::vector<detail::chunked_raw_array> &data, size_t nrow, size_t ncol,
		const NUMA_mapper &_mapper, const scalar_type &type): NUMA_matrix_store(
			nrow, ncol, type, mat_counter++), mapper(_mapper)
{
	this->data = data;
	std::vector<size_t> local_lens = mapper.cal_local_lengths(nrow);
	assert(data.size() == mapper.get_num_nodes());
	size_t block_bytes = mem_matrix_store::CHUNK_SIZE * ncol * get_entry_size();
	for (size_t node_id = 0; node_id < mapper.get_num_nodes(); node_id++) {
		assert((size_t) data[node_id].get_node_id() == node_id);
		assert(data[node_id].get_num_bytes()
				>= local_lens[node_id] * ncol * get_entry_size());
		assert(data[node_id].get_contig_block_size() % block_bytes == 0);
	}
}

NUMA_row_tall_matrix_store::NUMA_row_tall_matrix_store(size_t nrow, size_t ncol,
		int num_nodes, const scalar_type &type): NUMA_matrix_store(nrow, ncol,
			type, mat_counter++), mapper(num_nodes, NUMA_range_size_log)
{
	data.resize(num_nodes);
	std::vector<size_t> local_lens = mapper.cal_local_lengths(nrow);
	size_t block_bytes = mem_matrix_store::CHUNK_SIZE * ncol * get_entry_size();
	for (int node_id = 0; node_id < num_nodes; node_id++)
		data[node_id] = detail::chunked_raw_array(
				local_lens[node_id] * ncol * get_entry_size(), block_bytes,
				node_id);
}

char *NUMA_row_tall_matrix_store::get_row(size_t row_idx)
{
	auto phy_loc = mapper.map2physical(row_idx);
	return data[phy_loc.first].get_raw(
			phy_loc.second * get_num_cols() * get_entry_size());
}

const char *NUMA_row_tall_matrix_store::get_row(size_t row_idx) const
{
	auto phy_loc = mapper.map2physical(row_idx);
	return data[phy_loc.first].get_raw(
			phy_loc.second * get_num_cols() * get_entry_size());
}

matrix_store::const_ptr NUMA_row_tall_matrix_store::get_rows(
		const std::vector<off_t> &idxs) const
{
	// I assume we only need to get a small number of rows from the matrix.
	// So it's sufficient to output a SMP matrix to store the result.
	// TODO this might be slow
	mem_row_matrix_store::ptr ret = mem_row_matrix_store::create(
			idxs.size(), get_num_cols(), get_type());
	for (size_t i = 0; i < idxs.size(); i++)
		memcpy(ret->get_row(i), get_row(idxs[i]),
				get_entry_size() * get_num_cols());
	return ret;
}

matrix_store::const_ptr NUMA_row_tall_matrix_store::get_cols(
		const std::vector<off_t> &idxs) const
{
	NUMA_row_tall_matrix_store::ptr ret = NUMA_row_tall_matrix_store::create(
			get_num_rows(), idxs.size(), get_num_nodes(), get_type());
	std::vector<const char *> src_data(idxs.size());
	for (size_t i = 0; i < get_num_rows(); i++) {
		const char *src_row = get_row(i);
		for (size_t j = 0; j < src_data.size(); j++)
			src_data[j] = src_row + idxs[j] * get_entry_size();
		get_type().get_sg().gather(src_data, ret->get_row(i));
	}
	return ret;
}

NUMA_col_tall_matrix_store::NUMA_col_tall_matrix_store(size_t nrow,
		size_t ncol, int num_nodes, const scalar_type &type): NUMA_matrix_store(
			nrow, ncol, type, mat_counter++), mapper(num_nodes,
			NUMA_range_size_log)
{
	data.resize(num_nodes);
	std::vector<size_t> local_lens = mapper.cal_local_lengths(nrow);
	size_t block_bytes = mem_matrix_store::CHUNK_SIZE * ncol * get_entry_size();
	for (int node_id = 0; node_id < num_nodes; node_id++)
		data[node_id] = detail::chunked_raw_array(
				local_lens[node_id] * ncol * get_entry_size(), block_bytes,
				node_id);
}

matrix_store::const_ptr NUMA_row_tall_matrix_store::transpose() const
{
	return NUMA_col_wide_matrix_store::create_transpose(*this);
}

matrix_store::const_ptr NUMA_col_tall_matrix_store::transpose() const
{
	return NUMA_row_wide_matrix_store::create_transpose(*this);
}

matrix_store::const_ptr NUMA_row_wide_matrix_store::transpose() const
{
	return matrix_store::const_ptr(new NUMA_col_tall_matrix_store(store));
}

matrix_store::const_ptr NUMA_col_wide_matrix_store::transpose() const
{
	return matrix_store::const_ptr(new NUMA_row_tall_matrix_store(store));
}

bool NUMA_row_tall_matrix_store::get_portion_check(size_t start_row,
		size_t start_col, size_t num_rows, size_t num_cols) const
{
	if (start_row + num_rows > get_num_rows()) {
		BOOST_LOG_TRIVIAL(error) << "get a portion out of boundary";
		return false;
	}
	// We have to retrieve the entire rows.
	if (num_cols != get_num_cols() || start_col != 0) {
		BOOST_LOG_TRIVIAL(error)
			<< "has to get a portion with all elements in the short dimension";
		return false;
	}

	// The retrieved rows have to be stored contiguously.
	size_t range_size = mapper.get_range_size();
	if (ROUND(start_row, range_size)
			!= ROUND(start_row + num_rows - 1, range_size)) {
		BOOST_LOG_TRIVIAL(error) << "data isn't in the same range";
		return false;
	}
	return true;
}

local_matrix_store::const_ptr NUMA_row_tall_matrix_store::get_portion(
		size_t start_row, size_t start_col, size_t num_rows,
		size_t num_cols) const
{
	if (!get_portion_check(start_row, start_col, num_rows, num_cols))
		return local_matrix_store::const_ptr();

	// Let's only count read bytes from the const version of get_portion.
	detail::matrix_stats.inc_read_bytes(
			num_rows * num_cols * get_entry_size(), true);

	auto phy_loc_start = mapper.map2physical(start_row);
	auto phy_loc_end = mapper.map2physical(start_row + num_rows - 1);
	assert(phy_loc_start.first == phy_loc_end.first);
	const char * addr = data[phy_loc_start.first].get_raw(
			phy_loc_start.second * get_num_cols() * get_entry_size(),
			(phy_loc_end.second + 1) * get_num_cols() * get_entry_size());
	// If all rows are stored in contiguous memory, we can return them immediately.
	if (addr)
		return local_matrix_store::const_ptr(new local_cref_contig_row_matrix_store(
					addr, start_row, start_col, num_rows, num_cols,
					get_type(), phy_loc_start.first));

	// If the required data isn't stored in contiguous memory, let's copy them
	// to a piece of contiguous memory. This doesn't happen very frequently.
	// Normally, it only happens when we perform operations on matrices with
	// different portion sizes. so it should be fine.
	// TODO In the future, we should avoid memory copy here.
	size_t num_sub_portions
		= ceil(((double) num_rows) / mem_matrix_store::CHUNK_SIZE);
	size_t local_start_row = 0;
	local_buf_row_matrix_store::ptr ret(new local_buf_row_matrix_store(
				start_row, start_col, num_rows, num_cols, get_type(), -1));
	for (size_t i = 0; i < num_sub_portions; i++) {
		addr = get_row(start_row + local_start_row);
		size_t local_num_rows = std::min(mem_matrix_store::CHUNK_SIZE,
				num_rows - local_start_row);
		local_cref_contig_row_matrix_store sub_portion(addr,
				start_row + local_start_row, start_col, local_num_rows,
				num_cols, get_type(), -1);
		ret->resize(local_start_row, 0, local_num_rows, ret->get_num_cols());
		ret->copy_from(sub_portion);
		local_start_row += local_num_rows;
	}
	ret->reset_size();
	return ret;
}

local_matrix_store::ptr NUMA_row_tall_matrix_store::get_portion(
		size_t start_row, size_t start_col, size_t num_rows,
		size_t num_cols)
{
	if (!get_portion_check(start_row, start_col, num_rows, num_cols))
		return local_matrix_store::ptr();

	auto phy_loc_start = mapper.map2physical(start_row);
	auto phy_loc_end = mapper.map2physical(start_row + num_rows - 1);
	assert(phy_loc_start.first == phy_loc_end.first);
	char * addr = data[phy_loc_start.first].get_raw(
			phy_loc_start.second * get_num_cols() * get_entry_size(),
			(phy_loc_end.second + 1) * get_num_cols() * get_entry_size());
	// TODO let's hope all data is stored in contiguous memory.
	assert(addr);
	return local_matrix_store::ptr(new local_ref_contig_row_matrix_store(
				addr, start_row, start_col, num_rows, num_cols,
				get_type(), phy_loc_start.first));
}

local_matrix_store::const_ptr NUMA_row_tall_matrix_store::get_portion(
		size_t id) const
{
	size_t chunk_size = get_portion_size().first;
	size_t start_row = id * chunk_size;
	size_t start_col = 0;
	size_t num_rows = std::min(get_num_rows() - start_row, chunk_size);
	size_t num_cols = get_num_cols();
	auto phy_loc = mapper.map2physical(start_row);

	// Let's only count read bytes from the const version of get_portion.
	detail::matrix_stats.inc_read_bytes(
			num_rows * num_cols * get_entry_size(), true);
	return local_matrix_store::const_ptr(new local_cref_contig_row_matrix_store(
				get_row(start_row), start_row, start_col, num_rows, num_cols,
				get_type(), phy_loc.first));
}

local_matrix_store::ptr NUMA_row_tall_matrix_store::get_portion(size_t id)
{
	size_t chunk_size = get_portion_size().first;
	size_t start_row = id * chunk_size;
	size_t start_col = 0;
	size_t num_rows = std::min(get_num_rows() - start_row, chunk_size);
	size_t num_cols = get_num_cols();
	auto phy_loc = mapper.map2physical(start_row);
	return local_matrix_store::ptr(new local_ref_contig_row_matrix_store(
				get_row(start_row), start_row, start_col, num_rows, num_cols,
				get_type(), phy_loc.first));
}

int NUMA_row_tall_matrix_store::get_portion_node_id(size_t id) const
{
	size_t chunk_size = get_portion_size().first;
	size_t start_row = id * chunk_size;
	auto phy_loc = mapper.map2physical(start_row);
	return phy_loc.first;
}

matrix_store::const_ptr NUMA_col_tall_matrix_store::get_rows(
		const std::vector<off_t> &idxs) const
{
	// I assume we only need to get a small number of rows from the matrix.
	// So it's sufficient to output a SMP matrix to store the result.
	mem_col_matrix_store::ptr res = mem_col_matrix_store::create(
			idxs.size(), get_num_cols(), get_type());
	std::vector<const char *> src_data(idxs.size());
	for (size_t i = 0; i < get_num_cols(); i++) {
		// TODO this is too expensive.
		for (size_t j = 0; j < idxs.size(); j++)
			src_data[j] = get(idxs[j], i);
		char *dst_col = res->get_col(i);
		get_type().get_sg().gather(src_data, dst_col);
	}
	return res;
}

local_matrix_store::const_ptr NUMA_col_tall_matrix_store::get_portion(
		size_t start_row, size_t start_col, size_t num_rows,
		size_t num_cols) const
{
	if (start_row + num_rows > get_num_rows()) {
		BOOST_LOG_TRIVIAL(error) << "get a portion out of boundary";
		return local_matrix_store::const_ptr();
	}

	// Let's only count read bytes from the const version of get_portion.
	detail::matrix_stats.inc_read_bytes(
			num_rows * num_cols * get_entry_size(), true);

	size_t chunk_size = get_portion_size().first;
	size_t portion_id = start_row / chunk_size;
	// We can't get a portion stored across multiple chunks.
	if ((start_row + num_rows - 1) / chunk_size != portion_id) {
		BOOST_LOG_TRIVIAL(error) << "cannot get data across multiple chunks";
		return local_matrix_store::ptr();
	}

	local_col_matrix_store::const_ptr ret
		= std::static_pointer_cast<const local_col_matrix_store>(
				get_portion(portion_id));
	// If we want all rows in the portion.
	if (start_row % chunk_size == 0 && (num_rows == chunk_size
				|| start_row + num_rows == get_num_rows())) {
		// If we want the entire portion. It should be a common case.
		if (start_col == 0 && num_cols == get_num_cols())
			return ret;
		else
			return local_matrix_store::const_ptr(
					new local_cref_contig_col_matrix_store(
						ret->get_col(start_col), start_row, start_col,
						num_rows, num_cols, get_type(), ret->get_node_id()));
	}
	else {
		std::vector<const char *> cols(num_cols);
		off_t off
			= (start_row - ret->get_global_start_row()) * ret->get_entry_size();
		for (size_t i = 0; i < num_cols; i++)
			cols[i] = ret->get_col(i + start_col) + off;
		return local_matrix_store::const_ptr(new local_cref_col_matrix_store(
					cols, start_row, start_col, num_rows, num_cols,
					get_type(), ret->get_node_id()));
	}
}

local_matrix_store::ptr NUMA_col_tall_matrix_store::get_portion(
		size_t start_row, size_t start_col, size_t num_rows,
		size_t num_cols)
{
	if (start_row + num_rows > get_num_rows()) {
		BOOST_LOG_TRIVIAL(error) << "get a portion out of boundary";
		return local_matrix_store::ptr();
	}

	size_t chunk_size = get_portion_size().first;
	size_t portion_id = start_row / chunk_size;
	// We can't get a portion stored across multiple chunks.
	if ((start_row + num_rows - 1) / chunk_size != portion_id) {
		BOOST_LOG_TRIVIAL(error) << "cannot get data across multiple chunks";
		return local_matrix_store::ptr();
	}

	local_col_matrix_store::ptr ret
		= std::static_pointer_cast<local_col_matrix_store>(
				get_portion(portion_id));
	// If we want all rows in the portion.
	if (start_row % chunk_size == 0 && (num_rows == chunk_size
				|| start_row + num_rows == get_num_rows())) {
		// If we want the entire portion. It should be a common case.
		if (start_col == 0 && num_cols == get_num_cols())
			return ret;
		else
			return local_matrix_store::ptr(
					new local_ref_contig_col_matrix_store(
						ret->get_col(start_col), start_row, start_col,
						num_rows, num_cols, get_type(), ret->get_node_id()));
	}
	else {
		std::vector<char *> cols(num_cols);
		off_t off
			= (start_row - ret->get_global_start_row()) * ret->get_entry_size();
		for (size_t i = 0; i < num_cols; i++)
			cols[i] = ret->get_col(i + start_col) + off;
		return local_matrix_store::ptr(new local_ref_col_matrix_store(
					cols, start_row, start_col, num_rows, num_cols,
					get_type(), ret->get_node_id()));
	}
}

local_matrix_store::const_ptr NUMA_col_tall_matrix_store::get_portion(
		size_t id) const
{
	size_t chunk_size = get_portion_size().first;
	size_t start_row = chunk_size * id;
	size_t start_col = 0;
	size_t num_rows = std::min(get_num_rows() - start_row, chunk_size);
	size_t num_cols = get_num_cols();
	// Let's only count read bytes from the const version of get_portion.
	detail::matrix_stats.inc_read_bytes(
			num_rows * num_cols * get_entry_size(), true);

	auto phy_loc_start = mapper.map2physical(start_row);
	off_t last_row = start_row + num_rows - 1;
	auto phy_loc_end = mapper.map2physical(last_row);
	assert(phy_loc_start.first == phy_loc_end.first);
	const char * addr = data[phy_loc_start.first].get_raw(
			phy_loc_start.second * get_num_cols() * get_entry_size(),
			(phy_loc_end.second + 1) * get_num_cols() * get_entry_size());
	assert(addr);
	return local_matrix_store::const_ptr(new local_cref_contig_col_matrix_store(
				addr, start_row, start_col, num_rows, num_cols,
				get_type(), phy_loc_start.first));
}

local_matrix_store::ptr NUMA_col_tall_matrix_store::get_portion(size_t id)
{
	size_t chunk_size = get_portion_size().first;
	size_t start_row = chunk_size * id;
	size_t start_col = 0;
	size_t num_rows = std::min(get_num_rows() - start_row, chunk_size);
	size_t num_cols = get_num_cols();
	auto phy_loc_start = mapper.map2physical(start_row);
	off_t last_row = start_row + num_rows - 1;
	auto phy_loc_end = mapper.map2physical(last_row);
	assert(phy_loc_start.first == phy_loc_end.first);
	char * addr = data[phy_loc_start.first].get_raw(
			phy_loc_start.second * get_num_cols() * get_entry_size(),
			(phy_loc_end.second + 1) * get_num_cols() * get_entry_size());
	assert(addr);
	return local_matrix_store::ptr(new local_ref_contig_col_matrix_store(
				addr, start_row, start_col, num_rows, num_cols,
				get_type(), phy_loc_start.first));
}

int NUMA_col_tall_matrix_store::get_portion_node_id(size_t id) const
{
	size_t chunk_size = get_portion_size().first;
	size_t start_row = id * chunk_size;
	auto phy_loc_start = mapper.map2physical(start_row);
	return phy_loc_start.first;
}

local_matrix_store::const_ptr NUMA_row_wide_matrix_store::get_portion(
		size_t start_row, size_t start_col, size_t num_rows, size_t num_cols) const
{
	return local_matrix_store::cast(store.get_portion(start_col, start_row,
				num_cols, num_rows)->transpose());
}

local_matrix_store::ptr NUMA_row_wide_matrix_store::get_portion(
		size_t start_row, size_t start_col, size_t num_rows, size_t num_cols)
{
	return local_matrix_store::cast(store.get_portion(start_col, start_row,
				num_cols, num_rows)->transpose());
}

local_matrix_store::const_ptr NUMA_col_wide_matrix_store::get_portion(
		size_t start_row, size_t start_col, size_t num_rows, size_t num_cols) const
{
	return local_matrix_store::cast(store.get_portion(start_col, start_row,
				num_cols, num_rows)->transpose());
}

local_matrix_store::ptr NUMA_col_wide_matrix_store::get_portion(
		size_t start_row, size_t start_col, size_t num_rows, size_t num_cols)
{
	return local_matrix_store::cast(store.get_portion(start_col, start_row,
				num_cols, num_rows)->transpose());
}

local_matrix_store::const_ptr NUMA_row_wide_matrix_store::get_portion(
		size_t id) const
{
	return local_matrix_store::cast(store.get_portion(id)->transpose());
}

local_matrix_store::ptr NUMA_row_wide_matrix_store::get_portion(size_t id)
{
	return local_matrix_store::cast(store.get_portion(id)->transpose());
}

local_matrix_store::const_ptr NUMA_col_wide_matrix_store::get_portion(
		size_t id) const
{
	return local_matrix_store::cast(store.get_portion(id)->transpose());
}

local_matrix_store::ptr NUMA_col_wide_matrix_store::get_portion(size_t id)
{
	return local_matrix_store::cast(store.get_portion(id)->transpose());
}

char *NUMA_col_tall_matrix_store::get(size_t row_idx, size_t col_idx)
{
	size_t chunk_size = get_portion_size().first;
	size_t portion_start_row = ROUND(row_idx, chunk_size);
	size_t portion_num_rows = portion_start_row + chunk_size
		> get_num_rows() ? get_num_rows() - portion_start_row : chunk_size;
	auto phy_loc_start = mapper.map2physical(portion_start_row);
	char *addr = data[phy_loc_start.first].get_raw(
			phy_loc_start.second * get_num_cols() * get_entry_size());
	assert(addr);
	off_t local_row_idx = row_idx - portion_start_row;
	return addr + (portion_num_rows * col_idx + local_row_idx) * get_entry_size();
}

const char *NUMA_col_tall_matrix_store::get(size_t row_idx, size_t col_idx) const
{
	size_t chunk_size = get_portion_size().first;
	size_t portion_start_row = ROUND(row_idx, chunk_size);
	size_t portion_num_rows = portion_start_row + chunk_size
		> get_num_rows() ? get_num_rows() - portion_start_row : chunk_size;
	auto phy_loc_start = mapper.map2physical(portion_start_row);
	const char * addr = data[phy_loc_start.first].get_raw(
			phy_loc_start.second * get_num_cols() * get_entry_size());
	assert(addr);
	off_t local_row_idx = row_idx - portion_start_row;
	return addr + (portion_num_rows * col_idx + local_row_idx) * get_entry_size();
}

matrix_store::const_ptr NUMA_col_tall_matrix_store::get_cols(
		const std::vector<off_t> &idxs) const
{
	NUMA_col_tall_matrix_store::const_ptr copy(new NUMA_col_tall_matrix_store(*this));
	return matrix_store::const_ptr(new sub_col_matrix_store(idxs, copy));
}

}

}
