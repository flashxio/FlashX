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

#include "log.h"
#include "thread.h"

#include "mem_matrix_store.h"
#include "local_matrix_store.h"
#include "NUMA_dense_matrix.h"
#include "mem_worker_thread.h"

namespace fm
{

namespace detail
{

const size_t mem_matrix_store::CHUNK_SIZE = 64 * 1024;

mem_matrix_store::mem_matrix_store(size_t nrow, size_t ncol,
		const scalar_type &type): matrix_store(nrow, ncol, true,
			type), mat_id(mat_counter++)
{
}

vec_store::const_ptr mem_col_matrix_store::get_row_vec(off_t row) const
{
	assert(data.get_num_bytes()
			== get_num_rows() * get_num_cols() * get_entry_size());
	if (get_num_rows() == 1)
		return detail::smp_vec_store::create(data, get_type());
	else {
		BOOST_LOG_TRIVIAL(error)
			<< "Can't get a row from a column matrix with multiple rows";
		return vec_store::const_ptr();
	}
}

vec_store::const_ptr mem_col_matrix_store::get_col_vec(off_t col) const
{
	assert(data.get_num_bytes()
			== get_num_rows() * get_num_cols() * get_entry_size());
	detail::smp_vec_store::ptr ret = detail::smp_vec_store::create(data,
			get_type());
	ret->expose_sub_vec(col * get_num_rows(), get_num_rows());
	return ret;
}

vec_store::const_ptr mem_sub_col_matrix_store::get_row_vec(off_t row) const
{
	BOOST_LOG_TRIVIAL(error)
		<< "Can't get a row from a column sub matrix";
	return vec_store::const_ptr();
}

vec_store::const_ptr mem_sub_col_matrix_store::get_col_vec(off_t col) const
{
	// The original column matrix has at least this many columns.
	size_t orig_num_cols = orig_col_idxs[col] + 1;
	assert(get_data().get_num_bytes()
			>= get_num_rows() * orig_num_cols * get_entry_size());
	detail::smp_vec_store::ptr ret = detail::smp_vec_store::create(get_data(),
			get_type());
	ret->expose_sub_vec(orig_col_idxs[col] * get_num_rows(), get_num_rows());
	return ret;
}

vec_store::const_ptr mem_row_matrix_store::get_row_vec(off_t row) const
{
	assert(data.get_num_bytes()
			== get_num_rows() * get_num_cols() * get_entry_size());
	detail::smp_vec_store::ptr ret = detail::smp_vec_store::create(data,
			get_type());
	ret->expose_sub_vec(row * get_num_cols(), get_num_cols());
	return ret;
}

vec_store::const_ptr mem_row_matrix_store::get_col_vec(off_t col) const
{
	assert(data.get_num_bytes()
			== get_num_rows() * get_num_cols() * get_entry_size());
	if (get_num_cols() == 1)
		return detail::smp_vec_store::create(data, get_type());
	else {
		BOOST_LOG_TRIVIAL(error)
			<< "Can't get a column from a row matrix with multiple columns";
		return vec_store::const_ptr();
	}
}

vec_store::const_ptr mem_sub_row_matrix_store::get_row_vec(off_t row) const
{
	// The original row matrix has at least this many rows.
	size_t orig_num_rows = orig_row_idxs[row] + 1;
	assert(get_data().get_num_bytes()
			>= get_num_cols() * orig_num_rows * get_entry_size());
	detail::smp_vec_store::ptr ret = detail::smp_vec_store::create(get_data(),
			get_type());
	ret->expose_sub_vec(orig_row_idxs[row] * get_num_cols(), get_num_cols());
	return ret;
}

vec_store::const_ptr mem_sub_row_matrix_store::get_col_vec(off_t row) const
{
	BOOST_LOG_TRIVIAL(error)
		<< "Can't get a column from a row sub matrix";
	return vec_store::const_ptr();
}

mem_matrix_store::ptr mem_matrix_store::create(size_t nrow, size_t ncol,
		matrix_layout_t layout, const scalar_type &type, int num_nodes)
{
	if (num_nodes < 0) {
		if (layout == matrix_layout_t::L_ROW)
			return detail::mem_row_matrix_store::create(nrow, ncol, type);
		else
			return detail::mem_col_matrix_store::create(nrow, ncol, type);
	}
	else
		return detail::NUMA_matrix_store::create(nrow, ncol, num_nodes,
				layout, type);
}

local_matrix_store::const_ptr mem_col_matrix_store::get_portion(
		size_t start_row, size_t start_col, size_t num_rows,
		size_t num_cols) const
{
	if (start_row + num_rows > get_num_rows()
			|| start_col + num_cols > get_num_cols()) {
		BOOST_LOG_TRIVIAL(error) << "it's out of bounds";
		return local_matrix_store::const_ptr();
	}
	int node_id = -1;
	// For a wide matrix
	if (start_row == 0 && num_rows == get_num_rows())
		return local_matrix_store::const_ptr(new local_cref_contig_col_matrix_store(
					get_col(start_col), start_row, start_col,
					num_rows, num_cols, get_type(), node_id));
	// For a tall matrix
	else {
		std::vector<const char *> cols(num_cols);
		for (size_t i = 0; i < num_cols; i++)
			cols[i] = get_col(i + start_col) + start_row * get_entry_size();
		return local_matrix_store::const_ptr(new local_cref_col_matrix_store(
					cols, start_row, start_col, num_rows, num_cols, get_type(),
					node_id));
	}
}

local_matrix_store::ptr mem_col_matrix_store::get_portion(
		size_t start_row, size_t start_col, size_t num_rows,
		size_t num_cols)
{
	if (start_row + num_rows > get_num_rows()
			|| start_col + num_cols > get_num_cols()) {
		BOOST_LOG_TRIVIAL(error) << "it's out of bounds";
		return local_matrix_store::ptr();
	}
	int node_id = -1;
	// For a wide matrix
	if (start_row == 0 && num_rows == get_num_rows())
		return local_matrix_store::ptr(new local_ref_contig_col_matrix_store(
					get_col(start_col), start_row, start_col,
					num_rows, num_cols, get_type(), node_id));
	// For a tall matrix
	else {
		std::vector<char *> cols(num_cols);
		for (size_t i = 0; i < num_cols; i++)
			cols[i] = get_col(i + start_col) + start_row * get_entry_size();
		return local_matrix_store::ptr(new local_ref_col_matrix_store(
					cols, start_row, start_col, num_rows, num_cols, get_type(),
					node_id));
	}
}

local_matrix_store::const_ptr mem_row_matrix_store::get_portion(
		size_t start_row, size_t start_col, size_t num_rows,
		size_t num_cols) const
{
	if (start_row + num_rows > get_num_rows()
			|| start_col + num_cols > get_num_cols()) {
		BOOST_LOG_TRIVIAL(error) << "it's out of bounds";
		return local_matrix_store::const_ptr();
	}
	int node_id = -1;
	// For a tall matrix
	if (start_col == 0 && num_cols == get_num_cols())
		return local_matrix_store::const_ptr(new local_cref_contig_row_matrix_store(
					get_row(start_row), start_row, start_col,
					num_rows, num_cols, get_type(), node_id));
	// For a wide matrix
	else {
		std::vector<const char *> rows(num_rows);
		for (size_t i = 0; i < num_rows; i++)
			rows[i] = get_row(i + start_row) + start_col * get_entry_size();
		return local_matrix_store::const_ptr(new local_cref_row_matrix_store(
					rows, start_row, start_col, num_rows, num_cols, get_type(),
					node_id));
	}
}

local_matrix_store::ptr mem_row_matrix_store::get_portion(
		size_t start_row, size_t start_col, size_t num_rows,
		size_t num_cols)
{
	if (start_row + num_rows > get_num_rows()
			|| start_col + num_cols > get_num_cols()) {
		BOOST_LOG_TRIVIAL(error) << "it's out of bounds";
		return local_matrix_store::ptr();
	}
	int node_id = -1;
	// For a tall matrix
	if (start_col == 0 && num_cols == get_num_cols())
		return local_matrix_store::ptr(new local_ref_contig_row_matrix_store(
					get_row(start_row), start_row, start_col,
					num_rows, num_cols, get_type(), node_id));
	// For a wide matrix
	else {
		std::vector<char *> rows(num_rows);
		for (size_t i = 0; i < num_rows; i++)
			rows[i] = get_row(i + start_row) + start_col * get_entry_size();
		return local_matrix_store::ptr(new local_ref_row_matrix_store(
					rows, start_row, start_col, num_rows, num_cols, get_type(),
					node_id));
	}
}

local_matrix_store::ptr mem_sub_col_matrix_store::get_portion(
		size_t start_row, size_t start_col, size_t num_rows,
		size_t num_cols)
{
	if (start_row + num_rows > get_num_rows()
			|| start_col + num_cols > get_num_cols()) {
		BOOST_LOG_TRIVIAL(error) << "it's out of bounds";
		return local_matrix_store::ptr();
	}

	std::vector<char *> cols(num_cols);
	for (size_t i = 0; i < num_cols; i++)
		cols[i] = get_col(i + start_col) + start_row * get_entry_size();
	int node_id = -1;
	return local_matrix_store::ptr(new local_ref_col_matrix_store(
				cols, start_row, start_col, num_rows, num_cols, get_type(), node_id));
}

local_matrix_store::const_ptr mem_sub_col_matrix_store::get_portion(
		size_t start_row, size_t start_col, size_t num_rows,
		size_t num_cols) const
{
	if (start_row + num_rows > get_num_rows()
			|| start_col + num_cols > get_num_cols()) {
		BOOST_LOG_TRIVIAL(error) << "it's out of bounds";
		return local_matrix_store::ptr();
	}

	std::vector<const char *> cols(num_cols);
	for (size_t i = 0; i < num_cols; i++)
		cols[i] = get_col(i + start_col) + start_row * get_entry_size();
	int node_id = -1;
	return local_matrix_store::const_ptr(new local_cref_col_matrix_store(
				cols, start_row, start_col, num_rows, num_cols, get_type(),
				node_id));
}

local_matrix_store::const_ptr mem_sub_row_matrix_store::get_portion(
			size_t start_row, size_t start_col, size_t num_rows,
			size_t num_cols) const
{
	if (start_row + num_rows > get_num_rows()
			|| start_col + num_cols > get_num_cols()) {
		BOOST_LOG_TRIVIAL(error) << "it's out of bounds";
		return local_matrix_store::ptr();
	}

	std::vector<const char *> rows(num_rows);
	for (size_t i = 0; i < num_rows; i++)
		rows[i] = get_row(i + start_row) + start_col * get_entry_size();
	int node_id = -1;
	return local_matrix_store::const_ptr(new local_cref_row_matrix_store(
				rows, start_row, start_col, num_rows, num_cols, get_type(),
				node_id));
}

local_matrix_store::ptr mem_sub_row_matrix_store::get_portion(
			size_t start_row, size_t start_col, size_t num_rows,
			size_t num_cols)
{
	if (start_row + num_rows > get_num_rows()
			|| start_col + num_cols > get_num_cols()) {
		BOOST_LOG_TRIVIAL(error) << "it's out of bounds";
		return local_matrix_store::ptr();
	}

	std::vector<char *> rows(num_rows);
	for (size_t i = 0; i < num_rows; i++)
		rows[i] = get_row(i + start_row) + start_col * get_entry_size();
	int node_id = -1;
	return local_matrix_store::ptr(new local_ref_row_matrix_store(
				rows, start_row, start_col, num_rows, num_cols, get_type(),
				node_id));
}

matrix_store::const_ptr mem_col_matrix_store::transpose() const
{
	return mem_row_matrix_store::create(data, get_num_cols(), get_num_rows(),
			get_type());
}

matrix_store::const_ptr mem_row_matrix_store::transpose() const
{
	return mem_col_matrix_store::create(data, get_num_cols(), get_num_rows(),
			get_type());
}

matrix_store::const_ptr mem_sub_col_matrix_store::transpose() const
{
	return mem_sub_row_matrix_store::create(get_data(), orig_col_idxs,
			get_num_rows(), get_type());
}

matrix_store::const_ptr mem_sub_row_matrix_store::transpose() const
{
	return mem_sub_col_matrix_store::create(get_data(), orig_row_idxs,
			get_num_cols(), get_type());
}

matrix_store::const_ptr mem_col_matrix_store::get_cols(
		const std::vector<off_t> &idxs) const
{
	return mem_sub_col_matrix_store::create(*this, idxs);
}

matrix_store::const_ptr mem_row_matrix_store::get_rows(
		const std::vector<off_t> &idxs) const
{
	return mem_sub_row_matrix_store::create(*this, idxs);
}

matrix_store::const_ptr mem_sub_col_matrix_store::get_cols(
		const std::vector<off_t> &idxs) const
{
	std::vector<off_t> direct_idxs(idxs.size());
	for (size_t i = 0; i < idxs.size(); i++) {
		if ((size_t) idxs[i] >= get_num_cols()) {
			BOOST_LOG_TRIVIAL(error)
				<< "a column index is out of bounds";
			return matrix_store::ptr();
		}
		direct_idxs[i] = orig_col_idxs[idxs[i]];
	}

	return mem_sub_col_matrix_store::create(*this, direct_idxs);
}

matrix_store::const_ptr mem_sub_row_matrix_store::get_rows(
		const std::vector<off_t> &idxs) const
{
	std::vector<off_t> direct_idxs(idxs.size());
	for (size_t i = 0; i < idxs.size(); i++) {
		if ((size_t) idxs[i] >= get_num_rows()) {
			BOOST_LOG_TRIVIAL(error)
				<< "a row index is out of bounds";
			return matrix_store::ptr();
		}
		direct_idxs[i] = orig_row_idxs[idxs[i]];
	}

	return mem_sub_row_matrix_store::create(*this, direct_idxs);
}

bool mem_col_matrix_store::write2file(const std::string &file_name) const
{
	FILE *f = fopen(file_name.c_str(), "w");
	if (f == NULL) {
		BOOST_LOG_TRIVIAL(error)
			<< boost::format("can't open %1%: %2%") % file_name % strerror(errno);
		return false;
	}
	if (!write_header(f))
		return false;

	size_t ncol = get_num_cols();
	size_t col_size = get_num_rows() * get_entry_size();
	for (size_t i = 0; i < ncol; i++) {
		const char *col = get_col(i);
		size_t ret = fwrite(col, col_size, 1, f);
		if (ret == 0) {
			BOOST_LOG_TRIVIAL(error)
				<< boost::format("can't write to %1%: %2%")
				% file_name % strerror(errno);
			return false;
		}
	}
	fclose(f);
	return true;
}

bool mem_row_matrix_store::write2file(const std::string &file_name) const
{
	FILE *f = fopen(file_name.c_str(), "w");
	if (f == NULL) {
		BOOST_LOG_TRIVIAL(error)
			<< boost::format("can't open %1%: %2%") % file_name % strerror(errno);
		return false;
	}
	if (!write_header(f))
		return false;

	size_t nrow = get_num_rows();
	size_t row_size = get_num_cols() * get_entry_size();
	for (size_t i = 0; i < nrow; i++) {
		const char *row = get_row(i);
		size_t ret = fwrite(row, row_size, 1, f);
		if (ret == 0) {
			BOOST_LOG_TRIVIAL(error)
				<< boost::format("can't write to %1%: %2%")
				% file_name % strerror(errno);
			fclose(f);
			return false;
		}
	}
	fclose(f);
	return true;
}

bool mem_matrix_store::write_header(FILE *f) const
{
	matrix_header header(DENSE, get_entry_size(), get_num_rows(),
			get_num_cols(), store_layout(), get_type().get_type());

	size_t ret = fwrite(&header, sizeof(header), 1, f);
	if (ret == 0) {
		BOOST_LOG_TRIVIAL(error)
			<< boost::format("can't write header: %2%") % strerror(errno);
		return false;
	}
	return true;
}

mem_matrix_store::ptr mem_matrix_store::load(const std::string &file_name)
{
	matrix_header header;

	FILE *f = fopen(file_name.c_str(), "r");
	if (f == NULL) {
		fclose(f);
		BOOST_LOG_TRIVIAL(error) << boost::format("can't open %1%: %2%")
			% file_name % strerror(errno);
		return mem_matrix_store::ptr();
	}

	size_t ret = fread(&header, sizeof(header), 1, f);
	if (ret == 0) {
		fclose(f);
		BOOST_LOG_TRIVIAL(error) << boost::format("can't read header from %1%: %2%")
			% file_name % strerror(errno);
		return mem_matrix_store::ptr();
	}

	header.verify();
	if (header.is_sparse()) {
		fclose(f);
		BOOST_LOG_TRIVIAL(error) << "The matrix to be loaded is sparse";
		return mem_matrix_store::ptr();
	}

	size_t nrow = header.get_num_rows();
	size_t ncol = header.get_num_cols();
	const scalar_type &type = get_scalar_type(header.get_data_type());
	size_t mat_size = nrow * ncol * type.get_size();
	detail::raw_data_array data(mat_size);
	ret = fread(data.get_raw(), mat_size, 1, f);
	if (ret == 0) {
		BOOST_LOG_TRIVIAL(error)
			<< boost::format("can't read %1% bytes from the file") % mat_size;
		return mem_matrix_store::ptr();
	}

	mem_matrix_store::ptr m;
	if (header.get_layout() == matrix_layout_t::L_ROW)
		m = mem_row_matrix_store::create(data, nrow, ncol, type);
	else if (header.get_layout() == matrix_layout_t::L_COL)
		m = mem_col_matrix_store::create(data, nrow, ncol, type);
	else
		BOOST_LOG_TRIVIAL(error) << "wrong matrix data layout";

	fclose(f);
	return m;
}

mem_matrix_store::ptr mem_matrix_store::cast(matrix_store::ptr store)
{
	if (!store->is_in_mem()) {
		BOOST_LOG_TRIVIAL(error)
			<< "cast to col matrix: the matrix isn't in memory";
		return mem_matrix_store::ptr();
	}
	return std::static_pointer_cast<mem_matrix_store>(store);
}

mem_matrix_store::const_ptr mem_matrix_store::cast(matrix_store::const_ptr store)
{
	if (!store->is_in_mem()) {
		BOOST_LOG_TRIVIAL(error)
			<< "cast to col matrix: the matrix isn't in memory";
		return mem_matrix_store::const_ptr();
	}
	return std::static_pointer_cast<const mem_matrix_store>(store);
}

mem_col_matrix_store::const_ptr mem_col_matrix_store::cast(matrix_store::const_ptr store)
{
	if (!store->is_in_mem()) {
		BOOST_LOG_TRIVIAL(error)
			<< "cast to col matrix: the matrix isn't in memory";
		return mem_col_matrix_store::const_ptr();
	}
	mem_matrix_store::const_ptr mem_store
		= std::static_pointer_cast<const mem_matrix_store>(store);
	if (mem_store->get_num_nodes() >= 0) {
		BOOST_LOG_TRIVIAL(error)
			<< "can't cast from a NUMA matrix store";
		return mem_col_matrix_store::const_ptr();
	}
	if (store->store_layout() != matrix_layout_t::L_COL) {
		BOOST_LOG_TRIVIAL(error)
			<< "cast to col matrix: the matrix has to be column major";
		return mem_col_matrix_store::const_ptr();
	}

	return std::static_pointer_cast<const mem_col_matrix_store>(store);
}

mem_col_matrix_store::ptr mem_col_matrix_store::cast(matrix_store::ptr store)
{
	if (!store->is_in_mem()) {
		BOOST_LOG_TRIVIAL(error)
			<< "cast to col matrix: the matrix isn't in memory";
		return mem_col_matrix_store::ptr();
	}
	mem_matrix_store::ptr mem_store
		= std::static_pointer_cast<mem_matrix_store>(store);
	if (mem_store->get_num_nodes() >= 0) {
		BOOST_LOG_TRIVIAL(error)
			<< "can't cast from a NUMA matrix store";
		return mem_col_matrix_store::ptr();
	}
	if (store->store_layout() != matrix_layout_t::L_COL) {
		BOOST_LOG_TRIVIAL(error)
			<< "cast to col matrix: the matrix has to be column major";
		return mem_col_matrix_store::ptr();
	}

	return std::static_pointer_cast<mem_col_matrix_store>(store);
}

mem_row_matrix_store::const_ptr mem_row_matrix_store::cast(matrix_store::const_ptr store)
{
	if (!store->is_in_mem()) {
		BOOST_LOG_TRIVIAL(error)
			<< "cast to col matrix: the matrix isn't in memory";
		return mem_row_matrix_store::const_ptr();
	}
	mem_matrix_store::const_ptr mem_store
		= std::static_pointer_cast<const mem_matrix_store>(store);
	if (mem_store->get_num_nodes() >= 0) {
		BOOST_LOG_TRIVIAL(error)
			<< "can't cast from a NUMA matrix store";
		return mem_row_matrix_store::const_ptr();
	}
	if (store->store_layout() != matrix_layout_t::L_ROW) {
		BOOST_LOG_TRIVIAL(error)
			<< "cast to row matrix: the matrix has to be row major";
		return mem_row_matrix_store::const_ptr();
	}

	return std::static_pointer_cast<const mem_row_matrix_store>(store);
}

mem_row_matrix_store::ptr mem_row_matrix_store::cast(matrix_store::ptr store)
{
	if (!store->is_in_mem()) {
		BOOST_LOG_TRIVIAL(error)
			<< "cast to col matrix: the matrix isn't in memory";
		return mem_row_matrix_store::ptr();
	}
	mem_matrix_store::ptr mem_store
		= std::static_pointer_cast<mem_matrix_store>(store);
	if (mem_store->get_num_nodes() >= 0) {
		BOOST_LOG_TRIVIAL(error)
			<< "can't cast from a NUMA matrix store";
		return mem_row_matrix_store::ptr();
	}
	if (store->store_layout() != matrix_layout_t::L_ROW) {
		BOOST_LOG_TRIVIAL(error)
			<< "cast to row matrix: the matrix has to be row major";
		return mem_row_matrix_store::ptr();
	}

	return std::static_pointer_cast<mem_row_matrix_store>(store);
}

namespace
{

class reset_task: public thread_task
{
	detail::local_matrix_store::ptr local_store;
public:
	reset_task(detail::local_matrix_store::ptr local_store) {
		this->local_store = local_store;
	}

	void run() {
		local_store->reset_data();
	}
};

class set_task: public thread_task
{
	detail::local_matrix_store::ptr local_store;
	const set_operate &op;
public:
	set_task(detail::local_matrix_store::ptr local_store,
			const set_operate &_op): op(_op) {
		this->local_store = local_store;
	}

	void run() {
		local_store->set_data(op);
	}
};

class copy_task: public thread_task
{
	detail::local_matrix_store::const_ptr src_store;
	detail::local_matrix_store::ptr dest_store;
public:
	copy_task(detail::local_matrix_store::const_ptr src_store,
			detail::local_matrix_store::ptr dest_store) {
		this->src_store = src_store;
		this->dest_store = dest_store;
	}

	void run() {
		dest_store->copy_from(*src_store);
	}
};

}

void mem_matrix_store::reset_data()
{
	size_t num_chunks = get_num_portions();
	detail::mem_thread_pool::ptr mem_threads
		= detail::mem_thread_pool::get_global_mem_threads();
	for (size_t i = 0; i < num_chunks; i++) {
		detail::local_matrix_store::ptr local_store = get_portion(i);

		int node_id = local_store->get_node_id();
		// If the local matrix portion is not assigned to any node, 
		// assign the tasks in round robin fashion.
		if (node_id < 0)
			node_id = i % mem_threads->get_num_nodes();
		mem_threads->process_task(node_id, new reset_task(local_store));
	}
	mem_threads->wait4complete();
}

void mem_matrix_store::set_data(const set_operate &op)
{
	size_t num_chunks = get_num_portions();
	detail::mem_thread_pool::ptr mem_threads
		= detail::mem_thread_pool::get_global_mem_threads();
	for (size_t i = 0; i < num_chunks; i++) {
		detail::local_matrix_store::ptr local_store = get_portion(i);

		int node_id = local_store->get_node_id();
		// If the local matrix portion is not assigned to any node, 
		// assign the tasks in round robin fashion.
		if (node_id < 0)
			node_id = i % mem_threads->get_num_nodes();
		mem_threads->process_task(node_id, new set_task(local_store, op));
	}
	mem_threads->wait4complete();
}

matrix_store::const_ptr mem_col_matrix_store::append_cols(
		const std::vector<matrix_store::const_ptr> &mats) const
{
	for (size_t i = 0; i < mats.size(); i++) {
		if (!mats[i]->is_in_mem()) {
			BOOST_LOG_TRIVIAL(error)
				<< "The columns aren't in memory";
			return matrix_store::const_ptr();
		}
		if (mats[i]->get_num_rows() != get_num_rows()) {
			BOOST_LOG_TRIVIAL(error)
				<< "can't append columns with different length";
			return matrix_store::const_ptr();
		}
		if (mats[i]->get_type() != get_type()) {
			BOOST_LOG_TRIVIAL(error)
				<< "can't append columns with different type";
			return matrix_store::const_ptr();
		}
	}
	assert(!is_wide());
	std::vector<const char *> src_cols;
	for (size_t i = 0; i < get_num_cols(); i++)
		src_cols.push_back(get_col(i));
	for (size_t i = 0; i < mats.size(); i++) {
		for (size_t j = 0; j < mats[i]->get_num_cols(); j++) {
			const mem_col_matrix_store &mem_mat
				= dynamic_cast<const mem_col_matrix_store &>(*mats[i]);
			src_cols.push_back(mem_mat.get_col(j));
		}
	}
	mem_col_matrix_store::ptr ret = mem_col_matrix_store::create(
			get_num_rows(), src_cols.size(), get_type());
#pragma omp parallel for
	for (size_t i = 0; i < src_cols.size(); i++)
		memcpy(ret->get_col(i), src_cols[i], get_num_rows() * get_entry_size());
	return ret;
}

matrix_store::const_ptr mem_row_matrix_store::append_cols(
		const std::vector<matrix_store::const_ptr> &mats) const
{
	throw unsupported_exception(
			"can't add columns to a row-major matrix");
}

}

}
