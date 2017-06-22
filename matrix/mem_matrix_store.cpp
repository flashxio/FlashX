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
#include "matrix_stats.h"
#include "mem_vec_store.h"
#include "dense_matrix.h"

namespace fm
{

namespace detail
{

const size_t mem_matrix_store::CHUNK_SIZE = 16 * 1024;

bool mem_matrix_store::symmetrize(bool upper2lower)
{
	if (get_num_rows() != get_num_cols())
		return false;

	const scatter_gather &sg = get_type().get_sg();
	if (upper2lower && store_layout() == matrix_layout_t::L_ROW) {
		std::vector<char *> non_contig(get_num_rows());
		for (size_t i = 0; i < get_num_rows(); i++) {
			non_contig[i] = get_row(i);
			if (non_contig[i] == NULL)
				return false;
		}

		for (size_t i = 0; i < get_num_rows(); i++) {
			char *row = get_row(i) + i * get_entry_size();
			sg.scatter(row, non_contig);

			std::vector<char *> tmp(non_contig.size() - 1);
			for (size_t j = 0; j < tmp.size(); j++)
				tmp[j] = non_contig[j + 1] + get_entry_size();
			non_contig = tmp;
		}
	}
	else if (upper2lower) {
		std::vector<const char *> non_contig(get_num_cols());
		for (size_t i = 0; i < get_num_cols(); i++) {
			non_contig[i] = get_col(i);
			if (non_contig[i] == NULL)
				return false;
		}

		for (size_t i = 0; i < get_num_cols(); i++) {
			char *col = get_col(i) + i * get_entry_size();
			sg.gather(non_contig, col);

			std::vector<const char *> tmp(non_contig.size() - 1);
			for (size_t j = 0; j < tmp.size(); j++)
				tmp[j] = non_contig[j + 1] + get_entry_size();
			non_contig = tmp;
		}
	}
	else if (store_layout() == matrix_layout_t::L_ROW) {
		std::vector<char *> non_contig(get_num_rows());
		for (size_t i = 0; i < get_num_rows(); i++) {
			if (get_row(i) == NULL)
				return false;
			non_contig[i] = get_row(i) + (get_num_cols() - 1) * get_entry_size();
		}

		for (size_t i = get_num_rows() - 1; i > 0; i--) {
			char *row = get_row(i);
			sg.scatter(row, non_contig);

			std::vector<char *> tmp(non_contig.size() - 1);
			for (size_t j = 0; j < tmp.size(); j++)
				tmp[j] = non_contig[j] - get_entry_size();
			non_contig = tmp;
		}
	}
	else {
		std::vector<const char *> non_contig(get_num_cols());
		for (size_t i = 0; i < get_num_cols(); i++) {
			if (get_col(i) == NULL)
				return false;
			non_contig[i] = get_col(i) + (get_num_rows() - 1) * get_entry_size();
		}

		for (size_t i = get_num_cols() - 1; i > 0; i--) {
			char *col = get_col(i);
			sg.gather(non_contig, col);

			std::vector<const char *> tmp(non_contig.size() - 1);
			for (size_t j = 0; j < tmp.size(); j++)
				tmp[j] = non_contig[j] - get_entry_size();
			non_contig = tmp;
		}
	}
	return true;
}

mem_matrix_store::mem_matrix_store(size_t nrow, size_t ncol,
		const scalar_type &type): matrix_store(nrow, ncol, true,
			type), mat_id(mat_counter++)
{
}

void mem_matrix_store::write_portion_async(local_matrix_store::const_ptr portion,
			off_t start_row, off_t start_col)
{
	if (is_wide()) {
		assert(start_row == 0);
		assert(portion->get_num_rows() == get_num_rows());
		local_matrix_store::ptr lstore = get_portion(start_row,
				start_col, portion->get_num_rows(), portion->get_num_cols());
		assert(lstore);
		lstore->copy_from(*portion);
	}
	else {
		assert(start_col == 0);
		assert(portion->get_num_cols() == get_num_cols());
		local_matrix_store::ptr lstore = get_portion(start_row, start_col,
				portion->get_num_rows(), portion->get_num_cols());
		assert(lstore);
		lstore->copy_from(*portion);
	}
}

mem_matrix_store::ptr mem_matrix_store::create(size_t nrow, size_t ncol,
		matrix_layout_t layout, const scalar_type &type, int num_nodes)
{
	// If the number of nodes aren't specified, or this isn't a very tall or
	// wide matrix, we use a simple way of storing the matrix.
	if (num_nodes < 0 || (nrow <= CHUNK_SIZE && ncol <= CHUNK_SIZE)) {
		if (layout == matrix_layout_t::L_ROW)
			return detail::mem_row_matrix_store::create(nrow, ncol, type);
		else
			return detail::mem_col_matrix_store::create(nrow, ncol, type);
	}
	else
		return detail::NUMA_matrix_store::create(nrow, ncol, num_nodes,
				layout, type);
}

std::string mem_matrix_store::get_name() const
{
	size_t id = get_data_id();
	if (id == INVALID_MAT_ID)
		id = mat_id;
	return (boost::format("mem_mat-%1%(%2%,%3%,%4%)") % id % get_num_rows()
			% get_num_cols()
			% (store_layout() == matrix_layout_t::L_ROW ? "row" : "col")).str();
}

bool mem_matrix_store::share_data(const matrix_store &store) const
{
	const mem_matrix_store *mem_store
		= dynamic_cast<const mem_matrix_store *>(&store);
	size_t store_len = store.get_num_rows() * store.get_num_cols();
	size_t this_len = get_num_rows() * get_num_cols();
	// If the two matrices store data in the same memory and have
	// the same number of elements.
	if (get_raw_arr() && mem_store)
		return mem_store->get_raw_arr() == get_raw_arr()
			&& store_len == this_len;

	matrix_store::const_ptr tstore;
	// If the other matrix might be a transpose of this matrix.
	if (store.get_num_rows() == get_num_cols()
			&& store.get_num_cols() == get_num_rows()
			&& store.store_layout() != store_layout()) {
		tstore = store.transpose();
		mem_store = dynamic_cast<const mem_matrix_store *>(tstore.get());
	}
	if (!mem_store)
		return false;

	if (mem_store->get_num_rows() != get_num_rows()
			|| mem_store->get_num_cols() != get_num_cols()
			|| mem_store->store_layout() != store_layout())
		return false;

	if (mem_store->store_layout() == matrix_layout_t::L_ROW) {
		for (size_t i = 0; i < get_num_rows(); i++)
			if (get_row(i) != mem_store->get_row(i))
				return false;
		return true;
	}
	else {
		for (size_t i = 0; i < get_num_cols(); i++)
			if (get_col(i) != mem_store->get_col(i))
				return false;
		return true;
	}
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

	// Let's only count read bytes from the const version of get_portion.
	detail::matrix_stats.inc_read_bytes(
			num_rows * num_cols * get_entry_size(), true);
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

vec_store::const_ptr mem_col_matrix_store::conv2vec() const
{
	return smp_vec_store::create(get_data(), get_type());
}

vec_store::const_ptr mem_row_matrix_store::conv2vec() const
{
	return smp_vec_store::create(get_data(), get_type());
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

	// Let's only count read bytes from the const version of get_portion.
	detail::matrix_stats.inc_read_bytes(
			num_rows * num_cols * get_entry_size(), true);
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

	// Let's only count read bytes from the const version of get_portion.
	detail::matrix_stats.inc_read_bytes(
			num_rows * num_cols * get_entry_size(), true);
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

	// Let's only count read bytes from the const version of get_portion.
	detail::matrix_stats.inc_read_bytes(
			num_rows * num_cols * get_entry_size(), true);
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

matrix_store::const_ptr mem_col_matrix_store::get_rows(
		const std::vector<off_t> &idxs) const
{
	mem_col_matrix_store::ptr res = mem_col_matrix_store::create(
			idxs.size(), get_num_cols(), get_type());
	std::vector<const char *> src_data(idxs.size());
	for (size_t i = 0; i < get_num_cols(); i++) {
		const char *src_col = get_col(i);
		for (size_t j = 0; j < idxs.size(); j++)
			src_data[j] = src_col + idxs[j] * get_entry_size();
		char *dst_col = res->get_col(i);
		get_type().get_sg().gather(src_data, dst_col);
	}
	return res;
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
			BOOST_LOG_TRIVIAL(error) << "a column index is out of bounds";
			return matrix_store::ptr();
		}
		direct_idxs[i] = orig_col_idxs->at(idxs[i]);
	}

	return mem_sub_col_matrix_store::create(*this, direct_idxs);
}

matrix_store::const_ptr mem_sub_row_matrix_store::get_rows(
		const std::vector<off_t> &idxs) const
{
	std::vector<off_t> direct_idxs(idxs.size());
	for (size_t i = 0; i < idxs.size(); i++) {
		if ((size_t) idxs[i] >= get_num_rows()) {
			BOOST_LOG_TRIVIAL(error) << "a row index is out of bounds";
			return matrix_store::ptr();
		}
		direct_idxs[i] = orig_row_idxs->at(idxs[i]);
	}

	return mem_sub_row_matrix_store::create(*this, direct_idxs);
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

namespace
{

struct empty_deleter {
	void operator()(const matrix_store *addr) {
	}
};

}

static bool write_tall_portions(FILE *f,
		const std::vector<local_matrix_store::const_ptr> &stores)
{
	size_t entry_size = stores[0]->get_type().get_size();
	if (stores[0]->store_layout() == matrix_layout_t::L_ROW) {
		for (size_t i = 0; i < stores.size(); i++) {
			local_row_matrix_store::const_ptr row_store
				= std::dynamic_pointer_cast<const local_row_matrix_store>(
						stores[i]);
			for (size_t j = 0; j < row_store->get_num_rows(); j++) {
				size_t ret = fwrite(row_store->get_row(j),
						row_store->get_num_cols() * entry_size, 1, f);
				if (ret == 0) {
					BOOST_LOG_TRIVIAL(error) << boost::format("can't write: %2%")
						% strerror(errno);
					return false;
				}
			}
		}
	}
	else {
		std::vector<local_col_matrix_store::const_ptr> col_stores(stores.size());
		for (size_t i = 0; i < stores.size(); i++)
			col_stores[i] = std::dynamic_pointer_cast<const local_col_matrix_store>(
						stores[i]);
		size_t num_cols = stores[0]->get_num_cols();
		for (size_t j = 0; j < num_cols; j++) {
			for (size_t i = 0; i < col_stores.size(); i++) {
				size_t ret = fwrite(col_stores[i]->get_col(j),
						col_stores[i]->get_num_rows() * entry_size, 1, f);
				if (ret == 0) {
					BOOST_LOG_TRIVIAL(error) << boost::format("can't write: %2%")
						% strerror(errno);
					return false;
				}
			}
		}
	}
	return true;
}

static bool write_wide_portions(FILE *f,
		const std::vector<local_matrix_store::const_ptr> &stores)
{
	std::vector<local_matrix_store::const_ptr> tstores(stores.size());
	for (size_t i = 0; i < tstores.size(); i++) {
		auto tstore = stores[i]->transpose();
		assert(tstore);
		tstores[i] = std::dynamic_pointer_cast<const local_matrix_store>(tstore);
	}
	return write_tall_portions(f, tstores);
}

static bool write_portions_text(FILE *f,
		const std::vector<local_matrix_store::const_ptr> &stores,
		const std::string &sep)
{
	const scalar_type &type = stores[0]->get_type();
	for (size_t i = 0; i < stores.size(); i++) {
		auto row_store
			= std::dynamic_pointer_cast<const local_row_matrix_store>(stores[i]);
		assert(row_store);
		for (size_t j = 0; j < row_store->get_num_rows(); j++)
			fprintf(f, "%s\n", type.conv2str(row_store->get_row(j),
						row_store->get_num_cols(), sep).c_str());
	}
	return true;
}

bool mem_matrix_store::write2file(const std::string &file_name, bool text,
		std::string sep) const
{
	FILE *f = fopen(file_name.c_str(), "w");
	if (f == NULL) {
		BOOST_LOG_TRIVIAL(error)
			<< boost::format("can't open %1%: %2%") % file_name % strerror(errno);
		return false;
	}
	bool ret;
	if (!text) {
		if (!write_header(f)) {
			fclose(f);
			return false;
		}

		std::vector<local_matrix_store::const_ptr> lstores(get_num_portions());
		for (size_t i = 0; i < lstores.size(); i++)
			lstores[i] = get_portion(i);
		if (is_wide())
			ret = write_wide_portions(f, lstores);
		else
			ret = write_tall_portions(f, lstores);
	}
	else {
		if (is_wide()) {
			fclose(f);
			BOOST_LOG_TRIVIAL(error) << "can't write a wide matrix";
			return false;
		}

		dense_matrix::ptr mat = dense_matrix::create(
				detail::matrix_store::const_ptr(this, empty_deleter()));
		mat = mat->conv2(matrix_layout_t::L_ROW);
		matrix_store::const_ptr row_store = mat->get_raw_store();
		std::vector<local_matrix_store::const_ptr> lstores(get_num_portions());
		for (size_t i = 0; i < lstores.size(); i++) {
			lstores[i] = row_store->get_portion(i);
			assert(lstores[i]->store_layout() == matrix_layout_t::L_ROW);
		}
		ret = write_portions_text(f, lstores, sep);
	}
	fclose(f);
	return ret;
}

mem_matrix_store::const_ptr mem_matrix_store::load(const std::string &file_name,
		int num_nodes)
{
	file_io::ptr io = file_io::create_local(file_name);
	if (io == NULL)
		return mem_matrix_store::ptr();
	return load(io, num_nodes);
}

mem_matrix_store::const_ptr mem_matrix_store::load(file_io::ptr io, int num_nodes)
{
	size_t read_bytes = 0;
	auto data = io->read_bytes(sizeof(matrix_header), read_bytes);
	if (read_bytes != sizeof(matrix_header)) {
		BOOST_LOG_TRIVIAL(error)
			<< boost::format("%1% doesn't contain a header") % io->get_name();
		return mem_matrix_store::ptr();
	}
	const matrix_header *header
		= reinterpret_cast<const matrix_header *>(data.get());
	if (!header->is_matrix_file()) {
		BOOST_LOG_TRIVIAL(error)
			<< boost::format("%1% doesn't contain a valid header") % io->get_name();
		return mem_matrix_store::ptr();
	}
	if (header->is_sparse()) {
		BOOST_LOG_TRIVIAL(error) << "The matrix to be loaded is sparse";
		return mem_matrix_store::ptr();
	}

	size_t nrow = header->get_num_rows();
	size_t ncol = header->get_num_cols();
	const scalar_type &type = header->get_data_type();
	return load_raw(io, nrow, ncol, header->get_layout(), type, num_nodes);
}

mem_matrix_store::const_ptr mem_matrix_store::load_raw(
		const std::string &file_name, size_t nrow, size_t ncol,
		matrix_layout_t layout, const scalar_type &type, int num_nodes)
{
	file_io::ptr io = file_io::create_local(file_name);
	if (io == NULL)
		return mem_matrix_store::ptr();
	return load_raw(io, nrow, ncol, layout, type, num_nodes);
}

mem_matrix_store::const_ptr mem_matrix_store::load_raw(file_io::ptr io,
		size_t nrow, size_t ncol, matrix_layout_t layout,
		const scalar_type &type, int num_nodes)
{
	mem_matrix_store::ptr m;
	// In these two cases, we can read portion by portion.
	if ((ncol > nrow && layout == matrix_layout_t::L_COL)
			|| (nrow > ncol && layout == matrix_layout_t::L_ROW)) {
		m = mem_matrix_store::create(nrow, ncol, layout, type, num_nodes);
		for (size_t i = 0; i < m->get_num_portions(); i++) {
			local_matrix_store::ptr portion = m->get_portion(i);
			size_t psize = portion->get_num_rows() * portion->get_num_cols();
			psize *= portion->get_entry_size();
			size_t read_bytes = 0;
			auto data = io->read_bytes(psize, read_bytes);
			if (read_bytes != psize) {
				BOOST_LOG_TRIVIAL(error)
					<< boost::format("try to read %1% bytes and get %2% bytes")
					% psize % read_bytes;
				return mem_matrix_store::ptr();
			}
			memcpy(portion->get_raw_arr(), data.get(), psize);
		}
	}
	else {
		// Here the matrix can be tall col matrix or wide row matrix.
		// In these two cases, the data in the entire cols or rows are stored
		// contiguously. so we can assume the data in the file is a tall col
		// matrix and transpose the matrix if the data in the file is a wide
		// row matrix.
		m = mem_matrix_store::create(std::max(nrow, ncol), std::min(nrow, ncol),
				matrix_layout_t::L_COL, type, num_nodes);
		std::vector<local_matrix_store::ptr> lstores(m->get_num_portions());
		for (size_t i = 0; i < lstores.size(); i++)
			lstores[i] = m->get_portion(i);
		for (size_t i = 0; i < m->get_num_cols(); i++) {
			for (size_t j = 0; j < lstores.size(); j++) {
				local_col_matrix_store::ptr lstore
					= std::static_pointer_cast<local_col_matrix_store>(lstores[j]);
				size_t size = lstore->get_num_rows() * lstore->get_entry_size();
				size_t read_bytes = 0;
				auto data = io->read_bytes(size, read_bytes);
				if (size != read_bytes) {
					BOOST_LOG_TRIVIAL(error)
						<< boost::format("try to read %1% bytes and get %2% bytes")
						% size % read_bytes;
					return mem_matrix_store::ptr();
				}
				memcpy(lstore->get_col(i), data.get(), size);
			}
		}
	}

	if (m->store_layout() != layout)
		return std::static_pointer_cast<const mem_matrix_store>(m->transpose());
	else
		return m;
}

mem_matrix_store::ptr mem_matrix_store::cast(matrix_store::ptr store)
{
	if (!store->is_in_mem()) {
		BOOST_LOG_TRIVIAL(error)
			<< "cast to col matrix: the matrix isn't in memory";
		return mem_matrix_store::ptr();
	}
	return std::dynamic_pointer_cast<mem_matrix_store>(store);
}

mem_matrix_store::const_ptr mem_matrix_store::cast(matrix_store::const_ptr store)
{
	if (!store->is_in_mem()) {
		BOOST_LOG_TRIVIAL(error)
			<< "cast to col matrix: the matrix isn't in memory";
		return mem_matrix_store::const_ptr();
	}
	return std::dynamic_pointer_cast<const mem_matrix_store>(store);
}

mem_col_matrix_store::const_ptr mem_col_matrix_store::cast(matrix_store::const_ptr store)
{
	if (!store->is_in_mem()) {
		BOOST_LOG_TRIVIAL(error)
			<< "cast to col matrix: the matrix isn't in memory";
		return mem_col_matrix_store::const_ptr();
	}
	mem_matrix_store::const_ptr mem_store
		= std::dynamic_pointer_cast<const mem_matrix_store>(store);
	assert(mem_store);
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

	return std::dynamic_pointer_cast<const mem_col_matrix_store>(store);
}

mem_col_matrix_store::ptr mem_col_matrix_store::cast(matrix_store::ptr store)
{
	if (!store->is_in_mem()) {
		BOOST_LOG_TRIVIAL(error)
			<< "cast to col matrix: the matrix isn't in memory";
		return mem_col_matrix_store::ptr();
	}
	mem_matrix_store::ptr mem_store
		= std::dynamic_pointer_cast<mem_matrix_store>(store);
	assert(mem_store);
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

	return std::dynamic_pointer_cast<mem_col_matrix_store>(store);
}

mem_row_matrix_store::const_ptr mem_row_matrix_store::cast(matrix_store::const_ptr store)
{
	if (!store->is_in_mem()) {
		BOOST_LOG_TRIVIAL(error)
			<< "cast to col matrix: the matrix isn't in memory";
		return mem_row_matrix_store::const_ptr();
	}
	mem_matrix_store::const_ptr mem_store
		= std::dynamic_pointer_cast<const mem_matrix_store>(store);
	assert(mem_store);
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

	return std::dynamic_pointer_cast<const mem_row_matrix_store>(store);
}

mem_row_matrix_store::ptr mem_row_matrix_store::cast(matrix_store::ptr store)
{
	if (!store->is_in_mem()) {
		BOOST_LOG_TRIVIAL(error)
			<< "cast to col matrix: the matrix isn't in memory";
		return mem_row_matrix_store::ptr();
	}
	mem_matrix_store::ptr mem_store
		= std::dynamic_pointer_cast<mem_matrix_store>(store);
	assert(mem_store);
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

	return std::dynamic_pointer_cast<mem_row_matrix_store>(store);
}

size_t mem_col_matrix_store::get_data_id() const
{
	// TODO
	return INVALID_MAT_ID;
}

size_t mem_row_matrix_store::get_data_id() const
{
	// TODO
	return INVALID_MAT_ID;
}

size_t mem_sub_col_matrix_store::get_data_id() const
{
	// TODO
	return INVALID_MAT_ID;
}

size_t mem_sub_row_matrix_store::get_data_id() const
{
	// TODO
	return INVALID_MAT_ID;
}

bool mem_row_matrix_store::resize(size_t num_rows, size_t num_cols)
{
	if (num_rows > get_num_rows() || num_cols > get_num_cols()) {
		BOOST_LOG_TRIVIAL(error)
			<< "can't resize a matrix to a larger one";
		return false;
	}
	if (num_rows == get_num_rows() && num_cols == get_num_cols())
		return true;

	// If we only need to reduce the number of rows, we can simply adjust
	// the number of rows and columns.
	if (num_rows < get_num_rows() && num_cols == get_num_cols())
		return matrix_store::resize(num_rows, num_cols);

	// Otherwise, we need to copy the data before adjusting the number of
	// rows and columns.
	mem_row_matrix_store::ptr tmp = mem_row_matrix_store::create(num_rows,
			num_cols, get_type());
	for (size_t i = 0; i < num_rows; i++)
		memcpy(tmp->get_row(i), get_row(i), num_cols * get_entry_size());
	data = tmp->data;
	return matrix_store::resize(num_rows, num_cols);
}

bool mem_col_matrix_store::resize(size_t num_rows, size_t num_cols)
{
	if (num_rows > get_num_rows() || num_cols > get_num_cols()) {
		BOOST_LOG_TRIVIAL(error)
			<< "can't resize a matrix to a larger one";
		return false;
	}
	if (num_rows == get_num_rows() && num_cols == get_num_cols())
		return true;

	// If we only need to reduce the number of cols, we can simply adjust
	// the number of rows and columns.
	if (num_cols < get_num_cols() && num_rows == get_num_rows())
		return matrix_store::resize(num_rows, num_cols);

	// Otherwise, we need to copy the data before adjusting the number of
	// rows and columns.
	mem_col_matrix_store::ptr tmp = mem_col_matrix_store::create(num_rows,
			num_cols, get_type());
	for (size_t i = 0; i < num_cols; i++)
		memcpy(tmp->get_col(i), get_col(i), num_rows * get_entry_size());
	data = tmp->data;
	return matrix_store::resize(num_rows, num_cols);
}

}

}
