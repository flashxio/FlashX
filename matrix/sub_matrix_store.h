#ifndef __FM_SUB_MATRIX_STORE_H__
#define __FM_SUB_MATRIX_STORE_H__

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

#include <boost/format.hpp>

#include "matrix_store.h"
#include "EM_object.h"
#include "local_matrix_store.h"

namespace fm
{

namespace detail
{

class sub_matrix_store: public matrix_store, public EM_object
{
	matrix_store::const_ptr orig;
protected:
	const matrix_store &get_orig() const {
		return *orig;
	}
public:
	sub_matrix_store(matrix_store::const_ptr orig, size_t nrows,
			size_t ncols): matrix_store(nrows, ncols, orig->is_in_mem(),
				orig->get_type()) {
		this->orig = orig;
		assert(orig->get_num_rows() == nrows || orig->get_num_cols() == ncols);
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
		return orig->get_portion_node_id(id);
	}
	virtual int get_num_nodes() const {
		return orig->get_num_nodes();
	}

	virtual std::shared_ptr<local_matrix_store> get_portion(
			size_t start_row, size_t start_col, size_t num_rows,
			size_t num_cols) {
		// This doesn't need to be used. Changing the data in the local portion
		// doesn't affect the data in the disks.
		assert(0);
		return local_matrix_store::ptr();
	}
	virtual async_cres_t get_mem_portion_async(
			size_t start_row, size_t start_col, size_t num_rows,
			size_t num_cols, portion_compute::ptr compute) const = 0;
	virtual async_cres_t get_EM_portion_async(
			size_t start_row, size_t start_col, size_t num_rows,
			size_t num_cols, portion_compute::ptr compute) const = 0;
	virtual async_cres_t get_portion_async(
			size_t start_row, size_t start_col, size_t num_rows,
			size_t num_cols, portion_compute::ptr compute) const;
	virtual std::shared_ptr<const local_matrix_store> get_portion(
			size_t start_row, size_t start_col, size_t num_rows,
			size_t num_cols) const;

	virtual void write_portion_async(local_matrix_store::const_ptr portion,
			off_t start_row, off_t start_col) {
		BOOST_LOG_TRIVIAL(error)
			<< "Don't support write_portion_async in a sub EM matrix";
	}

	virtual std::vector<safs::io_interface::ptr> create_ios() const;
};

/*
 * This matrix store accesses a set of columns of a tall matrix.
 */
class sub_col_matrix_store: public sub_matrix_store
{
	const size_t mat_id;
	const data_id_t::ptr data_id;
	std::vector<off_t> rc_idxs;

public:
	sub_col_matrix_store(const std::vector<off_t> idxs,
			matrix_store::const_ptr store, data_id_t::ptr _data_id): sub_matrix_store(
				store, store->get_num_rows(), idxs.size()),
			mat_id(mat_counter++), data_id(_data_id) {
		this->rc_idxs = idxs;
		assert(!store->is_wide());
		assert(idxs.size() > 0);
	}

	sub_col_matrix_store(const std::vector<off_t> idxs,
			matrix_store::const_ptr store): sub_matrix_store(store,
				store->get_num_rows(), idxs.size()), mat_id(
				// The sub matrix has a different matrix ID and data ID from
				// its parent matrix.
				mat_counter++), data_id(data_id_t::create(mat_id)) {
		this->rc_idxs = idxs;
		assert(!store->is_wide());
		assert(idxs.size() > 0);
	}

	virtual void inc_dag_ref(size_t id) {
		data_id->inc_ref(id);
	}
	virtual void reset_dag_ref() {
		data_id->reset_ref();
	}
	virtual size_t get_dag_ref() const {
		return data_id->get_ref();
	}
	virtual size_t get_data_id() const {
		return data_id->get_id();
	}

	virtual matrix_store::const_ptr transpose() const;
	virtual std::string get_name() const {
		return (boost::format("sub_mat-%1%(%2%,%3%)") % mat_id
				% get_num_rows() % get_num_cols()).str();
	}

	virtual async_cres_t get_mem_portion_async(
			size_t start_row, size_t start_col, size_t num_rows,
			size_t num_cols, portion_compute::ptr compute) const;
	virtual async_cres_t get_EM_portion_async(
			size_t start_row, size_t start_col, size_t num_rows,
			size_t num_cols, portion_compute::ptr compute) const;

	virtual matrix_store::const_ptr get_cols(
			const std::vector<off_t> &idxs) const {
		assert(!is_wide());
		std::vector<off_t> orig_idxs(idxs.size());
		for (size_t i = 0; i < idxs.size(); i++)
			orig_idxs[i] = rc_idxs[idxs[i]];
		return get_orig().get_cols(orig_idxs);
	}
	virtual matrix_store::const_ptr get_rows(
			const std::vector<off_t> &idxs) const {
		assert(!is_wide());
		matrix_store::const_ptr rows = get_orig().get_rows(idxs);
		return rows->get_cols(rc_idxs);
	}
};

class sub_row_matrix_store: public sub_matrix_store
{
	const size_t mat_id;
	const data_id_t::ptr data_id;
	std::vector<off_t> rc_idxs;

public:
	sub_row_matrix_store(const std::vector<off_t> idxs,
			matrix_store::const_ptr store, data_id_t::ptr _data_id): sub_matrix_store(
				store, idxs.size(), store->get_num_cols()),
			mat_id(mat_counter++), data_id(_data_id) {
		this->rc_idxs = idxs;
		assert(store->is_wide());
		assert(idxs.size() > 0);
	}

	sub_row_matrix_store(const std::vector<off_t> idxs,
			matrix_store::const_ptr store): sub_matrix_store(store,
				idxs.size(), store->get_num_cols()), mat_id(
				// The sub matrix has a different matrix ID and data ID from
				// its parent matrix.
				mat_counter++), data_id(data_id_t::create(mat_id)) {
		this->rc_idxs = idxs;
		assert(store->is_wide());
		assert(idxs.size() > 0);
	}

	virtual void inc_dag_ref(size_t id) {
		data_id->inc_ref(id);
	}
	virtual void reset_dag_ref() {
		data_id->reset_ref();
	}
	virtual size_t get_dag_ref() const {
		return data_id->get_ref();
	}
	virtual size_t get_data_id() const {
		return data_id->get_id();
	}

	virtual matrix_store::const_ptr transpose() const;
	virtual std::string get_name() const {
		return (boost::format("sub_mat-%1%(%2%,%3%)") % mat_id
				% get_num_rows() % get_num_cols()).str();
	}

	virtual async_cres_t get_mem_portion_async(
			size_t start_row, size_t start_col, size_t num_rows,
			size_t num_cols, portion_compute::ptr compute) const;
	virtual async_cres_t get_EM_portion_async(
			size_t start_row, size_t start_col, size_t num_rows,
			size_t num_cols, portion_compute::ptr compute) const;

	virtual matrix_store::const_ptr get_cols(
			const std::vector<off_t> &idxs) const {
		assert(is_wide());
		matrix_store::const_ptr cols = get_orig().get_cols(idxs);
		return cols->get_rows(rc_idxs);
	}
	virtual matrix_store::const_ptr get_rows(
			const std::vector<off_t> &idxs) const {
		assert(is_wide());
		std::vector<off_t> orig_idxs(idxs.size());
		for (size_t i = 0; i < idxs.size(); i++)
			orig_idxs[i] = rc_idxs[idxs[i]];
		return get_orig().get_rows(orig_idxs);
	}
};

}

}

#endif
