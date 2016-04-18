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

#include "io_interface.h"

#include "sink_matrix.h"
#include "local_matrix_store.h"
#include "EM_object.h"
#include "dense_matrix.h"

namespace fm
{

namespace detail
{

namespace
{

class lmaterialize_col_matrix_store: public lvirtual_col_matrix_store
{
	std::vector<local_matrix_store::const_ptr> parts;
public:
	lmaterialize_col_matrix_store(
			const std::vector<local_matrix_store::const_ptr> &parts): lvirtual_col_matrix_store(
				parts[0]->get_global_start_row(), parts[0]->get_global_start_col(),
				parts[0]->get_num_rows(), parts[0]->get_num_cols(),
				parts[0]->get_type(), parts[0]->get_node_id()) {
		this->parts = parts;
	}

	virtual bool resize(off_t local_start_row, off_t local_start_col,
			size_t local_num_rows, size_t local_num_cols) {
		for (size_t i = 0; i < parts.size(); i++)
			const_cast<local_matrix_store &>(*parts[i]).resize(local_start_row,
					local_start_col, local_num_rows, local_num_cols);
		return local_matrix_store::resize(local_start_row, local_start_col,
				local_num_rows, local_num_cols);
	}
	virtual void reset_size() {
		for (size_t i = 0; i < parts.size(); i++)
			const_cast<local_matrix_store &>(*parts[i]).reset_size();
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
		for (size_t i = 0; i < parts.size(); i++)
			parts[i]->materialize_self();
	}
};

class lmaterialize_row_matrix_store: public lvirtual_row_matrix_store
{
	std::vector<local_matrix_store::const_ptr> parts;
public:
	lmaterialize_row_matrix_store(
			const std::vector<local_matrix_store::const_ptr> &parts): lvirtual_row_matrix_store(
				parts[0]->get_global_start_row(), parts[0]->get_global_start_col(),
				parts[0]->get_num_rows(), parts[0]->get_num_cols(),
				parts[0]->get_type(), parts[0]->get_node_id()) {
		this->parts = parts;
	}

	virtual bool resize(off_t local_start_row, off_t local_start_col,
			size_t local_num_rows, size_t local_num_cols) {
		for (size_t i = 0; i < parts.size(); i++)
			const_cast<local_matrix_store &>(*parts[i]).resize(local_start_row,
					local_start_col, local_num_rows, local_num_cols);
		return local_matrix_store::resize(local_start_row, local_start_col,
				local_num_rows, local_num_cols);
	}
	virtual void reset_size() {
		for (size_t i = 0; i < parts.size(); i++)
			const_cast<local_matrix_store &>(*parts[i]).reset_size();
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
		for (size_t i = 0; i < parts.size(); i++)
			parts[i]->materialize_self();
	}
};

/*
 * This class is used internally for materializing a group of blocks
 * that share the same underlying matrices. Once an object of the class
 * is created, we'll pass it to __mapply_portion. Thus, we only need to
 * implement a subset of its methods.
 */
class block_group: public virtual_matrix_store, public EM_object
{
	typedef std::vector<detail::matrix_store::const_ptr> store_vec_t;
	store_vec_t stores;
	const EM_object *obj;
public:
	block_group(const store_vec_t &stores): virtual_matrix_store(
			stores[0]->get_num_rows(), stores[0]->get_num_cols(),
			stores[0]->is_in_mem(), stores[0]->get_type()) {
		this->stores = stores;
		obj = dynamic_cast<const EM_object *>(stores[0].get());
		assert(obj);
	}

	/*
	 * These methods don't need to be implemented.
	 */

	virtual void materialize_self() const {
		assert(0);
	}

	virtual matrix_store::const_ptr materialize(bool in_mem,
		int num_nodes) const {
		assert(0);
		return matrix_store::const_ptr();
	}

	virtual vec_store::const_ptr get_col_vec(off_t idx) const {
		assert(0);
		return vec_store::const_ptr();
	}
	virtual vec_store::const_ptr get_row_vec(off_t idx) const {
		assert(0);
		return vec_store::const_ptr();
	}
	virtual matrix_store::const_ptr get_cols(const std::vector<off_t> &idxs) const {
		assert(0);
		return matrix_store::const_ptr();
	}
	virtual matrix_store::const_ptr get_rows(const std::vector<off_t> &idxs) const {
		assert(0);
		return matrix_store::const_ptr();
	}

	virtual matrix_store::const_ptr transpose() const {
		assert(0);
		return matrix_store::const_ptr();
	}

	/*
	 * These methods only indicate the information of the matrix.
	 */

	virtual int get_portion_node_id(size_t id) const {
		return stores[0]->get_portion_node_id(id);
	}

	virtual std::pair<size_t, size_t> get_portion_size() const {
		return stores[0]->get_portion_size();
	}

	virtual int get_num_nodes() const {
		return stores[0]->get_num_nodes();
	}

	virtual matrix_layout_t store_layout() const {
		return stores[0]->store_layout();
	}

	virtual std::string get_name() const {
		return stores[0]->get_name();
	}

	virtual std::unordered_map<size_t, size_t> get_underlying_mats() const {
		return stores[0]->get_underlying_mats();
	}

	virtual std::vector<safs::io_interface::ptr> create_ios() const {
		return obj->create_ios();
	}

	/*
	 * These are the main methods that we need to take care of.
	 */

	using virtual_matrix_store::get_portion;
	virtual local_matrix_store::const_ptr get_portion(size_t start_row,
			size_t start_col, size_t num_rows, size_t num_cols) const;
	virtual local_matrix_store::const_ptr get_portion(size_t id) const;
	using virtual_matrix_store::get_portion_async;
	virtual async_cres_t get_portion_async(
			size_t start_row, size_t start_col, size_t num_rows,
			size_t num_cols, std::shared_ptr<portion_compute> compute) const;
};

local_matrix_store::const_ptr create_local_store(matrix_layout_t store_layout,
		const std::vector<local_matrix_store::const_ptr> &portions)
{
	if (store_layout == matrix_layout_t::L_ROW)
		return local_matrix_store::const_ptr(
				new lmaterialize_row_matrix_store(portions));
	else
		return local_matrix_store::const_ptr(
				new lmaterialize_col_matrix_store(portions));
}

local_matrix_store::const_ptr block_group::get_portion(size_t start_row,
		size_t start_col, size_t num_rows, size_t num_cols) const
{
	std::vector<local_matrix_store::const_ptr> portions(stores.size());
	for (size_t i = 0; i < stores.size(); i++)
		portions[i] = stores[i]->get_portion(start_row, start_col,
				num_rows, num_cols);
	return create_local_store(store_layout(), portions);
}

local_matrix_store::const_ptr block_group::get_portion(size_t id) const
{
	std::vector<local_matrix_store::const_ptr> portions(stores.size());
	for (size_t i = 0; i < stores.size(); i++)
		portions[i] = stores[i]->get_portion(id);
	return create_local_store(store_layout(), portions);
}

/*
 * When a portion is read from diks, this portion compute is invoked.
 * It will be invoked multiple times. We need to make sure user's portion
 * compute is invoked only once.
 */
class collect_portion_compute: public portion_compute
{
	size_t num_expected;
	size_t num_reads;
	portion_compute::ptr orig_compute;
public:
	collect_portion_compute(portion_compute::ptr orig_compute,
			size_t num_expected) {
		this->num_expected = num_expected;
		this->num_reads = 0;
		this->orig_compute = orig_compute;
	}

	virtual void run(char *buf, size_t size) {
		num_reads++;
		if (num_reads == num_expected) {
			orig_compute->run(NULL, 0);
			// This only runs once.
			// Let's remove all user's portion compute to indicate that it has
			// been invoked.
			orig_compute = NULL;
		}
	}
};

async_cres_t block_group::get_portion_async(size_t start_row, size_t start_col,
		size_t num_rows, size_t num_cols, portion_compute::ptr compute) const
{
	std::vector<local_matrix_store::const_ptr> portions(stores.size());
	bool avail = false;
	portion_compute::ptr collect_compute(
			new collect_portion_compute(compute, stores.size()));
	for (size_t i = 0; i < stores.size(); i++) {
		auto ret = stores[i]->get_portion_async(start_row, start_col,
				num_rows, num_cols, collect_compute);
		portions[i] = ret.second;
		// These portions are all from the same underlying matrices, they
		// should all be available or unavailable at the same time.
		if (i == 0)
			avail = ret.first;
		else
			assert(avail == ret.first);
	}
	return async_cres_t(avail, create_local_store(store_layout(), portions));
}

}

std::vector<matrix_store::const_ptr> reorg_block_sinks(
		const std::vector<block_sink_store::const_ptr> &sinks)
{
	size_t num_blocks = sinks[0]->get_num_blocks();
	std::vector<detail::matrix_store::const_ptr> block_groups(num_blocks);
	for (size_t i = 0; i < num_blocks; i++) {
		// block i from all sink matrices.
		std::vector<detail::matrix_store::const_ptr> blocks(sinks.size());
		for (size_t j = 0; j < sinks.size(); j++)
			blocks[j] = sinks[j]->get_block(i);
		for (size_t j = 1; j < blocks.size(); j++) {
			// The blocks are from the block matrices. They should all be
			// the same.
			assert(blocks[j]->get_num_rows() == blocks[0]->get_num_rows());
			assert(blocks[j]->get_num_cols() == blocks[0]->get_num_cols());
			assert(blocks[j]->is_in_mem() == blocks[0]->is_in_mem());
		}
		block_groups[i] = detail::matrix_store::const_ptr(new block_group(blocks));
	}
	return block_groups;
}

typedef std::vector<block_sink_store::const_ptr> sink_vec_t;

std::vector<sink_vec_t> group_block_sinks(const sink_vec_t &sinks)
{
	// I assume the number of block sink matrices is small.
	// Maybe we can use hashing to reduce the number of groups we should
	// search through. For now we only need scan all groups to figure out
	// which group a sink matrix belongs to.
	std::vector<sink_vec_t> ret;
	for (size_t i = 0; i < sinks.size(); i++) {
		block_sink_store::const_ptr store = sinks[i];
		// Let's search through all groups and see which group does this sink
		// matrix belongs to.
		for (size_t j = 0; j < ret.size(); j++) {
			if (ret[j].front()->match(*store)) {
				ret[j].push_back(store);
				store = NULL;
			}
		}

		// This sink matrix doesn't belong to any group, let's create
		// a new group.
		if (store)
			ret.push_back(sink_vec_t(1, store));
	}

	return ret;
}

static size_t get_num_rows(const std::vector<matrix_store::const_ptr> &stores)
{
	if (stores[0]->is_wide()) {
		size_t num_rows = 0;
		for (size_t i = 0; i < stores.size(); i++)
			num_rows += stores[i]->get_num_rows();
		return num_rows;
	}
	else
		return stores[0]->get_num_rows();
}

static size_t get_num_cols(const std::vector<matrix_store::const_ptr> &stores)
{
	if (stores[0]->is_wide())
		return stores[0]->get_num_cols();
	else {
		size_t num_cols = 0;
		for (size_t i = 0; i < stores.size(); i++)
			num_cols += stores[i]->get_num_cols();
		return num_cols;
	}
}

block_sink_store::block_sink_store(
		const std::vector<matrix_store::const_ptr> &stores): virtual_matrix_store(
			detail::get_num_rows(stores), detail::get_num_cols(stores),
			stores[0]->is_in_mem(), stores[0]->get_type())
{
	this->stores = stores;
	under_mats.resize(stores.size());
	// Collect the underlying matrices for each input matrix.
	for (size_t i = 0; i < stores.size(); i++) {
		auto ret = stores[i]->get_underlying_mats();
		for (auto it = ret.begin(); it != ret.end(); it++)
			under_mats[i].push_back(it->first);
	}
}

/*
 * In order to have two block sink matrices materialized together,
 * all of the following conditions should be met:
 * * both matrices should the same number of blocks.
 * * each block should have the same underlying matrices.
 */
bool block_sink_store::match(const block_sink_store &store) const
{
	if (under_mats.size() != store.under_mats.size())
		return false;

	for (size_t i = 0; i < under_mats.size(); i++) {
		if (under_mats[i].size() != store.under_mats[i].size())
			return false;
		for (size_t j = 0; j < under_mats[i].size(); j++)
			if (under_mats[i][j] != store.under_mats[i][j])
				return false;
	}
	return true;
}

std::vector<matrix_store::const_ptr> block_sink_store::get_materialized_blocks() const
{
	std::vector<dense_matrix::ptr> mats(stores.size());
	for (size_t i = 0; i < mats.size(); i++)
		mats[i] = dense_matrix::create(stores[i]);
	// TODO we may need to disable caching on the matrices.
	bool ret = fm::materialize(mats, false);
	if (!ret)
		return std::vector<matrix_store::const_ptr>();

	std::vector<matrix_store::const_ptr> stores(mats.size());
	for (size_t i = 0; i < stores.size(); i++)
		stores[i] = mats[i]->get_raw_store();
	return stores;
}

vec_store::const_ptr block_sink_store::get_col_vec(off_t idx) const
{
	assert(0);
	return vec_store::const_ptr();
}

vec_store::const_ptr block_sink_store::get_row_vec(off_t idx) const
{
	assert(0);
	return vec_store::const_ptr();
}

matrix_store::const_ptr block_sink_store::get_cols(
		const std::vector<off_t> &idxs) const
{
	assert(0);
	return matrix_store::const_ptr();
}

matrix_store::const_ptr block_sink_store::get_rows(
		const std::vector<off_t> &idxs) const
{
	assert(0);
	return matrix_store::const_ptr();
}

local_matrix_store::const_ptr block_sink_store::get_portion(
			size_t start_row, size_t start_col, size_t num_rows,
			size_t num_cols) const
{
	assert(0);
	return local_matrix_store::const_ptr();
}

local_matrix_store::const_ptr block_sink_store::get_portion(size_t id) const
{
	assert(0);
	return local_matrix_store::const_ptr();
}

async_cres_t block_sink_store::get_portion_async(size_t start_row,
		size_t start_col, size_t num_rows, size_t num_cols,
		std::shared_ptr<portion_compute> compute) const
{
	assert(0);
	return async_cres_t();
}

matrix_store::const_ptr block_sink_store::transpose() const
{
	assert(0);
	return matrix_store::const_ptr();
}

int block_sink_store::get_portion_node_id(size_t id) const
{
	assert(0);
	return -1;
}

std::pair<size_t, size_t> block_sink_store::get_portion_size() const
{
	assert(0);
	return std::pair<size_t, size_t>();
}

int block_sink_store::get_num_nodes() const
{
	return stores[0]->get_num_nodes();
}

matrix_layout_t block_sink_store::store_layout() const
{
	// The matrix either contains a vector or a single element.
	// In either case, it should be stored in column-major order.
	return matrix_layout_t::L_COL;
}

std::string block_sink_store::get_name() const
{
	assert(0);
	return std::string();
}

std::unordered_map<size_t, size_t> block_sink_store::get_underlying_mats() const
{
	assert(0);
	return std::unordered_map<size_t, size_t>();
}

}

}
