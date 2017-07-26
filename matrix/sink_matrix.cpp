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
#include "materialize.h"
#include "mem_matrix_store.h"

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

std::unordered_map<size_t, sink_store::const_ptr> sink_store::sinks;

void sink_store::register_sink_matrices(sink_store::const_ptr store)
{
	if (!store->has_materialized())
		sinks.insert(std::pair<size_t, sink_store::const_ptr>(
					store->get_data_id(), store));
}

static inline bool share_io(const std::unordered_map<size_t, size_t> &underlying1,
		const std::vector<size_t> &underlying2)
{
	for (size_t i = 0; i < underlying2.size(); i++) {
		auto it = underlying1.find(underlying2[i]);
		if (INVALID_MAT_ID != underlying2[i] && it != underlying1.end())
			return true;
	}
	return false;
}

void sink_store::materialize_matrices(virtual_matrix_store::const_ptr store)
{
	// Get the underlying matrices from the virtual matrix.
	auto underlying = store->get_underlying_mats();
	std::vector<size_t> underlying_ids;
	for (auto it = underlying.begin(); it != underlying.end(); it++)
		underlying_ids.push_back(it->first);

	// All sink matrices in `sinks' are individual sink matrices, instead of
	// block sink matrices.
	std::unordered_map<size_t, sink_store::const_ptr> share_io_sinks;
	std::vector<sink_store::const_ptr> to_remove;
	for (auto it = sinks.begin(); it != sinks.end(); it++) {
		sink_store::const_ptr sink = it->second;
		// We should remove the matrices that have been materialized.
		bool materialized = sink->has_materialized();
		if (materialized)
			to_remove.push_back(sink);
		// We only select the matrices that haven't been materialized.
		if (!materialized && share_io(sink->get_underlying_mats(),
					underlying_ids)) {
			share_io_sinks.insert(std::pair<size_t, sink_store::const_ptr>(
						sink->get_data_id(), sink));
			to_remove.push_back(sink);
		}
	}
	for (auto it = to_remove.begin(); it != to_remove.end(); it++)
		sinks.erase((*it)->get_data_id());
	to_remove.clear();

	if (share_io_sinks.empty())
		store->materialize_self();
	else {
		std::vector<dense_matrix::ptr> mats;
		for (auto it = share_io_sinks.begin(); it != share_io_sinks.end(); it++) {
			if (!store->share_data(*it->second) && !it->second->has_materialized())
				mats.push_back(dense_matrix::create(it->second));
		}

		// The input matrix might be a block sink matrix.
		// If it is, we need to search among the individual blocks in
		// the block sink matrix and avoid materializing an individual
		// sink matrix multiple times.
		block_sink_store::const_ptr block_sink
			= std::dynamic_pointer_cast<const block_sink_store>(store);
		if (block_sink == NULL)
			mats.push_back(dense_matrix::create(store));
		else {
			const std::vector<sink_store::const_ptr> &blocks
				= block_sink->get_stores();
			for (size_t i = 0; i < blocks.size(); i++) {
				if (blocks[i] == NULL)
					continue;
				// We only need to materialize the individual matrices
				// that don't exist in `share_io_sinks'.
				auto it = share_io_sinks.find(blocks[i]->get_data_id());
				if (it == share_io_sinks.end())
					mats.push_back(dense_matrix::create(blocks[i]));
			}
		}

		fm::materialize(mats, true, false);
	}

	for (auto it = share_io_sinks.begin(); it != share_io_sinks.end(); it++)
		assert(it->second->has_materialized());
}

static size_t get_num_rows(const std::vector<sink_store::const_ptr> &stores,
		size_t num_block_rows, size_t num_block_cols)
{
	size_t num_rows = 0;
	for (size_t i = 0; i < num_block_rows; i++)
		num_rows += stores[i * num_block_cols]->get_num_rows();
	return num_rows;
}

static size_t get_num_cols(const std::vector<sink_store::const_ptr> &stores,
		size_t num_block_rows, size_t num_block_cols)
{
	size_t num_cols = 0;
	for (size_t i = 0; i < num_block_cols; i++)
		num_cols += stores[i]->get_num_cols();
	return num_cols;
}

block_sink_store::ptr block_sink_store::create(
		const std::vector<matrix_store::const_ptr> &stores,
		size_t num_block_rows, size_t num_block_cols)
{
	std::vector<sink_store::const_ptr> sink_stores(stores.size());
	size_t num_stores = 0;
	for (size_t i = 0; i < stores.size(); i++) {
		if (stores[i] == NULL)
			continue;

		num_stores++;
		sink_stores[i] = std::dynamic_pointer_cast<const sink_store>(stores[i]);
		// The input matrices have to be sink matrices.
		assert(sink_stores[i]);
	}
	if (num_stores != num_block_rows * num_block_cols)
		return block_sink_store::ptr();

	// We don't need to register a block sink matrix to be materialized
	// whenever possible because all individual sink matrices have been
	// registered for materialization.
	return block_sink_store::ptr(new block_sink_store(sink_stores,
				num_block_rows, num_block_cols, false));
}

static inline size_t sum(const std::vector<size_t> &v)
{
	return std::accumulate(v.begin(), v.end(), 0);
}

block_sink_store::ptr block_sink_store::create(
		const std::vector<size_t> &nrow_in_blocks,
		const std::vector<size_t> &ncol_in_blocks,
		bool in_mem, const scalar_type &type, bool is_sym)
{
	if (is_sym && sum(nrow_in_blocks) != sum(ncol_in_blocks))
		return block_sink_store::ptr();

	// We don't need to register a block sink matrix to be materialized
	// whenever possible because all individual sink matrices have been
	// registered for materialization.
	return block_sink_store::ptr(new block_sink_store(nrow_in_blocks,
				ncol_in_blocks, in_mem, type, is_sym));
}

block_sink_store::block_sink_store(const std::vector<size_t> &nrow_in_blocks,
		const std::vector<size_t> &ncol_in_blocks, bool in_mem,
		const scalar_type &type, bool is_sym): sink_store(sum(nrow_in_blocks),
			sum(ncol_in_blocks), in_mem, type)
{
	stores.resize(nrow_in_blocks.size() * ncol_in_blocks.size());
	this->nrow_in_blocks = nrow_in_blocks;
	this->ncol_in_blocks = ncol_in_blocks;
	this->is_sym = is_sym;
	this->underlying = get_underlying_mats();
}

block_sink_store::block_sink_store(
		// I assume all matrices are kept in row-major order.
		const std::vector<sink_store::const_ptr> &stores,
		size_t num_block_rows, size_t num_block_cols, bool is_sym): sink_store(
			detail::get_num_rows(stores, num_block_rows,
				num_block_cols), detail::get_num_cols(stores, num_block_rows,
				num_block_cols), stores[0]->is_in_mem(), stores[0]->get_type())
{
	// If we have all sink matrices, we don't care about this flag.
	this->is_sym = is_sym;
	this->stores = stores;
	nrow_in_blocks.resize(num_block_rows);
	ncol_in_blocks.resize(num_block_cols);

	// all matrices in the block row should have the same number of rows.
	for (size_t i = 0; i < num_block_rows; i++) {
		size_t num_rows = stores[get_idx(i, 0)]->get_num_rows();
		nrow_in_blocks[i] = num_rows;
		for (size_t j = 1; j < num_block_cols; j++)
			assert(stores[get_idx(i, j)]
					&& stores[get_idx(i, j)]->get_num_rows() == num_rows);
	}
	// all matrices in the block col should have the same number of cols.
	for (size_t i = 0; i < num_block_cols; i++) {
		size_t num_cols = stores[i]->get_num_cols();
		ncol_in_blocks[i] = num_cols;
		for (size_t j = 1; j < num_block_rows; j++)
			assert(stores[get_idx(j, i)]
					&& stores[get_idx(j, i)]->get_num_cols() == num_cols);
	}
	assert(num_block_rows * num_block_cols == stores.size());
	this->underlying = get_underlying_mats();
}

matrix_store::const_ptr block_sink_store::transpose() const
{
	if (result)
		return result->transpose();

	// If this is a symmetric matrix, the transpose will be the same
	// as the original matrix.
	if (is_sym) {
		block_sink_store *new_store = new block_sink_store(*this);
		return matrix_store::const_ptr(new_store);
	}
	else {
		block_sink_store::ptr new_sink = block_sink_store::create(ncol_in_blocks,
				nrow_in_blocks, is_in_mem(), get_type(), is_sym);
		for (size_t j = 0; j < get_num_block_cols(); j++)
			for (size_t i = 0; i < get_num_block_rows(); i++) {
				auto store = stores[get_idx(i, j)];
				if (store == NULL)
					continue;

				matrix_store::const_ptr t = store->transpose();
				sink_store::const_ptr sink
					= std::dynamic_pointer_cast<const sink_store>(t);
				assert(sink);
				new_sink->stores[new_sink->get_idx(j, i)] = sink;
			}
		return new_sink;
	}
}

std::string block_sink_store::get_name() const
{
	std::string str = "block_sink(";
	for (size_t i = 0; i < stores.size(); i++) {
		if (stores[i]) {
			str += stores[i]->get_name();
			if (i < stores.size() - 1)
				str += ", ";
		}
	}
	str += ")";
	return str;
}

std::unordered_map<size_t, size_t> block_sink_store::get_underlying_mats() const
{
	if (has_materialized())
		return std::unordered_map<size_t, size_t>();
	if (!this->underlying.empty())
		return this->underlying;

	std::unordered_map<size_t, size_t> underlying;
	for (size_t i = 0; i < stores.size(); i++) {
		// Some of the stores might be empty. For example, crossprod with
		// itself only needs to store half of the crossprod result on
		// individual matrices.
		if (stores[i] == NULL)
			continue;

		auto tmp = stores[i]->get_underlying_mats();
		underlying.insert(tmp.begin(), tmp.end());
	}
	return underlying;
}

bool block_sink_store::has_materialized() const
{
	if (result)
		return true;

	for (size_t i = 0; i < stores.size(); i++) {
		if (stores[i] && !stores[i]->has_materialized())
			return false;
	}
	return true;
}

matrix_store::const_ptr block_sink_store::get_result() const
{
	if (result == NULL)
		materialize_self();
	return result;
}

std::vector<virtual_matrix_store::const_ptr> block_sink_store::get_compute_matrices() const
{
	std::vector<virtual_matrix_store::const_ptr> ret;
	for (size_t i = 0; i < stores.size(); i++) {
		if (stores[i]) {
			auto tmp = stores[i]->get_compute_matrices();
			ret.insert(ret.end(), tmp.begin(), tmp.end());
		}
	}
	return ret;
}

void block_sink_store::materialize_self() const
{
	if (result)
		return;

	std::vector<dense_matrix::ptr> mat_vec;
	std::vector<dense_matrix::ptr> mat_map(stores.size());
	for (size_t i = 0; i < stores.size(); i++) {
		if (stores[i]) {
			dense_matrix::ptr mat = dense_matrix::create(stores[i]);
			if (!stores[i]->has_materialized())
				mat_vec.push_back(mat);
			mat_map[i] = mat;
		}
	}
	// We don't want to materialize all block matrices together because it
	// might consume a lot of memory.
	bool ret =fm::materialize(mat_vec, false);
	assert(ret);

	mem_matrix_store::ptr res = mem_matrix_store::create(get_num_rows(),
			get_num_cols(), store_layout(), get_type(), -1);
	off_t start_row = 0;
	for (size_t i = 0; i < get_num_block_rows(); i++) {
		off_t start_col = 0;
		for (size_t j = 0; j < get_num_block_cols(); j++) {
			auto mat = mat_map[get_idx(i, j)];
			if (mat == NULL) {
				mat = mat_map[get_idx(j, i)];
				assert(mat);
				mat = mat->transpose();
			}
			auto res_store = mat->get_raw_store();
			local_matrix_store::const_ptr tmp_portion
				= res_store->get_portion(0);
			assert(tmp_portion->get_num_rows() == res_store->get_num_rows());
			assert(tmp_portion->get_num_cols() == res_store->get_num_cols());

			local_matrix_store::ptr res_portion = res->get_portion(start_row,
					start_col, res_store->get_num_rows(),
					res_store->get_num_cols());
			res_portion->copy_from(*tmp_portion);
			start_col += res_store->get_num_cols();
		}
		start_row += nrow_in_blocks[i];
	}
	const_cast<block_sink_store *>(this)->result = res;
}

matrix_store::const_ptr block_sink_store::materialize(bool in_mem,
		int num_nodes) const
{
	materialize_self();
	return result;
}

void block_sink_store::set_store(size_t i, size_t j, matrix_store::const_ptr store)
{
	if (nrow_in_blocks[i] == 0)
		nrow_in_blocks[i] = store->get_num_rows();
	if (ncol_in_blocks[j] == 0)
		ncol_in_blocks[j] = store->get_num_cols();
	assert(nrow_in_blocks[i] == store->get_num_rows());
	assert(ncol_in_blocks[j] == store->get_num_cols());
	assert(store->is_in_mem() == is_in_mem());
	stores[get_idx(i, j)]
		= std::dynamic_pointer_cast<const sink_store>(store);
	assert(stores[get_idx(i, j)]);
}

matrix_store::const_ptr block_sink_store::get_store(size_t i, size_t j) const
{
	matrix_store::const_ptr ret = stores[get_idx(i, j)];
	if (ret)
		return ret;

	ret = stores[get_idx(j, i)];
	if (is_sym && ret)
		return ret->transpose();
	else
		return matrix_store::const_ptr();
}

static bool is_all_in_mem(
		const std::vector<detail::matrix_store::const_ptr> &stores) {
	for (size_t i = 0; i < stores.size(); i++)
		if (!stores[i]->is_in_mem())
			return false;
	return true;
}

mapply_sink_store::mapply_sink_store(
		const std::vector<matrix_store::const_ptr> &stores,
		portion_mapply_op::const_ptr op): sink_store(
			stores[0]->get_num_rows(), stores[0]->get_num_cols(),
			is_all_in_mem(stores), op->get_output_type())
{
	this->stores = stores;
	this->op = op;
}

mapply_sink_store::ptr mapply_sink_store::create(
		const std::vector<matrix_store::const_ptr> &stores,
		portion_mapply_op::const_ptr op)
{
	for (size_t i = 1; i < stores.size(); i++) {
		if (stores[0]->get_num_rows() != stores[i]->get_num_rows()
				|| stores[0]->get_num_cols() != stores[i]->get_num_cols()) {
			BOOST_LOG_TRIVIAL(error)
				<< "The input matrices don't have the same shape";
			return mapply_sink_store::ptr();
		}
	}
	return ptr(new mapply_sink_store(stores, op));
}

bool mapply_sink_store::has_materialized() const
{
	for (size_t i = 0; i < stores.size(); i++) {
		if (!stores[i]->is_virtual())
			continue;
		auto vmat = std::dynamic_pointer_cast<const virtual_matrix_store>(stores[i]);
		assert(vmat);
		if (!vmat->has_materialized())
			return false;
	}
	return true;
}

matrix_store::const_ptr mapply_sink_store::get_result() const
{
	if (result)
		return result;
	materialize_self();
	return result;
}

std::vector<virtual_matrix_store::const_ptr> mapply_sink_store::get_compute_matrices() const
{
	std::vector<virtual_matrix_store::const_ptr> ret;
	for (size_t i = 0; i < stores.size(); i++) {
		if (stores[i]->is_virtual()) {
			auto smat = std::dynamic_pointer_cast<const sink_store>(stores[i]);
			if (smat == NULL)
				continue;
			auto tmp = smat->get_compute_matrices();
			ret.insert(ret.end(), tmp.begin(), tmp.end());
		}
	}
	return ret;
}

void mapply_sink_store::materialize_self() const
{
	if (result)
		return;

	mem_matrix_store::ptr ret = mem_matrix_store::create(
			get_num_rows(), get_num_cols(), store_layout(), get_type(), -1);
	std::vector<local_matrix_store::const_ptr> ins(stores.size());
	for (size_t i = 0; i < ins.size(); i++) {
		ins[i] = stores[i]->get_portion(0);
		assert(ins[i]->get_num_rows() == stores[i]->get_num_rows());
		assert(ins[i]->get_num_cols() == stores[i]->get_num_cols());
	}
	local_matrix_store::ptr out = ret->get_portion(0);
	assert(out->get_num_rows() == ret->get_num_rows());
	assert(out->get_num_cols() == ret->get_num_cols());
	op->run(ins, *out);
	const_cast<mapply_sink_store *>(this)->result = ret;
}

matrix_store::const_ptr mapply_sink_store::materialize(bool in_mem, int num_nodes) const
{
	materialize_self();
	return result;
}

matrix_store::const_ptr mapply_sink_store::transpose() const
{
	if (result)
		return result->transpose();

	std::vector<matrix_store::const_ptr> tstores(stores.size());
	for (size_t i = 0; i < stores.size(); i++)
		tstores[i] = stores[i]->transpose();
	portion_mapply_op::const_ptr top = op->transpose();
	return mapply_sink_store::create(tstores, top);
}

std::unordered_map<size_t, size_t> mapply_sink_store::get_underlying_mats() const
{
	auto ret = stores[0]->get_underlying_mats();
	for (size_t i = 1; i < stores.size(); i++) {
		auto tmp = stores[i]->get_underlying_mats();
		ret.insert(tmp.begin(), tmp.end());
	}
	return ret;
}

void block_sink_store::inc_dag_ref(size_t id)
{
	for (size_t i = 0; i < stores.size(); i++)
		if (stores[i])
			const_cast<sink_store &>(*stores[i]).inc_dag_ref(get_data_id());
	// We don't need to increase the ref count of a sink matrix
	// because we never get a portion from a sink matrix.
}

void block_sink_store::reset_dag_ref()
{
	for (size_t i = 0; i < stores.size(); i++)
		if (stores[i])
			const_cast<sink_store &>(*stores[i]).reset_dag_ref();
}

void mapply_sink_store::inc_dag_ref(size_t id)
{
	for (size_t i = 0; i < stores.size(); i++)
		if (stores[i])
			const_cast<matrix_store &>(*stores[i]).inc_dag_ref(get_data_id());
	// We don't need to increase the ref count of a sink matrix
	// because we never get a portion from a sink matrix.
}

void mapply_sink_store::reset_dag_ref()
{
	for (size_t i = 0; i < stores.size(); i++)
		if (stores[i])
			const_cast<matrix_store &>(*stores[i]).reset_dag_ref();
}

}

}
