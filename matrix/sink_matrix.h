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

#ifndef __FM_SINK_MATRIX_H__
#define __FM_SINK_MATRIX_H__

#include "virtual_matrix_store.h"

namespace fm
{

namespace detail
{

class sink_store: public virtual_matrix_store
{
public:
	typedef std::shared_ptr<sink_store> ptr;
	typedef std::shared_ptr<const sink_store> const_ptr;

	sink_store(size_t nrow, size_t ncol, bool in_mem,
			const scalar_type &type): virtual_matrix_store(nrow, ncol,
				in_mem, type) {
	}

	// This returns the materialized result of the sink matrix.
	virtual matrix_store::const_ptr get_result() const = 0;
	// This returns a computation matrix that can be used for materialization
	// with other matrices.
	virtual virtual_matrix_store::const_ptr get_compute_matrix() const = 0;

	virtual matrix_store::const_ptr get_cols(const std::vector<off_t> &idxs) const {
		materialize_self();
		return get_result()->get_cols(idxs);
	}
	virtual matrix_store::const_ptr get_rows(const std::vector<off_t> &idxs) const {
		materialize_self();
		return get_result()->get_rows(idxs);
	}

	using virtual_matrix_store::get_portion;
	virtual std::shared_ptr<const local_matrix_store> get_portion(
			size_t start_row, size_t start_col, size_t num_rows,
			size_t num_cols) const {
		materialize_self();
		return get_result()->get_portion(start_row, start_col, num_rows,
				num_cols);
	}
	virtual std::shared_ptr<const local_matrix_store> get_portion(
			size_t id) const {
		materialize_self();
		return get_result()->get_portion(id);
	}
	using virtual_matrix_store::get_portion_async;
	virtual async_cres_t get_portion_async(
			size_t start_row, size_t start_col, size_t num_rows,
			size_t num_cols, std::shared_ptr<portion_compute> compute) const {
		materialize_self();
		return get_result()->get_portion_async(start_row, start_col,
				num_rows, num_cols, compute);
	}

	virtual int get_portion_node_id(size_t id) const {
		return -1;
	}

	virtual std::pair<size_t, size_t> get_portion_size() const {
		return std::pair<size_t, size_t>(mem_matrix_store::CHUNK_SIZE,
				mem_matrix_store::CHUNK_SIZE);
	}

	virtual int get_num_nodes() const {
		return -1;
	}

	virtual std::vector<safs::io_interface::ptr> create_ios() const {
		return std::vector<safs::io_interface::ptr>();
	}
};

/*
 * This represents the computation in a sink matrix.
 * This class is used as a proxy so that we can materialize this sink matrix
 * with other matrices together. Only the methods related to getting portions
 * of the input matrices matter.
 */
class sink_compute_store: public virtual_matrix_store
{
public:
	sink_compute_store(size_t nrow, size_t ncol, bool in_mem,
			const scalar_type &type): virtual_matrix_store(nrow, ncol,
				in_mem, type) {
	}

	bool has_materialized() const {
		return false;
	}

	virtual void materialize_self() const {
	}

	virtual matrix_store::const_ptr materialize(bool in_mem,
		int num_nodes) const {
		assert(0);
		return matrix_store::const_ptr();
	}

	virtual matrix_store::const_ptr transpose() const {
		assert(0);
		return matrix_store::const_ptr();
	}

	virtual std::string get_name() const {
		return "sink_compute_store";
	}
};

class block_sink_store: public virtual_matrix_store
{
	std::vector<matrix_store::const_ptr> stores;
	std::vector<std::vector<size_t> > under_mats;
public:
	typedef std::shared_ptr<block_sink_store> ptr;
	typedef std::shared_ptr<const block_sink_store> const_ptr;

	block_sink_store(const std::vector<matrix_store::const_ptr> &stores);
	/*
	 * Test if the two block sink matrices can be assigned to the same group
	 * for materialization.
	 */
	bool match(const block_sink_store &store) const;

	size_t get_num_blocks() const {
		return stores.size();
	}

	matrix_store::const_ptr get_block(size_t idx) const {
		return stores[idx];
	}

	std::vector<matrix_store::const_ptr> get_materialized_blocks() const;

	virtual matrix_store::const_ptr get_cols(const std::vector<off_t> &idxs) const;
	virtual matrix_store::const_ptr get_rows(const std::vector<off_t> &idxs) const;

	using virtual_matrix_store::get_portion;
	virtual std::shared_ptr<const local_matrix_store> get_portion(
			size_t start_row, size_t start_col, size_t num_rows,
			size_t num_cols) const;
	virtual std::shared_ptr<const local_matrix_store> get_portion(
			size_t id) const;
	using virtual_matrix_store::get_portion_async;
	virtual async_cres_t get_portion_async(
			size_t start_row, size_t start_col, size_t num_rows,
			size_t num_cols, std::shared_ptr<portion_compute> compute) const;

	virtual matrix_store::const_ptr transpose() const;

	virtual int get_portion_node_id(size_t id) const;
	virtual std::pair<size_t, size_t> get_portion_size() const;
	virtual int get_num_nodes() const;

	virtual matrix_layout_t store_layout() const;

	virtual std::string get_name() const;
	virtual std::unordered_map<size_t, size_t> get_underlying_mats() const;
};

/*
 * Here we assume all input block sink matrices run on the same block matrix.
 * As such, each block i from all sink matrices shares the same underlying
 * matrices. We reorganize the blocks and wrap block i from all matrices in
 * a virtual matrix. The idea is that when we materialize a virtual matrix,
 * all blocks that shared the same underlying matrices are materialized
 * together to reduce data movement between SSDs and main memory as well as
 * between main memory and CPU cache.
 */
std::vector<matrix_store::const_ptr> reorg_block_sinks(
		const std::vector<block_sink_store::const_ptr> &sinks);

/*
 * This function searches among all the block sink matrices and clusters them
 * based on their underlying matrices and the number of blocks in each of
 * the block sink matrices.
 */
std::vector<std::vector<block_sink_store::const_ptr> > group_block_sinks(
		const std::vector<block_sink_store::const_ptr> &sinks);

}

}

#endif
