#ifndef __FM_CACHED_MATRIX_STORE_H__
#define __FM_CACHED_MATRIX_STORE_H__

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

#include "matrix_store.h"
#include "EM_object.h"
#include "EM_dense_matrix.h"

namespace fm
{

namespace detail
{

/*
 * This matrix caches part of the underlying matrix.
 * It always caches the first few columns of a tall matrix and
 * the first few rows of a wide matrix.
 */
class cached_matrix_store: public matrix_store, public EM_object
{
	// The in-mem part.
	matrix_store::const_ptr cached;
	// The hybrid matrix store (in-mem + external-memory)
	matrix_store::const_ptr mixed;
	// The external-memory matrix.
	EM_matrix_store::const_ptr em_store;

	// The members below are only used for writing data.
	// It caches the first few vectors.
	matrix_store::ptr cached_buf;
	// The external-memory matrix store
	EM_matrix_store::ptr em_buf;

	cached_matrix_store(size_t num_rows, size_t num_cols,
			const scalar_type &type);
	cached_matrix_store(size_t num_rows, size_t num_cols, int num_nodes,
			const scalar_type &type, size_t num_cached_vecs,
			matrix_layout_t cached_layout, matrix_layout_t em_layout);
public:
	typedef std::shared_ptr<cached_matrix_store> ptr;
	typedef std::shared_ptr<const cached_matrix_store> const_ptr;

	static ptr create(size_t num_rows, size_t num_cols, int num_nodes,
			const scalar_type &type, size_t num_cached_vecs,
			matrix_layout_t cached_layout,
			matrix_layout_t em_layout = matrix_layout_t::L_NONE);

	static void drop_all_cache();

	virtual size_t get_data_id() const {
		return em_store->get_data_id();
	}

	virtual std::unordered_map<size_t, size_t> get_underlying_mats() const {
		return em_store->get_underlying_mats();
	}

	virtual matrix_layout_t store_layout() const {
		return em_store->store_layout();
	}

	virtual std::pair<size_t, size_t> get_portion_size() const {
		return em_store->get_portion_size();
	}

	virtual std::string get_name() const;

	virtual matrix_store::const_ptr transpose() const;

	virtual async_cres_t get_portion_async(size_t start_row, size_t start_col,
			size_t num_rows, size_t num_cols,
			std::shared_ptr<portion_compute> compute) const {
		return mixed->get_portion_async(start_row, start_col, num_rows,
				num_cols, compute);
	}
	virtual std::shared_ptr<const local_matrix_store> get_portion(
			size_t start_row, size_t start_col, size_t num_rows,
			size_t num_cols) const {
		return mixed->get_portion(start_row, start_col, num_rows, num_cols);
	}
	virtual std::shared_ptr<local_matrix_store> get_portion(
			size_t start_row, size_t start_col, size_t num_rows,
			size_t num_cols) {
		return std::shared_ptr<local_matrix_store>();
	}
	virtual int get_portion_node_id(size_t id) const {
		// The in-memory data may not be in the same NUMA node.
		return -1;
	}
	virtual void write_portion_async(
			std::shared_ptr<const local_matrix_store> portion,
			off_t start_row, off_t start_col);

	virtual matrix_store::const_ptr get_rows(
			const std::vector<off_t> &idxs) const;

	virtual void set_prefetches(size_t num, std::pair<size_t, size_t> range);
	virtual void set_cache_portion(bool cache_portion);
	virtual bool is_cache_portion() const;

	size_t get_num_cached_vecs() const {
		if (cached == NULL)
			return 0;
		else if (is_wide())
			return cached->get_num_rows();
		else
			return cached->get_num_cols();
	}

	void drop_cache() {
		cached = NULL;
		mixed = em_store;
	}

	EM_matrix_store::const_ptr get_underlying() const {
		return em_store;
	}

	virtual std::vector<safs::io_interface::ptr> create_ios() const {
		return em_store->create_ios();
	}
};

}

}

#endif
