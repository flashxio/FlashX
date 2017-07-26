#ifndef __MAPPLY_MATRIX_STORE_H__
#define __MAPPLY_MATRIX_STORE_H__

/*
 * Copyright 2015 Open Connectome Project (http://openconnecto.me)
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

#include "virtual_matrix_store.h"
#include "mem_matrix_store.h"
#include "EM_object.h"
#include "materialize.h"

namespace fm
{

namespace detail
{

class portion_mapply_op;
class materialized_mapply_tall_store;

/*
 * This class represents the result matrix of mapply operations.
 * It partially materializes a portion of the matrix when the portion
 * is needed. The underlying matrices that mapply runs on can be both
 * stored in memory and on disks. Therefore, this matrix should also
 * expose the EM object interface.
 */
class mapply_matrix_store: public virtual_matrix_store, public EM_object
{
	// This identifies the data in a matrix.
	// So when a matrix is transposed, it should share the same data id.
	const data_id_t::ptr data_id;

	/*
	 * This indicates whether the input matrices are accessed in parallel
	 * when the matrix is materialized.
	 */
	bool par_access;

	int num_nodes;

	matrix_layout_t layout;
	const std::vector<matrix_store::const_ptr> in_mats;
	portion_mapply_op::const_ptr op;

	std::shared_ptr<materialized_mapply_tall_store> res;
public:
	typedef std::shared_ptr<const mapply_matrix_store> const_ptr;

	// All input matrices should have the same shape, i.e., either all are
	// tall matrices or wide matrices. The output matrix should also have
	// the same shape as the input matrices.
	mapply_matrix_store(
			const std::vector<matrix_store::const_ptr> &in_mats,
			portion_mapply_op::const_ptr op,
			matrix_layout_t layout,
			data_id_t::ptr data_id = data_id_t::create(mat_counter++));

	virtual void inc_dag_ref(size_t data_id);
	virtual void reset_dag_ref();
	virtual size_t get_dag_ref() const {
		return data_id->get_ref();
	}

	bool is_wide() const {
		// We handle matrices with different shapes differently
		// We rely on this method to tell if a matrix is wide or tall.
		// This method doesn't work well on a square matrix. Instead of
		// judging the shape of the result matrix, we can judge
		// one of the input matrices.
		return in_mats.front()->is_wide();
	}

	const std::vector<matrix_store::const_ptr> get_input_mats() const {
		return in_mats;
	}

	virtual bool has_materialized() const;
	virtual size_t get_data_id() const {
		return data_id->get_id();
	}
	data_id_t::ptr get_id_ptr() const {
		return data_id;
	}

	virtual void set_prefetches(size_t num, std::pair<size_t, size_t> range);
	virtual void set_cache_portion(bool cache_portion);

	void set_par_access(bool par_access) {
		this->par_access = par_access;
	}

	virtual void set_materialize_level(materialize_level level,
			detail::matrix_store::ptr materialize_buf);

	virtual void materialize_self() const;

	virtual matrix_store::const_ptr materialize(bool in_mem,
		int num_nodes) const;

	virtual matrix_store::const_ptr get_rows(const std::vector<off_t> &idxs) const;

	using virtual_matrix_store::get_portion;
	virtual std::shared_ptr<const local_matrix_store> get_portion(
			size_t start_row, size_t start_col, size_t num_rows,
			size_t num_cols) const;
	virtual std::shared_ptr<const local_matrix_store> get_portion(
			size_t id) const;
	virtual int get_portion_node_id(size_t id) const;
	using virtual_matrix_store::get_portion_async;
	virtual async_cres_t get_portion_async(
			size_t start_row, size_t start_col, size_t num_rows,
			size_t num_cols, std::shared_ptr<portion_compute> compute) const;
	virtual std::pair<size_t, size_t> get_portion_size() const;
	virtual int get_num_nodes() const {
		return num_nodes;
	}

	virtual matrix_store::const_ptr transpose() const;

	virtual matrix_layout_t store_layout() const {
		return layout;
	}

	virtual std::vector<safs::io_interface::ptr> create_ios() const;

	virtual std::string get_name() const;
	virtual std::unordered_map<size_t, size_t> get_underlying_mats() const;
};

}

}

#endif
