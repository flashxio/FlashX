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

#ifndef __FM_AGG_MATRIX_STORE_H__
#define __FM_AGG_MATRIX_STORE_H__

#include "io_interface.h"

#include "sink_matrix.h"
#include "bulk_operate.h"
#include "bulk_operate_ext.h"
#include "materialize.h"
#include "EM_object.h"

namespace fm
{

namespace detail
{

class matrix_long_agg_op;

/*
 * This matrix store is to enable lazy evaluation on the aggregation on
 * the long dimension of a dense matrix or the entire dense matrix. 
 */
class agg_matrix_store: public sink_store
{
	std::shared_ptr<const matrix_long_agg_op> agg_op;
	matrix_store::const_ptr data;
	// We often need to check the underlying matrices of a sink matrix
	// when materializing a virtual matrix.
	// This caches the underlying matrices that are used to compute
	// this sink matrix. We rarely materialize any non-sink matrices,
	// so the underlying matrices usually remain the same.
	std::unordered_map<size_t, size_t> underlying;

	matrix_store::const_ptr get_agg_res() const;
	agg_matrix_store(matrix_store::const_ptr data, matrix_margin margin,
			agg_operate::const_ptr op);
	agg_matrix_store(matrix_store::const_ptr data,
			std::shared_ptr<const matrix_long_agg_op> portion_op);
public:
	typedef std::shared_ptr<const agg_matrix_store> const_ptr;

	static ptr create(matrix_store::const_ptr data, matrix_margin margin,
			agg_operate::const_ptr op) {
		ptr ret(new agg_matrix_store(data, margin, op));
		sink_store::register_sink_matrices(ret);
		return ret;
	}

	virtual void inc_dag_ref(size_t id) {
		const_cast<matrix_store &>(*data).inc_dag_ref(get_data_id());
		// We don't need to increase the ref count of a sink matrix
		// because we never get a portion from a sink matrix.
	}
	virtual void reset_dag_ref() {
		const_cast<matrix_store &>(*data).reset_dag_ref();
	}

	virtual bool has_materialized() const;

	virtual void materialize_self() const;

	virtual std::vector<virtual_matrix_store::const_ptr> get_compute_matrices() const;
	virtual matrix_store::const_ptr get_result() const {
		if (has_materialized())
			return get_agg_res();
		else
			return matrix_store::const_ptr();
	}
	virtual size_t get_data_id() const;

	virtual matrix_store::const_ptr materialize(bool in_mem,
		int num_nodes) const;

	virtual matrix_store::const_ptr transpose() const;

	virtual matrix_layout_t store_layout() const {
		if (is_wide())
			return matrix_layout_t::L_ROW;
		else
			return matrix_layout_t::L_COL;
	}

	virtual std::string get_name() const;

	virtual std::unordered_map<size_t, size_t> get_underlying_mats() const {
		if (has_materialized())
			return std::unordered_map<size_t, size_t>();
		else if (!underlying.empty())
			return underlying;
		else
			return data->get_underlying_mats();
	}
};

}

}

#endif
