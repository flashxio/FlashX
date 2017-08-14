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

#ifndef __FM_GROUPBY_MATRIX_H__
#define __FM_GROUPBY_MATRIX_H__

#include "sink_matrix.h"
#include "EM_object.h"
#include "factor.h"

namespace fm
{

class factor_col_vector;
class agg_operate;

namespace detail
{

class portion_mapply_op;

class groupby_matrix_store: public sink_store
{
	std::shared_ptr<portion_mapply_op> portion_op;
	std::shared_ptr<const agg_operate> agg_op;
	matrix_store::const_ptr data;
	matrix_store::const_ptr label_store;
	matrix_margin margin;
	factor f;
	// We often need to check the underlying matrices of a sink matrix
	// when materializing a virtual matrix.
	// This caches the underlying matrices that are used to compute
	// this sink matrix. We rarely materialize any non-sink matrices,
	// so the underlying matrices usually remain the same.
	std::unordered_map<size_t, size_t> underlying;

	matrix_store::ptr get_agg_res() const;

	groupby_matrix_store(matrix_store::const_ptr data,
			matrix_store::const_ptr label_store, const factor &f,
			matrix_margin margin, agg_operate::const_ptr op);

	groupby_matrix_store(matrix_store::const_ptr data,
			std::shared_ptr<const factor_col_vector> labels,
			matrix_margin margin, agg_operate::const_ptr op);
public:
	typedef std::shared_ptr<const agg_matrix_store> const_ptr;

	static ptr create(matrix_store::const_ptr data,
			std::shared_ptr<const factor_col_vector> labels,
			matrix_margin margin, agg_operate::const_ptr op) {
		ptr ret(new groupby_matrix_store(data, labels, margin, op));
		sink_store::register_sink_matrices(ret);
		return ret;
	}
	virtual size_t get_data_id() const;

	virtual void inc_dag_ref(size_t id) {
		const_cast<matrix_store &>(*data).inc_dag_ref(get_data_id());
		const_cast<matrix_store &>(*label_store).inc_dag_ref(get_data_id());
		// We don't need to increase the ref count of a sink matrix
		// because we never get a portion from a sink matrix.
	}
	virtual void reset_dag_ref() {
		const_cast<matrix_store &>(*data).reset_dag_ref();
		const_cast<matrix_store &>(*label_store).reset_dag_ref();
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
	virtual matrix_store::const_ptr materialize(bool in_mem,
		int num_nodes) const;

	virtual matrix_store::const_ptr transpose() const;

	virtual matrix_layout_t store_layout() const {
		// TODO what is the right layout?
		return data->store_layout();
	}

	virtual std::unordered_map<size_t, size_t> get_underlying_mats() const {
		if (has_materialized())
			return std::unordered_map<size_t, size_t>();
		else if (!underlying.empty())
			return underlying;
		else
			return data->get_underlying_mats();
	}

	virtual std::string get_name() const {
		std::vector<matrix_store::const_ptr> mats(1);
		mats[0] = data;
		return portion_op->to_string(mats);
	}
};

}

}

#endif
