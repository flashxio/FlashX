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

class portion_mapply_op;

/*
 * This matrix store is to enable lazy evaluation on the aggregation on
 * the long dimension of a dense matrix or the entire dense matrix. 
 */
class agg_matrix_store: public sink_store
{
	std::shared_ptr<portion_mapply_op> portion_op;
	matrix_store::const_ptr data;

	matrix_store::const_ptr get_agg_res() const;
public:
	typedef std::shared_ptr<const agg_matrix_store> const_ptr;

	agg_matrix_store(matrix_store::const_ptr data, matrix_margin margin,
			agg_operate::const_ptr op);

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

	virtual std::string get_name() const {
		std::vector<matrix_store::const_ptr> mats(1);
		mats[0] = data;
		return portion_op->to_string(mats);
	}
	virtual std::unordered_map<size_t, size_t> get_underlying_mats() const {
		return data->get_underlying_mats();
	}
};

}

}

#endif
