#ifndef __FM_CUM_MATRIX_H__
#define __FM_CUM_MATRIX_H__

/*
 * Copyright 2017
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
#include "local_vec_store.h"
#include "materialize.h"
#include "bulk_operate_ext.h"

namespace fm
{

namespace detail
{

class cum_long_dim_op: public portion_mapply_op
{
	struct cum_res_t {
		local_vec_store::ptr data;
		// This indidates whether this accumulated result is global or local.
		bool is_global;
		// This indicates whether this data has been used. Each accumulated
		// result only needs to be used once.
		bool is_used;
		cum_res_t(local_vec_store::ptr data, bool global) {
			this->data = data;
			this->is_global = global;
			this->is_used = false;
		}
	};
	pthread_mutex_t lock;
	pthread_cond_t cond;
	// The key is the current cum location
	// (start_row * num_cols + start_col * num_rows). We assume that start_row
	// or start_col is 0.
	std::map<off_t, cum_res_t> curr_cums;

	matrix_margin margin;
	agg_operate::const_ptr op;
	// num_rows * num_cols of a portion.
	size_t portion_size;
public:
	cum_long_dim_op(matrix_margin margin, agg_operate::const_ptr op,
			size_t portion_size, size_t num_rows, size_t num_cols);

	virtual portion_mapply_op::const_ptr transpose() const {
		matrix_margin new_margin = margin
			== matrix_margin::MAR_ROW ? matrix_margin::MAR_COL : matrix_margin::MAR_ROW;
		return portion_mapply_op::const_ptr(new cum_long_dim_op(new_margin, op,
					portion_size, get_out_num_cols(), get_out_num_rows()));
	}

	virtual void run(
			const std::vector<std::shared_ptr<const local_matrix_store> > &ins,
			local_matrix_store &out) const;

	virtual std::string to_string(
			const std::vector<matrix_store::const_ptr> &mats) const {
		return std::string("cum(") + mats[0]->get_name() + ")";
	}

	/*
	 * We prefer a local matrix can't be resized. In this way, it'll be much
	 * easier to manage the accumulated results.
	 */
	virtual bool is_resizable(size_t local_start_row, size_t local_start_col,
			size_t local_num_rows, size_t local_num_cols) const {
		return false;
	}

	// A portion operation may fail. If so, the portion operation should
	// notify the main thread that the operation failed by override this method.
	virtual bool is_success() const {
		return true;
	}
};

}

}

#endif
