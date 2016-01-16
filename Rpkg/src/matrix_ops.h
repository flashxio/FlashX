#ifndef __MATRIX_OPS_H__
#define __MATRIX_OPS_H__

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
#include <vector>

#include <Rcpp.h>

#include "bulk_operate.h"
#include "bulk_operate_ext.h"

namespace fmr
{

/*
 * Register a binary UDF.
 * A user has to provide UDFs for all different types.
 */
void register_udf(const std::vector<fm::bulk_operate::const_ptr> &ops,
		const std::string &name);
/*
 * Register a unary UDF.
 * A user has to provide UDFs for all different types.
 */
void register_udf(const std::vector<fm::bulk_uoperate::const_ptr> &ops,
		const std::string &name);
/*
 * This registers more UDFs provided by FlashMatrix, which aren't basic UDFs.
 */
void init_udf_ext();

/*
 * Get a binary operator.
 * Even though a binary operator takes two inputs, FlashR always first casts
 * the type of one input to match the other. So we assume that all binary
 * operators take inputs of the same type.
 */
fm::bulk_operate::const_ptr get_op(SEXP pfun, fm::prim_type type);
/* Get a unary operator. */
fm::bulk_uoperate::const_ptr get_uop(SEXP pfun, fm::prim_type type);

typedef int op_id_t;

/* Get the binary operator Id given a name. */
op_id_t get_op_id(const std::string &name);
/* Get the unary operator Id given a name. */
op_id_t get_uop_id(const std::string &name);

/* This construct an aggregation operator from binary operators. */
fm::agg_operate::const_ptr get_agg_op(SEXP pfun, const fm::scalar_type &mat_type);

}

#endif
