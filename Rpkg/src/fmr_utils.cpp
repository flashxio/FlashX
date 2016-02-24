/*
 * Copyright 2014 Open Connectome Project (http://openconnecto.me)
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

#include <stdio.h>

#include "data_frame.h"
#include "mem_matrix_store.h"
#include "sparse_matrix.h"
#include "factor.h"

#include "fmr_utils.h"
#include "rutils.h"

using namespace fm;

/*
 * Clean up a sparse matrix.
 */
static void fm_clean_SpM(SEXP p)
{
	object_ref<sparse_matrix> *ref
		= (object_ref<sparse_matrix> *) R_ExternalPtrAddr(p);
	delete ref;
}

/*
 * Clean up a dense matrix
 */
static void fm_clean_DM(SEXP p)
{
	object_ref<dense_matrix> *ref
		= (object_ref<dense_matrix> *) R_ExternalPtrAddr(p);
	delete ref;
}

SEXP create_FMR_matrix(sparse_matrix::ptr m, const std::string &name)
{
	Rcpp::List ret;
	ret["name"] = Rcpp::String(name);
	ret["type"] = Rcpp::String("sparse");
	if (m->is_type<int>())
		ret["ele_type"] = Rcpp::String("integer");
	else if (m->is_type<double>())
		ret["ele_type"] = Rcpp::String("double");
	else if (m->is_type<bool>())
		ret["ele_type"] = Rcpp::String("logical");
	else
		ret["ele_type"] = Rcpp::String("unknown");

	object_ref<sparse_matrix> *ref = new object_ref<sparse_matrix>(m);
	SEXP pointer = R_MakeExternalPtr(ref, R_NilValue, R_NilValue);
	R_RegisterCFinalizerEx(pointer, fm_clean_SpM, TRUE);
	ret["pointer"] = pointer;

	Rcpp::LogicalVector sym(1);
	sym[0] = m->is_symmetric();
	ret["sym"] = sym;

	Rcpp::NumericVector nrow(1);
	nrow[0] = m->get_num_rows();
	ret["nrow"] = nrow;

	Rcpp::NumericVector ncol(1);
	ncol[0] = m->get_num_cols();
	ret["ncol"] = ncol;

	return ret;
}

SEXP create_FMR_matrix(dense_matrix::ptr m, const std::string &name)
{
	Rcpp::List ret;
	ret["name"] = Rcpp::String(name);
	ret["type"] = Rcpp::String("dense");
	if (m->is_type<int>())
		ret["ele_type"] = Rcpp::String("integer");
	else if (m->is_type<double>())
		ret["ele_type"] = Rcpp::String("double");
	else if (m->is_type<bool>()) {
		m = m->cast_ele_type(get_scalar_type<int>());
		ret["ele_type"] = Rcpp::String("logical");
	}
	else
		ret["ele_type"] = Rcpp::String("unknown");

	object_ref<dense_matrix> *ref = new object_ref<dense_matrix>(m);
	SEXP pointer = R_MakeExternalPtr(ref, R_NilValue, R_NilValue);
	R_RegisterCFinalizerEx(pointer, fm_clean_DM, TRUE);
	ret["pointer"] = pointer;

	Rcpp::NumericVector nrow(1);
	nrow[0] = m->get_num_rows();
	ret["nrow"] = nrow;

	Rcpp::NumericVector ncol(1);
	ncol[0] = m->get_num_cols();
	ret["ncol"] = ncol;

	return ret;
}

SEXP create_FMR_vector(detail::vec_store::const_ptr vec, const std::string &name)
{
	detail::matrix_store::const_ptr mat = vec->conv2mat(vec->get_length(),
			1, false);
	return create_FMR_vector(dense_matrix::create(
				detail::mem_matrix_store::cast(mat)), name);
}

SEXP create_FMR_vector(dense_matrix::ptr m, const std::string &name)
{
	Rcpp::List ret;
	ret["name"] = Rcpp::String(name);
	ret["type"] = Rcpp::String("vector");
	if (m->is_type<int>())
		ret["ele_type"] = Rcpp::String("integer");
	else if (m->is_type<double>())
		ret["ele_type"] = Rcpp::String("double");
	else if (m->is_type<bool>()) {
		m = m->cast_ele_type(get_scalar_type<int>());
		ret["ele_type"] = Rcpp::String("logical");
	}
	else
		ret["ele_type"] = Rcpp::String("unknown");

	object_ref<dense_matrix> *ref = new object_ref<dense_matrix>(m);
	SEXP pointer = R_MakeExternalPtr(ref, R_NilValue, R_NilValue);
	R_RegisterCFinalizerEx(pointer, fm_clean_DM, TRUE);
	ret["pointer"] = pointer;

	Rcpp::NumericVector len(1);
	if (m->get_num_cols() == 1)
		len[0] = m->get_num_rows();
	else
		len[0] = m->get_num_cols();
	ret["len"] = len;

	return ret;
}

SEXP create_FMR_factor_vector(dense_matrix::ptr m, int num_levels,
		const std::string &name)
{
	Rcpp::List ret = create_FMR_vector(m, name);
	Rcpp::NumericVector levels(1);
	levels[0] = num_levels;
	ret["levels"] = num_levels;
	return ret;
}

vector::ptr get_vector(const Rcpp::S4 &vec)
{
	if (!is_vector(vec)) {
		fprintf(stderr, "The S4 object isn't a vector\n");
		return vector::ptr();
	}
	object_ref<dense_matrix> *ref
		= (object_ref<dense_matrix> *) R_ExternalPtrAddr(vec.slot("pointer"));
	dense_matrix::ptr mat = ref->get_object();
	// This should be a column matrix.
	assert(mat->store_layout == matrix_layout_t::L_COL
			&& mat->get_num_cols() == 1);
	detail::vec_store::const_ptr store = mat->get_data().get_col_vec(0);
	if (store == NULL) {
		fprintf(stderr, "can't convert a matrix to a vector");
		return vector::ptr();
	}
	return vector::create(store);
}

factor_col_vector::ptr get_factor_vector(const Rcpp::S4 &vec)
{
	if (!is_factor_vector(vec)) {
		fprintf(stderr, "The S4 object isn't a factor vector\n");
		return factor_col_vector::ptr();
	}
	object_ref<dense_matrix> *ref
		= (object_ref<dense_matrix> *) R_ExternalPtrAddr(vec.slot("pointer"));
	dense_matrix::ptr mat = ref->get_object();
	// This should be a column matrix.
	assert(mat->get_num_cols() == 1);
	size_t num_levels = vec.slot("num.levels");
	return factor_col_vector::create(factor(num_levels), mat);
}

SEXP create_FMR_data_frame(data_frame::ptr df, const std::string &name)
{
	Rcpp::List ret;
	for (size_t i = 0; i < df->get_num_vecs(); i++) {
		std::string vec_name = df->get_vec_name(i);
		ret[vec_name] = create_FMR_vector(df->get_vec(i), vec_name);
	}
	return ret;
}

SEXP create_FMR_sinkV(dense_matrix::ptr m, size_t len, const std::string &name)
{
	Rcpp::List ret;
	ret["name"] = Rcpp::String(name);
	ret["type"] = Rcpp::String("vector");
	if (m->is_type<int>())
		ret["ele_type"] = Rcpp::String("integer");
	else if (m->is_type<double>())
		ret["ele_type"] = Rcpp::String("double");
	else if (m->is_type<bool>()) {
		ret["ele_type"] = Rcpp::String("logical");
	}
	else
		ret["ele_type"] = Rcpp::String("unknown");

	object_ref<dense_matrix> *ref = new object_ref<dense_matrix>(m);
	SEXP pointer = R_MakeExternalPtr(ref, R_NilValue, R_NilValue);
	R_RegisterCFinalizerEx(pointer, fm_clean_DM, TRUE);
	ret["pointer"] = pointer;

	Rcpp::NumericVector nrow(1);
	nrow[0] = len;
	ret["nrow"] = nrow;

	Rcpp::NumericVector ncol(1);
	ncol[0] = 1;
	ret["ncol"] = ncol;

	return ret;
}
