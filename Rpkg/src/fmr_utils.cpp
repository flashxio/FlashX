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

static inline Rcpp::String trans_RType2Str(R_type type)
{
	if (type == R_type::R_INT)
		return Rcpp::String("integer");
	else if (type == R_type::R_REAL)
		return Rcpp::String("double");
	else if (type == R_type::R_LOGICAL)
		return Rcpp::String("logical");
	else
		return Rcpp::String("unknown");
}

SEXP create_FMR_matrix(sparse_matrix::ptr m, R_type type, const std::string &name)
{
	if (m == NULL) {
		fprintf(stderr, "can't create an empty matrix\n");
		return R_NilValue;
	}

	Rcpp::List ret;
	ret["name"] = Rcpp::String(name);
	ret["type"] = Rcpp::String("sparse");
	ret["ele_type"] = trans_RType2Str(type);

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

SEXP create_FMR_matrix(dense_matrix::ptr m, R_type type, const std::string &name)
{
	if (m == NULL) {
		fprintf(stderr, "can't create an empty matrix\n");
		return R_NilValue;
	}

	Rcpp::List ret;
	ret["name"] = Rcpp::String(name);
	ret["type"] = Rcpp::String("dense");
	ret["ele_type"] = trans_RType2Str(type);

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

SEXP create_FMR_vector(detail::vec_store::const_ptr vec, R_type type,
		const std::string &name)
{
	detail::matrix_store::const_ptr mat = vec->conv2mat(vec->get_length(),
			1, false);
	if (mat == NULL)
		return R_NilValue;
	return create_FMR_vector(dense_matrix::create(mat), type, name);
}

SEXP create_FMR_vector(dense_matrix::ptr m, R_type type, const std::string &name)
{
	if (m == NULL) {
		fprintf(stderr, "can't create a vector from an empty matrix\n");
		return R_NilValue;
	}

	if (m->get_num_cols() > 1)
		m = m->transpose();
	if (m->get_num_cols() > 1) {
		fprintf(stderr,
				"can't create a vector with a matrix with more than one col\n");
		return R_NilValue;
	}

	Rcpp::List ret;
	ret["name"] = Rcpp::String(name);
	ret["type"] = Rcpp::String("vector");
	ret["ele_type"] = trans_RType2Str(type);

	object_ref<dense_matrix> *ref = new object_ref<dense_matrix>(m);
	SEXP pointer = R_MakeExternalPtr(ref, R_NilValue, R_NilValue);
	R_RegisterCFinalizerEx(pointer, fm_clean_DM, TRUE);
	ret["pointer"] = pointer;

	Rcpp::NumericVector len(1);
	len[0] = m->get_num_rows();
	ret["len"] = len;

	return ret;
}

SEXP create_FMR_vector(fm::factor_col_vector::ptr v, R_type type,
		const std::string &name)
{
	if (v == NULL) {
		fprintf(stderr, "can't create a factor vector\n");
		return R_NilValue;
	}

	Rcpp::List ret;
	ret["name"] = Rcpp::String(name);
	ret["type"] = Rcpp::String("vector");
	ret["ele_type"] = Rcpp::String("integer");

	// The element type of a factor vector is factor_value_t in FlashMatrix.
	// FlashR can only handle int or double, so let's cast it into int.
	dense_matrix::ptr mat = v->cast_ele_type(get_scalar_type<int>(), true);
	object_ref<dense_matrix> *ref = new object_ref<dense_matrix>(mat);
	SEXP pointer = R_MakeExternalPtr(ref, R_NilValue, R_NilValue);
	R_RegisterCFinalizerEx(pointer, fm_clean_DM, TRUE);
	ret["pointer"] = pointer;

	Rcpp::NumericVector len(1);
	len[0] = v->get_length();
	ret["len"] = len;

	Rcpp::IntegerVector nlevels(1);
	nlevels[0] = v->get_num_levels();
	ret["num.levels"] = nlevels;

	auto vals = v->get_uniq_vals();
	if (vals == NULL)
		ret["vals"] = R_NilValue;
	else
		ret["vals"] = create_FMR_vector(vals, type, "");
	auto cnts = v->get_counts();
	if (cnts == NULL)
		ret["cnts"] = R_NilValue;
	else {
		// TODO we might want to use double floating-point to count
		// more elements.
		auto cnt_vec = vector::create(cnts);
		auto cnt_mat = cnt_vec->conv2mat(cnt_vec->get_length(), 1, false);
		cnt_mat = cnt_mat->cast_ele_type(get_scalar_type<int>());
		ret["cnts"] = create_FMR_vector(cnt_mat, R_type::R_INT, "");
	}

	return ret;
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

col_vec::ptr get_vector(const Rcpp::S4 &vec)
{
	dense_matrix::ptr mat = get_matrix<dense_matrix>(vec);
	if (mat->get_num_rows() > 1 && mat->get_num_cols() > 1) {
		fprintf(stderr, "The input object is a matrix");
		return col_vec::ptr();
	}
	return col_vec::create(mat);
}

SEXP create_FMR_data_frame(data_frame::ptr df,
		const std::vector<R_type> &col_types, const std::string &name)
{
	if (col_types.size() != df->get_num_vecs()) {
		fprintf(stderr,
				"# col types doesn't match with # cols in the data frame\n");
		return R_NilValue;
	}

	Rcpp::List ret;
	for (size_t i = 0; i < df->get_num_vecs(); i++) {
		std::string vec_name = df->get_vec_name(i);
		ret[vec_name] = create_FMR_vector(df->get_vec(i), col_types[i],
				vec_name);
	}
	return ret;
}
