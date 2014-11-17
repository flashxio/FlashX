/**
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

#include <unordered_map>
#include <boost/filesystem.hpp>
#include <Rcpp.h>

#include "log.h"
#include "safs_file.h"

#include "FGlib.h"
#include "sparse_matrix.h"

class matrix_ref
{
	sparse_matrix::ptr m;
public:
	matrix_ref(sparse_matrix::ptr m) {
		this->m = m;
	}

	sparse_matrix::ptr get_matrix() {
		return m;
	}
};

static void fm_clean_matrix(SEXP p)
{
	matrix_ref *ref = (matrix_ref *) R_ExternalPtrAddr(p);
	delete ref;
}

static SEXP create_FMR_object(sparse_matrix::ptr m, const std::string &name)
{
	Rcpp::List ret;
	ret["name"] = Rcpp::String(name);

	matrix_ref *ref = new matrix_ref(m);
	SEXP pointer = R_MakeExternalPtr(ref, R_NilValue, R_NilValue);
	R_RegisterCFinalizerEx(pointer, fm_clean_matrix, TRUE);
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

static sparse_matrix::ptr get_matrix(SEXP pmatrix)
{
	Rcpp::List matrix(pmatrix);
	matrix_ref *ref = (matrix_ref *) R_ExternalPtrAddr(matrix["pointer"]);
	return ref->get_matrix();
}

FG_graph::ptr R_FG_get_graph(SEXP pgraph);

RcppExport SEXP R_FM_get_matrix_fg(SEXP pgraph)
{
	Rcpp::List graph = Rcpp::List(pgraph);
	Rcpp::LogicalVector res(1);
	FG_graph::ptr fg = R_FG_get_graph(pgraph);
	sparse_matrix::ptr m = sparse_matrix::create(fg);
	std::string name = graph["name"];
	return create_FMR_object(m, name);
}

RcppExport SEXP R_FM_multiply_v(SEXP pmatrix, SEXP pvec)
{
	Rcpp::NumericVector vec(pvec);
	size_t length = vec.size();
	FG_vector<double>::ptr in_vec = FG_vector<double>::create(length);
	for (size_t i = 0; i < length; i++) {
		in_vec->get_data()[i] = vec[i];
	}
	sparse_matrix::ptr matrix = get_matrix(pmatrix);
	assert(matrix->get_num_rows() == length);
	assert(matrix->get_num_cols() == length);

	FG_vector<double>::ptr out_vec = matrix->multiply<double>(in_vec);
	Rcpp::NumericVector ret(out_vec->get_data(), out_vec->get_data() + length);
	return ret;
}
