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

#include <unordered_map>
#include <boost/filesystem.hpp>
#include <Rcpp.h>

#include "log.h"
#include "safs_file.h"

#include "FGlib.h"
#include "sparse_matrix.h"
#include "dense_matrix.h"
#include "bulk_operate.h"
#include "generic_type.h"

#include "rutils.h"

using namespace fm;

template<class ObjectType>
class object_ref
{
	typename ObjectType::ptr o;
public:
	object_ref(typename ObjectType::ptr o) {
		this->o = o;
	}

	typename ObjectType::ptr get_object() {
		return o;
	}
};

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

static SEXP create_FMR_matrix(sparse_matrix::ptr m, const std::string &name)
{
	Rcpp::List ret;
	ret["name"] = Rcpp::String(name);
	ret["type"] = Rcpp::String("sparse");

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

static SEXP create_FMR_matrix(dense_matrix::ptr m, const std::string &name)
{
	Rcpp::List ret;
	ret["name"] = Rcpp::String(name);
	ret["type"] = Rcpp::String("dense");

	object_ref<dense_matrix> *ref = new object_ref<dense_matrix>(m);
	SEXP pointer = R_MakeExternalPtr(ref, R_NilValue, R_NilValue);
	R_RegisterCFinalizerEx(pointer, fm_clean_SpM, TRUE);
	ret["pointer"] = pointer;

	Rcpp::NumericVector nrow(1);
	nrow[0] = m->get_num_rows();
	ret["nrow"] = nrow;

	Rcpp::NumericVector ncol(1);
	ncol[0] = m->get_num_cols();
	ret["ncol"] = ncol;

	return ret;
}

static SEXP create_FMR_vector(dense_matrix::ptr m, const std::string &name)
{
	Rcpp::List ret;
	ret["name"] = Rcpp::String(name);
	ret["type"] = Rcpp::String("vector");

	object_ref<dense_matrix> *ref = new object_ref<dense_matrix>(m);
	SEXP pointer = R_MakeExternalPtr(ref, R_NilValue, R_NilValue);
	R_RegisterCFinalizerEx(pointer, fm_clean_DM, TRUE);
	ret["pointer"] = pointer;

	Rcpp::NumericVector len(1);
	// TODO I assume the vector is stored as a nx1 matrix.
	len[0] = m->get_num_rows();
	ret["len"] = len;
	return ret;
}

/*
 * Test if this is a sparse matrix.
 */
static bool is_sparse(const Rcpp::List &matrix)
{
	return matrix["type"] == "sparse";
}

static bool is_vector(const Rcpp::List &matrix)
{
	return matrix["type"] == "vector";
}

template<class MatrixType>
static typename MatrixType::ptr get_matrix(const Rcpp::List &matrix)
{
	object_ref<MatrixType> *ref
		= (object_ref<MatrixType> *) R_ExternalPtrAddr(matrix["pointer"]);
	return ref->get_object();
}

fg::FG_graph::ptr R_FG_get_graph(SEXP pgraph);

template<class EntryType>
class set_const_operate: public set_operate
{
	EntryType v;
public:
	set_const_operate(EntryType v) {
		this->v = v;
	}

	virtual void set(void *arr, size_t num_eles, off_t row_idx,
			off_t col_idx) const {
		EntryType *ele_p = (EntryType *) arr;
		for (size_t i = 0; i < num_eles; i++)
			ele_p[i] = v;
	}

	virtual size_t entry_size() const {
		return sizeof(EntryType);
	}
};

template<class EntryType>
dense_matrix::ptr create_dense_matrix(size_t nrow, size_t ncol,
		matrix_layout_t layout, EntryType initv)
{
	dense_matrix::ptr m = dense_matrix::create(nrow, ncol, sizeof(EntryType),
			// TODO let's just use in-memory dense matrix first.
			layout, true);
	m->set_data(set_const_operate<EntryType>(initv));
	return m;
}

RcppExport SEXP R_FM_create_vector(SEXP plen, SEXP pinitv)
{
	size_t len = REAL(plen)[0];

	dense_matrix::ptr m;
	if (R_is_real(pinitv))
		m = create_dense_matrix<double>(len, 1, matrix_layout_t::L_COL,
				REAL(pinitv)[0]);
	else if (R_is_integer(pinitv))
		m = create_dense_matrix<int>(len, 1, matrix_layout_t::L_COL,
				INTEGER(pinitv)[0]);
	else {
		fprintf(stderr, "The initial value has unsupported type\n");
		return Rcpp::List();
	}

	return create_FMR_vector(m, "");
}

RcppExport SEXP R_FM_get_matrix_fg(SEXP pgraph)
{
	Rcpp::List graph = Rcpp::List(pgraph);
	Rcpp::LogicalVector res(1);
	fg::FG_graph::ptr fg = R_FG_get_graph(pgraph);
	sparse_matrix::ptr m = sparse_matrix::create(fg);
	std::string name = graph["name"];
	return create_FMR_matrix(m, name);
}

/*
 * R has only two data types in matrix multiplication: integer and numeric.
 * So we only need to predefine a small number of basic operations with
 * different types.
 */

static basic_ops_impl<int, int, int> R_basic_ops_II;
static basic_ops_impl<double, int, double> R_basic_ops_DI;
static basic_ops_impl<int, double, double> R_basic_ops_ID;
static basic_ops_impl<double, double, double> R_basic_ops_DD;

static basic_ops &get_inner_prod_left_ops(const dense_matrix &left,
		const dense_matrix &right)
{
	if (left.get_entry_size() == sizeof(int)
			&& right.get_entry_size() == sizeof(int))
		return R_basic_ops_II;
	else if (left.get_entry_size() == sizeof(double)
			&& right.get_entry_size() == sizeof(int))
		return R_basic_ops_DI;
	else if (left.get_entry_size() == sizeof(int)
			&& right.get_entry_size() == sizeof(double))
		return R_basic_ops_ID;
	else if (left.get_entry_size() == sizeof(double)
			&& right.get_entry_size() == sizeof(double))
		return R_basic_ops_DD;
	else {
		fprintf(stderr, "the matrix has a wrong type\n");
		abort();
	}
}

static basic_ops &get_inner_prod_right_ops(const bulk_operate &left_ops)
{
	if (left_ops.output_entry_size() == 4)
		return R_basic_ops_II;
	else if (left_ops.output_entry_size() == 8)
		return R_basic_ops_DD;
	else {
		fprintf(stderr,
				"the left operator of inner product has a wrong output type\n");
		abort();
	}
}

static SEXP SpMV(sparse_matrix::ptr matrix, dense_matrix::ptr right_mat)
{
	if (right_mat->is_type<double>()) {
		mem_vector<double>::ptr in_vec = mem_vector<double>::create(
				mem_dense_matrix::cast(right_mat));
		mem_vector<double>::ptr out_vec = matrix->multiply<double>(in_vec);
		return create_FMR_vector(out_vec->get_data(), "");
	}
	else if (right_mat->is_type<int>()) {
		mem_vector<int>::ptr in_vec = mem_vector<int>::create(
				mem_dense_matrix::cast(right_mat));
		mem_vector<int>::ptr out_vec = matrix->multiply<int>(in_vec);
		return create_FMR_vector(out_vec->get_data(), "");
	}
	else {
		fprintf(stderr, "the input vector has an unsupported type in SpMV\n");
		return R_NilValue;
	}
}

static SEXP SpMM(sparse_matrix::ptr matrix, dense_matrix::ptr right_mat)
{
	fprintf(stderr,
			"Sparse matrix dense matrix multiplication isn't supported\n");
	return R_NilValue;
}

static bool is_vector(const dense_matrix &mat)
{
	// If the matrix has one row or one column, we consider it as a vector.
	return mat.get_num_rows() == 1 || mat.get_num_cols() == 1;
}

RcppExport SEXP R_FM_multiply_sparse(SEXP pmatrix, SEXP pmat)
{
	dense_matrix::ptr right_mat = get_matrix<dense_matrix>(pmat);
	if (!right_mat->is_in_mem()) {
		fprintf(stderr, "we now only supports in-mem vector for SpMV\n");
		return R_NilValue;
	}
	sparse_matrix::ptr matrix = get_matrix<sparse_matrix>(pmatrix);
	if (is_vector(*right_mat))
		return SpMV(matrix, right_mat);
	else
		return SpMM(matrix, right_mat);
}

RcppExport SEXP R_FM_multiply_dense(SEXP pmatrix, SEXP pmat)
{
	dense_matrix::ptr right_mat = get_matrix<dense_matrix>(pmat);
	dense_matrix::ptr matrix = get_matrix<dense_matrix>(pmatrix);
	const bulk_operate &left_op = get_inner_prod_left_ops(*matrix,
			*right_mat).get_multiply();
	const bulk_operate &right_op = get_inner_prod_right_ops(left_op).get_add();
	dense_matrix::ptr prod_vec = matrix->inner_prod(*right_mat, left_op,
			right_op);
	return create_FMR_vector(prod_vec, "");
}

template<class T>
T matrix_sum(const dense_matrix &mat)
{
	scalar_type_impl<T> res;
	basic_ops_impl<T, T, T> ops;
	mat.aggregate(ops.get_add(), res);
	return res.get();
}

RcppExport SEXP R_FM_matrix_sum(SEXP pmat)
{
	dense_matrix::ptr mat = get_matrix<dense_matrix>(pmat);
	if (mat->is_type<double>()) {
		Rcpp::NumericVector ret(1);
		ret[0] = matrix_sum<double>(*mat);
		return ret;
	}
	else if (mat->is_type<int>()) {
		Rcpp::NumericVector ret(1);
		ret[0] = matrix_sum<int>(*mat);
		return ret;
	}
	else {
		fprintf(stderr, "The matrix has an unsupported type for sum\n");
		return R_NilValue;
	}
}

RcppExport SEXP R_FM_conv_matrix(SEXP pmat, SEXP pnrow, SEXP pncol, SEXP pbyrow)
{
	Rcpp::List matrix_obj(pmat);
	if (is_sparse(matrix_obj)) {
		fprintf(stderr, "We can't change the dimension of a sparse matrix\n");
		return R_NilValue;
	}

	size_t nrow = REAL(pnrow)[0];
	size_t ncol = REAL(pncol)[0];
	bool byrow = LOGICAL(pbyrow)[0];
	dense_matrix::ptr mat = get_matrix<dense_matrix>(pmat);
	return create_FMR_matrix(mat->conv2(nrow, ncol, byrow), "");
}

template<class T, class RContainerType>
void copy_FM2Rvector(const mem_vector<T> &vec, RContainerType &r_arr)
{
	printf("copy a FM vector to R vector\n");
	size_t length = vec.get_length();
	for (size_t i = 0; i < length; i++) {
		r_arr[i] = vec.get(i);
	}
}

template<class T, class RContainerType>
void copy_FM2Rmatrix(const type_mem_dense_matrix<T> &mat,
		RContainerType &r_mat)
{
	// TODO this is going to be slow. But I don't care about performance
	// for now.
	size_t nrow = mat.get_num_rows();
	size_t ncol = mat.get_num_cols();
	for (size_t i = 0; i < nrow; i++) {
		for (size_t j = 0; j < ncol; j++) {
			r_mat(i, j) = mat.get(i, j);
		}
	}
}

RcppExport SEXP R_FM_conv_FM2R(SEXP pobj)
{
	Rcpp::List matrix_obj(pobj);
	if (is_sparse(matrix_obj)) {
		fprintf(stderr, "We can't convert a sparse matrix to an R object\n");
		return R_NilValue;
	}

	dense_matrix::ptr mat = get_matrix<dense_matrix>(pobj);
	if (!mat->is_in_mem()) {
		fprintf(stderr, "We only support in-memory matrix right now\n");
		return R_NilValue;
	}

	mem_dense_matrix::ptr mem_mat = mem_dense_matrix::cast(mat);
	if (mem_mat->is_type<double>()) {
		if (is_vector(pobj)) {
			mem_vector<double>::ptr mem_vec = mem_vector<double>::create(mem_mat);
			Rcpp::NumericVector ret(mem_vec->get_length());
			copy_FM2Rvector<double, Rcpp::NumericVector>(*mem_vec, ret);
			return ret;
		}
		else {
			Rcpp::NumericMatrix ret(mem_mat->get_num_rows(),
					mem_mat->get_num_cols());
			copy_FM2Rmatrix<double, Rcpp::NumericMatrix>(
					*type_mem_dense_matrix<double>::create(mem_mat), ret);
			return ret;
		}
	}
	else if (mem_mat->is_type<int>()) {
		if (is_vector(pobj)) {
			mem_vector<int>::ptr mem_vec = mem_vector<int>::create(mem_mat);
			Rcpp::IntegerVector ret(mem_vec->get_length());
			copy_FM2Rvector<int, Rcpp::IntegerVector>(*mem_vec, ret);
			return ret;
		}
		else {
			Rcpp::IntegerMatrix ret(mem_mat->get_num_rows(),
					mem_mat->get_num_cols());
			copy_FM2Rmatrix<int, Rcpp::IntegerMatrix>(
					*type_mem_dense_matrix<int>::create(mem_mat), ret);
			return ret;
		}
	}
	else {
		fprintf(stderr, "the dense matrix doesn't have a right type\n");
		return R_NilValue;
	}
}

RcppExport SEXP R_FM_conv_RVec2FM(SEXP pobj)
{
	if (R_is_real(pobj)) {
		Rcpp::NumericVector vec(pobj);
		mem_vector<double>::ptr fm_vec = mem_vector<double>::create(vec.size());
		for (size_t i = 0; i < fm_vec->get_length(); i++)
			fm_vec->set(i, vec[i]);
		return create_FMR_vector(fm_vec->get_data(), "");
	}
	else if (R_is_integer(pobj)) {
		Rcpp::IntegerVector vec(pobj);
		mem_vector<int>::ptr fm_vec = mem_vector<int>::create(vec.size());
		for (size_t i = 0; i < fm_vec->get_length(); i++)
			fm_vec->set(i, vec[i]);
		return create_FMR_vector(fm_vec->get_data(), "");
	}
	else {
		fprintf(stderr, "The R vector has an unsupported type\n");
		return R_NilValue;
	}
}

RcppExport SEXP R_FM_conv_RMat2FM(SEXP pobj, SEXP pbyrow)
{
	bool byrow = LOGICAL(pbyrow)[0];
	if (R_is_real(pobj)) {
		Rcpp::NumericMatrix mat(pobj);
		size_t nrow = mat.nrow();
		size_t ncol = mat.ncol();
		type_mem_dense_matrix<double>::ptr fm_mat
			= type_mem_dense_matrix<double>::create(nrow, ncol,
					byrow ? matrix_layout_t::L_ROW : matrix_layout_t::L_COL);
		for (size_t i = 0; i < nrow; i++)
			for (size_t j = 0; j < ncol; j++)
				fm_mat->set(i, j, mat(i, j));
		return create_FMR_matrix(fm_mat->get_matrix(), "");
	}
	else if (R_is_integer(pobj)) {
		Rcpp::IntegerMatrix mat(pobj);
		size_t nrow = mat.nrow();
		size_t ncol = mat.ncol();
		type_mem_dense_matrix<int>::ptr fm_mat
			= type_mem_dense_matrix<int>::create(nrow, ncol,
					byrow ? matrix_layout_t::L_ROW : matrix_layout_t::L_COL);
		for (size_t i = 0; i < nrow; i++)
			for (size_t j = 0; j < ncol; j++)
				fm_mat->set(i, j, mat(i, j));
		return create_FMR_matrix(fm_mat->get_matrix(), "");
	}
	else {
		fprintf(stderr, "The R vector has an unsupported type\n");
		return R_NilValue;
	}
}
