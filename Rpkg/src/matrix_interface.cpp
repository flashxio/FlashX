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
#include "bulk_operate.h"
#include "generic_type.h"
#include "eigensolver/eigensolver.h"

#include "rutils.h"
#include "fm_utils.h"

using namespace fm;

fg::FG_graph::ptr R_FG_get_graph(SEXP pgraph);

int num_nodes = -1;

template<class EntryType>
dense_matrix::ptr create_dense_matrix(size_t nrow, size_t ncol,
		matrix_layout_t layout, EntryType initv)
{
	return dense_matrix::create_const(initv, nrow, ncol, layout);
}

RcppExport SEXP R_FM_create_vector(SEXP plen, SEXP pinitv)
{
	size_t len = REAL(plen)[0];

	vector::ptr vec;
	if (R_is_real(pinitv))
		vec = create_vector<double>(len, REAL(pinitv)[0]);
	else if (R_is_integer(pinitv))
		vec = create_vector<int>(len, INTEGER(pinitv)[0]);
	else {
		fprintf(stderr, "The initial value has unsupported type\n");
		return Rcpp::List();
	}

	return create_FMR_vector(vec->get_raw_store(), "");
}

template<class T>
class rand_set_operate: public type_set_vec_operate<T>
{
	const T min;
	const T max;

	T gen_rand() const {
		// We need to rescale and shift the random number accordingly.
		return unif_rand() * (max - min) + min;
	}
public:
	rand_set_operate(T _min, T _max): min(_min), max(_max) {
	}

	virtual void set(T *arr, size_t num_eles, off_t start_idx) const {
		for (size_t i = 0; i < num_eles; i++) {
			arr[i] = gen_rand();
		}
	}
};

RcppExport SEXP R_FM_create_rand(SEXP pn, SEXP pmin, SEXP pmax)
{
	size_t n;
	double min, max;
	bool ret1, ret2, ret3;
	ret1 = R_get_number<size_t>(pn, n);
	ret2 = R_get_number<double>(pmin, min);
	ret3 = R_get_number<double>(pmax, max);
	if (!ret1 || !ret2 || !ret3) {
		fprintf(stderr, "the arguments aren't of the supported type\n");
		return R_NilValue;
	}

	// TODO let's just use in-memory dense matrix first.
	GetRNGstate();
	vector::ptr v = vector::create(n, get_scalar_type<double>(), true,
			rand_set_operate<double>(min, max));
	PutRNGstate();
	return create_FMR_vector(v->get_raw_store(), "");
}

RcppExport SEXP R_FM_create_seq(SEXP pfrom, SEXP pto, SEXP pby)
{
	// This function always generates a sequence of real numbers.
	double from, to, by;
	bool ret1, ret2, ret3;
	ret1 = R_get_number<double>(pfrom, from);
	ret2 = R_get_number<double>(pto, to);
	ret3 = R_get_number<double>(pby, by);
	if (!ret1 || !ret2 || !ret3) {
		fprintf(stderr, "the arguments aren't of the supported type\n");
		return R_NilValue;
	}

	vector::ptr vec = create_vector<double>(from, to, by);
	return create_FMR_vector(vec->get_raw_store(), "");
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

RcppExport SEXP R_FM_load_matrix(SEXP pmat_file, SEXP pindex_file)
{
	std::string mat_file = CHAR(STRING_ELT(pmat_file, 0));
	std::string index_file = CHAR(STRING_ELT(pindex_file, 0));

	SpM_2d_index::ptr index;
	safs::safs_file index_f(safs::get_sys_RAID_conf(), index_file);
	if (index_f.exist())
		index = SpM_2d_index::safs_load(index_file);
	else
		index = SpM_2d_index::load(index_file);
	if (index == NULL) {
		fprintf(stderr, "can't load index\n");
		return R_NilValue;
	}

	SpM_2d_storage::ptr store;
	safs::safs_file mat_f(safs::get_sys_RAID_conf(), mat_file);
	if (mat_f.exist())
		store = SpM_2d_storage::safs_load(mat_file, index);
	else
		store = SpM_2d_storage::load(mat_file, index);
	if (store == NULL) {
		fprintf(stderr, "can't load matrix file\n");
		return R_NilValue;
	}
	sparse_matrix::ptr mat = sparse_matrix::create(index, store);
	return create_FMR_matrix(mat, "mat_file");
}

/*
 * R has only two data types in matrix multiplication: integer and numeric.
 * So we only need to predefine a small number of basic operations with
 * different types.
 */

static basic_ops_impl<int, int, int> R_basic_ops_II;
// This is a special version, used by multiplication in R.
static basic_ops_impl<int, int, double> R_basic_ops_IID;
static basic_ops_impl<double, int, double> R_basic_ops_DI;
static basic_ops_impl<int, double, double> R_basic_ops_ID;
static basic_ops_impl<double, double, double> R_basic_ops_DD;

static basic_uops_impl<int, int> R_basic_uops_I;
static basic_uops_impl<double, double> R_basic_uops_D;
static basic_uops_impl<bool, bool> R_basic_uops_B;

static SEXP SpMV(sparse_matrix::ptr matrix, vector::ptr vec)
{
	const detail::mem_vec_store &in_vec
		= dynamic_cast<const detail::mem_vec_store &>(vec->get_data());
	detail::mem_vec_store::ptr out_vec = detail::mem_vec_store::create(
			matrix->get_num_rows(), in_vec.get_num_nodes(),
			in_vec.get_type());
	if (vec->is_type<double>()) {
		matrix->multiply<double>(in_vec, *out_vec);
		return create_FMR_vector(out_vec, "");
	}
	else if (vec->is_type<int>()) {
		matrix->multiply<int>(in_vec, *out_vec);
		return create_FMR_vector(out_vec, "");
	}
	else {
		fprintf(stderr, "the input vector has an unsupported type in SpMV\n");
		return R_NilValue;
	}

}

static SEXP SpMM(sparse_matrix::ptr matrix, dense_matrix::ptr right_mat)
{
	if (right_mat->store_layout() != matrix_layout_t::L_ROW) {
		right_mat = right_mat->conv2(matrix_layout_t::L_ROW);
	}

	if (right_mat->is_type<double>()) {
		const detail::mem_matrix_store &in_mat
			= static_cast<const detail::mem_matrix_store &>(right_mat->get_data());
		detail::mem_matrix_store::ptr out_mat = detail::mem_matrix_store::create(
				matrix->get_num_rows(), right_mat->get_num_cols(),
				matrix_layout_t::L_ROW, right_mat->get_type(),
				in_mat.get_num_nodes());
		matrix->multiply<double>(right_mat->get_data(), *out_mat);
		return create_FMR_matrix(dense_matrix::create(out_mat), "");
	}
	else if (right_mat->is_type<int>()) {
		const detail::mem_matrix_store &in_mat
			= static_cast<const detail::mem_matrix_store &>(right_mat->get_data());
		detail::mem_matrix_store::ptr out_mat = detail::mem_matrix_store::create(
				matrix->get_num_rows(), right_mat->get_num_cols(),
				matrix_layout_t::L_ROW, right_mat->get_type(),
				in_mat.get_num_nodes());
		matrix->multiply<int>(right_mat->get_data(), *out_mat);
		return create_FMR_matrix(dense_matrix::create(out_mat), "");
	}
	else {
		fprintf(stderr, "the right matrix has an unsupported type in SpMM\n");
		return R_NilValue;
	}
}

RcppExport SEXP R_FM_multiply_sparse(SEXP pmatrix, SEXP pmat)
{
	sparse_matrix::ptr matrix = get_matrix<sparse_matrix>(pmatrix);
	if (is_vector(pmat)) {
		vector::ptr vec = get_vector(pmat);
		if (!vec->is_in_mem()) {
			fprintf(stderr, "we now only supports in-mem vector for SpMV\n");
			return R_NilValue;
		}
		return SpMV(matrix, vec);
	}
	else {
		dense_matrix::ptr right_mat = get_matrix<dense_matrix>(pmat);
		if (!right_mat->is_in_mem()) {
			fprintf(stderr, "we now only supports in-mem matrix for SpMM\n");
			return R_NilValue;
		}
		// We need to make sure the dense matrix is materialized before SpMM.
		right_mat->materialize_self();
		return SpMM(matrix, right_mat);
	}
}

RcppExport SEXP R_FM_multiply_dense(SEXP pmatrix, SEXP pmat)
{
	dense_matrix::ptr matrix = get_matrix<dense_matrix>(pmatrix);
	dense_matrix::ptr right_mat = get_matrix<dense_matrix>(pmat);
	if (matrix->is_type<int>() && right_mat->is_type<double>())
		matrix = matrix->cast_ele_type(get_scalar_type<double>());
	if (matrix->is_type<double>() && right_mat->is_type<int>())
		right_mat = right_mat->cast_ele_type(get_scalar_type<double>());
	dense_matrix::ptr res = matrix->multiply(*right_mat,
			matrix->store_layout(), true);
	if (res == NULL)
		return R_NilValue;

	if (!res->is_type<double>())
		res = res->cast_ele_type(get_scalar_type<double>());

	bool is_vec = is_vector(pmat);
	if (res && is_vec) {
		return create_FMR_vector(res, "");
	}
	else if (res && !is_vec) {
		return create_FMR_matrix(res, "");
	}
	else
		return R_NilValue;
}

RcppExport SEXP R_FM_conv_matrix(SEXP pvec, SEXP pnrow, SEXP pncol, SEXP pbyrow)
{
	Rcpp::List vec_obj(pvec);
	if (!is_vector(vec_obj)) {
		fprintf(stderr, "The input object isn't a vector\n");
		return R_NilValue;
	}

	size_t nrow = REAL(pnrow)[0];
	size_t ncol = REAL(pncol)[0];
	bool byrow = LOGICAL(pbyrow)[0];
	vector::ptr vec = get_vector(pvec);
	return create_FMR_matrix(vec->conv2mat(nrow, ncol, byrow), "");
}

template<class T, class RType>
void copy_FM2Rmatrix(const dense_matrix &mat, RType *r_vec)
{
	dense_matrix::ptr mem_mat;
	if (!mat.is_in_mem())
		mem_mat = mat.conv_store(true, -1);
	else {
		mem_mat = mat.clone();
		mem_mat->materialize_self();
	}
	const detail::mem_matrix_store &mem_store
		= dynamic_cast<const detail::mem_matrix_store &>(mem_mat->get_data());
	// TODO this is going to be slow. But I don't care about performance
	// for now.
	size_t nrow = mat.get_num_rows();
	size_t ncol = mat.get_num_cols();
	for (size_t i = 0; i < nrow; i++)
		for (size_t j = 0; j < ncol; j++)
			r_vec[i + j * nrow] = mem_store.get<T>(i, j);
}

template<class T, class RType>
void copy_FM2R_mem(dense_matrix::ptr mem_mat, bool is_vec, RType *ret)
{
	if (is_vec) {
		vector::ptr mem_vec = mem_mat->get_col(0);
		if (sizeof(T) == sizeof(RType))
			mem_vec->get_data().copy_to((char *) ret, mem_vec->get_length());
		else {
			std::unique_ptr<T[]> tmp(new T[mem_vec->get_length()]);
			mem_vec->get_data().copy_to((char *) tmp.get(), mem_vec->get_length());
			for (size_t i = 0; i < mem_vec->get_length(); i++)
				ret[i] = tmp[i];
		}
	}
	else
		copy_FM2Rmatrix<T>(*mem_mat, ret);
}

RcppExport SEXP R_FM_copy_FM2R(SEXP pobj, SEXP pRmat)
{
	Rcpp::LogicalVector ret(1);
	if (is_sparse(pobj)) {
		fprintf(stderr, "We can't copy a sparse matrix to an R object\n");
		ret[0] = false;
		return ret;
	}

	dense_matrix::ptr mat = get_matrix<dense_matrix>(pobj);
	assert(mat);
	if (!mat->is_in_mem()) {
		fprintf(stderr, "We only support in-memory matrix right now\n");
		ret[0] = false;
		return ret;
	}

	bool is_vec = is_vector(pobj);
	if (mat->is_type<double>()) {
		copy_FM2R_mem<double, double>(mat, is_vec, REAL(pRmat));
		ret[0] = true;
	}
	else if (mat->is_type<int>()) {
		copy_FM2R_mem<int, int>(mat, is_vec, INTEGER(pRmat));
		ret[0] = true;
	}
	else if (mat->is_type<bool>()) {
		copy_FM2R_mem<bool, int>(mat, is_vec, LOGICAL(pRmat));
		ret[0] = true;
	}
	else {
		fprintf(stderr, "the dense matrix doesn't have a right type\n");
		ret[0] = false;
	}

	return ret;
}

#if 0
template<class T, int SEXPType>
SEXP conv_FM2R_mem(mem_dense_matrix::ptr mem_mat, bool is_vec)
{
	if (is_vec) {
		typename type_mem_vector<T>::ptr mem_vec = type_mem_vector<T>::create(mem_mat);
		Rcpp::Vector<SEXPType> ret(mem_vec->get_length());
		copy_FM2Rvector<T, Rcpp::Vector<SEXPType> >(*mem_vec, ret);
		return ret;
	}
	else {
		Rcpp::Matrix<SEXPType> ret(mem_mat->get_num_rows(),
				mem_mat->get_num_cols());
		copy_FM2Rmatrix<T, Rcpp::Matrix<SEXPType>>(
				*type_mem_dense_matrix<T>::create(mem_mat), ret);
		return ret;
	}
}

RcppExport SEXP R_FM_conv_FM2R(SEXP pobj)
{
	if (is_sparse(pobj)) {
		fprintf(stderr, "We can't convert a sparse matrix to an R object\n");
		return R_NilValue;
	}

	dense_matrix::ptr mat = get_matrix<dense_matrix>(pobj);
	if (!mat->is_in_mem()) {
		fprintf(stderr, "We only support in-memory matrix right now\n");
		return R_NilValue;
	}

	mem_dense_matrix::ptr mem_mat = mem_dense_matrix::cast(mat);
	bool is_vec = is_vector(pobj);
	if (mem_mat->is_type<double>())
		return conv_FM2R_mem<double, REALSXP>(mem_mat, is_vec);
	else if (mem_mat->is_type<int>())
		return conv_FM2R_mem<int, INTSXP>(mem_mat, is_vec);
	else if (mem_mat->is_type<bool>())
		return conv_FM2R_mem<bool, LGLSXP>(mem_mat, is_vec);
	else {
		fprintf(stderr, "the dense matrix doesn't have a right type\n");
		return R_NilValue;
	}
}
#endif

RcppExport SEXP R_FM_conv_RVec2FM(SEXP pobj)
{
	if (R_is_real(pobj)) {
		Rcpp::NumericVector vec(pobj);
		// TODO Is there a way of avoiding the extra memory copy?
		std::unique_ptr<double[]> tmp(new double[vec.size()]);
		for (int i = 0; i < vec.size(); i++)
			tmp[i] = vec[i];

		detail::mem_vec_store::ptr fm_vec = detail::mem_vec_store::create(
				vec.size(), num_nodes, get_scalar_type<double>());
		fm_vec->copy_from((char *) tmp.get(),
				vec.size() * fm_vec->get_entry_size());
		return create_FMR_vector(fm_vec, "");
	}
	else if (R_is_integer(pobj)) {
		Rcpp::IntegerVector vec(pobj);
		// TODO Is there a way of avoiding the extra memory copy?
		std::unique_ptr<int[]> tmp(new int[vec.size()]);
		for (int i = 0; i < vec.size(); i++)
			tmp[i] = vec[i];

		detail::mem_vec_store::ptr fm_vec = detail::mem_vec_store::create(
				vec.size(), num_nodes, get_scalar_type<int>());
		fm_vec->copy_from((char *) tmp.get(),
				vec.size() * fm_vec->get_entry_size());
		return create_FMR_vector(fm_vec, "");
	}
	// TODO handle more types.
	else {
		fprintf(stderr, "The R vector has an unsupported type\n");
		return R_NilValue;
	}
}

RcppExport SEXP R_FM_conv_RMat2FM(SEXP pobj, SEXP pbyrow)
{
	bool byrow = LOGICAL(pbyrow)[0];
	matrix_layout_t layout
		= byrow ? matrix_layout_t::L_ROW : matrix_layout_t::L_COL;
	if (R_is_real(pobj)) {
		Rcpp::NumericMatrix mat(pobj);
		size_t nrow = mat.nrow();
		size_t ncol = mat.ncol();
		detail::mem_matrix_store::ptr fm_mat
			= detail::mem_matrix_store::create(nrow, ncol, layout,
					get_scalar_type<double>(), num_nodes);
		for (size_t i = 0; i < nrow; i++)
			for (size_t j = 0; j < ncol; j++)
				fm_mat->set<double>(i, j, mat(i, j));
		return create_FMR_matrix(dense_matrix::create(fm_mat), "");
	}
	else if (R_is_integer(pobj)) {
		Rcpp::IntegerMatrix mat(pobj);
		size_t nrow = mat.nrow();
		size_t ncol = mat.ncol();
		detail::mem_matrix_store::ptr fm_mat
			= detail::mem_matrix_store::create(nrow, ncol, layout,
					get_scalar_type<int>(), num_nodes);
		for (size_t i = 0; i < nrow; i++)
			for (size_t j = 0; j < ncol; j++)
				fm_mat->set<int>(i, j, mat(i, j));
		return create_FMR_matrix(dense_matrix::create(fm_mat), "");
	}
	// TODO handle more types.
	else {
		fprintf(stderr, "The R vector has an unsupported type\n");
		return R_NilValue;
	}
}

RcppExport SEXP R_FM_transpose(SEXP pmat)
{
	Rcpp::List matrix_obj(pmat);
	if (is_sparse(matrix_obj)) {
		fprintf(stderr, "We don't support transpose a sparse matrix yet\n");
		return R_NilValue;
	}

	dense_matrix::ptr m = get_matrix<dense_matrix>(matrix_obj);
	dense_matrix::ptr tm = m->transpose();
	return create_FMR_matrix(tm, "");
}

RcppExport SEXP R_FM_get_basic_op(SEXP pname)
{
	std::string name = CHAR(STRING_ELT(pname, 0));

	basic_ops::op_idx idx;
	if (name == "add")
		idx = basic_ops::op_idx::ADD;
	else if (name == "sub")
		idx = basic_ops::op_idx::SUB;
	else if (name == "mul")
		idx = basic_ops::op_idx::MUL;
	else if (name == "div")
		idx = basic_ops::op_idx::DIV;
	else if (name == "min")
		idx = basic_ops::op_idx::MIN;
	else if (name == "max")
		idx = basic_ops::op_idx::MAX;
	else if (name == "pow")
		idx = basic_ops::op_idx::POW;
	else if (name == "eq")
		idx = basic_ops::op_idx::EQ;
	else if (name == "gt")
		idx = basic_ops::op_idx::GT;
	else if (name == "ge")
		idx = basic_ops::op_idx::GE;
	else {
		fprintf(stderr, "Unsupported basic operator: %s\n", name.c_str());
		return R_NilValue;
	}

	Rcpp::List ret;
	Rcpp::IntegerVector r_info(1);
	// The index
	r_info[0] = idx;
	// The number of operands
	r_info[1] = 2;
	ret["info"] = r_info;
	ret["name"] = pname;
	ret.attr("class") = "fm.bo";
	return ret;
}

RcppExport SEXP R_FM_get_basic_uop(SEXP pname)
{
	std::string name = CHAR(STRING_ELT(pname, 0));

	basic_uops::op_idx idx;
	if (name == "neg")
		idx = basic_uops::op_idx::NEG;
	else if (name == "sqrt")
		idx = basic_uops::op_idx::SQRT;
	else if (name == "abs")
		idx = basic_uops::op_idx::ABS;
	else if (name == "not")
		idx = basic_uops::op_idx::NOT;
	else {
		fprintf(stderr, "Unsupported basic operator: %s\n", name.c_str());
		return R_NilValue;
	}

	Rcpp::List ret;
	Rcpp::IntegerVector r_info(1);
	// The index
	r_info[0] = idx;
	// The number of operands
	r_info[1] = 1;
	ret["info"] = r_info;
	ret["name"] = pname;
	ret.attr("class") = "fm.bo";
	return ret;
}

static int get_op_idx(const Rcpp::List &fun_obj)
{
	Rcpp::IntegerVector info = fun_obj["info"];
	return info[0];
}

static int get_op_nop(const Rcpp::List &fun_obj)
{
	Rcpp::IntegerVector info = fun_obj["info"];
	return info[1];
}

/*
 * Get a binary operator.
 */
static const bulk_operate *get_op(SEXP pfun, prim_type type1, prim_type type2)
{
	Rcpp::List fun_obj(pfun);
	basic_ops::op_idx bo_idx = (basic_ops::op_idx) get_op_idx(fun_obj);
	int noperands = get_op_nop(fun_obj);
	if (noperands != 2) {
		fprintf(stderr, "This isn't a binary operator\n");
		return NULL;
	}

	basic_ops *ops = NULL;
	if (type1 == prim_type::P_DOUBLE && type2 == prim_type::P_DOUBLE)
		ops = &R_basic_ops_DD;
	else if (type1 == prim_type::P_DOUBLE && type2 == prim_type::P_INTEGER)
		ops = &R_basic_ops_DI;
	else if (type1 == prim_type::P_INTEGER && type2 == prim_type::P_DOUBLE)
		ops = &R_basic_ops_ID;
	else if (type1 == prim_type::P_INTEGER && type2 == prim_type::P_INTEGER)
		ops = &R_basic_ops_II;
	else {
		fprintf(stderr, "wrong type\n");
		return NULL;
	}

	const bulk_operate *op = ops->get_op(bo_idx);
	if (op == NULL) {
		fprintf(stderr, "invalid basic binary operator\n");
		return NULL;
	}
	return op;
}

/*
 * Get a unary operator.
 */
static const bulk_uoperate *get_uop(SEXP pfun, prim_type type)
{
	Rcpp::List fun_obj(pfun);
	basic_uops::op_idx bo_idx = (basic_uops::op_idx) get_op_idx(fun_obj);
	int noperands = get_op_nop(fun_obj);
	if (noperands != 1) {
		fprintf(stderr, "This isn't a unary operator\n");
		return NULL;
	}

	basic_uops *ops = NULL;
	if (type == prim_type::P_DOUBLE)
		ops = &R_basic_uops_D;
	else if (type == prim_type::P_INTEGER)
		ops = &R_basic_uops_I;
	else if (type == prim_type::P_BOOL)
		ops = &R_basic_uops_B;
	else {
		fprintf(stderr, "wrong type\n");
		return NULL;
	}

	const bulk_uoperate *op = ops->get_op(bo_idx);
	if (op == NULL) {
		fprintf(stderr, "invalid basic unary operator\n");
		return NULL;
	}
	return op;
}

static prim_type get_prim_type(SEXP obj)
{
	if (R_is_integer(obj))
		return prim_type::P_INTEGER;
	else if (R_is_real(obj))
		return prim_type::P_DOUBLE;
	else
		return prim_type::NUM_TYPES;
}

RcppExport SEXP R_FM_mapply2(SEXP pfun, SEXP po1, SEXP po2)
{
	Rcpp::List obj1(po1);
	Rcpp::List obj2(po2);
	if (is_sparse(obj1) || is_sparse(obj2)) {
		fprintf(stderr, "mapply2 doesn't support sparse matrix\n");
		return R_NilValue;
	}


	// We only need to test on one vector.
	bool is_vec = is_vector(obj1);
	dense_matrix::ptr m1 = get_matrix<dense_matrix>(obj1);
	dense_matrix::ptr m2 = get_matrix<dense_matrix>(obj2);
	const bulk_operate *op = get_op(pfun, m1->get_type().get_type(),
			m2->get_type().get_type());
	if (op == NULL)
		return R_NilValue;

	dense_matrix::ptr out = m1->mapply2(*m2, bulk_operate::conv2ptr(*op));
	if (out == NULL)
		return R_NilValue;
	else if (is_vec)
		return create_FMR_vector(out, "");
	else
		return create_FMR_matrix(out, "");
}

/*
 * A wrapper class that perform array-element operation.
 * This class converts this binary operation into a unary operation.
 */
template<class T>
class AE_operator: public bulk_uoperate
{
	const bulk_operate &op;
	T v;
public:
	AE_operator(const bulk_operate &_op, T v): op(_op) {
		this->v = v;
		assert(sizeof(v) == op.right_entry_size());
	}

	virtual void runA(size_t num_eles, const void *in_arr,
			void *out_arr) const {
		op.runAE(num_eles, in_arr, &v, out_arr);
	}

	virtual const scalar_type &get_input_type() const {
		return op.get_left_type();
	}

	virtual const scalar_type &get_output_type() const {
		return op.get_output_type();
	}
};

RcppExport SEXP R_FM_mapply2_AE(SEXP pfun, SEXP po1, SEXP po2)
{
	Rcpp::List obj1(po1);
	if (is_sparse(obj1)) {
		fprintf(stderr, "mapply2 doesn't support sparse matrix\n");
		return R_NilValue;
	}

	bool is_vec = is_vector(obj1);
	dense_matrix::ptr m1 = get_matrix<dense_matrix>(obj1);

	const bulk_operate *op = get_op(pfun, m1->get_type().get_type(),
			get_prim_type(po2));
	if (op == NULL)
		return R_NilValue;

	dense_matrix::ptr out;
	if (R_is_real(po2)) {
		double res;
		R_get_number<double>(po2, res);
		out = m1->sapply(std::shared_ptr<bulk_uoperate>(
					new AE_operator<double>(*op, res)));
	}
	else if (R_is_integer(po2)) {
		int res;
		R_get_number<int>(po2, res);
		out = m1->sapply(std::shared_ptr<bulk_uoperate>(
					new AE_operator<int>(*op, res)));
	}
	else {
		fprintf(stderr, "wrong type of the right input\n");
		return R_NilValue;
	}

	if (out == NULL)
		return R_NilValue;
	else if (is_vec)
		return create_FMR_vector(out, "");
	else
		return create_FMR_matrix(out, "");
}

/*
 * A wrapper class that perform element-array operation.
 * This class converts this binary operation into a unary operation.
 */
template<class T>
class EA_operator: public bulk_uoperate
{
	const bulk_operate &op;
	T v;
public:
	EA_operator(const bulk_operate &_op, T v): op(_op) {
		this->v = v;
		assert(sizeof(v) == op.left_entry_size());
	}

	virtual void runA(size_t num_eles, const void *in_arr,
			void *out_arr) const {
		op.runEA(num_eles, &v, in_arr, out_arr);
	}

	virtual const scalar_type &get_input_type() const {
		return op.get_right_type();
	}

	virtual const scalar_type &get_output_type() const {
		return op.get_output_type();
	}
};

RcppExport SEXP R_FM_mapply2_EA(SEXP pfun, SEXP po1, SEXP po2)
{
	Rcpp::List obj2(po2);
	if (is_sparse(obj2)) {
		fprintf(stderr, "mapply2 doesn't support sparse matrix\n");
		return R_NilValue;
	}

	bool is_vec = is_vector(obj2);
	dense_matrix::ptr m2 = get_matrix<dense_matrix>(obj2);

	const bulk_operate *op = get_op(pfun, get_prim_type(po1),
			m2->get_type().get_type());
	if (op == NULL)
		return R_NilValue;

	dense_matrix::ptr out;
	if (R_is_real(po1)) {
		double res;
		R_get_number<double>(po1, res);
		out = m2->sapply(std::shared_ptr<bulk_uoperate>(
					new EA_operator<double>(*op, res)));
	}
	else if (R_is_integer(po1)) {
		int res;
		R_get_number<int>(po1, res);
		out = m2->sapply(std::shared_ptr<bulk_uoperate>(
					new EA_operator<int>(*op, res)));
	}
	else {
		fprintf(stderr, "wrong type of the left input\n");
		return R_NilValue;
	}

	if (out == NULL)
		return R_NilValue;
	else if (is_vec)
		return create_FMR_vector(out, "");
	else
		return create_FMR_matrix(out, "");
}

RcppExport SEXP R_FM_sapply(SEXP pfun, SEXP pobj)
{
	Rcpp::List obj(pobj);
	if (is_sparse(obj)) {
		fprintf(stderr, "sapply doesn't support sparse matrix\n");
		return R_NilValue;
	}

	// We only need to test on one vector.
	bool is_vec = is_vector(obj);
	dense_matrix::ptr m = get_matrix<dense_matrix>(obj);

	const bulk_uoperate *op = get_uop(pfun, m->get_type().get_type());
	if (op == NULL)
		return R_NilValue;

	dense_matrix::ptr out = m->sapply(bulk_uoperate::conv2ptr(*op));
	if (out == NULL)
		return R_NilValue;
	else if (is_vec)
		return create_FMR_vector(out, "");
	else
		return create_FMR_matrix(out, "");
}

template<class T, class ReturnType>
ReturnType matrix_agg(const dense_matrix &mat, const bulk_operate &op)
{
	ReturnType ret(1);
	scalar_variable::ptr res = mat.aggregate(op);
	assert(res->get_type() == get_scalar_type<T>());
	if (res != NULL) {
		ret[0] = *(const T *) res->get_raw();
		return ret;
	}
	else {
		fprintf(stderr, "fail to perform aggregation on the matrix\n");
		return R_NilValue;
	}
}

RcppExport SEXP R_FM_agg(SEXP pfun, SEXP pobj)
{
	Rcpp::List obj1(pobj);
	if (is_sparse(obj1)) {
		fprintf(stderr, "agg doesn't support sparse matrix\n");
		return R_NilValue;
	}

	dense_matrix::ptr m = get_matrix<dense_matrix>(obj1);
	// For aggregation, the left and right operands have the same type.
	const bulk_operate *op = get_op(pfun, m->get_type().get_type(),
			m->get_type().get_type());
	if (op == NULL)
		return R_NilValue;

	if (m->is_type<double>())
		return matrix_agg<double, Rcpp::NumericVector>(*m, *op);
	else if (m->is_type<int>())
		return matrix_agg<int, Rcpp::IntegerVector>(*m, *op);
	else {
		fprintf(stderr, "The matrix has an unsupported type for aggregation\n");
		return R_NilValue;
	}
}

RcppExport SEXP R_FM_matrix_layout(SEXP pmat)
{
	Rcpp::StringVector ret(1);
	if (is_sparse(pmat)) {
		ret[0] = Rcpp::String("adj");
	}
	else {
		dense_matrix::ptr mat = get_matrix<dense_matrix>(pmat);
		if (mat->store_layout() == matrix_layout_t::L_COL)
			ret[0] = Rcpp::String("col");
		else if (mat->store_layout() == matrix_layout_t::L_ROW)
			ret[0] = Rcpp::String("row");
		else
			ret[0] = Rcpp::String("unknown");
	}
	return ret;
}

RcppExport SEXP R_FM_typeof(SEXP pmat)
{
	Rcpp::StringVector ret(1);
	if (is_sparse(pmat)) {
		fprintf(stderr, "Don't support sparse matrix\n");
		return R_NilValue;
	}
	else {
		dense_matrix::ptr mat = get_matrix<dense_matrix>(pmat);
		switch(mat->get_type().get_type()) {
			case prim_type::P_BOOL:
				ret[0] = Rcpp::String("logical");
				break;
			case prim_type::P_INTEGER:
				ret[0] = Rcpp::String("integer");
				break;
			case prim_type::P_DOUBLE:
				ret[0] = Rcpp::String("double");
				break;
			default:
				ret[0] = Rcpp::String("unknown");
		}
	}
	return ret;
}

#if 0
RcppExport SEXP R_FM_set_cols(SEXP pmat, SEXP pidxs, SEXP pvs)
{
	Rcpp::LogicalVector ret(1);
	if (is_sparse(pmat)) {
		fprintf(stderr, "can't write columns to a sparse matrix\n");
		ret[0] = false;
		return ret;
	}

	dense_matrix::ptr mat = get_matrix<dense_matrix>(pmat);
	mem_col_dense_matrix::ptr col_m = mem_col_dense_matrix::cast(mat);
	if (col_m == NULL) {
		ret[0] = false;
		return ret;
	}

	dense_matrix::ptr vs = get_matrix<dense_matrix>(pvs);
	mem_col_dense_matrix::ptr mem_vs = mem_col_dense_matrix::cast(vs);
	if (mem_vs == NULL) {
		ret[0] = false;
		return ret;
	}

	Rcpp::IntegerVector r_idxs(pidxs);
	std::vector<off_t> c_idxs(r_idxs.size());
	for (size_t i = 0; i < c_idxs.size(); i++)
		// R is 1-based indexing, and C/C++ is 0-based.
		c_idxs[i] = r_idxs[i] - 1;

	ret[0] = col_m->set_cols(*mem_vs, c_idxs);;
	return ret;
}
#endif

RcppExport SEXP R_FM_get_cols(SEXP pmat, SEXP pidxs)
{
	if (is_sparse(pmat)) {
		fprintf(stderr, "can't get columns from a sparse matrix\n");
		return R_NilValue;
	}

	Rcpp::IntegerVector r_idxs(pidxs);
	std::vector<off_t> c_idxs(r_idxs.size());
	for (size_t i = 0; i < c_idxs.size(); i++)
		// R is 1-based indexing, and C/C++ is 0-based.
		c_idxs[i] = r_idxs[i] - 1;

	dense_matrix::ptr mat = get_matrix<dense_matrix>(pmat);
	dense_matrix::ptr sub_m = mat->get_cols(c_idxs);
	if (sub_m == NULL)
		return R_NilValue;
	else
		return create_FMR_matrix(sub_m, "");
}

RcppExport SEXP R_FM_as_vector(SEXP pmat)
{
	if (is_sparse(pmat)) {
		fprintf(stderr, "can't a sparse matrix to a vector\n");
		return R_NilValue;
	}

	dense_matrix::ptr mat = get_matrix<dense_matrix>(pmat);
	if (mat->get_num_rows() == 1 || mat->get_num_cols() == 1)
		return create_FMR_vector(mat, "");
	else
		return R_NilValue;
}

RcppExport SEXP R_FM_write_obj(SEXP pmat, SEXP pfile)
{
	if (is_sparse(pmat)) {
		fprintf(stderr, "Doesn't support write a sparse matrix to a file\n");
		return R_NilValue;
	}

	dense_matrix::ptr mat = get_matrix<dense_matrix>(pmat);
	std::string file_name = CHAR(STRING_ELT(pfile, 0));
	Rcpp::LogicalVector ret(1);
	ret[0] = dynamic_cast<const detail::mem_matrix_store &>(
			mat->get_data()).write2file(file_name);
	return ret;
}

RcppExport SEXP R_FM_read_obj(SEXP pfile)
{
	std::string file_name = CHAR(STRING_ELT(pfile, 0));
	detail::mem_matrix_store::ptr mat = detail::mem_matrix_store::load(file_name);
	if (mat == NULL)
		return R_NilValue;
	else
		return create_FMR_matrix(dense_matrix::create(mat), "");
}

class R_spm_function: public eigen::spm_function
{
	Rcpp::Function fun;
	SEXP pextra;
	SEXP penv;
	size_t n;
public:
	R_spm_function(SEXP pfun, SEXP pextra, SEXP penv, size_t n): fun(pfun) {
		this->pextra = pextra;
		this->penv = penv;
		this->n = n;
	}

	virtual dense_matrix::ptr run(dense_matrix::ptr x) const {
		// Force R to perform garbage collection. The R garbage collector
		// doesn't work frequently, so it won't clean up the existing
		// dense matrices and use a lot of memory.
		R_gc();
		SEXP r_mat = create_FMR_matrix(x, "x");
		SEXP pret = fun(r_mat, pextra);
		return get_matrix<dense_matrix>(pret);
	}

	virtual size_t get_num_cols() const {
		return n;
	}
	virtual size_t get_num_rows() const {
		return n;
	}
};

template<class T>
T get_scalar(SEXP val)
{
	if (R_is_integer(val))
		return INTEGER(val)[0];
	else
		return REAL(val)[0];
}

RcppExport SEXP R_FM_eigen(SEXP pfunc, SEXP pextra, SEXP psym, SEXP poptions,
		SEXP penv)
{
	Rcpp::LogicalVector sym(psym);
	Rcpp::List options(poptions);

	size_t nev = 1;
	if (options.containsElementNamed("nev"))
		nev = get_scalar<size_t>(options["nev"]);
	std::string solver = "KrylovSchur";
	if (options.containsElementNamed("solver"))
		solver = CHAR(STRING_ELT(options["solver"], 0));

	eigen::eigen_options opts(nev, solver);
	if (options.containsElementNamed("tol"))
		opts.tol = REAL(options["tol"])[0];
	if (options.containsElementNamed("num_blocks"))
		opts.num_blocks = get_scalar<int>(options["num_blocks"]);
	if (options.containsElementNamed("max_restarts"))
		opts.max_restarts = get_scalar<int>(options["max_restarts"]);
	if (options.containsElementNamed("max_iters"))
		opts.max_iters = get_scalar<int>(options["max_iters"]);
	if (options.containsElementNamed("block_size"))
		opts.block_size = get_scalar<int>(options["block_size"]);
	if (options.containsElementNamed("which"))
		opts.which = CHAR(STRING_ELT(options["which"], 0));

	if (!options.containsElementNamed("n")) {
		fprintf(stderr,
				"User needs to specify `n' (the size of the eigenproblem)\n");
		return R_NilValue;
	}
	size_t n = get_scalar<size_t>(options["n"]);
	eigen::eigen_res res = eigen::compute_eigen(new R_spm_function(pfunc,
				pextra, penv, n), sym[0], opts);

	// Return the options.

	Rcpp::IntegerVector nev_vec(1);
	nev_vec[0] = opts.nev;
	options["nev"] = nev_vec;

	Rcpp::NumericVector tol_vec(1);
	tol_vec[0] = opts.tol;
	options["tol"] = tol_vec;

	Rcpp::IntegerVector nblocks_vec(1);
	nblocks_vec[0] = opts.num_blocks;
	options["num_blocks"] = nblocks_vec;

	Rcpp::IntegerVector max_restarts_vec(1);
	max_restarts_vec[0] = opts.max_restarts;
	options["max_restarts"] = max_restarts_vec;

	Rcpp::IntegerVector max_iters_vec(1);
	max_iters_vec[0] = opts.max_iters;
	options["max_iters"] = max_iters_vec;

	Rcpp::IntegerVector block_size_vec(1);
	block_size_vec[0] = opts.block_size;
	options["block_size"] = block_size_vec;
	
	Rcpp::StringVector solver_str(1);
	solver_str[0] = Rcpp::String(opts.solver);
	options["solver"] = solver_str;

	Rcpp::StringVector which_str(1);
	which_str[0] = Rcpp::String(opts.which);
	options["which"] = which_str;

	Rcpp::List ret;
	Rcpp::NumericVector vals(res.vals.begin(), res.vals.end());
	ret["vals"] = vals;
	ret["vecs"] = create_FMR_matrix(res.vecs, "evecs");
	ret["options"] = options;
	return ret;
}

RcppExport SEXP R_FM_scale(SEXP pmat, SEXP pvec, SEXP pbyrow)
{
	bool byrow = LOGICAL(pbyrow)[0];
	if (is_sparse(pmat)) {
		fprintf(stderr, "Doesn't support scale rows/cols of a sparse matrix\n");
		return R_NilValue;
	}
	if (!is_vector(pvec)) {
		fprintf(stderr, "The second argument should be a vector\n");
		return R_NilValue;
	}

	dense_matrix::ptr mat = get_matrix<dense_matrix>(pmat);
	vector::ptr vec = get_vector(pvec);
	if (vec == NULL) {
		return R_NilValue;
	}
	dense_matrix::ptr res;
	if (byrow) {
		if (mat->get_num_rows() != vec->get_length()) {
			fprintf(stderr,
					"The length of the vector doesn't match the rows of the matrix\n");
			return R_NilValue;
		}
		res = mat->scale_rows(vec);
	}
	else {
		if (mat->get_num_cols() != vec->get_length()) {
			fprintf(stderr,
					"The length of the vector doesn't match the columns of the matrix\n");
			return R_NilValue;
		}
		res = mat->scale_cols(vec);
	}
	return create_FMR_matrix(res, "scale");
}
