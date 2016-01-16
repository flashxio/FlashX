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
#include "data_frame.h"
#include "sparse_matrix.h"
#include "bulk_operate.h"
#include "bulk_operate_ext.h"
#include "generic_type.h"
#ifdef ENABLE_TRILINOS
#include "eigensolver/eigensolver.h"
#endif
#include "factor.h"
#include "EM_dense_matrix.h"
#include "combined_matrix_store.h"

#include "rutils.h"
#include "fmr_utils.h"
#include "matrix_ops.h"

using namespace fm;

fg::FG_graph::ptr R_FG_get_graph(SEXP pgraph);

static inline bool is_supported_type(const scalar_type &type)
{
	return type == get_scalar_type<int>()
		|| type == get_scalar_type<double>();
}

static inline const scalar_type &get_common_type(const scalar_type &left,
		const scalar_type &right)
{
	if (left == right)
		return left;
	else if (left == get_scalar_type<int>())
		return right;
	else if (right == get_scalar_type<int>())
		return left;
	else {
		fprintf(stderr, "left type: %d, right type: %d\n", left.get_type(),
				right.get_type());
		return get_scalar_type<double>();
	}
}

static inline prim_type get_prim_type(SEXP obj)
{
	if (R_is_integer(obj))
		return prim_type::P_INTEGER;
	else if (R_is_real(obj))
		return prim_type::P_DOUBLE;
	else
		return prim_type::NUM_TYPES;
}

static inline scalar_variable::ptr get_scalar(SEXP po)
{
	if (R_is_real(po)) {
		return scalar_variable::ptr(
				new scalar_variable_impl<double>(REAL(po)[0]));
	}
	else if (R_is_integer(po)) {
		return scalar_variable::ptr(
				new scalar_variable_impl<int>(INTEGER(po)[0]));
	}
	else {
		fprintf(stderr, "The R variable has unsupported type\n");
		return scalar_variable::ptr();
	}
}

template<class EntryType>
dense_matrix::ptr create_dense_matrix(size_t nrow, size_t ncol,
		matrix_layout_t layout, EntryType initv)
{
	return dense_matrix::create_const(initv, nrow, ncol, layout);
}

RcppExport SEXP R_FM_create_vector(SEXP plen, SEXP pinitv)
{
	size_t len = REAL(plen)[0];

	int num_nodes = matrix_conf.get_num_nodes();
	// When there is only one NUMA node, it's better to use SMP vector.
	if (num_nodes == 1)
		num_nodes = -1;
	vector::ptr vec;
	if (R_is_real(pinitv))
		vec = create_rep_vector<double>(len, REAL(pinitv)[0], num_nodes, true);
	else if (R_is_integer(pinitv))
		vec = create_rep_vector<int>(len, INTEGER(pinitv)[0], num_nodes, true);
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
	int num_nodes = matrix_conf.get_num_nodes();
	// When there is only one NUMA node, it's better to use SMP vector.
	if (num_nodes == 1)
		num_nodes = -1;
	vector::ptr v = vector::create(n, get_scalar_type<double>(), num_nodes,
			true, rand_set_operate<double>(min, max));
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

	int num_nodes = matrix_conf.get_num_nodes();
	// When there is only one NUMA node, it's better to use SMP vector.
	if (num_nodes == 1)
		num_nodes = -1;
	vector::ptr vec = create_seq_vector<double>(from, to, by, num_nodes, true);
	return create_FMR_vector(vec->get_raw_store(), "");
}

RcppExport SEXP R_FM_get_matrix_fg(SEXP pgraph)
{
	Rcpp::List graph = Rcpp::List(pgraph);
	Rcpp::LogicalVector res(1);
	fg::FG_graph::ptr fg = R_FG_get_graph(pgraph);
	// TODO does this work if this isn't a binary matrix?
	sparse_matrix::ptr m = sparse_matrix::create(fg, NULL);
	std::string name = graph["name"];
	return create_FMR_matrix(m, name);
}

RcppExport SEXP R_FM_get_dense_matrix(SEXP pname)
{
	std::string mat_file = CHAR(STRING_ELT(pname, 0));
	if (!safs::exist_safs_file(mat_file)) {
		fprintf(stderr, "The dense matrix doesn't exist\n");
		return R_NilValue;
	}
	detail::EM_matrix_store::ptr store = detail::EM_matrix_store::create(
			mat_file);
	if (store == NULL)
		return R_NilValue;
	return create_FMR_matrix(dense_matrix::create(store), "");
}

RcppExport SEXP R_FM_load_matrix_sym(SEXP pmat_file, SEXP pindex_file, SEXP pin_mem)
{
	std::string mat_file = CHAR(STRING_ELT(pmat_file, 0));
	std::string index_file = CHAR(STRING_ELT(pindex_file, 0));
	bool in_mem = LOGICAL(pin_mem)[0];

	SpM_2d_index::ptr index;
	try {
		safs::native_file index_f(index_file);
		if (index_f.exist())
			index = SpM_2d_index::load(index_file);
		else
			index = SpM_2d_index::safs_load(index_file);
	} catch (std::exception &e) {
		fprintf(stderr, "load index: %s\n", e.what());
		return R_NilValue;
	}

	sparse_matrix::ptr mat;
	try {
		if (!safs::exist_safs_file(mat_file)) {
			SpM_2d_storage::ptr store = SpM_2d_storage::load(mat_file, index);
			if (store)
				mat = sparse_matrix::create(index, store);
		}
		else if (in_mem) {
			SpM_2d_storage::ptr store = SpM_2d_storage::safs_load(mat_file, index);
			if (store)
				mat = sparse_matrix::create(index, store);
		}
		else
			mat = sparse_matrix::create(index, safs::create_io_factory(
						mat_file, safs::REMOTE_ACCESS));
	} catch (std::exception &e) {
		fprintf(stderr, "load matrix: %s\n", e.what());
		return R_NilValue;
	}
	return create_FMR_matrix(mat, "mat_file");
}

RcppExport SEXP R_FM_load_matrix_asym(SEXP pmat_file, SEXP pindex_file,
		SEXP ptmat_file, SEXP ptindex_file, SEXP pin_mem)
{
	std::string mat_file = CHAR(STRING_ELT(pmat_file, 0));
	std::string index_file = CHAR(STRING_ELT(pindex_file, 0));
	std::string tmat_file = CHAR(STRING_ELT(ptmat_file, 0));
	std::string tindex_file = CHAR(STRING_ELT(ptindex_file, 0));
	bool in_mem = LOGICAL(pin_mem)[0];

	SpM_2d_index::ptr index;
	SpM_2d_index::ptr tindex;

	try {
		safs::native_file index_f(index_file);
		if (index_f.exist())
			index = SpM_2d_index::load(index_file);
		else
			index = SpM_2d_index::safs_load(index_file);

		safs::native_file tindex_f(tindex_file);
		if (tindex_f.exist())
			tindex = SpM_2d_index::load(tindex_file);
		else
			tindex = SpM_2d_index::safs_load(tindex_file);
	} catch (std::exception &e) {
		fprintf(stderr, "load index: %s\n", e.what());
		return R_NilValue;
	}

	sparse_matrix::ptr mat;
	// If one of the data matrices doesn't exist in SAFS or the user wants
	// to load the sparse matrix to memory.
	if (!safs::exist_safs_file(mat_file) || !safs::exist_safs_file(tmat_file)
			|| in_mem) {
		SpM_2d_storage::ptr store;
		SpM_2d_storage::ptr tstore;
		try {
			if (!safs::exist_safs_file(mat_file))
				store = SpM_2d_storage::load(mat_file, index);
			else
				store = SpM_2d_storage::safs_load(mat_file, index);

			if (!safs::exist_safs_file(tmat_file))
				tstore = SpM_2d_storage::load(tmat_file, tindex);
			else
				tstore = SpM_2d_storage::safs_load(tmat_file, tindex);
		} catch (std::exception &e) {
			fprintf(stderr, "load matrix: %s\n", e.what());
			return R_NilValue;
		}
		mat = sparse_matrix::create(index, store, tindex, tstore);
	}
	// Here both data matrices exist in SAFS.
	else {
		try {
			safs::file_io_factory::shared_ptr io_fac
				= safs::create_io_factory(mat_file, safs::REMOTE_ACCESS);
			safs::file_io_factory::shared_ptr tio_fac
				= safs::create_io_factory(tmat_file, safs::REMOTE_ACCESS);
			mat = sparse_matrix::create(index, io_fac, tindex, tio_fac);
		} catch (std::exception &e) {
			fprintf(stderr, "load matrix: %s\n", e.what());
			return R_NilValue;
		}
	}
	return create_FMR_matrix(mat, "mat_file");
}

static SEXP SpMV(sparse_matrix::ptr matrix, vector::ptr vec)
{
	detail::mem_vec_store::const_ptr in_vec
		= detail::mem_vec_store::cast(vec->get_raw_store());
	detail::vec_store::ptr out_vec = detail::mem_vec_store::create(
			matrix->get_num_rows(), in_vec->get_num_nodes(),
			in_vec->get_type());
	// TODO it only supports a binary matrix right now.
	assert(matrix->get_entry_size() == 0);
	if (vec->is_type<double>()) {
		matrix->multiply<double, bool>(in_vec, out_vec);
		return create_FMR_vector(out_vec, "");
	}
	else if (vec->is_type<int>()) {
		matrix->multiply<int, bool>(in_vec, out_vec);
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
	// The input matrix in the right operand might be a virtual matrix originally.
	// When we convert its data layout, it's definitely a virtual matrix.
	right_mat->materialize_self();

	// TODO it only supports a binary matrix right now.
	assert(matrix->get_entry_size() == 0);
	if (right_mat->is_type<double>()) {
		detail::mem_matrix_store::const_ptr in_mat
			= detail::mem_matrix_store::cast(right_mat->get_raw_store());
		detail::matrix_store::ptr out_mat = detail::mem_matrix_store::create(
				matrix->get_num_rows(), right_mat->get_num_cols(),
				matrix_layout_t::L_ROW, right_mat->get_type(),
				in_mat->get_num_nodes());
		matrix->multiply<double, bool>(in_mat, out_mat);
		return create_FMR_matrix(dense_matrix::create(out_mat), "");
	}
	else if (right_mat->is_type<int>()) {
		detail::mem_matrix_store::const_ptr in_mat
			= detail::mem_matrix_store::cast(right_mat->get_raw_store());
		detail::matrix_store::ptr out_mat = detail::mem_matrix_store::create(
				matrix->get_num_rows(), right_mat->get_num_cols(),
				matrix_layout_t::L_ROW, right_mat->get_type(),
				in_mat->get_num_nodes());
		matrix->multiply<int, bool>(in_mat, out_mat);
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

RcppExport SEXP R_FM_inner_prod_dense(SEXP pmatrix, SEXP pmat,
		SEXP pfun1, SEXP pfun2)
{
	dense_matrix::ptr matrix = get_matrix<dense_matrix>(pmatrix);
	dense_matrix::ptr right_mat = get_matrix<dense_matrix>(pmat);
	if (!is_supported_type(matrix->get_type())
			|| !is_supported_type(right_mat->get_type())) {
		fprintf(stderr, "The input matrices have unsupported type\n");
		return R_NilValue;
	}
	const scalar_type &common_type = get_common_type(matrix->get_type(),
			right_mat->get_type());
	if (common_type != matrix->get_type())
		matrix = matrix->cast_ele_type(common_type);
	if (common_type != right_mat->get_type())
		right_mat = right_mat->cast_ele_type(common_type);

	bulk_operate::const_ptr op1 = fmr::get_op(pfun1,
			matrix->get_type().get_type());
	if (op1 == NULL) {
		fprintf(stderr, "can't find a right form for the left operator\n");
		return R_NilValue;
	}
	bulk_operate::const_ptr op2 = fmr::get_op(pfun2,
			op1->get_output_type().get_type());
	if (op2 == NULL) {
		fprintf(stderr, "can't find a right form for the right operator\n");
		return R_NilValue;
	}

	dense_matrix::ptr res = matrix->inner_prod(*right_mat, op1, op2);
	if (res == NULL)
		return R_NilValue;

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
	Rcpp::S4 vec_obj(pvec);
	if (!is_vector(vec_obj)) {
		fprintf(stderr, "The input object isn't a vector\n");
		return R_NilValue;
	}

	size_t nrow = REAL(pnrow)[0];
	size_t ncol = REAL(pncol)[0];
	bool byrow = LOGICAL(pbyrow)[0];
	vector::ptr vec = get_vector(vec_obj);
	if (vec == NULL) {
		fprintf(stderr, "Can't get the vector\n");
		return R_NilValue;
	}
	dense_matrix::ptr mat = vec->conv2mat(nrow, ncol, byrow);
	if (mat == NULL) {
		fprintf(stderr, "can't convert a vector to a matrix\n");
		return R_NilValue;
	}
	return create_FMR_matrix(mat, "");
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
	int num_nodes = matrix_conf.get_num_nodes();
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
	int num_nodes = matrix_conf.get_num_nodes();
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
	Rcpp::S4 matrix_obj(pmat);
	if (is_sparse(matrix_obj)) {
		sparse_matrix::ptr m = get_matrix<sparse_matrix>(matrix_obj);
		return create_FMR_matrix(m->transpose(), "");
	}
	else {
		dense_matrix::ptr m = get_matrix<dense_matrix>(matrix_obj);
		dense_matrix::ptr tm = m->transpose();
		return create_FMR_matrix(tm, "");
	}
}

RcppExport SEXP R_FM_get_basic_op(SEXP pname)
{
	std::string name = CHAR(STRING_ELT(pname, 0));

	fmr::op_id_t id = fmr::get_op_id(name);
	if (id < 0) {
		fprintf(stderr, "Unsupported basic operator: %s\n", name.c_str());
		return R_NilValue;
	}

	Rcpp::List ret;
	Rcpp::IntegerVector r_info(2);
	// The index
	r_info[0] = id;
	// The number of operands
	r_info[1] = 2;
	ret["info"] = r_info;
	ret["name"] = pname;
	return ret;
}

RcppExport SEXP R_FM_get_basic_uop(SEXP pname)
{
	std::string name = CHAR(STRING_ELT(pname, 0));

	fmr::op_id_t id = fmr::get_uop_id(name);
	if (id < 0) {
		fprintf(stderr, "Unsupported basic operator: %s\n", name.c_str());
		return R_NilValue;
	}

	Rcpp::List ret;
	Rcpp::IntegerVector r_info(2);
	// The index
	r_info[0] = id;
	// The number of operands
	r_info[1] = 1;
	ret["info"] = r_info;
	ret["name"] = pname;
	return ret;
}

RcppExport SEXP R_FM_mapply2(SEXP pfun, SEXP po1, SEXP po2)
{
	Rcpp::S4 obj1(po1);
	Rcpp::S4 obj2(po2);
	if (is_sparse(obj1) || is_sparse(obj2)) {
		fprintf(stderr, "mapply2 doesn't support sparse matrix\n");
		return R_NilValue;
	}


	// We only need to test on one vector.
	bool is_vec = is_vector(obj1);
	dense_matrix::ptr m1 = get_matrix<dense_matrix>(obj1);
	dense_matrix::ptr m2 = get_matrix<dense_matrix>(obj2);
	if (!is_supported_type(m1->get_type())
			|| !is_supported_type(m2->get_type())) {
		fprintf(stderr, "The input matrices have unsupported type\n");
		return R_NilValue;
	}
	const scalar_type &common_type = get_common_type(m1->get_type(),
			m2->get_type());
	if (common_type != m1->get_type())
		m1 = m1->cast_ele_type(common_type);
	if (common_type != m2->get_type())
		m2 = m2->cast_ele_type(common_type);

	bulk_operate::const_ptr op = fmr::get_op(pfun, m1->get_type().get_type());
	if (op == NULL)
		return R_NilValue;

	dense_matrix::ptr out = m1->mapply2(*m2, op);
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
	bulk_operate::const_ptr op;
	T v;
public:
	AE_operator(bulk_operate::const_ptr op, T v) {
		this->op = op;
		this->v = v;
		assert(sizeof(v) == op->right_entry_size());
	}

	virtual void runA(size_t num_eles, const void *in_arr,
			void *out_arr) const {
		op->runAE(num_eles, in_arr, &v, out_arr);
	}

	virtual const scalar_type &get_input_type() const {
		return op->get_left_type();
	}

	virtual const scalar_type &get_output_type() const {
		return op->get_output_type();
	}
};

RcppExport SEXP R_FM_mapply2_AE(SEXP pfun, SEXP po1, SEXP po2)
{
	Rcpp::S4 obj1(po1);
	if (is_sparse(obj1)) {
		fprintf(stderr, "mapply2 doesn't support sparse matrix\n");
		return R_NilValue;
	}

	bool is_vec = is_vector(obj1);
	dense_matrix::ptr m1 = get_matrix<dense_matrix>(obj1);
	if (!is_supported_type(m1->get_type())) {
		fprintf(stderr, "The input matrix have unsupported type\n");
		return R_NilValue;
	}

	// Get the scalar.
	scalar_variable::ptr o2 = get_scalar(po2);
	if (o2 == NULL)
		return R_NilValue;

	const scalar_type &common_type = get_common_type(m1->get_type(),
			o2->get_type());
	if (common_type != m1->get_type())
		m1 = m1->cast_ele_type(common_type);
	if (common_type != o2->get_type())
		o2 = o2->cast_type(common_type);

	bulk_operate::const_ptr op = fmr::get_op(pfun, m1->get_type().get_type());
	if (op == NULL)
		return R_NilValue;

	dense_matrix::ptr out;
	if (m1->get_type() == get_scalar_type<double>()) {
		double val = scalar_variable::get_val<double>(*o2);
		out = m1->sapply(std::shared_ptr<bulk_uoperate>(
					new AE_operator<double>(op, val)));
	}
	else if (m1->get_type() == get_scalar_type<int>()) {
		int val = scalar_variable::get_val<int>(*o2);
		out = m1->sapply(std::shared_ptr<bulk_uoperate>(
					new AE_operator<int>(op, val)));
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
	bulk_operate::const_ptr op;
	T v;
public:
	EA_operator(bulk_operate::const_ptr op, T v) {
		this->op = op;
		this->v = v;
		assert(sizeof(v) == op->left_entry_size());
	}

	virtual void runA(size_t num_eles, const void *in_arr,
			void *out_arr) const {
		op->runEA(num_eles, &v, in_arr, out_arr);
	}

	virtual const scalar_type &get_input_type() const {
		return op->get_right_type();
	}

	virtual const scalar_type &get_output_type() const {
		return op->get_output_type();
	}
};

RcppExport SEXP R_FM_mapply2_EA(SEXP pfun, SEXP po1, SEXP po2)
{
	Rcpp::S4 obj2(po2);
	if (is_sparse(obj2)) {
		fprintf(stderr, "mapply2 doesn't support sparse matrix\n");
		return R_NilValue;
	}

	bool is_vec = is_vector(obj2);
	dense_matrix::ptr m2 = get_matrix<dense_matrix>(obj2);
	if (!is_supported_type(m2->get_type())) {
		fprintf(stderr, "The input matrix have unsupported type\n");
		return R_NilValue;
	}

	scalar_variable::ptr o1 = get_scalar(po1);;
	if (o1 == NULL)
		return R_NilValue;

	const scalar_type &common_type = get_common_type(m2->get_type(),
			o1->get_type());
	if (common_type != m2->get_type())
		m2 = m2->cast_ele_type(common_type);
	if (common_type != o1->get_type())
		o1 = o1->cast_type(common_type);

	bulk_operate::const_ptr op = fmr::get_op(pfun, m2->get_type().get_type());
	if (op == NULL)
		return R_NilValue;

	dense_matrix::ptr out;
	if (m2->get_type() == get_scalar_type<double>()) {
		double val = scalar_variable::get_val<double>(*o1);
		out = m2->sapply(std::shared_ptr<bulk_uoperate>(
					new EA_operator<double>(op, val)));
	}
	else if (m2->get_type() == get_scalar_type<int>()) {
		int val = scalar_variable::get_val<int>(*o1);
		out = m2->sapply(std::shared_ptr<bulk_uoperate>(
					new EA_operator<int>(op, val)));
	}

	if (out == NULL)
		return R_NilValue;
	else if (is_vec)
		return create_FMR_vector(out, "");
	else
		return create_FMR_matrix(out, "");
}

RcppExport SEXP R_FM_mapply2_MV(SEXP po1, SEXP po2, SEXP pmargin, SEXP pfun)
{
	if (is_sparse(po1)) {
		fprintf(stderr, "mapply2_MV doesn't support sparse matrix\n");
		return R_NilValue;
	}
	if (!is_vector(po2)) {
		fprintf(stderr, "the second argument must be a vector\n");
		return R_NilValue;
	}
	dense_matrix::ptr m = get_matrix<dense_matrix>(po1);
	dense_matrix::ptr m2 = get_matrix<dense_matrix>(po2);
	if (!is_supported_type(m->get_type())
			|| !is_supported_type(m2->get_type())) {
		fprintf(stderr, "The input matrices have unsupported type\n");
		return R_NilValue;
	}
	const scalar_type &common_type = get_common_type(m2->get_type(),
			m->get_type());
	if (common_type != m->get_type())
		m = m->cast_ele_type(common_type);
	if (common_type != m2->get_type())
		m2 = m2->cast_ele_type(common_type);

	vector::ptr v = m2->get_col(0);
	if (v == NULL) {
		fprintf(stderr, "can't get the vector\n");
		return R_NilValue;
	}

	int margin = INTEGER(pmargin)[0];
	bulk_operate::const_ptr op = fmr::get_op(pfun, m->get_type().get_type());
	if (op == NULL)
		return R_NilValue;
	dense_matrix::ptr res;
	if (margin == matrix_margin::MAR_ROW)
		res = m->mapply_rows(v, op);
	else if (margin == matrix_margin::MAR_COL)
		res = m->mapply_cols(v, op);
	else {
		fprintf(stderr, "a wrong margin\n");
		return R_NilValue;
	}

	if (res != NULL)
		return create_FMR_matrix(res, "");
	else
		return R_NilValue;
}

RcppExport SEXP R_FM_sapply(SEXP pfun, SEXP pobj)
{
	Rcpp::S4 obj(pobj);
	if (is_sparse(obj)) {
		fprintf(stderr, "sapply doesn't support sparse matrix\n");
		return R_NilValue;
	}

	// We only need to test on one vector.
	bool is_vec = is_vector(obj);
	dense_matrix::ptr m = get_matrix<dense_matrix>(obj);

	bulk_uoperate::const_ptr op = fmr::get_uop(pfun, m->get_type().get_type());
	if (op == NULL)
		return R_NilValue;

	dense_matrix::ptr out = m->sapply(op);
	if (out == NULL)
		return R_NilValue;
	else if (is_vec)
		return create_FMR_vector(out, "");
	else
		return create_FMR_matrix(out, "");
}

template<class T, class ReturnType>
ReturnType matrix_agg(const dense_matrix &mat, agg_operate::const_ptr op)
{
	ReturnType ret(1);
	scalar_variable::ptr res = mat.aggregate(op);
	if (res == NULL) {
		fprintf(stderr, "can't aggregate on the matrix\n");
		return R_NilValue;
	}
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

RcppExport SEXP R_FM_agg(SEXP pobj, SEXP pfun)
{
	Rcpp::S4 obj1(pobj);
	if (is_sparse(obj1)) {
		fprintf(stderr, "agg doesn't support sparse matrix\n");
		return R_NilValue;
	}

	dense_matrix::ptr m = get_matrix<dense_matrix>(obj1);
	if (m->is_type<bool>())
		m = m->cast_ele_type(get_scalar_type<int>());
	agg_operate::const_ptr op = fmr::get_agg_op(pfun, m->get_type());
	if (op == NULL)
		return R_NilValue;

	if (m->is_type<double>())
		return matrix_agg<double, Rcpp::NumericVector>(*m, op);
	else if (m->is_type<int>())
		return matrix_agg<int, Rcpp::IntegerVector>(*m, op);
	else {
		fprintf(stderr, "The matrix has an unsupported type for aggregation\n");
		return R_NilValue;
	}
}

RcppExport SEXP R_FM_agg_mat(SEXP pobj, SEXP pmargin, SEXP pfun)
{
	Rcpp::S4 obj1(pobj);
	if (is_sparse(obj1)) {
		fprintf(stderr, "agg_mat doesn't support sparse matrix\n");
		return R_NilValue;
	}

	dense_matrix::ptr m = get_matrix<dense_matrix>(obj1);
	if (m->is_type<bool>())
		m = m->cast_ele_type(get_scalar_type<int>());
	agg_operate::const_ptr op = fmr::get_agg_op(pfun, m->get_type());
	if (op == NULL)
		return R_NilValue;

	int margin = INTEGER(pmargin)[0];
	if (margin != matrix_margin::MAR_ROW && margin != matrix_margin::MAR_COL) {
		fprintf(stderr, "unknown margin\n");
		return R_NilValue;
	}
	dense_matrix::ptr res = m->aggregate((matrix_margin) margin, op);
	if (res == NULL) {
		fprintf(stderr, "can't aggregate on the matrix\n");
		return R_NilValue;
	}
	else
		return create_FMR_vector(res, "");
}

RcppExport SEXP R_FM_sgroupby(SEXP pvec, SEXP pfun)
{
	if (!is_vector(pvec)) {
		fprintf(stderr, "Doesn't support sgroupby on a matrix\n");
		return R_NilValue;
	}
	vector::ptr vec = get_vector(pvec);
	agg_operate::const_ptr op = fmr::get_agg_op(pfun, vec->get_type());
	data_frame::ptr groupby_res = vec->groupby(op, true);
	return create_FMR_data_frame(groupby_res, "");
}

RcppExport SEXP R_FM_groupby(SEXP pmat, SEXP pmargin, SEXP pfactor, SEXP pfun)
{
	if (is_vector(pmat)) {
		fprintf(stderr, "Doesn't support groupby on a vector\n");
		return R_NilValue;
	}
	if (is_sparse(pmat)) {
		fprintf(stderr, "Doesn't support groupby on a sparse matrix\n");
		return R_NilValue;
	}
	dense_matrix::ptr mat = get_matrix<dense_matrix>(pmat);

	int margin = INTEGER(pmargin)[0];
	if (margin != matrix_margin::MAR_ROW && margin != matrix_margin::MAR_COL) {
		fprintf(stderr, "invalid margin in groupby\n");
		return R_NilValue;
	}

	factor_col_vector::ptr factor = get_factor_vector(pfactor);
	if (factor == NULL) {
		fprintf(stderr, "groupby needs a factor vector\n");
		return R_NilValue;
	}
	if (margin == matrix_margin::MAR_ROW
			&& factor->get_length() != mat->get_num_cols()) {
		fprintf(stderr,
				"the factor vector needs to have the length as #columns");
		return R_NilValue;
	}
	else if (margin == matrix_margin::MAR_COL
			&& factor->get_length() != mat->get_num_rows()) {
		fprintf(stderr,
				"the factor vector needs to have the length as #rows");
		return R_NilValue;
	}

	if (margin == matrix_margin::MAR_ROW) {
		fprintf(stderr, "doesn't support grouping columns\n");
		return R_NilValue;
	}
	agg_operate::const_ptr op = fmr::get_agg_op(pfun, mat->get_type());
	dense_matrix::ptr groupby_res = mat->groupby_row(factor, op);
	if (groupby_res == NULL)
		return R_NilValue;
	else
		return create_FMR_matrix(groupby_res, "");
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

RcppExport SEXP R_FM_get_submat(SEXP pmat, SEXP pmargin, SEXP pidxs)
{
	if (is_sparse(pmat)) {
		fprintf(stderr, "can't get a submatrix from a sparse matrix\n");
		return R_NilValue;
	}

	int margin = INTEGER(pmargin)[0];
	if (margin != matrix_margin::MAR_ROW && margin != matrix_margin::MAR_COL) {
		fprintf(stderr, "the margin has invalid value\n");
		return R_NilValue;
	}

	Rcpp::IntegerVector r_idxs(pidxs);
	std::vector<off_t> c_idxs(r_idxs.size());
	for (size_t i = 0; i < c_idxs.size(); i++)
		// R is 1-based indexing, and C/C++ is 0-based.
		c_idxs[i] = r_idxs[i] - 1;

	dense_matrix::ptr mat = get_matrix<dense_matrix>(pmat);
	dense_matrix::ptr sub_m;
	if (margin == matrix_margin::MAR_COL)
		sub_m = mat->get_cols(c_idxs);
	else
		sub_m = mat->get_rows(c_idxs);
	if (sub_m == NULL) {
		fprintf(stderr, "can't get a submatrix from the matrix\n");
		return R_NilValue;
	}
	else
		return create_FMR_matrix(sub_m, "");
}

RcppExport SEXP R_FM_get_rows(SEXP pmat, SEXP pidxs)
{
	if (is_sparse(pmat)) {
		fprintf(stderr, "We don't support get rows from a sparse matrix yet\n");
		return R_NilValue;
	}
	dense_matrix::ptr mat = get_matrix<dense_matrix>(pmat);
	Rcpp::IntegerVector tmp(pidxs);
	std::vector<off_t> idxs(tmp.begin(), tmp.end());
	dense_matrix::ptr res = mat->get_rows(idxs);
	if (res == NULL) {
		fprintf(stderr, "can't get rows from the matrix\n");
		return R_NilValue;
	}
	return create_FMR_matrix(res, "");
}

RcppExport SEXP R_FM_as_vector(SEXP pmat)
{
	if (is_sparse(pmat)) {
		fprintf(stderr, "can't a sparse matrix to a vector\n");
		return R_NilValue;
	}

	dense_matrix::ptr mat = get_matrix<dense_matrix>(pmat);
	if (mat->get_num_cols() == 1)
		return create_FMR_vector(mat, "");
	else if (mat->get_num_rows() == 1)
		return create_FMR_vector(mat->transpose(), "");
	else
		return R_NilValue;
}

RcppExport SEXP R_FM_as_factor_vector(SEXP pmat, SEXP plevels)
{
	if (is_sparse(pmat)) {
		fprintf(stderr, "can't a sparse matrix to a vector\n");
		return R_NilValue;
	}

	dense_matrix::ptr mat = get_matrix<dense_matrix>(pmat);
	if (mat->get_num_rows() > 1 && mat->get_num_cols() > 1) {
		fprintf(stderr, "can't convert a matrix to a factor vector\n");
		return R_NilValue;
	}
	// TODO now I assume the elements in a factor vector are integers.
	if (!mat->is_type<int>()) {
		fprintf(stderr,
				"can't convert a non-integer vector to a factor vector\n");
		return R_NilValue;
	}
	Rcpp::IntegerVector num_levels(plevels);
	int num_levels1 = num_levels[0];
	if (num_levels1 <= 0) {
		scalar_variable::ptr tmp = mat->max();
		num_levels1 = *(int *) tmp->get_raw();
	}
	assert(num_levels1 > 0);
	return create_FMR_factor_vector(mat, num_levels1, "");
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
	detail::matrix_store::ptr store = detail::mem_matrix_store::load(file_name);
	if (store == NULL)
		return R_NilValue;
	else {
		int num_nodes = matrix_conf.get_num_nodes();
		dense_matrix::ptr mat = dense_matrix::create(store);
		if (num_nodes > 1)
			mat = mat->conv_store(true, num_nodes);
		return create_FMR_matrix(mat, "");
	}
}

template<class T>
T get_scalar(SEXP val)
{
	if (R_is_integer(val))
		return INTEGER(val)[0];
	else
		return REAL(val)[0];
}

#ifdef ENABLE_TRILINOS
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

	virtual dense_matrix::ptr run(dense_matrix::ptr &x) const {
		// Force R to perform garbage collection. The R garbage collector
		// doesn't work frequently, so it won't clean up the existing
		// dense matrices and use a lot of memory.
		R_gc();
		SEXP s4_mat = R_create_s4fm(create_FMR_matrix(x, "x"));
		SEXP pret;
		bool success;
		try {
			pret = fun(s4_mat, pextra);
			success = true;
		} catch (Rcpp::eval_error e) {
			std::cerr << "can't eval the multiply function" << std::endl;
			std::cerr << e.what() << std::endl;
			success = false;
		}
		UNPROTECT(2);
		if (success)
			return get_matrix<dense_matrix>(pret);
		else
			return dense_matrix::ptr();
	}

	virtual size_t get_num_cols() const {
		return n;
	}
	virtual size_t get_num_rows() const {
		return n;
	}
};

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

	eigen::eigen_options opts;
	if (!opts.init(nev, solver)) {
		fprintf(stderr, "can't init eigen options\n");
		return R_NilValue;
	}
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
	if (res.vals.empty())
		return R_NilValue;

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
#endif

RcppExport SEXP R_FM_set_materialize_level(SEXP pmat, SEXP plevel)
{
	Rcpp::LogicalVector res(1);
	if (is_sparse(pmat)) {
		fprintf(stderr, "Doesn't support materializing a sparse matrix\n");
		res[0] = false;
		return res;
	}
	Rcpp::IntegerVector level_vec(plevel);
	int level = level_vec[0];
	if (level < materialize_level::MATER_CPU
			|| level > materialize_level::MATER_FULL) {
		fprintf(stderr, "unknown materialization level: %d\n", level);
		res[0] = false;
		return res;
	}
	dense_matrix::ptr mat = get_matrix<dense_matrix>(pmat);
	mat->set_materialize_level((materialize_level) level);
	res[0] = true;
	return res;
}

RcppExport SEXP R_FM_materialize(SEXP pmat)
{
	if (is_sparse(pmat)) {
		fprintf(stderr, "Doesn't support materializing a sparse matrix\n");
		return R_NilValue;
	}
	dense_matrix::ptr mat = get_matrix<dense_matrix>(pmat);
	// I think it's OK to materialize on the original matrix.
	mat->materialize_self();
	return create_FMR_matrix(mat, "");
}

RcppExport SEXP R_FM_conv_layout(SEXP pmat, SEXP pbyrow)
{
	bool byrow = LOGICAL(pbyrow)[0];
	if (is_sparse(pmat)) {
		fprintf(stderr,
				"Doesn't support convert the layout of a sparse matrix\n");
		return R_NilValue;
	}
	dense_matrix::ptr mat = get_matrix<dense_matrix>(pmat);
	dense_matrix::ptr ret;
	if (byrow)
		ret = mat->conv2(matrix_layout_t::L_ROW);
	else
		ret = mat->conv2(matrix_layout_t::L_COL);
	return create_FMR_matrix(ret, "");
}

RcppExport SEXP R_FM_get_layout(SEXP pmat)
{
	Rcpp::StringVector ret(1);
	if (is_sparse(pmat))
		ret[0] = Rcpp::String("row-oriented");
	else {
		dense_matrix::ptr mat = get_matrix<dense_matrix>(pmat);
		if (mat->store_layout() == matrix_layout_t::L_ROW)
			ret[0] = Rcpp::String("row-oriented");
		else
			ret[0] = Rcpp::String("col-oriented");
	}
	return ret;
}

RcppExport SEXP R_FM_is_sym(SEXP pmat)
{
	Rcpp::LogicalVector res(1);
	if (is_sparse(pmat)) {
		sparse_matrix::ptr mat = get_matrix<sparse_matrix>(pmat);
		res[0] = mat->is_symmetric();
	}
	else
		res[0] = false;
	return res;
}

RcppExport SEXP R_FM_rbind(SEXP pmats)
{
	Rcpp::List mats(pmats);
	std::vector<detail::matrix_store::const_ptr> stores(mats.size());
	for (int i = 0; i < mats.size(); i++) {
		if (is_sparse(mats[i])) {
			fprintf(stderr, "can't bind sparse matrix\n");
			return R_NilValue;
		}
		dense_matrix::ptr mat = get_matrix<dense_matrix>(mats[i]);
		if (!mat->is_wide()) {
			fprintf(stderr, "can't rbind two tall matrix\n");
			return R_NilValue;
		}
		stores[i] = mat->get_raw_store();
	}
	// TODO does it matter what layout is here?
	detail::combined_matrix_store::ptr combined
		= detail::combined_matrix_store::create(stores, matrix_layout_t::L_ROW);
	return create_FMR_matrix(dense_matrix::create(combined), "");
}

RcppExport SEXP R_FM_cbind(SEXP pmats)
{
	Rcpp::List mats(pmats);
	std::vector<detail::matrix_store::const_ptr> stores(mats.size());
	for (int i = 0; i < mats.size(); i++) {
		if (is_sparse(mats[i])) {
			fprintf(stderr, "can't bind sparse matrix\n");
			return R_NilValue;
		}
		dense_matrix::ptr mat = get_matrix<dense_matrix>(mats[i]);
		if (mat->is_wide()) {
			fprintf(stderr, "can't cbind two wide matrix\n");
			return R_NilValue;
		}
		stores[i] = mat->get_raw_store();
	}
	// TODO does it matter what layout is here?
	detail::combined_matrix_store::ptr combined
		= detail::combined_matrix_store::create(stores, matrix_layout_t::L_COL);
	return create_FMR_matrix(dense_matrix::create(combined), "");
}

template<class T>
class ifelse2_op: public bulk_operate
{
	T no;

	ifelse2_op(T no) {
		this->no = no;
	}
public:
	static const_ptr create(scalar_variable::ptr no) {
		T val = scalar_variable::get_val<T>(*no);
		return const_ptr(new ifelse2_op<T>(val));
	}

	/*
	 * This performs element-wise operation on two input arrays, and stores
	 * the result on the output array.
	 */
	virtual void runAA(size_t num_eles, const void *left_arr,
			const void *right_arr, void *output_arr) const {
		const bool *test = reinterpret_cast<const bool *>(left_arr);
		const T *yes = reinterpret_cast<const T *>(right_arr);
		T *output = reinterpret_cast<T *>(output_arr);
		for (size_t i = 0; i < num_eles; i++) {
			if (test[i])
				output[i] = yes[i];
			else
				output[i] = no;
		}
	}
	/*
	 * This performs operations on the left input array and the right element,
	 * and stores the result on the output array.
	 */
	virtual void runAE(size_t num_eles, const void *left_arr,
			const void *right, void *output_arr) const {
		throw unsupported_exception("ifelse2_op doesn't support runAE");
	}
	/*
	 * This performs operations on the left element array and the right array,
	 * and stores the result on the output array.
	 */
	virtual void runEA(size_t num_eles, const void *left,
			const void *right_arr, void *output_arr) const {
		throw unsupported_exception("ifelse2_op doesn't support runEA");
	}

	/*
	 * This performs aggregation on the input array, combines the agg result
	 * with the original agg result and stores the result on output.
	 */
	virtual void runAgg(size_t num_eles, const void *left_arr, const void *orig,
			void *output) const {
		throw unsupported_exception("ifelse2_op doesn't support runAgg");
	}

	virtual const scalar_type &get_left_type() const {
		return get_scalar_type<bool>();
	}
	virtual const scalar_type &get_right_type() const {
		return get_scalar_type<T>();
	}
	virtual const scalar_type &get_output_type() const {
		return get_scalar_type<T>();
	}

	virtual std::string get_name() const {
		return std::string("ifelse2_op");
	}
};

/*
 * This version of ifelse only requires test and yes to be FlashMatrix matrices.
 */
RcppExport SEXP R_FM_ifelse2(SEXP ptest, SEXP pyes, SEXP pno)
{
	if (is_sparse(ptest) || is_sparse(pyes)) {
		fprintf(stderr, "ifelse doesn't support sparse matrices\n");
		return R_NilValue;
	}
	dense_matrix::ptr test = get_matrix<dense_matrix>(ptest);
	dense_matrix::ptr yes = get_matrix<dense_matrix>(pyes);
	if (test->get_num_rows() != yes->get_num_rows()
			|| test->get_num_cols() != yes->get_num_cols()) {
		fprintf(stderr, "the size of test and yes has to be the same\n");
		return R_NilValue;
	}

	scalar_variable::ptr no;
	if (R_is_integer(pno)) {
		int val = INTEGER(pno)[0];
		no = scalar_variable::ptr(new scalar_variable_impl<int>(val));
	}
	else if (R_is_real(pno)) {
		double val = REAL(pno)[0];
		no = scalar_variable::ptr(new scalar_variable_impl<double>(val));
	}
	else {
		fprintf(stderr,
				"ifelse2 only works with integer and float currently\n");
		return R_NilValue;
	}
	if (test->get_type() != get_scalar_type<bool>()) {
		fprintf(stderr, "test must be boolean\n");
		return R_NilValue;
	}

	// TODO we should cast type if they are different.
	if (yes->get_type() != no->get_type()) {
		fprintf(stderr,
				"ifelse2 doesn't support yes and no of different types\n");
		return R_NilValue;
	}

	// We need to cast type so that the type of yes and no matches.
	bulk_operate::const_ptr op;
	if (yes->get_type() == get_scalar_type<int>())
		op = ifelse2_op<int>::create(no);
	else if (yes->get_type() == get_scalar_type<double>())
		op = ifelse2_op<double>::create(no);
	else {
		fprintf(stderr, "unsupported type in ifelse2\n");
		return R_NilValue;
	}

	dense_matrix::ptr ret = test->mapply2(*yes, op);
	if (ret == NULL)
		return R_NilValue;
	else if (is_vector(ptest))
		return create_FMR_vector(ret, "");
	else
		return create_FMR_matrix(ret, "");
}

class double_isna_op: public bulk_uoperate
{
public:
	virtual void runA(size_t num_eles, const void *in_arr,
			void *out_arr) const {
		const double *in = reinterpret_cast<const double *>(in_arr);
		bool *out = reinterpret_cast<bool *>(out_arr);
		// is.na in R returns true for both NA and NaN.
		// we should do the same thing.
		for (size_t i = 0; i < num_eles; i++)
			out[i] = ISNAN(in[i]);
	}

	virtual const scalar_type &get_input_type() const {
		return get_scalar_type<double>();
	}

	virtual const scalar_type &get_output_type() const {
		return get_scalar_type<bool>();
	}
};

RcppExport SEXP R_FM_isna(SEXP px)
{
	if (is_sparse(px)) {
		fprintf(stderr, "isna doesn't support sparse matrices\n");
		return R_NilValue;
	}
	dense_matrix::ptr x = get_matrix<dense_matrix>(px);
	if (x->get_type() != get_scalar_type<double>()) {
		fprintf(stderr, "isna only works on float-point matrices\n");
		return R_NilValue;
	}
	dense_matrix::ptr ret
		= x->sapply(bulk_uoperate::const_ptr(new double_isna_op()));
	if (ret == NULL)
		return R_NilValue;
	else if (is_vector(px))
		return create_FMR_vector(ret, "");
	else
		return create_FMR_matrix(ret, "");
}

class double_isnan_op: public bulk_uoperate
{
public:
	virtual void runA(size_t num_eles, const void *in_arr,
			void *out_arr) const {
		const double *in = reinterpret_cast<const double *>(in_arr);
		bool *out = reinterpret_cast<bool *>(out_arr);
		for (size_t i = 0; i < num_eles; i++)
			out[i] = R_IsNaN(in[i]);
	}

	virtual const scalar_type &get_input_type() const {
		return get_scalar_type<double>();
	}

	virtual const scalar_type &get_output_type() const {
		return get_scalar_type<bool>();
	}
};

RcppExport SEXP R_FM_isnan(SEXP px)
{
	if (is_sparse(px)) {
		fprintf(stderr, "isnan doesn't support sparse matrices\n");
		return R_NilValue;
	}
	dense_matrix::ptr x = get_matrix<dense_matrix>(px);
	if (x->get_type() != get_scalar_type<double>()) {
		fprintf(stderr, "isnan only works on float-point matrices\n");
		return R_NilValue;
	}
	dense_matrix::ptr ret
		= x->sapply(bulk_uoperate::const_ptr(new double_isnan_op()));
	if (ret == NULL)
		return R_NilValue;
	else if (is_vector(px))
		return create_FMR_vector(ret, "");
	else
		return create_FMR_matrix(ret, "");
}

void init_flashmatrixr()
{
	fmr::init_udf_ext();
}

RcppExport SEXP R_FM_print_conf()
{
	matrix_conf.print();
	return R_NilValue;
}

RcppExport SEXP R_SAFS_print_conf()
{
	safs::params.print();
	return R_NilValue;
}
