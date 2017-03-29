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
#ifdef USE_PROFILER
#include <gperftools/profiler.h>
#endif
#include <unordered_map>
#include <Rcpp.h>
#include <Rmath.h>
#include <fmr_isna.h>

#include "log.h"
#include "safs_file.h"
#include "io_interface.h"

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
#include "EM_vector.h"
#include "combined_matrix_store.h"
#include "block_matrix.h"
#include "col_vec.h"
#include "project_matrix_store.h"
#include "fm_utils.h"

#include "rutils.h"
#include "fmr_utils.h"
#include "matrix_ops.h"
#include "data_io.h"

using namespace fm;

static inline bool is_supported_type(const scalar_type &type)
{
	return type == get_scalar_type<int>()
		|| type == get_scalar_type<double>();
}

static matrix_layout_t determine_layout(size_t nrow, size_t ncol)
{
	return nrow > ncol ? matrix_layout_t::L_COL : matrix_layout_t::L_ROW;
}

static inline const scalar_type &get_common_type(const scalar_type &left,
		const scalar_type &right)
{
	if (left == right)
		return left;
	else if (left == get_scalar_type<int>() && right == get_scalar_type<double>())
		return right;
	else if (right == get_scalar_type<int>() && left == get_scalar_type<double>())
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
	// boolean values in R are stored as integers.
	else if (R_is_logical(obj))
		return prim_type::P_INTEGER;
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
	else if (R_is_logical(po)) {
		// boolean values in R are stored as integers.
		return scalar_variable::ptr(
				new scalar_variable_impl<int>(LOGICAL(po)[0]));
	}
	else {
		fprintf(stderr, "The R variable has unsupported type\n");
		return scalar_variable::ptr();
	}
}

R_type FM_get_Rtype(Rcpp::S4 obj)
{
	Rcpp::String type = obj.slot("ele_type");
	if (type == "logical")
		return R_type::R_LOGICAL;
	else if (type == "integer")
		return R_type::R_INT;
	else if (type == "double")
		return R_type::R_REAL;
	else
		return R_type::R_NTYPES;
}

R_type trans_FM2R(const scalar_type &type)
{
	if (type == get_scalar_type<int>())
		return R_type::R_INT;
	else if (type == get_scalar_type<double>())
		return R_type::R_REAL;
	else if (type == get_scalar_type<bool>())
		return R_type::R_LOGICAL;
	else {
		fprintf(stderr, "unknown output type\n");
		return R_type::R_NTYPES;
	}
}

RcppExport SEXP R_FM_create_vector(SEXP plen, SEXP pinitv)
{
	size_t len = REAL(plen)[0];

	int num_nodes = matrix_conf.get_num_nodes();
	// When there is only one NUMA node, it's better to use SMP vector.
	if (num_nodes == 1)
		num_nodes = -1;
	vector::ptr vec;
	if (R_is_real(pinitv)) {
		dense_matrix::ptr vec = dense_matrix::create_const<double>(
				REAL(pinitv)[0], len, 1, matrix_layout_t::L_COL);
		return create_FMR_vector(vec, R_type::R_REAL, "");
	}
	else if (R_is_integer(pinitv)) {
		dense_matrix::ptr vec = dense_matrix::create_const<int>(
				INTEGER(pinitv)[0], len, 1, matrix_layout_t::L_COL);
		return create_FMR_vector(vec, R_type::R_INT, "");
	}
	else if (R_is_logical(pinitv)) {
		dense_matrix::ptr vec = dense_matrix::create_const<int>(
				LOGICAL(pinitv)[0], len, 1, matrix_layout_t::L_COL);
		return create_FMR_vector(vec, R_type::R_LOGICAL, "");
	}
	else {
		fprintf(stderr, "The initial value has unsupported type\n");
		return R_NilValue;
	}
}

RcppExport SEXP R_FM_create_randmat(SEXP ptype, SEXP pnrow, SEXP pncol,
		SEXP pin_mem, SEXP pname, SEXP pparams)
{
	std::string type = CHAR(STRING_ELT(ptype, 0));
	std::string name = CHAR(STRING_ELT(pname, 0));
	size_t nrow;
	size_t ncol;
	bool ret1 = R_get_number<size_t>(pnrow, nrow);
	bool ret2 = R_get_number<size_t>(pncol, ncol);
	if (!ret1 || !ret2) {
		fprintf(stderr, "the arguments aren't of the supported type\n");
		return R_NilValue;
	}
	bool in_mem = LOGICAL(pin_mem)[0];
	if (!in_mem && !safs::is_safs_init()) {
		fprintf(stderr, "can't create ext-mem matrix when SAFS is disabled\n");
		return R_NilValue;
	}

	int num_nodes = matrix_conf.get_num_nodes();
	// When there is only one NUMA node, it's better to use SMP vector.
	if (num_nodes == 1)
		num_nodes = -1;

	matrix_layout_t layout = determine_layout(nrow, ncol);
	dense_matrix::ptr mat;
	if (type == "uniform") {
		Rcpp::List params(pparams);
		double min, max;
		bool ret2 = R_get_number<double>(params["min"], min);
		bool ret3 = R_get_number<double>(params["max"], max);
		if (!ret2 || !ret3)
			fprintf(stderr, "min/max aren't of the supported type\n");
		else
			mat = dense_matrix::create_randu<double>(min, max, nrow, ncol,
					layout, num_nodes, in_mem);
	}
	else if (type == "norm") {
		Rcpp::List params(pparams);
		double mu, sigma;
		bool ret2 = R_get_number<double>(params["mu"], mu);
		bool ret3 = R_get_number<double>(params["sigma"], sigma);
		if (!ret2 || !ret3)
			fprintf(stderr, "mu/sigma aren't of the supported type\n");
		else
			mat = dense_matrix::create_randn<double>(mu, sigma, nrow, ncol,
					layout, num_nodes, in_mem);
	}
	else
		fprintf(stderr, "unsupported type\n");

	if (mat) {
		detail::EM_object::const_ptr em_mat
			= std::dynamic_pointer_cast<const detail::EM_object>(
					mat->get_raw_store());
		if (em_mat && !name.empty()) {
			bool ret = em_mat->set_persistent(name);
			if (!ret)
				fprintf(stderr, "Can't set matrix %s persistent\n", name.c_str());
		}
		return create_FMR_matrix(mat, R_type::R_REAL, "");
	}
	else
		return R_NilValue;
}

RcppExport SEXP R_FM_create_seq(SEXP pfrom, SEXP pto, SEXP pby)
{
	int num_nodes = matrix_conf.get_num_nodes();
	// When there is only one NUMA node, it's better to use SMP vector.
	if (num_nodes == 1)
		num_nodes = -1;

	bool ret1, ret2, ret3;
	bool any_double = R_is_real(pfrom) || R_is_real(pto) || R_is_real(pby);
	// If any of the arguments is floating-point, we output a floating-point
	// matrix.
	if (any_double) {
		double from, to, by;
		ret1 = R_get_number<double>(pfrom, from);
		ret2 = R_get_number<double>(pto, to);
		ret3 = R_get_number<double>(pby, by);
		if (!ret1 || !ret2 || !ret3) {
			fprintf(stderr, "the arguments aren't of the supported type\n");
			return R_NilValue;
		}
		vector::ptr vec = create_seq_vector<double>(from, to, by, num_nodes,
				true);
		return create_FMR_vector(vec->get_raw_store(), R_type::R_REAL, "");
	}
	else {
		int from, to, by;
		ret1 = R_get_number<int>(pfrom, from);
		ret2 = R_get_number<int>(pto, to);
		ret3 = R_get_number<int>(pby, by);
		if (!ret1 || !ret2 || !ret3) {
			fprintf(stderr, "the arguments aren't of the supported type\n");
			return R_NilValue;
		}
		vector::ptr vec = create_seq_vector<int>(from, to, by, num_nodes,
				true);
		return create_FMR_vector(vec->get_raw_store(), R_type::R_INT, "");
	}

}

RcppExport SEXP R_FM_create_seq_matrix(SEXP pfrom, SEXP pto, SEXP pnrow,
		SEXP pncol, SEXP pbyrow)
{
	size_t nrow, ncol;
	bool byrow = false;
	bool ret1, ret2, ret3, ret4, ret5;

	int num_nodes = matrix_conf.get_num_nodes();
	// When there is only one NUMA node, it's better to use SMP vector.
	if (num_nodes == 1)
		num_nodes = -1;

	ret3 = R_get_number<size_t>(pnrow, nrow);
	ret4 = R_get_number<size_t>(pncol, ncol);
	ret5 = R_get_number<bool>(pbyrow, byrow);

	// If from or to is floating-point, we output a floating-point
	// matrix.
	bool any_double = R_is_real(pfrom) || R_is_real(pto);
	if (any_double) {
		double from, to;
		ret1 = R_get_number<double>(pfrom, from);
		ret2 = R_get_number<double>(pto, to);
		if (!ret1 || !ret2 || !ret3 || !ret4 || !ret5) {
			fprintf(stderr, "the arguments aren't of the supported type\n");
			return R_NilValue;
		}

		double by = (to - from) / (nrow * ncol - 1);
		dense_matrix::ptr mat = dense_matrix::create_seq<double>(from, by,
				nrow, ncol, determine_layout(nrow, ncol), byrow, num_nodes);
		return create_FMR_matrix(mat, R_type::R_REAL, "");
	}
	else {
		int from, to;
		ret1 = R_get_number<int>(pfrom, from);
		ret2 = R_get_number<int>(pto, to);
		if (!ret1 || !ret2 || !ret3 || !ret4 || !ret5) {
			fprintf(stderr, "the arguments aren't of the supported type\n");
			return R_NilValue;
		}

		int by = (to - from) / (nrow * ncol - 1);
		dense_matrix::ptr mat = dense_matrix::create_seq<int>(from, by,
				nrow, ncol, determine_layout(nrow, ncol), byrow, num_nodes);
		return create_FMR_matrix(mat, R_type::R_INT, "");
	}
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
	return create_FMR_matrix(dense_matrix::create(store),
			// TODO how do we determine the matrix element type?
			trans_FM2R(store->get_type()), "");
}

RcppExport SEXP R_FM_load_dense_matrix(SEXP pname, SEXP pin_mem,
		SEXP pele_type, SEXP pdelim, SEXP pncol, SEXP pmat_name)
{
	Rcpp::StringVector rcpp_mats(pname);
	std::vector<std::string> mat_files(rcpp_mats.begin(), rcpp_mats.end());
	bool in_mem = LOGICAL(pin_mem)[0];
	std::string ele_type = CHAR(STRING_ELT(pele_type, 0));
	std::string delim = CHAR(STRING_ELT(pdelim, 0));
	std::string mat_name = CHAR(STRING_ELT(pmat_name, 0));

	if (!in_mem && !safs::is_safs_init()) {
		fprintf(stderr,
				"SAFS isn't init, can't store a matrix on SAFS\n");
		return R_NilValue;
	}

	dense_matrix::ptr mat;
	if (TYPEOF(pncol) == INTSXP) {
		int ncol = INTEGER(pncol)[0];
		// R uses 4-bytes integers.
		// When it's the maximal integer value, we need to convert it to
		// the maximal value of 8-byte integer to indicate automatic discovery
		// of the number of columns.
		if (ncol == std::numeric_limits<int>::max())
			ncol = std::numeric_limits<size_t>::max();
		mat = read_matrix(mat_files, in_mem, true, ele_type, delim, ncol);
	}
	else {
		std::string cols = CHAR(STRING_ELT(pncol, 0));
		mat = read_matrix(mat_files, in_mem, true, ele_type, delim, cols);
	}
	if (mat == NULL)
		return R_NilValue;

	// If a user provides the name for the matrix and the matrix is stored
	// on disks, let's make it persistent on disks.
	if (!mat_name.empty() && !mat->is_in_mem()) {
		const detail::EM_object *obj
			= dynamic_cast<const detail::EM_object *>(mat->get_raw_store().get());
		if (obj)
			obj->set_persistent(mat_name);
	}

	return create_FMR_matrix(mat, trans_FM2R(mat->get_type()), mat_name);
}

RcppExport SEXP R_FM_load_dense_matrix_bin(SEXP pname, SEXP pin_mem,
		SEXP pnrow, SEXP pncol, SEXP pbyrow, SEXP pele_type, SEXP pmat_name)
{
	std::string file_name = CHAR(STRING_ELT(pname, 0));
	bool in_mem = LOGICAL(pin_mem)[0];
	size_t nrow = REAL(pnrow)[0];
	size_t ncol = REAL(pncol)[0];
	bool byrow = LOGICAL(pbyrow)[0];
	std::string ele_type = CHAR(STRING_ELT(pele_type, 0));
	std::string mat_name = CHAR(STRING_ELT(pmat_name, 0));

	if (!in_mem && !safs::is_safs_init()) {
		fprintf(stderr,
				"SAFS isn't init, can't store a matrix on SAFS\n");
		return R_NilValue;
	}

	matrix_layout_t layout
		= byrow ? matrix_layout_t::L_ROW : matrix_layout_t::L_COL;
	if (!valid_ele_type(ele_type)) {
		fprintf(stderr, "wrong element type\n");
		return R_NilValue;
	}
	const scalar_type &type = get_ele_type(ele_type);

	safs::native_file ext_f(file_name);
	if (!ext_f.exist()) {
		fprintf(stderr, "%s doesn't exist\n", file_name.c_str());
		return R_NilValue;
	}
	if (ext_f.get_size() < (ssize_t) (nrow * ncol * type.get_size())) {
		fprintf(stderr, "%s doesn't contain enough data for the matrix\n",
				file_name.c_str());
		return R_NilValue;
	}

	dense_matrix::ptr mat;
	// Load the matrix to SAFS.
	if (!in_mem) {
		// If a user provides a matrix name, we need to make the matrix
		// persistent on SAFS.
		bool temp = true;
		if (!mat_name.empty())
			temp = false;

		detail::EM_matrix_store::ptr store = detail::EM_matrix_store::load(
				file_name, nrow, ncol, layout, type);
		if (store == NULL) {
			fprintf(stderr, "can't load %s to SAFS\n", file_name.c_str());
			return R_NilValue;
		}

		if (!temp) {
			const detail::EM_object *obj
				= dynamic_cast<const detail::EM_object *>(store.get());
			if (obj)
				obj->set_persistent(mat_name);
		}
		mat = dense_matrix::create(store);
	}
	else {
		detail::mem_matrix_store::ptr store = detail::mem_matrix_store::create(
				nrow, ncol, layout, type, -1);
		FILE *f = fopen(file_name.c_str(), "r");
		if (f == NULL) {
			fprintf(stderr, "can't open %s: %s\n", file_name.c_str(),
					strerror(errno));
			return R_NilValue;
		}
		size_t ret = fread(store->get_raw_arr(), nrow * ncol * type.get_size(),
				1, f);
		if (ret != 1) {
			fprintf(stderr, "can't read %s: %s\n", file_name.c_str(),
					strerror(errno));
			return R_NilValue;
		}
		mat = dense_matrix::create(store);
	}
	if (mat) {
		if (mat->get_type() == get_scalar_type<float>())
			mat = mat->cast_ele_type(get_scalar_type<double>());
		return create_FMR_matrix(mat, trans_FM2R(mat->get_type()), mat_name);
	}
	else
		return R_NilValue;
}

RcppExport SEXP R_FM_load_list_vecs(SEXP pname, SEXP pin_mem, SEXP pele_types,
		SEXP pdelim)
{
	Rcpp::StringVector rcpp_mats(pname);
	std::vector<std::string> mat_files(rcpp_mats.begin(), rcpp_mats.end());
	bool in_mem = LOGICAL(pin_mem)[0];
	Rcpp::StringVector rcpp_types(pele_types);
	std::vector<std::string> ele_types(rcpp_types.begin(), rcpp_types.end());
	std::string delim = CHAR(STRING_ELT(pdelim, 0));

	if (!in_mem && !safs::is_safs_init()) {
		fprintf(stderr,
				"SAFS isn't init, can't store a matrix on SAFS\n");
		return R_NilValue;
	}

	std::vector<ele_parser::const_ptr> parsers;
	for (size_t i = 0; i < ele_types.size(); i++) {
		ele_parser::const_ptr parser = get_ele_parser(ele_types[i]);
		if (parser == NULL) {
			fprintf(stderr, "cannot get a right parser\n");
			return R_NilValue;
		}
		parsers.push_back(parser);
	}
	data_frame::ptr df = read_data_frame(mat_files, in_mem, true, delim, parsers);
	if (df == NULL)
		return R_NilValue;

	Rcpp::List ret;
	for (size_t i = 0; i < df->get_num_vecs(); i++) {
		std::string i_str = itoa(i);
		auto vec = df->get_vec(i);
		ret[i_str] = create_FMR_vector(vec, trans_FM2R(vec->get_type()), "");
	}
	return ret;
}

RcppExport SEXP R_FM_load_spm(SEXP pfile, SEXP pin_mem, SEXP pis_sym,
		SEXP pele_type, SEXP pdelim, SEXP pname)
{
	std::string file = CHAR(STRING_ELT(pfile, 0));
	bool in_mem = LOGICAL(pin_mem)[0];
	bool is_sym = LOGICAL(pis_sym)[0];
	std::string ele_type = CHAR(STRING_ELT(pele_type, 0));
	std::string delim = CHAR(STRING_ELT(pdelim, 0));
	std::string mat_name = CHAR(STRING_ELT(pname, 0));
	const scalar_type *type_p = &get_ele_type(ele_type);

	if (!in_mem && !safs::is_safs_init()) {
		fprintf(stderr,
				"SAFS isn't init, can't store a matrix on SAFS\n");
		return R_NilValue;
	}

	std::vector<ele_parser::const_ptr> parsers;
	parsers.push_back(ele_parser::const_ptr(new int_parser<ele_idx_t>()));
	parsers.push_back(ele_parser::const_ptr(new int_parser<ele_idx_t>()));
	ele_parser::const_ptr attr_parser = get_ele_parser(ele_type);
	if (attr_parser != NULL)
		parsers.push_back(attr_parser);
	std::vector<std::string> files(1, file);
	std::vector<off_t> dup_idxs;
	if (is_sym) {
		for (size_t i = 0; i < parsers.size(); i++)
			dup_idxs.push_back(i);
		dup_idxs[0] = 1;
		dup_idxs[1] = 0;
	}
	// We are going to construct a sparse matrix from the edge list and
	// we will sort the edges any way, so we can read data out of order.
	data_frame::ptr df = read_data_frame(files, in_mem, false, delim, parsers,
			dup_idxs);
	if (df == NULL)
		return R_NilValue;

	sparse_matrix::ptr spm = create_2d_matrix(df,
			block_2d_size(16 * 1024, 16 * 1024), type_p, is_sym, mat_name);
	return create_FMR_matrix(spm, trans_FM2R(spm->get_type()), mat_name);
}

RcppExport SEXP R_FM_load_spm_bin_sym(SEXP pmat_file, SEXP pindex_file, SEXP pin_mem)
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
	return create_FMR_matrix(mat, trans_FM2R(mat->get_type()), "mat_file");
}

RcppExport SEXP R_FM_load_spm_bin_asym(SEXP pmat_file, SEXP pindex_file,
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
	return create_FMR_matrix(mat, trans_FM2R(mat->get_type()), "mat_file");
}

static dense_matrix::ptr SpMM(sparse_matrix::ptr matrix,
		dense_matrix::ptr right_mat)
{
	detail::matrix_store::ptr out_mat = detail::mem_matrix_store::create(
			matrix->get_num_rows(), right_mat->get_num_cols(),
			matrix_layout_t::L_ROW, right_mat->get_type(),
			right_mat->get_raw_store()->get_num_nodes());
	matrix->multiply(right_mat->get_raw_store(), out_mat);
	return dense_matrix::create(out_mat);
}

RcppExport SEXP R_FM_multiply_sparse(SEXP pmatrix, SEXP pmat)
{
	sparse_matrix::ptr matrix = get_matrix<sparse_matrix>(pmatrix);
	if (is_sparse(pmat)) {
		fprintf(stderr, "the right matrix can't be sparse\n");
		return R_NilValue;
	}

	dense_matrix::ptr right_mat = get_matrix<dense_matrix>(pmat);
	if (!is_supported_type(right_mat->get_type())) {
		fprintf(stderr, "multiply doesn't support the type\n");
		return R_NilValue;
	}
	if (!right_mat->is_in_mem()) {
		fprintf(stderr, "we now only supports in-mem matrix for SpMM\n");
		return R_NilValue;
	}
	dense_matrix::ptr ret = SpMM(matrix, right_mat);
	if (ret == NULL)
		return R_NilValue;

	if (is_vector(pmat))
		return create_FMR_vector(ret, trans_FM2R(ret->get_type()), "");
	else
		return create_FMR_matrix(ret, trans_FM2R(ret->get_type()), "");
}

RcppExport SEXP R_FM_multiply_dense(SEXP pmatrix, SEXP pmat)
{
	dense_matrix::ptr matrix = get_matrix<dense_matrix>(pmatrix);
	dense_matrix::ptr right_mat = get_matrix<dense_matrix>(pmat);
	if (!is_supported_type(matrix->get_type())
			|| !is_supported_type(right_mat->get_type())) {
		fprintf(stderr, "Input (%s, %s) of multiply have unsupported type\n",
				matrix->get_type().get_name().c_str(),
				right_mat->get_type().get_name().c_str());
		return R_NilValue;
	}
	if (matrix->is_type<int>() && right_mat->is_type<double>())
		matrix = fmr::cast_Rtype(matrix, FM_get_Rtype(pmatrix), R_type::R_REAL);
	if (matrix->is_type<double>() && right_mat->is_type<int>())
		right_mat = fmr::cast_Rtype(right_mat, FM_get_Rtype(pmat),
				R_type::R_REAL);
	dense_matrix::ptr res = matrix->multiply(*right_mat);
	if (res == NULL)
		return R_NilValue;

	if (!res->is_type<double>()) {
		// TODO this is really unnecessary. But we can't run sapply on an IPW
		// matrix.
		bool ret = res->materialize_self();
		if (!ret) {
			fprintf(stderr, "can't materialize the result matrix\n");
			return R_NilValue;
		}
		res = res->cast_ele_type(get_scalar_type<double>());
	}

	bool is_vec = is_vector(pmat);
	if (res && is_vec) {
		return create_FMR_vector(res, trans_FM2R(res->get_type()), "");
	}
	else if (res && !is_vec) {
		return create_FMR_matrix(res, trans_FM2R(res->get_type()), "");
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
		fprintf(stderr,
				"Input (%s, %s) of inner product have unsupported type\n",
				matrix->get_type().get_name().c_str(),
				right_mat->get_type().get_name().c_str());
		return R_NilValue;
	}
	R_type left_type = FM_get_Rtype(pmatrix);
	R_type right_type = FM_get_Rtype(pmat);
	R_type common_Rtype = get_common_Rtype(left_type, right_type);
	if (common_Rtype != left_type)
		matrix = fmr::cast_Rtype(matrix, left_type, common_Rtype);
	if (common_Rtype != right_type)
		right_mat = fmr::cast_Rtype(right_mat, right_type, common_Rtype);

	auto op1_res = fmr::get_op(pfun1, common_Rtype);
	bulk_operate::const_ptr op1 = op1_res.first;
	if (op1 == NULL) {
		fprintf(stderr, "can't find a right form for the left operator\n");
		return R_NilValue;
	}
	R_type common_out_type;
	if (op1->get_output_type() == get_scalar_type<int>())
		common_out_type = R_type::R_INT;
	else if (op1->get_output_type() == get_scalar_type<double>())
		common_out_type = R_type::R_REAL;
	else {
		fprintf(stderr, "a wrong output type in inner product\n");
		return R_NilValue;
	}
	auto op2_res = fmr::get_op(pfun2, common_out_type);
	bulk_operate::const_ptr op2 = op2_res.first;
	if (op2 == NULL) {
		fprintf(stderr, "can't find a right form for the right operator\n");
		return R_NilValue;
	}

	dense_matrix::ptr res = matrix->inner_prod(*right_mat, op1, op2);
	if (res == NULL)
		return R_NilValue;

	bool is_vec = is_vector(pmat);
	if (res && is_vec) {
		return create_FMR_vector(res, op2_res.second, "");
	}
	else if (res && !is_vec) {
		return create_FMR_matrix(res, op2_res.second, "");
	}
	else
		return R_NilValue;
}

RcppExport SEXP R_FM_create_rep_matrix(SEXP pvec, SEXP pnrow, SEXP pncol,
		SEXP pbyrow)
{
	size_t nrow = REAL(pnrow)[0];
	size_t ncol = REAL(pncol)[0];
	bool byrow = LOGICAL(pbyrow)[0];

	dense_matrix::ptr ret;
	R_type type;
	if (R_is_real(pvec)) {
		ret = dense_matrix::create_const<double>(REAL(pvec)[0], nrow, ncol,
				determine_layout(nrow, ncol));
		type = R_type::R_REAL;
	}
	else if (R_is_integer(pvec)) {
		ret = dense_matrix::create_const<int>(INTEGER(pvec)[0], nrow, ncol,
				determine_layout(nrow, ncol));
		type = R_type::R_INT;
	}
	else if (R_is_logical(pvec)) {
		ret = dense_matrix::create_const<int>(LOGICAL(pvec)[0], nrow, ncol,
				determine_layout(nrow, ncol));
		type = R_type::R_LOGICAL;
	}
	else {
		dense_matrix::ptr mat = get_matrix<dense_matrix>(pvec);
		if (mat == NULL) {
			fprintf(stderr, "we only accept FlashR vectors or R vectors\n");
			return R_NilValue;
		}
		if (mat->get_num_rows() > 1 && mat->get_num_cols() > 1) {
			fprintf(stderr, "we only accept vectors\n");
			return R_NilValue;
		}
		col_vec::ptr vec = col_vec::create(mat);
		vec->move_store(true, -1);
		ret = dense_matrix::create_repeat(vec, nrow, ncol,
				determine_layout(nrow, ncol), byrow, matrix_conf.get_num_nodes());
		type = FM_get_Rtype(pvec);
	}

	if (ret == NULL)
		return R_NilValue;
	else
		return create_FMR_matrix(ret, type, "");
}

template<class T, class RType>
class FM2R_portion_op: public detail::portion_mapply_op
{
	RType *r_vec;
	size_t global_nrow;
public:
	FM2R_portion_op(RType *r_vec, size_t global_nrow): detail::portion_mapply_op(
			0, 0, get_scalar_type<int>()) {
		this->r_vec = r_vec;
		this->global_nrow = global_nrow;
	}

	virtual detail::portion_mapply_op::const_ptr transpose() const {
		fprintf(stderr, "FM2R portion operator doesn't support transpose\n");
		return detail::portion_mapply_op::const_ptr();
	}

	virtual void run(
			const std::vector<detail::local_matrix_store::const_ptr> &ins) const {
		size_t nrow = ins[0]->get_num_rows();
		size_t ncol = ins[0]->get_num_cols();
		const detail::local_col_matrix_store &col_in
			= dynamic_cast<const detail::local_col_matrix_store &>(*ins[0]);
		for (size_t j = 0; j < ncol; j++) {
			const RType *src_col = reinterpret_cast<const RType *>(
					col_in.get_col(j));
			for (size_t i = 0; i < nrow; i++) {
				off_t global_row = i + ins[0]->get_global_start_row();
				off_t global_col = j + ins[0]->get_global_start_col();
				r_vec[global_row + global_col * global_nrow] = src_col[i];
			}
		}
	}

	virtual std::string to_string(
			const std::vector<detail::matrix_store::const_ptr> &mats) const {
		return "FM2R_portion_op";
	}

	virtual bool is_agg() const {
		return false;
	}
};

template<class T, class RType>
bool copy_FM2Rmatrix(dense_matrix::ptr mat, RType *r_vec)
{
	size_t chunk_size = detail::mem_matrix_store::CHUNK_SIZE;
	// If the matrix is in memory and is small, we can copy it directly.
	if ((mat->is_in_mem() && mat->get_num_rows() < chunk_size
			&& mat->get_num_cols() < chunk_size)
			|| mat->get_raw_store()->is_sink()) {
		bool ret = mat->materialize_self();
		if (!ret) {
			fprintf(stderr, "can't materialize the matrix\n");
			return R_NilValue;
		}
		if (mat->store_layout() == matrix_layout_t::L_ROW)
			mat = mat->conv2(matrix_layout_t::L_COL);

		detail::local_col_matrix_store::const_ptr col_lstore
			= std::dynamic_pointer_cast<const detail::local_col_matrix_store>(
					mat->get_raw_store()->get_portion(0));
		size_t nrow = mat->get_num_rows();
		size_t ncol = mat->get_num_cols();
		for (size_t j = 0; j < ncol; j++) {
			const RType *src_col = reinterpret_cast<const RType *>(
					col_lstore->get_col(j));
			for (size_t i = 0; i < nrow; i++)
				r_vec[i + j * nrow] = src_col[i];
		}
	}
	else {
		// R stores data in col-major order
		mat = mat->conv2(matrix_layout_t::L_COL);
		std::vector<detail::matrix_store::const_ptr> mats(1, mat->get_raw_store());
		detail::portion_mapply_op::const_ptr op(
				new FM2R_portion_op<T, RType>(r_vec, mat->get_num_rows()));
		detail::__mapply_portion(mats, op, matrix_layout_t::L_ROW);
	}
	return true;
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
	if (mat == NULL) {
		fprintf(stderr, "cannot get a dense matrix.\n");
		return R_NilValue;
	}
	if (mat->is_type<double>())
		ret[0] = copy_FM2Rmatrix<double,double>(mat, REAL(pRmat));
	else if (mat->is_type<int>())
		ret[0] = copy_FM2Rmatrix<int, int>(mat, INTEGER(pRmat));
	else {
		fprintf(stderr, "the dense matrix doesn't have a right type\n");
		ret[0] = false;
	}

	return ret;
}

template<class T>
T *get_Rdata(SEXP pobj)
{
	fprintf(stderr, "Wrong type for getting R object data");
	return NULL;
}

template<>
int *get_Rdata<int>(SEXP pobj)
{
	return INTEGER(pobj);
}

template<>
double *get_Rdata<double>(SEXP pobj)
{
	return REAL(pobj);
}

template<class T>
dense_matrix::ptr RVec2FM(SEXP pobj)
{
	size_t len = get_length(pobj);
	T *data = get_Rdata<T>(pobj);
	// FlashR stores vectors in matrices.
	// We can assume we only convert small R vectors. so we should store data
	// in a SMP matrix.
	detail::mem_matrix_store::ptr fm = detail::mem_matrix_store::create(
			len, 1, matrix_layout_t::L_COL, get_scalar_type<T>(), -1);
	memcpy(fm->get_raw_arr(), data, len * sizeof(data[0]));
	return dense_matrix::create(fm);
}

template<class T>
dense_matrix::ptr RMat2FM(SEXP pobj)
{
	T *data = get_Rdata<T>(pobj);
	size_t nrow = get_nrows(pobj);
	size_t ncol = get_ncols(pobj);
	size_t len = nrow * ncol;
	// We can assume we only convert small R matrices, so we should store data
	// in a SMP matrix.
	detail::mem_matrix_store::ptr fm = detail::mem_matrix_store::create(
			nrow, ncol, matrix_layout_t::L_COL, get_scalar_type<T>(), -1);
	memcpy(fm->get_raw_arr(), data, len * sizeof(data[0]));
	return dense_matrix::create(fm);
}

RcppExport SEXP R_FM_conv_RVec2FM(SEXP pobj)
{
	if (R_is_real(pobj)) {
		auto fm = RVec2FM<double>(pobj);
		return create_FMR_vector(fm, R_type::R_REAL, "");
	}
	else if (R_is_integer(pobj)) {
		auto fm = RVec2FM<int>(pobj);
		return create_FMR_vector(fm, R_type::R_INT, "");
	}
	else if (R_is_logical(pobj)) {
		auto fm = RVec2FM<int>(pobj);
		return create_FMR_vector(fm, R_type::R_LOGICAL, "");
	}
	// TODO handle more types.
	else {
		fprintf(stderr, "The R vector has an unsupported type\n");
		return R_NilValue;
	}
}

RcppExport SEXP R_FM_conv_RMat2FM(SEXP pobj)
{
	if (R_is_real(pobj)) {
		auto fm = RMat2FM<double>(pobj);
		return create_FMR_matrix(fm, R_type::R_REAL, "");
	}
	else if (R_is_integer(pobj)) {
		auto fm = RMat2FM<int>(pobj);
		return create_FMR_matrix(fm, R_type::R_INT, "");
	}
	else if (R_is_logical(pobj)) {
		auto fm = RMat2FM<int>(pobj);
		return create_FMR_matrix(fm, R_type::R_LOGICAL, "");
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
	Rcpp::List ret;
	if (is_sparse(matrix_obj)) {
		sparse_matrix::ptr m = get_matrix<sparse_matrix>(matrix_obj);
		ret = create_FMR_matrix(m->transpose(), FM_get_Rtype(pmat), "");
	}
	else {
		dense_matrix::ptr m = get_matrix<dense_matrix>(matrix_obj);
		dense_matrix::ptr tm = m->transpose();
		ret = create_FMR_matrix(tm, FM_get_Rtype(pmat), "");
	}
	return ret;
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
	static const_ptr create(bulk_operate::const_ptr op,
			detail::mem_matrix_store::const_ptr v) {
		return const_ptr(new AE_operator<T>(op,
					*reinterpret_cast<const T *>(v->get_raw_arr())));
	}
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
	virtual std::string get_name() const {
		return "mapply_AE";
	}
};

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
	static const_ptr create(bulk_operate::const_ptr op,
			detail::mem_matrix_store::const_ptr v) {
		return const_ptr(new EA_operator<T>(op,
					*reinterpret_cast<const T *>(v->get_raw_arr())));
	}
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
	virtual std::string get_name() const {
		return "mapply_EA";
	}
};

RcppExport SEXP R_FM_mapply2(SEXP pfun, SEXP po1, SEXP po2)
{
	Rcpp::S4 obj1(po1);
	Rcpp::S4 obj2(po2);
	if (is_sparse(obj1) || is_sparse(obj2)) {
		fprintf(stderr, "mapply2 doesn't support sparse matrix\n");
		return R_NilValue;
	}


	// We output a vector only if the two input objects are vectors.
	bool is_vec = is_vector(obj1) && is_vector(obj2);
	dense_matrix::ptr m1 = get_matrix<dense_matrix>(obj1);
	dense_matrix::ptr m2 = get_matrix<dense_matrix>(obj2);
	if (!is_supported_type(m1->get_type())
			|| !is_supported_type(m2->get_type())) {
		fprintf(stderr, "Input (%s, %s) of mapply2 have unsupported type\n",
				m1->get_type().get_name().c_str(),
				m2->get_type().get_name().c_str());
		return R_NilValue;
	}
	R_type left_type = FM_get_Rtype(obj1);
	R_type right_type = FM_get_Rtype(obj2);
	R_type common_Rtype = get_common_Rtype(left_type, right_type);
	if (common_Rtype != left_type)
		m1 = fmr::cast_Rtype(m1, left_type, common_Rtype);
	if (common_Rtype != right_type)
		m2 = fmr::cast_Rtype(m2, right_type, common_Rtype);

	auto op_res = fmr::get_op(pfun, common_Rtype);
	bulk_operate::const_ptr op = op_res.first;
	if (op == NULL)
		return R_NilValue;

	dense_matrix::ptr out;
	// If the input matrices have the same size.
	if (m1->get_num_rows() == m2->get_num_rows()
			&& m1->get_num_cols() == m2->get_num_cols())
		out = m1->mapply2(*m2, op);
	// If the left matrix is actually a scalar.
	else if (m1->get_num_rows() * m1->get_num_cols() == 1) {
		if (m1->is_in_mem())
			m1->materialize_self();
		else
			m1 = m1->conv_store(true, -1);
		auto var = std::static_pointer_cast<const detail::mem_matrix_store>(
				m1->get_raw_store());
		if (m1->get_type() == get_scalar_type<int>())
			out = m2->sapply(EA_operator<int>::create(op, var));
		else if (m1->get_type() == get_scalar_type<double>())
			out = m2->sapply(EA_operator<double>::create(op, var));
	}
	// If the right matrix is actually a scalar
	else if (m2->get_num_rows() * m2->get_num_cols() == 1) {
		if (m2->is_in_mem())
			m2->materialize_self();
		else
			m2 = m2->conv_store(true, -1);
		auto var = std::static_pointer_cast<const detail::mem_matrix_store>(
				m2->get_raw_store());
		if (m1->get_type() == get_scalar_type<int>())
			out = m1->sapply(AE_operator<int>::create(op, var));
		else if (m1->get_type() == get_scalar_type<double>())
			out = m1->sapply(AE_operator<double>::create(op, var));
	}
	// If the left matrix is actually a vector.
	else if (m1->get_num_rows() == m2->get_num_rows()
			&& m1->get_num_cols() == 1) {
		m1 = dense_matrix::create_repeat(col_vec::create(m1),
				m2->get_num_rows(), m2->get_num_cols(),
				matrix_layout_t::L_COL, false,
				m1->get_raw_store()->get_num_nodes());
		out = m1->mapply2(*m2, op);
	}
	// If the right matrix is actually a vector.
	else if (m1->get_num_rows() == m2->get_num_rows()
			&& m2->get_num_cols() == 1) {
		out = m1->mapply_cols(col_vec::create(m2), op);
	}
	else {
		fprintf(stderr, "The shape of the input matrices doesn't match");
		return R_NilValue;
	}

	if (out == NULL)
		return R_NilValue;
	else if (is_vec)
		return create_FMR_vector(out, op_res.second, "");
	else
		return create_FMR_matrix(out, op_res.second, "");
}

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
	R_type common_Rtype = get_common_Rtype(FM_get_Rtype(obj1),
			R_get_type(po2));

	auto op_res = fmr::get_op(pfun, common_Rtype);
	bulk_operate::const_ptr op = op_res.first;
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
		return create_FMR_vector(out, op_res.second, "");
	else
		return create_FMR_matrix(out, op_res.second, "");
}

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
	R_type common_Rtype = get_common_Rtype(R_get_type(po1),
			FM_get_Rtype(obj2));

	auto op_res = fmr::get_op(pfun, common_Rtype);
	bulk_operate::const_ptr op = op_res.first;
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
		return create_FMR_vector(out, op_res.second, "");
	else
		return create_FMR_matrix(out, op_res.second, "");
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
	dense_matrix::ptr v = get_matrix<dense_matrix>(po2);
	if (!is_supported_type(m->get_type())
			|| !is_supported_type(v->get_type())) {
		fprintf(stderr, "Input (%s, %s) of mapply2_MV have unsupported type\n",
				m->get_type().get_name().c_str(),
				v->get_type().get_name().c_str());
		return R_NilValue;
	}
	if (v->get_num_cols() > 1) {
		fprintf(stderr, "The second argument should be a vector\n");
		return R_NilValue;
	}
	R_type left_type = FM_get_Rtype(po1);
	R_type right_type = FM_get_Rtype(po2);
	R_type common_Rtype = get_common_Rtype(left_type, right_type);
	if (common_Rtype != left_type)
		m = fmr::cast_Rtype(m, left_type, common_Rtype);
	if (common_Rtype != right_type)
		v = fmr::cast_Rtype(v, right_type, common_Rtype);

	int margin = INTEGER(pmargin)[0];
	auto op_res = fmr::get_op(pfun, common_Rtype);
	bulk_operate::const_ptr op = op_res.first;
	if (op == NULL)
		return R_NilValue;
	dense_matrix::ptr res;
	if (margin == matrix_margin::MAR_ROW)
		res = m->mapply_rows(col_vec::create(v), op);
	else if (margin == matrix_margin::MAR_COL)
		res = m->mapply_cols(col_vec::create(v), op);
	else {
		fprintf(stderr, "a wrong margin\n");
		return R_NilValue;
	}

	if (res != NULL)
		return create_FMR_matrix(res, op_res.second, "");
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
	if (!is_supported_type(m->get_type())) {
		fprintf(stderr, "The input matrix has unsupported type\n");
		return R_NilValue;
	}

	auto op_res = fmr::get_uop(pfun, FM_get_Rtype(obj));
	bulk_uoperate::const_ptr op = op_res.first;
	if (op == NULL)
		return R_NilValue;

	dense_matrix::ptr out = m->sapply(op);
	if (out == NULL)
		return R_NilValue;
	else if (is_vec)
		return create_FMR_vector(out, op_res.second, "");
	else
		return create_FMR_matrix(out, op_res.second, "");
}

RcppExport SEXP R_FM_apply(SEXP pfun, SEXP pmargin, SEXP pobj)
{
	Rcpp::S4 obj(pobj);
	if (is_sparse(obj)) {
		fprintf(stderr, "apply doesn't support sparse matrix\n");
		return R_NilValue;
	}

	// We only need to test on one vector.
	bool is_vec = is_vector(obj);
	dense_matrix::ptr m = get_matrix<dense_matrix>(obj);
	if (!is_supported_type(m->get_type())) {
		fprintf(stderr, "The input matrix has unsupported type\n");
		return R_NilValue;
	}

	auto op_res = fmr::get_apply_op(pfun, FM_get_Rtype(obj));
	fm::arr_apply_operate::const_ptr op = op_res.first;
	if (op == NULL)
		return R_NilValue;

	int margin = INTEGER(pmargin)[0];
	if (margin != matrix_margin::MAR_ROW && margin != matrix_margin::MAR_COL) {
		fprintf(stderr, "unknown margin\n");
		return R_NilValue;
	}

	dense_matrix::ptr out = m->apply((fm::matrix_margin) margin, op);
	if (out == NULL)
		return R_NilValue;
	else if (is_vec)
		return create_FMR_vector(out, op_res.second, "");
	else
		return create_FMR_matrix(out, op_res.second, "");
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
	if (res != NULL) {
		ret[0] = *(const T *) res->get_raw();
		return ret;
	}
	else {
		fprintf(stderr, "fail to perform aggregation on the matrix\n");
		return R_NilValue;
	}
}

RcppExport SEXP R_FM_agg_lazy(SEXP pobj, SEXP pfun)
{
	Rcpp::S4 obj1(pobj);
	if (is_sparse(obj1)) {
		fprintf(stderr, "agg doesn't support sparse matrix\n");
		return R_NilValue;
	}

	dense_matrix::ptr m = get_matrix<dense_matrix>(obj1);
	if (!is_supported_type(m->get_type())) {
		fprintf(stderr, "The input matrix has unsupported type\n");
		return R_NilValue;
	}
	auto op_res = fmr::get_agg_op(pfun, FM_get_Rtype(obj1));
	agg_operate::const_ptr op = op_res.first;
	if (op == NULL)
		return R_NilValue;

	dense_matrix::ptr res = m->aggregate(matrix_margin::BOTH, op);
	return create_FMR_vector(res, op_res.second, "");
}

RcppExport SEXP R_FM_agg_mat_lazy(SEXP pobj, SEXP pmargin, SEXP pfun)
{
	Rcpp::S4 obj1(pobj);
	if (is_sparse(obj1)) {
		fprintf(stderr, "agg_mat doesn't support sparse matrix\n");
		return R_NilValue;
	}

	dense_matrix::ptr m = get_matrix<dense_matrix>(obj1);
	if (!is_supported_type(m->get_type())) {
		fprintf(stderr, "The input matrix has unsupported type\n");
		return R_NilValue;
	}
	auto op_res = fmr::get_agg_op(pfun, FM_get_Rtype(obj1));
	agg_operate::const_ptr op = op_res.first;
	if (op == NULL)
		return R_NilValue;

	int margin = INTEGER(pmargin)[0];
	if (margin != matrix_margin::MAR_ROW && margin != matrix_margin::MAR_COL) {
		fprintf(stderr, "unknown margin\n");
		return R_NilValue;
	}

	dense_matrix::ptr res = m->aggregate((matrix_margin) margin, op);
	return create_FMR_vector(res, op_res.second, "");
}

RcppExport SEXP R_FM_sgroupby(SEXP pvec, SEXP pfun)
{
	if (!is_vector(pvec)) {
		fprintf(stderr, "Doesn't support sgroupby on a matrix\n");
		return R_NilValue;
	}
	col_vec::ptr vec = get_vector(pvec);
	if (!is_supported_type(vec->get_type())) {
		fprintf(stderr, "The input vector has unsupported type\n");
		return R_NilValue;
	}
	auto op_res = fmr::get_agg_op(pfun, FM_get_Rtype(pvec));
	agg_operate::const_ptr op = op_res.first;
	data_frame::ptr groupby_res = vec->groupby(op, true);
	std::vector<R_type> col_types(2);
	col_types[0] = FM_get_Rtype(pvec);
	col_types[1] = op_res.second;
	return create_FMR_data_frame(groupby_res, col_types, "");
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
	if (!is_supported_type(mat->get_type())) {
		fprintf(stderr, "The input matrix has unsupported type\n");
		return R_NilValue;
	}

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
				"the factor vector needs to have the length as #columns\n");
		return R_NilValue;
	}
	else if (margin == matrix_margin::MAR_COL
			&& factor->get_length() != mat->get_num_rows()) {
		fprintf(stderr,
				"the factor vector needs to have the length as #rows\n");
		return R_NilValue;
	}

	if (margin == matrix_margin::MAR_ROW) {
		fprintf(stderr, "doesn't support grouping columns\n");
		return R_NilValue;
	}
	auto op_res = fmr::get_agg_op(pfun, FM_get_Rtype(pmat));
	agg_operate::const_ptr op = op_res.first;
	dense_matrix::ptr groupby_res = mat->groupby_row(factor, op);
	if (groupby_res == NULL)
		return R_NilValue;
	else
		return create_FMR_matrix(groupby_res, op_res.second, "");
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

RcppExport SEXP R_FM_is_inmem(SEXP pmat)
{
	Rcpp::LogicalVector ret(1);
	if (is_sparse(pmat)) {
		sparse_matrix::ptr mat = get_matrix<sparse_matrix>(pmat);
		// TODO let's assume it's always on SSDs first.
		ret[0] = false;
	}
	else {
		dense_matrix::ptr mat = get_matrix<dense_matrix>(pmat);
		ret[0] = mat->is_in_mem();
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

	dense_matrix::ptr mat = get_matrix<dense_matrix>(pmat);
	dense_matrix::ptr sub_m;

	if (R_is_real(pidxs)) {
		Rcpp::NumericVector r_idxs(pidxs);
		std::vector<off_t> c_idxs(r_idxs.size());
		for (size_t i = 0; i < c_idxs.size(); i++)
			// R is 1-based indexing, and C/C++ is 0-based.
			c_idxs[i] = r_idxs[i] - 1;

		sub_m = margin == matrix_margin::MAR_COL
			? mat->get_cols(c_idxs) : mat->get_rows(c_idxs);
	}
	else {
		dense_matrix::ptr idxs = get_matrix<dense_matrix>(pidxs);
		if (idxs->get_num_rows() > 1 && idxs->get_num_cols() > 1) {
			fprintf(stderr, "the index vector is a matrix\n");
			return R_NilValue;
		}
		if (FM_get_Rtype(pidxs) == R_type::R_LOGICAL)
			idxs = idxs->cast_ele_type(get_scalar_type<bool>());
		col_vec::ptr idx_vec = col_vec::create(idxs);
		sub_m = margin == matrix_margin::MAR_COL
			? mat->get_cols(idx_vec) : mat->get_rows(idx_vec);
	}

	if (sub_m == NULL) {
		fprintf(stderr, "can't get a submatrix from the matrix\n");
		return R_NilValue;
	}
	else
		return create_FMR_matrix(sub_m, FM_get_Rtype(pmat), "");
}

RcppExport SEXP R_FM_set_submat(SEXP pmat, SEXP pmargin, SEXP pidxs, SEXP pdata)
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
	dense_matrix::ptr mat = get_matrix<dense_matrix>(pmat);
	dense_matrix::ptr data = get_matrix<dense_matrix>(pdata);;

	Rcpp::NumericVector r_idxs(pidxs);
	std::vector<off_t> c_idxs(r_idxs.size());
	for (size_t i = 0; i < c_idxs.size(); i++)
		// R is 1-based indexing, and C/C++ is 0-based.
		c_idxs[i] = r_idxs[i] - 1;
	if (margin == matrix_margin::MAR_COL
			&& (c_idxs.size() != data->get_num_cols()
				|| data->get_num_rows() != mat->get_num_rows())) {
		fprintf(stderr, "The new data doesn't have the right dimensions\n");
		return R_NilValue;
	}
	if (margin == matrix_margin::MAR_ROW
			&& (c_idxs.size() != data->get_num_rows()
				|| data->get_num_cols() != mat->get_num_cols())) {
		fprintf(stderr, "The new data doesn't have the right dimensions\n");
		return R_NilValue;
	}

	dense_matrix::ptr new_mat;
	if (margin == matrix_margin::MAR_COL)
		new_mat = mat->set_cols(c_idxs, data);
	else
		new_mat = mat->set_rows(c_idxs, data);

	if (new_mat == NULL)
		return R_NilValue;
	else {
		printf("res: %s\n", new_mat->get_raw_store()->get_name().c_str());
		set_matrix<dense_matrix>(pmat, new_mat);
		return create_FMR_matrix(new_mat, FM_get_Rtype(pmat), "");
	}
}

RcppExport SEXP R_FM_get_vec_eles(SEXP pvec, SEXP pidxs)
{
	Rcpp::NumericVector r_idxs(pidxs);
	std::vector<off_t> c_idxs(r_idxs.size());
	for (size_t i = 0; i < c_idxs.size(); i++)
		// R is 1-based indexing, and C/C++ is 0-based.
		c_idxs[i] = r_idxs[i] - 1;

	dense_matrix::ptr mat = get_matrix<dense_matrix>(pvec);
	dense_matrix::ptr sub_m = mat->get_rows(c_idxs);
	if (sub_m == NULL) {
		fprintf(stderr, "can't get elements from the vector\n");
		return R_NilValue;
	}
	else
		return create_FMR_vector(sub_m, FM_get_Rtype(pvec), "");
}

RcppExport SEXP R_FM_as_vector(SEXP pmat)
{
	if (is_sparse(pmat)) {
		fprintf(stderr, "can't a sparse matrix to a vector\n");
		return R_NilValue;
	}

	Rcpp::S4 rcpp_mat(pmat);
	dense_matrix::ptr mat = get_matrix<dense_matrix>(pmat);
	if (mat->get_num_cols() == 1)
		return create_FMR_vector(mat, FM_get_Rtype(pmat), "");
	else if (mat->get_num_rows() == 1)
		return create_FMR_vector(mat->transpose(), FM_get_Rtype(pmat), "");
	else
		return R_NilValue;
}

RcppExport SEXP R_FM_write_obj(SEXP pmat, SEXP pfile, SEXP ptext, SEXP psep)
{
	if (is_sparse(pmat)) {
		fprintf(stderr, "Doesn't support write a sparse matrix to a file\n");
		return R_NilValue;
	}

	dense_matrix::ptr mat = get_matrix<dense_matrix>(pmat);
	// The input matrix might be a block matrix.
	mat = dense_matrix::create(mat->get_raw_store());
	// To write data to a Linux file, we need to make sure data is stored
	// in memory.
	if (!mat->is_in_mem() || mat->is_virtual())
		mat = mat->conv_store(true, -1);

	std::string sep = CHAR(STRING_ELT(psep, 0));
	std::string file_name = CHAR(STRING_ELT(pfile, 0));
	bool text = LOGICAL(ptext)[0];
	Rcpp::LogicalVector ret(1);
	ret[0] = dynamic_cast<const detail::mem_matrix_store &>(
			mat->get_data()).write2file(file_name, text, sep);
	return ret;
}

RcppExport SEXP R_FM_read_obj(SEXP pfile)
{
	std::string file_name = CHAR(STRING_ELT(pfile, 0));
	detail::matrix_store::const_ptr store = detail::mem_matrix_store::load(
			file_name, matrix_conf.get_num_nodes());
	if (store == NULL)
		return R_NilValue;
	else
		return create_FMR_matrix(dense_matrix::create(store),
				trans_FM2R(store->get_type()), "");
}

template<class T>
T get_scalar(SEXP val)
{
	if (R_is_integer(val))
		return INTEGER(val)[0];
	else if (R_is_logical(val))
		return LOGICAL(val)[0];
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
		SEXP s4_mat = R_create_s4fm(create_FMR_matrix(x,
					trans_FM2R(x->get_type()), "x"));
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
	if (options.containsElementNamed("in.mem"))
		opts.in_mem = get_scalar<bool>(options["in.mem"]);

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

	Rcpp::IntegerVector iter_vec(1);
	iter_vec[0] = res.status.num_iters;
	options["iter"] = iter_vec;

	Rcpp::IntegerVector numop_vec(1);
	numop_vec[0] = res.status.num_ops;
	options["numop"] = numop_vec;

	Rcpp::List ret;
	Rcpp::NumericVector vals(res.vals.begin(), res.vals.end());
	ret["values"] = vals;
	ret["vectors"] = create_FMR_matrix(res.vecs,
			trans_FM2R(res.vecs->get_type()), "evecs");
	ret["options"] = options;
	return ret;
}
#endif

RcppExport SEXP R_FM_set_materialize_level(SEXP pmat, SEXP pcached, SEXP pin_mem)
{
	Rcpp::LogicalVector res(1);
	if (is_sparse(pmat)) {
		fprintf(stderr, "Doesn't support materializing a sparse matrix\n");
		res[0] = false;
		return res;
	}

	bool cached = LOGICAL(pcached)[0];
	materialize_level level
		= cached ? materialize_level::MATER_FULL : materialize_level::MATER_CPU;

	bool in_mem = LOGICAL(pin_mem)[0];
	if (!in_mem && !safs::is_safs_init()) {
		fprintf(stderr,
				"can't materialize a matrix on SAFS when SAFS isn't init\n");
		return R_NilValue;
	}

	dense_matrix::ptr mat = get_matrix<dense_matrix>(pmat);
	if (in_mem == mat->is_in_mem())
		mat->set_materialize_level(level);
	else {
		// The store buffer has to be a tall matrix.
		size_t nrow = std::max(mat->get_num_rows(), mat->get_num_cols());
		size_t ncol = std::min(mat->get_num_rows(), mat->get_num_cols());
		detail::matrix_store::ptr store = detail::matrix_store::create(
				nrow, ncol, mat->store_layout(), mat->get_type(),
				matrix_conf.get_num_nodes(), in_mem);
		mat->set_materialize_level((materialize_level) level, store);
	}
	res[0] = true;
	return res;
}

static SEXP materialize_sparse(const SEXP &pmat)
{
	// don't do anything for a sparse matrix.
	sparse_matrix::ptr mat = get_matrix<sparse_matrix>(pmat);
	Rcpp::S4 rcpp_mat(pmat);
	Rcpp::String name = rcpp_mat.slot("name");
	return create_FMR_matrix(mat, FM_get_Rtype(pmat), name);
}

RcppExport SEXP R_FM_materialize(SEXP pmat)
{
	if (is_sparse(pmat))
		return materialize_sparse(pmat);

	dense_matrix::ptr mat = get_matrix<dense_matrix>(pmat);
	// I think it's OK to materialize on the original matrix.
	bool mater_ret = mat->materialize_self();
	if (!mater_ret) {
		fprintf(stderr, "can't materialize the matrix\n");
		return R_NilValue;
	}

	Rcpp::List ret;
	Rcpp::S4 rcpp_mat(pmat);
	Rcpp::String name = rcpp_mat.slot("name");
	if (is_vector(pmat))
		ret = create_FMR_vector(mat, FM_get_Rtype(rcpp_mat), name);
	else
		ret = create_FMR_matrix(mat, FM_get_Rtype(rcpp_mat), name);
	return ret;
}

RcppExport SEXP R_FM_materialize_list(SEXP plist)
{
	Rcpp::List ret_list;
	Rcpp::List in_list(plist);
	std::vector<int> dense_mat_idxs;
	std::vector<dense_matrix::ptr> dense_mats;
	for (int i = 0; i < in_list.size(); i++) {
		SEXP pmat = in_list[i];
		if (is_sparse(pmat))
			ret_list.push_back(materialize_sparse(pmat));
		else {
			// We collect the dense matrices for materialization.
			dense_matrix::ptr mat = get_matrix<dense_matrix>(pmat);
			dense_mats.push_back(mat);
			dense_mat_idxs.push_back(i);
			// We should have NULL to hold a place in the return list.
			ret_list.push_back(R_NilValue);
		}
	}
	bool ret = materialize(dense_mats);
	if (!ret)
		return R_NilValue;

	// We need to add the materialized matrix to the return list.
	for (size_t i = 0; i < dense_mats.size(); i++) {
		int orig_idx = dense_mat_idxs[i];
		dense_matrix::ptr mat = dense_mats[i];
		SEXP pmat = in_list[orig_idx];

		Rcpp::S4 rcpp_mat(pmat);
		Rcpp::String name = rcpp_mat.slot("name");
		Rcpp::List ret;
		if (is_vector(pmat))
			ret = create_FMR_vector(mat, FM_get_Rtype(rcpp_mat), name);
		else
			ret = create_FMR_matrix(mat, FM_get_Rtype(rcpp_mat), name);
		ret_list[orig_idx] = ret;
	}
	return ret_list;
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
	dense_matrix::ptr ret_mat;
	if (byrow)
		ret_mat = mat->conv2(matrix_layout_t::L_ROW);
	else
		ret_mat = mat->conv2(matrix_layout_t::L_COL);
	return create_FMR_matrix(ret_mat, FM_get_Rtype(pmat), "");
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

RcppExport SEXP R_FM_is_sink(SEXP pmat)
{
	Rcpp::LogicalVector res(1);
	if (is_sparse(pmat))
		res[0] = false;
	else {
		dense_matrix::ptr mat = get_matrix<dense_matrix>(pmat);
		res[0] = mat->get_raw_store()->is_sink();
	}
	return res;
}

const scalar_type *get_castup_type(const std::vector<dense_matrix::ptr> &mats)
{
	size_t i = 1;
	for (; i < mats.size(); i++)
		if (mats[i]->get_type() != mats[0]->get_type())
			break;
	// All matrices have the same type.
	if (i == mats.size())
		return NULL;

	// There are only two types in a dense matrix: int or double.
	for (size_t i = 0; i < mats.size(); i++) {
		assert(mats[i]->get_type() == get_scalar_type<double>()
				|| mats[i]->get_type() == get_scalar_type<int>());
		if (mats[i]->get_type() == get_scalar_type<double>())
			return &get_scalar_type<double>();
	}
	return &get_scalar_type<int>();
}

SEXP fm_bind(SEXP pmats, bool byrow)
{
	Rcpp::List rcpp_mats(pmats);
	std::vector<dense_matrix::ptr> mats(rcpp_mats.size());
	for (int i = 0; i < rcpp_mats.size(); i++) {
		if (is_sparse(rcpp_mats[i])) {
			fprintf(stderr, "can't bind sparse matrix\n");
			return R_NilValue;
		}
		mats[i] = get_matrix<dense_matrix>(rcpp_mats[i]);
	}

	const scalar_type *type = get_castup_type(mats);
	// If some of the matrices have different types, we should cast them
	// first.
	if (type) {
		for (size_t i = 0; i < mats.size(); i++)
			if (mats[i]->get_type() != *type)
				mats[i] = mats[i]->cast_ele_type(*type);
	}
	dense_matrix::ptr combined;
	if (byrow)
		combined = dense_matrix::rbind(mats);
	else
		combined = dense_matrix::cbind(mats);
	if (combined == NULL)
		return R_NilValue;

	return create_FMR_matrix(combined, FM_get_Rtype(rcpp_mats[0]), "");
}

RcppExport SEXP R_FM_rbind(SEXP pmats)
{
	return fm_bind(pmats, true);
}

RcppExport SEXP R_FM_cbind(SEXP pmats)
{
	return fm_bind(pmats, false);
}

template<class BoolType, class T>
class ifelse2_op: public bulk_operate
{
public:
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
	virtual void runAgg(size_t num_eles, const void *left_arr, void *output) const {
		throw unsupported_exception("ifelse2_op doesn't support runAgg");
	}

	virtual const scalar_type &get_left_type() const {
		return get_scalar_type<BoolType>();
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

template<class BoolType, class T>
class ifelse_no_op: public ifelse2_op<BoolType, T>
{
	T no;

	ifelse_no_op(T no) {
		this->no = no;
	}
public:
	static bulk_operate::const_ptr create(scalar_variable::ptr no) {
		T val = scalar_variable::get_val<T>(*no);
		return bulk_operate::const_ptr(new ifelse_no_op<BoolType, T>(val));
	}

	/*
	 * This performs element-wise operation on two input arrays, and stores
	 * the result on the output array.
	 */
	virtual void runAA(size_t num_eles, const void *left_arr,
			const void *right_arr, void *output_arr) const {
		const BoolType *test = reinterpret_cast<const BoolType *>(left_arr);
		const T *yes = reinterpret_cast<const T *>(right_arr);
		T *output = reinterpret_cast<T *>(output_arr);
		for (size_t i = 0; i < num_eles; i++) {
			if (test[i])
				output[i] = yes[i];
			else
				output[i] = no;
		}
	}
};

template<class BoolType, class T>
class ifelse_yes_op: public ifelse2_op<BoolType, T>
{
	T yes;

	ifelse_yes_op(T yes) {
		this->yes = yes;
	}
public:
	static bulk_operate::const_ptr create(scalar_variable::ptr yes) {
		T val = scalar_variable::get_val<T>(*yes);
		return bulk_operate::const_ptr(new ifelse_yes_op<BoolType, T>(val));
	}

	/*
	 * This performs element-wise operation on two input arrays, and stores
	 * the result on the output array.
	 */
	virtual void runAA(size_t num_eles, const void *left_arr,
			const void *right_arr, void *output_arr) const {
		const BoolType *test = reinterpret_cast<const BoolType *>(left_arr);
		const T *no = reinterpret_cast<const T *>(right_arr);
		T *output = reinterpret_cast<T *>(output_arr);
		for (size_t i = 0; i < num_eles; i++) {
			if (test[i])
				output[i] = yes;
			else
				output[i] = no[i];
		}
	}
};

static scalar_variable::ptr conv_R2scalar(SEXP pobj)
{
	if (R_is_logical(pobj)) {
		int val = LOGICAL(pobj)[0];
		return scalar_variable::ptr(new scalar_variable_impl<int>(val));
	}
	else if (R_is_integer(pobj)) {
		int val = INTEGER(pobj)[0];
		return scalar_variable::ptr(new scalar_variable_impl<int>(val));
	}
	else if (R_is_real(pobj)) {
		double val = REAL(pobj)[0];
		return scalar_variable::ptr(new scalar_variable_impl<double>(val));
	}
	else {
		fprintf(stderr,
				"ifelse2 only works with boolean, integer or float currently\n");
		return NULL;
	}
}

/*
 * This version of ifelse only requires test and yes to be FlashMatrix matrices.
 */
RcppExport SEXP R_FM_ifelse_no(SEXP ptest, SEXP pyes, SEXP pno)
{
	if (is_sparse(ptest) || is_sparse(pyes)) {
		fprintf(stderr, "ifelse doesn't support sparse matrices\n");
		return R_NilValue;
	}
	dense_matrix::ptr test = get_matrix<dense_matrix>(ptest);
	dense_matrix::ptr yes = get_matrix<dense_matrix>(pyes);
	if (test->get_num_rows() != yes->get_num_rows()
			|| test->get_num_cols() != yes->get_num_cols()) {
		fprintf(stderr, "test, yes and no don't have the same size\n");
		return R_NilValue;
	}

	scalar_variable::ptr no = conv_R2scalar(pno);
	if (no == NULL)
		return R_NilValue;

	if (test->get_type() != get_scalar_type<int>()) {
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
	test = test->cast_ele_type(get_scalar_type<int>());
	if (yes->get_type() == get_scalar_type<int>())
		op = ifelse_no_op<int, int>::create(no);
	else if (yes->get_type() == get_scalar_type<double>())
		op = ifelse_no_op<int, double>::create(no);
	else {
		fprintf(stderr, "unsupported type in ifelse2\n");
		return R_NilValue;
	}

	dense_matrix::ptr ret = test->mapply2(*yes, op);
	if (ret == NULL)
		return R_NilValue;
	else if (is_vector(ptest))
		return create_FMR_vector(ret, FM_get_Rtype(pyes), "");
	else
		return create_FMR_matrix(ret, FM_get_Rtype(pyes), "");
}

/*
 * This version of ifelse only requires test and no to be FlashMatrix matrices.
 */
RcppExport SEXP R_FM_ifelse_yes(SEXP ptest, SEXP pyes, SEXP pno)
{
	if (is_sparse(ptest) || is_sparse(pno)) {
		fprintf(stderr, "ifelse doesn't support sparse matrices\n");
		return R_NilValue;
	}
	dense_matrix::ptr test = get_matrix<dense_matrix>(ptest);
	dense_matrix::ptr no = get_matrix<dense_matrix>(pno);
	if (test->get_num_rows() != no->get_num_rows()
			|| test->get_num_cols() != no->get_num_cols()) {
		fprintf(stderr, "the size of test and no has to be the same\n");
		return R_NilValue;
	}

	scalar_variable::ptr yes = conv_R2scalar(pyes);
	if (yes == NULL)
		return R_NilValue;

	if (test->get_type() != get_scalar_type<int>()) {
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
	test = test->cast_ele_type(get_scalar_type<int>());
	if (no->get_type() == get_scalar_type<int>())
		op = ifelse_yes_op<int, int>::create(yes);
	else if (no->get_type() == get_scalar_type<double>())
		op = ifelse_yes_op<int, double>::create(yes);
	else {
		fprintf(stderr, "unsupported type in ifelse2\n");
		return R_NilValue;
	}

	dense_matrix::ptr ret = test->mapply2(*no, op);
	if (ret == NULL)
		return R_NilValue;
	else if (is_vector(ptest))
		return create_FMR_vector(ret, FM_get_Rtype(pno), "");
	else
		return create_FMR_matrix(ret, FM_get_Rtype(pno), "");
}

template<class BoolType, class T>
class ifelse_portion_op: public detail::portion_mapply_op
{
public:
	ifelse_portion_op(size_t out_num_rows, size_t out_num_cols,
			const scalar_type &_type): detail::portion_mapply_op(out_num_rows,
				out_num_cols, _type) {
	}

	virtual portion_mapply_op::const_ptr transpose() const {
		return portion_mapply_op::const_ptr(new ifelse_portion_op(
					get_out_num_cols(), get_out_num_rows(), get_output_type()));
	}

	virtual void run(
			const std::vector<detail::local_matrix_store::const_ptr> &ins,
			detail::local_matrix_store &out) const;

	virtual std::string to_string(
			const std::vector<detail::matrix_store::const_ptr> &mats) const {
		assert(mats.size() == 3);
		return std::string("ifelse(") + mats[0]->get_name() + ","
			+ mats[1]->get_name() + "," + mats[2]->get_name() + ")";
	}
};

template<class BoolType, class T>
void ifelse_portion_op<BoolType, T>::run(
		const std::vector<detail::local_matrix_store::const_ptr> &ins,
		detail::local_matrix_store &out) const
{
	assert(ins.size() == 3);
	assert(ins[0]->get_type() == get_scalar_type<BoolType>());
	assert(ins[1]->get_type() == get_scalar_type<T>());
	assert(ins[2]->get_type() == get_scalar_type<T>());
	assert(out.get_type() == get_scalar_type<T>());
	assert(ins[0]->store_layout() == ins[1]->store_layout());
	assert(ins[0]->store_layout() == ins[2]->store_layout());
	assert(ins[0]->store_layout() == out.store_layout());

	if (ins[0]->get_raw_arr() && ins[1]->get_raw_arr()
			&& ins[2]->get_raw_arr() && out.get_raw_arr()) {
		const BoolType *test = reinterpret_cast<const BoolType *>(
				ins[0]->get_raw_arr());
		const T *yes = reinterpret_cast<const T *>(ins[1]->get_raw_arr());
		const T *no = reinterpret_cast<const T *>(ins[2]->get_raw_arr());
		T *res = reinterpret_cast<T *>(out.get_raw_arr());
		size_t len = out.get_num_rows() * out.get_num_cols();
		for (size_t i = 0; i < len; i++) {
			if (test[i])
				res[i] = yes[i];
			else
				res[i] = no[i];
		}
	}
	else if (ins[0]->store_layout() == matrix_layout_t::L_ROW) {
		detail::local_row_matrix_store::const_ptr row_in0
			= std::static_pointer_cast<const detail::local_row_matrix_store>(
					ins[0]);
		detail::local_row_matrix_store::const_ptr row_in1
			= std::static_pointer_cast<const detail::local_row_matrix_store>(
					ins[1]);
		detail::local_row_matrix_store::const_ptr row_in2
			= std::static_pointer_cast<const detail::local_row_matrix_store>(
					ins[2]);
		detail::local_row_matrix_store &row_out
			= static_cast<detail::local_row_matrix_store &>(out);
		for (size_t i = 0; i < ins[0]->get_num_rows(); i++) {
			const BoolType *test = reinterpret_cast<const BoolType *>(
					row_in0->get_row(i));
			const T *yes = reinterpret_cast<const T *>(row_in1->get_row(i));
			const T *no = reinterpret_cast<const T *>(row_in2->get_row(i));
			T *res = reinterpret_cast<T *>(row_out.get_row(i));
			assert(test && yes && no && res);
			for (size_t j = 0; j < row_in0->get_num_cols(); j++) {
				if (test[j])
					res[j] = yes[j];
				else
					res[j] = no[j];
			}
		}
	}
	else {
		detail::local_col_matrix_store::const_ptr col_in0
			= std::static_pointer_cast<const detail::local_col_matrix_store>(
					ins[0]);
		detail::local_col_matrix_store::const_ptr col_in1
			= std::static_pointer_cast<const detail::local_col_matrix_store>(
					ins[1]);
		detail::local_col_matrix_store::const_ptr col_in2
			= std::static_pointer_cast<const detail::local_col_matrix_store>(
					ins[2]);
		detail::local_col_matrix_store &col_out
			= static_cast<detail::local_col_matrix_store &>(out);
		for (size_t i = 0; i < ins[0]->get_num_cols(); i++) {
			const BoolType *test = reinterpret_cast<const BoolType *>(
					col_in0->get_col(i));
			const T *yes = reinterpret_cast<const T *>(col_in1->get_col(i));
			const T *no = reinterpret_cast<const T *>(col_in2->get_col(i));
			T *res = reinterpret_cast<T *>(col_out.get_col(i));
			assert(test && yes && no && res);
			for (size_t j = 0; j < col_in0->get_num_rows(); j++) {
				if (test[j])
					res[j] = yes[j];
				else
					res[j] = no[j];
			}
		}
	}
}

RcppExport SEXP R_FM_ifelse(SEXP ptest, SEXP pyes, SEXP pno)
{
	if (is_sparse(ptest) || is_sparse(pno) || is_sparse(pyes)) {
		fprintf(stderr, "ifelse doesn't support sparse matrices\n");
		return R_NilValue;
	}
	dense_matrix::ptr test = get_matrix<dense_matrix>(ptest);
	dense_matrix::ptr yes = get_matrix<dense_matrix>(pyes);
	dense_matrix::ptr no = get_matrix<dense_matrix>(pno);
	if (test->get_num_rows() != no->get_num_rows()
			|| test->get_num_cols() != no->get_num_cols()
			|| test->get_num_rows() != yes->get_num_rows()
			|| test->get_num_cols() != yes->get_num_cols()) {
		fprintf(stderr, "the size of test, yes and no has to be the same\n");
		return R_NilValue;
	}

	// TODO we should cast type if they are different.
	if (yes->get_type() != no->get_type()) {
		fprintf(stderr,
				"ifelse2 doesn't support yes and no of different types\n");
		return R_NilValue;
	}

	test = test->cast_ele_type(get_scalar_type<int>());
	detail::portion_mapply_op::const_ptr op;
	if (no->get_type() == get_scalar_type<int>())
		op = detail::portion_mapply_op::const_ptr(new ifelse_portion_op<int, int>(
					test->get_num_rows(), test->get_num_cols(), yes->get_type()));
	else if (no->get_type() == get_scalar_type<double>())
		op = detail::portion_mapply_op::const_ptr(new ifelse_portion_op<int, double>(
					test->get_num_rows(), test->get_num_cols(), yes->get_type()));
	else {
		fprintf(stderr, "unsupported type in ifelse2\n");
		return R_NilValue;
	}

	std::vector<dense_matrix::const_ptr> mats(3);
	mats[0] = test;
	mats[1] = yes;
	mats[2] = no;
	dense_matrix::ptr ret = mapply_ele(mats, op, test->store_layout());

	Rcpp::List ret_obj;
	if (ret == NULL)
		return R_NilValue;
	else if (is_vector(ptest))
		return create_FMR_vector(ret, FM_get_Rtype(pyes), "");
	else
		return create_FMR_matrix(ret, FM_get_Rtype(pyes), "");
}

/*
 * This template is a little different from the ones in matrix_ops.cpp.
 * We want to return true for both NA and NaN.
 */

template<class T, bool is_logical>
bool R_is_na(T val)
{
	fprintf(stderr, "unknown type for NA\n");
	return false;
}

template<>
bool R_is_na<int, true>(int val)
{
	return val == NA_LOGICAL;
}

template<>
bool R_is_na<int, false>(int val)
{
	return val == NA_INTEGER;
}

template<>
bool R_is_na<double, true>(double val)
{
	return ISNAN(val);
}

template<>
bool R_is_na<double, false>(double val)
{
	return ISNAN(val);
}

template<class T, int is_logical>
class isna_op: public bulk_uoperate
{
public:
	virtual void runA(size_t num_eles, const void *in_arr,
			void *out_arr) const {
		const T *in = reinterpret_cast<const T *>(in_arr);
		int *out = reinterpret_cast<int *>(out_arr);
		// is.na in R returns true for both NA and NaN.
		// we should do the same thing.
		for (size_t i = 0; i < num_eles; i++)
			out[i] = R_is_na<T, is_logical>(in[i]);
	}

	virtual const scalar_type &get_input_type() const {
		return get_scalar_type<T>();
	}

	virtual const scalar_type &get_output_type() const {
		return get_scalar_type<int>();
	}
	virtual std::string get_name() const {
		return "isna";
	}
};

// This is true only for NA.
class double_isna_only_op: public bulk_uoperate
{
public:
	virtual void runA(size_t num_eles, const void *in_arr,
			void *out_arr) const {
		int *out = reinterpret_cast<int *>(out_arr);
#ifdef RCPP_HAS_LONG_LONG_TYPES
		const rcpp_ulong_long_type *in
			= reinterpret_cast<const rcpp_ulong_long_type *>(in_arr);
		for (size_t i = 0; i < num_eles; i++)
			out[i] = FM_IsNA(in + i);
#else
		const double *in = reinterpret_cast<const double *>(in_arr);
		for (size_t i = 0; i < num_eles; i++)
			out[i] = R_IsNA(in[i]);
#endif
	}

	virtual const scalar_type &get_input_type() const {
		return get_scalar_type<double>();
	}

	virtual const scalar_type &get_output_type() const {
		return get_scalar_type<int>();
	}
	virtual std::string get_name() const {
		return "isna_only";
	}
};

RcppExport SEXP R_FM_isna(SEXP px, SEXP ponly)
{
	if (is_sparse(px)) {
		fprintf(stderr, "isna doesn't support sparse matrices\n");
		return R_NilValue;
	}
	dense_matrix::ptr x = get_matrix<dense_matrix>(px);
	R_type type = FM_get_Rtype(px);

	bool na_only = LOGICAL(ponly)[0];
	dense_matrix::ptr ret;
	if (type == R_type::R_REAL && na_only)
		ret = x->sapply(bulk_uoperate::const_ptr(new double_isna_only_op()));
	else if (type == R_type::R_REAL)
		ret = x->sapply(bulk_uoperate::const_ptr(new isna_op<double, false>()));
	else if (type == R_type::R_INT)
		ret = x->sapply(bulk_uoperate::const_ptr(new isna_op<int, false>()));
	else if (type == R_type::R_LOGICAL)
		ret = x->sapply(bulk_uoperate::const_ptr(new isna_op<int, true>()));
	if (ret == NULL)
		return R_NilValue;
	else if (is_vector(px))
		return create_FMR_vector(ret, R_type::R_LOGICAL, "");
	else
		return create_FMR_matrix(ret, R_type::R_LOGICAL, "");
}

class double_isnan_op: public bulk_uoperate
{
public:
	virtual void runA(size_t num_eles, const void *in_arr,
			void *out_arr) const {
		const double *in = reinterpret_cast<const double *>(in_arr);
		int *out = reinterpret_cast<int *>(out_arr);
		for (size_t i = 0; i < num_eles; i++)
			out[i] = R_IsNaN(in[i]);
	}

	virtual const scalar_type &get_input_type() const {
		return get_scalar_type<double>();
	}

	virtual const scalar_type &get_output_type() const {
		return get_scalar_type<int>();
	}
	virtual std::string get_name() const {
		return "isnan";
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
		return create_FMR_vector(ret, R_type::R_LOGICAL, "");
	else
		return create_FMR_matrix(ret, R_type::R_LOGICAL, "");
}

RcppExport SEXP R_FM_init(SEXP pconf)
{
	set_log_level(c_log_level::warning);
	std::string conf_file;
	if (!R_is_null(pconf) && R_is_string(pconf))
		conf_file = CHAR(STRING_ELT(pconf, 0));

	config_map::ptr configs;
	if (!conf_file.empty() && safs::file_exist(conf_file)) {
		configs = config_map::create(conf_file);
		configs->add_options("writable=1");
	}
	else if (!conf_file.empty()) {
		fprintf(stderr, "conf file %s doesn't exist.\n", conf_file.c_str());
		configs = config_map::create();
	}
	// If there isn't a conf file, we just use the default settings.
	else
		configs = config_map::create();

	bool safs_success;
	bool standalone = true;
	try {
		safs::init_io_system(configs);
		standalone = false;
		safs_success = true;
	} catch (safs::init_error &e) {
		if (!conf_file.empty())
			fprintf(stderr, "init SAFS: %s\n", e.what());
		safs_success = true;
	} catch (std::exception &e) {
		fprintf(stderr, "exception in init: %s\n", e.what());
		safs_success = false;
	}

	bool fm_success;
	try {
		fm::init_flash_matrix(configs);
		fm_success = true;
	} catch (std::exception &e) {
		fprintf(stderr, "exception in init: %s\n", e.what());
		fm_success = false;
	}

	Rcpp::LogicalVector res(1);
	res[0] = safs_success && fm_success;
	if (standalone)
		printf("Run FlashR in standalone mode\n");
	else if (safs::is_safs_init())
		printf("Run FlashR\n");
	else {
		fprintf(stderr, "Can't enable the SAFS mode of FlashR\n");
		res[0] = false;
	}
	fmr::init_udf_ext();
	fmr::init_apply_ops();
	return res;
}

RcppExport SEXP R_FM_set_conf(SEXP pconf)
{
	fm::destroy_flash_matrix();
	safs::destroy_io_system();
	return R_FM_init(pconf);
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

RcppExport SEXP R_FM_print_mat_info(SEXP pmat)
{
	if (is_sparse(pmat)) {
		sparse_matrix::ptr mat = get_matrix<sparse_matrix>(pmat);
		printf("sparse matrix of %ld rows and %ld cols\n", mat->get_num_rows(),
				mat->get_num_cols());
		printf("sparse matrix element type: %s\n",
				mat->get_type().get_name().c_str());
		block_2d_size bsize = mat->get_block_size();
		printf("sparse matrix block size: %ld, %ld\n", bsize.get_num_rows(),
				bsize.get_num_cols());
	}
	else {
		dense_matrix::ptr mat = get_matrix<dense_matrix>(pmat);
		printf("dense matrix with %ld rows and %ld cols in %s-major order\n",
				mat->get_num_rows(), mat->get_num_cols(),
				mat->store_layout() == matrix_layout_t::L_COL ? "col" : "row");
		block_matrix::ptr block_mat = std::dynamic_pointer_cast<block_matrix>(mat);
		if (block_mat)
			printf("the matrix has %ld blocks of size %ld\n",
					block_mat->get_num_blocks(), block_mat->get_block_size());
		if (!mat->is_in_mem())
			printf("dense matrix is stored on disks\n");
		else if (mat->get_data().get_num_nodes() > 0)
			printf("dense matrix is stored on %d NUMA nodes\n",
					mat->get_data().get_num_nodes());
		else
			printf("dense matrix is stored on SMP\n");
		std::string name = mat->get_data().get_name();
		printf("matrix store: %s\n", name.c_str());
	}
	return R_NilValue;
}

RcppExport SEXP R_FM_conv_store(SEXP pmat, SEXP pin_mem, SEXP pname)
{
	if (is_sparse(pmat)) {
		fprintf(stderr, "we can't convert the store of a sparse matrix\n");
		return R_NilValue;
	}

	bool in_mem = LOGICAL(pin_mem)[0];
	if (!in_mem && !safs::is_safs_init()) {
		fprintf(stderr,
				"can't convert it to ext-mem matrix when SAFS is disabled\n");
		return R_NilValue;
	}

	dense_matrix::ptr mat = get_matrix<dense_matrix>(pmat);
	std::string name = CHAR(STRING_ELT(pname, 0));
	mat = mat->conv_store(in_mem, matrix_conf.get_num_nodes());
	bool ret = mat->materialize_self();
	if (!ret) {
		fprintf(stderr, "can't materialize the matrix\n");
		return R_NilValue;
	}
	if (!name.empty() && !in_mem) {
		detail::EM_matrix_store::const_ptr store
			= detail::EM_matrix_store::cast(mat->get_raw_store());
		bool ret = store->set_persistent(name);
		if (!ret)
			return R_NilValue;
	}
	if (is_vector(pmat))
		return create_FMR_vector(mat, FM_get_Rtype(pmat), name);
	else
		return create_FMR_matrix(mat, FM_get_Rtype(pmat), name);
}

#ifdef USE_PROFILER
RcppExport SEXP R_start_profiler(SEXP pfile)
{
	std::string file = CHAR(STRING_ELT(pfile, 0));
	ProfilerStart(file.c_str());
	return R_NilValue;
}

RcppExport SEXP R_stop_profiler()
{
	ProfilerStop();
	return R_NilValue;
}
#endif

RcppExport SEXP R_FM_rand_sparse_proj(SEXP pnrow, SEXP pncol, SEXP pdensity)
{
	size_t nrow;
	size_t ncol;
	double density;
	bool ret1 = R_get_number<size_t>(pnrow, nrow);
	bool ret2 = R_get_number<size_t>(pncol, ncol);
	bool ret3 = R_get_number<double>(pdensity, density);
	if (!ret1 || !ret2 || !ret3) {
		fprintf(stderr, "the arguments aren't of the supported type\n");
		return R_NilValue;
	}

	detail::sparse_project_matrix_store::ptr store;
	matrix_layout_t layout
		= nrow > ncol ? matrix_layout_t::L_ROW : matrix_layout_t::L_COL;
	store = detail::sparse_project_matrix_store::create_sparse_rand(nrow,
			ncol, layout, get_scalar_type<double>(), density);
	return create_FMR_matrix(dense_matrix::create(store),
			trans_FM2R(store->get_type()), "");
}

RcppExport SEXP R_FM_print_features()
{
	std::string features = safs::get_supported_features();
	printf("SAFS: %s\n", features.c_str());
	features = fm::get_supported_features();
	printf("FM: %s\n", features.c_str());
	return R_NilValue;
}

RcppExport SEXP R_FM_create_factor(SEXP pmat, SEXP pnum_levels)
{
	dense_matrix::ptr mat = get_matrix<dense_matrix>(pmat);
	int num_levels = INTEGER(pnum_levels)[0];
	factor_col_vector::ptr fvec;
	if (num_levels < 0)
		fvec = factor_col_vector::create(mat);
	else
		fvec = factor_col_vector::create(factor(num_levels), mat);
	return create_FMR_vector(fvec, FM_get_Rtype(pmat), "");
}

namespace fmr
{
void set_use_na_op(bool val);
}

RcppExport SEXP R_FM_set_test_NA(SEXP pval)
{
	fmr::set_use_na_op(LOGICAL(pval)[0]);
	return R_NilValue;
}
