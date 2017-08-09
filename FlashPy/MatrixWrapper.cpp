/*
 * Copyright 2017 Open Connectome Project (http://openconnecto.me)
 * Written by Da Zheng (zhengda1936@gmail.com)
 *
 * This file is part of FlashPy.
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

#include "mem_matrix_store.h"
#include "col_vec.h"
#include "data_frame.h"
#include "sparse_matrix.h"

#include "MatrixWrapper.h"

using namespace fm;

namespace flashpy
{

class invalid_operation: public std::exception
{
	std::string msg;
public:
	invalid_operation(std::string msg) {
		this->msg = msg;
	}

	~invalid_operation() throw() {
	}

	virtual const char *what() const throw() {
		return msg.c_str();
	}
};

struct py_type_info
{
	enum NPY_TYPES type;
	std::string name;
	py_type_info() {
		type = NPY_BYTE;
		name = "b";
	}
	py_type_info(enum NPY_TYPES type, std::string name) {
		this->type = type;
		this->name = name;
	}
};

static std::unordered_map<std::string, const fm::scalar_type *> py2fm;
// The index is prim_type and the value is the type name
static std::vector<py_type_info> fm2py;

static void init_map()
{
	// TODO support unsigned char.
	py2fm.insert(std::pair<std::string, const fm::scalar_type *>("?",
				&fm::get_scalar_type<bool>()));
	py2fm.insert(std::pair<std::string, const fm::scalar_type *>("b",
				&fm::get_scalar_type<char>()));
	py2fm.insert(std::pair<std::string, const fm::scalar_type *>("h",
				&fm::get_scalar_type<short>()));
	py2fm.insert(std::pair<std::string, const fm::scalar_type *>("i",
				&fm::get_scalar_type<int>()));
	py2fm.insert(std::pair<std::string, const fm::scalar_type *>("l",
				&fm::get_scalar_type<long>()));

#if 0
	py2fm.insert(std::pair<std::string, const fm::scalar_type *>("B",
				&fm::get_scalar_type<unsigned char>()));
#endif
	py2fm.insert(std::pair<std::string, const fm::scalar_type *>("H",
				&fm::get_scalar_type<unsigned short>()));
	py2fm.insert(std::pair<std::string, const fm::scalar_type *>("I",
				&fm::get_scalar_type<unsigned int>()));
	py2fm.insert(std::pair<std::string, const fm::scalar_type *>("L",
				&fm::get_scalar_type<unsigned long>()));

	py2fm.insert(std::pair<std::string, const fm::scalar_type *>("f",
				&fm::get_scalar_type<float>()));
	py2fm.insert(std::pair<std::string, const fm::scalar_type *>("d",
				&fm::get_scalar_type<double>()));
	py2fm.insert(std::pair<std::string, const fm::scalar_type *>("g",
				&fm::get_scalar_type<long double>()));

	fm2py.resize(fm::prim_type::NUM_TYPES);
	fm2py[fm::prim_type::P_BOOL] = py_type_info(NPY_BOOL, "?");
	fm2py[fm::prim_type::P_CHAR] = py_type_info(NPY_BYTE, "b");
	fm2py[fm::prim_type::P_SHORT] = py_type_info(NPY_SHORT, "h");
	fm2py[fm::prim_type::P_USHORT] = py_type_info(NPY_USHORT, "H");
	fm2py[fm::prim_type::P_INTEGER] = py_type_info(NPY_INT, "i");
	fm2py[fm::prim_type::P_UINT] = py_type_info(NPY_UINT, "I");
	fm2py[fm::prim_type::P_LONG] = py_type_info(NPY_LONG, "l");
	fm2py[fm::prim_type::P_ULONG] = py_type_info(NPY_ULONG, "L");
	fm2py[fm::prim_type::P_FLOAT] = py_type_info(NPY_FLOAT, "f");
	fm2py[fm::prim_type::P_DOUBLE] = py_type_info(NPY_DOUBLE, "d");
	fm2py[fm::prim_type::P_LDOUBLE] = py_type_info(NPY_LONGDOUBLE, "g");
}

const fm::scalar_type &convT_py2fm(const std::string &t)
{
	if (py2fm.empty())
		init_map();

	auto it = py2fm.find(t);
	if (it == py2fm.end())
		throw std::invalid_argument("invalid type");
	else
		return *it->second;
}

fm::scalar_variable::ptr scalar_wrapper::get_var(const fm::scalar_type *type) const
{
	if (type == NULL || *type == var->get_type())
		return var;
	else
		return var->cast_type(*type);
}

static size_t small_portion = 16 * 1024;

static int get_num_nodes(size_t nrow, size_t ncol)
{
	int num_nodes = matrix_conf.get_num_nodes();
	// When there is only one NUMA node, it's better to use SMP vector.
	if (num_nodes == 1)
		num_nodes = -1;
	// If the matrix is small, we should put it in non-NUMA matrix.
	if (nrow < small_portion && ncol < small_portion)
		num_nodes = -1;
	return num_nodes;
}

matrix_wrapper::matrix_wrapper(intptr_t data_ptr, size_t length,
		const std::string &t)
{
	detail::mem_matrix_store::ptr store = detail::mem_matrix_store::create(
			length, 1, matrix_layout_t::L_COL, convT_py2fm(t), -1);
	memcpy(store->get_raw_arr(), (void *) data_ptr,
			length * store->get_type().get_size());
	int num_nodes = get_num_nodes(length, 1);
	if (num_nodes < 0)
		this->mat = col_vec::create(store);
	else {
		auto mat = dense_matrix::create(store);
		mat = mat->conv_store(true, num_nodes);
		this->mat = col_vec::create(mat);
	}
}

matrix_wrapper::matrix_wrapper(intptr_t data_ptr, size_t nrow, size_t ncol,
		const std::string &t, const std::string layout)
{
	detail::mem_matrix_store::ptr store = detail::mem_matrix_store::create(
			nrow, ncol, get_layout(layout), convT_py2fm(t), -1);
	memcpy(store->get_raw_arr(), (void *) data_ptr,
			nrow * ncol * store->get_type().get_size());
	int num_nodes = get_num_nodes(nrow, ncol);
	if (num_nodes < 0)
		this->mat = dense_matrix::create(store);
	else {
		auto mat = dense_matrix::create(store);
		this->mat = mat->conv_store(true, num_nodes);
	}
}

matrix_wrapper::matrix_wrapper(size_t length, std::string &t)
{
	auto data = fm::dense_matrix::create(length, 1,
			fm::matrix_layout_t::L_COL, convT_py2fm(t),
			get_num_nodes(length, 1));
	mat = fm::col_vec::create(data);
}

matrix_wrapper::matrix_wrapper(size_t nrow, size_t ncol, const std::string &t,
		const std::string layout)
{
	mat = fm::dense_matrix::create(nrow, ncol, get_layout(layout),
			convT_py2fm(t), get_num_nodes(nrow, ncol));
}

bool matrix_wrapper::is_vector() const
{
	// TODO this might be too expensive.
	auto vec = std::dynamic_pointer_cast<col_vec>(mat);
	return vec != NULL;
}

std::string matrix_wrapper::get_type_str() const
{
	check_mat();
	if (fm2py.empty())
		init_map();
	return fm2py[mat->get_type().get_type()].name;
}

enum NPY_TYPES matrix_wrapper::get_type_py() const
{
	check_mat();
	if (fm2py.empty())
		init_map();
	return fm2py[mat->get_type().get_type()].type;
}

void matrix_wrapper::init_const_float(double val)
{
	check_mat();
	bool is_vec = is_vector();
	if (mat->get_type() == get_scalar_type<double>())
		mat = fm::dense_matrix::create_const<double>(val,
				mat->get_num_rows(), mat->get_num_cols(), mat->store_layout(),
				mat->get_raw_store()->get_num_nodes(), mat->is_in_mem());
	else if (mat->get_type() == get_scalar_type<float>())
		mat = fm::dense_matrix::create_const<float>(val,
				mat->get_num_rows(), mat->get_num_cols(), mat->store_layout(),
				mat->get_raw_store()->get_num_nodes(), mat->is_in_mem());
	else if (mat->get_type() == get_scalar_type<long double>())
		mat = fm::dense_matrix::create_const<long double>(val,
				mat->get_num_rows(), mat->get_num_cols(), mat->store_layout(),
				mat->get_raw_store()->get_num_nodes(), mat->is_in_mem());
	else
		throw invalid_operation("can't init as floating point");
	if (is_vec)
		mat = col_vec::create(mat);
}

void matrix_wrapper::init_const_int(long val)
{
	check_mat();
	bool is_vec = is_vector();
	if (mat->get_type() == get_scalar_type<char>())
		mat = fm::dense_matrix::create_const<char>(val,
				mat->get_num_rows(), mat->get_num_cols(), mat->store_layout(),
				mat->get_raw_store()->get_num_nodes(), mat->is_in_mem());
	else if (mat->get_type() == get_scalar_type<short>())
		mat = fm::dense_matrix::create_const<short>(val,
				mat->get_num_rows(), mat->get_num_cols(), mat->store_layout(),
				mat->get_raw_store()->get_num_nodes(), mat->is_in_mem());
	else if (mat->get_type() == get_scalar_type<int>())
		mat = fm::dense_matrix::create_const<int>(val,
				mat->get_num_rows(), mat->get_num_cols(), mat->store_layout(),
				mat->get_raw_store()->get_num_nodes(), mat->is_in_mem());
	else if (mat->get_type() == get_scalar_type<long>())
		mat = fm::dense_matrix::create_const<long>(val,
				mat->get_num_rows(), mat->get_num_cols(), mat->store_layout(),
				mat->get_raw_store()->get_num_nodes(), mat->is_in_mem());
	else if (mat->get_type() == get_scalar_type<unsigned short>())
		mat = fm::dense_matrix::create_const<unsigned short>(val,
				mat->get_num_rows(), mat->get_num_cols(), mat->store_layout(),
				mat->get_raw_store()->get_num_nodes(), mat->is_in_mem());
	else if (mat->get_type() == get_scalar_type<unsigned int>())
		mat = fm::dense_matrix::create_const<unsigned int>(val,
				mat->get_num_rows(), mat->get_num_cols(), mat->store_layout(),
				mat->get_raw_store()->get_num_nodes(), mat->is_in_mem());
	else if (mat->get_type() == get_scalar_type<unsigned long>())
		mat = fm::dense_matrix::create_const<unsigned long>(val,
				mat->get_num_rows(), mat->get_num_cols(), mat->store_layout(),
				mat->get_raw_store()->get_num_nodes(), mat->is_in_mem());
	else
		throw invalid_operation("can't init as integer");
	if (is_vec)
		mat = col_vec::create(mat);
}

matrix_wrapper matrix_wrapper::cast_ele_type(std::string dtype) const {
	check_mat();
	auto res = mat->cast_ele_type(convT_py2fm(dtype));
	if (is_vector())
		return matrix_wrapper(fm::col_vec::create(res));
	else
		return matrix_wrapper(res);
}

matrix_wrapper matrix_wrapper::as_factor(int num_levels) const
{
	check_mat();
	if (mat->get_num_rows() > 1 && mat->get_num_cols() > 1)
		throw invalid_operation("can't cast a matrix to a vector.");
	if (num_levels < 0)
		return matrix_wrapper(fm::factor_col_vector::create(mat));
	else
		return matrix_wrapper(fm::factor_col_vector::create(fm::factor(num_levels), mat));
}

matrix_wrapper matrix_wrapper::as_vector() const
{
	check_mat();
	if (mat->get_num_rows() > 1 && mat->get_num_cols() > 1)
		throw invalid_operation("can't cast a matrix to a vector.");
	return matrix_wrapper(fm::col_vec::create(mat));
}

matrix_wrapper matrix_wrapper::as_matrix() const
{
	check_mat();
	if (is_vector())
		return matrix_wrapper(fm::dense_matrix::create(mat->get_raw_store()));
	else
		return *this;
}

matrix_wrapper matrix_wrapper::inner_prod(matrix_wrapper m, bulk_op_idx_t left_op,
		bulk_op_idx_t right_op) const
{
	check_mat();
	m.check_mat();
	if (mat->get_type() == m.mat->get_type())
		return matrix_wrapper(mat->inner_prod(*m.mat,
					get_op(mat->get_type(), left_op),
					get_op(mat->get_type(), right_op)));
	else {
		const scalar_type &common_type = get_larger_type(mat->get_type(),
				m.mat->get_type());
		dense_matrix::ptr left = mat->cast_ele_type(common_type);
		dense_matrix::ptr right = m.mat->cast_ele_type(common_type);
		return matrix_wrapper(left->inner_prod(*right, get_op(common_type, left_op),
					get_op(common_type, right_op)));
	}
}

matrix_wrapper matrix_wrapper::multiply(matrix_wrapper m) const
{
	check_mat();
	m.check_mat();
	auto res = mat->multiply(*m.mat);
	if (m.is_vector())
		return matrix_wrapper(fm::col_vec::create(res));
	else
		return matrix_wrapper(res);
}

matrix_wrapper matrix_wrapper::mapply_cols(matrix_wrapper vals,
		bulk_op_idx_t op) const
{
	check_mat();
	vals.check_mat();
	if (mat->get_type() == vals.mat->get_type())
		return matrix_wrapper(mat->mapply_cols(get_vec(vals.mat),
					get_op(mat->get_type(), op)));
	else {
		const scalar_type &common_type = get_larger_type(mat->get_type(),
				vals.mat->get_type());
		dense_matrix::ptr left = mat->cast_ele_type(common_type);
		col_vec::ptr right = col_vec::create(vals.mat->cast_ele_type(common_type));
		return matrix_wrapper(left->mapply_cols(right, get_op(common_type, op)));
	}
}

matrix_wrapper matrix_wrapper::mapply_rows(matrix_wrapper vals,
		bulk_op_idx_t op) const
{
	check_mat();
	vals.check_mat();
	if (mat->get_type() == vals.mat->get_type())
		return matrix_wrapper(mat->mapply_rows(get_vec(vals.mat),
					get_op(mat->get_type(), op)));
	else {
		const scalar_type &common_type = get_larger_type(mat->get_type(),
				vals.mat->get_type());
		dense_matrix::ptr left = mat->cast_ele_type(common_type);
		col_vec::ptr right = col_vec::create(vals.mat->cast_ele_type(common_type));
		return matrix_wrapper(left->mapply_rows(right, get_op(common_type, op)));
	}
}

matrix_wrapper matrix_wrapper::apply_scalar(scalar_wrapper var,
		bulk_op_idx_t op) const
{
	check_mat();
	dense_matrix::ptr res;
	if (mat->get_type() == var.get_type())
		res = mat->apply_scalar(var.get_var(), get_op(mat->get_type(), op));
	else {
		const scalar_type &common_type = get_larger_type(mat->get_type(),
				var.get_type());
		dense_matrix::ptr left = mat->cast_ele_type(common_type);
		res = left->apply_scalar(var.get_var(&common_type),
				get_op(left->get_type(), op));
	}
	if (is_vector())
		return matrix_wrapper(fm::col_vec::create(res));
	else
		return matrix_wrapper(res);
}

matrix_wrapper matrix_wrapper::mapply2(matrix_wrapper m, bulk_op_idx_t op) const
{
	check_mat();
	m.check_mat();
	fm::dense_matrix::ptr res;
	dense_matrix::ptr left, right;
	if (mat->get_type() == m.mat->get_type()) {
		left = mat;
		right = m.mat;
	}
	else {
		const scalar_type &common_type = get_larger_type(mat->get_type(),
				m.mat->get_type());
		left = mat->cast_ele_type(common_type);
		right = m.mat->cast_ele_type(common_type);
	}

	// If the left and right one have the same shape.
	if (left->get_num_rows() == right->get_num_rows()
			&& left->get_num_cols() == right->get_num_cols())
		res = left->mapply2(*right, get_op(left->get_type(), op));
	else if (right->get_num_rows() == 1 && right->get_num_cols() == 1) {
		auto val = right->get(0, 0);
		res = left->apply_scalar(val, get_op(left->get_type(), op));
	}
	else if (left->get_num_rows() == 1 && left->get_num_cols() == 1) {
		auto val = left->get(0, 0);
		left = dense_matrix::create_const(val, right->get_num_rows(),
				right->get_num_cols(), right->store_layout(),
				right->get_data().get_num_nodes(), right->is_in_mem());
		res = left->mapply2(*right, get_op(left->get_type(), op));
	}
	// The left one is a matrix.
	else if (left->get_num_rows() > 1 && left->get_num_cols() > 1
			// The right one is a vector.
			&& m.is_vector()
			&& left->get_num_cols() == right->get_num_rows())
		res = left->mapply_rows(get_vec(right), get_op(left->get_type(), op));
	// If the left one is a vector.
	else if (is_vector() && left->get_num_rows() == right->get_num_cols()
			// and the right one is a matrix.
			&& right->get_num_rows() > 1 && right->get_num_cols() > 1) {
		auto left_vec = fm::col_vec::create(left);
		auto left_mat = fm::dense_matrix::create_repeat(left_vec,
				right->get_num_rows(), right->get_num_cols(),
				matrix_layout_t::L_ROW, true,
				right->get_data().get_num_nodes());
		res = left_mat->mapply2(*right, get_op(left->get_type(), op));
	}
	else if (left->get_num_rows() > 1 && left->get_num_cols() > 1
			// The right one is a one-col matrix.
			&& (right->get_num_rows() > 1 && right->get_num_cols() == 1)
			&& left->get_num_rows() == right->get_num_rows()) {
		auto right_vec = fm::col_vec::create(right);
		res = left->mapply_cols(right_vec, get_op(left->get_type(), op));
	}
	else if (left->get_num_rows() > 1 && left->get_num_cols() > 1
			// The right one is a one-row matrix.
			&& (right->get_num_rows() == 1 && right->get_num_cols() > 1)
			&& left->get_num_cols() == right->get_num_cols()) {
		auto right_vec = fm::col_vec::create(right);
		res = left->mapply_rows(right_vec, get_op(left->get_type(), op));
	}
	// If the left one is a one-col matrix.
	else if ((left->get_num_rows() > 1 && left->get_num_cols() == 1)
			&& left->get_num_rows() == right->get_num_rows()
			// and the right one is a matrix.
			&& right->get_num_rows() > 1 && right->get_num_cols() > 1) {
		auto left_vec = fm::col_vec::create(left);
		auto left_mat = fm::dense_matrix::create_repeat(left_vec,
				right->get_num_rows(), right->get_num_cols(),
				matrix_layout_t::L_COL, false,
				right->get_data().get_num_nodes());
		res = left_mat->mapply2(*right, get_op(left->get_type(), op));
	}
	// If the left one is a one-row matrix.
	else if ((left->get_num_rows() == 1 && left->get_num_cols() > 1)
			&& left->get_num_cols() == right->get_num_cols()
			// and the right one is a matrix.
			&& right->get_num_rows() > 1 && right->get_num_cols() > 1) {
		auto left_vec = fm::col_vec::create(left);
		auto left_mat = fm::dense_matrix::create_repeat(left_vec,
				right->get_num_rows(), right->get_num_cols(),
				matrix_layout_t::L_ROW, true,
				right->get_data().get_num_nodes());
		res = left_mat->mapply2(*right, get_op(left->get_type(), op));
	}
	else {
		throw std::invalid_argument(
				"The shape of the two matrices doesn't match");
	}

	if (m.is_vector() && is_vector())
		return matrix_wrapper(fm::col_vec::create(res));
	else
		return matrix_wrapper(res);
}

bool matrix_wrapper::copy_rows_to(char *arr, size_t len) const
{
	check_mat();
	dense_matrix::ptr data = mat;
	if (get_num_rows() > 1 && get_num_cols() > 1)
		data = mat->conv2(matrix_layout_t::L_ROW);
	data->move_store(true, -1);
	auto store = std::dynamic_pointer_cast<const fm::detail::mem_matrix_store>(
			data->get_raw_store());
	if (get_num_rows() * get_num_cols() * get_entry_size() > len) {
		fprintf(stderr, "the array is too small for the matrix\n");
		return false;
	}

	const char *src = store->get_raw_arr();
	if (src)
		memcpy(arr, src, get_num_rows() * get_num_cols() * get_entry_size());
	else {
		auto row_store = std::dynamic_pointer_cast<const fm::detail::mem_row_matrix_store>(
				data->get_raw_store());
		assert(store);
		for (size_t i = 0; i < get_num_rows(); i++)
			memcpy(arr + i * get_num_cols() * get_entry_size(),
					row_store->get_row(i), get_num_cols() * get_entry_size());
	}
	return true;
}

void matrix_wrapper::set_cached(bool cached)
{
	materialize_level level
		= cached ? materialize_level::MATER_FULL : materialize_level::MATER_CPU;
	mat->set_materialize_level(level);
}

bool init_flashpy_c(const std::string &conf_file)
{
	try {
		if (conf_file.empty())
			init_flash_matrix(config_map::create());
		else {
			config_map::ptr configs = config_map::create(conf_file);
			configs->add_options("writable=1");
			fm::init_flash_matrix(configs);
		}
		return true;
	} catch (std::exception &e) {
		fprintf(stderr, "exception in init: %s\n", e.what());
		return false;
	}
}

matrix_wrapper matrix_wrapper::ifelse(matrix_wrapper x, matrix_wrapper y) const
{
	check_mat();
	dense_matrix::ptr left, right;
	if (x.mat->get_type() == y.mat->get_type()) {
		left = x.mat;
		right = y.mat;
	}
	else {
		const scalar_type &common_type = get_larger_type(x.mat->get_type(),
				y.mat->get_type());
		left = x.mat->cast_ele_type(common_type);
		right = y.mat->cast_ele_type(common_type);
	}
	auto ret = mat->ifelse(*left, *right);
	if (is_vector())
		return matrix_wrapper(col_vec::create(ret));
	else
		return matrix_wrapper(ret);
}

static inline std::shared_ptr<fm::factor_col_vector> get_factor(
		fm::dense_matrix::ptr mat)
{
	auto ret = std::dynamic_pointer_cast<fm::factor_col_vector>(mat);
	if (ret == NULL)
		throw std::invalid_argument("invalid matrix, want a factor col vector");
	else
		return ret;
}

matrix_wrapper matrix_wrapper::groupby_row(matrix_wrapper labels,
		agg_op_idx_t op) const
{
	check_mat();
	labels.check_mat();
	return matrix_wrapper(mat->groupby_row(get_factor(labels.mat),
				get_agg(mat->get_type(), op)));
}

std::pair<matrix_wrapper, matrix_wrapper> matrix_wrapper::groupby(
		agg_op_idx_t op, bool with_val) const
{
	check_mat();
	auto ret = mat->groupby(get_agg(mat->get_type(), op), with_val);
	if (with_val) {
		matrix_wrapper agg(col_vec::create(vector::create(ret->get_vec("agg"))));
		matrix_wrapper val(col_vec::create(vector::create(ret->get_vec("val"))));
		return std::pair<matrix_wrapper, matrix_wrapper>(agg, val);
	}
	else {
		matrix_wrapper agg(col_vec::create(vector::create(ret->get_vec("agg"))));
		return std::pair<matrix_wrapper, matrix_wrapper>(agg, matrix_wrapper());
	}
}

matrix_wrapper matrix_wrapper::get_col(long idx) const
{
	check_mat();
	auto vec = mat->get_col(idx);
	if (vec == NULL)
		throw std::invalid_argument("can't get a col");
	return matrix_wrapper(col_vec::create(vec));
}

matrix_wrapper matrix_wrapper::get_row(long idx) const
{
	check_mat();
	auto vec = mat->get_row(idx);
	if (vec == NULL)
		throw std::invalid_argument("can't get a row");
	return matrix_wrapper(col_vec::create(vec));
}

matrix_wrapper matrix_wrapper::get_cols(const std::vector<off_t> &idxs) const
{
	check_mat();
	auto ret = mat->get_cols(idxs);
	if (ret == NULL)
		throw std::invalid_argument("can't get cols");
	if (is_vector())
		return matrix_wrapper(col_vec::create(ret));
	else
		return matrix_wrapper(ret);
}

matrix_wrapper matrix_wrapper::get_rows(const std::vector<off_t> &idxs) const
{
	check_mat();
	auto ret = mat->get_rows(idxs);
	if (ret == NULL)
		throw std::invalid_argument("can't get rows");
	if (is_vector())
		return matrix_wrapper(col_vec::create(ret));
	else
		return matrix_wrapper(ret);
}

matrix_wrapper matrix_wrapper::get_cols(matrix_wrapper idxs) const
{
	check_mat();
	idxs.check_mat();
	auto ret = mat->get_cols(get_vec(idxs.mat));
	if (ret == NULL)
		throw std::invalid_argument("can't get cols");
	if (is_vector())
		return matrix_wrapper(col_vec::create(ret));
	else
		return matrix_wrapper(ret);
}

matrix_wrapper matrix_wrapper::get_rows(matrix_wrapper idxs) const
{
	check_mat();
	idxs.check_mat();
	auto ret = mat->get_rows(get_vec(idxs.mat));
	if (ret == NULL)
		throw std::invalid_argument("can't get rows");
	if (is_vector())
		return matrix_wrapper(col_vec::create(ret));
	else
		return matrix_wrapper(ret);
}

matrix_wrapper matrix_wrapper::get_cols(size_t start, size_t end, long step) const
{
	check_mat();
	auto ret = mat->get_cols(start, end, step);
	if (ret == NULL)
		throw std::invalid_argument("can't get cols");
	if (is_vector())
		return matrix_wrapper(col_vec::create(ret));
	else
		return matrix_wrapper(ret);
}

matrix_wrapper matrix_wrapper::get_rows(size_t start, size_t end, long step) const
{
	check_mat();
	auto ret = mat->get_rows(start, end, step);
	if (ret == NULL)
		throw std::invalid_argument("can't get rows");
	if (is_vector())
		return matrix_wrapper(col_vec::create(ret));
	else
		return matrix_wrapper(ret);
}

matrix_wrapper matrix_wrapper::get_eles(matrix_wrapper idxs) const
{
	check_mat();
	auto ret = mat->get_eles(idxs.mat);
	if (ret == NULL)
		throw std::invalid_argument("can't get eles");
	return matrix_wrapper(ret);
}

matrix_wrapper matrix_wrapper::set_cols(const std::vector<off_t> &idxs,
		matrix_wrapper cols)
{
	check_mat();
	auto ret = mat->set_cols(idxs, cols.mat);
	if (ret == NULL)
		throw std::invalid_argument("can't set cols");
	if (is_vector())
		return matrix_wrapper(col_vec::create(ret));
	else
		return matrix_wrapper(ret);
}

matrix_wrapper matrix_wrapper::set_rows(const std::vector<off_t> &idxs,
		matrix_wrapper rows)
{
	check_mat();
	auto ret = mat->set_rows(idxs, rows.mat);
	if (ret == NULL)
		throw std::invalid_argument("can't set rows");
	if (is_vector())
		return matrix_wrapper(col_vec::create(ret));
	else
		return matrix_wrapper(ret);
}

matrix_wrapper matrix_wrapper::cbind(const std::vector<matrix_wrapper> &mats)
{
	std::vector<dense_matrix::ptr> tmp(mats.size());
	for (size_t i = 0; i < tmp.size(); i++)
		tmp[i] = mats[i].mat;
	auto ret = dense_matrix::cbind(tmp);
	if (ret == NULL)
		throw std::invalid_argument("can't cbind the matrices");
	return matrix_wrapper(ret);
}

matrix_wrapper matrix_wrapper::rbind(const std::vector<matrix_wrapper> &mats)
{
	std::vector<dense_matrix::ptr> tmp(mats.size());
	bool is_vec = mats[0].is_vector();
	for (size_t i = 0; i < tmp.size(); i++) {
		tmp[i] = mats[i].mat;
		is_vec = is_vec && mats[i].is_vector();
	}
	dense_matrix::ptr ret = dense_matrix::rbind(tmp);
	if (ret == NULL)
		throw std::invalid_argument("can't rbind the matrices");
	if (is_vec)
		return matrix_wrapper(col_vec::create(ret));
	else
		return matrix_wrapper(ret);
}

matrix_wrapper matrix_wrapper::cum_row(agg_op_idx_t op) const
{
	check_mat();
	auto ret = mat->cum(fm::matrix_margin::MAR_ROW,
			get_agg(mat->get_type(), op));
	if (is_vector())
		return matrix_wrapper(col_vec::create(ret));
	else
		return matrix_wrapper(ret);
}

matrix_wrapper matrix_wrapper::cum_col(agg_op_idx_t op) const
{
	check_mat();
	auto ret = mat->cum(fm::matrix_margin::MAR_COL,
			get_agg(mat->get_type(), op));
	if (is_vector())
		return matrix_wrapper(col_vec::create(ret));
	else
		return matrix_wrapper(ret);
}

}
