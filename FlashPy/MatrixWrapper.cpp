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
	// TODO support boolean and unsigned char.
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

matrix_wrapper::matrix_wrapper(intptr_t data_ptr, size_t length,
		const std::string &t)
{
	detail::mem_matrix_store::ptr store = detail::mem_matrix_store::create(
			length, 1, matrix_layout_t::L_COL, convT_py2fm(t), -1);
	memcpy(store->get_raw_arr(), (void *) data_ptr,
			length * store->get_type().get_size());
	this->mat = col_vec::create(store);
}

matrix_wrapper::matrix_wrapper(intptr_t data_ptr, size_t nrow, size_t ncol,
		const std::string &t, const std::string layout)
{
	detail::mem_matrix_store::ptr store = detail::mem_matrix_store::create(
			nrow, ncol, get_layout(layout), convT_py2fm(t), -1);
	memcpy(store->get_raw_arr(), (void *) data_ptr,
			nrow * ncol * store->get_type().get_size());
	this->mat = dense_matrix::create(store);
}

bool matrix_wrapper::is_vector() const
{
	// TODO this might be too expensive.
	auto vec = std::dynamic_pointer_cast<col_vec>(mat);
	return vec != NULL;
}

std::string matrix_wrapper::get_type_str() const
{
	if (fm2py.empty())
		init_map();
	return fm2py[mat->get_type().get_type()].name;
}

enum NPY_TYPES matrix_wrapper::get_type_py() const
{
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

matrix_wrapper matrix_wrapper::mapply2(matrix_wrapper m, bulk_op_idx_t op) const
{
	check_mat();
	m.check_mat();
	fm::dense_matrix::ptr res;
	// If the left and right one have the same shape.
	if (mat->get_num_rows() == m.mat->get_num_rows()
			&& mat->get_num_cols() == m.mat->get_num_cols())
		res = mat->mapply2(*m.mat, get_op(mat->get_type(), op));
	// The left one is a matrix.
	else if (mat->get_num_rows() > 1 && mat->get_num_cols() > 1
			// The right one is a vector.
			&& (m.mat->get_num_rows() > 1 && m.mat->get_num_cols() == 1)
			&& mat->get_num_cols() == m.mat->get_num_rows())
		res = mat->mapply_rows(get_vec(m.mat), get_op(mat->get_type(), op));
	// If the left one is a vector.
	else if ((mat->get_num_rows() > 1 && mat->get_num_cols() == 1)
			&& mat->get_num_rows() == m.mat->get_num_cols()
			// and the right one is a matrix.
			&& m.mat->get_num_rows() > 1 && m.mat->get_num_cols() > 1) {
		auto left = fm::col_vec::create(mat);
		auto left_mat = fm::dense_matrix::create_repeat(left,
				m.mat->get_num_rows(), m.mat->get_num_cols(),
				matrix_layout_t::L_ROW, true,
				m.mat->get_data().get_num_nodes());
		res = left_mat->mapply2(*m.mat, get_op(mat->get_type(), op));
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
			init_flash_matrix(NULL);
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

}
