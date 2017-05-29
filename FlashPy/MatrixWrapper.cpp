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
}

void matrix_wrapper::init_const_int(long val)
{
	check_mat();
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
	return matrix_wrapper(fm::dense_matrix::create(mat->get_raw_store()));
}

}
