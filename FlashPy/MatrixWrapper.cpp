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

static std::unordered_map<std::string, const fm::scalar_type *> py2fm;
// The index is prim_type and the value is the type name
static std::vector<enum NPY_TYPES> fm2py;

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

	py2fm.insert(std::pair<std::string, const fm::scalar_type *>("B",
				&fm::get_scalar_type<unsigned char>()));
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
	fm2py[fm::prim_type::P_CHAR] = NPY_BYTE;
	fm2py[fm::prim_type::P_SHORT] = NPY_SHORT;
	fm2py[fm::prim_type::P_USHORT] = NPY_USHORT;
	fm2py[fm::prim_type::P_INTEGER] = NPY_INT;
	fm2py[fm::prim_type::P_UINT] = NPY_UINT;
	fm2py[fm::prim_type::P_LONG] = NPY_LONG;
	fm2py[fm::prim_type::P_ULONG] = NPY_ULONG;
	fm2py[fm::prim_type::P_FLOAT] = NPY_FLOAT;
	fm2py[fm::prim_type::P_DOUBLE] = NPY_DOUBLE;
	fm2py[fm::prim_type::P_LDOUBLE] = NPY_LONGDOUBLE;
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
		const std::string &t, const std::string layout_str)
{
	matrix_layout_t layout;
	if (layout_str == "c")
		layout = matrix_layout_t::L_COL;
	else if (layout_str == "r")
		layout = matrix_layout_t::L_ROW;
	else
		throw std::invalid_argument("invalid type");

	detail::mem_matrix_store::ptr store = detail::mem_matrix_store::create(
			nrow, ncol, layout, convT_py2fm(t), -1);
	memcpy(store->get_raw_arr(), (void *) data_ptr,
			nrow * ncol * store->get_type().get_size());
	this->mat = dense_matrix::create(store);
}

bool matrix_wrapper::is_vector() const
{
	auto vec = std::dynamic_pointer_cast<col_vec>(mat);
	return vec != NULL;
}

std::string matrix_wrapper::get_type_str() const
{
	return mat->get_type().get_name();
}

enum NPY_TYPES matrix_wrapper::get_type_py() const
{
	if (fm2py.empty())
		init_map();
	return fm2py[mat->get_type().get_type()];
}

}
