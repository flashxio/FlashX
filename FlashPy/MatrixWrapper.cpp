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
