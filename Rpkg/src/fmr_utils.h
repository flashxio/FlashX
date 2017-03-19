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

#ifndef __FMR_UTILS_H__
#define __FMR_UTILS_H__

#include <Rcpp.h>
#include <memory>

#include "rutils.h"

namespace fm
{
	class col_vec;
	class factor_col_vector;
	class dense_matrix;
	class sparse_matrix;

namespace detail
{
	class vec_store;
}
}

template<class ObjectType>
class object_ref
{
	typename ObjectType::ptr o;
public:
	object_ref(typename ObjectType::ptr o) {
		this->o = o;
	}

	typename ObjectType::ptr get_object() const {
		return o;
	}

	void set_object(typename ObjectType::ptr obj) {
		this->o = obj;
	}
};

template<class MatrixType>
typename MatrixType::ptr get_matrix(const Rcpp::S4 &matrix)
{
	// TODO I should test if the pointer slot does exist.
	object_ref<MatrixType> *ref
		= (object_ref<MatrixType> *) R_ExternalPtrAddr(matrix.slot("pointer"));
	return ref->get_object();
}
template<class MatrixType>
void set_matrix(const Rcpp::S4 &matrix, typename MatrixType::ptr mat)
{
	// TODO I should test if the pointer slot does exist.
	object_ref<MatrixType> *ref
		= (object_ref<MatrixType> *) R_ExternalPtrAddr(matrix.slot("pointer"));
	ref->set_object(mat);
}

std::shared_ptr<fm::col_vec> get_vector(const Rcpp::S4 &vec);
std::shared_ptr<fm::factor_col_vector> get_factor_vector(const Rcpp::S4 &vec);

SEXP create_FMR_vector(std::shared_ptr<const fm::detail::vec_store> vec,
		R_type type, const std::string &name);
SEXP create_FMR_vector(std::shared_ptr<fm::dense_matrix> m,
		R_type type, const std::string &name);
SEXP create_FMR_vector(std::shared_ptr<fm::factor_col_vector> v,
		R_type type, const std::string &name);
SEXP create_FMR_matrix(std::shared_ptr<fm::dense_matrix> m,
		R_type type, const std::string &name);
SEXP create_FMR_matrix(std::shared_ptr<fm::sparse_matrix> m,
		R_type type, const std::string &name);
SEXP create_FMR_data_frame(std::shared_ptr<fm::data_frame> df,
		const std::vector<R_type> &type, const std::string &name);

#endif
