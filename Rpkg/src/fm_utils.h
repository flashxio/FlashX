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

#ifndef __FM_UTILS_H__
#define __FM_UTILS_H__

#include <Rcpp.h>
#include <memory>

namespace fm
{
	class vector;
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

	typename ObjectType::ptr get_object() {
		return o;
	}
};

template<class MatrixType>
typename MatrixType::ptr get_matrix(const Rcpp::List &matrix)
{
	object_ref<MatrixType> *ref
		= (object_ref<MatrixType> *) R_ExternalPtrAddr(matrix["pointer"]);
	return ref->get_object();
}

std::shared_ptr<fm::vector> get_vector(const Rcpp::List &vec);

SEXP create_FMR_vector(std::shared_ptr<const fm::detail::vec_store> vec, const std::string &name);
SEXP create_FMR_vector(std::shared_ptr<fm::dense_matrix> m, const std::string &name);
SEXP create_FMR_matrix(std::shared_ptr<fm::dense_matrix> m, const std::string &name);
SEXP create_FMR_matrix(std::shared_ptr<fm::sparse_matrix> m, const std::string &name);

#endif
