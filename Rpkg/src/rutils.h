#ifndef __R_UTILS_H__
#define __R_UTILS_H__

/*
 * Copyright 2014 Open Connectome Project (http://openconnecto.me)
 * Written by Da Zheng (zhengda1936@gmail.com)
 *
 * This file is part of FlashR.
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

#include <memory>

#include<Rcpp.h>

bool R_is_real(SEXP v);
bool R_is_integer(SEXP v);

template<class T>
bool R_get_number(SEXP v, T &ret) {
	if (R_is_real(v)) {
		ret = REAL(v)[0];
		return true;
	}
	else if (R_is_integer(v)) {
		ret = INTEGER(v)[0];
		return true;
	}
	else
		return false;
}

/*
 * Test if this is a sparse matrix.
 */
static inline bool is_sparse(const Rcpp::List &matrix)
{
	std::string type = matrix["type"];
	return type == "sparse";
}

static inline bool is_vector(const Rcpp::List &matrix)
{
	std::string type = matrix["type"];
	return type == "vector";
}

void R_gc();

#endif
