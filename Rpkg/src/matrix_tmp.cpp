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

/*
 * This file contains the code that should be replaced by R apply functions.
 * So all of the implementations here are temporary.
 */

#include "rutils.h"
#include "fm_utils.h"
#include "dense_matrix.h"
#include "bulk_operate.h"

using namespace fm;

template<class T>
class norm2: public apply_operate
{
public:
	// norm always outputs only one double float-point.
	norm2(size_t num_in_eles): apply_operate(num_in_eles) {
	}

	virtual void run(const void *input, size_t num_in_eles, void *output) const {
		const T *t_in = (const T *) input;
		double *d_out = (double *) output;

		T sum = 0;
		for (size_t i = 0; i < num_in_eles; i++) {
			T v = t_in[i];
			sum += v * v;
		}
		if (sum == 0) {
			for (size_t i = 0; i < num_in_eles; i++)
				d_out[i] = 0;
		}
		else {
			double norm = std::sqrt(sum);
			for (size_t i = 0; i < num_in_eles; i++)
				d_out[i] = t_in[i] / norm;
		}
	}

	virtual const scalar_type &get_input_type() const {
		static scalar_type_impl<T> t;
		return t;
	}

	virtual const scalar_type &get_output_type() const {
		static scalar_type_impl<double> t;
		return t;
	}
};

RcppExport SEXP R_FM_norm_matrix(SEXP pmat, SEXP pmargin, SEXP ptype)
{
	if (is_sparse(pmat)) {
		fprintf(stderr, "can't compute norm on a sparse matrix\n");
		return R_NilValue;
	}

	apply_margin margin = (apply_margin) INTEGER(pmargin)[0];
	if (margin != apply_margin::MAR_ROW && margin != apply_margin::MAR_COL) {
		fprintf(stderr, "unsupported margin\n");
		return R_NilValue;
	}
	std::string type = CHAR(STRING_ELT(ptype, 0));
	if (type != "2") {
		fprintf(stderr, "only norm2 is supported\n");
		return R_NilValue;
	}

	dense_matrix::ptr mat = get_matrix<dense_matrix>(pmat);
	dense_matrix::ptr ret;
	size_t num_in_eles;
	if (margin == apply_margin::MAR_ROW)
		num_in_eles = mat->get_num_cols();
	else
		num_in_eles = mat->get_num_rows();

	if (mat->is_type<double>())
		ret = mat->apply(margin, norm2<double>(num_in_eles));
	else if (mat->is_type<int>())
		ret = mat->apply(margin, norm2<int>(num_in_eles));
	else {
		fprintf(stderr, "unsupported data type\n");
		return R_NilValue;
	}

	if (is_vector(pmat))
		return create_FMR_vector(ret, "");
	else
		return create_FMR_matrix(ret, "");
}
