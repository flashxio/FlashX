/*
 * Copyright 2014 Open Connectome Project (http://openconnecto.me)
 * Written by Da Zheng (zhengda1936@gmail.com)
 *
 * This file is part of FlashGraph.
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

#include "FG_sparse_matrix.h"
#include "matrix_eigensolver.h"

namespace fg
{

class AcD_SPMV: public SPMV
{
	FG_adj_matrix::ptr A;
	FG_vector<double>::ptr cd;
public:
	AcD_SPMV(FG_adj_matrix::ptr A, double c, FG_vector<vsize_t>::ptr d) {
		this->A = A;
		cd = d->multiply<double, double>(c);
	}

	virtual void compute(const FG_vector<ev_float_t> &input,
			FG_vector<ev_float_t> &output) {
		A->multiply(input, output);
		FG_vector<ev_float_t>::ptr cdx = input.multiply<double, ev_float_t>(cd);
		output.add_in_place<ev_float_t>(cdx);
	}

	virtual size_t get_vector_size() {
		return A->get_num_rows();
	}
};

void compute_AcD_uw(FG_graph::ptr fg, double c, int ncv, int nev,
		const std::string &which, std::vector<eigen_pair_t> &eigen_pairs)
{
	if (fg->get_graph_header().is_directed_graph()) {
		fprintf(stderr, "AcD doesn't support a directed graph\n");
		return;
	}
	FG_adj_matrix::ptr matrix = FG_adj_matrix::create(fg);
	AcD_SPMV spmv(matrix, c, get_degree(fg, edge_type::BOTH_EDGES));
	eigen_solver(spmv, ncv, nev, which, eigen_pairs);
}

}
