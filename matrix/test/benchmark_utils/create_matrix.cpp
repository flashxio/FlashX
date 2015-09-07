/*
 * Copyright 2014 Open Connectome Project (http://openconnecto.me)
 * Written by Disa Mhembere (disa@jhu.edu)
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

// Create a test matrix to check the conversion to other formats for testing
//	kmeans

#include <fstream>

#include <math.h>
#include "log.h"

#include "../generic_type.h"
#include "../mem_matrix_store.h"
#include "convert_util.h"

using namespace fm;
using namespace fm::detail;

class set_col_operate: public type_set_operate<double>
{
	size_t num_cols;
	public:
	set_col_operate(size_t num_cols) {
		this->num_cols = num_cols;
	}

	void set(double *arr, size_t num_eles, off_t row_idx, off_t col_idx) const {
		for (size_t i = 0; i < num_eles; i++) {
			arr[i] = (row_idx + i) * num_cols + col_idx;
		}
	}
};

int main (int argc, char* argv[]) {
	if (argc < 2) {
		fprintf(stderr, "usage: ./create_matrix [row/col]");
		exit(911);
	}

	const size_t nrow = 5;
	const size_t ncol = 7;
	std::string outfn = "data/tiny/testmat";

	if (std::string(argv[1]) == "row") {
		const scalar_variable_impl<double> mean(.5);
		const set_col_operate sco(ncol);

		// L_COL because I expect to get this layout form eigensolver
		dense_matrix::ptr dmat = dense_matrix::create(nrow, ncol, matrix_layout_t::L_COL, 
				get_scalar_type<double>(), sco, -1, true);

		BOOST_LOG_TRIVIAL(info) << "TRANSPOSING matrix to row major ...";
		dmat = dmat->transpose();

#if 1
		BOOST_LOG_TRIVIAL(info) << "Printing row-wise matrix for verification!\n";
		print_dmat(dmat);
#endif

		BOOST_LOG_TRIVIAL(info) << "Writing the matrix";
		// convert to mem_matrix_store
		mem_matrix_store::const_ptr mms = mem_matrix_store::cast(dmat->get_raw_store());
		mms->write2file(outfn+"_dm_rw.bin");

#if 1
		BOOST_LOG_TRIVIAL(info) << "Test read of matrix";
		dense_matrix::ptr rdmat = dense_matrix::create(mem_matrix_store::load(outfn+"_dm_rw.bin"));
		BOOST_LOG_TRIVIAL(info) << "Test print of matrix";
		print_dmat(rdmat);
#endif
	} else if (std::string(argv[1]) == "col") {
		int min = 1; int max = 50;

		double* dmat = new double [nrow*ncol];
		for (size_t i = 0; i < nrow*ncol; i++) {
			dmat[i] = ceil(min + ((double)random() / (double)RAND_MAX * (max - min)));
		}

#if 1
		BOOST_LOG_TRIVIAL(info) << "Printing col-wise matrix for verification";
		for (size_t row = 0; row < nrow; row++) {
			for (size_t col = 0; col < ncol; col++) {
				std::cout << " | " << dmat[row + (nrow*col)];
			}
			std::cout << " |\n";
		}

#endif
#if 0
		BOOST_LOG_TRIVIAL(info) << "Printing col-wise matrix as vector for verification";
		std::cout << "[ ";
		for (size_t i = 0; i < nrow*ncol; i++) {
			std::cout << dmat[i] << " ";	
		}
		std::cout << "]\n";
#endif

		BOOST_LOG_TRIVIAL(info) << "Writing the matrix";
		std::ofstream outfile;
		outfile.open(outfn+"_cw.bin", std::ios::binary | std::ios::trunc | std::ios::out);
		outfile.write((char*)&dmat[0], (sizeof(double)*nrow*ncol));
		outfile.close();
		delete [] dmat;

		BOOST_LOG_TRIVIAL(info) << "Test read of matrix";
	    read_fg(outfn+"_cw.bin", L_COL, nrow, ncol); 
	} else {
		fprintf(stderr, "Unknown matrix type '%s'", argv[1]);
	}

	return EXIT_SUCCESS;
}
