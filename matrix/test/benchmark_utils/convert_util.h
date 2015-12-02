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

#include <fstream>

#include "../dense_matrix.h"
#include "../vector.h"

using namespace fm;
using namespace fm::detail;

enum conv_layout{ROW, COL, RAWROW, RAWCOL};

template <typename T>
static void print_vector(typename std::vector<T>& v)
{
	std::cout << "[";

	typename std::vector<T>::iterator itr = v.begin();
	for (; itr != v.end(); itr++) {
		std::cout << " "<< *itr;
	}
	std::cout <<  " ]\n";
}

template <typename T>
static void print_mat(T* matrix, const unsigned rows, const unsigned cols,
		conv_layout lay=RAWROW) {
	for (unsigned row = 0; row < rows; row++) {
		if (lay == RAWROW) { std::cout << "["; }
		for (unsigned col = 0; col < cols; col++) {
			if (lay == RAWROW) { std::cout << " " << matrix[row*cols + col]; }
			else { std::cout << " | " << matrix[row + (rows*col)]; }
		}
		if (lay == RAWROW) { std::cout <<  " ]\n"; }
		else { std::cout << " |\n"; }
	}
}

static void print_dmat(dense_matrix::ptr dmat) {
	for (size_t row = 0; row < dmat->get_num_rows(); row++) {
		std::shared_ptr<vector> curr_row = dmat->get_row(row);
		std::vector<double> stdvec = curr_row->conv2std<double>();
		print_vector<double>(stdvec);
	}
}

// read_nrow_ncol - does the data have the nrow, ncol first before it?
// No return since outmat is malloced and populated
static double* read_fg(std::string filename, conv_layout lay,
		size_t NUM_ROWS=0, size_t NUM_COLS=0) {
	std::ifstream infile;
	infile.open(filename, std::ios::in | std::ios::binary);
	double* outmat;
	if (infile.is_open()) {
		BOOST_LOG_TRIVIAL(info) << "Beginning read ...";

		if (NUM_ROWS == 0 || NUM_COLS == 0) {
			BOOST_LOG_TRIVIAL(info) << "Reading matrix dims from file ...";
			infile.read((char*)&NUM_ROWS, sizeof(size_t));
			infile.read((char*)&NUM_COLS, sizeof(size_t));
		}

		assert(NUM_ROWS > 0 && NUM_COLS > 0);

		if (lay == RAWROW) {
		} else if (lay == RAWCOL) {
			// Swap dimensions
			size_t tmp = NUM_ROWS;
			NUM_ROWS = NUM_COLS;
			NUM_COLS = tmp;
		} else {
			assert (0);
		}

		BOOST_LOG_TRIVIAL(info) << "Number of rows = " << NUM_ROWS;
		BOOST_LOG_TRIVIAL(info) << "Number of cols = " << NUM_COLS;
		outmat = new double[NUM_ROWS*NUM_COLS];
		infile.read((char*)&outmat[0], sizeof(double)*(NUM_ROWS)*(NUM_COLS));

		infile.close();
	}
#if 0
		BOOST_LOG_TRIVIAL(info) << "Print read matrix ...";
		print_mat(outmat, NUM_ROWS, NUM_COLS);
#endif
#if 0
	BOOST_LOG_TRIVIAL(info) << "Printing read matrix as vector ...";
	std::cout << "[ ";
	for (size_t i = 0; i < NUM_ROWS*NUM_COLS; i++) {
		std::cout << outmat[i] << " ";
	}
	std::cout << "]\n";
#endif
	return outmat;
}

