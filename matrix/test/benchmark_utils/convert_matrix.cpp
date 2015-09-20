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

#include "log.h"

#include "../mem_matrix_store.h"
#include "convert_util.h"

using namespace fm;
using namespace fm::detail;

#if 0
size_t g_nrow = 5; // FIXME
size_t g_ncol = 7; // FIXME

size_t g_nrow = 8; // FIXME
size_t g_ncol = 65608366; // FIXME
#endif

size_t g_nrow, g_ncol;

static dense_matrix::ptr get_mat(std::string filename) {
	mem_matrix_store::ptr mat = mem_matrix_store::load(filename);
	BOOST_LOG_TRIVIAL(info) << "Loaded matrix into mem";

	dense_matrix::ptr dmat = dense_matrix::create(mat);
	BOOST_LOG_TRIVIAL(info) << "Converted matrix into dense";
	BOOST_LOG_TRIVIAL(info) << "The matrix layout is (0=L_COL, 1=L_ROW): " << dmat->store_layout();

	if (dmat->store_layout() != L_ROW) {
		BOOST_LOG_TRIVIAL(info) << "Converting matrix to row major ...";
		dmat = dmat->transpose();
		BOOST_LOG_TRIVIAL(info) << "The matrix layout is now: " << dmat->store_layout();
		BOOST_LOG_TRIVIAL(info) << "Dim (" << dmat->get_num_rows() << ", " << dmat->get_num_cols()
			<< "), with: " << dmat->get_entry_size() << " entries";
	}
	return dmat;
}

// COMMON for reading contiguous double data
static double* bin_read(std::string fname) {
	std::ifstream _if;
	_if.open(fname, std::ios::binary | std::ios::in);
	double* outmat = new double[g_nrow*g_ncol];


	BOOST_LOG_TRIVIAL(info) << "Binary read";
	_if.read((char*)&outmat[0], sizeof(double)*g_nrow*g_ncol);

	_if.close();
	BOOST_LOG_TRIVIAL(info) << "Read completed!";
	return outmat;
}

// SPARK
static void to_spark(dense_matrix::ptr dmat, std::ofstream& of) {
	for (size_t row=0; row < dmat->get_num_rows(); row++) {
		std::shared_ptr<vector> curr_row = dmat->get_row(row);
		std::vector<double> stdvec = curr_row->conv2std<double>();

		for (size_t col=0; col < stdvec.size(); col++) {
			if (col < stdvec.size()-1) {
				of << stdvec[col] << " ";
			} else {
				of << stdvec[col];
			}
		}
		of << "\n";
	}
}

static void to_spark(double* outmat, std::ofstream& of, const conv_layout lay) {
	BOOST_LOG_TRIVIAL(info) << "Writing matrix";

	if (lay == RAWCOL) {
		// Write it
		for (size_t row=0; row < g_ncol; row++) {
			for (size_t col=0; col < g_nrow; col++) {
				if (col < g_nrow-1) {
					of << outmat[row*g_nrow+col] << " ";
				} else {
					of << outmat[row*g_nrow+col];
				}
			}
			of << "\n";
		}
	} else if (lay == RAWROW) {
		for (size_t row = 0; row < g_nrow; row++) {
			for (size_t col = 0; col < g_ncol; col++) {
				if (col < g_ncol-1) {
					of << outmat[row*g_ncol+col] << " ";
				} else {
					of << outmat[row*g_ncol+col];
				}
			}
			of << "\n";
		}
	} else { assert(0); }
	delete [] outmat;
}

// KMEANS_PAR
static void to_kmeans_par(dense_matrix::ptr dmat, std::ofstream& of) {
	for (size_t row=0; row < dmat->get_num_rows(); row++) {
		std::shared_ptr<vector> curr_row = dmat->get_row(row);
		std::vector<double> stdvec = curr_row->conv2std<double>();

		of << row+1 << " ";
		for (size_t col=0; col < stdvec.size(); col++) {
			if (col < stdvec.size()-1) {
				of << stdvec[col] << " ";
			} else {
				of << stdvec[col];
			}
		}
		of << "\n";
	}
}

static void to_kmeans_par(double* outmat, std::ofstream& of, const conv_layout lay) {
	BOOST_LOG_TRIVIAL(info) << "Writing matrix";
	// Write it
	if (lay == RAWCOL) {
		for (size_t row = 0; row < g_ncol; row++) {
			of << row+1 << " ";
			for (size_t col = 0; col < g_nrow; col++) {
				if (col < g_nrow-1) {
					of << outmat[row*g_nrow+col] << " ";
				} else {
					of << outmat[row*g_nrow+col];
				}
			}
			of << "\n";
		}
	} else if (lay == RAWROW) {
		for (size_t row = 0; row < g_nrow; row++) {
			of << row+1 << " ";
			for (size_t col = 0; col < g_ncol; col++) {
				if (col < g_ncol-1) {
					of << outmat[row*g_ncol+col] << " ";
				} else {
					of << outmat[row*g_ncol+col];
				}
			}
			of << "\n";
		}
	} else { assert(0); }
	delete [] outmat;
}

// Write to disk size_t bytes with nrow, size_t bytes with ncol, then all data row-wise
// dmat - the dense matrix to be converted
// of - the output filestream with which to write the converted matrix
static void to_fg(dense_matrix::ptr dmat, std::ofstream& of) {
	const size_t NUM_ROWS = dmat->get_num_rows();
	const size_t NUM_COLS = dmat->get_row(0)->get_length();
	BOOST_LOG_TRIVIAL(info) << "nrow = " << NUM_ROWS << ", ncol = " << NUM_COLS;

	of.write((char*)&NUM_ROWS, sizeof(size_t)); // size_t rows
	of.write((char*)&NUM_COLS, sizeof(size_t)) ; // size_t cols

	for (size_t row=0; row < dmat->get_num_rows(); row++) {
		std::shared_ptr<vector> curr_row = dmat->get_row(row);
		std::vector<double> stdvec = curr_row->conv2std<double>();

		of.write((char*)(&stdvec[0]), sizeof(stdvec[0])*NUM_COLS);
	}
}

static void to_fg(double* mat, std::ofstream& of, const conv_layout lay) {
	const size_t NUM_ROWS = g_nrow;
	const size_t NUM_COLS = g_ncol;
	BOOST_LOG_TRIVIAL(info) << "nrow = " << NUM_ROWS << ", ncol = " << NUM_COLS;

	of.write((char*)&NUM_ROWS, sizeof(size_t)); // size_t rows
	of.write((char*)&NUM_COLS, sizeof(size_t)); // size_t cols
	of.write((char*)&mat[0], sizeof(double)*g_ncol*g_nrow);
	delete [] mat;
}

static void to_h2o(dense_matrix::ptr dmat, std::ofstream& of) {
	BOOST_LOG_TRIVIAL(info) << "Writing matrix ...";
	for (size_t row=0; row < dmat->get_num_rows(); row++) {
		std::shared_ptr<vector> curr_row = dmat->get_row(row);
		std::vector<double> stdvec = curr_row->conv2std<double>();
		of << row + 1;
		for (size_t col=0; col < stdvec.size(); col++) {
			of << "," << stdvec[col];
		}
		of << "\n";
	}
}

static void to_h2o(double* outmat, std::ofstream& of, const conv_layout lay) {
	BOOST_LOG_TRIVIAL(info) << "Writing matrix ...";
	if (lay == RAWCOL) {
		// Write it
		for (size_t row = 0; row < g_ncol; row++) {
			of << row + 1;
			for (size_t col = 0; col < g_nrow; col++) {
				of << "," << outmat[row*g_nrow+col];
			}
			of << "\n";
		}
	} else if (lay == RAWROW) {
		for (size_t row = 0; row < g_nrow; row++) {
			of << row+1;
			for (size_t col = 0; col < g_ncol; col++) {
				of << "," << outmat[row*g_nrow+col];
			}
			of << "\n";
		}
	} else { assert(0); }
	delete [] outmat;
}

int main(int argc, char* argv[]) {
	if (argc < 5) {
		fprintf(stderr, "usage: ./convert_matrix in_filename to_format"
			   " out_filename layout[row/col/rrow/rcol] [nrow] [ncol]");
		exit(-1);
	}

	std::string infile = (argv[1]);
	std::string to_format = argv[2];
	std::string out_filename = argv[3];
	std::string argv4 = std::string(argv[4]);
	conv_layout lay = argv4 == "row" ? ROW 
		: argv4 == "col" ? COL
		: argv4 == "rrow" ? RAWROW
		: RAWCOL;

	if ((lay == RAWROW || lay == RAWCOL) && argc != 7) {
		fprintf(stderr, "Must provide 2 more args for col-wise");
		exit(-1);
	}

	if (argc == 7) {
		g_nrow = atol(argv[5]);
		g_ncol = atol(argv[6]);
	}

	std::ofstream out_file;

	if (to_format == "h2o") {
		BOOST_LOG_TRIVIAL(info) << "Converting to " << to_format << " format ...";
		out_file.open(out_filename, std::ios::out);
		if (out_file.is_open()) {
			if (lay == ROW || lay == COL) {
				to_h2o(get_mat(infile), out_file);
			} else {
				to_h2o(bin_read(infile), out_file, lay);
			}
			out_file.close();
		} else { 
			BOOST_LOG_TRIVIAL(info) << "Failed to open " << out_filename;
			exit(911);
		}
	} else	if (to_format == "spark" || to_format == "dato") {
		BOOST_LOG_TRIVIAL(info) << "Converting to " << to_format << " format ...";
		out_file.open(out_filename, std::ios::out);
		if (out_file.is_open()) {
			if (lay == ROW || lay == COL) {
				printf("Shouldnt be here!\n"); exit(-1);
				to_spark(get_mat(infile), out_file);
			} else {
				to_spark(bin_read(infile), out_file, lay);
			}
			out_file.close();
		} else {
			BOOST_LOG_TRIVIAL(info) << "Failed to open " << out_filename;
			exit(911);
		}
	} else if (to_format == "kmp") {
		BOOST_LOG_TRIVIAL(info) << "Converting to " << to_format << " format ...";
		out_file.open(out_filename, std::ios::out);
		if (out_file.is_open()) {
			fprintf(stderr, "Failed to open file\n");
			if (lay == ROW || lay == COL) {
				to_kmeans_par(get_mat(infile), out_file);
			} else {
				to_kmeans_par(bin_read(infile), out_file, lay);
			}
			out_file.close();
		} else {
			BOOST_LOG_TRIVIAL(info) << "Failed to open " << out_filename;
			exit(911);
		}
	} else if (to_format == "fg") {
		BOOST_LOG_TRIVIAL(info) << "Converting to " << to_format << " format ...";
		out_file.open(out_filename, std::ios::binary | std::ios::trunc | std::ios::out);
		if (out_file.is_open()) {
			if (lay == ROW || lay == COL) {
				to_fg(get_mat(infile), out_file);
			} else {
				to_fg(bin_read(infile), out_file, lay);
			}
			BOOST_LOG_TRIVIAL(info) << "Conversion to fg complete";
			out_file.close();
		} else { 
			BOOST_LOG_TRIVIAL(info) << "Failed to open " << out_filename;
			exit(911);
		}

#if 0
		BOOST_LOG_TRIVIAL(info) << "Reading back fg";
		double* out_mat;
		out_mat = read_fg(out_filename, lay);
		delete [] out_mat;
#endif
	} else {
		fprintf(stderr, "Unknown format '%s'\n", to_format.c_str());
	}

	BOOST_LOG_TRIVIAL(info) << "Conversion complete! file is '" << argv[3] << "'\n";
	return EXIT_SUCCESS;
}
