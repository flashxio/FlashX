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

#include "convert_util.h"
#include "libgraph-algs/sem_kmeans_util.h"

static size_t g_nrow, g_ncol;

// SPARK
static void to_spark(const std::string fn, std::ofstream& of, const conv_layout lay) {
    unsigned long size = g_nrow*g_ncol;
    BOOST_LOG_TRIVIAL(info) << "Malloc-ing matrix with size: " << size;

    double* outmat = new double [size];
    BOOST_LOG_TRIVIAL(info) << "Reading " << fn << ", with r:" << g_nrow
        << ", c: " << g_ncol;
    bin_reader<double> b(fn, g_nrow, g_ncol);
    b.read(outmat);

    BOOST_LOG_TRIVIAL(info) << "Writing matrix ...";
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
static void to_kmeans_par(const std::string fn, std::ofstream& of, const conv_layout lay) {
    unsigned long size = g_nrow*g_ncol;
    BOOST_LOG_TRIVIAL(info) << "Malloc-ing matrix with size: " << size;

    double* outmat = new double [size];
    BOOST_LOG_TRIVIAL(info) << "Reading " << fn << ", with r:" << g_nrow
        << ", c: " << g_ncol;
    bin_reader<double> b(fn, g_nrow, g_ncol);
    b.read(outmat);

    BOOST_LOG_TRIVIAL(info) << "Writing matrix ...";

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

static void to_fg(const std::string fn, std::ofstream& of, const conv_layout lay) {
    unsigned long size = g_nrow*g_ncol;
    BOOST_LOG_TRIVIAL(info) << "Malloc-ing matrix with size: " << size;

    double* outmat = new double [size];
    BOOST_LOG_TRIVIAL(info) << "Reading " << fn << ", with r:" << g_nrow
        << ", c: " << g_ncol;
    bin_reader<double> b(fn, g_nrow, g_ncol);
    b.read(outmat);

    BOOST_LOG_TRIVIAL(info) << "Writing matrix ...";

	const size_t NUM_ROWS = g_nrow;
	const size_t NUM_COLS = g_ncol;
	BOOST_LOG_TRIVIAL(info) << "nrow = " << NUM_ROWS << ", ncol = " << NUM_COLS;

	of.write((char*)&NUM_ROWS, sizeof(size_t)); // size_t rows
	of.write((char*)&NUM_COLS, sizeof(size_t)); // size_t cols
	of.write((char*)&outmat[0], sizeof(double)*size);
	delete [] outmat;
}

static void to_h2o(const std::string fn, std::ofstream& of,
        const conv_layout lay) {

    unsigned long size = g_nrow*g_ncol;
    BOOST_LOG_TRIVIAL(info) << "Malloc-ing matrix with size: " << size;

    double* outmat = new double [size];
    BOOST_LOG_TRIVIAL(info) << "Reading " << fn << ", with r:" << g_nrow
        << ", c: " << g_ncol;
    bin_reader<double> b(fn, g_nrow, g_ncol);
    b.read(outmat);

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
			of << row + 1;
			for (size_t col = 0; col < g_ncol; col++) {
				of << "," << outmat[row*g_ncol+col];
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

	std::string infile = argv[1];
	std::string to_format = argv[2];
	std::string out_filename = argv[3];
	std::string argv4 = std::string(argv[4]);
	conv_layout lay = argv4 == "rrow" ? RAWROW
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
			to_h2o(infile, out_file, lay);
			out_file.close();
		} else { 
			BOOST_LOG_TRIVIAL(info) << "Failed to open " << out_filename;
			exit(911);
		}
	} else	if (to_format == "spark" || to_format == "dato") {
		BOOST_LOG_TRIVIAL(info) << "Converting to " << to_format << " format ...";
		out_file.open(out_filename, std::ios::out);
		if (out_file.is_open()) {
			to_spark(infile, out_file, lay);
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
			to_kmeans_par(infile, out_file, lay);
			out_file.close();
		} else {
			BOOST_LOG_TRIVIAL(info) << "Failed to open " << out_filename;
			exit(911);
		}
	} else if (to_format == "fg") {
		BOOST_LOG_TRIVIAL(info) << "Converting to " << to_format << " format ...";
		out_file.open(out_filename, std::ios::binary | std::ios::trunc | std::ios::out);
		if (out_file.is_open()) {
			to_fg(infile, out_file, lay);
			BOOST_LOG_TRIVIAL(info) << "Conversion to fg complete";
			out_file.close();
		} else { 
			BOOST_LOG_TRIVIAL(info) << "Failed to open " << out_filename;
			exit(911);
		}

	} else {
		fprintf(stderr, "Unknown format '%s'\n", to_format.c_str());
	}

	BOOST_LOG_TRIVIAL(info) << "Conversion complete! file is '" << argv[3] << "'\n";
	return EXIT_SUCCESS;
}
