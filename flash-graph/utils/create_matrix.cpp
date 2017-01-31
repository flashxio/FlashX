/*
 * Copyright 2016 neurodata (http://neurodata.io/)
 * Written by Disa Mhembere (disa@jhu.edu)
 *
 * This file is part of k-par-means
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY CURRENT_KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

// Create a test matrix to check the conversion to other formats

#include <math.h>
#include <assert.h>

#include <fstream>
#include <string>
#include <iostream>
#include <random>

#include "types.hpp"
#include "knors_index.h"

namespace kpmutil = kpmeans::utils;
constexpr unsigned rand_min = 1;
constexpr unsigned rand_max = 5;

template <typename T>
static void append_bin(const size_t nrow, const size_t ncol,
        std::ofstream& of) {
    assert(of.is_open());
    std::default_random_engine generator;
    std::uniform_real_distribution<T> distribution(rand_min, rand_max);

    std::cout << "Writing the matrix" << std::endl;

    for (size_t i = 0; i < nrow*ncol; i++) {
        T val = distribution(generator);
        of.write(reinterpret_cast<char*>(&val), sizeof(T));
    }
}

template <typename T>
static void append_text(const size_t nrow, const size_t ncol,
        std::ofstream& of) {

    assert(of.is_open());
    std::default_random_engine generator;
    std::uniform_real_distribution<T> distribution(rand_min, rand_max);

    std::cout << "Writing TEXT for nrow: " << nrow <<
        ", ncol: " << ncol << std::endl;

    for (size_t row = 0; row < nrow; row++) {
        if (row > 0)
            of << "\n";
        for (size_t col = 0; col < ncol; col++) {
            T val = distribution(generator);
            if (col > 0)
                of << " " << val;
            else
                of << val;
        }
    }
}

int main (int argc, char* argv[]) {
	if (argc < 4) {
		fprintf(stderr, "usage: ./create_matrix nrow ncol [knori/knord/knors/text]\n");
		exit(911);
	}

	const size_t nrow = atol(argv[1]);
	const size_t ncol = atol(argv[2]);
	std::string outfn = "matrix_r"+ std::to_string(nrow)+"_c"+std::to_string(ncol);
	std::string in_format = std::string(argv[3]);

    kpmutil::layout lay =
        in_format == "knori" ? kpmutil::BIN_RM :
        in_format == "knord" ? kpmutil::BIN_RM :
        in_format == "knors" ? kpmutil::SEM :
        in_format == "text" ? kpmutil::TEXT : kpmutil::INVALID;

    std::ofstream outfile;

    switch (lay) {
        case kpmutil::layout::BIN_RM:
            outfile.open(outfn+".dat", std::ios::binary | std::ios::out);
            append_bin<double>(nrow, ncol, outfile);
            break;
        case kpmutil::layout::SEM:
            {
                outfile.open(outfn+".adj", std::ios::binary | std::ios::out);
                graph_header header = make_graph_header(nrow, ncol);
                outfile.write(reinterpret_cast<char*>(&header), sizeof(header));
                append_bin<double>(nrow, ncol, outfile);
            }
            break;
        case kpmutil::layout::TEXT:
            outfile.open(outfn+".txt", std::ios::out | std::ios::app);
            append_text<double>(nrow, ncol, outfile);
            outfile.close();
            break;
        case kpmutil::layout::INVALID:
        fprintf(stderr, "Unknown format/layout '%s\n', in_format",
                in_format.c_str());
            break;
        default:
            exit(1);
    }

	return EXIT_SUCCESS;
}
