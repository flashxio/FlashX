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

#include "convert.hpp"

namespace kpmutil = knor::utils;

int main(int argc, char* argv[]) {
    if (argc < 7) {
        fprintf(stderr, "usage: ./convert_matrix in_filename"
                " in_format [text/knori/knord/knors]\\\n"
                " \tout_filename out_format[text/knori/knord/knors]"
                " nsamples dim\n");
        exit(-1);
    }

    std::string infile = argv[1];
    std::string in_format = argv[2];
    std::string out_filename = argv[3];
    std::string out_format = argv[4];
    size_t nrow = atol(argv[5]);
    size_t ncol = atol(argv[6]);

    kpmutil::layout in_layout =
        in_format == "knori" ? kpmutil::BIN_RM :
        in_format == "knord" ? kpmutil::BIN_RM :
        in_format == "knors" ? kpmutil::SEM :
        in_format == "text" ? kpmutil::TEXT : kpmutil::INVALID;

    kpmutil::layout out_layout =
        out_format == "knori" ? kpmutil::BIN_RM :
        out_format == "knord" ? kpmutil::BIN_RM :
        out_format == "knors" ? kpmutil::SEM :
        out_format == "text" ? kpmutil::TEXT : kpmutil::INVALID;

    if (in_layout != kpmutil::TEXT && argc != 7) {
        fprintf(stderr, "Must provide `nrow` and `ncol` for "
                " format %s\n", in_format.c_str());
        exit(-1);
    }

    kpmutil::format_converter fc (infile, in_layout, nrow, ncol);
    fc.write(out_filename, out_layout);
    return EXIT_SUCCESS;
}