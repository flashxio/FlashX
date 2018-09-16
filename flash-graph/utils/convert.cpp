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
#include "../../../../libkcommon/io.hpp"
#include "knors_index.h"

namespace kpmbase = knor::base;
namespace knor { namespace utils {

format_converter::format_converter(const std::string infile,
        const layout inlayout, size_t nrow, size_t ncol) {
    this->infile = infile;
    this->inlayout = inlayout;
    this->nrow = nrow;
    this->ncol = ncol;
}

format_converter::format_converter(double* data,
        size_t nrow, size_t ncol) {
    this->data = data;
    this->nrow = nrow;
    this->ncol = ncol;
}

void format_converter::write(const std::string outfile,
        const layout outlayout) {

    std::vector<double> row(this->get_ncol());
    std::ofstream of;

    switch (inlayout) {
        case TEXT:
            {
                kpmbase::text_reader<double> rdr(infile);
                rdr.set_ncol(this->get_ncol());
                if (outlayout == BIN_RM) {
                    std::cout << "Reading text, writing binary" << std::endl;

                    of.open(outfile, std::ios::out | std::ios::binary);
                    while (rdr.readline(row)) {
                        of.write(reinterpret_cast<char*>(&row[0]),
                                row.size()*sizeof(double));
                    }
                } else if (outlayout == SEM) {
                    std::cout << "Reading text, writing SEM" << std::endl;

                    of.open(outfile, std::ios::binary | std::ios::out);
                    graph_header header =
                        make_graph_header(get_nrow(), get_ncol());
                    of.write(reinterpret_cast<char*>(&header), sizeof(header));

                    while (rdr.readline(row)) {
                        of.write(reinterpret_cast<char*>(&row[0]),
                                row.size()*sizeof(double));
                    }
                } else {
                    BOOST_ASSERT_MSG(false, "Unknown output format\n");
                }
            }
            break;
        case SEM:
        case BIN_RM:
            {
                kpmbase::bin_rm_reader<double> rdr(infile);
                rdr.set_ncol(this->get_ncol());

                if (inlayout == SEM) {
                    std::cout << "Reading SEM, ";
                    printf("Seeking %d bytes\n", fg::graph_header::HEADER_SIZE);
                    rdr.seek(fg::graph_header::HEADER_SIZE);
                } else {
                    std::cout << "Reading binary, ";
                }

                if (outlayout == TEXT) {
                    std::cout << "writing text" << std::endl;

                    of.open(outfile, std::ios::out);
                    while (rdr.readline(row)) {
                        for (size_t col = 0; col < get_ncol(); col++)
                            of << row[col] << " ";
                        of << "\n";
                    }
                } else if (outlayout == SEM || outlayout == BIN_RM) {

                    if (outlayout == SEM) {
                        std::cout << "writing SEM" << std::endl;

                        of.open(outfile, std::ios::binary | std::ios::out);
                        graph_header header =
                            make_graph_header(get_nrow(), get_ncol());
                        of.write(reinterpret_cast<char*>(&header),
                                sizeof(header));
                    } else {
                        std::cout << "writing binary" << std::endl;
                    }

                    while (rdr.readline(row)) {
                        of.write(reinterpret_cast<char*>(&row[0]),
                                row.size()*sizeof(double));
                    }
                } else {
                    BOOST_ASSERT_MSG(false, "Unknown output format\n");
                }
            }
            break;
        case BIN_CM:
            BOOST_ASSERT_MSG(false, "Not implemented yet!\n");
            break;
        case INVALID:
            BOOST_ASSERT_MSG(false, "!\n");
            break;
    }
    of.close();
}
} } // End namespace knor::util