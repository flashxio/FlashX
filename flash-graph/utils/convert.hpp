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

#ifndef __KPM_UTILS_CONVERT_HPP__
#define __KPM_UTILS_CONVERT_HPP__

//#include <boost/log/trivial.hpp>
#include <fstream>
#include <string>
#include "types.hpp"

namespace kpmeans { namespace utils {

class format_converter {

private:
    double* data;
    std::string infile;
    layout inlayout;
    size_t nrow, ncol;

public:
    format_converter(const std::string infile, const layout inlayout,
            size_t nrow=0, size_t ncol=0);
    format_converter(double* data, size_t nrow=0, size_t ncol=0);

    void write(const std::string outfile, const layout outlayout);

    const size_t get_nrow() const {
        return nrow;
    }

    const size_t get_ncol() const {
        return ncol;
    }

    void set_nrow(const size_t nrow) {
        this->nrow = nrow;
    }

    void set_ncol(const size_t col) {
        this->ncol = ncol;
    }
};
} }
#endif
