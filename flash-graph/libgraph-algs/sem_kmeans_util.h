/*
 * Copyright 2015 Open Connectome Project (http://openconnecto.me)
 * Written by Disa Mhembere (disa@jhu.edu)
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
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY CURRENT_KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */
#ifndef __SEM_KMEANS_UTIL_H__
#define __SEM_KMEANS_UTIL_H__

#include <stdlib.h>
#include <stdio.h>

#include <vector>
#include <iostream>
#include <boost/assert.hpp>

using namespace fg;

namespace {

    // A very C-style binary data reader
    template <typename T>
        class bin_reader {
            private:
            FILE* f;
            unsigned nrow, ncol;

            void cat(const T* arr) {
                std::cout << "[ ";
                for (unsigned i = 0; i < ncol; i++) {
                    std::cout << arr[i] << " ";
                }
                std::cout << "]\n";
            }

            public:
            bin_reader(std::string fn, unsigned nrow, unsigned ncol) {
                f = fopen(fn.c_str(), "rb");
                BOOST_VERIFY(NULL != f);
                this->nrow = nrow;
                this->ncol = ncol;
            }

            // Read data and cat in a viewer friendly fashion
            void read_cat() {
                T arr [ncol];
                for (unsigned i = 0; i < nrow; i++) {
                    size_t num_read = fread(&arr[0], sizeof(T)*ncol, 1, f);
                    BOOST_ASSERT_MSG(num_read == 1, "Error reading file!\n");
                    cat(arr);
                }
            }

            std::vector<double> readline() {
                std::vector<double> v;
                v.resize(ncol);
                size_t num_read = fread(&v[0], sizeof(v[0])*ncol, 1, f);
                BOOST_ASSERT_MSG(num_read == 1, "Error reading file!\n");
                return v;
            }

            ~bin_reader() {
                fclose(f);
            }
        };
}
#endif
