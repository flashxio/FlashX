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

#include "clusters.h"

//using namespace fg;

namespace {
    template <typename T>
        static double const eucl_dist(const T* lhs, const T* rhs,
                const unsigned size) {
            double dist = 0;
            //BOOST_VERIFY(lhs->size() == rhs->size());

            for (unsigned col = 0; col < size; col++) {
                double diff = lhs[col] - rhs[col];
                dist += diff * diff;
            }

            BOOST_VERIFY(dist >= 0);
            return sqrt(dist); // TODO: rm sqrt
        }

    template<typename T>
        double const cos_dist(const T* lhs, const T* rhs,
                const unsigned size) {
            T numr, ldenom, rdenom;
            numr = ldenom = rdenom = 0;

            for (unsigned col = 0; col < size; col++) {
                T a = lhs[col];
                T b = rhs[col];

                numr += a*b;
                ldenom += a*a;
                rdenom += b*b;
            }
            return  1 - (numr / ((sqrt(ldenom)*sqrt(rdenom))));
        }

    template <typename T>
        static void print_vector(typename std::vector<T> v, unsigned max_print=100) {
            unsigned print_len = v.size() > max_print ? max_print : v.size();

            std::cout << "[";
            typename std::vector<T>::iterator itr = v.begin();
            for (; itr != v.begin()+print_len; itr++) {
                std::cout << " "<< *itr;
            }

            if (v.size() > print_len) std::cout << " ...";
            std::cout <<  " ]\n";
        }

    // A very C-style binary data reader
    template <typename T>
        class bin_reader {
            private:
                FILE* f;
                size_t nrow, ncol;

                void cat(const T* arr) {
                    std::cout << "[ ";
                    for (size_t i = 0; i < ncol; i++) {
                        std::cout << arr[i] << " ";
                    }
                    std::cout << "]\n";
                }

            public:
                bin_reader(std::string fn, size_t nrow, size_t ncol) {
                    f = fopen(fn.c_str(), "rb");
                    BOOST_VERIFY(NULL != f);
                    this->nrow = nrow;
                    this->ncol = ncol;
                }

                // Read data and cat in a viewer friendly fashion
                void read_cat() {
                    T arr [ncol];
                    for (size_t i = 0; i < nrow; i++) {
                        size_t num_read = fread(&arr[0], sizeof(T)*ncol, 1, f);
                        BOOST_ASSERT_MSG(num_read == 1, "Error reading file!\n");
                        cat(arr);
                    }
                }

                std::vector<T> readline() {
                    std::vector<T> v;
                    v.resize(ncol);
                    size_t num_read = fread(&v[0], sizeof(T)*ncol, 1, f);
                    BOOST_ASSERT_MSG(num_read == 1, "Error reading file!\n");
                    return v;
                }

                void readline(T* v) {
                    size_t num_read = fread(&v[0], sizeof(T)*ncol, 1, f);
                    BOOST_ASSERT_MSG(num_read == 1, "Error reading file!\n");
                }

                // Read all the data!
                void read(std::vector<T>* v) {
                    size_t num_read = fread(&((*v)[0]), sizeof(T)*ncol*nrow, 1, f);
                    BOOST_ASSERT_MSG(num_read == 1, "Error reading file!\n");
                }

                // Read all the data!
                void read(T* v) {
                    size_t num_read = fread(&v[0], sizeof(T)*ncol*nrow, 1, f);
                    BOOST_ASSERT_MSG(num_read == 1, "Error reading file!\n");
                }

                ~bin_reader() {
                    fclose(f);
                }
        };
}
#endif
