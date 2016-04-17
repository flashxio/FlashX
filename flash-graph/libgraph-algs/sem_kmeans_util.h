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
#include <math.h>
//#include <cmath>

#include <vector>
#include <iostream>

#include <boost/assert.hpp>
#include "log.h"
#include "kmeans_types.h"

namespace {

    /*\Internal
     * \brief print a col wise matrix of type double / double.
     * Used for testing only.
     * \param matrix The col wise matrix.
     * \param rows The number of rows in the mat
     * \param cols The number of cols in the mat
     */
    template <typename T>
        static void print_mat(T* matrix, const unsigned rows, const unsigned cols) {
            for (unsigned row = 0; row < rows; row++) {
                std::cout << "[";
                for (unsigned col = 0; col < cols; col++) {
                    std::cout << " " << matrix[row*cols + col];
                }
                std::cout <<  " ]\n";
            }
        }

    template <typename T>
        static void print_arr(const T* arr, const unsigned len) {
            printf("[ ");
            for (unsigned i = 0; i < len; i++) {
                std::cout << arr[i] << " ";
            }
            printf("]\n");
        }

    double get_bic(const std::vector<double>& dist_v, const unsigned nrow,
            const unsigned ncol, const unsigned k) {
            double bic = 0;
#pragma omp parallel for reduction(+:bic) shared (dist_v)
        for (unsigned i = 0; i < dist_v.size(); i++) {
            bic += (dist_v[i] );
        }
        printf("Distance sum: %f\n", bic);

        return 2*bic + log(nrow)*ncol*k;
    }

    void spherical_projection(double* data, const unsigned nrow, const unsigned ncol) {
#pragma omp parallel for shared (data)
        for (unsigned row = 0; row < nrow; row++) {
            double norm2 = 0;
            for (unsigned col = 0; col < ncol; col++)
                norm2 += (data[row]*data[row]);
            sqrt(norm2);
            for (unsigned col = 0; col < ncol; col++)
                data[col] = data[col]/norm2;
        }
    }

    // For experiments
    void store_cluster(const unsigned id, const double* data,
            const unsigned numel, const unsigned* cluster_assignments,
            const unsigned nrow, const unsigned ncol, const std::string dir) {
        BOOST_LOG_TRIVIAL(info) <<
            "Storing cluster " << id;

        FILE* f;
        std::string fn = dir+"cluster_"+std::to_string(id)+
            "_r"+std::to_string(numel)+"_c"+std::to_string(ncol)+".bin";
        BOOST_VERIFY(f = fopen(fn.c_str(), "wb"));
        BOOST_LOG_TRIVIAL(info) << "[Warning]: Writing cluster file '" <<
            fn << "'";
        unsigned count = 0;

        for(unsigned i = 0; i < nrow; i++) {
            if (count == numel) { break; }
            if (cluster_assignments[i] == id) {
                BOOST_VERIFY(fwrite(&data[i*ncol],
                            (ncol*sizeof(double)), 1, f));
                count++;
            }
        }
        BOOST_VERIFY(count == numel);
        fclose(f);
    }


    // Vector equal function
    template <typename T>
        static bool v_eq(const T& lhs, const T& rhs) {
            return std::equal(lhs.begin(), lhs.end(), rhs.begin());
        }


    template <typename T>
        static const bool v_eq_const(const std::vector<T>& v, const T var) {
            for (unsigned i=0; i < v.size(); i++) {
                if (v[i] != var) return false;
            }
            return true;
        }

    template <typename T>
        bool eq_all(const T* v1, const T* v2, const unsigned len) {
            return (std::equal(&v1[0], &(v1[len-1]), &v2[0]));
        }

    template <typename T>
        static const double eucl_dist(const T* lhs, const T* rhs,
                const unsigned size) {
            double dist = 0;
            double diff;

            for (unsigned col = 0; col < size; col++) {
                diff = lhs[col] - rhs[col];
                dist += diff * diff;
            }
            return sqrt(dist);
        }

    template<typename T>
        const double cos_dist(const T* lhs, const T* rhs,
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
                bin_reader(const std::string fn, const size_t nrow, const size_t ncol) {
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
