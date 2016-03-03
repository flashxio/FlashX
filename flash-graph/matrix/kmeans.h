#ifndef __FG_KMEANS_H__
#define __FG_KMEANS_H__

/*
 * Copyright 2014 Open Connectome Project (http://openconnecto.me)
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

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <omp.h>
#include <string.h>
#include <sys/time.h>
#ifdef PROFILER
#include <gperftools/profiler.h>
#endif
#include <vector>
#include <limits>
#include <iostream>
#include <algorithm>
#include <boost/assert.hpp>

#include "log.h"
#include "common.h"
#include "libgraph-algs/kmeans_types.h"
#include "libgraph-algs/sem_kmeans_util.h"
#include "libgraph-algs/clusters.h"

namespace {
    /**
     * \brief Print an arry of some length `len`.
     *	 \param len The length of the array.
     */
    template <typename T>
        static void print_arr(T* arr, unsigned len) {
            printf("[ ");
            for (unsigned i = 0; i < len; i++) {
                std::cout << arr[i] << " ";
            }
            printf("]\n");
        }

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

    static dist_type_t g_dist_type;

    /** /brief Choose the correct distance function and return it
     * /param arg0 A pointer to data
     * /param arg1 Another pointer to data
     * /param len The number of elements used in the comparison
     * /return the distance based on the chosen distance metric
     */
    double dist_comp_raw(const double* arg0, const double* arg1,
            const unsigned len) {
        if (g_dist_type == EUCL)
            return eucl_dist(arg0, arg1, len);
        else if (g_dist_type == COS)
            return cos_dist(arg0, arg1, len);
        else
            BOOST_ASSERT_MSG(false, "Unknown distance metric!");
        exit(EXIT_FAILURE);
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
}

namespace fg {
/**
 * \brief Compute kmeans on matrix of features
 * \param matrix The matrix who's row IDs are being clustered.
 * \param clusters The cluster centers (means).
 * \param cluster_assignments Which cluster each sample falls into.
 * \param cluster_assignment_counts How many members each cluster has.
 * \param num_rows The number of rows in `matrix`.
 * \param nev The number of eigenvalues / number of columns in `matrix`.
 * \param k The number of clusters required.
 * \param max_iters The maximum number of iterations of K-means to perform.
 * \param init The type of initilization ["random", "forgy", "kmeanspp"]
 **/
unsigned compute_kmeans(const double* matrix, double* clusters,
		unsigned* cluster_assignments, unsigned* cluster_assignment_counts,
		const unsigned num_rows, const unsigned num_cols, const unsigned k,
		const unsigned MAX_ITERS, const int max_threads,
        const std::string init="kmeanspp", const double tolerance=-1,
        const std::string dist_type="eucl");

/** See `compute_kmeans` for argument list */
unsigned compute_min_kmeans(const double* matrix, double* clusters_ptr,
        unsigned* cluster_assignments, unsigned* cluster_assignment_counts,
		const unsigned num_rows, const unsigned num_cols, const unsigned k,
        const unsigned MAX_ITERS, const int max_threads,
        const std::string init="kmeanspp", const double tolerance=-1,
        const std::string dist_type="eucl");
}
#endif
