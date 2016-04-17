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

namespace {
    static const unsigned INVALID_CLUSTER_ID = std::numeric_limits<unsigned>::max();

    /** /brief Choose the correct distance function and return it
     * /param arg0 A pointer to data
     * /param arg1 Another pointer to data
     * /param len The number of elements used in the comparison
     * /return the distance based on the chosen distance metric
     */
    double dist_comp_raw(const double* arg0, const double* arg1,
            const unsigned len, km::dist_type_t dt) {
        if (dt == km::dist_type_t::EUCL)
            return eucl_dist(arg0, arg1, len);
        else if (dt == km::dist_type_t::COS)
            return cos_dist(arg0, arg1, len);
        else
            BOOST_ASSERT_MSG(false, "Unknown distance metric!");
        exit(EXIT_FAILURE);
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
