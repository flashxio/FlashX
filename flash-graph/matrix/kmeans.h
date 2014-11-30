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
#include <math.h>
#ifdef PROFILER
#include <gperftools/profiler.h>
#endif
#include <vector>
#include <limits> 

#include "FGlib.h"
#include "graph.h"
#include "graph_engine.h"
#include "graph_config.h"

#include "matrix/FG_sparse_matrix.h"
#include "matrix/FG_dense_matrix.h"
#include "matrix/matrix_eigensolver.h"

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
fg::vsize_t compute_kmeans(const double* matrix, double* clusters, 
		fg::vsize_t* cluster_assignments, fg::vsize_t* cluster_assignment_counts,
		const fg::vsize_t num_rows, const fg::vsize_t nev, const size_t k, 
		const fg::vsize_t MAX_ITERS, const std::string init="random");

#endif
