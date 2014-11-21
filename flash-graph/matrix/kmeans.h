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

// #include "omp.h"

#include "FGlib.h"
#include "graph.h"
#include "graph_engine.h"
#include "graph_config.h"

#include "matrix/FG_sparse_matrix.h"
#include "matrix/FG_dense_matrix.h"
#include "matrix/matrix_eigensolver.h"

static const char * PROF_FILE_LOC = 
"/home/disa/graph-engine/flash-graph/clustering/profile_out.log";

// TODO: Doc
template <typename T>
void print_arr(T* arr, vsize_t len);

/* 
 * \Internal
 * \brief Simple helper used to print a vector.
 * \param v The vector to print.
 */
template <typename T>
void print_vector(typename std::vector<T>& v);

/**
 * \brief Get the squared distance given two values.
 *	\param arg1 the first value.
 *	\param arg2 the second value.
 *  \return the squared distance.
 */
ev_float_t dist_sq(ev_float_t arg1, ev_float_t arg2);

/* See: http://en.wikipedia.org/wiki/K-means_clustering#Initialization_methods */

/**
  * \brief This initializes clusters by randomly choosing sample 
  *		membership in a cluster.
  *	\param cluster_assignments Which cluster each sample falls into.
  */
void random_partition_init(vsize_t* cluster_assignments);

/**
  * \brief Forgy init takes `K` random samples from the matrix
  *		and uses them as cluster centers.
  * \param matrix the flattened matrix who's rows are being clustered.
  * \param clusters The cluster centers (means) flattened matrix.
  */
void forgy_init(const ev_float_t* matrix, ev_float_t* clusters);

/**
 * \brief A parallel version of the kmeans++ initialization alg.
 */
void kmeanspp_init();


/*\Internal
* \brief print a col wise matrix of type ev_float_t / double.
* Used for testing only.
* \param matrix The col wise matrix.
* \param rows The number of rows in the mat
* \param cols The number of cols in the mat
*/
template <typename T>
void print_mat(const T* matrix, const vsize_t rows, const vsize_t cols); 

/**
  * \brief Update the cluster assignments while recomputing distance matrix.
  * \param matrix The flattened matrix who's rows are being clustered.
  * \param clusters The cluster centers (means) flattened matrix.
  *	\param cluster_assignments Which cluster each sample falls into.
  */
void E_step(const ev_float_t* matrix, ev_float_t* clusters, vsize_t* cluster_assignments);

/**
 * \brief Update the cluster means
 * \param matrix The matrix who's rows are being clustered.
 * \param clusters The cluster centers (means).
 * \param cluster_assignment_counts How many members each cluster has.
 * \param cluster_assignments Which cluster each sample falls into.
 */
void M_step(const ev_float_t* matrix, ev_float_t* clusters, 
		vsize_t* cluster_assignment_counts, const vsize_t* cluster_assignments);

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
vsize_t compute_kmeans(const ev_float_t* matrix, ev_float_t* clusters, 
		vsize_t* cluster_assignments, vsize_t* cluster_assignment_counts,
		const vsize_t num_rows, const vsize_t nev, const size_t k, const vsize_t MAX_ITERS, 
		const std::string init="random");

/**
  * \brief Solely used to return the result of k means to R binding
  */
class FG_kmeans_ret {
	private:
		FG_kmeans_ret() {
		}
	public:
		FG_vector<vsize_t>::ptr assignments;
		FG_col_wise_matrix<ev_float_t>::ptr clusters;
		std::vector<vsize_t> sizes;
		vsize_t iters;

		FG_kmeans_ret(FG_vector<vsize_t>::ptr asgn, 
				FG_col_wise_matrix<ev_float_t>::ptr clusts, std::vector<vsize_t>& sizes, vsize_t iters)
		{
			this->assignments = asgn;
			this->clusters = clusts;
			this->sizes = sizes;
			this->iters = iters;
		}
};
#endif
