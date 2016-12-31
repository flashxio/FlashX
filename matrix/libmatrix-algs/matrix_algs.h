#ifndef __MATRIX_ALGS_H__
#define __MATRIX_ALGS_H__

/*
 * Copyright 2015 Open Connectome Project (http://openconnecto.me)
 * Written by Da Zheng (zhengda1936@gmail.com)
 *
 * This file is part of FlashMatrix.
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

#include "sparse_matrix.h"
#include "dense_matrix.h"

namespace fm
{

namespace alg
{

/**
  * \brief Compute PageRank.
  *
  * \param mat The adjacency matrix of a sparse graph.
  * \param max_niters The maximal number of iterations.
  * \param damping_factor The damping factor. Originally .85.
  * \param num_in_mem The number of vectors that allow to be in memory.
  * \return A one-column dense matrix with a PageRank value for each vertex
  * in the graph.
  */
dense_matrix::ptr PageRank(sparse_matrix::ptr mat, size_t max_niters,
		float damping_factor, size_t num_in_mem);

std::pair<dense_matrix::ptr, dense_matrix::ptr> NMF(sparse_matrix::ptr mat,
		size_t k, size_t max_niters, size_t num_in_mem);
}

}

#endif
