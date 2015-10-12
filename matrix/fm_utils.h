#ifndef __FM_UTILS_H__
#define __FM_UTILS_H__

/*
 * Copyright 2014 Open Connectome Project (http://openconnecto.me)
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
#include "FGlib.h"

#include "sparse_matrix_format.h"
#include "data_frame.h"

namespace fm
{

/*
 * This function creates a row-major matrix from a data frame that
 * contains the locations of the non-zero entries in the matrix.
 * The input data frame needs to contain at least two vectors: the first
 * vector contains the row indexes and the second vector contains the column
 * indexes.
 * The order of the rows in the input data frame may be changed.
 * Each row of the output matrix is stored in the FlashGraph vertex format.
 * We can easily convert this to the FlashGraph format.
 */
vector_vector::ptr create_1d_matrix(data_frame::ptr df);

/*
 * This function creates an edge list stored in the data frame and converts
 * it into the FlashGraph format stored in memory.
 */
fg::FG_graph::ptr create_fg_graph(const std::string &graph_name,
		data_frame::ptr df, bool directed);

/*
 * This function creates a 2D-partitioned matrix from a data frame that
 * contains the locations of the non-zero entries in the matrix.
 * The matrix and its index are kept in memory.
 */
std::pair<SpM_2d_index::ptr, SpM_2d_storage::ptr> create_2d_matrix(
		data_frame::ptr df, const block_2d_size &block_size, size_t entry_size);
/*
 * This function creates a 2D-partitioned matrix from a vector of adjacency
 * lists. The matrix and its index are kept in memory.
 */
std::pair<SpM_2d_index::ptr, SpM_2d_storage::ptr> create_2d_matrix(
		vector_vector::ptr adjs, const block_2d_size &block_size,
		size_t entry_size);
/*
 * This function creates a 2D-partitioned matrix from a vector of adjacency
 * lists and stores the matrix and its index in files.
 */
void export_2d_matrix(vector_vector::ptr adjs, const block_2d_size &block_size,
		size_t entry_size, const std::string &mat_file,
		const std::string &mat_idx_file, bool to_safs);

/*
 * A set of functions that change the parameters in the matrix generator.
 */

void set_deduplicate(bool v);
void set_remove_self_edge(bool v);
};

#endif
