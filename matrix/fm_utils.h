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
 * The data frame that stores the edge list of a graph has at least two columns:
 * the source vertices are stored in column "source";
 * the destination vertices are stored in column "dest";
 * if attributes exist, the attributes are stored in column "attr".
 */
class edge_list
{
	bool directed;
	data_frame::const_ptr df;

	edge_list(data_frame::const_ptr df, bool directed) {
		this->df = df;
		this->directed = directed;
	}
public:
	typedef std::shared_ptr<edge_list> ptr;
	typedef std::shared_ptr<const edge_list> const_ptr;

	static ptr create(data_frame::ptr df, bool directed);

	size_t get_num_vecs() const {
		return df->get_num_vecs();
	}
	bool is_in_mem() const {
		return df->get_vec(0)->is_in_mem();
	}
	bool is_directed() const {
		return directed;
	}
	bool has_attr() const {
		return df->get_num_vecs() > 2;
	}
	const scalar_type &get_attr_type() const {
		assert(has_attr());
		return df->get_vec(2)->get_type();
	}
	size_t get_attr_size() const;

	edge_list::ptr sort_source() const;
	vector_vector::ptr groupby_source(
			const gr_apply_operate<sub_data_frame> &op) const;
	edge_list::ptr reverse_edge() const;
};

/*
 * This function creates a row-major matrix from a data frame, which is stored
 * in a vector of vectors. Each row is stored in a byte vector with the FlashGraph
 * vertex format (ext_mem_undirected_vertex).
 *
 * This function outputs a vector of vectors that contains the sparse matrix
 * and a scalar that indicates the number of columns in the sparse matrix.
 */
std::pair<vector_vector::ptr, size_t> create_1d_matrix(edge_list::ptr el);

/*
 * This function creates an edge list stored in the data frame and converts
 * it into the FlashGraph format stored in memory.
 */
fg::FG_graph::ptr create_fg_graph(const std::string &graph_name,
		edge_list::ptr el);

/*
 * This function creates a 2D-partitioned matrix from a data frame that
 * contains the locations of the non-zero entries in the matrix.
 * The matrix and its index are kept in memory.
 */
std::pair<SpM_2d_index::ptr, SpM_2d_storage::ptr> create_2d_matrix(
		edge_list::ptr el, const block_2d_size &block_size,
		const scalar_type *entry_type);
/*
 * This function creates a 2D-partitioned matrix from a vector of adjacency
 * lists. The matrix and its index are kept in memory.
 * We can't detect the number of columns in the matrix easily, so users
 * should indicate the number of columns in the matrix.
 */
std::pair<SpM_2d_index::ptr, SpM_2d_storage::ptr> create_2d_matrix(
		vector_vector::ptr adjs, size_t num_cols,
		const block_2d_size &block_size, const scalar_type *entry_type);
/*
 * This function creates a 2D-partitioned matrix from a vector of adjacency
 * lists and stores the matrix and its index in files.
 * We can't detect the number of columns in the matrix easily, so users
 * should indicate the number of columns in the matrix.
 */
void export_2d_matrix(vector_vector::ptr adjs, size_t num_cols,
		const block_2d_size &block_size, const scalar_type *entry_type,
		const std::string &mat_file, const std::string &mat_idx_file,
		bool to_safs);

/*
 * A set of functions that change the parameters in the matrix generator.
 */

void set_deduplicate(bool v);
void set_remove_self_edge(bool v);
};

#endif
