#ifndef __FG_UTILS_H__
#define __FG_UTILS_H__

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
#include "vertex_index.h"
#include "FGlib.h"

#include "sparse_matrix_format.h"
#include "data_frame.h"

namespace fg
{

/*
 * Get the total size of the in-edge lists.
 */
static inline size_t get_in_size(fg::vertex_index::ptr vindex)
{
	assert(vindex->get_graph_header().is_directed_graph());
	return vindex->get_out_part_loc() - vindex->get_header_size();
}

/*
 * Get the offset of the first out-edge list in the graph image.
 */
static inline size_t get_out_off(fg::vertex_index::ptr vindex)
{
	if (vindex->get_graph_header().is_directed_graph())
		return vindex->get_out_part_loc();
	else
		return vindex->get_header_size();
}

/*
 * Get the offset of the first in-edge list in the graph image.
 */
static inline size_t get_in_off(fg::vertex_index::ptr vindex)
{
	assert(vindex->get_graph_header().is_directed_graph());
	return vindex->get_header_size();
}

/*
 * Get the total size of the out-edge lists.
 */
size_t get_out_size(fg::vertex_index::ptr vindex);
/*
 * Get the offsets of all out-edge lists in the graph image.
 */
void init_out_offs(fg::vertex_index::ptr vindex, std::vector<off_t> &out_offs);
/*
 * Get the offsets of all in-edge lists in the graph image.
 */
void init_in_offs(fg::vertex_index::ptr vindex, std::vector<off_t> &in_offs);

/*
 * The data frame that stores the edge list of a graph has at least two columns:
 * the source vertices are stored in column "source";
 * the destination vertices are stored in column "dest";
 * if attributes exist, the attributes are stored in column "attr".
 */
class edge_list
{
	bool directed;
	fm::data_frame::const_ptr df;

	edge_list(fm::data_frame::const_ptr df, bool directed) {
		this->df = df;
		this->directed = directed;
	}
public:
	typedef std::shared_ptr<edge_list> ptr;
	typedef std::shared_ptr<const edge_list> const_ptr;

	static ptr create(fm::data_frame::ptr df, bool directed);

	size_t get_num_vecs() const {
		return df->get_num_vecs();
	}
	size_t get_num_edges() const {
		return df->get_num_entries();
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
	const fm::scalar_type &get_attr_type() const {
		assert(has_attr());
		return df->get_vec(2)->get_type();
	}
	size_t get_attr_size() const;

	fm::detail::vec_store::const_ptr get_source() const {
		return df->get_vec(0);
	}
	fm::detail::vec_store::const_ptr get_dest() const {
		return df->get_vec(1);
	}
	fm::detail::vec_store::const_ptr get_attr() const {
		if (df->get_num_vecs() > 2)
			return df->get_vec(2);
		else
			return fm::detail::vec_store::const_ptr();
	}

	edge_list::ptr sort_source() const;
	fm::vector_vector::ptr groupby_source(
			const fm::gr_apply_operate<fm::sub_data_frame> &op) const;
	edge_list::ptr reverse_edge() const;
};

fg::FG_graph::ptr construct_FG_graph(
		const std::pair<fg::vertex_index::ptr, fm::detail::vec_store::ptr> &g,
		const std::string &graph_name);

/*
 * This function creates a row-major matrix from a data frame, which is stored
 * in a vector of vectors. Each row is stored in a byte vector with the FlashGraph
 * vertex format (ext_mem_undirected_vertex).
 *
 * This function outputs a vector of vectors that contains the sparse matrix
 * and a scalar that indicates the number of columns in the sparse matrix.
 */
std::pair<fm::vector_vector::ptr, size_t> create_1d_matrix(edge_list::ptr el);

/*
 * This function creates an edge list stored in the data frame and converts
 * it into the FlashGraph format stored in memory.
 */
fg::FG_graph::ptr create_fg_graph(const std::string &graph_name,
		edge_list::ptr el);

/*
 * This prints a graph into an edge list format.
 */
void print_graph_el(fg::FG_graph::ptr, const std::string &delim,
		const std::string &edge_attr_type, FILE *f);

/*
 * Fetch the subgraph that contains the specified vertices from the graph.
 * If `compact' is false, the constructed subgraph contains the same number
 * of vertices as the original graph but the unspecified vertices contains
 * no neighbors.
 * If `compact' is true, the constructed subgraph contains only the specified
 * vertices.
 */
fg::FG_graph::ptr fetch_subgraph(fg::FG_graph::ptr graph,
		const std::vector<fg::vertex_id_t> &vertices,
		const std::string &graph_name, bool compact);

/*
 * This function creates a 2D-partitioned matrix from a data frame that
 * contains the locations of the non-zero entries in the matrix.
 * The matrix and its index are kept in memory.
 */
std::pair<fm::SpM_2d_index::ptr, fm::SpM_2d_storage::ptr> create_2d_matrix(
		edge_list::ptr el, const fm::block_2d_size &block_size,
		const fm::scalar_type *entry_type);
/*
 * This function creates a 2D-partitioned matrix from a vector of adjacency
 * lists. The matrix and its index are kept in memory.
 * We can't detect the number of columns in the matrix easily, so users
 * should indicate the number of columns in the matrix.
 */
std::pair<fm::SpM_2d_index::ptr, fm::SpM_2d_storage::ptr> create_2d_matrix(
		fm::vector_vector::ptr adjs, size_t num_cols,
		const fm::block_2d_size &block_size, const fm::scalar_type *entry_type);
/*
 * This function creates a 2D-partitioned matrix from a vector of adjacency
 * lists and stores the matrix and its index in files.
 * We can't detect the number of columns in the matrix easily, so users
 * should indicate the number of columns in the matrix.
 */
void export_2d_matrix(fm::vector_vector::ptr adjs, size_t num_cols,
		const fm::block_2d_size &block_size, const fm::scalar_type *entry_type,
		const std::string &mat_file, const std::string &mat_idx_file,
		bool to_safs);

/*
 * A set of functions that change the parameters in the matrix generator.
 */

void set_deduplicate(bool v);
void set_remove_self_edge(bool v);
}

#endif
