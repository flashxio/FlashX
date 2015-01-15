#ifndef __FG_SPARSE_MATRIX_H__
#define __FG_SPARSE_MATRIX_H__

/**
 * Copyright 2014 Open Connectome Project (http://openconnecto.me)
 * Written by Da Zheng (zhengda1936@gmail.com)
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
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

#include "graph_engine.h"
#include "FGlib.h"

namespace fg
{

class matrix_vertex: public compute_vertex
{
public:
	matrix_vertex(vertex_id_t id): compute_vertex(id) {
	}

	void run(vertex_program &prog) {
		vertex_id_t id = prog.get_vertex_id(*this);
		request_vertices(&id, 1);
	}

	void run(vertex_program &prog, const page_vertex &vertex) {
	}

	void run_on_message(vertex_program &prog, const vertex_message &msg) {
	}
};

/*
 * The vertex program for sparse matrix vector multiplication
 * on the general sparse matrix.
 */
template<class ResType, class MatEntryType>
class SPMV_vertex_program: public vertex_program_impl<matrix_vertex>
{
	edge_type type;
	const FG_vector<ResType> &input;
	FG_vector<ResType> &output;

	static safs::page_byte_array::seq_const_iterator<MatEntryType> get_data_iterator(
			const page_vertex &v, edge_type type) {
		if (v.is_directed())
			return ((const page_directed_vertex &) v).get_data_seq_it<MatEntryType>(type);
		else
			return ((const page_undirected_vertex &) v).get_data_seq_it<MatEntryType>();
	}
public:
	SPMV_vertex_program(edge_type type, const FG_vector<ResType> &_input,
			FG_vector<ResType> &_output): input(_input), output(_output) {
		this->type = type;
	}

	virtual void run(compute_vertex &, const page_vertex &vertex) {
		edge_seq_iterator it = vertex.get_neigh_seq_it(type, 0,
				vertex.get_num_edges(type));
		safs::page_byte_array::seq_const_iterator<MatEntryType> d_it
			= get_data_iterator(vertex, type);
		ResType w = 0;
		while (it.has_next()) {
			vertex_id_t neigh_id = it.next();
			MatEntryType val = d_it.next();
			w += input.get(neigh_id) * val;
		}
		output.set(vertex.get_id(), w);
	}
};

/*
 * The vertex program for sparse matrix vector multiplication
 * on the adjacency matrix.
 */
template<class ResType>
class SPMV_vertex_program<ResType, empty_data>: public vertex_program_impl<matrix_vertex>
{
	edge_type type;
	const FG_vector<ResType> &input;
	FG_vector<ResType> &output;
public:
	SPMV_vertex_program(edge_type type, const FG_vector<ResType> &_input,
			FG_vector<ResType> &_output): input(_input), output(_output) {
		this->type = type;
	}

	virtual void run(compute_vertex &, const page_vertex &vertex) {
		edge_seq_iterator it = vertex.get_neigh_seq_it(type, 0,
				vertex.get_num_edges(type));
		ResType w = 0;
		while (it.has_next()) {
			vertex_id_t neigh_id = it.next();
			w += input.get(neigh_id);
		}
		output.set(vertex.get_id(), w);
	}
};

template<class ResType, class MatEntryType>
class SPMV_vertex_program_creater: public vertex_program_creater
{
	const FG_vector<ResType> &input;
	FG_vector<ResType> &output;
	edge_type etype;
public:
	SPMV_vertex_program_creater(edge_type etype, const FG_vector<ResType> &_input,
			FG_vector<ResType> &_output): input(_input), output(_output) {
		this->etype = etype;
	}

	vertex_program::ptr create() const {
		return vertex_program::ptr(
				new SPMV_vertex_program<ResType, MatEntryType>(etype,
					input, output));
	}
};

template<class EntryType>
class FG_sparse_matrix
{
	size_t nrow;
	size_t ncol;

	// The type of edges that specifies the rows of the matrix.
	// For a symmetric matrix, the edge type does not matter.
	edge_type etype;
	graph_engine::ptr graph;

	FG_sparse_matrix() {
		etype = edge_type::NONE;
		nrow = 0;
		ncol = 0;
	}

protected:
	FG_sparse_matrix(FG_graph::ptr fg) {
		graph_index::ptr index = NUMA_graph_index<matrix_vertex>::create(
				fg->get_graph_header());
		graph = fg->create_engine(index);
		etype = edge_type::OUT_EDGE;
		this->nrow = graph->get_num_vertices();
		this->ncol = graph->get_num_vertices();
	}
public:
	typedef std::shared_ptr<FG_sparse_matrix<EntryType> > ptr;

	static ptr create(FG_graph::ptr fg) {
		return ptr(new FG_sparse_matrix<EntryType>(fg));
	}

	void resize(size_t nrow, size_t ncol) {
		this->nrow = nrow;
		this->ncol = ncol;
		assert(nrow <= graph->get_num_vertices());
		assert(ncol <= graph->get_num_vertices());
		// TODO we need to check if we can resize like this.
	}

	template<class T>
	void multiply(const FG_vector<T> &input, FG_vector<T> &output) const {
		assert(input.get_size() == get_num_cols());
		assert(output.get_size() == get_num_rows());
		graph->start_all(vertex_initializer::ptr(),
				vertex_program_creater::ptr(
					new SPMV_vertex_program_creater<T, EntryType>(
						etype, input, output)));
		graph->wait4complete();
	}

	size_t get_num_rows() const {
		return nrow;
	}

	size_t get_num_cols() const {
		return ncol;
	}

	typename FG_sparse_matrix<EntryType>::ptr transpose() const {
		typename FG_sparse_matrix<EntryType>::ptr t
			= typename FG_sparse_matrix<EntryType>::ptr(
					new FG_sparse_matrix<EntryType>());
		if (this->etype == IN_EDGE)
			t->etype = OUT_EDGE;
		else if (this->etype == OUT_EDGE)
			t->etype = IN_EDGE;
		else
			assert(0);
		t->graph = this->graph;
		t->nrow = this->ncol;
		t->ncol = this->nrow;
		return t;
	}
};

typedef FG_sparse_matrix<empty_data> FG_adj_matrix;

}

#endif
