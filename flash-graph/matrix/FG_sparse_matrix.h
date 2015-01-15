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

class adj_get_edge_iter
{
public:
	typedef int value_type;

	class iterator {
		edge_seq_iterator it;
	public:
		iterator(const page_vertex &v, edge_type type): it(
				v.get_neigh_seq_it(type, 0, v.get_num_edges(type))) {
		}

		bool has_next() {
			return it.has_next();
		}

		std::pair<vertex_id_t, int> next() {
			return std::pair<vertex_id_t, int>(it.next(), 1);
		}
	};

	iterator operator()(const page_vertex &v, edge_type type) const {
		return iterator(v, type);
	}
};

template<class T>
class general_get_edge_iter
{
	static safs::page_byte_array::seq_const_iterator<T> get_data_iterator(
			const page_vertex &v, edge_type type) {
		if (v.is_directed())
			return ((const page_directed_vertex &) v).get_data_seq_it<T>(type);
		else
			return ((const page_undirected_vertex &) v).get_data_seq_it<T>();
	}
public:
	typedef T value_type;

	class iterator {
		edge_seq_iterator n_it;
		safs::page_byte_array::seq_const_iterator<T> d_it;
	public:
		iterator(const page_vertex &v, edge_type type): n_it(
				v.get_neigh_seq_it(type, 0, v.get_num_edges(type))), d_it(
				get_data_iterator(v, type)) {
		}

		bool has_next() {
			return n_it.has_next();
		}

		std::pair<vertex_id_t, T> next() {
			return std::pair<vertex_id_t, T>(n_it.next(), d_it.next());
		}
	};

	iterator operator()(const page_vertex &v, edge_type type) {
		return iterator(v, type);
	}
};

inline static edge_type reverse_dir(edge_type type)
{
	switch(type) {
		case IN_EDGE:
			return OUT_EDGE;
		case OUT_EDGE:
			return IN_EDGE;
		default:
			ABORT_MSG("wrong edge type");
	}
}

/**
 * The vertex program for sparse matrix vector multiplication
 * on the adjacency matrix.
 */
template<class ResType, class GetEdgeIterator>
class SPMV_vertex_program: public vertex_program_impl<matrix_vertex>
{
	edge_type type;
	const FG_vector<ResType> &input;
	FG_vector<ResType> &output;
	GetEdgeIterator get_edge_iterator;
public:
	SPMV_vertex_program(edge_type type, const FG_vector<ResType> &_input,
			FG_vector<ResType> &_output): input(_input), output(_output) {
		this->type = type;
	}

	virtual void run(compute_vertex &, const page_vertex &vertex) {
		typename GetEdgeIterator::iterator it = get_edge_iterator(vertex,
				type);
		ResType w = 0;
		while (it.has_next()) {
			std::pair<vertex_id_t,typename GetEdgeIterator::value_type> p
				= it.next();
			w += input.get(p.first) * p.second;
		}
		output.set(vertex.get_id(), w);
	}
};

template<class ResType, class GetEdgeIterator>
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
				new SPMV_vertex_program<ResType, GetEdgeIterator>(etype,
					input, output));
	}
};

template<class GetEdgeIterator>
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
	typedef std::shared_ptr<FG_sparse_matrix<GetEdgeIterator> > ptr;

	static ptr create(FG_graph::ptr fg) {
		return ptr(new FG_sparse_matrix<GetEdgeIterator>(fg));
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
					new SPMV_vertex_program_creater<T, GetEdgeIterator>(
						etype, input, output)));
		graph->wait4complete();
	}

	size_t get_num_rows() const {
		return nrow;
	}

	size_t get_num_cols() const {
		return ncol;
	}

	typename FG_sparse_matrix<GetEdgeIterator>::ptr transpose() const {
		typename FG_sparse_matrix<GetEdgeIterator>::ptr t
			= typename FG_sparse_matrix<GetEdgeIterator>::ptr(
					new FG_sparse_matrix<GetEdgeIterator>());
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

typedef FG_sparse_matrix<adj_get_edge_iter> FG_adj_matrix;

}

#endif
