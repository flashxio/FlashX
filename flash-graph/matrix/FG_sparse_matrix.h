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

enum FG_sparse_type
{
	ADJACENCY,
	LAPLACIAN,
	GENERAL,
	INVALID,
};

class SPMV_vertex: public compute_vertex
{
public:
	SPMV_vertex() {
	}

	SPMV_vertex(vertex_id_t id,
			const vertex_index &index): compute_vertex(id, index) {
	}

	void run(vertex_program &prog) {
		vertex_id_t id = get_id();
		request_vertices(&id, 1);
	}

	void run(vertex_program &prog, const page_vertex &vertex) {
	}

	void run_on_message(vertex_program &prog, const vertex_message &msg) {
	}
};

template<class T>
class SPMV_vertex_program: public vertex_program_impl<SPMV_vertex>
{
	edge_type type;
	const FG_vector<T> &input;
	FG_vector<T> &output;

	const T &get_input(off_t idx) const {
		return input.get(idx);
	}

	void set_output(off_t idx, const T &v) {
		output.set(idx, v);
	}

	edge_type get_edge_type() const {
		return type;
	}
public:
	typedef std::shared_ptr<SPMV_vertex_program<T> > ptr;

	static ptr cast2(vertex_program::ptr prog) {
		return std::static_pointer_cast<SPMV_vertex_program<T>,
			   vertex_program>(prog);
	}

	SPMV_vertex_program(edge_type type, const FG_vector<T> &_input,
			FG_vector<T> &_output): input(_input), output(_output) {
		this->type = type;
	}

	virtual void run(compute_vertex &, const page_vertex &vertex) {
		edge_seq_iterator it = vertex.get_neigh_seq_it(get_edge_type(),
				0, vertex.get_num_edges(get_edge_type()));
		T w = 0;
		PAGE_FOREACH(vertex_id_t, id, it) {
			w += get_input(id);
		} PAGE_FOREACH_END
		set_output(vertex.get_id(), w);
	}
};

template<class T>
class SPMV_vertex_program_creater: public vertex_program_creater
{
	const FG_vector<T> &input;
	FG_vector<T> &output;
	edge_type etype;
	FG_sparse_type mtype;
public:
	SPMV_vertex_program_creater(edge_type etype, FG_sparse_type mtype,
			const FG_vector<T> &_input, FG_vector<T> &_output): input(
				_input), output(_output) {
		this->etype = etype;
		this->mtype = mtype;
	}

	vertex_program::ptr create() const {
		switch(mtype) {
			case FG_sparse_type::ADJACENCY:
				return vertex_program::ptr(new SPMV_vertex_program<T>(
							etype, input, output));
			default:
				assert(0);
		}
	}
};

class FG_sparse_matrix
{
	// The type of matrix represented by the graph.
	FG_sparse_type mtype;
	// The type of edges that multiplication occurs.
	// For a symmetric matrix, the edge type does not matter.
	edge_type etype;
	graph_engine::ptr graph;

	FG_sparse_matrix() {
		etype = edge_type::NONE;
		mtype = FG_sparse_type::INVALID;
	}

protected:
	FG_sparse_matrix(FG_graph::ptr fg, FG_sparse_type mtype) {
		graph_index::ptr index = NUMA_graph_index<SPMV_vertex>::create(
				fg->get_index_file());
		graph = graph_engine::create(fg->get_graph_file(),
				index, fg->get_configs());
		etype = edge_type::OUT_EDGE;
		this->mtype = mtype;
	}
public:
	typedef std::shared_ptr<FG_sparse_matrix> ptr;

	template<class T>
	void multiply(const FG_vector<T> &input, FG_vector<T> &output) const {
		graph->start_all(vertex_initiator::ptr(),
				vertex_program_creater::ptr(new SPMV_vertex_program_creater<T>(
						etype, mtype, input, output)));
		graph->wait4complete();
	}

	size_t get_num_rows() const {
		return graph->get_num_vertices();
	}

	size_t get_num_cols() const {
		return graph->get_num_vertices();
	}

	FG_sparse_matrix::ptr transpose() const {
		FG_sparse_matrix::ptr t = FG_sparse_matrix::ptr(new FG_sparse_matrix());
		if (this->etype == IN_EDGE)
			t->etype = OUT_EDGE;
		else if (this->etype == OUT_EDGE)
			t->etype = IN_EDGE;
		else
			assert(0);
		t->graph = this->graph;
		t->mtype = this->mtype;
		return t;
	}
};

class FG_adj_matrix: public FG_sparse_matrix
{
	FG_adj_matrix(FG_graph::ptr fg): FG_sparse_matrix(fg,
			FG_sparse_type::ADJACENCY) {
	}
public:
	static ptr create(FG_graph::ptr fg) {
		return ptr(new FG_adj_matrix(fg));
	}
};

#endif
