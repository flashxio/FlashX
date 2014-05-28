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

	void run(vertex_program &prog, const page_vertex &vertex);

	void run_on_message(vertex_program &prog, const vertex_message &msg) {
	}
};

class SPMV_vertex_program: public vertex_program_impl<SPMV_vertex>
{
	edge_type type;
public:
	SPMV_vertex_program(edge_type type) {
		this->type = type;
	}

	edge_type get_edge_type() const {
		return type;
	}

	virtual void compute(vertex_id_t vertex_id, edge_seq_iterator &it) = 0;
};

inline void SPMV_vertex::run(vertex_program &prog, const page_vertex &vertex)
{
	SPMV_vertex_program &SPMV_vprog = (SPMV_vertex_program &) prog;
	edge_seq_iterator it = vertex.get_neigh_seq_it(SPMV_vprog.get_edge_type(),
			0, vertex.get_num_edges(SPMV_vprog.get_edge_type()));
	SPMV_vprog.compute(get_id(), it);
}

template<class T>
class T_SPMV_vertex_program: public SPMV_vertex_program
{
	const FG_vector<T> &input;
	FG_vector<T> &output;

	const T &get_input(off_t idx) const {
		return input.get(idx);
	}

	void set_output(off_t idx, const T &v) {
		output.set(idx, v);
	}
public:
	typedef std::shared_ptr<T_SPMV_vertex_program<T> > ptr;

	static ptr cast2(vertex_program::ptr prog) {
		return std::static_pointer_cast<T_SPMV_vertex_program<T>, vertex_program>(
				prog);
	}

	T_SPMV_vertex_program(edge_type type, const FG_vector<T> &_input,
			FG_vector<T> &_output): SPMV_vertex_program(type), input(
				_input), output(_output) {
	}

	void compute(vertex_id_t vertex_id, edge_seq_iterator &it) {
		T w = 0;
		PAGE_FOREACH(vertex_id_t, id, it) {
			w += get_input(id);
		} PAGE_FOREACH_END
		set_output(vertex_id, w);
	}
};

template<class T>
class SPMV_vertex_program_creater: public vertex_program_creater
{
	const FG_vector<T> &input;
	FG_vector<T> &output;
	edge_type type;
public:
	SPMV_vertex_program_creater(edge_type type, const FG_vector<T> &_input,
			FG_vector<T> &_output): input(_input), output(_output) {
		this->type = type;
	}

	vertex_program::ptr create() const {
		return vertex_program::ptr(new T_SPMV_vertex_program<T>(type, input, output));
	}
};

class FG_sym_adj_matrix
{
	graph_engine::ptr graph;

	FG_sym_adj_matrix(FG_graph::ptr fg) {
		graph_index::ptr index = NUMA_graph_index<SPMV_vertex>::create(
				fg->get_index_file());
		graph = graph_engine::create(fg->get_graph_file(),
				index, fg->get_configs());
	}
public:
	typedef std::shared_ptr<FG_sym_adj_matrix> ptr;

	static ptr create(FG_graph::ptr fg) {
		return ptr(new FG_sym_adj_matrix(fg));
	}

	template<class T>
	void multiply(const FG_vector<T> &input, FG_vector<T> &output) const {
		graph->start_all(vertex_initiator::ptr(),
				vertex_program_creater::ptr(new SPMV_vertex_program_creater<T>(
						edge_type::BOTH_EDGES, input, output)));
		graph->wait4complete();
	}

	size_t get_num_rows() const {
		return graph->get_num_vertices();
	}

	size_t get_num_cols() const {
		return graph->get_num_vertices();
	}
};

#endif
