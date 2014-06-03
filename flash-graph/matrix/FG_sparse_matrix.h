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
	// The adjacency matrix. Each element is a binary.
	ADJACENCY,
	LAPLACIAN,
	// The general sparse matrix. Each element is a real value.
	GENERAL,
	INVALID,
};

class matrix_vertex: public compute_vertex
{
public:
	matrix_vertex() {
	}

	matrix_vertex(vertex_id_t id,
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

/**
 * The default vertex program for sparse matrix vector multiplication.
 * It works on the adjacency matrix.
 */
template<class T>
class SPMV_vertex_program: public vertex_program_impl<matrix_vertex>
{
	edge_type type;
	const FG_vector<T> &input;
	FG_vector<T> &output;
public:
	SPMV_vertex_program(edge_type type, const FG_vector<T> &_input,
			FG_vector<T> &_output): input(_input), output(_output) {
		this->type = type;
	}

	const T &get_input(off_t idx) const {
		return input.get(idx);
	}

	void set_output(off_t idx, const T &v) {
		output.set(idx, v);
	}

	edge_type get_edge_type() const {
		return type;
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

/**
 * This vertex program implements sparse matrix vector multiplication
 * on a general sparse matrix.
 */
template<class T>
class general_SPMV_vertex_program: public SPMV_vertex_program<T>
{
public:
	general_SPMV_vertex_program(edge_type type, const FG_vector<T> &input,
			FG_vector<T> &output): SPMV_vertex_program<T>(type, input, output) {
		assert(this->get_graph().get_graph_header().get_graph_type()
				== graph_type::DIRECTED);
	}

	virtual void run(compute_vertex &, const page_vertex &vertex1) {
		const page_directed_vertex &vertex = (const page_directed_vertex &) vertex1;
		page_byte_array::const_iterator<vertex_id_t> n_it
			= vertex.get_neigh_begin(this->get_edge_type());
		page_byte_array::const_iterator<vertex_id_t> n_end
			= vertex.get_neigh_end(this->get_edge_type());
		page_byte_array::const_iterator<T> d_it
			= vertex.get_data_begin<T>(this->get_edge_type());
		T w = 0;
		for (; n_it != n_end; ++n_it, ++d_it) {
			vertex_id_t id = *n_it;
			T v = *d_it;
			w += this->get_input(id) * v;
		}
		this->set_output(vertex.get_id(), w);
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
			case FG_sparse_type::GENERAL:
				return vertex_program::ptr(new general_SPMV_vertex_program<T>(
							etype, input, output));
			default:
				assert(0);
		}
	}
};

/**
 * The default vertex program that groups rows or columns of a sparse matrix
 * and aggregates rows or columns in each group.
 * It works on the adjacency matrix.
 */
template<class T, class AggOp>
class groupby_vertex_program: public vertex_program_impl<matrix_vertex>
{
public:
	typedef std::map<int, typename FG_vector<AggOp>::ptr> agg_map_t;

private:
	edge_type row_type;
	bool row_wise;
	const FG_vector<int> &labels;
	agg_map_t &agg_results;

	static edge_type reverse_dir(edge_type type) {
		switch(type) {
			case IN_EDGE:
				return OUT_EDGE;
			case OUT_EDGE:
				return IN_EDGE;
			default:
				assert(0);
		}
	}
public:
	groupby_vertex_program(edge_type row_type, bool row_wise,
			const FG_vector<int> &_labels,
			agg_map_t &_agg_results): labels(_labels), agg_results(_agg_results) {
		this->row_type = row_type;
		this->row_wise = row_wise;
	}

	edge_type get_edge_type() const {
		if (row_wise)
			return row_type;
		else
			return reverse_dir(row_type);
	}

	void aggregate(vertex_id_t id, T v) {
		int label = labels.get(id);
		typename agg_map_t::const_iterator it = agg_results.find(label);
		assert(it != agg_results.end());
		it->second->get(id) += v;
	}

	virtual void run(compute_vertex &, const page_vertex &vertex) {
		edge_seq_iterator it = vertex.get_neigh_seq_it(get_edge_type(),
				0, vertex.get_num_edges(get_edge_type()));
		PAGE_FOREACH(vertex_id_t, id, it) {
			aggregate(id, 1);
		} PAGE_FOREACH_END
	}
};

/**
 * This vertex program works on a general sparse matrix.
 */
template<class T, class AggOp>
class general_groupby_vertex_program: public groupby_vertex_program<T, AggOp>
{
	typedef typename groupby_vertex_program<T, AggOp>::agg_map_t agg_map_t;
public:
	general_groupby_vertex_program(edge_type row_type, bool row_wise,
			const FG_vector<int> &labels,
			agg_map_t &agg_results): groupby_vertex_program<T, AggOp>(
				row_type, row_wise, labels, agg_results) {
		assert(this->get_graph().get_graph_header().get_graph_type()
				== graph_type::DIRECTED);
	}

	virtual void run(compute_vertex &, const page_vertex &vertex1) {
		const page_directed_vertex &vertex = (const page_directed_vertex &) vertex1;
		page_byte_array::const_iterator<vertex_id_t> n_it
			= vertex.get_neigh_begin(this->get_edge_type());
		page_byte_array::const_iterator<vertex_id_t> n_end
			= vertex.get_neigh_end(this->get_edge_type());
		page_byte_array::const_iterator<T> d_it
			= vertex.get_data_begin<T>(this->get_edge_type());
		for (; n_it != n_end; ++n_it, ++d_it) {
			this->aggregate(*n_it, *d_it);
		}
	}
};

template<class T, class AggOp>
class groupby_vertex_program_creater: public vertex_program_creater
{
	edge_type row_type;
	FG_sparse_type mtype;
	bool row_wise;
	const FG_vector<int> &labels;
	typename groupby_vertex_program<T, AggOp>::agg_map_t &agg_results;
public:
	groupby_vertex_program_creater(edge_type row_type, FG_sparse_type mtype,
			bool row_wise, const FG_vector<int> &_labels,
			typename groupby_vertex_program<T, AggOp>::agg_map_t &_agg_results): labels(
				_labels), agg_results(_agg_results) {
		this->row_type = row_type;
		this->mtype = mtype;
		this->row_wise = row_wise;
	}

	vertex_program::ptr create() const {
		switch(mtype) {
			case FG_sparse_type::ADJACENCY:
				return vertex_program::ptr(new groupby_vertex_program<T, AggOp>(
							row_type, row_wise, labels, agg_results));
			case FG_sparse_type::GENERAL:
				return vertex_program::ptr(
						new general_groupby_vertex_program<T, AggOp>(
							row_type, row_wise, labels, agg_results));
			default:
				assert(0);
		}
	}
};

class FG_sparse_matrix
{
	// The type of matrix represented by the graph.
	FG_sparse_type mtype;
	// The type of edges that specifies the rows of the matrix.
	// For a symmetric matrix, the edge type does not matter.
	edge_type etype;
	graph_engine::ptr graph;

	FG_sparse_matrix() {
		etype = edge_type::NONE;
		mtype = FG_sparse_type::INVALID;
	}

protected:
	FG_sparse_matrix(FG_graph::ptr fg, FG_sparse_type mtype) {
		graph_index::ptr index = NUMA_graph_index<matrix_vertex>::create(
				fg->get_index_file());
		graph = graph_engine::create(fg->get_graph_file(),
				index, fg->get_configs());
		etype = edge_type::OUT_EDGE;
		this->mtype = mtype;
	}
public:
	typedef std::shared_ptr<FG_sparse_matrix> ptr;

	static ptr create(FG_graph::ptr fg, FG_sparse_type mtype) {
		return ptr(new FG_sparse_matrix(fg, mtype));
	}

	template<class T>
	void multiply(const FG_vector<T> &input, FG_vector<T> &output) const {
		graph->start_all(vertex_initiator::ptr(),
				vertex_program_creater::ptr(new SPMV_vertex_program_creater<T>(
						etype, mtype, input, output)));
		graph->wait4complete();
	}

	/**
	 * Group rows or columns based on labels and compute aggregation info
	 * of each group in each column or row. It returns a Kxn dense matrix
	 * where K is the number of groups and n is the number of columns or
	 * rows.
	 */
	template<class T, class AggOp>
	void group_by(const FG_vector<int> &labels, bool row_wise,
			std::map<int, typename FG_vector<AggOp>::ptr> &agg_results) {
		std::set<int> set;
		labels.unique(set);
		vsize_t vec_size;
		if (row_wise)
			vec_size = get_num_cols();
		else
			vec_size = get_num_rows();
		BOOST_FOREACH(int label, set) {
			agg_results.insert(std::pair<int, typename FG_vector<AggOp>::ptr>(
						label, FG_vector<AggOp>::create(vec_size)));
		}
		graph->start_all(vertex_initiator::ptr(),
				vertex_program_creater::ptr(
					new groupby_vertex_program_creater<T, AggOp>(
						etype, mtype, row_wise, labels, agg_results)));
		graph->wait4complete();
	}

	void group_by_mean(const FG_vector<int> &labels, bool row_wise,
			std::map<int, FG_vector<double>::ptr> &agg_results);

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
