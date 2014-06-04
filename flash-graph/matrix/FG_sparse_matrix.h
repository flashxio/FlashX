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

class adj_get_edge_iter
{
public:
	typedef int value_type;

	class iterator {
		page_byte_array::seq_const_iterator<vertex_id_t> it;
	public:
		iterator(const page_vertex &v, edge_type type): it(
				v.get_neigh_seq_it(type, 0, v.get_num_edges(type))) {
		}

		bool has_next() {
			return it.has_next();
		}

		vertex_id_t get_curr_id() const {
			return it.curr();
		}

		int get_curr_value() const {
			return 1;
		}

		void next() {
			it.next();
		}
	};

	iterator operator()(const page_vertex &v, edge_type type) const {
		return iterator(v, type);
	}
};

template<class T>
class general_get_edge_iter
{
public:
	typedef T value_type;

	class iterator {
		page_byte_array::seq_const_iterator<vertex_id_t> n_it;
		page_byte_array::seq_const_iterator<T> d_it;
	public:
		iterator(const page_vertex &v, edge_type type): n_it(
				v.get_neigh_seq_it(type, 0, v.get_num_edges(type))), d_it(
				// TODO it should also work for an undirected vertex.
				((const page_directed_vertex &) v).get_data_seq_it<T>(type)) {
		}

		bool has_next() {
			return n_it.has_next();
		}

		vertex_id_t get_curr_id() const {
			return n_it.curr();
		}

		T get_curr_value() const {
			return d_it.curr();
		}

		void next() {
			n_it.next();
			d_it.next();
		}
	};

	iterator operator()(const page_vertex &v, edge_type type) {
		return iterator(v, type);
	}
};

/**
 * The vertex program for sparse matrix vector multiplication.
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

	const ResType &get_input(off_t idx) const {
		return input.get(idx);
	}

	void set_output(off_t idx, const ResType &v) {
		output.set(idx, v);
	}

	edge_type get_edge_type() const {
		return type;
	}

	virtual void run(compute_vertex &, const page_vertex &vertex) {
		typename GetEdgeIterator::iterator it
			= get_edge_iterator(vertex, get_edge_type());
		ResType w = 0;
		while (it.has_next()) {
			vertex_id_t id = it.get_curr_id();
			typename GetEdgeIterator::value_type v = it.get_curr_value();
			w += get_input(id) * v;
			it.next();
		}
		set_output(vertex.get_id(), w);
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
				new SPMV_vertex_program<ResType, GetEdgeIterator>(
					etype, input, output));
	}
};

/**
 * The vertex program that groups rows or columns of a sparse matrix
 * and aggregates rows or columns in each group.
 */
template<class AggOp, class GetEdgeIterator>
class groupby_vertex_program: public vertex_program_impl<matrix_vertex>
{
public:
	typedef std::map<int, typename FG_vector<AggOp>::ptr> agg_map_t;

private:
	edge_type row_type;
	bool row_wise;
	const FG_vector<int> &labels;
	agg_map_t &agg_results;
	GetEdgeIterator get_edge_iterator;

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

	template<class T>
	void aggregate(vertex_id_t id, T v) {
		int label = labels.get(id);
		typename agg_map_t::const_iterator it = agg_results.find(label);
		assert(it != agg_results.end());
		it->second->get(id) += v;
	}

	virtual void run(compute_vertex &, const page_vertex &vertex) {
		typename GetEdgeIterator::iterator it
			= get_edge_iterator(vertex, get_edge_type());
		while (it.has_next()) {
			vertex_id_t id = it.get_curr_id();
			typename GetEdgeIterator::value_type v = it.get_curr_value();
			aggregate(id, v);
			it.next();
		}
	}
};

template<class AggOp, class GetEdgeIterator>
class groupby_vertex_program_creater: public vertex_program_creater
{
	typedef typename groupby_vertex_program<AggOp, GetEdgeIterator>::agg_map_t agg_map_t;
	edge_type row_type;
	bool row_wise;
	const FG_vector<int> &labels;
	agg_map_t &agg_results;
public:
	groupby_vertex_program_creater(edge_type row_type, bool row_wise,
			const FG_vector<int> &_labels, agg_map_t &_agg_results): labels(
				_labels), agg_results(_agg_results) {
		this->row_type = row_type;
		this->row_wise = row_wise;
	}

	vertex_program::ptr create() const {
		return vertex_program::ptr(
				new groupby_vertex_program<AggOp, GetEdgeIterator>(
					row_type, row_wise, labels, agg_results));
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
				fg->get_index_file());
		graph = graph_engine::create(fg->get_graph_file(),
				index, fg->get_configs());
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
		graph->start_all(vertex_initiator::ptr(),
				vertex_program_creater::ptr(
					new SPMV_vertex_program_creater<T, GetEdgeIterator>(
						etype, input, output)));
		graph->wait4complete();
	}

	/**
	 * Group rows or columns based on labels and compute aggregation info
	 * of each group in each column or row. It returns a Kxn dense matrix
	 * where K is the number of groups and n is the number of columns or
	 * rows.
	 */
	template<class AggOp>
	void group_by(const FG_vector<int> &labels, bool row_wise,
			std::map<int, typename FG_vector<AggOp>::ptr> &agg_results) {
		std::set<int> set;
		labels.unique(set);
		vsize_t vec_size;
		if (row_wise) {
			assert(labels.get_size() == get_num_rows());
			vec_size = get_num_cols();
		}
		else {
			assert(labels.get_size() == get_num_cols());
			vec_size = get_num_rows();
		}
		assert(vec_size == graph->get_num_vertices());
		BOOST_FOREACH(int label, set) {
			agg_results.insert(std::pair<int, typename FG_vector<AggOp>::ptr>(
						label, FG_vector<AggOp>::create(vec_size)));
		}
		graph->start_all(vertex_initiator::ptr(),
				vertex_program_creater::ptr(
					new groupby_vertex_program_creater<AggOp, GetEdgeIterator>(
						etype, row_wise, labels, agg_results)));
		graph->wait4complete();
	}

	void group_by_mean(const FG_vector<int> &labels, bool row_wise,
			std::map<int, FG_vector<double>::ptr> &agg_results) {
		group_by<double>(labels, row_wise, agg_results);

		count_map<int> cmap;
		labels.count_unique(cmap);

		typedef std::map<int, typename FG_vector<double>::ptr> sum_map_t;
		BOOST_FOREACH(typename sum_map_t::value_type v, agg_results) {
			int label = v.first;
			int count = cmap.get(label);
			typename FG_vector<double>::ptr sum = v.second;
			for (size_t i = 0; i < sum->get_size(); i++)
				sum->set(i, sum->get(i) / count);
		}
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

class FG_adj_matrix: public FG_sparse_matrix<adj_get_edge_iter>
{
	FG_adj_matrix(FG_graph::ptr fg): FG_sparse_matrix<adj_get_edge_iter>(fg) {
	}
public:
	static ptr create(FG_graph::ptr fg) {
		return ptr(new FG_adj_matrix(fg));
	}
};

template<class T>
class FG_general_sparse_matrix: public FG_sparse_matrix<general_get_edge_iter<T> >
{
	FG_general_sparse_matrix(
			FG_graph::ptr fg): FG_sparse_matrix<general_get_edge_iter<T> >(fg) {
	}
public:
	static typename FG_sparse_matrix<general_get_edge_iter<T> >::ptr create(
			FG_graph::ptr fg) {
		return typename FG_sparse_matrix<general_get_edge_iter<T> >::ptr(
				new FG_general_sparse_matrix<T>(fg));
	}
};

#endif
