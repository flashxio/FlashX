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
#ifdef PROFILER
#include <gperftools/profiler.h>
#endif

#include "FG_vector.h"
#include "FGlib.h"

#include "scan_graph.h"

using namespace fg;

namespace {

scan_stage_t scan_stage;

class count_msg: public vertex_message
{
	size_t num;
public:
	count_msg(size_t num): vertex_message(sizeof(count_msg), false) {
		this->num = num;
	}

	size_t get_num() const {
		return num;
	}
};

class extended_neighbor_list: public neighbor_list
{
	std::vector<uint32_t> count_list;
public:
	extended_neighbor_list(const page_vertex &vertex,
			const std::vector<attributed_neighbor> &neighbors): neighbor_list(
				vertex, neighbors) {
		count_list.resize(this->size());
	}

	uint32_t get_count(size_t idx) const {
		return count_list[idx];
	}

	size_t count_edges(const page_vertex *v);

	off_t find_idx(vertex_id_t id) const {
		edge_set_t::const_iterator it = neighbor_set->find(id);
		if (it == neighbor_set->end())
			return -1;
		else {
			size_t idx = (*it).get_idx();
			assert(idx < id_list.size());
			return idx;
		}
	}

	attributed_neighbor find(vertex_id_t id) const {
		edge_set_t::const_iterator it = neighbor_set->find(id);
		if (it == neighbor_set->end())
			return attributed_neighbor();
		else {
			size_t idx = (*it).get_idx();
			assert(idx < id_list.size());
			return at(idx);
		}
	}

	attributed_neighbor at(size_t idx) const {
		return attributed_neighbor(id_list[idx], num_dup_list[idx]);
	}
};

size_t extended_neighbor_list::count_edges(const page_vertex *v)
{
	assert(!this->empty());
	if (v->get_num_edges(edge_type::BOTH_EDGES) == 0)
		return 0;

	std::vector<vertex_id_t> common_neighs1;
	std::vector<vertex_id_t> common_neighs2;
	size_t ret = neighbor_list::count_edges(v, edge_type::IN_EDGE, &common_neighs1)
		+ neighbor_list::count_edges(v, edge_type::OUT_EDGE, &common_neighs2);

	class skip_self {
	public:
		bool operator()(vertex_id_t id) {
			return false;
		}
	};

	class merge_edge {
	public:
		vertex_id_t operator()(vertex_id_t id1, vertex_id_t id2) {
			assert(id1 == id2);
			return id1;
		}
	};

	std::vector<vertex_id_t> common_neighs(common_neighs1.size()
			+ common_neighs2.size());
	size_t num_neighbors = unique_merge(
			common_neighs1.begin(), common_neighs1.end(),
			common_neighs2.begin(), common_neighs2.end(),
			skip_self(), merge_edge(), common_neighs.begin());
	common_neighs.resize(num_neighbors);

	// The number of duplicated edges between v and this vertex.
	off_t neigh_off = this->find_idx(v->get_id());
	assert(neigh_off >= 0);
	attributed_neighbor neigh = this->at(neigh_off);
	int num_v_dups = neigh.get_num_dups();
	assert(num_v_dups > 0);
	size_t num_edges = 0;
	for (std::vector<vertex_id_t>::const_iterator it = common_neighs.begin();
			it != common_neighs.end(); it++) {
		off_t n_off = this->find_idx(*it);
		assert(n_off >= 0);
		attributed_neighbor n = this->at(n_off);
		num_edges += n.get_num_dups();
		count_list[n_off] += num_v_dups;
	}
	if (num_edges > 0)
		count_list[neigh_off] += num_edges;
	return ret;
}

class local_scan_vertex: public compute_vertex
{
	vsize_t degree;
	scan_multi_func_value local_value;
public:
	local_scan_vertex(vertex_id_t id): compute_vertex(id) {
		degree = 0;
		local_value.set_real_local(0);
	}

	vsize_t get_degree() const {
		return degree;
	}

	bool has_local_scan() const {
		return local_value.has_real_local();
	}

	size_t get_local_scan() const {
		return local_value.get_real_local();
	}

	void run(vertex_program &prog) {
		vertex_id_t id = prog.get_vertex_id(*this);
		if (scan_stage == scan_stage_t::INIT)
			request_vertex_headers(&id, 1);
		else if (scan_stage == scan_stage_t::RUN)
			request_vertices(&id, 1);
	}

	void run(vertex_program &prog, const page_vertex &vertex) {
		if (vertex.get_id() == prog.get_vertex_id(*this))
			run_on_itself(prog, vertex);
		else
			run_on_neighbor(prog, vertex);
	}

	void run_on_itself(vertex_program &prog, const page_vertex &vertex);
	void run_on_neighbor(vertex_program &prog, const page_vertex &vertex);

	void run_on_message(vertex_program &prog, const vertex_message &msg1) {
		const count_msg &msg = (const count_msg &) msg1;
		if (local_value.has_runtime_data())
			local_value.get_runtime_data()->local_scan += msg.get_num();
		else
			local_value.inc_real_local(msg.get_num());
	}

	void run_on_vertex_header(vertex_program &prog, const vertex_header &header) {
		assert(prog.get_vertex_id(*this) == header.get_id());
		degree = header.get_num_edges();
	}

	void finding_triangles_end(vertex_program &prog, runtime_data_t *data) {
		extended_neighbor_list *neighbors
			= (extended_neighbor_list *) data->neighbors.get();
		// Inform all neighbors in the in-edges.
		for (size_t i = 0; i < data->neighbors->size(); i++) {
			size_t count = neighbors->get_count(i);
			if (count > 0) {
				count_msg msg(count);
				prog.send_msg(data->neighbors->get_neighbor_id(i), msg);
			}
		}
	}

	size_t get_result() const {
		return get_local_scan();
	}
};

runtime_data_t *create_runtime(graph_engine &graph, local_scan_vertex &scan_v,
		const page_vertex &vertex)
{
	class skip_larger {
		vsize_t degree;
		vertex_id_t id;
		graph_engine &graph;

		vsize_t get_degree(vertex_id_t id) const {
			return ((local_scan_vertex &) graph.get_vertex(id)).get_degree();
		}
	public:
		skip_larger(graph_engine &_graph, vertex_id_t id): graph(_graph) {
			this->degree = get_degree(id);
			this->id = id;
		}

		bool operator()(attributed_neighbor &e) {
			return operator()(e.get_id());
		}

		/**
		 * We are going to count edges on the vertices with the most edges.
		 * If two vertices have the same number of edges, we compute
		 * on the vertices with the largest Id.
		 */
		bool operator()(vertex_id_t id) {
			vsize_t other_degree = get_degree(id);
			if (other_degree == degree)
				return id >= this->id;
			return other_degree > degree;
		}
	};

	class merge_edge
	{
	public:
		attributed_neighbor operator()(const attributed_neighbor &e1,
				const attributed_neighbor &e2) {
			assert(e1.get_id() == e2.get_id());
			return attributed_neighbor(e1.get_id(),
					e1.get_num_dups() + e2.get_num_dups());
		}
	};

	merge_edge merge;
	std::vector<attributed_neighbor> neighbors(
			vertex.get_num_edges(edge_type::BOTH_EDGES));
	size_t num_neighbors = unique_merge(
			vertex.get_neigh_begin(edge_type::IN_EDGE),
			vertex.get_neigh_end(edge_type::IN_EDGE),
			vertex.get_neigh_begin(edge_type::OUT_EDGE),
			vertex.get_neigh_end(edge_type::OUT_EDGE),
			skip_larger(graph, vertex.get_id()), merge,
			neighbors.begin());
	neighbors.resize(num_neighbors);
	return new runtime_data_t(std::shared_ptr<neighbor_list>(
				new extended_neighbor_list(vertex, neighbors)),
			scan_v.get_local_scan());
}

void destroy_runtime(runtime_data_t *data)
{
	delete data;
}

void local_scan_vertex::run_on_itself(vertex_program &prog, const page_vertex &vertex)
{
	assert(!local_value.has_runtime_data());

	size_t num_local_edges = vertex.get_num_edges(edge_type::BOTH_EDGES);
	if (num_local_edges == 0)
		return;

	runtime_data_t *local_data = create_runtime(prog.get_graph(), *this, vertex);
	local_value.set_runtime_data(local_data);

	size_t tmp = 0;
	edge_iterator it = vertex.get_neigh_begin(edge_type::IN_EDGE);
	edge_iterator end = vertex.get_neigh_end(edge_type::IN_EDGE);
	for (; it != end; ++it) {
		vertex_id_t id = *it;
		// Ignore loops
		if (id != vertex.get_id())
			tmp++;
	}
	it = vertex.get_neigh_begin(edge_type::OUT_EDGE);
	end = vertex.get_neigh_end(edge_type::OUT_EDGE);
	for (; it != end; ++it) {
		vertex_id_t id = *it;
		// Ignore loops
		if (id != vertex.get_id())
			tmp++;
	}

	local_data->local_scan += tmp;

	if (local_data->neighbors->empty()) {
		local_value.set_real_local(local_data->local_scan);
		finding_triangles_end(prog, local_data);
		destroy_runtime(local_data);
		return;
	}

	std::vector<vertex_id_t> neighbors;
	local_data->neighbors->get_neighbors(neighbors);
	request_vertices(neighbors.data(), neighbors.size());
}

void local_scan_vertex::run_on_neighbor(vertex_program &prog, const page_vertex &vertex)
{
	assert(local_value.has_runtime_data());
	runtime_data_t *local_data = local_value.get_runtime_data();
	local_data->num_joined++;
	size_t ret = local_data->neighbors->count_edges(&vertex);
	if (ret > 0)
		local_data->local_scan += ret;

	// If we have seen all required neighbors, we have complete
	// the computation. We can release the memory now.
	if (local_data->num_joined == local_data->neighbors->size()) {
		local_value.set_real_local(local_data->local_scan);

		finding_triangles_end(prog, local_data);

		destroy_runtime(local_data);
	}
}

}

#include "save_result.h"

namespace fg
{

FG_vector<size_t>::ptr compute_local_scan(FG_graph::ptr fg)
{
	bool directed = fg->get_graph_header().is_directed_graph();
	if (!directed) {
		BOOST_LOG_TRIVIAL(error)
			<< "This algorithm current works on a directed graph";
		return FG_vector<size_t>::ptr();
	}

	graph_index::ptr index = NUMA_graph_index<local_scan_vertex>::create(
			fg->get_graph_header());
	graph_engine::ptr graph = fg->create_engine(index);

	BOOST_LOG_TRIVIAL(info) << "local scan starts";
	BOOST_LOG_TRIVIAL(info) << "prof_file: " << graph_conf.get_prof_file();
#ifdef PROFILER
	if (!graph_conf.get_prof_file().empty())
		ProfilerStart(graph_conf.get_prof_file().c_str());
#endif

	struct timeval start, end;
	gettimeofday(&start, NULL);
	scan_stage = scan_stage_t::INIT;
	graph->start_all();
	graph->wait4complete();

	scan_stage = scan_stage_t::RUN;
	graph->start_all();
	graph->wait4complete();
	gettimeofday(&end, NULL);

#ifdef PROFILER
	if (!graph_conf.get_prof_file().empty())
		ProfilerStop();
#endif
	BOOST_LOG_TRIVIAL(info)
		<< boost::format("It takes %1% seconds to compute all local scan")
		% time_diff(start, end);

	FG_vector<size_t>::ptr vec = FG_vector<size_t>::create(graph);
	graph->query_on_all(vertex_query::ptr(
				new save_query<size_t, local_scan_vertex>(vec)));
	return vec;
}

}
