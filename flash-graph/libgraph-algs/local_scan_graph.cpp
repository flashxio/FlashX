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
#include <google/profiler.h>
#endif

#include "FG_vector.h"
#include "FGlib.h"

#include "scan_graph.h"

namespace {

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

	off_t find_idx(vertex_id_t id) {
		edge_set_t::const_iterator it = neighbor_set->find(id);
		if (it == neighbor_set->end())
			return -1;
		else {
			size_t idx = (*it).get_idx();
			assert(idx < id_list.size());
			return idx;
		}
	}

	attributed_neighbor find(vertex_id_t id) {
		edge_set_t::const_iterator it = neighbor_set->find(id);
		if (it == neighbor_set->end())
			return attributed_neighbor();
		else {
			size_t idx = (*it).get_idx();
			assert(idx < id_list.size());
			return at(idx);
		}
	}

	attributed_neighbor at(size_t idx) {
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

#ifdef PV_STAT
	rand_jumps += common_neighs.size() + 1;
#endif
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

class local_scan_vertex: public scan_vertex
{
public:
	local_scan_vertex(vertex_id_t id): scan_vertex(id) {
		local_value.set_real_local(0);
	}

	using scan_vertex::run;

	void run(vertex_program &prog) {
		vertex_id_t id = get_id();
		request_vertices(&id, 1);
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

	void run_on_message(vertex_program &prog, const vertex_message &msg1) {
		const count_msg &msg = (const count_msg &) msg1;
		if (local_value.has_runtime_data())
			local_value.get_runtime_data()->local_scan += msg.get_num();
		else
			local_value.inc_real_local(msg.get_num());
	}
};

class skip_larger {
	vsize_t degree;
	vertex_id_t id;
	graph_engine &graph;

	vsize_t get_degree(vertex_id_t id) const {
		return ((scan_vertex &) graph.get_vertex(id)).get_degree();
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

runtime_data_t *ls_create_runtime(graph_engine &graph, scan_vertex &scan_v,
		const page_vertex &vertex)
{
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
	return new runtime_data_t(std::unique_ptr<neighbor_list>(
				new extended_neighbor_list(vertex, neighbors)),
			scan_v.get_local_scan());
}

void ls_finding_triangles_end(vertex_program &prog, scan_vertex &scan_v,
		runtime_data_t *data)
{
	local_scan_vertex &ls_v = (local_scan_vertex &) scan_v;
	ls_v.finding_triangles_end(prog, data);
}

class save_scan_query: public vertex_query
{
	FG_vector<size_t>::ptr vec;
public:
	save_scan_query(FG_vector<size_t>::ptr vec) {
		this->vec = vec;
	}

	virtual void run(graph_engine &graph, compute_vertex &v) {
		local_scan_vertex &lv = (local_scan_vertex &) v;
		vec->set(lv.get_id(), lv.get_local_scan());
	}

	virtual void merge(graph_engine &graph, vertex_query::ptr q) {
	}

	virtual ptr clone() {
		return vertex_query::ptr(new save_scan_query(vec));
	}
};

}

FG_vector<size_t>::ptr compute_local_scan(FG_graph::ptr fg)
{
	finding_triangles_end = ls_finding_triangles_end;
	create_runtime = ls_create_runtime;

	graph_index::ptr index = NUMA_graph_index<local_scan_vertex>::create(
			fg->get_index_file());
	graph_engine::ptr graph = graph_engine::create(fg->get_graph_file(),
			index, fg->get_configs());

	printf("local scan starts\n");
	printf("prof_file: %s\n", graph_conf.get_prof_file().c_str());
#ifdef PROFILER
	if (!graph_conf.get_prof_file().empty())
		ProfilerStart(graph_conf.get_prof_file().c_str());
#endif

	struct timeval start, end;
	gettimeofday(&start, NULL);
	graph->start_all();
	graph->wait4complete();
	gettimeofday(&end, NULL);

#ifdef PROFILER
	if (!graph_conf.get_prof_file().empty())
		ProfilerStop();
#endif
	printf("It takes %f seconds to compute all local scan\n",
			time_diff(start, end));

	FG_vector<size_t>::ptr vec = FG_vector<size_t>::create(graph);
	graph->query_on_all(vertex_query::ptr(new save_scan_query(vec)));
	return vec;
}
