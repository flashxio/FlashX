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

#include "FGlib.h"
#include "vertex.h"

/******************* Implementation of fetching clusters *********************/

namespace {

class cluster_subgraph_vertex: public compute_vertex
{
public:
	cluster_subgraph_vertex(vertex_id_t id): compute_vertex(id) {
	}

	void run(vertex_program &prog);
	void run(vertex_program &prog, const page_vertex &vertex);
	void run_on_message(vertex_program &prog, const vertex_message &msg) {
	}
};

void create_clusters(const std::set<vertex_id_t> &ids, graph_type type,
		std::map<vertex_id_t, graph::ptr> &clusters)
{
	BOOST_FOREACH(vertex_id_t id, ids) {
		graph::ptr g;
		switch(type) {
			case DIRECTED:
				g = directed_graph<>::create(false);
				break;
			case UNDIRECTED:
				g = undirected_graph<>::create(false);
				break;
			default:
				assert(0);
		}
		clusters.insert(std::pair<vertex_id_t, graph::ptr>(
					id, g));
	}
}

class cluster_subgraph_vertex_program: public vertex_program_impl<cluster_subgraph_vertex>
{
	graph_type type;
	FG_vector<vertex_id_t>::ptr cluster_ids;
	std::map<vertex_id_t, graph::ptr> wanted_clusters;
public:
	typedef std::shared_ptr<cluster_subgraph_vertex_program> ptr;

	static ptr cast2(vertex_program::ptr prog) {
		return std::static_pointer_cast<cluster_subgraph_vertex_program,
			   vertex_program>(prog);
	}

	cluster_subgraph_vertex_program(FG_vector<vertex_id_t>::ptr cluster_ids,
			const std::set<vertex_id_t> &wanted_clusters, graph_type type) {
		this->type = type;
		this->cluster_ids = cluster_ids;
		create_clusters(wanted_clusters, type, this->wanted_clusters);
	}

	void add_vertex(const page_vertex &vertex);

	const std::map<vertex_id_t, graph::ptr> &get_clusters() const {
		return wanted_clusters;
	}

	bool is_wanted_vertex(vertex_id_t vid) const {
		vertex_id_t cluster_id = cluster_ids->get(vid);
		std::map<vertex_id_t, graph::ptr>::const_iterator it
			= wanted_clusters.find(cluster_id);
		return it != wanted_clusters.end();
	}
};

class cluster_subgraph_vertex_program_creater: public vertex_program_creater
{
	FG_vector<vertex_id_t>::ptr cluster_ids;
	std::set<vertex_id_t> wanted_clusters;
	graph_type type;
public:
	cluster_subgraph_vertex_program_creater(
			FG_vector<vertex_id_t>::ptr cluster_ids,
			const std::set<vertex_id_t> &wanted_clusters,
			graph_type type) {
		this->cluster_ids = cluster_ids;
		this->wanted_clusters = wanted_clusters;
		this->type = type;
	}

	vertex_program::ptr create() const {
		return vertex_program::ptr(new cluster_subgraph_vertex_program(
					cluster_ids, wanted_clusters, type));
	}
};

void cluster_subgraph_vertex::run(vertex_program &prog)
{
	if (((cluster_subgraph_vertex_program &) prog).is_wanted_vertex(get_id())) {
		vertex_id_t id = get_id();
		request_vertices(&id, 1);
	}
}

void cluster_subgraph_vertex::run(vertex_program &prog, const page_vertex &vertex)
{
	((cluster_subgraph_vertex_program &) prog).add_vertex(vertex);
}

void cluster_subgraph_vertex_program::add_vertex(const page_vertex &vertex)
{
	vertex_id_t cluster_id = cluster_ids->get(vertex.get_id());
	std::map<vertex_id_t, graph::ptr>::iterator it = wanted_clusters.find(
			cluster_id);
	assert(it != wanted_clusters.end());

	graph::ptr g = it->second;
	if (type == UNDIRECTED) {
		in_mem_undirected_vertex<> in_v((const page_undirected_vertex &) vertex,
				false);
		g->add_vertex(in_v);
	}
	else if (type == DIRECTED) {
		in_mem_directed_vertex<> in_v((const page_directed_vertex &) vertex,
				false);
		g->add_vertex(in_v);
	}
	else
		assert(0);
}

}

void fetch_subgraphs(FG_graph::ptr fg, FG_vector<vertex_id_t>::ptr cluster_ids,
		const std::set<vertex_id_t> &wanted_clusters,
		std::map<vertex_id_t, graph::ptr> &clusters)
{
	graph_index::ptr index = NUMA_graph_index<cluster_subgraph_vertex>::create(
			fg->get_index_file());
	graph_engine::ptr graph = graph_engine::create(fg->get_graph_file(),
			index, fg->get_configs());

	graph_type type = graph->get_graph_header().get_graph_type();
	struct timeval start, end;
	gettimeofday(&start, NULL);
	graph->start_all(vertex_initializer::ptr(), vertex_program_creater::ptr(
				new cluster_subgraph_vertex_program_creater(cluster_ids,
					wanted_clusters, type)));
	graph->wait4complete();
	gettimeofday(&end, NULL);
	printf("fetching subgraphs takes %f seconds\n", time_diff(start, end));

	create_clusters(wanted_clusters, type, clusters);
	std::vector<vertex_program::ptr> vprogs;
	graph->get_vertex_programs(vprogs);
	BOOST_FOREACH(vertex_program::ptr vprog, vprogs) {
		cluster_subgraph_vertex_program::ptr subgraph_vprog
			= cluster_subgraph_vertex_program::cast2(vprog);
		typedef std::pair<vertex_id_t, graph::ptr> pair_t;
		BOOST_FOREACH(pair_t v, subgraph_vprog->get_clusters()) {
			vertex_id_t cluster_id = v.first;
			graph::ptr g = v.second;

			std::map<vertex_id_t, graph::ptr>::iterator it = clusters.find(
					cluster_id);
			assert(it != clusters.end());
			it->second->merge(g);
		}
	}
}

/**************** Implementation of computing cluster sizes ******************/

namespace {

typedef std::set<vertex_id_t> v_set_t;
typedef std::map<vertex_id_t, v_set_t> cluster_map_t;
typedef std::map<vertex_id_t, std::pair<size_t, size_t> > size_map_t;

class subgraph_size_vertex: public compute_vertex
{
public:
	subgraph_size_vertex(vertex_id_t id): compute_vertex(id) {
	}

	void run(vertex_program &prog);
	void run(vertex_program &prog, const page_vertex &vertex);
	void run_on_message(vertex_program &prog, const vertex_message &msg) {
	}
};

class subgraph_size_vertex_program: public vertex_program_impl<subgraph_size_vertex>
{
	FG_vector<vertex_id_t>::ptr cluster_ids;
	const cluster_map_t &wanted_clusters;
	size_map_t wanted_cluster_sizes;
public:
	typedef std::shared_ptr<subgraph_size_vertex_program> ptr;

	static ptr cast2(vertex_program::ptr prog) {
		return std::static_pointer_cast<subgraph_size_vertex_program,
			   vertex_program>(prog);
	}

	subgraph_size_vertex_program(FG_vector<vertex_id_t>::ptr cluster_ids,
			const cluster_map_t &_wanted_clusters): wanted_clusters(
				_wanted_clusters) {
		this->cluster_ids = cluster_ids;
		BOOST_FOREACH(cluster_map_t::value_type v, wanted_clusters) {
			vertex_id_t cluster_id = v.first;
			wanted_cluster_sizes.insert(size_map_t::value_type(cluster_id,
						std::pair<size_t, size_t>(0, 0)));
		}
	}

	void add_vertex(const page_vertex &vertex);

	bool is_wanted_vertex(vertex_id_t vid) const {
		vertex_id_t cluster_id = cluster_ids->get(vid);
		cluster_map_t::const_iterator it
			= wanted_clusters.find(cluster_id);
		return it != wanted_clusters.end();
	}

	const size_map_t &get_cluster_sizes() const {
		return wanted_cluster_sizes;
	}
};

void subgraph_size_vertex::run(vertex_program &prog)
{
	if (((subgraph_size_vertex_program &) prog).is_wanted_vertex(get_id())) {
		vertex_id_t id = get_id();
		request_vertices(&id, 1);
	}
}

void subgraph_size_vertex::run(vertex_program &prog, const page_vertex &vertex)
{
	((subgraph_size_vertex_program &) prog).add_vertex(vertex);
}

void subgraph_size_vertex_program::add_vertex(const page_vertex &vertex)
{
	vertex_id_t cluster_id = cluster_ids->get(vertex.get_id());
	cluster_map_t::const_iterator it = wanted_clusters.find(
			cluster_id);
	// If this vertex doesn't belong to a cluster we want.
	if (it == wanted_clusters.end())
		return;

	const v_set_t &cluster = it->second;
	vsize_t num_edges = vertex.get_num_edges(edge_type::BOTH_EDGES);
	// The number of edges that are inside the cluster.
	vsize_t num_inside_edges = 0;
	edge_seq_iterator edge_it = vertex.get_neigh_seq_it(
			edge_type::BOTH_EDGES, 0, num_edges);
	PAGE_FOREACH(vertex_id_t, id, edge_it) {
		v_set_t::const_iterator it = cluster.find(id);
		if (it != cluster.end())
			num_inside_edges++;
	} PAGE_FOREACH_END

	size_map_t::iterator size_it = wanted_cluster_sizes.find(cluster_id);
	assert(size_it != wanted_cluster_sizes.end());
	// Increase the number of vertices by 1.
	size_it->second.first++;
	// Increase the number of edges.
	size_it->second.second += num_inside_edges;
}

void merge_clusters(const cluster_map_t &clusters, cluster_map_t &aggs)
{
	BOOST_FOREACH(cluster_map_t::value_type v, clusters) {
		vertex_id_t cluster_id = v.first;

		cluster_map_t::iterator it = aggs.find(cluster_id);
		if (it == aggs.end())
			aggs.insert(cluster_map_t::value_type(cluster_id, v.second));
		else
			it->second.insert(v.second.begin(), v.second.end());
	}
}

class cluster_query: public vertex_query
{
	cluster_map_t clusters;
	FG_vector<vertex_id_t>::ptr cluster_ids;
	const v_set_t &wanted_cluster_ids;
public:
	cluster_query(FG_vector<vertex_id_t>::ptr cluster_ids,
			const v_set_t &_wanted_cluster_ids): wanted_cluster_ids(
				_wanted_cluster_ids) {
		this->cluster_ids = cluster_ids;
	}

	virtual void run(graph_engine &graph, compute_vertex &v) {
		vertex_id_t cluster_id = cluster_ids->get(v.get_id());
		v_set_t::const_iterator v_it = wanted_cluster_ids.find(cluster_id);
		// If the vertex isn't in the cluster we want.
		if (v_it == wanted_cluster_ids.end())
			return;

		cluster_map_t::iterator it = clusters.find(cluster_id);
		if (it == clusters.end()) {
			v_set_t v_set;
			v_set.insert(v.get_id());
			clusters.insert(cluster_map_t::value_type(cluster_id, v_set));
		}
		else {
			it->second.insert(v.get_id());
		}
	}

	virtual void merge(graph_engine &graph, vertex_query::ptr q) {
		cluster_query &cq = (cluster_query &) *q;
		merge_clusters(cq.clusters, this->clusters);
	}

	virtual ptr clone() {
		return vertex_query::ptr(new cluster_query(cluster_ids,
					wanted_cluster_ids));
	}

	const cluster_map_t &get_clusters() const {
		return clusters;
	}
};

void merge_clusters(const size_map_t &cluster_sizes, size_map_t &agg_sizes)
{
	BOOST_FOREACH(size_map_t::value_type v, cluster_sizes) {
		vertex_id_t cluster_id = v.first;
		size_t num_v = v.second.first;
		size_t num_e = v.second.second;

		size_map_t::iterator it = agg_sizes.find(cluster_id);
		if (it == agg_sizes.end()) {
			agg_sizes.insert(size_map_t::value_type(cluster_id,
						std::pair<size_t, size_t>(num_v, num_e)));
		}
		else {
			it->second.first += num_v;
			it->second.second += num_e;
		}
	}
}

class subgraph_size_vertex_program_creater: public vertex_program_creater
{
	FG_vector<vertex_id_t>::ptr cluster_ids;
	const cluster_map_t &wanted_clusters;
public:
	subgraph_size_vertex_program_creater(
			FG_vector<vertex_id_t>::ptr cluster_ids,
			const cluster_map_t &_wanted_clusters): wanted_clusters(
				_wanted_clusters) {
		this->cluster_ids = cluster_ids;
	}

	vertex_program::ptr create() const {
		return vertex_program::ptr(new subgraph_size_vertex_program(
					cluster_ids, wanted_clusters));
	}
};

}

void compute_subgraph_sizes(FG_graph::ptr fg,
		FG_vector<vertex_id_t>::ptr cluster_ids,
		const std::set<vertex_id_t> &wanted_cluster_ids, size_map_t &cluster_sizes)
{
	graph_index::ptr index = NUMA_graph_index<subgraph_size_vertex>::create(
			fg->get_index_file());
	graph_engine::ptr graph = graph_engine::create(fg->get_graph_file(),
			index, fg->get_configs());

	vertex_query::ptr q(new cluster_query(cluster_ids, wanted_cluster_ids));
	graph->query_on_all(q);
	const cluster_map_t &wanted_clusters = ((cluster_query &) *q).get_clusters();

	struct timeval start, end;
	gettimeofday(&start, NULL);
	graph->start_all(vertex_initializer::ptr(), vertex_program_creater::ptr(
				new subgraph_size_vertex_program_creater(cluster_ids,
					wanted_clusters)));
	graph->wait4complete();
	gettimeofday(&end, NULL);
	printf("computing subgraph size takes %f seconds\n", time_diff(start, end));

	std::vector<vertex_program::ptr> vprogs;
	graph->get_vertex_programs(vprogs);
	BOOST_FOREACH(vertex_program::ptr vprog, vprogs) {
		subgraph_size_vertex_program::ptr subgraph_vprog
			= subgraph_size_vertex_program::cast2(vprog);
		merge_clusters(subgraph_vprog->get_cluster_sizes(), cluster_sizes);
	}

	BOOST_FOREACH(size_map_t::value_type &v, cluster_sizes) {
		assert(v.second.second % 2 == 0);
		v.second.second /= 2;
	}
}

/********************** Get the degree of vertices ****************************/

namespace {

class degree_vertex: public compute_vertex
{
public:
	degree_vertex(vertex_id_t id): compute_vertex(id) {
	}

	virtual void run(vertex_program &prog) {
		vertex_id_t id = get_id();
		request_vertex_headers(&id, 1);
	}

	void run(vertex_program &prog, const page_vertex &vertex) {
	}

	void run_on_message(vertex_program &, const vertex_message &msg) {
	}

	void run_on_vertex_header(vertex_program &prog,
			const vertex_header &header);
};

class degree_vertex_program: public vertex_program_impl<degree_vertex>
{
	edge_type type;
	FG_vector<vsize_t>::ptr degree_vec;
public:
	degree_vertex_program(FG_vector<vsize_t>::ptr degree_vec, edge_type type) {
		this->degree_vec = degree_vec;
		this->type = type;
	}

	edge_type get_edge_type() const {
		return type;
	}

	void set_degree(vertex_id_t id, vsize_t degree) {
		degree_vec->set(id, degree);
	}
};

class degree_vertex_program_creater: public vertex_program_creater
{
	FG_vector<vertex_id_t>::ptr degree_vec;
	edge_type type;
public:
	degree_vertex_program_creater(
			FG_vector<vertex_id_t>::ptr degree_vec,
			edge_type type) {
		this->degree_vec = degree_vec;
		this->type = type;
	}

	vertex_program::ptr create() const {
		return vertex_program::ptr(new degree_vertex_program(
					degree_vec, type));
	}
};

void degree_vertex::run_on_vertex_header(vertex_program &prog,
		const vertex_header &header)
{
	degree_vertex_program &degree_vprog = (degree_vertex_program &) prog;
	if (prog.get_graph().is_directed()) {
		const directed_vertex_header &dheader
			= (const directed_vertex_header &) header;
		if (degree_vprog.get_edge_type() == IN_EDGE)
			degree_vprog.set_degree(get_id(), dheader.get_num_in_edges());
		else if (degree_vprog.get_edge_type() == OUT_EDGE)
			degree_vprog.set_degree(get_id(), dheader.get_num_out_edges());
		else
			degree_vprog.set_degree(get_id(), dheader.get_num_edges());
	}
	else
		degree_vprog.set_degree(get_id(), header.get_num_edges());
}

}

FG_vector<vsize_t>::ptr get_degree(FG_graph::ptr fg, edge_type type)
{
	graph_index::ptr index = NUMA_graph_index<degree_vertex>::create(
			fg->get_index_file());
	graph_engine::ptr graph = graph_engine::create(fg->get_graph_file(),
			index, fg->get_configs());

	FG_vector<vsize_t>::ptr degree_vec = FG_vector<vsize_t>::create(graph);
	struct timeval start, end;
	gettimeofday(&start, NULL);
	graph->start_all(vertex_initializer::ptr(), vertex_program_creater::ptr(
				new degree_vertex_program_creater(degree_vec, type)));
	graph->wait4complete();
	gettimeofday(&end, NULL);
	printf("computing subgraph size takes %f seconds\n", time_diff(start, end));
	return degree_vec;
}

/*************** Get the degree of vertices in a timestamp ********************/

namespace {

class ts_degree_vertex: public compute_vertex
{
public:
	ts_degree_vertex(vertex_id_t id): compute_vertex(id) {
	}

	virtual void run(vertex_program &prog) {
		vertex_id_t id = get_id();
		request_vertices(&id, 1);
	}

	void run(vertex_program &prog, const page_vertex &vertex);

	void run_on_message(vertex_program &, const vertex_message &msg) {
	}
};

class ts_degree_vertex_program: public vertex_program_impl<ts_degree_vertex>
{
	time_t start_time;
	time_t time_interval;
	edge_type type;
	FG_vector<vsize_t>::ptr degree_vec;
public:
	ts_degree_vertex_program(FG_vector<vsize_t>::ptr degree_vec, edge_type type,
			time_t start_time, time_t time_interval) {
		this->degree_vec = degree_vec;
		this->type = type;
		this->start_time = start_time;
		this->time_interval = time_interval;
	}

	time_t get_start_time() const {
		return start_time;
	}

	time_t get_time_interval() const {
		return time_interval;
	}

	edge_type get_edge_type() const {
		return type;
	}

	void set_degree(vertex_id_t id, vsize_t degree) {
		degree_vec->set(id, degree);
	}
};

class ts_degree_vertex_program_creater: public vertex_program_creater
{
	time_t start_time;
	time_t time_interval;
	FG_vector<vertex_id_t>::ptr degree_vec;
	edge_type type;
public:
	ts_degree_vertex_program_creater(
			FG_vector<vertex_id_t>::ptr degree_vec, edge_type type,
			time_t start_time, time_t time_interval) {
		this->degree_vec = degree_vec;
		this->type = type;
		this->start_time = start_time;
		this->time_interval = time_interval;
	}

	vertex_program::ptr create() const {
		return vertex_program::ptr(new ts_degree_vertex_program(
					degree_vec, type, start_time, time_interval));
	}
};

void ts_degree_vertex::run(vertex_program &prog, const page_vertex &vertex)
{
	ts_degree_vertex_program &degree_vprog = (ts_degree_vertex_program &) prog;
	edge_type type = degree_vprog.get_edge_type();
	time_t start_time = degree_vprog.get_start_time();
	time_t time_interval = degree_vprog.get_time_interval();

	if (prog.get_graph().is_directed()) {
		const page_directed_vertex &dv = (const page_directed_vertex &) vertex;

		page_byte_array::const_iterator<ts_edge_data> begin_it
			= dv.get_data_begin<ts_edge_data>(type);
		page_byte_array::const_iterator<ts_edge_data> end_it
			= dv.get_data_end<ts_edge_data>(type);
		page_byte_array::const_iterator<ts_edge_data> ts_it = std::lower_bound(
				begin_it, end_it, ts_edge_data(start_time));
		page_byte_array::const_iterator<ts_edge_data> ts_end_it = std::lower_bound(
				begin_it, end_it, start_time + time_interval);
		degree_vprog.set_degree(vertex.get_id(), ts_end_it - ts_it);
	}
	else {
		assert(0);
	}
}

}

FG_vector<vsize_t>::ptr get_ts_degree(FG_graph::ptr fg, edge_type type,
		time_t start_time, time_t time_interval)
{
	graph_index::ptr index = NUMA_graph_index<ts_degree_vertex>::create(
			fg->get_index_file());
	graph_engine::ptr graph = graph_engine::create(fg->get_graph_file(),
			index, fg->get_configs());
	assert(graph->get_graph_header().has_edge_data());

	FG_vector<vsize_t>::ptr degree_vec = FG_vector<vsize_t>::create(graph);
	graph->start_all(vertex_initializer::ptr(), vertex_program_creater::ptr(
				new ts_degree_vertex_program_creater(degree_vec, type,
					start_time, time_interval)));
	graph->wait4complete();
	return degree_vec;
}

graph_header get_graph_header(FG_graph::ptr fg)
{
	graph_header header;
	file_io_factory::shared_ptr index_factory = create_io_factory(
			fg->get_index_file(), GLOBAL_CACHE_ACCESS);
	io_interface::ptr io = index_factory->create_io(thread::get_curr_thread());
	io->access((char *) &header, 0, sizeof(header), READ);
	return header;
}
