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

class cluster_subgraph_vertex: public compute_vertex
{
public:
	cluster_subgraph_vertex() {
	}

	cluster_subgraph_vertex(vertex_id_t id,
			const vertex_index &index): compute_vertex(id, index) {
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
	graph->start_all(vertex_initiator::ptr(), vertex_program_creater::ptr(
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
