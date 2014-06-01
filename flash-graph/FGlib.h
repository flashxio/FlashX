#ifndef __FGLIB_H__
#define __FGLIB_H__

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
#include "graph.h"
#include "FG_vector.h"

class FG_graph
{
	std::string graph_file;
	std::string index_file;
	config_map configs;
	FG_graph(const std::string &graph_file,
			const std::string &index_file, const config_map &configs) {
		this->graph_file = graph_file;
		this->index_file = index_file;
		this->configs = configs;
	}
public:
	typedef std::shared_ptr<FG_graph> ptr;

	static ptr create(const std::string &graph_file,
			const std::string &index_file, const config_map &configs) {
		return ptr(new FG_graph(graph_file, index_file, configs));
	}

	const std::string &get_graph_file() const {
		return graph_file;
	}

	const std::string &get_index_file() const {
		return index_file;
	}

	const config_map &get_configs() const {
		return configs;
	}
};

enum directed_triangle_type
{
	CYCLE,
	ALL,
};

FG_vector<vertex_id_t>::ptr compute_wcc(FG_graph::ptr fg);
FG_vector<vertex_id_t>::ptr compute_scc(FG_graph::ptr fg);
FG_vector<size_t>::ptr compute_directed_triangles(FG_graph::ptr fg,
		directed_triangle_type type);
FG_vector<size_t>::ptr compute_undirected_triangles(FG_graph::ptr fg);
FG_vector<size_t>::ptr compute_local_scan(FG_graph::ptr);
FG_vector<size_t>::ptr compute_topK_scan(FG_graph::ptr, size_t topK);

/**
 * Fetch the clusters with the wanted cluster IDs.
 */
void fetch_subgraphs(FG_graph::ptr graph, FG_vector<vertex_id_t>::ptr comp_ids,
		const std::set<vertex_id_t> &wanted_clusters, std::map<vertex_id_t,
		graph::ptr> &clusters);

/**
 * Compute the size of each subgraph identified by cluster IDs.
 */
void compute_subgraph_sizes(FG_graph::ptr graph, FG_vector<vertex_id_t>::ptr comp_ids,
		const std::set<vertex_id_t> &wanted_clusters,
		std::map<vertex_id_t, std::pair<size_t, size_t> > &sizes);

#endif
