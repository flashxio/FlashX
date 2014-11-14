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

#include <boost/foreach.hpp>

#include "graph.h"
#include "in_mem_storage.h"
#include "utils.h"

/*
 * This class is to compress the vertex ID space in a graph.
 */
class graph_compressor
{
	// Old vertex Id <-> new vertex Id
	std::unordered_map<vertex_id_t, vertex_id_t> map;
public:
	graph_compressor(const std::vector<vertex_id_t> &vids) {
		assert(std::is_sorted(vids.begin(), vids.end()));
		for (vertex_id_t new_id = 0; new_id < vids.size(); new_id++)
			map.insert(std::pair<vertex_id_t, vertex_id_t>(
						vids[new_id], new_id));
	}

	std::shared_ptr<in_mem_vertex> compress(const in_mem_vertex &v) const {
		return v.create_remapped_vertex(map);
	}
};

std::pair<in_mem_graph::ptr, vertex_index::ptr> in_mem_subgraph::compress(
		const std::string &name) const
{
	std::vector<vertex_id_t> vertex_ids;
	get_all_vertices(vertex_ids);
	std::sort(vertex_ids.begin(), vertex_ids.end());
	if (vertex_ids.empty())
		return std::pair<in_mem_graph::ptr, vertex_index::ptr>(
				in_mem_graph::ptr(), vertex_index::ptr());

	graph_compressor compressor(vertex_ids);
	size_t edge_data_size = get_vertex(vertex_ids.front()).get_edge_data_size();
	mem_serial_graph::ptr serial_g = mem_serial_graph::create(is_directed(),
			edge_data_size);
	BOOST_FOREACH(vertex_id_t id, vertex_ids) {
		const in_mem_vertex &v = get_vertex(id);
		assert(v.get_edge_data_size() == edge_data_size);
		serial_g->add_vertex(*compressor.compress(v));
	}
	return std::pair<in_mem_graph::ptr, vertex_index::ptr>(
			serial_g->dump_graph(name), serial_g->dump_index(true));
}

in_mem_subgraph::ptr in_mem_subgraph::create(graph_type type, bool has_data)
{
	if (type == graph_type::DIRECTED)
		return in_mem_directed_subgraph<>::create(has_data);
	else if (type == graph_type::UNDIRECTED)
		return in_mem_undirected_subgraph<>::create(has_data);
	else
		ABORT_MSG("wrong graph type");
}
