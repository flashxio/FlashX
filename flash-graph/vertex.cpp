/*
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

#include "vertex.h"
#include "vertex_index.h"

namespace fg
{

empty_data edge<empty_data>::data;

size_t ext_mem_undirected_vertex::serialize(const in_mem_vertex &v, char *buf,
		size_t size, edge_type type)
{
	assert(ext_mem_undirected_vertex::get_header_size() <= size);
	ext_mem_undirected_vertex *ext_v = (ext_mem_undirected_vertex *) buf;
	ext_v->set_id(v.get_id());
	ext_v->num_edges = v.get_num_edges(type);
	if (v.has_edge_data())
		ext_v->edge_data_size = v.get_edge_data_size();
	else
		ext_v->edge_data_size = 0;
	size_t mem_size = ext_v->get_size();
	assert(mem_size <= MAX_VERTEX_SIZE);
	assert(size >= mem_size);
	v.serialize_edges(ext_v->neighbors, type);
	// serialize edge data
	if (v.has_edge_data()) {
		v.serialize_edge_data(ext_v->get_edge_data_addr(), type);
	}

	return mem_size;
}

}
