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

#include "ext_mem_vertex_iterator.h"

#if 0
void ext_mem_vertex_iterator::read_vertices()
{
	std::vector<in_mem_vertex_info> vinfos = overflow_vinfos;
	overflow_vinfos.clear();

	// Find the number of vertices that can be stored in the allocated buffer.
	size_t read_size = 0;
	for (size_t i = 0; i < vinfos.size(); i++)
		read_size += vinfos[i].get_ext_mem_size();
	assert(read_size <= MEM_SIZE);
	while (index_it->has_next()) {
		in_mem_vertex_info info = index_it->next();
		if (read_size + info.get_ext_mem_size() < MEM_SIZE) {
			vinfos.push_back(info);
			read_size += info.get_ext_mem_size();
		}
		else {
			overflow_vinfos.push_back(info);
			break;
		}
	}
	assert(!vinfos.empty());
	assert((size_t) vinfos.back().get_ext_mem_off()
			- vinfos.front().get_ext_mem_off()
			+ vinfos.back().get_ext_mem_size() == read_size);
	assert(ftell(adj_f) == vinfos.front().get_ext_mem_off());

	size_t ret = fread(adj_buf, read_size, 1, adj_f);
	assert(ret == 1);

	graph_type type = header.get_graph_type();
	vertices.clear();
	off_t start_off = vinfos.front().get_ext_mem_off();
	for (size_t i = 0; i < vinfos.size(); i++) {
		size_t size = vinfos[i].get_ext_mem_size();
		off_t off = vinfos[i].get_ext_mem_off() - start_off;
		vertices.push_back(adj_buf + off);

		size_t vsize;
		if (type == graph_type::DIRECTED) {
			ext_mem_directed_vertex *v
				= (ext_mem_directed_vertex *) vertices.back();
			vsize = v->get_size();
		}
		else if (type == graph_type::UNDIRECTED) {
			ext_mem_undirected_vertex *v
				= (ext_mem_undirected_vertex *) vertices.back();
			vsize = v->get_size();
		}
		else if (type == graph_type::TS_DIRECTED) {
			ts_ext_mem_directed_vertex *v
				= (ts_ext_mem_directed_vertex *) vertices.back();
			vsize = v->get_size();
		}
		else {
			assert(0);
		}
		assert(off + vsize <= read_size);
		assert(vsize == size);
	}
	vit = vertices.begin();
}
#endif
