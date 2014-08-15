#ifndef __EXT_MEM_VERTEX_ITERATOR_H__
#define __EXT_MEM_VERTEX_ITERATOR_H__

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

#include "vertex_index.h"
#include "vertex.h"

#if 0

class ext_mem_vertex_iterator
{
	static const size_t MEM_SIZE = 1024 * 1024 * 1024;
	graph_header header;
	vertex_index_iterator *index_it;
	FILE *adj_f;
	char *adj_buf;
	std::vector<in_mem_vertex_info> overflow_vinfos;
	std::vector<char *> vertices;
	std::vector<char *>::const_iterator vit;

	void read_vertices();
	ext_mem_vertex_iterator(const ext_mem_vertex_iterator &);
	ext_mem_vertex_iterator &operator=(const ext_mem_vertex_iterator &);
public:
	ext_mem_vertex_iterator() {
		index_it = NULL;
		adj_f = NULL;
		adj_buf = NULL;
	}

	ext_mem_vertex_iterator(const std::string &index_file,
			const std::string &adj_file) {
		init(index_file, adj_file);
	}

	~ext_mem_vertex_iterator() {
		if (is_valid()) {
			fclose(adj_f);
			graph_type type = header.get_graph_type();
			if (type == graph_type::DIRECTED)
				directed_vertex_index_iterator::destroy(
						(directed_vertex_index_iterator *) index_it);
			else
				default_vertex_index_iterator::destroy(
						(default_vertex_index_iterator *) index_it);
			delete [] adj_buf;
		}
	}

	bool is_valid() const {
		return index_it != NULL;
	}

	void init(const std::string &index_file, const std::string &adj_file) {
		adj_f = fopen(adj_file.c_str(), "r");
		assert(adj_f);
		size_t ret = fread(&header, sizeof(header), 1, adj_f);
		assert(ret == 1);

		graph_type type = header.get_graph_type();
		if (type == graph_type::DIRECTED)
			index_it = directed_vertex_index_iterator::create(index_file);
		else
			index_it = default_vertex_index_iterator::create(index_file);

		adj_buf = new char[MEM_SIZE];
		// It's possible the graph has no vertices at all.
		if (index_it->has_next())
			read_vertices();
	}

	const graph_header &get_graph_header() const {
		return header;
	}

	bool has_next() const {
		if (vit != vertices.end())
			return true;
		else if (!overflow_vinfos.empty())
			return true;
		else
			return index_it->has_next();
	}

	/*
	 * The iterator uses an internal memory buffer to keep the vertices
	 * read from the disks. Once next() or next_vertices() is invoked,
	 * all vertices returned by the previous next() and next_vertices()
	 * will become invalid.
	 */

	template<class vertex_type>
	vertex_type *next() {
		if (vit == vertices.end())
			read_vertices();

		vertex_type *v = (vertex_type *) *vit;
		vit++;
		return v;
	}

	template<class vertex_type>
	size_t next_vertices(std::vector<vertex_type *> &ret) {
		if (vit == vertices.end())
			read_vertices();

		while (vit != vertices.end()) {
			ret.push_back((vertex_type *) *vit);
			vit++;
		}
		return ret.size();
	}
};
#endif

#endif
