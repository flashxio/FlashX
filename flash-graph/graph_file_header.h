#ifndef __GRAPH_FILE_HEADER_H__
#define __GRAPH_FILE_HEADER_H__

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

#include <stdlib.h>
#include <string.h>

#include "parameters.h"

namespace fg
{

const int64_t MAGIC_NUMBER = 0x123456789ABCDEFL;
const int CURR_VERSION = 4;

enum graph_type {
	DIRECTED,
	UNDIRECTED,
	TS_DIRECTED,
	TS_UNDIRECTED,
};

struct graph_header_struct
{
	int64_t magic_number;
	int version_number;
	graph_type type;
	size_t num_vertices;
	size_t num_edges;
	int edge_data_size;
	// This is only used for time-series graphs.
	int max_num_timestamps;
};

/**
 * This is the graph header. It should be stored on the files
 * of both vertex index and adjacency list.
 */
class graph_header
{
	union {
		struct graph_header_struct data;
		char page[PAGE_SIZE];
	} h;
public:
	static int get_header_size() {
		return sizeof(graph_header);
	}

	static void init(struct graph_header_struct &data) {
		data.magic_number = MAGIC_NUMBER;
		data.version_number = CURR_VERSION;
		data.type = DIRECTED;
		data.num_vertices = 0;
		data.num_edges = 0;
		data.edge_data_size = 0;
		data.max_num_timestamps = 0;
	}

	graph_header() {
		assert(sizeof(*this) == PAGE_SIZE);
		memset(this, 0, sizeof(*this));
		init(h.data);
	}

	graph_header(graph_type type, size_t num_vertices, size_t num_edges,
			int edge_data_size, int max_num_timestamps = 0) {
		assert(sizeof(*this) == PAGE_SIZE);
		memset(this, 0, sizeof(*this));
		h.data.magic_number = MAGIC_NUMBER;
		h.data.version_number = CURR_VERSION;
		h.data.type = type;
		h.data.num_vertices = num_vertices;
		h.data.num_edges = num_edges;
		h.data.edge_data_size = edge_data_size;
		h.data.max_num_timestamps = max_num_timestamps;
	}

	bool is_graph_file() const {
		return h.data.magic_number == MAGIC_NUMBER;
	}

	bool is_right_version() const {
		return h.data.version_number == CURR_VERSION;
	}

	bool is_directed_graph() const {
		return h.data.type == graph_type::DIRECTED
			|| h.data.type == graph_type::TS_DIRECTED;
	}

	graph_type get_graph_type() const {
		return h.data.type;
	}

	size_t get_num_vertices() const {
		return h.data.num_vertices;
	}

	size_t get_num_edges() const {
		return h.data.num_edges;
	}

	bool has_edge_data() const {
		return h.data.edge_data_size > 0;
	}

	int get_edge_data_size() const {
		return h.data.edge_data_size;
	}

	int get_max_num_timestamps() const {
		return h.data.max_num_timestamps;
	}

	void verify() const {
		if (!is_graph_file()) {
			fprintf(stderr, "wrong magic number: %ld\n", h.data.magic_number);
		}
		if (!is_right_version()) {
			fprintf(stderr, "wrong version number: %d\n", h.data.version_number);
		}
		assert(is_graph_file());
		assert(is_right_version());
	}
};

}

#endif
