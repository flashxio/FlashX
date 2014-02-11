#ifndef __GRAPH_FILE_HEADER_H__
#define __GRAPH_FILE_HEADER_H__

/**
 * Copyright 2013 Da Zheng
 *
 * This file is part of SA-GraphLib.
 *
 * SA-GraphLib is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * SA-GraphLib is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with SA-GraphLib.  If not, see <http://www.gnu.org/licenses/>.
 */

#include <stdlib.h>

const int64_t MAGIC_NUMBER = 0x123456789ABCDEFL;
const int CURR_VERSION = 3;

enum graph_type {
	DIRECTED,
	UNDIRECTED,
	TS_DIRECTED,
	TS_UNDIRECTED,
};

/**
 * This is the graph header. It should be stored on the files
 * of both vertex index and adjacency list.
 */
class graph_header
{
	union {
		struct {
			int64_t magic_number;
			int version_number;
			graph_type type;
			size_t num_vertices;
			size_t num_edges;
			bool has_data;
			// This is only used for time-series graphs.
			int max_num_timestamps;
		} data;
		char page[PAGE_SIZE];
	} h;
public:
	graph_header() {
		assert(sizeof(*this) == PAGE_SIZE);
		memset(this, 0, sizeof(*this));
		h.data.magic_number = MAGIC_NUMBER;
		h.data.version_number = CURR_VERSION;
		h.data.type = DIRECTED;
		h.data.num_vertices = 0;
		h.data.num_edges = 0;
		h.data.has_data = false;
		h.data.max_num_timestamps = 0;
	}

	graph_header(graph_type type, size_t num_vertices, size_t num_edges,
			bool has_edge_data, int max_num_timestamps = 0) {
		assert(sizeof(*this) == PAGE_SIZE);
		memset(this, 0, sizeof(*this));
		h.data.magic_number = MAGIC_NUMBER;
		h.data.version_number = CURR_VERSION;
		h.data.type = type;
		h.data.num_vertices = num_vertices;
		h.data.num_edges = num_edges;
		h.data.has_data = has_edge_data;
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
		return h.data.has_data;
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

#endif
