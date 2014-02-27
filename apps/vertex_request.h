#ifndef __VERTEX_REQUEST_H__
#define __VERTEX_REQUEST_H__

/**
 * Copyright 2014 Da Zheng
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

#include "vertex.h"

class graph_engine;

class vertex_request
{
};

/**
 * This class contains the request of a time-series vertex
 * from the user application.
 */
class ts_vertex_request: public vertex_request
{
	vertex_id_t id;
	timestamp_pair range;
	edge_type type;
	bool require_all;
	graph_engine *graph;
public:
	ts_vertex_request(graph_engine *graph) {
		id = 0;
		range = timestamp_pair(INT_MAX, INT_MIN);
		type = edge_type::BOTH_EDGES;
		require_all = false;
		this->graph = graph;
	}

	void set_require_all(bool require_all) {
		this->require_all = require_all;
	}

	void set_vertex(vertex_id_t id);

	void add_timestamp(int timestamp) {
		if (!require_all) {
			if (range.second < timestamp)
				range.second = timestamp + 1;
			if (range.first > timestamp)
				range.first = timestamp;
		}
	}

	void set_edge_type(edge_type type) {
		this->type = type;
	}

	void clear() {
		id = 0;
		range = timestamp_pair(INT_MAX, INT_MIN);
		type = edge_type::BOTH_EDGES;
		require_all = false;
	}

	vertex_id_t get_id() const {
		return id;
	}

	const timestamp_pair &get_range() const {
		return range;
	}

	edge_type get_edge_type() const {
		return type;
	}

	bool is_require_all() const {
		return require_all;
	}
};

#endif
