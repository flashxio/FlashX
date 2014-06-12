#ifndef __VERTEX_REQUEST_H__
#define __VERTEX_REQUEST_H__

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

#include "vertex.h"

class graph_engine;

class vertex_request
{
	vertex_id_t id;
public:
	vertex_request() {
		id = -1;
	}

	vertex_request(vertex_id_t id) {
		this->id = id;
	}

	vertex_id_t get_id() const {
		return id;
	}
};

class directed_vertex_request: public vertex_request
{
	edge_type type;
public:
	directed_vertex_request() {
		type = edge_type::BOTH_EDGES;
	}

	directed_vertex_request(vertex_id_t id, edge_type type): vertex_request(id) {
		this->type = type;
	}

	edge_type get_type() const {
		return type;
	}
};

/**
 * This class contains the request of a time-series vertex
 * from the user application.
 */
class ts_vertex_request: public vertex_request
{
	timestamp_pair range;
public:
	ts_vertex_request() {
		range = timestamp_pair(INT_MAX, INT_MIN);
	}

	ts_vertex_request(vertex_id_t id): vertex_request(id) {
		range = timestamp_pair(INT_MAX, INT_MIN);
	}

	ts_vertex_request(vertex_id_t id, timestamp_pair range): vertex_request(
			id) {
		this->range = range;
	}

	void add_timestamp(int timestamp) {
		if (range.second < timestamp)
			range.second = timestamp + 1;
		if (range.first > timestamp)
			range.first = timestamp;
	}

	const timestamp_pair &get_range() const {
		return range;
	}

	bool is_require_all() const {
		return range.first == INT_MAX && range.second == INT_MIN;
	}
};

#endif
