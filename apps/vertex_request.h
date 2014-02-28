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
