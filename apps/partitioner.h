#ifndef __GRAPH_PARTITIONER_H__
#define __GRAPH_PARTITIONER_H__

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

#include <math.h>

#include "vertex.h"

class vertex_partitioner
{
	int num_parts_log;
	vertex_id_t mask;
public:
	vertex_partitioner(int num_parts) {
		this->num_parts_log = log2(num_parts);
		assert((1 << num_parts_log) == num_parts);
		mask = (1 << num_parts_log) - 1;
	}

	int map(vertex_id_t id) const {
		return id & mask;
	}

	void map2loc(vertex_id_t id, int &part_id, off_t &off) const {
		part_id = id & mask;
		off = id >> num_parts_log;
	}

	void loc2map(int part_id, off_t off, vertex_id_t &id) const {
		id = (off << num_parts_log) + part_id;
	}

	size_t get_all_vertices_in_part(int part_id, size_t tot_num_vertices,
			std::vector<vertex_id_t> &ids) const {
		int num_parts = 1 << num_parts_log;
		for (size_t i = part_id; i < tot_num_vertices; i += num_parts)
			ids.push_back(i);
		return ids.size();
	}
};


#endif
