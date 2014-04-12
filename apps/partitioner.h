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

#include <utility>

#include "vertex.h"

/**
 * This data structure represents the local Id of a vertex used
 * in its own partition.
 */
struct local_vid_t
{
	vertex_id_t id;

	local_vid_t() {
		id = INVALID_VERTEX_ID;
	}

	explicit local_vid_t(vertex_id_t id) {
		this->id = id;
	}
};

/**
 * The vertex location in a partition.
 * first: the partition ID.
 * second: the location in the partition.
 */
typedef std::pair<int, struct local_vid_t> vertex_loc_t;

class graph_partitioner
{
public:
	virtual int map(vertex_id_t id) const = 0;
	virtual void map2loc(vertex_id_t id, int &part_id, off_t &off) const = 0;
	virtual void map2loc(vertex_id_t ids[], int num,
			vertex_loc_t locs[]) const = 0;
	virtual void map2loc(edge_seq_iterator &, vertex_loc_t locs[]) const = 0;
	virtual void loc2map(int part_id, off_t off, vertex_id_t &id) const = 0;
	virtual size_t get_all_vertices_in_part(int part_id,
			size_t tot_num_vertices, std::vector<vertex_id_t> &ids) const = 0;
	virtual size_t get_part_size(int part_id, size_t num_vertices) const = 0;
};

class modulo_graph_partitioner: public graph_partitioner
{
	int num_parts_log;
	vertex_id_t mask;
public:
	modulo_graph_partitioner(int num_parts) {
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
	void map2loc(vertex_id_t ids[], int num, vertex_loc_t locs[]) const;
	void map2loc(edge_seq_iterator &, vertex_loc_t locs[]) const;

	void loc2map(int part_id, off_t off, vertex_id_t &id) const {
		id = (off << num_parts_log) + part_id;
	}

	size_t get_all_vertices_in_part(int part_id, size_t tot_num_vertices,
			std::vector<vertex_id_t> &ids) const;

	virtual size_t get_part_size(int part_id, size_t num_vertices) const {
		int num_parts = 1 << num_parts_log;
		return (num_vertices - part_id + num_parts - 1) / num_parts;
	}
};

class range_graph_partitioner: public graph_partitioner
{
	static const int RANGE_SIZE_LOG = 10;
	static const int RANGE_SIZE = 1 << RANGE_SIZE_LOG;
	static const unsigned long RANGE_MASK = (1UL << RANGE_SIZE_LOG) - 1;

	const int num_parts_log;
	const vertex_id_t mask;

	vertex_id_t get_part_end(int part_id, size_t num_vertices) const;
public:
	range_graph_partitioner(int num_parts): num_parts_log(
			log2(num_parts)), mask((1 << num_parts_log) - 1) {
		assert((1 << num_parts_log) == num_parts);
	}

	virtual int map(vertex_id_t id) const {
		return (id >> RANGE_SIZE_LOG) & mask;
	}

	virtual void map2loc(vertex_id_t id, int &part_id, off_t &off) const {
		vertex_id_t shifted_id = id >> RANGE_SIZE_LOG;
		part_id = shifted_id & mask;
		off = ((shifted_id >> num_parts_log) << RANGE_SIZE_LOG) + (id & RANGE_MASK);
	}

	virtual void map2loc(vertex_id_t ids[], int num,
			vertex_loc_t locs[]) const;
	virtual void map2loc(edge_seq_iterator &, vertex_loc_t locs[]) const;

	virtual void loc2map(int part_id, off_t off, vertex_id_t &id) const {
		off_t local_range_id = off >> RANGE_SIZE_LOG;
		off_t in_range_off = off & RANGE_MASK;
		id = (local_range_id << (num_parts_log + RANGE_SIZE_LOG))
			+ (part_id << RANGE_SIZE_LOG) + in_range_off;
	}

	virtual size_t get_all_vertices_in_part(int part_id,
			size_t tot_num_vertices, std::vector<vertex_id_t> &ids) const;

	virtual size_t get_part_size(int part_id, size_t num_vertices) const {
		vertex_id_t vid_end = get_part_end(part_id, num_vertices);
		size_t num_local_ranges = vid_end >> (RANGE_SIZE_LOG + num_parts_log);
		size_t num_vertices_last_range = vid_end
			- (num_local_ranges << (RANGE_SIZE_LOG + num_parts_log))
			- (part_id << RANGE_SIZE_LOG);
		num_vertices_last_range = max(0, num_vertices_last_range);
		return (num_local_ranges << RANGE_SIZE_LOG) + num_vertices_last_range;
	}
};

#endif
