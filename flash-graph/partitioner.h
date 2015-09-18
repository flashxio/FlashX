#ifndef __GRAPH_PARTITIONER_H__
#define __GRAPH_PARTITIONER_H__

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
	virtual int get_num_partitions() const = 0;
	virtual int map(vertex_id_t id) const = 0;
	virtual void map2loc(vertex_id_t id, int &part_id, off_t &off) const = 0;
	virtual void map2loc(vertex_id_t ids[], int num,
			std::vector<local_vid_t> locs[], int num_parts) const = 0;
	virtual size_t map2loc(edge_seq_iterator &,
			vertex_loc_t locs[], size_t size) const = 0;
	virtual void map2loc(edge_seq_iterator &, std::vector<local_vid_t> locs[],
			int num_parts) const = 0;
	virtual void loc2map(int part_id, off_t off, vertex_id_t &id) const = 0;
	virtual size_t get_all_vertices_in_part(int part_id,
			size_t tot_num_vertices, std::vector<vertex_id_t> &ids) const = 0;
	virtual size_t get_part_size(int part_id, size_t num_vertices) const = 0;
};

class modulo_graph_partitioner: public graph_partitioner
{
	int num_parts_log;
	vertex_id_t mask;
	int get_num_parts() const {
		return 1 << num_parts_log;
	}
public:
	modulo_graph_partitioner(int num_parts) {
		this->num_parts_log = log2(num_parts);
		assert((1 << num_parts_log) == num_parts);
		mask = (1 << num_parts_log) - 1;
	}

	int get_num_partitions() const {
		return 1 << num_parts_log;
	}

	int map(vertex_id_t id) const {
		return id & mask;
	}

	void map2loc(vertex_id_t id, int &part_id, off_t &off) const {
		part_id = id & mask;
		off = id >> num_parts_log;
	}
	void map2loc(vertex_id_t ids[], int num,
			std::vector<local_vid_t> locs[], int num_parts) const;
	void map2loc(edge_seq_iterator &, std::vector<local_vid_t> locs[],
			int num_parts) const;
	size_t map2loc(edge_seq_iterator &,
			vertex_loc_t locs[], size_t size) const;

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
	const int RANGE_SIZE_LOG;
	const int RANGE_SIZE;
	const unsigned long RANGE_MASK;

	const int num_parts_log;
	const vertex_id_t mask;

	int get_num_parts() const {
		return 1 << num_parts_log;
	}

	vertex_id_t get_part_end(int part_id, size_t num_vertices) const;
public:
	range_graph_partitioner(int num_parts);

	int get_num_partitions() const {
		return 1 << num_parts_log;
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
			std::vector<local_vid_t> locs[], int num_parts) const;
	virtual void map2loc(edge_seq_iterator &, std::vector<local_vid_t> locs[],
			int num_parts) const;
	virtual size_t map2loc(edge_seq_iterator &,
			vertex_loc_t locs[], size_t num) const;

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
