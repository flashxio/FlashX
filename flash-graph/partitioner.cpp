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

#include "partitioner.h"
#include "graph_config.h"

namespace fg
{

range_graph_partitioner::range_graph_partitioner(int num_parts): RANGE_SIZE_LOG(
		graph_conf.get_part_range_size_log()), RANGE_SIZE(1 << RANGE_SIZE_LOG),
	RANGE_MASK((1UL << RANGE_SIZE_LOG) - 1), num_parts_log(log2(num_parts)),
	mask((1 << num_parts_log) - 1) {
		assert((1 << num_parts_log) == num_parts);
}

size_t modulo_graph_partitioner::get_all_vertices_in_part(int part_id,
		size_t tot_num_vertices, std::vector<vertex_id_t> &ids) const
{
	int num_parts = 1 << num_parts_log;
	for (size_t i = part_id; i < tot_num_vertices; i += num_parts)
		ids.push_back(i);
	return ids.size();
}

size_t range_graph_partitioner::get_all_vertices_in_part(int part_id,
			size_t tot_num_vertices, std::vector<vertex_id_t> &ids) const
{
	vertex_id_t vid_end = get_part_end(part_id, tot_num_vertices);
	for (vertex_id_t vid_start = part_id * RANGE_SIZE; vid_start < vid_end;
			vid_start += 1 << (num_parts_log + RANGE_SIZE_LOG)) {
		int size = min(RANGE_SIZE, vid_end - vid_start);
		for (int i = 0; i < size; i++) {
			ids.push_back(vid_start + i);
		}
	}
	return ids.size();
}

vertex_id_t range_graph_partitioner::get_part_end(int part_id,
		size_t tot_num_vertices) const
{
	// The number of vertex ID ranges in the partition, excluding the range
	// in the last stripe if there is one.
	size_t num_ranges = tot_num_vertices >> (RANGE_SIZE_LOG + num_parts_log);
	assert((num_ranges << (RANGE_SIZE_LOG + num_parts_log))
			<= tot_num_vertices);
	assert(tot_num_vertices
			< ((num_ranges + 1) << (RANGE_SIZE_LOG + num_parts_log)));

	// The number of vertices in the last stripe
	size_t num = tot_num_vertices
		- (num_ranges << (RANGE_SIZE_LOG + num_parts_log));
	// The number of vertices in the last stripe in the partition
	for (int i = 0; i < part_id; i++)
		num -= RANGE_SIZE;

	vertex_id_t vid_end = part_id * RANGE_SIZE
		+ (num_ranges << (RANGE_SIZE_LOG + num_parts_log));
	if (num > 0)
		vid_end += min(num, RANGE_SIZE);
	return vid_end;
}

void modulo_graph_partitioner::map2loc(vertex_id_t ids[], int num,
		std::vector<local_vid_t> locs[], int num_parts) const
{
	assert(num_parts <= get_num_parts());
	for (int i = 0; i < num; i++) {
		int part_id = ids[i] & mask;
		locs[part_id].push_back(local_vid_t(ids[i] >> num_parts_log));
	}
}

void range_graph_partitioner::map2loc(vertex_id_t ids[], int num,
		std::vector<local_vid_t> locs[], int num_parts) const
{
	assert(num_parts <= get_num_parts());
	for (int i = 0; i < num; i++) {
		vertex_id_t id = ids[i];
		vertex_id_t shifted_id = id >> RANGE_SIZE_LOG;
		int part_id = shifted_id & mask;
		locs[part_id].push_back(local_vid_t(((shifted_id >> num_parts_log) << RANGE_SIZE_LOG)
					+ (id & RANGE_MASK)));
	}
}

void modulo_graph_partitioner::map2loc(edge_seq_iterator &it,
		std::vector<local_vid_t> locs[], int num_parts) const
{
	assert(num_parts <= get_num_parts());
	PAGE_FOREACH(vertex_id_t, id, it) {
		int part_id = id & mask;
		locs[part_id].push_back(local_vid_t(id >> num_parts_log));
	} PAGE_FOREACH_END
}

void range_graph_partitioner::map2loc(edge_seq_iterator &it,
		std::vector<local_vid_t> locs[], int num_parts) const
{
	assert(num_parts <= get_num_parts());
	PAGE_FOREACH(vertex_id_t, id, it) {
		vertex_id_t shifted_id = id >> RANGE_SIZE_LOG;
		int part_id = shifted_id & mask;
		locs[part_id].push_back(local_vid_t(((shifted_id >> num_parts_log) << RANGE_SIZE_LOG)
					+ (id & RANGE_MASK)));
	} PAGE_FOREACH_END
}

size_t modulo_graph_partitioner::map2loc(edge_seq_iterator &it,
		vertex_loc_t locs[], size_t num) const
{
	size_t ret = 0;
	PAGE_FOREACH(vertex_id_t, id, it) {
		if ((size_t) page_foreach_idx == num)
			break;
		vertex_loc_t loc(id & mask, local_vid_t(id >> num_parts_log));
		locs[page_foreach_idx] = loc;
		ret++;
	} PAGE_FOREACH_END
	return ret;
}

size_t range_graph_partitioner::map2loc(edge_seq_iterator &it,
		vertex_loc_t locs[], size_t num) const
{
	size_t ret = 0;
	PAGE_FOREACH(vertex_id_t, id, it) {
		if ((size_t) page_foreach_idx == num)
			break;
		vertex_id_t shifted_id = id >> RANGE_SIZE_LOG;
		vertex_loc_t loc(shifted_id & mask,
				local_vid_t(((shifted_id >> num_parts_log) << RANGE_SIZE_LOG)
					+ (id & RANGE_MASK)));
		locs[page_foreach_idx] = loc;
		ret++;
	} PAGE_FOREACH_END
	return ret;
}

}
