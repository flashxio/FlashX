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

#include "partitioner.h"

size_t modulo_vertex_partitioner::get_all_vertices_in_part(int part_id,
		size_t tot_num_vertices, std::vector<vertex_id_t> &ids) const
{
	int num_parts = 1 << num_parts_log;
	for (size_t i = part_id; i < tot_num_vertices; i += num_parts)
		ids.push_back(i);
	return ids.size();
}

size_t range_vertex_partitioner::get_all_vertices_in_part(int part_id,
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

vertex_id_t range_vertex_partitioner::get_part_end(int part_id,
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
