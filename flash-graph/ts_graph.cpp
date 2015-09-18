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

#include "ts_graph.h"

page_byte_array::seq_const_iterator<vertex_id_t> get_ts_iterator(
		const page_directed_vertex &v, edge_type type, time_t time_start,
		time_t time_interval)
{
	page_byte_array::const_iterator<ts_edge_data> begin_it
		= v.get_data_begin<ts_edge_data>(type);
	page_byte_array::const_iterator<ts_edge_data> end_it
		= v.get_data_end<ts_edge_data>(type);
	page_byte_array::const_iterator<ts_edge_data> ts_it = std::lower_bound(
			begin_it, end_it, ts_edge_data(time_start));
	// All timestamps are smaller than time_start
	if (ts_it == end_it)
		return v.get_neigh_seq_it(type, 0, 0);
	size_t start = ts_it - begin_it;

	page_byte_array::const_iterator<ts_edge_data> ts_end_it = std::lower_bound(
			begin_it, end_it, time_start + time_interval);
	size_t end = ts_end_it - begin_it;

	return v.get_neigh_seq_it(type, start, end);
}
