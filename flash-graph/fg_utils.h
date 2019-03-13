#ifndef __FG_UTILS_H__
#define __FG_UTILS_H__

/*
 * Copyright 2014 Open Connectome Project (http://openconnecto.me)
 * Written by Da Zheng (zhengda1936@gmail.com)
 *
 * This file is part of FlashMatrix.
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
#include "vertex_index.h"
#include "FGlib.h"

namespace fg
{

/*
 * Get the total size of the in-edge lists.
 */
static inline size_t get_in_size(fg::vertex_index::ptr vindex)
{
	assert(vindex->get_graph_header().is_directed_graph());
	return vindex->get_out_part_loc() - vindex->get_header_size();
}

/*
 * Get the offset of the first out-edge list in the graph image.
 */
static inline size_t get_out_off(fg::vertex_index::ptr vindex)
{
	if (vindex->get_graph_header().is_directed_graph())
		return vindex->get_out_part_loc();
	else
		return vindex->get_header_size();
}

/*
 * Get the offset of the first in-edge list in the graph image.
 */
static inline size_t get_in_off(fg::vertex_index::ptr vindex)
{
	assert(vindex->get_graph_header().is_directed_graph());
	return vindex->get_header_size();
}

/*
 * Get the total size of the out-edge lists.
 */
size_t get_out_size(fg::vertex_index::ptr vindex);
/*
 * Get the offsets of all out-edge lists in the graph image.
 */
void init_out_offs(fg::vertex_index::ptr vindex, std::vector<off_t> &out_offs);
/*
 * Get the offsets of all in-edge lists in the graph image.
 */
void init_in_offs(fg::vertex_index::ptr vindex, std::vector<off_t> &in_offs);

/*
 * Fetch the subgraph that contains the specified vertices from the graph.
 * If `compact' is false, the constructed subgraph contains the same number
 * of vertices as the original graph but the unspecified vertices contains
 * no neighbors.
 * If `compact' is true, the constructed subgraph contains only the specified
 * vertices.
 */
fg::FG_graph::ptr fetch_subgraph(fg::FG_graph::ptr graph,
		const std::vector<fg::vertex_id_t> &vertices,
		const std::string &graph_name, bool compact);

/*
 * A set of functions that change the parameters in the matrix generator.
 */

void set_deduplicate(bool v);
void set_remove_self_edge(bool v);
}

#endif
