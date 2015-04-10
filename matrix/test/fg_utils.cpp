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

#include "fg_utils.h"

size_t get_out_size(fg::vertex_index::ptr vindex)
{
	if (vindex->is_compressed()) {
		fg::vsize_t num_vertices = vindex->get_num_vertices();
		if (vindex->get_graph_header().is_directed_graph()) {
			fg::in_mem_cdirected_vertex_index::ptr dindex
				= fg::in_mem_cdirected_vertex_index::create(*vindex);
			fg::directed_vertex_entry dentry = dindex->get_vertex(num_vertices - 1);
			return dentry.get_out_off() + dindex->get_out_size(
					num_vertices - 1) - vindex->get_out_part_loc();
		}
		else {
			fg::in_mem_cundirected_vertex_index::ptr uindex
				= fg::in_mem_cundirected_vertex_index::create(*vindex);
			fg::vertex_offset off = uindex->get_vertex(num_vertices - 1);
			return off.get_off() + uindex->get_size(
					num_vertices - 1) - vindex->get_header_size();
		}
	}
	else {
		if (vindex->get_graph_header().is_directed_graph()) {
			fg::directed_vertex_index::ptr dindex
				= fg::directed_vertex_index::cast(vindex);
			return dindex->get_graph_size() - vindex->get_out_part_loc();
		}
		else {
			fg::undirected_vertex_index::ptr uindex
				= fg::undirected_vertex_index::cast(vindex);
			return uindex->get_graph_size() - vindex->get_header_size();
		}
	}
}

void init_out_offs(fg::vertex_index::ptr vindex, std::vector<off_t> &out_offs)
{
	size_t num_vertices = vindex->get_num_vertices();
	assert(num_vertices + 1 == out_offs.size());
	if (vindex->is_compressed()) {
		if (vindex->get_graph_header().is_directed_graph()) {
			fg::in_mem_cdirected_vertex_index::ptr dindex
				= fg::in_mem_cdirected_vertex_index::create(*vindex);
			out_offs[0] = 0;
			for (size_t i = 1; i <= num_vertices; i++)
				out_offs[i] = out_offs[i - 1] + dindex->get_out_size(i - 1);
		}
		else {
			fg::in_mem_cundirected_vertex_index::ptr uindex
				= fg::in_mem_cundirected_vertex_index::create(*vindex);
			out_offs[0] = 0;
			for (size_t i = 1; i <= num_vertices; i++)
				out_offs[i] = out_offs[i - 1] + uindex->get_size(i - 1);
		}
		assert((size_t) out_offs[num_vertices] == get_out_size(vindex));
	}
	else {
		if (vindex->get_graph_header().is_directed_graph()) {
			off_t out_part_loc = vindex->get_out_part_loc();
			fg::directed_vertex_index::ptr dindex
				= fg::directed_vertex_index::cast(vindex);
			for (size_t i = 0; i < num_vertices; i++)
				out_offs[i] = dindex->get_vertex(i).get_out_off()
					- out_part_loc;
			out_offs[num_vertices] = get_out_size(vindex);
		}
		else {
			fg::undirected_vertex_index::ptr uindex
				= fg::undirected_vertex_index::cast(vindex);
			for (size_t i = 0; i < num_vertices; i++)
				out_offs[i] = uindex->get_vertex(i).get_off()
					- vindex->get_header_size();
			out_offs[num_vertices] = get_out_size(vindex);
		}
	}
	for (size_t i = 1; i <= num_vertices; i++)
		assert(out_offs[i] > out_offs[i - 1]);
}

void init_in_offs(fg::vertex_index::ptr vindex, std::vector<off_t> &in_offs)
{
	size_t num_vertices = vindex->get_num_vertices();
	assert(num_vertices + 1 == in_offs.size());
	assert(vindex->get_graph_header().is_directed_graph());
	if (vindex->is_compressed()) {
		fg::in_mem_cdirected_vertex_index::ptr dindex
			= fg::in_mem_cdirected_vertex_index::create(*vindex);
		in_offs[0] = 0;
		for (size_t i = 1; i <= num_vertices; i++)
			in_offs[i] = in_offs[i - 1] + dindex->get_in_size(i - 1);
		assert((size_t) in_offs[num_vertices] == get_in_size(vindex));
	}
	else {
		fg::directed_vertex_index::ptr dindex
			= fg::directed_vertex_index::cast(vindex);
		for (size_t i = 0; i < num_vertices; i++)
			in_offs[i] = dindex->get_vertex(i).get_in_off()
				- vindex->get_header_size();
		in_offs[num_vertices] = get_in_size(vindex);
	}
	for (size_t i = 1; i <= num_vertices; i++)
		assert(in_offs[i] > in_offs[i - 1]);
}
