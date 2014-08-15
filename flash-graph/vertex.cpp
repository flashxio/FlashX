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
#include "vertex_index.h"

empty_data edge<empty_data>::data;

size_t ext_mem_undirected_vertex::serialize(const in_mem_vertex &v, char *buf,
		size_t size, edge_type type)
{
	assert(ext_mem_undirected_vertex::get_header_size() <= size);
	ext_mem_undirected_vertex *ext_v = (ext_mem_undirected_vertex *) buf;
	ext_v->set_id(v.get_id());
	ext_v->num_edges = v.get_num_edges(type);
	if (v.has_edge_data())
		ext_v->edge_data_size = v.get_edge_data_size();
	else
		ext_v->edge_data_size = 0;
	size_t mem_size = ext_v->get_size();
	assert(mem_size <= MAX_VERTEX_SIZE);
	assert(size >= mem_size);
	v.serialize_edges(ext_v->neighbors, type);
	// serialize edge data
	if (v.has_edge_data()) {
		v.serialize_edge_data(ext_v->get_edge_data_addr(), type);
	}

	return mem_size;
}

#if 0
ts_ext_mem_directed_vertex::unique_ptr ts_ext_mem_directed_vertex::merge(
		const std::vector<const ext_mem_directed_vertex *> &vertices)
{
	// Collect information from the array of vertices.
	assert(vertices.size() > 0);
	vertex_id_t id = INVALID_VERTEX_ID;
	size_t edge_data_size = 0;
	vsize_t num_edges = 0;
	int num_timestamps = 0;
	for (size_t i = 0; i < vertices.size(); i++) {
		const ext_mem_directed_vertex *v = vertices[i];
		if (v == NULL)
			continue;

		if (edge_data_size == 0)
			edge_data_size = v->get_edge_data_size();
		if (id == INVALID_VERTEX_ID)
			id = v->get_id();
		assert(edge_data_size == v->get_edge_data_size());
		assert(id == v->get_id());
		num_edges += v->get_num_in_edges() + v->get_num_out_edges();
		if (v->get_num_in_edges() + v->get_num_out_edges() > 0)
			num_timestamps++;
	}
	assert(id != INVALID_VERTEX_ID);

	// Create the time-series vertex.
	size_t size = ts_ext_mem_directed_vertex::get_vertex_size(num_timestamps,
			num_edges, edge_data_size);
	char *buf = new char[size];
	ts_ext_mem_directed_vertex *ts_v = new (buf) ts_ext_mem_directed_vertex(
			id, num_edges, num_timestamps, edge_data_size);
	assert(ts_v->get_size() <= MAX_VERTEX_SIZE);

	// Create the timestamp table.
	int ts_idx = 0;
	size_t off = 0;
	for (size_t i = 0; i < vertices.size(); i++) {
		const ext_mem_directed_vertex *v = vertices[i];
		if (v == NULL)
			continue;
		if (v->get_num_in_edges() + v->get_num_out_edges() == 0)
			continue;

		ts_v->get_timestamps_begin()[ts_idx] = i;
		ts_v->get_edge_off_begin()[ts_idx].in_off = off;
		off += v->get_num_in_edges();
		ts_v->get_edge_off_begin()[ts_idx].out_off = off;
		off += v->get_num_out_edges();
		ts_idx++;
	}

	// Create the edge list.
	char *edge_buf = (char *) ts_v->get_edge_list_begin();
	for (size_t i = 0; i < vertices.size(); i++) {
		const ext_mem_directed_vertex *v = vertices[i];
		if (v == NULL)
			continue;
		if (v->get_num_in_edges() + v->get_num_out_edges() == 0)
			continue;

		size_t v_edge_size = v->get_size()
			- ext_mem_directed_vertex::get_header_size();
		memcpy(edge_buf, v->neighbors, v_edge_size);
		edge_buf += v_edge_size;
	}

	size_t v_size = edge_buf - buf;
	assert(v_size == ts_v->get_size());
	return unique_ptr(ts_v);
}

void ts_ext_mem_directed_vertex::construct_header(
		const ts_ext_mem_directed_vertex &header,
		off_t edge_list_off, size_t edge_list_size)
{
	*this = header;
	int timestamp_idx_start = -1;
	int timestamp_idx_end = -1;
	off_t edge_list_end = edge_list_off + edge_list_size;
	for (int i = 0; i < num_timestamps; i++) {
		const_edge_list list = header.get_const_edge_list_idx(i);
		off_t off = ((char *) list.get_in_edge_begin()) - ((char *) &header);
		if (edge_list_off == off)
			timestamp_idx_start = i;
		if (edge_list_end == off)
			timestamp_idx_end = i;
	}
	if ((size_t) edge_list_end == header.get_size())
		timestamp_idx_end = num_timestamps;

	assert(timestamp_idx_start >= 0);
	assert(timestamp_idx_end >= 0);
	this->num_timestamps = timestamp_idx_end - timestamp_idx_start;
	memcpy(this->get_timestamps_begin(), header.get_timestamps_begin()
			+ timestamp_idx_start, sizeof(short) * this->num_timestamps);
	size_t num_prev_edges
		= header.get_edge_off_begin()[timestamp_idx_start].in_off;
	for (int i = 0; i < this->num_timestamps; i++) {
		this->get_edge_off_begin()[i].in_off
			= header.get_edge_off_begin()[timestamp_idx_start + i].in_off
			- num_prev_edges;
		this->get_edge_off_begin()[i].out_off
			= header.get_edge_off_begin()[timestamp_idx_start + i].out_off
			- num_prev_edges;
	}
	if (timestamp_idx_end == header.get_num_timestamps())
		this->num_edges = header.get_num_edges() - num_prev_edges;
	else
		this->num_edges = header.get_edge_off_begin()[timestamp_idx_end].in_off
			- num_prev_edges;
	assert(this->num_edges * (sizeof(vertex_id_t) + this->edge_data_size)
			== edge_list_size);
	assert(this->get_edge_off_begin()[0].in_off == 0);
}
#endif

offset_pair TS_page_directed_vertex::get_edge_list_offset(
		const timestamp_pair &range) const
{
	assert(range.first <= range.second);
	int i;
	// Find the timestamp range in the vertex that is inside the specified
	// timestamp range.
	timestamp_pair new_range(ext_v.get_last_timestamp() + 1,
			ext_v.get_last_timestamp() + 1);
	for (i = 0; i < ext_v.get_num_timestamps(); i++) {
		if (range.first <= ext_v.get_timestamps_begin()[i])
			break;
	}
	if (i < ext_v.get_num_timestamps()) {
		new_range.first = ext_v.get_timestamps_begin()[i];
		for (; i < ext_v.get_num_timestamps(); i++) {
			if (ext_v.get_timestamps_begin()[i] >= range.second)
				break;
		}
		if (i < ext_v.get_num_timestamps())
			new_range.second = ext_v.get_timestamps_begin()[i];
		else
			new_range.second = ext_v.get_last_timestamp() + 1;
	}

	off_t begin = ext_v.get_edge_list_offset(new_range.first,
			edge_type::IN_EDGE);
	// The beginning of the timestamp range has to exist.
	assert((size_t) begin != entire_vertex_size);
	off_t end;
	if (new_range.second > ext_v.get_last_timestamp())
		end = entire_vertex_size;
	else {
		end = ext_v.get_edge_list_offset(new_range.second,
				edge_type::IN_EDGE); 
		// The end of the timestamp range also has to exist.
		assert((size_t) end != entire_vertex_size);
	}
	return offset_pair(begin, end);
}
