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

#include "vertex.h"
#include "vertex_index.h"

in_mem_vertex_info::in_mem_vertex_info(vertex_id_t id,
		const vertex_index *index1)
{
	if (index1->get_graph_header().get_graph_type() == graph_type::DIRECTED) {
		directed_vertex_index *index = (directed_vertex_index *) index1;
		this->off = index->get_vertex_off(id);
		this->size = index->get_vertex_size(id);
	}
	else {
		default_vertex_index *index = (default_vertex_index *) index1;
		this->off = index->get_vertex_off(id);
		this->size = index->get_vertex_size(id);
	}
}

ext_mem_directed_vertex::unique_ptr ext_mem_directed_vertex::merge(
		const std::vector<const ext_mem_directed_vertex *> &vertices)
{
	std::vector<vertex_id_t> in_edges;
	std::vector<vertex_id_t> out_edges;
	assert(vertices.size() > 0);
	vertex_id_t id = INVALID_VERTEX_ID;
	for (size_t i = 0; i < vertices.size(); i++) {
		if (vertices[i] == NULL)
			continue;

		assert(!vertices[i]->has_edge_data());
		if (id == INVALID_VERTEX_ID)
			id = vertices[i]->get_id();
		assert(id == vertices[i]->get_id());
		in_edges.insert(in_edges.end(), vertices[i]->neighbors,
				vertices[i]->neighbors + vertices[i]->num_in_edges);
		out_edges.insert(out_edges.end(),
				&vertices[i]->neighbors[vertices[i]->num_in_edges],
				&vertices[i]->neighbors[vertices[i]->num_in_edges
				+ vertices[i]->num_out_edges]);
	}
	assert(id != INVALID_VERTEX_ID);
	std::sort(in_edges.begin(), in_edges.end());
	std::sort(out_edges.begin(), out_edges.end());

	size_t buf_size = get_header_size() + (in_edges.size()
				+ out_edges.size()) * sizeof(vertex_id_t);
	char *vertex_buf = new char[buf_size];
	ext_mem_directed_vertex *out_v = (ext_mem_directed_vertex *) vertex_buf;
	out_v->id = id;
	out_v->edge_data_size = 0;
	out_v->num_in_edges = in_edges.size();
	out_v->num_out_edges = out_edges.size();
	memcpy(out_v->neighbors, in_edges.data(),
			sizeof(vertex_id_t) * in_edges.size());
	memcpy(&out_v->neighbors[out_v->num_in_edges], out_edges.data(),
			sizeof(vertex_id_t) * out_edges.size());
	return unique_ptr(out_v);
}

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
			== (long) edge_list_size);
	assert(this->get_edge_off_begin()[0].in_off == 0);
}

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
