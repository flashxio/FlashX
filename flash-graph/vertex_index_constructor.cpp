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

#include "exception.h"

#include "vertex_index_constructor.h"

namespace fg
{

void cdefault_in_mem_vertex_index::finalize()
{
	if (last_entries.empty())
		return;

	compressed_undirected_vertex_entry centry(last_entries.data(),
			edge_data_size, last_entries.size());
	add_large_vertices(centries.size() * ENTRY_SIZE, centry, last_entries);
	centries.push_back(centry);
	last_entries.clear();
}

void cdefault_in_mem_vertex_index::add_large_vertices(vertex_id_t start_vertex_id,
		const compressed_undirected_vertex_entry entry,
		const std::vector<vertex_offset> &entries)
{
	int num_entries = entries.size() - 1;
	for (int i = 0;i < num_entries; i++) {
		if (entry.is_large_vertex(i)) {
			size_t size = entries[i + 1].get_off() - entries[i].get_off();
			large_vertices.push_back(large_vertex_t(start_vertex_id + i,
					ext_mem_undirected_vertex::vsize2num_edges(
						size, edge_data_size)));
		}
	}
}

void cdefault_in_mem_vertex_index::add_vertex(const in_mem_vertex &v)
{
	assert(get_num_vertices() == v.get_id());
	assert(last_entries.size() <= ENTRY_SIZE + 1);
	assert(!last_entries.empty());
	vertex_offset last_entry = last_entries.back();
	if (last_entries.size() == ENTRY_SIZE + 1) {
		compressed_undirected_vertex_entry centry(last_entries.data(),
				edge_data_size, last_entries.size());
		add_large_vertices(centries.size() * ENTRY_SIZE, centry, last_entries);
		centries.push_back(centry);
		last_entries.clear();
		last_entries.push_back(last_entry);
	}
	vertex_offset entry;
	entry.init(last_entry, v);
	last_entries.push_back(entry);
}

void cdefault_in_mem_vertex_index::dump(const std::string &file,
		const graph_header &header, bool compressed)
{
	throw safs::unsupported_exception();
}

vertex_index::ptr cdefault_in_mem_vertex_index::dump(const graph_header &header,
		bool compressed)
{
	finalize();
	return std::static_pointer_cast<vertex_index>(
			cundirected_vertex_index::construct(*this, header));
}

void cdirected_in_mem_vertex_index::add_large_vertices(
		vertex_id_t start_vertex_id, const compressed_directed_vertex_entry entry,
		const std::vector<directed_vertex_entry> &entries)
{
	int num_entries = entries.size() - 1;
	for (int i = 0;i < num_entries; i++) {
		if (entry.is_large_in_vertex(i)) {
			size_t size = entries[i + 1].get_in_off() - entries[i].get_in_off();
			large_in_vertices.push_back(large_vertex_t(start_vertex_id + i,
					ext_mem_undirected_vertex::vsize2num_edges(
						size, edge_data_size)));
		}
		if (entry.is_large_out_vertex(i)) {
			size_t size = entries[i + 1].get_out_off() - entries[i].get_out_off();
			large_out_vertices.push_back(large_vertex_t(start_vertex_id + i,
					ext_mem_undirected_vertex::vsize2num_edges(
						size, edge_data_size)));
		}
	}
}

void cdirected_in_mem_vertex_index::add_vertex(const in_mem_vertex &v)
{
	assert(get_num_vertices() == v.get_id());
	assert(last_entries.size() <= ENTRY_SIZE + 1);
	assert(!last_entries.empty());
	directed_vertex_entry last_entry = last_entries.back();
	if (last_entries.size() == ENTRY_SIZE + 1) {
		compressed_directed_vertex_entry centry(last_entries.data(),
				edge_data_size, last_entries.size());
		add_large_vertices(centries.size() * ENTRY_SIZE, centry, last_entries);
		centries.push_back(centry);
		last_entries.clear();
		last_entries.push_back(last_entry);
	}
	directed_vertex_entry entry;
	entry.init(last_entry, v);
	last_entries.push_back(entry);
}

void cdirected_in_mem_vertex_index::finalize()
{
	if (last_entries.empty())
		return;

	compressed_directed_vertex_entry centry(last_entries.data(),
			edge_data_size, last_entries.size());
	add_large_vertices(centries.size() * ENTRY_SIZE, centry, last_entries);
	centries.push_back(centry);
	directed_vertex_entry last_entry = last_entries.back();
	last_entries.clear();

	size_t in_part_size = last_entry.get_in_off();
	// If all out-edge lists have been moved to behind in-edge lists,
	// ignore it.
	assert((size_t) centries.front().get_start_out_off() < in_part_size);
	for (size_t i = 0; i < centries.size(); i++) {
		off_t start_in_off = centries[i].get_start_in_off();
		off_t start_out_off = centries[i].get_start_out_off();
		centries[i].reset_start_offs(start_in_off, start_out_off + in_part_size);
	}
}

void cdirected_in_mem_vertex_index::dump(const std::string &file,
		const graph_header &header, bool compressed)
{
	throw safs::unsupported_exception();
}

vertex_index::ptr cdirected_in_mem_vertex_index::dump(const graph_header &header,
		bool compressed)
{
	finalize();
	return std::static_pointer_cast<vertex_index>(
			cdirected_vertex_index::construct(*this, header));
}

}
