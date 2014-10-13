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

#include <boost/log/trivial.hpp>
#include <boost/format.hpp>

#include "io_interface.h"

#include "vertex_compute.h"
#include "vertex_index.h"

static void verify_index(vertex_index::ptr idx)
{
	idx->get_graph_header().verify();
	if (idx->get_graph_header().is_directed_graph()) {
		if (idx->is_compressed())
			cdirected_vertex_index::cast(idx)->verify();
		else
			directed_vertex_index::cast(idx)->verify();
	}
	else {
		assert(!idx->is_compressed());
		default_vertex_index::cast(idx)->verify();
	}
}

size_t vertex_index::get_index_size() const
{
	if (is_compressed() && get_graph_header().is_directed_graph()) {
		return ((cdirected_vertex_index *) this)->cal_index_size();
	}
	else if (!is_compressed() && get_graph_header().is_directed_graph()) {
		return ((directed_vertex_index *) this)->cal_index_size();
	}
	else if (!is_compressed() && !get_graph_header().is_directed_graph()) {
		return ((default_vertex_index *) this)->cal_index_size();
	}
	else
		assert(0);
}

vertex_index::ptr vertex_index::load(const std::string &index_file)
{
	native_file local_f(index_file);
	ssize_t size = local_f.get_size();
	assert(size > 0);
	assert((size_t) size >= sizeof(vertex_index));
	char *buf = (char *) malloc(size);
	assert(buf);
	FILE *fd = fopen(index_file.c_str(), "r");
	size_t ret = fread(buf, size, 1, fd);
	assert(ret == 1);
	fclose(fd);

	vertex_index::ptr idx((vertex_index *) buf, destroy_index());
	assert((size_t) size >= idx->get_index_size());
	verify_index(idx);
	BOOST_LOG_TRIVIAL(info)
		<< boost::format("load vertex index: file size: %1%, index size: %2%")
		% size % idx->get_index_size();

	return idx;
}

vertex_index::ptr vertex_index::safs_load(const std::string &index_file)
{
	const int INDEX_HEADER_SIZE = PAGE_SIZE * 2;
	const int READ_SIZE = 100 * 1024 * 1024;

	// Right now only the cached I/O can support async I/O
	file_io_factory::shared_ptr factory = create_io_factory(index_file,
			REMOTE_ACCESS);
	assert(factory->get_file_size() >= INDEX_HEADER_SIZE);
	io_interface::ptr io = factory->create_io(thread::get_curr_thread());

	// Get the header of the index.
	char *tmp = NULL;
	int ret = posix_memalign((void **) &tmp, PAGE_SIZE, INDEX_HEADER_SIZE);
	assert(ret == 0);
	data_loc_t loc(factory->get_file_id(), 0);
	io_request req(tmp, loc, INDEX_HEADER_SIZE, READ);
	io->access(&req, 1);
	io->wait4complete(1);
	vertex_index *index = (vertex_index *) tmp;
	index->get_graph_header().verify();

	// Initialize the buffer for containing the index.
	size_t index_size = index->get_index_size();
	assert((ssize_t) index_size <= factory->get_file_size());
	char *buf = NULL;
	BOOST_LOG_TRIVIAL(info)
		<< boost::format("allocate %1% bytes for vertex index") % index_size;
	ret = posix_memalign((void **) &buf, PAGE_SIZE,
			std::max(index_size, (size_t) INDEX_HEADER_SIZE));
	assert(ret == 0);
	off_t off = 0;
	memcpy(buf, tmp, INDEX_HEADER_SIZE);
	off += INDEX_HEADER_SIZE;
	free(tmp);

	// Read the index to the memory.
	size_t aligned_index_size = ROUND_PAGE(index_size);
	while ((size_t) off < aligned_index_size) {
		assert(off % PAGE_SIZE == 0);
		size_t size = min(READ_SIZE, aligned_index_size - off);
		data_loc_t loc(factory->get_file_id(), off);
		io_request req(buf + off, loc, size, READ);
		io->access(&req, 1);
		off += size;
		if (io->num_pending_ios() > 100)
			io->wait4complete(io->num_pending_ios() / 10);
	}
	io->wait4complete(io->num_pending_ios());

	// Read the last page.
	// The data may only occupy part of the page.
	if (aligned_index_size < index_size) {
		char *tmp = NULL;
		int ret = posix_memalign((void **) &tmp, PAGE_SIZE, PAGE_SIZE);
		assert(ret == 0);
		data_loc_t loc(factory->get_file_id(), aligned_index_size);
		io_request req(tmp, loc, PAGE_SIZE, READ);
		io->access(&req, 1);
		io->wait4complete(1);
		memcpy(buf + aligned_index_size, tmp, index_size - aligned_index_size);
		free(tmp);
	}

	vertex_index::ptr index_ptr((vertex_index *) buf, destroy_index());
	verify_index(index_ptr);
	return index_ptr;
}

const size_t compressed_directed_vertex_entry::ENTRY_SIZE;
compressed_directed_vertex_entry::compressed_directed_vertex_entry(
		const directed_vertex_entry offs[], size_t edge_data_size, size_t num)
{
	start_offs = offs[0];
	size_t num_vertices = std::min(num - 1, ENTRY_SIZE);
	for (size_t i = 0; i < num_vertices; i++) {
		vsize_t num_in_edges = ext_mem_undirected_vertex::vsize2num_edges(
					offs[i + 1].get_in_off() - offs[i].get_in_off(),
					edge_data_size);
		if (num_in_edges < LARGE_VERTEX_SIZE)
			edges[i].first = num_in_edges;
		else
			edges[i].first = LARGE_VERTEX_SIZE;

		vsize_t num_out_edges = ext_mem_undirected_vertex::vsize2num_edges(
				offs[i + 1].get_out_off() - offs[i].get_out_off(),
				edge_data_size);
		if (num_out_edges < LARGE_VERTEX_SIZE)
			edges[i].second = num_out_edges;
		else
			edges[i].second = LARGE_VERTEX_SIZE;
	}
	for (size_t i = num_vertices; i < ENTRY_SIZE; i++) {
		edges[i].first = edges[i].second = 0;
	}
}

void in_mem_cdirected_vertex_index::init(const directed_vertex_index &index)
{
	BOOST_LOG_TRIVIAL(info) << "init from a regular vertex index";
	index.verify();
	edge_data_size = index.get_graph_header().get_edge_data_size();
	size_t num_entries = index.get_num_entries();
	num_vertices = num_entries - 1;
	entries.resize(ROUNDUP(num_vertices, ENTRY_SIZE) / ENTRY_SIZE);
	for (size_t off = 0; off < num_vertices;
			off += ENTRY_SIZE) {
		off_t entry_idx = off / ENTRY_SIZE;
		entries[entry_idx] = compressed_directed_vertex_entry(
					index.get_data() + off, edge_data_size,
					std::min(ENTRY_SIZE + 1, num_entries - off));

		vertex_id_t id = off;
		for (size_t i = 0; i < ENTRY_SIZE; i++) {
			if (entries[entry_idx].is_large_in_vertex(i)) {
				ext_mem_vertex_info info = index.get_vertex_info_in(id + i);
				large_in_vmap.insert(vertex_map_t::value_type(id + i,
							ext_mem_undirected_vertex::vsize2num_edges(
								info.get_size(), edge_data_size)));
			}
			if (entries[entry_idx].is_large_out_vertex(i)) {
				ext_mem_vertex_info info = index.get_vertex_info_out(id + i);
				large_out_vmap.insert(vertex_map_t::value_type(id + i,
							ext_mem_undirected_vertex::vsize2num_edges(
								info.get_size(), edge_data_size)));
			}
		}
	}
}

void in_mem_cdirected_vertex_index::init(const cdirected_vertex_index &index)
{
	struct timeval start, end;
	gettimeofday(&start, NULL);
	BOOST_LOG_TRIVIAL(info) << "init from a compressed vertex index";
	index.verify();
	edge_data_size = index.get_graph_header().get_edge_data_size();
	num_vertices = index.get_graph_header().get_num_vertices();
	entries.insert(entries.end(), index.get_entries(),
			index.get_entries() + index.get_num_entries());

	const large_vertex_t *l_in_vertex_array = index.get_large_in_vertices();
	size_t num_large_in_vertices = index.get_num_large_in_vertices();
	const large_vertex_t *l_out_vertex_array = index.get_large_out_vertices();
	size_t num_large_out_vertices = index.get_num_large_out_vertices();

	for (size_t i = 0; i < num_large_in_vertices; i++) {
		large_in_vmap.insert(l_in_vertex_array[i]);
	}

	for (size_t i = 0; i < num_large_out_vertices; i++) {
		large_out_vmap.insert(l_out_vertex_array[i]);
	}
	BOOST_LOG_TRIVIAL(info)
		<< boost::format("There are %1% large in-vertices and %2% large out-vertices")
		% num_large_in_vertices % num_large_out_vertices;
	gettimeofday(&end, NULL);
	BOOST_LOG_TRIVIAL(info)
		<< boost::format("init in-mem compressed index takes %1% seconds")
		% time_diff(start, end);
}

in_mem_cdirected_vertex_index::in_mem_cdirected_vertex_index(
		vertex_index &index)
{
	if (index.is_compressed())
		init((const cdirected_vertex_index &) index);
	else
		init((const directed_vertex_index &) index);
}

const size_t in_mem_cdirected_vertex_index::ENTRY_SIZE;
directed_vertex_entry in_mem_cdirected_vertex_index::get_vertex(
		vertex_id_t id) const
{
	directed_vertex_entry e = entries[id / ENTRY_SIZE].get_start_offs();
	int off = id % ENTRY_SIZE;
	off_t in_off = e.get_in_off();
	off_t out_off = e.get_out_off();
	vertex_id_t start_id = id & (~ENTRY_MASK);
	for (int i = 0; i < off; i++) {
		in_off += get_in_size(start_id + i);
		out_off += get_out_size(start_id + i);
	}
	return directed_vertex_entry(in_off, out_off);
}

#include "vertex_index_reader.h"

void in_mem_cdirected_vertex_index::verify_against(
		directed_vertex_index &index)
{
	index.verify();
	id_range_t range(10, std::min(100UL, index.get_num_vertices()));
	compressed_directed_index_iterator it(*this, range);
	vertex_id_t id = range.first;
	while (it.has_next()) {
		ext_mem_vertex_info in_info = index.get_vertex_info_in(id);
		ext_mem_vertex_info out_info = index.get_vertex_info_out(id);
		assert(in_info.get_off() == it.get_curr_off());
		assert(out_info.get_off() == it.get_curr_out_off());
		id++;
		it.move_next();
	}
}

cdirected_vertex_index::ptr cdirected_vertex_index::construct(
		directed_vertex_index &index)
{
	size_t edge_data_size = index.get_graph_header().get_edge_data_size();
	size_t num_entries = index.get_num_entries();
	size_t num_vertices = num_entries - 1;
	std::vector<large_vertex_t> large_in_vertices;
	std::vector<large_vertex_t> large_out_vertices;
	std::vector<compressed_directed_vertex_entry> entries(
			ROUNDUP(num_vertices, ENTRY_SIZE) / ENTRY_SIZE);
	for (size_t off = 0; off < num_vertices; off += ENTRY_SIZE) {
		off_t entry_idx = off / ENTRY_SIZE;
		entries[entry_idx] = compressed_directed_vertex_entry(
					index.get_data() + off, edge_data_size,
					std::min(ENTRY_SIZE + 1, num_entries - off));

		vertex_id_t id = off;
		for (size_t i = 0; i < ENTRY_SIZE; i++) {
			if (entries[entry_idx].is_large_in_vertex(i)) {
				ext_mem_vertex_info info = index.get_vertex_info_in(id + i);
				large_in_vertices.push_back(large_vertex_t(id + i,
							ext_mem_undirected_vertex::vsize2num_edges(
								info.get_size(), edge_data_size)));
			}
			if (entries[entry_idx].is_large_out_vertex(i)) {
				ext_mem_vertex_info info = index.get_vertex_info_out(id + i);
				large_out_vertices.push_back(large_vertex_t(id + i,
							ext_mem_undirected_vertex::vsize2num_edges(
								info.get_size(), edge_data_size)));
			}
		}
	}

	size_t tot_size = sizeof(cdirected_vertex_index)
		+ sizeof(entries[0]) * entries.size()
		+ sizeof(large_in_vertices[0]) * large_in_vertices.size()
		+ sizeof(large_out_vertices[0]) * large_out_vertices.size();
	char *buf = (char *) malloc(tot_size);
	memcpy(buf, &index, vertex_index::get_header_size());
	cdirected_vertex_index *cindex = (cdirected_vertex_index *) buf;
	cindex->h.data.entry_size = sizeof(entries[0]);
	cindex->h.data.num_entries = entries.size();
	cindex->h.data.compressed = true;
	cindex->h.data.num_large_in_vertices = large_in_vertices.size();
	cindex->h.data.num_large_out_vertices = large_out_vertices.size();

	memcpy(cindex->entries, entries.data(), entries.size() * sizeof(entries[0]));
	memcpy(cindex->get_large_in_vertices(), large_in_vertices.data(),
			sizeof(large_in_vertices[0]) * large_in_vertices.size());
	memcpy(cindex->get_large_out_vertices(), large_out_vertices.data(),
			sizeof(large_out_vertices[0]) * large_out_vertices.size());
	return ptr((cdirected_vertex_index *) buf, destroy_index());
}
