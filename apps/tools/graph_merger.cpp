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

#include <unistd.h>

#include "graph.h"

bool print_graph = false;
bool check_graph = false;

const int MEM_BUF_SIZE = 1024 * 1024 * 100;

class vertex_buffer
{
	vertex_index *index;
	char *buf;
	// The offset of the data in the external memory
	off_t off;
	// THe size of the data.
	size_t size;
public:
	vertex_buffer(vertex_index *index, off_t off, char *buf, size_t size) {
		this->index = index;
		this->off = off;
		this->buf = buf;
		this->size = size;
	}

	const ext_mem_directed_vertex *get_directed_vertex(vertex_id_t id) const {
		assert(is_valid());
		if (id >= index->get_num_vertices())
			return NULL;

		off_t v_off = get_vertex_off(index, id);
		size_t v_size = get_vertex_size(index, id);
		assert(off <= v_off && (size_t) v_off < off + size);
		assert(v_off + v_size <= off + size);
		ext_mem_directed_vertex *v
			= (ext_mem_directed_vertex *) (buf + (v_off - off));
		assert(v->get_id() == id);
		assert(v->get_size() == v_size);
		return v;
	}

	vertex_buffer get_vertex(vertex_id_t id) const {
		assert(is_valid());
		if (id >= index->get_num_vertices())
			return vertex_buffer(NULL, -1, NULL, 0);

		off_t v_off = get_vertex_off(index, id);
		size_t v_size = get_vertex_size(index, id);
		assert(off <= v_off && (size_t) v_off < off + size);
		assert(v_off + v_size <= off + size);
		return vertex_buffer(index, v_off, buf + (v_off - off), v_size);
	}

	size_t get_size() const {
		return size;
	}

	bool is_valid() const {
		return index != NULL;
	}
};

class graph_reader
{
	char *buf;
	vertex_index *index;
	FILE *graph_fd;
public:
	graph_reader() {
		buf = NULL;
		index = NULL;
		graph_fd = NULL;
	}

	void init(const std::string &graph_file, vertex_index *index) {
		buf = new char[MEM_BUF_SIZE];
		this->index = index;
		graph_fd = fopen(graph_file.c_str(), "r");
		assert(graph_fd);
	}

	~graph_reader() {
		if (graph_fd)
			fclose(graph_fd);
		if (buf)
			delete [] buf;
	}

	/**
	 * Get the vertex interval that can be buffered in the pre-allocated
	 * memory.
	 * It takes the lower end of the interval and returns the higher end
	 * of the interval. The higher end is excluded from the interval.
	 * Once we reach the end of the graph, it returns the maximal vertex ID.
	 */
	vertex_id_t get_buffered_interval(vertex_id_t since);

	vertex_buffer read_interval(
			const std::pair<vertex_id_t, vertex_id_t> interval);
};

vertex_id_t graph_reader::get_buffered_interval(vertex_id_t since)
{
	size_t size = 0;
	vertex_id_t id = since;
	while (size <= (size_t) MEM_BUF_SIZE && id <= index->get_max_id()) {
		size += get_vertex_size(index, id);
		id++;
	}
	if (size > (size_t) MEM_BUF_SIZE) {
		// The upper boundary is exclusive.
		id -= 1;
	}
	else {
		id = MAX_VERTEX_ID;
	}
	assert(id > since);
	return id;
}

vertex_buffer graph_reader::read_interval(
		const std::pair<vertex_id_t, vertex_id_t> interval)
{
	off_t start_off = get_vertex_off(index, interval.first);
	off_t end_off;
	if (interval.second < index->get_num_vertices())
		end_off = get_vertex_off(index, interval.second);
	else
		end_off = index->get_graph_size();
	size_t size = end_off - start_off;
	assert(size <= (size_t) MEM_BUF_SIZE);
	int seek_ret = fseek(graph_fd, start_off, SEEK_SET);
	assert(seek_ret == 0);
	size_t ret = fread(buf, size, 1, graph_fd);
	assert(ret == 1);
	return vertex_buffer(index, start_off, buf, size);
}

class graph_stat
{
	size_t num_non_empty;
	size_t num_vertices;
	size_t num_edges;
public:
	graph_stat() {
		num_vertices = 0;
		num_edges = 0;
	}

	void inc_num_non_empty_vertices(size_t num) {
		num_non_empty += num;
	}

	void inc_num_vertices(size_t num) {
		num_vertices += num;
	}

	void inc_num_edges(size_t num) {
		num_edges += num;
	}

	size_t get_num_non_empty_vertices() const {
		return num_non_empty;
	}

	size_t get_num_vertices() const {
		return num_vertices;
	}

	size_t get_num_edges() const {
		return num_edges;
	}
};

size_t merge_directed_vertex(vertex_id_t id,
		const std::vector<vertex_buffer> &buffers,
		graph_stat &stat, char *vertex_buf, size_t buf_size)
{
	std::vector<const ext_mem_directed_vertex *> vertices;
	for (size_t i = 0; i < buffers.size(); i++) {
		if (buffers[i].is_valid())
			vertices.push_back(buffers[i].get_directed_vertex(id));
	}
	ext_mem_directed_vertex *v = ext_mem_directed_vertex::merge(
			vertices, vertex_buf, buf_size);
	int num_edges = v->get_num_in_edges() + v->get_num_out_edges();
	if (num_edges > 0) {
		stat.inc_num_non_empty_vertices(1);
		stat.inc_num_edges(num_edges);
	}
	size_t ret = v->get_size();
	assert(ret <= buf_size);
	return ret;
}

void merge_dump_part(const std::vector<vertex_buffer> &buffers,
		std::pair<vertex_id_t, vertex_id_t> interval,
		in_mem_vertex_index &in_mem_index, graph_type type,
		graph_stat &stat, FILE *adjacency_fd)
{
	for (vertex_id_t id = interval.first; id < interval.second; id++) {
		std::vector<vertex_buffer> vertices;
		size_t buf_size = 0;
		for (size_t i = 0; i < buffers.size(); i++) {
			vertex_buffer buf = buffers[i].get_vertex(id);
			buf_size += buf.get_size();
			vertices.push_back(buf);
		}
		char *merged_vertex_buf = new char[buf_size];
		stat.inc_num_vertices(1);
		size_t used_size = 0;
		switch(type) {
			case graph_type::DIRECTED:
				used_size = merge_directed_vertex(id, vertices, stat,
						merged_vertex_buf, buf_size);
				break;
			case graph_type::TS_DIRECTED:
			case graph_type::TS_UNDIRECTED:
			case graph_type::UNDIRECTED:
			default:
				assert(0);
		}
		assert(used_size <= buf_size);
		in_mem_index.add_vertex(merged_vertex_buf);
		size_t ret = fwrite(merged_vertex_buf, used_size, 1, adjacency_fd);
		assert(ret);
		delete [] merged_vertex_buf;
	}
}

void merge_dump(const std::vector<std::string> &graph_files,
		const std::vector<vertex_index *> &indices, graph_type type,
		const std::string &adjacency_list_file, const std::string &index_file)
{
	assert(graph_files.size() > 0);
	assert(graph_files.size() == indices.size());
	std::vector<graph_reader> graphs(graph_files.size());
	for (size_t i = 0; i < graph_files.size(); i++)
		graphs[i].init(graph_files[i], indices[i]);

	FILE *adjacency_fd = fopen(adjacency_list_file.c_str(), "w");
	assert(adjacency_fd);
	// We don't know enough information about the graph, let's just
	// write some dump data to occupy the space for the graph header.
	graph_header header;
	size_t ret = fwrite(&header, sizeof(header), 1, adjacency_fd);
	assert(ret == 1);

	in_mem_vertex_index *in_mem_index;
	if (type == graph_type::DIRECTED)
		in_mem_index = new directed_in_mem_vertex_index();
	else
		assert(0);

	// Find the number of vertices in the merged graph.
	size_t num_vertices = indices[0]->get_num_vertices();
	for (size_t i = 0; i < indices.size(); i++)
		num_vertices = max(num_vertices, indices[i]->get_num_vertices());

	graph_stat stat;
	vertex_id_t since = 0;
	while (true) {
		vertex_id_t higher = graphs[0].get_buffered_interval(since);
		for (size_t i = 1; i < graphs.size(); i++)
			higher = min(higher, graphs[i].get_buffered_interval(since));
		if (higher == MAX_VERTEX_ID)
			higher = num_vertices;
		// If no graphs have vertices left, we have merged all vertices
		// in the graphs.
		if (since == higher)
			break;

		std::vector<vertex_buffer> buffers;
		for (size_t i = 0; i < graphs.size(); i++)
			buffers.push_back(graphs[i].read_interval(
						std::pair<vertex_id_t, vertex_id_t>(since, higher)));
		std::pair<vertex_id_t, vertex_id_t> interval(since, higher);
		merge_dump_part(buffers, interval, *in_mem_index, type, stat,
				adjacency_fd);
		since = higher;
	}

	assert(indices.size() > 0);
	bool has_edge_data = indices[0]->get_graph_header().has_edge_data();
	for (size_t i = 1; i < indices.size(); i++)
		assert(indices[i]->get_graph_header().has_edge_data() == has_edge_data);

	int num_timestamps = 0;
	if (type == graph_type::TS_DIRECTED || type == graph_type::TS_UNDIRECTED)
		num_timestamps = graph_files.size();
	header = graph_header(type, stat.get_num_vertices(), stat.get_num_edges(),
			has_edge_data, num_timestamps);
	int seek_ret = fseek(adjacency_fd, 0, SEEK_SET);
	assert(seek_ret == 0);
	ret = fwrite(&header, sizeof(header), 1, adjacency_fd);
	assert(ret == 1);
	fclose(adjacency_fd);

	in_mem_index->dump(index_file, header);
	delete in_mem_index;

	printf("There are %ld vertices, %ld non-empty vertices and %ld edges\n",
			stat.get_num_vertices(), stat.get_num_non_empty_vertices(),
			stat.get_num_edges());
}

void print_usage()
{
	fprintf(stderr, "merge multiple graphs to a single graph\n");
	fprintf(stderr,
			"graph_merger [options] adj_list_file index_file graph_files index_files\n");
	fprintf(stderr, "-p: print adjacency list\n");
	fprintf(stderr, "-v: verify the created adjacency list\n");
	fprintf(stderr, "-t: merge graphs into a time-series graph. \n");
}

int main(int argc, char *argv[])
{
	int opt;
	int num_opts = 0;
	bool ts_merge = false;
	while ((opt = getopt(argc, argv, "pvt")) != -1) {
		num_opts++;
		switch (opt) {
			case 'p':
				print_graph = true;
				break;
			case 'v':
				check_graph = true;
				break;
			case 't':
				ts_merge = true;
				break;
			default:
				print_usage();
		}
	}
	argv += 1 + num_opts;
	argc -= 1 + num_opts;
	if (argc < 4 || argc % 2 != 0) {
		print_usage();
		exit(-1);
	}

	std::string out_graph_file = argv[0];
	out_graph_file += std::string("-v") + itoa(CURR_VERSION);
	std::string out_index_file = argv[1];
	out_index_file += std::string("-v") + itoa(CURR_VERSION);
	std::vector<std::string> in_graph_files;
	std::vector<std::string> in_index_files;
	int num_graphs = (argc - 2) / 2;
	for (int i = 2; i < num_graphs + 2; i++)
		in_graph_files.push_back(argv[i]);
	for (int i = num_graphs + 2; i < argc; i++)
		in_index_files.push_back(argv[i]);

	std::vector<vertex_index *> indices;
	for (size_t i = 0; i < in_index_files.size(); i++)
		indices.push_back(vertex_index::load(in_index_files[i]));
	bool directed = indices[0]->get_graph_header().is_directed_graph();
	graph_type in_type = indices[0]->get_graph_header().get_graph_type();
	// We don't want to merge two time-series graphs.
	assert(in_type == graph_type::DIRECTED
			|| in_type == graph_type::UNDIRECTED);
	for (size_t i = 0; i < in_index_files.size(); i++) {
		assert(indices[i]->get_graph_header().is_directed_graph() == directed);
		in_type = indices[i]->get_graph_header().get_graph_type();
		assert(in_type == graph_type::DIRECTED
				|| in_type == graph_type::UNDIRECTED);
	}

	graph_type out_type;
	if (ts_merge && directed)
		out_type = graph_type::TS_DIRECTED;
	else if (ts_merge && !directed)
		out_type = graph_type::TS_UNDIRECTED;
	else if (directed)
		out_type = graph_type::DIRECTED;
	else
		out_type = graph_type::UNDIRECTED;

	// If we get more than one graph, we need to merge them and will have
	// a time-series graph.
	if (indices.size() > 1) {
		struct timeval start, end;
		gettimeofday(&start, NULL);
		merge_dump(in_graph_files, indices, out_type, out_graph_file,
				out_index_file);
		gettimeofday(&end, NULL);
		printf("It takes %f seconds to merge graphs\n",
				time_diff(start, end));
	}
}
