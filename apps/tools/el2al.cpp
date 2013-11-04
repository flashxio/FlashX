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

#include "graph.h"

#ifdef MEMCHECK
#include <algorithm>
#else
#include <parallel/algorithm>
#endif

#include "thread.h"
#include "native_file.h"

const int NUM_THREADS = 32;
const int NUM_NODES = 4;
const int EDGE_LIST_BLOCK_SIZE = 1 * 1024 * 1024;

directed_graph *par_load_edge_list_text(const std::string &file);

static struct comp_edge {
	bool operator() (const edge &e1, const edge &e2) {
		if (e1.get_from() == e2.get_from())
			return e1.get_to() < e2.get_to();
		else
			return e1.get_from() < e2.get_from();
	}
} edge_comparator;

static struct comp_in_edge {
	bool operator() (const edge &e1, const edge &e2) {
		if (e1.get_to() == e2.get_to())
			return e1.get_from() < e2.get_from();
		else
			return e1.get_to() < e2.get_to();
	}
} in_edge_comparator;

class directed_edge_graph
{
	std::vector<edge> in_edges;
	std::vector<edge> out_edges;
	pthread_mutex_t lock;
public:
	directed_edge_graph() {
		pthread_mutex_init(&lock, NULL);
	}
	void sort_edges();
	void add_edges(std::vector<edge> &edges);
	directed_graph *create() const;
};

class graph_file_io
{
	FILE *f;
	ssize_t file_size;
public:
	graph_file_io(const std::string file) {
		f = fopen(file.c_str(), "r");
		if (f == NULL) {
			perror("fopen");
			assert(0);
		}
		native_file local_f(file);
		file_size = local_f.get_size();
	}

	~graph_file_io() {
		if (f)
			fclose(f);
	}

	char *read_edge_list_text(const size_t wanted_bytes, size_t &read_bytes);

	size_t get_num_remaining_bytes() const {
		off_t curr_off = ftell(f);
		return file_size - curr_off;
	}
};

/**
 * It read a text of an edge list roughly the size of the wanted bytes.
 * The returned text may be a little more than the wanted bytes, but
 * it's guaranteed that all lines are complete.
 * The returned string ends with '\0'.
 */
char *graph_file_io::read_edge_list_text(const size_t wanted_bytes,
		size_t &read_bytes)
{
	off_t curr_off = ftell(f);
	off_t off = curr_off + wanted_bytes;
	// After we just to the new location, we need to further read another
	// page to search for the end of a line. If there isn't enough data,
	// we can just read all remaining data.
	if (off + PAGE_SIZE < file_size) {
		int ret = fseek(f, off, SEEK_SET);
		if (ret < 0) {
			perror("fseek");
			return NULL;
		}

		char buf[PAGE_SIZE];
		ret = fread(buf, sizeof(buf), 1, f);
		if (ret != 1) {
			perror("fread");
			return NULL;
		}
		unsigned i;
		for (i = 0; i < sizeof(buf); i++)
			if (buf[i] == '\n')
				break;
		// A line shouldn't be longer than a page.
		assert(i != sizeof(buf));

		// We read a little more than asked to make sure that we read
		// the entire line.
		read_bytes = wanted_bytes + i + 1;

		// Go back to the original offset in the file.
		ret = fseek(f, curr_off, SEEK_SET);
		assert(ret == 0);
	}
	else {
		read_bytes = file_size - curr_off;
	}

	// The line buffer must end with '\0'.
	char *line_buf = new char[read_bytes + 1];
	int ret = fread(line_buf, read_bytes, 1, f);
	assert(ret == 1);
	line_buf[read_bytes] = 0;

	return line_buf;
}

int parse_edge_list_line(char *line, edge &e)
{
	if (line[0] == '#')
		return 0;
	char *second = strchr(line, '\t');
	if (second == NULL) {
		fprintf(stderr, "wrong format 1: %s\n", line);
		return -1;
	}
	*second = 0;
	second++;
	if (!isnumeric(line) || !isnumeric(second)) {
		fprintf(stderr, "wrong format 2: %s\t%s\n", line, second);
		return -1;
	}
	vertex_id_t from = atol(line);
	vertex_id_t to = atol(second);
	e = edge(from, to);

	return 1;
}

/**
 * Parse the edge list in the character buffer.
 * `size' doesn't include '\0'.
 */
size_t parse_edge_list_text(char *line_buf, size_t size,
		std::vector<edge> &edges)
{
	char *line_end;
	char *line = line_buf;
	size_t num_edges = 0;
	while ((line_end = strchr(line, '\n'))) {
		*line_end = 0;
		if (*(line_end - 1) == '\r')
			*(line_end - 1) = 0;
		edge e;
		int num = parse_edge_list_line(line, e);
		if (num > 0)
			edges.push_back(e);
		num_edges += num;
		line = line_end + 1;
		assert(line - line_end <= (ssize_t) size);
	}
	if (line - line_buf < (ssize_t) size) {
		edge e;
		int num = parse_edge_list_line(line, e);
		if (num > 0)
			edges.push_back(e);
		num_edges += num;
	}
	return num_edges;
}

class text_edge_task: public thread_task
{
	char *line_buf;
	size_t size;
	directed_edge_graph *graph;
public:
	text_edge_task(directed_edge_graph *graph, char *line_buf, size_t size) {
		this->graph = graph;
		this->line_buf = line_buf;
		this->size = size;
	}

	void run();
};

void text_edge_task::run()
{
	std::vector<edge> edges;
	parse_edge_list_text(line_buf, size, edges);
	std::vector<edge> *local_edge_buf
		= (std::vector<edge> *) thread::get_curr_thread()->get_user_data();
	local_edge_buf->insert(local_edge_buf->end(), edges.begin(), edges.end());
	delete [] line_buf;
}

void directed_edge_graph::add_edges(std::vector<edge> &edges)
{
	pthread_mutex_lock(&lock);
	in_edges.insert(in_edges.end(), edges.begin(), edges.end());
	out_edges.insert(out_edges.end(), edges.begin(), edges.end());
	pthread_mutex_unlock(&lock);
}

void directed_edge_graph::sort_edges()
{
#ifdef MEMCHECK
	std::sort(out_edges.begin(), out_edges.end(), edge_comparator);
	std::sort(in_edges.begin(), in_edges.end(), in_edge_comparator);
#else
	__gnu_parallel::sort(out_edges.begin(), out_edges.end(), edge_comparator);
	__gnu_parallel::sort(in_edges.begin(), in_edges.end(), in_edge_comparator);
#endif
}

directed_graph *directed_edge_graph::create() const
{
	directed_graph *g = new directed_graph();

	vertex_id_t curr = 0;
	in_mem_directed_vertex v(curr);
	size_t out_idx = 0;
	size_t in_idx = 0;
	size_t num_edges = in_edges.size();
	assert(in_edges.size() == out_edges.size());
	while (out_idx < num_edges && in_idx < num_edges) {
		while (out_edges[out_idx].get_from() == curr
				&& out_idx < num_edges) {
			v.add_out_edge(out_edges[out_idx++].get_to());
		}
		while (in_edges[in_idx].get_to() == curr
				&& in_idx < num_edges) {
			v.add_in_edge(in_edges[in_idx++].get_from());
		}
		g->add_vertex(v);
		vertex_id_t prev = curr + 1;
		if (out_idx < num_edges && in_idx < num_edges)
			curr = min(out_edges[out_idx].get_from(),
					in_edges[in_idx].get_to());
		else if (out_idx < num_edges)
			curr = out_edges[out_idx].get_from();
		else if (in_idx < num_edges)
			curr = in_edges[in_idx].get_to();
		else
			break;
		// The vertices without edges won't show up in the edge list,
		// but we need to fill the gap in the vertex Id space with empty
		// vertices.
		while (prev < curr) {
			v = in_mem_directed_vertex(prev);
			prev++;
			g->add_vertex(v);
		}
		v = in_mem_directed_vertex(curr);
	}

	assert(g->get_num_in_edges() == num_edges);
	assert(g->get_num_out_edges() == num_edges);
	return g;
}

directed_graph *par_load_edge_list_text(const std::string &file)
{
	struct timeval start, end;
	gettimeofday(&start, NULL);
	std::vector<task_thread *> threads(NUM_THREADS);
	for (int i = 0; i < NUM_THREADS; i++) {
		task_thread *t = new task_thread(std::string(
					"graph-task-thread") + itoa(i), i % NUM_NODES);
		t->set_user_data(new std::vector<edge>());
		t->start();
		threads[i] = t;
	}
	graph_file_io io(file);
	directed_edge_graph edge_g;
	int thread_no = 0;
	while (io.get_num_remaining_bytes() > 0) {
		size_t size = 0;
		char *line_buf = io.read_edge_list_text(EDGE_LIST_BLOCK_SIZE, size);
		assert(line_buf);
		thread_task *task = new text_edge_task(&edge_g, line_buf, size);
		threads[thread_no % NUM_THREADS]->add_task(task);
		thread_no++;
	}
	for (int i = 0; i < NUM_THREADS; i++)
		threads[i]->wait4complete();
	gettimeofday(&end, NULL);
	printf("It takes %f seconds to construct edge list\n", time_diff(start, end));

	start = end;
	for (int i = 0; i < NUM_THREADS; i++) {
		std::vector<edge> *local_edges = (std::vector<edge> *) threads[i]->get_user_data();
		edge_g.add_edges(*local_edges);
		delete local_edges;
	}
	gettimeofday(&end, NULL);
	printf("It takes %f seconds to combine edge list\n", time_diff(start, end));

	start = end;
	edge_g.sort_edges();
	gettimeofday(&end, NULL);
	printf("It takes %f seconds to sort edge list\n", time_diff(start, end));

	start = end;
	directed_graph *g = edge_g.create();
	gettimeofday(&end, NULL);
	printf("It takes %f seconds to construct the graph\n", time_diff(start, end));
	return g;
}

int main(int argc, char *argv[])
{
	if (argc < 5) {
		fprintf(stderr,
				"el2al edge_list_file adjacency_list_file index_file directed\n");
		exit(-1);
	}

	const std::string edge_list_file = argv[1];
	const std::string adjacency_list_file = argv[2];
	const std::string index_file = argv[3];
	bool directed = atoi(argv[4]) != 0;

	if (directed) {
		struct timeval start, end;
		directed_graph *g = par_load_edge_list_text(
				edge_list_file);
		gettimeofday(&start, NULL);
		vertex_index *index = g->create_vertex_index();
		gettimeofday(&end, NULL);
		printf("It takes %f seconds to create vertex index\n",
				time_diff(start, end));

		start = end;
		g->dump(adjacency_list_file);
		gettimeofday(&end, NULL);
		printf("It takes %f seconds to dump adjacency list\n",
				time_diff(start, end));

		start = end;
		index->dump(index_file);
		gettimeofday(&end, NULL);
		printf("It takes %f seconds to dump index\n",
				time_diff(start, end));
		printf("There are %ld vertices, %ld non-empty vertices and %ld edges\n",
				g->get_num_vertices(), g->get_num_non_empty_vertices(),
				g->get_num_in_edges());
		directed_graph::destroy(g);
	}
	else {
		undirected_graph *g = undirected_graph::load_edge_list_text(
				edge_list_file);
		vertex_index *index = g->create_vertex_index();
		g->dump(adjacency_list_file);
		index->dump(index_file);
		printf("There are %ld vertices, %ld non-empty vertices and %ld edges\n",
				g->get_num_vertices(), g->get_num_non_empty_vertices(),
				g->get_num_edges());
		undirected_graph::destroy(g);
	}
}
