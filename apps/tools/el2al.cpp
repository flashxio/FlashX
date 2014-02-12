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
#include <unistd.h>

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
const char *delimiter = "\t";

bool compress = false;
bool simplfy = false;
bool print_graph = false;
bool check_graph = false;

template<class edge_data_type>
struct comp_edge {
	bool operator() (const edge<edge_data_type> &e1, const edge<edge_data_type> &e2) {
		if (e1.get_from() == e2.get_from())
			return e1.get_to() < e2.get_to();
		else
			return e1.get_from() < e2.get_from();
	}
};

template<class edge_data_type>
struct comp_in_edge {
	bool operator() (const edge<edge_data_type> &e1, const edge<edge_data_type> &e2) {
		if (e1.get_to() == e2.get_to())
			return e1.get_from() < e2.get_from();
		else
			return e1.get_to() < e2.get_to();
	}
};

template<class edge_data_type>
class directed_edge_graph;

/**
 * This is a disk-backed directed graph.
 */
template<class edge_data_type>
class disk_directed_graph: public graph
{
	size_t num_edges;
	size_t num_vertices;
	size_t num_non_empty;
	FILE *f;
	const directed_edge_graph<edge_data_type> *g;
	directed_in_mem_vertex_index index;
public:
	disk_directed_graph(const directed_edge_graph<edge_data_type> *g) {
		num_edges = 0;
		num_vertices = 0;
		num_non_empty = 0;
		this->g = g;
		f = NULL;
	}

	~disk_directed_graph() {
		delete g;
	}

	void add_vertex(const in_mem_vertex &v1) {
		assert(f);
		const in_mem_directed_vertex<edge_data_type> &v
			= (const in_mem_directed_vertex<edge_data_type> &) v1;
		num_vertices++;
		// To get the total number of edges, I only accumulate on in-edges
		// or out-edges.
		num_edges += v.get_num_edges(edge_type::IN_EDGE);
		if (v.get_num_edges(edge_type::BOTH_EDGES) > 0)
			num_non_empty++;

		int mem_size = v.get_serialize_size();
		char *buf = new char[mem_size];
		ext_mem_directed_vertex::serialize<edge_data_type>(v,
				buf, mem_size);
		index.add_vertex(buf);
		ssize_t ret = fwrite(buf, mem_size, 1, f);
		delete [] buf;
		assert(ret == 1);
	}

	virtual size_t get_num_edges() const {
		return num_edges;
	}

	virtual size_t get_num_vertices() const {
		return num_vertices;
	}

	virtual bool has_edge_data() const;

	virtual size_t get_num_non_empty_vertices() const {
		return num_non_empty;
	}

	virtual void dump(const std::string &index_file,
			const std::string &graph_file);

	virtual void check_ext_graph(const std::string &index_file,
			const std::string &adj_file) const;

	virtual vertex_index *create_vertex_index() const {
		assert(0);
	}

	void get_all_vertices(std::vector<vertex_id_t> &ids) const {
		assert(0);
	}

	virtual void print() const {
		assert(0);
	}
};

/**
 * This represents a directed graph in the form of edge list.
 * It maintains a sorted list of out-edges (sorted on the from vertices)
 * and a sorted list of in-edges (sorted on the to vertices).
 */
template<class edge_data_type = empty_data>
class directed_edge_graph
{
	bool has_data;
	std::vector<edge<edge_data_type> > in_edges;
	std::vector<edge<edge_data_type> > out_edges;
	pthread_mutex_t lock;
public:
	directed_edge_graph(bool has_data) {
		this->has_data = has_data;
		pthread_mutex_init(&lock, NULL);
	}

	/**
	 * num_edges tells the edge graph that there will be num_edges
	 * edges added to the graph.
	 */
	directed_edge_graph(size_t num_edges, bool has_data) {
		this->has_data = has_data;
		pthread_mutex_init(&lock, NULL);
		in_edges.reserve(num_edges);
		out_edges.reserve(num_edges);
	}
	void sort_edges();
	directed_edge_graph<edge_count> *compress_edges() const;
	directed_edge_graph<edge_data_type> *simplify_edges() const;
	void construct_graph(graph *g) const;

	bool has_edge_data() const {
		return has_data;
	}

	directed_graph<edge_data_type> *create_in_mem_graph() const {
		directed_graph<edge_data_type> *g
			= new directed_graph<edge_data_type>(has_data);
		construct_graph(g);
		return g;
	}

	disk_directed_graph<edge_data_type> *create_disk_graph() const {
		return new disk_directed_graph<edge_data_type>(this);
	}

	void add_edge(const edge<edge_data_type> &e) {
		in_edges.push_back(e);
		out_edges.push_back(e);
	}

	void add_edges(std::vector<edge<edge_data_type> > &edges) {
		pthread_mutex_lock(&lock);
		in_edges.insert(in_edges.end(), edges.begin(), edges.end());
		out_edges.insert(out_edges.end(), edges.begin(), edges.end());
		printf("the edge graph has %ld edges and use %ld bytes\n",
				in_edges.size() + out_edges.size(), (in_edges.capacity()
					+ out_edges.capacity()) * sizeof(edge<edge_data_type>));
		pthread_mutex_unlock(&lock);
	}

	size_t get_num_edges() const {
		assert(in_edges.size() == out_edges.size());
		return in_edges.size();
	}
};

template<class edge_data_type>
void disk_directed_graph<edge_data_type>::check_ext_graph(
		const std::string &index_file, const std::string &adj_file) const
{
	directed_graph<edge_data_type> *dg = g->create_in_mem_graph();
	dg->check_ext_graph(index_file, adj_file);
	delete g;
}

template<class edge_data_type>
bool disk_directed_graph<edge_data_type>::has_edge_data() const
{
	return g->has_edge_data();
}

template<class edge_data_type>
void disk_directed_graph<edge_data_type>::dump(
		const std::string &index_file, const std::string &graph_file)
{
	assert(g);
	f = fopen(graph_file.c_str(), "w");
	graph_header header;
	// Write a dumb header to the file to occupy the space.
	ssize_t ret = fwrite(&header, sizeof(header), 1, f);
	assert(ret);

	// Write the adjacency lists to the graph file.
	g->construct_graph(this);
	assert(g->get_num_edges() == get_num_edges());

	// Write the real graph header.
	header = graph_header(graph_type::DIRECTED, num_vertices,
			num_edges, g->has_edge_data());
	int seek_ret = fseek(f, 0, SEEK_SET);
	assert(seek_ret == 0);
	ret = fwrite(&header, sizeof(header), 1, f);
	assert(ret == 1);
	fclose(f);
	f = NULL;

	assert(this->get_num_edges() == g->get_num_edges());
	index.dump(index_file, header);
}

template<class edge_data_type = empty_data>
directed_edge_graph<edge_data_type> *par_load_edge_list_text(
		const std::vector<std::string> &files, bool has_edge_data);

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

class ts_edge_data
{
	time_t timestamp;
	float weight;
public:
	ts_edge_data() {
		timestamp = 0;
		weight = 1;
	}

	ts_edge_data(time_t timestamp, float weight) {
		this->timestamp = timestamp;
		this->weight = weight;
	}

	time_t get_timestamp() const {
		return timestamp;
	}

	float get_weight() const {
		return weight;
	}
};

size_t parse_edge_list_line(char *line, edge<ts_edge_data> &e)
{
	int len = strlen(line);
	/*
	 * The format of a line should be
	 * from_vertex to_vertex "time" weight
	 * Fields are separated by tabs.
	 */
	if (line[len - 1] == '\n')
		line[len - 1] = 0;
	if (line[len - 2] == '\r')
		line[len - 2] = 0;
	if (line[0] == '#')
		return -1;
	char *second = strstr(line, delimiter);
	assert(second);
	*second = 0;
	second += strlen(delimiter);

	char *third = strstr(second, delimiter);
	assert(third);
	*third = 0;
	third += strlen(delimiter);
	if (*third == '"')
		third++;

	char *forth = strstr(third, delimiter);
	assert(forth);
	*forth = 0;
	if (*(forth - 1) == '"')
		*(forth - 1) = 0;

	if (!isnumeric(line) || !isnumeric(second)) {
		printf("%s\t%s\t%s\n", line, second, third);
		return -1;
	}
	long lfrom = atol(line);
	long lto = atol(second);
	assert(lfrom >= 0 && lfrom < MAX_VERTEX_ID);
	assert(lto >= 0 && lto < MAX_VERTEX_ID);
	vertex_id_t from = lfrom;
	vertex_id_t to = lto;
	struct tm tm;
	strptime(third, "%Y-%m-%d %H:%M:%S", &tm);
	time_t timestamp = mktime(&tm);
	// Let's ignore the weight on the edge first.
	ts_edge_data data(timestamp, 1);
	e = edge<ts_edge_data>(from, to, data);
	return 1;
}

int parse_edge_list_line(char *line, edge<edge_count> &e)
{
	if (line[0] == '#')
		return 0;
	char *second = strstr(line, delimiter);
	if (second == NULL) {
		fprintf(stderr, "wrong format 1: %s\n", line);
		return -1;
	}
	*second = 0;
	second += strlen(delimiter);
	char *third = strstr(second, delimiter);
	if (third == NULL) {
		fprintf(stderr, "wrong format 2: %s\n", second);
		return -1;
	}
	*third = 0;
	third += strlen(delimiter);
	if (!isnumeric(line) || !isnumeric(second) || !isnumeric(third)) {
		fprintf(stderr, "wrong format 3: %s\t%s\t%s\n", line, second, third);
		return -1;
	}
	long lfrom = atol(line);
	long lto = atol(second);
	assert(lfrom >= 0 && lfrom < MAX_VERTEX_ID);
	assert(lto >= 0 && lto < MAX_VERTEX_ID);
	vertex_id_t from = lfrom;
	vertex_id_t to = lto;
	edge_count c(atol(third));
	e = edge<edge_count>(from, to, c);

	return 1;
}

int parse_edge_list_line(char *line, edge<> &e)
{
	if (line[0] == '#')
		return 0;
	char *second = strstr(line, delimiter);
	if (second == NULL) {
		fprintf(stderr, "wrong format 1: %s\n", line);
		return -1;
	}
	*second = 0;
	second += strlen(delimiter);
	if (!isnumeric(line) || !isnumeric(second)) {
		fprintf(stderr, "wrong format 2: %s\t%s\n", line, second);
		return -1;
	}
	long lfrom = atol(line);
	long lto = atol(second);
	assert(lfrom >= 0 && lfrom < MAX_VERTEX_ID);
	assert(lto >= 0 && lto < MAX_VERTEX_ID);
	vertex_id_t from = lfrom;
	vertex_id_t to = lto;
	e = edge<>(from, to);

	return 1;
}

/**
 * Parse the edge list in the character buffer.
 * `size' doesn't include '\0'.
 */
template<class edge_data_type>
size_t parse_edge_list_text(char *line_buf, size_t size,
		std::vector<edge<edge_data_type> > &edges)
{
	char *line_end;
	char *line = line_buf;
	size_t num_edges = 0;
	while ((line_end = strchr(line, '\n'))) {
		*line_end = 0;
		if (*(line_end - 1) == '\r')
			*(line_end - 1) = 0;
		edge<edge_data_type> e;
		int num = parse_edge_list_line(line, e);
		if (num > 0)
			edges.push_back(e);
		num_edges += num;
		line = line_end + 1;
		assert(line - line_end <= (ssize_t) size);
	}
	if (line - line_buf < (ssize_t) size) {
		edge<edge_data_type> e;
		int num = parse_edge_list_line(line, e);
		if (num > 0)
			edges.push_back(e);
		num_edges += num;
	}
	return num_edges;
}

template<class edge_data_type>
class text_edge_task: public thread_task
{
	char *line_buf;
	size_t size;
public:
	text_edge_task(char *line_buf, size_t size) {
		this->line_buf = line_buf;
		this->size = size;
	}

	void run() {
		std::vector<edge<edge_data_type> > edges;
		parse_edge_list_text(line_buf, size, edges);
		std::vector<edge<edge_data_type> > *local_edge_buf
			= (std::vector<edge<edge_data_type> > *) thread::get_curr_thread()->get_user_data();
		local_edge_buf->insert(local_edge_buf->end(), edges.begin(), edges.end());
		delete [] line_buf;
	}
};

template<class edge_data_type>
directed_edge_graph<edge_count> *
directed_edge_graph<edge_data_type>::compress_edges() const
{
	directed_edge_graph<edge_count> *new_graph
		= new directed_edge_graph<edge_count>(true);
	printf("before: %ld edges\n", get_num_edges());
	if (!in_edges.empty()) {
		vertex_id_t from = in_edges[0].get_from();
		vertex_id_t to = in_edges[0].get_to();
		int num_duplicates = 1;
		for (size_t i = 1; i < in_edges.size(); i++) {
			if (in_edges[i].get_from() == from && in_edges[i].get_to() == to) {
				num_duplicates++;
			}
			else {
				edge_count c(num_duplicates);
				edge<edge_count> e(from, to, num_duplicates);
				new_graph->add_edge(e);

				num_duplicates = 1;
				from = in_edges[i].get_from();
				to = in_edges[i].get_to();
			}
		}
		edge_count c(num_duplicates);
		edge<edge_count> e(from, to, num_duplicates);
		new_graph->add_edge(e);
	}
	new_graph->sort_edges();

	printf("after: %ld edges\n", new_graph->get_num_edges());
	return new_graph;
}

template<class edge_data_type>
directed_edge_graph<edge_data_type> *
directed_edge_graph<edge_data_type>::simplify_edges() const
{
	directed_edge_graph<edge_data_type> *new_graph
		= new directed_edge_graph<edge_data_type>(false);
	printf("before: %ld in-edges and %ld out-edges\n", in_edges.size(),
			out_edges.size());
	if (!in_edges.empty()) {
		new_graph->in_edges.push_back(in_edges[0]);
		for (size_t i = 1; i < in_edges.size(); i++) {
			if (in_edges[i].get_from() != new_graph->in_edges.back().get_from()
					|| in_edges[i].get_to() != new_graph->in_edges.back().get_to())
				new_graph->in_edges.push_back(in_edges[i]);
		}
	}

	if (!out_edges.empty()) {
		new_graph->out_edges.push_back(out_edges[0]);
		for (size_t i = 1; i < out_edges.size(); i++) {
			if (out_edges[i].get_from() != new_graph->out_edges.back().get_from()
					|| out_edges[i].get_to() != new_graph->out_edges.back().get_to())
				new_graph->out_edges.push_back(out_edges[i]);
		}
	}
	printf("after: %ld in-edges and %ld out-edges\n", new_graph->in_edges.size(),
			new_graph->out_edges.size());
	return new_graph;
}

template<class edge_data_type>
void directed_edge_graph<edge_data_type>::sort_edges()
{
	comp_edge<edge_data_type> edge_comparator;
	comp_in_edge<edge_data_type> in_edge_comparator;
#ifdef MEMCHECK
	std::sort(out_edges.begin(), out_edges.end(), edge_comparator);
	std::sort(in_edges.begin(), in_edges.end(), in_edge_comparator);
#else
	__gnu_parallel::sort(out_edges.begin(), out_edges.end(), edge_comparator);
	__gnu_parallel::sort(in_edges.begin(), in_edges.end(), in_edge_comparator);
#endif
}

template<class edge_data_type>
void directed_edge_graph<edge_data_type>::construct_graph(graph *g) const
{
	vertex_id_t curr = 0;
	in_mem_directed_vertex<edge_data_type> v(curr, has_edge_data());
	size_t out_idx = 0;
	size_t in_idx = 0;
	size_t num_edges = in_edges.size();
	assert(in_edges.size() == out_edges.size());
	while (out_idx < num_edges && in_idx < num_edges) {
		while (out_edges[out_idx].get_from() == curr
				&& out_idx < num_edges) {
			v.add_out_edge(out_edges[out_idx++]);
		}
		while (in_edges[in_idx].get_to() == curr
				&& in_idx < num_edges) {
			v.add_in_edge(in_edges[in_idx++]);
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
			v = in_mem_directed_vertex<edge_data_type>(prev, has_edge_data());
			prev++;
			g->add_vertex(v);
		}
		v = in_mem_directed_vertex<edge_data_type>(curr, has_edge_data());
	}

	// Add remaining out-edges.
	while (out_idx < num_edges) {
		while (out_edges[out_idx].get_from() == curr
				&& out_idx < num_edges) {
			v.add_out_edge(out_edges[out_idx++]);
		}
		g->add_vertex(v);
		vertex_id_t prev = curr + 1;
		if (out_idx < num_edges)
			curr = out_edges[out_idx].get_from();
		else
			break;
		// The vertices without edges won't show up in the edge list,
		// but we need to fill the gap in the vertex Id space with empty
		// vertices.
		while (prev < curr) {
			v = in_mem_directed_vertex<edge_data_type>(prev, has_edge_data());
			prev++;
			g->add_vertex(v);
		}
		v = in_mem_directed_vertex<edge_data_type>(curr, has_edge_data());
	}

	// Add remaining in-edges
	while (in_idx < num_edges) {
		while (in_edges[in_idx].get_to() == curr
				&& in_idx < num_edges) {
			v.add_in_edge(in_edges[in_idx++]);
		}
		g->add_vertex(v);
		vertex_id_t prev = curr + 1;
		if (in_idx < num_edges)
			curr = in_edges[in_idx].get_to();
		else
			break;
		// The vertices without edges won't show up in the edge list,
		// but we need to fill the gap in the vertex Id space with empty
		// vertices.
		while (prev < curr) {
			v = in_mem_directed_vertex<edge_data_type>(prev, has_edge_data());
			prev++;
			g->add_vertex(v);
		}
		v = in_mem_directed_vertex<edge_data_type>(curr, has_edge_data());
	}
}

/**
 * This function loads edge lists from a tex file, parses them in parallel,
 * and convert the graph into the form of adjacency lists.
 */
template<class edge_data_type = empty_data>
directed_edge_graph<edge_data_type> *par_load_edge_list_text(
		const std::vector<std::string> &files, bool has_edge_data)
{
	struct timeval start, end;
	gettimeofday(&start, NULL);
	std::vector<task_thread *> threads(NUM_THREADS);
	for (int i = 0; i < NUM_THREADS; i++) {
		task_thread *t = new task_thread(std::string(
					"graph-task-thread") + itoa(i), i % NUM_NODES);
		t->set_user_data(new std::vector<edge<edge_data_type> >());
		t->start();
		threads[i] = t;
	}
	int thread_no = 0;
	printf("start to read the edge list\n");
	for (size_t i = 0; i < files.size(); i++) {
		graph_file_io io(files[i]);
		while (io.get_num_remaining_bytes() > 0) {
			size_t size = 0;
			char *line_buf = io.read_edge_list_text(EDGE_LIST_BLOCK_SIZE, size);
			assert(line_buf);
			thread_task *task = new text_edge_task<edge_data_type>(line_buf,
					size);
			threads[thread_no % NUM_THREADS]->add_task(task);
			thread_no++;
		}
	}
	for (int i = 0; i < NUM_THREADS; i++)
		threads[i]->wait4complete();
	gettimeofday(&end, NULL);
	printf("It takes %f seconds to construct edge list\n",
			time_diff(start, end));

	size_t mem_size = 0;
	size_t num_edges = 0;
	for (int i = 0; i < NUM_THREADS; i++) {
		std::vector<edge<edge_data_type> > *local_edges
			= (std::vector<edge<edge_data_type> > *) threads[i]->get_user_data();
		num_edges += local_edges->size();
		mem_size += local_edges->capacity() * sizeof(edge<edge_data_type>);
	}
	printf("There are %ld edges and use %ld bytes\n", num_edges, mem_size);

	directed_edge_graph<edge_data_type> *edge_g
		= new directed_edge_graph<edge_data_type>(num_edges, has_edge_data);
	start = end;
	for (int i = 0; i < NUM_THREADS; i++) {
		std::vector<edge<edge_data_type> > *local_edges
			= (std::vector<edge<edge_data_type> > *) threads[i]->get_user_data();
		edge_g->add_edges(*local_edges);
		delete local_edges;
	}
	gettimeofday(&end, NULL);
	printf("It takes %f seconds to combine edge list\n", time_diff(start, end));

	for (int i = 0; i < NUM_THREADS; i++) {
		threads[i]->stop();
		threads[i]->join();
		delete threads[i];
	}

	return edge_g;
}

static int get_file_id(const std::string &file)
{
	size_t begin = file.rfind('-');
	assert(begin != std::string::npos);
	size_t end = file.rfind('.');
	std::string sub = file.substr(begin + 1, end - begin - 1);
	return atoi(sub.c_str());
}

/**
 * I assume all files are named with numbers.
 * I need to sort the files in a numeric order.
 */
static void sort_edge_list_files(std::vector<std::string> &files)
{
	std::multimap<int, std::string> sorted_files;
	for (size_t i = 0; i < files.size(); i++) {
		int id = get_file_id(files[i]);
		sorted_files.insert(std::pair<int, std::string>(id, files[i]));
	}
	files.clear();
	for (std::map<int, std::string>::const_iterator it = sorted_files.begin();
			it != sorted_files.end(); it++)
		files.push_back(it->second);
}

template<class edge_data_type = empty_data>
graph *construct_directed_graph_compressed(
		const std::vector<std::string> &edge_list_files)
{
	struct timeval start, end;
	directed_edge_graph<edge_data_type> *edge_g
		= par_load_edge_list_text<edge_data_type>(edge_list_files,
				true);

	gettimeofday(&start, NULL);
	edge_g->sort_edges();
	gettimeofday(&end, NULL);
	printf("It takes %f seconds to sort edge list\n", time_diff(start, end));

	start = end;
	size_t orig_num_edges = edge_g->get_num_edges();
	directed_edge_graph<edge_count> *new_edge_g = edge_g->compress_edges();
	delete edge_g;
	gettimeofday(&end, NULL);
	printf("It takes %f seconds to compress edge list from %ld to %ld\n",
			time_diff(start, end), orig_num_edges,
			new_edge_g->get_num_edges());

	return new_edge_g->create_disk_graph();
}

template<class edge_data_type = empty_data>
graph *construct_directed_graph(
		const std::vector<std::string> &edge_list_files, bool has_edge_data)
{
	struct timeval start, end;
	directed_edge_graph<edge_data_type> *edge_g
		= par_load_edge_list_text<edge_data_type>(edge_list_files,
				has_edge_data);

	gettimeofday(&start, NULL);
	edge_g->sort_edges();
	gettimeofday(&end, NULL);
	printf("It takes %f seconds to sort edge list\n", time_diff(start, end));

	if (simplfy) {
		start = end;
		size_t orig_num_edges = edge_g->get_num_edges();
		directed_edge_graph<edge_data_type> *new_edge_g
			= edge_g->simplify_edges();
		delete edge_g;
		edge_g = new_edge_g;
		gettimeofday(&end, NULL);
		printf("It takes %f seconds to remove duplicated edges from %ld to %ld\n",
				time_diff(start, end), orig_num_edges, edge_g->get_num_edges());
	}

	return edge_g->create_disk_graph();
}

/**
 * The type of edge data.
 */
enum {
	DEFAULT_TYPE,
	EDGE_COUNT,
	EDGE_TIMESTAMP,
};

struct str2int_pair {
	std::string str;
	int number;
};
struct str2int_pair edge_type_map[] = {
	{"count", EDGE_COUNT},
	{"timestamp", EDGE_TIMESTAMP},
};
int type_map_size = sizeof(edge_type_map) / sizeof(edge_type_map[0]);

void print_usage()
{
	fprintf(stderr, "convert an edge list to adjacency lists\n");
	fprintf(stderr,
			"el2al [options] adj_list_file index_file edge_list_files (or directories)\n");
	fprintf(stderr, "-u: undirected graph\n");
	fprintf(stderr, "-d delimiter: the delimiter to seperate the input edge list\n");
	fprintf(stderr, "-c: compress the graph (remove duplicated edges)\n");
	fprintf(stderr, "-p: print adjacency list\n");
	fprintf(stderr, "-v: verify the created adjacency list\n");
	fprintf(stderr, "-t type: the type of edge data. Supported type: ");
	for (int i = 0; i < type_map_size; i++) {
		fprintf(stderr, "%s, ", edge_type_map[i].str.c_str());
	}
	fprintf(stderr, "\n");
	fprintf(stderr, "-m: merge multiple edge lists into a single graph. \n");
}

graph *construct_graph(const std::vector<std::string> &edge_list_files,
		int input_type)
{
	graph *g = NULL;
	switch(input_type) {
		case EDGE_COUNT:
			if (compress)
				g = construct_directed_graph_compressed<edge_count>(
						edge_list_files);
			else
				g = construct_directed_graph<edge_count>(
						edge_list_files, true);
			break;
		case EDGE_TIMESTAMP:
			if (compress)
				g = construct_directed_graph_compressed<ts_edge_data>(
						edge_list_files);
			else
				g = construct_directed_graph<ts_edge_data>(
						edge_list_files, true);
			break;
		default:
			if (compress)
				g = construct_directed_graph_compressed<>(edge_list_files);
			else
				g = construct_directed_graph<>(edge_list_files, false);
	}
	return g;
}

int main(int argc, char *argv[])
{
	int opt;
	bool directed = true;
	int num_opts = 0;
	char *type_str = NULL;
	bool merge_graph = false;
	while ((opt = getopt(argc, argv, "ud:cpvt:m")) != -1) {
		num_opts++;
		switch (opt) {
			case 'u':
				directed = false;
				break;
			case 'd':
				delimiter = optarg;
				num_opts++;
				break;
			case 'c':
				compress = true;
				break;
			case 'p':
				print_graph = true;
				break;
			case 'v':
				check_graph = true;
				break;
			case 't':
				type_str = optarg;
				num_opts++;
				break;
			case 'm':
				merge_graph = true;
				break;
			default:
				print_usage();
		}
	}
	argv += 1 + num_opts;
	argc -= 1 + num_opts;
	if (argc < 3) {
		print_usage();
		exit(-1);
	}

	int input_type = DEFAULT_TYPE;
	if (type_str) {
		for (int i = 0; i < type_map_size; i++) {
			if (edge_type_map[i].str == type_str) {
				input_type = edge_type_map[i].number;
				break;
			}
		}
	}

	std::string adjacency_list_file = argv[0];
	adjacency_list_file += std::string("-v") + itoa(CURR_VERSION);
	std::string index_file = argv[1];
	index_file += std::string("-v") + itoa(CURR_VERSION);
	std::vector<std::string> edge_list_files;
	for (int i = 2; i < argc; i++) {
		native_dir dir(argv[i]);
		if (dir.is_dir()) {
			std::vector<std::string> files;
			std::string dir_name = argv[i];
			dir.read_all_files(files);
			for (size_t i = 0; i < files.size(); i++)
				edge_list_files.push_back(dir_name + "/" + files[i]);
		}
		else
			edge_list_files.push_back(argv[i]);
	}
	sort_edge_list_files(edge_list_files);

	for (size_t i = 0; i < edge_list_files.size(); i++)
		printf("edge list file: %s\n", edge_list_files[i].c_str());

	assert(directed);
	if (merge_graph) {
		struct timeval start, end;
		gettimeofday(&start, NULL);
		graph *g = construct_graph(edge_list_files, input_type);

		// Write the constructed individual graph to a file.
		g->dump(index_file, adjacency_list_file);

		gettimeofday(&end, NULL);
		printf("It takes %f seconds to dump the graph to a file\n",
				time_diff(start, end));
		printf("There are %ld vertices, %ld non-empty vertices and %ld edges\n",
				g->get_num_vertices(), g->get_num_non_empty_vertices(),
				g->get_num_edges());
		if (print_graph)
			g->print();
		if (check_graph)
			g->check_ext_graph(index_file, adjacency_list_file);
		delete g;
	}
	else {
		std::vector<std::string> graph_files;
		std::vector<std::string> index_files;
		if (edge_list_files.size() > 1) {
			for (size_t i = 0; i < edge_list_files.size(); i++) {
				graph_files.push_back(adjacency_list_file + "-" + itoa(i));
				index_files.push_back(index_file + "-" + itoa(i));
			}
		}
		else {
			graph_files.push_back(adjacency_list_file);
			index_files.push_back(index_file);
		}

		struct timeval start, end;
		for (size_t i = 0; i < edge_list_files.size(); i++) {
			// construct individual graphs.
			gettimeofday(&start, NULL);
			std::vector<std::string> files(1);
			files[0] = edge_list_files[i];
			graph *g = construct_graph(files, input_type);

			// Write the constructed individual graph to a file.
			g->dump(index_files[i], graph_files[i]);

			gettimeofday(&end, NULL);
			printf("It takes %f seconds to dump the graph to a file\n",
					time_diff(start, end));
			printf("There are %ld vertices, %ld non-empty vertices and %ld edges\n",
					g->get_num_vertices(), g->get_num_non_empty_vertices(),
					g->get_num_edges());
			if (print_graph)
				g->print();
			if (check_graph)
				g->check_ext_graph(index_files[i], graph_files[i]);
			delete g;
		}
	}
}
