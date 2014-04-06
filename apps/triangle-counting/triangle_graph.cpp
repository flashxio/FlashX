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

#include <signal.h>
#include <google/profiler.h>

#include "thread.h"
#include "io_interface.h"

#include "graph_engine.h"
#include "graph_config.h"
#include "graphlab/cuckoo_set_pow2.hpp"

const double BIN_SEARCH_RATIO = 100;
const int HASH_SEARCH_RATIO = 16;

atomic_number<long> num_working_vertices;
atomic_number<long> num_completed_vertices;

int hash_threshold = 1000;

#if 0
class vertex_size_scheduler: public vertex_scheduler
{
public:
	void schedule(std::vector<compute_vertex *> &vertices);
};

void vertex_size_scheduler::schedule(std::vector<compute_vertex *> &vertices)
{
	class comp_size
	{
	public:
		bool operator()(const compute_vertex *v1, const compute_vertex *v2) {
			return v1->get_ext_mem_size() > v2->get_ext_mem_size();
		}
	};

	std::sort(vertices.begin(), vertices.end(), comp_size());
}
#endif

class count_msg: public vertex_message
{
	int num;
public:
	count_msg(int num): vertex_message(sizeof(count_msg), false) {
		this->num = num;
	}

	int get_num() const {
		return num;
	}
};

struct runtime_data_t
{
	class index_entry
	{
		vertex_id_t id;
		uint32_t idx;
	public:
		index_entry() {
			id = -1;
			idx = -1;
		}

		index_entry(vertex_id_t id) {
			this->id = id;
			this->idx = -1;
		}

		index_entry(vertex_id_t id, uint32_t idx) {
			this->id = id;
			this->idx = idx;
		}

		vertex_id_t get_id() const {
			return id;
		}

		uint32_t get_idx() const {
			return idx;
		}

		bool operator==(const index_entry &e) const {
			return id == e.get_id();
		}
	};

	class index_hash
	{
		boost::hash<vertex_id_t> id_hash;
	public:
		size_t operator()(const index_entry &e) const {
			return id_hash(e.get_id());
		}
	};
	typedef graphlab::cuckoo_set_pow2<index_entry, 3, size_t,
			index_hash> edge_set_t;
	// It contains part of in-edges.
	// We only use the neighbors whose ID is smaller than this vertex.
	std::vector<vertex_id_t> in_edges;
	// The vector contains the number of the neighbors' triangles shared
	// with this vertex. It only keeps the triangles of neighbors in the
	// in-edges.
	std::vector<int> triangles;
	// The number of vertices that have joined with the vertex.
	size_t num_joined;
	size_t num_required;
	size_t num_triangles;

	edge_set_t in_edge_set;
public:
	runtime_data_t(const std::vector<vertex_id_t> &in_edges,
			const std::vector<vertex_id_t> &out_edges, size_t num_triangles): in_edge_set(
				index_entry(), 0, 2 * in_edges.size()) {
		num_joined = 0;
		num_required = out_edges.size();
		this->in_edges = in_edges;
		triangles.resize(in_edges.size());
		this->num_triangles = num_triangles;

		// We only build a hash table on large vertices
		if (in_edges.size() > (size_t) hash_threshold)
			for (size_t i = 0; i < in_edges.size(); i++)
				in_edge_set.insert(index_entry(in_edges[i], i));
	}
};

enum multi_func_flags
{
	NUM_TRIANGLES,
	POINTER,
	NUM_FLAGS,
};

class multi_func_value
{
	static const int VALUE_BITS = sizeof(size_t) * 8 - NUM_FLAGS;
	static const size_t FLAGS_MASK = ((1UL << VALUE_BITS) - 1);
	size_t value;

	void set_flag(int flag) {
		value |= 1UL << (VALUE_BITS + flag);
	}

	bool has_flag(int flag) const {
		return value & (1UL << (VALUE_BITS + flag));
	}
public:
	multi_func_value() {
		value = 0;
		// By default, it stores the number of triangles.
		set_flag(NUM_TRIANGLES);
	}

	void set_num_triangles(size_t num) {
		value = num;
		set_flag(NUM_TRIANGLES);
	}

	bool has_num_triangles() const {
		return has_flag(NUM_TRIANGLES);
	}

	size_t get_num_triangles() const {
		assert(has_flag(NUM_TRIANGLES));
		return value & FLAGS_MASK;
	}

	void inc_num_triangles(size_t num) {
		assert(has_flag(NUM_TRIANGLES));
		value += num;
	}

	/**
	 * Pointer to the runtime data.
	 */

	void set_runtime_data(runtime_data_t *data) {
		value = (size_t) data;
		set_flag(POINTER);
	}

	bool has_runtime_data() const {
		return has_flag(POINTER);
	}

	runtime_data_t *get_runtime_data() const {
		assert(has_flag(POINTER));
		return (runtime_data_t *) (value & FLAGS_MASK);
	}
};

class triangle_vertex: public compute_directed_vertex
{
	multi_func_value local_value;

	void inc_num_triangles(size_t num) {
		if (local_value.has_num_triangles())
			local_value.inc_num_triangles(num);
		else
			local_value.get_runtime_data()->num_triangles += num;
	}
public:
	triangle_vertex() {
	}

	triangle_vertex(vertex_id_t id,
			const vertex_index *index): compute_directed_vertex(id, index) {
	}

	int count_triangles(const page_vertex *v) const;

	int get_num_triangles() const {
		return local_value.get_num_triangles();
	}

	void run(graph_engine &graph) {
		vertex_id_t id = get_id();
		request_vertices(&id, 1);
	}

	void run(graph_engine &graph, const page_vertex &vertex) {
		if (vertex.get_id() == get_id())
			run_on_itself(graph, vertex);
		else
			run_on_neighbor(graph, vertex);
	}

	void run_on_itself(graph_engine &graph, const page_vertex &vertex);
	void run_on_neighbor(graph_engine &graph, const page_vertex &vertex);

	void run_on_message(graph_engine &graph, const vertex_message &msg) {
		inc_num_triangles(((count_msg &) msg).get_num());
	}
};

int triangle_vertex::count_triangles(const page_vertex *v) const
{
	int num_local_triangles = 0;
	assert(v->get_id() != this->get_id());

	if (v->get_num_edges(edge_type::OUT_EDGE) == 0)
		return 0;

	/*
	 * We search for triangles with two different ways:
	 * binary search if two adjacency lists have very different sizes,
	 * scan otherwise.
	 *
	 * when binary search for multiple neighbors, we can reduce binary search
	 * overhead by using the new end in the search range. We can further reduce
	 * overhead by searching in a reverse order (start from the largest neighbor).
	 * Since vertices of smaller ID has more neighbors, it's more likely
	 * that a neighbor is in the beginning of the adjacency list, and
	 * the search range will be narrowed faster.
	 */

	runtime_data_t *data = local_value.get_runtime_data();
	if (data->in_edge_set.size() > 0
			&& data->in_edges.size() > HASH_SEARCH_RATIO * v->get_num_edges(
				edge_type::OUT_EDGE)) {
		page_byte_array::const_iterator<vertex_id_t> other_it
			= v->get_neigh_begin(edge_type::OUT_EDGE);
		page_byte_array::const_iterator<vertex_id_t> other_end
			= v->get_neigh_end(edge_type::OUT_EDGE);
		for (; other_it != other_end; ++other_it) {
			vertex_id_t neigh_neighbor = *other_it;
			runtime_data_t::edge_set_t::const_iterator it
				= data->in_edge_set.find(neigh_neighbor);
			if (it != data->in_edge_set.end()) {
				if (neigh_neighbor != v->get_id()
						&& neigh_neighbor != this->get_id()) {
					num_local_triangles++;
					int idx = (*it).get_idx();
					data->triangles[idx]++;
				}
			}
		}
	}
	// If the neighbor vertex has way more edges than this vertex.
	else if (v->get_num_edges(edge_type::OUT_EDGE) / data->in_edges.size(
				) > BIN_SEARCH_RATIO) {
		page_byte_array::const_iterator<vertex_id_t> other_it
			= v->get_neigh_begin(edge_type::OUT_EDGE);
		page_byte_array::const_iterator<vertex_id_t> other_end
			= v->get_neigh_end(edge_type::OUT_EDGE);
		for (int i = data->in_edges.size() - 1; i >= 0; i--) {
			vertex_id_t this_neighbor = data->in_edges.at(i);
			// We need to skip loops.
			if (this_neighbor != v->get_id()
					&& this_neighbor != this->get_id()) {
				page_byte_array::const_iterator<vertex_id_t> first
					= std::lower_bound(other_it, other_end, this_neighbor);
				if (first != other_end && this_neighbor == *first) {
					num_local_triangles++;
					data->triangles[i]++;
				}
				other_end = first;
			}
		}
	}
	else {
		std::vector<vertex_id_t>::const_iterator this_it = data->in_edges.begin();
		std::vector<int>::iterator count_it = data->triangles.begin();
		std::vector<vertex_id_t>::const_iterator this_end = data->in_edges.end();
		page_byte_array::seq_const_iterator<vertex_id_t> other_it
			= v->get_neigh_seq_it(edge_type::OUT_EDGE, 0,
					v->get_num_edges(edge_type::OUT_EDGE));
		while (this_it != this_end && other_it.has_next()) {
			vertex_id_t this_neighbor = *this_it;
			vertex_id_t neigh_neighbor = other_it.curr();
			if (this_neighbor == neigh_neighbor) {
				// skip loop
				if (neigh_neighbor != v->get_id()
						&& neigh_neighbor != this->get_id()) {
					num_local_triangles++;
					(*count_it)++;
				}
				++this_it;
				other_it.next();
				++count_it;
			}
			else if (this_neighbor < neigh_neighbor) {
				++this_it;
				++count_it;
			}
			else
				other_it.next();
		}
	}
	return num_local_triangles;
}

void triangle_vertex::run_on_itself(graph_engine &graph, const page_vertex &vertex)
{
	assert(!local_value.has_runtime_data());

	long ret = num_working_vertices.inc(1);
	if (ret % 100000 == 0)
		printf("%ld working vertices\n", ret);
	// A vertex has to have in-edges and out-edges in order to form
	// a triangle. so we can simply skip the vertices that don't have
	// either of them.
	if (vertex.get_num_edges(edge_type::OUT_EDGE) == 0
			|| vertex.get_num_edges(edge_type::IN_EDGE) == 0) {
		long ret = num_completed_vertices.inc(1);
		if (ret % 100000 == 0)
			printf("%ld completed vertices\n", ret);
		return;
	}

	std::vector<vertex_id_t> in_edges;
	std::vector<vertex_id_t> out_edges;

	page_byte_array::const_iterator<vertex_id_t> it
		= vertex.get_neigh_begin(edge_type::IN_EDGE);
	page_byte_array::const_iterator<vertex_id_t> end
		= vertex.get_neigh_end(edge_type::IN_EDGE);
	int num_local_edges = this->get_num_in_edges() + this->get_num_out_edges();
	for (; it != end; ++it) {
		vertex_id_t id = *it;
		triangle_vertex &v1 = (triangle_vertex &) graph.get_vertex(id);
		int num_local_edges1 = v1.get_num_in_edges() + v1.get_num_out_edges();
		if ((num_local_edges1 < num_local_edges && id != vertex.get_id())
				|| (num_local_edges1 == num_local_edges
					&& id < vertex.get_id())) {
			in_edges.push_back(id);
		}
	}

	it = vertex.get_neigh_begin(edge_type::OUT_EDGE);
	end = vertex.get_neigh_end(edge_type::OUT_EDGE);
	for (; it != end; ++it) {
		vertex_id_t id = *it;
		triangle_vertex &v1 = (triangle_vertex &) graph.get_vertex(id);
		int num_local_edges1 = v1.get_num_in_edges() + v1.get_num_out_edges();
		if ((num_local_edges1 < num_local_edges && id != vertex.get_id())
				|| (num_local_edges1 == num_local_edges
					&& id < vertex.get_id())) {
			out_edges.push_back(id);
		}
	}

	if (in_edges.empty() || out_edges.empty()) {
		long ret = num_completed_vertices.inc(1);
		if (ret % 100000 == 0)
			printf("%ld completed vertices\n", ret);
		return;
	}

	std::vector<directed_vertex_request> reqs(out_edges.size());
	for (size_t i = 0; i < out_edges.size(); i++) {
		vertex_id_t id = out_edges[i];
		reqs[i] = directed_vertex_request(id, edge_type::OUT_EDGE);
	}
	// We have to set runtime data before calling request_partial_vertices.
	// It's possible that the request to a partial vertex can be completed
	// immediately and run_on_neighbor is called in request_partial_vertices.
	// TODO Maybe I should avoid that.
	local_value.set_runtime_data(new runtime_data_t(in_edges, out_edges,
				local_value.get_num_triangles()));
	request_partial_vertices(reqs.data(), reqs.size());
}

void triangle_vertex::run_on_neighbor(graph_engine &graph,
		const page_vertex &vertex)
{
	assert(local_value.has_runtime_data());
	runtime_data_t *data = local_value.get_runtime_data();
	data->num_joined++;
	int ret = count_triangles(&vertex);
	// If we find triangles with the neighbor, notify the neighbor
	// as well.
	if (ret > 0) {
		inc_num_triangles(ret);
		count_msg msg(ret);
		graph.send_msg(vertex.get_id(), msg);
	}

	// If we have seen all required neighbors, we have complete
	// the computation. We can release the memory now.
	if (data->num_joined == data->num_required) {
		long ret = num_completed_vertices.inc(1);
		if (ret % 100000 == 0)
			printf("%ld completed vertices\n", ret);

		// Inform all neighbors in the in-edges.
		for (size_t i = 0; i < data->triangles.size(); i++) {
			// Inform the neighbor if they share triangles.
			if (data->triangles[i] > 0) {
				count_msg msg(data->triangles[i]);
				graph.send_msg(data->in_edges[i], msg);
			}
		}
		size_t num_curr_triangles = data->num_triangles;
		delete data;
		local_value.set_num_triangles(num_curr_triangles);
	}
}

void int_handler(int sig_num)
{
	if (!graph_conf.get_prof_file().empty())
		ProfilerStop();
	exit(0);
}

void print_usage()
{
	fprintf(stderr,
			"triangle-counting [options] conf_file graph_file index_file\n");
	fprintf(stderr, "-o file: output local scan of each vertex to a file\n");
	fprintf(stderr, "-c confs: add more configurations to the system\n");
	graph_conf.print_help();
	params.print_help();
}

int main(int argc, char *argv[])
{
	int opt;
	std::string output_file;
	std::string confs;
	int num_opts = 0;
	while ((opt = getopt(argc, argv, "o:c:")) != -1) {
		num_opts++;
		switch (opt) {
			case 'o':
				output_file = optarg;
				num_opts++;
				break;
			case 'c':
				confs = optarg;
				num_opts++;
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

	std::string conf_file = argv[0];
	std::string graph_file = argv[1];
	std::string index_file = argv[2];

	config_map configs(conf_file);
	configs.add_options(confs);
	graph_conf.init(configs);
	graph_conf.print();

	signal(SIGINT, int_handler);
	init_io_system(configs);

	graph_index *index = NUMA_graph_index<triangle_vertex>::create(
			index_file, graph_conf.get_num_threads(), params.get_num_nodes());
	graph_engine *graph = graph_engine::create(
			graph_conf.get_num_threads(), params.get_num_nodes(),
			graph_file, index);
#if 0
	// Let's schedule the order of processing activated vertices according
	// to the size of vertices. We start with processing vertices with higher
	// degrees in the hope we can find the max scan as early as possible,
	// so that we can simple ignore the rest of vertices.
	graph->set_vertex_scheduler(new vertex_size_scheduler());
#endif

	printf("triangle counting starts\n");
	printf("prof_file: %s\n", graph_conf.get_prof_file().c_str());
	if (!graph_conf.get_prof_file().empty())
		ProfilerStart(graph_conf.get_prof_file().c_str());

	struct timeval start, end;
	gettimeofday(&start, NULL);
	graph->start_all();
	gettimeofday(&end, NULL);
	float activate_time = time_diff(start, end);

	start = end;
	graph->wait4complete();
	gettimeofday(&end, NULL);

	// Count the total number of triangles in the graph.
	graph_index::const_iterator it = index->begin();
	graph_index::const_iterator end_it = index->end();
	long num_triangles = 0;
	for (; it != end_it; ++it) {
		const triangle_vertex &v = (const triangle_vertex &) *it;
		num_triangles += v.get_num_triangles();
	}

	if (!graph_conf.get_prof_file().empty())
		ProfilerStop();
	if (graph_conf.get_print_io_stat())
		print_io_thread_stat();
	graph_engine::destroy(graph);
	destroy_io_system();
	printf("There are %ld vertices\n", index->get_num_vertices());
	printf("process %ld vertices and complete %ld vertices\n",
			num_working_vertices.get(), num_completed_vertices.get());
	printf("there are %ld triangles.\n", num_triangles);
	printf("It takes %f seconds to activate all and %f seconds to finish\n",
			activate_time, time_diff(start, end));
}
