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

#include <vector>
#include <unordered_map>

#include "thread.h"
#include "io_interface.h"
#include "container.h"
#include "concurrency.h"

#include "vertex_index.h"
#include "graph_engine.h"
#include "graph_config.h"

class global_max
{
	volatile size_t value;
	pthread_spinlock_t lock;
public:
	global_max() {
		value = 0;
		pthread_spin_init(&lock, PTHREAD_PROCESS_PRIVATE);
	}

	global_max(size_t init) {
		value = init;
		pthread_spin_init(&lock, PTHREAD_PROCESS_PRIVATE);
	}

	bool update(size_t new_v) {
		if (new_v <= value)
			return false;

		bool ret = false;
		pthread_spin_lock(&lock);
		if (new_v > value) {
			value = new_v;
			ret = true;
		}
		pthread_spin_unlock(&lock);
		return ret;
	}

	size_t get() const {
		return value;
	}
} max_dist;

const int K = 10;
int num_bfs = 1;

class diameter_message: public vertex_message
{
	uint16_t dists[K];
public:
	diameter_message(uint16_t dists[], int num): vertex_message(
			sizeof(diameter_message), true) {
		memcpy(this->dists, dists, sizeof(dists[0]) * num);
	}

	uint16_t get_dist(int idx) const {
		return dists[idx];
	}
};

class diameter_vertex: public compute_directed_vertex
{
	uint16_t dists[K];
	bool updated;
public:
	diameter_vertex() {
		updated = false;
		for (int i = 0; i < K; i++)
			dists[i] = USHRT_MAX;
	}

	diameter_vertex(vertex_id_t id,
			const vertex_index *index): compute_directed_vertex(id, index) {
		updated = false;
		for (int i = 0; i < K; i++)
			dists[i] = USHRT_MAX;
	}

	void init(int bfs_id) {
		dists[bfs_id] = 0;
		updated = true;
		printf("BFS %d starts at v%u\n", bfs_id, get_id());
	}

	void reset() {
		updated = false;
		for (int i = 0; i < K; i++)
			dists[i] = USHRT_MAX;
	}

	uint16_t get_dist(int idx) const {
		return dists[idx];
	}

	void run(graph_engine &graph) {
		if (updated) {
			updated = false;
			directed_vertex_request req(get_id(), edge_type::OUT_EDGE);
			request_partial_vertices(&req, 1);
		}
	}

	void run(graph_engine &graph, const page_vertex &vertex);

	void run_on_message(graph_engine &, const vertex_message &msg) {
		const diameter_message &dmsg = (const diameter_message &) msg;
		for (int i = 0; i < num_bfs; i++) {
			if (dmsg.get_dist(i) + 1 < dists[i]) {
				dists[i] = dmsg.get_dist(i) + 1;
				updated = true;
				max_dist.update(dists[i]);
			}
		}
	}
};

void diameter_vertex::run(graph_engine &graph, const page_vertex &vertex)
{
	int num_dests = vertex.get_num_edges(OUT_EDGE);
	if (num_dests == 0)
		return;

	// We need to add the neighbors of the vertex to the queue of
	// the next level.
	edge_seq_iterator it = vertex.get_neigh_seq_it(OUT_EDGE, 0, num_dests);
	diameter_message msg(dists, num_bfs);
	graph.multicast_msg(it, msg);
}

class diameter_initiator: public vertex_initiator
{
	std::unordered_map<vertex_id_t, int> start_vertices;
public:
	diameter_initiator(std::vector<vertex_id_t> &vertices) {
		for (size_t i = 0; i < vertices.size(); i++) {
			start_vertices.insert(std::pair<vertex_id_t, int>(vertices[i], i));
		}
	}

	void init(compute_vertex &v) {
		diameter_vertex &dv = (diameter_vertex &) v;
		std::unordered_map<vertex_id_t, int>::const_iterator it
			= start_vertices.find(v.get_id());
		assert(it != start_vertices.end());
		dv.init(it->second);
	}
};

void int_handler(int sig_num)
{
	if (!graph_conf.get_prof_file().empty())
		ProfilerStop();
	exit(0);
}

void print_usage()
{
	fprintf(stderr,
			"diameter_estimator [options] conf_file graph_file index_file\n");
	fprintf(stderr, "-c confs: add more configurations to the system\n");
	fprintf(stderr, "-p: preload the graph\n");
	graph_conf.print_help();
	params.print_help();
}

int main(int argc, char *argv[])
{
	int opt;
	std::string confs;
	int num_opts = 0;
	bool preload = false;
	bool in_parallel = false;
	while ((opt = getopt(argc, argv, "c:pn:P")) != -1) {
		num_opts++;
		switch (opt) {
			case 'c':
				confs = optarg;
				num_opts++;
				break;
			case 'p':
				preload = true;
				break;
			case 'n':
				num_bfs = atoi(optarg);
				num_bfs = min(num_bfs, K);
				num_opts++;
				break;
			case 'P':
				in_parallel = true;
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

	graph_index *index = NUMA_graph_index<diameter_vertex>::create(index_file,
			graph_conf.get_num_threads(), params.get_num_nodes());

	graph_engine *graph = graph_engine::create(graph_conf.get_num_threads(),
			params.get_num_nodes(), graph_file, index);
	if (preload)
		graph->preload_graph();
	printf("BFS starts\n");
	printf("prof_file: %s\n", graph_conf.get_prof_file().c_str());
	if (!graph_conf.get_prof_file().empty())
		ProfilerStart(graph_conf.get_prof_file().c_str());

	std::vector<vertex_id_t> start_vertices;
	srandom(time(NULL));
	while (start_vertices.size() < (size_t) num_bfs) {
		vertex_id_t id = random() % graph->get_max_vertex_id();
		// We should skip the empty vertices.
		if (graph->get_vertex_edges(id) == 0)
			continue;

		start_vertices.push_back(id);
	}

	if (in_parallel) {
		struct timeval start, end;
		gettimeofday(&start, NULL);
		graph->start(start_vertices.data(), start_vertices.size(),
				vertex_initiator::ptr(new diameter_initiator(start_vertices)));
		graph->wait4complete();
		gettimeofday(&end, NULL);
		printf("It takes %f seconds\n", time_diff(start, end));

#if 0
		for (size_t i = 0; i < start_vertices.size(); i++) {
			graph_index::const_iterator it = index->begin();
			graph_index::const_iterator end_it = index->end();
			for (; it != end_it; ++it) {
				diameter_vertex &v = (diameter_vertex &) *it;
				if (v.get_dist(i) != USHRT_MAX)
					fprintf(stderr, "v%u: %d\n", v.get_id(), v.get_dist(i));
			}
		}
#endif
	}
	else {
		num_bfs = 1;
		for (size_t i = 0; i < start_vertices.size(); i++) {
			struct timeval start, end;
			gettimeofday(&start, NULL);
			std::vector<vertex_id_t> local_starts;
			local_starts.push_back(start_vertices[i]);
			graph->start(local_starts.data(), local_starts.size(),
					vertex_initiator::ptr(new diameter_initiator(local_starts)));
			graph->wait4complete();
			gettimeofday(&end, NULL);

			// Clean up all vertices.
			graph_index::const_iterator it = index->begin();
			graph_index::const_iterator end_it = index->end();
			for (; it != end_it; ++it) {
				diameter_vertex &v = (diameter_vertex &) *it;
#if 0
				if (v.get_dist(0) != USHRT_MAX)
					fprintf(stderr, "v%u: %d\n", v.get_id(), v.get_dist(0));
#endif
				v.reset();
			}
			printf("BFS %ld takes %f seconds\n", i, time_diff(start, end));
		}
	}

	if (!graph_conf.get_prof_file().empty())
		ProfilerStop();
	if (graph_conf.get_print_io_stat())
		print_io_thread_stat();
	graph_engine::destroy(graph);
	destroy_io_system();
	printf("The max dist: %ld\n", max_dist.get());
}
