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

#include <signal.h>
#include <google/profiler.h>

#include <atomic>
#include <vector>

#include "thread.h"
#include "io_interface.h"
#include "container.h"
#include "concurrency.h"

#include "vertex_index.h"
#include "graph_engine.h"
#include "graph_config.h"

/**
 * Measure the performance of message passing.
 * It can remove remote random memory access, but still generates
 * a lot of random access in the local memory.
 */

class test_message: public vertex_message
{
	int v;
public:
	test_message(int v): vertex_message(sizeof(test_message), false) {
		this->v = v;
	}

	int get_value() const {
		return v;
	}
};

class test_vertex: public compute_directed_vertex
{
	long sum;
	char stuffing[64 - sizeof(compute_directed_vertex) - sizeof(long)];
public:
	test_vertex() {
		sum = 0;
	}

	test_vertex(vertex_id_t id, const vertex_index &index): compute_directed_vertex(
			id, index) {
		sum = 0;
	}

	void run(vertex_program &prog) {
		vertex_id_t id = get_id();
		request_vertices(&id, 1);
	}

	void run(vertex_program &prog, const page_vertex &vertex);

	virtual void run_on_message(vertex_program &, const vertex_message &msg1) {
		const test_message &msg = (const test_message &) msg1;
		sum += msg.get_value();
	}
};

void test_vertex::run(vertex_program &prog, const page_vertex &vertex)
{
	int num_dests = vertex.get_num_edges(BOTH_EDGES);
	if (num_dests == 0)
		return;

	edge_seq_iterator it = vertex.get_neigh_seq_it(BOTH_EDGES, 0, num_dests);
	test_message msg(1);
	prog.multicast_msg(it, msg);
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
			"test [options] conf_file graph_file index_file\n");
	fprintf(stderr, "-c confs: add more configurations to the system\n");
	graph_conf.print_help();
	params.print_help();
}

int main(int argc, char *argv[])
{
	int opt;
	std::string confs;
	int num_opts = 0;
	while ((opt = getopt(argc, argv, "c:")) != -1) {
		num_opts++;
		switch (opt) {
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
	printf("The size of vertex state: %ld\n", sizeof(test_vertex));

	signal(SIGINT, int_handler);

	graph_index::ptr index = NUMA_graph_index<test_vertex>::create(index_file);
	graph_engine::ptr graph = graph_engine::create(graph_file, index, configs);
	graph->preload_graph();
	printf("test starts\n");
	printf("prof_file: %s\n", graph_conf.get_prof_file().c_str());
	if (!graph_conf.get_prof_file().empty())
		ProfilerStart(graph_conf.get_prof_file().c_str());

	struct timeval start, end;
	gettimeofday(&start, NULL);
	graph->start_all();
	graph->wait4complete();
	gettimeofday(&end, NULL);

	if (!graph_conf.get_prof_file().empty())
		ProfilerStop();
	if (graph_conf.get_print_io_stat())
		print_io_thread_stat();
	printf("It takes %f seconds\n", time_diff(start, end));
}
