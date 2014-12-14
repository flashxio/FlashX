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
#ifdef PROFILER
#include <gperftools/profiler.h>
#endif

#include <atomic>
#include <vector>

#include "thread.h"
#include "io_interface.h"
#include "container.h"
#include "concurrency.h"

#include "vertex_index.h"
#include "graph_engine.h"
#include "graph_config.h"
#include "FGlib.h"

using namespace safs;
using namespace fg;

/**
 * Measure the performance of gather data from neighbors.
 * It generates a lot of remote random memory access.
 */
class test_get_vertex: public compute_directed_vertex
{
	long sum;
	char stuffing[64 - sizeof(compute_directed_vertex) - sizeof(long)];
public:
	test_get_vertex(vertex_id_t id): compute_directed_vertex(id) {
		sum = 0;
	}

	void run(vertex_program &prog) {
		vertex_id_t id = prog.get_vertex_id(*this);
		request_vertices(&id, 1);
	}

	void run(vertex_program &prog, const page_vertex &vertex);

	virtual void run_on_message(vertex_program &, const vertex_message &msg) {
	}
};

void test_get_vertex::run(vertex_program &prog, const page_vertex &vertex)
{
	long local_sum = 0;
	int num_dests = vertex.get_num_edges(BOTH_EDGES);
	if (num_dests == 0)
		return;

	edge_seq_iterator it = vertex.get_neigh_seq_it(BOTH_EDGES, 0, num_dests);
	PAGE_FOREACH(vertex_id_t, id, it) {
		test_get_vertex &v = (test_get_vertex &) prog.get_graph().get_vertex(id);
		local_sum += v.sum;
	} PAGE_FOREACH_END
	this->sum = local_sum;
}

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

/**
 * Measure the performance of message passing.
 * It can remove remote random memory access, but still generates
 * a lot of random access in the local memory.
 */
class test_msg_vertex: public compute_directed_vertex
{
	long sum;
	char stuffing[64 - sizeof(compute_directed_vertex) - sizeof(long)];
public:
	test_msg_vertex(vertex_id_t id): compute_directed_vertex(id) {
		sum = 0;
	}

	void run(vertex_program &prog) {
		vertex_id_t id = prog.get_vertex_id(*this);
		request_vertices(&id, 1);
	}

	void run(vertex_program &prog, const page_vertex &vertex);

	virtual void run_on_message(vertex_program &, const vertex_message &msg1) {
		const test_message &msg = (const test_message &) msg1;
		sum += msg.get_value();
	}
};

void test_msg_vertex::run(vertex_program &prog, const page_vertex &vertex)
{
	int num_dests = vertex.get_num_edges(BOTH_EDGES);
	if (num_dests == 0)
		return;

	edge_seq_iterator it = vertex.get_neigh_seq_it(BOTH_EDGES, 0, num_dests);
	test_message msg(1);
	prog.multicast_msg(it, msg);
}

/**
 * Measure the performance of gather data from neighbors.
 * It generates a lot of remote random memory access.
 */
class test_push_vertex: public compute_directed_vertex
{
	std::atomic<long> sum;
public:
	test_push_vertex(vertex_id_t id): compute_directed_vertex(id) {
		sum = 0;
	}

	void run(vertex_program &prog) {
		vertex_id_t id = prog.get_vertex_id(*this);
		request_vertices(&id, 1);
	}

	void run(vertex_program &prog, const page_vertex &vertex);

	virtual void run_on_message(vertex_program &, const vertex_message &msg) {
	}
};

void test_push_vertex::run(vertex_program &prog, const page_vertex &vertex)
{
	int num_dests = vertex.get_num_edges(BOTH_EDGES);
	if (num_dests == 0)
		return;

	edge_seq_iterator it = vertex.get_neigh_seq_it(BOTH_EDGES, 0, num_dests);
	PAGE_FOREACH(vertex_id_t, id, it) {
		test_push_vertex &v = (test_push_vertex &) prog.get_graph().get_vertex(id);
		v.sum.fetch_add(1, std::memory_order_seq_cst);
	} PAGE_FOREACH_END
}

void int_handler(int sig_num)
{
#ifdef PROFILER
	if (!graph_conf.get_prof_file().empty())
		ProfilerStop();
#endif
	exit(0);
}

void print_usage()
{
	fprintf(stderr,
			"test [options] conf_file graph_file index_file test\n");
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

	if (argc < 4) {
		print_usage();
		exit(-1);
	}

	std::string conf_file = argv[0];
	std::string graph_file = argv[1];
	std::string index_file = argv[2];
	std::string test = argv[3];

	config_map::ptr configs = config_map::create(conf_file);
	configs->add_options(confs);
	if (test == "push")
		printf("The size of vertex state: %ld\n", sizeof(test_push_vertex));
	else if (test == "msg")
		printf("The size of vertex state: %ld\n", sizeof(test_msg_vertex));
	else if (test == "get")
		printf("The size of vertex state: %ld\n", sizeof(test_get_vertex));
	else {
		fprintf(stderr, "wrong test\n");
		exit(1);
	}

	signal(SIGINT, int_handler);

	FG_graph::ptr fg = FG_graph::create(graph_file, index_file, configs);
	graph_index::ptr index;
	if (test == "push")
		index = NUMA_graph_index<test_push_vertex>::create(
				fg->get_graph_header());
	else if (test == "msg")
		index = NUMA_graph_index<test_msg_vertex>::create(
				fg->get_graph_header());
	else if (test == "get")
		index = NUMA_graph_index<test_get_vertex>::create(
				fg->get_graph_header());
	graph_engine::ptr graph = fg->create_engine(index);
	printf("test starts\n");
	printf("prof_file: %s\n", graph_conf.get_prof_file().c_str());
#ifdef PROFILER
	if (!graph_conf.get_prof_file().empty())
		ProfilerStart(graph_conf.get_prof_file().c_str());
#endif

	struct timeval start, end;
	gettimeofday(&start, NULL);
	graph->start_all();
	graph->wait4complete();
	gettimeofday(&end, NULL);

#ifdef PROFILER
	if (!graph_conf.get_prof_file().empty())
		ProfilerStop();
#endif
	printf("It takes %f seconds\n", time_diff(start, end));
}
