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
#ifdef PROFILER
#include <google/profiler.h>
#endif

#include <vector>

#include "thread.h"
#include "io_interface.h"
#include "container.h"
#include "concurrency.h"

#include "vertex_index.h"
#include "graph_engine.h"
#include "graph_config.h"
#include "worker_thread.h"

std::vector<int> thread_delays;

class count_message: public vertex_message
{
	int count;
public:
	count_message(): vertex_message(sizeof(count_message), false) {
		count = 1;
	}

	int get_count() const {
		return count;
	}
};

class test_vertex: public compute_vertex
{
	size_t num_edges;
	size_t add1_count;
public:
	test_vertex(vertex_id_t id): compute_vertex(id) {
		num_edges = 0;
		add1_count = 0;
	}

	void init() {
		add1_count = 0;
	}

	void verify_result() const {
		assert(add1_count == num_edges * 2);
	}

	void run(vertex_program &prog) {
		vertex_id_t id = get_id();
		request_vertices(&id, 1);
	}

	void run(vertex_program &prog, const page_vertex &vertex);

	virtual void run_on_message(vertex_program &, const vertex_message &msg1) {
		const count_message &msg = (const count_message &) msg1;
		add1_count += msg.get_count();
	}
};

void test_vertex::run(vertex_program &prog, const page_vertex &vertex)
{
	worker_thread *t = (worker_thread *) thread::get_curr_thread();
	int worker_id = t->get_worker_id();
	if ((size_t) worker_id < thread_delays.size() / 2) {
		if (thread_delays[worker_id] == 0) {
			int slp_time = random() % 10;
			printf("worker %d sleeps %d seconds\n", worker_id, slp_time);
			sleep(slp_time);
			thread_delays[worker_id] = 1;
		}
	}
	// We need to add the neighbors of the vertex to the queue of
	// the next level.
	page_byte_array::const_iterator<vertex_id_t> end_it
		= vertex.get_neigh_end(BOTH_EDGES);
	num_edges = vertex.get_num_edges(BOTH_EDGES);
	stack_array<vertex_id_t, 1024> dest_buf(num_edges);
	int num_dests = 0;
	for (page_byte_array::const_iterator<vertex_id_t> it
			= vertex.get_neigh_begin(BOTH_EDGES); it != end_it; ++it) {
		vertex_id_t id = *it;
		dest_buf[num_dests++] = id;
		add1_count++;
	}
	count_message msg;
	prog.multicast_msg(dest_buf.data(), num_dests, msg);
}

class verify_vertex_query: public vertex_query
{
public:
	virtual void run(graph_engine &graph, compute_vertex &v) {
		const test_vertex &test_v = (const test_vertex &) v;
		test_v.verify_result();
	}

	virtual void merge(graph_engine &graph, vertex_query::ptr q) {
	}

	virtual ptr clone() {
		return vertex_query::ptr(new verify_vertex_query());
	}
};

int main(int argc, char *argv[])
{
	if (argc < 4) {
		fprintf(stderr, "test conf_file graph_file index_file");
		exit(-1);
	}

	std::string conf_file = argv[1];
	std::string graph_file = argv[2];
	std::string index_file = argv[3];

	config_map configs(conf_file);

	graph_index::ptr index = NUMA_graph_index<test_vertex>::create(index_file);
	graph_engine::ptr graph = graph_engine::create(graph_file, index, configs);
	printf("prof_file: %s\n", graph_conf.get_prof_file().c_str());
#ifdef PROFILER
	if (!graph_conf.get_prof_file().empty())
		ProfilerStart(graph_conf.get_prof_file().c_str());
#endif
	thread_delays.resize(graph->get_num_threads());

	for (int i = 0; i < 1000; i++) {
		printf("test %d\n", i);
		struct timeval start, end;
		gettimeofday(&start, NULL);
		graph->start_all();
		graph->wait4complete();
		gettimeofday(&end, NULL);
		printf("It takes %f seconds\n", time_diff(start, end));

		vertex_query::ptr vvq(new verify_vertex_query());
		graph->query_on_all(vvq);
	}

#ifdef PROFILER
	if (!graph_conf.get_prof_file().empty())
		ProfilerStop();
#endif
}
