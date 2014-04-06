/**
 * Copyright 2013 Disa Mhembere, Da Zheng
 *
 * This file is part of FlashGraph.
 *
 * FlashGraph is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * FlashGraph is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with FlashGraph.  If not, see <http://www.gnu.org/licenses/>.
 */

/**
 * This version of pagerank uses message passing only.
 * If the pagerank chagnes, a vertex push the new value to its neighbors.
 */

#include <signal.h>
#include <google/profiler.h>

#include <limits>

#include "thread.h"
#include "io_interface.h"
#include "container.h"
#include "concurrency.h"

#include "vertex_index.h"
#include "graph_engine.h"
#include "graph_config.h"

float DAMPING_FACTOR = 0.85;
float TOLERANCE = 1.0E-2; 

class pr_message: public vertex_message
{
	float delta;
public:
	pr_message(float delta): vertex_message(sizeof(pr_message),
			true) {
		this->delta = delta;
	}

	float get_delta() const {
		return delta;
	}
};

class pgrank_vertex: public compute_directed_vertex
{
	float new_pr;
	float curr_itr_pr; // Current iteration's page rank
public:
	pgrank_vertex() {
		this->curr_itr_pr = 1 - DAMPING_FACTOR; // Must be this
		this->new_pr = curr_itr_pr;
	}

	pgrank_vertex(vertex_id_t id,
			const vertex_index *index): compute_directed_vertex(id, index) {
		this->curr_itr_pr = 1 - DAMPING_FACTOR; // Must be this
		this->new_pr = curr_itr_pr;
	}

	float get_curr_itr_pr() const{
		return curr_itr_pr;
	}

	void run(graph_engine &graph) { 
		directed_vertex_request req(get_id(), edge_type::OUT_EDGE);
		request_partial_vertices(&req, 1);
	};

	void run(graph_engine &graph, const page_vertex &vertex);

	void run_on_message(graph_engine &, const vertex_message &msg1) {
		const pr_message &msg = (const pr_message &) msg1;
		new_pr += msg.get_delta();
	}
};

void pgrank_vertex::run(graph_engine &graph, const page_vertex &vertex)
{
	int num_dests = vertex.get_num_edges(OUT_EDGE);
	stack_array<vertex_id_t, 1024> dest_buf(num_dests);
	vertex.read_edges(OUT_EDGE, dest_buf.data(), num_dests);

	// If this is the first iteration.
	if (graph.get_curr_level() == 0) {
		pr_message msg(curr_itr_pr / num_dests * DAMPING_FACTOR);
		graph.multicast_msg(dest_buf.data(), num_dests, msg);
	}
	else if (std::fabs(new_pr - curr_itr_pr) > TOLERANCE) {
		pr_message msg((new_pr - curr_itr_pr) / num_dests * DAMPING_FACTOR);
		graph.multicast_msg(dest_buf.data(), num_dests, msg);
		curr_itr_pr = new_pr;
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
			"page-rank [options] conf_file graph_file index_file damping_factor\n");
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
	DAMPING_FACTOR = atof(argv[3]);

	if (DAMPING_FACTOR < 0 || DAMPING_FACTOR > 1) {
		fprintf(stderr, "Damping factor must be between 0 and 1 inclusive\n");
		exit(-1);
	}

	config_map configs(conf_file);
	configs.add_options(confs);
	graph_conf.init(configs);
	graph_conf.print();

	signal(SIGINT, int_handler);
	init_io_system(configs);

	graph_index *index = NUMA_graph_index<pgrank_vertex>::create(index_file,
			graph_conf.get_num_threads(), params.get_num_nodes());
	graph_engine *graph = graph_engine::create(graph_conf.get_num_threads(),
			params.get_num_nodes(), graph_file, index);
	printf("Pagerank starting\n");
	printf("prof_file: %s\n", graph_conf.get_prof_file().c_str());
	if (!graph_conf.get_prof_file().empty())
		ProfilerStart(graph_conf.get_prof_file().c_str());


	struct timeval start, end;
	gettimeofday(&start, NULL);
	graph->start_all(); 
	graph->wait4complete();
	gettimeofday(&end, NULL);

	NUMA_graph_index<pgrank_vertex>::const_iterator it
		= ((NUMA_graph_index<pgrank_vertex> *) index)->begin();
	NUMA_graph_index<pgrank_vertex>::const_iterator end_it
		= ((NUMA_graph_index<pgrank_vertex> *) index)->end();

#if 0
	for (; it != end_it; ++it) {
		const pgrank_vertex &v = (const pgrank_vertex &) *it;
		printf("%d:%f\n", v.get_id()+1, v.get_curr_itr_pr());
	}
#endif

#if 1
	float total = 0;
	vsize_t count = 0;

	for (; it != end_it; ++it) {
		const pgrank_vertex &v = (const pgrank_vertex &) *it;
		total += v.get_curr_itr_pr();
		count++;
	}

	if (!graph_conf.get_prof_file().empty())
		ProfilerStop();
	if (graph_conf.get_print_io_stat())
		print_io_thread_stat();
	graph_engine::destroy(graph);
	destroy_io_system();

	printf("The %d vertices have page rank sum: %f\n in %f seconds\n", 
			count, total, time_diff(start, end));
#endif
}
