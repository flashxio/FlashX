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
int num_iters = INT_MAX;

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
			const vertex_index &index): compute_directed_vertex(id, index) {
		this->curr_itr_pr = 1 - DAMPING_FACTOR; // Must be this
		this->new_pr = curr_itr_pr;
	}

	float get_curr_itr_pr() const{
		return new_pr;
	}

	void run(vertex_program &prog) { 
		// We perform pagerank for at most `num_iters' iterations.
		if (prog.get_graph().get_curr_level() >= num_iters)
			return;
		directed_vertex_request req(get_id(), edge_type::OUT_EDGE);
		request_partial_vertices(&req, 1);
	};

	void run(vertex_program &, const page_vertex &vertex);

	void run_on_message(vertex_program &, const vertex_message &msg1) {
		const pr_message &msg = (const pr_message &) msg1;
		new_pr += msg.get_delta();
	}
};

void pgrank_vertex::run(vertex_program &prog, const page_vertex &vertex)
{
	int num_dests = vertex.get_num_edges(OUT_EDGE);
	edge_seq_iterator it = vertex.get_neigh_seq_it(OUT_EDGE, 0, num_dests);

	// If this is the first iteration.
	if (prog.get_graph().get_curr_level() == 0) {
		pr_message msg(curr_itr_pr / num_dests * DAMPING_FACTOR);
		prog.multicast_msg(it, msg);
	}
	else if (std::fabs(new_pr - curr_itr_pr) > TOLERANCE) {
		pr_message msg((new_pr - curr_itr_pr) / num_dests * DAMPING_FACTOR);
		prog.multicast_msg(it, msg);
		curr_itr_pr = new_pr;
	}
}

class count_vertex_query: public vertex_query
{
	double num;
public:
	count_vertex_query() {
		num = 0;
	}

	virtual void run(graph_engine &graph, compute_vertex &v) {
		pgrank_vertex &pr_v = (pgrank_vertex &) v;
		num += pr_v.get_curr_itr_pr();
	}

	virtual void merge(graph_engine &graph, vertex_query::ptr q) {
		count_vertex_query *cvq = (count_vertex_query *) q.get();
		num += cvq->num;
	}

	virtual ptr clone() {
		return vertex_query::ptr(new count_vertex_query());
	}

	double get_num() const {
		return num;
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
			"page-rank [options] conf_file graph_file index_file damping_factor\n");
	fprintf(stderr, "-c confs: add more configurations to the system\n");
	fprintf(stderr, "-p: preload the graph\n");
	fprintf(stderr, "-i num: specify the maximal number of iterations\n");
	graph_conf.print_help();
	params.print_help();
}

int main(int argc, char *argv[])
{
	int opt;
	std::string confs;
	int num_opts = 0;
	bool preload = false;
	while ((opt = getopt(argc, argv, "c:pi:")) != -1) {
		num_opts++;
		switch (opt) {
			case 'c':
				confs = optarg;
				num_opts++;
				break;
			case 'p':
				preload = true;
				break;
			case 'i':
				num_iters = atoi(optarg);
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

	signal(SIGINT, int_handler);

	graph_index::ptr index = NUMA_graph_index<pgrank_vertex>::create(index_file);
	graph_engine::ptr graph = graph_engine::create(graph_file, index, configs);
	if (preload)
		graph->preload_graph();
	printf("Pagerank (at maximal %d iterations) starting\n", num_iters);
	printf("prof_file: %s\n", graph_conf.get_prof_file().c_str());
	if (!graph_conf.get_prof_file().empty())
		ProfilerStart(graph_conf.get_prof_file().c_str());


	struct timeval start, end;
	gettimeofday(&start, NULL);
	graph->start_all(); 
	graph->wait4complete();
	gettimeofday(&end, NULL);

	vertex_query::ptr cvq(new count_vertex_query());
	graph->query_on_all(cvq);
	double total = ((count_vertex_query *) cvq.get())->get_num();

	if (!graph_conf.get_prof_file().empty())
		ProfilerStop();
	if (graph_conf.get_print_io_stat())
		print_io_thread_stat();

	printf("The %ld vertices have page rank sum: %lf\n in %f seconds\n", 
			graph->get_num_vertices(), total, time_diff(start, end));
}
