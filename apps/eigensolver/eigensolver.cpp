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

#include <signal.h>
#include <google/profiler.h>

#include <vector>

#include "thread.h"
#include "io_interface.h"
#include "container.h"
#include "concurrency.h"

#include "vertex_index.h"
#include "graph_engine.h"
#include "graph_config.h"

float beta = 0;

/**
 * This eigensolver implements the Lanzcos algorithm.
 */

class eigen_vertex: public compute_vertex
{
	// v of the current iteration
	float v1;
	// v of the previous iteration
	float v0;
	float w;
public:
	eigen_vertex() {
		v0 = 0;
		v1 = 0;
		w = 0;
	}

	eigen_vertex(vertex_id_t id,
			const vertex_index &index): compute_vertex(id, index) {
		v0 = 0;
		v1 = sqrt(1.0 / index.get_num_vertices());
		w = 0;
	}

	void init_eigen() {
		v0 = v1;
		v1 = w / beta;
	}

	float adjust_w(float alpha) {
		w = w - alpha * v1 - beta * v0;
		return w;
	}

	float get_w() const {
		return w;
	}

	float get_v() const {
		return v1;
	}

	void run(vertex_program &prog) {
		vertex_id_t id = get_id();
		request_vertices(&id, 1);
	}

	void run(vertex_program &prog, const page_vertex &vertex);

	void run_on_message(vertex_program &prog, const vertex_message &msg) {
	}
};

class eigen_vertex_program: public vertex_program_impl<eigen_vertex>
{
	float alpha;
public:
	typedef std::shared_ptr<eigen_vertex_program> ptr;

	static ptr cast2(vertex_program::ptr prog) {
		return std::static_pointer_cast<eigen_vertex_program, vertex_program>(
				prog);
	}

	eigen_vertex_program() {
		alpha = 0;
	}

	void add_vertex(const eigen_vertex &v) {
		alpha += v.get_w() * v.get_v();
	}

	float get_alpha() const {
		return alpha;
	}
};

void eigen_vertex::run(vertex_program &prog, const page_vertex &vertex)
{
	edge_seq_iterator it = vertex.get_neigh_seq_it(edge_type::BOTH_EDGES, 0,
			vertex.get_num_edges(edge_type::BOTH_EDGES));
	w = 0;
	PAGE_FOREACH(vertex_id_t, id, it) {
		const eigen_vertex &eigen_v = (const eigen_vertex &) prog.get_graph().get_vertex(id);
		w += eigen_v.v1;
	} PAGE_FOREACH_END

	((eigen_vertex_program &) prog).add_vertex(*this);
}

class eigen_vertex_initiator: public vertex_initiator
{
public:
	void init(compute_vertex &v) {
		eigen_vertex &ev = (eigen_vertex &) v;
		ev.init_eigen();
	}
};

class eigen_vertex_program_creater: public vertex_program_creater
{
public:
	vertex_program::ptr create() const {
		return vertex_program::ptr(new eigen_vertex_program());
	}
};

class w_query: public vertex_query
{
	const float alpha;
	float w_sq_sum;
public:
	w_query(float _alpha): alpha(_alpha) {
		w_sq_sum = 0;
	}

	virtual void run(graph_engine &graph, compute_vertex &v) {
		eigen_vertex &eigen_v = (eigen_vertex &) v;
		float w = eigen_v.adjust_w(alpha);
		w_sq_sum += w * w;
	}

	virtual void merge(graph_engine &graph, vertex_query::ptr q) {
		w_sq_sum += ((w_query *) q.get())->w_sq_sum;
	}

	virtual ptr clone() {
		return vertex_query::ptr(new w_query(alpha));
	}

	float get_new_beta() {
		return sqrt(w_sq_sum);
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
			"eigensolver [options] conf_file graph_file index_file\n");
	fprintf(stderr, "-c confs: add more configurations to the system\n");
	fprintf(stderr, "-p: preload the graph\n");
	fprintf(stderr, "-m: the dimension of the tridiagonal matrix\n");
	graph_conf.print_help();
	params.print_help();
}

int main(int argc, char *argv[])
{
	int opt;
	std::string confs;
	int num_opts = 0;
	bool preload = false;
	int m = 0;
	while ((opt = getopt(argc, argv, "c:pm:")) != -1) {
		num_opts++;
		switch (opt) {
			case 'c':
				confs = optarg;
				num_opts++;
				break;
			case 'p':
				preload = true;
				break;
			case 'm':
				m = atoi(optarg);
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

	signal(SIGINT, int_handler);

	graph_index::ptr index = NUMA_graph_index<eigen_vertex>::create(index_file);
	graph_engine::ptr graph = graph_engine::create(graph_file, index, configs);
	if (preload)
		graph->preload_graph();
	printf("Eigensolver starts\n");
	printf("prof_file: %s\n", graph_conf.get_prof_file().c_str());
	if (!graph_conf.get_prof_file().empty())
		ProfilerStart(graph_conf.get_prof_file().c_str());

	printf("There are %d iterations\n", m);
	for (int i = 0; i < m; i++) {
		struct timeval start, end;
		gettimeofday(&start, NULL);
		if (beta)
			graph->start_all(vertex_initiator::ptr(new eigen_vertex_initiator()),
					vertex_program_creater::ptr(new eigen_vertex_program_creater()));
		else
			graph->start_all(vertex_initiator::ptr(),
					vertex_program_creater::ptr(new eigen_vertex_program_creater()));
		graph->wait4complete();

		// Compute alpha_j = w_j . v_j
		std::vector<vertex_program::ptr> vprogs;
		graph->get_vertex_programs(vprogs);
		float alpha = 0;
		BOOST_FOREACH(vertex_program::ptr vprog, vprogs) {
			eigen_vertex_program::ptr eigen_vprog = eigen_vertex_program::cast2(vprog);
			alpha += eigen_vprog->get_alpha();
		}

		// Compute w_j = w_j - alpha_j * v_j - beta_j * v_j-1
		// beta_j+1 = || w_j ||
		vertex_query::ptr wq(new w_query(alpha));
		graph->query_on_all(wq);
		beta = ((w_query *) wq.get())->get_new_beta();
		printf("a%d: %f, b%d: %f\n", i, alpha, i + 1, beta);

		gettimeofday(&end, NULL);
		printf("Iteration %d takes %f seconds\n", i, time_diff(start, end));
	}

	if (!graph_conf.get_prof_file().empty())
		ProfilerStop();
	if (graph_conf.get_print_io_stat())
		print_io_thread_stat();
}
