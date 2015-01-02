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

#define USE_EIGEN

#include "matrix/FG_sparse_matrix.h"
#include "matrix/matrix_eigensolver.h"
#include "matrix/ASE.h"

using namespace fg;

typedef FG_sparse_matrix<general_get_edge_iter<char> > FG_general_sparse_char_matrix;

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
			"cc [options] conf_file graph_file index_file\n");
	fprintf(stderr, "-c confs: add more configurations to the system\n");
	fprintf(stderr, "-m type: the type of matrix\n");
	fprintf(stderr, "-p num: the number of potential eigenvalues\n");
	fprintf(stderr, "-k num: the number of real eigenvalues\n");
	fprintf(stderr, "-w which: which side of eigenvalues\n");
	fprintf(stderr, "-t type: the type of eigenvlaues\n");
}

int main(int argc, char *argv[])
{
	int opt;
	std::string confs;
	std::string matrix_type = "adj";
	int nv = 1;
	int m = 2;
	std::string which = "LA";
	std::string type = "EV";
	
	int num_opts = 0;
	while ((opt = getopt(argc, argv, "c:m:k:p:w:t:")) != -1) {
		num_opts++;
		switch (opt) {
			case 'c':
				confs = optarg;
				num_opts++;
				break;
			case 'm':
				matrix_type = optarg;
				num_opts++;
				break;
			case 'k':
				nv = atoi(optarg);
				num_opts++;
				break;
			case 'p':
				m = atoi(optarg);
				num_opts++;
				break;
			case 'w':
				which = optarg;
				num_opts++;
				break;
			case 't':
				type = optarg;
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

	config_map::ptr configs = config_map::create(conf_file);
	assert(configs);
	configs->add_options(confs);

	signal(SIGINT, int_handler);

	printf("compute %d eigenvalues\n", nv);
	assert(nv < m);
	FG_graph::ptr graph = FG_graph::create(graph_file, index_file, configs);
	std::vector<eigen_pair_t> eigen_pairs;
	// adjacency matrix
	if (matrix_type == "adj") {
		FG_adj_matrix::ptr matrix = FG_adj_matrix::create(graph);
		if (type == "EV")
			compute_eigen<FG_adj_matrix>(matrix, m, nv, which, eigen_pairs);
		else if (type == "LS" || type == "RS")
			compute_SVD<FG_adj_matrix>(matrix, m, nv, which, type, eigen_pairs);
	}
	// general sparse matrix of characters
	else if (matrix_type == "gc") {
		FG_general_sparse_char_matrix::ptr matrix
			= FG_general_sparse_char_matrix::create(graph);
		if (type == "EV")
			compute_eigen<FG_general_sparse_char_matrix>(matrix, m, nv, which, eigen_pairs);
		else if (type == "LS" || type == "RS")
			compute_SVD<FG_general_sparse_char_matrix>(matrix, m, nv, which, type, eigen_pairs);
	}
	else if (matrix_type == "AcD") {
		size_t nvertices = graph->get_graph_header().get_num_vertices();
		compute_AcD_uw(graph, 1.0/nvertices, m, nv, which, eigen_pairs);
	}
	else
		assert(0);

	BOOST_FOREACH(eigen_pair_t &v, eigen_pairs) {
		printf("value: %f, ||vector||: %f\n", v.first, v.second->norm2());
	}
}
