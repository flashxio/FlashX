/*
 * Copyright 2014 Open Connectome Project (http://openconnecto.me)
 * Written by Disa Mhembere (disa@jhu.edu)
 *
 * This file is part of FlashMatrix.
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

#include "FGlib.h"
#include "matrix/matrix_eigensolver.h"
#include "matrix/FG_sparse_matrix.h"
#include "matrix/kmeans.h"
#include "ts_graph.h"

using namespace fg;

void print_usage();

void int_handler(int sig_num)
{
#ifdef PROFILER
	if (!graph_conf.get_prof_file().empty())
		ProfilerStop();
#endif
	exit(0);
}

template <typename T>
void pairs_to_p_mat(T* eig_matrix, std::vector<eigen_pair_t>& eigen_pairs) 
{
	unsigned vec_cnt = 0;
	BOOST_FOREACH(eigen_pair_t &v, eigen_pairs) {
		for (unsigned row = 0; row < v.second->get_size(); row++) { // 0...NUM_ROWS
			eig_matrix[vec_cnt*v.second->get_size() + row] = v.second->get(row);
		}
		vec_cnt++;
	}
}


/* For testing only */
void dummy_eigs(unsigned nev, size_t g_size, std::vector<eigen_pair_t>& eigen_pairs) {
	// Fake eigen computation for WIKI
	std::cout << "Creating " << nev << " Fake eigs & vecs of len: " << g_size << std::endl;
	for (unsigned i = 0; i < nev; i++) {
		double eigval = (double) random() / (double)RAND_MAX;
		FG_vector<double>::ptr eigvect = 
			FG_vector<double>::create(g_size);

		for (vertex_id_t ii = 0; ii < g_size; ii++) {
			eigvect->set(ii, (double) random() / (double)RAND_MAX); 
		} 

		eigen_pair_t ep(eigval, eigvect);
		eigen_pairs.push_back(ep); 
	}
}

void run_kmeans(FG_graph::ptr graph, int argc, char* argv[])
{
	int opt;
	std::string confs;
	std::string matrix_type = "adj";
	int nev = 1;
	int ncv = 2;
	unsigned k = 1;
	std::string outfile = "";

	unsigned max_iters=std::numeric_limits<unsigned>::max();
	std::string which = "LA";
	std::string init = "forgy";

	int num_opts = 0;

	while ((opt = getopt(argc, argv, "c:k:e:p:w:i:t:o:")) != -1) {
		num_opts++;
		switch (opt) {
			case 'c':
				confs = optarg;
				num_opts++;
				break;
			case 'k':
				k = atol(optarg);
				num_opts++;
				break;
			case 'e': 
				nev = atoi(optarg);
				ncv = 2*nev; // Chosen as a defualt
				num_opts++;
				break;
			case 'p': 
				ncv = atoi(optarg);
				num_opts++;
				break;
			case 'w':
				which = optarg;
				num_opts++;
				break;
			case 'i':
				max_iters = atol(optarg);
				num_opts++;
				break;
			case 't':
				init = optarg;
				num_opts++;
				break;
			case 'o':
				outfile = optarg;
				num_opts++;
				break;
			default:
				print_usage();
		}
	}
	argv += 1 + num_opts;
	argc -= 1 + num_opts;

	if (k == 1) {
		fprintf(stderr, "k must be > 1\n");
		exit(-1);
	}

	struct timeval start, end;

	if (matrix_type == "adj") {
		std::vector<eigen_pair_t> eigen_pairs;
		FG_adj_matrix::ptr matrix = FG_adj_matrix::create(graph);

#if 0
		compute_eigen<FG_adj_matrix>(matrix, ncv, nev, which, eigen_pairs);
#else
		dummy_eigs(nev, matrix->get_num_rows(), eigen_pairs);
#endif
		// Convert the eigen_pairs to a flattened matrix
		double* p_eig_matrix = new double[matrix->get_num_rows()*nev];
		pairs_to_p_mat(p_eig_matrix, eigen_pairs);

		/* Malloc */
		double* p_clusters = new double [k*nev];
		unsigned* p_clust_asgns = new unsigned [matrix->get_num_rows()];
		unsigned* p_clust_asgn_cnt = new unsigned [k];
		/* End Malloc */

		gettimeofday(&start, NULL);
		unsigned iters = compute_kmeans(p_eig_matrix, p_clusters, p_clust_asgns, 
				p_clust_asgn_cnt, matrix->get_num_rows(), nev, k, max_iters, init);


		gettimeofday(&end, NULL);

		printf("Kmeans took %u iters and %.3f sec\n", iters, time_diff(start, end));

		printf("Printing cluster assignment counts:\n");
		printf("[ ");
		for (unsigned i = 0; i < k; i++) {
			std::cout << p_clust_asgn_cnt[i] << " ";
		}
		printf("]\n");
		
		/* Reclamation */
		printf("Freeing p_eig_matrix\n");
		delete [] p_eig_matrix;

		printf("Freeing p_clust_asgns\n");
		delete [] p_clust_asgns;

		printf("Freeing p_clust_asgn_cnt\n");
		delete [] p_clust_asgn_cnt;

		printf("Freeing p_clusters\n");
		delete [] p_clusters;
	} else {
		assert (0);
	}
}

std::string supported_algs[] = {
	"kmeans",
};
int num_supported = sizeof(supported_algs) / sizeof(supported_algs[0]);

void print_usage()
{
	fprintf(stderr,
			"test_matrx conf_file graph_file index_file algorithm [alg-options]\n");
    fprintf(stderr, "-c confs: add more configurations to the system\n");
    fprintf(stderr, "-p num: the number of Lancos base vectors\n");
    fprintf(stderr, "-e num: the number of actual eigenvalues\n");
    fprintf(stderr, "-k num: the number of clusters desired\n");
    fprintf(stderr, "-w which: which side of eigenvalues\n");
    fprintf(stderr, "-t type: type of initialization for kmeans ['random', 'forgy', 'kmeanspp']\n");
    fprintf(stderr, "-i iters: maximum number of iterations\n");

	fprintf(stderr, "supported graph algorithms:\n");
	for (int i = 0; i < num_supported; i++)
		fprintf(stderr, "\t%s\n", supported_algs[i].c_str());
	graph_conf.print_help();
	safs::params.print_help();
}

int main(int argc, char *argv[])
{
	argv++;
	argc--;
	if (argc < 4) {
		print_usage();
		exit(-1);
	}

	std::string conf_file = argv[0];
	std::string graph_file = argv[1];
	std::string index_file = argv[2];
	std::string alg = argv[3];
	// We should increase by 3 instead of 4. getopt() ignores the first
	// argument in the list.
	argv += 3;
	argc -= 3;

	config_map::ptr configs = config_map::create(conf_file);
	if (configs == NULL)
		configs = config_map::ptr();
	signal(SIGINT, int_handler);

	FG_graph::ptr graph = FG_graph::create(graph_file, index_file, configs);

	if (alg == "kmeans") {
		run_kmeans(graph, argc, argv);
	}
}
