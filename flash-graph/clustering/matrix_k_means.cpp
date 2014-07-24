/*
 * Copyright 2014 Open Connectome Project (http://openconnecto.me)
 * Written by Disa Mhembere (disa@jhu.edu)
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
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY CURRENT_KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

#include <stdlib.h>
#include <math.h>
#ifdef PROFILER
#include <google/profiler.h>
#endif
#include <vector>
#include <limits> 

#include "FGlib.h"
#include "matrix/FG_sparse_matrix.h"
#include "matrix/FG_dense_matrix.h"
#include "matrix/matrix_eigensolver.h"

vsize_t g_MAX_ITERS = std::numeric_limits<vsize_t>::max();
std::vector<eigen_pair_t> g_eigen_pairs;
vsize_t K = 1;
vsize_t g_NV;
float g_MIN_EIGV = std::numeric_limits<float>::max();
float g_MAX_EIGV = std::numeric_limits<float>::min();

/**
 * \brief Get Euclidean distance given two vectors.
 */
ev_float_t L2Norm(FG_vector<ev_float_t>::ptr v1, FG_vector<ev_float_t>::ptr v2) 
{
	assert(v1->get_size() == v2->get_size());
	ev_float_t sum = 0;
	for (size_t idx = 0; idx < v1->get_size(); idx++) {
		ev_float_t diff = v1->get(idx) - v2->get(idx);
		sum += diff*diff;
	}
	if (sum == 0) { return 0; }   
	return sqrt(sum);
}

float gen_random_float() 
{
	float random = ((float) rand()) / (float) RAND_MAX;
	float diff = g_MAX_EIGV - g_MIN_EIGV;
	float r = random * diff;
	return g_MIN_EIGV + r;
}

std::vector<ev_float_t> dist_vector; // For memoization
// Serial
size_t phi(ev_float_t c) {
	ev_float_t min_dist = std::numeric_limits<ev_float_t>::max();
	for(size_t i = 0; i < g_eigen_pairs.size(); i++) {
		for (size_t j = 0; j < g_eigen_pairs[i].second->get_size(); j++) {
			ev_float_t eucl_dist = fabs(g_eigen_pairs[i].second->get(j) - c);
			if (eucl_dist < min_dist) {
				min_dist = eucl_dist;
			}
		}
	}	
	return min_dist;
}

// FIXME DM: Complete me
float get_prob(FG_vector<ev_float_t>::ptr C) {
	return 0.5;
}

// \param prob The probability of selection
// FIXME: Complete me
ev_float_t get_indep_sample(float prob) {
	return 0.0;
}

ev_float_t get_random_sample() 
{
	// Sample a point at random from eigs
	size_t eigs_size = g_eigen_pairs.size(); // To ensure compiler optimization for / and %
	size_t sample_size = eigs_size * (g_eigen_pairs[0].second->get_size() - 1);
	size_t rand_sample = rand() % sample_size;

	// g_eigen_pairs[row][col]
	return g_eigen_pairs[rand_sample / eigs_size].second->get(rand_sample % eigs_size);
}

void k_means_bar_bar_init_clusters(float oversampling_factor, FG_row_wise_matrix<ev_float_t>::ptr & clusters) 
{
	std::vector<ev_float_t> C;
	for (size_t clusterID = 0; clusterID < clusters->get_num_rows(); clusterID++) {
		// Sample a point at random from eigs
		C.push_back(get_random_sample());
		dist_vector.push_back(phi(C.back()));

		size_t iters = (size_t) ceil(log2(dist_vector.back()));

		for (size_t i = 0; i < iters; i++) {
			C.push_back(get_indep_sample(get_prob(clusters->get_row_ref(clusterID)))); // Line 4
			// FIXME: Complete me	
		}
	}
}

// Begin Helpers
	template <typename T>
void print_vector(typename std::vector<T> v) 
{
	std::cout << "[";

	typename std::vector<T>::iterator itr = v.begin();
	for (; itr != v.end(); itr++) {
		std::cout << " "<< *itr;
	}
	std::cout <<  " ]\n";
}

void print_usage()
{
	fprintf(stderr,
			"k-means [options] conf_file graph_file index_file K\n");
	fprintf(stderr, "-c confs: add more configurations to the system\n");
	fprintf(stderr, "-m type: the type of matrix\n");
	fprintf(stderr, "-p num: the number of potential eigenvalues\n");
	fprintf(stderr, "-k num: the number of real eigenvalues\n");
	fprintf(stderr, "-w which: which side of eigenvalues\n");
	fprintf(stderr, "-t type: the type of eigenvlaues\n"); 
	fprintf(stderr, "-i iters: maximum number of iterations\n"); 

	graph_conf.print_help();
	params.print_help();
}

int main(int argc, char *argv[])
{
	int opt;
	std::string confs;
	std::string matrix_type = "adj";
	g_NV = 1;
	int m = 2; m = m;
	std::string which = "LA";

	int num_opts = 0;

	while ((opt = getopt(argc, argv, "c:m:k:p:w:i:")) != -1) {
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
			case 'k': /**/
				g_NV = atoi(optarg);
				num_opts++;
				break;
			case 'p': /**/
				m = atoi(optarg);
				num_opts++;
				break;
			case 'w':
				which = optarg;
				num_opts++;
				break;
			case 'i':
				g_MAX_ITERS = atol(optarg);
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
	K = atol(argv[3]);

	config_map configs(conf_file);
	configs.add_options(confs);

	FG_graph::ptr _graph = FG_graph::create(graph_file, index_file, configs);

#if 1
	if (matrix_type == "adj") {
		FG_adj_matrix::ptr matrix = FG_adj_matrix::create(_graph);
		compute_eigen<FG_adj_matrix>(matrix, m, g_NV, which, g_eigen_pairs);
	}
	else
		assert (0);
#endif

# if 0
	// Fake eigen computation for WIKI
	std::cout << "Creating Fake eigs ..." << std::endl;
	for (vsize_t i = 0; i < g_NV; i++) {
		ev_float_t eigval = (ev_float_t) rand() / (ev_float_t)RAND_MAX;
		FG_vector<ev_float_t>::ptr eigvect = FG_vector<ev_float_t>::create(8298);

		for (vertex_id_t ii = 0; ii < 8298; ii++) {
			eigvect->set(ii, (ev_float_t) rand() / (ev_float_t)RAND_MAX); 
		} 

		eigen_pair_t ep(eigval, eigvect);
		g_eigen_pairs.push_back(ep); 
	}
# endif

	// Create row-wise dense matrix to hold eigs
	FG_row_wise_matrix<ev_float_t>::ptr eigs = 
		FG_row_wise_matrix<ev_float_t>::create(g_eigen_pairs[0].second->get_size(), g_NV); // (n, x)
	eigs->resize(g_eigen_pairs[0].second->get_size(), g_NV);

	// Make vector for cluster assignments
	FG_vector<vsize_t>::ptr cluster_assignments = 
		FG_vector<vsize_t>::create(eigs->get_num_rows());
	cluster_assignments->init(INVALID_VERTEX_ID);

	FG_vector<vsize_t>::ptr prev_cluster_assignments = 
		FG_vector<vsize_t>::create(eigs->get_num_rows());
	cluster_assignments->init(INVALID_VERTEX_ID);

#ifdef PROFILER
	ProfilerStart("/home/disa/graph-engine/flash-graph/clustering/profile_out.log");
#endif


	// Insert eigs into `eig matrix` and discard eigenvalues
	// FIXME: Slow .. I assume ...
	vsize_t col = 0;

	BOOST_FOREACH(eigen_pair_t &v, g_eigen_pairs) {
		for (vsize_t row = 0; row < v.second->get_size(); row++) {
			eigs->set(row, col, v.second->get(row));

			// For random init we need a baseline
			ev_float_t tmp = v.second->min();
			if (tmp < g_MIN_EIGV) { g_MIN_EIGV = tmp; }
			tmp = v.second->max();
			if (tmp > g_MAX_EIGV) { g_MAX_EIGV = tmp; }
		}
		col++;
	}

	// Create `clusters` dense matrix
	FG_row_wise_matrix<ev_float_t>::ptr clusters = 
		FG_row_wise_matrix<ev_float_t>::create(K, g_NV); // (K, x)
	clusters->resize(K, g_NV);

	// Figure out max and min so we can use it to evenly distribute cluster init
	printf("\nNew g_MAX_EIGV is %f ...\n", g_MAX_EIGV);
	printf("New g_MIN_EIGV is %f ...\n", g_MIN_EIGV);

	// Randomly initialize cluster centers
	for (vsize_t clusterID = 0; clusterID < K; clusterID++) {
		FG_vector<ev_float_t>::ptr fgv = FG_vector<ev_float_t>::create(g_NV);

#pragma omp parallel for 
		for (vsize_t j = 0; j < g_NV; j++) {
			fgv->set(j, gen_random_float());
		}

		clusters->set_row(clusterID, fgv);
	}

	// Create matrix to hold distances
	FG_row_wise_matrix<ev_float_t>::ptr distance_matrix = 
		FG_row_wise_matrix<ev_float_t>::create(eigs->get_num_rows(), K); // (n, K) 
	distance_matrix->resize(eigs->get_num_rows(), K);

	printf("Matrix K-means starting ... \n");

	struct timeval start, end;
	gettimeofday(&start, NULL);

	std::cout << "Computing " << g_MAX_ITERS << " iterations\n";
	for (vsize_t iter = 0; iter < g_MAX_ITERS; iter++) {
		// Hold cluster assignment counter
		std::vector<vsize_t> cluster_assignment_counts;
		cluster_assignment_counts.assign(K, 0);

		printf("E-step Iteration %d ... Computing distance matrix ...\n", iter);
		// Compute distances. 
		// TODO: Make coalesced OMP loop

		for (vsize_t vertexID = 0; vertexID < eigs->get_num_rows(); vertexID++) {
			for (vsize_t clusterID = 0; clusterID < clusters->get_num_rows(); clusterID++) {

				// printf("Cluster ID: %u has mean ", clusterID);
				// clusters->get_row_ref(clusterID)->print();

				distance_matrix->set(vertexID, clusterID, 
						L2Norm(eigs->get_row_ref(vertexID), clusters->get_row_ref(clusterID)));
			}

			// printf("Vertex %u has distance vector: ", vertexID); 
			// distance_matrix->get_row_ref(vertexID)->print();

			size_t assigned_cluster = distance_matrix->get_row_ref(vertexID)->argmin();

			// printf("Vertex %u with eig: %f, is assigned cluster %lu\n", vertexID, eigs->get(vertexID, 0), assigned_cluster);

			cluster_assignments->set(vertexID, assigned_cluster);
			cluster_assignment_counts[assigned_cluster]++;
		}

		// printf("Previous clusters: "); prev_cluster_assignments->print();
		// printf("\nNew clusters: "); cluster_assignments->print();

		if (cluster_assignments->eq_all(prev_cluster_assignments)) { // eq_all is not ||
			printf("K-means converged in %u iterations\n", iter);
			break;
		}

		printf("Cluster assignment counts: \n");
		print_vector(cluster_assignment_counts);

		// I clear the cluster mean matrix because I don't want to use more memory
		printf("Clearing cluster means ...\n");
		clusters->assign_all(0.0); // Not ||

		printf("M-step Updating cluster means ...\n");
		for (vsize_t vertexID = 0; vertexID < cluster_assignments->get_size(); vertexID++) {
			// printf("Vertex %u ...\n", vertexID);
			// plus_eq is ||
			(clusters->get_row_ref(cluster_assignments->get(vertexID)))->plus_eq( eigs->get_row_ref(vertexID) );
		}
		// Take the mean of all added means
		printf("Div by in place ...\n");
		for (vsize_t clusterID = 0; clusterID < clusters->get_num_rows(); clusterID++) {
			// div_by_in_place is ||
			clusters->get_row_ref(clusterID)->div_by_in_place
				(cluster_assignment_counts[clusterID]); 
		}

		// Set prev to current
		prev_cluster_assignments->shallow_copy(cluster_assignments); // shallow_copy is ||
	}

#ifdef PROFILER
	ProfilerStop();
#endif
	gettimeofday(&end, NULL);

	printf("\n******************************************\n"
			"K-means complete in %f seconds\n"
			"\n******************************************\n",
			time_diff(start, end));

	// printf("Printing clusters: \n");
	// cluster_assignments->print();
}
