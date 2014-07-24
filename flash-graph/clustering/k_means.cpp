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

#include <signal.h>
#include <google/profiler.h>
#include <stdlib.h>
#include <math.h>

#include <vector>
#include <limits> 

#include "thread.h"
#include "io_interface.h"
#include "container.h"
#include "concurrency.h"

#include "vertex_index.h"
#include "graph_engine.h"
#include "graph_config.h"
#include "matrix/FG_sparse_matrix.h"
#include "matrix/matrix_eigensolver.h"

vsize_t g_MAX_ITERS = std::numeric_limits<vsize_t>::max();
std::vector<eigen_pair_t> g_eigen_pairs;
vsize_t K = 1;
int g_NV;
ev_float_t g_MIN_EIGV = std::numeric_limits<ev_float_t>::max();
ev_float_t g_MAX_EIGV = std::numeric_limits<ev_float_t>::min();
bool g_changed_membership = false; // If any node changes membership

float gen_random_float() {
	float random = ((float) rand()) / (float) RAND_MAX;
	float diff = g_MAX_EIGV - g_MIN_EIGV;
	float r = random * diff;
	return g_MIN_EIGV + r;
}

// Helpers
template <typename T>
void print_vector(typename std::vector<T> v) {
	std::cout << "[";

	typename std::vector<T>::iterator itr = v.begin();
	for (; itr != v.end(); itr++) {
		std::cout << " "<< *itr;
	}
	std::cout <<  " ]\n";
}

class cluster 
{
	// Cluster_id implicit in its position in the vector<cluster>

	std::vector<float> curr_cluster_mean;
	std::vector<float> next_cluster_mean; // Cluster mean for the next round of K-means
	vsize_t num_members;
	pthread_spinlock_t lock;

	public:

	cluster() {
		num_members = 0;
		pthread_spin_init(&lock, PTHREAD_PROCESS_PRIVATE);
		// Random init mean
		for (int i = 0; i < g_NV; i++) {
			curr_cluster_mean.push_back(gen_random_float());
		}
	}

	void set_curr_cluster_mean(std::vector<float>& v) {
		curr_cluster_mean = v;
	}

	void add_member(vertex_id_t id) {
		pthread_spin_lock(&lock);

		if ( next_cluster_mean.empty() ) { next_cluster_mean.resize(g_NV); }

		for (vsize_t i = 0; i < next_cluster_mean.size(); i++) {
			next_cluster_mean[i] += g_eigen_pairs[i].second->get(id);
		}
		num_members++;
		pthread_spin_unlock(&lock);
	}

	void update_cluster_mean() {
		float l_num_members = (float)num_members; // Necessary for omp fp
		if (num_members > 0) {
#pragma omp parallel for firstprivate(l_num_members)
			for (vsize_t i=0; i < next_cluster_mean.size(); i++) {
				next_cluster_mean[i] /= l_num_members;
			}
		}

		printf("Curr cluster_mean: ");
		print_vector(curr_cluster_mean);

		printf("Next cluster mean: ");
		print_vector(next_cluster_mean);

		// Shallow copy vector
		curr_cluster_mean = next_cluster_mean;

		next_cluster_mean.assign(next_cluster_mean.size(), 0); 
		num_members = 0; // reset num members
	}

	std::vector<float>& get_mean() {
		return this->curr_cluster_mean;
	}
};

/* VarDecl */ 
std::vector<cluster> g_clusters;

float L2Norm(std::vector<float>& cluster_mean, vertex_id_t id) {
	float sum = 0;
	for (vsize_t i=0; i < cluster_mean.size(); i++) {
		float diff = cluster_mean[i] - g_eigen_pairs[i].second->get(id); 
		sum += diff * diff;
	}

	if (sum == 0) { return sum; }
	return sqrt(sum);
}

class kmeans_vertex: public compute_vertex
{
	uint32_t cluster_id; 

	public:
	kmeans_vertex() {
		cluster_id = rand() % K; // Every vertex is randomly assigned cluster initially
		g_clusters[cluster_id].add_member(get_id());
	}

	kmeans_vertex(vertex_id_t id, const vertex_index &index1):
		compute_vertex(id, index1) { }

	vsize_t compute_cluster() {
		float min_diff = std::numeric_limits<float>::max();
		uint32_t computed_cluster = std::numeric_limits<uint32_t>::max();

		for (uint32_t i = 0; i < g_clusters.size(); i++) {
			float diff = L2Norm(g_clusters[i].get_mean(), get_id());
			if (diff < min_diff ) { 
				min_diff = diff; 
				computed_cluster = i;
			}
		}

		if (computed_cluster != cluster_id) {
			g_changed_membership = true; // May cause race, but OK
		}
		return computed_cluster;
	}

	const vsize_t get_cluster() const {
		return cluster_id;
	}

	void run(vertex_program &prog);
	void run(vertex_program &prog, const page_vertex &vertex) { }
	void run_on_message(vertex_program &prog, const vertex_message &msg) { }
};

// M-step
void recompute_cluster_means() {
	// #pragma omp parallel for // NOTE: No need for omp here because update_cluster_mean is ||
	for (vsize_t idx = 0; idx < g_clusters.size(); idx++) {
		g_clusters[idx].update_cluster_mean();
	}
}

void kmeans_vertex::run(vertex_program &prog) {
	cluster_id = compute_cluster();
	g_clusters[cluster_id].add_member(get_id());
}

void assemble_clusters(graph_index::ptr& index, std::vector<float>& v) {

	// Print out
	graph_index::const_iterator it = index->begin();
	graph_index::const_iterator end_it = index->end();
	printf("Cluster by vertex:\n[ ");

	for (; it != end_it; ++it) { 
		const kmeans_vertex &v = (const kmeans_vertex &) *it;
		std::cout << v.get_cluster() << " ";
	}
	printf("]\n");

	for (vertex_id_t i=0; i < v.size(); i++) {
		v[i] = ((kmeans_vertex&)(index->get_vertex(i))).get_cluster();
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

	signal(SIGINT, int_handler);

# if 1
	FG_graph::ptr _graph = FG_graph::create(graph_file, index_file, configs);

	if (matrix_type == "adj") {
		FG_adj_matrix::ptr matrix = FG_adj_matrix::create(_graph);
		compute_eigen<FG_adj_matrix>(matrix, m, g_NV, which, g_eigen_pairs);
	}
	else
		assert (0);
# endif

# if 0
	// Fake eigen computation for WIKI
	std::cout << "Fake eigs" << std::endl;
	for (int i = 0; i < g_NV; i++) {
		ev_float_t eigval = (ev_float_t) rand() / (ev_float_t)RAND_MAX;
		FG_vector<ev_float_t>::ptr eigvect = FG_vector<ev_float_t>::create(8298);

		for (vertex_id_t ii = 0; ii < 8298; ii++) {
			eigvect->set(ii, (ev_float_t) rand() / (ev_float_t)RAND_MAX); 
		} 

		eigen_pair_t ep(eigval, eigvect);
		g_eigen_pairs.push_back(ep); 
	}
# endif

	// Compute eig mins, maxs
	BOOST_FOREACH(eigen_pair_t &v, g_eigen_pairs) {
		printf("value: %f, ||vector||: ", v.first);
		ev_float_t max = v.second->max();
		ev_float_t min = v.second->min();

		if (max > g_MAX_EIGV)
			g_MAX_EIGV = max;
		if (min < g_MIN_EIGV)
			g_MIN_EIGV = min;
		// v.second->print();
	}

	// Figure out max and min so we can use it to evenly distribute cluster init
	printf("\nNew g_MAX_EIGV is %f ...\n", g_MAX_EIGV);
	printf("\nNew g_MIN_EIGV is %f ...\n", g_MIN_EIGV);

	// Instantiate clusters
	for (vsize_t i = 0; i < K; i++) {
		cluster cl;
		g_clusters.push_back(cl);
	}

	graph_index::ptr index = NUMA_graph_index<kmeans_vertex>::create(index_file);
	graph_engine::ptr graph = graph_engine::create(graph_file, index, configs);
	printf("K-means starting\n");
	printf("prof_file: %s\n", graph_conf.get_prof_file().c_str());
	if (!graph_conf.get_prof_file().empty())
		ProfilerStart("/home/disa/graph-engine/flash-graph/clustering/profile_out.log");

	struct timeval start, end;
	gettimeofday(&start, NULL);


	printf("Computing %u iterations\n", g_MAX_ITERS);
	for (vsize_t iter = 0; iter < g_MAX_ITERS; iter++) {

		printf("\n**********Iteration %u **********\n", iter);
		graph->start_all();
		graph->wait4complete();
		if (!g_changed_membership && iter > 0) {
			printf("K-means converged!\n");
			break;
		}
		else {
			recompute_cluster_means();
			g_changed_membership = false; // reset
		}
	}

	gettimeofday(&end, NULL);

	std::vector<float> cluster_vector( graph->get_max_vertex_id() -
			graph->get_min_vertex_id() + 1, 0.0); 

	assemble_clusters(index, cluster_vector);

	// printf("Printing clusters: \n");
	// print_vector(cluster_vector);

	printf("\n*************************************\n"
			"K-means took %f seconds\n"
			"\n*************************************\n", time_diff(start, end));

	if (!graph_conf.get_prof_file().empty())
		ProfilerStop();
}
