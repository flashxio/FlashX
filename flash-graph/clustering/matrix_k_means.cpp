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

// #include <signal.h>
#include <google/profiler.h>
#include <stdlib.h>

#include <vector>
#include <limits> 

#include "thread.h"
#include "io_interface.h"
#include "container.h"
#include "concurrency.h"

#include "vertex_index.h"
#include "FGlib.h"
#include "graph_engine.h"
#include "graph.h"
#include "graph_config.h"
#include "math.h"
#include "matrix/FG_dense_matrix.h"
#include "matrix/matrix_eigensolver.h"

vsize_t g_MAX_ITERS = std::numeric_limits<vsize_t>::max();
std::vector<eigen_pair_t> g_eigen_pairs;
vsize_t K = 1;
vsize_t g_NV;
float g_MIN_EIGV = std::numeric_limits<float>::max();
float g_MAX_EIGV = std::numeric_limits<float>::min();

enum phase_t {
  e_step, /* Update cluster_id for each vertex */
  m_step, /* Update cluster_means for each cluster */
};

// TODO: Fill in
ev_float_t L2Norm(FG_vector<ev_float_t>::ptr, FG_vector<ev_float_t>::ptr) {
  return 0.0;
}

float gen_random_float() {
    float random = ((float) rand()) / (float) RAND_MAX;
    float diff = g_MAX_EIGV - g_MIN_EIGV;
    float r = random * diff;
    return g_MIN_EIGV + r;
}

void print_func(vertex_id_t i) {
  std::cout << " " << i;
}
void print_vector(std::vector<float> v) {
  std::cout << "[";
    for_each (v.begin(), v.end(), print_func);
  std::cout <<  " ]\n";
}

void print_FG_vector(FG_vector<ev_float_t>::ptr fgv) {
  std::cout << "[";
    for (vsize_t i=0; i < fgv->get_size(); i++) {
      std::cout << " " << fgv->get(i);
    }
  std::cout <<  " ]\n\n";
  
}
// End Helpers


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

  FG_graph::ptr _graph = FG_graph::create(graph_file, index_file, configs);

#if 0
  if (matrix_type == "adj") {
    FG_adj_matrix::ptr matrix = FG_adj_matrix::create(_graph);
    // compute_eigen<FG_adj_matrix>(matrix, m, g_NV, which, g_eigen_pairs);
  }
  else
    assert (0);
#endif

# if 1
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
  
  // Make vector for cluster assignments
  FG_vector<vsize_t>::ptr cluster_assignments = 
    FG_vector<vsize_t>::create(g_eigen_pairs.size());
  cluster_assignments->init(INVALID_VERTEX_ID);

  // Create row-wise dense matrix to hold eigs
  // FG_col_wise_matrix<ev_float_t>::ptr eigs = 
  //  FG_col_wise_matrix<ev_float_t>::create(g_eigen_pairs.size(), g_NV); // (n, x)
  FG_row_wise_matrix<ev_float_t>::ptr eigs = 
    FG_row_wise_matrix<ev_float_t>::create(g_eigen_pairs.size(), g_NV); // (n, x)

  // Insert eigs into `eig matrix` and discard eigenvalues
  // FIXME: Slow I assume ...
  vsize_t col = 0;
  BOOST_FOREACH(eigen_pair_t &v, g_eigen_pairs) {
    for (vsize_t row = 0; row < v.second->get_size(); row++) {
      eigs->set(row, col, v.second->get(row));
    }
    col++;
  }
  
  // Create `clusters` dense matrix
  FG_row_wise_matrix<ev_float_t>::ptr clusters = 
    FG_row_wise_matrix<ev_float_t>::create(K, g_NV); // (K, x)

  // Figure out max and min so we can use it to evenly distribute cluster init
  printf("\nNew g_MAX_EIGV is %f ...\n", g_MAX_EIGV);
  printf("\nNew g_MIN_EIGV is %f ...\n", g_MIN_EIGV);

  // Randomly initialize cluster centers
  for (vsize_t row = 0; row < K; row++) {
    FG_vector<ev_float_t>::ptr fgv = FG_vector<ev_float_t>::create(g_NV);

    #pragma omp parallel for 
    for (vsize_t j = 0; j < g_NV; j++) {
      fgv->set(j, gen_random_float());
    }

    clusters->set_row(row, *fgv);
  }
  
  // Create matrix to hold distances
  FG_row_wise_matrix<ev_float_t>::ptr distance_matrix = 
    FG_row_wise_matrix<ev_float_t>::create(eigs->get_num_rows(), K); // (n, K) 
  
	printf("Matrix K-means starting ... \n");
	printf("prof_file: %s\n", graph_conf.get_prof_file().c_str());
	if (!graph_conf.get_prof_file().empty())
		ProfilerStart(graph_conf.get_prof_file().c_str());

  struct timeval start, end;
  gettimeofday(&start, NULL);

  std::cout << "Computing " << g_MAX_ITERS << " iterations\n";
  for (vsize_t iter = 0; iter < g_MAX_ITERS; iter++) {
    printf("Iteration %d ... Computing distance matrix ...\n", iter);
    // Compute distances. TODO: Make coalesced OMP loop

    for (vsize_t vertexID = 0;  vertexID < eigs->get_num_rows(); vertexID++) {
      for (vsize_t clusterID = 0; clusterID < clusters->get_num_rows(); clusterID++) {
        distance_matrix->set(vertexID, clusterID, 
            L2Norm(eigs->get_row_ref(vertexID), clusters->get_row_ref(clusterID)));
      }
    cluster_assignments->set(vertexID, distance_matrix->get_row_ref(vertexID)->argmin());
    }
  
    printf("Clearing cluster means ...\n");
    for (vsize_t clusterID = 0; clusterID < clusters->get_num_rows(); clusterID++) {
      clusters->get_row_ref(clusterID)->assign(0); // TODO: Make row_wise_mat support this in ||
    }

    printf("Updating cluster means ...\n");
    for (vsize_t vertexID = 0; vertexID < cluster_assignments->get_size(); vertexID++) {
       clusters->get_row_ref(cluster_assignments->get(vertexID)) += eigs->get_row_ref(vertexID);
    }
    // TODO: Take the mean of all added means

  }
  gettimeofday(&end, NULL);

  printf("\n******************************************\n"
      "K-means converges in %f seconds\n"
      "\n******************************************\n",
      time_diff(start, end));

  // printf("Printing clusters: \n");
  // print_FG_vector(cluster_assignments);
  
	if (!graph_conf.get_prof_file().empty())
		ProfilerStop();
}
