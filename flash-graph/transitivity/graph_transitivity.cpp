/**
 * Copyright 2014 Open Connectome Project (http://openconnecto.me)
 * Written by Disa Mhembere (dias@jhu.edu)
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
#include <google/profiler.h>

#include "stat.h"
#include "FGlib.h"
#include "graph.h"

#ifdef _OPENMP
  #include "omp.h"
#endif

// TODO: Add includes for triangles

void int_handler(int sig_num)
{
	if (!graph_conf.get_prof_file().empty())
		ProfilerStop();
	exit(0);
}

void print_usage()
{
	fprintf(stderr,
			"cc [options] conf_file graph_file index_file\n");
	fprintf(stderr, "-c confs: add more configurations to the system\n");
	fprintf(stderr, "-s size: the output min component size\n");
	fprintf(stderr, "-p: preload the graph\n");
	graph_conf.print_help();
	params.print_help();
}

int main(int argc, char *argv[])
{
	int opt;
	std::string confs;
	std::string output_file;
	size_t min_comp_size = 0;
	int num_opts = 0;
	bool preload = false;
	while ((opt = getopt(argc, argv, "c:o:p:")) != -1) {
		num_opts++;
		switch (opt) {
			case 'c':
				confs = optarg;
				num_opts++;
				break;
			case 'o':
				output_file = optarg;
				num_opts++;
				break;
			case 'p':
				preload = true;
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

	FG_graph::ptr graph = FG_graph::create(graph_file, index_file, configs);
  FG_vector<vsize_t> vec = compute_transitivity(graph);
  // TODO: Print me

}
/**
* Transitivity = 2#triangles / (deg * (deg - 1))
**/
FG_vector<vsize_t>::ptr compute_transitivity(FG_graph::ptr fg)
{
  graph_index::ptr index = NUMA_graph_index<wcc_vertex>::create(
      fg->get_index_file());
  graph_engine::ptr graph = graph_engine::create(fg->get_graph_file(),
      index, fg->get_configs());

  struct timeval start, end;
  printf("transitivity starts\n");
  // TODO: Get triangles
  FG_vector<vsize_t> tri_v; // = fg->compute_triangles(); // Vector with triangles
 
  // Compute degree 
  FG_vector<vsize_t> degree_v; // really (deg * (deg - 1))
  vsize_t vid = graph->get_min_vertex_id();

(void) omp_set_num_threads(graph.get_num_threads()); // TODO: Will this work without starting? Prolly not

#pragma omp parallel for private(vid) 
    for (; vid < graph->get_max_vertex_id(); vid++) {
      vsize_t tmp =  graph->get_vertex_edges(vid);
      degree_v[vid] = tmp*(tmp-1); // May be zero 
    }
#pragma omp parallel for private(degree_v, tri_v) // Mutate degree_v to save space
    for (vsize_t i = 0; i < degree_v.size(); i++) {
      tri_v[i] > 0 ? degree_v[i] = 2*tri_v[i]/degree_v[i] : degree_v[i] = 0; 
    }

  printf("Transitivity takes %f seconds\n", time_diff(start, end));
  return degree_v;
}


