/**
 * Copyright 2013 Disa Mhembere
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

// #include <vector>
#include <sstream>
#include <fstream>
#include <limits>

#include "thread.h"
#include "io_interface.h"
#include "container.h"
#include "concurrency.h"

#include "vertex_index.h"
#include "graph_engine.h"
#include "graph_config.h"

float DAMPING_FACTOR;
// Same as powergraph
float TOLERANCE = 1.0E-2; 
// size_t ITERATIONS = 0;

// TODO: Maybe use a count of num_out_neighs who have sent messages
// to know when I can update the prev_itr_pr. So I Know I have all
// neighbors I need
vsize_t converged_count = 0;

class pgrank_vertex: public compute_directed_vertex
{
  float prev_itr_pr; // Previous iteration's page rank
  float curr_itr_pr; // Current iteration's page rank

public:
	pgrank_vertex() {
	}

  pgrank_vertex(vertex_id_t id, const vertex_index *index): compute_directed_vertex(id, index) {
    this->prev_itr_pr = std::numeric_limits<float>::min(); // Must be < -1
    this->curr_itr_pr = 1 - DAMPING_FACTOR; // Must be this
  }

  float get_prev_itr_pr() const{
    return prev_itr_pr;
  }

  float get_curr_itr_pr() const{
    return curr_itr_pr;
  }

  void run(graph_engine &graph) { 
    vertex_id_t id = get_id();
    request_vertices(&id, 1); // put my edgelist in page cache
  };

	void run(graph_engine &graph, const page_vertex &vertex);

	virtual void run_on_messages(graph_engine &,
			const vertex_message *msgs[], int num) { }; // Only serves to activate on the next iteration
};

// Alerts every out-neighbor that my page rank has changed
// and they should update their own page rank accordingly
class pgrank_message: public vertex_message
{
  int dummy; // TODO: see if rm this still works
  public:
  pgrank_message( ): 
        vertex_message(sizeof(pgrank_message), true) { // Always activate. Only place vertices are activated
        }
};

void pgrank_vertex::run(graph_engine &graph, const page_vertex &vertex) {

#if 0
  page_byte_array::const_iterator<vertex_id_t> end_it
    = vertex.get_neigh_end(OUT_EDGE);
  for (page_byte_array::const_iterator<vertex_id_t> it
      = vertex.get_neigh_begin(OUT_EDGE); it != end_it; ++it) {
    vertex_id_t id = *it;
    printf("%d %d\n", get_id(), id);
  } 
#endif

  // Gather
  float accum = 0;
  page_byte_array::const_iterator<vertex_id_t> end_it
    = vertex.get_neigh_end(IN_EDGE);
  
  for (page_byte_array::const_iterator<vertex_id_t> it
      = vertex.get_neigh_begin(IN_EDGE); it != end_it; ++it) {
    vertex_id_t id = *it;
    pgrank_vertex& v = (pgrank_vertex&) graph.get_vertex(id);
    accum += (v.get_curr_itr_pr()/v.get_num_out_edges())  ; // Notice I want this iterations pagerank
  }   

  // Apply
  if (get_num_in_edges() > 0) {
    float update_pr = ((1 - DAMPING_FACTOR)) + (DAMPING_FACTOR*(accum)); // FIXME: Many of these may be expensive
    curr_itr_pr = update_pr;
  }
  prev_itr_pr = curr_itr_pr; // Happens no matter what

  // Scatter (activate your out-neighbors ... if you have any :) 
  if ( std::fabs( curr_itr_pr - prev_itr_pr ) > TOLERANCE ) {
    page_byte_array::const_iterator<vertex_id_t> end_it
      = vertex.get_neigh_end(OUT_EDGE);
    stack_array<vertex_id_t, 1024> dest_buf(vertex.get_num_edges(OUT_EDGE));
    int num_dests = 0;
    for (page_byte_array::const_iterator<vertex_id_t> it
        = vertex.get_neigh_begin(OUT_EDGE); it != end_it; ++it) {
      vertex_id_t id = *it;
      dest_buf[num_dests++] = id; 
    }   

    if (num_dests > 0) {
      graph.multicast_msg(dest_buf.data(), num_dests, pgrank_message( )) ;
    }   
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
			"k-core [options] conf_file graph_file index_file damping_factor\n");
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
  
#if 1
	for (; it != end_it; ++it) {
		const pgrank_vertex &v = (const pgrank_vertex &) *it;
    printf("%d:%f\n", v.get_id()+1, v.get_curr_itr_pr());
	}
#endif

#if 0
  // Write pgrank to file for comparison
  std::ostringstream dict_str;
	for (; it != end_it; ++it) {
		const pgrank_vertex &v = (const pgrank_vertex &) *it;
    dict_str << v.get_id() << ":" << v.get_curr_itr_pr() << "\n";
	}
  std::ofstream outFile;
  outFile.open("/home/disa/graph-engine/apps/page-rank/compare/wiki-Vote-pagerank.edge");
  outFile << dict_str.rdbuf();
  outFile.close();
#endif

#if 0
	float mean_pgrank = 0;
  vsize_t count = 0;

	for (; it != end_it; ++it) {
		const pgrank_vertex &v = (const pgrank_vertex &) *it;
    mean_pgrank += v.get_curr_itr_pr();
    count++;
	}

	if (!graph_conf.get_prof_file().empty())
		ProfilerStop();
	if (graph_conf.get_print_io_stat())
		print_io_thread_stat();
	graph_engine::destroy(graph);
	destroy_io_system();
  printf("Mean pgrank: %f\n", mean_pgrank);
  printf("Count: %d\n", count);

  printf("The average page rank of %d vertices is: %f\n", count, (mean_pgrank/count));
#endif
}
