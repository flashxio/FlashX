/**
 * Copyright 2013 Disa Mhembere
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

#include <signal.h>
#include <google/profiler.h>

#include <vector>
#include <algorithm>

#include "thread.h"
#include "io_interface.h"
#include "container.h"
#include "concurrency.h"

#include "vertex_index.h"
#include "graph_engine.h"
#include "graph_config.h"

vsize_t CURRENT_K; // Min degree necessary to be part of the k-core graph

class global_variable
{
	volatile size_t value;
	pthread_spinlock_t lock;
public:
	global_variable() {
		value = 0;
		pthread_spin_init(&lock, PTHREAD_PROCESS_PRIVATE);
	}

	global_variable(size_t init) {
		value = init;
		pthread_spin_init(&lock, PTHREAD_PROCESS_PRIVATE);
	}

	void update(size_t new_v) {
		pthread_spin_lock(&lock);
    value = new_v;
		pthread_spin_unlock(&lock);
	}

  void minus_minus() {
    pthread_spin_lock(&lock);
    --value;
    pthread_spin_unlock(&lock);
  }

	size_t get() const {
		return value;
	}
} IN_CORE_VERTICES;

class kcore_vertex: public compute_directed_vertex
{
  bool deleted;
  vsize_t degree; 

public:
	kcore_vertex() {
	}

  kcore_vertex(vertex_id_t id, const vertex_index *index1):
    compute_directed_vertex(id, index1) {
    this->deleted = false;
	directed_vertex_index *index = (directed_vertex_index *) index1;
    this->degree = index->get_num_in_edges(id) + index->get_num_out_edges(id);
  }

  bool is_deleted() const {
    return deleted;
  }

  void _delete() {
    this->deleted = true;
    IN_CORE_VERTICES.minus_minus(); // Lower vertices still in CURRENT_K-core 
  }

  vsize_t get_degree() {
    return degree;
  }

  void run(graph_engine &graph) {
    if (degree > CURRENT_K) { return; }

    if (!is_deleted()) {
			vertex_id_t id = get_id();
			request_vertices(&id, 1); // put my edgelist in page cache
    }
  }

	void run(graph_engine &graph, const page_vertex &vertex);

	virtual void run_on_message(graph_engine &, const vertex_message &msg); 
};

// If I am to be deleted, multicast this message to all my neighbors
// and activate them
class deleted_message: public vertex_message
{
  public:
  deleted_message(): vertex_message(sizeof(deleted_message), true) {
  }
};

void multicast_delete_msg(graph_engine &graph, 
      const page_vertex &vertex, edge_type E)
{
    page_byte_array::const_iterator<vertex_id_t> end_it
      = vertex.get_neigh_end(E);
    stack_array<vertex_id_t, 1024> dest_buf(vertex.get_num_edges(E));
    int num_dests = 0;
    for (page_byte_array::const_iterator<vertex_id_t> it
        = vertex.get_neigh_begin(E); it != end_it; ++it) {
      vertex_id_t id = *it;
      dest_buf[num_dests++] = id;
    } 
    // Doesn't matter who sent it, just --degree on reception 
    if (num_dests > 0) {
      graph.multicast_msg(dest_buf.data(), num_dests, deleted_message());
    }
}

// This is only run by 1st iteration active vertices
void kcore_vertex::run(graph_engine &graph, const page_vertex &vertex) {
  if(is_deleted()) {
    return; // Nothing to be done here
  }

  if ( get_degree() < CURRENT_K ) {
    _delete();
   
    // Send two multicast messages - [IN_EDGE, OUT_EDGE] 
    multicast_delete_msg(graph, vertex, IN_EDGE);
    multicast_delete_msg(graph, vertex, OUT_EDGE);
    
  }
}

void kcore_vertex::run_on_message(graph_engine &, const vertex_message &msg) {
  if (is_deleted()) {
    return; // nothing to be done here
  }

  degree--;
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
			"k-core [options] conf_file graph_file index_file kmin [kmax] (Default=Max Degree)\n");
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
	vsize_t kmin = atol(argv[3]); // set kmin


	config_map configs(conf_file);
	configs.add_options(confs);
	graph_conf.init(configs);
	graph_conf.print();

	signal(SIGINT, int_handler);
	init_io_system(configs);

	graph_index *index = NUMA_graph_index<kcore_vertex>::create(index_file,
			graph_conf.get_num_threads(), params.get_num_nodes());
	graph_engine *graph = graph_engine::create(graph_conf.get_num_threads(),
			params.get_num_nodes(), graph_file, index);
	printf("K-core starting\n");
	printf("prof_file: %s\n", graph_conf.get_prof_file().c_str());
	if (!graph_conf.get_prof_file().empty())
		ProfilerStart(graph_conf.get_prof_file().c_str());
  
  vsize_t kmax; 
  if (argc > 4) {
    kmax = atol(argv[4]);
  }
  else { // compute largest degree and set it
    printf("Computing kmax as max_degree ...\n");
    kmax = 0; 

    NUMA_graph_index<kcore_vertex>::const_iterator it
      = ((NUMA_graph_index<kcore_vertex> *) index)->begin();
    NUMA_graph_index<kcore_vertex>::const_iterator end_it
      = ((NUMA_graph_index<kcore_vertex> *) index)->end();

    // Get max degree
    for (; it != end_it; ++it) {
      const kcore_vertex &v = (const kcore_vertex &) *it;
      if ((v.get_num_in_edges() + v.get_num_out_edges()) > kmax) 
        kmax++;
    } 
  } 
  printf("Setting kmax to %u ... \n", kmax);

  IN_CORE_VERTICES.update(graph->get_max_vertex_id()+1); // Set active vertices
  printf("IN_CORE_VERTICES starts at %lu\n", IN_CORE_VERTICES.get());
  struct timeval start, end;
  gettimeofday(&start, NULL);

  for (CURRENT_K = kmin; CURRENT_K <= kmax; CURRENT_K++) {
    graph->start_all(); 
    graph->wait4complete();

    printf("\n******************************************\n K = %u has %lu" 
        " vertices >= %u\n******************************************\n",
        CURRENT_K, IN_CORE_VERTICES.get(), CURRENT_K); 
  }
  gettimeofday(&end, NULL);

	if (!graph_conf.get_prof_file().empty())
		ProfilerStop();
	if (graph_conf.get_print_io_stat())
		print_io_thread_stat();
	graph_engine::destroy(graph);
	destroy_io_system();
  printf("K-core completed in %f seconds\n", time_diff(start, end));
}
