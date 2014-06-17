/**
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
#ifdef PROFILER
#include <google/profiler.h>
#endif

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
vsize_t PREVIOUS_K; 
bool all_greater_than_core = true;

class kcore_vertex: public compute_directed_vertex
{
  bool deleted;
  vsize_t core;
  vsize_t degree; 

public:
	kcore_vertex() {
	}

  kcore_vertex(vertex_id_t id, const vertex_index &index1):
    compute_directed_vertex(id, index1) {
	const directed_vertex_index &index = (const directed_vertex_index &) index1;
    this->degree = index.get_num_in_edges(id) + index.get_num_out_edges(id);
    this->core = degree == 0 ? 0 :
                degree == 1 ? 1 : -1; // Everyone between kmin < core > kmax will get this core
    this->deleted = degree == 0 ? true : false; // If your degree you're already deleted
  }

  bool is_deleted() const {
    return deleted;
  }

  void _delete() {
    this->deleted = true;
  }

  void set_core(vsize_t core) {
    this->core = core;
  }

  const vsize_t get_core() const {
    return this->core;
  }

  vsize_t get_degree() {
    return degree;
  }

  void run(vertex_program &prog);

	void run(vertex_program &prog, const page_vertex &vertex);

	void run_on_message(vertex_program &prog, const vertex_message &msg); 
};

// If I am to be deleted, multicast this message to all my neighbors
// and activate them
class deleted_message: public vertex_message
{
  public:
  deleted_message(): vertex_message(sizeof(deleted_message), true) {
  }
};

void multicast_delete_msg(vertex_program &prog, 
      const page_vertex &vertex, edge_type E)
{
  int num_dests = vertex.get_num_edges(E);
  edge_seq_iterator it = vertex.get_neigh_seq_it(E, 0, num_dests);

  // Doesn't matter who sent it, just --degree on reception 
  deleted_message msg;
  prog.multicast_msg(it, msg);
}

void kcore_vertex::run(vertex_program &prog) {
  if ( degree > CURRENT_K ) { 
    return; 
  }

  if (!is_deleted()) {
    vertex_id_t id = get_id();
    request_vertices(&id, 1); // put my edgelist in page cache

    if (all_greater_than_core) {
      all_greater_than_core = false;
    }
  }
}

void kcore_vertex::run(vertex_program &prog, const page_vertex &vertex) {
  if (is_deleted()) {
    return; // Nothing to be done here
  }

  if ( get_degree() < CURRENT_K ) {
    set_core(CURRENT_K - 1); // This is true because you must make it past CURRENT_K-1 to be here
    _delete();

    // Send two multicast messages - [IN_EDGE, OUT_EDGE] 
    multicast_delete_msg(prog, vertex, IN_EDGE);
    multicast_delete_msg(prog, vertex, OUT_EDGE);
  }
}

void kcore_vertex::run_on_message(vertex_program &prog, const vertex_message &msg) {
  if (is_deleted()) {
    return; // nothing to be done here
  }
  // else
  degree--;
}

class count_vertex_query: public vertex_query
{
	size_t num;
public:
	count_vertex_query() {
		num = 0;
	}

	virtual void run(graph_engine &graph, compute_vertex &v) {
		kcore_vertex &kcore_v = (kcore_vertex &) v;
		if (!kcore_v.is_deleted())
			num++;
	}

	virtual void merge(graph_engine &graph, vertex_query::ptr q) {
		count_vertex_query *cvq = (count_vertex_query *) q.get();
		num += cvq->num;
	}

	virtual ptr clone() {
		return vertex_query::ptr(new count_vertex_query());
	}

	size_t get_num() const {
		return num;
	}
};

// Max degree corresponds to the highest core
class max_degree_query: public vertex_query
{
  vsize_t max_degree;
public:
  max_degree_query() {
    max_degree = 0;
  }

  virtual void run(graph_engine &graph, compute_vertex &v) {
    if (graph.get_vertex_edges(v.get_id()) > max_degree) {
      max_degree = graph.get_vertex_edges(v.get_id());
    }
  }

  virtual void merge(graph_engine &graph, vertex_query::ptr q) {
    max_degree_query *mdq = (max_degree_query *) q.get();
    if (max_degree < mdq->max_degree) {
      max_degree = mdq->max_degree;
    }
  }

  virtual ptr clone() {
    return vertex_query::ptr(new max_degree_query());
  }

  vsize_t get_max_degree() const {
    return max_degree;
  }
};

// Figure out the lowest REMAINING degree in the graph
class min_degree_query: public vertex_query
{
  vsize_t min_degree;
public:
  min_degree_query() {
    min_degree = std::numeric_limits<vsize_t>::max();
  }

  virtual void run(graph_engine &graph, compute_vertex &v) {
    kcore_vertex kcore_v = (kcore_vertex&) v;
    if (!kcore_v.is_deleted()) {
      if (kcore_v.get_degree() < min_degree) {
        min_degree = kcore_v.get_degree();
      }
    }
  }

  virtual void merge(graph_engine &graph, vertex_query::ptr q) {
    min_degree_query *mdq = (min_degree_query *) q.get();
    if (min_degree > mdq->min_degree) {
      min_degree = mdq->min_degree;
    }
  }

  virtual ptr clone() {
    return vertex_query::ptr(new min_degree_query());
  }

  vsize_t get_min_degree() const {
    return min_degree;
  }
};

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
			"k-core [options] conf_file graph_file index_file kmin [kmax] (=Max Degree)\n");
	fprintf(stderr, "-c confs: add more configurations to the system\n");
	graph_conf.print_help();
	params.print_help();
}

// Helpers
void print_func(vertex_id_t i) {
  std::cout << " " << i;
}
void print_active(std::vector<vertex_id_t> v) {
  std::cout << "[";
    for_each (v.begin(), v.end(), print_func);
  std::cout <<  " ]\n";
}
// End Helpers

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
	CURRENT_K = atol(argv[3]); // Set kmin
  
  if (CURRENT_K < 2) {
   fprintf(stderr, "[Error]: kmin cannot be < 2\n");
   exit(-1);
  }

	config_map configs(conf_file);
	configs.add_options(confs);

	signal(SIGINT, int_handler);

	graph_index::ptr index = NUMA_graph_index<kcore_vertex>::create(index_file);
	graph_engine::ptr graph = graph_engine::create(graph_file, index, configs);
	printf("K-core starting\n");
	printf("prof_file: %s\n", graph_conf.get_prof_file().c_str());

#ifdef PROFILER
	if (!graph_conf.get_prof_file().empty())
		ProfilerStart(graph_conf.get_prof_file().c_str());
#endif

  // Set kmax
  vsize_t kmax;
  if (argc > 4) {
    kmax = atol(argv[4]);
  }
  else { // compute largest degree and set it
    printf("Computing kmax as max_degree ...\n");
    vertex_query::ptr mdq(new max_degree_query());
    graph->query_on_all(mdq); 
    kmax = ((max_degree_query *) mdq.get())->get_max_degree();
  }
  printf("Setting kmax to %u ... \n", kmax);


  class activate_k_filter: public vertex_filter {
    vsize_t min;
    public:
    activate_k_filter (vsize_t min) {
      this->min = min;
    }
    bool keep(compute_vertex &v) {
      kcore_vertex &kcore_v = (kcore_vertex &) v;
      return kcore_v.get_degree() < min;
    }
  };

  for (; CURRENT_K <= kmax; CURRENT_K++) {

    std::shared_ptr<vertex_filter> filter
      = std::shared_ptr<vertex_filter>(new activate_k_filter(CURRENT_K));

    struct timeval start, end;
    gettimeofday(&start, NULL);
    graph->start(filter, vertex_program_creater::ptr()); 
    graph->wait4complete();
    gettimeofday(&end, NULL);

    if (all_greater_than_core) { // There's a chance we can hop forward
      printf("***All are greater than K = %u\n", CURRENT_K);
      vertex_query::ptr mdq(new min_degree_query());
      graph->query_on_all(mdq);
      vsize_t min_degree_remaining = ((min_degree_query *) mdq.get())->get_min_degree();

      if (min_degree_remaining == std::numeric_limits<vsize_t>::max()) {
        printf("No more active vertices left!\n");
        break;
      }
      
      printf("\n\nThe graphs minimum degree remaining is %u\n\n", min_degree_remaining);
      // Effectively jumps us to the CURRENT_K + 1th core
      CURRENT_K = min_degree_remaining; // NOTE: Careful - messing with the loop variable :/
      
      if (CURRENT_K > kmax) 
        printf("\nTerminating computation at kmax\n");
        break;
    }
    all_greater_than_core = true;

#if 1
    vertex_query::ptr cvq(new count_vertex_query());
    graph->query_on_all(cvq);
    size_t in_k_core = ((count_vertex_query *) cvq.get())->get_num();
    printf("\n******************************************\n"
        "%d-core shows %ld vertices > %d degree in %f seconds\n"
        "\n******************************************\n",
        CURRENT_K, in_k_core, CURRENT_K, time_diff(start, end));
#endif
    PREVIOUS_K = CURRENT_K;

  }

  // Print out
  graph_index::const_iterator it = index->begin();
  graph_index::const_iterator end_it = index->end();
  printf("Cores by vertex:\n");
  for (; it != end_it; ++it) {
    const kcore_vertex &v = (const kcore_vertex &) *it;
    printf("%u, ", v.get_core());
  }
  printf("\n");

#ifdef PROFILER
	if (!graph_conf.get_prof_file().empty())
		ProfilerStop();
#endif

	if (graph_conf.get_print_io_stat())
		print_io_thread_stat();
}
