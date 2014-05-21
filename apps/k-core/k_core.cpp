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

class kcore_vertex: public compute_directed_vertex
{
  bool deleted;
  vsize_t degree; 

public:
	kcore_vertex() {
	}

  kcore_vertex(vertex_id_t id, const vertex_index &index1):
    compute_directed_vertex(id, index1) {
    this->deleted = false;
	const directed_vertex_index &index = (const directed_vertex_index &) index1;
    this->degree = index.get_num_in_edges(id) + index.get_num_out_edges(id);
  }

  bool is_deleted() const {
    return deleted;
  }

  void _delete() {
    this->deleted = true;
  }

  vsize_t get_degree() {
    return degree;
  }

  void run(vertex_program &prog) {
    if (degree > CURRENT_K) { return; }

    if (!is_deleted()) {
      std::cout << "In run Vertex id: " << get_id() << " has degree:" << get_degree() << std::endl;
			vertex_id_t id = get_id();
			request_vertices(&id, 1); // put my edgelist in page cache
    }
  }

	void run(vertex_program &prog, const page_vertex &vertex);

	void run_on_message(vertex_program &prog, const vertex_message &msg); 
};

// If I am to be deleted, multicast this message to all my neighbors
// and ACTIVATE them
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


// Per thread ...
class kcore_vertex_program: public vertex_program_impl<kcore_vertex>
{
	std::vector<vertex_id_t> next_engine_active_vert; // One for each thread
public:
  typedef std::shared_ptr<kcore_vertex_program> ptr;

	kcore_vertex_program() {
	}

	void activate_next_engine(vertex_id_t id) {
		next_engine_active_vert.push_back(id);
	}

  const std::vector<vertex_id_t>& get_next_engine_active_vert() const{
    return next_engine_active_vert;
  } 

  static ptr cast2(vertex_program::ptr prog) {
    return std::static_pointer_cast<kcore_vertex_program, vertex_program>(prog);
  }
};

class kcore_vertex_program_creater: public vertex_program_creater
{
public:
	vertex_program::ptr create() const {
		return vertex_program::ptr(new kcore_vertex_program());
	}
};

// This is only run ONCE per iteration ?
void kcore_vertex::run(vertex_program &prog, const page_vertex &vertex) {
  //std::cout << "Vertex " << get_id() <<" has degree " << degree << "\n";

  if(is_deleted()) {
    return; // Nothing to be done here
  }

  if (get_degree() < CURRENT_K) {
    _delete();
    std::cout << "Vertex " << get_id() << " being deleted\n"; 
   
    // Send two multicast messages - [IN_EDGE, OUT_EDGE] 
    multicast_delete_msg(prog, vertex, IN_EDGE);
    multicast_delete_msg(prog, vertex, OUT_EDGE);
  }
  
  // This is how to tell who to activiate next engine.
  // This SHOULD happen exactly one.
 /* else if (get_degree() == CURRENT_K + 1) {
    ((kcore_vertex_program&) prog).activate_next_engine(vertex.get_id());
 }*/
}

void kcore_vertex::run_on_message(vertex_program &, const vertex_message &msg) {
  if (is_deleted()) {
    assert(degree == 0);
    return; // nothing to be done here
  }
  
  std::cout << "In rom My degree is " << degree << std::endl;
  assert(degree > 0);
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
		if (kcore_v.is_deleted())
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
    vsize_t degree = ((kcore_vertex&) v).get_degree();
    if (degree < min_degree) {
      min_degree = degree;
    }
  }

  virtual void merge(graph_engine &graph, vertex_query::ptr q) {
    min_degree_query *mdq = (min_degree_query *) q.get();
    if (min_degree > mdq->min_degree) { // TODO: Verify
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
	if (!graph_conf.get_prof_file().empty())
		ProfilerStop();
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
	vsize_t CURRENT_K = atol(argv[3]); // Set kmin

	config_map configs(conf_file);
	configs.add_options(confs);

	signal(SIGINT, int_handler);

	graph_index::ptr index = NUMA_graph_index<kcore_vertex>::create(index_file);
	graph_engine::ptr graph = graph_engine::create(graph_file, index, configs);
	printf("K-core starting\n");
	printf("prof_file: %s\n", graph_conf.get_prof_file().c_str());
	if (!graph_conf.get_prof_file().empty())
		ProfilerStart(graph_conf.get_prof_file().c_str());

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
  
  // Filter by those whose degree < kmin
  class activate_k_filter: public vertex_filter {
    vsize_t min;
    public:
    activate_k_filter (vsize_t min) {
      this->min = min;
    }
    bool keep(compute_vertex &v) {
      kcore_vertex &kcore_v = (kcore_vertex &) v;
      std::cout << "In filter Vertex id: " << v.get_id() << " has degree:" << kcore_v.get_degree() << std::endl;
      return kcore_v.get_num_in_edges() + kcore_v.get_num_out_edges() < min;
    }
  };

  // First K only
  std::shared_ptr<vertex_filter> filter
    = std::shared_ptr<vertex_filter>(new activate_k_filter(CURRENT_K));

  
  struct timeval start, end;
  gettimeofday(&start, NULL);
  graph->start(filter); 
  graph->wait4complete();
  gettimeofday(&end, NULL);

  vertex_query::ptr cvq(new count_vertex_query());
  graph->query_on_all(cvq);
  size_t in_k_core = ((count_vertex_query *) cvq.get())->get_num();
  printf("\n******************************************\n"
      "%d-core shows %ld vertices <= %d degree in %f seconds\n"
      "\n******************************************\n",
      CURRENT_K, in_k_core, CURRENT_K, time_diff(start, end));

  printf("I got here!");
  exit(EXIT_FAILURE);

  // Subsequent K's
  if (CURRENT_K + 1 < kmax) { 
    for (; CURRENT_K <= kmax; CURRENT_K++) {
    std::vector<vertex_id_t> active_vertices;

      std::vector<vertex_program::ptr> vertex_progs_th;
      graph->get_vertex_programs(vertex_progs_th);

      BOOST_FOREACH(vertex_program::ptr vprog, vertex_progs_th) {
        kcore_vertex_program::ptr kcore_vprog = kcore_vertex_program::cast2(vprog);
        active_vertices.insert(active_vertices.begin(),
            kcore_vprog->get_next_engine_active_vert().begin(),
            kcore_vprog->get_next_engine_active_vert().end());
      }
      
      std::cout << "The following are active in the next engine: " << std::endl;
      print_active(active_vertices);
      
      // This means the CURRENT_K is an empty core 
      // so hop to next viable CURRENT_K 
      if (active_vertices.empty()) {
        vertex_query::ptr mdq(new min_degree_query());
        graph->query_on_all(mdq); 
        CURRENT_K = ((min_degree_query *) mdq.get())->get_min_degree();

        if (CURRENT_K > kmax) break;
      }

      struct timeval start, end;
      gettimeofday(&start, NULL);
      graph->start(&active_vertices[0], active_vertices.size()); // TODO: Verify if OK
      graph->wait4complete();
      gettimeofday(&end, NULL);

      vertex_query::ptr cvq(new count_vertex_query());
      graph->query_on_all(cvq);
      size_t in_k_core = ((count_vertex_query *) cvq.get())->get_num();
      printf("\n******************************************\n"
          "%d-core shows %ld vertices <= %d degree in %f seconds\n"
          "\n******************************************\n",
          CURRENT_K, in_k_core, CURRENT_K, time_diff(start, end));
    }
  }

	if (!graph_conf.get_prof_file().empty())
		ProfilerStop();
	if (graph_conf.get_print_io_stat())
		print_io_thread_stat();
}
