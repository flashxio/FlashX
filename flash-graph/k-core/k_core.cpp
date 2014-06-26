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
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
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

vsize_t K; // Min degree necessary to be part of the k-core graph

class kcore_vertex: public compute_directed_vertex
{
  bool deleted;
  vsize_t degree; 

public:
  kcore_vertex(vertex_id_t id): compute_directed_vertex(id) {
    this->deleted = false;
	assert(0);
    this->degree = get_num_edges();
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
    if (degree > K) { return; }

    if (!is_deleted()) {
			vertex_id_t id = get_id();
			request_vertices(&id, 1); // put my edgelist in page cache
    }
  }

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
		deleted_message msg;
      prog.multicast_msg(dest_buf.data(), num_dests, msg);
    }
}

// This is only run by 1st iteration active vertices
void kcore_vertex::run(vertex_program &prog, const page_vertex &vertex) {
  if(is_deleted()) {
    return; // Nothing to be done here
  }

  if ( get_degree() < K ) {
    _delete();
   
    // Send two multicast messages - [IN_EDGE, OUT_EDGE] 
    multicast_delete_msg(prog, vertex, IN_EDGE);
    multicast_delete_msg(prog, vertex, OUT_EDGE);
    
  }
}

void kcore_vertex::run_on_message(vertex_program &, const vertex_message &msg) {
  if (is_deleted()) {
    return; // nothing to be done here
  }

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
			"k-core [options] conf_file graph_file index_file K\n");
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
	K = atoi(argv[3]); // set K

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

  // Filter for activation first time around
  class activate_k_filter: public vertex_filter {
    vsize_t min;
    public:
    activate_k_filter (vsize_t min) {
      this->min = min;
    }
    bool keep(compute_vertex &v) {
      kcore_vertex &kcore_v = (kcore_vertex &) v;
      return kcore_v.get_num_in_edges() + kcore_v.get_num_out_edges() < min;
    }
  };

	std::shared_ptr<vertex_filter> filter
		= std::shared_ptr<vertex_filter>(new activate_k_filter(K));

	struct timeval start, end;
	gettimeofday(&start, NULL);
	graph->start(filter); 
  graph->wait4complete();
	gettimeofday(&end, NULL);

	vertex_query::ptr cvq(new count_vertex_query());
	graph->query_on_all(cvq);
	size_t in_k_core = ((count_vertex_query *) cvq.get())->get_num();

#ifdef PROFILER
	if (!graph_conf.get_prof_file().empty())
		ProfilerStop();
#endif
	printf("\n%d-core shows %ld vertices > %d degree in %f seconds\n",
			K, in_k_core, K, time_diff(start, end));
}
