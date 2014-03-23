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

#include <vector>

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
size_t ITERATIONS = 0;

vsize_t converged_count = 0;

class pgrank_vertex: public compute_directed_vertex
{
  bool first_itr; // Is it the first iteration
  bool converged; // Has the vertex page rank converged
  float prev_itr_pr; // Previous iteration's page rank
  float curr_itr_pr; // Current iteration's page rank
  vsize_t num_diverged_out_neighs; // num of out-edge vertices that haven't converged
  bool sent_converged_msg; // Tell in-neigh vertices I've converged 

public:
	pgrank_vertex() {
	}

  pgrank_vertex(vertex_id_t id, const vertex_index *index): compute_directed_vertex(id, index) {
    first_itr = true;
    converged = false;
    sent_converged_msg = false;
    this->prev_itr_pr = 1 - DAMPING_FACTOR;
    this->curr_itr_pr = 0; // Must be 0
    this->num_diverged_out_neighs = get_num_out_edges();
  }

  float get_prev_itr_pr() {
    return prev_itr_pr;
  }

  float get_curr_itr_pr() const{
    return curr_itr_pr;
  }

  /**
  * Vertex is only complete once it has converged AND
  * all its out-edge neighbors have also converged
  **/
  bool is_complete() {
    return (converged && (num_diverged_out_neighs == 0) && sent_converged_msg);
  }

  void run(graph_engine &graph); 
	void run(graph_engine &graph, const page_vertex &vertex);

	virtual void run_on_messages(graph_engine &,
			const vertex_message *msgs[], int num); 
};

// If my page rank changed since the last iteration I need to tell my neighbors
// Also used to tell out-edge neighbors you've converged
class pgrank_message: public vertex_message
{
  float norm_pr; // My new page rank
  bool converged;
  vsize_t sender_id; // FIXME: rm testing
  public:
  // FIXME: rm `sender_id' testing
  pgrank_message(float norm_pr, bool converged, vsize_t sender_id, bool activate): 
        vertex_message(sizeof(pgrank_message), activate) {
    this->norm_pr = norm_pr;
    this->converged = converged;
    this->sender_id = sender_id; // FIXME: rm testing
  }

  vsize_t get_sender_id() { return sender_id; } // FIXME: rm testing

  float* get_norm_pr() {
    return &norm_pr;
  }

  bool* has_converged() {
    return &converged;
  }
};

void multicast_pgrank_msg(graph_engine &graph, const page_vertex &vertex,
              float norm_pr, bool complete, edge_type E, bool activate)
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

    if (num_dests > 0) {
      graph.multicast_msg(dest_buf.data(), num_dests, pgrank_message(norm_pr, complete, vertex.get_id(), activate)); // FIXME: rm `vertex.get_id()' testing
    }
}

void pgrank_vertex::run(graph_engine &graph) {
  if (is_complete()) { return; }
  // else
  if (first_itr) {
    first_itr = false;
  }
  else {
    if ( !converged ) { 
      float new_prev_itr_pr = 1 - DAMPING_FACTOR + curr_itr_pr;
      if (std::fabs(new_prev_itr_pr - prev_itr_pr) < TOLERANCE) {
        printf("Vertex %d converged with page rank = %f\n", get_id(), new_prev_itr_pr);
        converged = true; // Only place convergence is tested
        printf("%d vertices converged!\n", ++converged_count);
      } else { printf("Vertex: %d changed its PR from %f to %f\n", get_id(), prev_itr_pr, new_prev_itr_pr); } // FIXME: rm testing

      prev_itr_pr = new_prev_itr_pr;
      curr_itr_pr = 0; // reset this value
    }
  }

  vertex_id_t id = get_id();
  request_vertices(&id, 1); // put my edgelist in page cache
}

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
  if (is_complete()) { return; }
  // else
  float norm_pr = get_prev_itr_pr()/get_num_out_edges(); // FIXME: Many of these may be expensive
  multicast_pgrank_msg(graph, vertex, norm_pr, false, OUT_EDGE, true);

  if (!(sent_converged_msg) && converged) {
    // Send message to all in-edge neighs to say decrease
    // your num_diverged_out_neighs by 1
    multicast_pgrank_msg(graph, vertex, 0, true, IN_EDGE, false);
    printf("Vertex %d has sent converged message\n", get_id());
    sent_converged_msg = true;
  }
}

void pgrank_vertex::run_on_messages(graph_engine &,
    const vertex_message *msgs[], int num) {
#if 0
  // FIXME: Testing
  if (get_id() == 3073) {
    printf("Vertex: %d received: %d messages\n", get_id(), num);
    for (int i = 0; i < num; i++) {
      pgrank_message* msg = (pgrank_message*) msgs[i];
      printf("Message sent by Vertex: %d, sender converged? %d \n", msg->get_sender_id(), ( *msg->has_converged() == true ? 1 : 0));
    }
    // printf("\n");

    assert(false);
  }
  // End Testing
#endif

  for (int i = 0; i < num; i++) {
    pgrank_message* msg = (pgrank_message*) msgs[i];
    this->curr_itr_pr += *(msg->get_norm_pr());
    if (msg->has_converged()) { 
      this->num_diverged_out_neighs--; // This is a message from one of my out-edge neighs
      // if (num_diverged_out_neighs < 0) { assert(false)}; // FIXME: rm testing
    }
  }
  this->curr_itr_pr *= DAMPING_FACTOR;
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
}
