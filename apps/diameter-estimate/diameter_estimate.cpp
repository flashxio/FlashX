/**
 * Copyright 2013 Disa Mhembere
 *
 * This file is part of FlashGraph
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
// This file computes an estimated diameter by running multiple BFSs 
// in parallel then taking the max to give a course estimate of 
// of the graphs diameter.

#include <signal.h>
#include <google/profiler.h>

#include <vector>
#include <algorithm>
#include <bitset>

#include "thread.h"
#include "io_interface.h"
#include "container.h"
#include "concurrency.h"

#include "vertex_index.h"
#include "graph_engine.h"
#include "graph_config.h"

uint32_t MAX_PARR_BFS = 20; // Change if you want more parallel bfs
const uint32_t NUM_PARR_SEARCH = 5;
uint32_t MAX_EST_DIAM = 0;
std::vector<vertex_id_t> START_VERTICES;

// helper
void print_arr(vertex_id_t arr[], uint32_t size) {
  printf("[ ");
  for (uint32_t i = 0; i < size; i++) {
    printf("%d ", arr[i]);
  }
  printf(" ]");
}

class depth_message: public vertex_message
{
  std::bitset<NUM_PARR_SEARCH> bfs_bits; // The index in the bitmap
  uint32_t depth;

  public:
  depth_message(std::bitset<NUM_PARR_SEARCH>& bfs_bits, uint32_t depth): 
              vertex_message(sizeof(depth_message), true) {
    this->bfs_bits = bfs_bits;
    this->depth = depth;
  }

  std::bitset<NUM_PARR_SEARCH> get_bits() {
    return this->bfs_bits;
  }

  uint32_t get_depth() {
    return this->depth;
  }
};

class diam_vertex: public compute_directed_vertex {

  // I need to know:
  //  1. which which messages to ignore if already sent
  //  2. when to stop activating neighbors 
  uint32_t depth;
  bool active;
  bool discovered; // Some bfs has seem me
  std::bitset<NUM_PARR_SEARCH> bfs_bits;
public:
	diam_vertex() {
    this->depth = 0; //  A vertex need only track its **maximum depth**
    this->active = true;

    if (is_start_vertex()) { 
      // I'm a start vertex, my depth = 0 and set my bit in bfs_bits
      std::vector<vertex_id_t>::iterator low = 
        std::lower_bound(START_VERTICES.begin(), START_VERTICES.end(), get_id());
      this->bfs_bits.set((uint32_t) (low-START_VERTICES.begin())); // really NUM_PARR_SEARCH - # but doesn't matter
    }

	}

	diam_vertex(vertex_id_t id,
			const vertex_index *index): compute_directed_vertex(id, index) {
	}
  
  // I need to know if I'm a start vertex or not. Only ever called by constructor
  bool is_start_vertex() {
    return (std::binary_search(START_VERTICES.begin(), START_VERTICES.end(), get_id()));
  }

  bool is_active() {
    return active;
  }

  void deactivate() {
    this->active = false;
  }

  uint32_t get_depth() const{
    return this->depth;
  }
 
  bool visit(const std::bitset<NUM_PARR_SEARCH> &msg_bits) const;

  // TODO: decide if I want to update MAX_DEPTHS on the fly or not
  void update_depth(uint32_t depth) {
    if (depth < this->depth)
      return;
    // else
    this->depth = depth;
  }

	void run(graph_engine &graph) {
		if (is_active()) {
			directed_vertex_request req(get_id(), edge_type::BOTH_EDGES);
			request_partial_vertices(&req, 1);
		}
	}

	void run(graph_engine &graph, const page_vertex &vertex);

	void run_on_message(graph_engine &graph, const vertex_message &msg);
};

// This makes sure I don't participate in the same bfs twice  
bool diam_vertex::visit(const std::bitset<NUM_PARR_SEARCH> &msg_bits) const {
  // Avoid some redundant messages
  if ( msg_bits == this->bfs_bits ) 
    return false;
  // TODO: More bit ops to reduce redundancy 
  else {
    diam_vertex* pthis = const_cast<diam_vertex*>(this);

    // std::bitset<NUM_PARR_SEARCH>* pbfs_bits = const_cast<std::bitset<NUM_PARR_SEARCH>* >(&(this->bfs_bits));
    //*pbfs_bits |= msg_bits;
    pthis->bfs_bits |= msg_bits;
    std::bitset<NUM_PARR_SEARCH> ones;
    ones.set();
    if ((bfs_bits & ones) == ones) { // I've received all the msgs I can
      pthis->deactivate();
    }
    return true;
  }
}

void diam_vertex::run(graph_engine &graph, const page_vertex &vertex) { 
  if (!is_active())
    return;
  // multicast my bfs_bits
  
  // Only multicast if I have something to send
  const std::bitset<NUM_PARR_SEARCH> zeros;
  if ((this->bfs_bits | zeros) != this->bfs_bits) {
    int num_dests = vertex.get_num_edges(BOTH_EDGES);
    if (num_dests == 0)
      return;
    edge_seq_iterator it = vertex.get_neigh_seq_it(BOTH_EDGES, 0, num_dests);
    depth_message msg(this->bfs_bits, this->depth);
    graph.multicast_msg(it, msg);
  }
}

void diam_vertex::run_on_message(graph_engine &graph, const vertex_message &msg) {
  if (!is_active())
    return;
  // If visit is true I haven't already received this bfs depth msg
  depth_message &dmsg = (depth_message&) msg;
  if (visit(dmsg.get_bits())) {
    // This means I haven't received this msg before so broadcast to my neighs this dfs
    update_depth(dmsg.get_depth()); // try to update max depth
  }
}

void int_handler(int sig_num) {
	if (!graph_conf.get_prof_file().empty())
		ProfilerStop();
	exit(0);
}

void print_usage() {
	fprintf(stderr,
			"graph-diameter [options] conf_file graph_file index_file [num_parralel_searches]\n");
	fprintf(stderr, "-c confs: add more configurations to the system\n");
	fprintf(stderr, "-p: preload the graph\n");
	graph_conf.print_help();
	params.print_help();
}

int main(int argc, char *argv[]) {
	int opt;
	std::string confs;
	int num_opts = 0;
	bool preload = false;
	while ((opt = getopt(argc, argv, "c:p")) != -1) {
		num_opts++;
		switch (opt) {
			case 'c':
				confs = optarg;
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
	graph_conf.init(configs);
	graph_conf.print();

	signal(SIGINT, int_handler);
	init_io_system(configs);

	graph_index *index = NUMA_graph_index<diam_vertex>::create(index_file,
			graph_conf.get_num_threads(), params.get_num_nodes());

	graph_engine *graph = graph_engine::create(graph_conf.get_num_threads(),
			params.get_num_nodes(), graph_file, index);
	if (preload)
		graph->preload_graph();
  
	if (argc > 3) {
    uint32_t* ptr = const_cast<uint32_t*>(&NUM_PARR_SEARCH);
     *ptr = (uint32_t) atoi(argv[3]);
     printf("The value of NUM_PARR_SEARCH is %d", NUM_PARR_SEARCH);
    if (NUM_PARR_SEARCH > graph->get_max_vertex_id()) {
      fprintf(stderr, 
          "Argument: 'num_parralel_searches' cannot be > than the number of vertices\n");
    }
	}

	printf("Diameter estimator starts\n");
	printf("prof_file: %s\n", graph_conf.get_prof_file().c_str());
	if (!graph_conf.get_prof_file().empty())
		ProfilerStart(graph_conf.get_prof_file().c_str());

	struct timeval start, end;
	gettimeofday(&start, NULL);
  
  // Choose the random vertices
  while (START_VERTICES.size() != NUM_PARR_SEARCH) {
    vertex_id_t id = (vertex_id_t) (rand() % graph->get_max_vertex_id());
		if (graph->get_vertex_edges(id) > 0)
      START_VERTICES.push_back(id);

    printf("Choosing vertex: %d to randomly start\n", START_VERTICES.back());
  }  

  //  Sort em for binary search
  std::sort(START_VERTICES.begin(), START_VERTICES.end());
   
  printf("Start vertices ...\n");
  print_arr(&START_VERTICES[0], NUM_PARR_SEARCH);
  exit(-1); 

	graph->start(&START_VERTICES[0], NUM_PARR_SEARCH);
	graph->wait4complete();
	gettimeofday(&end, NULL);

	NUMA_graph_index<diam_vertex>::const_iterator it
		= ((NUMA_graph_index<diam_vertex> *) index)->begin();
	NUMA_graph_index<diam_vertex>::const_iterator end_it
		= ((NUMA_graph_index<diam_vertex> *) index)->end();

	uint32_t est_diam = 0;
	for (; it != end_it; ++it) {
		const diam_vertex &v = (const diam_vertex &) *it;
		if (v.get_depth() > est_diam)
      est_diam = v.get_depth();
	}

	if (!graph_conf.get_prof_file().empty())
		ProfilerStop();
	if (graph_conf.get_print_io_stat())
		print_io_thread_stat();
	graph_engine::destroy(graph);
	destroy_io_system();
	printf("Estimated graph diameter is %d. It takes %f seconds\n",
			est_diam, time_diff(start, end)); // FIXME: stub
}
