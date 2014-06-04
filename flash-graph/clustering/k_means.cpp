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
#include <limits> 

#include "thread.h"
#include "io_interface.h"
#include "container.h"
#include "concurrency.h"

#include "vertex_index.h"
#include "graph_engine.h"
#include "graph_config.h"
#include "math.h"

int g_MAX_ITERS = 100;
// float g_TERM_FRAC_SWITCHES = .05; // If I have < than this then end the k-means routine

class km_vector
{
  std::vector<float> data;

  public:
    km_vector() { };

    km_vector(std::vector<float>& data) { 
      this->data = data;
    };

    km_vector(vsize_t n) {
      data.resize(n);
    }

    void resize(vsize_t n) {
      data.resize(n);
    }

    km_vector& operator+=(km_vector& other) {
      assert(this->data.size() == other.data.size());
      // #pragma omp parallel for
      for (vsize_t i=0; i < this->data.size(); i++) {
        this->data[i]+= other.data[i];
      }
      return *this;
    }

    km_vector& operator=(km_vector& other) {
     this->data = other.data; 
     return *this;
    }

    float L2Norm(km_vector& other_v) {
      // Assumption is a dense vector
      float sum = 0;
      for (vsize_t i=0; i< other_v.size(); i++) {
        float diff = get_elem(i) - other_v.get_elem(i);
        sum += diff * diff;
      }
      return sqrt(sum);
    } 

    float get_elem(vsize_t i) {
      return this->data[i];
    }

    void set_elem(vsize_t idx, float val) {
      this->data[idx] = val;
    }

    vsize_t size() {
      return this->data.size();
    }

    void push_back(float val) {
      this->data.push_back(val);
    }

    bool empty() {
      return this->data.empty();
    }

    void set_data(std::vector<float> v) {
      this->data = v;
    }

    std::vector<float>& get_data() {
      return this->data;
    }

    void assign(vsize_t size, int val) {
      this->data.assign(size, val);
    }

    void print() {
      std::cout << "[";
      for (vsize_t i = 0; i < size(); i++) {
        std::cout << " " << this->data[i];
      }
      std::cout <<  " ]\n";
    }
};

enum phase_t {
  e_step, /* Update cluster_id for each vertex */
  m_step, /* Update cluster_means for each cluster */
};

// Parallel functions
#if 0
km_vector& compute_mean(km_vector& kmv, float num_members) {
  
      // #pragma omp parallel for firstprivate(num_members, kmv)
      for (vsize_t i=0; i < kmv.size(); i++) {
        kmv[i] /= num_members;
      }
      return kmv;
}
#endif

class cluster 
{
  // Cluster_id implicit in its position in the vector<cluster>

  km_vector curr_cluster_mean;
  km_vector next_cluster_mean; // Cluster mean for the next round of K-means
  bool changed; 
  km_vector member_sum;
  vsize_t num_members;

  public:

  cluster() {
    changed = false;
  }

  void set_curr_cluster_mean(std::vector<float>& v) {
    curr_cluster_mean.set_data(v);
  }
  
  void add_member(km_vector& vec) {
    if (next_cluster_mean.empty()) { next_cluster_mean.resize(vec.size()); }
    if (!changed) { 
      changed = true;
      num_members = 0;
    }

    next_cluster_mean += vec; // TODO: Locking for no races
    num_members++;
  }

  void update_cluster_mean() {
    if (changed) {
    //  curr_cluster_mean = compute_mean(next_cluster_mean, (float)num_members);

      for (vsize_t i=0; i < next_cluster_mean.size(); i++) {
        next_cluster_mean.set_elem(i, (next_cluster_mean.get_elem(i)/(float)num_members));
      }
      curr_cluster_mean = next_cluster_mean;
    }
    next_cluster_mean.assign(next_cluster_mean.size(), 0); 
    changed = false; // reset the changed val
  }

  km_vector& get_mean() {
    return this->curr_cluster_mean;
  }
};

/* VarDecl */ 
std::vector<cluster> g_clusters;

class kmeans_vertex: public compute_vertex
{
  km_vector comp_vector; // comparison vector
  vsize_t cluster_id; // TODO: Use this later to reduce computation
  phase_t phase;

public:
	kmeans_vertex() {
    cluster_id = -1; // Every vertex is in an invalid cluster initially
    phase = e_step;
	}

  kmeans_vertex(vertex_id_t id, const vertex_index &index1):
    compute_vertex(id, index1) {
  }

  void set_comp_vector(km_vector& kmv) {
    this->comp_vector = kmv;
  }

  
  vsize_t compute_cluster() {
    float min_diff = std::numeric_limits<float>::max();
    int new_cluster = INT_MAX;

    for (uint32_t i = 0; i < g_clusters.size(); i++) {
      float diff = comp_vector.L2Norm(g_clusters[i].get_mean());
      if (diff < min_diff ) { 
        min_diff = diff; 
        new_cluster = i;
      }
    }
    return new_cluster;
  }

  void run(vertex_program &prog);

	void run(vertex_program &prog, const page_vertex &vertex) { };

	void run_on_message(vertex_program &prog, const vertex_message &msg) { };
};

void kmeans_vertex::run(vertex_program &prog) {
  switch (phase) {
    case e_step:
      if (get_id() == 0) {
        printf("Activating vertex 0 for the m_step ...\n");
        vertex_id_t id = get_id();
        prog.activate_vertices(&id, 1); // Only one is active next itr
      }

      cluster_id = compute_cluster();
      g_clusters[cluster_id].add_member(comp_vector);
      break;
    case m_step: // This will only happen to only vertex with id == 0
      printf("Activating all other vertices for the next iteration ...\n");
      
      // Activate all  
      vertex_id_t min_id = prog.get_graph().get_max_vertex_id();
      const static vertex_id_t range = prog.get_graph().get_max_vertex_id() -
                                prog.get_graph().get_min_vertex_id() + 1;
       
      vertex_id_t ids [range];

      // #pragma omp parallel for firstprivate(min_id, range)
      for (vertex_id_t idx = 0; idx <= range; idx++) {
        ids[idx] = idx + min_id;
      }

      if ( prog.get_graph().get_curr_level() <= g_MAX_ITERS ) { // TODO: set g_NUM_ITERS = 2*graph iterations
        prog.activate_vertices(&ids[0], prog.get_graph().get_num_vertices()); // Activate all vertices
      }
      break;
  }
}

// We need to iniate this with a comp vector derived from Eigs computation
class kmeans_initiator: public vertex_initiator
{
public:
  std::vector<std::vector<float> > comp_vectors;

  kmeans_initiator(std::vector<std::vector<float> >& comp_vectors) {
    this->comp_vectors = comp_vectors;
  }

	void init(compute_vertex &v) {
		kmeans_vertex &kmv = (kmeans_vertex &) v;
    km_vector comp_v_kmv(comp_vectors[v.get_id()]);
    kmv.set_comp_vector(comp_v_kmv);
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
			"k-means [options] conf_file graph_file index_file max_iters (default=50) \n");
	fprintf(stderr, "-c confs: add more configurations to the system\n");
	graph_conf.print_help();
	params.print_help();
}

// Helpers
// FIXME: Finish me!
std::vector<float> assemble_clusters(graph_engine::ptr graph) {

#if 0
  assigned_clusters(0, num_elem);

#pragma parallel for private(assigned_clusters)
  for (int i=0; i <= num_elem; i++) {
    assigned_clusters[i] = 
  }
#endif

  return std::vector<float>(); // stub 
}

void print_func(vertex_id_t i) {
  std::cout << " " << i;
}
void print_vector(std::vector<float> v) {
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
	g_MAX_ITERS = atoi(argv[3]); // Set kmin
  
	config_map configs(conf_file);
	configs.add_options(confs);

	signal(SIGINT, int_handler);

	graph_index::ptr index = NUMA_graph_index<kmeans_vertex>::create(index_file);
	graph_engine::ptr graph = graph_engine::create(graph_file, index, configs);
	printf("K-means starting\n");
	printf("prof_file: %s\n", graph_conf.get_prof_file().c_str());
	if (!graph_conf.get_prof_file().empty())
		ProfilerStart(graph_conf.get_prof_file().c_str());

  // TODO: Run Eigensolver return a vector of eigenvectors
  std::vector<std::vector<float> > eigs;

  struct timeval start, end;
  gettimeofday(&start, NULL);

  // TODO: Fix vertex initializer
#if 0
  graph->start_all(vertex_initiator::ptr(eigs), 
      vertex_program_creater::ptr()); // TODO: Fix this mess
#endif

  graph->wait4complete();
  gettimeofday(&end, NULL);

  std::vector<float> assigned_clusters = assemble_clusters(graph);
  printf("Printing clusters: \n");
  print_vector(assigned_clusters);
  
	if (!graph_conf.get_prof_file().empty())
		ProfilerStop();
}
