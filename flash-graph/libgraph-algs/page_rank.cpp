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
#ifdef PROFILER
#include <google/profiler.h>
#endif

#include <limits>

#include "thread.h"
#include "io_interface.h"
#include "container.h"
#include "concurrency.h"

#include "vertex_index.h"
#include "graph_engine.h"
#include "graph_config.h"
#include "FGlib.h"

namespace {

float DAMPING_FACTOR = 0.85;
float TOLERANCE = 1.0E-2; 
int num_iters = INT_MAX;

/*
 * pgrank_vertex needs to be initialized first.
 * Therefore, it has two stages.
 * pgrank_vertex2 doesn't need this process.
 */
enum pr_stage_t
{
	INIT,
	RUN,
};
pr_stage_t pr_stage;

class pgrank_vertex: public compute_directed_vertex
{
  float curr_itr_pr; // Current iteration's page rank
  vsize_t num_out_edges;

public:
  pgrank_vertex(vertex_id_t id): compute_directed_vertex(id) {
    this->curr_itr_pr = 1 - DAMPING_FACTOR; // Must be this
  }

  vsize_t get_num_out_edges() const {
	  return num_out_edges;
  }

  float get_curr_itr_pr() const{
    return curr_itr_pr;
  }

  float get_result() const {
	  return get_curr_itr_pr();
  }

  void run(vertex_program &prog);

	void run(vertex_program &prog, const page_vertex &vertex);

	void run_on_message(vertex_program &,
/* Only serves to activate on the next iteration */
			const vertex_message &msg) { }; 

	void run_on_vertex_header(vertex_program &, const vertex_header &header) {
		assert(get_id() == header.get_id());
		directed_vertex_header &dheader = (directed_vertex_header &) header;
		this->num_out_edges = dheader.get_num_out_edges();
	}
};

void pgrank_vertex::run(vertex_program &prog)
{
	vertex_id_t id = get_id();
	if (pr_stage == pr_stage_t::INIT) {
		request_vertex_headers(&id, 1);
	}
	else if (pr_stage == pr_stage_t::RUN) {
		// We perform pagerank for at most `num_iters' iterations.
		if (prog.get_graph().get_curr_level() >= num_iters)
			return;
		request_vertices(&id, 1); // put my edgelist in page cache
	}
};

void pgrank_vertex::run(vertex_program &prog, const page_vertex &vertex) {
  // Gather
  float accum = 0;
  page_byte_array::const_iterator<vertex_id_t> end_it
    = vertex.get_neigh_end(IN_EDGE);
  
  for (page_byte_array::const_iterator<vertex_id_t> it
      = vertex.get_neigh_begin(IN_EDGE); it != end_it; ++it) {
    vertex_id_t id = *it;
    pgrank_vertex& v = (pgrank_vertex&) prog.get_graph().get_vertex(id);
    // Notice I want this iteration's pagerank
    accum += (v.get_curr_itr_pr()/v.get_num_out_edges()); 
  }   

  // Apply
  float last_change = 0;
  if (vertex.get_num_edges(IN_EDGE) > 0) {
    float new_pr = ((1 - DAMPING_FACTOR)) + (DAMPING_FACTOR*(accum));
    last_change = new_pr - curr_itr_pr;
    curr_itr_pr = new_pr;
  }   
  
  // Scatter (activate your out-neighbors ... if you have any :) 
  if ( std::fabs( last_change ) > TOLERANCE ) {
	int num_dests = vertex.get_num_edges(OUT_EDGE);
    if (num_dests > 0) {
		edge_seq_iterator it = vertex.get_neigh_seq_it(OUT_EDGE, 0, num_dests);
		prog.activate_vertices(it) ;
    }   
  }
}

class pr_message: public vertex_message
{
	float delta;
public:
	pr_message(float delta): vertex_message(sizeof(pr_message),
			true) {
		this->delta = delta;
	}

	float get_delta() const {
		return delta;
	}
};

class pgrank_vertex2: public compute_directed_vertex
{
	float new_pr;
	float curr_itr_pr; // Current iteration's page rank
public:
	pgrank_vertex2(vertex_id_t id): compute_directed_vertex(id) {
		this->curr_itr_pr = 1 - DAMPING_FACTOR; // Must be this
		this->new_pr = curr_itr_pr;
	}

	float get_result() const{
		return new_pr;
	}

	void run(vertex_program &prog) { 
		// We perform pagerank for at most `num_iters' iterations.
		if (prog.get_graph().get_curr_level() >= num_iters)
			return;
		directed_vertex_request req(get_id(), edge_type::OUT_EDGE);
		request_partial_vertices(&req, 1);
	};

	void run(vertex_program &, const page_vertex &vertex);

	void run_on_message(vertex_program &, const vertex_message &msg1) {
		const pr_message &msg = (const pr_message &) msg1;
		new_pr += msg.get_delta();
	}
};

void pgrank_vertex2::run(vertex_program &prog, const page_vertex &vertex)
{
	int num_dests = vertex.get_num_edges(OUT_EDGE);
	edge_seq_iterator it = vertex.get_neigh_seq_it(OUT_EDGE, 0, num_dests);

	// If this is the first iteration.
	if (prog.get_graph().get_curr_level() == 0) {
		pr_message msg(curr_itr_pr / num_dests * DAMPING_FACTOR);
		prog.multicast_msg(it, msg);
	}
	else if (std::fabs(new_pr - curr_itr_pr) > TOLERANCE) {
		pr_message msg((new_pr - curr_itr_pr) / num_dests * DAMPING_FACTOR);
		prog.multicast_msg(it, msg);
		curr_itr_pr = new_pr;
	}
}

}

#include "save_result.h"
FG_vector<float>::ptr compute_pagerank(FG_graph::ptr fg, int num_iters,
		float damping_factor)
{
	DAMPING_FACTOR = damping_factor;
	if (DAMPING_FACTOR < 0 || DAMPING_FACTOR > 1) {
		fprintf(stderr, "Damping factor must be between 0 and 1 inclusive\n");
		exit(-1);
	}

	graph_index::ptr index = NUMA_graph_index<pgrank_vertex>::create(
			fg->get_index_file());
	graph_engine::ptr graph = graph_engine::create(fg->get_graph_file(),
			index, fg->get_configs());
	printf("Pagerank (at maximal %d iterations) starting\n", num_iters);
	printf("prof_file: %s\n", graph_conf.get_prof_file().c_str());
#ifdef PROFILER
	if (!graph_conf.get_prof_file().empty())
		ProfilerStart(graph_conf.get_prof_file().c_str());
#endif

	struct timeval start, end;
	gettimeofday(&start, NULL);
	pr_stage = pr_stage_t::INIT;
	graph->start_all(); 
	graph->wait4complete();
	pr_stage = pr_stage_t::RUN;
	graph->start_all(); 
	graph->wait4complete();
	gettimeofday(&end, NULL);

	FG_vector<float>::ptr ret = FG_vector<float>::create(
			graph->get_num_vertices());
	graph->query_on_all(vertex_query::ptr(
				new save_query<float, pgrank_vertex>(ret)));

#ifdef PROFILER
	if (!graph_conf.get_prof_file().empty())
		ProfilerStop();
#endif

	printf("It takes %f seconds\n", time_diff(start, end));
	return ret;
}

FG_vector<float>::ptr compute_pagerank2(FG_graph::ptr fg, int num_iters,
		float damping_factor)
{
	DAMPING_FACTOR = damping_factor;
	if (DAMPING_FACTOR < 0 || DAMPING_FACTOR > 1) {
		fprintf(stderr, "Damping factor must be between 0 and 1 inclusive\n");
		exit(-1);
	}

	graph_index::ptr index = NUMA_graph_index<pgrank_vertex2>::create(
			fg->get_index_file());
	graph_engine::ptr graph = graph_engine::create(fg->get_graph_file(),
			index, fg->get_configs());
	printf("Pagerank (at maximal %d iterations) starting\n", num_iters);
	printf("prof_file: %s\n", graph_conf.get_prof_file().c_str());
#ifdef PROFILER
	if (!graph_conf.get_prof_file().empty())
		ProfilerStart(graph_conf.get_prof_file().c_str());
#endif

	struct timeval start, end;
	gettimeofday(&start, NULL);
	graph->start_all(); 
	graph->wait4complete();
	gettimeofday(&end, NULL);

	FG_vector<float>::ptr ret = FG_vector<float>::create(
			graph->get_num_vertices());
	graph->query_on_all(vertex_query::ptr(
				new save_query<float, pgrank_vertex2>(ret)));

#ifdef PROFILER
	if (!graph_conf.get_prof_file().empty())
		ProfilerStop();
#endif

	printf("It takes %f seconds\n", time_diff(start, end));
	return ret;
}
