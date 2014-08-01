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
#include "FGlib.h"
#include "save_result.h"

vsize_t CURRENT_K; // Min degree necessary to be part of the k-core graph
vsize_t PREVIOUS_K; 
bool all_greater_than_core = true;

enum kcore_stage_t
{
	INIT_DEGREE,
	KCORE,
};
kcore_stage_t stage;

class kcore_vertex: public compute_vertex
{
	bool deleted;
	vsize_t core;
	vsize_t degree; 

	public:
	kcore_vertex(vertex_id_t id): compute_vertex(id) {
		this->deleted = false;
		this->core = -1;
		this->degree = 0;
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

	size_t get_result() const {
		return get_core();
	}

	void run(vertex_program &prog);

	void run(vertex_program &prog, const page_vertex &vertex);

	void run_on_message(vertex_program &prog, const vertex_message &msg); 

	void run_on_vertex_header(vertex_program &prog, const vertex_header &header) {
		degree = header.get_num_edges();
		this->core = degree == 0 ? 0 :
			degree == 1 ? 1 : -1; // Everyone between kmin < core > kmax will get this core
		this->deleted = degree == 0 ? true : false; // If your degree you're already deleted
	}
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
	if (stage == INIT_DEGREE) {
		vertex_id_t id = get_id();
		request_vertex_headers(&id, 1);
		return;
	}

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
		kcore_vertex& kv = (kcore_vertex&) v;
		if ((kv.get_degree()) > max_degree) {
			max_degree = kv.get_degree();
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

void set_kmax(graph_engine::ptr graph, size_t& kmax)
{
	printf("Computing kmax as max_degree ...\n");
	vertex_query::ptr mdq(new max_degree_query());
	graph->query_on_all(mdq); 
	kmax = ((max_degree_query *) mdq.get())->get_max_degree();
}

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

FG_vector<size_t>::ptr compute_kcore(FG_graph::ptr fg,
		size_t k, size_t kmax)
{
	graph_index::ptr index = NUMA_graph_index<kcore_vertex>::create(
			fg->get_index_file());
	graph_engine::ptr graph = graph_engine::create(fg->get_graph_file(),
			index, fg->get_configs());

	if (k > graph->get_max_vertex_id()) {
		fprintf(stderr, "'k' must be between 2 and the number of nodes in the graph\n");
		exit(-1);
	}

	struct timeval start, end;
	gettimeofday(&start, NULL);

	CURRENT_K = k;
	printf("Running the init degree stage ...\n");
	stage = INIT_DEGREE;
	graph->start_all(); 
	graph->wait4complete();
	stage = KCORE;

	if (kmax == 0 && k != 0) {
		set_kmax(graph, kmax);
	}
	printf("Setting kmax as %lu\n", kmax);


	for (; CURRENT_K <= kmax; CURRENT_K++) {

		std::shared_ptr<vertex_filter> filter
			= std::shared_ptr<vertex_filter>(new activate_k_filter(CURRENT_K));

		graph->start(filter, vertex_program_creater::ptr()); 
		graph->wait4complete();

		if (all_greater_than_core) { // There's a chance we can hop forward
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

			if (CURRENT_K > kmax) {
				printf("\nTerminating computation at kmax\n");
				break;
			}
		}
		all_greater_than_core = true;

#if 1
		vertex_query::ptr cvq(new count_vertex_query());
		graph->query_on_all(cvq);
		size_t in_k_core = ((count_vertex_query *) cvq.get())->get_num();
		printf("\n******************************************\n"
				"%d-core shows %ld vertices > %d degree\n"
				"\n******************************************\n",
				CURRENT_K, in_k_core, CURRENT_K);
#endif
		PREVIOUS_K = CURRENT_K;
	}

	gettimeofday(&end, NULL);
	printf("\nK-core took %f sec to complete\n", time_diff(start, end)); 

	FG_vector<size_t>::ptr ret = FG_vector<size_t>::create(
			graph->get_num_vertices());

	graph->query_on_all(vertex_query::ptr(
				new save_query<size_t, kcore_vertex>(ret)));

	return ret;
}

