/**
 * Copyright 2014 Open Connectome Project (http://openconnecto.me)
 * Written by Da Zheng (zhengda1936@gmail.com)
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

#include <set>
#include <vector>

#include "thread.h"
#include "io_interface.h"

#include "graph_engine.h"
#include "graph_config.h"

#include "FG_vector.h"
#include "FGlib.h"

#include "scan_graph.h"

namespace {

scan_stage_t scan_stage;

struct timeval graph_start;

class global_max
{
	volatile size_t value;
	pthread_spinlock_t lock;
public:
	global_max() {
		value = 0;
		pthread_spin_init(&lock, PTHREAD_PROCESS_PRIVATE);
	}

	global_max(size_t init) {
		value = init;
		pthread_spin_init(&lock, PTHREAD_PROCESS_PRIVATE);
	}

	bool update(size_t new_v) {
		if (new_v <= value)
			return false;

		bool ret = false;
		pthread_spin_lock(&lock);
		if (new_v > value) {
			value = new_v;
			ret = true;
		}
		pthread_spin_unlock(&lock);
		return ret;
	}

	size_t get() const {
		return value;
	}
} max_scan;

typedef std::pair<vertex_id_t, size_t> vertex_scan;

/**
 * This class maintains the local scan that have been computed.
 */
class scan_collection
{
	bool sorted;
	std::vector<vertex_scan> scans;
	pthread_spinlock_t lock;

	class greater {
	public:
		bool operator()(const vertex_scan &s1, const vertex_scan &s2) {
			return s1.second > s2.second;
		}
	};
public:
	scan_collection() {
		pthread_spin_init(&lock, PTHREAD_PROCESS_PRIVATE);
		sorted = false;
	}

	/**
	 * Get the ith largest scan.
	 */
	vertex_scan get(int idx) {
		pthread_spin_lock(&lock);
		if (!sorted) {
			// It needs to stored in the descending order on scan.
			std::sort(scans.begin(), scans.end(), greater());
			sorted = true;
		}
		vertex_scan ret = scans[idx];
		pthread_spin_unlock(&lock);
		return ret;
	}

	void add(vertex_id_t id, size_t scan) {
		pthread_spin_lock(&lock);
		sorted = false;
		scans.push_back(vertex_scan(id, scan));
		pthread_spin_unlock(&lock);
	}

	size_t get_size() {
		pthread_spin_lock(&lock);
		size_t ret = scans.size();
		pthread_spin_unlock(&lock);
		return ret;
	}
} known_scans;

class topK_scan_vertex: public compute_vertex
{
protected:
	vsize_t degree;
	multi_func_value local_value;
public:
	topK_scan_vertex(vertex_id_t id): compute_vertex(id) {
		degree = 0;
	}

	vsize_t get_degree() const {
		return degree;
	}

	bool has_local_scan() const {
		return local_value.has_real_local();
	}

	size_t get_local_scan() const {
		return local_value.get_real_local();
	}

	bool has_est_local() const {
		return local_value.has_est_local();
	}

	size_t get_est_local_scan() const {
		return local_value.get_est_local();
	}

	void run(vertex_program &prog);

	void run(vertex_program &prog, const page_vertex &vertex) {
		assert(scan_stage == scan_stage_t::RUN);
		if (vertex.get_id() == get_id())
			run_on_itself(prog, vertex);
		else
			run_on_neighbor(prog, vertex);
	}

	void run_on_itself(vertex_program &prog, const page_vertex &vertex);
	void run_on_neighbor(vertex_program &prog, const page_vertex &vertex);

	void run_on_message(vertex_program &prog, const vertex_message &msg) {
	}

	void run_on_vertex_header(vertex_program &prog, const vertex_header &header) {
		assert(get_id() == header.get_id());
		degree = header.get_num_edges();
	}

	void finding_triangles_end(vertex_program &prog, runtime_data_t *data) {
		if (max_scan.update(data->local_scan)) {
			struct timeval curr;
			gettimeofday(&curr, NULL);
			printf("%d: new max scan: %ld at v%u\n",
					(int) time_diff(graph_start, curr),
					data->local_scan, get_id());
		}
		known_scans.add(get_id(), data->local_scan);
	}
};

runtime_data_t *create_runtime(graph_engine &graph, topK_scan_vertex &scan_v,
		const page_vertex &pg_v)
{
	class skip_self {
		vertex_id_t id;
		graph_engine &graph;
	public:
		skip_self(graph_engine &_graph, vertex_id_t id): graph(_graph) {
			this->id = id;
		}

		bool operator()(attributed_neighbor &e) {
			return operator()(e.get_id());
		}

		bool operator()(vertex_id_t id) {
			return this->id == id;
		}
	};

	class merge_edge
	{
	public:
		attributed_neighbor operator()(const attributed_neighbor &e1,
				const attributed_neighbor &e2) {
			assert(e1.get_id() == e2.get_id());
			return attributed_neighbor(e1.get_id(),
					e1.get_num_dups() + e2.get_num_dups());
		}
	};

	merge_edge merge;
	std::vector<attributed_neighbor> neighbors(
			pg_v.get_num_edges(edge_type::BOTH_EDGES));
	size_t num_neighbors = unique_merge(
			pg_v.get_neigh_begin(edge_type::IN_EDGE),
			pg_v.get_neigh_end(edge_type::IN_EDGE),
			pg_v.get_neigh_begin(edge_type::OUT_EDGE),
			pg_v.get_neigh_end(edge_type::OUT_EDGE),
			skip_self(graph, pg_v.get_id()), merge,
			neighbors.begin());
	neighbors.resize(num_neighbors);
	return new runtime_data_t(std::unique_ptr<neighbor_list>(
				new neighbor_list(pg_v, neighbors)));
}

void destroy_runtime(runtime_data_t *data)
{
	delete data;
}

size_t get_est_local_scan(graph_engine &graph, topK_scan_vertex &scan_v,
		const page_vertex &page_v)
{
	// We have estimated the local scan of this vertex, return
	// the estimated one
	if (scan_v.has_est_local())
		return scan_v.get_est_local_scan();

	class skip_self {
		vertex_id_t id;
	public:
		skip_self(vertex_id_t id) {
			this->id = id;
		}

		bool operator()(vertex_id_t id) {
			return this->id == id;
		}
	};

	class merge_edge {
	public:
		vertex_id_t operator()(vertex_id_t e1, vertex_id_t e2) {
			assert(e1 == e2);
			return e1;
		}
	};

	std::vector<vertex_id_t> all_neighbors(
			page_v.get_num_edges(edge_type::BOTH_EDGES));
	size_t num_neighbors = unique_merge(
			page_v.get_neigh_begin(edge_type::IN_EDGE),
			page_v.get_neigh_end(edge_type::IN_EDGE),
			page_v.get_neigh_begin(edge_type::OUT_EDGE),
			page_v.get_neigh_end(edge_type::OUT_EDGE),
			skip_self(page_v.get_id()), merge_edge(),
			all_neighbors.begin());
	all_neighbors.resize(num_neighbors);

	size_t tot_edges = scan_v.get_degree();
	for (size_t i = 0; i < all_neighbors.size(); i++) {
		topK_scan_vertex &v = (topK_scan_vertex &) graph.get_vertex(all_neighbors[i]);
		// The max number of common neighbors should be smaller than all neighbors
		// in the neighborhood, assuming there aren't duplicated edges.
		tot_edges += min(v.get_degree(), num_neighbors * 2);
	}
	tot_edges /= 2;
	return tot_edges;
}

void topK_scan_vertex::run_on_itself(vertex_program &prog, const page_vertex &vertex)
{
	size_t num_local_edges = vertex.get_num_edges(edge_type::BOTH_EDGES);
	if (num_local_edges == 0)
		return;

	size_t tot_edges = ::get_est_local_scan(prog.get_graph(), *this, vertex);
	local_value.set_est_local(tot_edges);
	if (tot_edges < max_scan.get())
		return;

	assert(!local_value.has_runtime_data());

	runtime_data_t *local_data = create_runtime(prog.get_graph(), *this, vertex);
	local_value.set_runtime_data(local_data);

	page_byte_array::const_iterator<vertex_id_t> it = vertex.get_neigh_begin(
			edge_type::BOTH_EDGES);
	page_byte_array::const_iterator<vertex_id_t> end = vertex.get_neigh_end(
			edge_type::BOTH_EDGES);
	size_t tmp = 0;
	for (; it != end; ++it) {
		vertex_id_t id = *it;
		// Ignore loops
		if (id != vertex.get_id())
			tmp++;
	}
	local_data->local_scan += tmp;

	if (local_data->neighbors->empty()) {
		local_value.set_real_local(local_data->local_scan);
		finding_triangles_end(prog, local_data);
		destroy_runtime(local_data);
		return;
	}

	std::vector<vertex_id_t> neighbors;
	local_data->neighbors->get_neighbors(neighbors);
	request_vertices(neighbors.data(), neighbors.size());
}

void topK_scan_vertex::run_on_neighbor(vertex_program &prog, const page_vertex &vertex)
{
	assert(local_value.has_runtime_data());
	runtime_data_t *local_data = local_value.get_runtime_data();
	local_data->num_joined++;
	size_t ret = local_data->neighbors->count_edges(&vertex);
	if (ret > 0)
		local_data->local_scan += ret;

	// If we have seen all required neighbors, we have complete
	// the computation. We can release the memory now.
	if (local_data->num_joined == local_data->neighbors->size()) {
		local_value.set_real_local(local_data->local_scan);

		finding_triangles_end(prog, local_data);
		destroy_runtime(local_data);
	}
}

void topK_scan_vertex::run(vertex_program &prog)
{
	vertex_id_t id = get_id();
	if (scan_stage == scan_stage_t::INIT) {
		request_vertex_headers(&id, 1);
	}
	else if (scan_stage == scan_stage_t::RUN) {
		bool req_itself = false;
		// If we have computed local scan on the vertex, skip the vertex.
		if (has_local_scan())
			return;
		// If we have estimated the local scan, we should use the estimated one.
		else if (has_est_local())
			req_itself = get_est_local_scan() > max_scan.get();
		else {
			// If this is the first time to compute on the vertex, we can still
			// skip a lot of vertices with this condition.
			size_t num_local_edges = get_degree();
			req_itself = num_local_edges * num_local_edges >= max_scan.get();
		}
		if (req_itself)
			request_vertices(&id, 1);
	}
}

class vertex_size_scheduler: public vertex_scheduler
{
	graph_engine::ptr graph;
public:
	vertex_size_scheduler(graph_engine::ptr graph) {
		this->graph = graph;
	}
	void schedule(std::vector<vertex_id_t> &vertices);
};

void vertex_size_scheduler::schedule(std::vector<vertex_id_t> &ids)
{
	std::vector<compute_vertex *> vertices(ids.size());
	for (size_t i = 0; i < ids.size(); i++)
		vertices[i] = &graph->get_vertex(ids[i]);
	class comp_size
	{
	public:
		bool operator()(const compute_vertex *v1, const compute_vertex *v2) {
			return ((const topK_scan_vertex *) v1)->get_degree()
				> ((const topK_scan_vertex *) v2)->get_degree();
		}
	};

	std::sort(vertices.begin(), vertices.end(), comp_size());
	for (size_t i = 0; i < ids.size(); i++)
		ids[i] = vertices[i]->get_id();
}

}

FG_vector<std::pair<vertex_id_t, size_t> >::ptr compute_topK_scan(
		FG_graph::ptr fg, size_t topK)
{
	graph_index::ptr index = NUMA_graph_index<topK_scan_vertex>::create(
			fg->get_index_file());
	graph_engine::ptr graph = graph_engine::create(fg->get_graph_file(),
			index, fg->get_configs());

	printf("scan statistics starts\n");
	printf("prof_file: %s\n", graph_conf.get_prof_file().c_str());
#ifdef PROFILER
	if (!graph_conf.get_prof_file().empty())
		ProfilerStart(graph_conf.get_prof_file().c_str());
#endif

	scan_stage = scan_stage_t::INIT;
	graph->start_all();
	graph->wait4complete();
	scan_stage = scan_stage_t::RUN;

	// Let's schedule the order of processing activated vertices according
	// to the size of vertices. We start with processing vertices with higher
	// degrees in the hope we can find the max scan as early as possible,
	// so that we can simple ignore the rest of vertices.
	graph->set_vertex_scheduler(vertex_scheduler::ptr(
				new vertex_size_scheduler(graph)));
	graph->set_max_processing_vertices(3);

	class remove_small_filter: public vertex_filter
	{
		size_t min;
	public:
		remove_small_filter(size_t min) {
			this->min = min;
		}

		bool keep(compute_vertex &v) {
			topK_scan_vertex &scan_v = (topK_scan_vertex &) v;
			return scan_v.get_degree() >= min;
		}
	};

	size_t min_edges = 1000;
	std::shared_ptr<vertex_filter> filter
		= std::shared_ptr<vertex_filter>(new remove_small_filter(min_edges));
	struct timeval start, end;
	gettimeofday(&start, NULL);
	graph_start = start;
	printf("Computing local scan on at least %ld vertices\n", topK);
	while (known_scans.get_size() < topK) {
		gettimeofday(&start, NULL);
		graph->start(filter);
		graph->wait4complete();
		gettimeofday(&end, NULL);
		printf("It takes %f seconds\n", time_diff(start, end));
		printf("global max scan: %ld\n", max_scan.get());
		max_scan = global_max(0);
	}

	class remove_small_scan_filter: public vertex_filter
	{
		size_t min;
	public:
		remove_small_scan_filter(size_t min) {
			this->min = min;
		}

		bool keep(compute_vertex &v) {
			topK_scan_vertex &scan_v = (topK_scan_vertex &) v;
			size_t num_local_edges = scan_v.get_degree();
			return num_local_edges * num_local_edges >= min;
		}
	};

	printf("Compute local scan on %ld vertices\n", known_scans.get_size());
	printf("Looking for top %ld local scan\n", topK);
	size_t prev_topK_scan;
	do {
		prev_topK_scan = known_scans.get(topK - 1).second;
		// Let's use the topK as the max scan for unknown vertices
		// and see if we can find a new vertex that has larger local scan.
		max_scan = global_max(prev_topK_scan);

		gettimeofday(&start, NULL);
		graph->start(std::shared_ptr<vertex_filter>(
					new remove_small_scan_filter(prev_topK_scan)));
		graph->wait4complete();
		gettimeofday(&end, NULL);
		printf("It takes %f seconds\n", time_diff(start, end));
		printf("global max scan: %ld\n", max_scan.get());
		// If the previous topK is different from the current one,
		// it means we have found new local scans that are larger
		// than the previous topK. We should use the new topK and
		// try again.
	} while (prev_topK_scan != known_scans.get(topK - 1).second);
	assert(known_scans.get_size() >= topK);

#ifdef PROFILER
	if (!graph_conf.get_prof_file().empty())
		ProfilerStop();
#endif
	printf("It takes %f seconds for top %ld\n", time_diff(graph_start, end),
			topK);

	FG_vector<std::pair<vertex_id_t, size_t> >::ptr vec
		= FG_vector<std::pair<vertex_id_t, size_t> >::create(topK);
	for (size_t i = 0; i < topK; i++)
		vec->set(i, known_scans.get(i));
	return vec;
}
