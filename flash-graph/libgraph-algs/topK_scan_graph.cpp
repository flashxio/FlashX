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
#include <gperftools/profiler.h>
#endif

#include <set>
#include <vector>

#include "graph_engine.h"
#include "graph_config.h"

#include "FG_vector.h"
#include "FGlib.h"

#include "scan_graph.h"

class part_local_t
{
	int num_parts;
	size_t local_scan;
public:
	part_local_t(size_t part) {
		num_parts = 1;
		local_scan = part;
	}

	void add_part(size_t part) {
		num_parts++;
		this->local_scan += part;
	}

	size_t get_local_scan() const {
		return local_scan;
	}

	int get_num_parts() const {
		return num_parts;
	}
};

namespace {

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

class scan_msg: public vertex_message
{
	size_t num;
public:
	scan_msg(size_t num): vertex_message(sizeof(scan_msg), false) {
		this->num = num;
	}

	size_t get() const {
		return num;
	}
};

enum part_scan_type
{
	EST_LOCAL,
	NEIGH,
};

class part_scan_msg: public vertex_message
{
	part_scan_type type;
	size_t est_local;
public:
	part_scan_msg(part_scan_type type, size_t est_local): vertex_message(
			sizeof(part_scan_msg), false) {
		this->est_local = est_local;
		this->type = type;
	}

	size_t get() const {
		return est_local;
	}

	part_scan_type get_type() const {
		return type;
	}
};

class topK_scan_vertex: public compute_vertex
{
	multi_func_value local_value;
public:
	topK_scan_vertex(vertex_id_t id): compute_vertex(id) {
	}

	bool has_part_local_scan() const {
		return local_value.has_part_local();
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
		if (vertex.get_id() == prog.get_vertex_id(*this))
			run_on_itself(prog, vertex);
		else
			run_on_neighbor(prog, vertex);
	}

	void run_on_itself(vertex_program &prog, const page_vertex &vertex);
	void run_on_neighbor(vertex_program &prog, const page_vertex &vertex);
	void run_on_message(vertex_program &prog, const vertex_message &msg1);

	void finding_triangles_end(vertex_program &prog, size_t local_scan) {
		vertex_id_t id = prog.get_vertex_id(*this);
		if (max_scan.update(local_scan)) {
			BOOST_LOG_TRIVIAL(info)
				<< boost::format("new max scan: %1% at v%2%")
				% local_scan % id;
		}
		known_scans.add(id, local_scan);
	}
};

struct scan_runtime_data_t: public runtime_data_t
{
	vsize_t num_required;

	scan_runtime_data_t(std::shared_ptr<neighbor_list> neighbors): runtime_data_t(
			neighbors) {
		num_required = this->neighbors->size();
	}
};

scan_runtime_data_t *create_runtime(graph_engine &graph, const page_vertex &pg_v)
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
	return new scan_runtime_data_t(std::shared_ptr<neighbor_list>(
				new neighbor_list(pg_v, neighbors)));
}

void destroy_runtime(scan_runtime_data_t *data)
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

	size_t tot_edges = page_v.get_num_edges(edge_type::BOTH_EDGES);
	for (size_t i = 0; i < all_neighbors.size(); i++) {
		vsize_t degree = graph.get_num_edges(all_neighbors[i]);
		// The max number of common neighbors should be smaller than all neighbors
		// in the neighborhood, assuming there aren't duplicated edges.
		tot_edges += min(degree, num_neighbors * 2);
	}
	tot_edges /= 2;
	return tot_edges;
}

void topK_scan_vertex::run_on_message(vertex_program &prog,
		const vertex_message &msg1)
{
	const scan_msg &msg = (const scan_msg &) msg1;
	if (!local_value.has_part_local()) {
		local_value.set_part_local(new part_local_t(msg.get()));
		assert(local_value.has_part_local());
	}
	else {
		part_local_t *data = local_value.get_part_local();
		data->add_part(msg.get());
		// We send the local degree of the vertex and the edges
		// counted in each partition in separately messages.
		if (data->get_num_parts() == graph_conf.get_num_vparts() + 1) {
			size_t local_scan = data->get_local_scan();
			local_value.set_real_local(local_scan);
			delete data;
			finding_triangles_end(prog, local_scan);
		}
	}
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

	scan_runtime_data_t *local_data = create_runtime(prog.get_graph(), vertex);
	local_value.set_runtime_data(local_data);

	size_t tmp = 0;
	page_byte_array::seq_const_iterator<vertex_id_t> it = vertex.get_neigh_seq_it(
			edge_type::IN_EDGE);
	PAGE_FOREACH(vertex_id_t, id, it) {
		// Ignore loops
		if (id != vertex.get_id())
			tmp++;
	} PAGE_FOREACH_END
	it = vertex.get_neigh_seq_it(edge_type::OUT_EDGE);
	PAGE_FOREACH(vertex_id_t, id, it) {
		// Ignore loops
		if (id != vertex.get_id())
			tmp++;
	} PAGE_FOREACH_END
	local_data->local_scan += tmp;

	if (local_data->neighbors->empty()) {
		local_value.set_real_local(local_data->local_scan);
		finding_triangles_end(prog, local_data->local_scan);
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
	scan_runtime_data_t *local_data
		= (scan_runtime_data_t *) local_value.get_runtime_data();
	local_data->num_joined++;
	size_t ret = local_data->neighbors->count_edges(&vertex);
	if (ret > 0)
		local_data->local_scan += ret;

	// If we have seen all required neighbors, we have complete
	// the computation. We can release the memory now.
	if (local_data->num_joined == local_data->neighbors->size()) {
		local_value.set_real_local(local_data->local_scan);

		finding_triangles_end(prog, local_data->local_scan);
		destroy_runtime(local_data);
	}
}

void topK_scan_vertex::run(vertex_program &prog)
{
	vertex_id_t id = prog.get_vertex_id(*this);
	bool req_itself = false;
	assert(!has_part_local_scan());
	// If we have computed local scan on the vertex, skip the vertex.
	if (has_local_scan())
		return;
	// If we have estimated the local scan, we should use the estimated one.
	else if (has_est_local())
		req_itself = get_est_local_scan() > max_scan.get();
	else {
		// If this is the first time to compute on the vertex, we can still
		// skip a lot of vertices with this condition.
		size_t num_local_edges = prog.get_num_edges(id);
		req_itself = num_local_edges * num_local_edges >= max_scan.get();
	}
	if (req_itself)
		request_vertices(&id, 1);
}

class part_topK_scan_vertex: public part_compute_vertex
{
	size_t est_local;
	int64_t part_local;
	scan_runtime_data_t *local_data;
	pthread_spinlock_t lock;
public:
	part_topK_scan_vertex(vertex_id_t id, int part_id): part_compute_vertex(id,
			part_id) {
		est_local = 0;
		part_local = -1;
		local_data = NULL;
		pthread_spin_init(&lock, PTHREAD_PROCESS_PRIVATE);
	}

	void run(vertex_program &prog);

	void run(vertex_program &prog, const page_vertex &vertex) {
		if (vertex.get_id() == get_id())
			run_on_itself(prog, vertex);
		else
			run_on_neighbor(prog, vertex);
	}

	void run_on_itself(vertex_program &prog, const page_vertex &vertex);
	void run_on_neighbor(vertex_program &prog, const page_vertex &vertex);

	void run_on_message(vertex_program &prog, const vertex_message &msg1);
};

void part_topK_scan_vertex::run_on_message(vertex_program &prog,
		const vertex_message &msg1)
{
	const part_scan_msg &msg = (const part_scan_msg &) msg1;
	pthread_spin_lock(&lock);
	switch (msg.get_type()) {
		case part_scan_type::EST_LOCAL:
			est_local = msg.get();
			break;
		case part_scan_type::NEIGH:
			if (local_data == NULL) {
				scan_runtime_data_t *other_data = (scan_runtime_data_t *) msg.get();
				local_data = new scan_runtime_data_t(std::shared_ptr<neighbor_list>(
							other_data->neighbors));
			}
			break;
		default:
			ABORT_MSG("wrong message type");
	}
	pthread_spin_unlock(&lock);
}

void part_topK_scan_vertex::run_on_itself(vertex_program &prog, const page_vertex &vertex)
{
	size_t num_local_edges = vertex.get_num_edges(edge_type::BOTH_EDGES);
	if (num_local_edges == 0)
		return;

	topK_scan_vertex &self_v
		= (topK_scan_vertex &) prog.get_graph().get_vertex(get_id());
	pthread_spin_lock(&lock);
	if (est_local == 0)
		est_local = ::get_est_local_scan(prog.get_graph(), self_v, vertex);
	pthread_spin_unlock(&lock);
	broadcast_vpart(part_scan_msg(part_scan_type::EST_LOCAL, est_local));
	if (est_local < max_scan.get())
		return;

	pthread_spin_lock(&lock);
	if (local_data == NULL) {
		local_data = create_runtime(prog.get_graph(), vertex);
		assert(!local_data->neighbors->empty());
		assert(local_data->local_scan == 0);
	}
	else
		assert(local_data->local_scan == 0);
	pthread_spin_unlock(&lock);
	// TODO this is absolute hacking.
	// There is no guarantee that the message will be delivered.
	broadcast_vpart(part_scan_msg(part_scan_type::NEIGH, (size_t) local_data));

	size_t tmp = 0;
	page_byte_array::seq_const_iterator<vertex_id_t> it = vertex.get_neigh_seq_it(
			edge_type::IN_EDGE);
	PAGE_FOREACH(vertex_id_t, id, it) {
		// Ignore loops
		if (id != vertex.get_id())
			tmp++;
	} PAGE_FOREACH_END
	it = vertex.get_neigh_seq_it(edge_type::OUT_EDGE);
	PAGE_FOREACH(vertex_id_t, id, it) {
		// Ignore loops
		if (id != vertex.get_id())
			tmp++;
	} PAGE_FOREACH_END

	if (get_part_id() == 0) {
		scan_msg msg(tmp);
		msg.set_flush(true);
		prog.send_msg(get_id(), msg);
	}

	std::vector<vertex_id_t> neighbors;
	local_data->neighbors->get_neighbors(neighbors);

	size_t part_size
		= ceil(((double) neighbors.size()) / graph_conf.get_num_vparts());
	size_t start_off = std::min(part_size * get_part_id(), neighbors.size());
	size_t end_off = std::min(part_size * (get_part_id() + 1), neighbors.size());
	local_data->num_required = end_off - start_off;
	if (start_off == end_off) {
		// We need to make sure the main vertex gets enough messages.
		scan_msg msg(0);
		msg.set_flush(true);
		prog.send_msg(get_id(), msg);
		destroy_runtime(local_data);
	}
	else {
		request_vertices(neighbors.data() + start_off, end_off - start_off);
	}
}

void part_topK_scan_vertex::run_on_neighbor(vertex_program &prog,
		const page_vertex &vertex)
{
	assert(local_data);
	local_data->num_joined++;
	size_t ret = local_data->neighbors->count_edges(&vertex);
	if (ret > 0)
		local_data->local_scan += ret;

	// If we have seen all required neighbors, we have complete
	// the computation. We can release the memory now.
	if (local_data->num_joined == local_data->num_required) {
		part_local = local_data->local_scan;
		scan_msg msg(part_local);
		msg.set_flush(true);
		prog.send_msg(get_id(), msg);
		destroy_runtime(local_data);
	}
}

void part_topK_scan_vertex::run(vertex_program &prog)
{
	vertex_id_t id = prog.get_vertex_id(*this);
	bool req_itself = false;
	// If we have computed local scan on the vertex, skip the vertex.
	if (part_local >= 0)
		return;
	// If we have estimated the local scan, we should use the estimated one.
	else if (est_local > 0)
		req_itself = est_local > max_scan.get();
	else {
		// If this is the first time to compute on the vertex, we can still
		// skip a lot of vertices with this condition.
		size_t num_local_edges = prog.get_num_edges(id);
		req_itself = num_local_edges * num_local_edges >= max_scan.get();
	}
	if (req_itself)
		request_vertices(&id, 1);
}

class vertex_size_scheduler: public vertex_scheduler
{
public:
	void schedule(vertex_program &prog,
			std::vector<compute_vertex_pointer> &vertices);
};

void vertex_size_scheduler::schedule(vertex_program &prog,
		std::vector<compute_vertex_pointer> &vertices)
{
	class comp_size
	{
		vertex_program *prog;

		vsize_t get_degree(compute_vertex_pointer v) const {
			vertex_id_t id = prog->get_vertex_id(v);
			assert(id != INVALID_VERTEX_ID);
			return prog->get_num_edges(id);
		}
	public:
		comp_size(vertex_program &prog) {
			this->prog = &prog;
		}

		bool operator()(compute_vertex_pointer v1, compute_vertex_pointer v2) {
			return get_degree(v1) > get_degree(v2);
		}
	};

	std::sort(vertices.begin(), vertices.end(), comp_size(prog));
}

}

FG_vector<std::pair<vertex_id_t, size_t> >::ptr compute_topK_scan(
		FG_graph::ptr fg, size_t topK)
{
	struct timeval end;
	gettimeofday(&graph_start, NULL);
	graph_index::ptr index = NUMA_graph_index<topK_scan_vertex,
		part_topK_scan_vertex>::create(fg->get_graph_header());
	graph_engine::ptr graph = fg->create_engine(index);

	BOOST_LOG_TRIVIAL(info) << "scan statistics starts";
	BOOST_LOG_TRIVIAL(info) << "prof_file: " << graph_conf.get_prof_file();
#ifdef PROFILER
	if (!graph_conf.get_prof_file().empty())
		ProfilerStart(graph_conf.get_prof_file().c_str());
#endif

	// Let's schedule the order of processing activated vertices according
	// to the size of vertices. We start with processing vertices with higher
	// degrees in the hope we can find the max scan as early as possible,
	// so that we can simple ignore the rest of vertices.
	graph->set_vertex_scheduler(vertex_scheduler::ptr(
				new vertex_size_scheduler()));
	graph->set_max_processing_vertices(3);

	class remove_small_filter: public vertex_filter
	{
		size_t min;
	public:
		remove_small_filter(size_t min) {
			this->min = min;
		}

		bool keep(vertex_program &prog, compute_vertex &v) {
			vertex_id_t id = prog.get_vertex_id(v);
			return prog.get_num_edges(id) >= min;
		}
	};

	size_t min_edges = 1;
	std::shared_ptr<vertex_filter> filter
		= std::shared_ptr<vertex_filter>(new remove_small_filter(min_edges));
	BOOST_LOG_TRIVIAL(info)
		<< boost::format("Computing local scan on at least %1% vertices")
		% topK;
	while (known_scans.get_size() < topK) {
		graph->start(filter);
		graph->wait4complete();
		BOOST_LOG_TRIVIAL(info)
			<< boost::format("There are %1% computed vertices")
			% known_scans.get_size();
		BOOST_LOG_TRIVIAL(info) << "global max scan: " << max_scan.get();
		max_scan = global_max(0);
	}

	class remove_small_scan_filter: public vertex_filter
	{
		size_t min;
	public:
		remove_small_scan_filter(size_t min) {
			this->min = min;
		}

		bool keep(vertex_program &prog, compute_vertex &v) {
			vertex_id_t id = prog.get_vertex_id(v);
			size_t num_local_edges = prog.get_num_edges(id);
			return num_local_edges * num_local_edges >= min;
		}
	};

	size_t prev_topK_scan, curr_topK_scan;
	off_t prev_start_loc, curr_start_loc;
	do {
		prev_topK_scan = known_scans.get(topK - 1).second;
		for (prev_start_loc = topK - 1; prev_start_loc > 0
				&& known_scans.get(prev_start_loc).second == prev_topK_scan;
				prev_start_loc--);
		prev_start_loc++;
		assert(known_scans.get(prev_start_loc).second == prev_topK_scan);
		BOOST_LOG_TRIVIAL(info)
			<< boost::format("prev topK scan: %1%, prev loc: %2%")
			% prev_topK_scan % prev_start_loc;

		// Let's use the topK as the max scan for unknown vertices
		// and see if we can find a new vertex that has larger local scan.
		max_scan = global_max(prev_topK_scan);
		graph->start(std::shared_ptr<vertex_filter>(
					new remove_small_scan_filter(prev_topK_scan)));
		graph->wait4complete();
		BOOST_LOG_TRIVIAL(info)
			<< boost::format("There are %1% computed vertices")
			% known_scans.get_size();

		// If the previous topK is different from the current one,
		// it means we have found new local scans that are larger
		// than the previous topK. We should use the new topK and
		// try again.
		curr_topK_scan = known_scans.get(topK - 1).second;
		for (curr_start_loc = topK - 1; curr_start_loc > 0
				&& known_scans.get(curr_start_loc).second == curr_topK_scan;
				curr_start_loc--);
		curr_start_loc++;
		assert(known_scans.get(curr_start_loc).second == curr_topK_scan);
		BOOST_LOG_TRIVIAL(info)
			<< boost::format("global max scan: %1%, topK scan: %2%, start loc: %3%")
			% max_scan.get() % curr_topK_scan % curr_start_loc;
	} while (prev_topK_scan != curr_topK_scan || prev_start_loc != curr_start_loc);
	assert(known_scans.get_size() >= topK);

#ifdef PROFILER
	if (!graph_conf.get_prof_file().empty())
		ProfilerStop();
#endif
	gettimeofday(&end, NULL);
	BOOST_LOG_TRIVIAL(info)
		<< boost::format("It takes %1% seconds for top %2%")
		% time_diff(graph_start, end) % topK;

	FG_vector<std::pair<vertex_id_t, size_t> >::ptr vec
		= FG_vector<std::pair<vertex_id_t, size_t> >::create(topK);
	for (size_t i = 0; i < topK; i++)
		vec->set(i, known_scans.get(i));
	return vec;
}
