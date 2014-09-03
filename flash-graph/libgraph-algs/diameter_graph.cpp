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

#include <vector>
#include <set>
#include <unordered_map>

#include "thread.h"
#include "io_interface.h"
#include "container.h"
#include "concurrency.h"

#include "vertex_index.h"
#include "graph_engine.h"
#include "graph_config.h"
#include "FGlib.h"

namespace {

size_t num_bfs = 1;
edge_type traverse_edge = edge_type::OUT_EDGE;
int start_level;

template<class T>
class embedded_bitmap
{
	T map;
public:
	embedded_bitmap() {
		map = 0;
	}

	void set(int idx) {
		assert((size_t) idx < sizeof(T) * 8);
		map |= ((T) 1) << idx;
	}

	void clear() {
		map = 0;
	}

	bool merge(const embedded_bitmap<T> &map) {
		bool old = this->map;
		this->map |= map.map;
		return old != map.map;
	}

	bool operator!=(const embedded_bitmap<T> &map) {
		return this->map != map.map;
	}

	T get_value() const {
		return map;
	}
};

class diameter_message: public vertex_message
{
	embedded_bitmap<uint64_t> bfs_ids;
public:
	diameter_message(embedded_bitmap<uint64_t> bfs_ids): vertex_message(
			sizeof(diameter_message), true) {
		this->bfs_ids = bfs_ids;
	}

	embedded_bitmap<uint64_t> get_bfs_ids() const {
		return bfs_ids;
	}
};

class diameter_vertex: public compute_directed_vertex
{
	embedded_bitmap<uint64_t> bfs_ids;
	embedded_bitmap<uint64_t> new_bfs_ids;
	// The largest distance from a start vertex among the BFS.
	short max_dist;
	bool updated;
public:
	diameter_vertex(vertex_id_t id): compute_directed_vertex(id) {
		updated = false;
		max_dist = 0;
	}

	short get_max_dist() const {
		return max_dist;
	}

	void init(int bfs_id) {
		max_dist = 0;
		bfs_ids.clear();
		bfs_ids.set(bfs_id);
		new_bfs_ids = bfs_ids;
		updated = true;
		printf("BFS %d starts at v%u\n", bfs_id, get_id());
	}

	void reset() {
		max_dist = 0;
		updated = false;
		bfs_ids.clear();
		new_bfs_ids.clear();
	}

	void run(vertex_program &prog) {
		if (updated) {
			updated = false;
			directed_vertex_request req(get_id(), traverse_edge);
			request_partial_vertices(&req, 1);
		}
	}

	void run(vertex_program &prog, const page_vertex &vertex);

	void run_on_message(vertex_program &vprog, const vertex_message &msg) {
		const diameter_message &dmsg = (const diameter_message &) msg;
		if (new_bfs_ids.merge(dmsg.get_bfs_ids()))
			vprog.request_notify_iter_end(*this);
	}

	void notify_iteration_end(vertex_program &vprog);
};

typedef std::pair<vertex_id_t, int> vertex_dist_t;

template<class vertex_type>
class diameter_vertex_program: public vertex_program_impl<vertex_type>
{
	int curr_iter;
	std::vector<vertex_dist_t> curr_vertices;
	std::deque<vertex_dist_t> prev_vertices;
public:
	typedef std::shared_ptr<diameter_vertex_program<vertex_type> > ptr;

	static ptr cast2(vertex_program::ptr prog) {
		return std::static_pointer_cast<diameter_vertex_program<vertex_type>,
			   vertex_program>(prog);
	}

	diameter_vertex_program() {
		curr_iter = 0;
	}

	void set_max_dist(vertex_type &v, int iter_no) {
		if (curr_iter == iter_no) {
			if (curr_vertices.size() < num_bfs)
				curr_vertices.push_back(vertex_dist_t(v.get_id(), curr_iter));
		}
		else {
			assert(curr_iter < iter_no);
			curr_iter = iter_no;
			prev_vertices.insert(prev_vertices.end(),
					curr_vertices.begin(), curr_vertices.end());
			curr_vertices.clear();
			curr_vertices.push_back(vertex_dist_t(v.get_id(), curr_iter));

			// Remove vertices in the previous iterations.
			while (prev_vertices.size() > num_bfs)
				prev_vertices.pop_front();
		}
	}

	void get_max_dist_vertices(std::vector<vertex_dist_t> &vertices) const {
		vertices.insert(vertices.end(), curr_vertices.begin(),
				curr_vertices.end());
		for (std::deque<vertex_dist_t>::const_reverse_iterator it
				= prev_vertices.crbegin(); it != prev_vertices.crend(); it++)
			vertices.push_back(*it);
	}
};

void diameter_vertex::run(vertex_program &prog, const page_vertex &vertex)
{
	int num_dests = vertex.get_num_edges(traverse_edge);
	if (num_dests == 0)
		return;

	// We need to add the neighbors of the vertex to the queue of
	// the next level.
	diameter_message msg(bfs_ids);
	if (traverse_edge == BOTH_EDGES) {
		edge_seq_iterator it = vertex.get_neigh_seq_it(IN_EDGE);
		prog.multicast_msg(it, msg);
		it = vertex.get_neigh_seq_it(OUT_EDGE);
		prog.multicast_msg(it, msg);
	}
	else {
		edge_seq_iterator it = vertex.get_neigh_seq_it(traverse_edge);
		prog.multicast_msg(it, msg);
	}
}

void diameter_vertex::notify_iteration_end(vertex_program &vprog)
{
	if (bfs_ids != new_bfs_ids) {
		updated = true;
		int iter_no = vprog.get_graph().get_curr_level() + 1 - start_level;
		max_dist = max(iter_no, max_dist);
		((diameter_vertex_program<diameter_vertex> &) vprog).set_max_dist(
			*this, iter_no);
	}
	bfs_ids = new_bfs_ids;
}

class dist_compare
{
public:
	bool operator()(const vertex_dist_t &v1, const vertex_dist_t &v2) {
		return v1.second > v2.second;
	}
};

/*
 * This diameter estimation uses a single BFS for each sweep.
 * It assumes an unweighted graph.
 */
class simple_diameter_vertex: public compute_directed_vertex
{
	short max_dist;
public:
	simple_diameter_vertex(vertex_id_t id): compute_directed_vertex(id) {
		max_dist = -1;
	}

	short get_max_dist() const {
		return max_dist;
	}

	void init(int bfs_id) {
	}

	void reset() {
		max_dist = -1;
	}

	void run(vertex_program &prog) {
		diameter_vertex_program<simple_diameter_vertex> & diam_prog
			= (diameter_vertex_program<simple_diameter_vertex> &) prog;
		if (max_dist < 0) {
			int iter_no = prog.get_graph().get_curr_level() - start_level;
			diam_prog.set_max_dist(*this, iter_no);
			max_dist = prog.get_graph().get_curr_level();
			directed_vertex_request req(get_id(), traverse_edge);
			request_partial_vertices(&req, 1);
		}
	}

	void run(vertex_program &prog, const page_vertex &vertex) {
		int num_dests = vertex.get_num_edges(traverse_edge);
		if (num_dests == 0)
			return;

		// We need to add the neighbors of the vertex to the queue of
		// the next level.
		if (traverse_edge == BOTH_EDGES) {
			edge_seq_iterator it = vertex.get_neigh_seq_it(IN_EDGE);
			prog.activate_vertices(it);
			it = vertex.get_neigh_seq_it(OUT_EDGE);
			prog.activate_vertices(it);
		}
		else {
			edge_seq_iterator it = vertex.get_neigh_seq_it(traverse_edge);
			prog.activate_vertices(it);
		}
	}

	void run_on_message(vertex_program &vprog, const vertex_message &msg) {
		assert(0);
	}
};

template<class vertex_type>
class diameter_vertex_program_creater: public vertex_program_creater
{
public:
	vertex_program::ptr create() const {
		return vertex_program::ptr(new diameter_vertex_program<vertex_type>());
	}
};

template<class vertex_type>
class diameter_reset: public vertex_initializer
{
public:
	void init(compute_vertex &v) {
		vertex_type &dv = (vertex_type &) v;
		dv.reset();
	}
};

template<class vertex_type>
class diameter_initializer: public vertex_initializer
{
	std::unordered_map<vertex_id_t, int> start_vertices;
public:
	diameter_initializer(const std::vector<vertex_id_t> &vertices) {
		for (size_t i = 0; i < vertices.size(); i++) {
			start_vertices.insert(std::pair<vertex_id_t, int>(vertices[i], i));
		}
	}

	void init(compute_vertex &v) {
		vertex_type &dv = (vertex_type &) v;
		std::unordered_map<vertex_id_t, int>::const_iterator it
			= start_vertices.find(v.get_id());
		assert(it != start_vertices.end());
		dv.init(it->second);
	}
};

template<class vertex_type>
std::vector<vertex_dist_t> estimate_diameter_1sweep(graph_engine::ptr graph,
		const std::vector<vertex_id_t> &start_vertices)
{
	start_level = graph->get_curr_level();
	graph->init_all_vertices(vertex_initializer::ptr(
				new diameter_reset<vertex_type>()));
	graph->start(start_vertices.data(), start_vertices.size(),
			vertex_initializer::ptr(
				new diameter_initializer<vertex_type>(start_vertices)),
			vertex_program_creater::ptr(
				new diameter_vertex_program_creater<vertex_type>()));
	graph->wait4complete();

	std::vector<vertex_program::ptr> vprogs;
	graph->get_vertex_programs(vprogs);
	std::vector<vertex_dist_t> max_dist_vertices;
	BOOST_FOREACH(vertex_program::ptr vprog, vprogs) {
		typename diameter_vertex_program<vertex_type>::ptr diameter_vprog
			= diameter_vertex_program<vertex_type>::cast2(vprog);
		diameter_vprog->get_max_dist_vertices(max_dist_vertices);
	}
	return max_dist_vertices;
}

}

size_t estimate_diameter(FG_graph::ptr fg, int num_para_bfs,
		bool directed, int num_sweeps)
{
	num_bfs = num_para_bfs;
	if (!directed)
		traverse_edge = edge_type::BOTH_EDGES;

	graph_index::ptr index = NUMA_graph_index<diameter_vertex>::create(
			fg->get_index_file());
	graph_engine::ptr graph = graph_engine::create(fg->get_graph_file(),
			index, fg->get_configs());
	printf("diameter estimation starts\n");
	printf("#para BFS: %d, #sweeps: %d, directed: %d\n", num_para_bfs,
			num_sweeps, directed);
	printf("prof_file: %s\n", graph_conf.get_prof_file().c_str());
#ifdef PROFILER
	if (!graph_conf.get_prof_file().empty())
		ProfilerStart(graph_conf.get_prof_file().c_str());
#endif

	std::vector<vertex_id_t> start_vertices;
	while (start_vertices.size() < num_bfs) {
		vertex_id_t id = random() % graph->get_max_vertex_id();
		start_vertices.push_back(id);
	}

	short global_max = 0;
	for (int i = 0; i < num_sweeps; i++) {
		printf("Sweep %d starts on %ld vertices, traverse edge: %d\n",
			i, start_vertices.size(), traverse_edge);
		BOOST_FOREACH(vertex_id_t v, start_vertices) {
			printf("v%d\n", v);
		}

		struct timeval start, end;
		gettimeofday(&start, NULL);
		std::vector<vertex_dist_t> max_dist_vertices;
		if (start_vertices.size() == 1)
			max_dist_vertices = estimate_diameter_1sweep<simple_diameter_vertex>(
					graph, start_vertices);
		else
			max_dist_vertices = estimate_diameter_1sweep<diameter_vertex>(
					graph, start_vertices);
		gettimeofday(&end, NULL);

		if (max_dist_vertices.empty()) {
			size_t num_bfs = start_vertices.size();
			start_vertices.clear();
			while (start_vertices.size() < num_bfs) {
				vertex_id_t id = random() % graph->get_max_vertex_id();
				start_vertices.push_back(id);
			}
		}
		else {
			std::sort(max_dist_vertices.begin(), max_dist_vertices.end(),
					dist_compare());
			assert(max_dist_vertices.front().second
					>= max_dist_vertices.back().second);
			start_vertices.clear();
			std::set<vertex_id_t> start_set;
			for (size_t i = 0; start_set.size() < num_bfs
					&& i < max_dist_vertices.size(); i++) {
				start_set.insert(max_dist_vertices[i].first);
			}
			start_vertices.insert(start_vertices.begin(), start_set.begin(),
					start_set.end());
			short max_dist = max_dist_vertices.front().second;
			printf("It takes %f seconds. The current max dist: %d\n",
					time_diff(start, end), max_dist);
			global_max = max(global_max, max_dist);
			// We should switch the direction if we search for the longest
			// directed path
			if (traverse_edge == edge_type::IN_EDGE)
				traverse_edge = edge_type::OUT_EDGE;
			else if (traverse_edge == edge_type::OUT_EDGE)
				traverse_edge = edge_type::IN_EDGE;
		}
	}

#ifdef PROFILER
	if (!graph_conf.get_prof_file().empty())
		ProfilerStop();
#endif

	return global_max;
}
