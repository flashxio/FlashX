/*
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
#include <unordered_map>

#include "graph_engine.h"
#include "graph_config.h"
#include "FGlib.h"
#include "ts_graph.h"

using namespace fg;

namespace {

class component_message: public vertex_message
{
	int id;
public:
	component_message(vertex_id_t id): vertex_message(
			sizeof(component_message), true) {
		this->id = id;
	}

	vertex_id_t get_id() const {
		return id;
	}
};

class wcc_vertex: public compute_directed_vertex
{
protected:
	bool empty;
private:
	bool updated;
	vertex_id_t component_id;
public:
	wcc_vertex(vertex_id_t id): compute_directed_vertex(id) {
		component_id = id;
		updated = true;
		empty = true;
	}

	bool is_empty(graph_engine &graph) const {
		return empty;
	}

	vertex_id_t get_component_id() const {
		return component_id;
	}

	void run(vertex_program &prog) {
		if (updated) {
			vertex_id_t id = prog.get_vertex_id(*this);
			request_vertices(&id, 1);
			updated = false;
		}
	}

	void run(vertex_program &prog, const page_vertex &vertex);

	void run_on_message(vertex_program &, const vertex_message &msg1) {
		component_message &msg = (component_message &) msg1;
		if (msg.get_id() < component_id) {
			updated = true;
			component_id = msg.get_id();
		}
	}

	vertex_id_t get_result() const {
		if (!empty)
			return get_component_id();
		else
			return INVALID_VERTEX_ID;
	}
};

template<class vertex_type>
class wcc_vertex_program: public vertex_program_impl<vertex_type>
{
	std::vector<vertex_id_t> buf;
public:
	typedef std::shared_ptr<wcc_vertex_program<vertex_type> > ptr;

	static ptr cast2(vertex_program::ptr prog) {
		return std::static_pointer_cast<wcc_vertex_program<vertex_type>,
			   vertex_program>(prog);
	}

	std::vector<vertex_id_t> &get_buf() {
		return buf;
	}
};

template<class vertex_type>
class wcc_vertex_program_creater: public vertex_program_creater
{
public:
	vertex_program::ptr create() const {
		return vertex_program::ptr(new wcc_vertex_program<vertex_type>());
	}
};

/*
 * This function get all unique neighbors on two edge lists.
 */
size_t get_unique_neighbors(edge_seq_iterator &it1, edge_seq_iterator &it2,
		std::vector<vertex_id_t> &buf, size_t max_size)
{
	vertex_id_t v1 = INVALID_VERTEX_ID;
	vertex_id_t v2 = INVALID_VERTEX_ID;
	if (it1.has_next())
		v1 = it1.next();
	if (it2.has_next())
		v2 = it2.next();
	while (it1.has_next() && it2.has_next() && buf.size() < max_size) {
		if (v1 == v2) {
			buf.push_back(v1);
			v1 = it1.next();
			v2 = it2.next();
		}
		else if (v1 > v2) {
			buf.push_back(v2);
			v2 = it2.next();
		}
		else {
			buf.push_back(v1);
			v1 = it1.next();
		}
	}
	if (v1 != INVALID_VERTEX_ID)
		buf.push_back(v1);
	if (v2 != INVALID_VERTEX_ID)
		buf.push_back(v2);
	while (it1.has_next() && buf.size() < max_size)
		buf.push_back(it1.next());
	while (it2.has_next() && buf.size() < max_size)
		buf.push_back(it2.next());
	return buf.size();
}

void wcc_vertex::run(vertex_program &prog, const page_vertex &vertex)
{
	const size_t BUF_SIZE = 512 * 1024 / sizeof(vertex_id_t) - 1;
	component_message msg(component_id);
	const page_directed_vertex &dvertex = (const page_directed_vertex &) vertex;
	assert(dvertex.has_in_part());
	assert(dvertex.has_out_part());
	empty = (vertex.get_num_edges(BOTH_EDGES) == 0);
	wcc_vertex_program<wcc_vertex> &wcc_vprog
		= (wcc_vertex_program<wcc_vertex> &) prog;
	std::vector<vertex_id_t> &buf = wcc_vprog.get_buf();
	edge_seq_iterator in_it = dvertex.get_neigh_seq_it(IN_EDGE);
	edge_seq_iterator out_it = dvertex.get_neigh_seq_it(OUT_EDGE);
	do {
		buf.clear();
		get_unique_neighbors(in_it, out_it, buf, BUF_SIZE);
		prog.multicast_msg(buf.data(), buf.size(), msg);
	} while (in_it.has_next() || out_it.has_next());
}

class cc_vertex: public compute_vertex
{
	bool empty;
	bool updated;
	vertex_id_t component_id;
public:
	cc_vertex(vertex_id_t id): compute_vertex(id) {
		component_id = id;
		updated = true;
		empty = true;
	}

	bool is_empty(graph_engine &graph) const {
		return empty;
	}

	vertex_id_t get_component_id() const {
		return component_id;
	}

	void run(vertex_program &prog) {
		if (updated) {
			vertex_id_t id = prog.get_vertex_id(*this);
			request_vertices(&id, 1);
			updated = false;
		}
	}

	void run(vertex_program &prog, const page_vertex &vertex);

	void run_on_message(vertex_program &, const vertex_message &msg1) {
		component_message &msg = (component_message &) msg1;
		if (msg.get_id() < component_id) {
			updated = true;
			component_id = msg.get_id();
		}
	}

	vertex_id_t get_result() const {
		if (!empty)
			return get_component_id();
		else
			return INVALID_VERTEX_ID;
	}
};

template<class vertex_type>
class cc_vertex_program: public vertex_program_impl<vertex_type>
{
	std::vector<vertex_id_t> buf;
public:
	typedef std::shared_ptr<cc_vertex_program<vertex_type> > ptr;

	static ptr cast2(vertex_program::ptr prog) {
		return std::static_pointer_cast<cc_vertex_program<vertex_type>,
			   vertex_program>(prog);
	}

	std::vector<vertex_id_t> &get_buf() {
		return buf;
	}
};

template<class vertex_type>
class cc_vertex_program_creater: public vertex_program_creater
{
public:
	vertex_program::ptr create() const {
		return vertex_program::ptr(new cc_vertex_program<vertex_type>());
	}
};

void cc_vertex::run(vertex_program &prog, const page_vertex &vertex)
{
	component_message msg(component_id);
	empty = (vertex.get_num_edges(BOTH_EDGES) == 0);
	edge_seq_iterator it = vertex.get_neigh_seq_it(OUT_EDGE);
	prog.multicast_msg(it, msg);
}

class sync_wcc_vertex: public compute_directed_vertex
{
	bool empty;
	bool updated;
	vertex_id_t component_id;
	vertex_id_t new_id;
public:
	sync_wcc_vertex(vertex_id_t id): compute_directed_vertex(id) {
		new_id = component_id = id;
		updated = true;
		empty = true;
	}

	bool is_empty(graph_engine &graph) const {
		return empty;
	}

	vertex_id_t get_component_id() const {
		return component_id;
	}

	void run(vertex_program &prog) {
		if (updated) {
			vertex_id_t id = prog.get_vertex_id(*this);
			request_vertices(&id, 1);
			updated = false;
		}
	}

	void run(vertex_program &prog, const page_vertex &vertex);

	void run_on_message(vertex_program &vprog, const vertex_message &msg1) {
		component_message &msg = (component_message &) msg1;
		if (msg.get_id() < new_id) {
			new_id = msg.get_id();
			vprog.request_notify_iter_end(*this);
		}
	}

	void notify_iteration_end(vertex_program &vprog) {
		if (new_id < component_id) {
			component_id = new_id;
			updated = true;
		}
	}

	vertex_id_t get_result() const {
		if (!empty)
			return get_component_id();
		else
			return INVALID_VERTEX_ID;
	}
};

void sync_wcc_vertex::run(vertex_program &prog, const page_vertex &vertex)
{
	const size_t BUF_SIZE = 512 * 1024 / sizeof(vertex_id_t) - 1;
	component_message msg(component_id);
	const page_directed_vertex &dvertex = (const page_directed_vertex &) vertex;
	assert(dvertex.has_in_part());
	assert(dvertex.has_out_part());
	empty = (vertex.get_num_edges(BOTH_EDGES) == 0);
	wcc_vertex_program<sync_wcc_vertex> &wcc_vprog
		= (wcc_vertex_program<sync_wcc_vertex> &) prog;
	std::vector<vertex_id_t> &buf = wcc_vprog.get_buf();
	edge_seq_iterator in_it = dvertex.get_neigh_seq_it(IN_EDGE);
	edge_seq_iterator out_it = dvertex.get_neigh_seq_it(OUT_EDGE);
	do {
		buf.clear();
		get_unique_neighbors(in_it, out_it, buf, BUF_SIZE);
		prog.multicast_msg(buf.data(), buf.size(), msg);
	} while (in_it.has_next() || out_it.has_next());
}

class ts_wcc_vertex: public wcc_vertex
{
public:
	ts_wcc_vertex(vertex_id_t id): wcc_vertex(id) {
	}

	void run(vertex_program &prog) {
		wcc_vertex::run(prog);
	}

	void run(vertex_program &prog, const page_vertex &vertex);
};

class ts_wcc_vertex_program: public vertex_program_impl<ts_wcc_vertex>
{
	time_t start_time;
	time_t time_interval;
public:
	ts_wcc_vertex_program(time_t start_time, time_t time_interval) {
		this->start_time = start_time;
		this->time_interval = time_interval;
	}

	time_t get_start_time() const {
		return start_time;
	}

	time_t get_time_interval() const {
		return time_interval;
	}
};

class ts_wcc_vertex_program_creater: public vertex_program_creater
{
	time_t start_time;
	time_t time_interval;
public:
	ts_wcc_vertex_program_creater(time_t start_time, time_t time_interval) {
		this->start_time = start_time;
		this->time_interval = time_interval;
	}

	vertex_program::ptr create() const {
		return vertex_program::ptr(new ts_wcc_vertex_program(
					start_time, time_interval));
	}
};

void ts_wcc_vertex::run(vertex_program &prog, const page_vertex &vertex)
{
	assert(prog.get_graph().is_directed());
	ts_wcc_vertex_program &wcc_vprog = (ts_wcc_vertex_program &) prog;
	const page_directed_vertex &dvertex = (const page_directed_vertex &) vertex;
	edge_seq_iterator in_it = get_ts_iterator(dvertex, edge_type::IN_EDGE,
			wcc_vprog.get_start_time(), wcc_vprog.get_time_interval());
	edge_seq_iterator out_it = get_ts_iterator(dvertex, edge_type::OUT_EDGE,
			wcc_vprog.get_start_time(), wcc_vprog.get_time_interval());
	empty = !in_it.has_next() && !out_it.has_next();
	component_message msg(get_component_id());
	prog.multicast_msg(in_it, msg);
	prog.multicast_msg(out_it, msg);
}

}

#include "save_result.h"

namespace fg
{

fm::vector::ptr compute_cc(FG_graph::ptr fg)
{
	bool directed = fg->get_graph_header().is_directed_graph();
	if (directed) {
		BOOST_LOG_TRIVIAL(error)
			<< "This algorithm works on an undirected graph";
		return fm::vector::ptr();
	}

	struct timeval start, end;
	gettimeofday(&start, NULL);
	graph_index::ptr index = NUMA_graph_index<cc_vertex>::create(
			fg->get_graph_header());
	graph_engine::ptr graph = fg->create_engine(index);
	BOOST_LOG_TRIVIAL(info) << "connected components starts";
#ifdef PROFILER
	if (!graph_conf.get_prof_file().empty())
		ProfilerStart(graph_conf.get_prof_file().c_str());
#endif
	if (graph->is_directed()) {
		fprintf(stderr, "cc has to run on an undirected graph\n");
		return fm::vector::ptr();
	}

	graph->start_all(vertex_initializer::ptr(),
			vertex_program_creater::ptr(new cc_vertex_program_creater<cc_vertex>()));
	graph->wait4complete();
	gettimeofday(&end, NULL);
	BOOST_LOG_TRIVIAL(info)
		<< boost::format("WCC takes %1% seconds in total")
		% time_diff(start, end);

#ifdef PROFILER
	if (!graph_conf.get_prof_file().empty())
		ProfilerStop();
#endif

	fm::detail::mem_vec_store::ptr res_store = fm::detail::mem_vec_store::create(
			fg->get_num_vertices(), safs::params.get_num_nodes(),
			fm::get_scalar_type<vertex_id_t>());
	graph->query_on_all(vertex_query::ptr(
				new save_query<vertex_id_t, cc_vertex>(res_store)));
	return fm::vector::create(res_store);
}

fm::vector::ptr compute_wcc(FG_graph::ptr fg)
{
	bool directed = fg->get_graph_header().is_directed_graph();
	if (!directed) {
		BOOST_LOG_TRIVIAL(error)
			<< "This algorithm works on a directed graph";
		return fm::vector::ptr();
	}

	graph_index::ptr index = NUMA_graph_index<wcc_vertex>::create(
			fg->get_graph_header());
	graph_engine::ptr graph = fg->create_engine(index);
	BOOST_LOG_TRIVIAL(info) << "weakly connected components starts";
#ifdef PROFILER
	if (!graph_conf.get_prof_file().empty())
		ProfilerStart(graph_conf.get_prof_file().c_str());
#endif
	if (!graph->is_directed()) {
		fprintf(stderr, "wcc has to run on a directed graph\n");
		return fm::vector::ptr();
	}

	struct timeval start, end;
	gettimeofday(&start, NULL);
	graph->start_all(vertex_initializer::ptr(),
			vertex_program_creater::ptr(new wcc_vertex_program_creater<wcc_vertex>()));
	graph->wait4complete();
	gettimeofday(&end, NULL);
	BOOST_LOG_TRIVIAL(info)
		<< boost::format("WCC takes %1% seconds in total")
		% time_diff(start, end);

#ifdef PROFILER
	if (!graph_conf.get_prof_file().empty())
		ProfilerStop();
#endif

	fm::detail::mem_vec_store::ptr res_store = fm::detail::mem_vec_store::create(
			fg->get_num_vertices(), safs::params.get_num_nodes(),
			fm::get_scalar_type<vertex_id_t>());
	graph->query_on_all(vertex_query::ptr(
				new save_query<vertex_id_t, wcc_vertex>(res_store)));
	return fm::vector::create(res_store);
}

fm::vector::ptr compute_sync_wcc(FG_graph::ptr fg)
{
	bool directed = fg->get_graph_header().is_directed_graph();
	if (!directed) {
		BOOST_LOG_TRIVIAL(error)
			<< "This algorithm works on a directed graph";
		return fm::vector::ptr();
	}

	graph_index::ptr index = NUMA_graph_index<sync_wcc_vertex>::create(
			fg->get_graph_header());
	graph_engine::ptr graph = fg->create_engine(index);
	BOOST_LOG_TRIVIAL(info) << "synchronous weakly connected components starts";
#ifdef PROFILER
	if (!graph_conf.get_prof_file().empty())
		ProfilerStart(graph_conf.get_prof_file().c_str());
#endif
	if (!graph->is_directed()) {
		fprintf(stderr, "wcc has to run on a directed graph\n");
		return fm::vector::ptr();
	}

	struct timeval start, end;
	gettimeofday(&start, NULL);
	graph->start_all(vertex_initializer::ptr(),
			vertex_program_creater::ptr(new wcc_vertex_program_creater<sync_wcc_vertex>()));
	graph->wait4complete();
	gettimeofday(&end, NULL);
	BOOST_LOG_TRIVIAL(info)
		<< boost::format("WCC takes %1% seconds in total")
		% time_diff(start, end);

#ifdef PROFILER
	if (!graph_conf.get_prof_file().empty())
		ProfilerStop();
#endif

	fm::detail::mem_vec_store::ptr res_store = fm::detail::mem_vec_store::create(
			fg->get_num_vertices(), safs::params.get_num_nodes(),
			fm::get_scalar_type<vertex_id_t>());
	graph->query_on_all(vertex_query::ptr(
				new save_query<vertex_id_t, sync_wcc_vertex>(res_store)));
	return fm::vector::create(res_store);
}

fm::vector::ptr compute_ts_wcc(FG_graph::ptr fg,
		time_t start_time, time_t time_interval)
{
	bool directed = fg->get_graph_header().is_directed_graph();
	if (!directed) {
		BOOST_LOG_TRIVIAL(error)
			<< "This algorithm works on a directed graph";
		return fm::vector::ptr();
	}

	graph_index::ptr index = NUMA_graph_index<ts_wcc_vertex>::create(
			fg->get_graph_header());
	graph_engine::ptr graph = fg->create_engine(index);
	assert(graph->get_graph_header().has_edge_data());
	BOOST_LOG_TRIVIAL(info) << "TS weakly connected components starts";
#ifdef PROFILER
	if (!graph_conf.get_prof_file().empty())
		ProfilerStart(graph_conf.get_prof_file().c_str());
#endif
	if (!graph->is_directed()) {
		fprintf(stderr, "wcc has to run on a directed graph\n");
		return fm::vector::ptr();
	}

	struct timeval start, end;
	gettimeofday(&start, NULL);
	graph->start_all(vertex_initializer::ptr(), vertex_program_creater::ptr(
				new ts_wcc_vertex_program_creater(start_time, time_interval)));
	graph->wait4complete();
	gettimeofday(&end, NULL);
	BOOST_LOG_TRIVIAL(info)
		<< boost::format("WCC takes %1% seconds in total")
		% time_diff(start, end);

#ifdef PROFILER
	if (!graph_conf.get_prof_file().empty())
		ProfilerStop();
#endif

	fm::detail::mem_vec_store::ptr res_store = fm::detail::mem_vec_store::create(
			fg->get_num_vertices(), safs::params.get_num_nodes(),
			fm::get_scalar_type<vertex_id_t>());
	graph->query_on_all(vertex_query::ptr(
				new save_query<vertex_id_t, wcc_vertex>(res_store)));
	return fm::vector::create(res_store);
}

}
