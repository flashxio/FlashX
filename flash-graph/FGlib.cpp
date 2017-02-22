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

#include <unordered_set>

#include "io_interface.h"

#include "FGlib.h"
#include "vertex.h"
#include "in_mem_storage.h"
#include "vertex_index.h"
#include "safs_file.h"

using namespace safs;

namespace fg
{

FG_graph::ptr FG_graph::create(const std::string &graph_file,
		const std::string &index_file, config_map::ptr configs)
{
	bool init_io = is_safs_init();

	/*
	 * The order of searching for the graph file and index file is memory,
	 * SAFS and then local Linux filesystem.
	 *
	 * If the graph file and index file have been loaded to memory, use
	 * the in-memory copy.
	 * If the graph file and index file are in SAFS and users want to run
	 * in-memory mode, read graph data and index from SAFS.
	 * If the graph file and index file are in the local Linux filesystem,
	 * read graph data and index from the local filesystem and maintain
	 * a in-memory copy.
	 */

	bool graph_in_safs = false;
	bool index_in_safs = false;
	if (init_io) {
		const RAID_config &raid_conf = get_sys_RAID_conf();
		safs_file safs_graph(raid_conf, graph_file);
		safs_file safs_index(raid_conf, index_file);
		graph_in_safs = safs_graph.exist();
		index_in_safs = safs_index.exist();
	}

	in_mem_graph::ptr graph_data;
	if (graph_conf.use_in_mem_graph() && graph_in_safs)
		graph_data = in_mem_graph::load_safs_graph(graph_file);
	else if (!graph_in_safs)
		// If we can't initialize SAFS, we assume the graph file is
		// in the local filesystem.
		graph_data = in_mem_graph::load_graph(graph_file);

	vertex_index::ptr index_data;
	if (index_in_safs)
		index_data = vertex_index::safs_load(index_file);
	else
		// If we can't initialize SAFS, we assume the index file is
		// in the local filesystem.
		index_data = vertex_index::load(index_file);

	if (graph_data)
		return ptr(new FG_graph(graph_data, index_data, graph_file, configs));
	else
		return ptr(new FG_graph(graph_file, index_data, configs));
}

FG_graph::FG_graph(const std::string &graph_file, vertex_index::ptr index_data,
		config_map::ptr configs)
{
	this->graph_file = graph_file;
	this->index_data = index_data;
	this->configs = configs;
	this->header = index_data->get_graph_header();
}

FG_graph::FG_graph(in_mem_graph::ptr graph_data, vertex_index::ptr index_data,
		const std::string &graph_name, config_map::ptr configs)
{
	this->graph_data = graph_data;
	this->index_data = index_data;
	this->configs = configs;
	graph_file = graph_name;
	header = index_data->get_graph_header();
}

graph_engine::ptr FG_graph::create_engine(graph_index::ptr index)
{
	return graph_engine::create(*this, index);
}

file_io_factory::shared_ptr FG_graph::get_graph_io_factory(int access_option)
{
	if (graph_data)
		return graph_data->create_io_factory();
	else
		// Right now only the cached I/O can support async I/O
		return create_io_factory(graph_file, access_option);
}

vertex_index::ptr FG_graph::get_index_data() const
{
	if (index_data)
		return index_data;
	else
		return vertex_index::safs_load(index_file);
}

in_mem_graph::ptr FG_graph::get_graph_data() const
{
	if (graph_data)
		return graph_data;
	else
		return in_mem_graph::load_safs_graph(graph_file);
}

/********************** Get the degree of vertices ****************************/

namespace {

class degree_vertex: public compute_vertex
{
public:
	degree_vertex(vertex_id_t id): compute_vertex(id) {
	}

	void run(vertex_program &prog);

	void run(vertex_program &prog, const page_vertex &vertex) {
	}

	void run_on_message(vertex_program &, const vertex_message &msg) {
	}
};

class degree_vertex_program: public vertex_program_impl<degree_vertex>
{
	edge_type type;
	fm::detail::mem_vec_store::ptr degree_vec;
public:
	degree_vertex_program(fm::detail::mem_vec_store::ptr degree_vec,
			edge_type type) {
		this->degree_vec = degree_vec;
		this->type = type;
	}

	edge_type get_edge_type() const {
		return type;
	}

	void set_degree(vertex_id_t id, vsize_t degree) {
		degree_vec->set(id, degree);
	}
};

class degree_vertex_program_creater: public vertex_program_creater
{
	fm::detail::mem_vec_store::ptr degree_vec;
	edge_type type;
public:
	degree_vertex_program_creater(
			fm::detail::mem_vec_store::ptr degree_vec, edge_type type) {
		this->degree_vec = degree_vec;
		this->type = type;
	}

	vertex_program::ptr create() const {
		return vertex_program::ptr(new degree_vertex_program(
					degree_vec, type));
	}
};

void degree_vertex::run(vertex_program &prog)
{
	vertex_id_t id = prog.get_vertex_id(*this);
	degree_vertex_program &degree_vprog = (degree_vertex_program &) prog;
	degree_vprog.set_degree(id, prog.get_graph().get_num_edges(id,
				degree_vprog.get_edge_type()));
}

}

fm::vector::ptr get_degree(FG_graph::ptr fg, edge_type type)
{
	graph_index::ptr index = NUMA_graph_index<degree_vertex>::create(
			fg->get_graph_header());
	graph_engine::ptr graph = fg->create_engine(index);

	fm::detail::mem_vec_store::ptr degree_vec = fm::detail::mem_vec_store::create(
			fg->get_num_vertices(), safs::params.get_num_nodes(),
			fm::get_scalar_type<vsize_t>());
	struct timeval start, end;
	gettimeofday(&start, NULL);
	graph->start_all(vertex_initializer::ptr(), vertex_program_creater::ptr(
				new degree_vertex_program_creater(degree_vec, type)));
	graph->wait4complete();
	gettimeofday(&end, NULL);
	return fm::vector::create(degree_vec);
}

/*************** Get the degree of vertices in a timestamp ********************/

namespace {

class ts_degree_vertex: public compute_vertex
{
public:
	ts_degree_vertex(vertex_id_t id): compute_vertex(id) {
	}

	void run(vertex_program &prog) {
		vertex_id_t id = prog.get_vertex_id(*this);
		request_vertices(&id, 1);
	}

	void run(vertex_program &prog, const page_vertex &vertex);

	void run_on_message(vertex_program &, const vertex_message &msg) {
	}
};

class ts_degree_vertex_program: public vertex_program_impl<ts_degree_vertex>
{
	time_t start_time;
	time_t time_interval;
	edge_type type;
	fm::detail::mem_vec_store::ptr degree_vec;
public:
	ts_degree_vertex_program(fm::detail::mem_vec_store::ptr degree_vec,
			edge_type type, time_t start_time, time_t time_interval) {
		this->degree_vec = degree_vec;
		this->type = type;
		this->start_time = start_time;
		this->time_interval = time_interval;
	}

	time_t get_start_time() const {
		return start_time;
	}

	time_t get_time_interval() const {
		return time_interval;
	}

	edge_type get_edge_type() const {
		return type;
	}

	void set_degree(vertex_id_t id, vsize_t degree) {
		degree_vec->set(id, degree);
	}
};

class ts_degree_vertex_program_creater: public vertex_program_creater
{
	time_t start_time;
	time_t time_interval;
	fm::detail::mem_vec_store::ptr degree_vec;
	edge_type type;
public:
	ts_degree_vertex_program_creater(
			fm::detail::mem_vec_store::ptr degree_vec, edge_type type,
			time_t start_time, time_t time_interval) {
		this->degree_vec = degree_vec;
		this->type = type;
		this->start_time = start_time;
		this->time_interval = time_interval;
	}

	vertex_program::ptr create() const {
		return vertex_program::ptr(new ts_degree_vertex_program(
					degree_vec, type, start_time, time_interval));
	}
};

void ts_degree_vertex::run(vertex_program &prog, const page_vertex &vertex)
{
	ts_degree_vertex_program &degree_vprog = (ts_degree_vertex_program &) prog;
	edge_type type = degree_vprog.get_edge_type();
	time_t start_time = degree_vprog.get_start_time();
	time_t time_interval = degree_vprog.get_time_interval();

	if (prog.get_graph().is_directed()) {
		const page_directed_vertex &dv = (const page_directed_vertex &) vertex;

		page_byte_array::const_iterator<ts_edge_data> begin_it
			= dv.get_data_begin<ts_edge_data>(type);
		page_byte_array::const_iterator<ts_edge_data> end_it
			= dv.get_data_end<ts_edge_data>(type);
		page_byte_array::const_iterator<ts_edge_data> ts_it = std::lower_bound(
				begin_it, end_it, ts_edge_data(start_time));
		page_byte_array::const_iterator<ts_edge_data> ts_end_it = std::lower_bound(
				begin_it, end_it, start_time + time_interval);
		degree_vprog.set_degree(vertex.get_id(), ts_end_it - ts_it);
	}
	else {
		throw unsupported_exception("undirected graph in TS graph");
	}
}

}

fm::vector::ptr get_ts_degree(FG_graph::ptr fg, edge_type type,
		time_t start_time, time_t time_interval)
{
	graph_index::ptr index = NUMA_graph_index<ts_degree_vertex>::create(
			fg->get_graph_header());
	graph_engine::ptr graph = fg->create_engine(index);
	assert(graph->get_graph_header().has_edge_data());

	fm::detail::mem_vec_store::ptr degree_vec = fm::detail::mem_vec_store::create(
			fg->get_num_vertices(), safs::params.get_num_nodes(),
			fm::get_scalar_type<vsize_t>());
	graph->start_all(vertex_initializer::ptr(), vertex_program_creater::ptr(
				new ts_degree_vertex_program_creater(degree_vec, type,
					start_time, time_interval)));
	graph->wait4complete();
	return fm::vector::create(degree_vec);
}

/************** Get the time range of the time-series graph *******************/

namespace {

class time_range_vertex: public compute_vertex
{
public:
	time_range_vertex(vertex_id_t id): compute_vertex(id) {
	}

	void run(vertex_program &prog) {
		vertex_id_t id = prog.get_vertex_id(*this);
		request_vertices(&id, 1);
	}

	void run(vertex_program &prog, const page_vertex &vertex);

	void run_on_message(vertex_program &, const vertex_message &msg) {
	}
};

class time_range_vertex_program: public vertex_program_impl<time_range_vertex>
{
	time_t start_time;
	time_t end_time;
public:
	time_range_vertex_program() {
		start_time = std::numeric_limits<time_t>::max();
		end_time = std::numeric_limits<time_t>::min();
	}

	time_t get_start_time() const {
		return start_time;
	}

	time_t get_end_time() const {
		return end_time;
	}

	void set_start_time(time_t start_time) {
		this->start_time = std::min(start_time, this->start_time);
	}

	void set_end_time(time_t end_time) {
		this->end_time = std::max(end_time, this->end_time);
	}
};

class time_range_vertex_program_creater: public vertex_program_creater
{
public:
	time_range_vertex_program_creater() {
	}

	vertex_program::ptr create() const {
		return vertex_program::ptr(new time_range_vertex_program());
	}
};

void time_range_vertex::run(vertex_program &prog, const page_vertex &vertex)
{
	time_range_vertex_program &range_vprog = (time_range_vertex_program &) prog;
	if (prog.get_graph().is_directed()) {
		const page_directed_vertex &dv = (const page_directed_vertex &) vertex;

		if (dv.get_num_edges(IN_EDGE) > 0) {
			page_byte_array::const_iterator<ts_edge_data> it
				= dv.get_data_begin<ts_edge_data>(IN_EDGE);
			range_vprog.set_start_time((*it).get_timestamp());
			it = it + (dv.get_num_edges(IN_EDGE) - 1);
			range_vprog.set_end_time((*it).get_timestamp());
		}

		if (dv.get_num_edges(OUT_EDGE) > 0) {
			page_byte_array::const_iterator<ts_edge_data> it
				= dv.get_data_begin<ts_edge_data>(OUT_EDGE);
			range_vprog.set_start_time((*it).get_timestamp());
			it = it + (dv.get_num_edges(OUT_EDGE) - 1);
			range_vprog.set_end_time((*it).get_timestamp());
		}
	}
	else {
		throw unsupported_exception("undirected graph in TS graph");
	}
}

}

std::pair<time_t, time_t> get_time_range(FG_graph::ptr fg)
{
	graph_index::ptr index = NUMA_graph_index<time_range_vertex>::create(
			fg->get_graph_header());
	graph_engine::ptr graph = fg->create_engine(index);
	assert(graph->get_graph_header().has_edge_data());

	graph->start_all(vertex_initializer::ptr(), vertex_program_creater::ptr(
				new time_range_vertex_program_creater()));
	graph->wait4complete();

	std::vector<vertex_program::ptr> vprogs;
	graph->get_vertex_programs(vprogs);
	time_t start_time = std::numeric_limits<time_t>::max();
	time_t end_time = std::numeric_limits<time_t>::min();
	BOOST_FOREACH(vertex_program::ptr prog, vprogs) {
		start_time = std::min(start_time,
				((time_range_vertex_program &)(*prog)).get_start_time());
		end_time = std::max(end_time,
				((time_range_vertex_program &)(*prog)).get_end_time());
	}
	assert(start_time <= end_time);
	return std::pair<time_t, time_t>(start_time, end_time);
}

}
