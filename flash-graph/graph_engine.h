#ifndef __GRAPH_ENGINE_H__
#define __GRAPH_ENGINE_H__

/**
 * Copyright 2013 Da Zheng
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

#include "thread.h"
#include "container.h"
#include "concurrency.h"
#include "slab_allocator.h"
#include "io_interface.h"

#include "vertex.h"
#include "vertex_index.h"
#include "trace_logger.h"
#include "messaging.h"
#include "vertex_interpreter.h"
#include "partitioner.h"
#include "graph_index.h"
#include "graph_config.h"
#include "vertex_request.h"
#include "vertex_program.h"

/**
 * The size of a message buffer used to pass vertex messages to other threads.
 */
const int GRAPH_MSG_BUF_SIZE = PAGE_SIZE * 4;

class graph_engine;
class vertex_request;

class compute_vertex
{
	vertex_id_t id;
public:
	compute_vertex(vertex_id_t id) {
		this->id = id;
	}

	void init_vertex(const vertex_index &index) {
	}

	vsize_t get_num_edges() const;

	/**
	 * This allows a vertex to request other vertices in the graph.
	 * @ids: the Ids of vertices.
	 */
	void request_vertices(vertex_id_t ids[], size_t num);

	vertex_id_t get_id() const {
		return id;
	}

	void notify_iteration_end(vertex_program &) {
	}
};

class compute_directed_vertex: public compute_vertex
{
	vsize_t num_in_edges;
public:
	compute_directed_vertex(vertex_id_t id): compute_vertex(id) {
		num_in_edges = 0;
	}

	void init_vertex(const vertex_index &index1) {
		assert(index1.get_graph_header().get_graph_type()
				== graph_type::DIRECTED);
		const directed_vertex_index &index
			= (const directed_vertex_index &) index1;
		num_in_edges = index.get_num_in_edges(get_id());
	}

	vsize_t get_num_in_edges() const {
		return num_in_edges;
	}

	vsize_t get_num_out_edges() const {
		return get_num_edges() - num_in_edges;
	}

	/**
	 * This allows a vertex to request partial vertices in the graph.
	 * @reqs: defines part of vertices..
	 */
	void request_partial_vertices(directed_vertex_request reqs[], size_t num);
};

class compute_ts_vertex: public compute_vertex
{
public:
	compute_ts_vertex(vertex_id_t id): compute_vertex(id) {
	}

	void init_vertex(const vertex_index &index) {
		assert(index.get_graph_header().get_graph_type()
				== graph_type::TS_DIRECTED
				|| index.get_graph_header().get_graph_type()
				== graph_type::TS_UNDIRECTED);
	}

	/**
	 * This allows a vertex to request partial vertices in the graph.
	 * @reqs: defines part of vertices..
	 */
	void request_partial_vertices(ts_vertex_request reqs[], size_t num);
};

class vertex_scheduler
{
public:
	typedef std::shared_ptr<vertex_scheduler> ptr;
	virtual void schedule(std::vector<vertex_id_t> &vertices) = 0;
};

/**
 * When the graph engine starts, a user can use this filter to decide
 * what vertices are activated for the first time.
 */
class vertex_filter
{
public:
	virtual bool keep(compute_vertex &v) = 0;
};

class vertex_initiator
{
public:
	typedef std::shared_ptr<vertex_initiator> ptr;
	virtual void init(compute_vertex &) = 0;
};

class graph_engine;
class vertex_query
{
public:
	typedef std::shared_ptr<vertex_query> ptr;
	virtual void run(graph_engine &, compute_vertex &v) = 0;
	virtual void merge(graph_engine &graph, vertex_query::ptr q) = 0;
	virtual ptr clone() = 0;
};

class worker_thread;

class graph_engine
{
	int vertex_header_size;
	graph_header header;
	graph_index::ptr vertices;
	std::unique_ptr<ext_mem_vertex_interpreter> interpreter;
	vertex_scheduler::ptr scheduler;

	// The number of activated vertices that haven't been processed
	// in the current level.
	atomic_number<size_t> num_remaining_vertices_in_level;
	atomic_integer level;
	volatile bool is_complete;

	// These are used for switching queues.
	pthread_mutex_t lock;
	pthread_barrier_t barrier1;
	pthread_barrier_t barrier2;

	int num_nodes;
	std::vector<worker_thread *> worker_threads;
	std::vector<vertex_program::ptr> vprograms;

	trace_logger::ptr logger;
	file_io_factory::shared_ptr factory;
	int max_processing_vertices;

	// The time when the current iteration starts.
	struct timeval start_time;

	void init_threads(vertex_program_creater::ptr creater);
protected:
	graph_engine(const std::string &graph_file, graph_index::ptr index,
			const config_map &configs);
public:
	typedef std::shared_ptr<graph_engine> ptr;
	static graph_engine::ptr create(const std::string &graph_file,
			graph_index::ptr index, const config_map &configs) {
		return graph_engine::ptr(new graph_engine(graph_file, index, configs));
	}

	~graph_engine();

	/*
	 * Query graph information.
	 */

	/**
	 * The following four variants of get_vertex return the same compute_vertex
	 * but with slightly lower overhead.
	 * These can only be used in a shared machine.
	 */

	compute_vertex &get_vertex(vertex_id_t id) {
		return vertices->get_vertex(id);
	}

	compute_vertex &get_vertex(int part_id, local_vid_t id) {
		return vertices->get_vertex(part_id, id);
	}

	size_t get_vertices(const vertex_id_t ids[], int num, compute_vertex *v_buf[]) {
		return vertices->get_vertices(ids, num, v_buf);
	}

	size_t get_vertices(int part_id, const local_vid_t ids[], int num,
			compute_vertex *v_buf[]) {
		return vertices->get_vertices(part_id, ids, num, v_buf);
	}

	/**
	 * This returns the location and the size of a vertex in the file.
	 */
	const in_mem_vertex_info get_vertex_info(vertex_id_t id) const {
		return vertices->get_vertex_info(id);
	}

	/**
	 * This returns the number of edges of a vertex.
	 */
	const vsize_t get_vertex_edges(vertex_id_t id) const {
		in_mem_vertex_info info = get_vertex_info(id);
		return (info.get_ext_mem_size() - vertex_header_size) / sizeof(vertex_id_t);
	}

	vertex_id_t get_max_vertex_id() const {
		return vertices->get_max_vertex_id();
	}

	vertex_id_t get_min_vertex_id() const {
		return vertices->get_min_vertex_id();
	}

	size_t get_num_vertices() const {
		return vertices->get_num_vertices();
	}

	bool is_directed() const {
		return header.is_directed_graph();
	}

	const graph_header &get_graph_header() const {
		return header;
	}

	void set_vertex_scheduler(vertex_scheduler::ptr scheduler);

	void start(std::shared_ptr<vertex_filter> filter,
			vertex_program_creater::ptr creater = vertex_program_creater::ptr());
	void start(vertex_id_t ids[], int num,
			vertex_initiator::ptr init = vertex_initiator::ptr(),
			vertex_program_creater::ptr creater = vertex_program_creater::ptr());
	void start_all(vertex_initiator::ptr init = vertex_initiator::ptr(),
			vertex_program_creater::ptr creater = vertex_program_creater::ptr());

	void wait4complete();

	/**
	 * This method preloads the entire graph to the page cache.
	 */
	void preload_graph();

	/**
	 * These two methods allow users to initialize vertices to certain state.
	 */
	void init_vertices(vertex_id_t ids[], int num, vertex_initiator::ptr init);
	void init_all_vertices(vertex_initiator::ptr init);

	/**
	 * This method allows users to query the information on the state of all vertices.
	 */
	void query_on_all(vertex_query::ptr query);

	void get_vertex_programs(std::vector<vertex_program::ptr> &programs) {
		programs = vprograms;
	}

	/**
	 * This returns the current iteration No.
	 */
	int get_curr_level() const {
		return level.get();
	}

	/**
	 * The methods below should be used internally.
	 */

	/**
	 * The algorithm progresses to the next level.
	 * It returns true if no more work can progress.
	 */
	bool progress_next_level();

	trace_logger::ptr get_logger() const {
		return logger;
	}

	/**
	 * Get the file id where the graph data is stored.
	 */
	int get_file_id() const {
		return factory->get_file_id();
	}

	ext_mem_vertex_interpreter &get_vertex_interpreter() const {
		return *interpreter;
	}

	const graph_partitioner *get_partitioner() const {
		return &vertices->get_partitioner();
	}

	int get_num_threads() const {
		return worker_threads.size();
	}

	worker_thread *get_thread(int idx) const {
		return worker_threads[idx];
	}

	/**
	 * The following two methods keep track of the number of active vertices
	 * globally in the current iteration.
	 */

	// We have processed the specified number of vertices.
	void process_vertices(int num) {
		num_remaining_vertices_in_level.dec(num);
	}

	// Get the number of activated vertices that still haven't been
	// processed in the current level.
	size_t get_num_remaining_vertices() const {
		return num_remaining_vertices_in_level.get();
	}
};

#endif
