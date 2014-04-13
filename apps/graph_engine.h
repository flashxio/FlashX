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
	compute_vertex() {
		id = INVALID_VERTEX_ID;
	}

	compute_vertex(vertex_id_t id, const vertex_index *index) {
		this->id = id;
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
};

class compute_directed_vertex: public compute_vertex
{
	vsize_t num_in_edges;
public:
	compute_directed_vertex() {
		num_in_edges = 0;
	}

	compute_directed_vertex(vertex_id_t id,
			const vertex_index *index1): compute_vertex(id, index1) {
		assert(index1->get_graph_header().get_graph_type()
				== graph_type::DIRECTED);
		const directed_vertex_index *index
			= (const directed_vertex_index *) index1;
		num_in_edges = index->get_num_in_edges(id);
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
	compute_ts_vertex() {
	}

	compute_ts_vertex(vertex_id_t id, const vertex_index *index): compute_vertex(id, index) {
		assert(index->get_graph_header().get_graph_type()
				== graph_type::TS_DIRECTED
				|| index->get_graph_header().get_graph_type()
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

class worker_thread;

class graph_engine
{
	int vertex_header_size;
	graph_header header;
	graph_index *vertices;
	std::unique_ptr<ext_mem_vertex_interpreter> interpreter;
	vertex_scheduler *scheduler;

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

	trace_logger *logger;
	file_io_factory::shared_ptr factory;
	int max_processing_vertices;

	void cleanup() {
		if (logger) {
			logger->close();
			logger = NULL;
		}
	}

	void init_threads(vertex_program::ptr prog);
protected:
	graph_engine(int num_threads, int num_nodes, const std::string &graph_file,
			graph_index *index);
public:
	static graph_engine *create(int num_threads, int num_nodes,
			const std::string &graph_file, graph_index *index) {
		return new graph_engine(num_threads, num_nodes, graph_file, index);
	}

	static void destroy(graph_engine *graph) {
		graph->cleanup();
		delete graph;
	}

	~graph_engine();

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

	const in_mem_vertex_info get_vertex_info(vertex_id_t id) const {
		return vertices->get_vertex_info(id);
	}

	const vsize_t get_vertex_edges(vertex_id_t id) const {
		in_mem_vertex_info info = get_vertex_info(id);
		int vertex_header_size = get_vertex_header_size();
		return (info.get_ext_mem_size() - vertex_header_size) / sizeof(vertex_id_t);
	}

	void start(std::shared_ptr<vertex_filter> filter,
			vertex_program::ptr prog = vertex_program::ptr());
	void start(vertex_id_t ids[], int num,
			vertex_program::ptr prog = vertex_program::ptr());
	void start_all(vertex_program::ptr prog = vertex_program::ptr());

	/**
	 * The algorithm progresses to the next level.
	 * It returns true if no more work can progress.
	 */
	bool progress_next_level();

	/**
	 * Activate vertices that may be processed in the next level.
	 */
	void activate_vertices(vertex_id_t ids[], int num);
	void activate_vertices(edge_seq_iterator &it);

	void activate_vertex(vertex_id_t vertex) {
		activate_vertices(&vertex, 1);
	}

	/**
	 * Get vertices to be processed in the current level.
	 */
	int get_curr_activated_vertices(vertex_id_t vertices[], int num);
	size_t get_num_curr_activated_vertices() const;

	vertex_id_t get_max_vertex_id() const {
		return vertices->get_max_vertex_id();
	}

	vertex_id_t get_min_vertex_id() const {
		return vertices->get_min_vertex_id();
	}

	size_t get_num_vertices() const {
		return vertices->get_num_vertices();
	}

	void wait4complete();

	int get_num_threads() const {
		return worker_threads.size();
	}

	bool is_directed() const {
		return header.is_directed_graph();
	}

	trace_logger *get_logger() const {
		return logger;
	}

	/**
	 * Get the file id where the graph data is stored.
	 */
	int get_file_id() const {
		return factory->get_file_id();
	}

	void multicast_msg(vertex_id_t ids[], int num, const vertex_message &msg);
	void multicast_msg(edge_seq_iterator &it, vertex_message &msg);

	void send_msg(vertex_id_t dest, vertex_message &msg);

	ext_mem_vertex_interpreter &get_vertex_interpreter() const {
		return *interpreter;
	}

	const graph_partitioner *get_partitioner() const {
		return &vertices->get_partitioner();
	}

	worker_thread *get_thread(int idx) const {
		return worker_threads[idx];
	}

	const graph_header &get_graph_header() const {
		return header;
	}

	void set_vertex_scheduler(vertex_scheduler *scheduler);

	// We have processed the specified number of vertices.
	void process_vertices(int num) {
		num_remaining_vertices_in_level.dec(num);
	}

	// Get the number of activated vertices that still haven't been
	// processed in the current level.
	size_t get_num_remaining_vertices() const {
		return num_remaining_vertices_in_level.get();
	}

	int get_max_processing_vertices() const {
		if (max_processing_vertices > 0)
			return max_processing_vertices;
		else
			return graph_conf.get_max_processing_vertices();
	}

	void set_max_processing_vertices(int max) {
		max_processing_vertices = max;
	}

	int get_curr_level() const {
		return level.get();
	}

	int get_vertex_header_size() const {
		return vertex_header_size;
	}

	void preload_graph();
};

#endif
