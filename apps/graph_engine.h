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

/**
 * The size of a message buffer used to pass vertex messages to other threads.
 */
const int GRAPH_MSG_BUF_SIZE = PAGE_SIZE * 4;

class graph_engine;
class vertex_request;

class compute_vertex: public in_mem_vertex_info
{
public:
	compute_vertex() {
	}

	compute_vertex(vertex_id_t id, const vertex_index *index): in_mem_vertex_info(
			id, index) {
	}

	virtual compute_allocator *create_part_compute_allocator(
			graph_engine *graph, thread *t) {
		return NULL;
	}

	virtual void init() {
	}

	/**
	 * This is a pre-run before users get any information of adjacency list
	 * of vertices.
	 * If it returns true, the graph engine will fetch the adjacency list
	 * of itself. So by default, the graph engine will fetch its adjacency
	 * list.
	 */
	virtual bool run(graph_engine &graph) {
		return true;
	}

	/**
	 * Run user's code when the adjacency list of the vertex is read
	 * from disks.
	 * It returns true if the vertex has completed the iteration.
	 */
	virtual bool run(graph_engine &graph, const page_vertex *vertex) = 0;

	/**
	 * Run user's code when the adjacency lists of its neighbors are read
	 * from disks.
	 * It returns true if the vertex has completed the iteration.
	 */
	virtual bool run_on_neighbors(graph_engine &graph,
			const page_vertex *vertices[], int num) = 0;

	/**
	 * Run user's code when the vertex receives messages from other.
	 */
	virtual void run_on_messages(graph_engine &,
			const vertex_message *msgs[], int num) = 0;
};

class ts_vertex_request;

class ts_compute_vertex: public compute_vertex
{
public:
	ts_compute_vertex() {
	}

	ts_compute_vertex(vertex_id_t id, const vertex_index *index): compute_vertex(
			id, index) {
	}

	virtual compute_allocator *create_part_compute_allocator(
			graph_engine *graph, thread *t);

	virtual bool run_on_neighbors(graph_engine &graph,
			const page_vertex *vertices[], int num);
	virtual bool run_on_neighbors(graph_engine &graph,
			const TS_page_vertex *vertices[], int num) = 0;
};

class vertex_scheduler
{
public:
	virtual void schedule(std::vector<vertex_id_t> &vertices) = 0;
};

class default_vertex_scheduler: public vertex_scheduler
{
public:
	void schedule(std::vector<vertex_id_t> &vertices) {
		std::sort(vertices.begin(), vertices.end());
	}
};
extern default_vertex_scheduler default_scheduler;

class worker_thread;

class graph_engine
{
	graph_header header;
	graph_index *vertices;
	ext_mem_vertex_interpreter *interpreter;
	vertex_partitioner *partitioner;
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

	edge_type required_neighbor_type;
	trace_logger *logger;
	file_io_factory *factory;
	int max_processing_vertices;

	void cleanup() {
		if (logger) {
			logger->close();
			logger = NULL;
		}
	}

	/**
	 * This method returns the message sender of the current thread that
	 * sends messages to the thread with the specified thread id.
	 */
	simple_msg_sender *get_msg_sender(int thread_id) const;
	multicast_msg_sender *get_multicast_sender(int thread_id) const;
	multicast_msg_sender *get_activate_sender(int thread_id) const;
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

	void start(vertex_id_t ids[], int num);
	void start_all();

	void set_required_neighbor_type(edge_type type) {
		required_neighbor_type = type;
	}

	edge_type get_required_neighbor_type() const {
		return required_neighbor_type;
	}

	/**
	 * The algorithm progresses to the next level.
	 * It returns true if no more work can progress.
	 */
	bool progress_next_level();

	/**
	 * Activate vertices that may be processed in the next level.
	 */
	void activate_vertices(vertex_id_t ids[], int num) {
		for (int i = 0; i < num; i++) {
			int part_id = partitioner->map(ids[i]);
			multicast_msg_sender *sender = get_activate_sender(part_id);
			bool ret = false;
			if (sender->has_msg()) {
				ret = sender->add_dest(ids[i]);
			}
			// If we can't add a destination vertex to the multicast msg,
			// or there isn't a msg in the sender.
			if (!ret) {
				vertex_message msg(sizeof(vertex_message), true);
				sender->init(msg);
				ret = sender->add_dest(ids[i]);
				assert(ret);
			}
		}
	}

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

	template<class T>
	void multicast_msg(vertex_id_t ids[], int num, const T &msg) {
		for (int i = 0; i < num; i++) {
			int part_id = partitioner->map(ids[i]);
			multicast_msg_sender *sender = get_multicast_sender(part_id);
			bool ret = false;
			if (sender->has_msg()) {
				ret = sender->add_dest(ids[i]);
			}
			// If we can't add a destination vertex to the multicast msg,
			// or there isn't a msg in the sender.
			if (!ret) {
				int retries = 0;
				do {
					retries++;
					sender->init(msg);
					ret = sender->add_dest(ids[i]);
				} while (!ret);
				// We shouldn't try it more than twice.
				assert(retries <= 2);
			}
		}
		// Now we have multicast the message, we need to notify all senders
		// of the end of multicast.
		for (int i = 0; i < this->get_num_threads(); i++) {
			multicast_msg_sender *sender = get_multicast_sender(i);
			// We only send the multicast on the sender that has received
			// the multicast message.
			if (sender->has_msg())
				sender->end_multicast();
		}
	}

	template<class T>
	void send_msg(vertex_id_t dest, T &msg) {
		vertex_id_t id = dest;
		simple_msg_sender *sender = get_msg_sender(partitioner->map(id));
		msg.set_dest(dest);
		sender->send_cached(msg);
	}

	/**
	 * This allows a vertex to request other vertices in the graph.
	 * @vertex: the vertex that sends the requests.
	 * @ids: the Ids of vertices requested by @vertex.
	 */
	void request_vertices(compute_vertex &vertex, vertex_id_t ids[], int num);
	/**
	 * This allows a vertex to request other vertices in the graph as above.
	 * The difference is that this interface allows a vertex to request
	 * part of other vertices.
	 * @vertex: the vertex that sends the requests.
	 * @reqs: defines part of vertices requested by @vertex.
	 */
	void request_partial_vertices(compute_vertex &vertex,
			vertex_request *reqs[], int num);

	ext_mem_vertex_interpreter *get_vertex_interpreter() const {
		return interpreter;
	}

	compute_allocator *create_part_compute_allocator(thread *t) {
		// We only need to get a vertex that exists.
		int min_id = vertices->get_min_vertex_id();
		return vertices->get_vertex(min_id).create_part_compute_allocator(
				this, t);
	}

	void destroy_part_compute_allocator(compute_allocator *alloc) {
		delete alloc;
	}

	const vertex_partitioner *get_partitioner() const {
		return partitioner;
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
};

#endif
