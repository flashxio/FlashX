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

#include "vertex.h"
#include "vertex_index.h"
#include "trace_logger.h"
#include "messaging.h"

/**
 * The size of a message buffer used to pass vertex messages to other threads.
 */
const int GRAPH_MSG_BUF_SIZE = PAGE_SIZE * 4;

class graph_engine;

class vertex_message
{
	vertex_id_t dest;
	int size;
public:
	static vertex_message *deserialize(char *buf, int size) {
		vertex_message *msg = (vertex_message *) buf;
		assert(msg->size <= size);
		return msg;
	}

	vertex_message(vertex_id_t dest) {
		this->dest = dest;
		this->size = sizeof(vertex_message);
	}

	vertex_message(vertex_id_t dest, int size) {
		this->dest = dest;
		this->size = size;
	}

	vertex_id_t get_dest() const {
		return dest;
	}

	int get_serialized_size() const {
		return size;
	}

	bool is_empty() const {
		return (size_t) size == sizeof(vertex_message);
	}

	int serialize(char *buf, int size) const {
		assert(this->size <= size);
		memcpy(buf, this, this->size);
		return this->size;
	}
};

class compute_vertex: public in_mem_vertex_info
{
public:
	compute_vertex(vertex_id_t id, off_t off, int size): in_mem_vertex_info(
			id, off, size) {
	}

	virtual compute_allocator *create_part_compute_allocator(
			graph_engine *graph, thread *t) {
		return NULL;
	}

	virtual bool has_required_vertices() const {
		return false;
	}

	/**
	 * This translate the required vertex to an I/O request to the file.
	 */
	virtual request_range get_next_request(graph_engine *graph);

	virtual vertex_id_t get_next_required_vertex() {
		assert(0);
		return -1;
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
	vertex_id_t get_next_required_vertex() {
		return -1;
	}
	virtual bool has_required_vertices() const {
		return has_required_ts_vertices();
	}
public:
	ts_compute_vertex(vertex_id_t id, off_t off, int size): compute_vertex(
			id, off, size) {
	}

	virtual compute_allocator *create_part_compute_allocator(
			graph_engine *graph, thread *t);

	virtual bool run_on_neighbors(graph_engine &graph,
			const page_vertex *vertices[], int num);

	virtual request_range get_next_request(graph_engine *graph);

	virtual void get_next_required_ts_vertex(ts_vertex_request &) = 0;
	virtual bool has_required_ts_vertices() const = 0;
	virtual bool run_on_neighbors(graph_engine &graph,
			const TS_page_vertex *vertices[], int num) = 0;
};

class graph_index
{
public:
	virtual compute_vertex &get_vertex(vertex_id_t id) = 0;

	virtual int is_empty_vertex(vertex_id_t) const = 0;

	virtual vertex_id_t get_max_vertex_id() const = 0;

	virtual vertex_id_t get_min_vertex_id() const = 0;

	virtual size_t get_num_vertices() const = 0;

	virtual size_t get_all_vertices(std::vector<vertex_id_t> &vec) const = 0;
};

template<class vertex_type>
class graph_index_impl: public graph_index
{
	// This contains the vertices with edges.
	std::vector<vertex_type> vertices;
	int min_vertex_size;
	
	graph_index_impl(const std::string &index_file, int min_vertex_size) {
		this->min_vertex_size = min_vertex_size;
		vertex_index *indices = vertex_index::load(index_file);
		size_t num_vertices = indices->get_num_vertices();
		size_t num_non_empty = 0;
		vertices.resize(num_vertices);
		for (size_t i = 0; i < num_vertices; i++) {
			off_t off = indices->get_vertex_off(i);
			int size = indices->get_vertex_size(i);
			vertices[i] = vertex_type(i, off, size);
			if (size > min_vertex_size) {
				num_non_empty++;
			}
		}
		printf("There are %ld vertices and %ld non-empty, vertex array capacity: %ld\n",
				num_vertices, num_non_empty, vertices.capacity());
		vertex_index::destroy(indices);
	}
public:
	static graph_index *create(const std::string &index_file,
			int min_vertex_size) {
		return new graph_index_impl<vertex_type>(index_file, min_vertex_size);
	}

	virtual compute_vertex &get_vertex(vertex_id_t id) {
		return vertices[id];
	}

	virtual size_t get_num_vertices() const {
		return vertices.size();
	}

	virtual size_t get_all_vertices(std::vector<vertex_id_t> &vec) const {
		vec.resize(vertices.size());
		for (size_t i = 0; i < vertices.size(); i++)
			vec[i] = vertices[i].get_id();
		return vec.size();
	}

	virtual vertex_id_t get_max_vertex_id() const {
		return vertices.back().get_id();
	}

	virtual vertex_id_t get_min_vertex_id() const {
		return vertices.front().get_id();
	}

	virtual int is_empty_vertex(vertex_id_t id) const {
		return vertices[id].get_ext_mem_size() == min_vertex_size;
	}
};

/**
 * This defines the interface of interpreting vertices in the external memory.
 */
class ext_mem_vertex_interpreter
{
public:
	/**
	 * Interpret the data in the page byte array, and construct the page vertex
	 * in the buffer.
	 */
	virtual page_vertex *interpret(page_byte_array &, char *buf,
			int size) const = 0;
	virtual page_vertex *interpret_part(const page_vertex *header,
			page_byte_array &, char *buf, int size) const = 0;
	/**
	 * The size of the vertex object.
	 */
	virtual int get_vertex_size() const = 0;
};

class ext_mem_directed_vertex_interpreter: public ext_mem_vertex_interpreter
{
public:
	virtual page_vertex *interpret(page_byte_array &array, char *buf,
			int size) const {
		assert(size >= (int) sizeof(page_directed_vertex));
		return new (buf) page_directed_vertex(array);
	}

	virtual page_vertex *interpret_part(const page_vertex *header,
			page_byte_array &, char *buf, int size) const {
		assert(0);
		return NULL;
	}

	virtual int get_vertex_size() const {
		return sizeof(page_directed_vertex);
	}
};

class ext_mem_undirected_vertex_interpreter: public ext_mem_vertex_interpreter
{
public:
	virtual page_vertex *interpret(page_byte_array &array, char *buf,
			int size) const {
		assert(size >= (int) sizeof(page_undirected_vertex));
		return new (buf) page_undirected_vertex(array);
	}

	virtual page_vertex *interpret_part(const page_vertex *header,
			page_byte_array &, char *buf, int size) const {
		assert(0);
		return NULL;
	}

	virtual int get_vertex_size() const {
		return sizeof(page_undirected_vertex);
	}
};

class ts_ext_mem_vertex_interpreter: public ext_mem_vertex_interpreter
{
	int num_timestamps;
public:
	ts_ext_mem_vertex_interpreter(int num_timestamps) {
		this->num_timestamps = num_timestamps;
	}

	virtual page_vertex *interpret(page_byte_array &array, char *buf,
			int size) const {
		assert(size >= TS_page_directed_vertex::get_size(num_timestamps));
		return TS_page_directed_vertex::create(array, buf, size);
	}

	virtual page_vertex *interpret_part(const page_vertex *header,
			page_byte_array &array, char *buf, int size) const {
		assert(size >= TS_page_directed_vertex::get_size(num_timestamps));
		return TS_page_directed_vertex::create(
				(const TS_page_directed_vertex *) header, array, buf, size);
	}

	virtual int get_vertex_size() const {
		return TS_page_directed_vertex::get_size(num_timestamps);
	}
};

class vertex_partitioner
{
	int num_parts_log;
	vertex_id_t mask;
public:
	vertex_partitioner(int num_parts) {
		this->num_parts_log = log2(num_parts);
		assert((1 << num_parts_log) == num_parts);
		mask = (1 << num_parts_log) - 1;
	}

	int map(vertex_id_t id) const {
		return id & mask;
	}

	void map2loc(vertex_id_t id, int &part_id, off_t &off) const {
		part_id = id & mask;
		off = id >> num_parts_log;
	}

	void loc2map(int part_id, off_t off, vertex_id_t &id) const {
		id = (off << num_parts_log) + part_id;
	}
};

class worker_thread;

class graph_engine
{
	graph_index *vertices;
	ext_mem_vertex_interpreter *interpreter;
	vertex_partitioner *partitioner;

	atomic_integer level;
	volatile bool is_complete;

	// These are used for switching queues.
	pthread_mutex_t lock;
	pthread_barrier_t barrier1;
	pthread_barrier_t barrier2;

	thread *first_thread;
	std::vector<worker_thread *> worker_threads;

	bool directed;
	edge_type required_neighbor_type;

	trace_logger *logger;

	int file_id;

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
protected:
	graph_engine(int num_threads, int num_nodes, const std::string &graph_file,
			graph_index *index, ext_mem_vertex_interpreter *interpreter,
			bool directed);
public:
	static graph_engine *create(int num_threads, int num_nodes,
			const std::string &graph_file, graph_index *index,
			ext_mem_vertex_interpreter *interpreter, bool directed) {
		return new graph_engine(num_threads, num_nodes, graph_file, index,
				interpreter, directed);
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
	void activate_vertices(vertex_id_t vertices[], int num) {
		for (int i = 0; i < num; i++) {
			vertex_message msg(vertices[i]);
			send_msg(msg);
		}
	}

	void activate_vertex(vertex_id_t vertex) {
		vertex_message msg(vertex);
		send_msg(msg);
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

	bool is_empty_vertex(vertex_id_t id) const {
		return vertices->is_empty_vertex(id);
	}

	void wait4complete();

	int get_num_threads() const {
		return worker_threads.size();
	}

	bool is_directed() const {
		return directed;
	}

	trace_logger *get_logger() const {
		return logger;
	}

	/**
	 * Get the file id where the graph data is stored.
	 */
	int get_file_id() const {
		return file_id;
	}

	template<class T>
	void send_msg(const T &msg) {
		vertex_id_t id = msg.get_dest();
		simple_msg_sender *sender = get_msg_sender(partitioner->map(id));
		sender->send_cached(msg);
	}

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
};

/**
 * This class contains the request of a time-series vertex
 * from the user application.
 */
class ts_vertex_request
{
	vertex_id_t id;
	timestamp_pair range;
	edge_type type;
	bool require_all;
	graph_engine *graph;
public:
	ts_vertex_request(graph_engine *graph) {
		id = 0;
		range = timestamp_pair(INT_MAX, INT_MIN);
		type = edge_type::BOTH_EDGES;
		require_all = false;
		this->graph = graph;
	}

	void set_require_all(bool require_all) {
		this->require_all = require_all;
	}

	void set_vertex(vertex_id_t id);

	void add_timestamp(int timestamp) {
		if (!require_all) {
			if (range.second < timestamp)
				range.second = timestamp + 1;
			if (range.first > timestamp)
				range.first = timestamp;
		}
	}

	void set_edge_type(edge_type type) {
		this->type = type;
	}

	void clear() {
		id = 0;
		range = timestamp_pair(INT_MAX, INT_MIN);
		type = edge_type::BOTH_EDGES;
		require_all = false;
	}

	vertex_id_t get_id() const {
		return id;
	}

	const timestamp_pair &get_range() const {
		return range;
	}

	edge_type get_edge_type() const {
		return type;
	}

	bool is_require_all() const {
		return require_all;
	}
};

#endif
