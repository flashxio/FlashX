#ifndef __WORKER_THREAD_H__
#define __WORKER_THREAD_H__

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

#include <pthread.h>

#include <vector>
#include <unordered_map>

#include "graph_engine.h"
#include "bitmap.h"
#include "scan_pointer.h"

class worker_thread;

/**
 * The queue for active vertices.
 */
class active_vertex_queue
{
public:
	virtual ~active_vertex_queue() {
	}

	virtual void init(const vertex_id_t buf[], size_t size, bool sorted) = 0;
	// This is the common case for iterations.
	virtual void init(worker_thread &) = 0;
	virtual int fetch(compute_vertex *vertices[], int num) = 0;
	virtual bool is_empty() = 0;
	virtual size_t get_num_vertices() = 0;

	void init(const std::vector<vertex_id_t> &vec, bool sorted) {
		init(vec.data(), vec.size(), sorted);
	}
};

/**
 * This vertex queue is sorted based on the vertex ID.
 */
class default_vertex_queue: public active_vertex_queue
{
	static const size_t VERTEX_BUF_SIZE = 10 * 1024;
	pthread_spinlock_t lock;
	// It contains the offset of the vertex in the local partition
	// instead of the real vertex Ids.
	std::vector<vertex_id_t> vertex_buf;
	std::unique_ptr<bitmap> active_bitmap;
	// The fetch index in the vertex buffer.
	scan_pointer buf_fetch_idx;
	// The fetech index in the active bitmap. It indicates the index of longs.
	scan_pointer bitmap_fetch_idx;
	graph_engine &graph;
	size_t num_active;
	int part_id;

	void fetch_from_map();
public:
	default_vertex_queue(graph_engine &_graph, int part_id,
			int node_id): buf_fetch_idx(0, true), bitmap_fetch_idx(0,
				true), graph(_graph) {
		pthread_spin_init(&lock, PTHREAD_PROCESS_PRIVATE);
		num_active = 0;
		this->part_id = part_id;
		this->active_bitmap = std::unique_ptr<bitmap>(new bitmap(
					_graph.get_partitioner()->get_part_size(
					part_id, _graph.get_num_vertices()), node_id));
	}

	virtual void init(const vertex_id_t buf[], size_t size, bool sorted);
	virtual void init(worker_thread &);
	virtual int fetch(compute_vertex *vertices[], int num);

	virtual bool is_empty() {
		return num_active == 0;
	}

	virtual size_t get_num_vertices() {
		return num_active;
	}
};

class customized_vertex_queue: public active_vertex_queue
{
	pthread_spinlock_t lock;
	std::vector<vertex_id_t> sorted_vertices;
	scan_pointer fetch_idx;
	vertex_scheduler::ptr scheduler;
	graph_engine &graph;
	int part_id;
public:
	customized_vertex_queue(graph_engine &_graph, vertex_scheduler::ptr scheduler,
			int part_id): fetch_idx(0, true), graph(_graph) {
		pthread_spin_init(&lock, PTHREAD_PROCESS_PRIVATE);
		this->scheduler = scheduler;
		this->part_id = part_id;
	}

	void init(const vertex_id_t buf[], size_t size, bool sorted) {
		pthread_spin_lock(&lock);
		bool forward = true;
		if (graph_conf.get_elevator_enabled())
			forward = graph.get_curr_level() % 2;
		this->fetch_idx = scan_pointer(size, forward);
		sorted_vertices.clear();
		sorted_vertices.insert(sorted_vertices.end(), buf, buf + size);
		if (!sorted)
			scheduler->schedule(sorted_vertices);
		pthread_spin_unlock(&lock);
	}

	void init(worker_thread &);

	int fetch(compute_vertex *vertices[], int num) {
		pthread_spin_lock(&lock);
		int num_fetches = min(num, fetch_idx.get_num_remaining());
		stack_array<vertex_id_t, 128> buf(num_fetches);
		if (num_fetches > 0) {
			size_t curr_loc = fetch_idx.get_curr_loc();
			size_t new_loc = fetch_idx.move(num_fetches);
			memcpy(buf.data(), sorted_vertices.data() + min(curr_loc, new_loc),
					num_fetches * sizeof(vertex_id_t));
		}
		pthread_spin_unlock(&lock);
		for (int i = 0; i < num_fetches; i++)
			vertices[i] = &graph.get_vertex(buf[i]);
		return num_fetches;
	}

	bool is_empty() {
		pthread_spin_lock(&lock);
		bool ret = fetch_idx.get_num_remaining() == 0;
		pthread_spin_unlock(&lock);
		return ret;
	}

	size_t get_num_vertices() {
		pthread_spin_lock(&lock);
		size_t num = fetch_idx.get_num_remaining();
		pthread_spin_unlock(&lock);
		return num;
	}
};

class vertex_compute;
class steal_state_t;
class message_processor;
class load_balancer;
class simple_index_reader;

class worker_thread: public thread
{
	int worker_id;
	file_io_factory::shared_ptr graph_factory;
	file_io_factory::shared_ptr index_factory;
	io_interface::ptr io;
	graph_engine *graph;
	compute_allocator *alloc;
	compute_allocator *part_alloc;
	vertex_program::ptr vprogram;
	std::shared_ptr<simple_index_reader> index_reader;

	// This buffers the I/O requests for adjacency lists.
	std::vector<io_request> adj_reqs;

	// When a thread process a vertex, the worker thread should keep
	// a vertex compute for the vertex. This is useful when a user-defined
	// compute vertex needs to reference its vertex compute.
	std::unordered_map<vertex_id_t, vertex_compute *> active_computes;
	// This references the vertex compute used/created by the current vertex
	// being processed.
	vertex_compute *curr_compute;

	/**
	 * A vertex is allowed to send messages to other vertices.
	 * The message passing scheme between vertices are implemented as follows:
	 * all non-empty vertices (with edges) have a message holder;
	 * vertices are partitioned and each worker thread is responsible for
	 * a certin number of vertices;
	 * when a vertex issues messages, the thread that processes the vertex will
	 * redirect them to the right threads;
	 * a thread that receives messages process them and place them in
	 * the right message holder.
	 */

	std::shared_ptr<slab_allocator> msg_alloc;
	std::unique_ptr<message_processor> msg_processor;
	std::unique_ptr<load_balancer> balancer;

	// This indicates the vertices that request the notification of the end
	// of an iteration.
	std::unique_ptr<bitmap> notify_vertices;
	// This is to collect vertices activated in the next level.
	std::unique_ptr<bitmap> next_activated_vertices;
	// This contains the vertices activated in the current level.
	std::unique_ptr<active_vertex_queue> curr_activated_vertices;
	vertex_scheduler::ptr scheduler;

	// Indicate that we need to start all vertices.
	bool start_all;
	std::vector<vertex_id_t> started_vertices;
	std::shared_ptr<vertex_filter> filter;
	vertex_initializer::ptr vinitializer;

	// The number of activated vertices processed in the current level.
	atomic_number<long> num_activated_vertices_in_level;
	// The number of vertices completed in the current level.
	atomic_number<long> num_completed_vertices_in_level;

	/**
	 * Get the number of vertices being processed in the current level.
	 */
	int get_num_vertices_processing() const {
		return num_activated_vertices_in_level.get()
			- num_completed_vertices_in_level.get();
	}
	int process_activated_vertices(int max);
public:
	worker_thread(graph_engine *graph, file_io_factory::shared_ptr graph_factory,
			file_io_factory::shared_ptr index_factory,
			vertex_program::ptr prog, int node_id, int worker_id,
			int num_threads, vertex_scheduler::ptr scheduler);

	~worker_thread();

	void init_messaging(const std::vector<worker_thread *> &threads);

	void run();
	void init();

	compute_allocator *get_part_compute_allocator() const {
		assert(part_alloc);
		return part_alloc;
	}

	/**
	 * When a vertex has been completed for the current iteration, this
	 * method is invoked to notify the worker thread, so that the worker
	 * thread can update its statistics on the number of completed vertices.
	 */
	void complete_vertex(const compute_vertex &v);

	size_t enter_next_level();

	void start_vertices(const std::vector<vertex_id_t> &vertices,
			vertex_initializer::ptr initializer) {
		this->vinitializer = initializer;
		started_vertices = vertices;
	}

	void start_all_vertices(vertex_initializer::ptr init) {
		start_all = true;
		this->vinitializer = init;
	}

	void start_vertices(std::shared_ptr<vertex_filter> filter) {
		this->filter = filter;
	}

	vertex_compute *get_vertex_compute(compute_vertex &v);
	vertex_compute *get_vertex_compute(vertex_id_t id) const {
		std::unordered_map<vertex_id_t, vertex_compute *>::const_iterator it
			= active_computes.find(id);
		assert(it != active_computes.end());
		return it->second;
	}

	/**
	 * Activate the vertex in its own partition for the next iteration.
	 */
	void activate_vertex(vertex_id_t id);
	void activate_vertex(local_vid_t id);

	void request_notify_iter_end(local_vid_t id) {
		notify_vertices->set(id.id);
	}

	int steal_activated_vertices(compute_vertex *vertices[], int num);
	void return_vertices(vertex_id_t ids[], int num);

	size_t get_num_local_vertices() const {
		return graph->get_partitioner()->get_part_size(worker_id,
					graph->get_num_vertices());
	}

	int get_worker_id() const {
		return worker_id;
	}

	vertex_program &get_vertex_program() {
		return *vprogram;
	}

	graph_engine &get_graph() {
		return *graph;
	}

	message_processor &get_msg_processor() {
		return *msg_processor;
	}

	simple_index_reader &get_index_reader() {
		return *index_reader;
	}

	void issue_io_request(io_request &req) {
		adj_reqs.push_back(req);
	}

	size_t get_activates() const {
		return curr_activated_vertices->get_num_vertices();
	}

	friend class load_balancer;
	friend class default_vertex_queue;
	friend class customized_vertex_queue;
};

#endif
