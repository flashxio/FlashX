#ifndef __WORKER_THREAD_H__
#define __WORKER_THREAD_H__

/**
 * Copyright 2014 Da Zheng
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

#include <pthread.h>

#include <vector>

#include "graph_engine.h"
#include "bitmap.h"
#include "scan_pointer.h"

class sorted_vertex_queue
{
	pthread_spinlock_t lock;
	std::vector<compute_vertex *> sorted_vertices;
	scan_pointer fetch_idx;
	vertex_scheduler *scheduler;
	graph_engine &graph;
public:
	sorted_vertex_queue(graph_engine &graph);

	void set_vertex_scheduler(vertex_scheduler *scheduler) {
		this->scheduler = scheduler;
	}

	void init(const vertex_id_t buf[], size_t size, bool sorted) {
		pthread_spin_lock(&lock);
		bool forward = true;
		if (graph_conf.get_elevator_enabled())
			forward = graph.get_curr_level() % 2;
		this->fetch_idx = scan_pointer(size, forward);
		sorted_vertices.clear();
		sorted_vertices.resize(size);
		for (size_t i = 0; i < size; i++)
			sorted_vertices[i] = &graph.get_vertex(buf[i]);
		if (!sorted)
			scheduler->schedule(sorted_vertices);
		pthread_spin_unlock(&lock);
	}

	void init(const std::vector<vertex_id_t> &vec, bool sorted) {
		init(vec.data(), vec.size(), sorted);
	}

	void init(const bitmap &map, int part_id,
			const graph_partitioner *partitioner);

	int fetch(compute_vertex *vertices[], int num) {
		pthread_spin_lock(&lock);
		int num_fetches = min(num, fetch_idx.get_num_remaining());
		if (num_fetches > 0) {
			size_t curr_loc = fetch_idx.get_curr_loc();
			size_t new_loc = fetch_idx.move(num_fetches);
			memcpy(vertices, sorted_vertices.data() + min(curr_loc, new_loc),
					num_fetches * sizeof(compute_vertex *));
		}
		pthread_spin_unlock(&lock);
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

class worker_thread: public thread
{
	int worker_id;
	file_io_factory::shared_ptr factory;
	io_interface *io;
	graph_engine *graph;
	compute_allocator *alloc;
	compute_allocator *part_alloc;
	vertex_program::ptr vprogram;

	// When a thread process a vertex, the worker thread should point to
	// a vertex compute. This is useful when a user-defined compute vertex
	// needs to reference its vertex compute.
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

	slab_allocator *msg_alloc;
	// The message senders to send messages to all other threads.
	// There are n senders, n is the total number of threads used by
	// the graph engine.
	std::vector<simple_msg_sender *> msg_senders;
	std::vector<multicast_msg_sender *> multicast_senders;
	std::vector<multicast_msg_sender *> activate_senders;
	std::unique_ptr<message_processor> msg_processor;

	std::unique_ptr<load_balancer> balancer;

	// This is to collect vertices activated in the next level.
	bitmap next_activated_vertices;
	// This contains the vertices activated in the current level.
	sorted_vertex_queue curr_activated_vertices;

	// Indicate that we need to start all vertices.
	bool start_all;
	std::shared_ptr<vertex_filter> filter;

	// The number of activated vertices processed in the current level.
	atomic_number<long> num_activated_vertices_in_level;
	// The number of vertices completed in the current level.
	atomic_number<long> num_completed_vertices_in_level;
public:
	worker_thread(graph_engine *graph, file_io_factory::shared_ptr factory,
			vertex_program::ptr prog, int node_id, int worker_id, int num_threads);

	~worker_thread();

	void init_messaging(const std::vector<worker_thread *> &threads);

	void run();
	void init() {
		io = factory->create_io(this);
		io->init();

		std::vector<vertex_id_t> local_ids;
		graph->get_partitioner()->get_all_vertices_in_part(worker_id,
				graph->get_num_vertices(), local_ids);

		std::vector<vertex_id_t> kept_ids;
		BOOST_FOREACH(vertex_id_t id, local_ids) {
			compute_vertex &v = graph->get_vertex(id);
			if (filter && filter->keep(v))
				kept_ids.push_back(id);
		}
		if (!kept_ids.empty()) {
			assert(curr_activated_vertices.is_empty());
			// Although we don't process the filtered vertices, we treat
			// them as if they were processed.
			graph->process_vertices(local_ids.size() - kept_ids.size());
			curr_activated_vertices.init(kept_ids, false);
			printf("worker %d has %ld vertices and activates %ld of them\n",
					worker_id, local_ids.size(), kept_ids.size());
		}
		// If a user wants to start all vertices.
		else if (start_all) {
			assert(curr_activated_vertices.is_empty());
			curr_activated_vertices.init(local_ids, false);
		}
		else if (filter) {
			// Although we don't process the filtered vertices, we treat
			// them as if they were processed.
			graph->process_vertices(local_ids.size());
		}
	}

	compute_allocator *get_part_compute_allocator() const {
		return part_alloc;
	}

	multicast_msg_sender *get_activate_sender(int thread_id) const {
		return activate_senders[thread_id];
	}

	multicast_msg_sender *get_multicast_sender(int thread_id) const {
		return multicast_senders[thread_id];
	}

	simple_msg_sender *get_msg_sender(int thread_id) const {
		return msg_senders[thread_id];
	}

	/**
	 * When a vertex has been completed for the current iteration, this
	 * method is invoked to notify the worker thread, so that the worker
	 * thread can update its statistics on the number of completed vertices.
	 */
	void complete_vertex(const compute_vertex &v);

	void flush_msgs() {
		for (size_t i = 0; i < msg_senders.size(); i++)
			msg_senders[i]->flush();
		for (size_t i = 0; i < multicast_senders.size(); i++)
			multicast_senders[i]->flush();
		for (size_t i = 0; i < activate_senders.size(); i++) {
			activate_senders[i]->flush();
			vertex_message msg(sizeof(vertex_message), true);
			activate_senders[i]->init(msg);
		}
	}

	int process_activated_vertices(int max);
	int enter_next_level();

	void start_vertices(const std::vector<vertex_id_t> &vertices) {
		assert(curr_activated_vertices.is_empty());
		curr_activated_vertices.init(vertices, false);
	}

	void start_all_vertices() {
		start_all = true;
	}

	void start_vertices(std::shared_ptr<vertex_filter> filter) {
		this->filter = filter;
	}

	void set_vertex_scheduler(vertex_scheduler *scheduler) {
		curr_activated_vertices.set_vertex_scheduler(scheduler);
	}

	/**
	 * Get the number of vertices being processed in the current level.
	 */
	int get_num_vertices_processing() const {
		return num_activated_vertices_in_level.get()
			- num_completed_vertices_in_level.get();
	}

	size_t get_num_activated_on_others() {
		return graph->get_num_remaining_vertices();
	}

	vertex_compute *get_curr_vertex_compute() const {
		return curr_compute;
	}

	vertex_compute *create_vertex_compute(compute_vertex *v);

	void set_curr_vertex_compute(vertex_compute *compute) {
		assert(this->curr_compute == NULL);
		curr_compute = compute;
	}

	void reset_curr_vertex_compute() {
		curr_compute = NULL;
	}

	/**
	 * Activate a vertex for the next iteration.
	 */
	void activate_vertex(vertex_id_t id);

	int steal_activated_vertices(compute_vertex *vertices[], int num);
	void return_vertices(vertex_id_t ids[], int num);

	size_t get_num_local_vertices() const {
		return next_activated_vertices.get_num_bits();
	}

	int get_worker_id() const {
		return worker_id;
	}

	slab_allocator *get_msg_allocator() const {
		return msg_alloc;
	}

	vertex_program &get_vertex_program() {
		return *vprogram;
	}

	graph_engine &get_graph() {
		return *graph;
	}

	friend class load_balancer;
};

#endif
