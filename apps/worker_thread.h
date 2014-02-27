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

class sorted_vertex_queue
{
	pthread_spinlock_t lock;
	std::vector<vertex_id_t> sorted_vertices;
	size_t fetch_idx;
	vertex_scheduler *scheduler;
public:
	sorted_vertex_queue();

	void set_vertex_scheduler(vertex_scheduler *scheduler) {
		this->scheduler = scheduler;
	}

	void init(vertex_id_t buf[], int size, bool sorted) {
		pthread_spin_lock(&lock);
		fetch_idx = 0;
		sorted_vertices.clear();
		sorted_vertices.assign(buf, buf + size);
		if (!sorted)
			scheduler->schedule(sorted_vertices);
		pthread_spin_unlock(&lock);
	}

	void init(const std::vector<vertex_id_t> &vec, bool sorted) {
		pthread_spin_lock(&lock);
		fetch_idx = 0;
		sorted_vertices.clear();
		sorted_vertices.assign(vec.begin(), vec.end());
		if (!sorted)
			scheduler->schedule(sorted_vertices);
		pthread_spin_unlock(&lock);
	}

	void init(const bitmap &map, int part_id,
			const vertex_partitioner *partitioner);

	int fetch(vertex_id_t vertices[], int num) {
		pthread_spin_lock(&lock);
		int num_fetches = min(num, sorted_vertices.size() - fetch_idx);
		memcpy(vertices, sorted_vertices.data() + fetch_idx,
				num_fetches * sizeof(vertex_id_t));
		fetch_idx += num_fetches;
		pthread_spin_unlock(&lock);
		return num_fetches;
	}

	bool is_empty() {
		pthread_spin_lock(&lock);
		bool ret = sorted_vertices.size() - fetch_idx == 0;
		pthread_spin_unlock(&lock);
		return ret;
	}

	size_t get_num_vertices() {
		pthread_spin_lock(&lock);
		size_t num = sorted_vertices.size() - fetch_idx;
		pthread_spin_unlock(&lock);
		return num;
	}

	void init_vertices(graph_engine &graph) {
		pthread_spin_lock(&lock);
		BOOST_FOREACH(vertex_id_t id, sorted_vertices) {
			compute_vertex &v = graph.get_vertex(id);
			v.init();
		}
		pthread_spin_unlock(&lock);
	}
};

class vertex_compute;

class worker_thread: public thread
{
	int worker_id;
	file_io_factory *factory;
	io_interface *io;
	graph_engine *graph;
	compute_allocator *alloc;
	compute_allocator *part_alloc;

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

	// The queue of messages sent from other threads.
	msg_queue msg_q;
	slab_allocator *msg_alloc;
	// The message senders to send messages to all other threads.
	// There are n senders, n is the total number of threads used by
	// the graph engine.
	std::vector<simple_msg_sender *> msg_senders;
	std::vector<multicast_msg_sender *> multicast_senders;
	std::vector<multicast_msg_sender *> activate_senders;
	// This is to collect vertices activated in the next level.
	bitmap next_activated_vertices;
	// This contains the vertices activated in the current level.
	sorted_vertex_queue curr_activated_vertices;
	// The thread where we should steal activated vertices from.
	int steal_thread_id;

	// Indicate that we need to start all vertices.
	bool start_all;

	// The number of activated vertices processed in the current level.
	atomic_number<long> num_activated_vertices_in_level;
	// The number of vertices completed in the current level.
	atomic_number<long> num_completed_vertices_in_level;
public:
	worker_thread(graph_engine *graph, file_io_factory *factory, int node_id,
			int worker_id, int num_threads);

	~worker_thread();

	void init_messaging(const std::vector<worker_thread *> &threads);

	void run();
	void init() {
		io = factory->create_io(this);
		io->init();

		// If a user wants to start all vertices.
		if (start_all) {
			std::vector<vertex_id_t> local_ids;
			graph->get_partitioner()->get_all_vertices_in_part(worker_id,
					graph->get_num_vertices(), local_ids);
			assert(curr_activated_vertices.is_empty());
			curr_activated_vertices.init(local_ids, false);
		}
		curr_activated_vertices.init_vertices(*graph);
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

	int process_activated_vertices(int max);

	void complete_vertex(const compute_vertex &v) {
		num_completed_vertices_in_level.inc(1);
	}

	void flush_msgs() {
		for (size_t i = 0; i < msg_senders.size(); i++)
			msg_senders[i]->flush();
		for (size_t i = 0; i < multicast_senders.size(); i++)
			multicast_senders[i]->flush();
		for (size_t i = 0; i < activate_senders.size(); i++)
			activate_senders[i]->flush();
	}

	void process_msgs();
	void process_msg(message &msg);
	void process_multicast_msg(multicast_message &mmsg);
	int enter_next_level();

	int steal_activated_vertices(vertex_id_t buf[], int num);

	void start_vertices(const std::vector<vertex_id_t> &vertices) {
		assert(curr_activated_vertices.is_empty());
		curr_activated_vertices.init(vertices, false);
	}

	void start_all_vertices() {
		start_all = true;
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

	vertex_compute *get_curr_vertex_compute() {
		if (curr_compute == NULL)
			curr_compute = (vertex_compute *) alloc->alloc();
		return curr_compute;
	}

	void set_curr_vertex_compute(vertex_compute *compute) {
		assert(this->curr_compute == NULL);
		curr_compute = compute;
	}

	void reset_curr_vertex_compute() {
		curr_compute = NULL;
	}
};

#endif
