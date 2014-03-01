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

#include "worker_thread.h"
#include "graph_engine.h"
#include "bitmap.h"
#include "vertex_compute.h"

const int MAX_STOLEN_VERTICES = 1024;

sorted_vertex_queue::sorted_vertex_queue()
{
	fetch_idx = 0;
	pthread_spin_init(&lock, PTHREAD_PROCESS_PRIVATE);
	this->scheduler = &default_scheduler;
}

void sorted_vertex_queue::init(const bitmap &map, int part_id,
		const vertex_partitioner *partitioner)
{
	pthread_spin_lock(&lock);
	fetch_idx = 0;
	sorted_vertices.clear();
	map.get_set_bits(sorted_vertices);
	// the bitmap only contains the locations of vertices in the bitmap.
	// We have to translate them back to vertex ids.
	for (size_t i = 0; i < sorted_vertices.size(); i++) {
		vertex_id_t id;
		partitioner->loc2map(part_id, sorted_vertices[i], id);
		sorted_vertices[i] = id;
	}

	if (scheduler != &default_scheduler)
		scheduler->schedule(sorted_vertices);
	pthread_spin_unlock(&lock);
}

worker_thread::worker_thread(graph_engine *graph, file_io_factory *factory,
		int node_id, int worker_id, int num_threads): thread("worker_thread",
			node_id), msg_q(get_node_id(), "graph_msg_queue", 16, INT_MAX),
		next_activated_vertices((size_t) ceil(((double) graph->get_max_vertex_id()
						+ 1) / num_threads))
{
	start_all = false;
	this->worker_id = worker_id;
	this->graph = graph;
	this->io = NULL;
	this->factory = factory;
	this->curr_compute = NULL;
	switch(graph->get_graph_header().get_graph_type()) {
		case graph_type::DIRECTED:
			alloc = new vertex_compute_allocator<vertex_compute>(graph, this);
			part_alloc = NULL;
			break;
		case graph_type::TS_DIRECTED:
			alloc = new vertex_compute_allocator<ts_vertex_compute>(graph, this);
			part_alloc = new vertex_compute_allocator<part_ts_vertex_compute>(
					graph, this);
			break;
		default:
			assert(0);

	}
}

worker_thread::~worker_thread()
{
	delete alloc;
	delete part_alloc;
	for (unsigned i = 0; i < msg_senders.size(); i++)
		simple_msg_sender::destroy(msg_senders[i]);
	for (unsigned i = 0; i < multicast_senders.size(); i++)
		multicast_msg_sender::destroy(multicast_senders[i]);
	for (unsigned i = 0; i < activate_senders.size(); i++)
		multicast_msg_sender::destroy(activate_senders[i]);
	delete msg_alloc;
	factory->destroy_io(io);
}

void worker_thread::init_messaging(const std::vector<worker_thread *> &threads)
{
	steal_thread_id = (worker_id + 1) % threads.size();
	// We increase the allocator by 1M each time.
	// It shouldn't need to allocate much memory.
	msg_alloc = new slab_allocator("graph-message-allocator",
			GRAPH_MSG_BUF_SIZE, 1024 * 1024, INT_MAX, get_node_id());

	int num_self = 0;
	for (unsigned i = 0; i < threads.size(); i++) {
		if (threads[i] == this)
			num_self++;
		msg_senders.push_back(simple_msg_sender::create(get_node_id(),
					msg_alloc, &threads[i]->msg_q));
		multicast_senders.push_back(multicast_msg_sender::create(msg_alloc,
					&threads[i]->msg_q));
		activate_senders.push_back(multicast_msg_sender::create(msg_alloc,
					&threads[i]->msg_q));
	}
	assert(num_self == 1);
}

/**
 * This steals vertices from other threads. It tries to steal more vertices
 * than it can process, and the remaining vertices will be placed in its
 * own activated vertex queue.
 */
int worker_thread::steal_activated_vertices(vertex_id_t vertex_buf[], int buf_size)
{
	if (steal_thread_id == this->worker_id)
		steal_thread_id = (steal_thread_id + 1) % graph->get_num_threads();
	int num_tries = 0;
	vertex_id_t *steal_buf = new vertex_id_t[MAX_STOLEN_VERTICES];
	int num;
	do {
		worker_thread *t = graph->get_thread(steal_thread_id);
		num_tries++;

		// We want to steal as much as possible, but we don't want
		// to overloaded by the stolen vertices.
		size_t num_steal = max(1,
				t->curr_activated_vertices.get_num_vertices() / graph->get_num_threads());
		num = t->curr_activated_vertices.fetch(steal_buf,
				min(MAX_STOLEN_VERTICES, num_steal));

		// If we can't steal vertices from the thread, we should move
		// to the next thread.
		if (num == 0)
			steal_thread_id = (steal_thread_id + 1) % graph->get_num_threads();
		// If we have tried to steal vertices from all threads.
	} while (num == 0 && num_tries < graph->get_num_threads());

	int ret = min(buf_size, num);
	memcpy(vertex_buf, steal_buf, sizeof(vertex_buf[0]) * ret);
	// We stole more vertices than we can process this time.
	// The vertices stolen from another thread will also be placed in
	// the queue for currently activated vertices.
	if (num - ret > 0)
		curr_activated_vertices.init(steal_buf + ret, num - ret, true);
	delete [] steal_buf;
	return ret;
}

/**
 * This is to process the activated vertices in the current iteration.
 */
int worker_thread::process_activated_vertices(int max)
{
	if (max <= 0)
		return 0;

	vertex_id_t vertex_buf[max];
	stack_array<io_request> reqs(max);
	int num = curr_activated_vertices.fetch(vertex_buf, max);
	if (num == 0)
		num = steal_activated_vertices(vertex_buf, max);
	if (num > 0) {
		num_activated_vertices_in_level.inc(num);
		graph->process_vertices(num);
	}

	int num_completed = 0;
	int num_to_process = 0;
	for (int i = 0; i < num; i++) {
		compute_vertex &info = graph->get_vertex(vertex_buf[i]);
		// We execute the pre-run to determine if the vertex has completed
		// in the current iteration.
		if (!info.run(*graph)) {
			data_loc_t loc(io->get_file_id(), info.get_ext_mem_off());
			reqs[num_to_process++] = io_request(alloc->alloc(), loc,
					// TODO I might need to set the node id.
					info.get_ext_mem_size(), READ, io, -1);
		}
		else
			num_completed++;
	}
	if (num_completed > 0)
		num_completed_vertices_in_level.inc(num_completed);
	if (graph->get_logger())
		graph->get_logger()->log(reqs.data(), num_to_process);
	io->access(reqs.data(), num_to_process);
	return num;
}

int worker_thread::enter_next_level()
{
	// We have to make sure all messages sent by other threads are processed.
	process_msgs();
	curr_activated_vertices.init(next_activated_vertices, worker_id,
			graph->get_partitioner());
	next_activated_vertices.clear();
	return curr_activated_vertices.get_num_vertices();
}

void worker_thread::process_multicast_msg(multicast_message &mmsg)
{
	int num_dests = mmsg.get_num_dests();
	multicast_dest_list dest_list = mmsg.get_dest_list();
	for (int i = 0; i < num_dests; i++) {
		vertex_id_t id = dest_list.get_dest(i);
		int part_id;
		off_t off;
		graph->get_partitioner()->map2loc(id, part_id, off);
		assert(part_id == worker_id);
		// TODO now the size is the entire message. Now the message
		// is considered as non-empty.
		if (!mmsg.is_empty()) {
			compute_vertex &info = graph->get_vertex(id);
			const vertex_message *msgs[] = {&mmsg};
			info.run_on_messages(*graph, msgs, 1);
		}
		if (mmsg.is_activate())
			next_activated_vertices.set(off);
	}
}

void worker_thread::process_msg(message &msg)
{
	const int VMSG_BUF_SIZE = 128;
	vertex_message *v_msgs[VMSG_BUF_SIZE];
	while (!msg.is_empty()) {
		int num = msg.get_next(v_msgs, VMSG_BUF_SIZE);
		for (int i = 0; i < num; i++) {
			if (v_msgs[i]->is_multicast()) {
				process_multicast_msg(*multicast_message::cast2multicast(
							v_msgs[i]));
				continue;
			}
			vertex_id_t id = v_msgs[i]->get_dest();
			int part_id;
			off_t off;
			graph->get_partitioner()->map2loc(id, part_id, off);
			assert(part_id == worker_id);
			if (!v_msgs[i]->is_empty()) {
				compute_vertex &info = graph->get_vertex(id);
				info.run_on_messages(*graph,
						(const vertex_message **) &v_msgs[i], 1);
			}
			if (v_msgs[i]->is_activate())
				next_activated_vertices.set(off);
		}
	}
}

void worker_thread::process_msgs()
{
	const int MSG_BUF_SIZE = 16;
	message msgs[MSG_BUF_SIZE];
	while (!msg_q.is_empty()) {
		int num_fetched = msg_q.fetch(msgs, MSG_BUF_SIZE);
		for (int i = 0; i < num_fetched; i++)
			process_msg(msgs[i]);
	}
}

/**
 * This method is the main function of the graph engine.
 */
void worker_thread::run()
{
	while (true) {
		int num_visited = 0;
		int num;
		do {
			num = process_activated_vertices(
					graph->get_max_processing_vertices()
					- io->num_pending_ios());
			num_visited += num;
			process_msgs();
			io->wait4complete(min(io->num_pending_ios() / 10, 2));
			// If there are vertices being processed, we need to call
			// wait4complete to complete processing them.
		} while (get_num_vertices_processing() > 0
				// We still have vertices remaining for processing
				|| !curr_activated_vertices.is_empty()
				// Even if we have processed all activated vertices belonging
				// to this thread, we still need to process vertices from
				// other threads in order to balance the load.
				|| get_num_activated_on_others() > 0);
		assert(curr_activated_vertices.is_empty());
		printf("thread %d visited %d vertices\n", this->get_id(), num_visited);

		// Now we have finished this level, we can progress to the next level.
		num_activated_vertices_in_level = atomic_number<long>(0);
		num_completed_vertices_in_level = atomic_number<long>(0);

		flush_msgs();
		bool completed = graph->progress_next_level();
		printf("thread %d finish in a level, completed? %d\n", get_id(), completed);
		if (completed)
			break;
	}
	stop();
	if (graph_conf.get_print_io_stat())
		io->print_stat(graph->get_num_threads());
}
