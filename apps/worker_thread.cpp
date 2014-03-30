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

#include <atomic>

#include "worker_thread.h"
#include "graph_engine.h"
#include "bitmap.h"
#include "vertex_compute.h"
#include "message_processor.h"
#include "load_balancer.h"
#include "steal_state.h"

sorted_vertex_queue::sorted_vertex_queue(graph_engine &_graph): fetch_idx(
		0, true), graph(_graph)
{
	pthread_spin_init(&lock, PTHREAD_PROCESS_PRIVATE);
	this->scheduler = &default_scheduler;
}

void sorted_vertex_queue::init(const bitmap &map, int part_id,
		const graph_partitioner *partitioner)
{
	pthread_spin_lock(&lock);
	sorted_vertices.clear();
	std::vector<vertex_id_t> local_ids;
	map.get_set_bits(local_ids);
	// the bitmap only contains the locations of vertices in the bitmap.
	// We have to translate them back to vertex ids.
	sorted_vertices.resize(local_ids.size());
	for (size_t i = 0; i < local_ids.size(); i++) {
		vertex_id_t id;
		partitioner->loc2map(part_id, local_ids[i], id);
		sorted_vertices[i] = &graph.get_vertex(id);
	}

	if (scheduler != &default_scheduler)
		scheduler->schedule(sorted_vertices);
	bool forward = true;
	if (graph_conf.get_elevator_enabled())
		forward = graph.get_curr_level() % 2;
	fetch_idx = scan_pointer(sorted_vertices.size(), forward);
	pthread_spin_unlock(&lock);
}

worker_thread::worker_thread(graph_engine *graph, file_io_factory::shared_ptr factory,
		vertex_program::ptr prog, int node_id, int worker_id,
		int num_threads): thread("worker_thread",
			node_id), next_activated_vertices(graph->get_partitioner(
					)->get_part_size(worker_id, graph->get_num_vertices()), node_id),
			curr_activated_vertices(*graph)
{
	this->vprogram = std::move(prog);
	start_all = false;
	this->worker_id = worker_id;
	this->graph = graph;
	this->io = NULL;
	this->factory = factory;
	this->curr_compute = NULL;
	msg_processor = std::unique_ptr<message_processor>(new message_processor(
				*graph, *this));
	balancer = std::unique_ptr<load_balancer>(new load_balancer(*graph, *this));
	switch(graph->get_graph_header().get_graph_type()) {
		case graph_type::DIRECTED:
			alloc = new vertex_compute_allocator<directed_vertex_compute>(graph, this);
			part_alloc = new vertex_compute_allocator<part_directed_vertex_compute>(
					graph, this);
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
	// We increase the allocator by 1M each time.
	// It shouldn't need to allocate much memory.
	msg_alloc = new slab_allocator("graph-message-allocator",
			GRAPH_MSG_BUF_SIZE, 1024 * 1024, INT_MAX, get_node_id());

	int num_self = 0;
	for (unsigned i = 0; i < threads.size(); i++) {
		if (threads[i] == this)
			num_self++;
		msg_senders.push_back(simple_msg_sender::create(get_node_id(),
					msg_alloc, &threads[i]->msg_processor->get_msg_queue()));
		multicast_senders.push_back(multicast_msg_sender::create(msg_alloc,
					&threads[i]->msg_processor->get_msg_queue()));
		multicast_msg_sender *activate_sender = multicast_msg_sender::create(
				msg_alloc, &threads[i]->msg_processor->get_msg_queue());
		// All activation messages are the same. We can initialize the sender
		// here.
		vertex_message msg(sizeof(vertex_message), true);
		activate_sender->init(msg);
		activate_senders.push_back(activate_sender);
	}
	assert(num_self == 1);
}

/**
 * This is to process the activated vertices in the current iteration.
 */
int worker_thread::process_activated_vertices(int max)
{
	if (max <= 0)
		return 0;

	compute_vertex *vertex_buf[max];
	stack_array<io_request> reqs(max);
	int num = curr_activated_vertices.fetch(vertex_buf, max);
	if (num == 0) {
		num = balancer->steal_activated_vertices(vertex_buf, max);
	}
	if (num > 0) {
		num_activated_vertices_in_level.inc(num);
		graph->process_vertices(num);
	}

	int num_to_process = 0;
	for (int i = 0; i < num; i++) {
		compute_vertex *info = vertex_buf[i];
		// We execute the pre-run to determine if the vertex has completed
		// in the current iteration.
		vertex_program &curr_vprog = get_vertex_program();
		assert(curr_compute == NULL);
		curr_vprog.run(*graph, *info);
		if (curr_compute) {
			// If the user code requests the vertices that are empty or whose
			// requested part is empty. These empty requests can be handled
			// immediately, so it's possible that the current vertex compute
			// may not have requests.
			if (curr_compute->has_requests()) {
				// It's mostly likely that it is requesting the adjacency list
				// of itself. But it doesn't really matter what the vertex
				// wants to request here.
				request_range range = curr_compute->get_next_request();
				reqs[num_to_process++] = io_request(range.get_compute(),
						range.get_loc(), range.get_size(),
						// TODO I might need to set the node id.
						range.get_access_method(), io, -1);
			}
			else {
				// The reason we reach here is that the vertex requests some
				// partial vertices and the request parts are empty.
				// The user compute is only referenced here. We need to delete
				// it.
				assert(curr_compute->get_ref() == 0);
				assert(curr_compute->has_completed());
				compute_allocator *alloc = curr_compute->get_allocator();
				alloc->free(curr_compute);
			}
		}
		else
			complete_vertex(*info);
		reset_curr_vertex_compute();
	}
	if (graph->get_logger())
		graph->get_logger()->log(reqs.data(), num_to_process);
	io->access(reqs.data(), num_to_process);
	return num;
}

int worker_thread::enter_next_level()
{
	// We have to make sure all messages sent by other threads are processed.
	msg_processor->process_msgs();
	curr_activated_vertices.init(next_activated_vertices, worker_id,
			graph->get_partitioner());
	next_activated_vertices.clear();
	balancer->reset();
	msg_processor->reset();
	return curr_activated_vertices.get_num_vertices();
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
			balancer->process_completed_stolen_vertices();
			num = process_activated_vertices(
					graph->get_max_processing_vertices()
					- io->num_pending_ios());
			num_visited += num;
			msg_processor->process_msgs();
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
		printf("worker %d visited %d vertices\n", worker_id, num_visited);
		assert(num_visited == num_activated_vertices_in_level.get());
		if (num_visited != num_completed_vertices_in_level.get()) {
			printf("worker %d: visits %d vertices and completes %ld\n",
					worker_id, num_visited, num_completed_vertices_in_level.get());
		}
		assert(num_visited == num_completed_vertices_in_level.get());

		// Now we have finished this level, we can progress to the next level.
		num_activated_vertices_in_level = atomic_number<long>(0);
		num_completed_vertices_in_level = atomic_number<long>(0);

		flush_msgs();
		// We have to make sure all stolen vertices are returned to their owner
		// threads.
		balancer->process_completed_stolen_vertices();
		balancer->reset();
		bool completed = graph->progress_next_level();
		printf("thread %d finish in a level, completed? %d\n", get_id(), completed);
		if (completed)
			break;
	}
	stop();
	if (graph_conf.get_print_io_stat())
		io->print_stat(graph->get_num_threads());
}

int worker_thread::steal_activated_vertices(compute_vertex *vertices[], int num)
{
	// We want to steal as much as possible, but we don't want
	// to overloaded by the stolen vertices.
	size_t num_steal = max(1,
			curr_activated_vertices.get_num_vertices() / graph->get_num_threads());
	num = curr_activated_vertices.fetch(vertices,
			min(num, num_steal));
	if (num > 0)
		// If the thread steals vertices from another thread successfully,
		// it needs to notify the thread of the stolen vertices.
		msg_processor->steal_vertices(vertices, num);
	return num;
}

void worker_thread::return_vertices(vertex_id_t ids[], int num)
{
	msg_processor->return_vertices(ids, num);
}

void worker_thread::complete_vertex(const compute_vertex &v)
{
	num_completed_vertices_in_level.inc(1);
	// The vertex might be stolen from another thread. Now we have
	// finished processing it, we should return it to its owner thread.
	int part_id = graph->get_partitioner()->map(v.get_id());
	if (part_id != worker_id) {
		vertex_id_t id = v.get_id();
		balancer->return_vertices(&id, 1);
	}
}

void worker_thread::activate_vertex(vertex_id_t id)
{
	int part_id;
	off_t off;
	graph->get_partitioner()->map2loc(id, part_id, off);
	assert(part_id == worker_id);
	next_activated_vertices.set(off);
}

vertex_compute *worker_thread::create_vertex_compute(compute_vertex *v)
{
	assert(curr_compute == NULL);
	curr_compute = (vertex_compute *) alloc->alloc();
	curr_compute->init(v);
	return curr_compute;
}
