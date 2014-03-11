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

#include <algorithm>
#ifndef MEMCHECK
#include <parallel/algorithm>
#endif

#include "io_interface.h"

#include "bitmap.h"
#include "graph_config.h"
#include "graph_engine.h"
#include "messaging.h"
#include "worker_thread.h"
#include "vertex_compute.h"
#include "vertex_request.h"

graph_config graph_conf;
default_vertex_scheduler default_scheduler;

/**
 * This I/O scheduler is to favor maximizing throughput.
 * Therefore, it processes all user tasks together to potentially increase
 * the page cache hit rate.
 */
class throughput_comp_io_scheduler: public comp_io_scheduler
{
	fifo_queue<io_request> req_buf;
public:
	throughput_comp_io_scheduler(int node_id): comp_io_scheduler(
			node_id), req_buf(node_id, 512, true) {
	}

	size_t get_requests(fifo_queue<io_request> &reqs);
};

class throughput_comp_io_sched_creater: public comp_io_sched_creater
{
public:
	comp_io_scheduler *create(int node_id) const {
		return new throughput_comp_io_scheduler(node_id);
	}
};

struct prio_compute
{
	user_compute *compute;
	io_request req;

	prio_compute(io_interface *io, user_compute *compute) {
		this->compute = compute;
		bool ret = compute->fetch_request(io, req);
		assert(ret);
	}
};

class comp_prio_compute
{
public:
	bool operator()(const prio_compute &c1, const prio_compute &c2) {
		// We want the priority queue returns requests with
		// the smallest offset first. If we use less, the priority queue
		// will return the request with the greatest offset.
		return c1.req.get_offset() > c2.req.get_offset();
	}
};

size_t throughput_comp_io_scheduler::get_requests(fifo_queue<io_request> &reqs)
{
	size_t num = 0;
	// Let's add the buffered requests first.
	if (!req_buf.is_empty()) {
		num = reqs.add(&req_buf);
	}

	if (!reqs.is_full()) {
		// Construct a priority queue on user tasks, ordered by the offset
		// of their next requests.
		std::priority_queue<prio_compute, std::vector<prio_compute>,
			comp_prio_compute> user_computes;
		compute_iterator it = this->get_begin();
		compute_iterator end = this->get_end();
		for (; it != end; ++it) {
			user_compute *compute = *it;
			// Skip the ones without user tasks.
			if (!compute->has_requests())
				continue;

			prio_compute prio_comp(get_io(), compute);
			user_computes.push(prio_comp);
		}

		// Add requests to the queue in a sorted order.
		off_t prev = 0;
		int num = 0;
		while (!reqs.is_full() && !user_computes.empty()) {
			num++;
			prio_compute prio_comp = user_computes.top();
			assert(prev <= prio_comp.req.get_offset());
			prev = prio_comp.req.get_offset();
			user_computes.pop();
			assert(prev <= user_computes.top().req.get_offset());
			reqs.push_back(prio_comp.req);
			num++;
			user_compute *compute = prio_comp.compute;
			if (compute->has_requests()) {
				prio_compute prio_comp(get_io(), compute);
				assert(prio_comp.req.get_offset() >= prev);
				user_computes.push(prio_comp);
			}
		}

		// We have got a request from each user task. We can't add them to
		// the queue this time. We need to buffer them.
		while (!user_computes.empty()) {
			prio_compute prio_comp = user_computes.top();
			user_computes.pop();
			if (req_buf.is_full())
				req_buf.expand_queue(req_buf.get_size() * 2);
			req_buf.push_back(prio_comp.req);
		}
	}
	return num;
}

graph_engine::graph_engine(int num_threads, int num_nodes,
		const std::string &graph_file, graph_index *index)
{
	max_processing_vertices = 0;
	this->scheduler = &default_scheduler;
	this->partitioner = new vertex_partitioner(num_threads);
	is_complete = false;
	this->vertices = index;

	pthread_mutex_init(&lock, NULL);
	pthread_barrier_init(&barrier1, NULL, num_threads);
	pthread_barrier_init(&barrier2, NULL, num_threads);

	// Right now only the cached I/O can support async I/O
	factory = create_io_factory(graph_file,
			GLOBAL_CACHE_ACCESS);
	factory->set_sched_creater(new throughput_comp_io_sched_creater());

	io_interface *io = factory->create_io(thread::get_curr_thread());
	io->access((char *) &header, 0, sizeof(header), READ);
	header.verify();
	factory->destroy_io(io);

	switch (header.get_graph_type()) {
		case graph_type::DIRECTED:
			interpreter = std::unique_ptr<ext_mem_vertex_interpreter>(
					new ext_mem_directed_vertex_interpreter());
			break;
		case graph_type::UNDIRECTED:
			interpreter = std::unique_ptr<ext_mem_vertex_interpreter>(
					new ext_mem_undirected_vertex_interpreter());
			break;
		case graph_type::TS_DIRECTED:
			interpreter = std::unique_ptr<ext_mem_vertex_interpreter>(
					new ts_ext_mem_vertex_interpreter(
					header.get_max_num_timestamps()));
			break;
		case graph_type::TS_UNDIRECTED:
			assert(0);
			break;
		default:
			assert(0);
	}

	this->num_nodes = num_nodes;
	assert(num_threads > 0 && num_nodes > 0);
	assert(num_threads % num_nodes == 0);
	worker_threads.resize(num_threads);

	if (graph_conf.get_trace_file().empty())
		logger = NULL;
	else
		logger = new trace_logger(graph_conf.get_trace_file());
}

graph_engine::~graph_engine()
{
	delete partitioner;
	for (unsigned i = 0; i < worker_threads.size(); i++)
		delete worker_threads[i];
	if (logger)
		delete logger;
}

multicast_msg_sender *graph_engine::get_activate_sender(int thread_id) const
{
	worker_thread *curr = (worker_thread *) thread::get_curr_thread();
	return curr->get_activate_sender(thread_id);
}

multicast_msg_sender *graph_engine::get_multicast_sender(int thread_id) const
{
	worker_thread *curr = (worker_thread *) thread::get_curr_thread();
	return curr->get_multicast_sender(thread_id);
}

simple_msg_sender *graph_engine::get_msg_sender(int thread_id) const
{
	worker_thread *curr = (worker_thread *) thread::get_curr_thread();
	return curr->get_msg_sender(thread_id);
}

void graph_engine::start(vertex_id_t ids[], int num)
{
	// Prepare the worker threads.
	int num_threads = get_num_threads();
	for (int i = 0; i < num_threads; i++) {
		worker_thread *t = new worker_thread(this, factory,
				i % num_nodes, i, num_threads);
		worker_threads[i] = t;
		worker_threads[i]->set_vertex_scheduler(scheduler);
	}
	for (int i = 0; i < num_threads; i++) {
		worker_threads[i]->init_messaging(worker_threads);
	}

	num_remaining_vertices_in_level.inc(num);
	std::vector<std::vector<vertex_id_t> > start_vertices(num_threads);
	for (int i = 0; i < num; i++) {
		int idx = get_partitioner()->map(ids[i]);
		start_vertices[idx].push_back(ids[i]);
	}

	for (int i = 0; i < num_threads; i++) {
		worker_threads[i]->start_vertices(start_vertices[i]);
		worker_threads[i]->start();
	}
}

void graph_engine::start_all()
{
	// Prepare the worker threads.
	int num_threads = get_num_threads();
	for (int i = 0; i < num_threads; i++) {
		worker_thread *t = new worker_thread(this, factory,
				i % num_nodes, i, num_threads);
		worker_threads[i] = t;
		worker_threads[i]->set_vertex_scheduler(scheduler);
	}
	for (int i = 0; i < num_threads; i++) {
		worker_threads[i]->init_messaging(worker_threads);
	}

	num_remaining_vertices_in_level.inc(get_num_vertices());
	BOOST_FOREACH(worker_thread *t, worker_threads) {
		t->start_all_vertices();
		t->start();
	}
}

bool graph_engine::progress_next_level()
{
	static atomic_number<long> tot_num_activates;
	static atomic_integer num_threads;
	// We have to make sure all threads have reach here, so we can switch
	// queues to progress to the next level.
	// If the queue of the next level is empty, the program can terminate.
	int rc = pthread_barrier_wait(&barrier1);
	if(rc != 0 && rc != PTHREAD_BARRIER_SERIAL_THREAD)
	{
		printf("Could not wait on barrier\n");
		exit(-1);
	}
	worker_thread *curr = (worker_thread *) thread::get_curr_thread();
	int num_activates = curr->enter_next_level();
	tot_num_activates.inc(num_activates);
	// If all threads have reached here.
	if (num_threads.inc(1) == get_num_threads()) {
		level.inc(1);
		printf("progress to level %d, there are %ld vertices in this level\n",
				level.get(), tot_num_activates.get());
		assert(num_remaining_vertices_in_level.get() == 0);
		num_remaining_vertices_in_level = atomic_number<size_t>(
				tot_num_activates.get());
		// If there aren't more activated vertices.
		is_complete = tot_num_activates.get() == 0;
		tot_num_activates = 0;
		num_threads = 0;
	}

	// We need to synchronize again. We have to make sure all threads see
	// the completion signal.
	rc = pthread_barrier_wait(&barrier2);
	if(rc != 0 && rc != PTHREAD_BARRIER_SERIAL_THREAD)
	{
		printf("Could not wait on barrier\n");
		exit(-1);
	}
	return is_complete;
}

void graph_engine::wait4complete()
{
	for (unsigned i = 0; i < worker_threads.size(); i++) {
		worker_threads[i]->join();
		delete worker_threads[i];
		worker_threads[i] = NULL;
	}
}

void graph_engine::set_vertex_scheduler(vertex_scheduler *scheduler)
{
	this->scheduler = scheduler;
}

void graph_engine::request_vertices(compute_vertex &vertex, vertex_id_t ids[],
		int num)
{
	worker_thread *curr = (worker_thread *) thread::get_curr_thread();
	vertex_compute *compute = curr->get_curr_vertex_compute();
	if (compute == NULL)
		compute = curr->create_vertex_compute();
	compute->request_vertices(ids, num);
}

void graph_engine::request_partial_vertices(compute_vertex &vertex,
		vertex_request *reqs[], int num)
{
	worker_thread *curr = (worker_thread *) thread::get_curr_thread();
	vertex_compute *compute = curr->get_curr_vertex_compute();
	if (compute == NULL)
		compute = curr->create_vertex_compute();
	compute->request_partial_vertices(reqs, num);
}

vertex_index *load_vertex_index(const std::string &index_file)
{
	const int INDEX_HEADER_SIZE = PAGE_SIZE * 2;
	const int READ_SIZE = 100 * 1024 * 1024;

	// Right now only the cached I/O can support async I/O
	file_io_factory *factory = create_io_factory(index_file,
			REMOTE_ACCESS);
	assert(factory->get_file_size() >= INDEX_HEADER_SIZE);
	io_interface *io = factory->create_io(thread::get_curr_thread());

	// Get the header of the index.
	char *tmp = NULL;
	int ret = posix_memalign((void **) &tmp, PAGE_SIZE, INDEX_HEADER_SIZE);
	assert(ret == 0);
	data_loc_t loc(factory->get_file_id(), 0);
	io_request req(tmp, loc, INDEX_HEADER_SIZE, READ, io, -1);
	io->access(&req, 1);
	io->wait4complete(1);
	vertex_index *index = (vertex_index *) tmp;
	index->get_graph_header().verify();

	// Initialize the buffer for containing the index.
	size_t index_size = index->get_index_size();
	assert((ssize_t) index_size <= factory->get_file_size());
	char *buf = NULL;
	ret = posix_memalign((void **) &buf, PAGE_SIZE, index_size);
	assert(ret == 0);
	off_t off = 0;
	memcpy(buf, tmp, INDEX_HEADER_SIZE);
	off += INDEX_HEADER_SIZE;
	free(tmp);

	// Read the index to the memory.
	size_t aligned_index_size = ROUND_PAGE(index_size);
	while ((size_t) off < aligned_index_size) {
		assert(off % PAGE_SIZE == 0);
		size_t size = min(READ_SIZE, aligned_index_size - off);
		data_loc_t loc(factory->get_file_id(), off);
		io_request req(buf + off, loc, size, READ, io, -1);
		io->access(&req, 1);
		off += size;
		if (io->num_pending_ios() > 100)
			io->wait4complete(io->num_pending_ios() / 10);
	}
	io->wait4complete(io->num_pending_ios());

	// Read the last page.
	// The data may only occupy part of the page.
	if (aligned_index_size < index_size) {
		char *tmp = NULL;
		int ret = posix_memalign((void **) &tmp, PAGE_SIZE, PAGE_SIZE);
		assert(ret == 0);
		data_loc_t loc(factory->get_file_id(), aligned_index_size);
		io_request req(tmp, loc, PAGE_SIZE, READ, io, -1);
		io->access(&req, 1);
		io->wait4complete(1);
		memcpy(buf + aligned_index_size, tmp, index_size - aligned_index_size);
		free(tmp);
	}
	factory->destroy_io(io);

	index = (vertex_index *) buf;
	if (index->get_graph_header().get_graph_type() == graph_type::DIRECTED)
		((directed_vertex_index *) index)->verify();
	else
		((default_vertex_index *) index)->verify();
	return index;
}
