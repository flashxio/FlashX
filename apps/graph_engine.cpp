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

class forward_comp_prio_compute
{
public:
	bool operator()(const prio_compute &c1, const prio_compute &c2) {
		// We want the priority queue returns requests with
		// the smallest offset first. If we use less, the priority queue
		// will return the request with the greatest offset.
		return c1.req.get_offset() > c2.req.get_offset();
	}
};

class backward_comp_prio_compute
{
public:
	bool operator()(const prio_compute &c1, const prio_compute &c2) {
		// We want the priority queue returns requests with
		// the smallest offset first. If we use less, the priority queue
		// will return the request with the greatest offset.
		return c1.req.get_offset() < c2.req.get_offset();
	}
};

template<class prio_queue_type>
class comp_io_schedule_queue
{
	bool forward;
	// Construct a priority queue on user tasks, ordered by the offset
	// of their next requests.
	prio_queue_type user_computes;
	comp_io_scheduler *scheduler;
public:
	comp_io_schedule_queue(comp_io_scheduler *scheduler, bool forward) {
		this->scheduler = scheduler;
		this->forward = forward;
	}

	size_t get_requests(fifo_queue<io_request> &reqs);

	bool is_empty() const {
		return user_computes.empty();
	}
};

template<class prio_queue_type>
size_t comp_io_schedule_queue<prio_queue_type>::get_requests(
		fifo_queue<io_request> &reqs)
{
	size_t num = 0;

	if (!reqs.is_full()) {
		// we don't have user computes, we can get some from the queue of
		// incomplete user computes in comp_io_scheduler.
		if (user_computes.empty()) {
			comp_io_scheduler::compute_iterator it = scheduler->get_begin();
			comp_io_scheduler::compute_iterator end = scheduler->get_end();
			for (; it != end; ++it) {
				user_compute *compute = *it;
				// Skip the ones without user tasks.
				if (!compute->has_requests())
					continue;

				// We have a reference to the user compute. Let's increase
				// its ref count. User computes should be in the queue of
				// comp_io_scheduler as long as they can generate more
				// requests. It might not be necessary to increase the ref
				// count, but it can work as a sanity check.
				compute->inc_ref();
				((vertex_compute *) compute)->set_scan_dir(forward);
				prio_compute prio_comp(scheduler->get_io(), compute);
				user_computes.push(prio_comp);
			}
		}

		// Add requests to the queue in a sorted order.
		int num = 0;
		while (!reqs.is_full() && !user_computes.empty()) {
			num++;
			prio_compute prio_comp = user_computes.top();
			user_computes.pop();
			reqs.push_back(prio_comp.req);
			num++;
			user_compute *compute = prio_comp.compute;
			if (compute->has_requests()) {
				prio_compute prio_comp(scheduler->get_io(), compute);
				user_computes.push(prio_comp);
			}
			else {
				// This user compute no longer needs to stay in the priority
				// queue. We can decrease its ref count.
				compute->dec_ref();
			}
		}
	}
	return num;
}

/**
 * This I/O scheduler is to favor maximizing throughput.
 * Therefore, it processes all user tasks together to potentially increase
 * the page cache hit rate.
 */
class throughput_comp_io_scheduler: public comp_io_scheduler
{
	// The scheduler process user computes in batches. It reads all
	// available user computes in the beginning of a batch and then processes
	// them. A batch ends if the scheduler processes all of the user computes.
	size_t batch_num;
	comp_io_schedule_queue<std::priority_queue<prio_compute,
		std::vector<prio_compute>, forward_comp_prio_compute> > forward_queue;
	comp_io_schedule_queue<std::priority_queue<prio_compute,
		std::vector<prio_compute>, backward_comp_prio_compute> > backward_queue;
public:
	throughput_comp_io_scheduler(int node_id): comp_io_scheduler(
			node_id), forward_queue(this, true), backward_queue(this, false) {
		batch_num = 0;
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

size_t throughput_comp_io_scheduler::get_requests(fifo_queue<io_request> &reqs)
{
	size_t ret;
	if (graph_conf.get_elevator_enabled()) {
		if (batch_num % 2 == 0) {
			ret = forward_queue.get_requests(reqs);
			if (forward_queue.is_empty())
				batch_num++;
		}
		else {
			ret = backward_queue.get_requests(reqs);
			if (backward_queue.is_empty())
				batch_num++;
		}
	}
	else
		ret = forward_queue.get_requests(reqs);
	return ret;
}

void compute_vertex::request_vertices(vertex_id_t ids[], size_t num)
{
	worker_thread *curr = (worker_thread *) thread::get_curr_thread();
	vertex_compute *compute = curr->get_curr_vertex_compute();
	if (compute == NULL)
		compute = curr->create_vertex_compute(this);
	compute->request_vertices(ids, num);
}

vsize_t compute_vertex::get_num_edges() const
{
	worker_thread *t = (worker_thread *) thread::get_curr_thread();
	in_mem_vertex_info info = t->get_graph().get_vertex_info(get_id());
	assert(!t->get_graph().get_graph_header().has_edge_data());
	int vertex_header_size = t->get_graph().get_vertex_header_size();
	assert(vertex_header_size > 0);
	return (info.get_ext_mem_size() - vertex_header_size) / sizeof(vertex_id_t);
}

void compute_directed_vertex::request_partial_vertices(
		directed_vertex_request reqs[], size_t num)
{
	worker_thread *curr = (worker_thread *) thread::get_curr_thread();
	directed_vertex_compute *compute
		= (directed_vertex_compute *) curr->get_curr_vertex_compute();
	if (compute == NULL)
		compute = (directed_vertex_compute *) curr->create_vertex_compute(this);
	compute->request_partial_vertices(reqs, num);
}

void compute_ts_vertex::request_partial_vertices(ts_vertex_request reqs[],
		size_t num)
{
	worker_thread *curr = (worker_thread *) thread::get_curr_thread();
	ts_vertex_compute *compute
		= (ts_vertex_compute *) curr->get_curr_vertex_compute();
	if (compute == NULL)
		compute = (ts_vertex_compute *) curr->create_vertex_compute(this);
	compute->request_partial_vertices(reqs, num);
}

graph_engine::graph_engine(int num_threads, int num_nodes,
		const std::string &graph_file, graph_index *index)
{
	max_processing_vertices = 0;
	this->scheduler = NULL;
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
			vertex_header_size = ext_mem_directed_vertex::get_header_size();
			break;
		case graph_type::UNDIRECTED:
			interpreter = std::unique_ptr<ext_mem_vertex_interpreter>(
					new ext_mem_undirected_vertex_interpreter());
			vertex_header_size = ext_mem_undirected_vertex::get_header_size();
			break;
		case graph_type::TS_DIRECTED:
			interpreter = std::unique_ptr<ext_mem_vertex_interpreter>(
					new ts_ext_mem_vertex_interpreter(
					header.get_max_num_timestamps()));
			vertex_header_size = -1;
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
	for (unsigned i = 0; i < worker_threads.size(); i++)
		delete worker_threads[i];
	if (logger)
		delete logger;
}

void graph_engine::init_threads(vertex_program::ptr prog)
{
	// Prepare the worker threads.
	int num_threads = get_num_threads();
	for (int i = 0; i < num_threads; i++) {
		vertex_program::ptr new_prog;
		if (prog)
			new_prog = prog->clone();
		else
			new_prog = vertices->create_def_vertex_program();
		worker_thread *t = new worker_thread(this, factory, std::move(new_prog),
				i % num_nodes, i, num_threads, scheduler);
		assert(worker_threads[i] == NULL);
		worker_threads[i] = t;
	}
	for (int i = 0; i < num_threads; i++) {
		worker_threads[i]->init_messaging(worker_threads);
	}
}

void graph_engine::start(vertex_id_t ids[], int num, vertex_program::ptr prog)
{
	init_threads(std::move(prog));
	num_remaining_vertices_in_level.inc(num);
	int num_threads = get_num_threads();
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

void graph_engine::start(std::shared_ptr<vertex_filter> filter,
		vertex_program::ptr prog)
{
	init_threads(std::move(prog));
	// Let's assume all vertices will be activated first.
	num_remaining_vertices_in_level.inc(get_num_vertices());
	BOOST_FOREACH(worker_thread *t, worker_threads) {
		t->start_vertices(filter);
		t->start();
	}
}

void graph_engine::start_all(vertex_program::ptr prog)
{
	init_threads(std::move(prog));
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

void graph_engine::preload_graph()
{
	const int BLOCK_SIZE = 1024 * 1024 * 32;
	io_interface *io = factory->create_io(thread::get_curr_thread());
	size_t preload_size = min(params.get_cache_size(), factory->get_file_size());
	printf("preload %ld bytes\n", preload_size);
	char *buf = new char[BLOCK_SIZE];
	for (size_t i = 0; i < preload_size; i += BLOCK_SIZE)
		io->access(buf, i, BLOCK_SIZE, READ);
	factory->destroy_io(io);
	printf("successfully preload\n");
}

void graph_engine::activate_vertices(vertex_id_t ids[], int num)
{
	worker_thread *curr = (worker_thread *) thread::get_curr_thread();
	curr->send_activation(ids, num);
}

void graph_engine::multicast_msg(vertex_id_t ids[], int num,
		const vertex_message &msg)
{
	worker_thread *curr = (worker_thread *) thread::get_curr_thread();
	curr->multicast_msg(ids, num, msg);
}

void graph_engine::send_msg(vertex_id_t dest, vertex_message &msg)
{
	worker_thread *curr = (worker_thread *) thread::get_curr_thread();
	curr->send_msg(dest, msg);
}

vertex_index *load_vertex_index(const std::string &index_file)
{
	const int INDEX_HEADER_SIZE = PAGE_SIZE * 2;
	const int READ_SIZE = 100 * 1024 * 1024;

	// Right now only the cached I/O can support async I/O
	file_io_factory::shared_ptr factory = create_io_factory(index_file,
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

size_t graph_get_vertices(graph_engine &graph, const worker_thread &t,
		const local_vid_t ids[], int num_ids, compute_vertex *v_buf[])
{
	return graph.get_vertices(t.get_worker_id(), ids, num_ids, v_buf);
}
