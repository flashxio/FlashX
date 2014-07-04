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
#include "vertex_index_reader.h"

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

	size_t get_requests(fifo_queue<io_request> &reqs, int max);

	bool is_empty() const {
		return user_computes.empty();
	}
};

template<class prio_queue_type>
size_t comp_io_schedule_queue<prio_queue_type>::get_requests(
		fifo_queue<io_request> &reqs, int max)
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
				compute->set_scan_dir(forward);
				prio_compute prio_comp(scheduler->get_io(), compute);
				user_computes.push(prio_comp);
			}
		}

		// Add requests to the queue in a sorted order.
		int num = 0;
		while (!reqs.is_full() && !user_computes.empty() && num < max) {
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

	size_t get_requests(fifo_queue<io_request> &reqs, size_t max);
};

class throughput_comp_io_sched_creater: public comp_io_sched_creater
{
public:
	comp_io_scheduler *create(int node_id) const {
		return new throughput_comp_io_scheduler(node_id);
	}
};

size_t throughput_comp_io_scheduler::get_requests(fifo_queue<io_request> &reqs,
		size_t max)
{
	size_t ret;
	if (graph_conf.get_elevator_enabled()) {
		if (batch_num % 2 == 0) {
			ret = forward_queue.get_requests(reqs, max);
			if (forward_queue.is_empty())
				batch_num++;
		}
		else {
			ret = backward_queue.get_requests(reqs, max);
			if (backward_queue.is_empty())
				batch_num++;
		}
	}
	else
		ret = forward_queue.get_requests(reqs, max);
	assert(ret <= max);
	return ret;
}

void compute_vertex::request_vertices(vertex_id_t ids[], size_t num)
{
	worker_thread *curr = (worker_thread *) thread::get_curr_thread();
	vertex_compute *compute = curr->get_vertex_compute(*this);
	compute->request_vertices(ids, num);
}

void compute_vertex::request_num_edges(vertex_id_t ids[], size_t num)
{
	worker_thread *curr = (worker_thread *) thread::get_curr_thread();
	vertex_compute *compute = curr->get_vertex_compute(*this);
	compute->request_num_edges(ids, num);
}

void compute_directed_vertex::request_partial_vertices(
		directed_vertex_request reqs[], size_t num)
{
	worker_thread *curr = (worker_thread *) thread::get_curr_thread();
	directed_vertex_compute *compute
		= (directed_vertex_compute *) curr->get_vertex_compute(*this);
	compute->request_partial_vertices(reqs, num);
}

#if 0
void compute_ts_vertex::request_partial_vertices(ts_vertex_request reqs[],
		size_t num)
{
	worker_thread *curr = (worker_thread *) thread::get_curr_thread();
	ts_vertex_compute *compute
		= (ts_vertex_compute *) curr->get_vertex_compute(*this);
	compute->request_partial_vertices(reqs, num);
}
#endif

graph_engine::graph_engine(const std::string &graph_file,
		graph_index::ptr index, const config_map &configs)
{
	graph_conf.init(configs);
	graph_conf.print();
	init_io_system(configs);
	int num_threads = graph_conf.get_num_threads();
	this->num_nodes = params.get_num_nodes();
	index->init(num_threads, num_nodes);

	max_processing_vertices = 0;
	is_complete = false;
	this->vertices = index;

	pthread_mutex_init(&lock, NULL);
	pthread_barrier_init(&barrier1, NULL, num_threads);
	pthread_barrier_init(&barrier2, NULL, num_threads);

	// Right now only the cached I/O can support async I/O
	graph_factory = create_io_factory(graph_file, GLOBAL_CACHE_ACCESS);
	graph_factory->set_sched_creater(new throughput_comp_io_sched_creater());
	index_factory = create_io_factory(index->get_index_file(), GLOBAL_CACHE_ACCESS);

	io_interface::ptr io = graph_factory->create_io(thread::get_curr_thread());
	io->access((char *) &header, 0, sizeof(header), READ);
	header.verify();

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

	assert(num_threads > 0 && num_nodes > 0);
	assert(num_threads % num_nodes == 0);
	worker_threads.resize(num_threads);
	vprograms.resize(num_threads);

	if (!graph_conf.get_trace_file().empty())
		logger = trace_logger::ptr(new trace_logger(graph_conf.get_trace_file()));

	if (graph_conf.preload())
		preload_graph();
}

graph_engine::~graph_engine()
{
	for (unsigned i = 0; i < worker_threads.size(); i++)
		delete worker_threads[i];
	graph_factory = file_io_factory::shared_ptr();
	index_factory = file_io_factory::shared_ptr();
	destroy_io_system();
}

void graph_engine::init_threads(vertex_program_creater::ptr creater)
{
	// Prepare the worker threads.
	int num_threads = get_num_threads();
	for (int i = 0; i < num_threads; i++) {
		vertex_program::ptr new_prog;
		if (creater)
			new_prog = creater->create();
		else
			new_prog = vertices->create_def_vertex_program();
		worker_thread *t = new worker_thread(this, graph_factory, index_factory,
				new_prog, i % num_nodes, i, num_threads, scheduler);
		assert(worker_threads[i] == NULL);
		worker_threads[i] = t;
		vprograms[i] = new_prog;
	}
	for (int i = 0; i < num_threads; i++) {
		worker_threads[i]->init_messaging(worker_threads);
	}
}

void graph_engine::start(vertex_id_t ids[], int num,
		vertex_initializer::ptr init, vertex_program_creater::ptr creater)
{
	init_threads(std::move(creater));
	num_remaining_vertices_in_level.inc(num);
	int num_threads = get_num_threads();
	std::vector<std::vector<vertex_id_t> > start_vertices(num_threads);
	for (int i = 0; i < num; i++) {
		int idx = get_partitioner()->map(ids[i]);
		start_vertices[idx].push_back(ids[i]);
	}

	for (int i = 0; i < num_threads; i++) {
		worker_threads[i]->start_vertices(start_vertices[i], init);
		worker_threads[i]->start();
	}
	gettimeofday(&start_time, NULL);
}

void graph_engine::start(std::shared_ptr<vertex_filter> filter,
		vertex_program_creater::ptr creater)
{
	init_threads(std::move(creater));
	// Let's assume all vertices will be activated first.
	num_remaining_vertices_in_level.inc(get_num_vertices());
	BOOST_FOREACH(worker_thread *t, worker_threads) {
		t->start_vertices(filter);
		t->start();
	}
	gettimeofday(&start_time, NULL);
}

void graph_engine::start_all(vertex_initializer::ptr init,
		vertex_program_creater::ptr creater)
{
	init_threads(std::move(creater));
	num_remaining_vertices_in_level.inc(get_num_vertices());
	BOOST_FOREACH(worker_thread *t, worker_threads) {
		t->start_all_vertices(init);
		t->start();
	}
	gettimeofday(&start_time, NULL);
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
		struct timeval curr;
		gettimeofday(&curr, NULL);
		printf("Iter %d takes %.2f seconds, and %ld vertices are in iter %d\n",
				level.get() - 1, time_diff(start_time, curr),
				tot_num_activates.get(), level.get());
		start_time = curr;
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

void graph_engine::set_vertex_scheduler(vertex_scheduler::ptr scheduler)
{
	this->scheduler = scheduler;
}

void graph_engine::preload_graph()
{
	const int BLOCK_SIZE = 1024 * 1024 * 32;
	std::unique_ptr<char[]> buf = std::unique_ptr<char[]>(new char[BLOCK_SIZE]);

	size_t cache_size = params.get_cache_size();
	io_interface::ptr io = index_factory->create_io(thread::get_curr_thread());
	size_t preload_size = min(cache_size,
			index_factory->get_file_size());
	cache_size -= preload_size;
	printf("preload %ld bytes of the index file\n", preload_size);
	for (size_t i = 0; i < preload_size; i += BLOCK_SIZE)
		io->access(buf.get(), i, BLOCK_SIZE, READ);

	io = graph_factory->create_io(thread::get_curr_thread());
	preload_size = min(cache_size, graph_factory->get_file_size());
	printf("preload %ld bytes of the graph file\n", preload_size);
	for (size_t i = 0; i < preload_size; i += BLOCK_SIZE)
		io->access(buf.get(), i, BLOCK_SIZE, READ);
	printf("successfully preload\n");
}

void graph_engine::init_vertices(vertex_id_t ids[], int num,
		vertex_initializer::ptr init)
{
#pragma omp parallel for
	for (int i = 0; i < num; i++) {
		init->init(get_vertex(ids[i]));
	}
}

void graph_engine::init_all_vertices(vertex_initializer::ptr init)
{
	vertex_id_t max_id = get_max_vertex_id();
#pragma omp parallel for
	for (vertex_id_t id = 0; id <= max_id; id++) {
		init->init(get_vertex(id));
	}
}

class query_thread: public thread
{
	graph_engine &graph;
	vertex_query::ptr query;
	int part_id;
public:
	query_thread(graph_engine &_graph, vertex_query::ptr query, int part_id,
			int node_id): thread("query_thread", node_id), graph(_graph) {
		this->query = query;
		this->part_id = part_id;
	}

	void run();

	vertex_query::ptr get_query() {
		return query;
	}
};

void query_thread::run()
{
	size_t part_size = graph.get_partitioner()->get_part_size(part_id,
			graph.get_num_vertices());
	// We only iterate over the vertices in the local partition.
	for (vertex_id_t id = 0; id < part_size; id++) {
		local_vid_t local_id(id);
		compute_vertex &v = graph.get_vertex(part_id, local_id);
		query->run(graph, v);
	}
	stop();
}

void graph_engine::query_on_all(vertex_query::ptr query)
{
	std::vector<query_thread *> threads(get_num_threads());
	for (size_t i = 0; i < threads.size(); i++) {
		threads[i] = new query_thread(*this, query->clone(),
				i, i % num_nodes);
		threads[i]->start();
	}
	for (size_t i = 0; i < threads.size(); i++) {
		threads[i]->join();
		query->merge(*this, threads[i]->get_query());
		delete threads[i];
	}
}

size_t graph_get_vertices(graph_engine &graph, const worker_thread &t,
		const local_vid_t ids[], int num_ids, compute_vertex *v_buf[])
{
	return graph.get_vertices(t.get_worker_id(), ids, num_ids, v_buf);
}
