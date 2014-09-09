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

#include "io_interface.h"

#include "bitmap.h"
#include "graph_config.h"
#include "graph_engine.h"
#include "messaging.h"
#include "worker_thread.h"
#include "vertex_compute.h"
#include "vertex_request.h"
#include "vertex_index_reader.h"
#include "in_mem_storage.h"

/**
 * The size of a message buffer used to pass vertex messages to other threads.
 */
const int GRAPH_MSG_BUF_SIZE = PAGE_SIZE * 4;
const int MAX_FLUSH_MSG_SIZE = 256;

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
	int num = 0;

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
		while (!reqs.is_full() && !user_computes.empty() && num < max) {
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

bool request_self(vertex_id_t ids[], size_t num, vertex_id_t self)
{
	return num == 1 && ids[0] == self;
}

bool request_self(directed_vertex_request reqs[], size_t num, vertex_id_t self)
{
	if (num == 1)
		return reqs[0].get_id() == self;
	// If there are two requests, we should make sure that they request
	// different parts.
	else if (num == 2)
		return reqs[0].get_id() == self && reqs[1].get_id() == self
			&& reqs[0].get_type() != reqs[1].get_type();
	else
		return false;
}

void compute_vertex::request_vertices(vertex_id_t ids[], size_t num)
{
	worker_thread *curr = (worker_thread *) thread::get_curr_thread();
	curr->request_on_vertex(get_id());
	if (request_self(ids, num, get_id())) {
		if (curr->get_graph().is_directed()) {
			directed_vertex_request req(ids[0], BOTH_EDGES);
			curr->get_index_reader().request_vertex(req);
		}
		else
			curr->get_index_reader().request_vertex(ids[0]);
	}
	else {
		compute_vertex_pointer curr_vertex = curr->get_curr_vertex();
		assert(curr_vertex.is_valid());
		assert(curr_vertex->get_id() == this->get_id());
		vertex_compute *compute = curr->get_vertex_compute(curr_vertex);
		compute->request_vertices(ids, num);
	}
}

void compute_vertex::request_vertex_headers(vertex_id_t ids[], size_t num)
{
	worker_thread *curr = (worker_thread *) thread::get_curr_thread();
	curr->request_on_vertex(get_id());
	compute_vertex_pointer curr_vertex = curr->get_curr_vertex();
	assert(curr_vertex.is_valid());
	assert(curr_vertex->get_id() == this->get_id());
	vertex_compute *compute = curr->get_vertex_compute(curr_vertex);
	compute->request_num_edges(ids, num);
}

void compute_directed_vertex::request_partial_vertices(
		directed_vertex_request reqs[], size_t num)
{
	worker_thread *curr = (worker_thread *) thread::get_curr_thread();
	curr->request_on_vertex(get_id());
	if (request_self(reqs, num, get_id())) {
		for (size_t i = 0; i < num; i++)
			curr->get_index_reader().request_vertex(reqs[i]);
	}
	else {
		compute_vertex_pointer curr_vertex = curr->get_curr_vertex();
		assert(curr_vertex.is_valid());
		assert(curr_vertex->get_id() == this->get_id());
		directed_vertex_compute *compute
			= (directed_vertex_compute *) curr->get_vertex_compute(curr_vertex);
		compute->request_partial_vertices(reqs, num);
	}
}

void part_compute_vertex::broadcast_vpart(const vertex_message &msg)
{
	worker_thread *curr = (worker_thread *) thread::get_curr_thread();
	std::vector<compute_vertex_pointer> ps(graph_conf.get_num_vparts());
	int ret = curr->get_graph().get_graph_index()->get_vpart_vertices(get_id(),
			ps.data(), ps.size());
	for (int i = 0; i < ret; i++) {
		compute_vertex_pointer v = ps[i];
		curr->get_vertex_program(true).run_on_message(*v, msg);
	}
}

void part_compute_vertex::request_vertices(vertex_id_t ids[], size_t num)
{
	// The trick for self request doesn't work for part compute vertex.
	worker_thread *curr = (worker_thread *) thread::get_curr_thread();
	curr->request_on_vertex(get_id());
	compute_vertex_pointer curr_vertex = curr->get_curr_vertex();
	assert(curr_vertex.is_valid());
	assert(curr_vertex->get_id() == this->get_id());
	vertex_compute *compute = curr->get_vertex_compute(curr_vertex);
	compute->request_vertices(ids, num);
}

void part_compute_directed_vertex::request_vertices(vertex_id_t ids[], size_t num)
{
	// The trick for self request doesn't work for part compute vertex.
	worker_thread *curr = (worker_thread *) thread::get_curr_thread();
	curr->request_on_vertex(get_id());
	compute_vertex_pointer curr_vertex = curr->get_curr_vertex();
	assert(curr_vertex.is_valid());
	assert(curr_vertex->get_id() == this->get_id());
	vertex_compute *compute = curr->get_vertex_compute(curr_vertex);
	compute->request_vertices(ids, num);
}

void part_compute_directed_vertex::request_partial_vertices(
		directed_vertex_request reqs[], size_t num)
{
	// The trick for self request doesn't work for part compute vertex.
	worker_thread *curr = (worker_thread *) thread::get_curr_thread();
	curr->request_on_vertex(get_id());
	compute_vertex_pointer curr_vertex = curr->get_curr_vertex();
	assert(curr_vertex.is_valid());
	assert(curr_vertex->get_id() == this->get_id());
	directed_vertex_compute *compute
		= (directed_vertex_compute *) curr->get_vertex_compute(curr_vertex);
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

namespace {

class vpart_index_compute: public index_compute
{
	std::vector<vertex_id_t> *large_degree_ids;
	vsize_t num_gets;
	graph_engine &graph;
public:
	vpart_index_compute(index_comp_allocator &_alloc,
			std::vector<vertex_id_t> *_ids,
			graph_engine &_graph): index_compute(_alloc), graph(_graph) {
		this->large_degree_ids = _ids;
		this->num_gets = 0;
	}

	void init(vertex_id_t start, vertex_id_t end) {
		index_compute::init(start);
		add_vertex(end - 1);
	}

	vsize_t get_num_vertices() const {
		return this->get_last_vertex() - this->get_first_vertex() + 1;
	}

	virtual bool run(vertex_id_t start_vid, index_iterator &it);
};

bool vpart_index_compute::run(vertex_id_t start_vid, index_iterator &it)
{
	vertex_id_t vid = start_vid;
	while (it.has_next()) {
		if (graph.is_directed()) {
			vsize_t num_edges = graph.cal_num_edges(it.get_curr_size())
				+ graph.cal_num_edges(it.get_curr_out_size());
			if (num_edges >= (vsize_t) graph_conf.get_min_vpart_degree())
				large_degree_ids->push_back(vid);
		}
		else {
			vsize_t num_edges = graph.cal_num_edges(it.get_curr_size());
			if (num_edges >= (vsize_t) graph_conf.get_min_vpart_degree())
				large_degree_ids->push_back(vid);
		}

		num_gets++;
		vid++;
		it.move_next();
	}
	return num_gets == get_num_vertices();
}

class index_comp_allocator_impl: public index_comp_allocator
{
	class compute_initiator: public obj_initiator<vpart_index_compute>
	{
		index_comp_allocator_impl &alloc;
	public:
		compute_initiator(
				index_comp_allocator_impl &_alloc): alloc(_alloc) {
		}

		virtual void init(vpart_index_compute *obj) {
			new (obj) vpart_index_compute(alloc, alloc.large_degree_ids,
					alloc.graph);
		}
	};

	class compute_destructor: public obj_destructor<vpart_index_compute>
	{
	public:
		void destroy(vpart_index_compute *obj) {
			obj->~vpart_index_compute();
		}
	};

	obj_allocator<vpart_index_compute> allocator;
	std::vector<vertex_id_t> *large_degree_ids;
	graph_engine &graph;
public:
	index_comp_allocator_impl(thread *t, std::vector<vertex_id_t> &ids,
			graph_engine &_graph): allocator("index-compute-allocator",
				t->get_node_id(), false, 1024 * 1024, params.get_max_obj_alloc_size(),
			obj_initiator<vpart_index_compute>::ptr(new compute_initiator(*this)),
			obj_destructor<vpart_index_compute>::ptr(new compute_destructor())),
			graph(_graph) {
		this->large_degree_ids = &ids;
	}

	virtual index_compute *alloc() {
		return allocator.alloc_obj();
	}

	virtual void free(index_compute *obj) {
		allocator.free((vpart_index_compute *) obj);
	}
};

class init_vpart_thread: public thread
{
	graph_index::ptr index;
	graph_engine &graph;
	int hpart_id;
	file_io_factory::shared_ptr io_factory;
public:
	init_vpart_thread(graph_index::ptr index, file_io_factory::shared_ptr io_factory,
			graph_engine &_graph, int hpart_id, int node_id): thread(
				"index-init-thread", node_id), graph(_graph) {
		this->io_factory = io_factory;
		this->index = index;
		this->hpart_id = hpart_id;
	}

	virtual void run();
};

class id_range_t
{
	vertex_id_t start;
	vertex_id_t end;
public:
	id_range_t() {
		start = INVALID_VERTEX_ID;
		end = INVALID_VERTEX_ID;
	}

	vertex_id_t get_start() const {
		return start;
	}

	vertex_id_t get_end() const {
		return end;
	}

	bool add_vertex(vertex_id_t id) {
		if (start == INVALID_VERTEX_ID) {
			start = id;
			end = id + 1;
			return true;
		}
		else if (end == id) {
			end++;
			return true;
		}
		else
			return false;
	}
};

void init_vpart_thread::run()
{
	vertex_index_reader::ptr index_reader;
	if (graph_conf.use_in_mem_index())
		index_reader = vertex_index_reader::create(graph.get_in_mem_index(),
				graph.is_directed());
	else
		index_reader = vertex_index_reader::create(
				io_factory->create_io(this), graph.is_directed());
	std::vector<vertex_id_t> large_degree_ids;
	std::unique_ptr<index_comp_allocator_impl> alloc
		= std::unique_ptr<index_comp_allocator_impl>(
				new index_comp_allocator_impl(this, large_degree_ids, graph));
	const graph_partitioner *partitioner = graph.get_partitioner();
	for (vertex_id_t id = 0; id < index->get_num_vertices(); ) {
		id_range_t range;
		while (id < index->get_num_vertices()) {
			if (partitioner->map(id) == hpart_id) {
				if (range.add_vertex(id))
					id++;
				// If we can't add the vertex, let's stop here
				// and try it in the next iteration of the for loop.
				else
					break;
			}
			else
				id++;
		}
		if (range.get_end() - range.get_start() > 0) {
			vpart_index_compute *compute = (vpart_index_compute *) alloc->alloc();
			compute->init(range.get_start(), range.get_end());
			index_reader->request_index(compute);
			index_reader->wait4complete(1);
		}
		index_reader->wait4complete(0);
		if (index_reader->get_num_pending_tasks() > 100)
			index_reader->wait4complete(1);
	}
	while (index_reader->get_num_pending_tasks() > 0)
		index_reader->wait4complete(1);
	index->init_vparts(hpart_id, graph_conf.get_num_vparts(), large_degree_ids);
	this->stop();
}

}

graph_engine::graph_engine(const std::string &graph_file,
		graph_index::ptr index, const config_map &configs)
{
	struct timeval init_start, init_end;
	gettimeofday(&init_start, NULL);

	init_flash_graph(configs);
	set_file_weight(index->get_index_file(), graph_conf.get_index_file_weight());
	int num_threads = graph_conf.get_num_threads();
	this->num_nodes = params.get_num_nodes();
	index->init(num_threads, num_nodes);

	if (graph_conf.use_in_mem_index()) {
		vertex_index::ptr raw_vindex = vertex_index::safs_load(index->get_index_file());
		assert(raw_vindex->get_graph_header().is_directed_graph());
		vindex = compressed_directed_vertex_index::create(
				(directed_vertex_index &) *raw_vindex);
	}

	max_processing_vertices = graph_conf.get_max_processing_vertices();
	is_complete = false;
	this->vertices = index;

	pthread_mutex_init(&lock, NULL);
	pthread_barrier_init(&barrier1, NULL, num_threads);
	pthread_barrier_init(&barrier2, NULL, num_threads);

	// Right now only the cached I/O can support async I/O
	if (graph_conf.use_in_mem_graph()) {
		graph_data = in_mem_graph::load_graph(graph_file);
		graph_factory = graph_data->create_io_factory();
	}
	else
		graph_factory = create_io_factory(graph_file, GLOBAL_CACHE_ACCESS);
	graph_factory->set_sched_creater(new throughput_comp_io_sched_creater());
	index_factory = create_io_factory(index->get_index_file(), GLOBAL_CACHE_ACCESS);

	io_interface::ptr io = index_factory->create_io(thread::get_curr_thread());
	io->access((char *) &header, 0, sizeof(header), READ);
	header.verify();
	out_part_off = 0;
	if (header.is_directed_graph()) {
		assert(sizeof(vertex_index) == sizeof(header));
		vertex_index *idx = (vertex_index *) &header;
		out_part_off = idx->get_out_part_loc();
	}

	assert(num_threads > 0 && num_nodes > 0);
	assert(num_threads % num_nodes == 0);
	worker_threads.resize(num_threads);
	vprograms.resize(num_threads);

	if (!graph_conf.get_trace_file().empty())
		logger = trace_logger::ptr(new trace_logger(graph_conf.get_trace_file()));

	if (graph_conf.preload())
		preload_graph();

	// If we need to perform vertical partitioning on the graph.
	if (graph_conf.get_num_vparts() > 1) {
		std::vector<init_vpart_thread *> threads(num_threads);
		for (int i = 0; i < num_threads; i++) {
			threads[i] = new init_vpart_thread(index, index_factory, *this,
					i, i % num_nodes);
			threads[i]->start();
		}
		for (int i = 0; i < num_threads; i++) {
			threads[i]->join();
			delete threads[i];
		}
	}

	gettimeofday(&init_end, NULL);
	printf("It takes %f seconds to initialize the graph engine\n",
			time_diff(init_start, init_end));
}

graph_engine::~graph_engine()
{
	graph_factory->print_statistics();
	index_factory->print_statistics();
	for (unsigned i = 0; i < worker_threads.size(); i++)
		delete worker_threads[i];
	graph_factory = file_io_factory::shared_ptr();
	index_factory = file_io_factory::shared_ptr();
	destroy_flash_graph();
}

static inline int get_node_id(int thread_idx, int num_nodes)
{
	return thread_idx % num_nodes;
}

void graph_engine::init_threads(vertex_program_creater::ptr creater)
{
	std::vector<std::shared_ptr<slab_allocator> > msg_allocs(num_nodes);
	std::vector<std::shared_ptr<slab_allocator> > flush_msg_allocs(num_nodes);
	// It turns out that it's important to respect the NUMA effect here.
	// When a worker thread uses the message allocator from the same NUMA node,
	// we can get noticeably better performance.
	for (int node_id = 0; node_id < num_nodes; node_id++) {
		// We increase the allocator by 1M each time.
		// It shouldn't need to allocate much memory.
		msg_allocs[node_id] = std::shared_ptr<slab_allocator>(
				new slab_allocator("graph-message-allocator",
					GRAPH_MSG_BUF_SIZE, 1024 * 1024, INT_MAX, node_id,
					false /* init */, false /* pinned */, 5 /* local_buf_size*/));
		flush_msg_allocs[node_id] = std::shared_ptr<slab_allocator>(
				new slab_allocator("graph-message-allocator",
					MAX_FLUSH_MSG_SIZE, 1024 * 1024, INT_MAX, node_id,
					false /* init */, false /* pinned */, 20 /* local_buf_size*/));
	}
	// Prepare the worker threads.
	int num_threads = get_num_threads();
	for (int i = 0; i < num_threads; i++) {
		vertex_program::ptr new_prog;
		if (creater)
			new_prog = creater->create();
		else
			new_prog = vertices->create_def_vertex_program();
		worker_thread *t = new worker_thread(this, graph_factory, index_factory,
				new_prog, vertices->create_def_part_vertex_program(),
				get_node_id(i, num_nodes), i, num_threads, scheduler,
				msg_allocs[get_node_id(i, num_nodes)]);
		assert(worker_threads[i] == NULL);
		worker_threads[i] = t;
		vprograms[i] = new_prog;
	}
	for (int i = 0; i < num_threads; i++) {
		worker_threads[i]->init_messaging(worker_threads,
				msg_allocs[get_node_id(i, num_nodes)],
				flush_msg_allocs[get_node_id(i, num_nodes)]);
	}
}

void graph_engine::start(const vertex_id_t ids[], int num,
		vertex_initializer::ptr init, vertex_program_creater::ptr creater)
{
	gettimeofday(&start_time, NULL);
	init_threads(std::move(creater));
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
	iter_start = start_time;
}

void graph_engine::start(std::shared_ptr<vertex_filter> filter,
		vertex_program_creater::ptr creater)
{
	gettimeofday(&start_time, NULL);
	init_threads(std::move(creater));
	// Let's assume all vertices will be activated first.
	BOOST_FOREACH(worker_thread *t, worker_threads) {
		t->start_vertices(filter);
		t->start();
	}
	iter_start = start_time;
}

void graph_engine::start_all(vertex_initializer::ptr init,
		vertex_program_creater::ptr creater)
{
	gettimeofday(&start_time, NULL);
	init_threads(std::move(creater));
	BOOST_FOREACH(worker_thread *t, worker_threads) {
		t->start_all_vertices(init);
		t->start();
	}
	iter_start = start_time;
}

bool graph_engine::progress_first_level()
{
	static atomic_number<long> tot_num_activates;
	static atomic_integer num_threads;

	worker_thread *curr = (worker_thread *) thread::get_curr_thread();
	int num_activates = curr->get_activates();
	tot_num_activates.inc(num_activates);
	// If all threads have reached here.
	if (num_threads.inc(1) == get_num_threads()) {
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
	int rc = pthread_barrier_wait(&barrier2);
	if(rc != 0 && rc != PTHREAD_BARRIER_SERIAL_THREAD)
	{
		printf("Could not wait on barrier\n");
		exit(-1);
	}
	return is_complete;
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
				level.get() - 1, time_diff(iter_start, curr),
				tot_num_activates.get(), level.get());
		iter_start = curr;
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
	struct timeval curr;
	gettimeofday(&curr, NULL);
	printf("The graph engine takes %f seconds to complete\n",
			time_diff(start_time, curr));
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

std::atomic<size_t> graph_engine::init_count;

void graph_engine::init_flash_graph(const config_map &configs)
{
	size_t count = init_count.fetch_add(1);
	if (count == 0) {
		graph_conf.init(configs);
		graph_conf.print();
		init_io_system(configs);
	}
}

void graph_engine::destroy_flash_graph()
{
	size_t count = init_count.fetch_sub(1);
	if (count == 1) {
		destroy_io_system();
	}
}
