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

const int MAX_IO_PEND_VERTICES = 2000;

graph_config graph_conf;

class worker_thread;

request_range compute_vertex::get_next_request(graph_engine *graph)
{
	vertex_id_t id = get_next_required_vertex();
	compute_vertex &info = graph->get_vertex(id);
	data_loc_t loc(graph->get_file_id(), info.get_ext_mem_off());
	return request_range(loc, info.get_ext_mem_size(), READ, NULL);
}

void ts_vertex_request::set_vertex(vertex_id_t id)
{
	this->id = id;
	compute_vertex &info = graph->get_vertex(id);
	// There is some overhead to fetch part of a vertex, so we should
	// minize the number of vertices fetched partially.
	// If a vertex is small enough (stored on <= 3 pages), we fetch the entire
	// vertex.
	off_t start_pg = ROUND_PAGE(info.get_ext_mem_off());
	off_t end_pg = ROUNDUP_PAGE(info.get_ext_mem_off() + info.get_ext_mem_size());
	if (end_pg - start_pg <= PAGE_SIZE * 3)
		require_all = true;
}

/**
 * This callback is to process a vertex.
 */
class vertex_compute: public user_compute
{
	// The number of requested vertices that will be read in the user compute.
	int num_complete_issues;
	// The number of vertices read by the user compute.
	int num_complete_fetched;

	graph_engine *graph;
	compute_vertex *v;
	// The thread that creates the vertex compute.
	worker_thread *issue_thread;
public:
	vertex_compute(graph_engine *graph,
			compute_allocator *alloc): user_compute(alloc) {
		this->graph = graph;
		v = NULL;
		issue_thread = (worker_thread *) thread::get_curr_thread();
		num_complete_issues = 0;
		num_complete_fetched = 0;
	}

	virtual int serialize(char *buf, int size) const {
		return 0;
	}

	virtual int get_serialized_size() const {
		return 0;
	}

	virtual int has_requests() const {
		if (v == NULL)
			return false;
		else
			return v->has_required_vertices();
	}

	virtual request_range get_next_request() {
		assert(v);
		request_range range = v->get_next_request(graph);
		if (range.get_compute() == NULL) {
			num_complete_issues++;
			range.set_compute(this);
		}
		return range;
	}

	virtual void run(page_byte_array &);

	virtual bool has_completed() const {
		// If the user compute has got all requested data and it has
		// no more requests to issue, we can consider the user compute
		// has been completed.
		// NOTE: it's possible that requested data may not be passed to
		// this user compute, so we only count the requests that are going
		// to be passed to this user compute.
		return num_complete_issues == num_complete_fetched && !has_requests();
	}
};

/**
 * Sometimes a vertex only needs to read part of its neighbors.
 * This class is to read part of a neighbor and pass the neighbor to
 * the specified vertex to perform computation.
 * An instance of the class only reads one neighbor.
 */
class part_ts_vertex_compute: public user_compute
{
	graph_engine *graph;
	// The vertex where computation should perform.
	compute_vertex *comp_v;
	const TS_page_vertex *required_vertex_header;
	// The part of the vertex will be read and passed to
	// the computation vertex.
	ts_vertex_request required_part;
	// The thread that creates the vertex compute.
	worker_thread *issue_thread;
	int num_issued;
	int num_fetched;
public:
	part_ts_vertex_compute(graph_engine *graph,
			compute_allocator *alloc): user_compute(alloc), required_part(graph) {
		this->graph = graph;
		comp_v = NULL;
		issue_thread = (worker_thread *) thread::get_curr_thread();
		required_vertex_header = NULL;
		num_issued = 0;
		num_fetched = 0;
	}

	void init(compute_vertex *v, const ts_vertex_request &req) {
		comp_v = v;
		required_part = req;
	}

	virtual int serialize(char *buf, int size) const {
		return 0;
	}

	virtual int get_serialized_size() const {
		return 0;
	}

	virtual int has_requests() const {
		return num_issued == 0;
	}

	virtual request_range get_next_request();

	virtual void run(page_byte_array &);

	virtual bool has_completed() const {
		return num_fetched > 0;
	}
};

template<class compute_type>
class vertex_compute_allocator: public compute_allocator
{
	class compute_initiator: public obj_initiator<compute_type>
	{
		graph_engine *graph;
		vertex_compute_allocator<compute_type> *alloc;
	public:
		compute_initiator(graph_engine *graph,
				vertex_compute_allocator<compute_type> *alloc) {
			this->graph = graph;
			this->alloc = alloc;
		}

		virtual void init(compute_type *obj) {
			new (obj) compute_type(graph, alloc);
		}
	};

	obj_allocator<compute_type> allocator;
public:
	vertex_compute_allocator(graph_engine *graph, thread *t): allocator(
			"vertex-compute-allocator", t->get_node_id(), 1024 * 1024,
			params.get_max_obj_alloc_size(),
			// TODO memory leak here
			new compute_initiator(graph, this)) {
	}

	virtual user_compute *alloc() {
		return allocator.alloc_obj();
	}

	virtual void free(user_compute *obj) {
		allocator.free((compute_type *) obj);
	}
};

class worker_thread: public thread
{
	file_io_factory *factory;
	io_interface *io;
	graph_engine *graph;
	compute_allocator *alloc;
	compute_allocator *part_alloc;
	// This is to collect vertices activated in the next level.
	bitmap activated_vertices;

	// The number of activated vertices processed in the current level.
	atomic_number<long> num_activated_vertices_in_level;
	// The number of vertices completed in the current level.
	atomic_number<long> num_completed_vertices_in_level;
public:
	worker_thread(graph_engine *graph, file_io_factory *factory,
			int node_id): thread("worker_thread", node_id), activated_vertices(
				graph->get_max_vertex_id() + 1) {
		this->graph = graph;
		this->io = NULL;
		this->factory = factory;
		alloc = new vertex_compute_allocator<vertex_compute>(graph, this);
		// TODO we need to fix this.
		part_alloc = graph->create_part_compute_allocator(this);
	}

	~worker_thread() {
		delete alloc;
		graph->destroy_part_compute_allocator(part_alloc);
		factory->destroy_io(io);
	}

	void run();
	void init() {
		io = factory->create_io(this);
		io->init();
	}

	compute_allocator *get_part_compute_allocator() const {
		return part_alloc;
	}

	bitmap &get_activated_vertices() {
		return activated_vertices;
	}

	int process_activated_vertices(int max);

	void complete_vertex(const compute_vertex &v) {
		num_completed_vertices_in_level.inc(1);
	}
};

request_range ts_compute_vertex::get_next_request(graph_engine *graph)
{
	ts_vertex_request ts_req(graph);
	get_next_required_ts_vertex(ts_req);
	assert(ts_req.get_edge_type() == edge_type::BOTH_EDGES);

	compute_vertex &info = graph->get_vertex(ts_req.get_id());
	data_loc_t loc(graph->get_file_id(), info.get_ext_mem_off());
	if (ts_req.is_require_all()) {
		return request_range(loc, info.get_ext_mem_size(), READ, NULL);
	}
	else {
		worker_thread *t = (worker_thread *) thread::get_curr_thread();
		compute_allocator *alloc = t->get_part_compute_allocator();
		assert(alloc);
		part_ts_vertex_compute *comp = (part_ts_vertex_compute *) alloc->alloc();
		comp->init(this, ts_req);
		// We assume the header of a ts-vertex is never larger than a page.
		return request_range(loc, PAGE_SIZE, READ, comp);
	}
}

bool ts_compute_vertex::run_on_neighbors(graph_engine &graph,
		const page_vertex *vertices[], int num)
{
	const TS_page_vertex **ts_vertices = (const TS_page_vertex **) vertices;
	return run_on_neighbors(graph, ts_vertices, num);
}

compute_allocator *ts_compute_vertex::create_part_compute_allocator(
		graph_engine *graph, thread *t)
{
	return new vertex_compute_allocator<part_ts_vertex_compute>(
			graph, t);
}

request_range part_ts_vertex_compute::get_next_request()
{
	assert(required_vertex_header);
	assert(num_issued == 0);
	num_issued++;

	compute_vertex &info = graph->get_vertex(
			required_vertex_header->get_id());
	offset_pair rel_offsets = required_vertex_header->get_edge_list_offset(
			required_part.get_range());
	data_loc_t loc(graph->get_file_id(), rel_offsets.first
			+ info.get_ext_mem_off());
	return request_range(loc, rel_offsets.second - rel_offsets.first,
			READ, this);
}

void part_ts_vertex_compute::run(page_byte_array &array)
{
	assert(!has_completed());
	bool completed = false;
	if (required_vertex_header == NULL) {
		ext_mem_vertex_interpreter *interpreter = graph->get_vertex_interpreter();
		char *buf = new char[interpreter->get_vertex_size()];
		required_vertex_header = (const TS_page_vertex *) interpreter->interpret(
				array, buf, interpreter->get_vertex_size());
		assert(!required_vertex_header->is_complete());
	}
	else {
		ext_mem_vertex_interpreter *interpreter = graph->get_vertex_interpreter();
		stack_array<char, 64> buf(interpreter->get_vertex_size());
		const page_vertex *ext_v = interpreter->interpret_part(
				required_vertex_header, array, buf.data(),
				interpreter->get_vertex_size());

		num_fetched++;
		assert(comp_v);
		completed = comp_v->run_on_neighbors(*graph, &ext_v, 1);

		char *tmp = (char *) required_vertex_header;
		delete [] tmp;
	}
	// We need to notify the thread that initiate processing the vertex
	// of the completion of the vertex.
	if (completed)
		issue_thread->complete_vertex(*comp_v);
}

void vertex_compute::run(page_byte_array &array)
{
	ext_mem_vertex_interpreter *interpreter = graph->get_vertex_interpreter();
	stack_array<char, 64> buf(interpreter->get_vertex_size());
	const page_vertex *ext_v = interpreter->interpret(array, buf.data(),
			interpreter->get_vertex_size());
	bool completed = false;
	// If the algorithm doesn't need to get the full information
	// of their neighbors
	if (graph->get_required_neighbor_type() == edge_type::NONE
			// Or we haven't perform computation on the vertex yet.
			|| v == NULL) {
		v = &graph->get_vertex(ext_v->get_id());
		completed = v->run(*graph, ext_v);
	}
	else {
		num_complete_fetched++;
		completed = v->run_on_neighbors(*graph, &ext_v, 1);
	}
	// We need to notify the thread that initiate processing the vertex
	// of the completion of the vertex.
	if (completed) {
		// If the vertex has completed computation, the user compute should
		// also finish its computation.
		assert(has_completed());
		issue_thread->complete_vertex(*v);
	}
}

class sorted_vertex_queue
{
	pthread_spinlock_t lock;
	std::vector<vertex_id_t> sorted_vertices;
	size_t fetch_idx;
public:
	sorted_vertex_queue() {
		fetch_idx = 0;
		pthread_spin_init(&lock, PTHREAD_PROCESS_PRIVATE);
	}

	void init(const std::vector<vertex_id_t> &vec, bool sorted) {
		fetch_idx = 0;
		sorted_vertices.clear();
		sorted_vertices.assign(vec.begin(), vec.end());
		if (!sorted)
			std::sort(sorted_vertices.begin(), sorted_vertices.end());
	}

	void init(bitmap *vecs[], int num_vecs) {
		fetch_idx = 0;
		sorted_vertices.clear();
		if (num_vecs == 0)
			return;

		size_t size = vecs[0]->get_num_bits();
		for (int i = 1; i < num_vecs; i++)
			assert(size == vecs[i]->get_num_bits());

		bitmap tmp(size);
		for (int i = 0; i < num_vecs; i++) {
			tmp.merge(*vecs[i]);
		}
		tmp.get_set_bits(sorted_vertices);
	}

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
};

/**
 * This class collects vertices added by each thread.
 */
class vertex_collection
{
	std::vector<thread *> threads;
	vertex_id_t max_id;
public:
	vertex_collection(const std::vector<thread *> &threads) {
		this->threads = threads;
		this->max_id = max_id;
	}

	void add(vertex_id_t vertices[], int num) {
		bitmap &vec = ((worker_thread *) thread::get_curr_thread(
					))->get_activated_vertices();
		for (int i = 0; i < num; i++)
			vec.set(vertices[i]);
	}

	void sort(sorted_vertex_queue &sorted_vertices) const {
		bitmap *vecs[threads.size()];
		for (unsigned i = 0; i < threads.size(); i++) {
			vecs[i] = &((worker_thread *) threads[i])->get_activated_vertices();
		}
		sorted_vertices.init(vecs, threads.size());
	}

	void clear() {
		for (unsigned i = 0; i < threads.size(); i++)
			((worker_thread *) threads[i])->get_activated_vertices().clear();
	}
};

/**
 * This is to process the activated vertices in the current iteration.
 */
int worker_thread::process_activated_vertices(int max)
{
	if (max <= 0)
		return 0;

	vertex_id_t vertex_buf[max];
	stack_array<io_request> reqs(max);
	int num = graph->get_curr_activated_vertices(vertex_buf, max);
	num_activated_vertices_in_level.inc(num);
	int num_completed = 0;
	int num_to_process = 0;
	for (int i = 0; i < num; i++) {
		compute_vertex &info = graph->get_vertex(vertex_buf[i]);
		// We need to execute the pre-run to determine if we should
		// fetch the adjacency list of itself.
		if (info.run(*graph)) {
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

/**
 * This method is the main function of the graph engine.
 */
void worker_thread::run()
{
	while (true) {
		int num_visited = 0;
		int num = process_activated_vertices(MAX_IO_PEND_VERTICES);
		num_visited += num;
		while (graph->get_num_curr_activated_vertices() > 0) {
			num = process_activated_vertices(
					MAX_IO_PEND_VERTICES - io->num_pending_ios());
			num_visited += num;
			io->wait4complete(io->num_pending_ios() / 10 + max<int>(0,
						io->num_pending_ios() - MAX_IO_PEND_VERTICES));
		}
		// We have completed processing the activated vertices in this iteration.
		while (num_activated_vertices_in_level.get()
				- num_completed_vertices_in_level.get() > 0) {
			io->wait4complete(1);
		}
		printf("thread %d visited %d vertices\n", this->get_id(), num_visited);

		// Now we have finished this level, we can progress to the next level.
		num_activated_vertices_in_level = atomic_number<long>(0);
		num_completed_vertices_in_level = atomic_number<long>(0);
		bool completed = graph->progress_next_level();
		printf("thread %d finish in a level, completed? %d\n", get_id(), completed);
		if (completed)
			break;
	}
	stop();
	if (graph_conf.get_print_io_stat())
		io->print_stat(graph->get_num_threads());
}

graph_engine::graph_engine(int num_threads, int num_nodes,
		const std::string &graph_file, graph_index *index,
		ext_mem_vertex_interpreter *interpreter, bool directed)
{
	this->interpreter = interpreter;
	this->required_neighbor_type = edge_type::NONE;
	this->directed = directed;
	is_complete = false;
	this->vertices = index;

	pthread_mutex_init(&lock, NULL);
	pthread_barrier_init(&barrier1, NULL, num_threads);
	pthread_barrier_init(&barrier2, NULL, num_threads);

	file_io_factory *factory = create_io_factory(graph_file,
			GLOBAL_CACHE_ACCESS);

	graph_header header;
	io_interface *io = factory->create_io(thread::get_curr_thread());
	io->access((char *) &header, 0, sizeof(header), READ);
	header.verify();
	factory->destroy_io(io);

	file_id = factory->get_file_id();
	assert(num_threads > 0 && num_nodes > 0);
	assert(num_threads % num_nodes == 0);
	for (int i = 0; i < num_threads; i++) {
		worker_thread *t = new worker_thread(this, factory,
				i % num_nodes);
		worker_threads.push_back(t);
	}
	first_thread = worker_threads[0];

	activated_vertices = new sorted_vertex_queue();
	activated_vertex_buf = new vertex_collection(worker_threads);

	if (graph_conf.get_trace_file().empty())
		logger = NULL;
	else
		logger = new trace_logger(graph_conf.get_trace_file());
}

graph_engine::~graph_engine()
{
	delete activated_vertices;
	delete activated_vertex_buf;
	for (unsigned i = 0; i < worker_threads.size(); i++)
		delete worker_threads[i];
	if (logger)
		delete logger;
}

void graph_engine::start(vertex_id_t ids[], int num)
{
	std::vector<vertex_id_t> starts;
	starts.assign(ids, ids + num);
	activated_vertices->init(starts, false);
	for (unsigned i = 0; i < worker_threads.size(); i++)
		worker_threads[i]->start();
}

void graph_engine::start_all()
{
	std::vector<vertex_id_t> all_vertices;
	vertices->get_all_vertices(all_vertices);
	printf("there are %ld activated vertices\n", all_vertices.size());
	activated_vertices->init(all_vertices, true);
	for (unsigned i = 0; i < worker_threads.size(); i++)
		worker_threads[i]->start();
}

void graph_engine::activate_vertices(vertex_id_t vertices[], int num)
{
	activated_vertex_buf->add(vertices, num);
}

int graph_engine::get_curr_activated_vertices(vertex_id_t vertices[], int num)
{
	return activated_vertices->fetch(vertices, num);
}

size_t graph_engine::get_num_curr_activated_vertices() const
{
	return activated_vertices->get_num_vertices();
}

bool graph_engine::progress_next_level()
{
	// We have to make sure all threads have reach here, so we can switch
	// queues to progress to the next level.
	// If the queue of the next level is empty, the program can terminate.
	int rc = pthread_barrier_wait(&barrier1);
	if(rc != 0 && rc != PTHREAD_BARRIER_SERIAL_THREAD)
	{
		printf("Could not wait on barrier\n");
		exit(-1);
	}
	pthread_mutex_lock(&lock);
	if (thread::get_curr_thread() == first_thread) {
		assert(activated_vertices->is_empty());
		activated_vertex_buf->sort(*activated_vertices);
		activated_vertex_buf->clear();
		level.inc(1);
		printf("progress to level %d, there are %ld vertices in this level\n",
				level.get(), activated_vertices->get_num_vertices());
		is_complete = activated_vertices->is_empty();
	}
	pthread_mutex_unlock(&lock);

	// We need to synchronize again. Only one thread can switch the queues.
	// We have to make sure that a thread checks the queue of the current level
	// after the first thread has switched the queues.
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
	}
}
