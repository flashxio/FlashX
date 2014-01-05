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
#include "slab_allocator.h"

#include "graph_config.h"
#include "graph_engine.h"

const int MAX_IO_PEND_VERTICES = 1000;

graph_config graph_conf;

/**
 * This callback is to process a vertex.
 */
class vertex_compute: public user_compute
{
	graph_engine *graph;
	compute_vertex *v;
public:
	vertex_compute(graph_engine *graph,
			compute_allocator *alloc): user_compute(alloc) {
		this->graph = graph;
		v = NULL;
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
		vertex_id_t id = v->get_next_required_vertex();
		compute_vertex &info = graph->get_vertex(id);
		data_loc_t loc(graph->get_file_id(), info.get_ext_mem_off());
		return request_range(loc, info.get_ext_mem_size(), READ);
	}

	virtual bool run(page_byte_array &);
};

class vertex_compute_allocator: public compute_allocator
{
	class compute_initiator: public obj_initiator<vertex_compute>
	{
		graph_engine *graph;
		vertex_compute_allocator *alloc;
	public:
		compute_initiator(graph_engine *graph, vertex_compute_allocator *alloc) {
			this->graph = graph;
			this->alloc = alloc;
		}

		virtual void init(vertex_compute *obj) {
			new (obj) vertex_compute(graph, alloc);
		}
	};

	obj_allocator<vertex_compute> allocator;
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
		allocator.free((vertex_compute *) obj);
	}
};

class worker_thread: public thread
{
	file_io_factory *factory;
	io_interface *io;
	graph_engine *graph;
	compute_allocator *alloc;
	std::vector<vertex_id_t> activated_vertices;
public:
	worker_thread(graph_engine *graph, file_io_factory *factory,
			int node_id): thread("worker_thread", node_id) {
		this->graph = graph;
		this->io = NULL;
		this->factory = factory;
		alloc = new vertex_compute_allocator(graph, this);
	}

	~worker_thread() {
		delete alloc;
	}

	void run();
	void init();

	std::vector<vertex_id_t> &get_activated_vertices() {
		return activated_vertices;
	}

	int process_activated_vertices(int max);
};

bool vertex_compute::run(page_byte_array &array)
{
	char buf[STACK_PAGE_VERTEX_SIZE];
	const page_vertex *ext_v;
	if (graph->is_directed())
		ext_v = new (buf) page_directed_vertex(array);
	else
		ext_v = new (buf) page_undirected_vertex(array);
	// If the algorithm doesn't need to get the full information
	// of their neighbors
	if (graph->get_required_neighbor_type() == edge_type::NONE
			// Or we haven't perform computation on the vertex yet.
			|| v == NULL) {
		v = &graph->get_vertex(ext_v->get_id());
		v->run(*graph, ext_v);
	}
	else {
		v->run_on_neighbors(*graph, &ext_v, 1);
	}
	return true;
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

	void init(std::vector<vertex_id_t> *vecs[], int num_vecs) {
		// all vertices have been sorted in each vector, we only need to
		// merge them.
		fetch_idx = 0;
		sorted_vertices.clear();
#ifndef MEMCHECK
		std::vector<std::pair<std::vector<vertex_id_t>::iterator,
			std::vector<vertex_id_t>::iterator> > seqs;
		size_t tot_length = 0;
		for (int i = 0; i < num_vecs; i++) {
			seqs.push_back(std::make_pair<std::vector<vertex_id_t>::iterator,
					std::vector<vertex_id_t>::iterator>(
						vecs[i]->begin(), vecs[i]->end()));
			tot_length += vecs[i]->size();
		}
		sorted_vertices.resize(tot_length);
		__gnu_parallel::multiway_merge(seqs.begin(), seqs.end(),
				sorted_vertices.begin(), tot_length, std::less<int>());
#else
		for (int i = 0; i < num_vecs; i++)
			sorted_vertices.insert(sorted_vertices.end(), vecs[i]->begin(),
					vecs[i]->end());
		std::sort(sorted_vertices.begin(), sorted_vertices.end());
#endif
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
public:
	vertex_collection(const std::vector<thread *> &threads) {
		this->threads = threads;
	}

	void add(vertex_id_t vertices[], int num) {
		std::vector<vertex_id_t> &vec = ((worker_thread *) thread::get_curr_thread(
					))->get_activated_vertices();
		vec.insert(vec.end(), vertices, vertices + num);
	}

	void local_sort() {
		std::vector<vertex_id_t> &vec = ((worker_thread *) thread::get_curr_thread(
					))->get_activated_vertices();
		std::sort(vec.begin(), vec.end());
	}

	void sort(sorted_vertex_queue &sorted_vertices) const {
		std::vector<vertex_id_t> *vecs[threads.size()];
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

void worker_thread::init()
{
	io = factory->create_io(this);
	io->init();
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
	int num = graph->get_curr_activated_vertices(vertex_buf, max);
	for (int i = 0; i < num; i++) {
		compute_vertex &info = graph->get_vertex(vertex_buf[i]);
		data_loc_t loc(io->get_file_id(), info.get_ext_mem_off());
		reqs[i] = io_request(alloc->alloc(), loc,
				// TODO I might need to set the node id.
				info.get_ext_mem_size(), READ, io, -1);
	}
	if (graph->get_logger())
		graph->get_logger()->log(reqs.data(), num);
	io->access(reqs.data(), num);
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
		printf("thread %d process %d activated vertices\n", get_id(), num);
		num_visited += num;
		while (graph->get_num_curr_activated_vertices() > 0) {
			num = process_activated_vertices(
					MAX_IO_PEND_VERTICES - io->num_pending_ios());
			num_visited += num;
			io->wait4complete(io->num_pending_ios() / 10 + max<int>(0,
						io->num_pending_ios() - MAX_IO_PEND_VERTICES));
		}
		// We have completed processing the activated vertices in this iteration.
		while (io->num_pending_ios() > 0) {
			io->wait4complete(1);
		}
		printf("thread %d visited %d vertices\n", this->get_id(), num_visited);

		// Now we have finished this level, we can progress to the next level.
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
		const std::string &graph_file, graph_index *index, bool directed)
{
	this->required_neighbor_type = edge_type::NONE;
	this->directed = directed;
	is_complete = false;
	this->vertices = index;

	pthread_mutex_init(&lock, NULL);
	pthread_barrier_init(&barrier1, NULL, num_threads);
	pthread_barrier_init(&barrier2, NULL, num_threads);

	file_io_factory *factory = create_io_factory(graph_file,
			GLOBAL_CACHE_ACCESS);
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
	vertex_id_t to_add[num];
	int num_to_add = 0;

	for (int i = 0; i < num; i++) {
		compute_vertex &v = get_vertex(vertices[i]);
		// When a vertex is added to the queue, we mark it as visited.
		// Therefore, a vertex can only be added to the queue once and
		// can only be visited once.
		if (v.activate_in(level.get()))
			to_add[num_to_add++] = vertices[i];
	}
	activated_vertex_buf->add(to_add, num_to_add);
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
	activated_vertex_buf->local_sort();
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
