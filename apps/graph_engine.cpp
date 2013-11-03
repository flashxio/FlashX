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

#include "graph_config.h"
#include "graph_engine.h"

const int MAX_IO_PEND_VERTICES = 1000;

graph_config graph_conf;

/**
 * This callback is to process a vertex.
 */
class vertex_callback: public callback
{
	graph_engine *graph;
	io_interface *io;
public:
	vertex_callback(graph_engine *graph, io_interface *io) {
		this->graph = graph;
		this->io = io;
	}

	int invoke(io_request *reqs[], int num);
};

class pending_vertex: public ext_mem_vertex
{
	// The number of neighbors that have been fetched from disks.
	int num_completed_neighbors;
	edge_type required_neighbor_type;

	pending_vertex(const ext_mem_vertex &v,
			edge_type required_neighbor_type): ext_mem_vertex(v) {
		num_completed_neighbors = 0;
		this->required_neighbor_type = required_neighbor_type;
	}

	pending_vertex(char *buf, int size, bool directed,
			edge_type required_neighbor_type): ext_mem_vertex(buf, size, directed) {
		num_completed_neighbors = 0;
		this->required_neighbor_type = required_neighbor_type;
	}
public:
	static pending_vertex *create(char *buf, int size, graph_engine *graph) {
		return new pending_vertex(buf, size, graph->is_directed(),
				graph->get_required_neighbor_type());
	}
	static pending_vertex *create(const ext_mem_vertex &v, graph_engine *graph) {
		return new pending_vertex(v, graph->get_required_neighbor_type());
	}

	static void destroy(pending_vertex *v) {
		delete [] v->get_buf();
		delete v;
	}

	void complete_neighbor() {
		num_completed_neighbors++;
	}

	bool is_complete() const {
		return num_completed_neighbors == get_num_edges(required_neighbor_type);
	}

	int get_num_neighbors() const {
		return get_num_edges(required_neighbor_type);
	}

	int fetch_neighbors(fifo_queue<vertex_id_t> &neighbors) const {
		return ext_mem_vertex::get_neighbors(required_neighbor_type, neighbors);
	}
};

class pending_neighbor_collection
{
	fifo_queue<vertex_id_t> neighbors;
	pending_vertex *vertex;
public:
	pending_neighbor_collection(): neighbors(-1, 1024, true) {
		this->vertex = NULL;
	}

	void init(pending_vertex *vertex) {
		assert(neighbors.is_empty());
		this->vertex = vertex;
		if (neighbors.get_size() < vertex->get_num_neighbors())
			neighbors.expand_queue(vertex->get_num_neighbors());
		assert(neighbors.get_size() >= vertex->get_num_neighbors());
		vertex->fetch_neighbors(neighbors);
	}

	bool is_empty() {
		return get_num_neighbors() == 0;
	}

	int get_num_neighbors() {
		return neighbors.get_num_entries();
	}

	int fetch_neighbors(vertex_id_t neighbors[], int num) {
		return this->neighbors.fetch(neighbors, num);
	}

	pending_vertex *get_pending_vertex() const {
		return vertex;
	}
};

class worker_thread: public thread
{
	file_io_factory *factory;
	io_interface *io;
	graph_engine *graph;
	fifo_queue<ext_mem_vertex> pending_vertices;

	pending_neighbor_collection curr_pending;

	std::vector<vertex_id_t> activated_vertices;
public:
	worker_thread(graph_engine *graph, file_io_factory *factory,
			int node_id): thread("worker_thread", node_id), pending_vertices(
				-1, 1024, true) {
		this->graph = graph;
		this->io = NULL;
		this->factory = factory;
	}

	void run();
	void init();

	std::vector<vertex_id_t> &get_activated_vertices() {
		return activated_vertices;
	}

	void add_pending_vertex(ext_mem_vertex v) {
		if (pending_vertices.is_full())
			pending_vertices.expand_queue(pending_vertices.get_size() * 2);
		pending_vertices.push_back(v);
	}

	int process_pending_vertex(int max);
	int process_activated_vertices(int max);
	int process_pending_vertices(int max);
};

int vertex_callback::invoke(io_request *reqs[], int num)
{
	worker_thread *curr_thread = (worker_thread *) thread::get_curr_thread();
	for (int i = 0; i < num; i++) {
		char *req_buf = reqs[i]->get_buf();
		size_t req_size = reqs[i]->get_size();
		assert(this->io == reqs[i]->get_io());
		ext_mem_vertex ext_v(req_buf, req_size, graph->is_directed());
		// If the algorithm doesn't need to get the full information
		// of their neighbors
		if (graph->get_required_neighbor_type() == edge_type::NONE
				// Or the vertex doesn't have neighbors.
				|| (reqs[i]->get_user_data() == NULL
					&& ext_v.get_num_edges(graph->get_required_neighbor_type()) == 0)) {
			// We can run user's code immediately on the vertex.
			compute_vertex &v = graph->get_vertex(ext_v.get_id());
			v.materialize(ext_v);
			v.run(*graph, NULL, 0);
			v.dematerialize();
			delete [] req_buf;
		}
		// We just fetched a vertex, we need to fetch its neighbors to
		// perform computation.
		else if (reqs[i]->get_user_data() == NULL) {
			curr_thread->add_pending_vertex(ext_v);
		}
		else {
			// Now a neighbor has been fetched, now we can do some computation
			// between the original pending vertex and its neighbor.
			pending_vertex *pending
				= (pending_vertex *) reqs[i]->get_user_data();
			compute_vertex &v = graph->get_vertex(
					pending->get_id());
			// We materialize the vertex and perform computation.
			// The callback function is guaranteed to be called in the thread
			// where a request is issued. Since all requests of fetching
			// neighbors are issued by one thread, we don't need to use a lock
			// to protect the pending vertex from concurrent access.
			v.materialize(*pending);
			v.run(*graph, &ext_v, 1);
			v.dematerialize();
			// The buffer contains the info of the neighbor, no we don't need
			// it any more.
			delete [] req_buf;
			pending->complete_neighbor();
			// Once we perform computation on all neighbors. We can destroy
			// the pending vertex.
			if (pending->is_complete())
				pending_vertex::destroy(pending);
		}
	}
	return 0;
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
	io->set_callback(new vertex_callback(graph, io));
}

int worker_thread::process_pending_vertex(int max)
{
	int num_neighbors = min(max, curr_pending.get_num_neighbors());
	if (num_neighbors <= 0)
		return 0;

	stack_array<io_request> reqs(num_neighbors);
	stack_array<vertex_id_t> remain_neighs(num_neighbors);

	int ret = curr_pending.fetch_neighbors(remain_neighs.data(),
			num_neighbors);
	assert(ret == num_neighbors);
	for (int j = 0; j < num_neighbors; j++) {
		vertex_id_t neighbor = remain_neighs[j];
		compute_vertex &info = graph->get_vertex(neighbor);
		reqs[j].init(new char[info.get_ext_mem_size()],
				info.get_ext_mem_off(),
				// TODO I might need to set the node id.
				info.get_ext_mem_size(), READ, io, -1);
		reqs[j].set_user_data(curr_pending.get_pending_vertex());
	}
	if (graph->get_logger())
		graph->get_logger()->log(reqs.data(), num_neighbors);
	io->access(reqs.data(), num_neighbors);
	return num_neighbors;
}

int worker_thread::process_pending_vertices(int max)
{
	if (max <= 0)
		return 0;

	int num_processed = 0;
	while (max > 0 && (!pending_vertices.is_empty() || !curr_pending.is_empty())) {
		if (!curr_pending.is_empty()) {
			int ret = process_pending_vertex(max);
			num_processed += ret;
			max -= ret;
			assert(max >= 0);
		}
		else {
			ext_mem_vertex ext_v = pending_vertices.pop_front();
			int num_neighbors = ext_v.get_num_edges(
					graph->get_required_neighbor_type());
			assert(num_neighbors > 0);
			curr_pending.init(pending_vertex::create(ext_v, graph));
			int ret = process_pending_vertex(max);
			assert(ret > 0);
			num_processed += ret;
			max -= ret;
			assert(max >= 0);
		}
	}
	return num_processed;
}

int worker_thread::process_activated_vertices(int max)
{
	if (max <= 0)
		return 0;

	vertex_id_t vertex_buf[max];
	stack_array<io_request> reqs(max);
	int num = graph->get_curr_activated_vertices(vertex_buf, max);
	for (int i = 0; i < num; i++) {
		compute_vertex &info = graph->get_vertex(vertex_buf[i]);
		reqs[i].init(new char[info.get_ext_mem_size()],
				info.get_ext_mem_off(),
				// TODO I might need to set the node id.
				info.get_ext_mem_size(), READ, io, -1);
	}
	if (graph->get_logger())
		graph->get_logger()->log(reqs.data(), num);
	io->access(reqs.data(), num);
	return num;
}

void worker_thread::run()
{
	while (true) {
		int num_visited = 0;
		int num = process_activated_vertices(MAX_IO_PEND_VERTICES);
		printf("thread %d process %d activated vertices\n", get_id(), num);
		num_visited += num;
		while (graph->get_num_curr_activated_vertices() > 0) {
			while (process_pending_vertices(
						MAX_IO_PEND_VERTICES - io->num_pending_ios()) > 0)
				io->wait4complete(io->num_pending_ios() / 10);
			num = process_activated_vertices(
					MAX_IO_PEND_VERTICES - io->num_pending_ios());
			num_visited += num;
			io->wait4complete(io->num_pending_ios() / 10 + max<int>(0,
						io->num_pending_ios() - MAX_IO_PEND_VERTICES));
		}
		// We have completed processing the activated vertices in this iteration.
		while (io->num_pending_ios() > 0) {
			io->wait4complete(1);
			process_pending_vertices(MAX_IO_PEND_VERTICES - io->num_pending_ios());
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
	assert(num_threads > 0 && num_nodes > 0);
	assert(num_threads % num_nodes == 0);
	for (int i = 0; i < num_threads; i++) {
		worker_thread *t = new worker_thread(this, factory,
				num_threads % num_nodes);
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
