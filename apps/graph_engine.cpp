#include <algorithm>
#include <parallel/algorithm>

#include "io_interface.h"

#include "graph_engine.h"

const int VERTEX_BUF_SIZE = 1024;

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
	int num_completed_neighbors;
	edge_type required_neighbor_type;

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
};

int vertex_callback::invoke(io_request *reqs[], int num)
{
	for (int i = 0; i < num; i++) {
		char *req_buf = reqs[i]->get_buf();
		size_t req_size = reqs[i]->get_size();
		assert(this->io == reqs[i]->get_io());
		ext_mem_vertex ext_v(req_buf, req_size, graph->is_directed());
		if (graph->get_required_neighbor_type() == edge_type::NONE) {
			compute_vertex &v = graph->get_vertex(ext_v.get_id());
			v.materialize(ext_v);
			v.run(*graph, NULL, 0);
			v.dematerialize();
			delete [] req_buf;
		}
		// We just fetched a vertex, we need to fetch its neighbors to
		// perform computation.
		else if (reqs[i]->get_user_data() == NULL) {
			int num_neighbors = ext_v.get_num_edges(
					graph->get_required_neighbor_type());
			io_request reqs[num_neighbors];
			pending_vertex *pending = pending_vertex::create(req_buf, req_size,
					graph);
			for (int j = 0; j < num_neighbors; j++) {
				vertex_id_t neighbor = ext_v.get_neighbor(
						graph->get_required_neighbor_type(), j);
				compute_vertex &info = graph->get_vertex(neighbor);
				reqs[j].init(new char[info.get_ext_mem_size()],
						info.get_ext_mem_off(),
						// TODO I might need to set the node id.
						info.get_ext_mem_size(), READ, io, -1);
				reqs[j].set_user_data(pending);
			}
			io->access(reqs, num);
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

class worker_thread: public thread
{
	file_io_factory *factory;
	io_interface *io;
	graph_engine *graph;

	std::vector<vertex_id_t> activated_vertices;
public:
	worker_thread(graph_engine *graph, file_io_factory *factory,
			int node_id): thread("worker_thread", node_id) {
		this->graph = graph;
		this->io = NULL;
		this->factory = factory;
	}

	void run();
	void init();

	std::vector<vertex_id_t> &get_activated_vertices() {
		return activated_vertices;
	}
};

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

void worker_thread::run()
{
	vertex_id_t vertex_buf[VERTEX_BUF_SIZE];
	io_request reqs[VERTEX_BUF_SIZE];

	while (true) {
		int num_visited = 0;
		int num = graph->get_curr_activated_vertices(vertex_buf, VERTEX_BUF_SIZE);
		num_visited += num;
		while (num > 0) {
			for (int i = 0; i < num; i++) {
				compute_vertex &info = graph->get_vertex(vertex_buf[i]);
				reqs[i].init(new char[info.get_ext_mem_size()],
						info.get_ext_mem_off(),
						// TODO I might need to set the node id.
						info.get_ext_mem_size(), READ, io, -1);
			}
			io->access(reqs, num);
			int num2wait = num / 10;
			if (num2wait == 0)
				num2wait = 1;
			io->wait4complete(num2wait);
			num = graph->get_curr_activated_vertices(vertex_buf, VERTEX_BUF_SIZE);
			num_visited += num;
		}
		io->wait4complete(io->num_pending_ios());
		printf("thread %d visited %d vertices\n", this->get_id(), num_visited);

		// Now we have finished this level, we can progress to the next level.
		bool completed = graph->progress_next_level();
		printf("thread %d finish in a level, completed? %d\n", get_id(), completed);
		if (completed)
			break;
	}
	stop();
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
		t->start();
		worker_threads.push_back(t);
	}
	first_thread = worker_threads[0];

	activated_vertices = new sorted_vertex_queue();
	activated_vertex_buf = new vertex_collection(worker_threads);
}

void graph_engine::start(vertex_id_t ids[], int num)
{
	std::vector<vertex_id_t> starts;
	starts.assign(ids, ids + num);
	activated_vertices->init(starts, false);
	worker_threads[0]->activate();
}

void graph_engine::start_all()
{
	std::vector<vertex_id_t> all_vertices;
	vertices->get_all_vertices(all_vertices);
	activated_vertices->init(all_vertices, true);
	for (unsigned i = 0; i < worker_threads.size(); i++)
		worker_threads[i]->activate();
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
