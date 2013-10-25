#include "thread.h"
#include "io_interface.h"

#include "bfs_graph.h"

const int VERTEX_BUF_SIZE = 1024;

/**
 * This callback is to process a vertex.
 */
class bfs_callback: public callback
{
	bfs_graph *graph;
public:
	bfs_callback(bfs_graph *graph) {
		this->graph = graph;
	}

	int invoke(io_request *reqs[], int num);
};

int bfs_callback::invoke(io_request *reqs[], int num)
{
	std::vector<vertex_id_t> neighbors;
	for (int i = 0; i < num; i++) {
		char *buf = reqs[i]->get_buf();
		ext_mem_undirected_vertex *v = ext_mem_undirected_vertex::deserialize(
				buf, reqs[i]->get_size());
		// We need to add the neighbors of the vertex to the queue of
		// the next level.
		for (int j = 0; j < v->get_num_edges(); j++) {
			vertex_id_t id = v->get_edge(j).get_to();
			bfs_vertex &info = graph->get_vertex(id);
			// If the vertex has been visited, we can skip it.
			if (info.has_visited())
				continue;
			neighbors.push_back(id);
		}
		delete [] buf;
	}
	graph->add_next_vertices(neighbors.data(), neighbors.size());
//	printf("add %ld vertices to the next level\n", neighbors.size());
	return 0;
}

class worker_thread: public thread
{
	file_io_factory *factory;
	io_interface *io;
	bfs_graph *graph;
public:
	worker_thread(bfs_graph *graph, file_io_factory *factory,
			int node_id): thread("worker_thread", node_id) {
		this->graph = graph;
		this->io = NULL;
		this->factory = factory;
	}

	void run();
	void init();
};

void worker_thread::init()
{
	io = factory->create_io(this);
	io->init();
	io->set_callback(new bfs_callback(graph));
}

void worker_thread::run()
{
	vertex_id_t vertex_buf[VERTEX_BUF_SIZE];
	io_request reqs[VERTEX_BUF_SIZE];

	while (true) {
		int num_visited = 0;
		int num = graph->get_curr_vertices(vertex_buf, VERTEX_BUF_SIZE);
		num_visited += num;
		while (num > 0) {
			for (int i = 0; i < num; i++) {
				bfs_vertex &info = graph->get_vertex(vertex_buf[i]);
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
			num = graph->get_curr_vertices(vertex_buf, VERTEX_BUF_SIZE);
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
}

bfs_graph::bfs_graph(int num_threads, int num_nodes,
		const std::string &graph_file, const std::string &index_file)
{
	vertex_index *indices = vertex_index::load(index_file);
	vertices.resize(indices->get_num_vertices());
	for (size_t i = 0; i < vertices.size(); i++) {
		off_t off = indices->get_vertex_off(i);
		int size = indices->get_vertex_size(i);
		vertices[i] = bfs_vertex(off, size);
	}
	vertex_index::destroy(indices);

	queue = new thread_safe_FIFO_queue<vertex_id_t>(-1, PAGE_SIZE, INT_MAX);
	next_queue = new thread_safe_FIFO_queue<vertex_id_t>(-1, PAGE_SIZE, INT_MAX);

	pthread_mutex_init(&lock, NULL);
	pthread_barrier_init(&barrier1, NULL, num_threads);
	pthread_barrier_init(&barrier2, NULL, num_threads);

	file_io_factory *factory = create_io_factory(graph_file,
			GLOBAL_CACHE_ACCESS);
	assert(num_threads % num_nodes == 0);
	for (int i = 0; i < num_threads; i++) {
		worker_thread *t = new worker_thread(this, factory,
				num_threads % num_nodes);
		t->start();
		worker_threads.push_back(t);
	}
	first_thread = worker_threads[0];
}

void bfs_graph::start(vertex_id_t id)
{
	printf("start on vertex %ld\n", id);
	queue->add(&id, 1);
	num_visited_vertices.inc(1);
	worker_threads[0]->activate();
}

void bfs_graph::add_next_vertices(vertex_id_t vertices[], int num)
{
	num_visited_vertices.inc(num);

	vertex_id_t to_add[num];
	int num_to_add = 0;

	for (int i = 0; i < num; i++) {
		bfs_vertex &v = get_vertex(vertices[i]);
		// When a vertex is added to the queue, we mark it as visited.
		// Therefore, a vertex can only be added to the queue once and
		// can only be visited once.
		if (!v.set_visited(true))
			to_add[num_to_add++] = vertices[i];
	}
	int num_added = next_queue->add(to_add, num_to_add);
	assert(num_added == num_to_add);
}

bool bfs_graph::progress_next_level()
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
		assert(queue->is_empty());
		thread_safe_FIFO_queue<vertex_id_t> *tmp = queue;
		queue = next_queue;
		next_queue = tmp;
		level.inc(1);
		printf("progress to level %d, there are %d vertices in this level\n",
				level.get(), queue->get_num_entries());
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
	return queue->is_empty();
}

void bfs_graph::wait4complete()
{
	for (unsigned i = 0; i < worker_threads.size(); i++) {
		worker_threads[i]->join();
	}
}

class graph_config
{
	int num_threads;
	std::string root_conf_file;
public:
	graph_config() {
		num_threads = 4;
	}

	void print_help();
	void print();

	void init(const config_map &map);

	int get_num_threads() const {
		return num_threads;
	}

	const std::string &get_root_conf_file() const {
		return root_conf_file;
	}
} graph_conf;

void graph_config::print_help()
{
}

void graph_config::print()
{
}

void graph_config::init(const config_map &map)
{
}

int main(int argc, char *argv[])
{
	if (argc < 5) {
		fprintf(stderr, "bfs conf_file graph_file index_file start_vertex\n");
		graph_conf.print_help();
		params.print_help();
		exit(-1);
	}

	std::string conf_file = argv[1];
	std::string graph_file = argv[2];
	std::string index_file = argv[3];
	vertex_id_t start_vertex = atoi(argv[4]);

	config_map configs(conf_file);
	configs.add_options(argv + 4, argc - 4);
	graph_conf.init(configs);
	graph_conf.print();

	init_io_system(configs);

	bfs_graph *graph = bfs_graph::create(graph_conf.get_num_threads(),
			params.get_num_nodes(), graph_file, index_file);
	graph->start(start_vertex);
	graph->wait4complete();
	printf("BFS from vertex 0 visits %d vertices\n", graph->get_num_visited_vertices());
}
