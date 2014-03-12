#include <google/profiler.h>

#include <vector>

#include "thread.h"
#include "io_interface.h"
#include "container.h"
#include "concurrency.h"

#include "vertex_index.h"
#include "graph_engine.h"
#include "graph_config.h"
#include "worker_thread.h"

std::vector<int> thread_delays;

class count_message: public vertex_message
{
	int count;
public:
	count_message(): vertex_message(sizeof(count_message), false) {
		count = 1;
	}

	int get_count() const {
		return count;
	}
};

class test_vertex: public compute_vertex
{
	size_t num_edges;
	size_t add1_count;
public:
	test_vertex() {
		num_edges = 0;
		add1_count = 0;
	}

	test_vertex(vertex_id_t id, const vertex_index *index): compute_vertex(
			id, index) {
		num_edges = 0;
		add1_count = 0;
	}

	void init() {
		add1_count = 0;
	}

	void verify_result() const {
		assert(add1_count == num_edges * 2);
	}

	void run(graph_engine &graph) {
		vertex_id_t id = get_id();
		graph.request_vertices(*this, &id, 1);
	}

	void run(graph_engine &graph, const page_vertex &vertex);

	virtual void run_on_messages(graph_engine &,
			const vertex_message *msgs[], int num) {
		for (int i = 0; i < num; i++) {
			const count_message *msg = (const count_message *) msgs[i];
			add1_count += msg->get_count();
		}
	}
};

void test_vertex::run(graph_engine &graph, const page_vertex &vertex)
{
	worker_thread *t = (worker_thread *) thread::get_curr_thread();
	int worker_id = t->get_worker_id();
	if ((size_t) worker_id < thread_delays.size() / 2) {
		if (thread_delays[worker_id] == 0) {
			int slp_time = random() % 10;
			printf("worker %d sleeps %d seconds\n", worker_id, slp_time);
			sleep(slp_time);
			thread_delays[worker_id] = 1;
		}
	}
	// We need to add the neighbors of the vertex to the queue of
	// the next level.
	page_byte_array::const_iterator<vertex_id_t> end_it
		= vertex.get_neigh_end(BOTH_EDGES);
	num_edges = vertex.get_num_edges(BOTH_EDGES);
	stack_array<vertex_id_t, 1024> dest_buf(num_edges);
	int num_dests = 0;
	for (page_byte_array::const_iterator<vertex_id_t> it
			= vertex.get_neigh_begin(BOTH_EDGES); it != end_it; ++it) {
		vertex_id_t id = *it;
		dest_buf[num_dests++] = id;
		add1_count++;
	}
	count_message msg;
	graph.multicast_msg(dest_buf.data(), num_dests, msg);
}

int main(int argc, char *argv[])
{
	if (argc < 4) {
		fprintf(stderr, "test conf_file graph_file index_file");
		exit(-1);
	}

	std::string conf_file = argv[1];
	std::string graph_file = argv[2];
	std::string index_file = argv[3];

	config_map configs(conf_file);
	graph_conf.init(configs);
	graph_conf.print();

	init_io_system(configs);

	graph_index *index = NUMA_graph_index<test_vertex>::create(index_file,
			graph_conf.get_num_threads(), params.get_num_nodes());

	graph_engine *graph = graph_engine::create(graph_conf.get_num_threads(),
			params.get_num_nodes(), graph_file, index);
	printf("prof_file: %s\n", graph_conf.get_prof_file().c_str());
	if (!graph_conf.get_prof_file().empty())
		ProfilerStart(graph_conf.get_prof_file().c_str());
	thread_delays.resize(graph->get_num_threads());

	for (int i = 0; i < 1000; i++) {
		printf("test %d\n", i);
		struct timeval start, end;
		gettimeofday(&start, NULL);
		graph->start_all();
		graph->wait4complete();
		gettimeofday(&end, NULL);

		NUMA_graph_index<test_vertex>::const_iterator it
			= ((NUMA_graph_index<test_vertex> *) index)->begin();
		NUMA_graph_index<test_vertex>::const_iterator end_it
			= ((NUMA_graph_index<test_vertex> *) index)->end();
		for (; it != end_it; ++it) {
			const test_vertex &v = (const test_vertex &) *it;
			v.verify_result();
		}
		printf("It takes %f seconds\n", time_diff(start, end));
	}

	if (!graph_conf.get_prof_file().empty())
		ProfilerStop();
	if (graph_conf.get_print_io_stat())
		print_io_thread_stat();
	graph_engine::destroy(graph);
	destroy_io_system();
}
