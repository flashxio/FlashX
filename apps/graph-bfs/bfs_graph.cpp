#include <signal.h>
#include <google/profiler.h>

#include "thread.h"
#include "io_interface.h"

#include "bfs_graph.h"

atomic_integer num_visited_vertices;

void bfs_vertex::run(graph_engine &graph, ext_mem_vertex &v,
		std::vector<vertex_id_t> &activated_vertices)
{
	vertex_id_t max_id = graph.get_max_vertex_id();
	vertex_id_t min_id = graph.get_min_vertex_id();

	// We need to add the neighbors of the vertex to the queue of
	// the next level.
	ext_mem_directed_vertex *directed_v = v.deserialize2directed();
	int num_activated = 0;
	for (int j = 0; j < directed_v->get_num_out_edges(); j++) {
		vertex_id_t id = directed_v->get_out_edge(j).get_to();
		assert(id >= min_id && id <= max_id);
		bfs_vertex &info = (bfs_vertex &) graph.get_vertex(id);
		// If the vertex has been visited, we can skip it.
		if (info.has_visited())
			continue;
		if (info.set_visited(true))
			continue;
		num_activated++;
		activated_vertices.push_back(id);
	}
	if (num_activated > 0)
		num_visited_vertices.inc(num_activated);
}

class graph_config
{
	int num_threads;
	std::string prof_file;
public:
	graph_config() {
		num_threads = 4;
	}

	void print_help();
	void print();

	void init(const config_map &map);

	const std::string &get_prof_file() const {
		return prof_file;
	}

	int get_num_threads() const {
		return num_threads;
	}
} graph_conf;

void graph_config::print_help()
{
	printf("Configuration parameters in graph algorithm.\n");
	printf("\tthreads: the number of threads processing the graph\n");
}

void graph_config::print()
{
	printf("Configuration parameters in graph algorithm.\n");
	printf("\tthreads: %d\n", num_threads);
}

void graph_config::init(const config_map &map)
{
	map.read_option_int("threads", num_threads);
	map.read_option("prof_file", prof_file);
}

void int_handler(int sig_num)
{
	if (!graph_conf.get_prof_file().empty())
		ProfilerStop();
	exit(0);
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

	signal(SIGINT, int_handler);
	init_io_system(configs);

	bfs_graph_index *index = bfs_graph_index::create(index_file);
	bfs_graph *graph = bfs_graph::create(graph_conf.get_num_threads(),
			params.get_num_nodes(), graph_file, index);
	printf("BFS starts\n");
	printf("prof_file: %s\n", graph_conf.get_prof_file().c_str());
	if (!graph_conf.get_prof_file().empty())
		ProfilerStart(graph_conf.get_prof_file().c_str());

	struct timeval start, end;
	gettimeofday(&start, NULL);
	graph->start(start_vertex);
	graph->wait4complete();
	gettimeofday(&end, NULL);

	if (!graph_conf.get_prof_file().empty())
		ProfilerStop();
	print_io_thread_stat();
	printf("BFS from vertex %ld visits %d vertices. It takes %f seconds\n",
			start_vertex, num_visited_vertices.get(),
			time_diff(start, end));
}
