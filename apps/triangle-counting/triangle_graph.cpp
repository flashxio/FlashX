#include <signal.h>
#include <google/profiler.h>

#include "thread.h"
#include "io_interface.h"

#include "graph_engine.h"
#include "graph_config.h"

atomic_number<long> num_triangles;

class triangle_vertex: public compute_vertex
{
public:
	triangle_vertex(): compute_vertex(-1, -1, 0) {
	}

	triangle_vertex(vertex_id_t id, off_t off, int size): compute_vertex(
			id, off, size) {
	}

	int count_triangle(const ext_mem_vertex &neighbor) const;

	void run(graph_engine &graph, const ext_mem_vertex vertices[],
			int num);
};

void triangle_vertex::run(graph_engine &graph, const ext_mem_vertex vertices[],
			int num)
{
	// If there aren't neighbors passed to the vertex's user code,
	// it's because the vertex doesn't have neighbors.
	if (num == 0) {
		assert(this->get_num_edges(graph.get_required_neighbor_type()) == 0);
		return;
	}

	int num_local_triangles = 0;
	for (int i = 0; i < num; i++) {
		const ext_mem_vertex *v = vertices + i;
		int k, l;
		for (k = 0, l = 0; k < this->get_num_edges(edge_type::IN_EDGE)
				&& l < v->get_num_edges(edge_type::OUT_EDGE); k++, l++) {
			vertex_id_t this_neighbor = this->get_neighbor(
					edge_type::IN_EDGE, k);
			vertex_id_t neigh_neighbor = v->get_neighbor(
					edge_type::OUT_EDGE, l);
			if (this_neighbor == neigh_neighbor) {
				num_local_triangles++;
				k++;
				l++;
			}
			else if (this_neighbor < neigh_neighbor)
				k++;
			else
				l++;
		}
	}
	if (num_local_triangles > 0)
		num_triangles.inc(num_local_triangles);
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
		fprintf(stderr, "triangle-counting conf_file graph_file index_file directed\n");
		graph_conf.print_help();
		params.print_help();
		exit(-1);
	}

	std::string conf_file = argv[1];
	std::string graph_file = argv[2];
	std::string index_file = argv[3];
	bool directed = atoi(argv[4]);

	config_map configs(conf_file);
	configs.add_options(argv + 5, argc - 5);
	graph_conf.init(configs);
	graph_conf.print();

	signal(SIGINT, int_handler);
	init_io_system(configs);

	graph_index *index = graph_index_impl<triangle_vertex>::create(index_file);
	graph_engine *graph = graph_engine::create(graph_conf.get_num_threads(),
			params.get_num_nodes(), graph_file, index, directed);
	graph->set_required_neighbor_type(edge_type::OUT_EDGE);
	printf("triangle counting starts\n");
	printf("prof_file: %s\n", graph_conf.get_prof_file().c_str());
	if (!graph_conf.get_prof_file().empty())
		ProfilerStart(graph_conf.get_prof_file().c_str());

	struct timeval start, end;
	gettimeofday(&start, NULL);
	graph->start_all();
	graph->wait4complete();
	gettimeofday(&end, NULL);

	if (!graph_conf.get_prof_file().empty())
		ProfilerStop();
	if (graph_conf.get_print_io_stat())
		print_io_thread_stat();
	printf("there are %ld triangles. It takes %f seconds\n",
			num_triangles.get(), time_diff(start, end));
}
