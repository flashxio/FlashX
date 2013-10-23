#include "graph.h"

int main(int argc, char *argv[])
{
	if (argc < 5) {
		fprintf(stderr,
				"el2al edge_list_file adjacency_list_file index_file directed\n");
		exit(-1);
	}

	const std::string edge_list_file = argv[1];
	const std::string adjacency_list_file = argv[2];
	const std::string index_file = argv[3];
	bool directed = atoi(argv[4]) != 0;

	if (directed) {
		directed_graph *g = directed_graph::load_edge_list_text(
				edge_list_file);
		vertex_index *index = g->create_vertex_index();
		g->dump(adjacency_list_file);
		index->dump(index_file);
		printf("There are %ld vertices, %ld non-empty vertices and %ld edges\n",
				g->get_num_vertices(), g->get_num_non_empty_vertices(),
				g->get_num_in_edges());
		directed_graph::destroy(g);
	}
	else {
		undirected_graph *g = undirected_graph::load_edge_list_text(
				edge_list_file);
		vertex_index *index = g->create_vertex_index();
		g->dump(adjacency_list_file);
		index->dump(index_file);
		printf("There are %ld vertices, %ld non-empty vertices and %ld edges\n",
				g->get_num_vertices(), g->get_num_non_empty_vertices(),
				g->get_num_edges());
		undirected_graph::destroy(g);
	}
}
