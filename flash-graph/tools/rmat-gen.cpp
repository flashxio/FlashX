#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/rmat_graph_generator.hpp>
#include <boost/random/linear_congruential.hpp>
#include <boost/graph/graph_traits.hpp>
#include <iostream>
#include <ctime>
#include <stdio.h>

using namespace boost;
typedef adjacency_list<> Graph;
typedef rmat_iterator<minstd_rand, Graph> RMATGen;
typedef graph_traits<Graph>::vertex_iterator vertex_iter;
typedef property_map<Graph, vertex_index_t>::type IndexMap;

int main(int argc, char* argv[])
{
	if (argc < 3) {
		fprintf(stderr, "usage: make_graph #vertices, #numedges [output]\n");
		return EXIT_FAILURE;
	}
	typedef boost::graph_traits<Graph>::vertices_size_type vertices_size_type;
	typedef boost::graph_traits<Graph>::edges_size_type edges_size_type;
	vertices_size_type n = atol(argv[1]);
	edges_size_type m = atol(argv[2]);

	FILE *f = stdout;
	std::string out_file;
	if (argc >= 4) {
		out_file = argv[3];
		f = fopen(argv[3], "w");
		assert(f);
	}

	fprintf(stderr, "Vertices = %ld\n", n);
	fprintf(stderr, "Edges = %ld\n", m);

	std::clock_t start;
	start = std::clock();
	minstd_rand gen;
	RMATGen gen_it(gen, n, m, 0.57, 0.19, 0.19, 0.05, true);
	RMATGen gen_end;
	for (; gen_it != gen_end; ++gen_it) {
		fprintf(f, "%ld %ld\n", gen_it->first, gen_it->second);
	}
#if 0
	// Create graph with 100 nodes and 400 edges
	Graph g(RMATGen(gen, n, m, 0.57, 0.19, 0.19, 0.05, true), RMATGen(), n);

	IndexMap index = get(vertex_index, g);

	// Get vertex set
#if 0
	std::pair<vertex_iter, vertex_iter> vp;
	for (vp = vertices(g); vp.first != vp.second; ++vp.first)
		std::cout << index[*vp.first] <<  " ";
	std::cout << std::endl;
#endif

	// Get edge set
	graph_traits<Graph>::edge_iterator ei, ei_end;
	for (tie(ei, ei_end) = edges(g); ei != ei_end; ++ei)
		std::cout << index[source(*ei, g)]<< " " << index[target(*ei, g)] << "\n";
#endif
	if (f != stdout) {
		printf("close %s\n", out_file.c_str());
		fclose(f);
	}

	std::cerr << "The time to build the graph = " << (std::clock()-start)/(double)CLOCKS_PER_SEC << std::endl;
	return 0;
}
