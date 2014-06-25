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
		fprintf(stderr, "usage: make_graph #vertices, #numedges\n");
		return EXIT_FAILURE;
	}
	std::size_t n = atol(argv[1]);
	std::size_t m = atol(argv[2]);

	fprintf(stderr, "Vertices = %ld\n", n);
	fprintf(stderr, "Edges = %ld\n", m);

	std::clock_t start;
	start = std::clock();
	minstd_rand gen;
	RMATGen gen_it(gen, n, m, 0.57, 0.19, 0.19, 0.05, true);
	RMATGen gen_end;
	for (; gen_it != gen_end; ++gen_it) {
		std::cout << gen_it->first << " " << gen_it->second << std::endl;
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

	std::cerr << "The time to build the graph = " << (std::clock()-start)/(double)CLOCKS_PER_SEC << std::endl;
	return 0;
}
