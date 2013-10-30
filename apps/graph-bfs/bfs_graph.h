#ifndef __BFS_GRAPH_H__
#define __BFS_GRAPH_H__

#include <vector>

#include "container.h"
#include "concurrency.h"
#include "io_interface.h"

#include "vertex_index.h"
#include "graph_engine.h"

class bfs_vertex: public compute_vertex
{
	enum {
		VISITED,
	};

	atomic_flags<int> flags;
public:
	bfs_vertex(): compute_vertex(-1, -1, 0) {
	}

	bfs_vertex(vertex_id_t id, off_t off, int size): compute_vertex(
			id, off, size) {
	}

	bool has_visited() const {
		return flags.test_flag(VISITED);
	}

	bool set_visited(bool visited) {
		if (visited)
			return flags.set_flag(VISITED);
		else
			return flags.clear_flag(VISITED);
	}

	void run(graph_engine &graph, const ext_mem_vertex vertices[],
			int num);
};

class bfs_graph: public graph_engine
{
	bfs_graph(int num_threads, int num_nodes, const std::string &graph_file,
			graph_index *index, bool directed): graph_engine(num_threads, num_nodes,
				graph_file, index, directed) {
	}
public:
	static bfs_graph *create(int num_threads, int num_nodes,
			const std::string &graph_file, graph_index *index, bool directed) {
		return new bfs_graph(num_threads, num_nodes, graph_file, index, directed);
	}
};

#endif
