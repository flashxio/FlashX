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

	void run(graph_engine &graph);
};

class bfs_graph_index: public graph_index
{
	std::vector<bfs_vertex> vertices;
	
	bfs_graph_index(const std::string &index_file) {
		vertex_index *indices = vertex_index::load(index_file);
		vertices.resize(indices->get_num_vertices());
		for (size_t i = 0; i < vertices.size(); i++) {
			off_t off = indices->get_vertex_off(i);
			int size = indices->get_vertex_size(i);
			vertices[i] = bfs_vertex(i, off, size);
		}
		vertex_index::destroy(indices);
	}
public:
	static bfs_graph_index *create(const std::string &index_file) {
		return new bfs_graph_index(index_file);
	}

	virtual compute_vertex &get_vertex(vertex_id_t id) {
		return vertices[id];
	}

	virtual vertex_id_t get_max_vertex_id() const {
		return vertices.back().get_id();
	}

	virtual vertex_id_t get_min_vertex_id() const {
		return vertices.front().get_id();
	}
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
