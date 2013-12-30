#ifndef __BFS_GRAPH_H__
#define __BFS_GRAPH_H__

/**
 * Copyright 2013 Da Zheng
 *
 * This file is part of SA-GraphLib.
 *
 * SA-GraphLib is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * SA-GraphLib is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with SA-GraphLib.  If not, see <http://www.gnu.org/licenses/>.
 */

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

	void run_on_neighbors(graph_engine &graph, const page_vertex *vertices[],
			int num) {
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
