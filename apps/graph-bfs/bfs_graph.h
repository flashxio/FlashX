#ifndef __BFS_GRAPH_H__
#define __BFS_GRAPH_H__

#include <vector>

#include "container.h"
#include "concurrency.h"
#include "io_interface.h"

#include "vertex_index.h"

class bfs_vertex: public in_mem_vertex_info
{
	enum {
		VISITED,
		ADD2Q,
	};

	atomic_flags<int> flags;
public:
	bfs_vertex() {
	}

	bfs_vertex(off_t off, int size): in_mem_vertex_info(off, size) {
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
};

class worker_thread;

class bfs_graph
{
	std::vector<bfs_vertex> vertices;
	atomic_integer num_visited_vertices;

	// These are two global queues. One contains the vertices that are being
	// processed in the current level. The other contains the vertices that
	// will be processed in the next level.
	// The queue for the current level.
	thread_safe_FIFO_queue<vertex_id_t> *queue;
	// The queue for the next level.
	thread_safe_FIFO_queue<vertex_id_t> *next_queue;
	atomic_integer level;
	volatile bool is_complete;

	// These are used for switching queues.
	pthread_mutex_t lock;
	pthread_barrier_t barrier1;
	pthread_barrier_t barrier2;

	worker_thread *first_thread;
	std::vector<worker_thread *> worker_threads;

	bfs_graph(int num_threads, int num_nodes, const std::string &graph_file,
			const std::string &index_file);
public:
	static bfs_graph *create(int num_threads, int num_nodes,
			const std::string &graph_file, const std::string &index_file) {
		return new bfs_graph(num_threads, num_nodes, graph_file, index_file);
	}

	bfs_vertex &get_vertex(vertex_id_t id) {
		return vertices[id];
	}

	void start(vertex_id_t id);

	/**
	 * The algorithm progresses to the next level.
	 * It returns true if no more work can progress.
	 */
	bool progress_next_level();

	/**
	 * Add vertices that may be processed in the next level.
	 */
	void add_next_vertices(vertex_id_t vertices[], int num);

	/**
	 * Get vertices to be processed in the current level.
	 */
	int get_curr_vertices(vertex_id_t vertices[], int num) {
		return queue->fetch(vertices, num);
	}

	int get_num_visited_vertices() const {
		return num_visited_vertices.get();
	}

	void wait4complete();
};

#endif
