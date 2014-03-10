#include "load_balancer.h"
#include "worker_thread.h"
#include "graph_engine.h"

const int MAX_STOLEN_VERTICES = 1024;

load_balancer::load_balancer(graph_engine &_graph,
		worker_thread &_owner): owner(_owner), graph(_graph)
{
	steal_thread_id = (owner.get_worker_id() + 1) % graph.get_num_threads();
	// TODO can I have a better way to do it?
	completed_stolen_vertices = (fifo_queue<vertex_id_t> *) malloc(
			graph.get_num_threads() * sizeof(fifo_queue<vertex_id_t>));
	for (int i = 0; i < graph.get_num_threads(); i++) {
		new (completed_stolen_vertices + i) fifo_queue<vertex_id_t>(
				_owner.get_node_id(), PAGE_SIZE, true);
	}
	num_completed_stolen_vertices = 0;
}

load_balancer::~load_balancer()
{
	for (int i = 0; i < graph.get_num_threads(); i++)
		completed_stolen_vertices[i].~fifo_queue<vertex_id_t>();
	free(completed_stolen_vertices);
}

/**
 * This steals vertices from other threads. It tries to steal more vertices
 * than it can process, and the remaining vertices will be placed in its
 * own activated vertex queue.
 */
int load_balancer::steal_activated_vertices(vertex_id_t vertex_buf[], int buf_size)
{
	if (steal_thread_id == owner.get_worker_id())
		steal_thread_id = (steal_thread_id + 1) % graph.get_num_threads();
	int num_tries = 0;
	vertex_id_t *steal_buf = new vertex_id_t[MAX_STOLEN_VERTICES];
	int num;
	do {
		worker_thread *t = graph.get_thread(steal_thread_id);
		num_tries++;

		num = t->steal_activated_vertices(steal_buf, MAX_STOLEN_VERTICES);
		// If we can't steal vertices from the thread, we should move
		// to the next thread.
		if (num == 0)
			steal_thread_id = (steal_thread_id + 1) % graph.get_num_threads();
		// If we have tried to steal vertices from all threads.
	} while (num == 0 && num_tries < graph.get_num_threads());

	int ret = min(buf_size, num);
	memcpy(vertex_buf, steal_buf, sizeof(vertex_buf[0]) * ret);
	// We stole more vertices than we can process this time.
	// The vertices stolen from another thread will also be placed in
	// the queue for currently activated vertices.
	if (num - ret > 0)
		owner.curr_activated_vertices.init(steal_buf + ret, num - ret, true);
	delete [] steal_buf;
	return ret;
}

void load_balancer::process_completed_stolen_vertices()
{
	if (num_completed_stolen_vertices == 0)
		return;

	int num_tot = 0;
	for (int i = 0; i < graph.get_num_threads(); i++) {
		fifo_queue<vertex_id_t> &q = completed_stolen_vertices[i];
		if (!q.is_empty()) {
			worker_thread *t = graph.get_thread(i);
			int num_completed = q.get_num_entries();
			num_tot += num_completed;
			stack_array<vertex_id_t> buf(num_completed);
			int ret = q.fetch(buf.data(), num_completed);
			assert(ret == num_completed);
			t->return_vertices(buf.data(), num_completed);
		}
	}
	assert(num_tot == num_completed_stolen_vertices);
	num_completed_stolen_vertices = 0;
}

void load_balancer::return_vertices(vertex_id_t ids[], int num)
{
	for (int i = 0; i < num; i++) {
		int part_id = graph.get_partitioner()->map(ids[i]);
		if (completed_stolen_vertices[part_id].is_full()) {
			completed_stolen_vertices[part_id].expand_queue(
					completed_stolen_vertices[part_id].get_size() * 2);
		}
		completed_stolen_vertices[part_id].push_back(ids[i]);
	}
	num_completed_stolen_vertices += num;
}

void load_balancer::reset()
{
	for (int i = 0; i < graph.get_num_threads(); i++)
		assert(completed_stolen_vertices[i].is_empty());
	assert(num_completed_stolen_vertices == 0);
}
