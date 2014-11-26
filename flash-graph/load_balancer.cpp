/**
 * Copyright 2014 Open Connectome Project (http://openconnecto.me)
 * Written by Da Zheng (zhengda1936@gmail.com)
 *
 * This file is part of FlashGraph.
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

#include "load_balancer.h"
#include "worker_thread.h"
#include "graph_engine.h"

namespace fg
{

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
int load_balancer::steal_activated_vertices(compute_vertex_pointer vertex_buf[],
		int buf_size)
{
	if (steal_thread_id == owner.get_worker_id())
		steal_thread_id = (steal_thread_id + 1) % graph.get_num_threads();
	int num_tries = 0;
	int num;
	do {
		worker_thread *t = graph.get_thread(steal_thread_id);
		num_tries++;

		num = t->steal_activated_vertices(vertex_buf, buf_size);
		// If we can't steal vertices from the thread, we should move
		// to the next thread.
		if (num == 0)
			steal_thread_id = (steal_thread_id + 1) % graph.get_num_threads();
		// If we have tried to steal vertices from all threads.
	} while (num == 0 && num_tries < graph.get_num_threads());

	// Record the owner thread of the stolen vertices.
	for (int i = 0; i < num; i++)
		stolen_vertex_map.insert(vertex_map_t::value_type(
					vertex_buf[i].get(), steal_thread_id));

	return num;
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
			BOOST_VERIFY(q.fetch(buf.data(), num_completed)
					== num_completed);
			t->return_vertices(buf.data(), num_completed);
		}
	}
	assert(num_tot == num_completed_stolen_vertices);
	num_completed_stolen_vertices = 0;
}

void load_balancer::return_vertices(const compute_vertex_pointer vs[], int num)
{
	for (int i = 0; i < num; i++) {
		compute_vertex_pointer v = vs[i];
		vertex_map_t::iterator it = stolen_vertex_map.find(v.get());
		assert(it != stolen_vertex_map.end());
		int part_id = it->second;
		// We don't need to return verticalled partitioned vertices to their
		// owner because messages are processed in the main vertices and the
		// main vertices cannot be stolen by other threads.
		if (!v.is_part()) {
			if (completed_stolen_vertices[part_id].is_full()) {
				completed_stolen_vertices[part_id].expand_queue(
						completed_stolen_vertices[part_id].get_size() * 2);
			}
			// TODO we can return compute_vertex_pointer and so we don't
			// map it back to local_vid_t.
			vertex_id_t id = graph.get_graph_index().get_vertex_id(part_id,
					*v.get());
			completed_stolen_vertices[part_id].push_back(id);
			num_completed_stolen_vertices++;
		}
		stolen_vertex_map.erase(it);
	}
}

void load_balancer::reset()
{
	for (int i = 0; i < graph.get_num_threads(); i++)
		assert(completed_stolen_vertices[i].is_empty());
	assert(num_completed_stolen_vertices == 0);
}

int load_balancer::get_stolen_vertex_part(const compute_vertex &v) const
{
	vertex_map_t::const_iterator it = stolen_vertex_map.find(&v);
	if (it != stolen_vertex_map.end())
		return it->second;
	else
		return -1;
}

}
