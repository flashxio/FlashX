#ifndef __VERTEX_PROGRAM_H__
#define __VERTEX_PROGRAM_H__

/**
 * Copyright 2014 Da Zheng
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

#include <memory>

#include "container.h"
#include "vertex.h"
#include "messaging.h"

class graph_engine;
class compute_vertex;
class page_vertex;
class vertex_message;
class worker_thread;

class vertex_program
{
	worker_thread *t;
	graph_engine *graph;
public:
	typedef std::shared_ptr<vertex_program> ptr;

	void init(graph_engine *graph, worker_thread *t) {
		this->t = t;
		this->graph = graph;
	}

	/**
	 * This is a pre-run before users get any information of adjacency list
	 * of vertices.
	 */
	virtual void run(compute_vertex &) = 0;

	/**
	 * Run user's code when the adjacency list of the vertex is read
	 * from disks.
	 */
	virtual void run(compute_vertex &, const page_vertex &vertex) = 0;

	/**
	 * Run user's code when the vertex receives messages from other.
	 */
	virtual void run_on_message(compute_vertex &, const vertex_message &) = 0;

	virtual void run_on_messages(const vertex_message *v_msgs[], int num) = 0;

	virtual void run_on_multicast_message(multicast_message &mmsg) = 0;

	const worker_thread &get_thread() const {
		return *t;
	}

	graph_engine &get_graph() {
		return *graph;
	}

	void multicast_msg(vertex_id_t ids[], int num, vertex_message &msg);
	void multicast_msg(edge_seq_iterator &it, vertex_message &msg);

	void send_msg(vertex_id_t dest, vertex_message &msg);

	/**
	 * Activate vertices that may be processed in the next level.
	 */
	void activate_vertices(vertex_id_t ids[], int num);
	void activate_vertices(edge_seq_iterator &it);

	void activate_vertex(vertex_id_t vertex) {
		activate_vertices(&vertex, 1);
	}
};

class vertex_program_creater
{
public:
	typedef std::unique_ptr<vertex_program_creater> ptr;

	virtual vertex_program::ptr create() const = 0;
};

size_t graph_get_vertices(graph_engine &graph, const worker_thread &,
		const local_vid_t ids[], int num_ids, compute_vertex *v_buf[]);

template<class vertex_type>
class vertex_program_impl: public vertex_program
{
	embedded_array<compute_vertex *, 1024> vertex_buf;
	embedded_array<local_vid_t, 1024> id_buf;
public:
	/**
	 * This is a pre-run before users get any information of adjacency list
	 * of vertices.
	 */
	virtual void run(compute_vertex &comp_v) {
		((vertex_type &) comp_v).run(*this);
	}

	/**
	 * Run user's code when the adjacency list of the vertex is read
	 * from disks.
	 */
	virtual void run(compute_vertex &comp_v, const page_vertex &vertex) {
		((vertex_type &) comp_v).run(*this, vertex);
	}

	/**
	 * Run user's code when the vertex receives messages from other.
	 */
	virtual void run_on_message(compute_vertex &comp_v,
			const vertex_message &msg) {
		((vertex_type &) comp_v).run_on_message(*this, msg);
	}

	virtual void run_on_messages(const vertex_message *v_msgs[],
			int num) {
		vertex_buf.resize(num);
		id_buf.resize(num);
		for (int i = 0; i < num; i++)
			id_buf[i] = v_msgs[i]->get_dest();
		graph_get_vertices(get_graph(), get_thread(), id_buf.data(), num,
				vertex_buf.data());
		for (int i = 0; i < num; i++) {
			assert(!v_msgs[i]->is_multicast());
			vertex_type *v = (vertex_type *) vertex_buf[i];
			v->run_on_message(*this, *v_msgs[i]);
		}
	}

	virtual void run_on_multicast_message(multicast_message &mmsg) {
		int num_dests = mmsg.get_num_dests();
		multicast_dest_list dest_list = mmsg.get_dest_list();

		vertex_buf.resize(num_dests);
		id_buf.resize(num_dests);
		for (int i = 0; i < num_dests; i++)
			id_buf[i] = dest_list.get_dest(i);
		graph_get_vertices(get_graph(), get_thread(), id_buf.data(), num_dests,
				vertex_buf.data());

		for (int i = 0; i < num_dests; i++) {
			vertex_type *v = (vertex_type *) vertex_buf[i];
			v->run_on_message(*this, mmsg);
		}
	}
};

#endif
