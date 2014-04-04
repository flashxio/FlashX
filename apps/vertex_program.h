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

class graph_engine;
class compute_vertex;
class page_vertex;
class vertex_message;

class vertex_program
{
public:
	typedef std::unique_ptr<vertex_program> ptr;
	/**
	 * This is a pre-run before users get any information of adjacency list
	 * of vertices.
	 */
	virtual void run(graph_engine &graph, compute_vertex &) = 0;

	/**
	 * Run user's code when the adjacency list of the vertex is read
	 * from disks.
	 */
	virtual void run(graph_engine &graph, compute_vertex &,
			const page_vertex &vertex) = 0;

	/**
	 * Run user's code when the vertex receives messages from other.
	 */
	virtual void run_on_message(graph_engine &, compute_vertex &,
			const vertex_message &msg) = 0;

	virtual void run_on_messages(graph_engine &graph, const vertex_message *v_msgs[],
			int num) = 0;

	virtual void run_on_multicast_message(graph_engine &graph,
			multicast_message &mmsg) = 0;

	virtual vertex_program::ptr clone() const = 0;
};

compute_vertex &graph_get_vertex(graph_engine &, vertex_id_t id);

template<class vertex_type>
class vertex_program_impl: public vertex_program
{
public:
	/**
	 * This is a pre-run before users get any information of adjacency list
	 * of vertices.
	 */
	virtual void run(graph_engine &graph, compute_vertex &comp_v) {
		((vertex_type &) comp_v).run(graph);
	}

	/**
	 * Run user's code when the adjacency list of the vertex is read
	 * from disks.
	 */
	virtual void run(graph_engine &graph, compute_vertex &comp_v,
			const page_vertex &vertex) {
		((vertex_type &) comp_v).run(graph, vertex);
	}

	/**
	 * Run user's code when the vertex receives messages from other.
	 */
	virtual void run_on_message(graph_engine &graph, compute_vertex &comp_v,
			const vertex_message &msg) {
		const vertex_message *msgs[1] = {&msg};
		((vertex_type &) comp_v).run_on_messages(graph, msgs, 1);
	}

	virtual vertex_program::ptr clone() const {
		return vertex_program::ptr(new vertex_program_impl<vertex_type>());
	}

	virtual void run_on_messages(graph_engine &graph, const vertex_message *v_msgs[],
			int num) {
		for (int i = 0; i < num; i++) {
			assert(!v_msgs[i]->is_multicast());
			vertex_id_t id = v_msgs[i]->get_dest();
			vertex_type &v = (vertex_type &) graph_get_vertex(graph, id);
			v.run_on_messages(graph, &v_msgs[i], 1);
		}
	}

	virtual void run_on_multicast_message(graph_engine &graph,
			multicast_message &mmsg) {
		int num_dests = mmsg.get_num_dests();
		multicast_dest_list dest_list = mmsg.get_dest_list();

		const vertex_message *msgs[1] = {&mmsg};
		for (int i = 0; i < num_dests; i++) {
			vertex_id_t id = dest_list.get_dest(i);
			// TODO now the size is the entire message. Now the message
			// is considered as non-empty.
			if (!mmsg.is_empty()) {
				vertex_type &v = (vertex_type &) graph_get_vertex(graph, id);
				v.run_on_messages(graph, msgs, 1);
			}
		}
	}
};

#endif
