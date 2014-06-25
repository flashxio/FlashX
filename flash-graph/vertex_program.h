#ifndef __VERTEX_PROGRAM_H__
#define __VERTEX_PROGRAM_H__

/*
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

#include <memory>

#include "container.h"
#include "vertex.h"
#include "messaging.h"

class graph_engine;
class compute_vertex;
class page_vertex;
class vertex_message;
class worker_thread;

/**
 *  This class allows users to customize the default `vertex_program`.
 *  For instance extending this class can allow a user to easily create
 *  & manage per-thread data structures as opposed to per-vertex which is
 *  possible by extending any flavor of `compute_vertex`.
 */
class vertex_program
{
	worker_thread *t;
	graph_engine *graph;

	std::unique_ptr<std::vector<local_vid_t>[]> vid_bufs;
	std::unique_ptr<vertex_loc_t[]> vertex_locs;
	size_t vloc_size;

	// The message senders to send messages to all other threads.
	// There are n senders, n is the total number of threads used by
	// the graph engine.
	std::vector<simple_msg_sender *> msg_senders;
	std::vector<multicast_msg_sender *> multicast_senders;
	std::vector<multicast_msg_sender *> activate_senders;
    
	multicast_msg_sender &get_activate_sender(int thread_id) const {
		return *activate_senders[thread_id];
	}

	multicast_msg_sender &get_multicast_sender(int thread_id) const {
		return *multicast_senders[thread_id];
	}

	simple_msg_sender &get_msg_sender(int thread_id) const {
		return *msg_senders[thread_id];
	}
public:
	typedef std::shared_ptr<vertex_program> ptr; /**Smart pointer by which `vertex_program`s should be accessed.*/
    
    /** \brief Destructor */
	virtual ~vertex_program();
    
    /**
	 * \internal
	 * \brief Initiate function a.k.a constructor.
     *  \param graph The graph engine pointer.
     *  \param t The array of worker threads running the program.
     */
	void init(graph_engine *graph, worker_thread *t) {
		this->t = t;
		this->graph = graph;
	}
    
    /* Internal */
	void init_messaging(const std::vector<worker_thread *> &threads,
			std::shared_ptr<slab_allocator> msg_alloc);

	/**
	 * \brief This is a pre-run before users get any information of adjacency list
	 * of vertices. This is commonly where a user would issue a request for the vertex
     * data it needs on disk for the two other run methods.
     *
     *  \param vertex A vertex.
	 */
	virtual void run(compute_vertex &vertex) = 0;

	/**
	 * \brief Run user's code ideally/generally when the adjacency list
     *      of a vertex is read from disks.
     * \param comp_v A `compute_vertex`.
     * \param vertex A `page vertex`.
	 */
	virtual void run(compute_vertex &comp_v, const page_vertex &vertex) = 0;

	/**
	 * \brief Run user's code when the vertex receives messages from others.
     *  \param c_vertex A `compute_vertex`.
     *  \param vertex_m A *single* `vertex_message` received from a vertex in the current iteration.
	 */
	virtual void run_on_message(compute_vertex &c_vertex, const vertex_message &vertex_m) = 0;
    
    /**
	 * \brief Run user's code when the vertex receives messages from others.
     *  \param v_msgs `vertex_message`s received from one or more vertices in the current iteration.
     *  \param num The number of messages received.
	 */
	virtual void run_on_messages(const vertex_message *v_msgs[], int num) = 0;
    
    /**
     * \brief Run user's code when a multicast message is received.
     *  \param mmsg The message(s) received by the vertex.
     */
	virtual void run_on_multicast_message(multicast_message &mmsg) = 0;
    
    /**
     * \brief Perform some user defined action on a vertex when the current iteration comes to an end.
     *  \param cv A `compute_vertex`.
     */
	virtual void notify_iteration_end(compute_vertex &cv) = 0;
    
    /* Internal */
	const worker_thread &get_thread() const {
		return *t;
	}
    
    /**
     * \brief Get a pointer to the `graph_engine`.
     *  \return A pointer to the `graph_engine`.
     */
	graph_engine &get_graph() {
		return *graph;
	}
    
    /** 
     * \brief Multicast the same message to several other vertices. If the number of vertices
     *      receiving the message is too small the graph engine will automatically alter the
     *      message type to point-to-point.
     *  \param ids The vertex IDs a user wants to send the message to.
     *  \param num The number of vertices a user wants to send the message to.
     *  \param msg The message intended for recepients.
     */
	void multicast_msg(vertex_id_t ids[], int num, vertex_message &msg);
    
    /**
     * \brief Multicast the same message to several other vertices. If the number of vertices
     *      receiving the message is too small the graph engine will automatically alter the
     *      message type to point-to-point.
     *  \param it An `edge_seq_iterator` defining which vertices to send the message to.
     *  \param msg The message intended for recepients.
     */
	void multicast_msg(edge_seq_iterator &it, vertex_message &msg);
    
    /**
     *  \brief Send a point-to-point message from one vertex to the next.
     *  \param dest The destination ID of the vertex a user is sending to.
     *  \param msg The message intended for the recepient.
     */
	void send_msg(vertex_id_t dest, vertex_message &msg);

	/**
	 * \brief Activate vertices to be processed in the next level (iteration).
     *  \param ids The unique IDs of the vertices to be activated.
     *  \param num The number of all the vertices to be activated.
	 */
	void activate_vertices(vertex_id_t ids[], int num);
    
	/**
	 * \brief Activate vertices to be processed in the next level (iteration).
     *  \param it An `edge_seq_iterator` defining which vertices to activate in the next iteration.
	 */
	void activate_vertices(edge_seq_iterator &it);
    
    /**
	 * \brief Activate a singel vertex to be processed in the next level (iteration).
     *  \param ids The unique ID of the vertex to be activated.
	 */
	void activate_vertex(vertex_id_t vertex) {
		activate_vertices(&vertex, 1);
	}

    /* Internal */
	void flush_msgs();

	/**
	 * \brief A vertex requests the end of an iteration.
	 * `notify_iteration_end' of the vertex will be invoked at the end
	 * of an iteration.
     *  \param v The vertex that will receive the notification.
	 */
	void request_notify_iter_end(const compute_vertex &v);
};

/**
 * \brief Extend/Override when defining a custom vertex program.
 *        The graph engine uses this to construct vertex programs for
 *        each worker thread.
 */
class vertex_program_creater
{
public:
	typedef std::unique_ptr<vertex_program_creater> ptr; /** Pointer defining object access. */
    
    /**
     *  \brief Much like a constructor -- implement this in lieu of that.
     */
	virtual vertex_program::ptr create() const = 0;
};

size_t graph_get_vertices(graph_engine &graph, const worker_thread &,
		const local_vid_t ids[], int num_ids, compute_vertex *v_buf[]);

/**
 * \brief The default implementation of a vertex program in the graph engine.
 */
template<class vertex_type>
class vertex_program_impl: public vertex_program
{
	embedded_array<compute_vertex *, 1024> vertex_buf;
	embedded_array<local_vid_t, 1024> id_buf;
public:
	
	virtual void run(compute_vertex &comp_v) {
		((vertex_type &) comp_v).run(*this);
	}

	virtual void run(compute_vertex &comp_v, const page_vertex &vertex) {
		((vertex_type &) comp_v).run(*this, vertex);
	}

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
    
	virtual void notify_iteration_end(compute_vertex &comp_v) {
		((vertex_type &) comp_v).notify_iteration_end(*this);
	}
};

#endif
