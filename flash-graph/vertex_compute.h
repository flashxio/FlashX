#ifndef __VERTEX_COMPUTE_H__
#define __VERTEX_COMPUTE_H__

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

#include <algorithm>
#include <deque>
#include <unordered_map>

#include "io_interface.h"
#include "slab_allocator.h"

#include "vertex_request.h"
#include "scan_pointer.h"
#include "graph_index.h"

class worker_thread;
class graph_engine;
class compute_vertex;
class compute_directed_vertex;

/**
 * This data structure represents an active vertex that is being processed
 * in a worker thread. It is used to handle two types of asynchronous
 * requests: for the adjacency list and the number of edges.
 */
class vertex_compute: public user_compute
{
	struct vertex_info_comp
	{
		bool operator()(const ext_mem_vertex_info &info1,
				const ext_mem_vertex_info &info2) {
			return info1.get_off() > info2.get_off();
		}
	};

	// TODO use the embedded array as the container.
	std::priority_queue<ext_mem_vertex_info, std::vector<ext_mem_vertex_info>,
		vertex_info_comp> requested_vertices;
protected:
	graph_engine *graph;

	// The thread that creates the vertex compute.
	worker_thread *issue_thread;
	compute_vertex_pointer v;

	/*
	 * These two variables keep track of the number of completed requests
	 * for adjacency lists.
	 */

	// The number of requested vertices that will be read in the user compute.
	size_t num_requested;
	// The number of issued requests.
	size_t num_issued;
	// The number of vertices read by the user compute.
	size_t num_complete_fetched;

	/*
	 * These two variables keep track of the number of completed requests
	 * for the number of edges.
	 */

	size_t num_edge_requests;
	size_t num_edge_completed;

	size_t get_num_pending_ios() const {
		assert(num_issued >= num_complete_fetched);
		return num_issued - num_complete_fetched;
	}

	bool issued_to_io() const {
		// When the vertex_compute is created, it has one reference.
		// If the vertex_compute has been issued to SAFS, its reference count
		// should be larger than 1.
		return get_ref() > 1 || get_num_pending_ios() > 0;
	}

	void start_run();
	void finish_run();
public:
	vertex_compute(graph_engine *graph,
			compute_allocator *alloc): user_compute(alloc) {
		this->graph = graph;
		issue_thread = (worker_thread *) thread::get_curr_thread();
		num_requested = 0;
		num_complete_fetched = 0;
		num_issued = 0;
		num_edge_requests = 0;
		num_edge_completed = 0;
	}

	void init(compute_vertex_pointer v) {
		this->v = v;
	}

	/*
	 * The methods below are from user_compute.
	 */

	virtual int serialize(char *buf, int size) const {
		return 0;
	}

	virtual int get_serialized_size() const {
		return 0;
	}

	virtual void set_scan_dir(bool forward) {
	}

	virtual int has_requests() {
		return requested_vertices.size() > 0;
	}

	virtual request_range get_next_request();

	virtual void run(page_byte_array &);

	virtual bool has_completed() {
		// If the user compute has got all requested data and it has
		// no more requests to issue, we can consider the user compute
		// has been completed.
		// NOTE: it's possible that requested data may not be passed to
		// this user compute, so we only count the requests that are going
		// to be passed to this user compute.
		return num_requested == num_complete_fetched && !has_requests();
	}

	/*
	 * The methods below deal with requesting the adjacency list of vertices.
	 */

	/*
	 * The method accepts the requests from graph applications and issues
	 * the request for the adjacency lists of vertices. It has to get
	 * the location and the size of the vertices from the vertex index
	 * before issuing the real requests to SAFS.
	 */
	virtual void request_vertices(vertex_id_t ids[], size_t num);

	/*
	 * This is a callback function. When the location and the size of
	 * a vertex is ready, the vertex index notifies the vertex compute
	 * of the information.
	 */
	void issue_io_request(const ext_mem_vertex_info &info);

	/*
	 * The methods below deal with requesting # edges of vertices.
	 */

	/*
	 * This is a callback function. When the vertex index gets the vertex size,
	 * it notifies the vertex_compute of this information.
	 */
	void run_on_vertex_size(vertex_id_t id, vsize_t size);

	/*
	 * This methods accepts the requests from graph applications and issues
	 * the requests to the vertex index.
	 */
	virtual void request_num_edges(vertex_id_t ids[], size_t num);

	/*
	 * The methods are used for both cases (requesting the adjacency list
	 * and #edges).
	 */

	/*
	 * This indicates the total number of pending requests for both adjacency
	 * lists and #edges.
	 */
	virtual size_t get_num_pending() const {
		return (num_edge_requests - num_edge_completed)
			+ (num_requested - num_complete_fetched);
	}

	vertex_id_t get_id() const;

	graph_engine &get_graph() {
		return *graph;
	}
};

class directed_vertex_compute: public vertex_compute
{
	typedef std::unordered_map<vertex_id_t, page_byte_array *> combine_map_t;
	combine_map_t combine_map;

	void run_on_page_vertex(page_directed_vertex &);
public:
	directed_vertex_compute(graph_engine *graph,
			compute_allocator *alloc): vertex_compute(graph, alloc) {
	}

	virtual void set_scan_dir(bool forward) {
		vertex_compute::set_scan_dir(forward);
	}

	virtual void run(page_byte_array &);

	/*
	 * These two methods accept the requests from graph applications and issue
	 * the request for the adjacency lists of vertices. It has to get
	 * the location and the size of the vertices from the vertex index
	 * before issuing the real requests to SAFS.
	 * `request_vertices' requests both edge lists of vertices.
	 * `request_partial_vertices' requests one type of edge lists of vertices.
	 */
	virtual void request_vertices(vertex_id_t ids[], size_t num);
	void request_partial_vertices(directed_vertex_request reqs[], size_t num);

	/*
	 * This is a callback function. When the vertex index gets the vertex size,
	 * it notifies the vertex_compute of this information.
	 */
	void run_on_vertex_size(vertex_id_t id, size_t in_size, size_t out_size);

	using vertex_compute::issue_io_request;
	void issue_io_request(const ext_mem_vertex_info &in_info,
			const ext_mem_vertex_info &out_info);

	/*
	 * This methods accepts the requests from graph applications and issues
	 * the requests to the vertex index.
	 */
	void request_num_edges(vertex_id_t ids[], size_t num);
};

/*
 * This class is to compute on the undirected vertices requested
 * by a single I/O request.
 */
class merged_vertex_compute: public user_compute
{
	vertex_id_t start_id;
	int num_vertices;
	bool complete;
	graph_engine *graph;
protected:
	worker_thread *issue_thread;

	void start_run(compute_vertex_pointer v);
	void finish_run(compute_vertex_pointer v);
public:
	merged_vertex_compute(graph_engine *graph,
			compute_allocator *alloc): user_compute(alloc) {
		this->graph = graph;
		start_id = INVALID_VERTEX_ID;
		num_vertices = 0;
		complete = false;
		issue_thread = (worker_thread *) thread::get_curr_thread();
	}

	graph_engine &get_graph() {
		return *graph;
	}

	vertex_id_t get_start_id() const {
		return start_id;
	}

	int get_num_vertices() const {
		return num_vertices;
	}

	virtual void init(vertex_id_t start_id, int num_vertices, edge_type type) {
		this->start_id = start_id;
		this->num_vertices = num_vertices;
	}

	virtual int serialize(char *buf, int size) const {
		return 0;
	}

	virtual int get_serialized_size() const {
		return 0;
	}

	virtual void set_scan_dir(bool forward) {
	}

	virtual int has_requests() {
		return false;
	}

	virtual request_range get_next_request() {
		assert(0);
	}

	virtual void run(page_byte_array &arr) {
		// TODO
		assert(0);
		complete = true;
	}

	virtual bool has_completed() {
		return complete;
	}
};

/*
 * This class is to compute on the directed vertices requested
 * by a single I/O request.
 */
class merged_directed_vertex_compute: public merged_vertex_compute
{
	edge_type type;
	int num_fetched_arrs;
	int num_required_arrs;
	page_byte_array *buffered_arr;

	void run_on_array(page_byte_array &arr);
	void run_on_arrays(page_byte_array &in_arr, page_byte_array &out_arr);
public:
	merged_directed_vertex_compute(graph_engine *graph,
			compute_allocator *alloc): merged_vertex_compute(graph, alloc) {
		type = edge_type::NONE;
		num_fetched_arrs = 0;
		num_required_arrs = 0;
		buffered_arr = NULL;
	}

	void init(vertex_id_t start_id, int num_vertices, edge_type type) {
		merged_vertex_compute::init(start_id, num_vertices, type);
		this->type = type;
		this->num_fetched_arrs = 0;
		switch(type) {
			case IN_EDGE:
			case OUT_EDGE:
				this->num_required_arrs = 1;
				break;
			case BOTH_EDGES:
				this->num_required_arrs = 2;
				break;
			default:
				assert(0);
		}
	}

	virtual void run(page_byte_array &arr);

	virtual bool has_completed() {
		return num_fetched_arrs == num_required_arrs;
	}
};

template<class compute_type>
class vertex_compute_allocator: public compute_allocator
{
	class compute_initiator: public obj_initiator<compute_type>
	{
		graph_engine *graph;
		vertex_compute_allocator<compute_type> *alloc;
	public:
		compute_initiator(graph_engine *graph,
				vertex_compute_allocator<compute_type> *alloc) {
			this->graph = graph;
			this->alloc = alloc;
		}

		virtual void init(compute_type *obj) {
			new (obj) compute_type(graph, alloc);
		}
	};

	class compute_destructor: public obj_destructor<compute_type>
	{
	public:
		void destroy(compute_type *obj) {
			obj->~compute_type();
		}
	};

	obj_allocator<compute_type> allocator;
public:
	vertex_compute_allocator(graph_engine *graph, thread *t): allocator(
			"vertex-compute-allocator", t->get_node_id(), false, 1024 * 1024,
			params.get_max_obj_alloc_size(),
			typename obj_initiator<compute_type>::ptr(new compute_initiator(graph, this)),
			typename obj_destructor<compute_type>::ptr(new compute_destructor())) {
	}

	virtual user_compute *alloc() {
		return allocator.alloc_obj();
	}

	virtual void free(user_compute *obj) {
		allocator.free((compute_type *) obj);
	}
};

#endif
