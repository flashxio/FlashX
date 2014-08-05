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
		bool operator()(const in_mem_vertex_info &info1,
				const in_mem_vertex_info &info2) {
			return info1.get_ext_mem_off() > info2.get_ext_mem_off();
		}
	};

	// TODO use the embedded array as the container.
	std::priority_queue<in_mem_vertex_info, std::vector<in_mem_vertex_info>,
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
	void issue_io_request(const in_mem_vertex_info &info);

	/*
	 * Complete a request for the adjacency list.
	 */
	void complete_request();

	/*
	 * The methods below deal with requesting # edges of vertices.
	 */

	/*
	 * This is a callback function. When the vertex index gets the number of
	 * edges, it notifies the vertex_compute of this information.
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

class part_directed_vertex_compute;

class directed_vertex_compute: public vertex_compute
{
	class part_request_info
	{
		directed_vertex_request req;
		in_mem_directed_vertex_info info;
	public:
		part_request_info() {
		}

		part_request_info(const directed_vertex_request &req,
				const in_mem_directed_vertex_info &info) {
			this->req = req;
			this->info = info;
		}

		void init(const directed_vertex_request &req,
				const in_mem_directed_vertex_info &info) {
			this->req = req;
			this->info = info;
		}

		void reset() {
			req = directed_vertex_request();
			info = in_mem_directed_vertex_info();
		}

		const directed_vertex_request &get_request() const {
			return req;
		}

		const in_mem_directed_vertex_info &get_info() const {
			return info;
		}

		bool is_valid() const {
			return req.is_valid();
		}

		off_t get_offset() const {
			if (req.get_type() == edge_type::IN_EDGE)
				return info.get_ext_mem_off()
					+ ext_mem_directed_vertex::get_header_size();
			else
				return info.get_ext_mem_off()
					+ ext_mem_directed_vertex::get_header_size()
					+ info.get_num_in_edges() * sizeof(vertex_id_t);
		}

		size_t get_size() const {
			if (req.get_type() == edge_type::IN_EDGE)
				return info.get_num_in_edges() * sizeof(vertex_id_t);
			else
				return info.get_num_out_edges() * sizeof(vertex_id_t);
		}

		bool operator<(const part_request_info &info) const {
			return this->req.get_id() > info.req.get_id();
		}
	};

	std::priority_queue<part_request_info> reqs;
	part_request_info part_req;

	request_range generate_request(const directed_vertex_request &req,
			const in_mem_directed_vertex_info &info);
public:
	directed_vertex_compute(graph_engine *graph,
			compute_allocator *alloc): vertex_compute(graph, alloc) {
	}

	virtual void set_scan_dir(bool forward) {
		vertex_compute::set_scan_dir(forward);
	}

	virtual int has_requests();

	virtual request_range get_next_request();

	virtual void run(page_byte_array &);

	/*
	 * The requested part of the vertex may have no edges.
	 * We can notify the user immediately.
	 */
	void complete_empty_part(const directed_vertex_request &req);

	void request_partial_vertices(directed_vertex_request reqs[], size_t num);

	void issue_io_request(const directed_vertex_request &req,
			const in_mem_directed_vertex_info &info);

	/*
	 * This is a callback function. When the vertex index gets the number of
	 * edges, it notifies the vertex_compute of this information.
	 */
	void run_on_num_edges(vertex_id_t id, vsize_t num_in_edges,
			vsize_t num_out_edges);

	/*
	 * This methods accepts the requests from graph applications and issues
	 * the requests to the vertex index.
	 */
	void request_num_edges(vertex_id_t ids[], size_t num);
};

class part_directed_vertex_compute: public user_compute
{
	graph_engine *graph;
	compute_vertex_pointer comp_v;
	directed_vertex_compute *compute;
	directed_vertex_request req;
	int num_fetched;
public:
	part_directed_vertex_compute(graph_engine *graph,
			compute_allocator *alloc): user_compute(alloc) {
		this->graph = graph;
		num_fetched = 0;
		compute = NULL;
	}

	void init(compute_vertex_pointer v, directed_vertex_compute *compute,
			const directed_vertex_request &req) {
		this->comp_v = v;
		this->compute = compute;
		this->req = req;
	}

	virtual int serialize(char *buf, int size) const {
		return 0;
	}

	virtual int get_serialized_size() const {
		return 0;
	}

	virtual int has_requests() {
		return false;
	}

	virtual request_range get_next_request() {
		// It shouldn't be invoked.
		assert(0);
	}

	virtual void run(page_byte_array &);

	virtual bool has_completed() {
		return num_fetched > 0;
	}
};

#if 0
class part_ts_vertex_compute;

class ts_vertex_compute: public vertex_compute
{
	std::priority_queue<ts_vertex_request> reqs;
public:
	ts_vertex_compute(graph_engine *graph,
			compute_allocator *alloc): vertex_compute(graph, alloc) {
	}

	virtual void set_scan_dir(bool forward) {
		vertex_compute::set_scan_dir(forward);
	}

	virtual int has_requests() const {
		return vertex_compute::has_requests() || reqs.size() > 0;
	}

	virtual request_range get_next_request();

	void request_partial_vertices(ts_vertex_request reqs[], size_t num);
};

/**
 * Sometimes a vertex only needs to read part of its neighbors.
 * This class is to read part of a neighbor and pass the neighbor to
 * the specified vertex to perform computation.
 * An instance of the class only reads one neighbor.
 */
class part_ts_vertex_compute: public user_compute
{
	graph_engine *graph;
	// The vertex where computation should perform.
	compute_vertex *comp_v;
	ts_vertex_compute *ts_compute;
	const TS_page_vertex *required_vertex_header;
	// The part of the vertex will be read and passed to
	// the computation vertex.
	ts_vertex_request required_part;
	int num_issued;
	int num_fetched;
public:
	part_ts_vertex_compute(graph_engine *graph,
			compute_allocator *alloc): user_compute(alloc) {
		this->graph = graph;
		comp_v = NULL;
		ts_compute = NULL;
		required_vertex_header = NULL;
		num_issued = 0;
		num_fetched = 0;
	}

	void init(compute_vertex *v, ts_vertex_compute *compute,
			const ts_vertex_request &req) {
		comp_v = v;
		ts_compute = compute;
		required_part = req;
	}

	virtual int serialize(char *buf, int size) const {
		return 0;
	}

	virtual int get_serialized_size() const {
		return 0;
	}

	virtual int has_requests() const {
		return num_issued == 0;
	}

	virtual request_range get_next_request();

	virtual void run(page_byte_array &);

	virtual bool has_completed() const {
		return num_fetched > 0;
	}
};
#endif

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
