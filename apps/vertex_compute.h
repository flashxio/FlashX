#ifndef __VERTEX_COMPUTE_H__
#define __VERTEX_COMPUTE_H__

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
#include <algorithm>
#include <deque>

#include "io_interface.h"
#include "slab_allocator.h"

#include "vertex_request.h"

class worker_thread;
class graph_engine;
class compute_vertex;
class compute_directed_vertex;

/**
 * This callback is to process a vertex.
 */
class vertex_compute: public user_compute
{
	std::vector<vertex_id_t> requested_vertices;
	size_t fetch_idx;

protected:
	graph_engine *graph;

	// The thread that creates the vertex compute.
	worker_thread *issue_thread;
	compute_vertex *v;
	// The number of requested vertices that will be read in the user compute.
	size_t num_complete_issues;
	// The number of vertices read by the user compute.
	size_t num_complete_fetched;
public:
	vertex_compute(graph_engine *graph,
			compute_allocator *alloc): user_compute(alloc) {
		this->graph = graph;
		v = NULL;
		issue_thread = (worker_thread *) thread::get_curr_thread();
		num_complete_issues = 0;
		num_complete_fetched = 0;
		fetch_idx = 0;
	}

	void init(compute_vertex *v) {
		this->v = v;
	}

	virtual int serialize(char *buf, int size) const {
		return 0;
	}

	virtual int get_serialized_size() const {
		return 0;
	}

	virtual int has_requests() const {
		return fetch_idx < requested_vertices.size();
	}

	virtual request_range get_next_request();

	virtual void run(page_byte_array &);

	virtual bool has_completed() const {
		// If the user compute has got all requested data and it has
		// no more requests to issue, we can consider the user compute
		// has been completed.
		// NOTE: it's possible that requested data may not be passed to
		// this user compute, so we only count the requests that are going
		// to be passed to this user compute.
		return num_complete_issues == num_complete_fetched && !has_requests();
	}

	virtual void request_vertices(vertex_id_t ids[], size_t num);

	graph_engine &get_graph() {
		return *graph;
	}

	void complete_request();
};

class part_directed_vertex_compute;

class directed_vertex_compute: public vertex_compute
{
	std::vector<directed_vertex_request> reqs;
	size_t fetch_idx;
public:
	directed_vertex_compute(graph_engine *graph,
			compute_allocator *alloc): vertex_compute(graph, alloc) {
		fetch_idx = 0;
	}

	virtual int has_requests() const {
		return vertex_compute::has_requests() || fetch_idx < reqs.size();
	}

	virtual request_range get_next_request();

	void request_partial_vertices(directed_vertex_request reqs[], size_t num);
};

class part_directed_vertex_compute: public user_compute
{
	graph_engine *graph;
	compute_directed_vertex *comp_v;
	directed_vertex_compute *compute;
	directed_vertex_request req;
	int num_fetched;
public:
	part_directed_vertex_compute(graph_engine *graph,
			compute_allocator *alloc): user_compute(alloc) {
		this->graph = graph;
		num_fetched = 0;
		comp_v = NULL;
		compute = NULL;
	}

	void init(compute_directed_vertex *v, directed_vertex_compute *compute,
			const directed_vertex_request &req) {
		this->comp_v = v;
		this->compute = compute;
		compute->inc_ref();
		this->req = req;
	}

	virtual int serialize(char *buf, int size) const {
		return 0;
	}

	virtual int get_serialized_size() const {
		return 0;
	}

	virtual int has_requests() const {
		return false;
	}

	virtual request_range get_next_request() {
		// It shouldn't be invoked.
		assert(0);
	}

	virtual void run(page_byte_array &);

	virtual bool has_completed() const {
		return num_fetched > 0;
	}
};

class part_ts_vertex_compute;

class ts_vertex_compute: public vertex_compute
{
	std::vector<ts_vertex_request> reqs;
	size_t fetch_idx;
public:
	ts_vertex_compute(graph_engine *graph,
			compute_allocator *alloc): vertex_compute(graph, alloc) {
		fetch_idx = 0;
	}

	virtual int has_requests() const {
		return vertex_compute::has_requests() || fetch_idx < reqs.size();
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
		compute->inc_ref();
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
			"vertex-compute-allocator", t->get_node_id(), 1024 * 1024,
			params.get_max_obj_alloc_size(),
			// TODO memory leak here
			new compute_initiator(graph, this), new compute_destructor()) {
	}

	virtual user_compute *alloc() {
		return allocator.alloc_obj();
	}

	virtual void free(user_compute *obj) {
		allocator.free((compute_type *) obj);
	}
};

#endif
