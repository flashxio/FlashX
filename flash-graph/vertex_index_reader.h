#ifndef __VERTEX_INDEX_READER__
#define __VERTEX_INDEX_READER__

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

#include "FG_basic_types.h"
#include "vertex.h"
#include "vertex_index.h"

class directed_vertex_request;
class index_comp_allocator;

/*
 * This iterates on the entries of the vertex index.
 * This can iterate on both undirected and directed vertex index.
 */
class index_iterator
{
	char *p;
	char *p_end;
	int entry_size;	// in bytes.
public:
	index_iterator(char *p, char *end, int entry_size) {
		this->p = p;
		this->p_end = end;
		this->entry_size = entry_size;
		assert(end > p && (end - p) % entry_size == 0);
	}

	bool has_next() const {
		// TO iterate over an array of n entries, the actual array
		// has n + 1 entries because we need to have the next entry
		// to compute the vertex size.
		return p < p_end - entry_size;
	}

	void move_next() {
		p += entry_size;
	}

	off_t get_curr_off() const {
		return ((vertex_offset *) p)->get_off();
	}

	vsize_t get_curr_vertex_size() const {
		return ((vertex_offset *) (p + entry_size))->get_off()
			- ((vertex_offset *) p)->get_off();
	}

	vsize_t get_curr_num_in_edges() const {
		return ((directed_vertex_entry *) p)->get_num_in_edges();
	}

	vsize_t get_curr_num_out_edges() const {
		return ((directed_vertex_entry *) p)->get_num_out_edges();
	}
};

/**
 * This interface defines the method invoked in vertex_index_reader.
 * It is designed to compute on multiple index entries. If we require
 * vertices that are adjacent in vertex ids, we should merge all of
 * these requests in a single index_compute.
 */
class index_compute
{
public:
	typedef std::pair<vertex_id_t, vertex_id_t> id_range_t;
private:
	index_comp_allocator &alloc;
	id_range_t id_range;
public:
	index_compute(index_comp_allocator &_alloc): alloc(_alloc) {
		this->id_range.first = INVALID_VERTEX_ID;
		this->id_range.second = INVALID_VERTEX_ID;
	}

	virtual ~index_compute() {
	}

	bool add_vertex(vertex_id_t id) {
		bool ret;
		if (id_range.first == INVALID_VERTEX_ID) {
			id_range.first = id;
			id_range.second = id + 1;
			ret = true;
		}
		else if (id_range.second == id) {
			id_range.second++;
			ret = true;
		}
		else
			ret = false;
		return ret;
	}

	int get_num_vertices() const {
		return id_range.second - id_range.first;
	}

	vertex_id_t get_first_vertex() const {
		return id_range.first;
	}

	const id_range_t &get_range() const {
		return id_range;
	}

	virtual bool run(vertex_id_t start_vid, index_iterator &it) = 0;

	virtual index_comp_allocator &get_allocator() const {
		return alloc;
	}
};

class index_comp_allocator
{
public:
	virtual ~index_comp_allocator() {
	}
	virtual index_compute *alloc() = 0;
	virtual void free(index_compute *compute) = 0;
};

/*
 * This interface reads vertex index from SSDs.
 * It accepts the requests of reading vertices or partial vertices as well
 * as other vertex information in the vertex index such as the number of edges.
 * This is a per-thread data structure.
 */
class vertex_index_reader
{
public:
	typedef std::shared_ptr<vertex_index_reader> ptr;

	static ptr create(vertex_index::ptr index, bool directed);
	static ptr create(io_interface::ptr io, bool directed);

	virtual ~vertex_index_reader() {
	}

	virtual void request_index(index_compute *compute) = 0;
	virtual void wait4complete(int num) = 0;
	virtual size_t get_num_pending_tasks() const = 0;
};

/*
 * The classes below implements 4 index_compute to request:
 *	vertex,
 *	part of vertex,
 *	# edges,
 *	# directed edges.
 */

class base_vertex_compute: public index_compute
{
protected:
	static const int MAX_COMPUTE_SIZE = 512;
	embedded_array<vertex_compute *> computes;
	int num_gets;
public:
	base_vertex_compute(index_comp_allocator &_alloc): index_compute(_alloc) {
		num_gets = 0;
	}

	vertex_compute *get_compute(vertex_id_t id) const {
		assert(id >= get_first_vertex());
		return computes[id - get_first_vertex()];
	}

	bool add_vertex(vertex_id_t id, vertex_compute *compute) {
		// We don't want to have too many requests in a compute.
		// Otherwise, we need to use use malloc to allocate memory.
		if (computes.get_capacity() <= get_num_vertices()) {
			if (computes.get_capacity() >= MAX_COMPUTE_SIZE)
				return false;
			else
				computes.resize(computes.get_capacity() * 2);
		}
		bool ret = index_compute::add_vertex(id);
		if (ret)
			computes[get_num_vertices() - 1] = compute;
		return ret;
	}
};

/*
 * This requests an entire vertex.
 */
class req_vertex_compute: public base_vertex_compute
{
public:
	req_vertex_compute(index_comp_allocator &_alloc): base_vertex_compute(_alloc) {
	}

	virtual bool run(vertex_id_t vid, index_iterator &it);
};

/*
 * This requests part of a vertex.
 */
class req_part_vertex_compute: public base_vertex_compute
{
	edge_type type;
public:
	req_part_vertex_compute(index_comp_allocator &_alloc): base_vertex_compute(_alloc) {
		this->type = edge_type::NONE;
	}

	void set_type(edge_type type) {
		this->type = type;
	}

	virtual bool run(vertex_id_t vid, index_iterator &it);
};

/*
 * This requests # edges of a vertex.
 */
class req_edge_compute: public base_vertex_compute
{
public:
	req_edge_compute(index_comp_allocator &_alloc): base_vertex_compute(_alloc) {
	}

	virtual bool run(vertex_id_t vid, index_iterator &it);
};

/*
 * This requests # directed edges of a vertex.
 */
class req_directed_edge_compute: public base_vertex_compute
{
public:
	req_directed_edge_compute(index_comp_allocator &_alloc): base_vertex_compute(
			_alloc) {
	}

	virtual bool run(vertex_id_t vid, index_iterator &it);
};

template<class compute_type>
class index_comp_allocator_impl: public index_comp_allocator
{
	class compute_initiator: public obj_initiator<compute_type>
	{
		index_comp_allocator_impl<compute_type> &alloc;
	public:
		compute_initiator(
				index_comp_allocator_impl<compute_type> &_alloc): alloc(_alloc) {
		}

		virtual void init(compute_type *obj) {
			new (obj) compute_type(alloc);
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
	index_comp_allocator_impl(thread *t): allocator(
			"index-compute-allocator", t->get_node_id(), 1024 * 1024,
			params.get_max_obj_alloc_size(),
			typename obj_initiator<compute_type>::ptr(new compute_initiator(*this)),
			typename obj_destructor<compute_type>::ptr(new compute_destructor())) {
	}

	virtual index_compute *alloc() {
		return allocator.alloc_obj();
	}

	virtual void free(index_compute *obj) {
		allocator.free((compute_type *) obj);
	}
};

/*
 * This is a helper class that requests
 *	vertex,
 *	part of vertex,
 *	# edges,
 *	# directed edges.
 */
class simple_index_reader
{
	req_vertex_compute *whole_compute;
	req_part_vertex_compute *part_computes[edge_type::NUM_TYPES];
	req_edge_compute *edge_compute;
	req_directed_edge_compute *directed_edge_compute;

	index_comp_allocator_impl<req_vertex_compute> req_vertex_comp_alloc;
	index_comp_allocator_impl<req_part_vertex_compute> req_part_vertex_comp_alloc;
	index_comp_allocator_impl<req_edge_compute> req_edge_comp_alloc;
	index_comp_allocator_impl<req_directed_edge_compute> req_directed_edge_comp_alloc;

	vertex_index_reader::ptr index_reader;

	void init() {
		whole_compute = (req_vertex_compute *) req_vertex_comp_alloc.alloc();
		part_computes[edge_type::NONE] = NULL;
		part_computes[edge_type::IN_EDGE]
			= (req_part_vertex_compute *) req_part_vertex_comp_alloc.alloc();
		part_computes[edge_type::IN_EDGE]->set_type(edge_type::IN_EDGE);
		part_computes[edge_type::OUT_EDGE]
			= (req_part_vertex_compute *) req_part_vertex_comp_alloc.alloc();
		part_computes[edge_type::OUT_EDGE]->set_type(edge_type::OUT_EDGE);
		part_computes[edge_type::BOTH_EDGES]
			= (req_part_vertex_compute *) req_part_vertex_comp_alloc.alloc();
		part_computes[edge_type::BOTH_EDGES]->set_type(edge_type::BOTH_EDGES);
		edge_compute = (req_edge_compute *) req_edge_comp_alloc.alloc();
		directed_edge_compute
			= (req_directed_edge_compute *) req_directed_edge_comp_alloc.alloc();
	}

	void flush_computes() {
		if (whole_compute->get_num_vertices() > 0) {
			index_reader->request_index(whole_compute);
			whole_compute = (req_vertex_compute *) req_vertex_comp_alloc.alloc();
		}
		for (int type = edge_type::IN_EDGE; type < edge_type::NUM_TYPES; type++) {
			if (part_computes[type]->get_num_vertices() > 0) {
				index_reader->request_index(part_computes[type]);
				part_computes[type]
					= (req_part_vertex_compute *) req_part_vertex_comp_alloc.alloc();
				edge_type etype = (edge_type) type;
				part_computes[type]->set_type(etype);
			}
		}
		if (edge_compute->get_num_vertices() > 0) {
			index_reader->request_index(edge_compute);
			edge_compute = (req_edge_compute *) req_edge_comp_alloc.alloc();
		}
		if (directed_edge_compute->get_num_vertices() > 0) {
			index_reader->request_index(directed_edge_compute);
			directed_edge_compute
				= (req_directed_edge_compute *) req_directed_edge_comp_alloc.alloc();
		}
	}

	simple_index_reader(vertex_index::ptr index, bool directed,
			thread *t): req_vertex_comp_alloc(t), req_part_vertex_comp_alloc(
			t), req_edge_comp_alloc(t), req_directed_edge_comp_alloc(t) {
		init();
		index_reader = vertex_index_reader::create(index, directed);
	}

	simple_index_reader(io_interface::ptr io, bool directed,
			thread *t): req_vertex_comp_alloc(t), req_part_vertex_comp_alloc(
			t), req_edge_comp_alloc(t), req_directed_edge_comp_alloc(t) {
		init();
		index_reader = vertex_index_reader::create(io, directed);
	}
public:
	typedef std::shared_ptr<simple_index_reader> ptr;

	static ptr create(io_interface::ptr io, bool directed, thread *t) {
		return ptr(new simple_index_reader(io, directed, t));
	}

	static ptr create(vertex_index::ptr index, bool directed, thread *t) {
		return ptr(new simple_index_reader(index, directed, t));
	}

	~simple_index_reader() {
		assert(get_num_pending_tasks() == 0);
	}

	void request_vertex(vertex_id_t id, vertex_compute &compute) {
		if (!whole_compute->add_vertex(id, &compute)) {
			index_reader->request_index(whole_compute);
			whole_compute = (req_vertex_compute *) req_vertex_comp_alloc.alloc();
			bool ret = whole_compute->add_vertex(id, &compute);
			assert(ret);
		}
	}

	void request_vertex(const directed_vertex_request &req,
			directed_vertex_compute &compute) {
		assert(req.get_type() != edge_type::NONE);
		if (!part_computes[req.get_type()]->add_vertex(req.get_id(), &compute)) {
			index_reader->request_index(part_computes[req.get_type()]);
			part_computes[req.get_type()]
				= (req_part_vertex_compute *) req_part_vertex_comp_alloc.alloc();
			part_computes[req.get_type()]->set_type(req.get_type());
			bool ret = part_computes[req.get_type()]->add_vertex(req.get_id(),
					&compute);
			assert(ret);
		}
	}

	void request_num_edges(vertex_id_t id, vertex_compute &compute) {
		if (!edge_compute->add_vertex(id, &compute)) {
			index_reader->request_index(edge_compute);
			edge_compute = (req_edge_compute *) req_edge_comp_alloc.alloc();
			bool ret = edge_compute->add_vertex(id, &compute);
			assert(ret);
		}
	}

	void request_num_directed_edges(vertex_id_t id,
			directed_vertex_compute &compute) {
		if (!directed_edge_compute->add_vertex(id, &compute)) {
			index_reader->request_index(directed_edge_compute);
			directed_edge_compute
				= (req_directed_edge_compute *) req_directed_edge_comp_alloc.alloc();
			bool ret = directed_edge_compute->add_vertex(id, &compute);
			assert(ret);
		}
	}

	void request_vertices(vertex_id_t ids[], int num, vertex_compute &compute) {
		for (int i = 0; i < num; i++)
			request_vertex(ids[i], compute);
	}

	void request_vertices(const directed_vertex_request reqs[], int num,
			directed_vertex_compute &compute) {
		for (int i = 0; i < num; i++)
			request_vertex(reqs[i], compute);
	}

	void request_num_edges(vertex_id_t ids[], int num, vertex_compute &compute) {
		for (int i = 0; i < num; i++)
			request_num_edges(ids[i], compute);
	}

	void request_num_directed_edges(vertex_id_t ids[], int num,
			directed_vertex_compute &compute) {
		for (int i = 0; i < num; i++)
			request_num_directed_edges(ids[i], compute);
	}

	void wait4complete(int num) {
		flush_computes();
		index_reader->wait4complete(num);
	}

	size_t get_num_pending_tasks() const {
		return index_reader->get_num_pending_tasks();
	}
};

#endif
