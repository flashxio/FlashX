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
	static const int BUF_SIZE = sizeof(directed_vertex_entry);
protected:
	char curr_buf[BUF_SIZE];
	char next_buf[BUF_SIZE];
	bool _has_next;
public:
	bool has_next() const {
		return _has_next;
	}

	virtual void move_next() = 0;
	virtual bool move_to(int idx) = 0;
	virtual int get_num_entries() const = 0;

	off_t get_curr_off() const {
		return ((vertex_offset *) curr_buf)->get_off();
	}

	vsize_t get_curr_vertex_size() const {
		return ((vertex_offset *) next_buf)->get_off()
			- ((vertex_offset *) curr_buf)->get_off();
	}

	vsize_t get_curr_num_in_edges() const {
		return ((directed_vertex_entry *) curr_buf)->get_num_in_edges();
	}

	vsize_t get_curr_num_out_edges() const {
		return ((directed_vertex_entry *) curr_buf)->get_num_out_edges();
	}
};

template<class EntryType>
class page_index_iterator_impl: public index_iterator
{
	page_byte_array::seq_const_iterator<EntryType> it;
public:
	page_index_iterator_impl(page_byte_array::seq_const_iterator<EntryType> &_it): it(_it) {
		assert(it.get_num_tot_entries() >= 2);
		assert(it.has_next());
		*(EntryType *) curr_buf = it.next();
		assert(it.has_next());
		*(EntryType *) next_buf = it.next();
		_has_next = true;
	}

	virtual void move_next() {
		_has_next = it.has_next();
		if (_has_next) {
			*(EntryType *) curr_buf = *(EntryType *) next_buf;
			*(EntryType *) next_buf = it.next();
		}
	}

	virtual bool move_to(int idx) {
		bool ret = it.move_to(idx);
		if (!ret) {
			_has_next = false;
			return false;
		}
		*(EntryType *) curr_buf = it.next();
		if (it.has_next()) {
			_has_next = true;
			*(EntryType *) next_buf = it.next();
		}
		else
			_has_next = false;
		return _has_next;
	}

	virtual int get_num_entries() const {
		return it.get_num_tot_entries();
	}
};

template<class EntryType>
class array_index_iterator_impl: public index_iterator
{
	EntryType *start;
	EntryType *p;
	EntryType *end;
public:
	array_index_iterator_impl(EntryType *start, EntryType *end) {
		this->start = start;
		this->p = start;
		this->end = end;
		assert(end - p >= 2);
		*(EntryType *) curr_buf = *p;
		p++;
		*(EntryType *) next_buf = *p;
		_has_next = true;
	}

	virtual void move_next() {
		*(EntryType *) curr_buf = *(EntryType *) next_buf;
		p++;
		_has_next = p < end;
		if (_has_next)
			*(EntryType *) next_buf = *p;
	}

	virtual bool move_to(int idx) {
		p = start + idx;
		if (p + 1 < end) {
			*(EntryType *) curr_buf = *p;
			*(EntryType *) next_buf = *(p + 1);
			_has_next = true;
		}
		else
			_has_next = false;
		return _has_next;
	}

	virtual int get_num_entries() const {
		return end - p;
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

	void clear() {
		this->id_range.first = INVALID_VERTEX_ID;
		this->id_range.second = INVALID_VERTEX_ID;
	}

	void init(vertex_id_t id) {
		assert(id_range.first == INVALID_VERTEX_ID);
		id_range.first = id;
		id_range.second = id + 1;
	}

	void add_vertex(vertex_id_t id) {
		assert(id_range.second - 1 <= id);
		id_range.second = id + 1;
	}

	bool empty() const {
		return id_range.first == INVALID_VERTEX_ID;
	}

	vertex_id_t get_first_vertex() const {
		return id_range.first;
	}

	vertex_id_t get_last_vertex() const {
		return id_range.second - 1;
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

/*
 * The function class of requesting a vertex.
 */
struct req_vertex_func
{
	void operator()(vertex_id_t vid, vertex_compute &compute,
			const index_iterator &it) const {
		in_mem_vertex_info info(vid, it.get_curr_off(),
				it.get_curr_vertex_size());
		compute.issue_io_request(info);
	}
};

/*
 * The function class of requesting part of a vertex.
 */
class req_part_vertex_func
{
	edge_type type;
public:
	req_part_vertex_func() {
		type = edge_type::NONE;
	}

	req_part_vertex_func(edge_type type) {
		this->type = type;
	}

	void operator()(vertex_id_t vid, vertex_compute &compute,
			const index_iterator &it) const {
		directed_vertex_compute &dcompute = (directed_vertex_compute &) compute;
		in_mem_directed_vertex_info info(vid, it.get_curr_off(),
				it.get_curr_vertex_size(), it.get_curr_num_in_edges(),
				it.get_curr_num_out_edges());
		directed_vertex_request req(vid, type);
		dcompute.issue_io_request(req, info);
	}
};

/*
 * The function class of requesting the number of edges of a vertex.
 */
struct req_edge_func
{
	void operator()(vertex_id_t vid, vertex_compute &compute,
			const index_iterator &it) const {
		compute.run_on_vertex_size(vid, it.get_curr_vertex_size());
	}
};

/*
 * The function class of requesting the number of directed edges of a vertex.
 */
struct req_directed_edge_func
{
	void operator()(vertex_id_t vid, vertex_compute &compute,
			const index_iterator &it) const {
		directed_vertex_compute &dcompute
			= (directed_vertex_compute &) compute;
		dcompute.run_on_num_edges(vid, it.get_curr_num_in_edges(),
				it.get_curr_num_out_edges());
	}
};

/*
 * These are the implementation for dense index computes.
 * The dense index computes are adjacent to each other. The next index compute
 * runs on the next vertex in the vertex ID space. There are no repeated vertex
 * IDs or gaps between vertices.
 */
class dense_vertex_compute: public index_compute
{
	static const int MAX_INDEX_COMPUTE_SIZE = 512;
protected:
	embedded_array<vertex_compute *> computes;
	int num_gets;
public:
	dense_vertex_compute(index_comp_allocator &_alloc): index_compute(_alloc) {
		num_gets = 0;
	}

	vertex_compute *get_compute(vertex_id_t id) const {
		assert(id >= get_first_vertex());
		return computes[id - get_first_vertex()];
	}

	vertex_compute *get_first_compute() const {
		assert(get_num_vertices() >= 1);
		return computes[0];
	}

	int get_num_vertices() const {
		return get_last_vertex() - get_first_vertex() + 1;
	}

	void init(vertex_id_t id, vertex_compute *compute) {
		index_compute::init(id);
		computes[0] = compute;
	}

	bool add_vertex(vertex_id_t id, vertex_compute *compute) {
		// We don't want to have too many requests in a compute.
		// Otherwise, we need to use use malloc to allocate memory.
		if (computes.get_capacity() <= get_num_vertices()) {
			if (computes.get_capacity() >= MAX_INDEX_COMPUTE_SIZE)
				return false;
			else
				computes.resize(computes.get_capacity() * 2);
		}

		// The next vertex has to be the ID of the previous vertex + 1.
		// There should be no space.
		if (get_last_vertex() + 1 == id) {
			index_compute::add_vertex(id);
			computes[get_num_vertices() - 1] = compute;
			return true;
		}
		else
			return false;
	}

	template<class Func>
	bool run_temp(vertex_id_t vid, index_iterator &it, Func &func) {
		while (it.has_next()) {
			num_gets++;
			func(vid, *get_compute(vid), it);
			vid++;
			it.move_next();
		}
		return num_gets == get_num_vertices();
	}
};

/*
 * This requests an entire vertex.
 */
class req_vertex_compute: public dense_vertex_compute
{
public:
	req_vertex_compute(index_comp_allocator &_alloc): dense_vertex_compute(_alloc) {
	}

	virtual bool run(vertex_id_t vid, index_iterator &it) {
		req_vertex_func func;
		return run_temp<req_vertex_func>(vid, it, func);
	}

	req_vertex_func get_func() const {
		return req_vertex_func();
	}
};

/*
 * This requests part of a vertex.
 */
class req_part_vertex_compute: public dense_vertex_compute
{
	edge_type type;
public:
	req_part_vertex_compute(
			index_comp_allocator &_alloc): dense_vertex_compute(_alloc) {
		this->type = edge_type::NONE;
	}

	edge_type get_type() const {
		return type;
	}

	void init(const directed_vertex_request &req, vertex_compute *compute) {
		dense_vertex_compute::init(req.get_id(), compute);
		this->type = req.get_type();
	}

	bool add_vertex(const directed_vertex_request &req, vertex_compute *compute) {
		assert(this->type == req.get_type());
		return dense_vertex_compute::add_vertex(req.get_id(), compute);
	}

	virtual bool run(vertex_id_t vid, index_iterator &it) {
		req_part_vertex_func func(type);
		return run_temp<req_part_vertex_func>(vid, it, func);
	}

	req_part_vertex_func get_func() const {
		return req_part_vertex_func(type);
	}
};

/*
 * This requests # edges of a vertex.
 */
class req_edge_compute: public dense_vertex_compute
{
public:
	req_edge_compute(index_comp_allocator &_alloc): dense_vertex_compute(_alloc) {
	}

	virtual bool run(vertex_id_t vid, index_iterator &it) {
		req_edge_func func;
		return run_temp<req_edge_func>(vid, it, func);
	}

	req_edge_func get_func() const {
		return req_edge_func();
	}
};

/*
 * This requests # directed edges of a vertex.
 */
class req_directed_edge_compute: public dense_vertex_compute
{
public:
	req_directed_edge_compute(index_comp_allocator &_alloc): dense_vertex_compute(
			_alloc) {
	}

	virtual bool run(vertex_id_t vid, index_iterator &it) {
		req_directed_edge_func func;
		return run_temp<req_directed_edge_func>(vid, it, func);
	}

	req_directed_edge_func get_func() const {
		return req_directed_edge_func();
	}
};

/*
 * These index computes handle the general cases. That is, there might be a gap
 * between vertices as well as repeated vertices.
 */
class general_vertex_compute: public index_compute
{
	static const int MAX_INDEX_COMPUTE_SIZE = 512;
	int num_vertices;
	embedded_array<vertex_id_t> ids;
	embedded_array<vertex_compute *> computes;
	int entry_size_log;

	off_t get_page(vertex_id_t id) const {
		return ROUND_PAGE(id << entry_size_log);
	}

	off_t get_last_page() const {
		return get_page(get_last_vertex());
	}
protected:
	int num_gets;

	int get_start_idx(vertex_id_t start_vid) const {
		if (start_vid == get_id(0))
			return 0;
		else {
			const vertex_id_t *loc = std::lower_bound(ids.data(),
					ids.data() + num_vertices, start_vid);
			assert(loc != ids.data() + num_vertices);
			return loc - ids.data();
		}
	}
public:
	general_vertex_compute(index_comp_allocator &_alloc,
			int entry_size_log): index_compute(_alloc) {
		num_vertices = 0;
		num_gets = 0;
		this->entry_size_log = entry_size_log;
	}

	vertex_compute *get_compute(int idx) const {
		return computes[idx];
	}

	vertex_compute *get_first_compute() const {
		assert(get_num_vertices() >= 1);
		return computes[0];
	}

	vertex_id_t get_id(int idx) const {
		return ids[idx];
	}

	int get_num_vertices() const {
		return num_vertices;
	}

	void init(vertex_id_t id, vertex_compute *compute) {
		index_compute::init(id);
		ids[0] = id;
		computes[0] = compute;
		num_vertices++;
	}

	void clear() {
		index_compute::clear();
		num_vertices = 0;
	}

	bool add_vertex(vertex_id_t id, vertex_compute *compute) {
		// We don't want to have too many requests in a compute.
		// Otherwise, we need to use use malloc to allocate memory.
		if (computes.get_capacity() <= get_num_vertices()) {
			if (computes.get_capacity() >= MAX_INDEX_COMPUTE_SIZE)
				return false;
			else {
				ids.resize(ids.get_capacity() * 2);
				computes.resize(computes.get_capacity() * 2);
			}
		}

		if (get_last_page() == get_page(id)
				|| get_last_page() + 1 == get_page(id)) {
			index_compute::add_vertex(id);
			computes[get_num_vertices()] = compute;
			ids[get_num_vertices()] = id;
			num_vertices++;
			return true;
		}
		else
			return false;
	}

	template<class Func>
	bool run_temp(vertex_id_t start_vid, index_iterator &it, Func &func) {
		int num_vertices = get_num_vertices();
		for (int i = get_start_idx(start_vid); i < num_vertices; i++) {
			vertex_id_t id = get_id(i);
			vertex_compute *compute = get_compute(i);
			assert(id >= start_vid);
			if (!it.move_to(id - start_vid)) {
				break;
			}
			func(id, *compute, it);
			num_gets++;
		}
		return num_gets == get_num_vertices();
	}
};

class genrq_vertex_compute: public general_vertex_compute
{
public:
	genrq_vertex_compute(index_comp_allocator &_alloc,
			int entry_size_log): general_vertex_compute(_alloc, entry_size_log) {
	}

	virtual bool run(vertex_id_t vid, index_iterator &it) {
		req_vertex_func func;
		return run_temp<req_vertex_func>(vid, it, func);
	}

	req_vertex_func get_func() const {
		return req_vertex_func();
	}
};

class genrq_part_vertex_compute: public general_vertex_compute
{
	edge_type type;
public:
	genrq_part_vertex_compute(index_comp_allocator &_alloc,
			int entry_size_log): general_vertex_compute(_alloc, entry_size_log) {
		this->type = edge_type::NONE;
	}

	edge_type get_type() const {
		return type;
	}

	void init(const directed_vertex_request &req, vertex_compute *compute) {
		general_vertex_compute::init(req.get_id(), compute);
		this->type = req.get_type();
	}

	bool add_vertex(const directed_vertex_request &req, vertex_compute *compute) {
		assert(this->type == req.get_type());
		return general_vertex_compute::add_vertex(req.get_id(), compute);
	}

	virtual bool run(vertex_id_t vid, index_iterator &it) {
		req_part_vertex_func func(type);
		return run_temp<req_part_vertex_func>(vid, it, func);
	}

	req_part_vertex_func get_func() const {
		return req_part_vertex_func(type);
	}
};

class genrq_edge_compute: public general_vertex_compute
{
public:
	genrq_edge_compute(index_comp_allocator &_alloc,
			int entry_size_log): general_vertex_compute(_alloc, entry_size_log) {
	}

	virtual bool run(vertex_id_t vid, index_iterator &it) {
		req_edge_func func;
		return run_temp<req_edge_func>(vid, it, func);
	}

	req_edge_func get_func() const {
		return req_edge_func();
	}
};

class genrq_directed_edge_compute: public general_vertex_compute
{
public:
	genrq_directed_edge_compute(index_comp_allocator &_alloc,
			int entry_size_log): general_vertex_compute(_alloc, entry_size_log) {
	}

	virtual bool run(vertex_id_t vid, index_iterator &it) {
		req_directed_edge_func func;
		return run_temp<req_directed_edge_func>(vid, it, func);
	}

	req_directed_edge_func get_func() const {
		return req_directed_edge_func();
	}
};

template<class Func>
class single_index_compute: public index_compute
{
	vertex_compute *compute;
	Func func;
public:
	single_index_compute(
			index_comp_allocator &alloc): index_compute(alloc) {
		compute = NULL;
	}

	template<class MultIndexCompute>
	void init(const MultIndexCompute &compute) {
		assert(compute.get_num_vertices() == 1);
		index_compute::init(compute.get_first_vertex());
		assert(this->compute == NULL);
		this->compute = compute.get_first_compute();
		this->func = compute.get_func();
	}

	virtual bool run(vertex_id_t start_vid, index_iterator &it) {
		assert(start_vid == get_first_vertex());
		assert(it.has_next());
		func(start_vid, *compute, it);
		return true;
	}
};

typedef single_index_compute<req_vertex_func> single_vertex_compute;
typedef single_index_compute<req_part_vertex_func> single_part_vertex_compute;
typedef single_index_compute<req_edge_func> single_edge_compute;
typedef single_index_compute<req_directed_edge_func> single_directed_edge_compute;

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
			"index-compute-allocator", t->get_node_id(), false, 1024 * 1024,
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

template<class compute_type>
class general_index_comp_allocator_impl: public index_comp_allocator
{
	class compute_initiator: public obj_initiator<compute_type>
	{
		general_index_comp_allocator_impl<compute_type> &alloc;
		int entry_size_log;
	public:
		compute_initiator(general_index_comp_allocator_impl<compute_type> &_alloc,
				int entry_size_log): alloc(_alloc) {
			this->entry_size_log = entry_size_log;
		}

		virtual void init(compute_type *obj) {
			new (obj) compute_type(alloc, entry_size_log);
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
	general_index_comp_allocator_impl(thread *t, int entry_size_log): allocator(
			"sparse-index-compute-allocator", t->get_node_id(), false, 1024 * 1024,
			params.get_max_obj_alloc_size(),
			typename obj_initiator<compute_type>::ptr(new compute_initiator(*this,
					entry_size_log)),
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
	index_comp_allocator_impl<req_vertex_compute> *req_vertex_comp_alloc;
	index_comp_allocator_impl<req_part_vertex_compute> *req_part_vertex_comp_alloc;
	index_comp_allocator_impl<req_edge_compute> *req_edge_comp_alloc;
	index_comp_allocator_impl<req_directed_edge_compute> *req_directed_edge_comp_alloc;

	general_index_comp_allocator_impl<genrq_vertex_compute> *genrq_vertex_comp_alloc;
	general_index_comp_allocator_impl<genrq_part_vertex_compute> *genrq_part_vertex_comp_alloc;
	general_index_comp_allocator_impl<genrq_edge_compute> *genrq_edge_comp_alloc;
	general_index_comp_allocator_impl<genrq_directed_edge_compute> *genrq_directed_edge_comp_alloc;

	index_comp_allocator_impl<single_vertex_compute> *single_vertex_comp_alloc;
	index_comp_allocator_impl<single_part_vertex_compute> *single_part_vertex_comp_alloc;
	index_comp_allocator_impl<single_edge_compute> *single_edge_comp_alloc;
	index_comp_allocator_impl<single_directed_edge_compute> *single_directed_edge_comp_alloc;

	typedef std::pair<vertex_id_t, vertex_compute *> id_compute_t;
	typedef std::pair<directed_vertex_request, directed_vertex_compute *> directed_compute_t;
	std::vector<id_compute_t> vertex_comps;
	std::vector<directed_compute_t> part_vertex_comps[edge_type::NUM_TYPES];
	std::vector<id_compute_t> edge_comps;
	std::vector<id_compute_t> directed_edge_comps;

	struct id_compute_lesseq
	{
		bool operator()(const id_compute_t &comp1, const id_compute_t &comp2) {
			return comp1.first <= comp2.first;
		}
	};

	struct id_compute_less
	{
		bool operator()(const id_compute_t &comp1, const id_compute_t &comp2) {
			return comp1.first < comp2.first;
		}
	};

	struct directed_compute_lesseq
	{
		bool operator()(const directed_compute_t &comp1, const directed_compute_t &comp2) {
			return comp1.first.get_id() <= comp2.first.get_id();
		}
	};

	struct directed_compute_less
	{
		bool operator()(const directed_compute_t &comp1, const directed_compute_t &comp2) {
			return comp1.first.get_id() < comp2.first.get_id();
		}
	};

	vertex_index_reader::ptr index_reader;

	void flush_computes();

	static int get_index_entry_size_log(bool directed) {
		return log2(get_index_entry_size(
					directed ? graph_type::DIRECTED : graph_type::UNDIRECTED));
	}

	void init(thread *t, bool directed) {
		req_vertex_comp_alloc
			= new index_comp_allocator_impl<req_vertex_compute>(t);
		req_part_vertex_comp_alloc
			= new index_comp_allocator_impl<req_part_vertex_compute>(t);
		req_edge_comp_alloc
			= new index_comp_allocator_impl<req_edge_compute>(t);
		req_directed_edge_comp_alloc
			= new index_comp_allocator_impl<req_directed_edge_compute>(t);

		genrq_vertex_comp_alloc
			= new general_index_comp_allocator_impl<genrq_vertex_compute>(
					t, get_index_entry_size_log(directed));
		genrq_part_vertex_comp_alloc
			= new general_index_comp_allocator_impl<genrq_part_vertex_compute>(
					t, get_index_entry_size_log(directed));
		genrq_edge_comp_alloc
			= new general_index_comp_allocator_impl<genrq_edge_compute>(
					t, get_index_entry_size_log(directed));
		genrq_directed_edge_comp_alloc
			= new general_index_comp_allocator_impl<genrq_directed_edge_compute>(
					t, get_index_entry_size_log(directed));

		single_vertex_comp_alloc
			= new index_comp_allocator_impl<single_vertex_compute>(t);
		single_part_vertex_comp_alloc
			= new index_comp_allocator_impl<single_part_vertex_compute>(t);
		single_edge_comp_alloc
			= new index_comp_allocator_impl<single_edge_compute>(t);
		single_directed_edge_comp_alloc
			= new index_comp_allocator_impl<single_directed_edge_compute>(t);
	}

	simple_index_reader(vertex_index::ptr index, bool directed, thread *t) {
		init(t, directed);
		index_reader = vertex_index_reader::create(index, directed);
	}

	simple_index_reader(io_interface::ptr io, bool directed, thread *t) {
		init(t, directed);
		index_reader = vertex_index_reader::create(io, directed);
	}

	/*
	 * This gets the number of partitions in the vector. Each partition
	 * has all vertices adjacent to each other. There are no gaps or repeated
	 * vertices in a partition.
	 */
	template<class EntryType, class GetId>
	static int get_num_parts(const std::vector<EntryType> &vec)
	{
		GetId get_id;
		if (get_id(vec.back().first) - get_id(vec.front().first) + 1
				== vec.size())
			return 1;

		vertex_id_t id = get_id(vec[0].first);
		int num_parts = 1;
		for (size_t i = 1; i < vec.size(); i++) {
			if (id + 1 == get_id(vec[i].first))
				id++;
			else {
				num_parts++;
				id = get_id(vec[i].first);
			}
		}
		return num_parts;
	}

	static bool is_dense(const std::vector<id_compute_t> &vec)
	{
		if (vec.size() <= 1)
			return true;
		struct get_id {
			vertex_id_t operator()(vertex_id_t id) {
				return id;
			}
		};
		int num_parts = get_num_parts<id_compute_t, get_id>(vec);
		if (num_parts == 1)
			return true;
		return vec.size() / num_parts >= 16;
	}

	static bool is_dense(const std::vector<directed_compute_t> &vec)
	{
		if (vec.size() <= 1)
			return true;
		struct get_id {
			vertex_id_t operator()(const directed_vertex_request &req) {
				return req.get_id();
			}
		};
		int num_parts = get_num_parts<directed_compute_t, get_id>(vec);
		if (num_parts == 1)
			return true;
		return vec.size() / num_parts >= 16;
	}

#if 0
	void request_vertex(vertex_id_t id, vertex_compute *compute) {
		if (!whole_compute->add_vertex(id, compute)) {
			index_reader->request_index(whole_compute);
			whole_compute = (req_vertex_compute *) req_vertex_comp_alloc.alloc();
			whole_compute->init(id, compute);
		}
	}

	void request_vertex(const directed_vertex_request &req,
			directed_vertex_compute *compute) {
		assert(req.get_type() != edge_type::NONE);
		if (!part_computes[req.get_type()]->add_vertex(req.get_id(), compute)) {
			index_reader->request_index(part_computes[req.get_type()]);
			part_computes[req.get_type()]
				= (req_part_vertex_compute *) req_part_vertex_comp_alloc.alloc();
			part_computes[req.get_type()]->init(req, compute);
		}
	}

	void request_num_edges(vertex_id_t id, vertex_compute *compute) {
		if (!edge_compute->add_vertex(id, compute)) {
			index_reader->request_index(edge_compute);
			edge_compute = (req_edge_compute *) req_edge_comp_alloc.alloc();
			edge_compute->init(id, compute);
		}
	}

	void request_num_directed_edges(vertex_id_t id,
			directed_vertex_compute *compute) {
		if (!directed_edge_compute->add_vertex(id, compute)) {
			if (directed_edge_compute->get_num_vertices() == 1) {
				single_directed_edge_compute *comp
					= (single_directed_edge_compute *) single_directed_edge_comp_alloc.alloc();
				comp->init(*directed_edge_compute);
				index_reader->request_index(comp);
				directed_edge_compute->clear();
			}
			else {
				index_reader->request_index(directed_edge_compute);
				directed_edge_compute
					= (req_directed_edge_compute *) req_directed_edge_comp_alloc.alloc();
			}
			directed_edge_compute->init(id, compute);
		}
	}
#endif
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
		delete req_vertex_comp_alloc;
		delete req_part_vertex_comp_alloc;
		delete req_edge_comp_alloc;
		delete req_directed_edge_comp_alloc;

		delete genrq_vertex_comp_alloc;
		delete genrq_part_vertex_comp_alloc;
		delete genrq_edge_comp_alloc;
		delete genrq_directed_edge_comp_alloc;

		delete single_vertex_comp_alloc;
		delete single_part_vertex_comp_alloc;
		delete single_edge_comp_alloc;
		delete single_directed_edge_comp_alloc;
	}

	void request_vertices(vertex_id_t ids[], int num, vertex_compute &compute) {
		for (int i = 0; i < num; i++)
			vertex_comps.push_back(id_compute_t(ids[i], &compute));
	}

	void request_vertices(const directed_vertex_request reqs[], int num,
			directed_vertex_compute &compute) {
		for (int i = 0; i < num; i++)
			part_vertex_comps[reqs[i].get_type()].push_back(
					directed_compute_t(reqs[i], &compute));
	}

	void request_num_edges(vertex_id_t ids[], int num, vertex_compute &compute) {
		for (int i = 0; i < num; i++)
			edge_comps.push_back(id_compute_t(ids[i], &compute));
	}

	void request_num_directed_edges(vertex_id_t ids[], int num,
			directed_vertex_compute &compute) {
		for (int i = 0; i < num; i++)
			directed_edge_comps.push_back(id_compute_t(ids[i], &compute));
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
