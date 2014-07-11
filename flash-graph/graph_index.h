#ifndef __GRAPH_INDEX_H__
#define __GRAPH_INDEX_H__

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

#include <boost/foreach.hpp>

#include "vertex.h"
#include "partitioner.h"
#include "vertex_program.h"
#include "graph_file_header.h"

class compute_vertex;

/**
 * \brief This file contains a set of graph index implementation.
 * The graph index is the in-memory container to keep all vertices of a graph. <br>
 *
 * It's different from `vertex_index`, which is an index to the adjacency
 * list of vertices in the external memory. **NOTE: This is class is mostly used internally**
 */
class graph_index
{
	graph_index(const graph_index &);
	graph_index &operator=(const graph_index &);
public:
	typedef std::shared_ptr<graph_index> ptr; /** Type used to access object*/
    
	class const_iterator
	{
		vertex_id_t id;
		graph_index *index;
	public:
		const_iterator(graph_index *index, vertex_id_t id) {
			this->index = index;
			this->id = id;
		}

		const compute_vertex &operator*() const {
			return index->get_vertex(id);
		}

		const_iterator &operator++() {
			id++;
			return *this;
		}

		bool operator==(const const_iterator &it) const {
			return id == it.id;
		}

		bool operator!=(const const_iterator &it) const {
			return id != it.id;
		}
	};

	graph_index() {
	}

	virtual ~graph_index() {
	}

	virtual void init(int num_threads, int num_nodes) {
	}

	virtual size_t get_vertices(const vertex_id_t ids[], int num,
			compute_vertex *v_buf[]) const = 0;
	virtual size_t get_vertices(int part_id, const local_vid_t ids[], int num,
			compute_vertex *v_buf[]) const = 0;
	virtual compute_vertex &get_vertex(vertex_id_t id) = 0;
	virtual compute_vertex &get_vertex(int part_id, local_vid_t id) = 0;

	virtual vertex_id_t get_max_vertex_id() const = 0;

	virtual vertex_id_t get_min_vertex_id() const = 0;

	virtual size_t get_num_vertices() const = 0;

	virtual const graph_partitioner &get_partitioner() const = 0;

	virtual vertex_program::ptr create_def_vertex_program() const = 0;

	const_iterator begin() const {
		return const_iterator((graph_index *) this, this->get_min_vertex_id());
	}

	const_iterator end() const {
		return const_iterator((graph_index *) this,
				this->get_max_vertex_id() + 1);
	}

	virtual std::string get_index_file() const = 0;
};

template<class vertex_type>
class NUMA_graph_index;

template<class vertex_type>
class NUMA_local_graph_index: public graph_index
{
	int node_id;
	int part_id;
	size_t tot_num_vertices;
	size_t num_vertices;
	vertex_type *vertex_arr;
	const graph_partitioner &partitioner;

	NUMA_local_graph_index(const graph_partitioner &_partitioner, int part_id,
			int node_id, size_t tot_num_vertices): partitioner(_partitioner) {
		this->part_id = part_id;
		this->tot_num_vertices = tot_num_vertices;
		this->num_vertices = partitioner.get_part_size(part_id,
				tot_num_vertices);
		this->node_id = node_id;
		vertex_arr = NULL;
	}
public:
	~NUMA_local_graph_index() {
		if (vertex_arr)
			free_large(vertex_arr, sizeof(vertex_arr[0]) * num_vertices);
	}

	void init() {
		vertex_arr = (vertex_type *) malloc_large(
				sizeof(vertex_arr[0]) * num_vertices);
		assert(vertex_arr);
		std::vector<vertex_id_t> local_ids;
		local_ids.reserve(num_vertices);
		partitioner.get_all_vertices_in_part(this->part_id, tot_num_vertices,
				local_ids);
		assert(local_ids.size() == num_vertices);
		BOOST_FOREACH(vertex_id_t vid, local_ids) {
			int part_id;
			off_t part_off;
			partitioner.map2loc(vid, part_id, part_off);
			assert(this->part_id == part_id);
			new (vertex_arr + part_off) vertex_type(vid);
		}
		if (local_ids.size() < num_vertices) {
			assert(local_ids.size() == num_vertices - 1);
			new (vertex_arr + local_ids.size()) vertex_type(INVALID_VERTEX_ID);
		}
	}

	virtual size_t get_vertices(const vertex_id_t ids[], int num,
			compute_vertex *v_buf[]) const {
		assert(0);
	}

	virtual compute_vertex &get_vertex(vertex_id_t id) {
		int part_id;
		off_t part_off;
		partitioner.map2loc(id, part_id, part_off);
		return vertex_arr[part_off];
	}

	virtual size_t get_vertices(int part_id, const local_vid_t ids[], int num,
			compute_vertex *v_buf[]) const {
		assert(0);
	}

	virtual compute_vertex &get_vertex(int part_id, local_vid_t id) {
		assert(0);
	}

	virtual vertex_id_t get_max_vertex_id() const {
		return vertex_arr[num_vertices - 1].get_id();
	}

	virtual vertex_id_t get_min_vertex_id() const {
		return vertex_arr[0].get_id();
	}

	virtual size_t get_num_vertices() const {
		return num_vertices;
	}

	int get_node_id() const {
		return node_id;
	}

	virtual const graph_partitioner &get_partitioner() const {
		return partitioner;
	}

	virtual vertex_program::ptr create_def_vertex_program() const {
		assert(0);
		return vertex_program::ptr();
	}

	virtual std::string get_index_file() const {
		assert(0);
	}

	friend class NUMA_graph_index<vertex_type>;
};

static inline size_t get_min_ext_mem_vertex_size(graph_type type)
{
	switch(type) {
		case graph_type::DIRECTED:
			return sizeof(ext_mem_directed_vertex);
		case graph_type::UNDIRECTED:
			return sizeof(ext_mem_undirected_vertex);
		case graph_type::TS_DIRECTED:
			return sizeof(ts_ext_mem_directed_vertex);
		case graph_type::TS_UNDIRECTED:
			assert(0);
			break;
		default:
			assert(0);
	}
}

template<class vertex_type>
class NUMA_graph_index: public graph_index
{
	std::string index_file;
	graph_header header;

	vertex_id_t max_vertex_id;
	vertex_id_t min_vertex_id;
	std::unique_ptr<range_graph_partitioner> partitioner;
	// A graph index per thread
	std::vector<std::unique_ptr<NUMA_local_graph_index<vertex_type> > > index_arr;

	class init_thread: public thread
	{
		NUMA_local_graph_index<vertex_type> *index;
	public:
		init_thread(NUMA_local_graph_index<vertex_type> *index): thread(
				"index-init-thread", index->get_node_id()) {
			this->index = index;
		}

		virtual void run() {
			index->init();
			this->stop();
		}
	};

	NUMA_graph_index(const std::string &index_file) {
		this->index_file = index_file;
	}
public:
	static graph_index::ptr create(const std::string &index_file) {
		return graph_index::ptr(new NUMA_graph_index<vertex_type>(index_file));
	}

	void init(int num_threads, int num_nodes) {
		partitioner = std::unique_ptr<range_graph_partitioner>(
				new range_graph_partitioner(num_threads));

		file_io_factory::shared_ptr index_factory = create_io_factory(
				index_file, GLOBAL_CACHE_ACCESS);
		io_interface::ptr io = index_factory->create_io(thread::get_curr_thread());
		io->access((char *) &header, 0, sizeof(header), READ);

		// Construct the indices.
		for (int i = 0; i < num_threads; i++) {
			index_arr.emplace_back(new NUMA_local_graph_index<vertex_type>(
						// The partitions are assigned to worker threads.
						// The memory used to store the partitions should
						// be on the same NUMA as the worker threads.
						*partitioner, i, i % num_nodes,
						header.get_num_vertices()));
		}

		std::vector<init_thread *> threads(num_threads);
		for (int i = 0; i < num_threads; i++) {
			threads[i] = new init_thread(index_arr[i].get());
			threads[i]->start();
		}
		for (int i = 0; i < num_threads; i++) {
			threads[i]->join();
			delete threads[i];
		}

		min_vertex_id = 0;
		max_vertex_id = header.get_num_vertices() - 1;

		printf("There are %ld vertices\n", header.get_num_vertices());
	}

	virtual size_t get_vertices(const vertex_id_t ids[], int num,
			compute_vertex *v_buf[]) const {
		for (int i = 0; i < num; i++) {
			vertex_id_t id = ids[i];
			int part_id;
			off_t part_off;
			partitioner->map2loc(id, part_id, part_off);
			v_buf[i] = &index_arr[part_id]->vertex_arr[part_off];
		}
		return num;
	}

	virtual compute_vertex &get_vertex(vertex_id_t id) {
		int part_id;
		off_t part_off;
		partitioner->map2loc(id, part_id, part_off);
		return index_arr[part_id]->vertex_arr[part_off];
	}

	virtual size_t get_vertices(int part_id, const local_vid_t ids[], int num,
			compute_vertex *v_buf[]) const {
		for (int i = 0; i < num; i++)
			v_buf[i] = &index_arr[part_id]->vertex_arr[ids[i].id];
		return num;
	}

	virtual compute_vertex &get_vertex(int part_id, local_vid_t id) {
		return index_arr[part_id]->vertex_arr[id.id];
	}

	virtual vertex_id_t get_max_vertex_id() const {
		return max_vertex_id;
	}

	virtual vertex_id_t get_min_vertex_id() const {
		return min_vertex_id;
	}

	virtual size_t get_num_vertices() const {
		return header.get_num_vertices();
	}

	virtual const graph_partitioner &get_partitioner() const {
		return *partitioner;
	}

	virtual vertex_program::ptr create_def_vertex_program(
			) const {
		return vertex_program::ptr(new vertex_program_impl<vertex_type>());
	}

	virtual std::string get_index_file() const {
		return index_file;
	}
};

#if 0
template<class vertex_type>
class graph_index_impl: public graph_index
{
	// This contains the vertices with edges.
	std::vector<vertex_type> vertices;
	size_t min_vertex_size;
	
	graph_index_impl(const std::string &index_file) {
		vertex_index *index = load_vertex_index(index_file);
		this->min_vertex_size = get_min_ext_mem_vertex_size(
				index->get_graph_header().get_graph_type());
		size_t num_vertices = index->get_num_vertices();
		size_t num_non_empty = 0;
		vertices.resize(num_vertices);
		for (size_t i = 0; i < num_vertices; i++) {
			vertices[i] = vertex_type(i, index);
			if (vertices[i].get_ext_mem_size() > min_vertex_size) {
				num_non_empty++;
			}
		}
		printf("There are %ld vertices and %ld non-empty, vertex array capacity: %ld\n",
				num_vertices, num_non_empty, vertices.capacity());
		vertex_index::destroy(index);
	}
public:
	static graph_index *create(const std::string &index_file) {
		return new graph_index_impl<vertex_type>(index_file);
	}

	virtual compute_vertex &get_vertex(vertex_id_t id) {
		return vertices[id];
	}

	virtual size_t get_num_vertices() const {
		return vertices.size();
	}

	virtual vertex_id_t get_max_vertex_id() const {
		return vertices.back().get_id();
	}

	virtual vertex_id_t get_min_vertex_id() const {
		return vertices.front().get_id();
	}
};
#endif

#endif
