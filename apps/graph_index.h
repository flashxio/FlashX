#ifndef __GRAPH_INDEX_H__
#define __GRAPH_INDEX_H__

/**
 * Copyright 2013 Da Zheng
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

#include <boost/foreach.hpp>

#include "vertex.h"
#include "partitioner.h"
#include "vertex_program.h"

class compute_vertex;

vertex_index *load_vertex_index(const std::string &index_file);

/**
 * This file contains a set of graph index implementation.
 * The graph index not only indexes vertices in a graph, but also is
 * the in-memory container to keep all vertices of a graph.
 * It's different from vertex_index, which is an index to the adjacency
 * list of vertices in the external memory.
 */

class graph_index
{
	graph_index(const graph_index &);
	graph_index &operator=(const graph_index &);
public:
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

	virtual compute_vertex &get_vertex(vertex_id_t id) = 0;
	virtual const in_mem_vertex_info &get_vertex_info(vertex_id_t id) const = 0;

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
};

template<class vertex_type>
class NUMA_graph_index;

template<class vertex_type>
class NUMA_local_graph_index: public graph_index
{
	size_t min_vertex_size;
	int node_id;
	int part_id;
	size_t tot_num_vertices;
	size_t num_vertices;
	size_t num_non_empty;
	vertex_type *vertex_arr;
	in_mem_vertex_info *info_arr;
	graph_partitioner *partitioner;
	const vertex_index *index;

	NUMA_local_graph_index(const vertex_index *index,
			graph_partitioner *partitioner, int part_id, int node_id,
			size_t tot_num_vertices, size_t min_vertex_size) {
		this->part_id = part_id;
		this->index = index;
		this->tot_num_vertices = tot_num_vertices;
		this->num_vertices = partitioner->get_part_size(part_id,
				tot_num_vertices);
		this->num_non_empty = 0;
		this->partitioner = partitioner;
		this->min_vertex_size = min_vertex_size;
		this->node_id = node_id;
		vertex_arr = (vertex_type *) numa_alloc_onnode(
				sizeof(vertex_arr[0]) * num_vertices, node_id);
		info_arr = (in_mem_vertex_info *) numa_alloc_onnode(
				sizeof(info_arr[0]) * num_vertices, node_id);
	}
public:
	~NUMA_local_graph_index() {
		if (vertex_arr)
			numa_free(vertex_arr, sizeof(vertex_arr[0]) * num_vertices);
		if (info_arr)
			numa_free(info_arr, sizeof(info_arr[0]) * num_vertices);
	}

	void init() {
		std::vector<vertex_id_t> local_ids;
		local_ids.reserve(num_vertices);
		partitioner->get_all_vertices_in_part(this->part_id, tot_num_vertices,
				local_ids);
		assert(local_ids.size() == num_vertices);
		BOOST_FOREACH(vertex_id_t vid, local_ids) {
			int part_id;
			off_t part_off;
			partitioner->map2loc(vid, part_id, part_off);
			assert(this->part_id == part_id);
			new (vertex_arr + part_off) vertex_type(vid, index);
			new (info_arr + part_off) in_mem_vertex_info(vid, index);
			if (info_arr[part_off].get_ext_mem_size() > min_vertex_size) {
				num_non_empty++;
			}
		}
		if (local_ids.size() < num_vertices) {
			assert(local_ids.size() == num_vertices - 1);
			new (vertex_arr + local_ids.size()) vertex_type();
			new (info_arr + local_ids.size()) in_mem_vertex_info();
		}
	}

	virtual const in_mem_vertex_info &get_vertex_info(vertex_id_t id) const {
		int part_id;
		off_t part_off;
		partitioner->map2loc(id, part_id, part_off);
		return info_arr[part_off];
	}

	virtual compute_vertex &get_vertex(vertex_id_t id) {
		int part_id;
		off_t part_off;
		partitioner->map2loc(id, part_id, part_off);
		return vertex_arr[part_off];
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

	size_t get_num_non_empty() const {
		return num_non_empty;
	}

	int get_node_id() const {
		return node_id;
	}

	virtual const graph_partitioner &get_partitioner() const {
		return *partitioner;
	}

	virtual vertex_program::ptr create_def_vertex_program() const {
		assert(0);
		return vertex_program::ptr();
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
	vertex_id_t max_vertex_id;
	vertex_id_t min_vertex_id;
	size_t num_vertices;
	range_graph_partitioner partitioner;
	// A graph index per thread
	std::vector<NUMA_local_graph_index<vertex_type> *> index_arr;

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

	NUMA_graph_index(const std::string &index_file,
			int num_threads, int num_nodes): partitioner(num_threads) {
		vertex_index *index = load_vertex_index(index_file);
		size_t min_vertex_size = get_min_ext_mem_vertex_size(
				index->get_graph_header().get_graph_type());
		num_vertices = index->get_num_vertices();
		// Construct the indices.
		for (int i = 0; i < num_threads; i++) {
			index_arr.push_back(new NUMA_local_graph_index<vertex_type>(
						// The partitions are assigned to worker threads.
						// The memory used to store the partitions should
						// be on the same NUMA as the worker threads.
						index, &partitioner, i, i % num_nodes,
						num_vertices, min_vertex_size));
		}

		std::vector<init_thread *> threads(num_threads);
		for (int i = 0; i < num_threads; i++) {
			threads[i] = new init_thread(index_arr[i]);
			threads[i]->start();
		}
		size_t num_non_empty = 0;
		for (int i = 0; i < num_threads; i++) {
			threads[i]->join();
			delete threads[i];
			num_non_empty += index_arr[i]->get_num_non_empty();
		}

		min_vertex_id = 0;
		max_vertex_id = num_vertices - 1;

		printf("There are %ld vertices and %ld non-empty\n", num_vertices,
				num_non_empty);
		vertex_index::destroy(index);
	}
public:
	static graph_index *create(const std::string &index_file,
			int num_threads, int num_nodes) {
		return new NUMA_graph_index<vertex_type>(index_file,
				num_threads, num_nodes);
	}

	virtual compute_vertex &get_vertex(vertex_id_t id) {
		int part_id;
		off_t part_off;
		partitioner.map2loc(id, part_id, part_off);
		return index_arr[part_id]->vertex_arr[part_off];
	}

	virtual const in_mem_vertex_info &get_vertex_info(vertex_id_t id) const {
		int part_id;
		off_t part_off;
		partitioner.map2loc(id, part_id, part_off);
		return index_arr[part_id]->info_arr[part_off];
	}

	virtual vertex_id_t get_max_vertex_id() const {
		return max_vertex_id;
	}

	virtual vertex_id_t get_min_vertex_id() const {
		return min_vertex_id;
	}

	virtual size_t get_num_vertices() const {
		return num_vertices;
	}

	virtual const graph_partitioner &get_partitioner() const {
		return partitioner;
	}

	virtual vertex_program::ptr create_def_vertex_program(
			) const {
		return vertex_program::ptr(new vertex_program_impl<vertex_type>());
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
