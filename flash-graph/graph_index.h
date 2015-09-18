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

#include <algorithm>
#include <boost/foreach.hpp>
#include <boost/format.hpp>

#include "log.h"
#include "vertex.h"
#include "partitioner.h"
#include "vertex_program.h"
#include "graph_file_header.h"
#include "vertex_pointer.h"

class compute_vertex;
class part_compute_vertex;

/*
 * A pointer to the vertically partitioned vertices in the graph index.
 */
class vpart_vertex_pointer
{
	vertex_id_t id;
	vsize_t off;
public:
	vpart_vertex_pointer() {
		id = INVALID_VERTEX_ID;
		off = -1;
	}

	vpart_vertex_pointer(vertex_id_t id, vsize_t off) {
		this->id = id;
		this->off = off;
	}

	vertex_id_t get_vertex_id() const {
		return id;
	}

	vsize_t get_off() const {
		return off;
	}
};

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
    
	graph_index() {
	}

	virtual ~graph_index() {
	}

	virtual void init(int num_threads, int num_nodes) {
	}
	virtual void init_vparts(int hpart_id, int num_vparts,
			std::vector<vertex_id_t> &ids) = 0;

	virtual size_t get_vertices(const vertex_id_t ids[], int num,
			compute_vertex *v_buf[]) const = 0;
	virtual size_t get_vertices(int part_id, const local_vid_t ids[], int num,
			compute_vertex *v_buf[]) const = 0;
	virtual compute_vertex &get_vertex(vertex_id_t id) = 0;
	virtual compute_vertex &get_vertex(int part_id, local_vid_t id) = 0;

	/*
	 * The interface of getting vertically partitioned vertices.
	 */
	virtual size_t get_num_vpart_vertices(int hpart_id) const = 0;
	virtual size_t get_vpart_vertex_pointers(int hpart_id,
			vpart_vertex_pointer ps[], int num) const = 0;
	virtual size_t get_vpart_vertices(int hpart_id, int vpart_id,
			vpart_vertex_pointer ps[], int num,
			compute_vertex_pointer vertices[]) const = 0;
	virtual size_t get_vpart_vertices(vertex_id_t id,
			compute_vertex_pointer vertices[], int num) const = 0;

	virtual vertex_id_t get_max_vertex_id() const = 0;

	virtual vertex_id_t get_min_vertex_id() const = 0;

	virtual size_t get_num_vertices() const = 0;

	virtual const graph_partitioner &get_partitioner() const = 0;

	virtual vertex_program::ptr create_def_vertex_program() const = 0;
	virtual vertex_program::ptr create_def_part_vertex_program() const = 0;

	virtual local_vid_t get_local_id(int part_id, const compute_vertex &v) const = 0;
	virtual vertex_id_t get_vertex_id(int part_id, const compute_vertex &v) const = 0;
	virtual vertex_id_t get_vertex_id(int part_id, compute_vertex_pointer v) const = 0;
	virtual vertex_id_t get_vertex_id(const compute_vertex &v) const = 0;
	virtual bool belong2part(const compute_vertex &v, int part_id) const = 0;
};

template<class vertex_type, class part_vertex_type>
class NUMA_graph_index;

template<class vertex_type, class part_vertex_type>
class graph_local_partition
{
	typedef std::pair<size_t, part_vertex_type *> part_vertex_array;

	struct vpart_vertex_comp
	{
		bool operator()(const part_vertex_type &v1,
				const vertex_id_t v2) const {
			return v1.get_id() < v2;
		}
	};

	int node_id;
	int part_id;
	size_t tot_num_vertices;
	size_t num_vertices;
	vertex_type *vertex_arr;
	std::vector<part_vertex_array> part_vertex_arrs;
	const graph_partitioner &partitioner;

	graph_local_partition(const graph_partitioner &_partitioner, int part_id,
			int node_id, size_t tot_num_vertices): partitioner(_partitioner) {
		this->part_id = part_id;
		this->tot_num_vertices = tot_num_vertices;
		this->num_vertices = partitioner.get_part_size(part_id,
				tot_num_vertices);
		this->node_id = node_id;
		vertex_arr = NULL;
	}
public:
	~graph_local_partition() {
		if (vertex_arr)
			free_large(vertex_arr, sizeof(vertex_arr[0]) * num_vertices);

		BOOST_FOREACH(part_vertex_array arr, part_vertex_arrs)
			free_large(arr.second, sizeof(arr.second[0]) * arr.first);
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

	/*
	 * This method is to initialize the vertical partitions inside
	 * the horizontal partition.
	 */
	void init_vparts(int num_parts, std::vector<vertex_id_t> &ids) {
		assert(std::is_sorted(ids.begin(), ids.end()));
		assert(num_parts > 1);
		part_vertex_arrs.resize(num_parts);
		for (int i = 0; i < num_parts; i++) {
			part_vertex_arrs[i].first = ids.size();
			part_vertex_arrs[i].second = (part_vertex_type *) malloc_large(
					sizeof(part_vertex_arrs[i].second[0]) * ids.size());
		}
		for (size_t i = 0; i < ids.size(); i++) {
			vertex_id_t id = ids[i];
			// horizontal part id.
			BOOST_VERIFY(this->part_id == partitioner.map(id));
			// vertical parts.
			for (int vpart_id = 0; vpart_id < num_parts; vpart_id++)
				new (part_vertex_arrs[vpart_id].second + i) part_vertex_type(id, vpart_id);
		}
	}

	local_vid_t get_local_id(const compute_vertex &v) const {
		vertex_type *addr1 = (vertex_type *) &v;
		// If the compute vertex is a main vertex in this partition.
		if (vertex_arr <= addr1 && addr1 < vertex_arr + num_vertices)
			return local_vid_t(addr1 - vertex_arr);
		else
			return local_vid_t(INVALID_VERTEX_ID);
	}

	vertex_id_t get_vertex_id(const compute_vertex &v) const {
		vertex_type *addr1 = (vertex_type *) &v;
		// If the compute vertex is a main vertex in this partition.
		if (vertex_arr <= addr1 && addr1 < vertex_arr + num_vertices) {
			vertex_id_t local_id = addr1 - vertex_arr;
			vertex_id_t id;
			partitioner.loc2map(part_id, local_id, id);
			return id;
		}
		else {
			// If not, is it one of the vertically partitioned vertices?
			part_vertex_type *addr2 = (part_vertex_type *) &v;
			for (size_t i = 0; i < part_vertex_arrs.size(); i++) {
				part_vertex_type *arr_start = part_vertex_arrs[i].second;
				part_vertex_type *arr_end
					= part_vertex_arrs[i].second + part_vertex_arrs[i].first;
				if (arr_start <= addr2 && addr2 < arr_end)
					return addr2->get_id();
			}
			// If not, the vertex doesn't belong to this partition.
			return INVALID_VERTEX_ID;
		}
	}

	compute_vertex &get_vertex(vertex_id_t id) {
		return vertex_arr[id];
	}

	size_t get_num_vpart_vertices() const {
		if (part_vertex_arrs.empty())
			return 0;
		else
			return part_vertex_arrs[0].first;
	}

	size_t get_vpart_vertex_pointers(vpart_vertex_pointer ps[],
			int num) const {
		int act_num = min(num, get_num_vpart_vertices());
		for (int i = 0; i < act_num; i++)
			ps[i] = vpart_vertex_pointer(part_vertex_arrs[0].second[i].get_id(), i);
		return act_num;
	}

	size_t get_vpart_vertices(int vpart_id, vpart_vertex_pointer ps[], int num,
			compute_vertex_pointer vertices[]) const {
		assert((size_t) vpart_id < part_vertex_arrs.size());
		int act_num = min(num, get_num_vpart_vertices());
		for (int i = 0; i < act_num; i++)
			vertices[i] = compute_vertex_pointer(
					&part_vertex_arrs[vpart_id].second[ps[i].get_off()], true);
		return act_num;
	}

	size_t get_vpart_vertices(vertex_id_t id,
			compute_vertex_pointer vertices[], int num) const {
		assert(!part_vertex_arrs.empty());
		// All vpart vertices are sorted on vertex Id.
		part_vertex_type *entry = std::lower_bound(part_vertex_arrs[0].second,
				part_vertex_arrs[0].second + part_vertex_arrs[0].first, id,
				vpart_vertex_comp());
		if (entry == part_vertex_arrs[0].second + part_vertex_arrs[0].first)
			return 0;
		off_t off = entry - part_vertex_arrs[0].second;
		int act_num = std::min((size_t) num, part_vertex_arrs.size());
		for (int i = 0; i < act_num; i++)
			vertices[i] = compute_vertex_pointer(
					part_vertex_arrs[i].second + off, true);
		return act_num;
	}

	size_t get_num_vertices() const {
		return num_vertices;
	}

	int get_node_id() const {
		return node_id;
	}

	friend class NUMA_graph_index<vertex_type, part_vertex_type>;
};

template<class vertex_type, class part_vertex_type>
class NUMA_graph_index: public graph_index
{
	graph_header header;
	vertex_id_t max_vertex_id;
	vertex_id_t min_vertex_id;
	std::unique_ptr<range_graph_partitioner> partitioner;
	// A graph index per thread
	std::vector<std::unique_ptr<graph_local_partition<vertex_type, part_vertex_type> > > index_arr;

	class init_thread: public thread
	{
		graph_local_partition<vertex_type, part_vertex_type> *index;
	public:
		init_thread(graph_local_partition<vertex_type, part_vertex_type> *index): thread(
				"index-init-thread", index->get_node_id()) {
			this->index = index;
		}

		virtual void run() {
			index->init();
			this->stop();
		}
	};
public:
	static graph_index::ptr create(const graph_header &header) {
		NUMA_graph_index<vertex_type, part_vertex_type> *index
			= new NUMA_graph_index<vertex_type, part_vertex_type>();
		index->header = header;
		return graph_index::ptr(index);
	}

	static graph_index::ptr create(const std::string &index_file) {
		NUMA_graph_index<vertex_type, part_vertex_type> *index
			= new NUMA_graph_index<vertex_type, part_vertex_type>();
		file_io_factory::shared_ptr factory = create_io_factory(index_file,
				GLOBAL_CACHE_ACCESS);
		io_interface::ptr io = factory->create_io(thread::get_curr_thread());
		io->access((char *) &index->header, 0, sizeof(header), READ);
		return graph_index::ptr(index);
	}

	void init(int num_threads, int num_nodes) {
		partitioner = std::unique_ptr<range_graph_partitioner>(
				new range_graph_partitioner(num_threads));

		// Construct the indices.
		for (int i = 0; i < num_threads; i++) {
			index_arr.emplace_back(new graph_local_partition<vertex_type, part_vertex_type>(
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

		BOOST_LOG_TRIVIAL(info) << boost::format("There are %1% vertices")
			% header.get_num_vertices();
	}

	virtual void init_vparts(int hpart_id, int num_vparts,
			std::vector<vertex_id_t> &ids) {
		index_arr[hpart_id]->init_vparts(num_vparts, ids);
	}

	virtual size_t get_vertices(const vertex_id_t ids[], int num,
			compute_vertex *v_buf[]) const {
		for (int i = 0; i < num; i++) {
			vertex_id_t id = ids[i];
			int part_id;
			off_t part_off;
			partitioner->map2loc(id, part_id, part_off);
			v_buf[i] = &index_arr[part_id]->get_vertex(part_off);
		}
		return num;
	}

	virtual compute_vertex &get_vertex(vertex_id_t id) {
		int part_id;
		off_t part_off;
		partitioner->map2loc(id, part_id, part_off);
		return index_arr[part_id]->get_vertex(part_off);
	}

	virtual size_t get_vertices(int part_id, const local_vid_t ids[], int num,
			compute_vertex *v_buf[]) const {
		for (int i = 0; i < num; i++)
			v_buf[i] = &index_arr[part_id]->get_vertex(ids[i].id);
		return num;
	}

	virtual compute_vertex &get_vertex(int part_id, local_vid_t id) {
		return index_arr[part_id]->get_vertex(id.id);
	}

	virtual size_t get_num_vpart_vertices(int hpart_id) const {
		return index_arr[hpart_id]->get_num_vpart_vertices();
	}

	virtual size_t get_vpart_vertex_pointers(int hpart_id,
			vpart_vertex_pointer ps[], int num) const {
		return index_arr[hpart_id]->get_vpart_vertex_pointers(ps, num);
	}

	virtual size_t get_vpart_vertices(int hpart_id, int vpart_id,
			vpart_vertex_pointer ps[], int num,
			compute_vertex_pointer vertices[]) const {
		return index_arr[hpart_id]->get_vpart_vertices(vpart_id, ps, num,
				vertices);
	}

	virtual size_t get_vpart_vertices(vertex_id_t id,
			compute_vertex_pointer vertices[], int num) const {
		int part_id = partitioner->map(id);
		return index_arr[part_id]->get_vpart_vertices(id, vertices, num);
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

	virtual vertex_program::ptr create_def_part_vertex_program(
			) const {
		return vertex_program::ptr(new vertex_program_impl<part_vertex_type>());
	}

	virtual vertex_id_t get_vertex_id(const compute_vertex &v) const {
		for (size_t i = 0; i < index_arr.size(); i++) {
			vertex_id_t id = index_arr[i]->get_vertex_id(v);
			if (id != INVALID_VERTEX_ID)
				return id;
		}
		ABORT_MSG("can't find vertex ID");
	}

	virtual vertex_id_t get_vertex_id(int part_id, compute_vertex_pointer v) const {
		if (v.is_part())
			return ((part_vertex_type &) *v).get_id();
		else {
			local_vid_t local_id = index_arr[part_id]->get_local_id(*v);
			vertex_id_t id;
			partitioner->loc2map(part_id, local_id.id, id);
			return id;
		}
	}

	virtual vertex_id_t get_vertex_id(int part_id, const compute_vertex &v) const {
		return index_arr[part_id]->get_vertex_id(v);
	}

	virtual local_vid_t get_local_id(int part_id, const compute_vertex &v) const {
		return index_arr[part_id]->get_local_id(v);
	}

	virtual bool belong2part(const compute_vertex &v, int part_id) const {
		// TODO there might be a more light-weight implementation.
		return get_vertex_id(part_id, v) != INVALID_VERTEX_ID;
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
