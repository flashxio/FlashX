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

#include "thread.h"
#include "io_interface.h"

#include "graph_engine.h"
#include "graph_config.h"
#include "graphlab/cuckoo_set_pow2.hpp"
#include "FG_vector.h"
#include "FGlib.h"

/**
 * This contains the data structures shared directed triangle counting
 * and undirected triangle counting.
 */

const double BIN_SEARCH_RATIO = 100;
const int HASH_SEARCH_RATIO = 16;
const int hash_threshold = 1000;

static atomic_number<long> num_working_vertices;
static atomic_number<long> num_completed_vertices;

class count_msg: public vertex_message
{
	int num;
public:
	count_msg(int num): vertex_message(sizeof(count_msg), false) {
		this->num = num;
	}

	int get_num() const {
		return num;
	}
};

struct runtime_data_t
{
	class index_entry
	{
		vertex_id_t id;
		uint32_t idx;
	public:
		index_entry() {
			id = -1;
			idx = -1;
		}

		index_entry(vertex_id_t id) {
			this->id = id;
			this->idx = -1;
		}

		index_entry(vertex_id_t id, uint32_t idx) {
			this->id = id;
			this->idx = idx;
		}

		vertex_id_t get_id() const {
			return id;
		}

		uint32_t get_idx() const {
			return idx;
		}

		bool operator==(const index_entry &e) const {
			return id == e.get_id();
		}
	};

	class index_hash
	{
		boost::hash<vertex_id_t> id_hash;
	public:
		size_t operator()(const index_entry &e) const {
			return id_hash(e.get_id());
		}
	};
	typedef graphlab::cuckoo_set_pow2<index_entry, 3, size_t,
			index_hash> edge_set_t;
	// It contains part of the edge list.
	// We only use the neighbors whose ID is smaller than this vertex.
	std::vector<vertex_id_t> edges;
	// The vector contains the number of the neighbors' triangles shared
	// with this vertex. It only keeps the triangles of neighbors in the
	// in-edges.
	std::vector<int> triangles;
	// The number of vertices that have joined with the vertex.
	size_t num_joined;
	size_t num_required;
	size_t num_triangles;

	edge_set_t edge_set;
public:
	runtime_data_t(size_t num_edges, size_t num_triangles): edge_set(
				index_entry(), 0, 2 * num_edges) {
		num_joined = 0;
		this->num_required = 0;
		this->num_triangles = num_triangles;
	}

	void finalize_init() {
		// We only build a hash table on large vertices
		if (edges.size() > (size_t) hash_threshold)
			for (size_t i = 0; i < edges.size(); i++)
				edge_set.insert(index_entry(edges[i], i));
		triangles.resize(edges.size());
	}
};

enum multi_func_flags
{
	NUM_TRIANGLES,
	POINTER,
	NUM_FLAGS,
};

class multi_func_value
{
	static const int VALUE_BITS = sizeof(size_t) * 8 - NUM_FLAGS;
	static const size_t FLAGS_MASK = ((1UL << VALUE_BITS) - 1);
	size_t value;

	void set_flag(int flag) {
		value |= 1UL << (VALUE_BITS + flag);
	}

	bool has_flag(int flag) const {
		return value & (1UL << (VALUE_BITS + flag));
	}
public:
	multi_func_value() {
		value = 0;
		// By default, it stores the number of triangles.
		set_flag(NUM_TRIANGLES);
	}

	void set_num_triangles(size_t num) {
		value = num;
		set_flag(NUM_TRIANGLES);
	}

	bool has_num_triangles() const {
		return has_flag(NUM_TRIANGLES);
	}

	size_t get_num_triangles() const {
		assert(has_flag(NUM_TRIANGLES));
		return value & FLAGS_MASK;
	}

	void inc_num_triangles(size_t num) {
		assert(has_flag(NUM_TRIANGLES));
		value += num;
	}

	/**
	 * Pointer to the runtime data.
	 */

	void set_runtime_data(runtime_data_t *data) {
		value = (size_t) data;
		set_flag(POINTER);
	}

	bool has_runtime_data() const {
		return has_flag(POINTER);
	}

	runtime_data_t *get_runtime_data() const {
		assert(has_flag(POINTER));
		return (runtime_data_t *) (value & FLAGS_MASK);
	}
};
