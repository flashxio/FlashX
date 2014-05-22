#ifndef __FG_VECTOR_H__
#define __FG_VECTOR_H__

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

#include <memory>

#include "graph_engine.h"
#include "stat.h"

template<class T>
class FG_vector
{
	// TODO I might need to split the vector into partitions.
	std::vector<T> eles;
	graph_engine::ptr graph;

	FG_vector(graph_engine::ptr graph) {
		this->graph = graph;
		eles.resize(graph->get_num_vertices());
	}

public:
	typedef typename std::shared_ptr<FG_vector<T> > ptr;

	static ptr create(graph_engine::ptr graph) {
		return ptr(new FG_vector<T>(graph));
	}

	void count_unique(count_map<T> &map) {
		// TODO we need a parallel implementation.
		BOOST_FOREACH(T v, eles) {
			map.add(v);
		}
	}

	size_t get_size() const {
		return eles.size();
	}

	// TODO these interfaces assume shared memory.

	void set(vertex_id_t id, const T &v) {
		eles[id] = v;
	}

	const T &get(vertex_id_t id) const {
		return eles[id];
	}
};

#endif
