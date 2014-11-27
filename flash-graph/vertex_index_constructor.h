#ifndef __VERTEX_INDEX_CONSTRUCTOR_H__
#define __VERTEX_INDEX_CONSTRUCTOR_H__

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
#include "vertex_index.h"

namespace fg
{

/*
 * These are the in-mem counterparts of vertex index above.
 * These in-memory data structures are used to construct vertex indices.
 */

class in_mem_vertex_index
{
public:
	typedef std::shared_ptr<in_mem_vertex_index> ptr;

	static ptr create(bool directed);
	static ptr create_compressed(bool directed, size_t edge_data_size);

	virtual ~in_mem_vertex_index() {
	}

	virtual void add_vertex(const in_mem_vertex &) = 0;
	virtual void dump(const std::string &file, const graph_header &header,
			bool compressed) = 0;
	virtual vertex_index::ptr dump(const graph_header &header, bool compressed) = 0;
};

}

#endif
