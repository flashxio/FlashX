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

class vertex_compute;
class directed_vertex_compute;
class directed_vertex_request;

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

	virtual void request_vertices(vertex_id_t ids[], size_t num,
			vertex_compute &compute) = 0;
	virtual void request_num_edges(vertex_id_t vertices[], size_t num,
			vertex_compute &compute) = 0;
	virtual void request_num_directed_edges(vertex_id_t ids[], size_t num,
			directed_vertex_compute &compute) = 0;
	virtual void request_vertices(directed_vertex_request reqs[], size_t num,
			directed_vertex_compute &compute) = 0;
	virtual void wait4complete(int num) = 0;
	virtual size_t get_num_pending_tasks() const = 0;
};

#endif
