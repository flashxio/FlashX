#ifndef __VERTEX_INTERPRETER_H__
#define __VERTEX_INTERPRETER_H__

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

#include "cache.h"

#include "vertex.h"

/**
 * This defines the interface of interpreting vertices in the external memory.
 */
class ext_mem_vertex_interpreter
{
public:
	/**
	 * Interpret the data in the page byte array, and construct the page vertex
	 * in the buffer.
	 */
	virtual page_vertex *interpret(page_byte_array &, char *buf,
			size_t size) const = 0;
	virtual page_vertex *interpret_part(const page_vertex *header,
			page_byte_array &, char *buf, size_t size) const = 0;
	/**
	 * The size of the vertex object.
	 */
	virtual size_t get_vertex_size() const = 0;
};

class ext_mem_directed_vertex_interpreter: public ext_mem_vertex_interpreter
{
public:
	virtual page_vertex *interpret(page_byte_array &array, char *buf,
			size_t size) const {
		assert(size >= sizeof(page_directed_vertex));
		return new (buf) page_directed_vertex(array);
	}

	virtual page_vertex *interpret_part(const page_vertex *header,
			page_byte_array &, char *buf, size_t size) const {
		assert(0);
		return NULL;
	}

	virtual size_t get_vertex_size() const {
		return sizeof(page_directed_vertex);
	}
};

class ext_mem_undirected_vertex_interpreter: public ext_mem_vertex_interpreter
{
public:
	virtual page_vertex *interpret(page_byte_array &array, char *buf,
			size_t size) const {
		assert(size >= sizeof(page_undirected_vertex));
		return new (buf) page_undirected_vertex(array);
	}

	virtual page_vertex *interpret_part(const page_vertex *header,
			page_byte_array &, char *buf, size_t size) const {
		assert(0);
		return NULL;
	}

	virtual size_t get_vertex_size() const {
		return sizeof(page_undirected_vertex);
	}
};

class ts_ext_mem_vertex_interpreter: public ext_mem_vertex_interpreter
{
	int num_timestamps;
public:
	ts_ext_mem_vertex_interpreter(int num_timestamps) {
		this->num_timestamps = num_timestamps;
	}

	virtual page_vertex *interpret(page_byte_array &array, char *buf,
			size_t size) const {
		assert(size >= TS_page_directed_vertex::get_size(num_timestamps));
		return TS_page_directed_vertex::create(array, buf, size);
	}

	virtual page_vertex *interpret_part(const page_vertex *header,
			page_byte_array &array, char *buf, size_t size) const {
		assert(size >= TS_page_directed_vertex::get_size(num_timestamps));
		return TS_page_directed_vertex::create(
				(const TS_page_directed_vertex *) header, array, buf, size);
	}

	virtual size_t get_vertex_size() const {
		return TS_page_directed_vertex::get_size(num_timestamps);
	}
};

#endif
