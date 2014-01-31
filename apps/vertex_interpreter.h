#ifndef __VERTEX_INTERPRETER_H__
#define __VERTEX_INTERPRETER_H__

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
			int size) const = 0;
	virtual page_vertex *interpret_part(const page_vertex *header,
			page_byte_array &, char *buf, int size) const = 0;
	/**
	 * The size of the vertex object.
	 */
	virtual int get_vertex_size() const = 0;
};

class ext_mem_directed_vertex_interpreter: public ext_mem_vertex_interpreter
{
public:
	virtual page_vertex *interpret(page_byte_array &array, char *buf,
			int size) const {
		assert(size >= (int) sizeof(page_directed_vertex));
		return new (buf) page_directed_vertex(array);
	}

	virtual page_vertex *interpret_part(const page_vertex *header,
			page_byte_array &, char *buf, int size) const {
		assert(0);
		return NULL;
	}

	virtual int get_vertex_size() const {
		return sizeof(page_directed_vertex);
	}
};

class ext_mem_undirected_vertex_interpreter: public ext_mem_vertex_interpreter
{
public:
	virtual page_vertex *interpret(page_byte_array &array, char *buf,
			int size) const {
		assert(size >= (int) sizeof(page_undirected_vertex));
		return new (buf) page_undirected_vertex(array);
	}

	virtual page_vertex *interpret_part(const page_vertex *header,
			page_byte_array &, char *buf, int size) const {
		assert(0);
		return NULL;
	}

	virtual int get_vertex_size() const {
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
			int size) const {
		assert(size >= TS_page_directed_vertex::get_size(num_timestamps));
		return TS_page_directed_vertex::create(array, buf, size);
	}

	virtual page_vertex *interpret_part(const page_vertex *header,
			page_byte_array &array, char *buf, int size) const {
		assert(size >= TS_page_directed_vertex::get_size(num_timestamps));
		return TS_page_directed_vertex::create(
				(const TS_page_directed_vertex *) header, array, buf, size);
	}

	virtual int get_vertex_size() const {
		return TS_page_directed_vertex::get_size(num_timestamps);
	}
};

#endif
