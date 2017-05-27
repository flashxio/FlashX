#ifndef __DATA_IO_H__
#define __DATA_IO_H__

/*
 * Copyright 2014 Open Connectome Project (http://openconnecto.me)
 * Written by Da Zheng (zhengda1936@gmail.com)
 *
 * This file is part of FlashMatrix.
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

/*
 * This file contains the code that imports data from other sources and
 * exports data.
 */

#include <stdlib.h>

#include <memory>
#include <vector>
#include <limits>

#include "generic_type.h"

namespace fm
{

class data_frame;
class scalar_type;
class dense_matrix;

/*
 * This class provides a common interface to read data from I/O.
 * Data might be in a local file, in the network or compressed.
 * The data is read from an I/O stream. Once data is read from
 * the stream, we can't read it again. If the data is going to be
 * reused, we need to buffer it.
 */
class file_io
{
public:
	static std::shared_ptr<char> alloc_io_buf(size_t num_bytes);

	typedef std::shared_ptr<file_io> ptr;

	static ptr create_local(const std::string name);
	static ptr create_gz(const std::string name);

	virtual std::shared_ptr<char> read_bytes(size_t wanted_bytes,
			size_t &read_bytes) = 0;
	virtual bool eof() const = 0;

	virtual std::string get_name() const = 0;
};

/*
 * This class is to read text lines from streaming I/O.
 * It allows users to peek data in the stream. Users can peek
 * the same data as many times as they want.
 */
class text_io
{
	// This points to some location that contains data to be read.
	const char *data_ptr;
	// The number of bytes for the data to be read.
	size_t data_size;
	// The buffer contains some data from the I/O stream.
	std::shared_ptr<char> data_buf;
	file_io::ptr io;
public:
	typedef std::shared_ptr<text_io> ptr;

	static ptr create(const std::string &file_name);
	static ptr create(file_io::ptr io) {
		return ptr(new text_io(io));
	}

	text_io() {
		data_ptr = NULL;
		data_size = 0;
	}

	text_io(file_io::ptr io) {
		data_ptr = NULL;
		data_size = 0;
		this->io = io;
	}

	virtual ~text_io() {
	}

	/*
	 * Read some bytes and the current location to read will be moved
	 * accordingly.
	 */
	virtual std::shared_ptr<char> read_lines(size_t wanted_bytes,
			size_t &read_bytes);
	/*
	 * Read some bytes and the current location to read doesn't move.
	 */
	virtual std::shared_ptr<char> peek(size_t wanted_bytes,
			size_t &read_bytes);

	/*
	 * Test if we have reached the end of the file.
	 */
	virtual bool eof() const {
		if (data_buf)
			return false;
		else
			return io->eof();
	}
};

class line_parser
{
public:
	virtual size_t parse(const std::vector<std::string> &lines,
			data_frame &df) const = 0;
	virtual size_t get_num_cols() const = 0;
	virtual const scalar_type &get_col_type(off_t idx) const = 0;
	virtual std::string get_col_name(off_t idx) const = 0;
};

/*
 * This converts a string to an element.
 */
class ele_parser
{
public:
	typedef std::shared_ptr<const ele_parser> const_ptr;

	virtual void set_zero(void *buf) const = 0;
	virtual void parse(const std::string &str, void *buf) const = 0;
	virtual const scalar_type &get_type() const = 0;
};

/*
 * Convert a string of decimal to an integer.
 */
template<class T>
class int_parser: public ele_parser
{
	int base;
public:
	int_parser() {
		base = 10;
	}

	int_parser(int base) {
		this->base = base;
	}

	virtual void set_zero(void *buf) const {
		T *val = reinterpret_cast<T *>(buf);
		*val = 0;
	}

	virtual void parse(const std::string &str, void *buf) const {
		T *val = reinterpret_cast<T *>(buf);
		val[0] = strtol(str.c_str(), NULL, base);
	}
	virtual const scalar_type &get_type() const {
		return get_scalar_type<T>();
	}
};

template<class T>
class float_parser: public ele_parser
{
public:
	virtual void parse(const std::string &str, void *buf) const {
		T *val = reinterpret_cast<T *>(buf);
		val[0] = atof(str.c_str());
	}

	virtual void set_zero(void *buf) const {
		T *val = reinterpret_cast<T *>(buf);
		*val = 0;
	}
	virtual const scalar_type &get_type() const {
		return get_scalar_type<T>();
	}
};

ele_parser::const_ptr get_ele_parser(const std::string &type);
const scalar_type &get_ele_type(const std::string &type);
bool valid_ele_type(const std::string &type);

std::shared_ptr<data_frame> read_data_frame(const std::vector<std::string> &files,
		bool in_mem, bool sequential, const std::string &delim,
		const std::vector<ele_parser::const_ptr> &parsers,
		const std::vector<off_t> &dup_col_idxs = std::vector<off_t>());
std::shared_ptr<dense_matrix> read_matrix(const std::vector<std::string> &files,
		bool in_mem, bool sequential, const std::string &ele_type,
		const std::string &delim,
		size_t num_cols = std::numeric_limits<size_t>::max());
std::shared_ptr<dense_matrix> read_matrix(const std::vector<std::string> &files,
		bool in_mem, bool sequential, const std::string &ele_type,
		const std::string &delim, const std::string &col_indicator);

std::shared_ptr<data_frame> read_data_frame(
		const std::vector<text_io::ptr> &file_ios,
		bool in_mem, bool sequential, const std::string &delim,
		const std::vector<ele_parser::const_ptr> &parsers,
		const std::vector<off_t> &dup_col_idxs = std::vector<off_t>());
std::shared_ptr<dense_matrix> read_matrix(
		const std::vector<text_io::ptr> &file_ios,
		bool in_mem, bool sequential, const std::string &ele_type,
		const std::string &delim,
		size_t num_cols = std::numeric_limits<size_t>::max());
std::shared_ptr<dense_matrix> read_matrix(
		const std::vector<text_io::ptr> &file_ios,
		bool in_mem, bool sequential, const std::string &ele_type,
		const std::string &delim, const std::string &col_indicator);

}

#endif
