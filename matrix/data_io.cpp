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

#include <boost/format.hpp>
#include <boost/foreach.hpp>

#include "log.h"
#include "native_file.h"
#include "thread.h"

#include "data_io.h"
#include "generic_type.h"
#include "matrix_config.h"
#include "data_frame.h"
#include "mem_vector.h"
#include "mem_worker_thread.h"

namespace fm
{

static const int LINE_BLOCK_SIZE = 16 * 1024 * 1024;

class file_io
{
public:
	typedef std::shared_ptr<file_io> ptr;

	virtual ~file_io() {
	}

	virtual std::unique_ptr<char[]> read_lines(size_t wanted_bytes,
			size_t &read_bytes) = 0;

	virtual bool eof() const = 0;
};

class text_file_io: public file_io
{
	FILE *f;
	ssize_t file_size;

	text_file_io(FILE *f, const std::string file) {
		this->f = f;
		safs::native_file local_f(file);
		file_size = local_f.get_size();
	}
public:
	static ptr create(const std::string file);

	~text_file_io() {
		if (f)
			fclose(f);
	}

	std::unique_ptr<char[]> read_lines(size_t wanted_bytes,
			size_t &read_bytes);

	bool eof() const {
		off_t curr_off = ftell(f);
		return file_size - curr_off == 0;
	}
};

file_io::ptr text_file_io::create(const std::string file)
{
	FILE *f = fopen(file.c_str(), "r");
	if (f == NULL) {
		BOOST_LOG_TRIVIAL(error)
			<< boost::format("fail to open %1%: %2%") % file % strerror(errno);
		return ptr();
	}
	return ptr(new text_file_io(f, file));
}

std::unique_ptr<char[]> text_file_io::read_lines(
		size_t wanted_bytes, size_t &read_bytes)
{
	off_t curr_off = ftell(f);
	off_t off = curr_off + wanted_bytes;
	// After we just to the new location, we need to further read another
	// page to search for the end of a line. If there isn't enough data,
	// we can just read all remaining data.
	if (off + PAGE_SIZE < file_size) {
		int ret = fseek(f, off, SEEK_SET);
		if (ret < 0) {
			perror("fseek");
			return NULL;
		}

		char buf[PAGE_SIZE];
		ret = fread(buf, sizeof(buf), 1, f);
		if (ret != 1) {
			perror("fread");
			return NULL;
		}
		unsigned i;
		for (i = 0; i < sizeof(buf); i++)
			if (buf[i] == '\n')
				break;
		// A line shouldn't be longer than a page.
		assert(i != sizeof(buf));

		// We read a little more than asked to make sure that we read
		// the entire line.
		read_bytes = wanted_bytes + i + 1;

		// Go back to the original offset in the file.
		ret = fseek(f, curr_off, SEEK_SET);
		assert(ret == 0);
	}
	else {
		read_bytes = file_size - curr_off;
	}

	// The line buffer must end with '\0'.
	char *line_buf = new char[read_bytes + 1];
	BOOST_VERIFY(fread(line_buf, read_bytes, 1, f) == 1);
	line_buf[read_bytes] = 0;

	return std::unique_ptr<char[]>(line_buf);
}

/*
 * Parse the lines in the character buffer.
 * `size' doesn't include '\0'.
 */
static size_t parse_lines(std::unique_ptr<char[]> line_buf, size_t size,
		const line_parser &parser, data_frame &df)
{
	char *line_end;
	char *line = line_buf.get();
	char *buf_start = line_buf.get();
	// TODO I need to be careful here. It potentially needs to allocate
	// a lot of small pieces of memory.
	std::vector<std::string> lines;
	while ((line_end = strchr(line, '\n'))) {
		assert(line_end - buf_start <= (ssize_t) size);
		*line_end = 0;
		if (*(line_end - 1) == '\r')
			*(line_end - 1) = 0;
		lines.push_back(std::string(line));
		line = line_end + 1;
	}
	if (line - buf_start < (ssize_t) size)
		lines.push_back(std::string(line));
	
	return parser.parse(lines, df);
}

namespace
{

class data_frame_set
{
	std::vector<data_frame::ptr> dfs;
	pthread_spinlock_t lock;
public:
	data_frame_set() {
		pthread_spin_init(&lock, PTHREAD_PROCESS_PRIVATE);
	}

	void add(data_frame::ptr df) {
		pthread_spin_lock(&lock);
		dfs.push_back(df);
		pthread_spin_unlock(&lock);
	}

	const std::vector<data_frame::ptr> &get_data_frames() const {
		return dfs;
	}
};

class parse_task: public thread_task
{
	std::unique_ptr<char[]> lines;
	size_t size;
	const line_parser &parser;
	data_frame_set &dfs;
public:
	parse_task(std::unique_ptr<char[]> _lines, size_t size,
			const line_parser &_parser, data_frame_set &_dfs): parser(
				_parser), dfs(_dfs) {
		this->lines = std::move(_lines);
		this->size = size;
	}

	void run() {
		data_frame::ptr df = data_frame::create();
		df->add_vec(parser.get_col_name(0),
				detail::smp_vec_store::create(0, parser.get_col_type(0)));
		df->add_vec(parser.get_col_name(1),
				detail::smp_vec_store::create(0, parser.get_col_type(1)));
		parse_lines(std::move(lines), size, parser, *df);
		dfs.add(df);
	}
};

}

static bool read_lines(const std::string &file, const line_parser &parser,
		data_frame &df)
{
	file_io::ptr io = text_file_io::create(file);
	if (io == NULL)
		return false;

	printf("parse edge list\n");
	detail::mem_thread_pool::ptr mem_threads
		= detail::mem_thread_pool::get_global_mem_threads();
	data_frame_set dfs;
	size_t i = 0;
	while (!io->eof()) {
		size_t size = 0;
		std::unique_ptr<char[]> lines = io->read_lines(LINE_BLOCK_SIZE, size);

		int node_id = i % mem_threads->get_num_nodes();
		mem_threads->process_task(node_id,
				new parse_task(std::move(lines), size, parser, dfs));
		i++;
	}
	mem_threads->wait4complete();

	printf("merge parse results\n");
	df.append(dfs.get_data_frames().begin(), dfs.get_data_frames().end());
	return true;
}

data_frame::ptr read_lines(const std::vector<std::string> &files,
		const line_parser &parser)
{
	data_frame::ptr df = data_frame::create();
	df->add_vec(parser.get_col_name(0),
			detail::smp_vec_store::create(0, parser.get_col_type(0)));
	df->add_vec(parser.get_col_name(1),
			detail::smp_vec_store::create(0, parser.get_col_type(1)));
	BOOST_FOREACH(std::string file, files) {
		bool ret = read_lines(file, parser, *df);
		if (!ret)
			return data_frame::ptr();
	}
	return df;
}

}
