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

#ifdef USE_GZIP
#include <zlib.h>
#endif

#include <boost/format.hpp>
#include <boost/foreach.hpp>

#include "log.h"
#include "native_file.h"
#include "thread.h"

#include "data_io.h"
#include "generic_type.h"
#include "matrix_config.h"
#include "data_frame.h"
#include "mem_worker_thread.h"

namespace fm
{

static const int LINE_BLOCK_SIZE = 16 * 1024 * 1024;

class file_io
{
public:
	typedef std::shared_ptr<file_io> ptr;

	static ptr create(const std::string &file_name);

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

#ifdef USE_GZIP
class gz_file_io: public file_io
{
	std::vector<char> prev_buf;
	size_t prev_buf_bytes;

	gzFile f;

	gz_file_io(gzFile f) {
		this->f = f;
		prev_buf_bytes = 0;
		prev_buf.resize(PAGE_SIZE);
	}
public:
	static ptr create(const std::string &file);

	~gz_file_io() {
		gzclose(f);
	}

	std::unique_ptr<char[]> read_lines(const size_t wanted_bytes,
			size_t &read_bytes);

	bool eof() const {
		return gzeof(f) && prev_buf_bytes == 0;
	}
};

std::unique_ptr<char[]> gz_file_io::read_lines(
		const size_t wanted_bytes1, size_t &read_bytes)
{
	read_bytes = 0;
	size_t wanted_bytes = wanted_bytes1;
	size_t buf_size = wanted_bytes + PAGE_SIZE;
	char *buf = new char[buf_size];
	std::unique_ptr<char[]> ret_buf(buf);
	if (prev_buf_bytes > 0) {
		memcpy(buf, prev_buf.data(), prev_buf_bytes);
		buf += prev_buf_bytes;
		read_bytes += prev_buf_bytes;
		wanted_bytes -= prev_buf_bytes;
		prev_buf_bytes = 0;
	}

	if (!gzeof(f)) {
		int ret = gzread(f, buf, wanted_bytes + PAGE_SIZE);
		if (ret <= 0) {
			if (ret < 0 || !gzeof(f)) {
				BOOST_LOG_TRIVIAL(fatal) << gzerror(f, &ret);
				exit(1);
			}
		}

		if ((size_t) ret > wanted_bytes) {
			int i = 0;
			int over_read = ret - wanted_bytes;
			for (; i < over_read; i++) {
				if (buf[wanted_bytes + i] == '\n') {
					i++;
					break;
				}
			}
			read_bytes += wanted_bytes + i;
			buf += wanted_bytes + i;

			prev_buf_bytes = over_read - i;
			assert(prev_buf_bytes <= (size_t) PAGE_SIZE);
			memcpy(prev_buf.data(), buf, prev_buf_bytes);
		}
		else
			read_bytes += ret;
	}
	// The line buffer must end with '\0'.
	assert(read_bytes < buf_size);
	ret_buf[read_bytes] = 0;
	return ret_buf;
}

file_io::ptr gz_file_io::create(const std::string &file)
{
	gzFile f = gzopen(file.c_str(), "rb");
	if (f == Z_NULL) {
		BOOST_LOG_TRIVIAL(error)
			<< boost::format("fail to open gz file %1%: %2%")
			% file % strerror(errno);
		return ptr();
	}
	return ptr(new gz_file_io(f));
}

#endif

file_io::ptr file_io::create(const std::string &file_name)
{
#ifdef USE_GZIP
	size_t loc = file_name.rfind(".gz");
	// If the file name ends up with ".gz", we consider it as a gzip file.
	if (loc != std::string::npos && loc + 3 == file_name.length())
		return gz_file_io::create(file_name);
	else
#endif
		return text_file_io::create(file_name);
}

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
	std::atomic<size_t> num_dfs;
	std::vector<data_frame::ptr> dfs;
	pthread_mutex_t lock;
	pthread_cond_t fetch_cond;
	pthread_cond_t add_cond;
	size_t max_queue_size;
	bool wait_for_fetch;
	bool wait_for_add;
public:
	data_frame_set(size_t max_queue_size) {
		this->max_queue_size = max_queue_size;
		pthread_mutex_init(&lock, NULL);
		pthread_cond_init(&fetch_cond, NULL);
		pthread_cond_init(&add_cond, NULL);
		num_dfs = 0;
		wait_for_fetch = false;
		wait_for_add = false;
	}
	~data_frame_set() {
		pthread_mutex_destroy(&lock);
		pthread_cond_destroy(&fetch_cond);
		pthread_cond_destroy(&add_cond);
	}

	void add(data_frame::ptr df) {
		pthread_mutex_lock(&lock);
		while (dfs.size() >= max_queue_size) {
			// If the consumer thread is wait for adding new data frames, we
			// should wake them up before going to sleep. There is only
			// one thread waiting for fetching data frames, so we only need
			// to signal one thread.
			if (wait_for_fetch)
				pthread_cond_signal(&fetch_cond);
			wait_for_add = true;
			pthread_cond_wait(&add_cond, &lock);
			wait_for_add = false;
		}
		dfs.push_back(df);
		num_dfs++;
		pthread_mutex_unlock(&lock);
		pthread_cond_signal(&fetch_cond);
	}

	std::vector<data_frame::ptr> fetch_data_frames() {
		std::vector<data_frame::ptr> ret;
		pthread_mutex_lock(&lock);
		while (dfs.empty()) {
			// If some threads are wait for adding new data frames, we should
			// wake them up before going to sleep. Potentially, there are
			// multiple threads waiting at the same time, we should wake
			// all of them up.
			if (wait_for_add)
				pthread_cond_broadcast(&add_cond);
			wait_for_fetch = true;
			pthread_cond_wait(&fetch_cond, &lock);
			wait_for_fetch = false;
		}
		ret = dfs;
		dfs.clear();
		num_dfs = 0;
		pthread_mutex_unlock(&lock);
		pthread_cond_broadcast(&add_cond);
		return ret;
	}

	size_t get_num_dfs() const {
		return num_dfs;
	}
};

static data_frame::ptr create_data_frame(const line_parser &parser)
{
	data_frame::ptr df = data_frame::create();
	for (size_t i = 0; i < parser.get_num_cols(); i++)
		df->add_vec(parser.get_col_name(i),
				detail::smp_vec_store::create(0, parser.get_col_type(i)));
	return df;
}

static data_frame::ptr create_data_frame(const line_parser &parser, bool in_mem)
{
	data_frame::ptr df = data_frame::create();
	for (size_t i = 0; i < parser.get_num_cols(); i++)
		df->add_vec(parser.get_col_name(i),
				detail::vec_store::create(0, parser.get_col_type(i), -1, in_mem));
	return df;
}

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
		data_frame::ptr df = create_data_frame(parser);
		parse_lines(std::move(lines), size, parser, *df);
		dfs.add(df);
	}
};

class file_parse_task: public thread_task
{
	file_io::ptr io;
	const line_parser &parser;
	data_frame_set &dfs;
public:
	file_parse_task(file_io::ptr io, const line_parser &_parser,
			data_frame_set &_dfs): parser(_parser), dfs(_dfs) {
		this->io = io;
	}

	void run();
};

void file_parse_task::run()
{
	while (!io->eof()) {
		size_t size = 0;
		std::unique_ptr<char[]> lines = io->read_lines(LINE_BLOCK_SIZE, size);
		assert(size > 0);
		data_frame::ptr df = create_data_frame(parser);
		parse_lines(std::move(lines), size, parser, *df);
		dfs.add(df);
	}
}

}

data_frame::ptr read_lines(const std::string &file, const line_parser &parser,
		bool in_mem)
{
	data_frame::ptr df = create_data_frame(parser, in_mem);
	file_io::ptr io = file_io::create(file);
	if (io == NULL)
		return data_frame::ptr();

	printf("parse edge list\n");
	detail::mem_thread_pool::ptr mem_threads
		= detail::mem_thread_pool::get_global_mem_threads();
	const size_t MAX_PENDING = mem_threads->get_num_threads() * 3;
	data_frame_set dfs(MAX_PENDING);

	while (!io->eof()) {
		size_t num_tasks = MAX_PENDING - mem_threads->get_num_pending();
		for (size_t i = 0; i < num_tasks && !io->eof(); i++) {
			size_t size = 0;
			std::unique_ptr<char[]> lines = io->read_lines(LINE_BLOCK_SIZE, size);

			mem_threads->process_task(-1,
					new parse_task(std::move(lines), size, parser, dfs));
		}
		if (dfs.get_num_dfs() > 0) {
			std::vector<data_frame::ptr> tmp_dfs = dfs.fetch_data_frames();
			if (!tmp_dfs.empty())
				df->append(tmp_dfs.begin(), tmp_dfs.end());
		}
	}
	mem_threads->wait4complete();
	std::vector<data_frame::ptr> tmp_dfs = dfs.fetch_data_frames();
	if (!tmp_dfs.empty())
		df->append(tmp_dfs.begin(), tmp_dfs.end());

	return df;
}

data_frame::ptr read_lines(const std::vector<std::string> &files,
		const line_parser &parser, bool in_mem)
{
	if (files.size() == 1)
		return read_lines(files[0], parser, in_mem);

	data_frame::ptr df = data_frame::create();
	// TODO should I make these NUMA vectors?
	df->add_vec(parser.get_col_name(0),
			detail::vec_store::create(0, parser.get_col_type(0), -1, in_mem));
	df->add_vec(parser.get_col_name(1),
			detail::vec_store::create(0, parser.get_col_type(1), -1, in_mem));

	detail::mem_thread_pool::ptr mem_threads
		= detail::mem_thread_pool::get_global_mem_threads();
	const size_t MAX_PENDING = mem_threads->get_num_threads() * 3;
	data_frame_set dfs(MAX_PENDING);
	/*
	 * We assign a thread to each file. This works better if there are
	 * many small input files. If the input files are compressed, this
	 * approach also parallelizes decompression.
	 *
	 * TODO it may not work so well if there are a small number of large
	 * input files.
	 */
	auto file_it = files.begin();
	while (file_it != files.end()) {
		size_t num_tasks = MAX_PENDING - mem_threads->get_num_pending();
		for (size_t i = 0; i < num_tasks && file_it != files.end(); i++) {
			file_io::ptr io = file_io::create(*file_it);
			file_it++;
			mem_threads->process_task(-1,
					new file_parse_task(io, parser, dfs));
		}
		// This is the only thread that can fetch data frames from the queue.
		// If there are pending tasks in the thread pool, it's guaranteed
		// that we can fetch data frames from the queue.
		if (mem_threads->get_num_pending() > 0) {
			std::vector<data_frame::ptr> tmp_dfs = dfs.fetch_data_frames();
			if (!tmp_dfs.empty())
				df->append(tmp_dfs.begin(), tmp_dfs.end());
		}
	}
	// TODO It might be expensive to calculate the number of pending
	// tasks every time.
	while (mem_threads->get_num_pending() > 0) {
		std::vector<data_frame::ptr> tmp_dfs = dfs.fetch_data_frames();
		if (!tmp_dfs.empty())
			df->append(tmp_dfs.begin(), tmp_dfs.end());
	}
	mem_threads->wait4complete();
	// At this point, all threads have stoped working. If there are
	// data frames in the queue, they are the very last ones.
	if (dfs.get_num_dfs() > 0) {
		std::vector<data_frame::ptr> tmp_dfs = dfs.fetch_data_frames();
		if (!tmp_dfs.empty())
			df->append(tmp_dfs.begin(), tmp_dfs.end());
	}

	return df;
}

}
