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
#include <boost/algorithm/string.hpp>

#include "log.h"
#include "native_file.h"
#include "thread.h"

#include "data_io.h"
#include "generic_type.h"
#include "matrix_config.h"
#include "data_frame.h"
#include "mem_worker_thread.h"
#include "dense_matrix.h"

namespace fm
{

static const int LINE_BLOCK_SIZE = 16 * 1024 * 1024;

namespace
{

struct del_arr {
	void operator()(char *buf) {
		delete [] buf;
	}
};

}

std::shared_ptr<char> file_io::alloc_io_buf(size_t num_bytes)
{
	return std::shared_ptr<char>(new char[num_bytes], del_arr());
}

static size_t get_last_newline(std::shared_ptr<char> data, size_t size)
{
	off_t off = size - 1;
	const char *buf = data.get();
	for (; buf[off] != '\n' && off > 0; off--);
	if (off == 0 && buf[0] != '\n')
		return size;
	else
		return off + 1;
};

std::shared_ptr<char> text_io::read_lines(size_t wanted_bytes,
		size_t &read_bytes)
{
	if (data_buf && wanted_bytes <= data_size) {
		std::shared_ptr<char> new_data = file_io::alloc_io_buf(wanted_bytes + 1);
		memcpy(new_data.get(), data_ptr, wanted_bytes);
		read_bytes = get_last_newline(new_data, wanted_bytes);
		new_data.get()[read_bytes] = 0;
		// Update the data buffer.
		data_ptr += read_bytes;
		data_size -= read_bytes;
		if (data_size == 0) {
			data_ptr = NULL;
			data_buf = NULL;
		}
		return new_data;
	}
	else if (data_buf && wanted_bytes > data_size) {
		// Read the remaining data.
		size_t new_read_size = 0;
		std::shared_ptr<char> remain = io->read_bytes(wanted_bytes - data_size,
				new_read_size);

		// Copy data to be returned.
		size_t tot_data_size = new_read_size + data_size;
		std::shared_ptr<char> new_data = file_io::alloc_io_buf(tot_data_size + 1);
		memcpy(new_data.get(), data_ptr, data_size);
		if (remain && new_read_size > 0)
			memcpy(new_data.get() + data_size, remain.get(), new_read_size);
		read_bytes = get_last_newline(new_data, tot_data_size);

		// Update the data buffer.
		if (tot_data_size - read_bytes > 0) {
			data_size = tot_data_size - read_bytes;
			data_buf = file_io::alloc_io_buf(data_size);
			memcpy(data_buf.get(), new_data.get() + read_bytes, data_size);
			data_ptr = data_buf.get();
		}
		else {
			data_size = 0;
			data_ptr = NULL;
			data_buf = NULL;
		}

		// We should change the data after we save the remaining data.
		new_data.get()[read_bytes] = 0;

		return new_data;
	}
	else {
		size_t tot_data_size = 0;
		auto new_data = io->read_bytes(wanted_bytes, tot_data_size);
		read_bytes = get_last_newline(new_data, tot_data_size);

		// Update the data buffer.
		if (tot_data_size - read_bytes > 0) {
			data_size = tot_data_size - read_bytes;
			data_buf = file_io::alloc_io_buf(data_size);
			memcpy(data_buf.get(), new_data.get() + read_bytes, data_size);
			data_ptr = data_buf.get();
		}

		// We should change the data after we save the remaining data.
		new_data.get()[read_bytes] = 0;
		return new_data;
	}
}

std::shared_ptr<char> text_io::peek(size_t wanted_bytes,
		size_t &read_bytes)
{
	if (data_buf && wanted_bytes <= data_size) {
		std::shared_ptr<char> new_data = file_io::alloc_io_buf(wanted_bytes + 1);
		memcpy(new_data.get(), data_ptr, wanted_bytes);
		read_bytes = get_last_newline(new_data, wanted_bytes);
		new_data.get()[read_bytes] = 0;
		return new_data;
	}
	else if (data_buf && wanted_bytes > data_size) {
		// Read the remaining data.
		size_t new_read_size = 0;
		std::shared_ptr<char> remain = io->read_bytes(wanted_bytes - data_size,
				new_read_size);

		// Copy data to be returned.
		size_t tot_data_size = new_read_size + data_size;
		std::shared_ptr<char> new_data = file_io::alloc_io_buf(tot_data_size + 1);
		memcpy(new_data.get(), data_ptr, data_size);
		if (remain && new_read_size > 0)
			memcpy(new_data.get() + data_size, remain.get(), new_read_size);

		// Update the data buffer.
		if (new_read_size > 0) {
			data_buf = file_io::alloc_io_buf(tot_data_size);
			memcpy(data_buf.get(), new_data.get(), tot_data_size);
			data_size = tot_data_size;
			data_ptr = data_buf.get();
		}

		read_bytes = get_last_newline(new_data, tot_data_size);
		new_data.get()[read_bytes] = 0;

		return new_data;
	}
	else {
		size_t tot_read_bytes = 0;
		data_buf = io->read_bytes(wanted_bytes, tot_read_bytes);
		if (tot_read_bytes == 0) {
			read_bytes = 0;
			return std::shared_ptr<char>();
		}
		read_bytes = get_last_newline(data_buf, tot_read_bytes);
		std::shared_ptr<char> new_data = file_io::alloc_io_buf(read_bytes + 1);
		memcpy(new_data.get(), data_buf.get(), read_bytes);
		new_data.get()[read_bytes] = 0;
		// Update the data buffer.
		data_size = tot_read_bytes;
		data_ptr = data_buf.get();
		return new_data;
	}
}

class local_file_io: public file_io
{
	struct del_off_ptr {
		void *orig_addr;

		del_off_ptr(void *addr) {
			this->orig_addr = addr;
		}

		void operator()(char *buf) {
			free(orig_addr);
		}
	};

	int fd;
	off_t curr_off;
	ssize_t file_size;
	std::string name;

	local_file_io(int fd, const std::string file) {
		this->curr_off = 0;
		this->fd = fd;
		safs::native_file local_f(file);
		file_size = local_f.get_size();
		this->name = file;
	}
public:
	static ptr create(const std::string file);

	~local_file_io() {
		close(fd);
	}

	std::shared_ptr<char> read_bytes(size_t wanted_bytes,
			size_t &read_bytes);

	bool eof() const {
		return file_size - curr_off == 0;
	}

	std::string get_name() const {
		return name;
	}
};

#ifdef USE_GZIP

class gz_file_io: public file_io
{
	gzFile f;
	std::string name;

	gz_file_io(gzFile f, std::string name) {
		this->f = f;
		this->name = name;
	}
public:
	static ptr create(const std::string &file);

	~gz_file_io() {
		gzclose(f);
	}

	std::shared_ptr<char> read_bytes(const size_t wanted_bytes,
			size_t &read_bytes);

	bool eof() const {
		return gzeof(f);
	}

	std::string get_name() const {
		return name;
	}
};

std::shared_ptr<char> gz_file_io::read_bytes(
		const size_t wanted_bytes, size_t &read_bytes)
{
	read_bytes = 0;
	std::shared_ptr<char> ret_buf(new char[wanted_bytes], del_arr());
	if (!gzeof(f)) {
		int ret = gzread(f, ret_buf.get(), wanted_bytes);
		if (ret <= 0) {
			if (ret < 0 || !gzeof(f)) {
				BOOST_LOG_TRIVIAL(fatal) << gzerror(f, &ret);
				return std::unique_ptr<char[]>();
			}
		}
		read_bytes += ret;
	}
	assert(read_bytes <= wanted_bytes);
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
	return ptr(new gz_file_io(f, file));
}

#endif

file_io::ptr file_io::create_local(const std::string name)
{
	return local_file_io::create(name);
}

file_io::ptr file_io::create_gz(const std::string name)
{
#ifdef USE_GZIP
	return gz_file_io::create(name);
#else
	return file_io::ptr();
#endif
}

text_io::ptr text_io::create(const std::string &file_name)
{
	file_io::ptr io;
#ifdef USE_GZIP
	size_t loc = file_name.rfind(".gz");
	// If the file name ends up with ".gz", we consider it as a gzip file.
	if (loc != std::string::npos && loc + 3 == file_name.length())
		io = gz_file_io::create(file_name);
	else
#endif
		io = local_file_io::create(file_name);
	if (io == NULL)
		return text_io::ptr();
	else
		return text_io::ptr(new text_io(io));
}

file_io::ptr local_file_io::create(const std::string file)
{
	int fd = open(file.c_str(), O_RDONLY | O_DIRECT);
	if (fd < 0)
		fd = open(file.c_str(), O_RDONLY);
	if (fd < 0) {
		BOOST_LOG_TRIVIAL(error)
			<< boost::format("fail to open %1%: %2%") % file % strerror(errno);
		return ptr();
	}
	return ptr(new local_file_io(fd, file));
}

void read_complete(int fd, char *buf, size_t buf_size, size_t expected_size)
{
	assert(buf_size >= expected_size);
	while (expected_size > 0) {
		ssize_t ret = read(fd, buf, buf_size);
		assert(ret >= 0);
		buf += ret;
		buf_size -= ret;
		expected_size -= ret;
	}
}

std::shared_ptr<char> local_file_io::read_bytes(
		size_t wanted_bytes, size_t &read_bytes)
{
	off_t align_start = ROUND_PAGE(curr_off);
	off_t align_end = ROUNDUP_PAGE(curr_off + wanted_bytes);
	off_t local_off = curr_off - align_start;
	off_t seek_ret = lseek(fd, align_start, SEEK_SET);
	assert(seek_ret >= 0);

	size_t buf_size = align_end - align_start;
	void *addr = NULL;
	int alloc_ret = posix_memalign(&addr, PAGE_SIZE, buf_size);
	assert(alloc_ret == 0);

	assert(file_size > align_start);
	size_t expected_size = std::min(buf_size, (size_t) (file_size - align_start));
	read_complete(fd, (char *) addr, buf_size, expected_size);

	char *buf = ((char *) addr) + local_off;
	read_bytes = std::min(wanted_bytes, (size_t) (file_size - curr_off));
	curr_off += read_bytes;
	assert(curr_off <= file_size);
	return std::shared_ptr<char>(buf, del_off_ptr(addr));
}

/*
 * Parse the lines in the character buffer.
 * `size' doesn't include '\0'.
 */
static size_t parse_lines(std::shared_ptr<char> line_buf, size_t size,
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
	std::atomic<size_t> tot_num_dfs;
	std::atomic<size_t> num_dfs;
	// The key should be a sequence number to order the data frames.
	std::map<off_t, data_frame::ptr> dfs;
	pthread_mutex_t lock;
	pthread_cond_t fetch_cond;
	pthread_cond_t add_cond;
	size_t max_queue_size;
	// This is a producer-consumer model. There are multiple producers
	// but only one consumer. When the queue is full and the consumer
	// isn't waiting (wait_for_fetch=false), a producer needs to go to
	// sleep and increase wait_for_add by 1. When wait_for_add isn't 0,
	// it means there is at least one consumer is sleeping.
	size_t wait_for_add;
	bool wait_for_fetch;
public:
	data_frame_set(size_t max_queue_size) {
		this->max_queue_size = max_queue_size;
		pthread_mutex_init(&lock, NULL);
		pthread_cond_init(&fetch_cond, NULL);
		pthread_cond_init(&add_cond, NULL);
		tot_num_dfs = 0;
		num_dfs = 0;
		wait_for_fetch = false;
		wait_for_add = 0;
	}
	~data_frame_set() {
		pthread_mutex_destroy(&lock);
		pthread_cond_destroy(&fetch_cond);
		pthread_cond_destroy(&add_cond);
	}

	void add(off_t seq_id, data_frame::ptr df);

	/*
	 * If a user wants to fetch the data frame sequentially, we only fetch
	 * the ones with the expected sequence numbers.
	 */
	std::vector<data_frame::ptr> fetch_data_frames(bool sequential,
			off_t seq_start);

	size_t get_num_dfs() const {
		return num_dfs;
	}
};

void data_frame_set::add(off_t seq_id, data_frame::ptr df)
{
	pthread_mutex_lock(&lock);
	// If the consumer thread isn't sleeping and the queue is full,
	// the producer thread needs to sleep.
	while (!wait_for_fetch && dfs.size() >= max_queue_size) {
		wait_for_add++;
		pthread_cond_wait(&add_cond, &lock);
		wait_for_add--;
	}
	// If the consumer is sleeping, the producer adds the new data
	// to the queue regardless the status of the queue.
	if (seq_id < 0)
		seq_id = tot_num_dfs;
	dfs.insert(std::pair<off_t, data_frame::ptr>(seq_id, df));
	tot_num_dfs++;
	num_dfs++;
	// If the consumer thread is sleeping, we need to wake it up.
	pthread_mutex_unlock(&lock);
	pthread_cond_signal(&fetch_cond);
}

std::vector<data_frame::ptr> data_frame_set::fetch_data_frames(bool sequential,
		off_t seq_start)
{
	std::vector<data_frame::ptr> ret;
	pthread_mutex_lock(&lock);
	while (dfs.empty()) {
		// If some threads are wait for adding new data frames, we should
		// wake them up before going to sleep. Potentially, there are
		// multiple threads waiting at the same time, we should wake
		// all of them up.
		if (wait_for_add > 0)
			pthread_cond_broadcast(&add_cond);
		wait_for_fetch = true;
		pthread_cond_wait(&fetch_cond, &lock);
		wait_for_fetch = false;
	}
	// If we want to fetch the data frames sequentially according
	// to the sequence number.
	if (sequential) {
		auto it = dfs.begin();
		// We first need to make sure the first sequence number in
		// the map is the one we expect. If not, we have to wait.
		while (it->first != seq_start) {
			wait_for_fetch = true;
			if (wait_for_add > 0)
				pthread_cond_broadcast(&add_cond);
			pthread_cond_wait(&fetch_cond, &lock);
			wait_for_fetch = false;
			it = dfs.begin();
		}

		off_t expect_seq = seq_start;
		for (; it != dfs.end(); it++)
			if (it->first == expect_seq) {
				expect_seq++;
				ret.push_back(it->second);
			}
		dfs.erase(dfs.begin(), it);
	}
	else {
		for (auto it = dfs.begin(); it != dfs.end(); it++)
			ret.push_back(it->second);
		dfs.clear();
	}
	num_dfs = dfs.size();
	pthread_mutex_unlock(&lock);
	pthread_cond_broadcast(&add_cond);
	return ret;
}

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
	std::shared_ptr<char> lines;
	size_t size;
	const line_parser &parser;
	data_frame_set &dfs;
	off_t task_id;
public:
	parse_task(std::shared_ptr<char> _lines, size_t size,
			const line_parser &_parser, data_frame_set &_dfs,
			off_t task_id): parser(_parser), dfs(_dfs) {
		this->lines = _lines;
		this->size = size;
		this->task_id = task_id;
	}

	void run() {
		data_frame::ptr df = create_data_frame(parser);
		parse_lines(lines, size, parser, *df);
		// Each task processes a chunk of data, so we can use task Id as
		// the sequence number for the data frame.
		dfs.add(task_id, df);
	}
};

class file_parse_task: public thread_task
{
	text_io::ptr io;
	const line_parser &parser;
	data_frame_set &dfs;
public:
	file_parse_task(text_io::ptr io, const line_parser &_parser,
			data_frame_set &_dfs): parser(_parser), dfs(_dfs) {
		this->io = io;
	}

	void run();
};

void file_parse_task::run()
{
	while (!io->eof()) {
		size_t size = 0;
		std::shared_ptr<char> lines = io->read_lines(LINE_BLOCK_SIZE, size);
		if (size == 0)
			continue;

		data_frame::ptr df = create_data_frame(parser);
		parse_lines(lines, size, parser, *df);
		// In this case, we process multiple files simultaneously.
		// It's hard to determine the order of the data frames, so we just
		// ignore the orders.
		dfs.add(-1, df);
	}
}

}

data_frame::ptr read_lines(text_io::ptr io, const line_parser &parser,
		bool in_mem, bool sequential)
{
	data_frame::ptr df = create_data_frame(parser, in_mem);
	detail::mem_thread_pool::ptr mem_threads
		= detail::mem_thread_pool::get_global_mem_threads();
	const size_t MAX_PENDING = mem_threads->get_num_threads() * 3;
	data_frame_set dfs(MAX_PENDING);

	off_t task_id = 0;
	off_t seq_num = 0;
	while (!io->eof()) {
		size_t num_tasks = MAX_PENDING - mem_threads->get_num_pending();
		for (size_t i = 0; i < num_tasks && !io->eof(); i++) {
			size_t size = 0;
			std::shared_ptr<char> lines = io->read_lines(LINE_BLOCK_SIZE, size);

			if (size > 0)
				mem_threads->process_task(-1,
						new parse_task(lines, size, parser, dfs, task_id++));
		}
		if (dfs.get_num_dfs() > 0) {
			std::vector<data_frame::ptr> tmp_dfs = dfs.fetch_data_frames(
					sequential, seq_num);
			seq_num += tmp_dfs.size();
			if (!tmp_dfs.empty())
				df->append(tmp_dfs.begin(), tmp_dfs.end());
		}
	}
	mem_threads->wait4complete();
	std::vector<data_frame::ptr> tmp_dfs = dfs.fetch_data_frames(sequential,
			seq_num);
	seq_num += tmp_dfs.size();
	if (!tmp_dfs.empty())
		df->append(tmp_dfs.begin(), tmp_dfs.end());

	return df;
}

data_frame::ptr read_lines(const std::vector<text_io::ptr> &file_ios,
		const line_parser &parser, bool in_mem, bool sequential)
{
	if (file_ios.size() == 1)
		return read_lines(file_ios[0], parser, in_mem, sequential);

	data_frame::ptr df = create_data_frame(parser, in_mem);
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
	auto io_it = file_ios.begin();
	while (io_it != file_ios.end()) {
		size_t num_tasks = MAX_PENDING - mem_threads->get_num_pending();
		for (size_t i = 0; i < num_tasks && io_it != file_ios.end(); i++) {
			text_io::ptr io = *io_it;
			io_it++;
			mem_threads->process_task(-1,
					new file_parse_task(io, parser, dfs));
		}
		// This is the only thread that can fetch data frames from the queue.
		// If there are pending tasks in the thread pool, it's guaranteed
		// that we can fetch data frames from the queue.
		if (mem_threads->get_num_pending() > 0) {
			// If we process multiple files simultaneously, we don't need
			// to fetch the data frames sequentially.
			std::vector<data_frame::ptr> tmp_dfs = dfs.fetch_data_frames(
					false, -1);
			if (!tmp_dfs.empty())
				df->append(tmp_dfs.begin(), tmp_dfs.end());
		}
	}
	// TODO It might be expensive to calculate the number of pending
	// tasks every time.
	while (mem_threads->get_num_pending() > 0) {
		std::vector<data_frame::ptr> tmp_dfs = dfs.fetch_data_frames(false, -1);
		if (!tmp_dfs.empty())
			df->append(tmp_dfs.begin(), tmp_dfs.end());
	}
	mem_threads->wait4complete();
	// At this point, all threads have stoped working. If there are
	// data frames in the queue, they are the very last ones.
	if (dfs.get_num_dfs() > 0) {
		std::vector<data_frame::ptr> tmp_dfs = dfs.fetch_data_frames(false, -1);
		if (!tmp_dfs.empty())
			df->append(tmp_dfs.begin(), tmp_dfs.end());
	}

	return df;
}

/*
 * This parses one row of a dense matrix at a time.
 */
class row_parser: public line_parser
{
	const std::string delim;
	const size_t num_cols;
	std::vector<ele_parser::const_ptr> parsers;
	// If this vector contains elements, it contains how columns are duplicated.
	std::vector<off_t> dup_col_idxs;

	static std::string interpret_delim(const std::string &delim) {
		std::string new_delim = delim;
		if (delim == "\\t")
			new_delim = "\t";
		else if (delim == "\\n")
			new_delim = "\n";
		else if (delim == "\\r")
			new_delim = "\r";
		return new_delim;
	}
public:
	row_parser(const std::string &_delim,
			const std::vector<ele_parser::const_ptr> &_parsers,
			const std::vector<off_t> &dup_col_idxs): delim(
				interpret_delim(_delim)), num_cols(_parsers.size()) {
		this->parsers = _parsers;
		this->dup_col_idxs = dup_col_idxs;
		assert(dup_col_idxs.empty() || dup_col_idxs.size() == num_cols);
	}

	size_t parse(const std::vector<std::string> &lines, data_frame &df) const;
	size_t get_num_cols() const {
		return num_cols;
	}

	std::string get_col_name(off_t idx) const {
		return boost::str(boost::format("c%1%") % idx);
	}

	const scalar_type &get_col_type(off_t idx) const {
		return parsers[idx]->get_type();
	}
};

size_t row_parser::parse(const std::vector<std::string> &lines,
		data_frame &df) const
{
	std::vector<detail::smp_vec_store::ptr> cols(num_cols);
	for (size_t i = 0; i < cols.size(); i++)
		cols[i] = detail::smp_vec_store::create(lines.size(), get_col_type(i));
	size_t num_rows = 0;
	std::vector<std::string> strs;
	for (size_t i = 0; i < lines.size(); i++) {
		const char *line = lines[i].c_str();
		// Skip space
		for (; isspace(*line); line++);
		if (*line == '#')
			continue;

		// Split a string
		strs.clear();
		boost::split(strs, line, boost::is_any_of(delim));
		// If the line doesn't have enough values than expected, we fill
		// the remaining elements in the row with 0.
		while (strs.size() < num_cols)
			strs.push_back("0");

		// Parse each element.
		for (size_t j = 0; j < num_cols; j++) {
			// If the value is missing. We make it 0.
			if (strs[j].empty())
				parsers[j]->set_zero(cols[j]->get(num_rows));
			else
				parsers[j]->parse(strs[j], cols[j]->get(num_rows));
		}
		num_rows++;
	}
	for (size_t j = 0; j < num_cols; j++) {
		cols[j]->resize(num_rows);
		df.get_vec(j)->append(*cols[j]);
	}
	if (!dup_col_idxs.empty())
		for (size_t j = 0; j < num_cols; j++)
			df.get_vec(dup_col_idxs[j])->append(*cols[j]);
	return num_rows;
}

static std::shared_ptr<char> read_first_chunk(text_io::ptr io)
{
	// Read at max 1M
	long buf_size = 1024 * 1024;
	size_t read_bytes = 0;
	std::shared_ptr<char> buf = io->peek(buf_size, read_bytes);
	if (buf == NULL || read_bytes == 0)
		return std::shared_ptr<char>();
	buf.get()[read_bytes - 1] = 0;
	return buf;
}

static std::string detect_delim(std::shared_ptr<char> buf)
{
	char *res = strchr(buf.get(), '\t');
	if (res)
		return "\t";

	res = strchr(buf.get(), ' ');
	if (res)
		return " ";

	res = strchr(buf.get(), ',');
	if (res)
		return ",";

	BOOST_LOG_TRIVIAL(error) << "Cannot detect a known delimiter";
	return "";
}

static size_t detect_ncols(std::shared_ptr<char> buf, const std::string &delim)
{
	// Find the first line.
	char *res = strchr(buf.get(), '\n');
	// If the buffer doesn't have '\n'
	if (res == NULL) {
		BOOST_LOG_TRIVIAL(error)
			<< "read 1M data, can't find the end of the line";
		return std::numeric_limits<size_t>::max();
	}
	*res = 0;

	// Split a string
	std::string line = buf.get();
	std::vector<std::string> strs;
	boost::split(strs, line, boost::is_any_of(delim));
	return strs.size();
}

static std::string detect_delim(text_io::ptr io)
{
	auto buf = read_first_chunk(io);
	if (buf == NULL) {
		BOOST_LOG_TRIVIAL(error) << "can't read data to detect a delimiter";
		return "";
	}
	return detect_delim(buf);
}

static size_t detect_ncols(text_io::ptr io, const std::string &delim)
{
	auto buf = read_first_chunk(io);
	if (buf == NULL) {
		BOOST_LOG_TRIVIAL(error) << "can't read data to detect #cols";
		return std::numeric_limits<size_t>::max();
	}
	return detect_ncols(buf, delim);
}

static std::vector<text_io::ptr> files2ios(const std::vector<std::string> &files)
{
	std::vector<text_io::ptr> file_ios(files.size());
	for (size_t i = 0; i < file_ios.size(); i++)
		file_ios[i] = text_io::create(files[i]);
	return file_ios;
}

data_frame::ptr read_data_frame(const std::vector<std::string> &files,
		bool in_mem, bool sequential, const std::string &delim,
		const std::vector<ele_parser::const_ptr> &ele_parsers,
		const std::vector<off_t> &dup_col_idxs)
{
	auto file_ios = files2ios(files);
	if (file_ios.size() != files.size())
		return data_frame::ptr();
	return read_data_frame(file_ios, in_mem, sequential, delim, ele_parsers,
			dup_col_idxs);
}

data_frame::ptr read_data_frame(const std::vector<text_io::ptr> &file_ios,
		bool in_mem, bool sequential, const std::string &delim,
		const std::vector<ele_parser::const_ptr> &ele_parsers,
		const std::vector<off_t> &dup_col_idxs)
{
	std::string act_delim = delim;
	if (delim == "auto")
		act_delim = detect_delim(file_ios.front());
	if (act_delim.empty())
		return data_frame::ptr();

	std::shared_ptr<line_parser> parser = std::shared_ptr<line_parser>(
			new row_parser(act_delim, ele_parsers, dup_col_idxs));
	return read_lines(file_ios, *parser, in_mem, sequential);
}

dense_matrix::ptr read_matrix(const std::vector<std::string> &files,
		bool in_mem, bool sequential, const std::string &ele_type,
		const std::string &delim, size_t num_cols)
{
	auto file_ios = files2ios(files);
	if (file_ios.size() != files.size())
		return dense_matrix::ptr();
	return read_matrix(file_ios, in_mem, sequential, ele_type, delim,
			num_cols);
}

dense_matrix::ptr read_matrix(const std::vector<text_io::ptr> &file_ios,
		bool in_mem, bool sequential, const std::string &ele_type,
		const std::string &delim, size_t num_cols)
{
	// Detect the delimiter
	std::string act_delim = delim;
	if (delim == "auto")
		act_delim = detect_delim(file_ios.front());
	if (act_delim.empty())
		return dense_matrix::ptr();

	// Detect the number of columns.
	if (num_cols == std::numeric_limits<size_t>::max())
		num_cols = detect_ncols(file_ios.front(), act_delim);
	if (num_cols == std::numeric_limits<size_t>::max())
		return dense_matrix::ptr();

	std::shared_ptr<line_parser> parser;
	std::vector<ele_parser::const_ptr> ele_parsers(num_cols);
	for (size_t i = 0; i < num_cols; i++) {
		ele_parsers[i] = get_ele_parser(ele_type);
		if (ele_parsers[i] == NULL)
			return dense_matrix::ptr();
	}
	parser = std::shared_ptr<line_parser>(new row_parser(act_delim,
				ele_parsers, std::vector<off_t>()));
	data_frame::ptr df = read_lines(file_ios, *parser, in_mem, sequential);
	return dense_matrix::create(df);
}

dense_matrix::ptr read_matrix(const std::vector<std::string> &files,
		bool in_mem, bool sequential, const std::string &ele_type,
		const std::string &delim, const std::string &col_indicator)
{
	auto file_ios = files2ios(files);
	if (file_ios.size() != files.size())
		return dense_matrix::ptr();
	return read_matrix(file_ios, in_mem, sequential, ele_type, delim,
			col_indicator);
}

dense_matrix::ptr read_matrix(const std::vector<text_io::ptr> &file_ios,
		bool in_mem, bool sequential, const std::string &ele_type,
		const std::string &delim, const std::string &col_indicator)
{
	std::string act_delim = delim;
	if (delim == "auto")
		act_delim = detect_delim(file_ios.front());
	if (act_delim.empty())
		return dense_matrix::ptr();

	std::vector<std::string> strs;
	boost::split(strs, col_indicator, boost::is_any_of(" "));
	std::vector<ele_parser::const_ptr> ele_parsers(strs.size());
	assert(strs.size());
	for (size_t i = 0; i < ele_parsers.size(); i++) {
		ele_parsers[i] = get_ele_parser(strs[i]);
		if (ele_parsers[i] == NULL)
			return dense_matrix::ptr();
	}

	for (size_t i = 1; i < ele_parsers.size(); i++)
		if (ele_parsers[i]->get_type() != ele_parsers[0]->get_type()) {
			BOOST_LOG_TRIVIAL(error) << "element parsers output different types";
			return dense_matrix::ptr();
		}

	std::shared_ptr<line_parser> parser = std::shared_ptr<line_parser>(
			new row_parser(act_delim, ele_parsers, std::vector<off_t>()));
	data_frame::ptr df = read_lines(file_ios, *parser, in_mem, sequential);
	return dense_matrix::create(df);
}

ele_parser::const_ptr get_ele_parser(const std::string &type)
{
	if (type == "B" || type.empty())
		return ele_parser::const_ptr();
	else if (type == "I")
		return ele_parser::const_ptr(new int_parser<int>());
	else if (type == "L")
		return ele_parser::const_ptr(new int_parser<long>());
	else if (type == "F")
		return ele_parser::const_ptr(new float_parser<float>());
	else if (type == "D")
		return ele_parser::const_ptr(new float_parser<double>());
	else if (type == "H")
		return ele_parser::const_ptr(new int_parser<int>(16));
	else if (type == "LH")
		return ele_parser::const_ptr(new int_parser<long>(16));
	else {
		BOOST_LOG_TRIVIAL(error) << "unknown element parser: " << type;
		return ele_parser::const_ptr();
	}
}

const scalar_type &get_ele_type(const std::string &type_name)
{
	if (type_name == "B")
		return get_scalar_type<bool>();
	else if (type_name == "I")
		return get_scalar_type<int>();
	else if (type_name == "L")
		return get_scalar_type<long>();
	else if (type_name == "F")
		return get_scalar_type<float>();
	else if (type_name == "D")
		return get_scalar_type<double>();
	else if (type_name == "H")
		return get_scalar_type<int>();
	else if (type_name == "LH")
		return get_scalar_type<long>();
	else
		throw std::invalid_argument("unknown element type");
}

bool valid_ele_type(const std::string &type_name)
{
	return type_name == "B" || type_name == "I" || type_name == "L"
		|| type_name == "F" || type_name == "D" || type_name == "H"
		|| type_name == "LH";
}

}
