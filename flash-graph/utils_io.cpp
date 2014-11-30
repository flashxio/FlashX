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
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>

#include <unordered_map>
#include <boost/format.hpp>

#include "log.h"
#include "common.h"

#include "safs_file.h"
#include "io_interface.h"

#include "utils.h"

namespace fg
{

static size_t buf_cap = 128 * 1024 * 1024;

///////////////////////////large I/O for native filesystem//////////////////////

/*
 * This class is optimized for writing large amount of data to Linux filesystem.
 */
class native_large_writer: public large_writer
{
	int fd;
	char *write_buf;
	size_t write_bytes;
	size_t tot_write_bytes;
	std::string file_name;

	native_large_writer(const std::string &file);
	void close_file() {
		flush();
		close(fd);
		fd = -1;
		write_bytes = 0;
		tot_write_bytes = 0;
	}
public:
	~native_large_writer() {
		close_file();
		free(write_buf);
	}

	static large_writer::ptr create(const std::string &file) {
		return ptr(new native_large_writer(file));
	}

	off_t seek(off_t off, int whence) {
		if (fd < 0)
			return -1;
		// If there are data buffered, let's not move to another location.
		if (write_bytes > 0)
			flush();

		return lseek(fd, off, whence);
	}
	ssize_t flush();
	ssize_t write(const char *buf, size_t bytes);
	size_t get_write_bytes() const {
		return tot_write_bytes;
	}

	virtual int delete_file() {
		if (fd < 0)
			return -1;
		close_file();
		return unlink(file_name.c_str());
	}

	virtual int rename2(const std::string &new_name) {
		if (fd < 0)
			return -1;
		close_file();
		int ret = rename(file_name.c_str(), new_name.c_str());
		file_name = new_name;
		return ret;
	}
};

ssize_t native_large_writer::flush()
{
	if (fd < 0)
		return -1;


	if (write_bytes == 0)
		return 0;
	// TODO
	// We might write some blank data to disks. If the file contains data
	// in the location where the blank data is written to, we should read
	// the data from the file first.
	size_t bytes = ROUNDUP(write_bytes, 512);
	size_t ret = bytes;
	char *tmp = write_buf;
	do {
		ssize_t ret = ::write(fd, tmp, bytes);
		if (ret < 0) {
			BOOST_LOG_TRIVIAL(error) << boost::format(
					"fail to write %1% bytes: %2%") % bytes % strerror(errno);
			return ret;
		}
		tmp += ret;
		bytes -= ret;
	} while (bytes > 0);
	write_bytes = 0;
	return ret;
}

ssize_t native_large_writer::write(const char *buf, size_t bytes)
{
	if (fd < 0)
		return -1;

	ssize_t ret = 0;
	do {
		size_t remain = buf_cap - write_bytes;
		size_t copy_bytes = std::min(remain, bytes);
		memcpy(write_buf + write_bytes, buf, copy_bytes);
		buf += copy_bytes;
		bytes -= copy_bytes;
		write_bytes += copy_bytes;
		if (write_bytes == buf_cap) {
			ssize_t ret1 = flush();
			if (ret1 < 0)
				return ret1;
		}
		ret += copy_bytes;
	} while (bytes > 0);
	tot_write_bytes += ret;
	return ret;
}

native_large_writer::native_large_writer(const std::string &file)
{
	file_name = file;
	fd = open(file.c_str(), O_WRONLY | O_CREAT | O_DIRECT, S_IRUSR | S_IWUSR | S_IRGRP);
	if (fd < 0) {
		BOOST_LOG_TRIVIAL(error) << boost::format(
				"fail to open %1%: %2%") % file % strerror(errno);
		exit(1);
	}

	write_buf = (char *) valloc(buf_cap);
	write_bytes = 0;
	tot_write_bytes = 0;
}

class native_large_reader: public large_reader
{
	FILE *f;

	native_large_reader(const std::string &file) {
		f = fopen(file.c_str(), "r");
		if (f == NULL) {
			BOOST_LOG_TRIVIAL(error) << boost::format(
					"fail to open %1%: %2%") % file % strerror(errno);
			exit(1);
		}
	}
public:
	static large_reader::ptr create(const std::string &file) {
		return ptr(new native_large_reader(file));
	}

	virtual ~native_large_reader() {
		fclose(f);
	}
	virtual ssize_t read(char *buf, size_t bytes) {
		ssize_t ret = fread(buf, bytes, 1, f);
		if (ret == 1)
			return bytes;
		else
			return ret;
	}
	virtual off_t seek(off_t off, int whence) {
		if (fseek(f, off, whence) == 0)
			return ftell(f);
		else
			return -1;
	}

	virtual bool is_safs() {
		return false;
	}
};

class native_large_io_creator: public large_io_creator
{
	const std::string curr_dir;
public:
	native_large_io_creator(const std::string &_curr_dir): curr_dir(_curr_dir) {
	}

	virtual large_writer::ptr create_writer(const std::string &file) {
		std::string path = curr_dir + "/" + file;
		return native_large_writer::create(path);
	}

	virtual large_reader::ptr create_reader(const std::string &file) {
		std::string path = curr_dir + "/" + file;
		return native_large_reader::create(path);
	}
};

///////////////////////////////large I/O for SAFS///////////////////////////////

struct buf_deleter
{
	void operator()(char *buf) {
		free(buf);
	}
};
typedef std::unique_ptr<char[], buf_deleter> align_buf_ptr;

/*
 * This callback only needs to make sure the buffers are free'd when the data
 * is written to SAFS.
 */
class large_writer_callback: public safs::callback
{
	std::unordered_map<char *, align_buf_ptr> buf_map;
public:
	typedef std::shared_ptr<large_writer_callback> ptr;

	~large_writer_callback() {
		assert(buf_map.empty());
	}

	void add_buf(align_buf_ptr buf) {
		buf_map.insert(std::pair<char *, align_buf_ptr>(buf.get(),
					std::move(buf)));
	}

	virtual int invoke(safs::io_request *reqs[], int num) {
		for (int i = 0; i < num; i++)
			buf_map.erase(reqs[i]->get_buf());
		return 0;
	}
};

class safs_large_writer: public large_writer
{
	static const int MAX_PENDING_IOS = 8;
	align_buf_ptr write_buf;
	size_t write_bytes;
	size_t tot_write_bytes;

	// The location where the data is written to next time.
	off_t curr_off;

	safs::file_io_factory::shared_ptr factory;
	safs::io_interface::ptr io;
	large_writer_callback::ptr cb;

	safs_large_writer(const std::string &file) {
		factory = safs::create_io_factory(file, safs::REMOTE_ACCESS);
		write_bytes = 0;
		write_buf = align_buf_ptr((char *) valloc(buf_cap));
		tot_write_bytes = 0;
		curr_off = 0;
	}

	void open_file() {
		io = factory->create_io(thread::get_curr_thread());
		cb = large_writer_callback::ptr(new large_writer_callback());
		io->set_callback(std::static_pointer_cast<safs::callback>(cb));
	}

	void close_file() {
		flush();
		if (io) {
			io->wait4complete(io->num_pending_ios());
			io->cleanup();
			io = NULL;
		}
		cb = NULL;
		factory = NULL;
		write_bytes = 0;
		tot_write_bytes = 0;
		curr_off = 0;
	}
public:
	static large_writer::ptr create(const std::string &file) {
		safs::safs_file f(safs::get_sys_RAID_conf(), file);
		if (!f.exist()) {
			bool ret = f.create_file(0);
			if (!ret)
				return ptr();
		}
		return ptr(new safs_large_writer(file));
	}

	virtual ~safs_large_writer() {
		close_file();
	}

	virtual off_t seek(off_t off, int whence) {
		if (factory == NULL)
			return -1;

		// If there are data buffered, let's not move to another location.
		if (write_bytes > 0)
			flush();

		if (whence == SEEK_SET)
			this->curr_off = off;
		else if (whence == SEEK_CUR)
			this->curr_off += off;
		else
			return -1;
		assert(curr_off % 512 == 0);
		return curr_off;
	}
	virtual ssize_t flush();
	virtual ssize_t write(const char *buf, size_t bytes);

	virtual size_t get_write_bytes() const {
		return tot_write_bytes;
	}

	virtual int delete_file() {
		if (factory == NULL)
			return -1;
		std::string file_name = factory->get_name();
		close_file();
		write_buf = NULL;
		safs::safs_file f(safs::get_sys_RAID_conf(), file_name);
		if (f.delete_file())
			return 0;
		else
			return -1;
	}

	virtual int rename2(const std::string &new_name) {
		if (factory == NULL)
			return -1;
		std::string file_name = factory->get_name();
		close_file();
		safs::safs_file f(safs::get_sys_RAID_conf(), file_name);
		if (f.rename(new_name)) {
			file_name = new_name;
			return 0;
		}
		else
			return -1;
	}
};

ssize_t safs_large_writer::flush()
{
	if (factory == NULL)
		return -1;

	if (write_bytes == 0)
		return 0;
	if (io == NULL)
		open_file();
	// TODO
	// We might write some blank data to disks. If the file contains data
	// in the location where the blank data is written to, we should read
	// the data from the file first.
	size_t bytes = ROUNDUP(write_bytes, 512);
	safs::data_loc_t loc(io->get_file_id(), curr_off);
	safs::io_request req(write_buf.get(), loc, bytes, WRITE);
	io->access(&req, 1);
	io->flush_requests();
	cb->add_buf(std::move(write_buf));
	write_bytes = 0;
	curr_off += bytes;
	write_buf = align_buf_ptr((char *) valloc(buf_cap));
	if (io->num_pending_ios() > MAX_PENDING_IOS)
		io->wait4complete(1);
	return bytes;
}

ssize_t safs_large_writer::write(const char *buf, size_t bytes)
{
	if (factory == NULL)
		return -1;

	ssize_t ret = 0;
	do {
		size_t remain = buf_cap - write_bytes;
		size_t copy_bytes = std::min(remain, bytes);
		memcpy(write_buf.get() + write_bytes, buf, copy_bytes);
		buf += copy_bytes;
		bytes -= copy_bytes;
		write_bytes += copy_bytes;
		if (write_bytes == buf_cap) {
			ssize_t ret1 = flush();
			if (ret1 < 0)
				return ret1;
		}
		ret += copy_bytes;
	} while (bytes > 0);
	tot_write_bytes += ret;
	return ret;
}

class safs_large_reader: public large_reader
{
	off_t curr_off;
	safs::file_io_factory::shared_ptr factory;
	safs::io_interface::ptr io;
	const size_t max_req_size;

	safs_large_reader(const std::string &file): max_req_size(std::min(
				128UL * 1024 * 1024, safs::io_request::get_max_req_size())) {
		factory = safs::create_io_factory(file, safs::GLOBAL_CACHE_ACCESS);
		curr_off = 0;
	}

	void open_file() {
		io = factory->create_io(thread::get_curr_thread());
	}
public:
	static large_reader::ptr create(const std::string &file) {
		safs::safs_file f(safs::get_sys_RAID_conf(), file);
		if (!f.exist())
			return ptr();
		return ptr(new safs_large_reader(file));
	}

	virtual ~safs_large_reader() {
		io = NULL;
	}
	virtual ssize_t read(char *buf, size_t bytes) {
		if (io == NULL)
			open_file();

		ssize_t ret = bytes;
		while (bytes > 0) {
			size_t req_size = std::min(bytes, max_req_size);
			io->access(buf, curr_off, req_size, READ);
			bytes -= req_size;
			buf += req_size;
			curr_off += req_size;
		}
		return ret;
	}

	virtual off_t seek(off_t off, int whence) {
		if (whence == SEEK_SET)
			this->curr_off = off;
		else if (whence == SEEK_CUR)
			this->curr_off += off;
		else
			return -1;
		return curr_off;
	}

	virtual bool is_safs() {
		return true;
	}
};

class safs_large_io_creator: public large_io_creator
{
public:
	virtual large_writer::ptr create_writer(const std::string &file) {
		if (safs::params.is_writable())
			return safs_large_writer::create(file);
		else
			return large_writer::ptr();
	}
	virtual large_reader::ptr create_reader(const std::string &file) {
		return safs_large_reader::create(file);
	}
};

large_io_creator::ptr large_io_creator::create(bool safs, const std::string &curr_dir)
{
	if (safs) {
		if (safs::is_safs_init())
			return ptr(new safs_large_io_creator());
		else
			return ptr();
	}
	else
		return ptr(new native_large_io_creator(curr_dir));
}

void set_buf_cap(size_t new_cap)
{
	buf_cap = new_cap;
}

}
