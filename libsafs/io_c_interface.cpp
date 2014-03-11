/**
 * Copyright 2013 Da Zheng
 *
 * This file is part of SAFSlib.
 *
 * SAFSlib is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * SAFSlib is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with SAFSlib.  If not, see <http://www.gnu.org/licenses/>.
 */

#include <fcntl.h>
#include <stdlib.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>

#include <vector>
#include <tr1/unordered_map>

#include "slab_allocator.h"
#include "io_interface.h"
#include "concurrency.h"
#include "RAID_config.h"
#include "cache.h"
#include "cache_config.h"
extern "C" {
#include "io_c_interface.h"
}

static pthread_mutex_t mutex = PTHREAD_MUTEX_INITIALIZER;
static std::tr1::unordered_map<std::string, file_io_factory *> opened_files;

struct data_fill_struct
{
	pthread_t tid;
	int fd;
	off_t off;
	size_t tot_size;
};

struct ssd_file_desc
{
	io_interface *io;
	file_io_factory *factory;
};

static void *fill_space(void *arg)
{
	struct data_fill_struct *data = (struct data_fill_struct *) arg;
	size_t tot_size = data->tot_size;
	int fd = data->fd;
	off_t off = data->off;

	size_t block_size = 2 * 1024 * 1024;
	char *buf = (char *) valloc(block_size);
	memset(buf, 0, block_size);
	off_t ret = lseek(fd, off, SEEK_SET);
	if (ret < 0) {
		perror("lseek");
		assert(0);
	}
	size_t remaining_size = tot_size;
	while (remaining_size > 0) {
		size_t count = min(remaining_size, block_size);
		ssize_t written = write(fd, buf, count);
		if (written < 0) {
			perror("write");
			assert(0);
		}
		if (written == 0) {
			printf("WARNING! nothing was written\n");
		}
		remaining_size -= written;
	}
	free(buf);
	return NULL;
}

static size_t get_filesize(int fd)
{
	struct stat stat_buf;
	if (fstat(fd, &stat_buf) != 0) {
		perror("fstat");
		assert(0);
	}
	return stat_buf.st_size;
}

class user_callback: public callback
{
	ssd_callback_func_t cb;
	void *cb_data;
public:
	user_callback(ssd_callback_func_t cb, void *cb_data) {
		this->cb = cb;
		this->cb_data = cb_data;
	}

	int invoke(io_request *rqs[], int num) {
		for (int i = 0; i < num; i++) {
			cb(rqs[i]->get_offset(), rqs[i]->get_buf(), rqs[i]->get_size(),
					cb_data, 0);
		}
		return 0;
	}
};

static config_map c_global_map;

extern "C" {

void set_cache_size(long size)
{
	std::map<std::string, std::string> parameter;
	parameter.insert(std::pair<std::string, std::string>("cache_size", itoa(size)));
	params.init(parameter);
}

void set_cache_type(int type)
{
	std::map<std::string, std::string> parameter;
	parameter.insert(std::pair<std::string, std::string>("cache_type", itoa(type)));
	params.init(parameter);
}

void set_RAID_block_size(int num_pages)
{
	std::map<std::string, std::string> parameter;
	parameter.insert(std::pair<std::string, std::string>("RAID_block_size",
				itoa(num_pages * PAGE_SIZE)));
	params.init(parameter);
}

void ssd_init_io_system(const char *root_conf_file, int *node_ids, int num_nodes)
{
	// TODO I need to pass the options to the IO system.
	assert(0);
	init_io_system(c_global_map);
}

void ssd_file_io_init(const char *name, int flags, int num_threads, int num_nodes,
		int *suggested_nodes)
{
	pthread_mutex_lock(&mutex);
	if (opened_files.find(name) != opened_files.end()) {
		pthread_mutex_unlock(&mutex);
		return;
	}

	printf("Init SSDIO with %d threads and %d nodes\n",
			num_threads, num_nodes);

	// TODO the suggested nodes should be added to the system parameters
	// in order to pass to the IO system.
	assert(0);
#if 0
	std::vector<int> node_id_array;
	// Users can suggest nodes where the IO should be. It make sense for cached IO
	// because we are going to place cache on those nodes.
	if (suggested_nodes) {
		for (int i = 0; i < num_nodes; i++)
			node_id_array.push_back(suggested_nodes[i]);
	}
	else {
		for (int i = 0; i < num_nodes; i++)
			node_id_array.push_back(i);
	}
	printf("There are %ld nodes\n", node_id_array.size());
#endif

	int access_option = GLOBAL_CACHE_ACCESS;
	if (flags | O_DIRECT)
		access_option = REMOTE_ACCESS;

	file_io_factory::shared_ptr factory = create_io_factory(name, access_option);
	// The shared pointer is going to destroy the I/O factory.
	assert(0);
	opened_files.insert(std::pair<std::string, file_io_factory *>(name,
				factory.get()));

	pthread_mutex_unlock(&mutex);
}

int ssd_create(const char *name, size_t tot_size)
{
	printf("create %s, size: %ld\n", name, tot_size);
	std::vector<file_info> files;
	retrieve_data_files(name, files);
	size_t file_size = tot_size / files.size();
	struct data_fill_struct data[files.size()];
	for (unsigned i = 0; i < files.size(); i++) {
		int fd = open(files[i].name.c_str(), O_RDWR | O_CREAT | O_DIRECT,
				S_IRUSR | S_IWUSR | S_IRGRP | S_IROTH);
		if (fd < 0) {
			perror("open");
			assert(0);
		}

		size_t curr_size = get_filesize(fd);
		if (file_size > curr_size) {
			int err = posix_fallocate(fd, 0, file_size);
			if (err) {
				fprintf(stderr, "posix_fallocate error: %s\n", strerror(err));
				assert(0);
			}
			data[i].fd = fd;
			data[i].off = 0;
			data[i].tot_size = file_size;
			pthread_create(&data[i].tid, NULL, fill_space, &data[i]);
		}
		else {
			close(fd);
			data[i].fd = -1;
		}
	}
	for (unsigned i = 0; i < files.size(); i++) {
		if (data[i].fd >= 0) {
			pthread_join(data[i].tid, NULL);
			close(data[i].fd);
			printf("fill the file %s\n", files[i].name.c_str());
		}
	}
	return 0;
}

ssd_file_desc_t ssd_open(const char *name, int node_id, int flags)
{
	thread *t = thread::represent_thread(node_id);
	file_io_factory *factory = opened_files[name];
	io_interface *io = factory->create_io(t);
	assert(io);
	ssd_file_desc_t desc = new struct ssd_file_desc;
	desc->io = io;
	desc->factory = factory;
	return desc;
}

ssize_t ssd_read(ssd_file_desc_t fd, void *buf, size_t count, off_t off)
{
	io_interface *io = fd->io;
	io_status status = io->access((char *) buf, off, count, READ);
	assert(status == IO_OK);
	return status.get_priv_data();
}

ssize_t ssd_write(ssd_file_desc_t fd, void *buf, size_t count, off_t off)
{
	io_interface *io = fd->io;
	io_status status = io->access((char *) buf, off, count, WRITE);
	assert(status == IO_OK);
	return status.get_priv_data();
}

int ssd_close(ssd_file_desc_t fd)
{
	io_interface *io = fd->io;
	io->cleanup();
	fd->factory->destroy_io(io);
	// TODO I need to free the space of fd.
	return 0;
}

int ssd_delete(const char *name)
{
	std::vector<file_info> files;
	retrieve_data_files(name, files);
	for (unsigned i = 0; i < files.size(); i++) {
		int ret = unlink(files[i].name.c_str());
		if (ret < 0) {
			perror("unlink");
			assert(0);
		}
	}
	return 0;
}

size_t ssd_get_filesize(const char *name)
{
	std::vector<file_info> files;
	retrieve_data_files(name, files);
	size_t tot_size = 0;
	for (unsigned i = 0; i < files.size(); i++) {
		struct stat stat_buf;
		if (stat(files[i].name.c_str(), &stat_buf) != 0) {
			perror("stat");
			assert(0);
		}
		tot_size += stat_buf.st_size;
	}
	return tot_size;
}

void ssd_set_callback(ssd_file_desc_t fd, ssd_callback_func_t cb, void *cb_data)
{
	io_interface *io = fd->io;
	assert(io->support_aio());
	io->set_callback(new user_callback(cb, cb_data));
}

ssize_t ssd_aread(ssd_file_desc_t fd, void *buf, size_t count, off_t off)
{
	io_interface *io = fd->io;
	assert(io->support_aio());
	data_loc_t loc(io->get_file_id(), off);
	io_request req((char *) buf, loc, count, READ, io, io->get_node_id());
	io->access(&req, 1);
	return 0;
}

ssize_t ssd_awrite(ssd_file_desc_t fd, void *buf, size_t count, off_t off)
{
	io_interface *io = fd->io;
	assert(io->support_aio());
	data_loc_t loc(io->get_file_id(), off);
	io_request req((char *) buf, loc, count, WRITE, io, io->get_node_id());
	io->access(&req, 1);
	return 0;
}

int ssd_wait(ssd_file_desc_t fd, int num)
{
	io_interface *io = fd->io;
	assert(io->support_aio());
	return io->wait4complete(num);
}

int ssd_get_io_slots(ssd_file_desc_t fd)
{
	io_interface *io = fd->io;
	assert(io->support_aio());
	return io->get_remaining_io_slots();
}

int ssd_fd_node_id(ssd_file_desc_t fd)
{
	io_interface *io = fd->io;
	return io->get_node_id();
}

struct buf_pool
{
	slab_allocator *alloc;
	int entry_size;
};

struct buf_pool *create_buf_pool(int obj_size, long pool_size, int node_id)
{
	struct buf_pool *pool = new struct buf_pool;
	pool->alloc = new slab_allocator(std::string("buf_pool-") + itoa(node_id),
			obj_size, obj_size * 1000, pool_size, node_id);
	pool->entry_size = obj_size;
	return pool;
}

void *alloc_buf(struct buf_pool *pool)
{
	return pool->alloc->alloc();
}

void free_buf(struct buf_pool *pool, void *buf)
{
	pool->alloc->free((char *) buf);
}

void destroy_buf_pool(struct buf_pool *pool)
{
	delete pool->alloc;
	delete pool;
}

}
