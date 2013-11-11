#ifndef _GNU_SOURCE
#define _GNU_SOURCE
#endif

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

#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <pthread.h>
#include <signal.h>
#ifdef PROFILER
#include <google/profiler.h>
#endif

#include "io_c_interface.h"

enum {
	READ,
	WRITE,
};

off_t *offs;
int num_threads = 16;
char *prof_file = "test_c_interface.prof";
long data_size = 1024L * 1024 * 4096;
int block_size = 4096;
int num_nodes = 1;
int sync = 0;
int access = READ;
int file_per_proc = 0;

struct thread_data
{
	pthread_t tid;
	char file_name[128];
	int node_id;
	int off_start;
	int num;
	int idx;
};

void rand_permute_array(off_t arr[], int num)
{
	int i;
	for (i = num - 1; i >= 1; i--) {
		int j = random() % i;
		off_t tmp = arr[j];
		arr[j] = arr[i];
		arr[i] = tmp;
	}
}

struct callback_data
{
	char *buffer;
	volatile int *num_completes;
	struct buf_pool *buf_allocator;
};

void cb_func(off_t off, void *buf, int size, void *cb_data, int status)
{
	struct callback_data *data = cb_data;
	struct buf_pool *buf_allocator = data->buf_allocator;

	__sync_fetch_and_add(data->num_completes, 1);
	free_buf(buf_allocator, buf);
}

void *SyncThreadWriteOrRead(void *arg)
{
	struct thread_data *data = arg;
	bind2node_id(data->node_id);
	int node_id = data->node_id;
	int num = data->num;
	int off_start = data->off_start;
	int i;
	char *buffer = (char *) numa_alloc_onnode(block_size, node_id);

	ssd_file_desc_t fd = ssd_open(data->file_name, node_id, 0);

	printf("thread %d: access %d blocks\n", data->idx, num);
	for (i = 0; i < num; i++) {
		off_t offset = offs[off_start + i];
		if (access == READ)
			ssd_read(fd, (void *) buffer, block_size, offset);
		else
			ssd_write(fd, (void *) buffer, block_size, offset);
	}
	numa_free(buffer, block_size);

	ssd_close(fd);
	return NULL;
}

void *AsyncThreadWriteOrRead(void *arg)
{
	struct thread_data *data = arg;
	bind2node_id(data->node_id);
	int node_id = data->node_id;
	int num_completes = 0;
	int num = data->num;
	int off_start = data->off_start;
	int i;

	ssd_file_desc_t fd = ssd_open(data->file_name, node_id, 0);

	struct buf_pool *buf_allocator = create_buf_pool(block_size,
			block_size * 10000, node_id);
	struct callback_data *cb_data = (struct callback_data *) malloc(sizeof(*cb_data));
	assert(cb_data);
	cb_data->num_completes = &num_completes;
	cb_data->buf_allocator = buf_allocator;

	ssd_set_callback(fd, cb_func, cb_data);

	for (i = 0; i < num; i++) {
		char *buffer = alloc_buf(buf_allocator);
		assert(buffer);
		off_t offset = offs[off_start + i];
		if (access == READ)
			ssd_aread(fd, (void *) buffer, block_size, offset);
		else
			ssd_awrite(fd, (void *) buffer, block_size, offset);
		if (ssd_get_io_slots(fd) == 0) {
			ssd_wait(fd, 1);
		}
	}

	ssd_close(fd);
	return NULL;
}

void int_handler(int sig_num)
{
#ifdef PROFILER
	ProfilerStop();
#endif
	exit(0);
}

int main(int argc, char *argv[])
{
	int i;
	int num_offs;
	char *root_conf_file;
	char *data_file;

	if (argc < 3) {
		fprintf(stderr, "test root_conf data_file\n");
		return -1;
	}
	root_conf_file = argv[1];
	data_file = argv[2];

	if (file_per_proc)
		num_offs = data_size / block_size / num_threads;
	else
		num_offs = data_size / block_size;
	offs = (off_t *) malloc(num_offs * sizeof(off_t));
	for (i = 0; i < num_offs; i++)
		offs[i] = i * (long) block_size;
	rand_permute_array(offs, num_offs);
	printf("generate random offsets\n");

	signal(SIGINT, int_handler);
	int node_ids[num_nodes];
	for (i = 0; i < num_nodes; i++)
		node_ids[i] = i;
	ssd_init_io_system(root_conf_file, node_ids, num_nodes);
	if (file_per_proc)
		set_cache_size(512 * 1024 * 1024 / num_threads);
	printf("init IO system\n");

	char file_name_buf[128];
	char *file_name;

#ifdef PROFILER
	ProfilerStart(prof_file);
#endif
	struct thread_data data[num_threads];
	long start = get_curr_ms();
	for (i = 0; i < num_threads; i++) {
		data[i].idx = i;
		data[i].node_id = i % num_nodes;
		if (file_per_proc) {
			data[i].off_start = 0;
			data[i].num = num_offs;
		}
		else {
			data[i].num = num_offs / num_threads;
			data[i].off_start = i * data[i].num; 
		}

		if (file_per_proc) {
			sprintf(file_name_buf, "%s.%08d", data_file, data[i].idx);
			file_name = file_name_buf;
			int suggested_nodes[1] = {data->node_id};
			ssd_file_io_init(file_name, 0, 1, 1, suggested_nodes);
		}
		else {
			file_name = data_file;
			ssd_file_io_init(file_name, 0, num_threads, num_nodes, NULL);
		}
		printf("init file IO\n");

		printf("open file %s\n", file_name);
		strncpy(data[i].file_name, file_name, sizeof(data[i].file_name));
		int ret;
		if (sync)
			ret = pthread_create(&data[i].tid, NULL, SyncThreadWriteOrRead,
					(void *) &data[i]);
		else
			ret = pthread_create(&data[i].tid, NULL, AsyncThreadWriteOrRead,
					(void *) &data[i]);
		assert(ret == 0);
	}
	for (i = 0; i < num_threads; i++)
		pthread_join(data[i].tid, NULL);
	long end = get_curr_ms();
#ifdef PROFILER
	ProfilerStop();
#endif
	printf("It takes %ld ms\n", end - start);
}
