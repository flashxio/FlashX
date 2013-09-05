#ifndef _GNU_SOURCE
#define _GNU_SOURCE
#endif

#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <pthread.h>
#include <signal.h>
#include <google/profiler.h>

#include "io_c_interface.h"

enum {
	READ,
	WRITE,
};

off_t *offs;
int num_threads = 16;
char *file_name_root = "../conf/data_files.txt";
char *prof_file = "test_c_interface.prof";
long data_size = 1024L * 1024 * 4096 * 32;
int block_size = 4096;
int num_nodes = 1;
int sync = 0;
int access = READ;
int file_per_proc = 1;

struct thread_data
{
	pthread_t tid;
	ssd_file_desc_t fd;
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
	struct buf_pool *cb_allocator;
};

void cb_func(void *arg, int status)
{
	struct callback_data *data = arg;
	struct buf_pool *buf_allocator = data->buf_allocator;
	struct buf_pool *cb_allocator = data->cb_allocator;

	__sync_fetch_and_add(data->num_completes, 1);
	free_buf(buf_allocator, data->buffer);
	free_buf(cb_allocator, data);
}

void *SyncThreadWriteOrRead(void *arg)
{
	struct thread_data *data = arg;
	bind2node_id(data->node_id);
	ssd_file_desc_t fd = data->fd;
	int node_id = ssd_fd_node_id(fd);
	assert(node_id == data->node_id);
	int num = data->num;
	int off_start = data->off_start;
	int i;
	char *buffer = (char *) numa_alloc_onnode(block_size, node_id);

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
	ssd_file_desc_t fd = data->fd;
	int node_id = ssd_fd_node_id(fd);
	assert(node_id == data->node_id);
	int num_completes = 0;
	int num = data->num;
	int off_start = data->off_start;
	int i;

	struct buf_pool *buf_allocator = create_buf_pool(block_size,
			block_size * 10000, node_id);
	struct buf_pool *cb_allocator = create_buf_pool(sizeof(struct callback_data),
			sizeof(struct callback_data) * 10000, node_id);

	ssd_set_callback(fd, cb_func);

	for (i = 0; i < num; i++) {
		char *buffer = alloc_buf(buf_allocator);
		off_t offset = offs[off_start + i];
		struct callback_data *cb_data = alloc_buf(cb_allocator);
		cb_data->buffer = buffer;
		cb_data->num_completes = &num_completes;
		cb_data->buf_allocator = buf_allocator;
		cb_data->cb_allocator = cb_allocator;
		if (access == READ)
			ssd_aread(fd, (void *) buffer, block_size, offset, (void *) cb_data);
		else
			ssd_awrite(fd, (void *) buffer, block_size, offset, (void *) cb_data);
	}

	ssd_close(fd);
	return NULL;
}

void int_handler(int sig_num)
{
	ProfilerStop();
	exit(0);
}

int main()
{
	int i;
	int num_offs;

	if (file_per_proc)
		num_offs = data_size / block_size / num_threads;
	else
		num_offs = data_size / block_size;
	offs = (off_t *) malloc(num_offs * sizeof(off_t));
	for (i = 0; i < num_offs; i++)
		offs[i] = i * (long) block_size;
	rand_permute_array(offs, num_offs);

	signal(SIGINT, int_handler);
	int node_ids[num_nodes];
	for (i = 0; i < num_nodes; i++)
		node_ids[i] = i;
	ssd_init_io_system(file_name_root, node_ids, num_nodes);
	if (file_per_proc)
		set_cache_size(512 * 1024 * 1024 / num_threads);

	char file_name_buf[128];
	char *file_name;

	ProfilerStart(prof_file);
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
			sprintf(file_name_buf, "%s.%08d", file_name_root, data[i].idx);
			file_name = file_name_buf;
			int suggested_nodes[1] = {data->node_id};
			ssd_file_io_init(file_name, 0, 1, 1, suggested_nodes);
		}
		else {
			file_name = file_name_root;
			ssd_file_io_init(file_name, 0, num_threads, num_nodes, NULL);
		}

		data[i].fd = ssd_open(file_name, data->node_id, 0);
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
	ProfilerStop();
	printf("It takes %ld ms\n", end - start);
}
