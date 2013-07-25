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
int num_offs;
int num_threads = 2;
char *file_name = "../conf/data_files.txt";
char *prof_file = "test_c_interface.prof";
long data_size = 1024L * 10240 * 4096;
int block_size = 4096 * 4;
int node_id = 1;
int sync = 1;
int access = WRITE;

struct thread_data
{
	pthread_t tid;
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
	int fd = ssd_open(file_name, 0);
	ssd_set_callback(fd, cb_func);

	struct thread_data *data = arg;
	int num = num_offs / num_threads;
	int off_start = data->idx * num;
	int i;
	char buffer[block_size];

	for (i = 0; i < num; i++) {
		off_t offset = offs[off_start + i];
		if (access == READ)
			ssd_read(fd, (void *) buffer, block_size, offset);
		else
			ssd_write(fd, (void *) buffer, block_size, offset);
	}

	ssd_close(fd);
	return NULL;
}

void *AsyncThreadWriteOrRead(void *arg)
{
	int fd = ssd_open(file_name, 0);
	ssd_set_callback(fd, cb_func);
	int num_completes = 0;
	struct buf_pool *buf_allocator = create_buf_pool(block_size,
			block_size * 10000, node_id);
	struct buf_pool *cb_allocator = create_buf_pool(sizeof(struct callback_data),
			sizeof(struct callback_data) * 10000, node_id);

	struct thread_data *data = arg;
	int num = num_offs / num_threads;
	int off_start = data->idx * num;
	int i;

	for (i = 0; i < num; i++) {
		char *buffer = alloc_buf(buf_allocator);
		assert(off_start + i < num_offs);
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
	num_offs = data_size / block_size;
	offs = (off_t *) malloc(num_offs * sizeof(off_t));
	for (i = 0; i < num_offs; i++)
		offs[i] = i * (long) block_size;
	rand_permute_array(offs, num_offs);

	signal(SIGINT, int_handler);

	ProfilerStart(prof_file);
	ssd_io_init(file_name, 0, num_threads);
	struct thread_data data[num_threads];
	for (i = 0; i < num_threads; i++) {
		data[i].idx = i;
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
	ProfilerStop();
	printf("complete\n");
}
