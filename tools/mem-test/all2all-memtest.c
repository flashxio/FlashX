#include <stdio.h>
#include <unistd.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <sys/time.h>
#include <stdlib.h>
#include <sys/resource.h>
#include <string.h>
#include <numa.h>
#include <numaif.h>
#include <pthread.h>
#include <assert.h>

#define NUM_NODES 4
#define NUM_THREADS 8
#define NUM_COPY 16
#define PAGE_SIZE 4096
#define ENTRY_SIZE PAGE_SIZE
#define ARRAY_SIZE (128 * 1024 * 1024)

struct timeval global_start;

static inline void bind2node_id(int node_id)
{
	struct bitmask *bmp = numa_allocate_nodemask();
	numa_bitmask_setbit(bmp, node_id);
	numa_bind(bmp);
	numa_free_nodemask(bmp);
}

float time_diff(struct timeval time1, struct timeval time2)
{
	return time2.tv_sec - time1.tv_sec
			+ ((float)(time2.tv_usec - time1.tv_usec))/1000000;
}

/**
 * Make sure the pages are allocated in the memory.
 */
void materialize_buf(char *buf, int size)
{
	int j;

	assert(((long) buf) % PAGE_SIZE == 0);
	assert(size % PAGE_SIZE == 0);
	for (j = 0; j < size; j += PAGE_SIZE)
		buf[j] = 0;
}

struct data_buffer
{
	char *addr;
	int size;
	int node_id;
};

void init_buffer(struct data_buffer *buf)
{
	buf->addr = NULL;
	buf->size = 0;
	buf->node_id = -1;
}

void set_buffer(struct data_buffer *buf, char *addr, int size, int node_id)
{
	buf->addr = addr;
	buf->size = size;
	buf->node_id = node_id;
}

int is_valid_buffer(struct data_buffer *buf)
{
	return buf->addr != NULL;
}

struct buf_init_data
{
	int buf_size;
	int node_id;

	// The source data buffers for other threads to copy data from.
	struct data_buffer src_bufs[NUM_NODES][NUM_THREADS];
	// The local data buffers for memcpy threads to store data to.
	struct data_buffer local_bufs[NUM_NODES][NUM_THREADS];
};

/**
 * Each thread initializes the buffers in the local NUMA node for memory copy
 * in the next phase.
 */
void *buf_init_func(void *arg)
{
	int i, j;
	struct buf_init_data *data = (struct buf_init_data *) arg;

	bind2node_id(data->node_id);
	for (i = 0; i < NUM_NODES; i++) {
		for (j = 0; j < NUM_THREADS; j++) {
			if (i == data->node_id) {
				init_buffer(&data->src_bufs[i][j]);
				init_buffer(&data->local_bufs[i][j]);
			}
			else {
				char *buf;
				
				buf = (char *) numa_alloc_onnode(data->buf_size, data->node_id);
				materialize_buf(buf, data->buf_size);
				set_buffer(&data->src_bufs[i][j], buf, data->buf_size,
						data->node_id);

				buf = (char *) numa_alloc_onnode(data->buf_size, data->node_id);
				materialize_buf(buf, data->buf_size);
				set_buffer(&data->local_bufs[i][j], buf, data->buf_size,
						data->node_id);
			}
		}
	}
	return NULL;
}

struct buf_copy_data
{
	struct {
		// A data buffer in a remote NUMA node.
		struct data_buffer from;
		// the data buffer local to the memcpy thread.
		struct data_buffer to;
	} copy_entries[NUM_NODES - 1];
	int node_id;

	size_t copy_size;
};

/**
 * Each thread copies memory from a remote memory buffer.
 */
void *buf_copy_func(void *arg)
{
	int i, j;
	size_t size = 0;
	struct buf_copy_data *data = (struct buf_copy_data *) arg;

	bind2node_id(data->node_id);
	data->copy_size = 0;
	for (j = 0; j < NUM_COPY; j++)
		for (i = 0; i < NUM_NODES - 1; i++) {
			int size = data->copy_entries[i].to.size;
			assert(size == data->copy_entries[i].from.size);
			assert(data->copy_entries[i].to.node_id == data->node_id);
			assert(data->copy_entries[i].from.node_id != data->node_id);
			memcpy(data->copy_entries[i].to.addr, data->copy_entries[i].from.addr,
					size);
			data->copy_size += size;
		}
	return NULL;
}

int main(int argc, char *argv[])
{
	struct timeval start_time, end_time;
	ssize_t read_bytes = 0;
	/* the number of entries the array can contain. */
	int node;
	struct buf_init_data node_buf_data[NUM_NODES];

	{
		pthread_t threads[NUM_NODES];
		int ret;
		int i;

		printf("initialize the buffer for memcpy\n");
		gettimeofday(&start_time, NULL);
		for (i = 0; i < NUM_NODES; i++) {
			node_buf_data[i].buf_size = ARRAY_SIZE;
			node_buf_data[i].node_id = i;

			ret = pthread_create(&threads[i], NULL,
					buf_init_func, &node_buf_data[i]);
			if (ret) {
				perror("pthread_create");
				exit(1);
			}
		}

		for (i = 0; i < NUM_NODES; i++) {
			ret = pthread_join(threads[i], NULL);
			if (ret) {
				perror("pthread_join");
				exit(1);
			}
		}
		gettimeofday(&end_time, NULL);
		printf("buffer initialization: done. It takes %fs\n",
				time_diff(start_time, end_time));
	}

	{
		int ret;
		int local_node_id, thread_id, remote_node_id;
		pthread_t threads[NUM_NODES][NUM_THREADS];
		struct buf_copy_data copy_data[NUM_NODES][NUM_THREADS];
		size_t tot_size = 0;

		gettimeofday(&start_time, NULL);
		global_start = start_time;
		for (local_node_id = 0; local_node_id < NUM_NODES; local_node_id++) {
			for (thread_id = 0; thread_id < NUM_THREADS; thread_id++) {
				int idx = 0;
				copy_data[local_node_id][thread_id].node_id = local_node_id;
				for (remote_node_id = 0; remote_node_id < NUM_NODES;
						remote_node_id++) {
					struct data_buffer buf1, buf2;
					if (local_node_id == remote_node_id)
						continue;

					buf1 = node_buf_data[remote_node_id].src_bufs[local_node_id][thread_id];
					assert(is_valid_buffer(&buf1));
					copy_data[local_node_id][thread_id].copy_entries[idx].from = buf1;

					buf2 = node_buf_data[local_node_id].local_bufs[remote_node_id][thread_id];
					assert(is_valid_buffer(&buf2));
					copy_data[local_node_id][thread_id].copy_entries[idx].to = buf2;
					printf("%p\t%p\n", buf1.addr, buf2.addr);
					idx++;
				}

				ret = pthread_create(&threads[local_node_id][thread_id], NULL,
						buf_copy_func, &copy_data[local_node_id][thread_id]);
				if (ret) {
					perror("pthread_create");
					exit(1);
				}
			}
		}

		for (local_node_id = 0; local_node_id < NUM_NODES; local_node_id++) {
			for (thread_id = 0; thread_id < NUM_THREADS; thread_id++) {
				ret = pthread_join(threads[local_node_id][thread_id], NULL);
				if (ret) {
					perror("pthread_join");
					exit(1);
				}
				tot_size += copy_data[local_node_id][thread_id].copy_size;
			}
		}
		gettimeofday(&end_time, NULL);
		printf("memcpy %ld bytes, it takes %f seconds\n", tot_size,
				time_diff(start_time, end_time));
	}
}
