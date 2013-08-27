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

/**
 * This program tests memory speed in a NUMA machine.
 * It has 5 test modes:
 * For remote memory copy (use_remote flag is enabled), there are three cases:
 *	memcpy-pull: read data from a remote NUMA node and store it in the
 *		local memory with multiple threads.
 *	memcpy-push: read data from the local memory and store it to a
 *		remote NUMA node.
 *	memcpy-r2r: read data from a NUMA node and store it to another NUMA
 *		node.
 *	memread: read data from memory.
 *	memwrite: write data from memory.
 *
 * There is only one type of local memory copy.
 *
 * To avoid threads from copying data from the same node simultaneously,
 * the order of remote nodes visited by memcpy threads is carefully
 * taken, i.e., a thread always starts to copy data from the NUMA node
 * whose node id is right behind its own node id.
 * For example, the order of nodes visited by a memcpy thread on node 0 is
 * 1, 2, 3; the order by a thread on node 1 is 2, 3, 0; etc.
 */

#define NUM_NODES 4
#define NUM_THREADS 8
#define NUM_COPY 16
#define PAGE_SIZE 4096
#define ENTRY_SIZE PAGE_SIZE
#define ARRAY_SIZE (128 * 1024 * 1024)

int use_remote;

enum
{
	MEMCPY_PULL,
	MEMCPY_PUSH,
	MEMCPY_R2R,
	MEMREAD,
	MEMWRITE,
};

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
	int mode;

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
			/*
			 * For remote memory access, NUM_NODES * NUM_THREADS pieces of
			 * memory are allocated, even though only (NUM_NODES - 1) * NUM_THREADS
			 * pieces of memory are actually used.
			 * For local memory access, only NUM_THREADS pieces of memory
			 * are allocated.
			 */
			if (/*(i == data->node_id && use_remote)
					||*/ (i != data->node_id && !use_remote)) {
				init_buffer(&data->src_bufs[i][j]);
				init_buffer(&data->local_bufs[i][j]);
			}
			if ((/*i != data->node_id && */use_remote)
					|| (i == data->node_id && !use_remote)) {
				char *buf;
				
				if (data->mode == MEMCPY_PULL || data->mode == MEMCPY_PUSH
						|| data->mode == MEMCPY_R2R || data->mode == MEMREAD) {
					buf = (char *) numa_alloc_onnode(data->buf_size,
							data->node_id);
					materialize_buf(buf, data->buf_size);
					set_buffer(&data->src_bufs[i][j], buf, data->buf_size,
							data->node_id);
				}
				else
					init_buffer(&data->src_bufs[i][j]);

				if (data->mode == MEMCPY_PULL || data->mode == MEMCPY_PUSH
						|| data->mode == MEMCPY_R2R || data->mode == MEMWRITE) {
					buf = (char *) numa_alloc_onnode(data->buf_size,
							data->node_id);
					materialize_buf(buf, data->buf_size);
					set_buffer(&data->local_bufs[i][j], buf, data->buf_size,
							data->node_id);
				}
				else
					init_buffer(&data->local_bufs[i][j]);
			}
		}
	}
	return NULL;
}

struct buf_copy_data
{
	struct copy_entry {
		// A data buffer in a remote NUMA node.
		struct data_buffer from;
		// the data buffer local to the memcpy thread.
		struct data_buffer to;
	} copy_entries[NUM_NODES - 1];
	int node_id;
	int mode;

	size_t copy_size;
};

int sum_buf(char *buf, int size)
{
	long *lbuf = (long *) buf;
	int num = size / sizeof(long);
	int i;
	long sum = 0;

	for (i = 0; i < num; i++)
		sum += lbuf[i];
	return sum;
}

int memread(struct data_buffer buf)
{
	long sum = 0;
	char page[4096];
	char *addr = buf.addr;
	char *end = buf.addr + buf.size;
	for (; addr < end; addr += sizeof(page)) {
		memcpy(page, addr, sizeof(page));
		sum += sum_buf(page, sizeof(page));
	}
	// The optimization of a compiler may remove the code entirely
	// if we don't return the sum.
	return sum;
}

int memwrite(struct data_buffer buf)
{
	char page[4096];
	char *addr = buf.addr;
	char *end = buf.addr + buf.size;
	memset(page, 0, sizeof(page));
	for (; addr < end; addr += sizeof(page))
		memcpy(addr, page, sizeof(page));
	return buf.size;
}

/**
 * Each thread copies memory from a remote memory buffer.
 */
void *buf_copy_func(void *arg)
{
	int i, j;
	size_t size = 0;
	struct buf_copy_data *data = (struct buf_copy_data *) arg;
	long memread_sum = 0;
	int start_entry_idx = 0;

	bind2node_id(data->node_id);
	data->copy_size = 0;
	if (use_remote) {
		/*
		 * We make sure threads on different nodes start to access memory
		 * from different nodes. In the current implementation, a thread
		 * always starts to copy data from the NUMA node whose node id is
		 * right behind its own node id. For example, the order of nodes
		 * visited by a memcpy thread on node 0 is 1, 2, 3; the order by
		 * a thread on node 1 is 2, 3, 0; etc.
		 */
		for (i = 0; i < NUM_NODES - 1; i++) {
			if (data->mode == MEMCPY_PULL || data->mode == MEMCPY_R2R
					|| data->mode == MEMREAD) {
				if (data->copy_entries[i].from.node_id == (data->node_id + 1)
						% NUM_NODES)
					break;
			}
			else if (data->mode == MEMCPY_PUSH || data->mode == MEMWRITE) {
				if (data->copy_entries[i].to.node_id == (data->node_id + 1)
						% NUM_NODES)
					break;
			}
		}
		assert(i != NUM_NODES - 1);
		start_entry_idx = i;
	}
	for (j = 0; j < NUM_COPY; j++)
		for (i = 0; i < NUM_NODES - 1; i++) {
			struct buf_copy_data::copy_entry *entry = &data->copy_entries[(i
					+ start_entry_idx) % (NUM_NODES - 1)];
			int size = entry->to.size;
			if (data->mode == MEMCPY_PULL) {
				if (use_remote) {
					assert(entry->to.node_id == data->node_id);
					assert(entry->from.node_id != data->node_id);
				}
				else {
					assert(entry->to.node_id == data->node_id);
					assert(entry->from.node_id == data->node_id);
				}
				assert(size == entry->from.size);
				memcpy(entry->to.addr, entry->from.addr, size);
			}
			else if (data->mode == MEMCPY_PUSH) {
				if (use_remote) {
					assert(entry->to.node_id != data->node_id);
					assert(entry->from.node_id == data->node_id);
				}
				else {
					assert(entry->to.node_id == data->node_id);
					assert(entry->from.node_id == data->node_id);
				}
				assert(size == entry->from.size);
				memcpy(entry->to.addr, entry->from.addr, size);
			}
			else if (data->mode == MEMCPY_R2R) {
				if (use_remote) {
					assert(entry->to.node_id != data->node_id);
					assert(entry->from.node_id != data->node_id);
				}
				else {
					assert(entry->to.node_id == data->node_id);
					assert(entry->from.node_id == data->node_id);
				}
				assert(size == entry->from.size);
				memcpy(entry->to.addr, entry->from.addr, size);
			}
			else if (data->mode == MEMREAD) {
				assert(entry->to.addr == NULL);
				if (use_remote)
					assert(entry->from.node_id != data->node_id);
				else
					assert(entry->from.node_id == data->node_id);
				memread_sum += memread(entry->from);
			}
			else if (data->mode == MEMWRITE) {
				if (use_remote)
					assert(entry->to.node_id != data->node_id);
				else
					assert(entry->to.node_id == data->node_id);
				assert(entry->from.addr == NULL);
				memwrite(entry->to);
			}
			data->copy_size += size;
		}
	return (void *) memread_sum;
}

int main(int argc, char *argv[])
{
	struct timeval start_time, end_time;
	ssize_t read_bytes = 0;
	/* the number of entries the array can contain. */
	int node;
	int mode;
	struct buf_init_data node_buf_data[NUM_NODES];

	if (argc < 3) {
		fprintf(stderr, "memtest mode remote\n");
		fprintf(stderr, "mode: 0(memcpy_pull), 1(memcpy_push), 2(memcpy_r2r), 3(memread), 4(memwrite)\n");
		fprintf(stderr, "remote: 0(local) 1(remote)\n");
		return -1;
	}
	mode = atoi(argv[1]);
	assert(mode <= MEMWRITE);
	use_remote = atoi(argv[2]);

	{
		pthread_t threads[NUM_NODES];
		int ret;
		int i;

		printf("initialize the buffer for memcpy\n");
		gettimeofday(&start_time, NULL);
		for (i = 0; i < NUM_NODES; i++) {
			node_buf_data[i].buf_size = ARRAY_SIZE;
			node_buf_data[i].node_id = i;
			node_buf_data[i].mode = mode;

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
				copy_data[local_node_id][thread_id].mode = mode;
				for (remote_node_id = 0; remote_node_id < NUM_NODES;
						remote_node_id++) {
					struct data_buffer buf1, buf2;
					int node1, node2;

					if (use_remote) {
						if (local_node_id == remote_node_id)
							continue;
						if (mode == MEMCPY_R2R) {
							node1 = (local_node_id + 1) % NUM_NODES;
							node2 = remote_node_id;
						}
						else if (mode == MEMCPY_PULL || mode == MEMREAD) {
							node1 = remote_node_id;
							node2 = local_node_id;
						}
						else if (mode == MEMCPY_PUSH || mode == MEMWRITE) {
							node1 = local_node_id;
							node2 = remote_node_id;
						}
						else
							assert(0);
					}
					else {
						// We only need to run for (NUM_NODES -1) times, and
						// each time we copy data from the same piece of memory.
						if (remote_node_id == NUM_NODES - 1)
							break;
						node1 = local_node_id;
						node2 = local_node_id;
					}

					if (mode == MEMCPY_PULL || mode == MEMCPY_PUSH
							|| mode == MEMCPY_R2R) {
						buf1 = node_buf_data[node1].src_bufs[node2][thread_id];
						assert(is_valid_buffer(&buf1));
						copy_data[local_node_id][thread_id].copy_entries[idx].from = buf1;

						buf2 = node_buf_data[node2].local_bufs[node1][thread_id];
						assert(is_valid_buffer(&buf2));
						copy_data[local_node_id][thread_id].copy_entries[idx].to = buf2;
					}
					else if (mode == MEMREAD) {
						buf1 = node_buf_data[node1].src_bufs[node2][thread_id];
						assert(is_valid_buffer(&buf1));
						copy_data[local_node_id][thread_id].copy_entries[idx].from = buf1;
						init_buffer(&copy_data[local_node_id][
								thread_id].copy_entries[idx].to);
					}
					else if (mode == MEMWRITE) {
						buf2 = node_buf_data[node2].local_bufs[node1][thread_id];
						assert(is_valid_buffer(&buf2));
						copy_data[local_node_id][thread_id].copy_entries[idx].to = buf2;
						init_buffer(&copy_data[local_node_id][
								thread_id].copy_entries[idx].from);
					}
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
