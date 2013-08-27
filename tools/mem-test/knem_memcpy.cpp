#include <stdio.h>
#include <fcntl.h>
#include <errno.h>
#include <assert.h>
#include <unistd.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <sched.h>
#include <getopt.h>
#include <sys/ioctl.h>
#include <sys/mman.h>

#include <string>

#include "knem_io.h"
#include "common.h"

#define LEN (1024 * 1024 * 32)
#define NUM_THREADS 1
#define NUM_PROCESSORS 4

int node0 = 0;
int node1 = 0;

static knem_cookie_t
create_region(int fd, int write, struct knem_cmd_param_iovec *iovec, int nr)
{
	int err;
	struct knem_cmd_create_region create;
	create.iovec_array = (uintptr_t) iovec;
	create.iovec_nr = nr;
	create.flags = KNEM_FLAG_SINGLEUSE;
	create.protection = write ? PROT_WRITE : PROT_READ;
	err = ioctl(fd, KNEM_CMD_CREATE_REGION, &create);
	if (err < 0) {
		perror("ioctl (create region)");
		exit(1);
	}
	return create.cookie;
}

static void
destroy_region(int fd, knem_cookie_t cookie)
{
	int err = ioctl(fd, KNEM_CMD_DESTROY_REGION, &cookie);
	assert(!err);
}

static void
check_status(knem_status_t status, const std::string &msg)
{
	if (status == KNEM_STATUS_SUCCESS) {
//		printf("%ld: %s OK\n", gettid(), msg.c_str());
	}
	else if (status == KNEM_STATUS_FAILED) {
//		printf("%ld: %s Failed\n", gettid(), msg.c_str());
		assert(0);
	} else {
//		printf("%ld: %s Pending?!\n", gettid(), msg.c_str());
		assert(0);
	}
}

static void
copy(int fd, knem_cookie_t src_cookie, knem_cookie_t dst_cookie, long length,
		unsigned long flags, volatile knem_status_t *sarray, int sindex)
{
	int err;
	struct knem_cmd_copy_bounded copy;
	copy.src_cookie = src_cookie;
	copy.src_offset = 0;
	copy.dst_cookie = dst_cookie;
	copy.dst_offset = 0;
	copy.length = length;
	copy.async_status_index = sindex;
	copy.flags = flags;
	err = ioctl(fd, KNEM_CMD_COPY_BOUNDED, &copy);
	if (err < 0) {
		perror("ioctl (copy)");
		exit(1);
	}
	if (copy.current_status != KNEM_STATUS_PENDING) {
		check_status(copy.current_status, "sync ");
	} else {
		if (!sarray) {
			fprintf(stderr, "no async array given\n");
			exit(1);
		}
		while (sarray[sindex] == KNEM_STATUS_PENDING);
		check_status(sarray[sindex], "async ");
	}
}

void check_mem(char *arr, int len)
{
	for (int i = 0; i < len; i++)
		assert(arr[i] == 1);
}

void *dma_memcpy(void *arg)
{
	struct knem_cmd_param_iovec send_iovec[2], recv_iovec[2];
	knem_cookie_t cookie, cookie2;
	volatile knem_status_t *sarray;
	int fd;
	int i;
	int thread_id = (long) arg;
	int cpu_id = NUM_PROCESSORS * thread_id + node0;

	bind2cpu(cpu_id);
	printf("bind to CPU %d\n", cpu_id);

	bind2node_id(node0);
	send_iovec[0].base = (uintptr_t) numa_alloc_onnode(LEN, node0);;
	assert(send_iovec[0].base);
	send_iovec[0].len = LEN;
	send_iovec[1].base = (uintptr_t) numa_alloc_onnode(LEN, node0);;
	assert(send_iovec[1].base);
	send_iovec[1].len = LEN;

	bind2node_id(node1);
	recv_iovec[0].base = (uintptr_t) numa_alloc_onnode(LEN, node1);;
	assert(recv_iovec[0].base);
	recv_iovec[0].len = LEN;
	recv_iovec[1].base = (uintptr_t) numa_alloc_onnode(LEN, node1);;
	assert(recv_iovec[1].base);
	recv_iovec[1].len = LEN;

	bind2node_id(node0);
	fd = open("/dev/knem", O_RDWR);
	if (fd < 0) {
		perror("open");
		exit(1);
	}

	sarray = (knem_status_t *) mmap(NULL, sizeof(knem_status_t), PROT_READ|PROT_WRITE,
			MAP_SHARED, fd, KNEM_STATUS_ARRAY_FILE_OFFSET);
	if (sarray == MAP_FAILED) {
		perror("mmap status array");
		exit(1);
	}

	memset((void *) send_iovec[0].base, 1, LEN);
	memset((void *) send_iovec[1].base, 1, LEN);
	long duration = 0;
	long length = 0;
	long start = get_curr_ms();
	for (i = 0; i < 1024; i++) {
		memset((void *) recv_iovec[0].base, 0, LEN);
		memset((void *) recv_iovec[1].base, 0, LEN);
		cookie = create_region(fd, 0, send_iovec, 2);
		cookie2 = create_region(fd, 1, recv_iovec, 2);
		copy(fd, cookie, cookie2, LEN * 2, KNEM_FLAG_DMA
				/*| KNEM_FLAG_ASYNCDMACOMPLETE | KNEM_FLAG_DMATHREAD*/, sarray, 0);
		length += LEN * 2;
		check_mem((char *) recv_iovec[0].base, LEN);
		check_mem((char *) recv_iovec[1].base, LEN);
	}
	long end = get_curr_ms();
	duration += (end - start);
	printf("cpu %d: copying %ld bytes takes %ld ms\n", cpu_id, length, duration);
}

int main()
{
	int i;
	int ret;
	pthread_t threads[NUM_THREADS];

	for (i = 0; i < NUM_THREADS; i++) {
		ret = pthread_create(&threads[i], NULL,
				dma_memcpy, (void *) (long) i);
		if (ret) {
			perror("pthread_create");
			exit(1);
		}
	}

	for (i = 0; i < NUM_THREADS; i++) {
		ssize_t size;
		ret = pthread_join(threads[i], NULL);
		if (ret) {
			perror("pthread_join");
			exit(1);
		}
	}
}
