#define _GNU_SOURCE
#include <stdio.h>
#include <unistd.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <sys/resource.h>
#include <fcntl.h>
#include <sys/time.h>
#include <stdlib.h>
#include <string.h>

#define NUM 0xFFFFFFFL
#define NUM_THREADS 8

/**
 * this program tries to compare the performance difference
 * between updating global variables and CPU-local variables.
 * Updating global variables cause a lot of cache invalidation
 * among CPUs, which can cause significant performance overhead.
 */

volatile long global_counter;

void thread_worker(void *arg)
{
	int i;
	volatile long *counter = arg;

//#ifdef USE_GLOBAL
//	printf("thread_worker starts (using global variable)\n");
//#else
//	printf("thread_worker starts (using thread-private variable)\n");
//#endif
	for (i = 0; i < NUM; i++)
#ifdef USE_GLOBAL
//		global_counter++;
		__sync_fetch_and_add(&global_counter, 1);
#else
		(*counter)++;
#endif
//#ifdef USE_GLOBAL
//	pthread_exit((void *) global_counter);
//#else
//	pthread_exit((void *) *counter);
//#endif
}

int main(int argc, char *argv[])
{
	int ret;
	int i;
	struct timeval start_time, end_time;
	ssize_t read_bytes = 0;
	/* so each counter is in a separate cache line */
	long counters[NUM_THREADS * 16];
	pthread_t threads[NUM_THREADS];
	int nthreads;

	if (argc != 2) {
		fprintf(stderr, "%s num_threads\n", argv[0]);
		exit(1);
	}

	nthreads = atoi(argv[1]);
	if (nthreads > NUM_THREADS) {
		fprintf(stderr, "too many threads\n");
		exit(1);
	}

	ret = setpriority(PRIO_PROCESS, getpid(), -20);
	if (ret < 0) {
		perror("setpriority");
		exit(1);
	}

	gettimeofday(&start_time, NULL);
	
	memset(counters, 0, sizeof(counters));
	for (i = 0; i < nthreads; i++) {
		ret = pthread_create(&threads[i], NULL,
				thread_worker, &counters[i * 16]);
		if (ret) {
			perror("pthread_create");
			exit(1);
		}
	}

	for (i = 0; i < nthreads; i++) {
		ssize_t size;
		ret = pthread_join(threads[i], (void **) &size);
		if (ret) {
			perror("pthread_join");
			exit(1);
		}
		read_bytes += size;
	}
	gettimeofday(&end_time, NULL);
	printf("read %ld bytes, takes %f seconds\n",
			read_bytes, end_time.tv_sec - start_time.tv_sec
			+ ((float)(end_time.tv_usec - start_time.tv_usec))/1000000);
}
