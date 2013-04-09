#include <stdio.h>
#include <assert.h>
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
#include <google/profiler.h>

#define NUM_THREADS 64
#define PAGE_SIZE 4096
#define ENTRY_SIZE 32
#define ARRAY_SIZE 1073741824

unsigned int nentries;
int nthreads;
struct timeval global_start;
char *array;
int read_size;

float time_diff(struct timeval time1, struct timeval time2)
{
	return time2.tv_sec - time1.tv_sec
			+ ((float)(time2.tv_usec - time1.tv_usec))/1000000;
}

void rand_read(void *arg)
{
	int fd;
	ssize_t ret;
	int i, j, start_i, end_i;
	struct timeval start_time, end_time;
	long res = 0;

	start_i = (long) arg;
	end_i = start_i + nentries / nthreads;
	gettimeofday(&start_time, NULL);
	int num_entries = 0;
	for (j = 0; j < 8; j++) {
		struct timeval start_time1, end_time1;
		gettimeofday(&start_time1, NULL);
		for (i = start_i; i < end_i; i++) {
			res += *((long *) (array + i * ENTRY_SIZE));
			res += *((long *) (array + i * ENTRY_SIZE + 8));
			res += *((long *) (array + i * ENTRY_SIZE + 16));
			res += *((long *) (array + i * ENTRY_SIZE + 24));
			num_entries++;
		}
		gettimeofday(&end_time1, NULL);
		printf("iteration %d: takes %f seconds\n", j, time_diff(start_time1, end_time1));
	}
	gettimeofday(&end_time, NULL);
	printf("add %d entries from %p to %p, res: %ld, start at %f seconds, takes %f seconds\n",
			num_entries, array + start_i * ENTRY_SIZE, array + end_i * ENTRY_SIZE,
			res, time_diff(global_start, start_time),
			time_diff(start_time, end_time));
}

int main(int argc, char *argv[])
{
	int ret;
	int i;
	struct timeval start_time, end_time;
	ssize_t read_bytes = 0;
	pthread_t threads[NUM_THREADS];
	/* the number of entries the array can contain. */
	int node;

	if (argc != 3) {
		fprintf(stderr, "read node_id num_threads\n");
		exit(1);
	}

	nentries = ARRAY_SIZE / ENTRY_SIZE;
	node = atoi(argv[1]);
	nthreads = atoi(argv[2]);
	if (nthreads > NUM_THREADS) {
		fprintf(stderr, "too many threads\n");
		exit(1);
	}

#if 0
	int ncpus = numa_num_configured_cpus();
	printf("there are %d cores in the machine\n", ncpus);
	for (i = 0; i < ncpus; i++) {
		printf("cpu %d belongs to node %d\n",
			i, numa_node_of_cpu(i));
	}
#endif

	array = numa_alloc_onnode(ARRAY_SIZE, 0);
	/* we need to avoid the cost of page fault. */
	for (i = 0; i < ARRAY_SIZE; i += PAGE_SIZE)
		array[i] = 0;
	printf("allocate source array %p in node 0\n", array);

	printf("run on node %d\n", node);
	if (numa_run_on_node(node) < 0) {
		perror("numa_run_on_node");
		exit(1);
	}

//	ret = setpriority(PRIO_PROCESS, getpid(), -20);
//	if (ret < 0) {
//		perror("setpriority");
//		exit(1);
//	}

	ProfilerStart("rand-memcpy");
	gettimeofday(&start_time, NULL);
	global_start = start_time;
	for (i = 0; i < nthreads; i++) {
		ret = pthread_create(&threads[i], NULL,
				rand_read, (void *) (long) (nentries / nthreads * i));
		if (ret) {
			perror("pthread_create");
			exit(1);
		}
	}

	for (i = 0; i < nthreads; i++) {
		ssize_t size;
		ret = pthread_join(threads[i], NULL);
		if (ret) {
			perror("pthread_join");
			exit(1);
		}
	}
	gettimeofday(&end_time, NULL);
	ProfilerStop();
	printf("takes %f seconds\n",
			end_time.tv_sec - start_time.tv_sec
			+ ((float)(end_time.tv_usec - start_time.tv_usec))/1000000);
}
