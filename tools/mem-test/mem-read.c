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
#include <sys/types.h>
#include <sys/syscall.h>
#include <pthread.h>
//#include <hugetlbfs.h>
//#include <emmintrin.h>

#define gettid() syscall(__NR_gettid)

#define NUM_THREADS 64
#define PAGE_SIZE 4096
#define CACHE_LINE 64
#define ENTRY_SIZE 32
#define NUM_PROCESSORS 4
#define NUM_CORES (NUM_PROCESSORS * 8)
//#define ARRAY_SIZE (1024 * 1024 * 1024)

unsigned int nentries;
int nthreads;
struct timeval global_start;
char *array;
int read_size;
int node_id;

long sum1(char *array, int start_i, int end_i);
long sum2(char *array, int start_i, int end_i);

float time_diff(struct timeval time1, struct timeval time2)
{
	return time2.tv_sec - time1.tv_sec
			+ ((float)(time2.tv_usec - time1.tv_usec))/1000000;
}

static inline void bind2node_id(int node_id)
{
	struct bitmask *bmp = numa_allocate_nodemask();
	numa_bitmask_setbit(bmp, node_id);
	numa_bind(bmp);
	numa_free_nodemask(bmp);
}

void flush_cache(void *start, void *end)
{
//	void *p = start;
//	while (p < end) {
//		_mm_clflush(p);
//		p += CACHE_LINE;
//	}
}

void *rand_read(void *arg)
{
	int fd;
	ssize_t ret;
	int i, j, start_i, end_i;
	struct timeval start_time, end_time;
	long res;

	int thread_id = (long) arg;
	start_i = nentries / nthreads * thread_id;
	end_i = start_i + nentries / nthreads;
	assert(((long) array) % CACHE_LINE == 0);
	int cpu_id = (thread_id * NUM_PROCESSORS + node_id) % NUM_CORES;
	cpu_set_t cpu_mask;
	CPU_ZERO(&cpu_mask);
	CPU_SET(cpu_id, &cpu_mask);
	int len = sizeof(cpu_mask);
	if (sched_setaffinity(gettid(), len, &cpu_mask) < 0) {
		perror("sched_setaffinity");
		exit(1);
	}

	gettimeofday(&start_time, NULL);
	char local_page[PAGE_SIZE];
	for (j = 0; j < 8; j++) {
		struct timeval start_time1, end_time1;
		flush_cache(array + start_i, array + end_i);
//		gettimeofday(&start_time1, NULL);
		res = sum1(array, start_i, end_i);
//		gettimeofday(&end_time1, NULL);
//		printf("iteration %d from %d to %d: takes %f seconds, res: %ld\n",
//				0, start_i, end_i, time_diff(start_time1, end_time1), res);

		flush_cache(array + start_i, array + end_i);
//		gettimeofday(&start_time1, NULL);
		res += sum2(array, start_i, end_i);
//		gettimeofday(&end_time1, NULL);
//		printf("iteration %d from %d to %d: takes %f seconds, res: %ld\n",
//				1, start_i, end_i, time_diff(start_time1, end_time1), res);
	}
	gettimeofday(&end_time, NULL);
	printf("add %d cache lines, res: %ld, start at %f seconds, takes %f seconds\n",
			(end_i - start_i) * 2 * 8, res, time_diff(global_start, start_time),
			time_diff(start_time, end_time));
	return NULL;
}

int main(int argc, char *argv[])
{
	int ret;
	int i;
	struct timeval start_time, end_time;
	ssize_t read_bytes = 0;
	pthread_t threads[NUM_THREADS];

	long pagesize = 1024 * 1024 * 2;
//	long pagesize = gethugepagesize();
	if (argc != 3) {
		fprintf(stderr, "read node_id num_threads\n");
		exit(1);
	}
	printf("huge page size: %ld\n", pagesize);
	int array_size = pagesize * 512;

	gettimeofday(&start_time, NULL);
	global_start = start_time;

	nentries = array_size / ENTRY_SIZE;
	node_id = atoi(argv[1]);
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

	printf("process %d\n", getpid());
	numa_set_bind_policy(1);
	bind2node_id(0);
	array = (char *) numa_alloc_onnode(array_size, 0);
//	array = get_huge_pages(array_size, GHP_DEFAULT);
	if (array == NULL) {
		perror("get_huge_pages");
		exit(1);
	}
	printf("array: %p\n", array);
	/* we need to avoid the cost of page fault. */
	for (i = 0; i < array_size; i += ENTRY_SIZE)
		*(long *) (array + i) = i;
	printf("allocate source array %p in node 0\n", array);

	printf("run on node %d\n", node_id);
	bind2node_id(node_id);

	for (i = 0; i < nthreads; i++) {
		ret = pthread_create(&threads[i], NULL,
				rand_read, (void *) (long) i);
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

//	free_huge_pages(array);
	gettimeofday(&end_time, NULL);
	printf("takes %f seconds\n",
			end_time.tv_sec - start_time.tv_sec
			+ ((float)(end_time.tv_usec - start_time.tv_usec))/1000000);
}
