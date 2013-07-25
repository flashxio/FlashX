#ifndef __COMMON_C_H__
#define __COMMON_C_H__

#include <sys/time.h>
#include <stdlib.h>
#include <stdio.h>
#include <sys/types.h>
#include <sys/syscall.h>
#include <numa.h>
#include <assert.h>
#include <execinfo.h>

#define gettid() syscall(__NR_gettid)

#define ROUND(off, base) (((long) off) & (~((long) (base) - 1)))
#define ROUNDUP(off, base) (((long) off + (base) - 1) & (~((long) (base) - 1)))

#define ROUND_PAGE(off) (((long) off) & (~((long) PAGE_SIZE - 1)))
#define ROUNDUP_PAGE(off) (((long) off + PAGE_SIZE - 1) & (~((long) PAGE_SIZE - 1)))

#define PRINT_BACKTRACE()							\
	do {											\
		void *buf[100];								\
		char **strings;								\
		int nptrs = backtrace(buf, 100);			\
		strings = backtrace_symbols(buf, nptrs);	\
		if (strings == NULL) {						\
			perror("backtrace_symbols");			\
			exit(EXIT_FAILURE);						\
		}											\
		for (int i = 0; i < nptrs; i++)				\
			printf("%s\n", strings[i]);				\
		free(strings);								\
	} while (0)

#define ASSERT_EQ(x, y)								\
	if ((x) != (y))	{								\
		PRINT_BACKTRACE();							\
		assert(x == y);								\
	}

#define ASSERT_TRUE(x)								\
	if (!(x)) {										\
		PRINT_BACKTRACE();							\
		assert(x);									\
	}

enum {
	READ,
	WRITE
};

#ifdef __cplusplus
extern "C" {
#endif

inline float time_diff(struct timeval time1, struct timeval time2)
{
	return time2.tv_sec - time1.tv_sec
			+ ((float)(time2.tv_usec - time1.tv_usec))/1000000;
}

inline long time_diff_us(struct timeval time1, struct timeval time2)
{
	return (time2.tv_sec - time1.tv_sec) * 1000000
			+ (time2.tv_usec - time1.tv_usec);
}

inline static int min(int v1, int v2)
{
	return v1 > v2 ? v2 : v1;
}

inline static int max(int v1, int v2)
{
	return v1 < v2 ? v2 : v1;
}

inline static long get_curr_ms()
{
	struct timeval tv;
	gettimeofday(&tv, NULL);
	return ((long) tv.tv_sec) * 1000 + tv.tv_usec / 1000;
};

inline static long get_curr_us()
{
	struct timeval tv;
	gettimeofday(&tv, NULL);
	return ((long) tv.tv_sec) * 1000 * 1000 + tv.tv_usec;
}

const int CONST_A = 27644437;
const long CONST_P = 68718952447L;

static inline int universal_hash(off_t v, int modulo)
{
	return (v * CONST_A) % CONST_P % modulo;
}

ssize_t get_file_size(const char *file_name);

void permute_offsets(int num, int repeats, int stride, off_t start,
		off_t offsets[]);

/**
 * This returns the first node id where the process can allocate memory.
 */
static inline int numa_get_mem_node()
{
	struct bitmask *bmp = numa_get_membind();
	int nbytes = numa_bitmask_nbytes(bmp);
	int num_nodes = 0;
	int node_id = -1;
	int i;
	for (i = 0; i < nbytes * 8; i++)
		if (numa_bitmask_isbitset(bmp, i)) {
			num_nodes++;
			node_id = i;
		}
	assert(num_nodes == 1);
	return node_id;
}

static inline void bind2node_id(int node_id)
{
	struct bitmask *bmp = numa_allocate_nodemask();
	numa_bitmask_setbit(bmp, node_id);
	numa_bind(bmp);
	numa_free_nodemask(bmp);
}

static inline void bind_mem2node_id(int node_id)
{
	struct bitmask *bmp = numa_allocate_nodemask();
	numa_bitmask_setbit(bmp, node_id);
	numa_set_membind(bmp);
	numa_free_nodemask(bmp);
}


static inline void bind2cpu(int cpu_id)
{
	cpu_set_t cpu_mask;
	CPU_ZERO(&cpu_mask);
	CPU_SET(cpu_id, &cpu_mask);
	int len = sizeof(cpu_mask);
	if (sched_setaffinity(gettid(), len, &cpu_mask) < 0) {
		perror("sched_setaffinity");
		exit(1);
	}
}

#ifdef __cplusplus
}
#endif

#endif
