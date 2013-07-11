#ifndef __COMMON_H__
#define __COMMON_H__

#include <sys/time.h>
#include <stdlib.h>
#include <stdio.h>
#include <sys/types.h>
#include <sys/syscall.h>
#include <numa.h>
#include <assert.h>
#include <execinfo.h>

#include <string>
#include <vector>

#include "parameters.h"

#define gettid() syscall(__NR_gettid)

enum {
	READ,
	WRITE
};

extern int nthreads;
extern struct timeval global_start;

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

class thread_private;
extern thread_private *get_thread(int idx);

/**
 * Check if the integer is a power of two.
 */
bool align_check(size_t alignment);

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

struct file_info
{
	std::string name;
	// The NUMA node id where the disk is connected to.
	int node_id;
};

static inline std::string itoa(int n)
{
	char buf[32];
	snprintf(buf, sizeof(buf), "%d", n);
	return buf;
}

static inline int universal_hash(off_t v, int modulo)
{
	return (v * CONST_A) % CONST_P % modulo;
}

int retrieve_data_files(std::string file_file,
		std::vector<file_info> &data_files);
ssize_t get_file_size(const char *file_name);

#define ROUND(off, base) (((long) off) & (~((long) (base) - 1)))
#define ROUNDUP(off, base) (((long) off + (base) - 1) & (~((long) (base) - 1)))

#define ROUND_PAGE(off) (((long) off) & (~((long) PAGE_SIZE - 1)))
#define ROUNDUP_PAGE(off) (((long) off + PAGE_SIZE - 1) & (~((long) PAGE_SIZE - 1)))

void permute_offsets(int num, int repeats, int stride, off_t start,
		off_t offsets[]);

template<class T>
void rand_permute_array(T arr[], int num)
{
	for (int i = num - 1; i >= 1; i--) {
		int j = random() % i;
		T tmp = arr[j];
		arr[j] = arr[i];
		arr[i] = tmp;
	}
}

class rand_permute
{
	off_t *offset;
	long num;

public:
	/**
	 * @start: the index of the first entry.
	 */
	rand_permute(long num, int stride, long start) {
		offset = new off_t[num];
		for (int i = 0; i < num; i++) {
			offset[i] = ((off_t) i) * stride + start * stride;
		}

		rand_permute_array(offset, num);
	}

	~rand_permute() {
		delete [] offset;
	}

	off_t get_offset(long idx) const {
		return offset[idx];
	}
};

long str2size(std::string str);

class sys_parameters
{
	int RAID_block_size;
	int SA_min_cell_size;
	int test_hit_rate;
public:
	sys_parameters() {
		RAID_block_size = get_default_RAID_block_size();
		SA_min_cell_size = get_default_SA_min_cell_size();
	}

	void init(int RAID_block_size, int SA_min_cell_size, int test_hit_rate) {
		this->RAID_block_size = RAID_block_size;
		this->SA_min_cell_size = SA_min_cell_size;
		this->test_hit_rate = test_hit_rate;
	}

	int get_RAID_block_size() {
		return RAID_block_size;
	}

	int get_SA_min_cell_size() {
		return SA_min_cell_size;
	}

	int get_test_hit_rate() {
		return test_hit_rate;
	}

	static int get_default_SA_min_cell_size() {
		return 12;
	}

	static int get_default_RAID_block_size() {
		return 16;
	}
};

/**
 * This returns the first node id where the process can allocate memory.
 */
static inline int numa_get_mem_node()
{
	struct bitmask *bmp = numa_get_membind();
	int nbytes = numa_bitmask_nbytes(bmp);
	int num_nodes = 0;
	int node_id = -1;
	for (int i = 0; i < nbytes * 8; i++)
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

extern sys_parameters params;

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

#endif
