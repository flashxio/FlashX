#ifndef __COMMON_H__
#define __COMMON_H__

#include <sys/time.h>
#include <stdlib.h>
#include <stdio.h>
#include <sys/types.h>
#include <sys/syscall.h>

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

		for (int i = num - 1; i >= 1; i--) {
			int j = random() % i;
			off_t tmp = offset[j];
			offset[j] = offset[i];
			offset[i] = tmp;
		}
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
public:
	sys_parameters() {
		RAID_block_size = get_default_RAID_block_size();
		SA_min_cell_size = get_default_SA_min_cell_size();
	}

	void init(int RAID_block_size, int SA_min_cell_size) {
		this->RAID_block_size = RAID_block_size;
		this->SA_min_cell_size = SA_min_cell_size;
	}

	int get_RAID_block_size() {
		return RAID_block_size;
	}

	int get_SA_min_cell_size() {
		return SA_min_cell_size;
	}

	static int get_default_SA_min_cell_size() {
		return 12;
	}

	static int get_default_RAID_block_size() {
		return 16;
	}
};

extern sys_parameters params;

#endif
