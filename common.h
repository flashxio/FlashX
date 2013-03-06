#ifndef __COMMON_H__
#define __COMMON_H__

#include <sys/time.h>
#include <stdlib.h>
#include <stdio.h>
#include <sys/types.h>
#include <sys/syscall.h>

#include <string>

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

#endif
