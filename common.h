#ifndef __COMMON_H__
#define __COMMON_H__

#include <sys/time.h>
#include <stdlib.h>

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

class thread_private;
extern thread_private *get_thread(int idx);

/**
 * Check if the integer is a power of two.
 */
bool align_check(size_t alignment);

#endif
