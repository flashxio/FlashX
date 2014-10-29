#ifndef __COMMON_C_H__
#define __COMMON_C_H__

/**
 * Copyright 2014 Open Connectome Project (http://openconnecto.me)
 * Written by Da Zheng (zhengda1936@gmail.com)
 *
 * This file is part of SAFSlib.
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

#include <errno.h>
#include <sys/time.h>
#include <stdlib.h>
#include <stdio.h>
#include <sys/types.h>
#include <sys/syscall.h>
#include <numa.h>
#include <assert.h>
#include <execinfo.h>
#include <unistd.h>

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
		for (int i = 0; i < nptrs; i++)	{			\
			char syscom[256];						\
			printf("[bt] #%d %s\n", i, strings[i]);	\
			sprintf(syscom,"addr2line %p -e %s", buf[i], program_invocation_name);\
			int ret __attribute__((unused)) = system(syscom);	\
		}											\
		free(strings);								\
	} while (0)

#define ASSERT_TRUE(x)								\
	if (!(x)) {										\
		PRINT_BACKTRACE();							\
		assert(x);									\
	}

#ifdef __cplusplus
extern "C" {
#endif

inline static float time_diff(struct timeval time1, struct timeval time2)
{
	return time2.tv_sec - time1.tv_sec
			+ ((float)(time2.tv_usec - time1.tv_usec))/1000000;
}

inline static long time_diff_us(struct timeval time1, struct timeval time2)
{
	return (time2.tv_sec - time1.tv_sec) * 1000000
			+ (time2.tv_usec - time1.tv_usec);
}

inline static long min(long v1, long v2)
{
	return v1 > v2 ? v2 : v1;
}

inline static long max(long v1, long v2)
{
	return v1 < v2 ? v2 : v1;
}

inline static int power2(int num)
{
	if (num < 0)
		num = -num;
	if (num == 1)
		return 1;
	while (num > 1) {
		if (num % 2)
			return 0;
		num = num / 2;
	}
	return 1;
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

int isnumeric(char *str);

static const int CONST_A = 27644437;
static const long CONST_P = 68718952447L;

static inline int universal_hash(off_t v, int modulo)
{
	return (v * CONST_A) % CONST_P % modulo;
}

/*
 * These two functions should be used to allocate/free large chunks of memory.
 */
void set_use_huge_page(bool v);
void *malloc_large(size_t size);
void free_large(void *addr, size_t size);

#define likely(x)       __builtin_expect((x),1)
#define unlikely(x)     __builtin_expect((x),0)

#ifdef __cplusplus
}
#endif

#endif
