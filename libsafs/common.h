#ifndef __MY_COMMON_H__
#define __MY_COMMON_H__

/*
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

#include <string>
#include <vector>
#include <iostream>

#include "log.h"
#include "common_c.h"

#define ASSERT_EQ(x, y)								\
	if ((x) != (y))	{								\
		std::cerr << "x: " << x << ", y: " << y << std::endl;\
		PRINT_BACKTRACE();							\
		assert(x == y);								\
	}

#define ASSERT_LT(x, y)								\
	if ((x) <= (y))	{								\
		std::cerr << "x: " << x << ", y: " << y << std::endl;\
		PRINT_BACKTRACE();							\
		assert(x > y);								\
	}

#define ASSERT_LTEQ(x, y)							\
	if ((x) < (y))	{								\
		std::cerr << "x: " << x << ", y: " << y << std::endl;\
		PRINT_BACKTRACE();							\
		assert(x >= y);								\
	}

#define ASSERT_ST(x, y)								\
	if ((x) >= (y))	{								\
		std::cerr << "x: " << x << ", y: " << y << std::endl;\
		PRINT_BACKTRACE();							\
		assert(x < y);								\
	}

#define ASSERT_STEQ(x, y)							\
	if ((x) > (y))	{								\
		std::cerr << "x: " << x << ", y: " << y << std::endl;\
		PRINT_BACKTRACE();							\
		assert(x <= y);								\
	}

enum {
	READ,
	WRITE
};

template<class T>
inline static T min(T v1, T v2)
{
	return v1 > v2 ? v2 : v1;
}

template<class T>
inline static T max(T v1, T v2)
{
	return v1 < v2 ? v2 : v1;
}

/*
 * Check if the integer is a power of two.
 */
bool align_check(size_t alignment);

static inline std::string itoa(int n)
{
	char buf[32];
	snprintf(buf, sizeof(buf), "%d", n);
	return buf;
}

static inline std::string ltoa(long n)
{
	char buf[32];
	snprintf(buf, sizeof(buf), "%ld", n);
	return buf;
}

long str2size(std::string str);

int split_string(const std::string &str, char delim,
		std::vector<std::string> &strs);

bool check_read_content(char *buf, int size, off_t off, int file_id);
void create_write_data(char *buf, int size, off_t off, int file_id);
int get_num_omp_threads();
int get_omp_thread_num();

std::string gen_rand_name(size_t len);

template<class T>
T div_ceil(T v1, T v2)
{
	if (v1 == 0)
		return 0;
	assert(v2 > 0);
	return 1 + ((v1 - 1) / v2); // if v1 != 0
}

#endif
