#ifndef __MY_COMMON_H__
#define __MY_COMMON_H__

/**
 * Copyright 2013 Da Zheng
 *
 * This file is part of SAFSlib.
 *
 * SAFSlib is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * SAFSlib is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with SAFSlib.  If not, see <http://www.gnu.org/licenses/>.
 */

#include <string>
#include <vector>
#include <iostream>

#include "common_c.h"
#include "parameters.h"

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

struct file_info
{
	std::string name;
	// The NUMA node id where the disk is connected to.
	int node_id;
};

/**
 * Check if the integer is a power of two.
 */
bool align_check(size_t alignment);

int retrieve_data_files(std::string file_file,
		std::vector<file_info> &data_files);

static inline std::string itoa(int n)
{
	char buf[32];
	snprintf(buf, sizeof(buf), "%d", n);
	return buf;
}

long str2size(std::string str);

int split_string(const std::string &str, char delim,
		std::vector<std::string> &strs);

bool check_read_content(char *buf, int size, off_t off);
void create_write_data(char *buf, int size, off_t off);

#endif
