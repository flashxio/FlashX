#ifndef __MY_COMMON_H__
#define __MY_COMMON_H__

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

bool check_read_content(char *buf, int size, off_t off);
void create_write_data(char *buf, int size, off_t off);

#endif
