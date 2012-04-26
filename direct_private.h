#ifndef __DIRECT_PRIVATE_H__
#define __DIRECT_PRIVATE_H__

#include "read_private.h"

class direct_private: public read_private
{
public:
	direct_private(const char *names[], int num, long size, int idx,
			int entry_size): read_private(names, num, size, idx,
			entry_size, O_DIRECT | O_RDWR) {
	}

	ssize_t access(char *buf, off_t offset, ssize_t size, int access_method);
};

#endif
