#ifndef __DIRECT_PRIVATE_H__
#define __DIRECT_PRIVATE_H__

#include "read_private.h"

class direct_io: public buffered_io
{
public:
	direct_io(const char *names[], int num,
			long size): buffered_io(names, num, size, O_DIRECT | O_RDWR) {
	}

	ssize_t access(char *buf, off_t offset, ssize_t size, int access_method);
};

#endif
