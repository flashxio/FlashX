#ifndef __DIRECT_PRIVATE_H__
#define __DIRECT_PRIVATE_H__

#include "read_private.h"

class direct_private: public read_private
{
	char *pages;
	int buf_idx;
public:
	direct_private(const char *names[], int num, long size, int idx,
			int entry_size): read_private(names, num, size, idx,
			entry_size, O_DIRECT | O_RDONLY) {
		pages = (char *) valloc(PAGE_SIZE * 4096);
		buf_idx = 0;
	}

	ssize_t access(char *buf, off_t offset, ssize_t size, int access_method);
};

#endif
