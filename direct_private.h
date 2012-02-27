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
			entry_size, O_DIRECT | O_RDWR) {
		pages = (char *) valloc(PAGE_SIZE * 4096);
		buf_idx = 0;
	}

	ssize_t access(char *buf, off_t offset, ssize_t size, int access_method) {
		ssize_t ret;
		/* for simplicity, I assume all request sizes are smaller than a page size */
		assert(size <= PAGE_SIZE);
		if (ROUND_PAGE(offset) == offset
				&& (long) buf == ROUND_PAGE(buf)
				&& size == PAGE_SIZE) {
			ret = read_private::access(buf, offset, size, access_method);
		}
		else {
			assert(access_method == READ);
			buf_idx++;
			if (buf_idx == 4096)
				buf_idx = 0;
			char *page = pages + buf_idx * PAGE_SIZE;
			ret = read_private::access(page, ROUND_PAGE(offset), PAGE_SIZE, access_method);
			if (ret < 0)
				return ret;
			else
				memcpy(buf, page + (offset - ROUND_PAGE(offset)), size);
			ret = size;
		}
		return ret;
	}
};

#endif
