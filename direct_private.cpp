#include "direct_private.h"

ssize_t direct_private::access(char *buf, off_t offset,
		ssize_t size, int access_method) {
	ssize_t ret;
	/* for simplicity, I assume all request sizes are smaller than a page size */
	assert(size <= PAGE_SIZE);
	if (ROUND_PAGE(offset) == offset
			&& (long) buf == ROUND_PAGE(buf)
			&& size == PAGE_SIZE) {
		ret = read_private::access(buf, offset, size, access_method);
	}
	else {
		printf("read a buffer smaller than a page size\n");
		assert(access_method == READ);
		buf_idx++;
		if (buf_idx == 4096)
			buf_idx = 0;
		char *page = pages + buf_idx * PAGE_SIZE;
		ret = read_private::access(page,
				ROUND_PAGE(offset), PAGE_SIZE, access_method);
		if (ret < 0)
			return ret;
		else
			memcpy(buf, page + (offset - ROUND_PAGE(offset)), size);
		ret = size;
	}
	return ret;
}
