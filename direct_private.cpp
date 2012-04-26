#include "direct_private.h"

ssize_t direct_private::access(char *buf, off_t offset,
		ssize_t size, int access_method) {
	assert(size >= MIN_BLOCK_SIZE);
	assert(size % MIN_BLOCK_SIZE == 0);
	assert(offset % MIN_BLOCK_SIZE == 0);
	assert((long) buf % MIN_BLOCK_SIZE == 0);
	return read_private::access(buf, offset, size, access_method);
}
