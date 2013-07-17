#ifndef __DIRECT_PRIVATE_H__
#define __DIRECT_PRIVATE_H__

#include "read_private.h"

class direct_io: public buffered_io
{
public:
	direct_io(const logical_file_partition &partition, long size,
			int node_id): buffered_io(partition, size, node_id, O_DIRECT | O_RDWR) {
	}

	io_status access(char *buf, off_t offset, ssize_t size, int access_method);
};

#endif
