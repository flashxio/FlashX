#ifndef __IO_C_INTERFACE_H__
#define __IO_C_INTERFACE_H__

#include <stdlib.h>

int ssd_create(const char *name, size_t size);
int ssd_open(const char *name);
ssize_t ssd_read(int fd, void *buf, size_t count, off_t off);
ssize_t ssd_write(int fd, void *buf, size_t count, off_t off);
int ssd_close(int fd);
int ssd_delete(const char *name);

/**
 * Set the parameters of the system.
 */
size_t ssd_get_filesize(const char *name);
void set_cache_size(long size);
void set_cache_type(int type);
void set_access_option(int option);
void set_num_threads(int num);
void set_num_nodes(int num);
void set_RAID_mapping_option(int option);
void set_RAID_block_size(int num_pages);

#endif
