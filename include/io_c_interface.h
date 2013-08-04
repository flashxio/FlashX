#ifndef __IO_C_INTERFACE_H__
#define __IO_C_INTERFACE_H__

#include <stdlib.h>

#include "common_c.h"

struct ssd_file_desc;
typedef struct ssd_file_desc * ssd_file_desc_t;

void ssd_init_io_system(const char *name, int *node_ids, int num_nodes);
void ssd_file_io_init(const char *name, int flags, int num_threads, int num_ndoes,
		int *suggested_nodes);
int ssd_create(const char *name, size_t size);
ssd_file_desc_t ssd_open(const char *name, int node_id, int flags);
ssize_t ssd_read(ssd_file_desc_t fd, void *buf, size_t count, off_t off);
ssize_t ssd_write(ssd_file_desc_t fd, void *buf, size_t count, off_t off);
int ssd_close(ssd_file_desc_t fd);
int ssd_delete(const char *name);
int ssd_fd_node_id(ssd_file_desc_t fd);

/**
 * For async IO.
 */
typedef void (*ssd_callback_func_t)(void *, int);
void ssd_set_callback(ssd_file_desc_t fd, ssd_callback_func_t);
ssize_t ssd_aread(ssd_file_desc_t fd, void *buf, size_t count, off_t off,
		void *callback_data);
ssize_t ssd_awrite(ssd_file_desc_t fd, void *buf, size_t count, off_t off,
		void *callback_data);

struct buf_pool;
struct buf_pool *create_buf_pool(int obj_size, long pool_size, int node_id);
void *alloc_buf(struct buf_pool *pool);
void free_buf(struct buf_pool *pool, void *buf);
void destroy_buf_pool(struct buf_pool *pool);

/**
 * Set the parameters of the system.
 */
size_t ssd_get_filesize(const char *name);
void set_cache_size(long size);
void set_cache_type(int type);
void set_RAID_mapping_option(int option);
void set_RAID_block_size(int num_pages);

#endif
