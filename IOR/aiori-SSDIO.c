#ifndef _GNU_SOURCE
#define _GNU_SOURCE
#endif
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>

#include "io_c_interface.h"
#include "aiori.h"

IOR_offset_t GetFileSize(IOR_param_t *test);

/**
 * Create the specified number of node Ids.
 * One of them must be the suggested node Id.
 */
void get_nodes(int suggested_node, int *node_ids, int num)
{
	int i;
	int exist = 0;
	for (i = 0; i < num; i++) {
		node_ids[i] = i;
		if (i == suggested_node)
			exist = 1;
	}
	// If the suggested node doesn't exist, we make the first node id
	// the suggested node id.
	if (!exist)
		node_ids[0] = suggested_node;
}

void *IOR_Create_SSDIO(char *name, IOR_param_t *test)
{
	printf("SSDIO: create and open %s on node %d\n", name, test->nodeId);
	printf("use %d threads for IO\n", test->numThreads);
    IOR_offset_t fileSize = GetFileSize(test);
	ssd_create(name, fileSize);
	ssd_file_desc_t fdp;
	int flags = O_RDWR;
	if (test->useO_DIRECT == TRUE)
		flags |= O_DIRECT;
	int node_ids[test->numNodes];
	get_nodes(test->nodeId, node_ids, test->numNodes);
	ssd_init_io_system(name, node_ids, test->numNodes);
	if (test->filePerProc) {
		int suggested_nodes[1] = {test->nodeId};
		set_cache_size(512 * 1024 * 1024 / test->numThreads);
		ssd_file_io_init(name, flags, 1, 1, suggested_nodes);
	}
	else
		ssd_file_io_init(name, flags, test->numThreads, test->numNodes, NULL);
	fdp = ssd_open(name, test->nodeId, flags);
	return (void *) fdp;
}

void *IOR_Open_SSDIO(char *name, IOR_param_t *test)
{
	printf("SSDIO: open %s on node %d\n", name, test->nodeId);
	ssd_file_desc_t fdp;
	int flags = O_RDWR;
	if (test->useO_DIRECT == TRUE)
		flags |= O_DIRECT;
	int node_ids[test->numNodes];
	get_nodes(test->nodeId, node_ids, test->numNodes);
	ssd_init_io_system(name, node_ids, test->numNodes);
	if (test->filePerProc) {
		int suggested_nodes[1] = {test->nodeId};
		set_cache_size(512 * 1024 * 1024 / test->numThreads);
		ssd_file_io_init(name, flags, 1, 1, suggested_nodes);
	}
	else
		ssd_file_io_init(name, flags, test->numThreads, test->numNodes, NULL);
	fdp = ssd_open(name, test->nodeId, flags);
	assert(test->nodeId == ssd_fd_node_id(fdp));
	return (void *) fdp;
}

void IOR_SetAsyncCallback_SSDIO(void *file, AsyncCallbackFunc_t func)
{
	ssd_set_callback((ssd_file_desc_t) file, func);
}

int IOR_AsyncXfer_SSDIO(int access, void *file, IOR_size_t *buffer,
		IOR_offset_t length, IOR_offset_t offset, IOR_param_t *test,
		struct AsyncData *data)
{
	ssd_file_desc_t fdp = (ssd_file_desc_t) file;
	if (access == READ) {
		return ssd_aread(fdp, (void *) buffer, length, offset, (void *) data);
	}
	else {
		return ssd_awrite(fdp, (void *) buffer, length, offset, (void *) data);
	}
}

IOR_offset_t IOR_Xfer_SSDIO(int access, void *file, IOR_size_t *buffer,
		IOR_offset_t length, IOR_offset_t offset, IOR_param_t *test)
{
	ssd_file_desc_t fdp = (ssd_file_desc_t) file;
	if (access == READ) {
		return ssd_read(fdp, (void *) buffer, length, offset);
	}
	else {
		return ssd_write(fdp, (void *) buffer, length, offset);
	}
}

void IOR_Close_SSDIO(void *file, IOR_param_t *test)
{
	printf("SSDIO: close file\n");
	ssd_file_desc_t fdp = (ssd_file_desc_t) file;
	ssd_close(fdp);
	free(file);
}

void IOR_Delete_SSDIO(char *file, IOR_param_t *test)
{
	printf("SSDIO: delete %s\n", file);
	ssd_delete(file);
}

void IOR_SetVersion_SSDIO(IOR_param_t *test)
{
}

void IOR_Fsync_SSDIO(void *file, IOR_param_t *test)
{
}

IOR_offset_t IOR_GetFileSize_SSDIO(IOR_param_t *test, MPI_Comm comm, char *testFileName)
{
	IOR_offset_t size = ssd_get_filesize(testFileName);
	printf("file %s has %lld bytes\n", testFileName, size);
	return size;
}
