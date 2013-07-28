#ifndef _GNU_SOURCE
#define _GNU_SOURCE
#endif
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>

#include "io_c_interface.h"
#include "aiori.h"

void *IOR_Create_SSDIO(char *name, IOR_param_t *test)
{
	printf("SSDIO: create and open %s\n", name);
	printf("use %d threads for IO\n", test->numThreads);
    IOR_offset_t fileSize = test->blockSize * test->segmentCount;
    if (test->filePerProc == FALSE) {
        fileSize *= test->numTasks;
    }
	ssd_create(name, fileSize);
	int *fdp = (int *) malloc(sizeof(*fdp));
	int flags = O_RDWR;
	if (test->useO_DIRECT == TRUE)
		flags |= O_DIRECT;
	ssd_io_init(name, flags, test->numThreads, test->numNodes);
	*fdp = ssd_open(name, flags);
	return (void *) fdp;
}

void *IOR_Open_SSDIO(char *name, IOR_param_t *test)
{
	printf("SSDIO: open %s\n", name);
	int *fdp = (int *) malloc(sizeof(*fdp));
	int flags = O_RDWR;
	if (test->useO_DIRECT == TRUE)
		flags |= O_DIRECT;
	ssd_io_init(name, flags, test->numThreads, test->numNodes);
	*fdp = ssd_open(name, flags);
	return (void *) fdp;
}

void IOR_SetAsyncCallback_SSDIO(void *file, AsyncCallbackFunc_t func)
{
	int fd = *(int *) file;
	ssd_set_callback(fd, func);
}

int IOR_AsyncXfer_SSDIO(int access, void *file, IOR_size_t *buffer,
		IOR_offset_t length, IOR_offset_t offset, IOR_param_t *test,
		struct AsyncData *data)
{
	int fd = *(int *) file;
	if (access == READ) {
		return ssd_aread(fd, (void *) buffer, length, offset, (void *) data);
	}
	else {
		return ssd_awrite(fd, (void *) buffer, length, offset, (void *) data);
	}
}

IOR_offset_t IOR_Xfer_SSDIO(int access, void *file, IOR_size_t *buffer,
		IOR_offset_t length, IOR_offset_t offset, IOR_param_t *test)
{
	int fd = *(int *) file;
	if (access == READ) {
		return ssd_read(fd, (void *) buffer, length, offset);
	}
	else {
		return ssd_write(fd, (void *) buffer, length, offset);
	}
}

void IOR_Close_SSDIO(void *file, IOR_param_t *test)
{
	printf("SSDIO: close file\n");
	int fd = *(int *) file;
	ssd_close(fd);
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
