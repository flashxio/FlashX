#include "io_c_interface.h"
#include "aiori.h"

void *IOR_Create_SSDIO(char *name, IOR_param_t *test)
{
	printf("SSDIO: create and open %s\n", name);
    IOR_offset_t fileSize = test->blockSize * test->segmentCount;
    if (test->filePerProc == FALSE) {
        fileSize *= test->numTasks;
    }
	ssd_create(name, fileSize);
	int *fdp = (int *) malloc(sizeof(*fdp));
	*fdp = ssd_open(name);
	return (void *) fdp;
}

void *IOR_Open_SSDIO(char *name, IOR_param_t *test)
{
	printf("SSDIO: open %s\n", name);
	int *fdp = (int *) malloc(sizeof(*fdp));
	*fdp = ssd_open(name);
	return (void *) fdp;
}

IOR_offset_t IOR_Xfer_SSDIO(int access, void *file, IOR_size_t *buffer,
		IOR_offset_t length, IOR_offset_t offset, IOR_param_t *test)
{
	int fd = *(int *) file;
	if (access == READ) {
		ssd_read(fd, (void *) buffer, length, offset);
	}
	else {
		ssd_write(fd, (void *) buffer, length, offset);
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
