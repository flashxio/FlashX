#ifndef _WPAIO_H_
#define _WPAIO_H_

#ifdef ENABLE_AIO

# ifndef _GNU_SOURCE
#define _GNU_SOURCE
#endif

#include <unistd.h>
#include <stdio.h>
#include <sys/types.h>
#include <sys/param.h>
#include <fcntl.h>
#include <stdlib.h>


#include <libaio.h>

#define A_READ 0
#define A_WRITE 1

struct free_list_s
{
	struct iocb** array;
	int pos;
};

struct aio_ctx
{
	int max_aio;
	int busy_aio;
	io_context_t ctx;
	struct free_list_s* free_list;
};

struct io_callback_s
{
	char *buf;
	off_t offset;
	ssize_t size;
	void (*func) (io_context_t, struct iocb*,
			struct io_callback_s *, long, long);
};

//extern struct aio_ctx* a_ctx;
//extern struct f_info* file_info;

extern int f_fd;
extern long long f_offset;

struct aio_ctx* create_aio_ctx(int max_aio);
//void init_free_list();
struct iocb* make_io_request(struct aio_ctx* a_ctx, int fd, size_t iosize, long long offset,
							 void* buffer, int io_type, io_callback_s *cb);
void submit_io_request(struct aio_ctx* a_ctx, struct iocb* ioq[], int num);
int io_wait(struct aio_ctx* a_ctx, struct timespec* to);
int max_io_slot(struct aio_ctx* a_ctx);

#endif

#endif
