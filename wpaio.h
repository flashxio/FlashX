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

#include "parameters.h"

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

typedef void (*callback_t) (io_context_t, struct iocb*,
		void *, long, long);

struct io_callback_s
{
	callback_t func;
	char *buf;
	off_t offset;
	ssize_t size;
};

struct iovec_callback_s
{
	callback_t func;
	struct iovec vecs[MAX_NUM_IOVECS];
	int num_vecs;
	off_t offset;
};

//extern struct aio_ctx* a_ctx;
//extern struct f_info* file_info;

extern int f_fd;
extern long long f_offset;

struct aio_ctx* create_aio_ctx(int max_aio);
//void init_free_list();
struct iocb* make_io_request(struct aio_ctx* a_ctx, int fd, size_t iosize, long long offset,
							 void* buffer, int io_type, io_callback_s *cb);
struct iocb *make_io_request(struct aio_ctx *a_ctx, int fd,
		const struct iovec *iov, int count, long long offset,
		int io_type, iovec_callback_s *cb);
void submit_io_request(struct aio_ctx* a_ctx, struct iocb* ioq[], int num);
int io_wait(struct aio_ctx* a_ctx, struct timespec* to, int num = 1);
int max_io_slot(struct aio_ctx* a_ctx);

#endif

#endif
