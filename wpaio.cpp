#include <unistd.h>
#include <stdio.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <sys/param.h>
#include <fcntl.h>
#include <errno.h>
#include <stdlib.h>
#include <assert.h>
#include <sys/select.h>

//#include <libaio.h>
#include "wpaio.h"




#define INIT_CAPACITY 8

#define AIO_BLKSIZE	(4*1024)
#define AIO_MAXIO	32
#define FREE_LIST_SIZE 128
//static int aio_blksize = 0;
//static int aio_maxio = 0;

#ifdef ENABLE_AIO
int f_fd;
long long f_offset;


static void init_free_list(struct free_list_s* free_list, int size)
{
	int i;
	free_list->array = (struct iocb**)malloc(size * sizeof(struct iocb*));
	if (NULL == free_list->array)
	{
		perror("malloc free_list");
	}
	for (i = 0; i < size; i++)
	{
		free_list->array[i] = (struct iocb*)malloc(sizeof(struct iocb));
		if (NULL == free_list->array[i])
		{
			perror("malloc free_list[i]");
		}
	}
	free_list->pos = size;
}


struct aio_ctx* create_aio_ctx(int max_aio)
{
	struct aio_ctx* a_ctx = (struct aio_ctx*)malloc(sizeof(struct aio_ctx));
	if (a_ctx == NULL)
	{
		perror("malloc aio_ctx");
		exit(1);
	}
	a_ctx->max_aio = max_aio;
	a_ctx->busy_aio = 0;
	memset(&a_ctx->ctx, 0, sizeof(a_ctx->ctx));
	printf("size of queue: %d\n", max_aio);
	int ret = io_queue_init(a_ctx->max_aio, &a_ctx->ctx);
	if (ret < 0)
	{
		perror ("io_queue_init");
		exit (1);
	}
	a_ctx->free_list = (struct free_list_s*)malloc(sizeof(struct free_list_s));
	init_free_list(a_ctx->free_list, max_aio);
	
	return a_ctx;
}



inline struct iocb* get_iocb(struct aio_ctx* a_ctx)
{
	if (!a_ctx->free_list->pos)
	{
		printf("out of free_list\n");
		exit(1);
	}
	return a_ctx->free_list->array[--a_ctx->free_list->pos];
}

inline void put_iocb(struct aio_ctx* a_ctx, struct iocb* io)
{
	a_ctx->free_list->array[a_ctx->free_list->pos++] = io;
}

struct iocb *make_io_request(struct aio_ctx *a_ctx, int fd,
		const struct iovec iov[], int count, long long offset,
		int io_type, io_callback_s *cb)
{
	struct iocb* a_req = get_iocb(a_ctx);
	if (io_type == A_READ) {
		io_prep_preadv(a_req, fd, iov, count, offset);
		io_set_callback(a_req, (io_callback_t) cb);
	}
	else if (io_type == A_WRITE) {
		io_prep_pwritev(a_req, fd, iov, count, offset);
		io_set_callback(a_req, (io_callback_t) cb);
	}
	else {
		perror("unknown operation");
		exit(1);
	}
	return a_req;
}

struct iocb* make_io_request(struct aio_ctx* a_ctx, int fd, size_t iosize, long long offset,
							 void* buffer, int io_type, io_callback_s *cb)
{
	struct iocb* a_req = get_iocb(a_ctx);
  if (io_type == A_READ)
  {
    io_prep_pread(a_req, fd, buffer, iosize, offset);
    io_set_callback(a_req, (io_callback_t) cb);
  }
  else
  {
    if (io_type == A_WRITE)
    {
      io_prep_pwrite(a_req, fd, buffer, iosize, offset);
      io_set_callback(a_req, (io_callback_t) cb);
    }
    else
    {
      perror("unknown operation");
      exit(1);
    }
  }
  return a_req;
}

int io_wait(struct aio_ctx* a_ctx, struct timespec* to, int num)
{
  struct io_event events[a_ctx->max_aio];
  struct io_event* ep = events;
  int ret, n;
  do {
	  ret = n = io_getevents(a_ctx->ctx, num, a_ctx->max_aio, events, to);
  } while (ret == -EINTR);
  if (ret < 0)
  {
    fprintf(stderr, "io_wait: %s\n", strerror(-ret));
    //exit(1);
  }

  struct iocb *iocbs[n];
  long res[n];
  long res2[n];
  io_callback_s *cbs[n];
  callback_t cb_func = NULL;
  for (int i = 0; i < n; ep++, i++)
  {
	cbs[i] = (io_callback_s *)ep->data;
	if (cb_func == NULL)
		cb_func = cbs[i]->func;
	assert(cb_func == cbs[i]->func);
	iocbs[i] = ep->obj;
	res[i] = ep->res;
	res2[i] = ep->res2;
  }

  cb_func(a_ctx->ctx, iocbs, (void **) cbs, res, res2, n);

  a_ctx->busy_aio -= n;
  for (int i = 0; i < n; i++)
	  put_iocb(a_ctx, iocbs[i]);
  return ret;
}

void submit_io_request(struct aio_ctx* a_ctx, struct iocb* ioq[], int num)
{
  int rc;
  rc = io_submit(a_ctx->ctx, num, ioq);
  if (rc < 0)
  {
    fprintf(stderr, "io_submit: %s", strerror(-rc));
    exit(1);
  }
  a_ctx->busy_aio += num;
} 

int max_io_slot(struct aio_ctx* a_ctx)
{
	return a_ctx->max_aio - a_ctx->busy_aio;
}

#endif
