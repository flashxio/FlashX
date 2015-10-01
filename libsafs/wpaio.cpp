/*
 * Copyright 2014 Open Connectome Project (http://openconnecto.me)
 * Written by Da Zheng (zhengda1936@gmail.com)
 *
 * This file is part of SAFSlib.
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

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

#include "wpaio.h"
#include "virt_aio_ctx.h"
#include "parameters.h"

#define INIT_CAPACITY 8

#define AIO_BLKSIZE	(4*1024)
#define AIO_MAXIO	32
#define FREE_LIST_SIZE 128

namespace safs
{

aio_ctx::aio_ctx(int node_id, int max_aio): iocb_allocator(std::string(
			"iocb_allocator-") + itoa(node_id), node_id, true,
		sizeof(struct iocb) * max_aio, params.get_max_obj_alloc_size())
{
}

struct iocb *aio_ctx::make_iovec_request(int fd, const struct iovec iov[],
		int count, long long offset, int io_type, io_callback_s *cb)
{
	struct iocb* a_req = iocb_allocator.alloc_obj();
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

struct iocb* aio_ctx::make_io_request(int fd, size_t iosize, long long offset,
							 void* buffer, int io_type, io_callback_s *cb)
{
	struct iocb* a_req = iocb_allocator.alloc_obj();
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

int aio_ctx_impl::io_wait(struct timespec* to, int num)
{
  struct io_event events[max_aio];
  struct io_event* ep = events;
  int ret, n;
  do {
	  ret = n = io_getevents(ctx, num, max_aio, events, to);
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

  cb_func(ctx, iocbs, (void **) cbs, res, res2, n);

  busy_aio -= n;
  destroy_io_requests(iocbs, n);
  return ret;
}

void aio_ctx_impl::submit_io_request(struct iocb* ioq[], int num)
{
  int rc;
  rc = io_submit(ctx, num, ioq);
  if (rc < 0)
  {
    fprintf(stderr, "io_submit: %s", strerror(-rc));
    exit(1);
  }
  busy_aio += num;
} 

int aio_ctx_impl::max_io_slot()
{
	return max_aio - busy_aio;
}

}
