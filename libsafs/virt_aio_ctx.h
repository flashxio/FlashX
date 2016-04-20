#ifndef __VIRT_AIO_CTX_H__
#define __VIRT_AIO_CTX_H__

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

#include "container.h"
#include "wpaio.h"

namespace safs
{

struct req_entry {
	struct iocb *req;
	struct timeval issue_time;
};

/*
 * This class defines the data inside the virtual SSDs.
 */
class virt_data
{
public:
	virtual void create_data(int fd, void *data, int size, off_t off) = 0;
	virtual bool verify_data(int fd, void *data, int size, off_t off) = 0;
};

class ssd_perf_model
{
public:
	virtual long get_read_delay(off_t off, size_t size) = 0;
	virtual long get_write_delay(off_t off, size_t size) = 0;
};

/*
 * This emulates an SSD and provides the interface of an AIO context.
 * It is used for performance evaluation and debugging.
 * It accepts requests and returns them after a certain period of time.
 * The delay from the virtual AIO instance should be similar to the delay
 * from an SSD.
 */
class virt_aio_ctx: public aio_ctx
{
	int max_aio;
	fifo_queue<struct req_entry> pending_reqs;
	virt_data *data;
	ssd_perf_model *model;

	long read_bytes;
	long write_bytes;
	long read_bytes_ps;		// the bytes to read within a second.
	long write_bytes_ps;	// the bytes to write within a second.
	struct timeval prev_print_time;
public:
	virt_aio_ctx(virt_data *data, int node_id, int max_aio);

	virtual void submit_io_request(struct iocb* ioq[], int num);
	virtual int io_wait(struct timespec* to, int num);

	virtual int max_io_slot();

	virtual void print_stat() {
		printf("the virtual AIO context reads %ld bytes and writes %ld bytes\n",
				read_bytes, write_bytes);
	}
};

}

#endif
