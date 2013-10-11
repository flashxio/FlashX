#ifndef __VIRT_AIO_CTX_H__
#define __VIRT_AIO_CTX_H__

#include "container.h"
#include "wpaio.h"

struct req_entry {
	struct iocb *req;
	struct timeval issue_time;
};

/**
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

/**
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

#endif
