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

	long read_bytes;
	long write_bytes;
	long read_bytes_ps;		// the bytes to read within a second.
	long write_bytes_ps;	// the bytes to write within a second.
	struct timeval prev_print_time;
public:
	virt_aio_ctx(virt_data *data, int node_id, int max_aio): aio_ctx(node_id,
			max_aio), pending_reqs(node_id, max_aio) {
		this->max_aio = max_aio;
		this->data = data;

		read_bytes = 0;
		write_bytes = 0;
		read_bytes_ps = 0;
		write_bytes_ps = 0;
		memset(&prev_print_time, 0, sizeof(prev_print_time));
	}

	virtual void submit_io_request(struct iocb* ioq[], int num);
	virtual int io_wait(struct timespec* to, int num);

	virtual int max_io_slot();

	virtual void print_stat() {
		printf("the virtual AIO context reads %ld bytes and writes %ld bytes\n",
				read_bytes, write_bytes);
	}
};

#endif
