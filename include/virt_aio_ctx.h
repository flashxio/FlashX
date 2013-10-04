#ifndef __VIRT_AIO_CTX_H__
#define __VIRT_AIO_CTX_H__

#include "container.h"
#include "wpaio.h"

struct req_entry {
	struct iocb *req;
	struct timeval issue_time;
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

	virt_aio_ctx(int node_id, int max_aio): aio_ctx(node_id,
			max_aio), pending_reqs(node_id, max_aio) {
		this->max_aio = max_aio;
	}
public:
	virtual void submit_io_request(struct iocb* ioq[], int num);
	virtual int io_wait(struct timespec* to, int num);

	virtual int max_io_slot();

	friend aio_ctx* aio_ctx::create_aio_ctx(int node_id, int max_aio);
};

#endif
