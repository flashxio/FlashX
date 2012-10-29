#ifndef __BALANCED_AIO_PRIVATE_H__
#define __BALANCED_AIO_PRIVATE_H__

#include "aio_private.h"

/**
 * this class extends aio_private.
 * It's only different when it accesses multiple files.
 * SSDs may have different access rates, so this class
 * tries to fully utilize all SSDs. As a result, it
 * drops requests to the slow SSDs.
 */
class balanced_aio_private: public aio_private
{
	/*
	 * This is to buffer requests, so if the requests
	 * to a file are more than other files, they will
	 * be buffered here first.
	 * This is only needed if the underlying layer reads
	 * data from multiple files.
	 */
	std::deque<io_request> *reqs_array;
	int *outstanding_nreqs;
public:
	/**
	 * @names: the names of files to be accessed
	 * @num: the number of files
	 * @size: the size of data to be accessed in all files
	 * @idx: the thread index
	 * @entry_size: the size of an entry to be accessed.
	 */
	balanced_aio_private(const char *names[], int num, long size,
			int idx, int entry_size): aio_private(names, num,
				size, idx, entry_size) {
		reqs_array = new std::deque<io_request>[this->num_open_files()];
		outstanding_nreqs = new int[this->num_open_files()];
		memset(outstanding_nreqs, 0, sizeof(int) * this->num_open_files());
	}

	~balanced_aio_private() {
		delete [] reqs_array;
		delete [] outstanding_nreqs;
	}

	void buffer_reqs(io_request *requests, int num);

	virtual void req_complete(io_request *rq) {
		aio_private::req_complete(rq);
		outstanding_nreqs[get_fd_idx(rq->offset)]--;
	}

	void drop_req(io_request *req) {
		buf->free_entry(req->get_buf());
	}

	ssize_t access(io_request *requests, int num, int access_method);
	virtual void cleanup();
};

#endif
