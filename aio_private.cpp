#include <limits.h>

#include "aio_private.h"
#include "container.cpp"

template class blocking_FIFO_queue<thread_callback_s *>;

#define EVEN_DISTRIBUTE

const int MAX_BUF_REQS = 1024 * 3;

/* 
 * each file gets the same number of outstanding requests.
 */
#ifdef EVEN_DISTRIBUTE
#define MAX_OUTSTANDING_NREQS (AIO_DEPTH / num_open_files())
#define ALLOW_DROP
#endif

struct thread_callback_s
{
	struct io_callback_s cb;
	async_io *aio;
	callback *aio_callback;
	obj_allocator<thread_callback_s> *cb_allocator;
	io_request req;
};

void aio_callback(io_context_t ctx, struct iocb* iocb[],
		void *cbs[], long res[], long res2[], int num) {
	async_io *aio = NULL;
	thread_callback_s *tcbs[num];
	for (int i = 0; i < num; i++) {
		assert(res2[i] == 0);
		tcbs[i] = (thread_callback_s *) cbs[i];
		if (aio == NULL)
			aio = tcbs[i]->aio;
		// This is true when disks are only accessed by disk access threads.
		assert(aio == tcbs[i]->aio);
	}

	aio->return_cb(tcbs, num);
}

async_io::async_io(const logical_file_partition &partition,
		const std::tr1::unordered_map<int, aio_complete_thread *> &complete_threads,
		long size, int aio_depth_per_file, int node_id): buffered_io(partition, size,
			node_id, O_DIRECT | O_RDWR), AIO_DEPTH(aio_depth_per_file *
				partition.get_num_files()), allocator(PAGE_SIZE,
				AIO_DEPTH * PAGE_SIZE, INT_MAX, node_id), cb_allocator(
					AIO_DEPTH * sizeof(thread_callback_s))
{
	printf("aio is used\n");
	buf_idx = 0;
	ctx = create_aio_ctx(AIO_DEPTH);
	cb = NULL;
	num_iowait = 0;
	for (std::tr1::unordered_map<int, aio_complete_thread *>::const_iterator it
			= complete_threads.begin(); it != complete_threads.end(); it++) {
		complete_senders.insert(std::pair<int, aio_complete_sender *>(it->first,
					new aio_complete_sender(it->second->get_queue())));
		remote_tcbs.insert(std::pair<int, fifo_queue<thread_callback_s *> *>(
					it->first, new fifo_queue<thread_callback_s *>(AIO_DEPTH)));
	}
	num_completed_reqs = 0;
	num_local_alloc = 0;
}

void async_io::cleanup()
{
	int slot = max_io_slot(ctx);

	while (slot < AIO_DEPTH) {
		io_wait(ctx, NULL);
		slot = max_io_slot(ctx);
	}
	buffered_io::cleanup();
}

async_io::~async_io()
{
	cleanup();
}

struct iocb *async_io::construct_req(io_request &io_req, callback_t cb_func)
{
	thread_callback_s *tcb = cb_allocator.alloc_obj();
	io_callback_s *cb = (io_callback_s *) tcb;

	cb->func = cb_func;
	tcb->req = io_req;
	tcb->aio = this;
	tcb->aio_callback = this->get_callback();
	tcb->cb_allocator = &cb_allocator;

	assert(tcb->req.get_size() >= MIN_BLOCK_SIZE);
	assert(tcb->req.get_size() % MIN_BLOCK_SIZE == 0);
	assert(tcb->req.get_offset() % MIN_BLOCK_SIZE == 0);
	assert((long) tcb->req.get_buf() % MIN_BLOCK_SIZE == 0);
	int io_type = tcb->req.get_access_method() == READ ? A_READ : A_WRITE;
	block_identifier bid;
	get_partition().map(tcb->req.get_offset() / PAGE_SIZE, bid);
	if (tcb->req.get_num_bufs() == 1)
		return make_io_request(ctx, get_fd(tcb->req.get_offset()),
				tcb->req.get_size(), bid.off * PAGE_SIZE, tcb->req.get_buf(),
				io_type, cb);
	else {
		int num_bufs = tcb->req.get_num_bufs();
		for (int i = 0; i < num_bufs; i++) {
			assert((long) tcb->req.get_buf(i) % MIN_BLOCK_SIZE == 0);
			assert(tcb->req.get_buf_size(i) % MIN_BLOCK_SIZE == 0);
		}
		return make_io_request(ctx, get_fd(tcb->req.get_offset()),
				/* 
				 * iocb only contains a pointer to the io vector.
				 * the space for the IO vector is stored
				 * in the callback structure.
				 */
				tcb->req.get_vec(), num_bufs, bid.off * PAGE_SIZE,
				io_type, cb);
	}
}

ssize_t async_io::access(io_request *requests, int num)
{
	ssize_t ret = 0;

	while (num > 0) {
		int slot = max_io_slot(ctx);
		if (slot == 0) {
			/*
			 * To achieve the best performance, we need to submit requests
			 * as long as there is a slot available.
			 */
			num_iowait++;
			io_wait(ctx, NULL, 1);
			slot = max_io_slot(ctx);
		}
		struct iocb *reqs[slot];
		int min = slot > num ? num : slot;
		for (int i = 0; i < min; i++) {
			ret += requests->get_size();
			reqs[i] = construct_req(*requests, aio_callback);
			requests++;
		}
		submit_io_request(ctx, reqs, min);
		num -= min;
	}
	return ret;
}

void async_io::return_cb(thread_callback_s *tcbs[], int num)
{
	thread_callback_s *local_tcbs[num];
	// If there is a dedicated thread to process the completed requests,
	// send the requests to it.
	// But we process completed requests ourselves if there aren't
	// so many. It's very often that there are only few completed requests,
	// It seems many numbers work. 5 is just randomly picked.
	if (complete_senders.size() > 0) {
		int num_remote = 0;
		int num_local = 0;
		for (int i = 0; i < num; i++) {
			thread_callback_s *tcb = tcbs[i];
			// We have allocated a local buffer, and it is a read request,
			// we need to copy data back to the issuer processor.
			// Pushing data to remote memory is more expensive than pulling
			// data from remote memory, so we let the issuer processor pull
			// data.
			if  (tcb->req.get_node_id() != this->get_node_id()) {
				remote_tcbs[tcb->req.get_node_id()]->push_back(tcb);
				num_remote++;
			}
			else
				local_tcbs[num_local++] = tcb;
		}
		if (num_remote > 0) {
			thread_callback_s *tcbs1[num];
			for (std::tr1::unordered_map<int, aio_complete_sender *>::iterator it
					= complete_senders.begin(); it != complete_senders.end(); it++) {
				aio_complete_sender *sender = it->second;
				int ret = remote_tcbs[it->first]->fetch(tcbs1, num);
				assert(ret <= num);
				assert(remote_tcbs[it->first]->is_empty());
				if (ret == 0)
					continue;

				sender->send_cached(tcbs1, ret);
				int num_msg = sender->get_num_remaining();
				if (num_msg >= AIO_DEPTH_PER_FILE) {
					sender->flush(false);
					assert(sender->get_num_remaining() < num_msg);
				}
			}
		}
		tcbs = local_tcbs;
		num = num_local;
	}
	if (num == 0)
		return;

	// Otherwise, we process them ourselves.
	io_request *reqs[num];
	for (int i = 0; i < num; i++) {
		thread_callback_s *tcb = tcbs[i];
		reqs[i] = &tcb->req;
	}
	if (this->cb) {
		this->cb->invoke(reqs, num);
	}
	for (int i = 0; i < num; i++) {
		thread_callback_s *tcb = tcbs[i];
		cb_allocator.free(tcb);
	}
	num_completed_reqs += num;
}

const int AIO_NUM_PROCESS_REQS = AIO_DEPTH_PER_FILE * 16;

int aio_complete_queue::process(int max_num, bool blocking)
{
	if (queue.is_empty() && !blocking)
		return 0;
	thread_callback_s *tcbs[max_num];
	int num;
	if (blocking)
		num = queue.fetch(tcbs, max_num);
	else {
		num = queue.non_blocking_fetch(tcbs, max_num);
	}

	for (int i = 0; i < num; i++) {
		thread_callback_s *tcb = tcbs[i];
		io_request *reqs[1];
		reqs[0] = &tcb->req;
		tcb->aio_callback->invoke(reqs, 1);
		tcb->cb_allocator->free(tcb);
	}
	return num;
}

void aio_complete_thread::run()
{
	while(true) {
		num_completed_reqs += queue.process(AIO_NUM_PROCESS_REQS, true);
	}
}
