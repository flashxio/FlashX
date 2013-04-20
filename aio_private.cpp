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

class extended_io_request: public io_request
{
	char *orig_embedded_bufs[NUM_EMBEDDED_IOVECS];
	char **orig_bufs;
	slab_allocator *allocator;

	void assign_extended(extended_io_request &req) {
		this->allocator = req.allocator;
		if (req.orig_bufs == req.orig_embedded_bufs) {
			this->orig_bufs = this->orig_embedded_bufs;
			memcpy(this->orig_embedded_bufs, req.orig_embedded_bufs,
					sizeof(orig_embedded_bufs));
		}
		else {
			this->orig_bufs = req.orig_bufs;
			req.orig_bufs = NULL;
		}
	}
public:
	extended_io_request(): io_request(-1, NULL, READ, -1) {
		orig_bufs = NULL;
		allocator = NULL;
		memset(orig_embedded_bufs, 0, sizeof(orig_embedded_bufs));
	}

	/**
	 * The buffers in the IO request are allocated in a different NUMA node
	 * than the SSDs are connected to, and allocate buffers on the local node.
	 */
	extended_io_request(io_request &req,
			slab_allocator *allocator): io_request(-1, NULL, READ, -1) {
		init(req, allocator);
	}

	/**
	 * The buffers are allocated in the same NUMA node as the SSDs
	 * are connected to.
	 */
	extended_io_request(io_request &req): io_request(req) {
		allocator = NULL;
		orig_bufs = NULL;
		memset(orig_embedded_bufs, 0, sizeof(orig_embedded_bufs));
	}

	extended_io_request(extended_io_request &req): io_request(-1,
			NULL, READ, -1) {
		io_request::assign(req);
		assign_extended(req);
	}

	~extended_io_request() {
		reset();
	}

	extended_io_request &operator=(io_request &req) {
		io_request::assign(req);
		orig_bufs = NULL;
		allocator = NULL;
		memset(orig_embedded_bufs, 0, sizeof(orig_embedded_bufs));
		return *this;
	}

	extended_io_request &operator=(extended_io_request &req) {
		io_request::assign(req);
		assign_extended(req);
		return *this;
	}

	void reset();
	void init(io_request &req, slab_allocator *allocator);

	bool is_replaced() const {
		return orig_bufs != NULL;
	}

	void use_orig_bufs();
};

void extended_io_request::use_orig_bufs()
{
	for (int i = 0; i < get_num_bufs(); i++) {
		char *buf = get_buf(i);
		// This memory copy can significantly decrease the performance.
		// But it seems there isn't a better way to avoid it.
		if (this->get_access_method() == READ)
			memcpy(orig_bufs[i], buf, get_buf_size(i));
		set_buf(i, orig_bufs[i]);
		if (this->get_buf_size(i) <= PAGE_SIZE)
			allocator->free(&buf, 1);
		else
			free(buf);
	}
	// We have to reset orig_bufs because all original buffers
	// will be destroyed when the object is destructed.
	if (orig_bufs != orig_embedded_bufs)
		delete orig_bufs;
	orig_bufs = NULL;
}

void extended_io_request::init(io_request &req, slab_allocator *allocator)
{
	assert(this->get_access_method() == READ);
	io_request::assign(req);
	this->allocator = allocator;
	memset(orig_embedded_bufs, 0, sizeof(orig_embedded_bufs));
	if (this->get_num_bufs() > NUM_EMBEDDED_IOVECS)
		orig_bufs = new char *[this->get_num_bufs()];
	else
		orig_bufs = orig_embedded_bufs;
	for (int i = 0; i < this->get_num_bufs(); i++) {
		char *remote_buf = this->get_buf(i);
		char *local_buf;
		if (this->get_buf_size(i) <= PAGE_SIZE) {
			int ret = allocator->alloc(&local_buf, 1);
			assert(ret == 1);
		}
		else
			local_buf = (char *) valloc(this->get_buf_size(i));
		this->set_buf(i, local_buf);
		orig_bufs[i] = remote_buf;
	}
}

void extended_io_request::reset()
{
	if (orig_bufs) {
		char *local_pages[get_num_bufs()];
		for (int i = 0; i < get_num_bufs(); i++) {
			local_pages[i] = get_buf(i);
		}
		allocator->free(local_pages, get_num_bufs());
		if (orig_bufs != orig_embedded_bufs)
			delete orig_bufs;
		orig_bufs = NULL;
	}
}

struct thread_callback_s
{
	struct io_callback_s cb;
	async_io *aio;
	callback *aio_callback;
	obj_allocator<thread_callback_s> *cb_allocator;
	extended_io_request req;
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
		complete_queues.insert(std::pair<int, aio_complete_queue *>(it->first,
					it->second->get_queue()));
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
	if (get_node_id() == io_req.get_node_id() || get_node_id() == -1
			// It seems remote write requests can perform well.
			|| io_req.get_access_method() == WRITE)
		tcb->req = io_req;
	else {
		num_local_alloc++;
		tcb->req.init(io_req, &allocator);
	}
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
	if (complete_queues.size() > 0 && num > 5) {
		int num_remote = 0;
		int num_local = 0;
		for (int i = 0; i < num; i++) {
			thread_callback_s *tcb = tcbs[i];
			// We have allocated a local buffer, and it is a read request,
			// we need to copy data back to the issuer processor.
			// Pushing data to remote memory is more expensive than pulling
			// data from remote memory, so we let the issuer processor pull
			// data.
			if (tcb->req.is_replaced() && tcb->req.get_access_method() == READ) {
				remote_tcbs[tcb->req.get_node_id()]->push_back(tcb);
				num_remote++;
			}
			else
				local_tcbs[num_local++] = tcb;
		}
		if (num_remote > 0) {
			thread_callback_s *tcbs1[num];
			for (std::tr1::unordered_map<int, aio_complete_queue *>::iterator it
					= complete_queues.begin(); it != complete_queues.end(); it++) {
				aio_complete_queue *complete_queue = it->second;
				int ret = remote_tcbs[it->first]->fetch(tcbs1, num);
				assert(ret <= num);
				assert(remote_tcbs[it->first]->is_empty());

				int num_added = complete_queue->add(tcbs1, ret);
				// If we can't add the remote requests to the queue,
				// we have to process them locally.
				if (num_added < ret) {
					for (int i = num_added; i < ret; i++)
						local_tcbs[num_local++] = tcbs1[i];
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
		if (tcb->req.is_replaced())
			tcb->req.use_orig_bufs();
		reqs[i] = &tcb->req;
	}
	if (this->cb) {
		this->cb->invoke(reqs, num);
	}
	for (int i = 0; i < num; i++) {
		thread_callback_s *tcb = tcbs[i];
		tcb->req.reset();
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
		if (tcb->req.is_replaced())
			tcb->req.use_orig_bufs();
		io_request *reqs[1];
		reqs[0] = &tcb->req;
		tcb->aio_callback->invoke(reqs, 1);
		tcb->req.reset();
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
