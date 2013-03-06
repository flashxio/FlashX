#include <limits.h>

#include "aio_private.h"

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
	extended_io_request() {
		orig_bufs = NULL;
		allocator = NULL;
		memset(orig_embedded_bufs, 0, sizeof(orig_embedded_bufs));
	}

	/**
	 * The buffers in the IO request are allocated in a different NUMA node
	 * than the SSDs are connected to, and allocate buffers on the local node.
	 */
	extended_io_request(io_request &req, slab_allocator *allocator) {
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

	extended_io_request(extended_io_request &req) {
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
	char *bufs[get_num_bufs()];
	for (int i = 0; i < get_num_bufs(); i++) {
		bufs[i] = get_buf(i);
		// This memory copy can significantly decrease the performance.
		// But it seems there isn't a better way to avoid it.
		memcpy(orig_bufs[i], bufs[i], get_buf_size(i));
		set_buf(i, orig_bufs[i]);
	}
	allocator->free(bufs, get_num_bufs());
	// We have to reset orig_bufs because all original buffers
	// will be destroyed when the object is destructed.
	if (orig_bufs != orig_embedded_bufs)
		delete orig_bufs;
	orig_bufs = NULL;
}

void extended_io_request::init(io_request &req, slab_allocator *allocator)
{
	io_request::assign(req);
	this->allocator = allocator;
	memset(orig_embedded_bufs, 0, sizeof(orig_embedded_bufs));
	if (this->get_num_bufs() > NUM_EMBEDDED_IOVECS)
		orig_bufs = new char *[this->get_num_bufs()];
	else
		orig_bufs = orig_embedded_bufs;
	char *local_pages[this->get_num_bufs()];
	int npages = allocator->alloc(local_pages, this->get_num_bufs());
	assert(npages == this->get_num_bufs());
	for (int i = 0; i < this->get_num_bufs(); i++) {
		char *remote_buf = this->get_buf(i);
		char *local_buf = local_pages[i];
		assert(this->get_buf_size(i) == PAGE_SIZE);
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
	extended_io_request req;
};

void aio_callback(io_context_t ctx, struct iocb* iocb,
		void *cb, long res, long res2) {
	assert(res2 == 0);
	thread_callback_s *tcb = (thread_callback_s *) cb;

	tcb->aio->return_cb(tcb);
}

async_io::async_io(const char *names[], int num,
		long size, int aio_depth_per_file, int node_id): buffered_io(names,
			num, size, node_id, O_DIRECT | O_RDWR), AIO_DEPTH(aio_depth_per_file * num),
		allocator(PAGE_SIZE, AIO_DEPTH * PAGE_SIZE, INT_MAX, node_id)
{
	printf("aio is used\n");
	buf_idx = 0;
	ctx = create_aio_ctx(AIO_DEPTH);
	for (int i = 0; i < AIO_DEPTH * 5; i++) {
		cbs.push_back(new thread_callback_s());
	}
	cb = NULL;
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
	if (cbs.empty()) {
		fprintf(stderr, "no callback object left\n");
		return NULL;
	}

	thread_callback_s *tcb = cbs.front();
	io_callback_s *cb = (io_callback_s *) tcb;
	cbs.pop_front();

	cb->func = cb_func;
	if (get_node_id() == io_req.get_node_id() || get_node_id() == -1)
		tcb->req = io_req;
	else
		tcb->req.init(io_req, &allocator);
	tcb->aio = this;

	assert(tcb->req.get_size() >= MIN_BLOCK_SIZE);
	assert(tcb->req.get_size() % MIN_BLOCK_SIZE == 0);
	assert(tcb->req.get_offset() % MIN_BLOCK_SIZE == 0);
	assert((long) tcb->req.get_buf() % MIN_BLOCK_SIZE == 0);
	int io_type = tcb->req.get_access_method() == READ ? A_READ : A_WRITE;
	if (tcb->req.get_num_bufs() == 1)
		return make_io_request(ctx, get_fd(tcb->req.get_offset()),
				tcb->req.get_size(), tcb->req.get_offset(), tcb->req.get_buf(),
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
				tcb->req.get_vec(), num_bufs, tcb->req.get_offset(),
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

void async_io::return_cb(thread_callback_s *tcb)
{
	if (tcb->req.is_replaced())
		tcb->req.use_orig_bufs();
	if (this->cb) {
		this->cb->invoke(&tcb->req);
	}
	tcb->req.reset();
	cbs.push_back(tcb);
}
