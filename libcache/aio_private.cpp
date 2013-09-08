#include <limits.h>

#include "aio_private.h"
#include "messaging.h"
#include "read_private.h"
#include "file_partition.h"
#include "slab_allocator.h"

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

class aio_complete_sender: public simple_sender<thread_callback_s *>
{
public:
	aio_complete_sender(int node_id,
			aio_complete_queue *queue): simple_sender<thread_callback_s *>(
			node_id, queue->get_queue(), AIO_DEPTH_PER_FILE) {
	}
};

struct thread_callback_s
{
	struct io_callback_s cb;
	async_io *aio;
	callback_allocator *cb_allocator;
	io_request req;
};

/**
 * This slab allocator makes sure all requests in the callback structure
 * are extended requests.
 */
class callback_allocator: public obj_allocator<thread_callback_s>
{
	class callback_initiator: public obj_initiator<thread_callback_s>
	{
	public:
		void init(thread_callback_s *cb) {
			cb->req.init();
		}
	} initiator;
public:
	callback_allocator(int node_id, long increase_size,
			long max_size = MAX_SIZE): obj_allocator<thread_callback_s>(node_id,
				increase_size, max_size, &initiator) {
	}
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
		int aio_depth_per_file, int node_id): io_interface(node_id), AIO_DEPTH(
			aio_depth_per_file * partition.get_num_files())
{
	cb_allocator = new callback_allocator(node_id,
			AIO_DEPTH * sizeof(thread_callback_s));;
	buf_idx = 0;
	ctx = create_aio_ctx(AIO_DEPTH);
	cb = NULL;
	num_iowait = 0;
	for (std::tr1::unordered_map<int, aio_complete_thread *>::const_iterator it
			= complete_threads.begin(); it != complete_threads.end(); it++) {
		complete_senders.insert(std::pair<int, aio_complete_sender *>(it->first,
					new aio_complete_sender(node_id, it->second->get_queue())));
		remote_tcbs.insert(std::pair<int, fifo_queue<thread_callback_s *> *>(
					it->first, fifo_queue<thread_callback_s *>::create(node_id, AIO_DEPTH)));
	}
	num_completed_reqs = 0;
	if (partition.is_active()) {
		int file_id = partition.get_file_id();
		buffered_io *io = new buffered_io(partition, node_id,
				O_DIRECT | O_RDWR);
		default_io = io;
		open_files.insert(std::pair<int, buffered_io *>(file_id, io));
	}
}

void async_io::cleanup()
{
	int slot = max_io_slot(ctx);

	while (slot < AIO_DEPTH) {
		io_wait(ctx, NULL, 1);
		slot = max_io_slot(ctx);
	}
	for (std::tr1::unordered_map<int, buffered_io *>::iterator it
			= open_files.begin(); it != open_files.end(); it++) {
		buffered_io *io = it->second;
		// Files may have been closed.
		if (io)
			io->cleanup();
	}
}

async_io::~async_io()
{
	cleanup();
	destroy_aio_ctx(ctx);
	for (std::tr1::unordered_map<int, aio_complete_sender *>::const_iterator it
			= complete_senders.begin(); it != complete_senders.end(); it++) {
		aio_complete_sender *sender = it->second;
		assert(sender->get_num_remaining() == 0);
		delete sender;
	}
	for (std::tr1::unordered_map<int, fifo_queue<thread_callback_s *> *>::const_iterator it
			= remote_tcbs.begin(); it != remote_tcbs.end(); it++) {
		assert(it->second->get_num_entries() == 0);
		delete it->second;
	}
	for (std::tr1::unordered_map<int, buffered_io *>::const_iterator it
			= open_files.begin(); it != open_files.end(); it++)
		delete it->second;
	delete cb_allocator;
}

int async_io::get_file_id() const
{
	if (default_io)
		return default_io->get_file_id();
	else
		return -1;
}

struct iocb *async_io::construct_req(io_request &io_req, callback_t cb_func)
{
	thread_callback_s *tcb = cb_allocator->alloc_obj();
	io_callback_s *cb = (io_callback_s *) tcb;

	cb->func = cb_func;
	// init doesn't pass the ownership of an extension from one request
	// to another.
	tcb->req = io_req;
	tcb->aio = this;
	tcb->cb_allocator = cb_allocator;

	assert(tcb->req.get_size() >= MIN_BLOCK_SIZE);
	assert(tcb->req.get_size() % MIN_BLOCK_SIZE == 0);
	assert(tcb->req.get_offset() % MIN_BLOCK_SIZE == 0);
	assert((long) tcb->req.get_buf() % MIN_BLOCK_SIZE == 0);
	int io_type = tcb->req.get_access_method() == READ ? A_READ : A_WRITE;
	block_identifier bid;
	buffered_io *io;
	std::tr1::unordered_map<int, buffered_io *>::iterator it
		= open_files.find(io_req.get_file_id());
	if (it != open_files.end())
		io = it->second;
	else
		io = default_io;
	assert(io);
	io->get_partition().map(tcb->req.get_offset() / PAGE_SIZE, bid);
	if (tcb->req.get_num_bufs() == 1)
		return make_io_request(ctx, io->get_fd(tcb->req.get_offset()),
				tcb->req.get_size(), bid.off * PAGE_SIZE, tcb->req.get_buf(),
				io_type, cb);
	else {
		int num_bufs = tcb->req.get_num_bufs();
		for (int i = 0; i < num_bufs; i++) {
			assert((long) tcb->req.get_buf(i) % MIN_BLOCK_SIZE == 0);
			assert(tcb->req.get_buf_size(i) % MIN_BLOCK_SIZE == 0);
		}
		struct iovec vec[num_bufs];
		int ret = tcb->req.get_vec(vec, num_bufs);
		assert(ret == num_bufs);
		struct iocb *req = make_iovec_request(ctx, io->get_fd(tcb->req.get_offset()),
				/* 
				 * iocb only contains a pointer to the io vector.
				 * the space for the IO vector is stored
				 * in the callback structure.
				 */
				vec, num_bufs, bid.off * PAGE_SIZE,
				io_type, cb);
		// I need to submit the request immediately. The iovec array is
		// allocated in the stack.
		submit_io_request(ctx, &req, 1);
		return NULL;
	}
}

void async_io::access(io_request *requests, int num, io_status *status)
{
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
		int num_iocb = 0;
		for (int i = 0; i < min; i++) {
			struct iocb *req = construct_req(*requests, aio_callback);
			requests++;
			if (req)
				reqs[num_iocb++] = req;
		}
		if (num_iocb > 0)
			submit_io_request(ctx, reqs, num_iocb);
		num -= min;
	}
	if (status)
		for (int i = 0; i < num; i++)
			status[i] = IO_PENDING;
}

void async_io::return_cb(thread_callback_s *tcbs[], int num)
{
	num_completed_reqs += num;
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
			assert(tcb->req.get_node_id() >= 0);
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
				bool sync = false;
				for (int i = 0; i < ret; i++)
					if (tcbs1[i]->req.is_sync()) {
						sync = true;
						break;
					}
				// Some requests may be synchronous, we should send them back
				// as quickly as possible.
				if (sync)
					sender->flush(false);
				else {
					int num_msg = sender->get_num_remaining();
					if (num_msg >= AIO_DEPTH_PER_FILE) {
						sender->flush(false);
						assert(sender->get_num_remaining() < num_msg);
					}
				}
			}
		}
		if (num_local > 0) {
			aio_complete_sender *sender = complete_senders[this->get_node_id()];
			sender->send_cached(local_tcbs, num_local);
			sender->flush(false);
		}
		return;
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
		cb_allocator->free(tcb);
	}
}

int async_io::open_file(const logical_file_partition &partition)
{
	int file_id = partition.get_file_id();
	if (open_files.find(file_id) == open_files.end()) {
		buffered_io *io = new buffered_io(partition, get_node_id(),
				O_DIRECT | O_RDWR);
		open_files.insert(std::pair<int, buffered_io *>(file_id, io));
	}
	else {
		fprintf(stderr, "the file id has been used\n");
		abort();
	}
	return 0;
}

int async_io::close_file(int file_id)
{
	buffered_io *io = open_files[file_id];
	// TODO I don't delete the entry.
	open_files[file_id] = NULL;
	io->cleanup();
	delete io;
	return 0;
}

void async_io::flush_requests()
{
	// There is nothing we can flush for incoming requests,
	// but we can flush completed requests.
	for (std::tr1::unordered_map<int, aio_complete_sender *>::iterator it
			= complete_senders.begin(); it != complete_senders.end(); it++) {
		aio_complete_sender *sender = it->second;
		sender->flush(true);
	}
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

	// We should try to invoke for as many requests as possible,
	// so the upper layer has the opportunity to optimize the request completion.
	std::tr1::unordered_map<io_interface *, std::vector<io_request *> > map;
	for (int i = 0; i < num; i++) {
		thread_callback_s *tcb = tcbs[i];
		std::vector<io_request *> *v;
		std::tr1::unordered_map<io_interface *, std::vector<io_request *> >::iterator it;
		io_interface *io = tcb->req.get_io();
		if ((it = map.find(io)) == map.end()) {
			map.insert(std::pair<io_interface *, std::vector<io_request *> >(
						io, std::vector<io_request *>()));
			v = &map[io];
		}
		else
			v = &it->second;

		v->push_back(&tcb->req);
	}
	for (std::tr1::unordered_map<io_interface *, std::vector<io_request *> >::iterator it
			= map.begin(); it != map.end(); it++) {
		io_interface *io = it->first;
		std::vector<io_request *> *v = &it->second;
		io->notify_completion(v->data(), v->size());
	}
	for (int i = 0; i < num; i++) {
		tcbs[i]->cb_allocator->free(tcbs[i]);
	}
	return num;
}

void aio_complete_thread::run()
{
	while(true) {
		num_completed_reqs += queue.process(AIO_NUM_PROCESS_REQS, true);
	}
}
