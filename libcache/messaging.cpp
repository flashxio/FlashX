#include <stdio.h>

#include "messaging.h"
#include "container.cpp"
#include "io_interface.h"
#include "slab_allocator.h"
#include "cache.h"

void *io_buf::get_buf() const
{
	if (is_page)
		return u.p->get_data();
	else
		return u.buf;
}

void io_buf::init(thread_safe_page *p)
{
	assert(p->get_ref() > 0);
	u.p = p;
	size = PAGE_SIZE;
	is_page = 1;
}

void io_request::init(char *buf, off_t off, ssize_t size,
		int access_method, io_interface *io, int node_id)
{
	this->offset = off;
	if (is_extended_req()) {
		get_extension()->io = io;
		if (buf)
			add_buf(buf, size);
	}
	else {
		this->io_idx = io->get_io_idx();
		assert(io_idx <= MAX_IO_IDX);
		assert(size <= MAX_BUF_SIZE);
		this->buf_size = size;
		this->buf_addr = (long) buf;
	}
	this->access_method = access_method & 0x1;
	// by default, a request is of high priority.
	this->high_prio = 1;
	assert(node_id <= MAX_NODE_ID);
	this->node_id = node_id;
}

void io_req_extension::add_io_buf(const io_buf &buf)
{
	if (num_bufs >= vec_capacity) {
		if (vec_pointer == embedded_vecs) {
			vec_capacity = MIN_NUM_ALLOC_IOVECS;
			vec_pointer = new io_buf[vec_capacity];
			memcpy(vec_pointer, embedded_vecs,
					sizeof(embedded_vecs[0]) * NUM_EMBEDDED_IOVECS);
		}
		else {
			vec_capacity *= 2;
			io_buf *tmp = new io_buf[vec_capacity];
			memcpy(tmp, vec_pointer,
					sizeof(vec_pointer[0]) * vec_capacity / 2);
			delete [] vec_pointer;
			vec_pointer = tmp;
		}
	}
	assert(num_bufs < vec_capacity);
	vec_pointer[num_bufs] = buf;
	num_bufs++;
}

void io_req_extension::add_buf(char *buf, int size, bool is_page)
{
	io_buf tmp;
	tmp.init(buf, size, is_page);
	add_io_buf(tmp);
}

void io_req_extension::add_buf_front(char *buf, int size, bool is_page)
{
	if (num_bufs >= vec_capacity) {
		if (vec_pointer == embedded_vecs) {
			vec_capacity = MIN_NUM_ALLOC_IOVECS;
			vec_pointer = new io_buf[vec_capacity];
			memcpy(vec_pointer + 1, embedded_vecs,
					sizeof(embedded_vecs[0]) * NUM_EMBEDDED_IOVECS);
		}
		else {
			vec_capacity *= 2;
			io_buf *tmp = new io_buf[vec_capacity];
			memcpy(tmp + 1, vec_pointer,
					sizeof(vec_pointer[0]) * vec_capacity / 2);
			delete [] vec_pointer;
			vec_pointer = tmp;
		}
	}
	else {
		memmove(vec_pointer + 1, vec_pointer,
				sizeof(vec_pointer[0]) * num_bufs);
	}
	assert(num_bufs < vec_capacity);
	vec_pointer[0].init((void *) buf, size, is_page);
	num_bufs++;
}

template<class T>
msg_sender<T>::msg_sender(int buf_size, fifo_queue<T> **queues,
		int num_queues) {
	buf = (T *) numa_alloc_local(sizeof(T) * buf_size);
	this->buf_size = buf_size;
	num_current = 0;
	dest_queues = (fifo_queue<T> **) numa_alloc_local(
			sizeof(fifo_queue<T> *) * num_queues);
	memcpy(dest_queues, queues, sizeof(fifo_queue<T> *) * num_queues);
	this->num_queues = num_queues;
}

/**
 * flush the entries in the buffer to the queues.
 * A queue is randomly picked. If the queue is full, pick the next queue
 * until all queues are tried or all entries in the buffer is flushed.
 * return the number of entries that have been flushed.
 */
template<class T>
int msg_sender<T>::flush() {
	if (num_current == 0) {
		return 0;
	}

	int base_idx;
	if (num_queues == 1)
		base_idx = 0;
	else
		base_idx = random() % num_queues;
	int num_sent = 0;
	T *tmp = buf;
	for (int i = 0; num_current > 0 && i < num_queues; i++) {
		fifo_queue<T> *q = dest_queues[(base_idx + i) % num_queues];
		assert(q);

		// TODO the thread might be blocked if it's full.
		// it might hurt performance. We should try other
		// queues first before being blocked.
		int ret = q->add(tmp, num_current);
		tmp += ret;
		num_current -= ret;
		num_sent += ret;
	}

	/* move the remaining entries to the beginning of the buffer. */
	if (num_current && buf != tmp) {
		for (int i = 0; i < num_current; i++) {
			assert(tmp[i].get_offset() >= 0);
			buf[i] = tmp[i];
			assert(buf[i].get_offset() >= 0);
		}
	}

	return num_sent;
}

template<class T>
int msg_sender<T>::send_cached(T *msg) {
	/* 
	 * if the buffer is full, and we can't flush
	 * any messages, there is nothing we can do.
	 */
	if (num_current == buf_size && flush() == 0) {
		return 0;
	}

	buf[num_current++] = *msg;
	if (num_current == buf_size)
		flush();
	/* one message has been cached. */
	return 1;
}
/**
 * flush the entries in the buffer to the queues.
 * A queue is randomly picked. If the queue is full, pick the next queue
 * until all queues are tried or all entries in the buffer is flushed.
 * return the number of entries that have been flushed.
 */
template<class T>
int thread_safe_msg_sender<T>::flush() {
	int base_idx;
	if (dest_queues.size() == 1)
		base_idx = 0;
	else
		base_idx = random() % dest_queues.size();
	int num_sent = 0;
	for (size_t i = 0; !buf.is_empty() && i < dest_queues.size(); i++) {
		fifo_queue<T> *q = dest_queues[(base_idx + i) % dest_queues.size()];
		assert(q);

		// TODO the thread might be blocked if it's full.
		// it might hurt performance. We should try other
		// queues first before being blocked.
		int ret = q->add(&buf);
		num_sent += ret;
	}

	return num_sent;
}

template<class T>
int thread_safe_msg_sender<T>::send_cached(T *msg) {
	int ret = buf.add(msg, 1);
	if (ret == 1)
		return 1;
	// We expect the method is always successful.
	// so we try again and again until we succeed.
	do {
		// If the buffer is full, we should flush the buffer.
		flush();
		ret = buf.add(msg, 1);
	} while (ret == 0);
	return ret;
}

template<class T>
int thread_safe_msg_sender<T>::send_cached(T *msg, int num)
{
	int num_added = 0;
	// We expect the method is always successful.
	// so we try again and again until we succeed.
	while (true) {
		int ret = buf.add(msg, num);
		msg += ret;
		num -= ret;
		num_added += ret;
		if (num == 0)
			return num_added;
		// If the buffer is full, we should flush the buffer.
		flush();
	}
}

// Send msgs to the destinatiion queue directly without caching.
template<class T>
int thread_safe_msg_sender<T>::send(T *msg, int num)
{
	// We should flush the msgs in the cache first.
	// but it won't flush all msgs in the cache.
	flush();

	int num_sent = 0;
	while (num > 0) {
		int base_idx;
		if (dest_queues.size() == 1)
			base_idx = 0;
		else
			base_idx = random() % dest_queues.size();
		for (size_t i = 0; num > 0 && i < dest_queues.size(); i++) {
			fifo_queue<T> *q = dest_queues[(base_idx + i) % dest_queues.size()];
			assert(q);

			// TODO the thread might be blocked if it's full.
			// it might hurt performance. We should try other
			// queues first before being blocked.
			int ret = q->add(msg, num);
			msg += ret;
			num -= ret;
			num_sent += ret;
		}
	}

	return num_sent;
}

/**
 * these are to force to instantiate the templates
 * for io_request and io_reply.
 */
template class thread_safe_FIFO_queue<io_request>;
template class thread_safe_FIFO_queue<io_reply>;
template class blocking_FIFO_queue<io_request>;
template class blocking_FIFO_queue<io_reply>;
template class msg_sender<io_request>;
template class msg_sender<io_reply>;
template class thread_safe_msg_sender<io_reply>;

atomic_unsigned_integer io_req_extension::num_creates;
