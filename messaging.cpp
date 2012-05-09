#include <stdio.h>

#include "messaging.h"

inline int min(int v1, int v2)
{
	return v1 > v2 ? v2 : v1;
}

template<class T>
int bulk_queue<T>::fetch(T *entries, int num) {
	pthread_spin_lock(&_lock);
	int n = min(num, num_entries);
	for (int i = 0; i < n; i++) {
		entries[i] = buf[(start + i) % this->size];
	}
	start = (start + n) % this->size;
	num_entries -= n;
	pthread_spin_unlock(&_lock);
	return n;
}

/**
 * this is non-blocking. 
 * It adds entries to the queue as much as possible,
 * and returns the number of entries that have been
 * added.
 */
template<class T>
int bulk_queue<T>::add(T *entries, int num) {
	pthread_spin_lock(&_lock);
	int n = min(num, this->size - num_entries);
	int end = (start + num_entries) % this->size;
	for (int i = 0; i < n; i++) {
		buf[(end + i) % this->size] = entries[i];
	}
	num_entries += n;
	pthread_spin_unlock(&_lock);
	return n;
}

template<class T>
msg_sender<T>::msg_sender(int buf_size, bulk_queue<T> **queues,
		int num_queues) {
	buf = (T *) numa_alloc_local(sizeof(T) * buf_size);
	this->buf_size = buf_size;
	num_current = 0;
	dest_queues = (bulk_queue<T> **) numa_alloc_local(
			sizeof(bulk_queue<T> *) * num_queues);
	memcpy(dest_queues, queues, sizeof(bulk_queue<T> *) * num_queues);
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
	if (num_current == 0)
		return 0;

	int base_idx;
	if (num_queues == 1)
		base_idx = 0;
	else
		base_idx = random() % num_queues;
	int num_sent = 0;
	T *tmp = buf;
	for (int i = 0; num_current > 0 && i < num_queues; i++) {
		bulk_queue<T> *q = dest_queues[(base_idx + i) % num_queues];
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
	if (num_current)
		memmove(buf, tmp, num_current * sizeof(T));

	return num_sent;
}

template<class T>
int msg_sender<T>::send_cached(T *msg) {
	/* 
	 * if the buffer is full, and we can't flush
	 * any messages, there is nothing we can do.
	 */
	if (num_current == buf_size && flush() == 0)
		return 0;

	buf[num_current++] = *msg;
	if (num_current == buf_size) {
		int ret = flush();
		return ret;
	}
	else
		/* one message has been cached. */
		return 1;
}

/**
 * these are to force to instantiate the templates
 * for io_request and io_reply.
 */
template class bulk_queue<io_request>;
template class bulk_queue<io_reply>;
template class msg_sender<io_request>;
template class msg_sender<io_reply>;
