#ifndef __MESSAGING_H__
#define __MESSAGING_H__

#include <pthread.h>
#include <numa.h>

inline int min(int v1, int v2)
{
	return v1 > v2 ? v2 : v1;
}

/**
 * this is a thread-safe FIFO queue.
 * It supports bulk operations.
 */
template<class T>
class bulk_queue
{
	T *buf;
	volatile int size;
	int start;
	int num_entries;
	pthread_spinlock_t _lock;
public:
	bulk_queue(int size) {
		buf = new T[size];
		this->size = size;
		start = 0;
		num_entries = 0;
		pthread_spin_init(&_lock, PTHREAD_PROCESS_PRIVATE);
	}

	~bulk_queue() {
		pthread_spin_destroy(&_lock);
		delete [] buf;
	}

	int fetch(T *entries, int num) {
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

	int add(T *entries, int num) {
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

	int get_num_entries() {
		return num_entries;
	}

	bool is_full() {
		return num_entries == size;
	}

	bool is_empty() {
		return num_entries == 0;
	}
};

template<class T>
class msg_sender
{
	T *buf;
	int buf_size;		// the max number of messages that can be buffered
	int num_current;	// the current number of messages in the buffer.
	bulk_queue<T> **dest_queues;
	int num_queues;
public:
	/**
	 * buf_size: the number of messages that can be buffered in the sender.
	 */
	msg_sender(int buf_size, bulk_queue<T> **queues, int num_queues) {
		buf = (T *) numa_alloc_local(sizeof(T) * buf_size);
		this->buf_size = buf_size;
		num_current = 0;
		dest_queues = (bulk_queue<T> **) numa_alloc_local(sizeof(bulk_queue<T> *) * num_queues);
		memcpy(dest_queues, queues, sizeof(bulk_queue<T> *) * num_queues);
		this->num_queues = num_queues;
	}

	~msg_sender() {
		numa_free(buf, sizeof(T) * buf_size);
		numa_free(dest_queues, sizeof(bulk_queue<T> *) * num_queues);
	}

	/**
	 * flush the entries in the buffer to the queues.
	 * A queue is randomly picked. If the queue is full, pick the next queue
	 * until all queues are tried or all entries in the buffer is flushed.
	 * return the number of entries that have been flushed.
	 */
	int flush() {
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
			/* 
			 * is_full is pre-check, it can't guarantee
			 * the queue isn't full.
			 */
			if (!q->is_full()) {
				int ret = q->add(tmp, num_current);
				tmp += ret;
				num_current -= ret;
				num_sent += ret;
			}
		}

		/* move the remaining entries to the beginning of the buffer. */
		if (num_current)
			memmove(buf, tmp, num_current * sizeof(T));

		return num_sent;
	}

	int send_cached(T *msg) {
		/* 
		 * if the buffer is full, and we can't flush
		 * any messages, there is nothing we can do.
		 */
		if (num_current == buf_size && flush() == 0)
			return 0;

		buf[num_current++] = *msg;
		if (num_current == buf_size)
			return flush();
		else
			/* one message has been cached. */
			return 1;
	}
};

#endif
