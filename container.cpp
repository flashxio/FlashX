#include <limits.h>

#include "container.h"

template<class T>
bool fifo_queue<T>::expand_queue(int new_size)
{
	assert(resizable && size < new_size);
	T *tmp = new T[new_size];
	int num = fifo_queue<T>::get_num_entries();
	for (int i = 0; i < num; i++) {
		tmp[i] = buf[(start + i) % size];
	}
	delete [] buf;
	buf = tmp;
	size = new_size;
	start = 0;
	end = num;
	return true;
}

template<class T>
int thread_safe_FIFO_queue<T>::fetch(T *entries, int num)
{
	pthread_spin_lock(&_lock);
	long curr_fetch_offset = fetch_offset;
	int n = min(num, get_actual_num_entries());
	fetch_offset += n;

	for (int i = 0; i < n; i++) {
		entries[i] = ((T*)buf)[(curr_fetch_offset + i) % this->capacity];
	}
	fetched_offset = curr_fetch_offset + n;
	pthread_spin_unlock(&_lock);
	return n;
}

template<class T>
void thread_safe_FIFO_queue<T>::addByForce(T *entries, int num)
{
	int added = add(entries, num);
	assert(added == num);
	// TODO I should make the queue extensible.
}

/**
 * this is non-blocking. 
 * It adds entries to the queue as much as possible,
 * and returns the number of entries that have been
 * added.
 */
template<class T>
int thread_safe_FIFO_queue<T>::add(T *entries, int num)
{
	pthread_spin_lock(&_lock);
	int n = min(num, get_remaining_space());
	long curr_alloc_offset = alloc_offset;
	alloc_offset += n;

	for (int i = 0; i < n; i++) {
		((T *)buf)[(curr_alloc_offset + i) % this->capacity] = entries[i];
	}
	add_offset = curr_alloc_offset + n;
	pthread_spin_unlock(&_lock);
	return n;
}

template<class T>
int blocking_FIFO_queue<T>::non_blocking_fetch(T *entries, int num)
{
	pthread_mutex_lock(&mutex);
	bool full = this->is_full();
	int ret = fifo_queue<T>::fetch(entries, num);
	pthread_mutex_unlock(&mutex);

	/* wake up all threads to send more requests */
	if (full)
		pthread_cond_broadcast(&cond);

	return ret;
}

template<class T>
int blocking_FIFO_queue<T>::fetch(T *entries, int num) {
	/* we have to wait for coming requests. */
	pthread_mutex_lock(&mutex);
	while(this->is_empty()) {
		num_empty++;
//		printf("the blocking queue %s is empty, wait...\n", name.c_str());
		pthread_cond_wait(&cond, &mutex);
	}
	bool full = this->is_full();
	int ret = fifo_queue<T>::fetch(entries, num);
	pthread_mutex_unlock(&mutex);

	/* wake up all threads to send more requests */
	if (full)
		pthread_cond_broadcast(&cond);

	return ret;
}

/**
 * This is a blocking version.
 * It adds all entries to the queue. If the queue is full,
 * wait until it can add all entries.
 */
template<class T>
int blocking_FIFO_queue<T>::add(T *entries, int num) {
	int orig_num = num;

	while (num > 0) {
		pthread_mutex_lock(&mutex);
		bool empty = this->is_empty();
		if (this->get_size() - fifo_queue<T>::get_num_entries() < num
				&& this->get_size() < max_size) {
			int new_size = this->get_size() * 2;
			new_size = max(new_size, fifo_queue<T>::get_num_entries() + num);
			new_size = min(new_size, max_size);
			printf("try to expand queue %s to %d\n", name.c_str(), new_size);
			bool ret = fifo_queue<T>::expand_queue(new_size);
			printf("expand queue to %d: %d\n", new_size, ret);
		}
		int ret = fifo_queue<T>::add(entries, num);
		entries += ret;
		num -= ret;
		/* signal the thread of reading disk to wake up. */
		if (empty)
			pthread_cond_broadcast(&cond);

		while (this->is_full() && num > 0) {
			num_full++;
//			printf("the blocking queue %s is full, wait...\n", name.c_str());
			pthread_cond_wait(&cond, &mutex);
		}
		pthread_mutex_unlock(&mutex);
	}
	return orig_num;
}

template<class T>
int blocking_FIFO_queue<T>::non_blocking_add(T *entries, int num) {
	pthread_mutex_lock(&mutex);
	bool empty = this->is_empty();
	if (this->get_size() - fifo_queue<T>::get_num_entries() < num
			&& this->get_size() < max_size) {
		int new_size = this->get_size() * 2;
		new_size = max(new_size, fifo_queue<T>::get_num_entries() + num);
		new_size = min(new_size, max_size);
		printf("try to expand queue %s to %d\n", name.c_str(), new_size);
		bool ret = fifo_queue<T>::expand_queue(new_size);
		printf("expand queue to %d: %d\n", new_size, ret);
	}
	int ret = fifo_queue<T>::add(entries, num);
	entries += ret;
	num -= ret;
	pthread_mutex_unlock(&mutex);
	/* signal the thread of reading disk to wake up. */
	if (empty)
		pthread_cond_broadcast(&cond);
	return ret;
}

template<class T>
int blocking_FIFO_queue<T>::add(fifo_queue<T> *queue)
{
	return add_partial(queue, INT_MAX);
}

template<class T>
int blocking_FIFO_queue<T>::add_partial(fifo_queue<T> *queue, int min_added)
{
	int num_added = 0;
	while (!queue->is_empty() && num_added < min_added) {
		int num = queue->get_num_entries();
		pthread_mutex_lock(&mutex);
		bool empty = this->is_empty();
		if (this->get_size() - fifo_queue<T>::get_num_entries() < num
				&& this->get_size() < max_size) {
			int new_size = this->get_size() * 2;
			new_size = max(new_size, fifo_queue<T>::get_num_entries() + num);
			new_size = min(new_size, max_size);
			printf("try to expand queue %s to %d\n", name.c_str(), new_size);
			bool ret = fifo_queue<T>::expand_queue(new_size);
			printf("expand queue to %d: %d\n", new_size, ret);
		}
		int ret = fifo_queue<T>::add(queue);
		num_added += ret;
		/* signal the thread of reading disk to wake up. */
		if (empty)
			pthread_cond_broadcast(&cond);

		while (this->is_full() && !queue->is_empty()) {
			num_full++;
			pthread_cond_wait(&cond, &mutex);
		}
		pthread_mutex_unlock(&mutex);
	}
	return num_added;
}

template<class T>
int blocking_FIFO_queue<T>::non_blocking_add(fifo_queue<T> *queue)
{
	int num_added = 0;
	if (!queue->is_empty()) {
		int num = queue->get_num_entries();
		pthread_mutex_lock(&mutex);
		bool empty = this->is_empty();
		if (this->get_size() - fifo_queue<T>::get_num_entries() < num
				&& this->get_size() < max_size) {
			int new_size = this->get_size() * 2;
			new_size = max(new_size, fifo_queue<T>::get_num_entries() + num);
			new_size = min(new_size, max_size);
			printf("try to expand queue %s to %d\n", name.c_str(), new_size);
			bool ret = fifo_queue<T>::expand_queue(new_size);
			printf("expand queue to %d: %d\n", new_size, ret);
		}
		int ret = fifo_queue<T>::add(queue);
		pthread_mutex_unlock(&mutex);

		num_added += ret;
		/* signal the thread of reading disk to wake up. */
		if (empty)
			pthread_cond_broadcast(&cond);
	}
	return num_added;
}

#ifdef USE_SHADOW_PAGE

/*
 * remove the idx'th element in the queue.
 * idx is the logical position in the queue,
 * instead of the physical index in the buffer.
 */
template<class T, int SIZE>
void embedded_queue<T, SIZE>::remove(int idx)
{
	assert(idx < num);
	/* the first element in the queue. */
	if (idx == 0) {
		pop_front();
	}
	/* the last element in the queue. */
	else if (idx == num - 1){
		num--;
	}
	/*
	 * in the middle.
	 * now we need to move data.
	 */
	else {
		T tmp[num];
		T *p = tmp;
		/* if the end of the queue is physically behind the start */
		if (start + num <= SIZE) {
			/* copy elements in front of the removed element. */
			memcpy(p, &buf[start], sizeof(T) * idx);
			p += idx;
			/* copy elements behind the removed element. */
			memcpy(p, &buf[start + idx + 1], sizeof(T) * (num - idx - 1));
		}
		/* 
		 * the removed element is between the first element
		 * and the end of the buffer.
		 */
		else if (idx + start < SIZE) {
			/* copy elements in front of the removed element. */
			memcpy(p, &buf[start], sizeof(T) * idx);
			p += idx;
			/*
			 * copy elements behind the removed element
			 * and before the end of the buffer.
			 */
			memcpy(p, &buf[start + idx + 1], sizeof(T) * (SIZE - start - idx - 1));
			p += (SIZE - start - idx - 1);
			/* copy the remaining elements in the beginning of the buffer. */
			memcpy(p, buf, sizeof(T) * (num - (SIZE - start)));
		}
		/*
		 * the removed element is between the beginning of the buffer
		 * and the last element.
		 */
		else {
			/* copy elements between the first element and the end of the buffer. */
			memcpy(p, &buf[start], sizeof(T) * (SIZE - start));
			p += (SIZE - start);
			/* copy elements between the beginning of the buffer and the removed element. */
			idx = (idx + start) % SIZE;
			memcpy(p, buf, sizeof(T) * idx);
			p += idx;
			/* copy elements after the removed element and before the last element */
			memcpy(p, &buf[idx + 1], sizeof(T) * ((start + num) % SIZE - idx - 1));
		}
		memcpy(buf, tmp, sizeof(T) * (num - 1));
		start = 0;
		num--;
	}
}    

#endif
