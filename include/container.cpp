#include "container.h"

template<class T>
bool fifo_queue<T>::expand_queue(int new_size)
{
	int log_size = (int) ceil(log2(new_size));
	new_size = 1 << log_size;
	assert(resizable && get_size() < new_size);
	T *tmp = new T[new_size];
	int num = fifo_queue<T>::get_num_entries();
	for (int i = 0; i < num; i++) {
		tmp[i] = buf[loc_in_queue(start + i)];
	}
	delete [] buf;
	buf = tmp;
	size_mask = new_size - 1;
	start = 0;
	end = num;
	return true;
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
#ifdef DEBUG
			printf("try to expand queue %s to %d\n", name.c_str(), new_size);
#endif
			bool ret = fifo_queue<T>::expand_queue(new_size);
			assert(ret);
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
#ifdef DEBUG
		printf("try to expand queue %s to %d\n", name.c_str(), new_size);
#endif
		bool ret = fifo_queue<T>::expand_queue(new_size);
		assert(ret);
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
#ifdef DEBUG
			printf("try to expand queue %s to %d\n", name.c_str(), new_size);
#endif
			bool ret = fifo_queue<T>::expand_queue(new_size);
			assert(ret);
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
#ifdef DEBUG
			printf("try to expand queue %s to %d\n", name.c_str(), new_size);
#endif
			bool ret = fifo_queue<T>::expand_queue(new_size);
			assert(ret);
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

template<class T>
int blocking_FIFO_queue<T>::fetch(T *entries, int num, bool blocking,
		bool interruptible)
{
	/* we have to wait for coming requests. */
	pthread_mutex_lock(&mutex);
	if (blocking) {
		while(this->is_empty()) {
			num_empty++;
			if (interruptible && interrupted) {
				// We need to reset the interrupt signal.
				interrupted = false;
				break;
			}
			pthread_cond_wait(&cond, &mutex);
		}
	}
	bool full = this->is_full();
	int ret = fifo_queue<T>::fetch(entries, num);
	pthread_mutex_unlock(&mutex);

	/* wake up all threads to send more requests */
	if (full)
		pthread_cond_broadcast(&cond);

	return ret;
}

template<class T>
int blocking_FIFO_queue<T>::add(T *entries, int num, bool blocking,
		bool interruptible)
{
	// TODO
	return -1;
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
