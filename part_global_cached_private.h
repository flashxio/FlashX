#ifndef __PART_GLOBAL_CACHED_PRIVATE_H__
#define __PART_GLOBAL_CACHED_PRIVATE_H__

#define BUF_SIZE 100
#define REQ_QUEUE_SIZE 1000

#include "global_cached_private.h"

struct io_request
{
	char *buf;
	off_t offset;
	ssize_t size;
};

class part_global_cached_private;

struct thread_group
{
	int id;
	int nthreads;
	part_global_cached_private **threads;
};

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
		printf ("start: %d, num: %d\n", start, num_entries);
		pthread_spin_unlock(&_lock);
		return n;
	}

	int add(T *entries, int num) {
		pthread_spin_lock(&_lock);
		int n = min(num, this->size - num_entries);
		int end = (start + num_entries) % this->size;
		printf("add to the point %d\n", end);
		for (int i = 0; i < n; i++) {
			buf[(end + i) % this->size] = entries[i];
		}
		num_entries += n;
		pthread_spin_unlock(&_lock);
		printf ("start: %d, num: %d\n", start, num_entries);
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

class part_global_cached_private: public global_cached_private
{
	static page_cache **caches;
	static thread_group *groups;
	int num_groups;
	int group_idx;
	
	bulk_queue<io_request> queue;

	int hash_req(io_request *req)
	{
		return req->offset % num_groups;
	}

	int send(int node_id, io_request *reqs, int num) {
		int num_sent = 0;
		thread_group *group = &groups[node_id];
		int base = random() % group->nthreads;
		for (int i = 0; i < group->nthreads; i++) {
			bulk_queue<io_request> *q = &group->threads[(base
					+ i) % group->nthreads]->queue;
			/* 
			 * is_full is pre-check, it can't guarantee
			 * the queue isn't full.
			 */
			if (!q->is_full()) {
				int ret = q->add(reqs, num);
				reqs += ret;
				num -= ret;
				num_sent += ret;
			}
		}
		// TODO there should be a better way to handle this.
		assert (num == 0);
		return num_sent;
	}

	/* distribute requests to nodes. */
	void distribute_reqs(io_request *requests, int num) {
		io_request thread_reqs[num_groups][BUF_SIZE];
		int nreqs[num_groups];
		memset(nreqs, 0, sizeof(nreqs));

		for (int i = 0; i < num; i++) {
			int idx = hash_req(&requests[i]);
			assert (idx < num_groups);
			thread_reqs[idx][nreqs[idx]++] = requests[i];
			if (nreqs[idx] == BUF_SIZE) {
				send(idx, thread_reqs[idx], BUF_SIZE);
				nreqs[idx] = 0;
			}
		}
		for (int i = 0; i < num_groups; i++) {
			send(i, thread_reqs[i], nreqs[i]);
		}
	}

public:
	part_global_cached_private(int num_groups, const char *names[],
			int num, long size, int idx, long cache_size, int entry_size,
			int cache_type): global_cached_private(names, num,
				size, idx, entry_size), queue(REQ_QUEUE_SIZE) {
		this->num_groups = num_groups;
		this->group_idx = idx % num_groups;

		if (caches == NULL) {
			caches = new page_cache *[num_groups];
			for (int i = 0; i < num_groups; i++)
				caches[i] = NULL;
		}

		if (caches[group_idx] == NULL) {
			caches[group_idx] = global_cached_private::create_cache(cache_type,
					cache_size / num_groups);
		}

		if (groups == NULL) {
			groups = new thread_group[num_groups];
			int remaining = nthreads % num_groups;
			for (int i = 0; i < num_groups; i++) {
				groups[i].id = i;
				groups[i].nthreads = nthreads / num_groups + (remaining > 0);
				if (remaining > 0)
					remaining--;
				groups[i].threads = new part_global_cached_private*[groups[i].nthreads];
				for (int j = 0; j < groups[i].nthreads; j++)
					groups[i].threads[j] = NULL;
			}
		}

		thread_group *group = NULL;
		for (int i = 0; i < num_groups; i++) {
			if (groups[i].id == group_idx) {
				group = &groups[i];
				break;
			}
		}
		assert (group);
		for (int i = 0; i < group->nthreads; i++) {
			if (group->threads[i] == NULL) {
				group->threads[i] = this;
				break;
			}
		}
	}

	virtual page_cache *get_global_cache() {
		return caches[group_idx];
	}

	ssize_t access(io_request *requests, int num, int access_method) {
		distribute_reqs(requests, num);

		/* process the requests sent to this thread */
		int num_processed = 0;
		io_request local_reqs[BUF_SIZE];
		/*
		 * let's process twice as many requests as demanded by the upper layer,
		 * I hope this can help load balancing problem.
		 */
		while (!queue.is_empty() && num_processed < num * 2) {
			int num = queue.fetch(local_reqs, BUF_SIZE);
			for (int i = 0; i < num; i++) {
				io_request *req = &local_reqs[i];
				// TODO will it be better if I collect all data
				// and send them back to the initiator in blocks?
				global_cached_private::access(req->buf, req->offset, req->size, access_method);
			}
			num_processed += num;
		}

		return 0;
	}
};
thread_group *part_global_cached_private::groups;
page_cache **part_global_cached_private::caches;

#endif
