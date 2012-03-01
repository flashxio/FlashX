#ifndef __PART_GLOBAL_CACHED_PRIVATE_H__
#define __PART_GLOBAL_CACHED_PRIVATE_H__

#include <errno.h>

#define BUF_SIZE 1000
#define REQ_QUEUE_SIZE 100000
#define REPLY_QUEUE_SIZE REQ_QUEUE_SIZE

#include "garbage_collection.h"
#include "global_cached_private.h"

class io_reply
{
	char *buf;
	off_t offset;
	ssize_t size: 32;
	int success: 1;
	int status: 16;
	int access_method: 1;
	void init(char *buf, off_t off, ssize_t size, int success,
			int status, int access_method) {
		this->buf = buf;
		this->offset = off;
		this->size = size;
		this->success = success;
		this->status = status;
		this->access_method = access_method;
	}
public:
	io_reply() {
		init(NULL, 0, 0, 0, 0, READ);
	}

	io_reply(io_request *req, int success, int status) {
		init(req->get_buf(), req->get_offset(), req->get_size(),
					success, status, req->get_access_method());
	}

	int get_status() {
		return status;
	}

	bool is_success() {
		return success;
	}

	char *get_buf() {
		return buf;
	}

	off_t get_offset() {
		return offset;
	}

	ssize_t get_size() {
		return size;
	}

	int get_access_method() {
		return access_method;
	}
};

class part_global_cached_private;

struct thread_group
{
	int id;
	int nthreads;
	part_global_cached_private **threads;
	page_cache *cache;
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

class part_global_cached_private: public global_cached_private
{
	static thread_group *groups;
	/* this mutex just for helping initialize cache. */
	static pthread_mutex_t init_mutex;
	/* indicates the number of threads that finish initialization. */
	static int num_finish_init;
	static pthread_mutex_t wait_mutex;
	static pthread_cond_t cond;

	int num_groups;
	int group_idx;

	/* the size of the cache associated to the thread. */
	long cache_size;
	int cache_type;
	
	bulk_queue<io_request> *request_queue;
	bulk_queue<io_reply> *reply_queue;

	volatile int finished_threads;

	/* 
	 * we have request and reply buffer to distribute requests and replies
	 * in case we can't send them all.
	 */
	/* there is a request buffer for each group it sends to */
	io_request **thread_reqs;
	int *nreqs;
	/* there is a reply buffer for each thread it sends to */
	io_reply **thread_replies;
	int *nreplies;

	int hash_req(io_request *req)
	{
		return req->get_offset() / get_entry_size() % num_groups;
	}

	long processed_requests;

public:
	static int group_id(int thread_id, int num_groups) {
		/* number of threads in a group */
		int num_threads = nthreads / num_groups;
		return thread_id / num_threads;
	}
	/* get the location of a thread in the group. */
	static int thread_idx(int thread_id, int num_groups) {
		/* number of threads in a group */
		int num_threads = nthreads / num_groups;
		return thread_id % num_threads;
	}

	part_global_cached_private *id2thread(int thread_id) {
		return groups[group_id(thread_id, num_groups)].threads[thread_idx(thread_id, num_groups)];
	}

	~part_global_cached_private() {
		for (int i = 0; i < num_groups; i++) {
			numa_free(thread_reqs[i], sizeof(io_request) * BUF_SIZE);
			numa_free(thread_replies[i], sizeof(io_reply) * BUF_SIZE);
		}
		numa_free(thread_reqs, sizeof(thread_reqs[0]) * num_groups);
		numa_free(thread_replies, sizeof(thread_replies[0]) * num_groups);
		delete request_queue;
		delete reply_queue;
	}

	int thread_init() {
		/* let's bind the thread to a specific node first. */
		struct bitmask *nodemask = numa_allocate_cpumask();
		numa_bitmask_clearall(nodemask);
		printf("thread %d is associated to node %d\n", idx, get_group_id());
		numa_bitmask_setbit(nodemask, get_group_id());
		numa_bind(nodemask);
		numa_set_strict(1);
		numa_set_bind_policy(1);

		read_private::thread_init();
		request_queue = new bulk_queue<io_request>(REQ_QUEUE_SIZE);
		reply_queue = new bulk_queue<io_reply>(REPLY_QUEUE_SIZE);

		/* initialize the request and reply buffer. */
		thread_reqs = (io_request **) numa_alloc_local(sizeof(thread_reqs[0]) * num_groups);
		thread_replies = (io_reply **) numa_alloc_local(sizeof(thread_replies[0]) * nthreads);
		for (int i = 0; i < num_groups; i++) {
			thread_reqs[i] = (io_request *) numa_alloc_local(sizeof(io_request) * BUF_SIZE);
		}
		for (int i = 0; i < nthreads; i++) {
			thread_replies[i] = (io_reply *) numa_alloc_local(sizeof(io_reply) * BUF_SIZE);
		}
		nreqs = (int *) numa_alloc_local(sizeof(*nreqs) * num_groups);
		nreplies = (int *) numa_alloc_local(sizeof(*nreplies) * num_groups);
		memset(nreqs, 0, sizeof(*nreqs) * num_groups);
		memset(nreplies, 0, sizeof(nreplies) * num_groups);

		/* 
		 * there is a global lock for all threads.
		 * so this lock makes sure cache initialization is serialized
		 */
		pthread_mutex_lock(&init_mutex);
		thread_group *group = &groups[group_idx];
		if (group->cache == NULL) {
			/* this allocates all pages for the cache. */
			page::allocate_cache(cache_size);
			group->cache = global_cached_private::create_cache(cache_type, cache_size);
		}
		num_finish_init++;
		pthread_mutex_unlock(&init_mutex);

		pthread_mutex_lock(&wait_mutex);
		while (num_finish_init < nthreads) {
			pthread_cond_wait(&cond, &wait_mutex);
		}
		pthread_mutex_unlock(&wait_mutex);
		pthread_mutex_lock(&wait_mutex);
		pthread_cond_broadcast(&cond);
		pthread_mutex_unlock(&wait_mutex);
		printf("thread %d finishes initialization\n", idx);
		return 0;
	}

	part_global_cached_private(int num_groups, const char *names[],
			int num, long size, int idx, long cache_size, int entry_size,
			int cache_type): global_cached_private(names, num,
				size, idx, entry_size) {
		assert(nthreads % num_groups == 0);
		this->num_groups = num_groups;
		this->group_idx = group_id(idx, num_groups);
		this->cache_size = cache_size / num_groups;
		this->cache_type = cache_type;
		processed_requests = 0;
		finished_threads = 0;

		printf("cache is partitioned\n");
		printf("thread id: %d, group id: %d, num groups: %d\n", idx, group_idx, num_groups);

		if (groups == NULL) {
			pthread_mutex_init(&init_mutex, NULL);
			pthread_mutex_init(&wait_mutex, NULL);
			pthread_cond_init(&cond, NULL);
			num_finish_init = 0;
			groups = new thread_group[num_groups];
			for (int i = 0; i < num_groups; i++) {
				groups[i].id = i;
				groups[i].nthreads = nthreads / num_groups;
				groups[i].threads = new part_global_cached_private*[groups[i].nthreads];
				groups[i].cache = NULL;
				for (int j = 0; j < groups[i].nthreads; j++)
					groups[i].threads[j] = NULL;
			}
		}

		/* assign a thread to a group. */
		thread_group *group = NULL;
		for (int i = 0; i < num_groups; i++) {
			if (groups[i].id == group_idx) {
				group = &groups[i];
				break;
			}
		}
		assert (group);
		int i = thread_idx(idx, num_groups);
		assert (group->threads[i] == NULL);
		group->threads[i] = this;
	}

	virtual page_cache *get_global_cache() {
		return groups[group_idx].cache;
	}

	/**
	 * Send the requests to any core in the specified node.
	 * To achieve load balancing, I pick a random core
	 * in the node, and try to copy requests to its request
	 * queue. If its queue is full, try the next core
	 * and so on. 
	 * Ideally, all requests should be sent after we try 
	 * all cores in the node. If not, we return the rest of requests.
	 */
	io_request *send(int node_id, io_request *reqs, int num) {
		int num_sent = 0;
		thread_group *group = &groups[node_id];
		int base = random() % group->nthreads;
		for (int i = 0; num > 0 && i < group->nthreads; i++) {
			part_global_cached_private *thread = group->threads[(base + i) % group->nthreads];
			bulk_queue<io_request> *q = thread->request_queue;
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
		return reqs;
	}

	/**
	 * send replies to the thread that sent the requests.
	 */
	int reply(io_request *requests, io_reply *replies, int num) {
		for (int i = 0; i < num; i++) {
			part_global_cached_private *thread = (part_global_cached_private *) requests[i].get_thread();
			int thread_id = thread->idx;
			if (nreplies[thread_id] == BUF_SIZE) {
				// TODO the buffer is already full.
				// discard the request for now.
				printf("the reply buffer for thread %d is already full\n", thread_id);
				continue;
			}
			thread_replies[thread_id][nreplies[thread_id]++] = replies[i];
			assert(thread == id2thread(thread_id));
			if (nreplies[thread_id] == BUF_SIZE) {
				int ret = thread->reply_queue->add(thread_replies[thread_id], BUF_SIZE);
				if (ret != 0 && ret != BUF_SIZE) {
					memmove(thread_replies[thread_id],
							&thread_replies[thread_id][ret], BUF_SIZE - ret);
				}
				nreplies[thread_id] = BUF_SIZE - ret;
			}
		}
		for (int i = 0; i < nthreads; i++)
			if (nreplies[i] > 0) {
				part_global_cached_private *thread = id2thread(i);
				int ret = thread->reply_queue->add(thread_replies[i], nreplies[i]);
				if (ret != 0 && ret != nreplies[i]) {
					memmove(thread_replies[i], &thread_replies[i][ret], nreplies[i] - ret);
				}
				nreplies[i] -= ret;
			}
		return 0;
	}

	/* distribute requests to nodes. */
	void distribute_reqs(io_request *requests, int num) {
		for (int i = 0; i < num; i++) {
			int idx = hash_req(&requests[i]);
			// TODO if we fail to send the requests to the specific node,
			// we should rehash it and give it to another node.
			assert (idx < num_groups);
			if (nreqs[idx] == BUF_SIZE) {
				// TODO the buffer is already full.
				// discard the request for now.
				printf("the request buffer for group %d is already full\n", idx);
				continue;
			}
			thread_reqs[idx][nreqs[idx]++] = requests[i];
			if (nreqs[idx] == BUF_SIZE) {
				io_request *remaining = send(idx, thread_reqs[idx], BUF_SIZE);
				/*
				 * if we have some requests we can't send,
				 * move them to the beginning of the buffer.
				 */
				if (remaining != thread_reqs[idx]
						&& remaining != &thread_reqs[idx][BUF_SIZE]) {
					printf("there are %ld requests left\n",
							&thread_reqs[idx][BUF_SIZE] - remaining);
					memmove(thread_reqs[idx], remaining,
							((char *) &thread_reqs[idx][BUF_SIZE] - (char *) remaining));
				}
				nreqs[idx] = &thread_reqs[idx][BUF_SIZE] - remaining;
			}
		}
		for (int i = 0; i < num_groups; i++) {
			if (nreqs[i] > 0) {
				io_request *remaining = send(i, thread_reqs[i], nreqs[i]);
				/*
				 * if we have some requests we can't send,
				 * move them to the beginning of the buffer.
				 */
				if (remaining != thread_reqs[i]
						&& remaining != &thread_reqs[i][nreqs[i]]) {
					printf("there are %ld requests left\n",
							&thread_reqs[i][nreqs[i]] - remaining);
					memmove(thread_reqs[i], remaining,
							((char *) &thread_reqs[i][nreqs[i]] - (char *) remaining));
				}
				nreqs[i] = &thread_reqs[i][nreqs[i]] - remaining;
			}
		}
	}

	/* process the requests sent to this thread */
	int process_requests(int max_nreqs) {
		int num_processed = 0;
		io_request local_reqs[BUF_SIZE];
		io_reply local_replies[BUF_SIZE];
		while (!request_queue->is_empty() && num_processed < max_nreqs) {
			int num = request_queue->fetch(local_reqs, BUF_SIZE);
			for (int i = 0; i < num; i++) {
				io_request *req = &local_reqs[i];
				// TODO will it be better if I collect all data
				// and send them back to the initiator in blocks?
				int access_method = req->get_access_method();
				assert(req->get_offset() >= 0);
				int ret = global_cached_private::access(req->get_buf(),
						req->get_offset(), req->get_size(), access_method);
				local_replies[i] = io_reply(req, ret >= 0, errno);
			}
			num_processed += num;
			reply(local_reqs, local_replies, num);
		}
		processed_requests += num_processed;
		return num_processed;
	}

	/* 
	 * process the replies and return the number
	 * of bytes that have been accessed.
	 */
	int process_replies(int max_nreplies) {
		int num_processed = 0;
		io_reply local_replies[BUF_SIZE];
		int size = 0;
		while(!reply_queue->is_empty() && num_processed < max_nreplies) {
			int num = reply_queue->fetch(local_replies, BUF_SIZE);
			for (int i = 0; i < num; i++) {
				io_reply *reply = &local_replies[i];
				int ret = process_reply(reply);
				if (ret >= 0)
					size += ret;
			}
			num_processed += num;
		}
		return num_processed;
	}

	int process_reply(io_reply *reply) {
		int access_method = reply->get_access_method();
		int ret = -1;
		if (reply->is_success()) {
			read_bytes += reply->get_size();
			if (access_method == READ) {
				assert(*(unsigned long *) reply->get_buf()
						== reply->get_offset() / sizeof(long));
			}
			ret = reply->get_size();
			buf->free_entry(reply->get_buf());
		}
		else {
			fprintf(stderr, "access error: %s\n",
					strerror(reply->get_status()));
		}
		return ret;
	}

	ssize_t access(io_request *requests, int num, int access_method) {
		distribute_reqs(requests, num);
		/*
		 * let's process up to twice as many requests as demanded by the upper layer,
		 * I hope this can help load balancing problem.
		 */
		int num_recv = 0;
		/* we need to process them at least once. */
		process_requests(num * 2);
		num_recv += process_replies(num * 4);
		while (buf->is_full()) {
			process_requests(num * 2);
			num_recv += process_replies(num * 4);
		}
		return num_recv;
	}

	void cleanup() {
		int num = 0;
		printf("thread %d: start to clean up\n", idx);
		for (int i = 0; i < num_groups; i++) {
			for (int j = 0; j < groups[i].nthreads; j++)
				__sync_fetch_and_add(&groups[i].threads[j]->finished_threads, 1);
		}
		while (!request_queue->is_empty()
				|| !reply_queue->is_empty()
				/*
				 * if finished_threads == nthreads,
				 * then all threads have reached the point.
				 */
				|| finished_threads < nthreads) {
			process_requests(200);
			process_replies(200);
			num++;
		}
		printf("thread %d processed %ld requests\n", idx, processed_requests);
	}

	int get_group_id() {
		return group_idx;
	}

	bool support_bulk() {
		return true;
	}
};
thread_group *part_global_cached_private::groups;
pthread_mutex_t part_global_cached_private::init_mutex;
int part_global_cached_private::num_finish_init;
pthread_mutex_t part_global_cached_private::wait_mutex;
pthread_cond_t part_global_cached_private::cond;

#endif
