#ifndef __WORKLOAD_H__
#define __WORKLOAD_H__

#include <stdio.h>
#include <unistd.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <malloc.h>
#include <stdlib.h>

#include <string>
#include <deque>

#include "container.h"
#include "cache.h"

#define CHUNK_SLOTS 1024

const int WORKLOAD_BUF_SIZE = 20;

typedef struct workload_type
{
	off_t off;
	int size: 31;
	int read: 1;
} workload_t;

class workload_gen
{
	workload_t access;
	static int default_entry_size;
	static int default_access_method;
public:
	static void set_default_entry_size(int entry_size) {
		default_entry_size = entry_size;
	}
	static int get_default_entry_size() {
		return default_entry_size;
	}
	static void set_default_access_method(int access_method) {
		default_access_method = access_method;
	}
	static int get_default_access_method() {
		return default_access_method;
	}

	workload_gen() {
		memset(&access, 0, sizeof(access));
	}

	/**
	 * This is a wrapper to the original interface `next_offset'
	 * so more info of an access is provided.
	 */
	virtual const workload_t &next() {
		assert(default_access_method >= 0);
		access.off = next_offset();
		access.size = default_entry_size;
		access.read = default_access_method == READ;
		return access;
	}

	/**
	 * The enxt offset in bytes.
	 */
	virtual off_t next_offset() = 0;
	virtual bool has_next() = 0;
	virtual ~workload_gen() { }
	virtual void print_state() {
	}
};

class seq_workload: public workload_gen
{
	long start;
	long end;
	long cur;
	int entry_size;

public:
	seq_workload(long start, long end, int entry_size) {
		this->start = start;
		this->end = end;
		this->cur = start;
		this->entry_size = entry_size;
	}

	off_t next_offset() {
		off_t next = cur;
		cur++;
		return next * entry_size;
	}

	bool has_next() {
		return cur < end;
	}
};

class global_rand_permute_workload: public workload_gen
{
	static thread_safe_FIFO_queue<off_t> *permuted_offsets;
	int num_reads_in_100;
	int num_accesses;
	workload_t access;
	fifo_queue<off_t> local_buf;

public:
	/**
	 * @start: the index of the first entry.
	 * @end: the index of the last entry.
	 */
	global_rand_permute_workload(int stride, int length, int repeats,
			double read_ratio): local_buf(-1, WORKLOAD_BUF_SIZE) {
		if (permuted_offsets == NULL) {
			int tot_length = length * repeats;
			off_t *offsets = new off_t[tot_length];
			permute_offsets(length, repeats, stride, 0, offsets);
			permuted_offsets = thread_safe_FIFO_queue<off_t>::create(-1, tot_length);
			int ret = permuted_offsets->add(offsets, tot_length);
			assert(ret == tot_length);
			delete offsets;
		}
		this->num_reads_in_100 = (int) (read_ratio * 100);
		this->num_accesses = 0;
	}

	virtual ~global_rand_permute_workload() {
		if (permuted_offsets) {
			delete permuted_offsets;
			permuted_offsets = NULL;
		}
	}

	bool is_initialized() const {
		return permuted_offsets != NULL;
	}

	off_t next_offset() {
		assert(!local_buf.is_empty());
		return local_buf.pop_front();
	}

	bool has_next() {
		if (local_buf.is_empty()) {
			off_t offs[WORKLOAD_BUF_SIZE];
			int num = permuted_offsets->fetch(offs, WORKLOAD_BUF_SIZE);
			local_buf.add(offs, num);
		}
		return !local_buf.is_empty();
	}

	virtual const workload_t &next() {
		off_t off = next_offset();
		access.off = off;
		access.size = workload_gen::get_default_entry_size();
		if (num_accesses < num_reads_in_100)
			access.read = 1;
		else
			access.read = 0;
		num_accesses++;
		if (num_accesses >= 100)
			num_accesses = 0;
		return access;
	}
};

/**
 * In this workload generator, the expected cache hit ratio
 * can be defined by the user.
 */
class cache_hit_defined_workload: public global_rand_permute_workload
{
	long num_pages;
	double cache_hit_ratio;
	std::deque<off_t> cached_pages;
	long seq;				// the sequence number of accesses
	long cache_hit_seq;		// the sequence number of cache hits
public:
	cache_hit_defined_workload(int stride, int length, long cache_size,
			double hit_ratio, double read_ratio): global_rand_permute_workload(
				stride, length, 1, read_ratio) {
		// only to access the most recent pages.
		this->num_pages = cache_size / PAGE_SIZE / 100;
		cache_hit_ratio = hit_ratio;
		seq = 0;
		cache_hit_seq = 0;
	}

	off_t next_offset();
};

off_t *load_java_dump(const std::string &file, long &num_offsets);
workload_t *load_file_workload(const std::string &file, long &num);

/* this class reads workload from a file dumped by a Java program. */
class java_dump_workload: public workload_gen
{
	off_t *offsets;
	long curr;
	long end;

public:
	java_dump_workload(off_t offsets[], long length, long start, long end) {
		assert(length >= end);
		this->offsets = (off_t *) valloc((end - start) * sizeof(off_t));
		memcpy(this->offsets, &offsets[start], sizeof(off_t) * (end - start));
		this->end = end - start;
		this->curr = 0;
		printf("start at %ld end at %ld\n", curr, end);
	}

	virtual ~java_dump_workload() {
		free(offsets);
	}

	off_t next_offset() {
		/*
		 * the data in the file is generated by a Java program,
		 * its byte order is different from Intel architectures.
		 */
		return swap_bytesl(offsets[curr++]);
	}

	bool has_next() {
		return curr < end;
	}

	long static swap_bytesl(long num) {
		long res;

		char *src = (char *) &num;
		char *dst = (char *) &res;
		for (unsigned int i = 0; i < sizeof(long); i++) {
			dst[sizeof(long) - 1 - i] = src[i];
		}
		return res;
	}
};

/**
 * This workload generator statistically divides workloads between threads.
 */
class divided_file_workload: public workload_gen
{
	fifo_queue<workload_t> local_buf;
	workload_t curr;
public:
	divided_file_workload(workload_t workloads[], long length, int thread_id,
			int nthreads): local_buf(-1, length / nthreads) {
		long part_size = length / nthreads;
		local_buf.add(workloads + thread_id * part_size, part_size);
	}

	const workload_t &next() {
		assert(!local_buf.is_empty());
		curr = local_buf.pop_front();
		if (get_default_access_method() >= 0)
			curr.read = get_default_access_method() == READ;
		return curr;
	}

	off_t next_offset() {
		assert(!local_buf.is_empty());
		workload_t access = local_buf.pop_front();
		return access.off;
	}

	bool has_next() {
		return !local_buf.is_empty();
	}
};

/**
 * This workload generator dynamically assigns workloads to threads.
 */
class file_workload: public workload_gen
{
	static thread_safe_FIFO_queue<workload_t> *workload_queue;
	fifo_queue<workload_t> local_buf;
	workload_t curr;
public:
	file_workload(workload_t workloads[], long length, int thread_id,
			int nthreads, int read_percent = -1): local_buf(-1, WORKLOAD_BUF_SIZE) {
		if (workload_queue == NULL) {
			workload_queue = thread_safe_FIFO_queue<workload_t>::create(-1, length);
			if (read_percent >= 0) {
				for (int i = 0; i < length; i++) {
					if (random() % 100 < read_percent)
						workloads[i].read = 1;
					else
						workloads[i].read = 0;
				}
			}
			int ret = workload_queue->add(workloads, length);
			assert(ret == length);
		}
	}

	virtual ~file_workload() {
		if (workload_queue) {
			thread_safe_FIFO_queue<workload_t>::destroy(workload_queue);
			workload_queue = NULL;
		}
	}

	const workload_t &next() {
		assert(!local_buf.is_empty());
		curr = local_buf.pop_front();
		return curr;
	}

	off_t next_offset() {
		assert(!local_buf.is_empty());
		workload_t access = local_buf.pop_front();
		return access.off;
	}

	bool has_next() {
		if (local_buf.is_empty()) {
			workload_t buf[WORKLOAD_BUF_SIZE];
			int num = workload_queue->fetch(buf, WORKLOAD_BUF_SIZE);
			local_buf.add(buf, num);
		}
		return !local_buf.is_empty();
	}

	virtual void print_state() {
		printf("file workload has %d global works and %d local works\n",
				workload_queue->get_num_entries(), local_buf.get_num_entries());
	}
};

/**
 * This generates workloads mostly sequential but with random jumps.
 */
class rand_seq_workload: public workload_gen
{
	long start;
	long end;
	int entry_size;
	long seq_len;

	long num;
	long tot_num;
	long tot_num_seqs;		// The number of sequential ranges.
	long off_in_seq;		// The offset in the current sequential range.
	long seq_num;			// The sequential range ID.
public:
	rand_seq_workload(long start, long end, int entry_size, long seq_len) {
		this->start = start;
		this->end = end;
		this->entry_size = entry_size;
		this->seq_len = seq_len;

		this->num = 0;
		this->tot_num = end - start;
		this->tot_num_seqs = (end - start) * entry_size / seq_len;
		this->off_in_seq = 0;
		this->seq_num = random() % tot_num_seqs;
	}

	off_t next_offset() {
		if (off_in_seq + entry_size > seq_len) {
			seq_num = random() % tot_num_seqs;
			off_in_seq = 0;
		}
		num++;
		long tmp = off_in_seq;
		off_in_seq += entry_size;
		return seq_num * seq_len + start + tmp;
	}

	bool has_next() {
		return num < tot_num;
	}
};

class rand_workload: public workload_gen
{
	workload_t access;
	long start;
	long range;
	long num;
	long tot_accesses;
	off_t *offsets;
	bool *access_methods;
public:
	rand_workload(long start, long end, int stride, long tot_accesses,
			int read_percent) {
		this->start = start;
		this->range = end - start;
		this->tot_accesses = tot_accesses;
		num = 0;
		offsets = (off_t *) valloc(sizeof(*offsets) * tot_accesses);
		for (int i = 0; i < tot_accesses; i++) {
			offsets[i] = (start + random() % range) * stride;
		}
		access_methods = (bool *) valloc(sizeof(bool) * tot_accesses);
		for (int i = 0; i < tot_accesses; i++) {
			if (random() % 100 < read_percent)
				access_methods[i] = READ;
			else
				access_methods[i] = WRITE;
		}
	}

	virtual ~rand_workload() {
		free(offsets);
		free(access_methods);
	}

	off_t next_offset() {
		return offsets[num++];
	}

	bool has_next() {
		return num < tot_accesses;
	}

	virtual const workload_t &next() {
		access.off = offsets[num];
		access.size = get_default_entry_size();
		access.read = access_methods[num] == READ;
		num++;
		return access;
	}

	virtual void print_state() {
		printf("rand workload has %ld works left\n", tot_accesses - num);
	}
};

class workload_chunk
{
public:
	virtual bool get_workload(off_t *, int num) = 0;
	virtual ~workload_chunk() { }
};

class stride_workload_chunk: public workload_chunk
{
	long first;	// the first entry
	long last;	// the last entry but it's not included in the range
	long curr;	// the current location
	int stride;
	int entry_size;
	pthread_spinlock_t _lock;
public:
	stride_workload_chunk(long first, long last, int entry_size) {
		pthread_spin_init(&_lock, PTHREAD_PROCESS_PRIVATE);
		this->first = first;
		this->last = last;
		this->entry_size = entry_size;
		printf("first: %ld, last: %ld\n", first, last);
		curr = first;
		stride = PAGE_SIZE / entry_size;
	}

	virtual ~stride_workload_chunk() {
		pthread_spin_destroy(&_lock);
	}

	bool get_workload(off_t *offsets, int num);
};

class balanced_workload: public workload_gen
{
	off_t offsets[CHUNK_SLOTS];
	int curr;
	static workload_chunk *chunks;
public:
	balanced_workload(workload_chunk *chunks) {
		memset(offsets, 0, sizeof(offsets));
		curr = CHUNK_SLOTS;
		this->chunks = chunks;
	}

	virtual ~balanced_workload() {
		if (chunks) {
			delete chunks;
			chunks = NULL;
		}
	}

	off_t next_offset() {
		return offsets[curr++];
	}

	bool has_next() {
		if (curr < CHUNK_SLOTS)
			return true;
		else {
			bool ret = chunks->get_workload(offsets, CHUNK_SLOTS);
			curr = 0;
			return ret;
		}
	}
};

#endif
