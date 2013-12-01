#ifndef __CACHE_H__
#define __CACHE_H__

/**
 * Copyright 2013 Da Zheng
 *
 * This file is part of SAFSlib.
 *
 * SAFSlib is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * SAFSlib is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with SAFSlib.  If not, see <http://www.gnu.org/licenses/>.
 */

#include <unistd.h>
#include <malloc.h>
#include <string.h>
#include <stdio.h>
#include <pthread.h>
#include <assert.h>
#include <numa.h>

#include <map>

#include "common.h"
#include "concurrency.h"
#include "io_request.h"
#include "parameters.h"

const off_t PAGE_INVALID_OFFSET = ((off_t) -1) << LOG_PAGE_SIZE;

enum {
	/* All the 4 bites need to be protected by the lock bit. */
	DATA_READY_BIT = 0,
	IO_PENDING_BIT,
	/* The two bits exclude each other. */
	DIRTY_BIT,
	OLD_DIRTY_BIT,
	/*
	 * This flag indicates the IO request is in the queue for
	 * writing back.
	 */
	PREPARE_WRITEBACK,

	LOCK_BIT,

	/* 
	 * These bits don't need to be protected by the lock.
	 * They are used by LRU2Q, so it can be deleted in the future.
	 */
	ACTIVE_BIT,
	REFERENCED_BIT,
};

static inline bool page_set_flag(char &flags, int flag, bool v)
{
	char orig = flags;
	if (v)
		flags |= 0x1 << flag;
	else
		flags &= ~(0x1 << flag);
	return orig & (0x1 << flag);
}

typedef data_loc_t page_id_t;

class page
{
	/*
	 * the offset in the file in pages.
	 * it can cover a file of 8 TB.
	 */
	int offset;
	file_id_t file_id;

	/* offset in the file in bytes */
	void set_offset(off_t off) {
		offset = off >> LOG_PAGE_SIZE;
	}

	/*
	 * in pages.
	 */
protected:
	volatile void *data;
	volatile short refcnt;
	volatile char flags;
	volatile unsigned char hits;
	volatile char flush_score;
public:
	page(): data(NULL) {
		file_id = -1;
		offset = -1;
		refcnt = 0;
		flags = 0;
		hits = 0;
		flush_score = 0;
	}

	page(const page_id_t &pg_id, char *data) {
		set_offset(pg_id.get_offset());
		file_id = pg_id.get_file_id();
		this->data = data;
		refcnt = 0;
		flags = 0;
		hits = 0;
		flush_score = 0;
	}

	void set_id(const page_id_t &pg_id) {
		set_offset(pg_id.get_offset());
		file_id = pg_id.get_file_id();
	}

	bool is_valid() const {
		return offset != -1;
	}

	bool set_flag(int flag, bool v) {
		char orig = flags;
		if (v)
			flags |= 0x1 << flag;
		else
			flags &= ~(0x1 << flag);
		return orig & (0x1 << flag);
	}
	char test_flags(char flags) {
		return this->flags & flags;
	}

	void set_flush_score(int score) {
		this->flush_score = score;
	}

	int get_flush_score() const {
		return flush_score;
	}

	bool initialized() const {
		return offset != -1;
	}

	file_id_t get_file_id() const { return file_id; }
	off_t get_offset() const { return ((off_t) offset) << LOG_PAGE_SIZE; }
	void *get_data() const { return (void *) data; }

	bool data_ready() const { return flags & (0x1 << DATA_READY_BIT); }
	bool set_data_ready(bool ready) {
		return set_flag(DATA_READY_BIT, ready);
	}
	bool is_dirty() const { return flags & (0x1 << DIRTY_BIT); }
	bool set_dirty(bool dirty) {
		return set_flag(DIRTY_BIT, dirty);
	}

	bool is_old_dirty() const { return flags & (0x1 << OLD_DIRTY_BIT); }
	bool set_old_dirty(bool dirty) {
		return set_flag(OLD_DIRTY_BIT, dirty);
	}

	bool is_prepare_writeback() const {
		return flags & (0x1 << PREPARE_WRITEBACK);
	}
	bool set_prepare_writeback(bool writeback) {
		return set_flag(PREPARE_WRITEBACK, writeback);
	}

	bool set_referenced(bool referenced) {
		return set_flag(REFERENCED_BIT, referenced);
	}
	bool referenced() const {
		return flags & (0x1 << REFERENCED_BIT);
	}

	bool set_active(bool active) {
		return set_flag(ACTIVE_BIT, active);
	}
	bool active() const {
		return flags & (0x1 << ACTIVE_BIT);
	}

	void reset_hits() {
		hits = 0;
	}
	int get_hits() const {
		return hits;
	}
	void set_hits(int hits) {
		assert (hits <= 0xff);
		this->hits = hits;
	}
	/* the page is accessed */
	void hit() {
		hits++;
	}

	virtual void inc_ref() {
		refcnt++;
	}
	virtual void dec_ref() {
		refcnt--;
	}
	virtual short get_ref() const {
		return refcnt;
	}
};

class original_io_request;

class thread_safe_page: public page
{
#ifdef PTHREAD_WAIT
	pthread_cond_t ready_cond;
	pthread_cond_t dirty_cond;
	pthread_mutex_t mutex;
#endif

	bool set_flags_bit(int i, bool v) {
		char orig;
		if (v)
			orig = __sync_fetch_and_or(&flags, 0x1 << i);
		else
			orig = __sync_fetch_and_and(&flags, ~(0x1 << i));
		return orig & (0x1 << i);
	}

	bool get_flags_bit(int i) const {
		return flags & (0x1 << i);
	}

	original_io_request *reqs;
	int node_id;
	pthread_spinlock_t _lock;

public:
	thread_safe_page(): page() {
		pthread_spin_init(&_lock, PTHREAD_PROCESS_PRIVATE);
#ifdef PTHREAD_WAIT
		pthread_cond_init(&ready_cond, NULL);
		pthread_cond_init(&dirty_cond, NULL);
		pthread_mutex_init(&mutex, NULL);
#endif
		reqs = NULL;
		node_id = -1;
	}

	thread_safe_page(const page_id_t &pg_id, char *data,
			int node_id): page(pg_id, data) {
		pthread_spin_init(&_lock, PTHREAD_PROCESS_PRIVATE);
#ifdef PTHREAD_WAIT
		pthread_cond_init(&ready_cond, NULL);
		pthread_cond_init(&dirty_cond, NULL);
		pthread_mutex_init(&mutex, NULL);
#endif
		reqs = NULL;
		this->node_id = node_id;
	}

	~thread_safe_page() {
		pthread_spin_destroy(&_lock);
#ifdef PTHREAD_WAIT
		pthread_mutex_destroy(&mutex);
		pthread_cond_destroy(&ready_cond);
		pthread_cond_destroy(&dirty_cond);
#endif
	}

	int get_node_id() const {
		return node_id;
	}

	/* this is enough for x86 architecture */
	bool data_ready() const { return get_flags_bit(DATA_READY_BIT); }
	void wait_ready() {
#ifdef PTHREAD_WAIT
		printf("wait data to be ready\n");
		pthread_mutex_lock(&mutex);
#endif
		while (!data_ready()) {
#ifdef PTHREAD_WAIT
			pthread_cond_wait(&ready_cond, &mutex);
			printf("thread %ld wait for data ready\n", pthread_self());
#endif
		}
#ifdef PTHREAD_WAIT
		pthread_mutex_unlock(&mutex);
#endif
	}
	bool set_data_ready(bool ready) {
#ifdef PTHREAD_WAIT
		pthread_mutex_lock(&mutex);
#endif
		bool ret = set_flags_bit(DATA_READY_BIT, ready);
#ifdef PTHREAD_WAIT
		pthread_cond_signal(&ready_cond);
		pthread_mutex_unlock(&mutex);
#endif
		return ret;
	}

	bool set_dirty(bool dirty) {
#ifdef PTHREAD_WAIT
		pthread_mutex_lock(&mutex);
#endif
		bool ret = set_flags_bit(DIRTY_BIT, dirty);
#ifdef PTHREAD_WAIT
		pthread_cond_signal(&dirty_cond);
		pthread_mutex_unlock(&mutex);
#endif
		return ret;
	}
	void wait_cleaned() {
#ifdef PTHREAD_WAIT
		pthread_mutex_lock(&mutex);
#endif
		while (is_dirty()) {
#ifdef PTHREAD_WAIT
			pthread_cond_wait(&dirty_cond, &mutex);
			printf("thread %ld wait for the page to be written\n",
					pthread_self());
#endif
		}
#ifdef PTHREAD_WAIT
		pthread_mutex_unlock(&mutex);
#endif
	}

	bool is_io_pending() const {
		return get_flags_bit(IO_PENDING_BIT);
	}
	bool set_io_pending(bool pending) {
		return set_flags_bit(IO_PENDING_BIT, pending);
	}
	/* we set the status to io pending,
	 * and return the original status */
	bool test_and_set_io_pending() {
		int old = __sync_fetch_and_or(&flags, 0x1 << IO_PENDING_BIT);
		return old & (0x1 << IO_PENDING_BIT);
	}

	bool is_prepare_writeback() const {
		return get_flags_bit(PREPARE_WRITEBACK);
	}
	bool set_prepare_writeback(bool writeback) {
		return set_flags_bit(PREPARE_WRITEBACK, writeback);
	}

	void lock() {
		pthread_spin_lock(&_lock);
//		int old;
//		do {
//			old = __sync_fetch_and_or(&flags, 0x1 << LOCK_BIT);
//		} while (old & (0x1 << LOCK_BIT));
	}
	void unlock() {
		pthread_spin_unlock(&_lock);
//		set_flags_bit(LOCK_BIT, false);
	}

	void inc_ref() {
		__sync_fetch_and_add(&refcnt, 1);
	}
	void dec_ref() {
		__sync_fetch_and_sub(&refcnt, 1);
	}

	void wait_unused() {
		while(get_ref()) {
#ifdef DEBUG
			printf("thread %ld wait for used\n", pthread_self());
#endif
		}
	}

	void add_req(original_io_request *req);
	original_io_request *reset_reqs();
	original_io_request *get_io_req() const;
};

class page_filter
{
public:
	virtual int filter(const thread_safe_page *pages[], int num,
			const thread_safe_page *returned_pages[]) = 0;
};

class dirty_page_flusher;
class io_interface;
class page_filter;
class page_cache
{
public:
	/**
	 * This method searches for a page with the specified offset.
	 * It may evict a page if the specificed page doesn't exist.
	 * If the returned page is evicted, its original page offset is
	 * saved in `old_off'.
	 */
	virtual page *search(const page_id_t &pg_id, page_id_t &old_id) = 0;
	/**
	 * This method searches for a page with the specified offset.
	 * If the page doesn't exist, it returns NULL.
	 */
	virtual page *search(const page_id_t &pg_id) {
		return NULL;
	}
	/**
	 * The size of allocated pages in the cache in bytes.
	 */
	virtual long size() = 0;
	/* This method should be called within each thread. */
	virtual void init(io_interface *underlying) {
	}
	virtual bool shrink(int npages, char *pages[]) {
		return false;
	}
	virtual dirty_page_flusher *create_flusher(io_interface *io,
			page_cache *global_cache) {
		return NULL;
	}
	virtual void mark_dirty_pages(thread_safe_page *pages[], int num,
			io_interface *io) {
	}
	virtual int flush_dirty_pages(page_filter *filter, int max_num) {
		return 0;
	}
	virtual void flush_callback(io_request &req) {
	}
	virtual int get_node_id() const {
		return -1;
	}

	virtual void print_stat() const {
	}
};

/**
 * This is used to implement hashtable-based indexing.
 */
class frame: public thread_safe_page
{
	/* equivalent to the number of hits */
	atomic_integer wcount;
	/* equivalent to the number of references */
	atomic_integer pinning;

public:
	frame(): thread_safe_page() {
	}

	frame(const page_id_t &pg_id, char *data): thread_safe_page(pg_id, data, -1) {
	}

	/**
	 * To remove the warning of clang++, I need to define a virtual destructor.
	 * But I really don't need to clean up anything.
	 */
	virtual ~frame() {
	}

	void *volatileGetValue() {
		/*
		 * maybe a memory fence here,
		 * but it's not needed in the x86 machine.
		 */
		return get_data();
	}

	bool CASValue(void *expect,void *update) {
		return __sync_bool_compare_and_swap(&data, expect, update);
	}

	int getWC() const {
		return wcount.get();
	}

	void incrWC(int num = 1) {
		wcount.inc(num);
	}

	int decrWC(int num = 1) {
		return wcount.dec(num);
	}

	bool tryEvict() {
		return pinning.CAS(0, -1);
	}

	void evictUnshared() {
		pinning.CAS(1, -1);
	}

	int pinCount() const {
		return pinning.get();
	}

	bool pin() {
		int x;
		do {
			x = pinning.get();
			if (x <= -1)
				return false;
		} while (!pinning.CAS(x, x + 1));
		return true;
	}

	void dec_ref() {
		unpin();
	}

	void unpin() {
		pinning.dec(1);
	}
};

#endif
