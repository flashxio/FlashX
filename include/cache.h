#ifndef __CACHE_H__
#define __CACHE_H__

/**
 * Copyright 2014 Open Connectome Project (http://openconnecto.me)
 * Written by Da Zheng (zhengda1936@gmail.com)
 *
 * This file is part of SAFSlib.
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
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
	void hit();

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
	virtual ~page_cache() {
	}
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
	virtual void create_flusher(io_interface *io, page_cache *global_cache) {
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

	// For test
	virtual void print_stat() const {
	}

	virtual void sanity_check() const = 0;
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

/*
 * The interface of allocating a page_byte_array.
 */
class byte_array_allocator
{
public:
	virtual ~byte_array_allocator() {
	}

	virtual page_byte_array *alloc() = 0;
	virtual void free(page_byte_array *) = 0;
};

/**
 * This byte array helps users access data in non-contiguous pages
 * in the page cache. It provides a STL-style iterator and a Java-style
 * iterator to access elements in the byte array.
 */
class page_byte_array
{
	byte_array_allocator *alloc;
public:
	static void destroy(page_byte_array *arr) {
		assert(arr->alloc);
		arr->alloc->free(arr);
	}

	page_byte_array() {
		alloc = NULL;
	}

	page_byte_array(byte_array_allocator &alloc) {
		this->alloc = &alloc;
	}

	byte_array_allocator &get_allocator() {
		return *alloc;
	}

	virtual ~page_byte_array() {
	}

	/**
	 * This method locks the byte array.
	 * It's currently not implemented yet.
	 */
	virtual void lock() = 0;
	/**
	 * This method unlocks the byte array.
	 * It's currently not implemented yet.
	 */
	virtual void unlock() = 0;
	/**
	 * This method gets the size of the byte array.
	 * \return the size of the byte array.
	 */
	virtual size_t get_size() const = 0;
	/**
	 * This clones the byte array.
	 */
	virtual page_byte_array *clone() = 0;

	/**
	 * This method gets the location of the byte array in the SAFS file.
	 * \return the location of the byte array in the SAFS file.
	 */
	off_t get_offset() const {
		const thread_safe_page *p = get_page(0);
		return p->get_offset() + get_offset_in_first_page();
	}

	virtual off_t get_offset_in_first_page() const = 0;
	virtual thread_safe_page *get_page(int idx) const = 0;

	/**
	 * This is a STL-compatile iterator. Users can redefine the type of
	 * elements in the byte array and iterate the elements stored in
	 * the byte array.
	 * It overrides *, -, +, ++, ==, !=, += operators.
	 * It supports random access in the byte array.
	 */
	template<class T>
	class const_iterator: public std::iterator<std::random_access_iterator_tag, T>
	{
		const page_byte_array *arr;

		// The byte offset in the pages.
		off_t off;
		// The end byte offset in the pages.
		off_t end;
	public:
		typedef typename std::iterator<std::random_access_iterator_tag,
				T>::difference_type difference_type;

		const_iterator(const page_byte_array *arr, off_t byte_off, off_t byte_end) {
			this->arr = arr;
			off = arr->get_offset_in_first_page() + byte_off;
			end = arr->get_offset_in_first_page() + byte_end;
			assert((size_t) byte_end <= arr->get_size());
			// Either no elements cross the page boundary,
			assert((PAGE_SIZE - (off % PAGE_SIZE)) % sizeof(T) == 0
					// or the entire range is inside a page.
					|| end < PAGE_SIZE);
			assert((byte_end - byte_off) % sizeof(T) == 0);
		}

		/**
		 * This method gets the current element.
		 * \return the current element.
		 */
		T operator*() const {
			off_t pg_idx = off / PAGE_SIZE;
			off_t off_in_pg = off % PAGE_SIZE;
			char *data = (char *) arr->get_page(pg_idx)->get_data();
			return *(T *) (data + off_in_pg);
		}

		/**
		 * This method calculates the distance between this iterator
		 * and another iterator.
		 * \param it the other iterator.
		 * \return the distance between the two iterators.
		 */
		difference_type operator-(const const_iterator<T> &it) const {
			return (off - it.off) / sizeof(T);
		}

		/**
		 * This method creates a new iterator that points to a location
		 * the specified number of elements after the current iterator.
		 * The current iterator isn't changed.
		 * \param num the number of elements after the current location.
		 * \return the new iterator points to the new location.
		 */
		const_iterator<T> operator+(size_t num) const {
			const_iterator<T> ret = *this;
			ret.off += num * sizeof(T);
			assert(ret.end >= ret.off);
			return ret;
		}

		/**
		 * Move the current iterator forward by 1.
		 * This is prefix ++.
		 * \return the reference to the current iterator.
		 */
		const_iterator<T> &operator++() {
			off += sizeof(T);
			return *this;
		}

		/**
		 * Test whether the current iterator is the same as
		 * the other iterator.
		 * \param it the other iterator.
		 * \return true if they are the same.
		 */
		bool operator==(const const_iterator<T> &it) const {
			return off == it.off;
		}

		/**
		 * Test whether the current iterator isn't the same as
		 * the other iterator.
		 * \param it the other iterator.
		 * \return true if they aren't the same.
		 */
		bool operator!=(const const_iterator<T> &it) const {
			return !(*this == it);
		}

		/**
		 * This method moves the current iterator forward
		 * by the specified number.
		 * \param num the number of elements that the current iterator
		 * is moved.
		 * \return the reference to the current iterator.
		 */
		const_iterator<T> &operator+=(size_t num) {
			off += num * sizeof(T);
			assert(end >= off);
			return *this;
		}

		bool has_next() const {
			return end - off >= sizeof(T);
		}
	};

	template<class T>
	class seq_const_page_iterator
	{
		thread_safe_page *pg;
		T *data;
		T *data_end;
	public:
		seq_const_page_iterator() {
			pg = NULL;
			data = NULL;
			data_end = NULL;
		}

		seq_const_page_iterator(thread_safe_page *pg, off_t byte_off,
				off_t byte_end) {
			this->pg = pg;
			data = (T *) (((char *) pg->get_data()) + byte_off);
			data_end = (T *) (((char *) pg->get_data()) + byte_end);
		}

		int get_num_entries() const {
			return data_end - data;
		}

		bool has_next() const {
			return data < data_end;
		}

		T next() {
			T *old_data = data;
			data++;
			return *old_data;
		}

		T curr() const {
			return *data;
		}
	};

	/**
	 * This is a Java-style iterator that accesses elements in the byte
	 * array sequentially. Users can redefine the type of elements
	 * in the byte array and iterate the elements stored in the byte array.
	 */
	template<class T>
	class seq_const_iterator
	{
		const page_byte_array *arr;
		seq_const_page_iterator<T> curr_page_it;

		off_t start;
		// The byte offset in the pages.
		off_t off;
		// The end byte offset in the pages.
		off_t end;
	public:
		seq_const_iterator(const page_byte_array *arr, off_t byte_off,
				off_t byte_end) {
			this->arr = arr;
			assert((size_t) byte_end <= arr->get_size());
			assert(byte_off <= byte_end);

			start = arr->get_offset_in_first_page() + byte_off;
			off = arr->get_offset_in_first_page() + byte_off;
			end = arr->get_offset_in_first_page() + byte_end;

			if (byte_off == byte_end) {
				curr_page_it = seq_const_page_iterator<T>();
			}
			else {
				off_t pg_end;
				if (end - ROUND_PAGE(off) >= PAGE_SIZE)
					pg_end = PAGE_SIZE;
				else
					pg_end = end - ROUND_PAGE(off);
				curr_page_it = seq_const_page_iterator<T>(arr->get_page(
							off / PAGE_SIZE), off % PAGE_SIZE, pg_end);
			}
			// TODO remove the constraints later.
			assert((PAGE_SIZE - (off % PAGE_SIZE)) % sizeof(T) == 0
					// or the entire range is inside a page.
					|| end < PAGE_SIZE);
			assert((byte_end - byte_off) % sizeof(T) == 0);
		}

		/**
		 * This method tests whether the iterator can move to
		 * the next element in the byte array.
		 * \return true if there are more elements in the iterator.
		 */
		bool has_next() {
			if (curr_page_it.has_next())
				return true;
			else {
				off = ROUNDUP_PAGE(off + 1);
				if (off < end) {
					curr_page_it = seq_const_page_iterator<T>(arr->get_page(
								off / PAGE_SIZE), 0, min(PAGE_SIZE, end - off));
					return true;
				}
				else {
					return false;
				}
			}
		}

		/**
		 * This method gets the total number of remaining elements
		 * in the iterator.
		 * \return the total number of remaining elements.
		 */
		int get_num_tot_entries() const {
			return (end - off) / sizeof(T);
		}

		/**
		 * This method gets the number of remaining elements
		 * in the page where the current iterator is on.
		 * \return the number of remaining elements in the current page.
		 */
		int get_num_entries_in_page() const {
			return curr_page_it.get_num_entries();
		}

		/**
		 * This method moves to the next element.
		 * \return the current element.
		 */
		T next() {
			return curr_page_it.next();
		}

		bool move_to(size_t idx) {
			off = start + idx * sizeof(T);
			if (off >= end)
				return false;

			off_t pg_end;
			if (end - ROUND_PAGE(off) >= PAGE_SIZE)
				pg_end = PAGE_SIZE;
			else
				pg_end = end - ROUND_PAGE(off);
			curr_page_it = seq_const_page_iterator<T>(arr->get_page(
						off / PAGE_SIZE), off % PAGE_SIZE, pg_end);
			return true;
		}

		/**
		 * This method gets the current element.
		 * \return the current element.
		 */
		T curr() const {
			seq_const_iterator<T> *it = (seq_const_iterator<T> *) this;
			// TODO I need to fix this.
			assert(it->has_next());
			return curr_page_it.curr();
		}
	};

	/**
	 * This macro iterates over all elements in the iterator.
	 * This is the most light-weight method of iterating over all elements
	 * in an iterator.
	 * \param type is the type of element.
	 * \param v is the variable that stores the current element.
	 * \param it is the Java-style iterator.
	 */
#define PAGE_FOREACH(type, v, it)					\
	for (int page_foreach_idx = 0; (it).has_next();) {			\
		for (int page_foreach_i = 0; page_foreach_i < (it).get_num_entries_in_page();\
				page_foreach_i++, page_foreach_idx++) {	\
			type v = it.next();


#define PAGE_FOREACH_END }}

	template<class T>
	class iterator: public std::iterator<std::random_access_iterator_tag, T>
	{
		const page_byte_array *arr;

		// The byte offset in the page array.
		off_t off;
		// The end byte offset in the page array.
		off_t end;
	public:
		typedef typename std::iterator<std::random_access_iterator_tag,
				T>::difference_type difference_type;

		iterator(page_byte_array *arr, off_t byte_off, off_t byte_end) {
			this->arr = arr;
			off = arr->get_offset_in_first_page() + byte_off;
			end = arr->get_offset_in_first_page() + byte_end;
			assert((size_t) byte_end <= arr->get_size());
			// Either no elements cross the page boundary,
			assert((PAGE_SIZE - (off % PAGE_SIZE)) % sizeof(T) == 0
					// or the entire range is inside a page.
					|| end < PAGE_SIZE);
			assert((byte_end - byte_off) % sizeof(T) == 0);

			// This iterator can change data, so I set all pages that can be
			// accessed by the iterator dirty.
			off_t num_pages = ROUNDUP_PAGE(end) / PAGE_SIZE;
			for (off_t i = off / PAGE_SIZE; i < num_pages; i++)
				arr->get_page(i)->set_dirty(true);
		}

		T &operator*() {
			off_t pg_idx = off / PAGE_SIZE;
			off_t off_in_pg = off % PAGE_SIZE;
			char *data = (char *) arr->get_page(pg_idx)->get_data();
			return *(T *) (data + off_in_pg);
		}

		difference_type operator-(const iterator<T> &it) const {
			return (off - it.off) / sizeof(T);
		}

		// Prefix ++
		iterator<T> &operator++() {
			off += sizeof(T);
			return *this;
		}

		bool operator==(const iterator<T> &it) const {
			return off == it.off;
		}

		bool operator!=(const iterator<T> &it) const {
			return !(*this == it);
		}

		iterator<T> &operator+=(int num) {
			assert(num >= 0);
			off += num * sizeof(T);
			assert(end >= off);
			return *this;
		}

		bool has_next() const {
			return end - off >= sizeof(T);
		}
	};

	/**
	 * This method gets the STL-compatible iterator that points to
	 * the user-specified location in the byte array.
	 * \param byte_off the offset (in bytes) in the byte array.
	 * \return the iterator to the specified location in the byte array.
	 */
	template<class T>
	const_iterator<T> begin(off_t byte_off = 0) const {
		return const_iterator<T>(this, byte_off, this->get_size());
	}

	/**
	 * This method gets the STL-compatible iterator that points to the end
	 * of the byte array.
	 * \return the iterator to the end of the byte array.
	 */
	template<class T>
	const_iterator<T> end() const {
		return const_iterator<T>(this, this->get_size(), this->get_size());
	}

	/**
	 * This method returns a pair of iterators that iterate over the range
	 * on the byte array, specified by the user.
	 * \param byte_off the beginning offset (in bytes) in the byte array.
	 * \param byte_end the end offset (in bytes) in the byte array.
	 * \return a pair of iterators. The first iterator points to
	 * the beginning offset and the second one points to the end offset.
	 */
	template<class T>
	std::pair<const_iterator<T>, const_iterator<T> > get_iterator(
			off_t byte_off, off_t byte_end) const {
		return std::pair<const_iterator<T>, const_iterator<T> >(
				const_iterator<T>(this, byte_off, byte_end),
				const_iterator<T>(this, byte_end, byte_end));
	}

	/**
	 * This method gets the Java-style iterator that iterates elements
	 * in the specified range in the byte array.
	 * \param byte_off the beginning offset (in bytes) in the byte array.
	 * \param byte_end the end offset (in bytes) in the byte array.
	 * \return the Java-style iterator.
	 */
	template<class T>
	seq_const_iterator<T> get_seq_iterator(off_t byte_off, off_t byte_end) const {
		return seq_const_iterator<T>(this, byte_off, byte_end);
	}

	template<class T>
	iterator<T> begin(off_t byte_off = 0) {
		return iterator<T>(this, byte_off, this->get_size());
	}

	template<class T>
	iterator<T> end() {
		return iterator<T>(this, this->get_size(), this->get_size());
	}

	template<class T>
	std::pair<iterator<T>, iterator<T> > get_iterator(
			off_t byte_off, off_t byte_end) const {
		return std::pair<iterator<T>, iterator<T> >(
				iterator<T>(this, byte_off, byte_end),
				iterator<T>(this, byte_end, byte_end));
	}

	template<class T>
	T get(off_t idx) const {
		T ele;
		memcpy(sizeof(T) * idx, (char *) &ele, sizeof(T));
		return ele;
	}

	/**
	 * This method copies data to a buffer.
	 * \param rel_off the offset relative to the beginning of the array.
	 * \param buf the buffer
	 * \param size the size of data copied to the buffer.
	 */
	void memcpy(off_t rel_off, char buf[], size_t size) const;

	static const page_byte_array &const_cast_ref(page_byte_array &arr) {
		return (const page_byte_array &) arr;
	}
};

#endif
