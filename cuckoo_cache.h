#ifndef __CUCKOO_HASH_H__
#define __CUCKOO_HASH_H__

#include <stdlib.h>
#include <math.h>

#include "cache.h"

#define MAXLOOP 5

#ifdef STATISTICS
volatile int removed_indices;
#endif

class cuckoo_hash
{
	page **tables[2];
	int log_size;
	long a[2];

	int hash(off_t key, int a_idx) {
		int v = ((a[a_idx] * key) & 0xFFFFFFFF) >> (32 - log_size);
		assert (v < (1 << log_size));
		return v;
	}
public:
	cuckoo_hash(int size) {
		/* cuckoo hash needs tables to be half-empty in order to be efficient. */
		tables[0] = (page **) calloc(size * 2, sizeof(page *));
		tables[1] = (page **) calloc(size * 2, sizeof(page *));
		log_size = log2(size);
		a[0] = random();
		a[1] = random();
#ifdef STATISTICS
		removed_indices = 0;
#endif
	}

	~cuckoo_hash() {
		free(tables[0]);
		free(tables[1]);
	}

	void insert(off_t key, page *value) {
		for (int i = 0; i < MAXLOOP; i++) {
			page *tmp;
			// TODO I need to test if this atomic operation works
			tmp = tables[0][hash(key, 0)];
			tables[0][hash(key, 0)] = value;
			value = tmp;
//			value = __sync_lock_test_and_set(&tables[0][hash(key, 0)], value);
			if (value == NULL || value->get_offset() == key)
				return;
			key = value->get_offset();
			tmp = tables[1][hash(key, 1)];
			tables[1][hash(key, 1)] = value;
			value = tmp;
//			value = __sync_lock_test_and_set(&tables[1][hash(key, 1)], value);
			if (value == NULL || value->get_offset() == key)
				return;
			key = value->get_offset();
		}
		/* 
		 * we don't need to rehash the table.
		 * just keep silent. The worst case is that the page can't be indexed
		 * at the moment, and it will be read again from the file. 
		 * So it doesn't hurt the correctness.
		 */
#ifdef STATISTICS
		__sync_fetch_and_add(&removed_indices, 1);
#endif
	}

	void remove(off_t key) {
		page *v = tables[0][hash(key, 0)];
		if (v && v->get_offset() == key) {
			tables[0][hash(key, 0)] = NULL;
			return;
		}

		/*
		 * if it's not in the first table,
		 * try the second one.
		 */
		v = tables[1][hash(key, 1)];
		if (v && v->get_offset() == key)
			tables[1][hash(key, 1)] = NULL;
	}

	page *search(off_t key) {
		page *v = tables[0][hash(key, 0)];
		/* find it in the first table. */
		if (v && v->get_offset() == key)
			return v;

		/* search the second table. */
		v = tables[1][hash(key, 1)];
		if (v && v->get_offset() != key)
			v = NULL;
		return v;
	}
};

class cuckoo_cache: public page_cache
{
//	page_buffer *bufs;
	page_buffer<thread_safe_page> *buf;
	cuckoo_hash table;
	pthread_spinlock_t _lock;
public:
	cuckoo_cache(long cache_size): table(cache_size / PAGE_SIZE) {
		long npages = cache_size / PAGE_SIZE;

//		/* each thread has a page buffer, and page eviction is done in the local thread. */
//		bufs = new page_buffer[nthreads](npages / nthreads);
		buf = new page_buffer<thread_safe_page>(npages);
		pthread_spin_init(&_lock, PTHREAD_PROCESS_PRIVATE);
	}

	~cuckoo_cache() {
		delete buf;
//		delete [] bufs;
		pthread_spin_destroy(&_lock);
	}

	page *search(off_t offset) {
		pthread_spin_lock(&_lock);
		thread_safe_page *pg = (thread_safe_page *) table.search(offset);
		if (pg == NULL) {
			pg = buf->get_empty_page();
			while (pg->get_ref()) {
				pthread_spin_unlock(&_lock);
				pg->wait_unused();
				pthread_spin_lock(&_lock);
			}
			/*
			 * the page has been used,
			 * it must be in the hash table.
			 */
			if (pg->get_offset() >= 0)
				table.remove(pg->get_offset());
			pg->set_data_ready(false);
			pg->set_offset(offset);
			table.insert(offset, pg);
		}
		pg->inc_ref();
		pthread_spin_unlock(&_lock);
		return pg;
	}
};

#endif
