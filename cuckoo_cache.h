#ifndef __CUCKOO_HASH_H__
#define __CUCKOO_HASH_H__

#include <stdlib.h>
#include <math.h>

#include "cache.h"

#define MAXLOOP 2

#ifdef STATISTICS
volatile int removed_indices;
#endif

template <class T>
class lockable_pointer
{
	T *p;
public:
	lockable_pointer() {
		this->p = NULL;
	}

	lockable_pointer(T *p) {
		this->p = p;
		assert((((long) p) & 0x1L) == 0);
	}

	~lockable_pointer() {
	}

	T *operator->() {
		return (T *) (((long) p) & ~0x1L);
	}

	void lock() {
#ifdef DEBUG
		printf("thread %ld lock\n", pthread_self());
#endif
		while (__sync_fetch_and_or((long *) &p, 0x1L) & 0x1L) {}
	}

	void unlock() {
#ifdef DEBUG
		printf("thread %ld unlock\n", pthread_self());
#endif
		__sync_fetch_and_and((long *) &p, ~0x1L);
	}

	void set_pointer(T *p) {
		this->p = (T *) (((long) p) | (((long) this->p) & 0x1L));
	}
};

class cuckoo_hash
{
	lockable_pointer<thread_safe_page> *tables[2];
	int log_table_sizes[2];
	long a[2];

	int hash(off_t key, int idx) {
		int v = ((a[idx] * key) & 0xFFFFFFFF) >> (32 - log_table_sizes[idx]);
		assert (v < (1 << log_table_sizes[idx]));
		return v;
	}
public:
	cuckoo_hash(int size) {
		const int table_size0 = size * 4;
		const int table_size1 = size * 2;
		/* cuckoo hash needs tables to be half-empty in order to be efficient. */
		tables[0] = new lockable_pointer<thread_safe_page>[table_size0];
		tables[1] = new lockable_pointer<thread_safe_page>[table_size1];
		log_table_sizes[0] = log2(table_size0);
		log_table_sizes[1] = log2(table_size1);
		a[0] = random();
		a[1] = random();
#ifdef STATISTICS
		removed_indices = 0;
#endif
	}

	~cuckoo_hash() {
		delete [] tables[0];
		delete [] tables[1];
	}

	thread_safe_page *swap_entry(const int i, const off_t key,
			thread_safe_page *const value) {
		thread_safe_page *tmp;
		// TODO I need to test if this atomic operation works
		tables[i][hash(key, i)].lock();
		if(tables[i][hash(key, i)].operator->() == NULL) {
			tables[i][hash(key, i)].set_pointer(value);
			tables[i][hash(key, i)].unlock();
			return NULL;
		}
		else if(tables[i][hash(key, i)]->get_offset() == key) {
			tables[i][hash(key, i)].unlock();
			return NULL;
		}
		else {
			tmp = tables[i][hash(key, i)].operator->();
			tables[i][hash(key, i)].set_pointer(value);
			tables[i][hash(key, i)].unlock();
		}

		return tmp;
	}

	void insert(off_t key, thread_safe_page *value) {
		for (int i = 0; i < MAXLOOP; i++) {
			value = swap_entry(0, key, value);
			if (value == NULL)
				return;
			key = value->get_offset();
			value = swap_entry(1, key, value);
			if (value == NULL)
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

	bool remove_entry(const int i, const off_t key) {
		bool ret = false;
		tables[i][hash(key, i)].lock();
		thread_safe_page *v = tables[i][hash(key, i)].operator->();
		if (v && v->get_offset() == key) {
			tables[i][hash(key, i)].set_pointer(NULL);
			ret = true;
		}
		tables[i][hash(key, i)].unlock();
		return ret;
	}

	void remove(off_t key) {
		if(!remove_entry(0, key))
			remove_entry(1, key);
	}

	thread_safe_page *search_entry(const int i, const off_t key) {
		thread_safe_page *ret = NULL;
		tables[i][hash(key, i)].lock();
		thread_safe_page *v = tables[i][hash(key, i)].operator->();
		if (v && v->get_offset() == key) {
			ret = v;
			v->inc_ref();
		}
		tables[i][hash(key, i)].unlock();
		return ret;
	}

	thread_safe_page *search(off_t key) {
		thread_safe_page *v = search_entry(0, key);
		if (v)
			return v;
		return search_entry(1, key);
	}
};

class cuckoo_cache: public page_cache
{
//	page_buffer *bufs;
	page_buffer<thread_safe_page> *buf;
	cuckoo_hash table;
public:
	cuckoo_cache(long cache_size): table(cache_size / PAGE_SIZE) {
		long npages = cache_size / PAGE_SIZE;

//		/* each thread has a page buffer, and page eviction is done in the local thread. */
//		bufs = new page_buffer[nthreads](npages / nthreads);
		buf = new page_buffer<thread_safe_page>(npages, 0);
	}

	~cuckoo_cache() {
		delete buf;
//		delete [] bufs;
	}

	page *search(off_t offset) {
		thread_safe_page *pg = table.search(offset);
		if (pg == NULL) {
			pg = buf->get_empty_page();
			/*
			 * after this point, no one else can find the page
			 * in the table. but other threads might have the 
			 * reference to the page. therefore, we need to wait
			 * until all other threads have release their reference.
			 */
			table.remove(pg->get_offset());

			pg->wait_unused();
			// TODO I should put a barrier so that the thread
			// wait until the status chanage.
			/* at this point, no other threads are using the page. */
			pg->set_data_ready(false);
			pg->set_offset(offset);
			pg->inc_ref();
			table.insert(offset, pg);
		}
		return pg;
	}
};

#endif
