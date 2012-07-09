#ifndef __HASH_INDEX_CACHE_H__
#define __HASH_INDEX_CACHE_H__

#include <pthread.h>
#include <math.h>

/* the header files for the lock-free hashtable. */
extern "C" {
#include <map.h>
#include <hashtable.h>
#include <common.h>
}

#include "memory_manager.h"
#include "cache.h"
#include "concurrency.h"
#include "hashtable.h"
#include "SA_hash_table.h"
#include "SA_hash_table.cpp"

template<class KeyT, class ValueT>
class lock_free_hashtable: hashtable_interface<KeyT, ValueT>
{
	map_t *map;
public:
	lock_free_hashtable() {
		static const map_impl_t *map_type = { &MAP_IMPL_HT };
		map = map_alloc(map_type, NULL);
	}

	~lock_free_hashtable() {
		map_free(map);
	}

	ValueT get(KeyT key) {
		/*
		 * key == 0 isn't allowed,
		 * and 1 doesn't exist in the key space.
		 */
		if (key == 0)
			key = 1;
		return (ValueT) map_get(map, (map_key_t) key);
	}

	bool remove(KeyT key, ValueT value) {
		if (key == 0)
			key = 1;
		map_val_t ret = map_remove(map, (map_key_t) key);
		// TODO add the value back if 
		return ret == (map_val_t) value;;
	}

	ValueT putIfAbsent(KeyT key, ValueT value) {
		if (key == 0)
			key = 1;
		return (ValueT) map_cas(map, (map_key_t) key,
				0, (map_val_t) value);
	}

	bool replace(KeyT key, ValueT expect, ValueT new_value) {
		if (key == 0)
			key = 1;
		map_val_t ret = map_cas(map, (map_key_t) key,
				(map_val_t) expect, (map_val_t) new_value);
		return ret != DOES_NOT_EXIST;
	}
};

class frame: public thread_safe_page
{
	/* equivalent to the number of hits */
	atomic_integer wcount;
	/* equivalent to the number of references */
	atomic_integer pinning;
	frame *prev, *next;

public:
	frame(): thread_safe_page() {
		prev = next = this;
	}

	frame(off_t offset, char *data): thread_safe_page(offset, data) {
		prev = next = this;
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

	void incrWC() {
		wcount.inc(1);
	}

	int decrWC() {
		return wcount.dec(1);
	}

	bool tryEvict() {
		return pinning.CAS(1, -1);
	}

	void evictUnshared() {
		pinning.CAS(1, -1);
	}

	int pinCount() {
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

	/* for linked pages */
	void add_front(frame *pg) {
		frame *next = this->next;
		pg->next = next;
		pg->prev = this;
		this->next = pg;
		next->prev = pg;
	}

	void add_back(frame *pg) {
		frame *prev = this->prev;
		pg->next = this;
		pg->prev = prev;
		this->prev = pg;
		prev->next = pg;
	}

	void remove_from_list() {
		frame *prev = this->prev;
		frame *next = this->next;
		prev->next = next;
		next->prev = prev;
		this->next = this;
		this->prev = this;
	}

	bool is_empty() {
		return this->next == this;
	}

	frame *front() {
		return next;
	}

	frame *back() {
		return prev;
	}

};

class clock_buffer
{
	atomic_array<frame *> pool;
	atomic_integer free;
	atomic_integer clock_hand;
	int size;
public:
	clock_buffer(int size): pool(size), free(size) {
		this->size = size;
	}

	frame *add(frame *entry);
	frame *swap(frame *entry);

	void moveClockHand(int curr, int start) {
		int delta;
		if (curr < start)
			delta = curr + size - start + 1;
		else
			delta = curr - start + 1;
		clock_hand.inc(delta);
	}
};

class frame_allocator
{
	frame list_head;
	int size;
	pthread_spinlock_t lock;
	frame *p;
public:
	frame_allocator(int size) {
		pthread_spin_init(&lock, PTHREAD_PROCESS_PRIVATE);
		this->size = size;
		p = new frame[size];
		for (int i = 0; i < size; i++)
			list_head.add_front(&p[i]);
	}

	~frame_allocator() {
		delete [] p;
	}

	frame *alloc() {
		pthread_spin_lock(&lock);
		if (list_head.is_empty()) {
			pthread_spin_unlock(&lock);
			return NULL;
		}
		frame *p = list_head.front();
		p->remove_from_list();
		pthread_spin_unlock(&lock);
		return p;
	}

	void free(frame *p) {
		pthread_spin_lock(&lock);
		list_head.add_front(p);
		pthread_spin_unlock(&lock);
	}
};

class hash_index_cache: public page_cache
{
	// TODO these should be thread private.
	frame_allocator *allocator;
	memory_manager *manager;

	hashtable_interface<off_t, frame *> *hashtable;
	clock_buffer *clock_buf;
public:
	hash_index_cache(memory_manager *manager) {
		hashtable = new SA_hashtable<off_t, frame *>(1024);
//		hashtable = new lock_free_hashtable<off_t, frame *>();
		this->manager = manager;
		manager->register_cache(this);
		int max_npages = manager->get_max_size() / PAGE_SIZE;
		clock_buf = new clock_buffer(max_npages);
		/* we need more frames than the maximal number of pages. */
		allocator = new frame_allocator(max_npages * 2);
	}

	~hash_index_cache() {
		delete hashtable;
		delete clock_buf;
		manager->unregister_cache(this);
	}

	void purge_frame(frame *p) {
		char *pg = (char *) p->volatileGetValue();
		if (pg) {
			manager->free_pages(1, &pg);
		}
		allocator->free(p);
	}

	frame *addEntry(off_t offset, char *data);
	page *search(off_t offset, off_t &old_off);
};

#endif

