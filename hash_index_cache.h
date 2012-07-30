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
#include "gclock.h"

template<class KeyT, class ValueT>
class lock_free_hashtable: public hashtable_interface<KeyT, ValueT>
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

class frame_allocator
{
	linked_page_queue queue;
	int size;
	pthread_spinlock_t lock;
	frame *p;
public:
	frame_allocator(int size) {
		pthread_spin_init(&lock, PTHREAD_PROCESS_PRIVATE);
		this->size = size;
		p = new frame[size];
		for (int i = 0; i < size; i++)
			queue.push_back(&p[i]);
	}

	~frame_allocator() {
		delete [] p;
	}

	frame *alloc() {
		pthread_spin_lock(&lock);
		if (queue.empty()) {
			pthread_spin_unlock(&lock);
			return NULL;
		}
		frame *p = queue.front();
		queue.pop_front();
		pthread_spin_unlock(&lock);
		return p;
	}

	void free(frame *p) {
		pthread_spin_lock(&lock);
		queue.push_back(p);
		pthread_spin_unlock(&lock);
	}
};

class hash_index_cache: public page_cache
{
	// TODO these should be thread private.
	frame_allocator *allocator;
	memory_manager *manager;

	hashtable_interface<off_t, frame *> *hashtable;
	LF_gclock_buffer *clock_buf;
public:
	hash_index_cache(memory_manager *manager) {
		hashtable = new SA_hashtable<off_t, frame *>(1024);
//		hashtable = new lock_free_hashtable<off_t, frame *>();
		this->manager = manager;
		manager->register_cache(this);
		int max_npages = manager->get_max_size() / PAGE_SIZE;
		clock_buf = new LF_gclock_buffer(max_npages);
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

