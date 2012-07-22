#include <errno.h>

#include "hash_index_cache.h"

frame *hash_index_cache::addEntry(off_t key, char *data) {
	for (;;) {
		frame *new_entry = allocator->alloc();
		*new_entry = frame(key, data);
		/*
		 * all allocated frames will be added to the clock buffer.
		 * so if the frame can't be added to the hashtable, 
		 * it doesn't need to be freed. The frame will be freed
		 * when it is evicted from the clock buffer.
		 */
		frame *removed = clock_buf->add(new_entry);
		/* 
		 * the removed page can't be pinned any more,
		 * so other threads can't reference it except the thread
		 * that calls `putIfAbsent()'.
		 */
		if (removed) {
			/*
			 * this is a place that is different from the paper.
			 * remove may fail, but in any case, the page can be
			 * purged now.
			 * BTW, the only place that can cause remove to fail
			 * is replace() below.
			 */
			hashtable->remove(removed->get_offset() / PAGE_SIZE, removed);
			/* TODO There is a problem here. How can I guarantee that
			 * no other threads have a reference to this frame?
			 * It's possible a thread gets a reference to the page before
			 * the page is removed from the hashtable. */
			purge_frame(removed);
		}
		frame *prev_entry = hashtable->putIfAbsent(key / PAGE_SIZE, new_entry);
		if (prev_entry) {
			/* this happens if the page is just added to the hash table. */
			if (!prev_entry->pin()) {
				/*
				 * It's possible that a page has been evicted,
				 * but it couldn't be removed from the clock buffer.
				 * Therefore, I can't purge the page.
				 */
				if (hashtable->replace(key / PAGE_SIZE, prev_entry, new_entry)) {
					/*
					 * this place is different from the code
					 * in the paper. The value of the frame
					 * will be set later.
					 */
					return new_entry;
				}
				new_entry->evictUnshared();
				continue;
			}
			new_entry->evictUnshared();
			prev_entry->incrWC();
			return prev_entry;
		}
		return new_entry;
	}
}

page *hash_index_cache::search(off_t offset, off_t &old_off) {
	frame *entry = hashtable->get(offset / PAGE_SIZE);
	/*
	 * since the key of a frame nevers changes,
	 * the frame returned from the hash table doesn't
	 * have an old offset.
	 */
	old_off = 0;
	if (entry && entry->pin()) {
		entry->incrWC();
		return entry;
	}
	else {
		entry = addEntry(offset, NULL);
		/* 
		 * the frame might just have been added to the hashtable.
		 * In this case, its value should be NULL.
		 */
		if (entry->volatileGetValue() == NULL) {
			char *pg;
			manager->get_free_pages(1, &pg, this);
			/*
			 * if the value of the frame has been set
			 * by another thread, free the allocated page.
			 */
			if (!entry->CASValue(NULL, pg)) {
				manager->free_pages(1, &pg);
			}
		}
		return entry;
	}
}
