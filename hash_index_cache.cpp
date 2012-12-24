#include <errno.h>

#include "hash_index_cache.h"

frame *hash_index_cache::addEntry(off_t key, char *data) {
#ifdef PER_CPU
	gclock_buffer *gclock_buf = (gclock_buffer *) pthread_getspecific(gclock_key);
	frame_allocator *allocator = (frame_allocator *) pthread_getspecific(allocator_key);
#endif
	for (;;) {
		frame *new_entry = allocator->alloc();
		*new_entry = frame(key, data);
		/*
		 * all allocated frames will be added to the clock buffer.
		 * so if the frame can't be added to the hashtable, 
		 * it doesn't need to be freed. The frame will be freed
		 * when it is evicted from the clock buffer.
		 */
		frame *removed = gclock_buf->add(new_entry);
		/*
		 * If the gclock buffer doesn't return a page,
		 * it doesn't evict any pages, so we should increase
		 * the number of pages by one. It doesn't happen very
		 * often, so it shouldn't hurt performance.
		 */
		if (removed == NULL)
			num_pages.inc(1);
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
			char *pg = (char *) removed->volatileGetValue();
			if (pg) {
				/*
				 * We can give the page in the old frame to the new page.
				 * No one else can get access to the new frame, the operation
				 * below should be always successful.
				 */
				if (!new_entry->CASValue(NULL, pg)) {
					fprintf(stderr, "we can't set the value of the new frame!\n");
#ifdef PER_CPU
					memory_manager *manager
						= (memory_manager *) pthread_getspecific(manager_key);
#endif
					manager->free_pages(1, &pg);
				}
			}
			allocator->free(removed);
		}

		/* 
		 * We have to make sure the new entry gets a page before it is
		 * published in the index.
		 */
		if (new_entry->volatileGetValue() == NULL) {
			char *pg;
#ifdef PER_CPU
			memory_manager *manager = (memory_manager *) pthread_getspecific(manager_key);
#endif
			bool ret = manager->get_free_pages(1, &pg, this);
			if (!ret) {
				fprintf(stderr, "can't allocate a page from the memory manager.\n");
				fprintf(stderr, "there are %d pages allocated\n", num_pages.get());
				exit(1);
			}
			/*
			 * if the value of the frame has been set
			 * by another thread, free the allocated page.
			 */
			if (!new_entry->CASValue(NULL, pg)) {
				manager->free_pages(1, &pg);
			}
		}

		frame *prev_entry = hashtable->putIfAbsent(key / PAGE_SIZE, new_entry);
		if (prev_entry) {
			/* this happens if the page is just added to the hash table. */

			if (!prev_entry->pin()) {	// If we can't pin the old page
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
				/*
				 * This can happen if another thread happens to replace the old
				 * page successfully. In this case, we can't do anything and
				 * have to start over. The new page will be evicted by the 
				 * gclock algorithm later.
				 * This case should happen very rarely.
				 */
				new_entry->evictUnshared();
				continue;
			}
			/*
			 * Since we can pin the old page, we are just using the old one.
			 * So we should set the new page unused, and gclock algorithm
			 * will use evict later.
			 */
			new_entry->evictUnshared();
			prev_entry->incrWC();
			return prev_entry;
		}
		if (new_entry->pin()) {
			new_entry->incrWC();
			return new_entry;
		}
		else {
			// In a rare case, we can't pin the new frame.
			// we should remove it from the hash table and try again.
			hashtable->remove(key / PAGE_SIZE, new_entry);
		}
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
		return entry;
	}
}
