#include "cuckoo_cache.h"

#ifdef STATISTICS
volatile int removed_indices;
#endif

page_buffer<thread_safe_page> *lockable_pointer::buf;

thread_safe_page *cuckoo_hash::swap_entry(const int i, const off_t key,
		thread_safe_page *const value) {
	thread_safe_page *tmp;
	// TODO I need to test if this atomic operation works
	tables[i][hash(key, i)].lock();
	if(tables[i][hash(key, i)].get_pointer() == NULL) {
		tables[i][hash(key, i)].set_pointer(value);
		tables[i][hash(key, i)].unlock();
		return NULL;
	}
	else if(tables[i][hash(key, i)].get_pointer()->get_offset() == key) {
		tables[i][hash(key, i)].unlock();
		return NULL;
	}
	else {
		tmp = tables[i][hash(key, i)].get_pointer();
		tables[i][hash(key, i)].set_pointer(value);
		tables[i][hash(key, i)].unlock();
	}

	return tmp;
}

void cuckoo_hash::insert(off_t key, thread_safe_page *value) {
	for (int i = 0; i < MAXLOOP; i++) {
		value = swap_entry(0, key, value);
		if (value == NULL)
			return;
		value = swap_entry(1, key, value);
		if (value == NULL)
			return;
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

bool cuckoo_hash::remove_entry(const int i, const off_t key) {
	bool ret = false;
	tables[i][hash(key, i)].lock();
	thread_safe_page *v = tables[i][hash(key, i)].get_pointer();
	if (v && v->get_offset() == key) {
		tables[i][hash(key, i)].set_pointer(NULL);
		ret = true;
	}
	tables[i][hash(key, i)].unlock();
	return ret;
}

thread_safe_page *cuckoo_hash::search_entry(const int i, const off_t key) {
	thread_safe_page *ret = NULL;
	tables[i][hash(key, i)].lock();
	thread_safe_page *v = tables[i][hash(key, i)].get_pointer();
	if (v && v->get_offset() == key) {
		ret = v;
		v->inc_ref();
	}
	tables[i][hash(key, i)].unlock();
	return ret;
}

page *cuckoo_cache::search(off_t offset, off_t &old_off) {
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
