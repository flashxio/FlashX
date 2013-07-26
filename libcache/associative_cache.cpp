#include <errno.h>
#include <limits.h>

#include "io_interface.h"
#include "associative_cache.h"
#include "flush_thread.h"
#include "container.cpp"
#include "exception.h"

#ifdef STATISTICS
volatile int avail_cells;
volatile int num_wait_unused;
volatile int lock_contentions;
#endif

const long default_init_cache_size = 128 * 1024 * 1024;

template<class T>
void page_cell<T>::set_pages(char *pages[], int num, int node_id)
{
	assert(num <= CELL_SIZE);
	for (int i = 0; i < num; i++) {
		buf[i] = T(-1, pages[i], node_id);
	}
	idx = 0;
	num_pages = num;
	for (int i = 0; i < num; i++) {
		maps[i] = i;
	}
}

template<class T>
void page_cell<T>::rebuild_map()
{
	int j = 0;
	for (int i = 0; i < CELL_SIZE; i++) {
		if (buf[i].get_data())
			maps[j++] = i;
	}
	assert(j == num_pages);
}

template<class T>
void page_cell<T>::add_pages(char *pages[], int num, int node_id)
{
	int num_added = 0;
	assert(num_pages == get_num_used_pages());
	assert(num + num_pages <= CELL_SIZE);
	for (int i = 0; i < CELL_SIZE && num_added < num; i++) {
		if (buf[i].get_data() == NULL)
			buf[i] = T(-1, pages[num_added++], node_id);
	}
	num_pages += num;
	rebuild_map();
}

template<class T>
void page_cell<T>::inject_pages(T pages[], int npages)
{
	int num_copied = 0;
	for (int i = 0; i < CELL_SIZE && num_copied < npages; i++) {
		if (buf[i].get_data() == NULL) {
			buf[i] = pages[num_copied++];
		}
	}
	assert(num_copied == npages);
	num_pages += num_copied;
	rebuild_map();
}

template<class T>
void page_cell<T>::steal_pages(T pages[], int &npages)
{
	int num_copied = 0;
	for (int i = 0; i < CELL_SIZE && num_copied < npages; i++) {
		if (buf[i].get_data()) {
			// We have to make sure the page isn't being referenced.
			// TODO busy wait.
			while (buf[i].get_ref() > 0) {}
			pages[num_copied++] = buf[i];
			buf[i] = T();
		}
	}
	npages = num_copied;
	num_pages -= num_copied;
	if (num_pages > 0)
		rebuild_map();
	else
		memset(maps, 0, sizeof(maps));
}

template<class T>
void page_cell<T>::sanity_check() const
{
	assert(params.get_SA_min_cell_size() <= num_pages);
	int num_used_pages = 0;
	for (int i = 0; i < CELL_SIZE; i++)
		if (buf[i].get_data())
			num_used_pages++;
	assert(num_used_pages == num_pages);
	int prev_map = -1;
	for (int i = 0; i < num_pages; i++) {
		int map = maps[i];
		if (prev_map >= 0)
			assert(map > prev_map);
		assert(buf[map].get_data());
		prev_map = map;
	}
}

template<class T>
int page_cell<T>::get_num_used_pages() const
{
	int num = 0;
	for (int i = 0; i < CELL_SIZE; i++)
		if (buf[i].get_data())
			num++;
	return num;
}

hash_cell::hash_cell(associative_cache *cache, long hash, bool get_pages) {
	this->hash = hash;
	assert(hash < INT_MAX);
	pthread_spin_init(&_lock, PTHREAD_PROCESS_PRIVATE);
	this->table = cache;
	if (get_pages) {
		char *pages[CELL_SIZE];
		if (!table->get_manager()->get_free_pages(params.get_SA_min_cell_size(),
					pages, cache))
			throw oom_exception();
		buf.set_pages(pages, params.get_SA_min_cell_size(), table->get_node_id());
	}
	num_accesses = 0;
	num_evictions = 0;
}

void hash_cell::sanity_check()
{
	pthread_spin_lock(&_lock);
	buf.sanity_check();
	pthread_spin_unlock(&_lock);
}

void hash_cell::add_pages(char *pages[], int num)
{
	buf.add_pages(pages, num, table->get_node_id());
}

int hash_cell::add_pages_to_min(char *pages[], int num)
{
	int num_required = CELL_MIN_NUM_PAGES - buf.get_num_pages();
	if (num_required > 0) {
		num_required = min(num_required, num);
		buf.add_pages(pages, num_required, table->get_node_id());
		return num_required;
	}
	else
		return 0;
}

void hash_cell::merge(hash_cell *cell)
{
	pthread_spin_lock(&_lock);
	pthread_spin_lock(&cell->_lock);

	assert(cell->get_num_pages() + this->get_num_pages() <= CELL_SIZE);
	thread_safe_page pages[CELL_SIZE];
	int npages = CELL_SIZE;
	// TODO there may be busy wait in this method.
	cell->buf.steal_pages(pages, npages);
	buf.inject_pages(pages, npages);

	pthread_spin_unlock(&cell->_lock);
	pthread_spin_unlock(&_lock);
}

/**
 * rehash the pages in the current cell to the expanded cell.
 */
void hash_cell::rehash(hash_cell *expanded)
{
	pthread_spin_lock(&_lock);
	pthread_spin_lock(&expanded->_lock);
	thread_safe_page *exchanged_pages_pointers[CELL_SIZE];
	int num_exchanges = 0;
	for (unsigned int i = 0; i < buf.get_num_pages(); i++) {
		thread_safe_page *pg = buf.get_page(i);
		int hash1 = table->hash1_locked(pg->get_offset());
		/*
		 * It's possible that a page is in a wrong cell.
		 * It's likely because the page is added to the cell 
		 * right when `level' is increased.
		 * But the case is rare, so we can just simple ignore
		 * the case. It doesn't affect the correctness of 
		 * the implementation. The only penalty is that
		 * we might get a cache miss.
		 * Since the page is in a wrong cell, it won't be 
		 * accessed any more, so we should shorten the time
		 * it gets evicted by setting its hit to 1.
		 */
		if (hash1 != expanded->hash) {
			pg->set_hits(1);
			continue;
		}
		/* 
		 * if the two hash values don't match,
		 * it means the page is mapped to the expanded cell.
		 * we exchange the pages in the two cells.
		 */
		if (this->hash != hash1) {
			/* 
			 * we have to make sure no other threads are using them
			 * before we can exchange them.
			 * If the pages are in use, skip them.
			 */
			if (pg->get_ref())
				continue;

			exchanged_pages_pointers[num_exchanges++] = pg;
			// We can't steal pages while iterating them.
		}
	}
	if (num_exchanges > 0) {
		// We can only steal pages here.
		thread_safe_page exchanged_pages[CELL_SIZE];
		for (int i = 0; i < num_exchanges; i++) {
			exchanged_pages[i] = *exchanged_pages_pointers[i];
			buf.steal_page(exchanged_pages_pointers[i], false);
			*exchanged_pages_pointers[i] = thread_safe_page();
		}
		buf.rebuild_map();
		expanded->buf.inject_pages(exchanged_pages, num_exchanges);
	}

	// Move empty pages to the expanded cell if it doesn't have enough pages.
	int num_required = params.get_SA_min_cell_size() - expanded->buf.get_num_pages();
	int num_empty = 0;
	if (num_required > 0) {
		thread_safe_page *empty_pages_pointers[num_required];
		thread_safe_page *empty_pages = new thread_safe_page[params.get_SA_min_cell_size()];
		for (unsigned int i = 0; i < buf.get_num_pages()
				&& num_empty < num_required; i++) {
			thread_safe_page *pg = buf.get_page(i);
			if (!pg->initialized())
				empty_pages_pointers[num_empty++] = pg;
		}
		for (int i = 0; i < num_empty; i++) {
			// For the same reason, we can't steal pages
			// while iterating them.
			empty_pages[i] = *empty_pages_pointers[i];
			buf.steal_page(empty_pages_pointers[i], false);
			*empty_pages_pointers[i] = thread_safe_page();
		}
		buf.rebuild_map();
		expanded->buf.inject_pages(empty_pages, num_empty);
		delete [] empty_pages;
	}
	pthread_spin_unlock(&expanded->_lock);
	pthread_spin_unlock(&_lock);
}

void hash_cell::steal_pages(char *pages[], int &npages)
{
	int num_stolen = 0;
	while (num_stolen < npages) {
		thread_safe_page *pg = get_empty_page();
		if (pg == NULL)
			break;
		assert(!pg->is_dirty());
		pages[num_stolen++] = (char *) pg->get_data();
		*pg = thread_safe_page();
		buf.steal_page(pg, false);
	}
	buf.rebuild_map();
	npages = num_stolen;
}

void hash_cell::rebalance(hash_cell *cell)
{
	// TODO
}

page *hash_cell::search(off_t offset)
{
	pthread_spin_lock(&_lock);
	page *ret = NULL;
	for (unsigned int i = 0; i < buf.get_num_pages(); i++) {
		if (buf.get_page(i)->get_offset() == offset) {
			ret = buf.get_page(i);
			break;
		}
	}
	if (ret) {
		if (ret->get_hits() == 0xff)
			buf.scale_down_hits();
		ret->inc_ref();
		ret->hit();
	}
	pthread_spin_unlock(&_lock);
	return ret;
}

/**
 * search for a page with the offset.
 * If the page doesn't exist, return an empty page.
 */
page *hash_cell::search(off_t off, off_t &old_off) {
	thread_safe_page *ret = NULL;
#ifndef STATISTICS
	pthread_spin_lock(&_lock);
#else
	if (pthread_spin_trylock(&_lock) == EBUSY) {
		__sync_fetch_and_add(&lock_contentions, 1);
		pthread_spin_lock(&_lock);
	}
#endif
	num_accesses++;

	for (unsigned int i = 0; i < buf.get_num_pages(); i++) {
		if (buf.get_page(i)->get_offset() == off) {
			ret = buf.get_page(i);
			break;
		}
	}
	if (ret == NULL) {
		num_evictions++;
		ret = get_empty_page();
		// We need to clear flags here.
		ret->set_data_ready(false);
		assert(!ret->is_io_pending());
		ret->set_prepare_writeback(false);
		if (ret->is_dirty() && !ret->is_old_dirty()) {
			ret->set_dirty(false);
			ret->set_old_dirty(true);
		}
		old_off = ret->get_offset();
		if (old_off == -1)
			old_off = PAGE_INVALID_OFFSET;
		/*
		 * I have to change the offset in the spinlock,
		 * to make sure when the spinlock is unlocked, 
		 * the page can be seen by others even though
		 * it might not have data ready.
		 */
		ret->set_offset(off);
#ifdef USE_SHADOW_PAGE
		shadow_page shadow_pg = shadow.search(off);
		/*
		 * if the page has been seen before,
		 * we should set the hits info.
		 */
		if (shadow_pg.is_valid())
			ret->set_hits(shadow_pg.get_hits());
#endif
	}
	else
		policy.access_page(ret, buf);
	/* it's possible that the data in the page isn't ready */
	ret->inc_ref();
	if (ret->get_hits() == 0xff) {
		buf.scale_down_hits();
#ifdef USE_SHADOW_PAGE
		shadow.scale_down_hits();
#endif
	}
	ret->hit();
	pthread_spin_unlock(&_lock);
	return ret;
}

/* this function has to be called with lock held */
thread_safe_page *hash_cell::get_empty_page() {
	thread_safe_page *ret = NULL;

search_again:
	ret = policy.evict_page(buf);
	if (ret == NULL) {
#ifdef DEBUG
		printf("all pages in the cell were all referenced\n");
#endif
		/* 
		 * If all pages in the cell are referenced, there is
		 * nothing we can do but wait. However, before busy waiting,
		 * we should unlock the lock, so other threads may still
		 * search the cell.
		 */
		pthread_spin_unlock(&_lock);
		bool all_referenced = true;
		while (all_referenced) {
			for (unsigned int i = 0; i < buf.get_num_pages(); i++) {
				thread_safe_page *pg = buf.get_page(i);
				/* If a page isn't referenced. */
				if (!pg->get_ref()) {
					all_referenced = false;
					break;
				}
			}
		}
		pthread_spin_lock(&_lock);
		goto search_again;
	}

	/* we record the hit info of the page in the shadow cell. */
#ifdef USE_SHADOW_PAGE
	if (ret->get_hits() > 0)
		shadow.add(shadow_page(*ret));
#endif

	return ret;
}

/* 
 * the end of the vector points to the pages
 * that are most recently accessed.
 */
thread_safe_page *LRU_eviction_policy::evict_page(
		page_cell<thread_safe_page> &buf)
{
	int pos;
	if (pos_vec.size() < buf.get_num_pages()) {
		pos = pos_vec.size();
	}
	else {
		/* evict the first page */
		pos = pos_vec[0];
		pos_vec.erase(pos_vec.begin());
	}
	thread_safe_page *ret = buf.get_page(pos);
	while (ret->get_ref()) {}
	pos_vec.push_back(pos);
	ret->set_data_ready(false);
	return ret;
}

void LRU_eviction_policy::access_page(thread_safe_page *pg,
		page_cell<thread_safe_page> &buf)
{
	/* move the page to the end of the pos vector. */
	int pos = buf.get_idx(pg);
	for (std::vector<int>::iterator it = pos_vec.begin();
			it != pos_vec.end(); it++) {
		if (*it == pos) {
			pos_vec.erase(it);
			break;
		}
	}
	pos_vec.push_back(pos);
}

thread_safe_page *LFU_eviction_policy::evict_page(
		page_cell<thread_safe_page> &buf)
{
	thread_safe_page *ret = NULL;
	int min_hits = 0x7fffffff;
	do {
		unsigned int num_io_pending = 0;
		for (unsigned int i = 0; i < buf.get_num_pages(); i++) {
			thread_safe_page *pg = buf.get_page(i);
			if (pg->get_ref()) {
				if (pg->is_io_pending())
					num_io_pending++;
				continue;
			}

			/* 
			 * refcnt only increases within the lock of the cell,
			 * so if the page's refcnt is 0 above,
			 * it'll be always 0 within the lock.
			 */

			if (min_hits > pg->get_hits()) {
				min_hits = pg->get_hits();
				ret = pg;
			}

			/* 
			 * if a page hasn't been accessed before,
			 * it's a completely new page, just use it.
			 */
			if (min_hits == 0)
				break;
		}
		if (num_io_pending == buf.get_num_pages()) {
			printf("all pages are at io pending\n");
			// TODO do something...
			// maybe we should use pthread_wait
		}
		/* it happens when all pages in the cell is used currently. */
	} while (ret == NULL);
	ret->set_data_ready(false);
	ret->reset_hits();
	return ret;
}

thread_safe_page *FIFO_eviction_policy::evict_page(
		page_cell<thread_safe_page> &buf)
{
	thread_safe_page *ret = buf.get_empty_page();
	/*
	 * This happens a lot if we actually read pages from the disk.
	 * So basically, we shouldn't use this eviction policy for SSDs
	 * or magnetic hard drive..
	 */
	while (ret->get_ref()) {
		ret = buf.get_empty_page();
	}
	ret->set_data_ready(false);
	return ret;
}

thread_safe_page *gclock_eviction_policy::evict_page(
		page_cell<thread_safe_page> &buf)
{
	thread_safe_page *ret = NULL;
	unsigned int num_referenced = 0;
	unsigned int num_dirty = 0;
	bool avoid_dirty = true;
	do {
		thread_safe_page *pg = buf.get_page(clock_head % buf.get_num_pages());
		if (num_dirty + num_referenced >= buf.get_num_pages()) {
			num_dirty = 0;
			num_referenced = 0;
			avoid_dirty = false;
		}
		if (pg->get_ref()) {
			num_referenced++;
			clock_head++;
			/*
			 * If all pages in the cell are referenced, we should
			 * return NULL to notify the invoker.
			 */
			if (num_referenced >= buf.get_num_pages())
				return NULL;
			continue;
		}
		if (avoid_dirty && pg->is_dirty()) {
			num_dirty++;
			clock_head++;
			continue;
		}
		if (pg->get_hits() == 0) {
			ret = pg;
			break;
		}
		pg->set_hits(pg->get_hits() - 1);
		clock_head++;
	} while (ret == NULL);
	ret->set_data_ready(false);
	return ret;
}

int gclock_eviction_policy::predict_evicted_pages(
		page_cell<thread_safe_page> &buf, int num_pages, int set_flags,
		int clear_flags, std::map<off_t, thread_safe_page *> &pages)
{
	// We are just predicting. We don't actually evict any pages.
	// So we need to make a copy of the hits of each page.
	short hits[CELL_SIZE];
	for (int i = 0; i < (int) buf.get_num_pages(); i++) {
		hits[i] = buf.get_page(i)->get_hits();
	}
	// The function returns when we get the expected number of pages
	// or we get pages 
	while (true) {
		int num_with_hits = 0;
		for (int i = 0; i < (int) buf.get_num_pages(); i++) {
			int idx = (i + clock_head) % buf.get_num_pages();
			short *hit = &hits[idx];
			// The page is already in the page map.
			if (*hit < 0)
				continue;
			else if (*hit == 0) {
				*hit = -1;
				thread_safe_page *p = buf.get_page(idx);
				if (p->test_flags(set_flags) && !p->test_flags(clear_flags)) {
					pages.insert(std::pair<off_t, thread_safe_page *>(
								p->get_offset(), p));
					if ((int) pages.size() == num_pages)
						return pages.size();
				}
			}
			else {
				(*hit)--;
				num_with_hits++;
			}
		}
		// No pages have hits.
		if (num_with_hits == 0)
			return pages.size();
	}
}

thread_safe_page *clock_eviction_policy::evict_page(
		page_cell<thread_safe_page> &buf)
{
	thread_safe_page *ret = NULL;
	unsigned int num_referenced = 0;
	unsigned int num_dirty = 0;
	bool avoid_dirty = true;
	do {
		thread_safe_page *pg = buf.get_page(clock_head % buf.get_num_pages());
		if (num_dirty + num_referenced >= buf.get_num_pages()) {
			num_dirty = 0;
			num_referenced = 0;
			avoid_dirty = false;
		}
		if (pg->get_ref()) {
			num_referenced++;
			if (num_referenced >= buf.get_num_pages())
				return NULL;
			clock_head++;
			continue;
		}
		if (avoid_dirty && pg->is_dirty()) {
			num_dirty++;
			clock_head++;
			continue;
		}
		if (pg->get_hits() == 0) {
			ret = pg;
			break;
		}
		pg->reset_hits();
		clock_head++;
	} while (ret == NULL);
	ret->set_data_ready(false);
	ret->reset_hits();
	return ret;
}

bool associative_cache::shrink(int npages, char *pages[])
{
	if (flags.set_flag(TABLE_EXPANDING)) {
		/*
		 * if the flag has been set before,
		 * it means another thread is expanding the table,
		 */
		return false;
	}

	/* starting from this point, only one thred can be here. */

	int pg_idx = 0;
	int orig_ncells = get_num_cells();
	while (pg_idx < npages) {
		// The cell table isn't in the stage of splitting.
		if (split == 0) {
			hash_cell *cell = get_cell(expand_cell_idx);
			while (height >= params.get_SA_min_cell_size()) {
				int num = max(0, cell->get_num_pages() - height);
				num = min(npages - pg_idx, num);
				if (num > 0) {
					cell->steal_pages(&pages[pg_idx], num);
					pg_idx += num;
				}

				if (expand_cell_idx <= 0) {
					height--;
					expand_cell_idx = orig_ncells;
				}
				expand_cell_idx--;
				cell = get_cell(expand_cell_idx);
			}
			if (pg_idx == npages) {
				cache_npages.dec(npages);
				flags.clear_flag(TABLE_EXPANDING);
				return true;
			}
		}

		/* From here, we shrink the cell table. */

		// When the thread is within in the while loop, other threads can
		// hardly access the cells in the table.
		if (level == 0)
			break;
		int num_half = (1 << level) * init_ncells / 2;
		table_lock.write_lock();
		if (split == 0) {
			split = num_half - 1;
			level--;
		}
		table_lock.write_unlock();
		while (split > 0) {
			hash_cell *high_cell = get_cell(split + num_half);
			hash_cell *cell = get_cell(split);
			// At this point, the high cell and the low cell together
			// should have no more than CELL_MIN_NUM_PAGES pages.
			cell->merge(high_cell);
			table_lock.write_lock();
			split--;
			table_lock.write_unlock();
		}
		int orig_narrays = (1 << level);
		// It's impossible to access the arrays after `narrays' now.
		int narrays = orig_narrays / 2;
		for (int i = narrays; i < orig_narrays; i++) {
			delete [] cells_table[i];
			cells_table[i] = NULL;
		}
	}
	flags.clear_flag(TABLE_EXPANDING);
	cache_npages.dec(npages);
	return true;
}

/**
 * This method increases the cache size by `npages'.
 */
int associative_cache::expand(int npages)
{
	if (flags.set_flag(TABLE_EXPANDING)) {
		/*
		 * if the flag has been set before,
		 * it means another thread is expanding the table,
		 */
		return 0;
	}

	/* starting from this point, only one thred can be here. */

	char *pages[npages];
	if (!manager->get_free_pages(npages, pages, this)) {
		flags.clear_flag(TABLE_EXPANDING);
		fprintf(stderr, "expand: can't allocate %d pages\n", npages);
		return 0;
	}
	int pg_idx = 0;
	bool expand_over = false;
	while (pg_idx < npages && !expand_over) {
		// The cell table isn't in the stage of splitting.
		if (split == 0) {
			int orig_ncells = get_num_cells();
			/* We first try to add pages to the existing cells. */
			hash_cell *cell = get_cell(expand_cell_idx);
			while (height <= CELL_SIZE && pg_idx < npages) {
				int num_missing;

				assert(pages[pg_idx]);
				// We should skip the cells with more than `height'.
				if (cell->get_num_pages() >= height)
					goto next_cell;

				num_missing = height - cell->get_num_pages();
				num_missing = min(num_missing, npages - pg_idx);
				cell->add_pages(&pages[pg_idx], num_missing);
				pg_idx += num_missing;
next_cell:
				expand_cell_idx++;
				if (expand_cell_idx >= (unsigned) orig_ncells) {
					expand_cell_idx = 0;
					height++;
				}
				cell = get_cell(expand_cell_idx);
			}
			if (pg_idx == npages) {
				cache_npages.inc(npages);
				flags.clear_flag(TABLE_EXPANDING);
				return npages;
			}

			/* We have to expand the cell table in order to add more pages. */

			/* Double the size of the cell table. */
			/* create cells and put them in a temporary table. */
			std::vector<hash_cell *> table;
			int orig_narrays = (1 << level);
			for (int i = orig_narrays; i < orig_narrays * 2; i++) {
				hash_cell *cells = new hash_cell[init_ncells];
				printf("create %d cells: %p\n", init_ncells, cells);
				for (int j = 0; j < init_ncells; j++) {
					cells[j] = hash_cell(this, i * init_ncells + j, false);
				}
				table.push_back(cells);
			}
			/*
			 * here we need to hold the lock because other threads
			 * might be accessing the table. by using the write lock,
			 * we notify others the table has been changed.
			 */
			table_lock.write_lock();
			for (unsigned int i = 0; i < table.size(); i++) {
				cells_table[orig_narrays + i] = table[i];
			}
			table_lock.write_unlock();
		}
		height = params.get_SA_min_cell_size() + 1;

		// When the thread is within in the while loop, other threads
		// can hardly access the cells in the table.
		int num_half = (1 << level) * init_ncells;
		while (split < num_half) {
			hash_cell *expanded_cell = get_cell(split + num_half);
			hash_cell *cell = get_cell(split);
			cell->rehash(expanded_cell);

			/*
			 * After rehashing, there is no guarantee that two cells will have
			 * the same number of pages. We need to either add empty pages to
			 * the cell without enough pages or rebalance the two cells.
			 */

			/* Add pages to the cell without enough pages. */
			int num_required = max(expanded_cell->get_num_pages()
				- params.get_SA_min_cell_size(), 0);
			num_required += max(cell->get_num_pages() - params.get_SA_min_cell_size(), 0);
			if (num_required <= npages - pg_idx) {
				/* 
				 * Actually only one cell requires more pages, the other
				 * one will just take 0 pages.
				 */
				pg_idx += cell->add_pages_to_min(&pages[pg_idx],
						npages - pg_idx);
				pg_idx += expanded_cell->add_pages_to_min(&pages[pg_idx],
						npages - pg_idx);
			}

			if (expanded_cell->get_num_pages() < params.get_SA_min_cell_size()
					|| cell->get_num_pages() < params.get_SA_min_cell_size()) {
				// If we failed to split a cell, we should merge the two half
				cell->merge(expanded_cell);
				expand_over = true;
				fprintf(stderr, "A cell can't have enough pages, merge back\n");
				break;
			}

			table_lock.write_lock();
			split++;
			table_lock.write_unlock();
		}
		table_lock.write_lock();
		if (split == num_half) {
			split = 0;
			level++;
		}
		table_lock.write_unlock();
	}
	if (pg_idx < npages)
		manager->free_pages(npages - pg_idx, &pages[pg_idx]);
	flags.clear_flag(TABLE_EXPANDING);
	cache_npages.inc(npages);
	return npages - pg_idx;
}

page *associative_cache::search(off_t offset, off_t &old_off) {
	/*
	 * search might change the structure of the cell,
	 * and cause the cell table to expand.
	 * Thus, the page might not be placed in the cell
	 * we found before. Therefore, we need to research
	 * for the cell.
	 */
	do {
		return get_cell_offset(offset)->search(offset, old_off);
	} while (true);
}

page *associative_cache::search(off_t offset)
{
	do {
		return get_cell_offset(offset)->search(offset);
	} while (true);
}

int associative_cache::get_num_used_pages() const
{
	unsigned long count;
	int npages = 0;
	do {
		table_lock.read_lock(count);
		int ncells = get_num_cells();
		for (int i = 0; i < ncells; i++) {
			npages += get_cell(i)->get_num_pages();
		}
	} while (!table_lock.read_unlock(count));
	return npages;
}

void associative_cache::sanity_check() const
{
	unsigned long count;
	do {
		table_lock.read_lock(count);
		int ncells = get_num_cells();
		for (int i = 0; i < ncells; i++) {
			hash_cell *cell = get_cell(i);
			cell->sanity_check();
		}
	} while (!table_lock.read_unlock(count));
}

associative_cache::associative_cache(long cache_size, long max_cache_size,
		int node_id, int offset_factor, bool expandable)
{
	this->offset_factor = offset_factor;
	pthread_mutex_init(&init_mutex, NULL);
	printf("associative cache is created on node %d, cache size: %ld, min cell size: %d\n",
			node_id, cache_size, params.get_SA_min_cell_size());
	this->node_id = node_id;
	level = 0;
	split = 0;
	height = params.get_SA_min_cell_size();
	expand_cell_idx = 0;
	this->expandable = expandable;
	this->manager = new memory_manager(max_cache_size, node_id);
	manager->register_cache(this);
	long init_cache_size = default_init_cache_size;
	if (init_cache_size > cache_size
			// If the cache isn't expandable, let's just use the maximal
			// cache size at the beginning.
			|| !expandable)
		init_cache_size = cache_size;
	int min_cell_size = params.get_SA_min_cell_size();
	if (init_cache_size < min_cell_size * PAGE_SIZE)
		init_cache_size = min_cell_size * PAGE_SIZE;
	int npages = init_cache_size / PAGE_SIZE;
	init_ncells = npages / min_cell_size;
	hash_cell *cells = new hash_cell[init_ncells];
	printf("%d cells: %p\n", init_ncells, cells);
	int max_npages = manager->get_max_size() / PAGE_SIZE;
	try {
		for (int i = 0; i < init_ncells; i++)
			cells[i] = hash_cell(this, i, true);
	} catch (oom_exception e) {
		fprintf(stderr,
				"out of memory: max npages: %d, init npages: %d\n",
				max_npages, npages);
		exit(1);
	}

	cells_table.push_back(cells);

	int max_ncells = max_npages / min_cell_size;
	for (int i = 1; i < max_ncells / init_ncells; i++)
		cells_table.push_back(NULL);

	if (expandable && cache_size > init_cache_size)
		expand((cache_size - init_cache_size) / PAGE_SIZE);
}

/**
 * This defines an interface for implementing the policy of selecting
 * dirty pages for flushing.
 */
class select_dirty_pages_policy
{
public:
	// It select a specified number of pages from the page set.
	virtual int select(hash_cell *cell, int num_pages,
			std::map<off_t, thread_safe_page *> &pages) = 0;
};

/**
 * The policy selects dirty pages that are most likely to be evicted
 * by the eviction policy.
 */
class eviction_select_dirty_pages_policy: public select_dirty_pages_policy
{
public:
	int select(hash_cell *cell, int num_pages,
			std::map<off_t, thread_safe_page *> &pages) {
		char set_flags = 0;
		char clear_flags = 0;
		page_set_flag(set_flags, DIRTY_BIT, true);
		page_set_flag(clear_flags, IO_PENDING_BIT, true);
		page_set_flag(clear_flags, PREPARE_WRITEBACK, true);
		cell->predict_evicted_pages(num_pages, set_flags, clear_flags, pages);
		return pages.size();
	}
};

/**
 * This policy simply selects some dirty pages in a page set.
 */
class default_select_dirty_pages_policy: public select_dirty_pages_policy
{
public:
	int select(hash_cell *cell, int num_pages,
			std::map<off_t, thread_safe_page *> &pages) {
		char set_flags = 0;
		char clear_flags = 0;
		page_set_flag(set_flags, DIRTY_BIT, true);
		page_set_flag(clear_flags, IO_PENDING_BIT, true);
		page_set_flag(clear_flags, PREPARE_WRITEBACK, true);
		cell->get_pages(num_pages, set_flags, clear_flags, pages);
		return pages.size();
	}
};

class associative_flush_thread: public flush_thread
{
	// For the case of NUMA cache, cache and local_cache are different.
	page_cache *cache;
	associative_cache *local_cache;

	io_interface *io;
	thread_safe_FIFO_queue<hash_cell *> dirty_cells;
	select_dirty_pages_policy *policy;
public:
	associative_flush_thread(page_cache *cache, associative_cache *local_cache,
			io_interface *io, int node_id): flush_thread(node_id), dirty_cells(
				MAX_NUM_DIRTY_CELLS_IN_QUEUE) {
		this->cache = cache;
		this->local_cache = local_cache;
		if (this->cache == NULL)
			this->cache = local_cache;

		this->io = io;
		policy = new eviction_select_dirty_pages_policy();
	}

	void run();
	void request_callback(io_request &req);
	void dirty_pages(thread_safe_page *pages[], int num);
	int flush_cell(hash_cell *cell, io_request *req_array, int req_array_size);
};

void associative_flush_thread::request_callback(io_request &req)
{
	if (req.get_num_bufs() == 1) {
		thread_safe_page *p = (thread_safe_page *) req.get_page(0);
		p->lock();
		assert(p->is_dirty());
		p->set_dirty(false);
		p->set_io_pending(false);
		p->dec_ref();
		p->unlock();
	}
	else {
		off_t off = req.get_offset();
		for (int i = 0; i < req.get_num_bufs(); i++) {
			thread_safe_page *p = req.get_page(i);
			assert(p);
			p->lock();
			assert(p->is_dirty());
			p->set_dirty(false);
			p->set_io_pending(false);
			p->dec_ref();
			assert(p->get_ref() >= 0);
			p->unlock();
			off += PAGE_SIZE;
		}
	}
}

void merge_pages2req(io_request &req, page_cache *cache);

int associative_flush_thread::flush_cell(hash_cell *cell,
		io_request *req_array, int req_array_size)
{
	std::map<off_t, thread_safe_page *> dirty_pages;
	policy->select(cell, NUM_WRITEBACK_DIRTY_PAGES, dirty_pages);
	int num_init_reqs = 0;
	for (std::map<off_t, thread_safe_page *>::const_iterator it
			= dirty_pages.begin(); it != dirty_pages.end(); it++) {
		thread_safe_page *p = it->second;
		p->lock();
		assert(!p->is_old_dirty());
		assert(p->data_ready());

		assert(num_init_reqs < req_array_size);
		// Here we flush dirty pages with normal requests.
		if (!p->is_io_pending() && !p->is_prepare_writeback()
				// The page may have been cleaned.
				&& p->is_dirty()) {
			// TODO global_cached_io may delete the extension.
			// I'll fix it later.
			if (!req_array[num_init_reqs].is_extended_req()) {
				io_request tmp(true);
				req_array[num_init_reqs] = tmp;
			}
			req_array[num_init_reqs].init(p->get_offset(), io, WRITE,
					get_node_id(), NULL, cache, NULL);
			req_array[num_init_reqs].add_page(p);
			req_array[num_init_reqs].set_high_prio(true);
			p->set_io_pending(true);
			p->unlock();
			merge_pages2req(req_array[num_init_reqs], cache);
			num_init_reqs++;
		}
		else
			p->unlock();

		// The code blow flush dirty pages with low-priority requests.
#if 0
		if (!p->is_io_pending() && !p->is_prepare_writeback()
				// The page may have been cleaned.
				&& p->is_dirty()) {
			// TODO global_cached_io may delete the extension.
			// I'll fix it later.
			if (!req_array[num_init_reqs].is_extended_req()) {
				io_request tmp(true);
				req_array[num_init_reqs] = tmp;
			}
			req_array[num_init_reqs].init(p->get_offset(), io, WRITE,
					get_node_id(), NULL, cache, NULL);
			req_array[num_init_reqs].add_page(p);
			req_array[num_init_reqs].set_high_prio(false);
			num_init_reqs++;
			p->set_prepare_writeback(true);
		}
		// When a page is put in the queue for writing back,
		// the queue of the IO thread doesn't own the page, which
		// means that the page can be evicted.
		p->unlock();
		p->dec_ref();
#endif
	}
	return num_init_reqs;
}

void associative_flush_thread::run()
{
	int num_fetches;
	// We can't get more requests than the number of pages in a cell.
	const int req_array_size = NUM_WRITEBACK_DIRTY_PAGES * 100;
	io_request req_array[req_array_size];
	while ((num_fetches = dirty_cells.get_num_entries()) > 0) {
		hash_cell *cells[num_fetches];
		int ret = dirty_cells.fetch(cells, num_fetches);
		// This is the only place where we fetches entries in the queue,
		// and there is only one thread fetching entries, so we can be 
		// very sure we can fetch the number of entries we specify.
		assert(ret == num_fetches);

		int num_init_reqs = 0;
		for (int i = 0; i < num_fetches; i++) {
			int ret = flush_cell(cells[i], &req_array[num_init_reqs],
					req_array_size - num_init_reqs);
			num_init_reqs += ret;
			if (num_init_reqs >= req_array_size - NUM_WRITEBACK_DIRTY_PAGES) {
				io->access(req_array, num_init_reqs);
				num_init_reqs = 0;
			}
			// We can clear the in_queue flag now.
			// The cell won't be added to the queue for flush until its dirty pages
			// have been written back successfully.
			// A cell is added to the queue only when the number of dirty pages
			// that aren't being written back is larger than a threshold.
			cells[i]->set_in_queue(false);
		}

		io->access(req_array, num_init_reqs);
	}
}

flush_thread *associative_cache::create_flush_thread(io_interface *io,
		page_cache *global_cache)
{
	pthread_mutex_lock(&init_mutex);
	if (_flush_thread == NULL && io) {
		_flush_thread = new associative_flush_thread(global_cache, this,
				io->clone(), node_id);
		_flush_thread->start();
	}
	pthread_mutex_unlock(&init_mutex);
	return _flush_thread;
}

void associative_cache::mark_dirty_pages(thread_safe_page *pages[], int num)
{
	if (_flush_thread)
		_flush_thread->dirty_pages(pages, num);
}

void associative_cache::flush_callback(io_request &req)
{
	if (_flush_thread)
		_flush_thread->request_callback(req);
}

void associative_cache::init(io_interface *underlying)
{
	create_flush_thread(underlying, this);
}

hash_cell *associative_cache::get_prev_cell(hash_cell *cell) {
	long index = cell->get_hash();
	// The first cell in the hash table.
	if (index == 0)
		return NULL;
	// The cell is in the middle of a cell array.
	if (index % init_ncells)
		return cell - 1;
	else {
		unsigned i;
		for (i = 0; i < cells_table.size(); i++) {
			if (cell == cells_table[i]) {
				assert(i > 0);
				// return the last cell in the previous cell array.
				return (cells_table[i - 1] + init_ncells - 1);
			}
		}
		// we should reach here if the cell exists in the table.
		abort();
	}
}

hash_cell *associative_cache::get_next_cell(hash_cell *cell)
{
	long index = cell->get_hash();
	// If it's not the last cell in the cell array.
	if (index % init_ncells != init_ncells - 1)
		return cell + 1;
	else {
		unsigned i;
		hash_cell *first = cell + 1 - init_ncells;
		for (i = 0; i < cells_table.size(); i++) {
			if (first == cells_table[i]) {
				if (i == cells_table.size() - 1)
					return NULL;
				else
					return cells_table[i + 1];
			}
		}
		// We should reach here.
		abort();
	}
}

int hash_cell::num_pages(char set_flags, char clear_flags)
{
	int num = 0;
	pthread_spin_lock(&_lock);
	for (unsigned int i = 0; i < buf.get_num_pages(); i++) {
		thread_safe_page *p = buf.get_page(i);
		if (p->test_flags(set_flags) && !p->test_flags(clear_flags))
			num++;
	}
	pthread_spin_unlock(&_lock);
	return num;
}

void hash_cell::predict_evicted_pages(int num_pages, char set_flags,
		char clear_flags, std::map<off_t, thread_safe_page *> &pages)
{
	pthread_spin_lock(&_lock);
	policy.predict_evicted_pages(buf, num_pages, set_flags,
			clear_flags, pages);
	for (std::map<off_t, thread_safe_page *>::iterator it = pages.begin();
			it != pages.end(); it++) {
		it->second->inc_ref();
	}
	pthread_spin_unlock(&_lock);
}

void hash_cell::get_pages(int num_pages, char set_flags, char clear_flags,
		std::map<off_t, thread_safe_page *> &pages)
{
	pthread_spin_lock(&_lock);
	for (int i = 0; i < (int) buf.get_num_pages(); i++) {
		thread_safe_page *p = buf.get_page(i);
		if (p->test_flags(set_flags) && !p->test_flags(clear_flags)) {
			p->inc_ref();
			pages.insert(std::pair<off_t, thread_safe_page *>(
						p->get_offset(), p));
		}
	}
	pthread_spin_unlock(&_lock);
}

#define ENABLE_FLUSH_THREAD

void associative_flush_thread::dirty_pages(thread_safe_page *pages[], int num)
{
#ifdef ENABLE_FLUSH_THREAD
	hash_cell *cells[num];
	int num_queued_cells = 0;
	for (int i = 0; i < num; i++) {
		hash_cell *cell = local_cache->get_cell_offset(pages[i]->get_offset());
		if (!cell->is_in_queue()) {
			char dirty_flag = 0;
			char skip_flags = 0;
			page_set_flag(dirty_flag, DIRTY_BIT, true);
			// We should skip pages in IO pending or in a writeback queue.
			page_set_flag(skip_flags, IO_PENDING_BIT, true);
			page_set_flag(skip_flags, PREPARE_WRITEBACK, true);
			/*
			 * We only count the number of dirty pages without IO pending.
			 * If a page is dirty but has IO pending, it means the page
			 * is being written back, so we don't need to do anything with it.
			 */
			int n = cell->num_pages(dirty_flag, skip_flags);
			if (n > DIRTY_PAGES_THRESHOLD && !cell->set_in_queue(true))
				cells[num_queued_cells++] = cell;
		}
	}
	if (num_queued_cells > 0) {
		// TODO currently, there is only one flush thread. Adding dirty cells
		// requires to grab a spin lock. It may not work well on a NUMA machine.
		dirty_cells.add(cells, num_queued_cells);
		activate();
	}
#endif
}
