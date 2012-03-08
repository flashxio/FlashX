#ifndef __ASSOCIATIVE_CACHE_H__
#define __ASSOCIATIVE_CACHE_H__

#include <vector>

#include "cache.h"

#define CELL_SIZE 8

#ifdef STATISTICS
static volatile int avail_cells;
static volatile int num_wait_unused;
static volatile int lock_contentions;
#endif

#define USE_SHADOW_PAGE

const int CACHE_LINE = 128;

/**
 * This data structure is to implement LRU.
 */
template<class T, int BUF_SIZE>
class page_cell
{
	unsigned int idx;		// to the point where we can evict a page in the buffer
	T buf[BUF_SIZE];			// a circular buffer to keep pages.

public:
	/*
	 * @size: the size of the page buffer
	 * @page_buf: the offset of the page array in the global page cache.
	 */
	page_cell() {
		idx = 0;
	}

	void set_pages(long page_buf) {
		for (int i = 0; i < BUF_SIZE; i++) {
			buf[i] = T(-1, page_buf + i * PAGE_SIZE);
		}
		idx = 0;
	}

	/**
	 * return an empty page.
	 * I expected the page will be filled with data,
	 * so I change the begin and end index of the circular buffer.
	 */
	T *get_empty_page() {
		/* TODO I ignore the case of integer overflow */
		T *ret = &buf[idx % BUF_SIZE];
		idx++;
		return ret;
	}

	T *get_page(int i) {
		if (i >= BUF_SIZE)
			return NULL;
		return &buf[i];
	}

	int get_idx(T *page) {
		int idx = page - buf;
		assert (idx >= 0 && idx < BUF_SIZE);
		return idx;
	}

	void scale_down_hits() {
		for (int i = 0; i < BUF_SIZE; i++) {
			buf[i].set_hits(buf[i].get_hits() / 2);
		}
	}
};

class shadow_page
{
	int offset;
	unsigned char hits;
public:
	shadow_page() {
		offset = -1;
		hits = 0;
	}
	shadow_page(page &pg) {
		offset = pg.get_offset() >> LOG_PAGE_SIZE;
		hits = pg.get_hits();
	}

	off_t get_offset() const {
		return ((off_t) offset) << LOG_PAGE_SIZE;
	}

	int get_hits() {
		return hits;
	}

	void set_hits(int hits) {
		assert (hits <= 0xff);
		this->hits = hits;
	}

	bool is_valid() {
		return offset != -1;
	}
};

int middle_evicts = 0;
int end_evicts = 0;

template<class T, int SIZE>
class generic_queue
{
	unsigned short start;
	unsigned short num;
	/* the size of the buffer is specified by SIZE. */
	T buf[0];
public:
	generic_queue() {
		assert(SIZE < 0xffff);
		start = 0;
		num = 0;
	}

	void push_back(T v) {
		assert(num < SIZE);
		buf[(start + num) % SIZE] = v;
		num++;
	}

	void pop_front() {
		end_evicts++;
		assert(num > 0);
		start = (start + 1) % SIZE;
		num--;
	}

	/*
	 * remove the idx'th element in the queue.
	 * idx is the logical position in the queue,
	 * instead of the physical index in the buffer.
	 */
	void remove(int idx) {
		assert(idx < num);
		/* the first element in the queue. */
		if (idx == 0) {
			end_evicts++;
			pop_front();
		}
		/* the last element in the queue. */
		else if (idx == num - 1){
			end_evicts++;
			num--;
		}
		/*
		 * in the middle.
		 * now we need to move data.
		 */
		else {
			middle_evicts++;
			T tmp[num];
			T *p = tmp;
			/* if the end of the queue is physically behind the start */
			if (start + num <= SIZE) {
				/* copy elements in front of the removed element. */
				memcpy(p, &buf[start], sizeof(T) * idx);
				p += idx;
				/* copy elements behind the removed element. */
				memcpy(p, &buf[start + idx + 1], sizeof(T) * (num - idx - 1));
			}
			/* 
			 * the removed element is between the first element
			 * and the end of the buffer.
			 */
			else if (idx + start < SIZE) {
				/* copy elements in front of the removed element. */
				memcpy(p, &buf[start], sizeof(T) * idx);
				p += idx;
				/*
				 * copy elements behind the removed element
				 * and before the end of the buffer.
				 */
				memcpy(p, &buf[start + idx + 1], sizeof(T) * (SIZE - start - idx - 1));
				p += (SIZE - start - idx - 1);
				/* copy the remaining elements in the beginning of the buffer. */
				memcpy(p, buf, sizeof(T) * (num - (SIZE - start)));
			}
			/*
			 * the removed element is between the beginning of the buffer
			 * and the last element.
			 */
			else {
				/* copy elements between the first element and the end of the buffer. */
				memcpy(p, &buf[start], sizeof(T) * (SIZE - start));
				p += (SIZE - start);
				/* copy elements between the beginning of the buffer and the removed element. */
				idx = (idx + start) % SIZE;
				memcpy(p, buf, sizeof(T) * idx);
				p += idx;
				/* copy elements after the removed element and before the last element */
				memcpy(p, &buf[idx + 1], sizeof(T) * ((start + num) % SIZE - idx - 1));
			}
			memcpy(buf, tmp, sizeof(T) * (num - 1));
			start = 0;
			num--;
		}
	}

	bool is_empty() {
		return num == 0;
	}

	bool is_full() {
		return num == SIZE;
	}

	int size() {
		return num;
	}

	T &back() {
		assert(num > 0);
		return buf[(start + num - 1) % SIZE];
	}

	T &front() {
		assert(num > 0);
		return buf[start];
	}

	T &get(int idx) {
		assert(num > 0);
		return buf[(start + idx) % SIZE];
	}

	void print_state() {
		printf("start: %d, num: %d\n", start, num);
		for (int i = 0; i < this->size(); i++)
			printf("%ld\t", this->get(i));
		printf("\n");
	}
};

template<int SHADOW_SIZE>
class shadow_cell
{
	generic_queue<shadow_page, SHADOW_SIZE> queue;
public:
	shadow_cell() {
	}

	/*
	 * add a page to the cell, which is evicted from hash_cell.
	 * the only thing we need to record is the number of hits
	 * of this page.
	 */
	void add(shadow_page pg) {
		if (queue.is_full())
			queue.pop_front();
		queue.push_back(pg);
	}

	shadow_page search(off_t off) {
		for (int i = 0; i < queue.size(); i++) {
			shadow_page pg = queue.get(i);
			if (pg.get_offset() == off) {
				queue.remove(i);
				queue.push_back(pg);
				return pg;
			}
		}
		return shadow_page();
	}

	void scale_down_hits() {
		for (int i = 0; i < queue.size(); i++) {
			queue.get(i).set_hits(queue.get(i).get_hits() / 2);
		}
	}
};

// TODO the entire cell should be put in the same cache line
// so each access to the hash table has only one cache miss.
class hash_cell
{
	pthread_spinlock_t _lock;
	page_cell<thread_safe_page, CELL_SIZE> buf;
#ifdef USE_LRU
	std::vector<int> pos_vec;
#endif
#ifdef USE_SHADOW_PAGE
#define NUM_SHADOW_PAGES 37
	shadow_cell<NUM_SHADOW_PAGES> shadow;
#define STUFFING_SIZE (CACHE_LINE * 4 - sizeof(_lock) - sizeof(buf) - sizeof(shadow) - 8)
	char stuffing[STUFFING_SIZE];
#endif

	/* this function has to be called with lock held */
	thread_safe_page *get_empty_page() {
		thread_safe_page *ret = NULL;

		do {
			int min_hits = 0x7fffffff;
			for (int i = 0; i < CELL_SIZE; i++) {
				thread_safe_page *pg = buf.get_page(i);
				if (pg->get_ref())
					continue;

				/* 
				 * refcnt only increases within the lock of the cell,
				 * so if the page's refcnt is 0 above,
				 * it'll be always 0 within the lock.
				 */

				if (min_hits > pg->get_hits()) {
					min_hits = pg->get_hits();
					ret = pg;
				}
			}
			/* it happens when all pages in the cell is used currently. */
		} while (ret == NULL);

		/* we record the hit info of the page in the shadow cell. */
#ifdef USE_SHADOW_PAGE
		if (ret->get_hits() > 0)
			shadow.add(shadow_page(*ret));
#endif

		ret->reset_hits();
		ret->set_data_ready(false);
		return ret;
	}

#ifdef USE_LRU
	/* 
	 * the end of the vector points to the pages
	 * that are most recently accessed.
	 */
	thread_safe_page *get_empty_page() {
		int pos;
		if (pos_vec.size() < CELL_SIZE) {
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
#endif

#ifdef USE_FIFO
    /* this function has to be called with lock held */
	thread_safe_page *get_empty_page() {
		thread_safe_page *ret = buf.get_empty_page();
		// TODO I assume this situation is rare
		while (ret->get_ref()) {
			ret = buf.get_empty_page();
			printf("try another empty page.\n");
		}
		ret->set_data_ready(false);
		return ret;
	}
#endif

public:
	hash_cell() {
		pthread_spin_init(&_lock, PTHREAD_PROCESS_PRIVATE);
	}

	~hash_cell() {
		pthread_spin_destroy(&_lock);
	}

	void *operator new[](size_t size) {
		printf("allocate %ld bytes\n", size);
		void *addr = memalign(CACHE_LINE, size + CACHE_LINE);
		// TODO 8 might be architecture specific. It's 8 for 64-bit machines.
		return (void *) ((long) addr + CACHE_LINE - 8);
	}

	void operator delete[](void *p) {
		free((void *) ((long) p - (CACHE_LINE - 8)));
	}

	void set_pages(long page_buf) {
		buf.set_pages(page_buf);
	}

	/**
	 * search for a page with the offset.
	 * If the page doesn't exist, return an empty page.
	 */
	page *search(off_t off, off_t &old_off) {
		thread_safe_page *ret = NULL;
#ifndef STATISTICS
		pthread_spin_lock(&_lock);
#else
		if (pthread_spin_trylock(&_lock) == EBUSY) {
			__sync_fetch_and_add(&lock_contentions, 1);
			pthread_spin_lock(&_lock);
		}
#endif

		for (int i = 0; i < CELL_SIZE; i++) {
			if (buf.get_page(i)->get_offset() == off) {
				ret = buf.get_page(i);
				break;
			}
		}
		if (ret == NULL) {
			ret = get_empty_page();
			old_off = ret->get_offset();
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
#ifdef USE_LRU
		else {
			/* move the page to the end of the pos vector. */
			int pos = buf.get_idx(ret);
			for (std::vector<int>::iterator it = pos_vec.begin();
					it != pos_vec.end(); it++) {
				if (*it == pos) {
					pos_vec.erase(it);
					break;
				}
			}
			pos_vec.push_back(pos);
		}
#endif
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

	void print_cell() {
		for (int i = 0; i < CELL_SIZE; i++)
			printf("%lx\t", buf.get_page(i)->get_offset());
		printf("\n");
	}
};

class associative_cache: public page_cache
{
	hash_cell *cells;
	int ncells;
	// TODO maybe it's not a good hash function
	int hash(off_t offset) {
		return offset / PAGE_SIZE % ncells;
	}

public:
	associative_cache(long cache_size) {
		printf("associative cache is used\n");
		int npages = cache_size / PAGE_SIZE;
		ncells = npages / CELL_SIZE;
		cells = new hash_cell[ncells];
		printf("%d cells: %p\n", ncells, cells);
		for (int i = 0; i < ncells; i++)
			cells[i].set_pages(i * PAGE_SIZE * CELL_SIZE);
	}

	~associative_cache() {
		delete [] cells;
	}

	page *search(off_t offset, off_t &old_off) {
		hash_cell *cell = &cells[hash(offset)];
		return cell->search(offset, old_off);
	}

	void print_cell(off_t off) {
		cells[hash(off)].print_cell();
	}
};

#endif
