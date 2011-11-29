#include <iostream>

#include <assert.h>
#include <stdio.h>

#include "cache.h"
#include "tree_cache.h"
#include "cuckoo_cache.h"
#include "associative_cache.h"

int main()
{
	printf("cell size: %ld\n", sizeof(hash_cell));
#if 0
	page pg(0, 0x51 * 4096);
	printf("data: %p\n", pg.get_data());
	pg.set_data_ready(true);
	printf("data: %p, ready: %d\n", pg.get_data(), pg.data_ready());

	thread_safe_page tspg(0, 0x52 * 4096);
	printf("thread_safe_page data: %p, ready: %d, ref: %d\n",
			tspg.get_data(), tspg.data_ready(), tspg.get_ref());
	tspg.set_data_ready(true);
	printf("thread_safe_page data: %p, ready: %d, ref: %d\n",
			tspg.get_data(), tspg.data_ready(), tspg.get_ref());
	tspg.inc_ref();
	printf("thread_safe_page data: %p, ready: %d, ref: %d\n",
			tspg.get_data(), tspg.data_ready(), tspg.get_ref());
	tspg.set_data_ready(false);
	printf("thread_safe_page data: %p, ready: %d, ref: %d\n",
			tspg.get_data(), tspg.data_ready(), tspg.get_ref());
	tspg.dec_ref();
	printf("thread_safe_page data: %p, ready: %d, ref: %d\n",
			tspg.get_data(), tspg.data_ready(), tspg.get_ref());
	printf("thread_safe_page data: %p, ready: %d, ref: %d, io pending: %d\n",
			tspg.get_data(), tspg.data_ready(), tspg.get_ref(), tspg.test_and_set_io_pending());
	printf("thread_safe_page data: %p, ready: %d, ref: %d, io pending: %d\n",
			tspg.get_data(), tspg.data_ready(), tspg.get_ref(), tspg.test_and_set_io_pending());
#endif

#if 0
	// test lockable_pointer
	page pg;
	printf("orig: %p\n", &pg);
	lockable_pointer<page> pointer(&pg);
	pointer.lock();
	printf("lock: %p\n", pointer.operator->());
	pointer.unlock();
	printf("unlock: %p\n", pointer.operator->());
#endif

#if 0
	// test cuckoo cache
	cuckoo_hash hash(16);

	for (int i = 0; i < 16; i++) {
		hash.insert(i, new page(i, NULL));
	}
	for (int i = 0; i < 16; i++) {
		page *pg = hash.search(i);
		if (pg == NULL)
			printf("page %d doesn't exist in the hash table\n", i);
		else if (pg->get_offset() == i)
			printf("find page %d\n", i);
		else
			printf("wrong page %ld in the hash table\n", pg->get_offset());
	}
#endif

#if 0
	// test page
	page::allocate_cache(4096 * 10);
	page *buf = new page[10];
	for (int i = 0; i < 10; i++) {
		buf[i] = page(0xfffffff0L + i * PAGE_SIZE, i * PAGE_SIZE);
	}
	for (int i = 0; i < 10; i++) {
		printf("page %d: offset: %lx, data: %p, data_ready? %d\n",
				i, buf[i].get_offset(), buf[i].get_data(), buf[i].data_ready());
		buf[i].set_data_ready(true);
		printf("page %d: data ready? %d\n", i, buf[i].data_ready());
	}

	thread_safe_page *buf1 = new thread_safe_page[10];
	for (int i = 0; i < 10; i++) {
		buf1[i] = thread_safe_page(0xfffffff0L + i * PAGE_SIZE, i * PAGE_SIZE);
	}
	for (int i = 0; i < 10; i++) {
		printf("page %d: offset: %lx, data: %p, data_ready? %d, io_pending? %d\n",
				i, buf1[i].get_offset(), buf1[i].get_data(),
				buf1[i].data_ready(), buf1[i].test_and_set_io_pending());
		buf1[i].set_data_ready(true);
		printf("page %d: data ready? %d, io_pending? %d\n",
				i, buf1[i].data_ready(), buf1[i].test_and_set_io_pending());
	}

	page empty_page;
	printf("offset: %lx, data: %p\n", empty_page.get_offset(),
			empty_page.get_data());
	page new_page(1000, 0);
	printf("offset: %lx, data: %p\n", new_page.get_offset(),
			new_page.get_data());
	page copied_page;
	printf("offset: %lx, data: %p\n", copied_page.get_offset(),
			copied_page.get_data());
	copied_page = new_page;
	printf("offset: %lx, data: %p\n", copied_page.get_offset(),
			copied_page.get_data());

	// test page_buffer
	page_buffer<page> *buffer = new page_buffer<page> (10, 0);
	for (int i = 0; i < 15; i++) {
		bool is_full = buffer->is_full();
		struct page new_page(i * PAGE_SIZE, 0);
		struct page *ret = buffer->get_empty_page();
		/* the number of pages that can be contained in the buffer is
		 * 9, so we should see 9 true, and then 6 fales */
		printf("buffer full? %d, empty page: %lx, initialized? %d\n",
				is_full, ret->get_offset(), ret->initialized());
		*ret = new_page;
	}
	page_cache *cache = new tree_cache(10, 0);
#endif
}
