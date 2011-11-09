#include <iostream>

#include <assert.h>
#include <stdio.h>

#include "cache.h"
#include "tree_cache.h"

int main()
{
	page empty_page;
	printf("offset: %ld, data: %p\n", empty_page.get_offset(),
			empty_page.get_data());
	page new_page(1000, &empty_page);
	printf("offset: %ld, data: %p\n", new_page.get_offset(),
			new_page.get_data());
	page copied_page;
	printf("offset: %ld, data: %p\n", copied_page.get_offset(),
			copied_page.get_data());
	copied_page = new_page;
	printf("offset: %ld, data: %p\n", copied_page.get_offset(),
			copied_page.get_data());

	page_buffer<page> *buffer = new page_buffer<page> (10);
	for (int i = 0; i < 15; i++) {
		bool is_full = buffer->is_full();
		struct page new_page(i, NULL);
		struct page *ret = buffer->get_empty_page();
		/* the number of pages that can be contained in the buffer is
		 * 9, so we should see 9 true, and then 6 fales */
		std::cout<<"buffer full? "<<is_full
			<<", empty page: "<<ret->get_offset()<<std::endl;
		*ret = new_page;
	}
	page_cache *cache = new tree_cache(10);
}
