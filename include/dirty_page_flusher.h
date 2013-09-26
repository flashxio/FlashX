#ifndef __DIRTY_PAGE_FLUSHER_H__
#define __DIRTY_PAGE_FLUSHER_H__

#include "common.h"

class thread_safe_page;
class io_interface;

class dirty_page_flusher
{
public:
	/**
	 * This method flushes dirty pages with the specified io instance.
	 * The method has to be thread-safe, as it will be invoked in different
	 * thread contexts. It also shouldn't block threads.
	 */
	virtual void flush_dirty_pages(thread_safe_page *pages[], int num,
			io_interface *io) = 0;

	virtual int flush_dirty_pages(page_filter *filter, int max_num) = 0;
};

#endif
