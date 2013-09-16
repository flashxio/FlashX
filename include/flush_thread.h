#ifndef __FLUSH_THREAD_H__
#define __FLUSH_THREAD_H__

#include "common.h"
#include "thread.h"

class io_request;
class thread_safe_page;
class page_filter;
class flush_thread: public thread
{
public:
	flush_thread(int node_id): thread(std::string("flush_thread-")
			+ itoa(node_id), node_id) {
	}

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
