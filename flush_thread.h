#ifndef __FLUSH_THREAD_H__
#define __FLUSH_THREAD_H__

#include "common.h"
#include "thread.h"

class io_request;
class thread_safe_page;
class flush_thread: public thread
{
public:
	flush_thread(int node_id): thread(std::string("flush_thread-")
			+ itoa(node_id), node_id) {
	}
	virtual void request_callback(io_request &req) = 0;
	virtual void dirty_pages(thread_safe_page *pages[], int num) = 0;
};

#endif
