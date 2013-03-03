#ifndef __FLUSH_THREAD_H__
#define __FLUSH_THREAD_H__

#include "thread.h"

class io_request;
class thread_safe_page;
class flush_thread: public thread
{
public:
	flush_thread(int node_id): thread(node_id) {
	}
	virtual void request_callback(io_request &req) = 0;
	virtual void dirty_pages(thread_safe_page *pages[], int num) = 0;
};

#endif
