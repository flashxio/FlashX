#ifndef __FLUSH_THREAD_H__
#define __FLUSH_THREAD_H__

#include "thread.h"

class io_request;
class thread_safe_page;
class flush_thread: public thread
{
public:
	virtual void request_callback(io_request &req) = 0;
	virtual void dirty_pages(thread_safe_page *pages[], int num) = 0;
};

#endif
