#include "cache.h"
#include "messaging.h"

/**
 * Add an IO request to the page when the page is in IO pending state.
 *
 * The request can be added only when data in the page isn't ready.
 * If data is ready, return NULL to notify the invoker.
 */
io_request *thread_safe_page::add_io_req(const io_request &req)
{
	// TODO I need to my own allocator for this.
	io_request *orig = new io_request(req);
	orig->set_priv((void *) this);
	lock();
	bool data_ready = this->data_ready();
	// If data is ready, reqs may be NULL.
//	assert(data_ready || reqs);
	// If data isn't ready, we should add the request to the page.
	if (!data_ready) {
		orig->set_next_req(reqs);
		reqs = orig;
	}
	unlock();
	if (data_ready) {
		delete orig;
		return NULL;
	}
	else
		return orig;
}

/**
 * Add the first IO request to a page.
 * There are three possibilities:
 * 1. this is really the first IO request added to the page.
 * 2. the page already has an IO request, but its data isn't ready yet.
 * 3. data in the page has been ready for reading.
 * It returns NULL for the last case.
 */
io_request *thread_safe_page::add_first_io_req(const io_request &req,
		int &status)
{
	// TODO I need to my own allocator for this.
	io_request *orig = new io_request(req);
	orig->set_priv((void *) this);
	lock();
	bool data_ready = this->data_ready();
	// case 1.
	if (reqs == NULL && !data_ready) {
		reqs = orig;
		status = ADD_REQ_SUCCESS;
	}
	else if (!data_ready) {
		orig->set_next_req(reqs);
		reqs = orig;
		status = ADD_REQ_NOT_DATA_READY;
	}
	else
		status = ADD_REQ_DATA_READY;
	unlock();
	if (status != ADD_REQ_DATA_READY)
		return orig;
	else {
		delete orig;
		return NULL;
	}
}

void thread_safe_page::reset_reqs()
{
	lock();
	io_request *orig = reqs;
	reqs = NULL;
	unlock();

	// We have to make sure no other threads are referencing io requests
	while (orig) {
		io_request *next = orig->get_next_req();
		// TODO use my own deallocator.
		delete orig;
		orig = next;
	}
}
