/**
 * Copyright 2014 Open Connectome Project (http://openconnecto.me)
 * Written by Da Zheng (zhengda1936@gmail.com)
 *
 * This file is part of SAFSlib.
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

#include <boost/format.hpp>

#include "remote_access.h"
#include "parameters.h"
#include "slab_allocator.h"
#include "disk_read_thread.h"
#include "file_mapper.h"

namespace safs
{

static const int COMPLETE_QUEUE_SIZE = 10240;

/**
 * An IO request may be split into multiple requests.
 * This helper class represents the original I/O request issued by users.
 */
class remote_orig_io_request: public io_request
{
	atomic_number<ssize_t> completed_size;
public:
	static remote_orig_io_request *cast2original(io_request *req) {
		return (remote_orig_io_request *) req;
	}

	void init(const io_request &req) {
		assert(req.get_req_type() == io_request::BASIC_REQ);
		io_request::init(req);
		completed_size = atomic_number<ssize_t>();
	}

	bool complete_part(const io_request &part) {
		ssize_t ret = completed_size.inc(part.get_size());
		return ret == this->get_size();
	}

	bool is_complete() const {
		return completed_size.get() >= this->get_size();
	}
};

/*
 * This method is invoked in the I/O thread.
 * It queues the completed I/O requests in the queue and these I/O requests
 * will be processed in the application threads later.
 */
void remote_io::notify_completion(io_request *reqs[], int num)
{
	stack_array<io_request> req_copies(num);
	for (int i = 0; i < num; i++) {
		req_copies[i] = *reqs[i];
		assert(req_copies[i].get_io());
	}

	BOOST_VERIFY(complete_queue.add(req_copies.data(), num) == num);
	get_thread()->activate();
}

remote_io::remote_io(const std::vector<disk_io_thread::ptr> &remotes,
		slab_allocator &_msg_allocator, file_mapper *mapper, thread *t,
		const safs_header &header, int max_reqs): io_interface(t,
			header), max_disk_cached_reqs(max_reqs), complete_queue(std::string(
					"disk_complete_queue-") + itoa(t->get_node_id()), t->get_node_id(),
				COMPLETE_QUEUE_SIZE, std::numeric_limits<int>::max()),
			msg_allocator(_msg_allocator)
{
	int node_id = t->get_node_id();
	num_ios.inc(1);
	this->io_threads = remotes;
	senders.resize(remotes.size());
	low_prio_senders.resize(remotes.size());
	// create a msg sender for each disk read thread.
	for (unsigned i = 0; i < remotes.size(); i++) {
		senders[i] = request_sender::create(node_id, &msg_allocator,
				remotes[i]->get_queue());
		low_prio_senders[i] = request_sender::create(node_id, &msg_allocator,
				remotes[i]->get_low_prio_queue());
	}
	cb = NULL;
	this->block_mapper = mapper;
}

remote_io::~remote_io()
{
	cleanup();
	assert(senders.size() == low_prio_senders.size());
	int num_senders = senders.size();
	for (int i = 0; i < num_senders; i++) {
		request_sender::destroy(senders[i]);
		request_sender::destroy(low_prio_senders[i]);
	}
}

io_interface *remote_io::clone(thread *t) const
{
	// An IO may not be associated to any threads.
	ASSERT_TRUE(t);
	num_ios.inc(1);
	remote_io *copy = new remote_io(io_threads, msg_allocator,
			block_mapper, t, get_header(), this->max_disk_cached_reqs);
	copy->cb = this->cb;
	return copy;
}

void remote_io::cleanup()
{
	process_completed_requests(complete_queue.get_num_entries());

	num_ios.dec(1);
	for (unsigned i = 0; i < senders.size(); i++) {
		senders[i]->flush();
		low_prio_senders[i]->flush();
	}
	int num;
	do {
		num = 0;
		assert(senders.size() == low_prio_senders.size());
		for (unsigned i = 0; i < senders.size(); i++) {
			num += senders[i]->get_queue()->get_num_entries();
			num += low_prio_senders[i]->get_queue()->get_num_entries();
		}
		/* 
		 * if there are still messages in the queue, wait.
		 * this might be the best I can do right now
		 * unless the queues can notify me when they are
		 * empty.
		 */
		if (num > 0) {
			// Let's wake up all IO threads if there are still
			// some low-priority requests that need to be processed.
			for (unsigned i = 0; i < io_threads.size(); i++) {
				io_threads[i]->activate();
			}
			usleep(100000);
		}
	} while (num > 0);

	for (unsigned i = 0; i < io_threads.size(); i++)
		io_threads[i]->flush_requests();
}

void remote_io::access(io_request *requests, int num,
		io_status *status)
{
	ASSERT_EQ(get_thread(), thread::get_curr_thread());
	num_issued_reqs.inc(num);

	bool syncd = false;
	for (int i = 0; i < num; i++) {
		if (requests[i].get_io() == NULL) {
			requests[i].set_io(this);
			requests[i].set_node_id(this->get_node_id());
		}

		if (requests[i].get_access_method() == WRITE && !params.is_writable())
			throw io_exception((boost::format(
							"The I/O object can't write data. offset: %1%, size: %2%")
						% requests[i].get_offset() % requests[i].get_size()).str());
		if (requests[i].get_offset() % MIN_BLOCK_SIZE > 0)
			throw io_exception((boost::format(
						"The IO request offset isn't aligned. offset: %1%, size: %2%")
					% requests[i].get_offset() % requests[i].get_size()).str());
		if (requests[i].get_size() % MIN_BLOCK_SIZE > 0)
			throw io_exception((boost::format(
							"The IO request size isn't aligned. offset: %1%, size: %2%")
						% requests[i].get_offset() % requests[i].get_size()).str());
		if (requests[i].get_req_type() == io_request::USER_COMPUTE)
			throw io_exception("user compute isn't supported");

		if (requests[i].is_flush()) {
			syncd = true;
			num_completed_reqs.inc(1);
			continue;
		}
		else if (requests[i].is_sync()) {
			syncd = true;
		}

		// If the request accesses one RAID block, it's simple.
		if (requests[i].inside_RAID_block(get_block_size())) {
			off_t pg_off = requests[i].get_offset() / PAGE_SIZE;
			int idx = block_mapper->map2file(pg_off);
			// The cache inside a sender is extensible, so it can absorb
			// all requests.
			int ret;
			if (requests[i].is_high_prio())
				ret = senders[idx]->send_cached(&requests[i]);
			else {
				ret = low_prio_senders[idx]->send_cached(&requests[i]);
			}
			assert(ret == 1);
		}
		else {
			// If the request accesses multiple RAID blocks, we have to
			// split the request.
			// I still use the default memory allocator, but since it is used
			// when the request size is large, it should normally be OK.
			// TODO I can use slab allocators later.
			remote_orig_io_request *orig = new remote_orig_io_request();
			// global_cached_io doesn't issue requests across a block boundary.
			// It can only be application issued requst, so it shouldn't have
			// extension.
			assert(!requests[i].is_extended_req());
			orig->init(requests[i]);
			off_t end = orig->get_offset() + orig->get_size();
			const off_t block_size = get_block_size() * PAGE_SIZE;
			for (off_t begin = orig->get_offset(); begin < end;
					begin = ROUND(begin + block_size, block_size)) {
				io_req_extension *ext = new io_req_extension();
				ext->set_priv(orig);
				io_request req(ext, INVALID_DATA_LOC, 0, NULL, 0);
				int size = ROUND(begin + block_size, block_size) - begin;
				size = min(size, end - begin);
				// It only supports to extract a specified request from
				// a single-buffer request.
				orig->extract(begin, size, req);
				req.set_io(this);
				assert(req.inside_RAID_block(get_block_size()));

				// Send a request.
				off_t pg_off = req.get_offset() / PAGE_SIZE;
				int idx = block_mapper->map2file(pg_off);
				// The cache inside a sender is extensible, so it can absorb
				// all requests.
				int ret;
				if (req.is_high_prio())
					ret = senders[idx]->send_cached(&req);
				else {
					ret = low_prio_senders[idx]->send_cached(&req);
				}
				assert(ret == 1);
			}
		}
	}

	int num_remaining = 0;
	assert(senders.size() == low_prio_senders.size());
	for (unsigned i = 0; i < senders.size(); i++) {
		num_remaining += senders[i]->get_num_remaining();
		num_remaining += low_prio_senders[i]->get_num_remaining();
	}
	if (num_remaining > MAX_DISK_CACHED_REQS) {
		flush_requests(MAX_DISK_CACHED_REQS);
	}
	if (syncd)
		flush_requests();

	if (status)
		for (int i = 0; i < num; i++)
			status[i] = IO_PENDING;
}

void remote_io::flush_requests()
{
	flush_requests(0);
}

/*
 * This method is invoked in the application threads.
 * It processes the completed I/O requests returned by the I/O threads.
 */
int remote_io::process_completed_requests(int num)
{
	if (num > 0) {
		stack_array<io_request> reqs(num);
		int ret = complete_queue.fetch(reqs.data(), num);
		process_completed_requests(reqs.data(), ret);
		return ret;
	}
	else
		return 0;
}

int remote_io::process_completed_requests(io_request reqs[], int num)
{
	// There are a few cases for the incoming requests.
	//	the requests issued by the upper layer IO;
	//	the requests split by the current IO;
	//	the requests issued by an application.
	io_request *from_upper[num];
	io_request *from_app[num];
	io_interface *upper_io = NULL;
	int num_from_upper = 0;
	int num_from_app = 0;
	int num_part_reqs = 0;
	std::vector<remote_orig_io_request *> completes;
	for (int i = 0; i < num; i++) {
		assert(reqs[i].get_io());
		// The requests issued by the upper layer IO.
		if (reqs[i].get_io() != this) {
			if (upper_io == NULL)
				upper_io = reqs[i].get_io();
			else
				// They should be from the same upper layer IO.
				assert(upper_io == reqs[i].get_io());
			from_upper[num_from_upper++] = &reqs[i];
			continue;
		}
		if (reqs[i].get_io() == this && !reqs[i].is_extended_req()) {
			from_app[num_from_app++] = &reqs[i];
			continue;
		}

		// Handle large I/O requests that were split in remote I/O.
		remote_orig_io_request *orig = remote_orig_io_request::cast2original(
				(io_request *) reqs[i].get_priv());
		io_request *req = &reqs[i];
		if (orig->complete_part(*req))
			completes.push_back(remote_orig_io_request::cast2original(orig));
		else {
			num_part_reqs++;
		}
		delete req->get_extension();
	}
	if (num_from_upper > 0) {
		assert(upper_io);
		upper_io->notify_completion(from_upper, num_from_upper);
	}
	if (num_from_app > 0 && this->have_callback())
		this->get_callback().invoke(from_app, num_from_app);

	for (unsigned i = 0; i < completes.size(); i++) {
		remote_orig_io_request *orig = completes[i];
		io_request *req = (io_request *) orig;
		io_interface *io = orig->get_io();
		// It's from an application.
		if (io == this) {
			if (io->have_callback())
				io->get_callback().invoke(&req, 1);
		}
		else
			io->notify_completion(&req, 1);
		// Now we can delete it.
		delete orig;
	}

	num_completed_reqs.inc(num - num_part_reqs);
	return num - num_part_reqs;
}

void remote_io::flush_requests(int max_cached)
{
	// Now let's flush requests to the queues, but we first try to
	// flush requests non-blockingly.
	assert(senders.size() == low_prio_senders.size());
	int num_senders = senders.size();
	for (int i = 0; i < num_senders; i++) {
		bool has_data = false;
		if (senders[i]->get_num_remaining() > 0) {
			has_data = true;
			senders[i]->flush();
		}
		if (low_prio_senders[i]->get_num_remaining() > 0) {
			has_data = true;
			low_prio_senders[i]->flush();
		}
		if (has_data)
			io_threads[i]->activate();
	}
}

int remote_io::get_file_id() const
{
	return block_mapper->get_file_id();
}

/**
 * We wait for at least the specified number of requests to complete.
 */
int remote_io::wait4complete(int num_to_complete)
{
	flush_requests();
	int pending = num_pending_ios();
	num_to_complete = min(pending, num_to_complete);

	process_completed_requests(complete_queue.get_num_entries());
	int iters = 0;
	while (pending - num_pending_ios() < num_to_complete) {
		iters++;
		get_thread()->wait();
		process_completed_requests(complete_queue.get_num_entries());
	}
	return pending - num_pending_ios();
}

void remote_io::print_state()
{
	printf("remote_io %d has %d pending reqs, %d completed reqs\n",
			get_io_id(), num_pending_ios(), complete_queue.get_num_entries());
	for (unsigned i = 0; i < senders.size(); i++)
		printf("\tsender %d: remain %d reqs\n", i,
				senders[i]->get_num_remaining());
	for (unsigned i = 0; i < low_prio_senders.size(); i++)
		printf("\tlow-prio sender %d: remain %d reqs\n", i,
				low_prio_senders[i]->get_num_remaining());
}

atomic_integer remote_io::num_ios;

}
