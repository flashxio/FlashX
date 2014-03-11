#ifndef __PART_GLOBAL_CACHED_PRIVATE_H__
#define __PART_GLOBAL_CACHED_PRIVATE_H__

/**
 * Copyright 2013 Da Zheng
 *
 * This file is part of SAFSlib.
 *
 * SAFSlib is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * SAFSlib is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with SAFSlib.  If not, see <http://www.gnu.org/licenses/>.
 */

#include <errno.h>

#include <tr1/unordered_map>

#include "parameters.h"
#include "messaging.h"
#include "global_cached_private.h"

const int REPLY_BUF_SIZE = 1000;
const int REQ_BUF_SIZE = 1000;
const int MSG_BUF_SIZE = 128;

// The size of a reply >= sizeof(io_reply).
const int NUMA_REPLY_BUF_SIZE = NUMA_MSG_SIZE / sizeof(io_reply);

struct thread_group;
class part_io_process_table;
class disk_io_thread;
class file_mapper;
class group_request_sender;
class NUMA_cache;

/**
 * This provides interface for application threads issue IO requests
 * in the NUMA architecture. If the requests are sent to disks connected
 * to the local processor, they are processed directly. No message passing
 * is needed.
 */
class part_global_cached_io: public io_interface
{
	part_io_process_table *global_table;
	const struct thread_group *local_group;
	const cache_config *cache_conf;

	global_cached_io *underlying;

	msg_queue<io_reply> *reply_queue;
	// This reply message buffer is used when copying remote messages
	// to the local memory.
	message<io_reply> local_reply_msgs[MSG_BUF_SIZE];
	// This request buffer is used when distributing requests.
	io_request local_req_buf[REQ_BUF_SIZE];
	// This reply buffer is used when processing replies.
	io_reply local_reply_buf[NUMA_REPLY_BUF_SIZE];

	/*
	 * there is a sender for each node.
	 */
	// group id <-> msg sender
	std::tr1::unordered_map<int, group_request_sender *> req_senders;

	/*
	 * These reply senders are to send replies to this IO. They are made
	 * to be thread-safe, so all threads can use them. However, the remote
	 * access on a NUMA machine is slow, so each NUMA node has a copy of
	 * the reply sender to improve performance.
	 */
	thread_safe_msg_sender<io_reply> *reply_sender;

	// All these variables are updated in one thread, so it's fine without
	// any concurrency control.
	long processed_requests;
	long sent_requests;
	atomic_number<long> processed_replies;
	long remote_reads;

	// It's the callback from the user.
	callback *final_cb;

	int process_replies();
	void notify_upper(io_request *reqs[], int num);

	part_global_cached_io(io_interface *underlying, part_io_process_table *);
	~part_global_cached_io();
public:
	static part_io_process_table *init_subsystem(
			const std::vector<disk_io_thread *> &io_threads,
			file_mapper *mapper, NUMA_cache *cache);
	static int destroy_subsystem(part_io_process_table *table);

	static part_global_cached_io *create(io_interface *underlying,
			part_io_process_table *table) {
		int node_id = underlying->get_node_id();
		assert(node_id >= 0);
		return new part_global_cached_io(underlying, table);
	}

	static void destroy(part_global_cached_io *io) {
		delete io;
	}

	int init();

	virtual bool set_callback(callback *cb) {
		this->final_cb = cb;
		return true;
	}

	virtual callback *get_callback() {
		return final_cb;
	}

	int reply(io_request *requests, io_reply *replies, int num);

	int process_requests(int max_nreqs);

	int process_reply(io_reply *reply);

	virtual void notify_completion(io_request *reqs[], int num);
	void access(io_request *requests, int num, io_status *status);
	io_status access(char *, off_t, ssize_t, int) {
		return IO_UNSUPPORTED;
	}

	void flush_requests();

	void cleanup();
	int preload(off_t start, long size);

	bool support_aio() {
		return true;
	}

	virtual int get_file_id() const {
		return underlying->get_file_id();
	}
	virtual int wait4complete(int num);
	virtual int num_pending_ios() const {
		// the number of pending requests on the remote nodes.
		return sent_requests - processed_replies.get()
			// The number of pending requests in the local IO instance.
			+ underlying->num_pending_ios();
	}

	friend class node_cached_io;

#ifdef STATISTICS
	virtual void print_stat(int nthreads);
#endif
	virtual void print_state();
};

#endif
