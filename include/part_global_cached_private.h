#ifndef __PART_GLOBAL_CACHED_PRIVATE_H__
#define __PART_GLOBAL_CACHED_PRIVATE_H__

#include <errno.h>

#include <tr1/unordered_map>

#include "parameters.h"
#include "messaging.h"
#include "global_cached_private.h"

const int REPLY_BUF_SIZE = 1000;
const int REQ_BUF_SIZE = 1000;
const int MSG_BUF_SIZE = 128;

class part_global_cached_io;

struct thread_group;
class part_io_process_table;

/**
 * This provides interface for application threads issue IO requests
 * in the NUMA architecture. If the requests are sent to disks connected
 * to the local processor, they are processed directly. No message passing
 * is needed.
 */
class part_global_cached_io: public global_cached_io
{
	const part_io_process_table *global_table;
	const struct thread_group *local_group;
	const cache_config *cache_conf;

	msg_queue<io_reply> *reply_queue;
	// This reply message buffer is used when copying remote messages
	// to the local memory.
	message<io_reply> local_reply_msgs[MSG_BUF_SIZE];
	// This request buffer is used when distributing requests.
	io_request local_req_buf[REQ_BUF_SIZE];

	/*
	 * there is a sender for each node.
	 */
	// group id <-> msg sender
	std::tr1::unordered_map<int, request_sender *> req_senders;

	/*
	 * These reply senders are to send replies to this IO. They are made
	 * to be thread-safe, so all threads can use them. However, the remote
	 * access on a NUMA machine is slow, so each NUMA node has a copy of
	 * the reply sender to improve performance.
	 */
	std::vector<thread_safe_msg_sender<io_reply> *> reply_senders;

	// All these variables are updated in one thread, so it's fine without
	// any concurrency control.
	long processed_requests;
	long sent_requests;
	long processed_replies;
	long remote_reads;

	// It's the callback from the user.
	callback *final_cb;

	int process_replies();
	int distribute_reqs(io_request *requests, int num);

	part_global_cached_io(int node_id, part_io_process_table *);
	~part_global_cached_io();
public:
	static part_io_process_table *open_file(
			std::map<int, io_interface *> &underlyings,
			const cache_config *config, int num_ssds);
	static int close_file(part_io_process_table *table);

	static part_global_cached_io *create(int node_id,
			part_io_process_table *table) {
		assert(node_id >= 0);
		void *addr = numa_alloc_onnode(sizeof(part_global_cached_io), node_id);
		return new(addr) part_global_cached_io(node_id, table);
	}

	static void destroy(part_global_cached_io *io) {
		io->~part_global_cached_io();
		numa_free(io, sizeof(*io));
	}

	int init();

	thread_safe_msg_sender<io_reply> *get_reply_sender(int node_id) const {
		return reply_senders[node_id];
	}

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

	void flush_requests() {
		for (std::tr1::unordered_map<int, request_sender *>::const_iterator it
				= req_senders.begin(); it != req_senders.end(); it++) {
			it->second->flush(true);
		}
		global_cached_io::flush_requests();
	}

	void cleanup();
	int preload(off_t start, long size);

	bool support_aio() {
		return true;
	}

	friend class node_cached_io;

#ifdef STATISTICS
	virtual void print_stat(int nthreads);
#endif
};

#endif
