#ifndef __IO_REQUEST_H__
#define __IO_REQUEST_H__

/*
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

#include <string.h>
#include <sys/uio.h>
#include <sys/time.h>
#include <limits.h>

#include <algorithm>
#include <queue>

#include "common.h"
#include "concurrency.h"
#include "container.h"
#include "parameters.h"

namespace safs
{

static const int PAGE_SIZE = 4096;
static const int LOG_PAGE_SIZE = 12;

#define ROUND_PAGE(off) (((long) off) & (~((long) safs::PAGE_SIZE - 1)))
#define ROUNDUP_PAGE(off) (((long) off + safs::PAGE_SIZE - 1) & (~((long) safs::PAGE_SIZE - 1)))

class thread_safe_page;
class io_interface;

class io_buf
{
	union {
		thread_safe_page *p;
		void *buf;
	} u;
	unsigned int size: 31;
	unsigned int is_page: 1;
public:
	io_buf() {
		u.buf = NULL;
		size = 0;
		is_page = 0;
	}

	void init(void *p, int size, bool is_page) {
		u.buf = p;
		if (is_page)
			assert(size == PAGE_SIZE);
		this->size = size;
		this->is_page = is_page;
	}

	void init(thread_safe_page *p);
	void init(void *buf, int size) {
		u.buf = buf;
		this->size = size;
		is_page = 0;
	}

	void *get_buf() const;

	int get_size() const {
		return size;
	}

	thread_safe_page *get_page() const {
		assert(is_page);
		return u.p;
	}
};

class io_request;

class io_req_extension
{
	void *priv;

	int num_bufs: 16;
	int vec_capacity: 15;

	io_buf *vec_pointer;
	io_buf embedded_vecs[NUM_EMBEDDED_IOVECS];

	struct timeval issue_time;

public:
	io_req_extension() {
		vec_pointer = embedded_vecs;
		vec_capacity = NUM_EMBEDDED_IOVECS;
		init();
	}

	~io_req_extension() {
		if (vec_pointer != embedded_vecs)
			delete [] vec_pointer;
	}

	bool is_valid() const {
		// Valid extension must have vec_pointer points to embedded_vecs
		// or an array.
		return vec_pointer != NULL;
	}

	void init() {
		this->priv = NULL;
		this->num_bufs = 0;
		memset(vec_pointer, 0, vec_capacity * sizeof(io_buf));
		memset(&issue_time, 0, sizeof(issue_time));
	}

	void init(const io_req_extension &ext) {
		this->priv = ext.priv;
		this->num_bufs = ext.num_bufs;
		assert(this->vec_capacity >= ext.vec_capacity);
		memcpy(vec_pointer, ext.vec_pointer, num_bufs * sizeof(*vec_pointer));
		memset(&issue_time, 0, sizeof(issue_time));
	}

	void *get_priv() const {
		return priv;
	}

	void set_priv(void *priv) {
		this->priv = priv;
	}

	void set_timestamp() {
		gettimeofday(&issue_time, NULL);
	}

	struct timeval get_timestamp() {
		return issue_time;
	}

	void add_io_buf(const io_buf &buf);
	void add_buf(char *buf, int size, bool is_page);
	void add_buf_front(char *buf, int size, bool is_page);
	int get_num_bufs() const {
		return num_bufs;
	}

	const io_buf &get_buf(int idx) const {
		return vec_pointer[idx];
	}

	int get_size() const {
		ssize_t size = 0;
		for (int i = 0; i < num_bufs; i++)
			size += vec_pointer[i].get_size();
		return size;
	}
};

const int MAX_INLINE_SIZE=128;

/**
 * The data type for a SAFS file identifier.
 */
typedef int file_id_t;
/**
 * Invalid SAFS file identifier.
 */
const int INVALID_FILE_ID = -1;

/**
 * This class specifies a data location in SAFS.
 * It includes a SAFS file ID and the location in the file.
 */
class data_loc_t
{
	file_id_t file_id;
	off_t off;
public:
	data_loc_t() {
		file_id = -1;
		off = -1;
	}

	/**
	 * The constructor.
	 * \param file_id the SAFS file.
	 * \param off the location in the file.
	 */
	data_loc_t(file_id_t file_id, off_t off) {
		this->file_id = file_id;
		this->off = off;
	}

	/**
	 * This method gets the SAFS file ID.
	 * \return the SAFS file ID.
	 */
	file_id_t get_file_id() const {
		return file_id;
	}

	/**
	 * This method gets the location in the SAFS file.
	 * \return the location.
	 */
	off_t get_offset() const {
		return off;
	}
};

const data_loc_t INVALID_DATA_LOC;

class user_compute;

/**
 * The class defines a compact data structure for containing the info of
 * an I/O request issued by a user task.
 */
class request_range
{
	data_loc_t loc;
	unsigned long size: 63;
	unsigned long access_method: 1;
	user_compute *compute;
public:
	request_range() {
		size = 0;
		access_method = 0;
		compute = NULL;
	}

	/**
	 * The constructor.
	 * \param loc the data location in SAFS.
	 * \param size the data size of the request.
	 * \param access_method indicates whether to read or write.
	 * \param compute the user task associated with the I/O request.
	 */
	request_range(const data_loc_t &loc, size_t size, int access_method,
			user_compute *compute) {
		this->loc = loc;
		this->size = size;
		this->access_method = access_method & 0x1;
		this->compute = compute;
	}

	/**
	 * This method gets the data location in SAFS.
	 * \return the data location.
	 */
	const data_loc_t &get_loc() const {
		return loc;
	};

	/**
	 * This method gets the I/O request size.
	 * \return the request size.
	 */
	size_t get_size() const {
		return size;
	}

	/**
	 * This method indicates whether to read or write.
	 * \return whether to read or write.
	 */
	int get_access_method() const {
		return access_method & 0x1;
	}

	/**
	 * This method gets the user task associated with the I/O request.
	 * \return the user task.
	 */
	user_compute *get_compute() const {
		return compute;
	}

	/**
	 * This method sets the user task associated with the I/O request.
	 * \param compute the user task.
	 */
	void set_compute(user_compute *compute) {
		this->compute = compute;
	}
};

typedef fifo_queue<io_request> user_comp_req_queue;
class page_byte_array;
class compute_allocator;

/**
 * This class defines the interface of a user task assocaited with
 * an I/O request. The user task is executed in the page cache when
 * the I/O request is complete. It is executd as follows:
 * upon the completion of an I/O request, run() is invoked;
 * SAFS invokes has_requests() to check whether the user task has
 * more I/O requests;
 * if the user task has more requests, SAFS invokes get_next_request() 
 * to get more requests;
 * SAFS invokes has_completed() to check whether the user task has completed;
 * if a user task has been completed, SAFS destroys the user task.
 */
class user_compute
{
	compute_allocator *alloc;
	atomic_flags<int> flags;
	int num_refs;
public:
	enum {
		IN_QUEUE,
	};

	/**
	 * The constructor.
	 * \param alloc the object allocator that allocates the user task.
	 */
	user_compute(compute_allocator *alloc) {
		this->alloc = alloc;
	}

	/**
	 * This method gets the object allocator that allocates the user task.
	 * \return the object allocator.
	 */
	compute_allocator *get_allocator() const {
		return alloc;
	}

	virtual ~user_compute() {
	}

	/**
	 * This method serialize the user task to a buffer.
	 * It's currently not used.
	 * \param buf the data buffer where the user task is serialized to.
	 * \param size the buffer size.
	 */
	virtual int serialize(char *buf, int size) const = 0;
	/**
	 * This method gets the serialized size of the user task.
	 * It's currently not used.
	 * \return the serialized size of the user task.
	 */
	virtual int get_serialized_size() const = 0;
	/**
	 * This method executes the user task on the data read by the I/O request
	 * that the user task is associated with. The data read by the I/O
	 * request is stored in the page cache.
	 * \param arr the byte array that contains the data read by the I/O
	 * request and is stored in the page cache.
	 */
	virtual void run(page_byte_array &arr) = 0;

	/**
	 * This method indicates whether the user task has been completed.
	 * \return whether the user task has been completed.
	 */
	virtual bool has_completed() = 0;

	/**
	 * This method indicates whether the user task has more I/O requests to
	 * be issued.
	 * \return whether the user task has more I/O requests to be issued.
	 */
	virtual int has_requests() = 0;

	/**
	 * This method get the next I/O request generated by the user task.
	 * When the method is invoked, we have to call has_requests() to check
	 * that the user task has more I/O requests.
	 * \return the next I/O request.
	 */
	virtual request_range get_next_request() = 0;

	virtual void set_flag(int flag, bool value) {
		if (value)
			flags.set_flag(flag);
		else
			flags.clear_flag(flag);
	}

	virtual bool test_flag(int flag) const {
		return flags.test_flag(flag);
	}

	/**
	 * FIXME set the direction of iterating the I/O requests generated
	 * by the user task.
	 */
	virtual void set_scan_dir(bool forward) {
	}

	/**
	 * This method fetches an I/O request from the user task. This is
	 * a helper method that wraps on the user-defined get_next_request.
	 * \param io the I/O instance associated with the fetched I/O request.
	 * \param req the fetched I/O request.
	 * \return true if a user can fetch an I/O request.
	 */
	bool fetch_request(io_interface *io, io_request &req);

	/**
	 * This method fetches multiple I/O requests from the user task.
	 * This is also a helper method that wraps on the user-defined
	 * get_next_request.
	 * \param io the I/O instance associated with the fetched I/O request.
	 * \param reqs the array where the fetched I/O requests should be stored.
	 * \param max_fetch the maximal number of I/O requests should be fetched
	 * from the user task.
	 * \return the number of I/O requests fetched from the user task.
	 */
	int fetch_requests(io_interface *io, user_comp_req_queue &reqs,
			int max_fetch);

	/**
	 * This method increases the reference count of the object.
	 */
	void inc_ref() {
		num_refs++;
	}

	/**
	 * This method decreases the reference count of this object.
	 */
	void dec_ref() {
		num_refs--;
	}

	/**
	 * This method gets the reference count of this object.
	 * \return the reference count.
	 */
	int get_ref() const {
		return num_refs;
	}
};

/**
 * This class defines the interface of allocating customized user tasks.
 */
class compute_allocator
{
public:
	virtual ~compute_allocator() {
	}
	/**
	 * This method alloates a user task.
	 * \return the allocated user task.
	 */
	virtual user_compute *alloc() = 0;
	/**
	 * This method deallocates a user task.
	 * \param compute the user task to be deallocated.
	 */
	virtual void free(user_compute *compute) = 0;
};

/**
 * This class defines an I/O request from users.
 * There are three forms of I/O reuqests:
 *	simple form,
 *	user-task form,
 *	extended form.
 *
 * In the simple form, an I/O request contains a single data buffer where
 * data read from the disk needs to be stored or data needs to be written
 * to the disk.
 *
 * In the user-task form, an I/O request does not contains a data buffer.
 * Instead, it contains a user task to be executed on the data covered
 * by the I/O request. This form is currently only used for read requests.
 * The user task is executed once the completion of the I/O request.
 *
 * In the extended form, an I/O request can contain multiple data buffers
 * to contain data. This form is currently only used internally and is not
 * allowed to be used by users.
 */
class io_request
{
	static const off_t MAX_FILE_SIZE = LONG_MAX;
	static const int MAX_NODE_ID = (1 << 8) - 1;

	size_t buf_size;
	off_t offset;

	// These two flags decide how the payload is interpreted, so they are
	// initialized when the object is created and can't be changed manually.
	// The only case that payload_type is changed is when
	// another request object is copied to this request object.
	// What needs to be guaranteed is that the pointer in the payload always
	// points to the object of the right type.
	// However, data_inline is never changed. It is usually false. It is
	// used when to interpret the io requests in a message.
	unsigned int payload_type: 2;
	unsigned int data_inline: 1;

	unsigned int access_method: 1;
	// Is this synchronous IO?
	unsigned int sync: 1;
	unsigned int high_prio: 1;
	unsigned int low_latency: 1;
	unsigned int discarded: 1;
	unsigned int node_id: 8;
	int file_id;

	io_interface *io;
	void *user_data;

	union {
		/**
		 * These two can be used by users.
		 */
		void *buf_addr;
		user_compute *compute;

		/**
		 * This field should only be used by the library itself.
		 */
		io_req_extension *ext;

		/**
		 * This field is used to help address the memory following the request.
		 * It is used when a user wants to embed user data in the request.
		 */
		char buf[0];
	} payload;

	void use_default_flags() {
		sync = 0;
		high_prio = 1;
		low_latency = 0;
		discarded = 0;
	}

	void copy_flags(const io_request &req) {
		this->sync = req.sync;
		this->high_prio = req.high_prio;
		this->low_latency = req.low_latency;
	}

	void set_int_buf_size(size_t size) {
		buf_size = size;
	}

	size_t get_int_buf_size() const {
		return buf_size;
	}

public:
	enum {
		BASIC_REQ,
		EXT_REQ,
		USER_COMPUTE,
	};

	static size_t get_max_req_size() {
		return std::numeric_limits<size_t>::max();
	}

	// By default, a request is initialized as a flush request.
	explicit io_request(bool sync = false) {
		memset(this, 0, sizeof(*this));
		payload_type = BASIC_REQ;
		use_default_flags();
		this->sync = sync;
	}

	/**
	 * The constructor of an I/O request in the simple form.
	 * \param buf the data buffer.
	 * \param loc the location in SAFS.
	 * \param size the request size.
	 * \param access_method indicates whether to read or write.
	 * \param io the I/O instance associated with the I/O request.
	 * \param node_id the NUMA node where the I/O request is issued.
	 * \param sync
	 */
	io_request(char *buf, const data_loc_t &loc, ssize_t size,
			int access_method, io_interface *io = NULL,
			int node_id = MAX_NODE_ID, bool sync = false) {
		payload_type = BASIC_REQ;
		data_inline = 0;
		user_data = NULL;
		init(buf, loc, size, access_method, io, node_id);
		use_default_flags();
		this->sync = sync;
	}

	io_request(io_req_extension *ext, const data_loc_t &loc, int access_method,
			io_interface *io = NULL, int node_id = MAX_NODE_ID) {
		payload_type = EXT_REQ;
		data_inline = 0;
		payload.ext = ext;
		user_data = NULL;
		init(NULL, loc, 0, access_method, io, node_id);
		use_default_flags();
		this->sync = false;
	}

	/**
	 * The constructor of an I/O request in the user-task form.
	 * \param compute the user task.
	 * \param loc the location in SAFS.
	 * \param size the request size.
	 * \param access_method indicates whether to read or write.
	 * \param io the I/O instance associated with the I/O request.
	 * \param node_id the NUMA node where the I/O request is issued.
	 */
	io_request(user_compute *compute, const data_loc_t &loc, ssize_t size,
			int access_method, io_interface *io = NULL,
			int node_id = MAX_NODE_ID) {
		payload_type = USER_COMPUTE;
		data_inline = 0;
		user_data = NULL;
		init(NULL, loc, size, access_method, io, node_id);
		payload.compute = compute;
		use_default_flags();
		this->sync = false;
	}

	void init(const io_request &req) {
		data_loc_t loc(req.get_file_id(), req.get_offset());
		assert(!data_inline);
		if (req.payload_type == USER_COMPUTE
				|| this->payload_type == USER_COMPUTE) {
			assert(req.payload_type == USER_COMPUTE
					&& this->payload_type == USER_COMPUTE);
			// TODO
		}
		else if (!req.is_extended_req()) {
			this->init(req.get_buf(), loc, req.get_size(),
					req.get_access_method(), req.get_io(), req.get_node_id());
		}
		// Both requests have extensions.
		else if (this->is_extended_req()) {
			this->init(NULL, loc, 0, req.get_access_method(),
					req.get_io(), req.get_node_id());
			this->get_extension()->init(*req.get_extension());
		}
		else {
			// The last case is that this request doesn't have extension,
			// but the given request has extension. This request can't keep
			// all information in the given request.
			assert(!this->is_extended_req() && req.is_extended_req());
		}
		copy_flags(req);
		this->user_data = req.user_data;
	}

	void init() {
		data_inline = 0;
		if (is_extended_req()) {
			io_req_extension *ext = get_extension();
			assert(ext);
			ext->init();
		}
		else {
			payload_type = BASIC_REQ;
			payload.buf_addr = NULL;
		}
		file_id = 0;
		offset = 0;
		high_prio = 0;
		sync = 0;
		node_id = MAX_NODE_ID;
		io = NULL;
		access_method = 0;
		set_int_buf_size(0);
		user_data = NULL;
	}

	void init(char *buf, const data_loc_t &loc, ssize_t size,
			int access_method, io_interface *io, int node_id);
	void init(const data_loc_t &loc, int access_method, io_interface *io,
			int node_id) {
		init(NULL, loc, 0, access_method, io, node_id);
	}

	int get_req_type() const {
		return payload_type;
	}

	io_req_extension *get_extension() const {
		assert(is_extended_req() && payload.ext);
		if (data_inline)
			return (io_req_extension *) payload.buf;
		else
			return payload.ext;
	}

	file_id_t get_file_id() const;

	/*
	 * Test whether the request is a flush request.
	 * The difference of a sync'd request and a flush request is that
	 * a flush request isn't a valid request for accessing data.
	 */
	bool is_flush() const {
		return sync && high_prio && (payload.buf_addr == NULL);
	}

	bool is_sync() const {
		return sync;
	}

	bool is_extended_req() const {
		return payload_type == EXT_REQ;
	}

	off_t get_offset() const {
		return offset;
	}

	void set_data_loc(const data_loc_t &loc) {
		this->file_id = loc.get_file_id();
		this->offset = loc.get_offset();
	}

	int get_access_method() const {
		return access_method & 0x1;
	}

	void set_io(io_interface *io) {
		this->io = io;
	}

	io_interface *get_io() const {
		return io;
	}

	int get_node_id() const {
		return node_id;
	}

	void set_node_id(int node_id) {
		assert(node_id <= MAX_NODE_ID);
		this->node_id = node_id;
	}

	bool is_discarded() const {
		return (discarded & 0x1) == 1;
	}

	void set_discarded(bool discarded) {
		this->discarded = discarded;
	}

	bool is_high_prio() const {
		return (high_prio & 0x1) == 1;
	}

	void set_high_prio(bool high_prio) {
		this->high_prio = high_prio;
	}

	bool is_low_latency() const {
		return (low_latency & 0x1) == 1;
	}

	void set_low_latency(bool low_latency) {
		this->low_latency = low_latency;
	}

	/*
	 * The requested data is inside a page on the disk.
	 */
	bool within_1page() const {
		return get_offset() + get_size() <= ROUND_PAGE(get_offset()) + PAGE_SIZE;
	}

	bool contain_offset(off_t off) const {
		return get_offset() <= off && off < get_offset() + get_size();
	}

	/*
	 * Test if the request has overlap with the specified range.
	 */
	bool has_overlap(off_t off, ssize_t size) const {
		// the beginning of the range inside the input request
		return (off >= this->get_offset() && off < this->get_offset()
				+ this->get_size())
			// or the end of the range inside the input request
			|| (off + size >= this->get_offset()
					&& off + size < this->get_offset() + this->get_size())
			// or the input request is inside the range.
			|| (off <= this->get_offset()
					&& off + size >= this->get_offset() + this->get_size());
	}

	int get_overlap_size(thread_safe_page *pg) const;

	int get_num_covered_pages() const {
		off_t begin_pg = ROUND_PAGE(get_offset());
		off_t end_pg = ROUNDUP_PAGE(get_offset() + get_size());
		return (end_pg - begin_pg) / PAGE_SIZE;
	}

	bool inside_RAID_block(int block_size) const {
		int block_size_bytes = block_size * PAGE_SIZE;
		return ROUND(this->get_offset(), block_size_bytes)
			== ROUND(this->get_offset() + this->get_size() - 1, block_size_bytes);
	}

	void *get_user_data() const {
		return user_data;
	}

	void set_user_data(void *data) {
		this->user_data = data;
	}

	void *get_priv() const {
		return get_extension()->get_priv();
	}

	void set_priv(void *priv) {
		get_extension()->set_priv(priv);
	}

	bool is_empty() const {
		return get_extension()->get_num_bufs() == 0;
	}

	bool is_valid() const {
		return get_offset() != -1;
	}

	ssize_t get_size() const {
		if (!is_extended_req()) {
			return get_int_buf_size();
		}
		else {
			return get_extension()->get_size();
		}
	}

	user_compute *get_compute() const {
		assert(get_req_type() == io_request::USER_COMPUTE);
		return payload.compute;
	}

	/*
	 * By default, we get the first buffer. This makes sense
	 * for a single buffer request.
	 */
	char *get_buf(int idx = 0) const {
		if (!is_extended_req()) {
			return (char *) payload.buf_addr;
		}
		else {
			return (char *) get_extension()->get_buf(idx).get_buf();
		}
	}

	thread_safe_page *get_page(int idx) const {
		return get_extension()->get_buf(idx).get_page();
	}

	void add_buf(char *buf, int size) {
		get_extension()->add_buf(buf, size, 0);
	}

	void add_page(thread_safe_page *p) {
		get_extension()->add_buf((char *) p, PAGE_SIZE, 1);
	}

	void add_io_buf(const io_buf &buf) {
		get_extension()->add_io_buf(buf);
	}

	void add_buf_front(char *buf, int size) {
		get_extension()->add_buf_front(buf, size, 0);
	}

	void add_page_front(thread_safe_page *p) {
		get_extension()->add_buf_front((char *) p, PAGE_SIZE, 1);
	}

	int get_num_bufs() const {
		if (is_extended_req())
			return get_extension()->get_num_bufs();
		else
			return 1;
	}

	int get_buf_size(int idx) const {
		if (!is_extended_req()) {
			assert(idx == 0);
			return get_int_buf_size();
		}
		else
			return get_extension()->get_buf(idx).get_size();
	}

	const io_buf &get_io_buf(int idx) const {
		return get_extension()->get_buf(idx);
	}

	const struct iovec get(int idx) const {
		struct iovec ret;
		io_req_extension *ext = get_extension();
		ret.iov_base = ext->get_buf(idx).get_buf();
		ret.iov_len = ext->get_buf(idx).get_size();
		return ret;
	}

	const int get_vec(struct iovec *vec, int num) const {
		num = min(get_num_bufs(), num);
		for (int i = 0; i < num; i++) {
			io_req_extension *ext = get_extension();
			vec[i].iov_base = ext->get_buf(i).get_buf();
			vec[i].iov_len = ext->get_buf(i).get_size();
		}
		return num;
	}

	bool is_data_inline() const {
		return data_inline == 1;
	}

	void set_timestamp() {
		get_extension()->set_timestamp();
	}

	struct timeval get_timestamp() {
		return get_extension()->get_timestamp();
	}

	/*
	 * Extract a request from the input request.
	 * The extract request is within the range [off, off + size).
	 */
	void extract(off_t off, int size, io_request &extracted) const {
		off_t req_off;
		char *req_buf;
		ssize_t req_size;
		assert(get_num_bufs() == 1);
		// We have to make sure the extracted range has overlap with
		// the input request.
		bool check = has_overlap(off, size);
		if (!check)
			fprintf(stderr, "req %lx, size: %lx, page off: %lx\n",
					this->get_offset(), this->get_size(), off);
		assert(check);
		// this is the first page in the request.
		if (off <= this->get_offset()) {
			req_off = this->get_offset();
			req_buf = this->get_buf();
		}
		else {
			req_off = off;
			/* 
			 * We can't be sure if the request buffer is aligned
			 * with the page size.
			 */
			req_buf = this->get_buf() + (off - this->get_offset());
		}
		req_size = min(off + size - req_off,
				this->get_offset() + this->get_size() - req_off);
		data_loc_t loc(this->get_file_id(), req_off);
		extracted.init(req_buf, loc, req_size, this->get_access_method(),
				this->get_io(), this->get_node_id());
	}

	/*
	 * We need to serialize an io request to a buffer so it can be sent to
	 * another thread.
	 * @accept_inline: indicates whether the IO request can inline data
	 * to the buffer.
	 */
	int serialize(char *buf, int size, bool accept_inline) {
		int serialized_size;
		if (is_data_inline()) {
			assert(accept_inline);
			serialized_size = get_serialized_size();
			assert(serialized_size <= size);
			memcpy(buf, this, serialized_size);
		}
		else if (payload_type == EXT_REQ) {
			// We never serialize the io request extension to the message.
			serialized_size = sizeof(io_request);
			assert(serialized_size <= size);
			memcpy(buf, this, sizeof(*this));
		}
		else if (payload_type == BASIC_REQ) {
			// We only serialize the data buffer to the message for write
			// requests. The size of the data buffer has to be small.
			if (get_size() > MAX_INLINE_SIZE || access_method == READ
					|| !accept_inline) {
				serialized_size = sizeof(io_request);
				assert(serialized_size <= size);
				memcpy(buf, this, sizeof(*this));
			}
			else {
				serialized_size = sizeof(io_request) - sizeof(this->payload)
					+ get_size();
				assert(serialized_size <= size);
				memcpy(buf, this, sizeof(*this));
				io_request *p = (io_request *) buf;
				p->data_inline = 1;
				memcpy(p->payload.buf, (char *) this->payload.buf_addr,
						this->get_size());
			}
		}
		else {
			// The user compute object is always serialized to the message.
			serialized_size = sizeof(io_request);
			assert(serialized_size <= size && sizeof(io_request)
					<= (unsigned) size);
			memcpy(buf, this, serialized_size);
		}
		return serialized_size;
	}

	/*
	 * This method returns the size of an IO request after it is serialized.
	 */
	int get_serialized_size() const {
		if (payload_type == EXT_REQ)
			return sizeof(io_request);
		else if (payload_type == BASIC_REQ && (get_size() > MAX_INLINE_SIZE
					|| access_method == READ))
			return sizeof(io_request);
		else if (payload_type == BASIC_REQ)
			return sizeof(io_request) - sizeof(this->payload) + get_size();
		else {
			return sizeof(io_request);
		}
	}

	/*
	 * This method deserialize an request from the buffer.
	 * If the request data is inline in the buffer, instead of allocating
	 * memory for the extra objects of the IO request, the extra objects
	 * will stay in the buffer and the created request will point to the buffer.
	 */
	static void deserialize(io_request &req, char *buf, int size) {
		assert((unsigned) size >= sizeof(io_request));
		io_request *p = (io_request *) buf;
		memcpy(&req, buf, sizeof(req));
		if (req.is_data_inline()) {
			assert(req.payload_type != EXT_REQ);
			switch(req.payload_type) {
				case BASIC_REQ:
					req.payload.buf_addr = p->payload.buf;
					break;
				case USER_COMPUTE:
					req.payload.compute = (user_compute *) p->payload.buf;
					break;
				default:
					assert(0);
			}
			req.data_inline = 0;
		}
	}

	static io_request *deserialize(char *buf, int size) {
		assert((unsigned) size >= sizeof(io_request));
		io_request *ret = (io_request *) buf;
		assert(ret->get_serialized_size() <= size);
		return ret;
	}
};

template<class req_type, class get_io_func>
class comp_req_io
{
	get_io_func get_io;

public:
	comp_req_io(get_io_func func) {
		this->get_io = func;
	}

	bool operator() (const req_type req1, const req_type req2) {
		return (long) get_io(req1) < (long) get_io(req2);
	}
};

typedef void (*req_process_func_t)(io_interface *io, io_request *reqs[], int num);
/*
 * Perform the same function to the requests with the same IO instance.
 * It turns out it's a common function when delivering requests to
 * the upper layer.
 */
template<class req_type, class get_io_func, class process_req_func>
void process_reqs_on_io(req_type reqs[], int num,
		get_io_func func, process_req_func proc_func)
{
	comp_req_io<req_type, get_io_func> req_io_comparator(func);

	std::sort(reqs, reqs + num, req_io_comparator);
	io_interface *prev = func(reqs[0]);
	int begin_idx = 0;
	for (int end_idx = 1; end_idx < num; end_idx++) {
		if (func(reqs[end_idx]) != prev) {
			proc_func(prev, reqs + begin_idx, end_idx - begin_idx);
			begin_idx = end_idx;
			prev = func(reqs[end_idx]);
		}
	}
	assert(begin_idx < num);
	proc_func(prev, reqs + begin_idx, num - begin_idx);
}

static inline void process_reqs_on_io(io_request *reqs[],
		int num, req_process_func_t func)
{
	class get_io_func {
	public:
		io_interface *operator()(io_request *req) {
			return req->get_io();
		}
	} io_func;

	class process_req_func {
		req_process_func_t func;
	public:
		process_req_func(req_process_func_t func) {
			this->func = func;
		}

		void operator()(io_interface *io, io_request *reqs[], int num) {
			func(io, reqs, num);
		}
	} proc_func(func);

	process_reqs_on_io(reqs, num, io_func, proc_func);
}

}

#endif
