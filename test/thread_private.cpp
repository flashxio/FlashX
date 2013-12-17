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

#include <sys/time.h>

#include "thread_private.h"
#include "parameters.h"
#include "exception.h"

bool align_req = false;
int align_size = PAGE_SIZE;

extern struct timeval global_start;

static atomic_number<long> global_sum;

struct issued_workload_t
{
	workload_t work;
	struct timeval issue_time;
};

class cleanup_callback: public callback
{
	rand_buf *buf;
	ssize_t read_bytes;
	int thread_id;
	thread_private *thread;
	int file_id;
public:
#ifdef DEBUG
	std::tr1::unordered_map<char *, issued_workload_t> pending_reqs;
#endif

	cleanup_callback(rand_buf *buf, int idx, thread_private *thread, int file_id) {
		this->buf = buf;
		read_bytes = 0;
		this->thread_id = idx;
		this->thread = thread;
		this->file_id = file_id;
	}

	int invoke(io_request *rqs[], int num) {
		assert(thread == thread::get_curr_thread());
		for (int i = 0; i < num; i++) {
			io_request *rq = rqs[i];
			if (rq->get_access_method() == READ && params.is_verify_content()) {
				off_t off = rq->get_offset();
				for (int i = 0; i < rq->get_num_bufs(); i++) {
					assert(check_read_content(rq->get_buf(i),
								rq->get_buf_size(i), off, file_id));
					off += rq->get_buf_size(i);
				}
			}
#ifdef DEBUG
			assert(rq->get_num_bufs() == 1);
			pending_reqs.erase(rq->get_buf(0));
#endif
			for (int i = 0; i < rq->get_num_bufs(); i++)
				buf->free_entry(rq->get_buf(i));
			read_bytes += rq->get_size();
		}
#ifdef STATISTICS
		thread->num_completes.inc(num);
		int res = thread->num_pending.dec(num);
		assert(res >= 0);
#endif
		return 0;
	}

	ssize_t get_size() {
		return read_bytes;
	}
};

class sum_user_compute: public user_compute
{
public:
	sum_user_compute(compute_allocator *alloc): user_compute(alloc) {
	}

	virtual int serialize(char *buf, int size) const {
		assert(0);
		return 0;
	}

	virtual int get_serialized_size() const {
		assert(0);
		return 0;
	}

	virtual bool run(page_byte_array &array) {
		long sum = 0;
		const page_byte_array &const_array = page_byte_array::const_cast_ref(
				array);
		page_byte_array::const_iterator<long> end = const_array.end<long>();
		for (page_byte_array::const_iterator<long> it = const_array.begin<long>();
				it != end; ++it) {
			sum += *it;
		}

		global_sum.inc(sum);
		return true;
	}
};

class write_user_compute: public user_compute
{
	int file_id;
public:
	write_user_compute(compute_allocator *alloc): user_compute(
			alloc) {
		this->file_id = -1;
	}

	void init(int file_id) {
		this->file_id = file_id;
	}

	virtual int serialize(char *buf, int size) const {
		assert(0);
		return 0;
	}

	virtual int get_serialized_size() const {
		assert(0);
		return 0;
	}

	virtual bool run(page_byte_array &array) {
		page_byte_array::iterator<long> end = array.end<long>();
		off_t off = array.get_offset();
		for (page_byte_array::iterator<long> it = array.begin<long>();
				it != end; ++it) {
			*it = off / sizeof(off_t) + file_id;
			off += sizeof(off_t);
		}
		return true;
	}
};

class write_compute_allocator: public compute_allocator
{
	class compute_initiator: public obj_initiator<write_user_compute>
	{
		write_compute_allocator *alloc;
	public:
		compute_initiator(write_compute_allocator *alloc) {
			this->alloc = alloc;
		}

		virtual void init(write_user_compute *obj) {
			new (obj) write_user_compute(alloc);
		}
	};

	obj_allocator<write_user_compute> allocator;
public:
	write_compute_allocator(thread *t): allocator("write-compute-allocator",
			t->get_node_id(), PAGE_SIZE, INT_MAX,
			// TODO memory leak
			new compute_initiator(this)) {
	}

	virtual user_compute *alloc() {
		return allocator.alloc_obj();
	}

	virtual void free(user_compute *obj) {
		allocator.free((write_user_compute *) obj);
	}
};

class sum_compute_allocator: public compute_allocator
{
	class compute_initiator: public obj_initiator<sum_user_compute>
	{
		sum_compute_allocator *alloc;
	public:
		compute_initiator(sum_compute_allocator *alloc) {
			this->alloc = alloc;
		}

		virtual void init(sum_user_compute *obj) {
			new (obj) sum_user_compute(alloc);
		}
	};

	obj_allocator<sum_user_compute> allocator;
public:
	sum_compute_allocator(thread *t): allocator("sum-compute-allocator",
			t->get_node_id(), PAGE_SIZE, INT_MAX,
			// TODO memory leak
			new compute_initiator(this)) {
	}

	virtual user_compute *alloc() {
		return allocator.alloc_obj();
	}

	virtual void free(user_compute *obj) {
		allocator.free((sum_user_compute *) obj);
	}
};

ssize_t thread_private::get_read_bytes() {
	if (cb)
		return cb->get_size();
	else
		return read_bytes;
}

thread_private::thread_private(int node_id, int idx, int entry_size,
		file_io_factory *factory, workload_gen *gen): thread(
			std::string("test_thread") + itoa(idx), node_id)
{
	this->cb = NULL;
	this->node_id = node_id;
	this->idx = idx;
	buf = NULL;
	this->gen = gen;
	this->io = NULL;
	this->factory = factory;
	read_bytes = 0;
	num_accesses = 0;
	num_sampling = 0;
	tot_num_pending = 0;
	max_num_pending = 0;
	sum_alloc = new sum_compute_allocator(this);
	write_alloc = new write_compute_allocator(this);
}

void thread_private::init() {
	io = factory->create_io(this);
	io->set_max_num_pending_ios(params.get_aio_depth_per_file());
	io->init();

	rand_buf *buf = new rand_buf(config.get_buf_size(),
			config.get_entry_size(), node_id);
	this->buf = buf;
	if (io->support_aio()) {
		cb = new cleanup_callback(buf, idx, this, io->get_file_id());
		io->set_callback(cb);
	}
}

/**
 * This converts a workload to I/O access ranges.
 */
class work2req_range_converter
{
	int file_id;
	workload_t workload;
public:
	work2req_range_converter() {
		file_id = -1;
		workload.off = -1;
		workload.size = -1;
		workload.read = 0;
	}

	work2req_range_converter(int file_id, const workload_t &workload, int align_size) {
		this->file_id = file_id;
		this->workload = workload;
		if (align_size > 0) {
			this->workload.off = ROUND(workload.off, align_size);
			this->workload.size = ROUNDUP(workload.off
					+ workload.size, align_size)
				- ROUND(workload.off, align_size);
		}
	}

	bool has_complete() const {
		return workload.size <= 0;
	}

	request_range get_request(int buf_type) {
		int access_method = workload.read ? READ : WRITE;
		if (buf_type == MULTI_BUF) {
			throw unsupported_exception();
#if 0
			assert(off % PAGE_SIZE == 0);
			int num_vecs = size / PAGE_SIZE;
			reqs[i].init(off, io, access_method, node_id);
			assert(buf->get_entry_size() >= PAGE_SIZE);
			for (int k = 0; k < num_vecs; k++) {
				reqs[i].add_buf(buf->next_entry(PAGE_SIZE), PAGE_SIZE);
			}
			workload.off += size;
			workload.size = 0;
#endif
		}
		else if (buf_type == SINGLE_SMALL_BUF) {
			off_t off = workload.off;
			int size = workload.size;
			// Get to the next page.
			off_t next_off = ROUNDUP_PAGE(off + 1);
			if (next_off > off + size)
				next_off = off + size;
			size -= next_off - off;
			workload.off = next_off;
			workload.size = size;

			data_loc_t loc(file_id, off);
			return request_range(loc, next_off - off, access_method);
		}
		else {
			int size = workload.size;
			off_t off = workload.off;
			workload.off += size;
			workload.size = 0;

			data_loc_t loc(file_id, off);
			return request_range(loc, size, access_method);
		}
	}
};

/**
 * This converts workloads to IO requests.
 */
class work2req_converter
{
	work2req_range_converter workload;
	int align_size;
	rand_buf *buf;
	io_interface *io;
	sum_compute_allocator *sum_alloc;
	write_compute_allocator *write_alloc;
public:
	work2req_converter(io_interface *io, rand_buf * buf, int align_size,
			sum_compute_allocator *sum_alloc, write_compute_allocator *write_alloc) {
		this->align_size = align_size;
		this->buf = buf;
		this->io = io;
		this->sum_alloc = sum_alloc;
		this->write_alloc = write_alloc;
	}

	void init(const workload_t &workload) {
		this->workload = work2req_range_converter(io->get_file_id(),
				workload, align_size);
	}

	bool has_complete() const {
		return workload.has_complete();
	}

	int to_reqs(io_interface *io, int buf_type, int num, io_request reqs[]);
};

int work2req_converter::to_reqs(io_interface *io, int buf_type, int num,
		io_request reqs[])
{
	int node_id = io->get_node_id();


	int i = 0;
	while (!workload.has_complete() && i < num) {
		request_range range = workload.get_request(buf_type);

		if (config.is_user_compute()) {
			if (range.get_access_method() == READ) {
				sum_user_compute *compute = (sum_user_compute *) sum_alloc->alloc();
				reqs[i] = io_request(compute, range.get_loc(), range.get_size(),
						range.get_access_method(), io, node_id);
			}
			else {
				write_user_compute *compute = (write_user_compute *) write_alloc->alloc();
				compute->init(io->get_file_id());
				reqs[i] = io_request(compute, range.get_loc(), range.get_size(),
						range.get_access_method(), io, node_id);
			}
		}
		else {
			data_loc_t loc = range.get_loc();
			char *p = buf->next_entry(range.get_size());
			assert(p);
			if (range.get_access_method() == WRITE && params.is_verify_content())
				create_write_data(p, range.get_size(), loc.get_offset(),
						io->get_file_id());
			reqs[i].init(p, loc, range.get_size(), range.get_access_method(),
					io, node_id);
		}
		i++;
	}

	return i;
}

void thread_private::run()
{
	gettimeofday(&start_time, NULL);
	io_request reqs[NUM_REQS_BY_USER];
	char *entry = NULL;
	if (config.is_use_aio())
		assert(io->support_aio());
	if (!config.is_use_aio()) {
		entry = (char *) valloc(config.get_entry_size());
	}
	work2req_converter converter(io, buf, align_size, sum_alloc, write_alloc);
	while (gen->has_next()) {
		if (config.is_use_aio()) {
			int i;
			bool no_mem = false;
			int num_reqs_by_user = min(io->get_remaining_io_slots(), NUM_REQS_BY_USER);
			for (i = 0; i < num_reqs_by_user; ) {
				if (converter.has_complete() && gen->has_next()) {
					converter.init(gen->next());
				}
				if (converter.has_complete())
					break;
				int ret = converter.to_reqs(io, config.get_buf_type(),
						num_reqs_by_user - i, reqs + i);
				if (ret == 0) {
					no_mem = true;
					break;
				}
				i += ret;
			}
#ifdef STATISTICS
			/*
			 * cached IO may complete a request before calling wait4complte().
			 * so we need to increase num_pending before access() is called.
			 */
			int curr = num_pending.inc(i);
#endif
			if (i > 0) {
#ifdef DEBUG
				struct timeval curr;

				gettimeofday(&curr, NULL);
				for (int k = 0; k < i; k++) {
					workload_t work = {reqs[k].get_offset(),
						(int) reqs[k].get_size(), reqs[k].get_access_method() == READ};
					issued_workload_t issue_work;
					issue_work.work = work;
					issue_work.issue_time = curr;
					cb->pending_reqs.insert(std::pair<char *, issued_workload_t>(
								reqs[k].get_buf(0), issue_work));
				}
#endif
				io->access(reqs, i);
			}
#ifdef STATISTICS
			if (max_num_pending < curr)
				max_num_pending = curr;
			if (num_accesses % 100 == 0) {
				num_sampling++;
				tot_num_pending += curr;
			}
#endif
			// We wait if we don't have IO slots left or we can't issue
			// more requests due to the lack of memory.
			if (io->get_remaining_io_slots() <= 0 || no_mem) {
				int num_ios = io->get_max_num_pending_ios() / 10;
				if (num_ios == 0)
					num_ios = 1;
				io->wait4complete(num_ios);
			}
			num_accesses += i;
		}
		else {
			int ret = 0;
			workload_t workload = gen->next();
			off_t off = workload.off;
			int access_method = workload.read ? READ : WRITE;
			int entry_size = workload.size;
			if (align_req) {
				off = ROUND(off, align_size);
				entry_size = ROUNDUP(off + entry_size, align_size)
					- ROUND(off, align_size);
			}

			if (config.get_buf_type() == SINGLE_SMALL_BUF) {
				while (entry_size > 0) {
					/*
					 * generate the data for writing the file,
					 * so the data in the file isn't changed.
					 */
					if (access_method == WRITE && params.is_verify_content()) {
						create_write_data(entry, entry_size, off, io->get_file_id());
					}
					// There is at least one byte we need to access in the page.
					// By adding 1 and rounding up the offset, we'll get the next page
					// behind the current offset.
					off_t next_off = ROUNDUP_PAGE(off + 1);
					if (next_off > off + entry_size)
						next_off = off + entry_size;
					io_status status = io->access(entry, off, next_off - off,
							access_method);
					assert(!(status == IO_UNSUPPORTED));
					if (status == IO_OK) {
						num_accesses++;
						if (access_method == READ && params.is_verify_content()) {
							check_read_content(entry, next_off - off, off, io->get_file_id());
						}
						read_bytes += ret;
					}
					if (status == IO_FAIL) {
						perror("access");
						::exit(1);
					}
					entry_size -= next_off - off;
					off = next_off;
				}
			}
			else {
				if (access_method == WRITE && params.is_verify_content()) {
					create_write_data(entry, entry_size, off, io->get_file_id());
				}
				io_status status = io->access(entry, off, entry_size,
						access_method);
				assert(!(status == IO_UNSUPPORTED));
				if (status == IO_OK) {
					num_accesses++;
					if (access_method == READ && params.is_verify_content()) {
						check_read_content(entry, entry_size, off, io->get_file_id());
					}
					read_bytes += ret;
				}
				if (status == IO_FAIL) {
					perror("access");
					::exit(1);
				}
			}
		}
	}
	struct timeval curr;
	gettimeofday(&curr, NULL);
	printf("thread %d has issued all requests at %ld\n", idx,
			time_diff_us(start_time, curr));
	io->cleanup();
	gettimeofday(&end_time, NULL);

	// Stop itself.
	stop();
}

int thread_private::attach2cpu()
{
#if 0
	cpu_set_t cpuset;
	pthread_t thread = pthread_self();
	CPU_ZERO(&cpuset);
	int cpu_num = idx % NCPUS;
	CPU_SET(cpu_num, &cpuset);
	int ret = pthread_setaffinity_np(thread, sizeof(cpu_set_t), &cpuset);
	if (ret != 0) {
		perror("pthread_setaffinity_np");
		exit(1);
	}
	return ret;
#endif
	return 0;
}

#ifdef USE_PROCESS
static int process_create(pid_t *pid, void (*func)(void *), void *priv)
{
	pid_t id = fork();

	if (id < 0)
		return -1;

	if (id == 0) {	// child
		func(priv);
		exit(0);
	}

	if (id > 0)
		*pid = id;
	return 0;
}

static int process_join(pid_t pid)
{
	int status;
	pid_t ret = waitpid(pid, &status, 0);
	return ret < 0 ? ret : 0;
}
#endif

void thread_private::print_stat()
{
#ifdef STATISTICS
#ifdef DEBUG
	assert(num_pending.get() == (int) cb->pending_reqs.size());
	if (cb->pending_reqs.size() > 0) {
		for (std::tr1::unordered_map<char *, issued_workload_t>::const_iterator it
				= cb->pending_reqs.begin(); it != cb->pending_reqs.end(); it++) {
			workload_t work = it->second.work;
			printf("missing req %lx, size %d, read: %d, issued at %ld\n",
					work.off, work.size, work.read, time_diff_us(
						start_time, it->second.issue_time));
		}
	}
#endif
	io->print_stat(config.get_nthreads());
	int avg_num_pending = 0;
	if (num_sampling > 0)
		avg_num_pending = tot_num_pending / num_sampling;
	printf("access %ld bytes in %ld accesses (%d completes), avg pending: %d, max pending: %d, remaining pending: %d\n",
			get_read_bytes(), num_accesses, num_completes.get(),
			avg_num_pending, max_num_pending, num_pending.get());
#endif
	extern struct timeval global_start;
	printf("thread %d: start at %f seconds, takes %f seconds, access %ld bytes in %ld accesses\n", idx,
			time_diff(global_start, start_time), time_diff(start_time, end_time),
			get_read_bytes(), num_accesses);
}
