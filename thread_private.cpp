#include <sys/time.h>

#include "thread_private.h"
#include "parameters.h"
#include "exception.h"

bool align_req = false;
int align_size = PAGE_SIZE;

extern struct timeval global_start;

class cleanup_callback: public callback
{
	rand_buf *buf;
	ssize_t read_bytes;
	int thread_id;
	thread_private *thread;
public:
#ifdef DEBUG
	std::tr1::unordered_map<char *, workload_t> pending_reqs;
#endif

	cleanup_callback(rand_buf *buf, int idx, thread_private *thread) {
		this->buf = buf;
		read_bytes = 0;
		this->thread_id = idx;
		this->thread = thread;
	}

	int invoke(io_request *rqs[], int num) {
		assert(thread == thread::get_curr_thread());
		for (int i = 0; i < num; i++) {
			io_request *rq = rqs[i];
			if (rq->get_access_method() == READ && params.is_verify_content()) {
				off_t off = rq->get_offset();
				for (int i = 0; i < rq->get_num_bufs(); i++) {
					assert(check_read_content(rq->get_buf(i),
								rq->get_buf_size(i), off));
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

ssize_t thread_private::get_read_bytes() {
	if (cb)
		return cb->get_size();
	else
		return read_bytes;
}

void thread_private::init() {
	io = factory->create_io(this);
	io->set_max_num_pending_ios(params.get_aio_depth_per_file());
	io->init();

	rand_buf *buf = new rand_buf(config.get_buf_size(),
			config.get_entry_size(), node_id);
	this->buf = buf;
	if (io->support_aio()) {
		cb = new cleanup_callback(buf, idx, this);
		io->set_callback(cb);
	}
}

class work2req_converter
{
	workload_t workload;
	int align_size;
	rand_buf *buf;
	io_interface *io;
public:
	work2req_converter(io_interface *io, rand_buf * buf, int align_size) {
		workload.off = -1;
		workload.size = -1;
		workload.read = 0;
		this->align_size = align_size;
		this->buf = buf;
		this->io = io;
	}

	void init(const workload_t &workload) {
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

	int to_reqs(int buf_type, int num, io_request reqs[]);
};

int work2req_converter::to_reqs(int buf_type, int num, io_request reqs[])
{
	int node_id = io->get_node_id();
	int access_method = workload.read ? READ : WRITE;
	off_t off = workload.off;
	int size = workload.size;

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
		int i = 0;
		while (size > 0 && i < num) {
			off_t next_off = ROUNDUP_PAGE(off + 1);
			if (next_off > off + size)
				next_off = off + size;
			char *p = buf->next_entry(next_off - off);
			if (p == NULL)
				break;
			if (access_method == WRITE && params.is_verify_content())
				create_write_data(p, next_off - off, off);
			reqs[i].init(p, off, next_off - off, access_method,
					io, node_id);
			size -= next_off - off;
			off = next_off;
			i++;
		}
		workload.off = off;
		workload.size = size;
		return i;
	}
	else {
		char *p = buf->next_entry(size);
		if (p == NULL)
			return 0;
		if (access_method == WRITE && params.is_verify_content())
			create_write_data(p, size, off);
		reqs[0].init(p, off, size, access_method, io, node_id);
		workload.off += size;
		workload.size = 0;
		return 1;
	}
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
	work2req_converter converter(io, buf, align_size);
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
				int ret = converter.to_reqs(config.get_buf_type(),
						num_reqs_by_user - i, reqs + i);
				if (ret == 0) {
					no_mem = true;
					break;
				}
				i += ret;
			}
			if (i > 0) {
#ifdef DEBUG
				for (int k = 0; k < i; k++) {
					workload_t work = {reqs[k].get_offset(),
						(int) reqs[k].get_size(), reqs[k].get_access_method() == READ};
					cb->pending_reqs.insert(std::pair<char *, workload_t>(
								reqs[k].get_buf(0), work));
				}
#endif
				io->access(reqs, i);
			}
#ifdef STATISTICS
			int curr = num_pending.inc(i);
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
						create_write_data(entry, entry_size, off);
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
							check_read_content(entry, next_off - off, off);
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
					create_write_data(entry, entry_size, off);
				}
				io_status status = io->access(entry, off, entry_size,
						access_method);
				assert(!(status == IO_UNSUPPORTED));
				if (status == IO_OK) {
					num_accesses++;
					if (access_method == READ && params.is_verify_content()) {
						check_read_content(entry, entry_size, off);
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
	printf("thread %d has issued all requests\n", idx);
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
		for (std::tr1::unordered_map<char *, workload_t>::const_iterator it
				= cb->pending_reqs.begin(); it != cb->pending_reqs.end(); it++) {
			workload_t work = it->second;
			printf("missing req %lx, size %d, read: %d\n", work.off, work.size, work.read);
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
