#include <sys/time.h>

#include "thread_private.h"
#include "parameters.h"
#include "exception.h"

#define NUM_PAGES (4096 * nthreads)

bool align_req = false;
int align_size = PAGE_SIZE;
extern bool use_aio;

extern bool verify_read_content;
extern struct timeval global_start;
extern int nthreads;

void check_read_content(char *buf, int size, off_t off)
{
	// I assume the space in the buffer is larger than 8 bytes.
	off_t aligned_off = off & (~(sizeof(off_t) - 1));
	long data[2];
	data[0] = aligned_off / sizeof(off_t);
	data[1] = aligned_off / sizeof(off_t) + 1;
	long expected = 0;
	int copy_size = size < (int) sizeof(off_t) ? size : (int) sizeof(off_t);
	memcpy(&expected, ((char *) data) + (off - aligned_off), copy_size);
	long read_value = 0;
	memcpy(&read_value, buf, copy_size);
	if(read_value != expected)
		printf("%ld %ld\n", read_value, expected);
	assert(read_value == expected);
}

void create_write_data(char *buf, int size, off_t off)
{
	off_t aligned_start = off & (~(sizeof(off_t) - 1));
	off_t aligned_end = (off + size) & (~(sizeof(off_t) - 1));
	long start_data = aligned_start / sizeof(off_t);
	long end_data = aligned_end / sizeof(off_t);

	/* If all data is in one 8-byte word. */
	if (aligned_start == aligned_end) {
		memcpy(buf, ((char *) &start_data) + (off - aligned_start), size);
		return;
	}

	int first_size =  (int)(sizeof(off_t) - (off - aligned_start));
	int last_size = (int) (off + size - aligned_end);

	if (first_size == sizeof(off_t))
		first_size = 0;
	if (first_size)
		memcpy(buf, ((char *) &start_data) + (off - aligned_start),
				first_size);
	// Make each buffer written to SSDs different, so it's hard for SSDs
	// to do some tricks on it.
	struct timeval zero = {0, 0};
	long diff = time_diff_us(zero, global_start);
	for (int i = first_size; i < aligned_end - off; i += sizeof(off_t)) {
		*((long *) (buf + i)) = (off + i) / sizeof(off_t) + diff;
	}
	if (aligned_end > aligned_start
			|| (aligned_end == aligned_start && first_size == 0)) {
		if (last_size)
			memcpy(buf + (aligned_end - off), (char *) &end_data, last_size);
	}

	check_read_content(buf, size, off);
}

class cleanup_callback: public callback
{
	rand_buf *buf;
	ssize_t read_bytes;
	int thread_id;
	thread_private *thread;
public:
	cleanup_callback(rand_buf *buf, int idx, thread_private *thread) {
		this->buf = buf;
		read_bytes = 0;
		this->thread_id = idx;
		this->thread = thread;
	}

	int invoke(io_request *rqs[], int num) {
		for (int i = 0; i < num; i++) {
			io_request *rq = rqs[i];
			extern bool verify_read_content;
			if (rq->get_access_method() == READ && verify_read_content) {
				off_t off = rq->get_offset();
				for (int i = 0; i < rq->get_num_bufs(); i++) {
					check_read_content(rq->get_buf(i), rq->get_buf_size(i), off);
					off += rq->get_buf_size(i);
				}
			}
			for (int i = 0; i < rq->get_num_bufs(); i++)
				buf->free_entry(rq->get_buf(i));
			read_bytes += rq->get_size();
		}
#ifdef STATISTICS
		thread->num_completes.inc(num);
		thread->num_pending.dec(num);
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
	extern int io_depth_per_file;
	io = factory->create_io(this);
	io->set_max_num_pending_ios(io_depth_per_file);
	io->init();

	extern int buf_size;
	rand_buf *buf = new rand_buf(NUM_PAGES / (nthreads
				// TODO maybe I should set the right entry size for a buffer.
				// If each access size is irregular, I'll break each access
				// into pages so each access is no larger than a page, so it
				// should workl fine.
				/ NUM_NODES) * PAGE_SIZE, buf_size, node_id);
	this->buf = buf;
	if (io->support_aio()) {
		cb = new cleanup_callback(buf, idx, this);
		io->set_callback(cb);
	}
}

void thread_private::run()
{
	int node_id = io->get_node_id();
	gettimeofday(&start_time, NULL);
	io_request reqs[NUM_REQS_BY_USER];
	char *entry = NULL;
	if (!io->support_aio())
		use_aio = false;
	if (!use_aio) {
		extern int buf_size;
		entry = (char *) valloc(buf_size);
	}
	while (gen->has_next()) {
		if (use_aio) {
			int i;
			int num_reqs_by_user = min(io->get_remaining_io_slots(), NUM_REQS_BY_USER);
			for (i = 0; i < num_reqs_by_user && gen->has_next(); ) {
				workload_t workload = gen->next();
				int access_method = workload.read ? READ : WRITE;
				off_t off = workload.off;
				int size = workload.size;
				if (align_req) {
					off = ROUND(off, align_size);
					size = ROUNDUP(off + size, align_size)
						- ROUND(off, align_size);
				}
				/*
				 * If the size of the request is larger than a page size,
				 * and the user explicitly wants to use multibuf requests.
				 */
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
					i++;
#endif
				}
				else if (buf_type == SINGLE_SMALL_BUF) {
again:
					num_reqs_by_user = min(io->get_remaining_io_slots(), NUM_REQS_BY_USER);
					while (size > 0 && i < num_reqs_by_user) {
						off_t next_off = ROUNDUP_PAGE(off + 1);
						if (next_off > off + size)
							next_off = off + size;
						char *p = buf->next_entry(next_off - off);
						if (p == NULL)
							break;
						if (access_method == WRITE && verify_read_content)
							create_write_data(p, next_off - off, off);
						reqs[i].init(p, off, next_off - off, access_method,
								io, node_id);
						size -= next_off - off;
						off = next_off;
						i++;
					}
					if (size > 0) {
						io->access(reqs, i);
						if (io->get_remaining_io_slots() <= 0) {
							int num_ios = io->get_max_num_pending_ios() / 10;
							if (num_ios == 0)
								num_ios = 1;
							io->wait4complete(num_ios);
						}
#ifdef STATISTICS
						num_pending.inc(i);
#endif
						num_accesses += i;
						i = 0;
						goto again;
					}
				}
				else {
					char *p = buf->next_entry(size);
					if (p == NULL)
						break;
					if (access_method == WRITE && verify_read_content)
						create_write_data(p, size, off);
					reqs[i++].init(p, off, size, access_method, io, node_id);
				}
			}
			io->access(reqs, i);
			if (io->get_remaining_io_slots() <= 0) {
				int num_ios = io->get_max_num_pending_ios() / 10;
				if (num_ios == 0)
					num_ios = 1;
				io->wait4complete(num_ios);
			}
			num_accesses += i;
#ifdef STATISTICS
			int curr = num_pending.inc(i);
			if (max_num_pending < curr)
				max_num_pending = curr;
			if (num_accesses % 100 == 0) {
				num_sampling++;
				tot_num_pending += curr;
			}
#endif
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

			if (buf_type == SINGLE_SMALL_BUF) {
				while (entry_size > 0) {
					/*
					 * generate the data for writing the file,
					 * so the data in the file isn't changed.
					 */
					if (access_method == WRITE && verify_read_content) {
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
						if (access_method == READ && verify_read_content) {
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
				if (access_method == WRITE && verify_read_content) {
					create_write_data(entry, entry_size, off);
				}
				io_status status = io->access(entry, off, entry_size,
						access_method);
				assert(!(status == IO_UNSUPPORTED));
				if (status == IO_OK) {
					num_accesses++;
					if (access_method == READ && verify_read_content) {
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
#if NCPUS > 0
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
#else
	return -1;
#endif
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
