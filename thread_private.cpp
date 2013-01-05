#include <sys/time.h>
#include <sys/types.h>
#include <sys/syscall.h>
#define gettid() syscall(__NR_gettid)

#include "thread_private.h"
#include "parameters.h"

#define NUM_PAGES (40960 * nthreads)

extern bool verify_read_content;

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

	int first_size =  (int)(sizeof(off_t) - (off - aligned_start));
	if (first_size == sizeof(off_t))
		first_size = 0;
	if (first_size)
		memcpy(buf, ((char *) &start_data) + (off - aligned_start),
				first_size);
	for (int i = first_size; i < size; i += sizeof(off_t)) {
		*((long *) (buf + i)) = (off + i) / sizeof(off_t);
	}
	if (aligned_end > aligned_start) {
		int last_size = (int) (off + size - aligned_end);
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
public:
	cleanup_callback(rand_buf *buf, int idx) {
		this->buf = buf;
		read_bytes = 0;
		this->thread_id = idx;
	}

	int invoke(io_request *rq) {
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

int thread_private::thread_init() {
	attach2cpu();
	io->init();

	extern int buf_size;
	rand_buf *buf = new rand_buf(NUM_PAGES / (nthreads
				// TODO maybe I should set the right entry size for a buffer.
				// If each access size is irregular, I'll break each access
				// into pages so each access is no larger than a page, so it
				// should workl fine.
				/ NUM_NODES) * PAGE_SIZE, buf_size);
	this->buf = buf;
	if (io->support_aio()) {
		cb = new cleanup_callback(buf, idx);
		io->set_callback(cb);
	}
	return 0;
}

int thread_private::run()
{
	ssize_t ret = -1;
	gettimeofday(&start_time, NULL);
	io_request reqs[NUM_REQS_BY_USER];
	char *entry = NULL;
	if (!io->support_aio()) {
		extern int buf_size;
		entry = (char *) valloc(buf_size);
	}
	while (gen->has_next()) {
		if (io->support_aio()) {
			int i;
			for (i = 0; i < NUM_REQS_BY_USER && gen->has_next(); ) {
				workload_t workload = gen->next();
				int access_method = workload.read ? READ : WRITE;
				off_t off = workload.off;
				int size = workload.size;
				/*
				 * If the size of the request is larger than a page size,
				 * and the user explicitly wants to use multibuf requests.
				 */
				if (buf_type == MULTI_BUF) {
					assert(off % PAGE_SIZE == 0);
					int num_vecs = size / PAGE_SIZE;
					reqs[i].init(off, io, access_method);
					assert(buf->get_entry_size() >= PAGE_SIZE);
					for (int k = 0; k < num_vecs; k++) {
						reqs[i].add_buf(buf->next_entry(), PAGE_SIZE);
					}
					i++;
				}
				else if (buf_type == SINGLE_SMALL_BUF) {
again:
					while (size > 0 && i < NUM_REQS_BY_USER) {
						off_t next_off = ROUNDUP_PAGE(off + 1);
						if (next_off > off + size)
							next_off = off + size;
						char *p = buf->next_entry();
						if (access_method == WRITE)
							create_write_data(p, next_off - off, off);
						reqs[i].init(p, off, next_off - off, access_method, io);
						size -= next_off - off;
						off = next_off;
						i++;
					}
					if (size > 0) {
						ret = io->access(reqs, i);
						i = 0;
						goto again;
					}
				}
				else {
					assert(buf->get_entry_size() >= size);
					char *p = buf->next_entry();
					if (access_method == WRITE)
						create_write_data(p, size, off);
					reqs[i++].init(p, off, size, access_method, io);
				}
			}
			ret = io->access(reqs, i);
			if (ret < 0) {
				perror("access_vector");
				exit(1);
			}
		}
		else {
			workload_t workload = gen->next();
			off_t off = workload.off;
			// TODO let's just read data first.
			int access_method = workload.read ? READ : WRITE;
			int entry_size = workload.size;

			if (buf_type == SINGLE_SMALL_BUF) {
				while (entry_size > 0) {
					/*
					 * generate the data for writing the file,
					 * so the data in the file isn't changed.
					 */
					if (access_method == WRITE) {
						create_write_data(entry, entry_size, off);
					}
					// There is at least one byte we need to access in the page.
					// By adding 1 and rounding up the offset, we'll get the next page
					// behind the current offset.
					off_t next_off = ROUNDUP_PAGE(off + 1);
					if (next_off > off + entry_size)
						next_off = off + entry_size;
					ret = io->access(entry, off, next_off - off, access_method);
					if (ret > 0) {
						if (access_method == READ && verify_read_content) {
							check_read_content(entry, next_off - off, off);
						}
						read_bytes += ret;
					}
					if (ret < 0) {
						perror("access");
						exit(1);
					}
					entry_size -= next_off - off;
					off = next_off;
				}
			}
			else {
				if (access_method == WRITE) {
					create_write_data(entry, entry_size, off);
				}
				ret = io->access(entry, off, entry_size, access_method);
				if (ret > 0) {
					if (access_method == READ && verify_read_content) {
						check_read_content(entry, entry_size, off);
					}
					read_bytes += ret;
				}
				if (ret < 0) {
					perror("access");
					exit(1);
				}
			}
		}
	}
	io->cleanup();
	printf("thread %d exits\n", idx);
	gettimeofday(&end_time, NULL);
	return 0;
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
	printf("attach thread %d to CPU %d\n", idx, cpu_num);
	return ret;
#else
	return -1;
#endif
}

static void *rand_read(void *arg)
{
	thread_private *priv = (thread_private *) arg;

	printf("rand_read: pid: %d, tid: %ld\n", getpid(), gettid());
	priv->thread_init();
	priv->run();
	return NULL;
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

int thread_private::start_thread()
{
	int ret;
#ifdef USE_PROCESS
	ret = process_create(&id, rand_read, (void *) this);
#else
	ret = pthread_create(&id, NULL, rand_read, (void *) this);
#endif
	return ret;
}

int thread_private::wait_thread_end()
{
	int ret;
#ifdef USE_PROCESS
	ret = process_join(id);
#else
	ret = pthread_join(id, NULL);
#endif
	return ret;
}
