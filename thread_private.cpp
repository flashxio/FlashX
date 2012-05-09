#include <sys/time.h>
#include <sys/types.h>
#include <sys/syscall.h>
#define gettid() syscall(__NR_gettid)

#include "thread_private.h"

#define BULK_SIZE 1000

extern int access_method;
extern bool verify_read_content;

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
			if(*(unsigned long *) rq->get_buf() != rq->get_offset() / sizeof(long))
				printf("%ld %ld\n", *(unsigned long *) rq->get_buf(),
						rq->get_offset() / sizeof(long));
			assert(*(unsigned long *) rq->get_buf()
					== rq->get_offset() / sizeof(long));
		}
//		printf("thread %d: free %p\n", thread_id, rq->get_buf());
		buf->free_entry(rq->get_buf());
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

	rand_buf *buf = new rand_buf(NUM_PAGES / (nthreads
				/ NUM_NODES) * PAGE_SIZE, get_entry_size());
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
	while (gen->has_next()) {
		if (io->support_aio()) {
			io_request reqs[BULK_SIZE];
			int i;
//			io_request *reqs = gc->allocate_obj(BULK_SIZE);
			for (i = 0; i < BULK_SIZE
					&& gen->has_next() && !buf->is_full(); i++) {
				char *p = buf->next_entry();
//				printf("thread %d: allocate %p\n", idx, p);
				// TODO right now it only support read.
				reqs[i].init(p, gen->next_offset(), entry_size, READ, io);
			}
			ret = io->access(reqs, i);
			if (ret < 0) {
				perror("access_vector");
				exit(1);
			}
		}
		else {
			char *entry = buf->next_entry();
			off_t off = gen->next_offset();

			/*
			 * generate the data for writing the file,
			 * so the data in the file isn't changed.
			 */
			if (access_method == WRITE) {
				unsigned long *p = (unsigned long *) entry;
				long start = off / sizeof(long);
				for (unsigned int i = 0; i < entry_size / sizeof(*p); i++)
					p[i] = start++;
			}

			ret = io->access(entry, off, entry_size, access_method);
			if (ret > 0) {
//				assert(ret == buf->get_entry_size());
				if (access_method == READ && verify_read_content) {
					if (*(unsigned long *) entry != off / sizeof(long))
						printf("entry: %ld, off: %ld\n", *(unsigned long *) entry, off / sizeof(long));
					assert(*(unsigned long *) entry == off / sizeof(long));
				}
				if (ret > 0)
					read_bytes += ret;
				else
					break;
			}
			buf->free_entry(entry);
			if (ret < 0) {
				perror("access");
				exit(1);
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

	printf("pid: %d, tid: %ld\n", getpid(), gettid());
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
