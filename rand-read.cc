#include <stdio.h>
#include <unistd.h>
#include <sys/types.h>
#include <sys/wait.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <sys/time.h>
#include <stdlib.h>
#include <sys/resource.h>
#include <sys/mman.h>
#include <string.h>
#include <errno.h>
#include <pthread.h>
#include <assert.h>
#include <numa.h>
#include <numaif.h>
#ifdef PROFILER
#include <google/profiler.h>
#endif
#include <sys/syscall.h>
#define gettid() syscall(__NR_gettid)

#include <iostream>
#include <string>
#include <deque>

#include "wpaio.h"
#include "cache.h"
#include "tree_cache.h"
#include "associative_cache.h"
#include "cuckoo_cache.h"
#include "workload.h"
#include "LRU2Q.h"

//#define USE_PROCESS

#define NUM_PAGES 16384
#define NUM_THREADS 32
#define AIO_DEPTH 32
#define ENTRY_READ_SIZE 128

enum {
	NORMAL,
	DIRECT,
	MMAP,
};

enum {
	READ,
	WRITE
};

enum {
	TREE_CACHE,
	ASSOCIATIVE_CACHE,
	CUCKOO_CACHE,
	LRU2Q_CACHE,
};

long npages;
int nthreads = 1;
struct timeval global_start;
char static_buf[PAGE_SIZE * 8] __attribute__((aligned(PAGE_SIZE)));
volatile int first[NUM_THREADS];
int access_method = READ;

float time_diff(struct timeval time1, struct timeval time2)
{
	return time2.tv_sec - time1.tv_sec
			+ ((float)(time2.tv_usec - time1.tv_usec))/1000000;
}

class rand_buf
{
	/* where the data read from the disk is stored */
	char *buf;
	/* shows the locations in the array where data has be to stored.*/
	rand_permute buf_offset;
	int entry_size;
	int num_entries;

	int current;
public:
	rand_buf(int buf_size, int entry_size): buf_offset(buf_size / entry_size, entry_size) {
		this->entry_size = entry_size;
		num_entries = buf_size / entry_size;
		buf = (char *) valloc(buf_size);

		if (buf == NULL){
			fprintf(stderr, "can't allocate buffer\n");
			exit(1);
		}
		/* trigger page faults and bring pages to memory. */
		for (int i = 0; i < buf_size / PAGE_SIZE; i++)
			buf[i * PAGE_SIZE] = 0;

		current = 0;
	}

	~rand_buf() {
		free(buf);
	}

	char *next_entry() {
		int off = buf_offset.get_offset(current);
		current = (current + 1) % num_entries;;
		return &buf[off];
	}

	int get_entry_size() {
		return entry_size;
	}
};

/* this data structure stores the thread-private info. */
class thread_private
{
public:
#ifdef USE_PROCESS
	pid_t id;
#else
	pthread_t id;
#endif
	/* the location in the thread descriptor array. */
	int idx;
	rand_buf buf;
	workload_gen *gen;
	ssize_t read_bytes;
	struct timeval start_time;
	struct timeval end_time;

#ifdef STATISTICS
	int cache_hits;
#endif

	virtual ssize_t access(char *, off_t, ssize_t, int) = 0;
	virtual int thread_init() = 0;

	thread_private(int idx, int entry_size): buf(NUM_PAGES / nthreads * PAGE_SIZE, entry_size) {
		this->idx = idx;
		read_bytes = 0;
#ifdef STATISTICS
		cache_hits = 0;
#endif
	}

#ifdef STATISTICS
	virtual void print_stat() {
		static int tot_hits = 0;
		static int seen_threads = 0;
		seen_threads++;
		tot_hits += cache_hits;
		printf("read %ld bytes, start at %f seconds, takes %f seconds\n",
				read_bytes, time_diff(global_start, start_time),
				time_diff(start_time, end_time));
		if (seen_threads == nthreads)
			printf("there are %d cache hits\n", tot_hits);
	}
#endif
};

class read_private: public thread_private
{
	/* the array of files that it's going to access */
	const char **file_names;
	int *fds;
	/* the number of files */
	int num;
	/* the size of data it's going to access, and it'll be divided for each file */
	long size;

	int flags;
#ifdef STATISTICS
	long read_time; // in us
	long num_reads;
#endif
public:
	read_private(const char *names[], int num, long size, int idx, int entry_size,
			int flags = O_RDWR): thread_private(idx, entry_size) {
		this->flags = flags;
#ifdef STATISTICS
		read_time = 0;
		num_reads = 0;
#endif
		file_names = new const char *[num];
		for (int i = 0; i < num; i++)
			file_names[i] = names[i];
		fds = new int[num];
		this->num = num;
		this->size = size;
	}

	~read_private() {
		delete [] file_names;
		delete [] fds;
	}

	int thread_init() {
		int ret;

		for (int i = 0; i < num; i++) {
			fds[i] = open(file_names[i], flags);
			if (fds[i] < 0) {
				perror("open");
				exit (1);
			}
			ret = posix_fadvise(fds[i], 0, 0, POSIX_FADV_RANDOM);
			if (ret < 0) {
				perror("posix_fadvise");
				exit(1);
			}
		}
		return 0;
	}

	int thread_end() {
		for (int i = 0; i < num; i++)
			close(fds[i]);
		return 0;
	}

	ssize_t access(char *buf, off_t offset, ssize_t size, int access_method) {
		int fd_idx = offset / (this->size / num);
		if (fd_idx >= num) {
			printf("offset: %ld, fd_idx: %d, size: %ld, num: %d\n", offset, fd_idx, this->size, num);
		}
		assert (fd_idx < num);
		int fd = fds[fd_idx];
#ifdef STATISTICS
		if (access_method == READ)
			num_reads++;
		struct timeval start, end;
		gettimeofday(&start, NULL);
#endif
		assert(offset < 0x100000000L);
		ssize_t ret;
		if (access_method == WRITE)
			ret = pwrite(fd, buf, size, offset);
		else
			ret = pread(fd, buf, size, offset);
#ifdef STATISTICS
		gettimeofday(&end, NULL);
		read_time += ((long) end.tv_sec - start.tv_sec) * 1000000 + end.tv_usec - start.tv_usec;
#endif
		return ret;
	}

#ifdef STATISTICS
	virtual void print_stat() {
		thread_private::print_stat();
		static int seen_threads = 0;
		static long tot_nreads;
		static long tot_read_time;
		tot_nreads += num_reads;
		tot_read_time += read_time;
		seen_threads++;
		if (seen_threads == nthreads)
			printf("there are %ld reads and takes %ldus\n", tot_nreads, tot_read_time);
	}
#endif
};

class direct_private: public read_private
{
	char *pages;
	int buf_idx;
public:
	direct_private(const char *names[], int num, long size, int idx,
			int entry_size): read_private(names, num, size, idx,
			entry_size, O_DIRECT | O_RDWR) {
		pages = (char *) valloc(PAGE_SIZE * 4096);
		buf_idx = 0;
	}

	ssize_t access(char *buf, off_t offset, ssize_t size, int access_method) {
		ssize_t ret;
		/* for simplicity, I assume all request sizes are smaller than a page size */
		assert(size <= PAGE_SIZE);
		if (ROUND_PAGE(offset) == offset
				&& (long) buf == ROUND_PAGE(buf)
				&& size == PAGE_SIZE) {
			ret = read_private::access(buf, offset, size, access_method);
		}
		else {
			assert(access_method == READ);
			buf_idx++;
			if (buf_idx == 4096)
				buf_idx = 0;
			char *page = pages + buf_idx * PAGE_SIZE;
			ret = read_private::access(page, ROUND_PAGE(offset), PAGE_SIZE, access_method);
			if (ret < 0)
				return ret;
			else
				memcpy(buf, page + (offset - ROUND_PAGE(offset)), size);
			ret = size;
		}
		return ret;
	}
};

#ifdef ENABLE_AIO
class aio_private;
void aio_callback(io_context_t, struct iocb*,
		struct io_callback_s *, long, long);

struct thread_callback_s
{
	struct io_callback_s cb;
	aio_private *thread;
};

class aio_private: public read_private
{
	char *pages;
	int buf_idx;
	struct aio_ctx *ctx;
	std::deque<thread_callback_s *> cbs;

public:
	aio_private(const char *names[], int num, long size, int idx,
			int entry_size): read_private(names, num, size, idx,
			entry_size, O_DIRECT | O_RDWR) {
		printf("aio is used\n");
		pages = (char *) valloc(PAGE_SIZE * 4096);
		buf_idx = 0;
		ctx = create_aio_ctx(AIO_DEPTH);
		for (int i = 0; i < AIO_DEPTH * 5; i++) {
			cbs.push_back(new thread_callback_s());
		}
	}

	~aio_private() {
		int slot = max_io_slot(ctx);

		while (slot < AIO_DEPTH) {
			io_wait(ctx, NULL);
			slot = max_io_slot(ctx);
		}
	}

	ssize_t access(char *buf, off_t offset, ssize_t size, int access_method) {
		struct iocb *req;
		int slot = max_io_slot(ctx);

		assert(access_method == READ);
		if (slot == 0) {
			io_wait(ctx, NULL);
		}

		if (cbs.empty()) {
			fprintf(stderr, "no callback object left\n");
			return -1;
		}

		thread_callback_s *tcb = cbs.front();
		io_callback_s *cb = (io_callback_s *) tcb;
		cbs.pop_front();
		cb->buf = buf;
		cb->offset = offset;
		cb->size = size;
		cb->func = aio_callback;
		tcb->thread = this;

		/* for simplicity, I assume all request sizes are smaller than a page size */
		assert(size <= PAGE_SIZE);
		if (ROUND_PAGE(offset) == offset
				&& (long) buf == ROUND_PAGE(buf)
				&& size == PAGE_SIZE) {
			req = make_io_request(ctx, get_fd(), PAGE_SIZE, offset,
					buf, A_READ, cb);
		}
		else {
			buf_idx++;
			if (buf_idx == 4096)
				buf_idx = 0;
			char *page = pages + buf_idx * PAGE_SIZE;
			req = make_io_request(ctx, get_fd(), PAGE_SIZE, ROUND_PAGE(offset),
					page, A_READ, cb);
		}
		submit_io_request(ctx, &req, 1);
		return 0;
	}

	void return_cb(thread_callback_s *cb) {
		cbs.push_back(cb);
	}
};

void aio_callback(io_context_t ctx, struct iocb* iocb,
		struct io_callback_s *cb, long res, long res2) {
	thread_callback_s *tcb = (thread_callback_s *) cb;

	memcpy(cb->buf, ((char *) iocb->u.c.buf)
			+ (cb->offset - ROUND_PAGE(cb->offset)), cb->size);
	if(*(unsigned long *) cb->buf != cb->offset / sizeof(long))
		printf("%ld %ld\n", *(unsigned long *) cb->buf, cb->offset / sizeof(long));
	assert(*(unsigned long *) cb->buf == cb->offset / sizeof(long));
	tcb->thread->return_cb(tcb);
	tcb->thread->read_bytes += cb->size;
}
#endif

class mmap_private: public thread_private
{
	void *addr;

public:
	mmap_private(const char *new_name, int idx,
			int entry_size): thread_private(idx, entry_size) {
		static void *addr = NULL;
		static const char *file_name = NULL;
		/* if we are mapping to a different file, do the real mapping. */
		if (file_name == NULL || strcmp(file_name, new_name)) {
			int fd = open(new_name, O_RDWR);
			int ret;

			if (fd < 0) {
				perror("open");
				exit (1);
			}
			addr = mmap(NULL, ((ssize_t) npages) * PAGE_SIZE,
					PROT_READ, MAP_PRIVATE, fd, 0);
			if (addr == NULL) {
				perror("mmap");
				exit(1);
			}
			ret = madvise(addr, ((ssize_t) npages) * PAGE_SIZE, MADV_RANDOM);
			if (ret < 0) {
				perror("madvise");
				exit(1);
			}
			file_name = new_name;
		}
		this->addr = addr;
	}

	int thread_init() {
		return 0;
	}

	ssize_t access(char *buf, off_t offset, ssize_t size, int access_method) {
		char *page = (char *) addr + offset;
		/* I try to avoid gcc optimization eliminating the code below,
		 * and it works. */
		first[idx] = page[0];
		__asm__ __volatile__("" : : "g"(&first[idx]));
		return size;
	}
};

class part_cached_private: public direct_private
{
	page_cache *cache;

public:
	part_cached_private(const char *names[], int num, long size, int idx, long cache_size,
			int entry_size): direct_private(names, num, size, idx, entry_size) {
		/* all local cache has the same size */
		cache = new tree_cache(cache_size, idx * cache_size);
	}

	ssize_t access(char *buf, off_t offset, ssize_t size, int access_method) {
		ssize_t ret;
		off_t old_off = -1;
		assert(access_method == READ);
		page *p = cache->search(ROUND_PAGE(offset), old_off);

		assert(p->get_offset() == ROUND_PAGE(offset));
		if (!p->data_ready()) {
			ret = read_private::access((char *) p->get_data(),
					ROUND_PAGE(offset), PAGE_SIZE, access_method);
			if (ret < 0) {
				perror("read");
				exit(1);
			}
			p->set_data_ready(true);
		}
		else {
#ifdef STATISTICS
			cache_hits++;
#endif
		}

		offset -= ROUND_PAGE(offset);
		/* I assume the data I read never crosses the page boundary */
		memcpy(buf, (char *) p->get_data() + offset, size);
		ret = size;
		return ret;
	}
};

page_cache *global_cache;
class global_cached_private: public direct_private
{
	int num_waits;
	long cache_size;
//	static page_cache *global_cache;
public:
	global_cached_private(const char *names[], int num, long size, int idx, long cache_size,
			int entry_size, int cache_type): direct_private(names, num, size, idx, entry_size) {
		num_waits = 0;
		this->cache_size = cache_size;
		if (global_cache == NULL) {
			switch (cache_type) {
				case TREE_CACHE:
					global_cache = new tree_cache(cache_size, 0);
					break;
				case ASSOCIATIVE_CACHE:
					global_cache = new associative_cache(cache_size);
					break;
				case CUCKOO_CACHE:
					global_cache = new cuckoo_cache(cache_size);
					break;
				case LRU2Q_CACHE:
					global_cache = new LRU2Q_cache(cache_size);
					break;
				default:
					fprintf(stderr, "wrong cache type\n");
					exit(1);
			}
		}
	}

	int preload(off_t start, long size) {
		if (size > cache_size) {
			fprintf(stderr, "we can't preload data larger than the cache size\n");
			exit(1);
		}

		/* open the file. It's a hack, but it works for now. */
		thread_init();

		assert(ROUND_PAGE(start) == start);
		for (long offset = start; offset < start + size; offset += PAGE_SIZE) {
			off_t old_off = -1;
			thread_safe_page *p = (thread_safe_page *) (global_cache->search(ROUND_PAGE(offset), old_off));
			if (!p->data_ready()) {
				ssize_t ret = read_private::access((char *) p->get_data(),
						ROUND_PAGE(offset), PAGE_SIZE, READ);
				if (ret < 0) {
					perror("read");
					return ret;
				}
				p->set_io_pending(false);
				p->set_data_ready(true);
			}
		}
		/* close the file as it will be opened again in the real workload. */
		thread_end();
		return 0;
	}

	ssize_t access(char *buf, off_t offset, ssize_t size, int access_method) {
		ssize_t ret;
		off_t old_off = -1;
		thread_safe_page *p = (thread_safe_page *) (global_cache->search(ROUND_PAGE(offset), old_off));

		/*
		 * the page isn't in the cache,
		 * so the cache evict a page and return it to us.
		 */
		if (old_off != ROUND_PAGE(offset) && old_off != -1) {
			/* 
			 * if the new page we get is dirty,
			 * we need to write its data back to the file first.
			 */
			if (p->is_dirty()) {
				unsigned long *l = (unsigned long *) p->get_data();
				unsigned long start = old_off / sizeof(long);
				if (*l != start)
					printf("start: %ld, l: %ld\n", start, *l);
				p->lock();
				read_private::access((char *) p->get_data(),
						old_off, PAGE_SIZE, WRITE);
				p->set_dirty(false);
				p->unlock();
			}
		}

		if (!p->data_ready()) {
			/* if the page isn't io pending,
			 * set it io pending, and return
			 * original result. otherwise,
			 * just return the original value.
			 *
			 * This is an ugly hack, but with this
			 * atomic operation, I can avoid using
			 * locks. */
			if(!p->test_and_set_io_pending()) {
				p->lock();
				/* 
				 * we have to make sure the page is clean
				 * before reading data to the page.
				 */
				while(p->is_dirty()) {}
				ret = read_private::access((char *) p->get_data(),
						ROUND_PAGE(offset), PAGE_SIZE, READ);
				p->unlock();
				if (ret < 0) {
					perror("read");
					exit(1);
				}
				p->set_io_pending(false);
				p->set_data_ready(true);
			}
			else {
				num_waits++;
				p->wait_ready();
			}
		}
		else {
#ifdef STATISTICS
			cache_hits++;
#endif
		}
		int page_off = offset - ROUND_PAGE(offset);
		p->lock();
		if (access_method == WRITE) {
			unsigned long *l = (unsigned long *) ((char *) p->get_data() + page_off);
			unsigned long start = offset / sizeof(long);
			if (*l != start)
				printf("write: start: %ld, l: %ld, offset: %ld\n",
						start, *l, p->get_offset() / sizeof(long));
			memcpy((char *) p->get_data() + page_off, buf, size);
			p->set_dirty(true);
		}
		else 
			/* I assume the data I read never crosses the page boundary */
			memcpy(buf, (char *) p->get_data() + page_off, size);
		p->unlock();
		p->dec_ref();
		ret = size;
		return ret;
	}

#ifdef STATISTICS
	void print_stat() {
		direct_private::print_stat();
		printf("there are %d waits in thread %d\n", num_waits, idx);
	}
#endif
};

thread_private *threads[NUM_THREADS];

void *rand_read(void *arg)
{
	ssize_t ret = -1;
	thread_private *priv = threads[(long) arg];
	rand_buf *buf;

	printf("pid: %d, tid: %ld\n", getpid(), gettid());
	priv->thread_init();
	buf = &priv->buf;

	gettimeofday(&priv->start_time, NULL);
	while (priv->gen->has_next()) {
		char *entry = buf->next_entry();
		off_t off = priv->gen->next_offset();

		/*
		 * generate the data for writing the file,
		 * so the data in the file isn't changed.
		 */
		if (access_method == WRITE) {
			unsigned long *p = (unsigned long *) entry;
			long start = off / sizeof(long);
			unsigned int entry_size = buf->get_entry_size();
			for (unsigned int i = 0; i < entry_size / sizeof(*p); i++)
				p[i] = start++;
		}

		ret = priv->access(entry, off, ENTRY_READ_SIZE, access_method);
		if (ret > 0) {
			assert(ret == buf->get_entry_size());
			if (access_method == READ)
				assert(*(unsigned long *) entry == off / sizeof(long));
			if (ret > 0)
				priv->read_bytes += ret;
			else
				break;
		}
		if (ret < 0) {
			perror("access");
			exit(1);
		}
	}
	if (ret < 0) {
		perror("read");
		exit(1);
	}
	gettimeofday(&priv->end_time, NULL);
	
#ifdef USE_PROCESS
	exit(priv->read_bytes);
#else
	pthread_exit((void *) priv->read_bytes);
#endif
}

int process_create(pid_t *pid, void (*func)(void *), void *priv)
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

int process_join(pid_t pid)
{
	int status;
	pid_t ret = waitpid(pid, &status, 0);
	return ret < 0 ? ret : 0;
}

struct str2int {
	std::string name;
	int value;
};

enum {
	READ_ACCESS,
	DIRECT_ACCESS,
#ifdef ENABLE_AIO
	AIO_ACCESS,
#endif
	MMAP_ACCESS,
	LOCAL_CACHE_ACCESS,
	GLOBAL_CACHE_ACCESS
};

str2int access_methods[] = {
	{ "normal", READ_ACCESS },
	{ "direct", DIRECT_ACCESS },
#ifdef ENABLE_AIO
	{ "aio", AIO_ACCESS },
#endif
	{ "mmap", MMAP_ACCESS },
	{ "local_cache", LOCAL_CACHE_ACCESS },
	{ "global_cache", GLOBAL_CACHE_ACCESS },
};

str2int cache_types[] = {
	{ "tree", TREE_CACHE } ,
	{ "associative", ASSOCIATIVE_CACHE },
	{ "cuckoo", CUCKOO_CACHE },
	{ "lru2q", LRU2Q_CACHE },
};

enum {
	RAND_OFFSET,
	RAND_PERMUTE,
	STRIDE,
	BALANCED,
	USER_FILE_WORKLOAD = -1
};

str2int workloads[] = {
	{ "RAND", RAND_OFFSET },
	{ "RAND_PERMUTE", RAND_PERMUTE },
	{ "STRIDE", STRIDE },
	{ "BALANCED", BALANCED },
	{ "user_file", USER_FILE_WORKLOAD },
};

class str2int_map {
	str2int *maps;
	int num;
public:
	str2int_map(str2int *maps, int num) {
		this->maps = maps;
		this->num = num;
	}

	int map(const std::string &str) {
		for (int i = 0; i < num; i++) {
			if (maps[i].name.compare(str) == 0)
				return i;
		}
		return -1;
	}

	void print(const std::string &str) {
		std::cout<<str;
		for (int i = 0; i < num; i++) {
			std::cout<<maps[i].name;
			if (i < num - 1)
				std::cout<<", ";
		}
		std::cout<<std::endl;
	}
};

long str2size(std::string str)
{
	int len = str.length();
	long multiply = 1;
	if (str[len - 1] == 'M' || str[len - 1] == 'm') {
		multiply *= 1024 * 1024;
		str[len - 1] = 0;
	}
	else if (str[len - 1] == 'K' || str[len - 1] == 'k') {
		multiply *= 1024;
		str[len - 1] = 0;
	}
	else if (str[len - 1] == 'G' || str[len - 1] == 'g') {
		multiply *= 1024 * 1024 * 1024;
		str[len - 1] = 0;
	}
	return atol(str.c_str()) * multiply;
}

int main(int argc, char *argv[])
{
	bool preload = false;
	long cache_size = 512 * 1024 * 1024;
	int entry_size = 128;
	int access_option = -1;
	int ret;
	int i, j;
	struct timeval start_time, end_time;
	ssize_t read_bytes = 0;
	int num_files = 0;
	int cache_type = -1;
	std::string file_names[NUM_THREADS];
	int workload = RAND_OFFSET;
	std::string workload_file;
	str2int_map access_map(access_methods,
			sizeof(access_methods) / sizeof(access_methods[0]));
	str2int_map workload_map(workloads, 
			sizeof(workloads) / sizeof(workloads[0]));
	str2int_map cache_map(cache_types, 
			sizeof(cache_types) / sizeof(cache_types[0]));
#ifdef PROFILER
	std::string prof_file;
#endif

	if (argc < 5) {
		fprintf(stderr, "there are %d argments\n", argc);
		fprintf(stderr, "read files option pages threads cache_size entry_size preload workload cache_type\n");
		access_map.print("available access options: ");
		workload_map.print("available workloads: ");
		cache_map.print("available cache types: ");
		exit(1);
	}

	for (int i = 1; i < argc; i++) {
		std::string str = argv[i];
		size_t found = str.find("=");
		/* if there isn't `=', I assume it's a file name*/
		if (found == std::string::npos) {
			file_names[num_files++] = str;
			continue;
		}

		std::string value = str.substr(found + 1);
		std::string key = str.substr(0, found);
		if (key.compare("option") == 0) {
			access_option = access_map.map(value);
			if (access_option < 0) {
				fprintf(stderr, "can't find the right option\n");
				exit(1);
			}
		}
		else if(key.compare("cache_type") == 0) {
			cache_type = cache_map.map(value);
		}
		else if(key.compare("pages") == 0) {
			npages = atoi(value.c_str());
		}
		else if(key.compare("threads") == 0) {
			nthreads = atoi(value.c_str());
		}
		else if(key.compare("cache_size") == 0) {
			cache_size = str2size(value);
		}
		else if(key.compare("entry_size") == 0) {
			entry_size = (int) str2size(value);
		}
		else if(key.compare("workload") == 0) {
			workload = workload_map.map(value);
			if (workload == -1) {
				workload_file = value;
			}
		}
		else if(key.compare("access") == 0) {
			if(value.compare("read") == 0)
				access_method = READ;
			else if(value.compare("write") == 0)
				access_method = WRITE;
			else {
				fprintf(stderr, "wrong access method\n");
				exit(1);
			}
		}
		else if(key.compare("preload") == 0) {
			preload = true;
		}
#ifdef PROFILER
		else if(key.compare("prof") == 0) {
			prof_file = value;
		}
#endif
		else {
			fprintf(stderr, "wrong option\n");
			exit(1);
		}
	}
	printf("access: %d, npages: %ld, nthreads: %d, cache_size: %ld, cache_type: %d, entry_size: %d, workload: %d\n",
			access_option, npages, nthreads, cache_size, cache_type, entry_size, workload);

	int num_entries = npages * (PAGE_SIZE / entry_size);

	if (nthreads > NUM_THREADS) {
		fprintf(stderr, "too many threads\n");
		exit(1);
	}

	page::allocate_cache(cache_size);
	int remainings = npages % nthreads;
	int shift = 0;
	long start;
	long end = 0;
	const char *cnames[num_files];
	int num;
	for (int k = 0; k < num_files; k++)
		cnames[k] = file_names[k].c_str();
	num = num_files;
	/* initialize the threads' private data. */
	for (j = 0; j < nthreads; j++) {
		switch (access_option) {
			case READ_ACCESS:
				threads[j] = new read_private(cnames, num, npages * PAGE_SIZE, j, entry_size);
				break;
			case DIRECT_ACCESS:
				threads[j] = new direct_private(cnames, num, npages * PAGE_SIZE, j, entry_size);
				break;
#if ENABLE_AIO
			case AIO_ACCESS:
				threads[j] = new aio_private(cnames, num, npages * PAGE_SIZE, j, entry_size);
				break;
#endif
			case MMAP_ACCESS:
				/* TODO for now mmap doesn't support accessing multiple files. */
				threads[j] = new mmap_private(cnames[0], j, entry_size);
				break;
			case LOCAL_CACHE_ACCESS:
				threads[j] = new part_cached_private(cnames, num, npages * PAGE_SIZE,
						j, cache_size / nthreads, entry_size);
				break;
			case GLOBAL_CACHE_ACCESS:
				threads[j] = new global_cached_private(cnames, num, npages * PAGE_SIZE, j,
						cache_size, entry_size, cache_type);
				if (preload)
					((global_cached_private *) threads[j])->preload(0, npages * PAGE_SIZE);
				break;
			default:
				fprintf(stderr, "wrong access option\n");
				exit(1);
		}
		
		if (num_files > 1) {
			start = 0;
			end = npages * PAGE_SIZE / entry_size;
		}
		else {
			start = end;
			end = start + ((long) npages / nthreads + (shift < remainings)) * PAGE_SIZE / entry_size;
			if (remainings != shift)
				shift++;
		}
		printf("thread %d starts %ld ends %ld\n", j, start, end);

		workload_gen *gen;
		workload_chunk *chunk = NULL;
		switch (workload) {
			case RAND_OFFSET:
				gen = new rand_workload(start, end, entry_size);
				break;
			case RAND_PERMUTE:
				gen = new local_rand_permute_workload(num_entries,
						entry_size, start, end);
				break;
			case STRIDE:
				gen = new stride_workload(start, end, entry_size);
				break;
			case BALANCED:
				if (chunk == NULL) {
					chunk = new stride_workload_chunk(0,
							(long) npages * PAGE_SIZE / entry_size, entry_size);
				}
				gen = new balanced_workload(chunk);
				break;
			case -1:
				gen = new file_workload(workload_file, nthreads);
				break;
			default:
				fprintf(stderr, "unsupported workload\n");
				exit(1);
		}
		threads[j]->gen = gen;
	}

	ret = setpriority(PRIO_PROCESS, getpid(), -20);
	if (ret < 0) {
		perror("setpriority");
		exit(1);
	}

	gettimeofday(&start_time, NULL);
	global_start = start_time;
#ifdef PROFILER
	if (!prof_file.empty())
		ProfilerStart(prof_file.c_str());
#endif
	for (i = 0; i < nthreads; i++) {
#ifdef USE_PROCESS
		ret = process_create(&threads[i]->id, rand_read, (void *) i);
#else
		ret = pthread_create(&threads[i]->id, NULL, rand_read, (void *) i);
#endif
		if (ret) {
			perror("pthread_create");
			exit(1);
		}
	}

	for (i = 0; i < nthreads; i++) {
		ssize_t size;
#ifdef USE_PROCESS
		ret = process_join(threads[i]->id);
#else
		ret = pthread_join(threads[i]->id, (void **) &size);
#endif
		if (ret) {
			perror("pthread_join");
			exit(1);
		}
		read_bytes += size;
	}
#ifdef PROFILER
	if (!prof_file.empty())
		ProfilerStop();
#endif
	gettimeofday(&end_time, NULL);
	printf("read %ld bytes, takes %f seconds\n",
			read_bytes, end_time.tv_sec - start_time.tv_sec
			+ ((float)(end_time.tv_usec - start_time.tv_usec))/1000000);

#ifdef STATISTICS
	for (int i = 0; i < nthreads; i++) {
		threads[i]->print_stat();
	}
	printf("there are %d cells\n", avail_cells);
	printf("there are %d waits for unused\n", num_wait_unused);
	printf("%d keys are evicted from the hash table because of conflicts\n", removed_indices);
	printf("there are %d lock contentions\n", lock_contentions);
#endif
}

const rand_permute *local_rand_permute_workload::permute;
