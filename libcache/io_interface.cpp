#include <assert.h>
#include <signal.h>
#include <time.h>
#include <pthread.h>

#include <tr1/unordered_set>
#include <vector>

#include "RAID_config.h"
#include "io_interface.h"
#include "read_private.h"
#include "direct_private.h"
#include "aio_private.h"
#include "remote_access.h"
#include "global_cached_private.h"
#include "part_global_cached_private.h"
#include "cache_config.h"
#include "disk_read_thread.h"

#define DEBUG

/**
 * This global data collection is very static.
 * Once the data is initialized, no data needs to be changed.
 * The mutex is to used only at the initialization.
 * As long as all threads call init_io_system() first before using
 * the global data, they will all see the complete global data.
 */
struct global_data_collection
{
	RAID_config raid_conf;
	std::vector<disk_read_thread *> read_threads;
	pthread_mutex_t mutex;

#ifdef DEBUG
	std::tr1::unordered_set<io_interface *> ios;
	std::tr1::unordered_set<file_io_factory *> factories;
#endif

	global_data_collection() {
		pthread_mutex_init(&mutex, NULL);
	}

#ifdef DEBUG
	void register_io(io_interface *io) {
		ios.insert(io);
	}

	void register_io_factory(file_io_factory *factory) {
		factories.insert(factory);
	}
#endif
};

static global_data_collection global_data;

static bool enable_debug = false;

static void enable_debug_handler(int sig, siginfo_t *si, void *uc)
{
	enable_debug = true;
	printf("debug mode is enabled\n");
	for (unsigned i = 0; i < global_data.read_threads.size(); i++)
		global_data.read_threads[i]->print_state();
#ifdef DEBUG
	for (std::tr1::unordered_set<io_interface *>::iterator it
			= global_data.ios.begin(); it != global_data.ios.end(); it++)
		(*it)->print_state();
	for (std::tr1::unordered_set<file_io_factory *>::iterator it
			= global_data.factories.begin(); it != global_data.factories.end(); it++)
		(*it)->print_state();
#endif
}

bool is_debug_enabled()
{
	return enable_debug;
}

static void set_enable_debug_signal()
{
	struct sigaction sa;

	/* Establish handler for timer signal */

	sa.sa_flags = SA_SIGINFO;
	sa.sa_sigaction = enable_debug_handler;
	sigemptyset(&sa.sa_mask);
	if (sigaction(SIGUSR1, &sa, NULL) == -1) {
		perror("sigaction");
		exit(1);
	}
}

#if 0
/**
 * We flush notifications of IO completion buffered by I/O threads
 * periodically. By default, it's every 200ms.
 */
static void handler(int sig, siginfo_t *si, void *uc)
{
	// We can't call the flush requests to wake up I/O threads.
	// There might be a deadlock if the lock has been held by the thread
	// interrupted by the signal.
#if 0
	for (unsigned i = 0; i < global_data.read_threads.size(); i++) {
		global_data.read_threads[i]->flush_requests();
	}
#endif
}

static void set_completion_flush_timer()
{
	/**
	 * The code here is copied from the example code in the manual of
	 * timer_create.
	 */
	timer_t timerid;
	struct sigevent sev;
	struct itimerspec its;
	sigset_t mask;
	struct sigaction sa;

	/* Establish handler for timer signal */

	sa.sa_flags = SA_SIGINFO;
	sa.sa_sigaction = handler;
	sigemptyset(&sa.sa_mask);
	if (sigaction(SIGRTMIN, &sa, NULL) == -1) {
		perror("sigaction");
		exit(1);
	}

	/* Block timer signal temporarily */

	sigemptyset(&mask);
	sigaddset(&mask, SIGRTMIN);
	if (sigprocmask(SIG_SETMASK, &mask, NULL) == -1) {
		perror("sigprocmask");
		exit(1);
	}

	/* Create the timer */

	sev.sigev_notify = SIGEV_SIGNAL;
	sev.sigev_signo = SIGRTMIN;
	sev.sigev_value.sival_ptr = &timerid;
	if (timer_create(CLOCK_REALTIME, &sev, &timerid) == -1) {
		perror("timer_create");
		exit(1);
	}

	/* Start the timer */

	its.it_value.tv_sec = 0;
	its.it_value.tv_nsec = COMPLETION_FLUSH_INTERVAL;
	its.it_interval.tv_sec = its.it_value.tv_sec;
	its.it_interval.tv_nsec = its.it_value.tv_nsec;

	if (timer_settime(timerid, 0, &its, NULL) == -1) {
		perror("timer_settime");
		exit(1);
	}

	/* Unlock the timer signal, so that timer notification
	 * can be delivered */

	if (sigprocmask(SIG_UNBLOCK, &mask, NULL) == -1) {
		perror("sigprocmask");
		exit(1);
	}
}
#endif

void init_io_system(const std::string root_conf_file)
{
	numa_set_bind_policy(1);
	thread::thread_class_init();

	RAID_config raid_conf(root_conf_file, params.get_RAID_mapping_option(),
			params.get_RAID_block_size());
	int num_files = raid_conf.get_num_disks();
	global_data.raid_conf = raid_conf;

	/* 
	 * The mutex is enough to guarantee that all threads will see initialized
	 * global data. The first thread that enters the critical area will
	 * initialize the global data. If another thread tries to run the code,
	 * it will be blocked by the mutex. When a thread is returned from
	 * the function, they all can see the global data.
	 */
	pthread_mutex_lock(&global_data.mutex);
	// The global data hasn't been initialized.
	if (global_data.read_threads.size() == 0) {
		global_data.read_threads.resize(num_files);
		for (int k = 0; k < num_files; k++) {
			std::vector<int> indices(1, k);
			logical_file_partition partition(indices);
			// Create disk accessing threads.
			global_data.read_threads[k] = new disk_read_thread(partition,
					global_data.raid_conf.get_disk(k).node_id, NULL, k);
		}

		set_enable_debug_signal();
	}
	pthread_mutex_unlock(&global_data.mutex);
}

void destroy_io_system()
{
	// TODO
}

class posix_io_factory: public file_io_factory
{
	int access_option;
public:
	posix_io_factory(const std::string &file_name,
			int access_option): file_io_factory(file_name) {
		this->access_option = access_option;
	}

	virtual io_interface *create_io(thread *t);

	virtual void destroy_io(io_interface *io) {
	}
};

class aio_factory: public file_io_factory
{
public:
	aio_factory(const std::string &file_name): file_io_factory(file_name) {
	}

	virtual io_interface *create_io(thread *t);

	virtual void destroy_io(io_interface *io) {
	}
};

class remote_io_factory: public file_io_factory
{
protected:
	file_mapper *mapper;
public:
	remote_io_factory(const std::string &file_name);

	~remote_io_factory() {
		delete mapper;
	}

	virtual io_interface *create_io(thread *t);

	virtual void destroy_io(io_interface *io) {
	}
};

class global_cached_io_factory: public remote_io_factory
{
	const cache_config *cache_conf;
	static page_cache *global_cache;
public:
	global_cached_io_factory(const std::string &file_name,
			const cache_config *cache_conf): remote_io_factory(file_name) {
		this->cache_conf = cache_conf;
		if (global_cache == NULL) {
			global_cache = cache_conf->create_cache(MAX_NUM_FLUSHES_PER_FILE *
					global_data.raid_conf.get_num_disks());
			int num_files = global_data.read_threads.size();
			for (int k = 0; k < num_files; k++) {
				global_data.read_threads[k]->register_cache(global_cache);
			}
		}
	}

	~global_cached_io_factory() {
		cache_conf->destroy_cache(global_cache);
	}

	virtual io_interface *create_io(thread *t);

	virtual void destroy_io(io_interface *io) {
	}
};

class part_global_cached_io_factory: public remote_io_factory
{
	const cache_config *cache_conf;
	part_io_process_table *table;
	int num_nodes;
public:
	part_global_cached_io_factory(const std::string &file_name,
			const cache_config *cache_conf);

	~part_global_cached_io_factory() {
		part_global_cached_io::close_file(table);
	}

	virtual io_interface *create_io(thread *t);

	virtual void destroy_io(io_interface *io) {
	}
	virtual void print_state() {
		print_part_cached_io_state(table);
	}
};

io_interface *posix_io_factory::create_io(thread *t)
{
	file_mapper *mapper = global_data.raid_conf.create_file_mapper(get_name());
	int num_files = mapper->get_num_files();
	std::vector<int> indices;
	for (int i = 0; i < num_files; i++)
		indices.push_back(i);
	// The partition contains all files.
	logical_file_partition global_partition(indices, mapper);

	io_interface *io;
	switch (access_option) {
		case READ_ACCESS:
			io = new buffered_io(global_partition, t);
			break;
		case DIRECT_ACCESS:
			io = new direct_io(global_partition, t);
			break;
		default:
			fprintf(stderr, "a wrong posix access option\n");
			assert(0);
	}
#ifdef DEBUG
	global_data.register_io(io);
#endif
	return io;
}

io_interface *aio_factory::create_io(thread *t)
{
	file_mapper *mapper = global_data.raid_conf.create_file_mapper(get_name());
	int num_files = mapper->get_num_files();
	std::vector<int> indices;
	for (int i = 0; i < num_files; i++)
		indices.push_back(i);
	// The partition contains all files.
	logical_file_partition global_partition(indices, mapper);

	io_interface *io;
	io = new async_io(global_partition, params.get_aio_depth_per_file(), t);
#ifdef DEBUG
	global_data.register_io(io);
#endif
	return io;
}

remote_io_factory::remote_io_factory(const std::string &file_name): file_io_factory(
			file_name)
{
	mapper = global_data.raid_conf.create_file_mapper(get_name());
	int num_files = mapper->get_num_files();
	assert((int) global_data.read_threads.size() == num_files);

	for (int i = 0; i < num_files; i++) {
		global_data.read_threads[i]->open_file(mapper);
	}
}

io_interface *remote_io_factory::create_io(thread *t)
{
	io_interface *io = new remote_disk_access(global_data.read_threads,
			mapper, t);
#ifdef DEBUG
	global_data.register_io(io);
#endif
	return io;
}

io_interface *global_cached_io_factory::create_io(thread *t)
{
	io_interface *underlying = new remote_disk_access(
			global_data.read_threads, mapper, t);
	global_cached_io *io = new global_cached_io(t, underlying,
			global_cache);
#ifdef DEBUG
	global_data.register_io(io);
#endif
	return io;
}

part_global_cached_io_factory::part_global_cached_io_factory(
		const std::string &file_name,
		const cache_config *cache_conf): remote_io_factory(file_name)
{
	std::vector<int> node_id_array;
	cache_conf->get_node_ids(node_id_array);

	std::map<int, io_interface *> underlyings;
	table = part_global_cached_io::open_file(global_data.read_threads,
			mapper, cache_conf);
	this->cache_conf = cache_conf;
	this->num_nodes = node_id_array.size();

}

io_interface *part_global_cached_io_factory::create_io(thread *t)
{
	part_global_cached_io *io = part_global_cached_io::create(
			new remote_disk_access(global_data.read_threads, mapper, t), table);
#ifdef DEBUG
	global_data.register_io(io);
#endif
	return io;
}

file_io_factory *create_io_factory(const std::string &file_name,
		const int access_option, const cache_config *cache_conf)
{
	file_io_factory *factory;
	switch (access_option) {
		case READ_ACCESS:
		case DIRECT_ACCESS:
			factory = new posix_io_factory(file_name, access_option);
			break;
		case AIO_ACCESS:
			factory = new aio_factory(file_name);
			break;
		case REMOTE_ACCESS:
			factory = new remote_io_factory(file_name);
			break;
		case GLOBAL_CACHE_ACCESS:
			assert(cache_conf);
			factory = new global_cached_io_factory(file_name, cache_conf);
			break;
		case PART_GLOBAL_ACCESS:
			assert(cache_conf);
			factory = new part_global_cached_io_factory(file_name, cache_conf);
			break;
		default:
			fprintf(stderr, "a wrong access option\n");
			assert(0);
	}
	global_data.register_io_factory(factory);
	return factory;
}

void destroy_io_factory(file_io_factory *factory)
{
	delete factory;
}

void print_io_thread_stat()
{
	for (unsigned i = 0; i < global_data.read_threads.size(); i++) {
		disk_read_thread *t = global_data.read_threads[i];
		if (t)
			t->stop();
	}
	sleep(1);
	for (unsigned i = 0; i < global_data.read_threads.size(); i++) {
		disk_read_thread *t = global_data.read_threads[i];
		if (t)
			t->print_stat();
	}
}

atomic_integer io_interface::io_counter;
page_cache *global_cached_io_factory::global_cache;
