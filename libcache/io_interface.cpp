#include <assert.h>
#include <pthread.h>

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

class io_tracker
{
	const std::string file_name;
	io_interface *io;
	bool has_init;
	pthread_spinlock_t lock;
	bool taken;
public:
	io_tracker(const std::string &_file_name,
			io_interface *io): file_name(_file_name) {
		this->io = io;
		pthread_spin_init(&lock, PTHREAD_PROCESS_PRIVATE);
		taken = false;
		has_init = false;
	}

	const std::string &get_file_name() const {
		return file_name;
	}

	io_interface *take() {
		pthread_spin_lock(&lock);
		io_interface *ret = NULL;
		if (!taken)
			ret = io;
		taken = true;
		// If the IO is used for the first time, initialize it.
		if (!has_init && ret) {
			ret->init();
			has_init = true;
		}
		pthread_spin_unlock(&lock);
		return ret;
	}

	io_interface *get() {
//		assert(has_init && taken);
		return io;
	}

	void release() {
		pthread_spin_lock(&lock);
		taken = false;
		pthread_spin_unlock(&lock);
	}
};

class io_table
{
	pthread_spinlock_t lock;
	std::vector<io_tracker *> table;
public:
	io_table() {
		pthread_spin_init(&lock, PTHREAD_PROCESS_PRIVATE);
	}

	int size() {
		pthread_spin_lock(&lock);
		int ret = table.size();
		pthread_spin_unlock(&lock);
		return ret;
	}

	void register_io(const std::string &file_name, io_interface *io) {
		pthread_spin_lock(&lock);
		// If the IO table is too small, we need to resize the table.
		// TODO after a while, the IO index might be very large. I should use
		// a hashtable here.
		if ((int) table.size() <= io->get_io_idx())
			table.resize(io->get_io_idx() * 2);
		table[io->get_io_idx()] = new io_tracker(file_name, io);
		pthread_spin_unlock(&lock);
	}

	io_interface *get_io(int idx) {
		// TODO I need to make thread-safe here.
//		pthread_spin_lock(&lock);
		assert(idx >= 0);
		io_interface *io = table[idx]->get();
		assert(io);
		assert(io->get_io_idx() == idx);
//		pthread_spin_unlock(&lock);
		return io;
	}

	io_interface *allocate_io(const std::string &file_name, int node_id) {
		pthread_spin_lock(&lock);
		io_interface *io = NULL;
		for (unsigned i = 0; i < table.size(); i++) {
			io_interface *tmp = table[i]->get();
			if (tmp->get_node_id() == node_id
					&& table[i]->get_file_name() == file_name) {
				io = table[i]->take();
				if (io)
					break;
			}
		}
		pthread_spin_unlock(&lock);
		return io;
	}

	void release_io(io_interface *io) {
		pthread_spin_lock(&lock);
		table[io->get_io_idx()]->release();
		pthread_spin_unlock(&lock);
	}
};

static io_table ios;

void register_io(const std::string &file_name, io_interface *io)
{
	ios.register_io(file_name, io);
}

io_interface *get_io(int idx)
{
	return ios.get_io(idx);
}

io_interface *allocate_io(const std::string &file_name, int node_id)
{
	io_interface *ret = ios.allocate_io(file_name, node_id);
	assert(ret);
	return ret;
}

void release_io(io_interface *io)
{
	ios.release_io(io);
}

/**
 * This global data collection is very static.
 * Once the data is initialized, no data needs to be changed.
 * The mutex is to used only at the initialization.
 * As long as all threads call init_io_system() first before using
 * the global data, they will all see the complete global data.
 */
struct global_data_collection
{
	std::vector<disk_read_thread *> read_threads;
	std::tr1::unordered_map<int, aio_complete_thread *>  complete_threads;
	pthread_mutex_t mutex;

	global_data_collection() {
		pthread_mutex_init(&mutex, NULL);
	}
};

static global_data_collection global_data;

void init_io_system(const RAID_config &raid_conf,
		const std::vector<int> &node_id_array)
{
	numa_set_bind_policy(1);

	int num_nodes = node_id_array.size();
	file_mapper *mapper = raid_conf.create_file_mapper();
	int num_files = mapper->get_num_files();

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
		// Create threads for helping process completed AIO requests.
		for (int i = 0; i < num_nodes; i++) {
			int node_id = node_id_array[i];
			global_data.complete_threads.insert(
					std::pair<int, aio_complete_thread *>(node_id,
						new aio_complete_thread(node_id)));
		}
		for (int k = 0; k < num_files; k++) {
			std::vector<int> indices(1, k);
			logical_file_partition partition(indices);
			// Create disk accessing threads.
			global_data.read_threads[k] = new disk_read_thread(partition,
					global_data.complete_threads,
					raid_conf.get_file(k).node_id);
		}
	}
	pthread_mutex_unlock(&global_data.mutex);
	delete mapper;
}

void destroy_io_system()
{
	// TODO
}

class posix_io_factory: public file_io_factory
{
	int access_option;
	RAID_config raid_conf;
public:
	posix_io_factory(const RAID_config &_raid_conf,
			int access_option): file_io_factory(
				_raid_conf.get_conf_file()), raid_conf(_raid_conf) {
		this->access_option = access_option;
	}

	virtual io_interface *create_io(int node_id);

	virtual void destroy_io(io_interface *io) {
	}
};

class aio_factory: public file_io_factory
{
	int io_depth_per_file;
	RAID_config raid_conf;
public:
	aio_factory(const RAID_config &_raid_conf,
			int io_depth_per_file): file_io_factory(
				_raid_conf.get_conf_file()), raid_conf(_raid_conf) {
		this->io_depth_per_file = io_depth_per_file;
	}

	virtual io_interface *create_io(int node_id);

	virtual void destroy_io(io_interface *io) {
	}
};

class remote_io_factory: public file_io_factory
{
protected:
	file_mapper *mapper;
	RAID_config raid_conf;
public:
	remote_io_factory(const RAID_config &raid_conf,
			const std::vector<int> &node_id_array);

	~remote_io_factory() {
		delete mapper;
	}

	virtual io_interface *create_io(int node_id);

	virtual void destroy_io(io_interface *io) {
	}
};

class global_cached_io_factory: public remote_io_factory
{
	const cache_config *cache_conf;
	page_cache *global_cache;
public:
	global_cached_io_factory(const RAID_config &raid_conf,
			const std::vector<int> &node_id_array,
			const cache_config *cache_conf): remote_io_factory(
				raid_conf, node_id_array) {
		this->cache_conf = cache_conf;
		global_cache = cache_conf->create_cache();
	}

	~global_cached_io_factory() {
		cache_conf->destroy_cache(global_cache);
	}

	virtual io_interface *create_io(int node_id);

	virtual void destroy_io(io_interface *io) {
	}
};

class part_global_cached_io_factory: public remote_io_factory
{
	const cache_config *cache_conf;
	part_io_process_table *table;
	int num_nodes;
public:
	part_global_cached_io_factory(const RAID_config &raid_conf,
			const std::vector<int> &node_id_array,
			const cache_config *cache_conf);

	~part_global_cached_io_factory() {
		part_global_cached_io::close_file(table);
	}

	virtual io_interface *create_io(int node_id);

	virtual void destroy_io(io_interface *io) {
	}
};

io_interface *posix_io_factory::create_io(int node_id)
{
	file_mapper *mapper = raid_conf.create_file_mapper();
	int num_files = mapper->get_num_files();
	std::vector<int> indices;
	for (int i = 0; i < num_files; i++)
		indices.push_back(i);
	// The partition contains all files.
	logical_file_partition global_partition(indices, mapper);

	io_interface *io;
	switch (access_option) {
		case READ_ACCESS:
			io = new buffered_io(global_partition, node_id);
			break;
		case DIRECT_ACCESS:
			io = new direct_io(global_partition, node_id);
			break;
		default:
			fprintf(stderr, "a wrong posix access option\n");
			abort();
	}
	register_io(raid_conf.get_conf_file(), io);
	return io;
}

io_interface *aio_factory::create_io(int node_id)
{
	file_mapper *mapper = raid_conf.create_file_mapper();
	int num_files = mapper->get_num_files();
	std::vector<int> indices;
	for (int i = 0; i < num_files; i++)
		indices.push_back(i);
	// The partition contains all files.
	logical_file_partition global_partition(indices, mapper);

	io_interface *io;
	std::tr1::unordered_map<int, aio_complete_thread *> no_complete_threads;
	io = new async_io(global_partition, no_complete_threads,
			io_depth_per_file, node_id);
	register_io(raid_conf.get_conf_file(), io);
	return io;
}

remote_io_factory::remote_io_factory(const RAID_config &_raid_conf,
		const std::vector<int> &node_id_array): file_io_factory(
			_raid_conf.get_conf_file()), raid_conf(_raid_conf)
{
	init_io_system(raid_conf, node_id_array);

	mapper = raid_conf.create_file_mapper();
	int num_files = mapper->get_num_files();
	assert(global_data.read_threads.size() >= node_id_array.size());
	assert((int) global_data.read_threads.size() == num_files);

	for (int i = 0; i < num_files; i++) {
		global_data.read_threads[i]->open_file(mapper);
	}

	// We have to make sure the nodes where IOs are going to create
	// has I/O complete threads.
	for (unsigned i = 0; i < node_id_array.size(); i++) {
		int node_id = node_id_array[i];
		assert(global_data.complete_threads.find(node_id)
				!= global_data.complete_threads.end());
	}
}

io_interface *remote_io_factory::create_io(int node_id)
{
	io_interface *io = new remote_disk_access(global_data.read_threads,
			global_data.complete_threads[node_id], mapper, node_id);
	register_io(raid_conf.get_conf_file(), io);
	return io;
}

io_interface *global_cached_io_factory::create_io(int node_id)
{
	io_interface *underlying = new remote_disk_access(
			global_data.read_threads,
			global_data.complete_threads[node_id],
			mapper, node_id);
	global_cached_io *io = new global_cached_io(underlying,
			global_cache);
	register_io(raid_conf.get_conf_file(), io);
	return io;
}

part_global_cached_io_factory::part_global_cached_io_factory(
		const RAID_config &raid_conf, const std::vector<int> &node_id_array,
		const cache_config *cache_conf): remote_io_factory(raid_conf,
			node_id_array)
{
	std::map<int, io_interface *> underlyings;
	for (unsigned i = 0; i < node_id_array.size(); i++) {
		int node_id = node_id_array[i];
		underlyings.insert(std::pair<int, io_interface *>(node_id,
					new remote_disk_access(
						global_data.read_threads,
						global_data.complete_threads[node_id],
						mapper, node_id)));
	}
	table = part_global_cached_io::open_file(underlyings, cache_conf);
	this->cache_conf = cache_conf;
	this->num_nodes = node_id_array.size();

}

io_interface *part_global_cached_io_factory::create_io(int node_id)
{
	part_global_cached_io *io = part_global_cached_io::create(node_id, table);
	register_io(raid_conf.get_conf_file(), io);
	return io;
}

std::vector<io_interface *> create_ios(const RAID_config &raid_conf,
		cache_config *cache_conf, const std::vector<int> &node_id_array,
		int nthreads, int access_option, long size, bool preload)
{
	std::vector<io_interface *> ios;
	int io_depth = AIO_DEPTH_PER_FILE / nthreads;
	if (io_depth == 0)
		io_depth = 1;
	file_io_factory *factory = create_io_factory(raid_conf, node_id_array,
			access_option, io_depth, cache_conf);
	int num_nodes = node_id_array.size();
	int nthreads_per_node = nthreads / num_nodes;
	for (int i = 0, j = 0; i < num_nodes; i++) {
		int node_id = node_id_array[i];
		for (int k = 0; k < nthreads_per_node; k++, j++) {
			ios.push_back(factory->create_io(node_id));
		}
	}
	delete factory;

	return ios;
}

file_io_factory *create_io_factory(const RAID_config &raid_conf,
		const std::vector<int> &node_id_array, const int access_option,
		const int io_depth, const cache_config *cache_conf)
{
	switch (access_option) {
		case READ_ACCESS:
		case DIRECT_ACCESS:
			return new posix_io_factory(raid_conf, access_option);
		case AIO_ACCESS:
			assert(io_depth > 0);
			return new aio_factory(raid_conf, io_depth);
		case REMOTE_ACCESS:
			return new remote_io_factory(raid_conf, node_id_array);
		case GLOBAL_CACHE_ACCESS:
			assert(cache_conf);
			return new global_cached_io_factory(raid_conf, node_id_array,
					cache_conf);
		case PART_GLOBAL_ACCESS:
			assert(cache_conf);
			return new part_global_cached_io_factory(raid_conf,
					node_id_array, cache_conf);
		default:
			fprintf(stderr, "a wrong access option\n");
			abort();
	}
	return NULL;
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
			printf("queue on file %s wait for requests for %d times, is full for %d times, and %d accesses and %d io waits, complete %d reqs and %d low-prio reqs, ignore %d low-prio reqs\n",
					t->get_file_name().c_str(), t->get_queue()->get_num_empty(),
					t->get_queue()->get_num_full(), t->get_num_accesses(),
					t->get_num_iowait(), t->get_num_completed_reqs(),
					t->get_num_low_prio_accesses(), t->get_num_ignored_low_prio_accesses());
	}
}

atomic_integer io_interface::io_counter;
