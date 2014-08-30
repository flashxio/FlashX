/**
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

#include <assert.h>
#include <signal.h>
#include <time.h>
#include <pthread.h>

#include <unordered_set>
#include <unordered_map>
#include <vector>
#include <boost/foreach.hpp>

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
#include "debugger.h"
#include "mem_tracker.h"
#include "native_file.h"
#include "safs_file.h"

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
	std::vector<disk_io_thread *> read_threads;
	pthread_mutex_t mutex;
	cache_config *cache_conf;
	page_cache *global_cache;
#ifdef PART_IO
	// For part_global_cached_io
	part_io_process_table *table;
#endif

	global_data_collection() {
#ifdef PART_IO
		table = NULL;
#endif
		cache_conf = NULL;
		global_cache = NULL;
		pthread_mutex_init(&mutex, NULL);
	}
};

static global_data_collection global_data;

class file_mapper_set
{
	std::unordered_map<std::string, file_mapper *> map;
	spin_lock lock;
public:
	file_mapper &get(const std::string &name) {
		lock.lock();
		std::unordered_map<std::string, file_mapper *>::const_iterator it
			= map.find(name);
		file_mapper *mapper;
		if (it == map.end()) {
			mapper = global_data.raid_conf.create_file_mapper(name);
			map.insert(std::pair<std::string, file_mapper *>(name, mapper));
		}
		else
			mapper = it->second;
		lock.unlock();
		return *mapper;
	}
};
static file_mapper_set file_mappers;

class debug_global_data: public debug_task
{
public:
	void run();
};

void debug_global_data::run()
{
	for (unsigned i = 0; i < global_data.read_threads.size(); i++)
		global_data.read_threads[i]->print_state();
}

const RAID_config &get_sys_RAID_conf()
{
	return global_data.raid_conf;
}

static std::vector<int> file_weights;

void set_file_weight(const std::string &file_name, int weight)
{
	file_mapper &mapper = file_mappers.get(file_name);
	if ((size_t) mapper.get_file_id() >= file_weights.size())
		file_weights.resize(mapper.get_file_id() + 1);
	file_weights[mapper.get_file_id()] = weight;
	printf("%s: id: %d, weight: %d\n", file_name.c_str(),
			mapper.get_file_id(), weight);
}

void parse_file_weights(const std::string &str)
{
	std::vector<std::string> file_strs;
	split_string(str, ',', file_strs);
	if (file_strs.size() > 0)
		file_weights.resize(file_strs.size());
	BOOST_FOREACH(std::string s, file_strs) {
		std::vector<std::string> ss;
		split_string(s, ':', ss);
		if (ss.size() != 2) {
			printf("file weight in wrong format: %s\n", s.c_str());
			continue;
		}

		int weight = atoi(ss[1].c_str());
		set_file_weight(ss[0], weight);
	}
	// When we resize the vector, the files that are not assigned a weight
	// get 0 as weight. We need to make their weight 1.
	for (size_t i = 0; i < file_weights.size(); i++)
		if (file_weights[i] == 0)
			file_weights[i] = 1;
}

/*
 * This method returns user-defined weight for a SAFS file in
 * the configuration. If the weight isn't defined, return 1.
 * By using this mechanism, users can tweak the performance of
 * the page cache while dealing with multiple files.
 */
int get_file_weight(file_id_t file_id)
{
	if ((size_t) file_id < file_weights.size())
		return file_weights[file_id];
	else
		return 1;
}

void init_io_system(const config_map &configs, bool with_cache)
{
#ifdef ENABLE_MEM_TRACE
	init_mem_tracker();
#endif
	
	params.init(configs.get_options());
	params.print();

	numa_set_bind_policy(1);
	thread::thread_class_init();

	std::string root_conf_file = configs.get_option("root_conf");
	printf("The root conf file: %s\n", root_conf_file.c_str());
	RAID_config raid_conf(root_conf_file, params.get_RAID_mapping_option(),
			params.get_RAID_block_size());
	int num_files = raid_conf.get_num_disks();
	global_data.raid_conf = raid_conf;

	std::set<int> disk_node_set = raid_conf.get_node_ids();
	std::vector<int> disk_node_ids(disk_node_set.begin(), disk_node_set.end());
	printf("There are %ld nodes with disks\n", disk_node_ids.size());
	init_aio(disk_node_ids);

	file_mapper *mapper = raid_conf.create_file_mapper();
	if (configs.has_option("file_weights"))
		parse_file_weights(configs.get_option("file_weights"));
	/* 
	 * The mutex is enough to guarantee that all threads will see initialized
	 * global data. The first thread that enters the critical area will
	 * initialize the global data. If another thread tries to run the code,
	 * it will be blocked by the mutex. When a thread is returned from
	 * the function, they all can see the global data.
	 */
	pthread_mutex_lock(&global_data.mutex);
	int flags = O_RDONLY;
	if (params.is_writable())
		flags = O_RDWR;
	// The global data hasn't been initialized.
	if (global_data.read_threads.size() == 0) {
		global_data.read_threads.resize(num_files);
		for (int k = 0; k < num_files; k++) {
			std::vector<int> indices(1, k);
			logical_file_partition partition(indices, mapper);
			// Create disk accessing threads.
			global_data.read_threads[k] = new disk_io_thread(partition,
					global_data.raid_conf.get_disk(k).node_id, NULL, k, flags);
		}
		debug.register_task(new debug_global_data());
	}

	// Assign a thread object to the current thread.
	if (thread::get_curr_thread() == NULL)
		thread::represent_thread(0);

	if (global_data.global_cache == NULL && with_cache) {
		std::vector<int> node_id_array;
		for (int i = 0; i < params.get_num_nodes(); i++)
			node_id_array.push_back(i);

		global_data.cache_conf = new even_cache_config(params.get_cache_size(),
				params.get_cache_type(), node_id_array);
		global_data.global_cache = global_data.cache_conf->create_cache(
				MAX_NUM_FLUSHES_PER_FILE *
				global_data.raid_conf.get_num_disks());
		int num_files = global_data.read_threads.size();
		for (int k = 0; k < num_files; k++) {
			global_data.read_threads[k]->register_cache(
					global_data.global_cache);
		}

		// The remote IO will never be used. It's only used for creating
		// more remote IOs for flushing dirty pages, so it doesn't matter
		// what thread is used here.
		thread *curr = thread::get_curr_thread();
		assert(curr);
		io_interface::ptr underlying = io_interface::ptr(new remote_io(
					global_data.read_threads, mapper, curr));
		global_data.global_cache->init(underlying);
	}
#ifdef PART_IO
	if (global_data.table == NULL && with_cache) {
		if (params.get_num_nodes() > 1)
			global_data.table = part_global_cached_io::init_subsystem(
					global_data.read_threads, mapper,
					(NUMA_cache *) global_data.global_cache);
	}
#endif
	pthread_mutex_unlock(&global_data.mutex);
}

void destroy_io_system()
{
	if (global_data.global_cache)
		global_data.global_cache->sanity_check();
#ifdef PART_IO
	// TODO destroy part global cached io table.
	if (global_data.table) {
		part_global_cached_io::destroy_subsystem(global_data.table);
		global_data.table = NULL;
	}
#endif
	if (global_data.cache_conf) {
		global_data.cache_conf->destroy_cache(global_data.global_cache);
		global_data.global_cache = NULL;
		delete global_data.cache_conf;
		global_data.cache_conf = NULL;
	}
	BOOST_FOREACH(disk_io_thread *t, global_data.read_threads) {
		t->stop();
		t->join();
		delete t;
	}
	global_data.read_threads.resize(0);
	destroy_aio();

#ifdef ENABLE_MEM_TRACE
	printf("memleak: %ld objects and %ld bytes\n", get_alloc_objs(),
			get_alloc_bytes());
	printf("max: %ld objs and %ld bytes, max alloc %ld bytes\n",
			get_max_alloc_objs(), get_max_alloc_bytes(), get_max_alloc());
#endif
}

class posix_io_factory: public file_io_factory
{
	int access_option;
	// The number of existing IO instances.
	std::atomic<size_t> num_ios;
	file_mapper &mapper;
public:
	posix_io_factory(file_mapper &_mapper, int access_option): file_io_factory(
				_mapper.get_name()), mapper(_mapper) {
		this->access_option = access_option;
		num_ios = 0;
	}

	~posix_io_factory() {
		assert(num_ios == 0);
	}

	virtual io_interface::ptr create_io(thread *t);

	virtual void destroy_io(io_interface *io);

	virtual int get_file_id() const {
		assert(0);
		return -1;
	}
};

class aio_factory: public file_io_factory
{
	// The number of existing IO instances.
	std::atomic<size_t> num_ios;
	file_mapper &mapper;
public:
	aio_factory(file_mapper &_mapper): file_io_factory(
			_mapper.get_name()), mapper(_mapper) {
		num_ios = 0;
	}

	~aio_factory() {
		assert(num_ios == 0);
	}

	virtual io_interface::ptr create_io(thread *t);

	virtual void destroy_io(io_interface *io);

	virtual int get_file_id() const {
		assert(0);
		return -1;
	}
};

class remote_io_factory: public file_io_factory
{
	std::atomic_ulong tot_accesses;
protected:
	// The number of existing IO instances.
	std::atomic<size_t> num_ios;
	file_mapper &mapper;
public:
	remote_io_factory(file_mapper &_mapper);

	~remote_io_factory();

	virtual io_interface::ptr create_io(thread *t);

	virtual void destroy_io(io_interface *io);

	virtual int get_file_id() const {
		return mapper.get_file_id();
	}

	virtual void collect_stat(io_interface &io) {
		remote_io &rio = (remote_io &) io;
		tot_accesses += rio.get_num_reqs();
	}

	virtual void print_statistics() const {
		printf("%s gets %ld I/O accesses\n", mapper.get_name().c_str(),
				tot_accesses.load());
	}
};

class global_cached_io_factory: public remote_io_factory
{
	std::atomic_ulong tot_bytes;
	std::atomic_ulong tot_accesses;
	std::atomic_ulong tot_hits;
	std::atomic_ulong tot_fast_process;

	page_cache *global_cache;
public:
	global_cached_io_factory(file_mapper &_mapper,
			page_cache *cache): remote_io_factory(_mapper) {
		this->global_cache = cache;
		tot_bytes = 0;
		tot_accesses = 0;
		tot_hits = 0;
		tot_fast_process = 0;
	}

	virtual io_interface::ptr create_io(thread *t);

	virtual void destroy_io(io_interface *io);

	virtual void collect_stat(io_interface &io) {
		global_cached_io &gio = (global_cached_io &) io;

		tot_bytes += gio.get_num_bytes();
		tot_accesses += gio.get_num_pg_accesses();
		tot_hits += gio.get_cache_hits();
		tot_fast_process += gio.get_num_fast_process();
	}

	virtual void print_statistics() const {
		printf("%s gets %ld I/O accesses, %ld in bytes\n", mapper.get_name().c_str(),
				tot_accesses.load(), tot_bytes.load());
		printf("There are %ld cache hits, %ld of them are in the fast process\n",
				tot_hits.load(), tot_fast_process.load());
	}
};

#ifdef PART_IO
class part_global_cached_io_factory: public remote_io_factory
{
public:
	part_global_cached_io_factory(
			file_mapper &_mapper): remote_io_factory(_mapper) {
	}

	virtual io_interface::ptr create_io(thread *t);

	virtual void destroy_io(io_interface *io);
};
#endif

class io_deleter
{
	file_io_factory &factory;
public:
	io_deleter(file_io_factory &_factory): factory(_factory) {
	}

	void operator()(io_interface *io) {
		factory.collect_stat(*io);
		factory.destroy_io(io);
	}
};

io_interface::ptr posix_io_factory::create_io(thread *t)
{
	int num_files = mapper.get_num_files();
	std::vector<int> indices;
	for (int i = 0; i < num_files; i++)
		indices.push_back(i);
	// The partition contains all files.
	logical_file_partition global_partition(indices, &mapper);

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
	num_ios++;
	return io_interface::ptr(io, io_deleter(*this));
}

void posix_io_factory::destroy_io(io_interface *io)
{
	num_ios--;
	delete io;
}

io_interface::ptr aio_factory::create_io(thread *t)
{
	int num_files = mapper.get_num_files();
	std::vector<int> indices;
	for (int i = 0; i < num_files; i++)
		indices.push_back(i);
	// The partition contains all files.
	logical_file_partition global_partition(indices, &mapper);

	io_interface *io;
	io = new async_io(global_partition, params.get_aio_depth_per_file(),
			t, O_RDWR);
	num_ios++;
	return io_interface::ptr(io, io_deleter(*this));
}

void aio_factory::destroy_io(io_interface *io)
{
	num_ios--;
	delete io;
}

remote_io_factory::remote_io_factory(file_mapper &_mapper): file_io_factory(
			_mapper.get_name()), mapper(_mapper)
{
	tot_accesses = 0;
	num_ios = 0;
	int num_files = mapper.get_num_files();
	assert((int) global_data.read_threads.size() == num_files);

	for (int i = 0; i < num_files; i++) {
		global_data.read_threads[i]->open_file(&mapper);
	}
}

remote_io_factory::~remote_io_factory()
{
	assert(num_ios == 0);
	int num_files = mapper.get_num_files();
	for (int i = 0; i < num_files; i++)
		global_data.read_threads[i]->close_file(&mapper);
}

io_interface::ptr remote_io_factory::create_io(thread *t)
{
	num_ios++;
	io_interface *io = new remote_io(global_data.read_threads, &mapper, t);
	return io_interface::ptr(io, io_deleter(*this));
}

void remote_io_factory::destroy_io(io_interface *io)
{
	num_ios--;
	delete io;
}

io_interface::ptr global_cached_io_factory::create_io(thread *t)
{
	io_interface *underlying = new remote_io(
			global_data.read_threads, &mapper, t);
	comp_io_scheduler *scheduler = NULL;
	if (get_sched_creater())
		scheduler = get_sched_creater()->create(underlying->get_node_id());
	global_cached_io *io = new global_cached_io(t, underlying,
			global_cache, scheduler);
	num_ios++;
	return io_interface::ptr(io, io_deleter(*this));
}

void global_cached_io_factory::destroy_io(io_interface *io)
{
	num_ios--;
	// The underlying IO is deleted in global_cached_io's destructor.
	delete io;
}


#ifdef PART_IO
io_interface::ptr part_global_cached_io_factory::create_io(thread *t)
{
	part_global_cached_io *io = part_global_cached_io::create(
			new remote_io(global_data.read_threads, &mapper, t),
			global_data.table);
	num_ios++;
	return io_interface::ptr(io, io_deleter(*this));
}

void part_global_cached_io_factory::destroy_io(io_interface *io)
{
	num_ios--;
	part_global_cached_io::destroy((part_global_cached_io *) io);
}
#endif

class destroy_io_factory
{
public:
	void operator()(file_io_factory *factory) {
		delete factory;
	}
};

file_io_factory::shared_ptr create_io_factory(const std::string &file_name,
		const int access_option)
{
	for (int i = 0; i < global_data.raid_conf.get_num_disks(); i++) {
		std::string abs_path = global_data.raid_conf.get_disk(i).name
			+ "/" + file_name;
		native_file f(abs_path);
		if (!f.exist()) {
			fprintf(stderr, "the underlying file %s doesn't exist\n",
					abs_path.c_str());
			return file_io_factory::shared_ptr();
		}
	}

	file_mapper &mapper = file_mappers.get(file_name);
	file_io_factory *factory = NULL;
	switch (access_option) {
		case READ_ACCESS:
		case DIRECT_ACCESS:
			factory = new posix_io_factory(mapper, access_option);
			break;
		case AIO_ACCESS:
			factory = new aio_factory(mapper);
			break;
		case REMOTE_ACCESS:
			factory = new remote_io_factory(mapper);
			break;
		case GLOBAL_CACHE_ACCESS:
			if (global_data.global_cache)
				factory = new global_cached_io_factory(mapper,
						global_data.global_cache);
			break;
#ifdef PART_IO
		case PART_GLOBAL_ACCESS:
			if (global_data.global_cache)
				factory = new part_global_cached_io_factory(file_name);
			break;
#endif
		default:
			fprintf(stderr, "a wrong access option\n");
			assert(0);
	}
	if (factory)
		return file_io_factory::shared_ptr(factory, destroy_io_factory());
	else
		return file_io_factory::shared_ptr();
}

void print_io_thread_stat()
{
	sleep(1);
	for (unsigned i = 0; i < global_data.read_threads.size(); i++) {
		disk_io_thread *t = global_data.read_threads[i];
		if (t)
			t->print_stat();
	}
}

ssize_t file_io_factory::get_file_size() const
{
	safs_file f(global_data.raid_conf, name);
	return f.get_file_size();
}

atomic_integer io_interface::io_counter;
