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

#include <assert.h>
#include <signal.h>
#include <time.h>
#include <pthread.h>

#include <unordered_set>
#include <unordered_map>
#include <vector>
#include <boost/foreach.hpp>
#include <boost/format.hpp>

#include "log.h"
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
#include "safs_exception.h"
#include "direct_comp_access.h"

namespace safs
{

/*
 * This global data collection is very static.
 * Once the data is initialized, no data needs to be changed.
 * The mutex is to used only at the initialization.
 * As long as all threads call init_io_system() first before using
 * the global data, they will all see the complete global data.
 */
struct global_data_collection
{
	// Count the number of times init_io_system is executed successfully.
	std::atomic<long> init_count;
	RAID_config::ptr raid_conf;
	// This contains a set of unique I/O threads.
	std::set<disk_io_thread::ptr> read_thread_set;
	// Each element points to an I/O thread in the set above.
	// Multiple elements may point to the same I/O thread.
	std::vector<disk_io_thread::ptr> read_threads;
	pthread_mutex_t mutex;
	// TODO there is memory leak here.
	cache_config::ptr cache_conf;
	page_cache::ptr global_cache;
	std::vector<int> io_cpus;
#ifdef PART_IO
	// For part_global_cached_io
	part_io_process_table *table;
#endif

	global_data_collection() {
#ifdef PART_IO
		table = NULL;
#endif
		pthread_mutex_init(&mutex, NULL);
	}
};

static global_data_collection global_data;

class debug_global_data: public debug_task
{
public:
	void run();
};

void debug_global_data::run()
{
	for (auto it = global_data.read_thread_set.begin();
			it != global_data.read_thread_set.end(); it++)
		(*it)->print_state();
}

const RAID_config &get_sys_RAID_conf()
{
	assert(global_data.raid_conf);
	return *global_data.raid_conf;
}

void init_io_system(config_map::ptr configs, bool with_cache)
{
#ifdef ENABLE_MEM_TRACE
	init_mem_tracker();
#endif
	long count = global_data.init_count.fetch_add(1);
	// If the I/O system has been initialized before, even if the previous
	// init failed, we still require users to destroy the I/O system before
	// they can initialize it again.
	if (count > 0)
		return;

	// We should initialize threads even if there aren't configs.
	thread::thread_class_init();
	// Assign a thread object to the main thread.
	if (thread::get_curr_thread() == NULL)
		thread::represent_thread(-1);

	if (configs == NULL)
		throw init_error("config map doesn't contain any options");
	
	params.init(configs->get_options());

	// The I/O system has been initialized.
	if (is_safs_init()) {
		assert(!global_data.read_thread_set.empty());
		return;
	}

	// If we don't have libaio, we disable the initialization of SAFS.
#ifdef USE_LIBAIO
	if (!configs->has_option("root_conf"))
		throw init_error("RAID config file doesn't exist");
	std::string root_conf_file = configs->get_option("root_conf");
	BOOST_LOG_TRIVIAL(info) << "The root conf file: " << root_conf_file;
	RAID_config::ptr raid_conf = RAID_config::create(root_conf_file,
			params.get_RAID_mapping_option(), params.get_RAID_block_size());
	// If we can't initialize RAID, there is nothing we can do.
	if (raid_conf == NULL) {
		throw init_error("can't create RAID config");
	}

	int num_files = raid_conf->get_num_disks();
	global_data.raid_conf = raid_conf;

	std::set<int> disk_node_set = raid_conf->get_node_ids();
	std::vector<int> disk_node_ids(disk_node_set.begin(), disk_node_set.end());
	BOOST_LOG_TRIVIAL(info) << boost::format("There are %1% nodes with disks")
		% disk_node_ids.size();
	init_aio(disk_node_ids);

	file_mapper::ptr mapper = raid_conf->create_file_mapper();
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
		// Determine a map to indicate which NUMA nodes the disks are
		// attached to.
		std::map<int, std::vector<int> > indices;
		for (int i = 0; i < num_files; i++) {
			int node_id = mapper->get_file_node_id(i);
			auto it = indices.find(node_id);
			if (it == indices.end()) {
				std::vector<int> per_node_indices(1);
				per_node_indices[0] = i;
				indices.insert(std::pair<int, std::vector<int> >(node_id,
							per_node_indices));
			}
			else {
				it->second.push_back(i);
			}
		}
		// Iterate over the NUMA nodes with disks.
		size_t tot_num_threads = 0;
		for (auto it = indices.begin(); it != indices.end(); it++) {
			size_t num_io_threads = params.get_num_io_threads();
			// If the number of I/O threads isn't specified, we give an I/O
			// thread for each SSD.
			if (num_io_threads == 0)
				num_io_threads = it->second.size();
			// Create disk accessing threads.
			assert(it->second.size() % num_io_threads == 0);

			std::vector<disk_io_thread::ptr> ts(num_io_threads);
			tot_num_threads += num_io_threads;
			for (size_t i = 0; i < ts.size(); i++) {
				std::vector<int> disks(it->second.size() / ts.size());
				for (size_t j = 0; j < disks.size(); j++)
					disks[j] = it->second[i * disks.size() + j];
				logical_file_partition partition(disks, mapper);

#ifdef USE_HWLOC
				if (params.is_bind_io_thread()) {
					// If we bind an I/O thread to a specific CPU core, the CPU
					// core will be used by the thread exclusively.
					const CPU_core &core = cpus.get_node(it->first).get_core(i);
					std::vector<int> units = core.get_units();
					global_data.io_cpus.insert(global_data.io_cpus.end(),
							units.begin(), units.end());
					ts[i] = disk_io_thread::ptr(new disk_io_thread(partition,
								units[0], it->first, flags));
				}
				else
#endif
					ts[i] = disk_io_thread::ptr(new disk_io_thread(partition,
								it->first, flags));
				for (size_t j = 0; j < disks.size(); j++) {
					int file_idx = disks[j];
					global_data.read_threads[file_idx] = ts[i];
				}
			}
		}
		BOOST_LOG_TRIVIAL(info) << boost::format(
				"SAFS runs on %1% SSDs with %2% I/O threads") % num_files
			% tot_num_threads;
		global_data.read_thread_set.insert(global_data.read_threads.begin(),
				global_data.read_threads.end());
#if 0
		debug.register_task(new debug_global_data());
#endif
	}

	if (global_data.global_cache == NULL && with_cache
			&& params.get_cache_size() > 0) {
		std::vector<int> node_id_array;
		for (int i = 0; i < params.get_num_nodes(); i++)
			node_id_array.push_back(i);

		global_data.cache_conf = cache_config::ptr(new even_cache_config(
					params.get_cache_size(), params.get_cache_type(),
					node_id_array));
		global_data.global_cache = global_data.cache_conf->create_cache(
				MAX_NUM_FLUSHES_PER_FILE *
				global_data.raid_conf->get_num_disks());

		// The remote IO will never be used. It's only used for creating
		// more remote IOs for flushing dirty pages, so it doesn't matter
		// what thread is used here.
#if 0
		thread *curr = thread::get_curr_thread();
		assert(curr);
		io_interface::ptr underlying = io_interface::ptr(new remote_io(
					global_data.read_threads, mapper, curr));
		global_data.global_cache->init(underlying);
#endif
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
#else
	throw init_error("There isn't libaio. SAFS isn't initialized.");
#endif
}

void destroy_io_system()
{
	long count = global_data.init_count.fetch_sub(1);
	assert(count > 0);
	// We only destroy the I/O system when destroy_io_system is invoked
	// the same number of times that init_io_system is invoked.
	if (count > 1) {
		return;
	}

	BOOST_LOG_TRIVIAL(info) << "I/O system is destroyed";
	global_data.raid_conf.reset();
	if (global_data.global_cache)
		global_data.global_cache->sanity_check();
#ifdef PART_IO
	// TODO destroy part global cached io table.
	if (global_data.table) {
		part_global_cached_io::destroy_subsystem(global_data.table);
		global_data.table = NULL;
	}
#endif
	size_t num_reads = 0;
	size_t num_writes = 0;
	size_t num_read_bytes = 0;
	size_t num_write_bytes = 0;
	BOOST_FOREACH(disk_io_thread::ptr t, global_data.read_thread_set) {
		if (t) {
			num_reads += t->get_num_reads();
			num_writes += t->get_num_writes();
			num_read_bytes += t->get_num_read_bytes();
			num_write_bytes += t->get_num_write_bytes();
		}
	}
	global_data.read_threads.clear();
	global_data.read_thread_set.clear();
	destroy_aio();
	BOOST_LOG_TRIVIAL(info)
		<< boost::format("I/O threads get %1% reads (%2% bytes) and %3% writes (%4% bytes)")
		% num_reads % num_read_bytes % num_writes % num_write_bytes;

#ifdef ENABLE_MEM_TRACE
	BOOST_LOG_TRIVIAL(info) << boost::format("memleak: %1% objects and %2% bytes")
		% get_alloc_objs() % get_alloc_bytes();
	BOOST_LOG_TRIVIAL(info)
		<< boost::format("max: %1% objs and %2% bytes, max alloc %3% bytes")
		% get_max_alloc_objs() % get_max_alloc_bytes() % get_max_alloc();
#endif
}

class posix_io_factory: public file_io_factory
{
	int access_option;
	// The number of existing IO instances.
	std::atomic<size_t> num_ios;
	file_mapper::ptr mapper;
public:
	posix_io_factory(file_mapper::ptr mapper, int access_option): file_io_factory(
				mapper->get_name()) {
		this->mapper = mapper;
		this->access_option = access_option;
		num_ios = 0;
	}

	~posix_io_factory() {
		assert(num_ios == 0);
	}

	virtual io_interface::ptr create_io(thread *t);

	virtual void destroy_io(io_interface &io);

	virtual int get_file_id() const {
		throw unsupported_exception("get_file_id");
	}
};

class aio_factory: public file_io_factory
{
	// The number of existing IO instances.
	std::atomic<size_t> num_ios;
	file_mapper::ptr mapper;
public:
	aio_factory(file_mapper::ptr mapper): file_io_factory(
			mapper->get_name()) {
		this->mapper = mapper;
		num_ios = 0;
	}

	~aio_factory() {
		assert(num_ios == 0);
	}

	virtual io_interface::ptr create_io(thread *t);

	virtual void destroy_io(io_interface &io);

	virtual int get_file_id() const {
		throw unsupported_exception("get_file_id");
	}
};

class remote_io_factory: public file_io_factory
{
	std::shared_ptr<slab_allocator> unbind_msg_allocator;
	std::vector<std::shared_ptr<slab_allocator> > msg_allocators;
	std::atomic_ulong tot_accesses;
	// The number of existing IO instances.
	std::atomic<size_t> num_ios;
	file_mapper::ptr mapper;

	slab_allocator &get_msg_allocator(int node_id) {
		if (node_id < 0)
			return *unbind_msg_allocator;
		else
			return *msg_allocators[node_id];
	}
public:
	remote_io_factory(file_mapper::ptr mapper);

	~remote_io_factory();

	virtual io_interface::ptr create_io(thread *t);

	virtual void destroy_io(io_interface &io);

	virtual int get_file_id() const {
		return mapper->get_file_id();
	}

	virtual void collect_stat(io_interface &io) {
		remote_io &rio = (remote_io &) io;
		tot_accesses += rio.get_num_reqs();
	}

	virtual void print_statistics() const {
		BOOST_LOG_TRIVIAL(info) << boost::format("%1% gets %2% I/O accesses")
			% mapper->get_name() % tot_accesses.load();
	}
};

class global_cached_io_factory: public file_io_factory
{
	std::atomic_ulong tot_bytes;
	std::atomic_ulong tot_accesses;
	std::atomic_ulong tot_pg_accesses;
	std::atomic_ulong tot_hits;
	std::atomic_ulong tot_fast_process;

	page_cache::ptr global_cache;
	remote_io_factory::shared_ptr remote_factory;
public:
	global_cached_io_factory(file_mapper::ptr mapper,
			page_cache::ptr cache): file_io_factory(mapper->get_name()) {
		this->global_cache = cache;
		tot_bytes = 0;
		tot_accesses = 0;
		tot_pg_accesses = 0;
		tot_hits = 0;
		tot_fast_process = 0;
		remote_factory = remote_io_factory::shared_ptr(new remote_io_factory(mapper));
	}

	virtual io_interface::ptr create_io(thread *t);

	virtual void destroy_io(io_interface &io);

	virtual int get_file_id() const {
		return remote_factory->get_file_id();
	}

	virtual void collect_stat(io_interface &io) {
		global_cached_io &gio = (global_cached_io &) io;

		tot_bytes += gio.get_num_bytes();
		tot_accesses += gio.get_num_areqs();
		tot_pg_accesses += gio.get_num_pg_accesses();
		tot_hits += gio.get_cache_hits();
		tot_fast_process += gio.get_num_fast_process();
	}

	virtual void print_statistics() const {
		BOOST_LOG_TRIVIAL(info)
			<< boost::format("%1% gets %2% async I/O accesses, %3% in bytes")
			% get_name() % tot_accesses.load() % tot_bytes.load();
		BOOST_LOG_TRIVIAL(info)
			<< boost::format("There are %1% pages accessed, %2% cache hits, %3% of them are in the fast process")
			% tot_pg_accesses.load() % tot_hits.load() % tot_fast_process.load();
	}
};

class direct_comp_io_factory: public file_io_factory
{
	// The number of bytes accessed from the disks.
	std::atomic_ulong tot_disk_bytes;
	// The number of bytes requested by applications.
	std::atomic_ulong tot_req_bytes;
	std::atomic_ulong tot_accesses;

	remote_io_factory::shared_ptr remote_factory;
public:
	direct_comp_io_factory(file_mapper::ptr mapper): file_io_factory(
				mapper->get_name()) {
		tot_disk_bytes = 0;
		tot_req_bytes = 0;
		tot_accesses = 0;
		remote_factory = remote_io_factory::shared_ptr(new remote_io_factory(mapper));
	}

	virtual io_interface::ptr create_io(thread *t);

	virtual void destroy_io(io_interface &io);

	virtual int get_file_id() const {
		return remote_factory->get_file_id();
	}

	virtual void collect_stat(io_interface &io) {
		direct_comp_io &dio = (direct_comp_io &) io;

		tot_disk_bytes += dio.get_num_disk_bytes();
		tot_req_bytes += dio.get_num_req_bytes();
		tot_accesses += dio.get_num_reqs();
	}

	virtual void print_statistics() const {
		BOOST_LOG_TRIVIAL(info)
			<< boost::format("%1% gets %2% async I/O accesses, %3% req bytes, %4% disk bytes")
			% get_name() % tot_accesses.load() % tot_req_bytes.load()
			% tot_disk_bytes.load();
	}
};

#ifdef PART_IO
class part_global_cached_io_factory: public remote_io_factory
{
public:
	part_global_cached_io_factory(
			file_mapper::ptr mapper): remote_io_factory(mapper) {
	}

	virtual io_interface::ptr create_io(thread *t);

	virtual void destroy_io(io_interface &io);
};
#endif

io_interface::ptr posix_io_factory::create_io(thread *t)
{
	int num_files = mapper->get_num_files();
	std::vector<int> indices;
	for (int i = 0; i < num_files; i++)
		indices.push_back(i);
	// The partition contains all files.
	logical_file_partition global_partition(indices, mapper);

	io_interface *io;
	switch (access_option) {
		case READ_ACCESS:
			io = new buffered_io(global_partition, t, get_header());
			break;
		case DIRECT_ACCESS:
			io = new direct_io(global_partition, t, get_header());
			break;
		default:
			fprintf(stderr, "a wrong posix access option\n");
			return io_interface::ptr();
	}
	num_ios++;
	return io_interface::ptr(io);
}

void posix_io_factory::destroy_io(io_interface &io)
{
	num_ios--;
}

io_interface::ptr aio_factory::create_io(thread *t)
{
	int num_files = mapper->get_num_files();
	std::vector<int> indices;
	for (int i = 0; i < num_files; i++)
		indices.push_back(i);
	// The partition contains all files.
	logical_file_partition global_partition(indices, mapper);

	io_interface *io;
	io = new async_io(global_partition, params.get_aio_depth_per_file(),
			t, get_header(), O_RDWR);
	num_ios++;
	return io_interface::ptr(io);
}

void aio_factory::destroy_io(io_interface &io)
{
	num_ios--;
}

remote_io_factory::remote_io_factory(file_mapper::ptr mapper): file_io_factory(
			mapper->get_name())
{
	this->mapper = mapper;
	msg_allocators.resize(params.get_num_nodes());
	for (int i = 0; i < params.get_num_nodes(); i++) {
		msg_allocators[i] = std::shared_ptr<slab_allocator>(new slab_allocator(
					std::string("disk_msg_allocator-") + itoa(i),
					IO_MSG_SIZE * sizeof(io_request),
					IO_MSG_SIZE * sizeof(io_request) * 1024, INT_MAX, i));
	}
	unbind_msg_allocator = std::shared_ptr<slab_allocator>(new slab_allocator(
				std::string("disk_msg_allocator-unbind"),
				IO_MSG_SIZE * sizeof(io_request),
				IO_MSG_SIZE * sizeof(io_request) * 1024, INT_MAX, -1));
	tot_accesses = 0;
	num_ios = 0;
	int num_files = mapper->get_num_files();
	assert((int) global_data.read_threads.size() == num_files);

	for (auto it = global_data.read_thread_set.begin();
			it != global_data.read_thread_set.end(); it++)
		(*it)->open_file(mapper);
}

remote_io_factory::~remote_io_factory()
{
	assert(num_ios == 0);
	// If the I/O threads haven't been destroyed.
	for (auto it = global_data.read_thread_set.begin();
			it != global_data.read_thread_set.end(); it++)
		(*it)->close_file(mapper);
}

io_interface::ptr remote_io_factory::create_io(thread *t)
{
	if (t->get_node_id() >= (int) msg_allocators.size()) {
		fprintf(stderr, "thread %d are not in a right node (%d)\n",
				t->get_id(), t->get_node_id());
		return io_interface::ptr();
	}

	num_ios++;
	io_interface *io = new remote_io(global_data.read_threads,
			get_msg_allocator(t->get_node_id()), mapper, t, get_header());
	return io_interface::ptr(io);
}

void remote_io_factory::destroy_io(io_interface &io)
{
	num_ios--;
}

io_interface::ptr global_cached_io_factory::create_io(thread *t)
{
	io_interface::ptr underlying = safs::create_io(remote_factory, t);
	comp_io_scheduler::ptr scheduler;
	if (get_sched_creator())
		scheduler = get_sched_creator()->create(underlying->get_node_id());
	global_cached_io *io = new global_cached_io(t, underlying,
			global_cache, scheduler);
	return io_interface::ptr(io);
}

void global_cached_io_factory::destroy_io(io_interface &io)
{
	// num_ios is decreased by the underlying remote I/O instance.

	// The underlying IO is deleted in global_cached_io's destructor.
}

io_interface::ptr direct_comp_io_factory::create_io(thread *t)
{
	io_interface::ptr underlying = safs::create_io(remote_factory, t);
	direct_comp_io *io = new direct_comp_io(
			std::static_pointer_cast<remote_io>(underlying));
	return io_interface::ptr(io);
}

void direct_comp_io_factory::destroy_io(io_interface &io)
{
	// num_ios is decreased by the underlying remote I/O instance.

	// The underlying IO is deleted in global_cached_io's destructor.
}

#ifdef PART_IO
io_interface::ptr part_global_cached_io_factory::create_io(thread *t)
{
	part_global_cached_io *io = part_global_cached_io::create(
			new remote_io(global_data.read_threads,
				get_msg_allocator(t->get_node_id()), mapper, t),
			global_data.table);
	num_ios++;
	return io_interface::ptr(io);
}

void part_global_cached_io_factory::destroy_io(io_interface &io)
{
	num_ios--;
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
	if (!safs::is_safs_init())
		throw io_exception("safs isn't init");
	safs_file f(get_sys_RAID_conf(), file_name);
	if (!f.exist())
		throw io_exception(boost::str(
					boost::format("safs file %1% doesn't exist") % file_name));

	for (int i = 0; i < global_data.raid_conf->get_num_disks(); i++) {
		std::string abs_path = global_data.raid_conf->get_disk(i).get_file_name()
			+ "/" + file_name;
		native_file f(abs_path);
		if (!f.exist())
			throw io_exception((boost::format(
							"the underlying file %1% doesn't exist")
						% abs_path).str());
	}

	file_mapper::ptr mapper = global_data.raid_conf->create_file_mapper(file_name);
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
			else
				throw io_exception("There is no page cache for global cache IO");
			break;
		case DIRECT_COMP_ACCESS:
			factory = new direct_comp_io_factory(mapper);
			break;
#ifdef PART_IO
		case PART_GLOBAL_ACCESS:
			if (global_data.global_cache)
				factory = new part_global_cached_io_factory(file_name);
			else
				throw io_exception(
						"There is no page cache for part'ed global cache IO");
			break;
#endif
		default:
			throw io_exception("a wrong access option");
	}
	return file_io_factory::shared_ptr(factory, destroy_io_factory());
}

io_interface::ptr create_io(file_io_factory::shared_ptr factory, thread *t)
{
	io_interface::ptr io = factory->create_io(t);
	io->set_owner(factory);
	return io;
}

void print_io_thread_stat()
{
	sleep(1);
	BOOST_FOREACH(disk_io_thread::ptr t, global_data.read_thread_set) {
		if (t)
			t->print_stat();
	}
}

void print_io_summary()
{
	size_t num_reads = 0;
	size_t num_read_bytes = 0;
	size_t num_writes = 0;
	size_t num_write_bytes = 0;

	sleep(1);
	BOOST_FOREACH(disk_io_thread::ptr t, global_data.read_thread_set) {
		if (t) {
			num_reads += t->get_num_reads();
			num_read_bytes += t->get_num_read_bytes();
			num_writes += t->get_num_writes();
			num_write_bytes += t->get_num_write_bytes();
		}
	}
	printf("It reads %ld bytes (in %ld reqs) and writes %ld bytes (in %ld reqs)\n",
			num_read_bytes, num_reads, num_write_bytes, num_writes);
}

ssize_t file_io_factory::get_file_size() const
{
	safs_file f(*global_data.raid_conf, name);
	return f.get_size();
}

bool is_safs_init()
{
	return global_data.raid_conf != NULL;
}

atomic_integer io_interface::io_counter;

io_interface::~io_interface()
{
	if (io_factory) {
		io_factory->collect_stat(*this);
		io_factory->destroy_io(*this);
	}
}

file_io_factory::file_io_factory(const std::string _name): name(_name)
{
	// It's possible that SAFS hasn't been initialized.
	if (global_data.raid_conf) {
		safs_file f(*global_data.raid_conf, name);
		header = f.get_header();
	}
}

namespace
{

class empty_io_select: public io_select
{
public:
	virtual bool add_io(io_interface::ptr io) {
		return true;
	}
	virtual int num_pending_ios() const {
		return 0;
	}
	virtual int wait4complete(int num_to_complete) {
		return 0;
	}
};

}

io_select::ptr create_io_select(const std::vector<io_interface::ptr> &ios)
{
	if (ios.empty())
		return io_select::ptr();

	// Let's try to find a valid I/O select from the I/O objects.
	io_select::ptr select;
	for (size_t i = 0; i < ios.size(); i++) {
		select = ios[i]->create_io_select();
		if (select)
			break;
	}
	// If none of the I/O objects can create a valid I/O object, it means
	// none of them actually need to access data from disks. We only need
	// to create an empty I/O select.
	if (select == NULL)
		select = io_select::ptr(new empty_io_select());
	for (size_t i = 0; i < ios.size(); i++)
		select->add_io(ios[i]);
	return select;
}

size_t wait4ios(safs::io_select::ptr select, size_t max_pending_ios)
{
	size_t num_pending;
	do {
		num_pending = select->num_pending_ios();

		// Figure out how many I/O requests we have to wait for in
		// this iteration.
		int num_to_process;
		if (num_pending > max_pending_ios)
			num_to_process = num_pending - max_pending_ios;
		else
			num_to_process = 0;
		select->wait4complete(num_to_process);

		// Test if all I/O instances have pending I/O requests left.
		// When a portion of a matrix is ready in memory and being processed,
		// it may result in writing data to another matrix. Therefore, we
		// need to process all completed I/O requests (portions with data
		// in memory) first and then count the number of new pending I/Os.
		num_pending = select->num_pending_ios();
	} while (num_pending > max_pending_ios);
	return num_pending;
}

const std::vector<int> &get_io_cpus()
{
	return global_data.io_cpus;
}

std::string get_supported_features()
{
	std::string ret;
#ifdef USE_NUMA
	ret += "+NUMA ";
#else
	ret += "-NUMA ";
#endif

#ifdef USE_LIBAIO
	ret += "+libaio ";
#else
	ret += "-libaio ";
#endif
	return ret;
}

}
