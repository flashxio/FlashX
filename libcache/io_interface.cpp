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

class io_tracker
{
	io_interface *io;
	bool has_init;
	pthread_spinlock_t lock;
	bool taken;
public:
	io_tracker(io_interface *io) {
		this->io = io;
		pthread_spin_init(&lock, PTHREAD_PROCESS_PRIVATE);
		taken = false;
		has_init = false;
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
	/* node id <-> ios */
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

	void register_io(io_interface *io) {
		pthread_spin_lock(&lock);
		// Make sure the index hasn't been set.
		assert(io->get_io_idx() < 0);

		table.push_back(new io_tracker(io));
		int idx = table.size() - 1;
		io->set_io_idx(idx);
		pthread_spin_unlock(&lock);
	}

	io_interface *get_io(int idx) {
		pthread_spin_lock(&lock);
		assert(idx >= 0);
		io_interface *io = table[idx]->get();
		assert(io);
		assert(io->get_io_idx() == idx);
		pthread_spin_unlock(&lock);
		return io;
	}

	io_interface *allocate_io(int node_id) {
		pthread_spin_lock(&lock);
		io_interface *io = NULL;
		for (unsigned i = 0; i < table.size(); i++) {
			io_interface *tmp = table[i]->get();
			if (tmp->get_node_id() == node_id) {
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

int get_num_ios()
{
	return ios.size();
}

void register_io(io_interface *io)
{
	ios.register_io(io);
}

io_interface *get_io(int idx)
{
	return ios.get_io(idx);
}

io_interface *allocate_io(int node_id)
{
	return ios.allocate_io(node_id);
}

void release_io(io_interface *io)
{
	ios.release_io(io);
}

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

const std::vector<disk_read_thread *> &get_read_threads()
{
	return global_data.read_threads;
}

std::vector<io_interface *> create_ios(const RAID_config &raid_conf,
		cache_config *cache_conf, const std::vector<int> &node_id_array,
		int nthreads, int access_option, long size, bool preload)
{
	int num_nodes = node_id_array.size();
	assert(num_nodes <= nthreads && nthreads % num_nodes == 0);
	file_mapper *mapper = raid_conf.create_file_mapper();

	std::vector<int> indices;
	int num_files = mapper->get_num_files();
	for (int i = 0; i < num_files; i++)
		indices.push_back(i);
	// The partition contains all files.
	logical_file_partition global_partition(indices, mapper);

	pthread_mutex_lock(&global_data.mutex);
	if (access_option == GLOBAL_CACHE_ACCESS
			|| access_option == PART_GLOBAL_ACCESS
			|| access_option == REMOTE_ACCESS) {
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
		for (int i = 0; i < num_files; i++) {
			global_data.read_threads[i]->open_file(mapper);
		}
	}
	pthread_mutex_unlock(&global_data.mutex);

	// TODO the number of threads shouldn't be a global variable.
	if (access_option == PART_GLOBAL_ACCESS)
		part_global_cached_io::set_num_threads(nthreads);

	std::vector<io_interface *> ios;
	int nthreads_per_node = nthreads / num_nodes;
	for (int i = 0, j = 0; i < num_nodes; i++) {
		int node_id = node_id_array[i];
		for (int k = 0; k < nthreads_per_node; k++, j++) {
			io_interface *io;
			switch (access_option) {
				case READ_ACCESS:
					io = new buffered_io(global_partition, node_id);
					register_io(io);
					ios.push_back(io);
					break;
				case DIRECT_ACCESS:
					io = new direct_io(global_partition, node_id);
					register_io(io);
					ios.push_back(io);
					break;
				case AIO_ACCESS:
					{
						int depth_per_file = AIO_DEPTH_PER_FILE / nthreads;
						if (depth_per_file == 0)
							depth_per_file = 1;
						std::tr1::unordered_map<int, aio_complete_thread *> no_complete_threads;
						io = new async_io(global_partition, no_complete_threads,
								depth_per_file, node_id);
						register_io(io);
						ios.push_back(io);
					}
					break;
				case REMOTE_ACCESS:
					io = new remote_disk_access(global_data.read_threads.data(),
							global_data.complete_threads[node_id], num_files,
							mapper, node_id);
					register_io(io);
					ios.push_back(io);
					break;
				case GLOBAL_CACHE_ACCESS:
					{
						io_interface *underlying = new remote_disk_access(
								global_data.read_threads.data(),
								global_data.complete_threads[node_id],
								num_files, mapper, node_id);
						global_cached_io *io = new global_cached_io(underlying,
								cache_conf);
						if (preload && j == 0)
							io->preload(0, size);
						register_io(io);
						ios.push_back(io);
					}
					break;
				case PART_GLOBAL_ACCESS:
					{
						io_interface *underlying = new remote_disk_access(
								global_data.read_threads.data(),
								global_data.complete_threads[node_id],
								num_files, mapper, node_id);
						part_global_cached_io *io = new part_global_cached_io(
								num_nodes, underlying, j, cache_conf);
						if (preload && j == 0)
							io->preload(0, size);
						register_io(io);
						ios.push_back(io);
					}
					break;
				default:
					fprintf(stderr, "wrong access option\n");
					exit(1);
			}
		}
	}
	return ios;
}
