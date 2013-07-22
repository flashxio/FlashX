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
		assert(has_init && taken);
		return io;
	}

	void release() {
		pthread_spin_lock(&lock);
		taken = false;
		pthread_spin_unlock(&lock);
	}
};

static std::vector<io_tracker *> io_table;

int get_num_ios()
{
	return io_table.size();
}

void register_io(io_interface *io)
{
	// Make sure the index hasn't been set.
	assert(io->get_io_idx() < 0);

	io_table.push_back(new io_tracker(io));
	int idx = io_table.size() - 1;
	io->set_io_idx(idx);
}

io_interface *get_io(int idx)
{
	assert(idx >= 0);
	io_interface *io = io_table[idx]->get();
	assert(io);
	assert(io->get_io_idx() == idx);
	return io;
}

io_interface *allocate_io()
{
	for (unsigned i = 0; i < io_table.size(); i++) {
		io_interface *io = io_table[i]->take();
		if (io)
			return io;
	}
	return NULL;
}

void release_io(io_interface *io)
{
	io_table[io->get_io_idx()]->release();
}

std::vector<io_interface *> create_ios(const RAID_config &raid_conf,
		cache_config *cache_conf, const std::vector<int> &node_id_array,
		int nthreads, int access_option, long size, bool preload)
{
	int num_nodes = node_id_array.size();
	file_mapper *mapper = raid_conf.create_file_mapper();

	std::vector<int> indices;
	int num_files = mapper->get_num_files();
	for (int i = 0; i < num_files; i++)
		indices.push_back(i);
	// The partition contains all files.
	logical_file_partition global_partition(indices, mapper);

	std::vector<disk_read_thread *> read_threads(num_files);
	std::tr1::unordered_map<int, aio_complete_thread *> complete_threads;
	if (access_option == GLOBAL_CACHE_ACCESS
			|| access_option == PART_GLOBAL_ACCESS
			|| access_option == REMOTE_ACCESS) {
		// Create threads for helping process completed AIO requests.
		for (int i = 0; i < num_nodes; i++) {
			int node_id = node_id_array[i];
			complete_threads.insert(std::pair<int, aio_complete_thread *>(node_id,
						new aio_complete_thread(node_id)));
		}
		for (int k = 0; k < num_files; k++) {
			std::vector<int> indices(1, k);
			logical_file_partition partition(indices, mapper);
			// Create disk accessing threads.
			read_threads[k] = new disk_read_thread(partition, complete_threads,
					raid_conf.get_file(k).node_id);
		}
	}

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
#if ENABLE_AIO
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
#endif
				case REMOTE_ACCESS:
					io = new remote_disk_access(read_threads.data(),
							complete_threads[node_id], num_files,
							mapper, node_id);
					register_io(io);
					ios.push_back(io);
					break;
				case GLOBAL_CACHE_ACCESS:
					{
						io_interface *underlying = new remote_disk_access(
								read_threads.data(), complete_threads[node_id],
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
								read_threads.data(), complete_threads[node_id],
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
