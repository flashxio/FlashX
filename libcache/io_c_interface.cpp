#include <fcntl.h>
#include <stdlib.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>

#include <vector>

#include "io_interface.h"
#include "concurrency.h"
#include "RAID_config.h"
#include "cache.h"
#include "cache_config.h"

long cache_size;
int cache_type;
int access_option;
int num_threads;
int num_nodes;
int RAID_mapping_option;
int RAID_block_size;

void ssd_io_init(const char *name)
{
	static atomic_unsigned_integer has_init;

	// this is the first time it is called.
	if (has_init.inc(1) == 1) {
		// Init RAID configuration.
		std::vector<file_info> files;
		int num_files = retrieve_data_files(name, files);
		printf("There are %d data files\n", num_files);
		RAID_config raid_conf(files, RAID_mapping_option, RAID_block_size);

		// Init node id array.
		std::set<int> node_ids = raid_conf.get_node_ids();
		// In this way, we can guarantee that the cache is created
		// on the nodes with the data files.
		for (int i = 0; i < num_nodes
				&& node_ids.size() < (unsigned) num_nodes; i++)
			node_ids.insert(i);
		std::vector<int> node_id_array;
		// We only get a specified number of nodes.
		for (std::set<int>::const_iterator it = node_ids.begin();
				it != node_ids.end() && (int) node_id_array.size() < num_nodes; it++)
			node_id_array.push_back(*it);
		printf("There are %ld nodes\n", node_id_array.size());

		// Init cache configuration.
		cache_config *cache_conf = NULL;
		if (access_option == GLOBAL_CACHE_ACCESS)
			cache_conf = new even_cache_config(cache_size, cache_type,
					node_id_array);
		else if (access_option == PART_GLOBAL_ACCESS) {
			assert(num_nodes == 4);
			cache_conf = new test_cache_config(cache_size, cache_type,
					node_id_array);
		}

		create_ios(raid_conf, cache_conf, node_id_array, access_option,
				0, false);
	}
}

extern "C" {
int ssd_create(const char *name, size_t tot_size)
{
	std::vector<file_info> files;
	retrieve_data_files(name, files);
	size_t file_size = tot_size / files.size();
	for (unsigned i = 0; i < files.size(); i++) {
		int fd = open(files[i].name.c_str(), O_RDWR | O_CREAT,
				S_IRUSR | S_IWUSR | S_IRGRP | S_IROTH);
		if (fd < 0) {
			perror("open");
			abort();
		}
		int err = posix_fallocate(fd, 0, file_size);
		if (err) {
			fprintf(stderr, "posix_fallocate error: %s\n", strerror(err));
			abort();
		}
	}
	return 0;
}

int ssd_open(const char *name)
{
	static atomic_unsigned_integer open_id;
	ssd_io_init(name);
	int id = open_id.inc(1) - 1;
	printf("open id: %d\n", id);
	assert(id < get_num_ios());

	io_interface *io = get_io(id);
	io->init();
	return id;
}

ssize_t ssd_read(int fd, void *buf, size_t count, off_t off)
{
	io_interface *io = get_io(fd);
	io_status status = io->access((char *) buf, off, count, READ);
	assert(status == IO_OK);
	return status.get_priv_data();
}

ssize_t ssd_write(int fd, void *buf, size_t count, off_t off)
{
	io_interface *io = get_io(fd);
	io_status status = io->access((char *) buf, off, count, WRITE);
	assert(status == IO_OK);
	return status.get_priv_data();
}

int ssd_close(int fd)
{
	io_interface *io = get_io(fd);
	io->cleanup();
	return 0;
}

int ssd_delete(const char *name)
{
	std::vector<file_info> files;
	retrieve_data_files(name, files);
	for (unsigned i = 0; i < files.size(); i++) {
		int ret = unlink(files[i].name.c_str());
		if (ret < 0) {
			perror("unlink");
			abort();
		}
	}
	return 0;
}

size_t ssd_get_filesize(const char *name)
{
	std::vector<file_info> files;
	retrieve_data_files(name, files);
	size_t tot_size = 0;
	for (unsigned i = 0; i < files.size(); i++) {
		struct stat stat_buf;
		if (stat(files[i].name.c_str(), &stat_buf) != 0) {
			perror("stat");
			abort();
		}
		tot_size += stat_buf.st_size;
	}
	return tot_size;
}

}
