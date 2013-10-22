#include <sys/types.h>
#include <sys/stat.h>
#include <stdio.h>
#include <string.h>
#include <fcntl.h>
#include <stdlib.h>
#include <unistd.h>
#include <numaif.h>
#include <numa.h>
#include <sys/time.h>
#include <assert.h>

#include <vector>

#include "thread.h"
#include "file_mapper.h"

class write_thread: public thread
{
	int fd;
	std::vector<off_t> block_ids;
	file_mapper *mapper;
	int block_size;
	int thread_id;

public:
	write_thread(int id, const file_info &file, file_mapper *mapper,
			int block_size): thread("write-thread", file.node_id) {
		this->thread_id = id;
		this->block_size = block_size;
		this->mapper = mapper;
		std::string file_name = file.name + "/test";
		printf("open file %s for write in thread %d\n", file_name.c_str(), id);
		fd = open(file_name.c_str(), O_DIRECT | O_RDWR | O_CREAT, 00644);
		if (fd < 0) {
			perror("open");
			::exit(1);
		}
	}

	void add_block(off_t id) {
		block_ids.push_back(id);
	}

	int get_num_blocks() const {
		return (int) block_ids.size();
	}

	void run();
};

void write_buf(int fd, char *buf, int size, off_t off)
{
	while (size > 0) {
		ssize_t ret = pwrite(fd, buf, size, off);
		if (ret < 0) {
			perror("write");
			exit(1);
		}
		size -= ret;
		buf += ret;
		off += ret;
	}
}

void init_block(char *buf, int size, off_t block_id)
{
	off_t *long_pointer = (off_t *) buf;
	int num_longs = size / sizeof(off_t);
	long value = block_id * num_longs;
	for (int i = 0; i < num_longs; i++, value++) {
		long_pointer[i] = value;
	}
}

void write_thread::run()
{
	printf("thread %d starts to run\n", thread_id);
	char *buf = (char *) valloc(block_size * PAGE_SIZE);
	for (size_t i = 0; i < block_ids.size(); i++) {
		off_t block_off = block_ids[i] * block_size;
		block_identifier bid;
		mapper->map(block_off, bid);
		assert(bid.idx == thread_id);
		off_t off_in_file = bid.off * PAGE_SIZE;
		init_block(buf, block_size * PAGE_SIZE, block_ids[i]);
		write_buf(fd, buf, block_size * PAGE_SIZE, bid.off * PAGE_SIZE);
	}
	stop();
}

int main(int argc, char *argv[])
{
	if (argc < 5) {
		fprintf(stderr, "create_file file_file size RAID_mapping block_size\n");
		exit(1);
	}

	char *file_file = argv[1];
	long size = str2size(argv[2]);
	const int block_size = atoi(argv[4]);		// number of pages.
	printf("create a file of %ld bytes\n", size);
	std::vector<file_info> data_files;
	retrieve_data_files(file_file, data_files);
	std::string mapping = argv[3];
	file_mapper *mapper;
	if (mapping == "RAID0") {
		mapper = new RAID0_mapper(data_files, block_size);
	}
	else if (mapping == "RAID5") {
		mapper = new RAID5_mapper(data_files, block_size);
	}
	else if (mapping == "HASH") {
		mapper = new hash_mapper(data_files, block_size);
	}
	else {
		fprintf(stderr, "wrong mapping\n");
		exit(1);
	}

	std::vector<write_thread *> threads(data_files.size(), NULL);
	for (size_t i = 0; i < data_files.size(); i++) {
		threads[i] = new write_thread(i, data_files[i], mapper, block_size);
	}
	const int num_blocks = size / PAGE_SIZE / block_size;
	for (int i = 0; i < num_blocks; i++) {
		off_t block_off = i * block_size;
		int idx = mapper->map2file(block_off);
		threads[idx]->add_block(i);
	}

	/* write data of `size' bytes to the file. */
	struct timeval start_time, end_time;
	gettimeofday(&start_time, NULL);
	for (size_t i = 0; i < threads.size(); i++) {
		printf("thread %ld gets %d blocks\n", i, threads[i]->get_num_blocks());
		threads[i]->start();
	}
	for (size_t i = 0; i < threads.size(); i++) {
		threads[i]->join();
	}
	gettimeofday(&end_time, NULL);
	long tot_size = ((long) num_blocks) * PAGE_SIZE * block_size;
	printf("write %ld bytes, takes %f seconds\n",
			tot_size, end_time.tv_sec - start_time.tv_sec
			+ ((float)(end_time.tv_usec - start_time.tv_usec))/1000000);

	return 0;
}
