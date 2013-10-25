#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>

#include <string>

#include "io_interface.h"

const int BUF_SIZE = 512 * PAGE_SIZE;

ssize_t complete_read(int fd, char *buf, size_t count)
{
	ssize_t bytes = 0;
	do {
		ssize_t ret = read(fd, buf, count);
		if (ret < 0)
			return ret;
		if (ret == 0)
			return bytes;
		bytes += ret;
		count -= ret;
		buf += ret;
	} while (count > 0);
	return bytes;
}

void comm_load_file2fs(int argc, char *argv[])
{
	if (argc < 2) {
		fprintf(stderr, "load ext_file file_name\n");
		fprintf(stderr, "ext_file is the file in the external file system\n");
		fprintf(stderr, "file_name is the file name in the SA-FS file system\n");
		exit(-1);
	}

	std::string ext_file = argv[0];
	std::string int_file_name = argv[1];

	file_io_factory *factory = create_io_factory(int_file_name, REMOTE_ACCESS);
	assert(factory->get_file_size() >= get_file_size(ext_file.c_str()));
	thread *curr_thread = thread::represent_thread(0);
	io_interface *io = factory->create_io(curr_thread);

	int fd = open(ext_file.c_str(), O_RDONLY);
	if (fd < 0) {
		perror("open");
		exit(-1);
	}

	char *buf = (char *) valloc(BUF_SIZE);
	off_t off = 0;

	while (true) {
		ssize_t ret = complete_read(fd, buf, BUF_SIZE);
		if (ret < 0) {
			perror("complete_read");
			exit(-1);
		}
		if (ret == 0)
			break;

		ssize_t write_bytes = ROUNDUP(ret, 512);
		io_request req(buf, off, write_bytes, WRITE, io, 0);
		io->access(&req, 1);
		io->wait4complete(1);
		off += write_bytes;
	}
	printf("write all data\n");
	
	close(fd);
	io->cleanup();
	factory->destroy_io(io);
}

/**
 * This is a utility tool for the SA-FS.
 */

int main(int argc, char *argv[])
{
	if (argc < 3) {
		fprintf(stderr, "util conf_file command ...\n");
		exit(-1);
	}

	std::string conf_file = argv[1];
	std::string command = argv[2];

	config_map configs(conf_file);
	init_io_system(configs);

	if (command == "load") {
		comm_load_file2fs(argc - 3, argv + 3);
	}
}
