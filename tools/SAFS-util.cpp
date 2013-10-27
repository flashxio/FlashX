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

int comp_bytes(const char *bs1, const char *bs2, int size)
{
	for (int i = 0; i < size; i++) {
		if (bs1[i] != bs2[i])
			return bs1[i] - bs2[i];
	}
	return 0;
}

class verify_callback: public callback
{
	char *orig_buf;
	int fd;
	size_t file_size;
public:
	verify_callback(const std::string &ext_file) {
		fd = open(ext_file.c_str(), O_RDONLY);
		if (fd < 0) {
			perror("open");
			exit(-1);
		}
		file_size = get_file_size(ext_file.c_str());
		orig_buf = (char *) malloc(BUF_SIZE);
	}

	~verify_callback() {
		close(fd);
		free(orig_buf);
	}

	int invoke(io_request *rqs[], int num) {
		int read_bytes = min(rqs[0]->get_size(),
				file_size - rqs[0]->get_offset());
		assert(read_bytes > 0);
		ssize_t ret = complete_read(fd, orig_buf, read_bytes);
		if (ret < 0) {
			perror("complete_read");
			exit(-1);
		}
		fprintf(stderr, "verify block %lx of %d bytes\n", rqs[0]->get_offset(), read_bytes);
		assert(ret == read_bytes);
		assert(comp_bytes(rqs[0]->get_buf(), orig_buf, read_bytes) == 0);
	}
};

void comm_verify_file(int argc, char *argv[])
{
	if (argc < 2) {
		fprintf(stderr, "verify ext_file file_name");
		fprintf(stderr, "ext_file is the file in the external file system\n");
		fprintf(stderr, "file_name is the file name in the SA-FS file system\n");
		exit(-1);
	}

	std::string ext_file = argv[0];
	std::string int_file_name = argv[1];

	file_io_factory *factory = create_io_factory(int_file_name, REMOTE_ACCESS);
	size_t file_size = get_file_size(ext_file.c_str());
	assert(factory->get_file_size() >= file_size);
	thread *curr_thread = thread::represent_thread(0);
	io_interface *io = factory->create_io(curr_thread);
	io->set_callback(new verify_callback(ext_file));

	file_size = ROUNDUP(file_size, BUF_SIZE);
	char *buf = (char *) valloc(BUF_SIZE);
	for (off_t off = 0; off < file_size; off += BUF_SIZE) {
		io_request req(buf, off, BUF_SIZE, READ, io, 0);
		io->access(&req, 1);
		io->wait4complete(1);
	}
	printf("verify all data\n");
	
	io->cleanup();
	delete io->get_callback();
	factory->destroy_io(io);
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
	else if (command == "verify") {
		comm_verify_file(argc - 3, argv + 3);
	}
}
