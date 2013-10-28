#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>

#include <string>

#include "io_interface.h"

const int BUF_SIZE = 512 * PAGE_SIZE;
const size_t DATA_SIZE = 10L * 1024 * 1024 * 1024;

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

int verify_bytes(const char *bs1, const char *bs2, int size)
{
	for (int i = 0; i < size; i++) {
		assert(bs1[i] == bs2[i]);
	}
	return 0;
}

class data_source
{
public:
	virtual ssize_t get_data(off_t off, size_t size, char *buf) const = 0;
	virtual size_t get_size() const = 0;
};

class file_data_source: public data_source
{
	int fd;
	size_t file_size;
public:
	file_data_source(const std::string &ext_file) {
		fd = open(ext_file.c_str(), O_RDONLY);
		if (fd < 0) {
			perror("open");
			exit(-1);
		}
		file_size = get_file_size(ext_file.c_str());
	}

	virtual ssize_t get_data(off_t off, size_t size, char *buf) const {
		long new_off = lseek(fd, off, SEEK_SET);
		assert(new_off == off);
		ssize_t ret = complete_read(fd, buf, size);
		if (ret < 0) {
			perror("complete_read");
			exit(-1);
		}
		return ret;
	}

	virtual size_t get_size() const {
		return file_size;
	}
};

class synthetic_data_source: public data_source
{
	size_t size;
public:
	synthetic_data_source(size_t size) {
		this->size = size;
	}

	virtual ssize_t get_data(off_t off, size_t size, char *buf) const {
		off_t *long_pointer = (off_t *) buf;
		int num_longs = size / sizeof(off_t);
		long value = off / sizeof(off_t);
		for (int i = 0; i < num_longs; i++, value++) {
			long_pointer[i] = value;
		}
		return size;
	}

	virtual size_t get_size() const {
		return size;
	}
};

class verify_callback: public callback
{
	char *orig_buf;
	data_source *source;
public:
	verify_callback(data_source *source) {
		this->source = source;
		orig_buf = (char *) malloc(BUF_SIZE);
	}

	~verify_callback() {
		free(orig_buf);
	}

	int invoke(io_request *rqs[], int num) {
		size_t read_bytes = min<size_t>(rqs[0]->get_size(),
				source->get_size() - rqs[0]->get_offset());
		assert(read_bytes > 0);
		size_t ret = source->get_data(rqs[0]->get_offset(), read_bytes, orig_buf);
		fprintf(stderr, "verify block %lx of %ld bytes\n", rqs[0]->get_offset(), read_bytes);
		assert(ret == read_bytes);
		verify_bytes(rqs[0]->get_buf(), orig_buf, read_bytes);
	}
};

void comm_verify_file(int argc, char *argv[])
{
	if (argc < 1) {
		fprintf(stderr, "verify file_name [ext_file]");
		fprintf(stderr, "file_name is the file name in the SA-FS file system\n");
		fprintf(stderr, "ext_file is the file in the external file system\n");
		exit(-1);
	}

	std::string int_file_name = argv[0];
	std::string ext_file;
	if (argc >= 2)
		ext_file = argv[1];

	file_io_factory *factory = create_io_factory(int_file_name, REMOTE_ACCESS);
	thread *curr_thread = thread::represent_thread(0);
	io_interface *io = factory->create_io(curr_thread);
	data_source *source;
	if (ext_file.empty())
		source = new synthetic_data_source(DATA_SIZE);
	else
		source = new file_data_source(ext_file);
	io->set_callback(new verify_callback(source));

	size_t file_size = source->get_size();
	assert(factory->get_file_size() >= file_size);
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
	if (argc < 1) {
		fprintf(stderr, "load file_name [ext_file]\n");
		fprintf(stderr, "file_name is the file name in the SA-FS file system\n");
		fprintf(stderr, "ext_file is the file in the external file system\n");
		exit(-1);
	}

	std::string int_file_name = argv[0];
	std::string ext_file;
	if (argc >= 2)
		ext_file = argv[1];

	file_io_factory *factory = create_io_factory(int_file_name, REMOTE_ACCESS);
	thread *curr_thread = thread::represent_thread(0);
	io_interface *io = factory->create_io(curr_thread);

	data_source *source;
	if (ext_file.empty()) {
		printf("use synthetic data\n");
		source = new synthetic_data_source(DATA_SIZE);
	}
	else {
		printf("use file %s\n", ext_file.c_str());
		source = new file_data_source(ext_file);
	}
	assert(factory->get_file_size() >= source->get_size());
	printf("source size: %ld\n", source->get_size());

	char *buf = (char *) valloc(BUF_SIZE);
	off_t off = 0;

	while (off < source->get_size()) {
		size_t size = min<size_t>(BUF_SIZE, source->get_size() - off);
		size_t ret = source->get_data(off, size, buf);
		assert(ret == size);
		ssize_t write_bytes = ROUNDUP(ret, 512);
		io_request req(buf, off, write_bytes, WRITE, io, 0);
		io->access(&req, 1);
		io->wait4complete(1);
		off += write_bytes;
	}
	printf("write all data\n");
	
	io->cleanup();
	factory->destroy_io(io);
}

void comm_create_file(int argc, char *argv[])
{
	if (argc < 2) {
		fprintf(stderr, "create file_name size\n");
		fprintf(stderr, "file_name is the file name in the SA-FS file system\n");
		exit(-1);
	}

	std::string file_name = argv[0];
	size_t file_size = str2size(argv[1]);
	safs_file file(file_name);
	file.create_file(file_size);
	printf("create file %s of %ld bytes\n", file_name.c_str(),
			file.get_file_size());
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
	else if (command == "create") {
		comm_create_file(argc - 3, argv + 3);
	}
}
