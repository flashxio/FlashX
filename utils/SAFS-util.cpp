/**
 * Copyright 2013 Da Zheng
 *
 * This file is part of SAFSlib.
 *
 * SAFSlib is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * SAFSlib is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with SAFSlib.  If not, see <http://www.gnu.org/licenses/>.
 */

#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>

#include <string>

#include "io_interface.h"
#include "native_file.h"
#include "safs_file.h"

const int BUF_SIZE = 1024 * 64 * PAGE_SIZE;
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
		native_file f(ext_file);
		file_size = f.get_size();
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
	size_t verified_bytes;
public:
	verify_callback(data_source *source) {
		this->source = source;
		orig_buf = (char *) malloc(BUF_SIZE);
		verified_bytes = 0;
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
		verified_bytes += read_bytes;
		verify_bytes(rqs[0]->get_buf(), orig_buf, read_bytes);
	}

	size_t get_verified_bytes() const {
		return verified_bytes;
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
	thread *curr_thread = thread::get_curr_thread();
	assert(curr_thread);
	io_interface *io = factory->create_io(curr_thread);
	data_source *source;
	if (ext_file.empty())
		source = new synthetic_data_source(DATA_SIZE);
	else
		source = new file_data_source(ext_file);
	verify_callback *cb = new verify_callback(source);
	io->set_callback(cb);

	size_t file_size = source->get_size();
	assert(factory->get_file_size() >= file_size);
	file_size = ROUNDUP(file_size, BUF_SIZE);
	char *buf = (char *) valloc(BUF_SIZE);
	for (off_t off = 0; off < file_size; off += BUF_SIZE) {
		data_loc_t loc(io->get_file_id(), off);
		io_request req(buf, loc, BUF_SIZE, READ, io, 0);
		io->access(&req, 1);
		io->wait4complete(1);
	}
	printf("verify %ld bytes\n", cb->get_verified_bytes());
	
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

	safs_file file(get_sys_RAID_conf(), int_file_name);
	if (!file.exist()) {
		fprintf(stderr, "%s doesn't exist\n", int_file_name.c_str());
		return;
	}

	file_io_factory *factory = create_io_factory(int_file_name, REMOTE_ACCESS);
	assert(factory);
	thread *curr_thread = thread::get_curr_thread();
	assert(curr_thread);
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
		data_loc_t loc(io->get_file_id(), off);
		io_request req(buf, loc, write_bytes, WRITE, io, 0);
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
	safs_file file(get_sys_RAID_conf(), file_name);
	file.create_file(file_size);
	printf("create file %s of %ld bytes\n", file_name.c_str(),
			file.get_file_size());
}

void print_help();

void comm_help(int argc, char *argv[])
{
	print_help();
}

void comm_list(int argc, char *argv[])
{
	std::set<std::string> files;
	const RAID_config &conf = get_sys_RAID_conf();

	// First find all individual file names in the root directories.
	for (int i = 0; i < conf.get_num_disks(); i++) {
		std::string dir_name = conf.get_disk(i).name;
		native_dir dir(dir_name);
		std::vector<std::string> file_names;
		dir.read_all_files(file_names);
		files.insert(file_names.begin(), file_names.end());
	}

	for (std::set<std::string>::const_iterator it = files.begin();
			it != files.end(); it++) {
		safs_file file(conf, *it);
		if (file.exist()) {
			printf("%s: %ld bytes\n", file.get_name().c_str(),
					file.get_file_size());
		}
		else {
			printf("%s is corrupted\n", file.get_name().c_str());
		}
	}
}

void comm_delete_file(int argc, char *argv[])
{
	if (argc < 1) {
		fprintf(stderr, "delete file_name\n");
		return;
	}

	std::string file_name = argv[0];
	safs_file file(get_sys_RAID_conf(), file_name);
	if (!file.exist()) {
		fprintf(stderr, "%s doesn't exist\n", file_name.c_str());
		return;
	}
	file.delete_file();
}

typedef void (*command_func_t)(int argc, char *argv[]);

struct command
{
	std::string name;
	command_func_t func;
	std::string help_info;
};

struct command commands[] = {
	{"create", comm_create_file,
		"create file_name size: create a file with the specified size"},
	{"delete", comm_delete_file,
		"delete file_name: delete a file"},
	{"help", comm_help,
		"help: print the help info"},
	{"list", comm_list, "list: list existing files in SAFS"},
	{"load", comm_load_file2fs,
		"load file_name [ext_file]: load data to the file"},
	{"verify", comm_verify_file,
		"verify file_name [ext_file]: verify data in the file"},
};

int get_num_commands()
{
	return sizeof(commands) / sizeof(commands[0]);
}

const command *get_command(const std::string &name)
{
	int num_comms = get_num_commands();
	for (int i = 0; i < num_comms; i++) {
		if (commands[i].name == name)
			return &commands[i];
	}
	return NULL;
}

void print_help()
{
	printf("SAFS-util conf_file command ...\n");
	printf("The supported commands are\n");
	int num_commands = get_num_commands();
	for (int i =0; i < num_commands; i++) {
		printf("\t%s\n", commands[i].help_info.c_str());
	}
}

/**
 * This is a utility tool for the SA-FS.
 */

int main(int argc, char *argv[])
{
	if (argc < 3) {
		print_help();
		exit(-1);
	}

	std::string conf_file = argv[1];
	std::string command = argv[2];

	config_map configs(conf_file);
	init_io_system(configs);

	const struct command *comm = get_command(command);
	if (comm == NULL) {
		fprintf(stderr, "wrong command %s\n", comm->name.c_str());
		print_help();
		return -1;
	}
	comm->func(argc - 3, argv + 3);
}
