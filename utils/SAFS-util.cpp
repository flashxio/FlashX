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

#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>

#include <string>
#include <boost/format.hpp>

#include "io_interface.h"
#include "native_file.h"
#include "safs_file.h"
#include "file_mapper.h"
#include "RAID_config.h"

const int BUF_SIZE = 1024 * 64 * PAGE_SIZE;

config_map::ptr configs;

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
		BOOST_VERIFY(new_off == off);
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
	const file_mapper *fmapper;
public:
	verify_callback(data_source *source, const file_mapper *fmapper) {
		this->source = source;
		this->fmapper = fmapper;
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
		BOOST_VERIFY(ret == read_bytes);
		verified_bytes += read_bytes;
		for (size_t i = 0; i < read_bytes; i++) {
			if (rqs[0]->get_buf()[i] != orig_buf[i]) {
				struct block_identifier bid;
				fmapper->map((rqs[0]->get_offset() + i) / PAGE_SIZE, bid);
				ABORT_MSG(boost::format(
							"bytes at %1% (in partition %2%) doesn't match")
						% (rqs[0]->get_offset() + i) % bid.idx);
			}
		}
		return 0;
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

	init_io_system(configs, false);
	file_io_factory::shared_ptr factory = create_io_factory(int_file_name,
			REMOTE_ACCESS);
	thread *curr_thread = thread::get_curr_thread();
	assert(curr_thread);
	io_interface::ptr io = factory->create_io(curr_thread);
	data_source *source;
	if (ext_file.empty())
		source = new synthetic_data_source(factory->get_file_size());
	else
		source = new file_data_source(ext_file);
	const RAID_config &conf = get_sys_RAID_conf();
	verify_callback *cb = new verify_callback(source, conf.create_file_mapper());
	io->set_callback(cb);

	ssize_t file_size = source->get_size();
	printf("verify %ld bytes\n", file_size);
	assert(factory->get_file_size() >= file_size);
	file_size = ROUNDUP(file_size, BUF_SIZE);
	char *buf = (char *) valloc(BUF_SIZE);
	for (off_t off = 0; off < file_size; off += BUF_SIZE) {
		data_loc_t loc(io->get_file_id(), off);
		io_request req(buf, loc, BUF_SIZE, READ);
		io->access(&req, 1);
		io->wait4complete(1);
	}
	printf("verify %ld bytes\n", cb->get_verified_bytes());
	
	io->cleanup();
	delete io->get_callback();
}

void comm_load_file2fs(int argc, char *argv[])
{
	if (argc < 2) {
		fprintf(stderr, "load file_name ext_file\n");
		fprintf(stderr, "file_name is the file name in the SA-FS file system\n");
		fprintf(stderr, "ext_file is the file in the external file system\n");
		exit(-1);
	}

	std::string int_file_name = argv[0];
	std::string ext_file = argv[1];

	configs->add_options("writable=1");
	init_io_system(configs, false);
	data_source *source = new file_data_source(ext_file);

	safs_file file(get_sys_RAID_conf(), int_file_name);
	// If the file in SAFS doesn't exist, create a new one.
	if (!file.exist()) {
		safs_file file(get_sys_RAID_conf(), int_file_name);
		file.create_file(source->get_size());
		printf("create file %s of %ld bytes\n", int_file_name.c_str(),
				file.get_file_size());
	}

	file_io_factory::shared_ptr factory = create_io_factory(int_file_name,
			REMOTE_ACCESS);
	assert(factory);
	assert((size_t) factory->get_file_size() >= source->get_size());
	printf("source size: %ld\n", source->get_size());

	thread *curr_thread = thread::get_curr_thread();
	assert(curr_thread);
	io_interface::ptr io = factory->create_io(curr_thread);

	char *buf = (char *) valloc(BUF_SIZE);
	off_t off = 0;

	while (off < (off_t) source->get_size()) {
		size_t size = min<size_t>(BUF_SIZE, source->get_size() - off);
		size_t ret = source->get_data(off, size, buf);
		assert(ret == size);
		ssize_t write_bytes = ROUNDUP(ret, 512);
		data_loc_t loc(io->get_file_id(), off);
		io_request req(buf, loc, write_bytes, WRITE);
		io->access(&req, 1);
		io->wait4complete(1);
		off += write_bytes;
	}
	printf("write all data\n");
	
	io->cleanup();
}

void comm_load_part_file2fs(int argc, char *argv[])
{
	if (argc < 3) {
		fprintf(stderr, "load_part file_name ext_file part_id\n");
		fprintf(stderr, "file_name is the file name in the SA-FS file system\n");
		fprintf(stderr, "ext_file is the file in the external file system\n");
		fprintf(stderr, "part_id is the partition of the file loaded to SAFS\n");
		exit(-1);
	}

	std::string int_file_name = argv[0];
	std::string ext_file = argv[1];
	int part_id = atoi(argv[2]);

	configs->add_options("writable=1");
	init_io_system(configs, false);
	const RAID_config &conf = get_sys_RAID_conf();
	file_mapper *fmapper = conf.create_file_mapper();
	std::string part_path = fmapper->get_file_name(part_id) + "/"
		+ int_file_name;

	// Create the directory in the native filesystem for the partition
	// of the specified file.
	native_dir dir(part_path);
	assert(!dir.exist());
	bool ret = dir.create_dir(false);
	BOOST_VERIFY(ret);

	std::string file_path = part_path + "/" + std::string(argv[2]);
	FILE *out_f = fopen(file_path.c_str(), "w");
	assert(out_f);
	FILE *in_f = fopen(ext_file.c_str(), "r");
	assert(in_f);
	native_file nfile(ext_file);
	size_t ext_file_size = nfile.get_size();
	// The last write location in the partition file in SAFS.
	off_t expected_write_pos = 0;
	const size_t block_size = fmapper->STRIPE_BLOCK_SIZE * PAGE_SIZE;
	std::unique_ptr<char[]> buf = std::unique_ptr<char[]>(new char[block_size]);
	for (size_t off = 0; off < ext_file_size; off += block_size) {
		struct block_identifier bid;
		fmapper->map(off / PAGE_SIZE, bid);
		// If the block doesn't belong to the specified partition, skip it.
		if (bid.idx != part_id)
			continue;
		int ret = fseek(in_f, off, SEEK_SET);
		if (ret < 0) {
			perror("fseek");
			exit(1);
		}
		size_t remain_size = ext_file_size - off;
		size_t read_size = min(block_size, remain_size);
		size_t rret = fread(buf.get(), read_size, 1, in_f);
		assert(rret == 1);
		assert(expected_write_pos == bid.off * PAGE_SIZE);
		rret = fwrite(buf.get(), read_size, 1, out_f);
		assert(rret == 1);
		expected_write_pos = bid.off * PAGE_SIZE + read_size;
	}

	fclose(out_f);
	fclose(in_f);
}

void comm_create_file(int argc, char *argv[])
{
	if (argc < 2) {
		fprintf(stderr, "create file_name size\n");
		fprintf(stderr, "file_name is the file name in the SA-FS file system\n");
		exit(-1);
	}

	configs->add_options("writable=1");
	init_io_system(configs, false);
	std::string file_name = argv[0];
	size_t file_size = str2size(argv[1]);
	safs_file file(get_sys_RAID_conf(), file_name);
	file.create_file(file_size);
	printf("create file %s of %ld bytes\n", file_name.c_str(),
			file.get_file_size());

	file_io_factory::shared_ptr factory = create_io_factory(file_name,
			REMOTE_ACCESS);
	assert(factory);
	data_source *source = new synthetic_data_source(factory->get_file_size());
	assert((size_t) factory->get_file_size() >= source->get_size());
	printf("source size: %ld\n", source->get_size());

	thread *curr_thread = thread::get_curr_thread();
	assert(curr_thread);
	io_interface::ptr io = factory->create_io(curr_thread);

	char *buf = (char *) valloc(BUF_SIZE);
	off_t off = 0;

	while (off < (off_t) source->get_size()) {
		size_t size = min<size_t>(BUF_SIZE, source->get_size() - off);
		size_t ret = source->get_data(off, size, buf);
		assert(ret == size);
		ssize_t write_bytes = ROUNDUP(ret, 512);
		data_loc_t loc(io->get_file_id(), off);
		io_request req(buf, loc, write_bytes, WRITE);
		io->access(&req, 1);
		io->wait4complete(1);
		off += write_bytes;
	}
	printf("write all data\n");
	
	io->cleanup();
}

void print_help();

void comm_help(int argc, char *argv[])
{
	print_help();
}

void comm_list(int argc, char *argv[])
{
	init_io_system(configs, false);
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

	configs->add_options("writable=1");
	init_io_system(configs, false);
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
	{"load_part", comm_load_part_file2fs,
		"load_part file_name ext_file part_id: load part of the file to SAFS"},
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

	configs = config_map::create(conf_file);
	const struct command *comm = get_command(command);
	if (comm == NULL) {
		fprintf(stderr, "wrong command %s\n", command.c_str());
		print_help();
		return -1;
	}
	comm->func(argc - 3, argv + 3);
}
