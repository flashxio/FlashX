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

#include <boost/format.hpp>

#include <limits.h>

#include "log.h"
#include "native_file.h"
#include "safs_file.h"
#include "RAID_config.h"
#include "io_interface.h"
#include "file_mapper.h"

namespace safs
{

static std::vector<int> shuffle_disks(int num_disks)
{
	std::vector<int> permute(num_disks);
	for (size_t i = 0; i < permute.size(); i++)
		permute[i] = i;
	random_shuffle(permute.begin(), permute.end());
	return permute;
}

safs_file::safs_file(const RAID_config &conf, const std::string &file_name)
{
	sys_mapping_option = conf.get_mapping_option();
	sys_block_size = conf.get_block_size();
	native_dirs = conf.get_disks();
	for (unsigned i = 0; i < native_dirs.size(); i++)
		native_dirs[i] = part_file_info(
				native_dirs[i].get_file_name() + "/" + file_name,
				native_dirs[i].get_disk_id(), native_dirs[i].get_node_id());
	this->name = file_name;
}

std::vector<std::string> safs_file::erase_header_file(
		const std::vector<std::string> &files)
{
	std::vector<std::string> ret;
	for (auto it = files.begin(); it != files.end(); it++)
		if (*it != "header")
			ret.push_back(*it);
	return ret;
}

bool safs_file::exist() const
{
	std::set<int> part_ids;
	for (unsigned i = 0; i < native_dirs.size(); i++) {
		native_dir dir(native_dirs[i].get_file_name());
		if (!dir.exist())
			return false;
		std::vector<std::string> files;
		dir.read_all_files(files);
		if (files.size() > 1)
			files = erase_header_file(files);
		if (files.size() != 1) {
			fprintf(stderr, "%s doesn't have exactly one file\n",
					dir.get_name().c_str());
			return false;
		}
		part_ids.insert(atoi(files[0].c_str()));
	}
	if (part_ids.size() < native_dirs.size()) {
		fprintf(stderr, "there are duplicated partition ids in %s.\n",
				name.c_str());
		return false;
	}
	return true;
}

/*
 * This returns the absolute path of the data files.
 */
std::vector<std::string> safs_file::get_data_files() const
{
	std::vector<std::string> files;
	for (unsigned i = 0; i < native_dirs.size(); i++) {
		native_dir dir(native_dirs[i].get_file_name());
		std::vector<std::string> local_files;
		dir.read_all_files(local_files);
		if (local_files.size() > 1)
			local_files = erase_header_file(local_files);
		assert(local_files.size() == 1);
		files.push_back(dir.get_name() + "/" + local_files[0]);
	}
	return files;
}

ssize_t safs_file::get_size() const
{
	if (!exist())
		return -1;
	size_t ret = 0;
	std::vector<std::string> data_files = get_data_files();
	for (unsigned i = 0; i < data_files.size(); i++) {
		native_file f(data_files[i]);
		ret += f.get_size();
	}
	return ret;
}

bool safs_file::resize(size_t new_size)
{
	if (!exist()) {
		fprintf(stderr, "%s doesn't exist\n", name.c_str());
		return false;
	}

	// TODO right now we can only extend the file size.
	// otherwise, the system on top of it doesn't work correctly.
	// Why?
	ssize_t orig_size = get_size();
	assert(orig_size >= 0);
	if ((size_t) orig_size < new_size) {
		std::vector<std::string> data_files = get_data_files();
		std::vector<size_t> sizes_per_disk = get_size_per_disk(new_size);
		for (size_t i = 0; i < data_files.size(); i++) {
			native_file f(data_files[i]);
			bool ret = f.resize(sizes_per_disk[i]);
			if (!ret)
				return false;
		}
	}

	// Save the new file size to the header of the SAFS file.
	std::string header_file = get_header_file();
	BOOST_LOG_TRIVIAL(info) << "header file: " << header_file;
	if (!file_exist(header_file))
		return false;
	FILE *f = fopen(header_file.c_str(), "r+");
	if (f == NULL) {
		fprintf(stderr, "fopen %s: %s\n", header_file.c_str(), strerror(errno));
		return false;
	}

	safs_header header;
	size_t num_reads = fread(&header, sizeof(header), 1, f);
	if (num_reads != 1) {
		perror("fread");
		fclose(f);
		return false;
	}
	int seek_ret = fseek(f, 0, SEEK_SET);
	if (seek_ret < 0) {
		perror("seek");
		fclose(f);
		return false;
	}

	safs_header new_header(header.get_block_size(),
			header.get_mapping_option(), header.is_writable(),
			new_size);
	size_t num_writes = fwrite(&new_header, sizeof(new_header), 1, f);
	if (num_writes != 1) {
		perror("fwrite");
		fclose(f);
		return false;
	}
	fclose(f);
	return true;
}

bool safs_file::rename(const std::string &new_name)
{
	if (!exist()) {
		fprintf(stderr, "can't rename: the new name %s exists\n",
				new_name.c_str());
		return false;
	}

	for (unsigned i = 0; i < native_dirs.size(); i++) {
		native_file f(native_dirs[i].get_file_name());
		if (!f.rename(f.get_dir_name() + "/" + new_name))
			return false;
	}
	name = new_name;
	for (unsigned i = 0; i < native_dirs.size(); i++) {
		native_file f(native_dirs[i].get_file_name());
		native_dirs[i] = part_file_info(f.get_dir_name() + "/" + new_name,
				native_dirs[i].get_disk_id(), native_dirs[i].get_node_id());
	}
	return true;
}

std::vector<size_t> safs_file::get_size_per_disk(size_t file_size) const
{
	auto header = get_header();
	if (!header.is_valid())
		header = safs_header(sys_block_size, sys_mapping_option, false, file_size);
	file_mapper::ptr map = file_mapper::create(header, native_dirs, name);
	assert(map);
	std::vector<size_t> ret = map->get_size_per_disk(div_ceil<size_t>(file_size,
				PAGE_SIZE));
	for (size_t i = 0; i < ret.size(); i++)
		ret[i] = ret[i] * PAGE_SIZE;
	return ret;
}

bool safs_file::create_file(size_t file_size, int block_size,
		int mapping_option, safs_file_group::ptr group)
{
	// We use the random index to reorder the native directories.
	// So different files map their data chunks to disks in different order.
	// The benefit is that when we access data in the same location but from
	// different files, the data is likely fetched from different disks.
	// Thus, this leads to better I/O utilization.
	std::vector<int> dir_idxs;
	if (group == NULL)
		dir_idxs = shuffle_disks(native_dirs.size());
	else
		dir_idxs = group->add_file(*this);

	safs_header header(block_size, mapping_option, true, file_size);
	std::vector<size_t> sizes_per_disk = get_size_per_disk(file_size);
	for (unsigned i = 0; i < native_dirs.size(); i++) {
		native_dir dir(native_dirs[dir_idxs[i]].get_file_name());
		bool ret = dir.create_dir(true);
		if (!ret)
			return false;
		// We store the metadata of the SAFS in the directory that
		// stores the first part.
		if (i == 0) {
			header_file = dir.get_name() + "/header";
			FILE *f = fopen(header_file.c_str(), "w");
			if (f == NULL) {
				fprintf(stderr, "fopen %s: %s\n", header_file.c_str(),
						strerror(errno));
				return false;
			}
			size_t num_writes = fwrite(&header, sizeof(header), 1, f);
			if (num_writes != 1) {
				perror("fwrite");
				return false;
			}
			int ret = fclose(f);
			assert(ret == 0);
		}
		native_file f(dir.get_name() + "/" + itoa(i));
		ret = f.create_file(sizes_per_disk[i]);
		if (!ret)
			return false;
	}
	assert(!header_file.empty());
	return true;
}

bool safs_file::delete_file()
{
	for (unsigned i = 0; i < native_dirs.size(); i++) {
		native_dir dir(native_dirs[i].get_file_name());
		bool ret = dir.delete_dir(true);
		if (!ret)
			return false;
	}
	return true;
}

std::string safs_file::get_header_file() const
{
	if (!header_file.empty())
		return header_file;

	for (size_t i = 0; i < native_dirs.size(); i++) {
		std::string dir_str = native_dirs[i].get_file_name();
		if (file_exist(dir_str) && file_exist(dir_str + "/header")) {
			const_cast<safs_file *>(this)->header_file = dir_str + "/header";
			break;
		}
	}
	return header_file;
}

safs_header safs_file::get_header() const
{
	std::string header_file = get_header_file();
	if (!file_exist(header_file))
		return safs_header();
	FILE *f = fopen(header_file.c_str(), "r");
	if (f == NULL) {
		fprintf(stderr, "fopen %s: %s\n", header_file.c_str(), strerror(errno));
		return safs_header();
	}
	safs_header header;
	size_t num_reads = fread(&header, sizeof(header), 1, f);
	if (num_reads != 1) {
		perror("fread");
		return safs_header();
	}
	int ret = fclose(f);
	assert(ret == 0);
	return header;
}

bool safs_file::set_user_metadata(const std::vector<char> &data)
{
	std::string header_file = get_header_file();
	native_file native_f(header_file);
	assert(native_f.exist());

	FILE *f = fopen(header_file.c_str(), "r+");
	if (f == NULL) {
		fprintf(stderr, "fopen %s: %s\n", header_file.c_str(), strerror(errno));
		return false;
	}
	int ret = fseek(f, safs_header::get_header_size(), SEEK_SET);
	if (ret != 0) {
		perror("fseek");
		return false;
	}
	size_t num_writes = fwrite(data.data(), data.size(), 1, f);
	if (num_writes != 1) {
		perror("fwrite");
		return false;
	}
	ret = fclose(f);
	assert(ret == 0);
	return true;
}

std::vector<char> safs_file::get_user_metadata() const
{
	std::string header_file = get_header_file();
	native_file native_f(header_file);
	assert(native_f.exist());
	size_t file_size = native_f.get_size();
	assert(file_size >= safs_header::get_header_size());
	if (file_size == safs_header::get_header_size())
		return std::vector<char>();

	FILE *f = fopen(header_file.c_str(), "r");
	if (f == NULL) {
		fprintf(stderr, "fopen %s: %s\n", header_file.c_str(), strerror(errno));
		return std::vector<char>();
	}
	int ret = fseek(f, safs_header::get_header_size(), SEEK_SET);
	if (ret != 0) {
		perror("fseek");
		return std::vector<char>();
	}
	std::vector<char> data(file_size - safs_header::get_header_size());
	size_t num_reads = fread(data.data(), data.size(), 1, f);
	if (num_reads != 1) {
		perror("fread");
		return std::vector<char>();
	}
	ret = fclose(f);
	assert(ret == 0);
	return data;
}

namespace
{

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

	file_data_source(int fd, size_t file_size) {
		this->fd = fd;
		this->file_size = file_size;
	}
public:
	// Create a data source that contains `size' bytes.
	static std::shared_ptr<file_data_source> create(const std::string &ext_file,
			size_t size) {
		int fd = open(ext_file.c_str(), O_RDONLY);
		if (fd < 0) {
			fprintf(stderr, "can't open %s: %s\n", ext_file.c_str(), strerror(errno));
			return std::shared_ptr<file_data_source>();
		}
		native_file f(ext_file);
		size_t file_size = f.get_size();
		if (size > file_size)
			size = file_size;
		return std::shared_ptr<file_data_source>(new file_data_source(fd, size));
	}

	virtual ssize_t get_data(off_t off, size_t size, char *buf) const {
		long new_off = lseek(fd, off, SEEK_SET);
		BOOST_VERIFY(new_off == off);
		ssize_t ret = complete_read(fd, buf, size);
		if (ret < 0)
			perror("complete_read");
		return ret;
	}

	virtual size_t get_size() const {
		return file_size;
	}
};

const int BUF_SIZE = 1024 * 64 * 4096;

}

bool safs_file::load_data(const std::string &ext_file, size_t block_size)
{
	std::shared_ptr<data_source> source = file_data_source::create(ext_file,
			get_size());
	if (source == NULL)
		return false;

	// If the file in SAFS doesn't exist, create a new one.
	if (!exist())
		create_file(source->get_size(), block_size);

	file_io_factory::shared_ptr factory = create_io_factory(name,
			REMOTE_ACCESS);
	if (factory == NULL) {
		fprintf(stderr, "can't create I/O factory\n");
		return false;
	}
	assert((size_t) factory->get_file_size() >= source->get_size());

	thread *curr_thread = thread::get_curr_thread();
	assert(curr_thread);
	io_interface::ptr io = create_io(factory, curr_thread);

	char *buf = (char *) valloc(BUF_SIZE);
	off_t off = 0;

	while (off < (off_t) source->get_size()) {
		size_t size = min<size_t>(BUF_SIZE, source->get_size() - off);
		size_t ret = source->get_data(off, size, buf);
		assert(ret == size);
		ssize_t write_bytes = ROUNDUP(ret, 512);
		memset(buf + size, 0, write_bytes - size);
		data_loc_t loc(io->get_file_id(), off);
		io_request req(buf, loc, write_bytes, WRITE);
		io->access(&req, 1);
		io->wait4complete(1);
		off += write_bytes;
	}
	io->cleanup();
	free(buf);
	return true;
}

size_t get_all_safs_files(std::set<std::string> &files)
{
	std::set<std::string> all_files;
	const RAID_config &conf = get_sys_RAID_conf();

	// First find all individual file names in the root directories.
	for (int i = 0; i < conf.get_num_disks(); i++) {
		std::string dir_name = conf.get_disk(i).get_file_name();
		native_dir dir(dir_name);
		std::vector<std::string> file_names;
		dir.read_all_files(file_names);
		all_files.insert(file_names.begin(), file_names.end());
	}

	for (std::set<std::string>::const_iterator it = all_files.begin();
			it != all_files.end(); it++) {
		safs_file file(conf, *it);
		if (file.exist()) {
			files.insert(*it);
		}
		else {
			BOOST_LOG_TRIVIAL(error) << boost::format("%1% is corrupted")
				% file.get_name();
		}
	}
	return 0;
}

class naive_file_group: public safs_file_group
{
	size_t num_files;
	size_t num_disks;
public:
	naive_file_group(const RAID_config &conf) {
		num_files = 0;
		num_disks = conf.get_num_disks();
	}
	std::vector<int> add_file(safs_file &file) {
		std::vector<int> ret(num_disks);
		for (size_t i = 0; i < num_disks; i++)
			ret[i] = i;
		num_files++;
		return ret;
	}
	std::string get_name() const {
		return "naive";
	}
};

class rotate_file_group: public safs_file_group
{
	size_t num_files;
	size_t num_disks;
public:
	rotate_file_group(const RAID_config &conf) {
		num_files = 0;
		num_disks = conf.get_num_disks();
	}
	std::vector<int> add_file(safs_file &file) {
		std::vector<int> ret(num_disks);
		for (size_t i = 0; i < num_disks; i++)
			ret[i] = (num_files + i) % num_disks;
		num_files++;
		return ret;
	}
	std::string get_name() const {
		return "rotate";
	}
};

class rand_rotate_file_group: public safs_file_group
{
	// The base permute is the permutation that other permutations are based
	// on. Other permutations just rotate the base permutation from a random
	// location. Every #disks files share the same base permutation.
	std::vector<std::vector<int> > base_permutes;
	std::vector<std::vector<int> > rand_rotates;
	size_t num_files;
public:
	rand_rotate_file_group(const RAID_config &conf);
	std::vector<int> add_file(safs_file &file);
	std::string get_name() const {
		return "rand_rotate";
	}
};

rand_rotate_file_group::rand_rotate_file_group(const RAID_config &conf)
{
	int num_disks = conf.get_num_disks();
	num_files = 0;
	base_permutes.push_back(shuffle_disks(num_disks));
	// Every #disks files share the same base permutation.
	rand_rotates.push_back(shuffle_disks(num_disks));
}

std::vector<int> rand_rotate_file_group::add_file(safs_file &file)
{
	size_t num_disks = base_permutes.front().size();
	size_t base_idx = num_files / num_disks;
	if (base_idx >= base_permutes.size()) {
		base_permutes.push_back(shuffle_disks(num_disks));
		rand_rotates.push_back(shuffle_disks(num_disks));
	}
	assert(base_permutes.size() > base_idx);

	std::vector<int> base = base_permutes[base_idx];
	std::vector<int> ret(num_disks);
	size_t rotate = rand_rotates[base_idx][num_files % num_disks];
	for (size_t i = 0; i < ret.size(); i++)
		ret[i] = base[(rotate + i) % num_disks];
	num_files++;
	return ret;
}

safs_file_group::ptr safs_file_group::create(const RAID_config &conf,
		group_t type)
{
	switch (type) {
		case group_t::NAIVE:
			return safs_file_group::ptr(new naive_file_group(conf));
		case group_t::ROTATE:
			return safs_file_group::ptr(new rotate_file_group(conf));
		case group_t::RAND_ROTATE:
			return safs_file_group::ptr(new rand_rotate_file_group(conf));
		default:
			fprintf(stderr, "unknow group type: %d\n", type);
			return safs_file_group::ptr();
	}
}

bool exist_safs_file(const std::string &name)
{
	if (!is_safs_init())
		return false;

	safs::safs_file mat_f(safs::get_sys_RAID_conf(), name);
	return mat_f.exist();
}

ssize_t get_safs_size(const std::string &name)
{
	if (!is_safs_init())
		return -1;

	safs::safs_file mat_f(safs::get_sys_RAID_conf(), name);
	return mat_f.get_size();
}

}
