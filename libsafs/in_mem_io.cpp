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

#include "in_mem_io.h"
#include "slab_allocator.h"
#include "native_file.h"
#include "safs_file.h"
#include "io_interface.h"

namespace safs
{

namespace
{
class numa_delete
{
	size_t size;
public:
	numa_delete(size_t size) {
		this->size = size;
	}

	void operator()(char *buf) const {
		numa_free(buf, size);
	}
};

}

NUMA_buffer::NUMA_buffer(std::shared_ptr<char> data, size_t length,
		const NUMA_mapper &_mapper): mapper(_mapper)
{
	assert(mapper.get_num_nodes() == 1);
	this->length = length;
	bufs.resize(1);
	buf_lens.resize(1);
	bufs[0] = data;
	buf_lens[0] = length;
}

NUMA_buffer::NUMA_buffer(size_t length,
		const NUMA_mapper &_mapper): mapper(_mapper)
{
	length = ROUNDUP(length, PAGE_SIZE);
	this->length = length;
	bufs.resize(mapper.get_num_nodes());
	buf_lens.resize(bufs.size());

	// Calculate the data size in each NUMA node.
	size_t last_off = length - 1;
	for (size_t i = 0; i < buf_lens.size(); i++) {
		auto loc = mapper.map2physical(last_off);
		// Data ranges are assigned to NUMA nodes in a round-robin fashion.
		// So when we determine the data size in a NUMA node, it shouldn't
		// have been decided before.
		assert(buf_lens[loc.first] == 0);
		// We look at the last byte in a range, the length should increase by 1.
		buf_lens[loc.first] = loc.second + 1;
		// If the data in the buffer is small, not all NUMA nodes get data.
		if (last_off < mapper.get_range_size())
			break;
		last_off = ROUND(last_off, mapper.get_range_size()) - 1;
	}

	// Allocate memory for each NUMA node.
	for (size_t i = 0; i < bufs.size(); i++) {
		if (buf_lens[i] > 0) {
			bufs[i] = std::shared_ptr<char>(
					(char *) numa_alloc_onnode(buf_lens[i], i),
					numa_delete(buf_lens[i]));
			assert(bufs[i]);
		}
	}
}

NUMA_buffer::data_loc_info NUMA_buffer::get_data_loc(off_t off,
		size_t size) const
{
	if (off + size > length) {
		fprintf(stderr, "data of %ld bytes in %ld exceeds length %ld\n",
				size, off, length);
		return data_loc_info(-1, -1, 0);
	}

	std::pair<int, size_t> loc = mapper.map2physical(off);
	if (buf_lens[loc.first] <= loc.second) {
		fprintf(stderr, "local off %ld on node %d exceeds local size %ld\n",
				loc.second, loc.first, buf_lens[loc.first]);
		return data_loc_info(-1, -1, 0);
	}

	size_t local_size;
	// If it's at the beginning of the range, we have the entire range.
	if (loc.second % mapper.get_range_size() == 0)
		local_size = mapper.get_range_size();
	else
		local_size = ROUNDUP(loc.second, mapper.get_range_size()) - loc.second;
	// This is the amount of data we should get from the range.
	local_size = std::min(local_size, size);
	// If the amount of required data is larger than the data existing
	// in the NUMA node.
	if (local_size > buf_lens[loc.first] - loc.second) {
		fprintf(stderr, "local end off %ld on node %d exceeds local size %ld\n",
				loc.second + local_size, loc.first, buf_lens[loc.first]);
		return data_loc_info(-1, -1, 0);
	}

	return data_loc_info(loc.first, loc.second, local_size);
}

NUMA_buffer::cdata_info NUMA_buffer::get_data(off_t off, size_t size) const
{
	data_loc_info loc = get_data_loc(off, size);
	if (loc.node_id < 0)
		return cdata_info(NULL, 0);

	return cdata_info(bufs[loc.node_id].get() + loc.local_off, loc.local_size);
}

NUMA_buffer::data_info NUMA_buffer::get_data(off_t off, size_t size)
{
	data_loc_info loc = get_data_loc(off, size);
	if (loc.node_id < 0)
		return data_info(NULL, 0);

	return data_info(bufs[loc.node_id].get() + loc.local_off, loc.local_size);
}

void NUMA_buffer::copy_from(const char *buf, size_t size, off_t off)
{
	// The required data may not be stored in contiguous memory.
	while (size > 0) {
		auto info = get_data(off, size);
		assert(info.first);
		memcpy(info.first, buf, info.second);

		size -= info.second;
		off += info.second;
		buf += info.second;
		// If the required data isn't stored in contiguous memory,
		// the next portion should be stored at the beginning of a data range.
		if (size > 0)
			assert(off % mapper.get_range_size() == 0);
	}
}

void NUMA_buffer::copy_to(char *buf, size_t size, off_t off) const
{
	// The required data may not be stored in contiguous memory.
	while (size > 0) {
		auto info = get_data(off, size);
		assert(info.first);
		memcpy(buf, info.first, info.second);

		size -= info.second;
		off += info.second;
		buf += info.second;
		// If the required data isn't stored in contiguous memory,
		// the next portion should be stored at the beginning of a data range.
		if (size > 0)
			assert(off % mapper.get_range_size() == 0);
	}
}

NUMA_buffer::ptr NUMA_buffer::load(const std::string &file_name,
		const NUMA_mapper &mapper)
{
	native_file local_f(file_name);
	if (!local_f.exist())
		throw io_exception(file_name + std::string(" doesn't exist"));

	ssize_t file_size = local_f.get_size();
	assert(file_size > 0);

	NUMA_buffer::ptr numa_buf(new NUMA_buffer(file_size, mapper));
	FILE *fd = fopen(file_name.c_str(), "r");
	if (fd == NULL) {
		int err = errno;
		fclose(fd);
		throw io_exception(boost::str(boost::format("can't open %1%: %2%")
					% file_name % strerror(err)));
	}
	for (off_t off = 0; off < file_size;) {
		size_t load_size = std::min(mapper.get_range_size(),
				(size_t) (file_size - off));
		data_info data = numa_buf->get_data(off, load_size);
		// The total size of all physical buffers may be larger than
		// the NUMA buffer size.
		size_t size = std::min(data.second, (size_t) (file_size - off));
		if (fread(data.first, size, 1, fd) != 1) {
			int err = errno;
			fclose(fd);
			throw io_exception(boost::str(boost::format(
							"can't read from %1%: %2%")
						% file_name % strerror(err)));
		}
		off += load_size;
		if (off < file_size)
			assert(off % mapper.get_range_size() == 0);
	}
	fclose(fd);

	return numa_buf;
}

void NUMA_buffer::dump(const std::string &file_name)
{
	FILE *f = fopen(file_name.c_str(), "w");
	if (f == NULL)
		throw io_exception(boost::str(boost::format("can't open %1%: %2%")
					% file_name % strerror(errno)));

	for (off_t off = 0; (size_t) off < get_length();) {
		size_t write_size = std::min(mapper.get_range_size(),
				get_length() - off);
		data_info data = get_data(off, write_size);
		// The total size of all physical buffers may be larger than
		// the NUMA buffer size.
		size_t size = std::min(data.second, get_length() - off);
		if (fwrite(data.first, size, 1, f) != 1) {
			int err = errno;
			fclose(f);
			throw io_exception(boost::str(boost::format(
							"can't write to %1%: %2%")
						% file_name % strerror(err)));
		}
		off += write_size;
		if ((size_t) off < get_length())
			assert(off % mapper.get_range_size() == 0);
	}
	fclose(f);
}

NUMA_buffer::ptr NUMA_buffer::load_safs(const std::string &file_name,
		const NUMA_mapper &mapper)
{
	file_io_factory::shared_ptr io_factory = create_io_factory(file_name,
			REMOTE_ACCESS);
	if (io_factory == NULL)
		throw io_exception(std::string("can't create io factory for ")
				+ file_name);

	size_t file_size = io_factory->get_file_size();
	NUMA_buffer::ptr numa_buf(new NUMA_buffer(file_size, mapper));
	io_interface::ptr io = create_io(io_factory, thread::get_curr_thread());
	if (io == NULL)
		throw io_exception(std::string("can't create io instance for ")
				+ file_name);
	const size_t MAX_IO_SIZE = std::min(256UL * 1024 * 1024,
			mapper.get_range_size());
	for (off_t off = 0; (size_t) off < file_size; ) {
		data_loc_t loc(io_factory->get_file_id(), off);
		size_t load_size = std::min(MAX_IO_SIZE, file_size - off);
		data_info data = numa_buf->get_data(off, load_size);
		size_t req_size = std::min(data.second, file_size - off);
		io_request req(data.first, loc, req_size, READ);
		io->access(&req, 1);
		io->wait4complete(1);
		off += load_size;
		if ((size_t) off < file_size)
			assert(off % MAX_IO_SIZE == 0);
	}
	return numa_buf;
}

NUMA_buffer::ptr NUMA_buffer::create(std::shared_ptr<char> data, size_t length,
		const NUMA_mapper &mapper)
{
	if (mapper.get_num_nodes() != 1) {
		throw io_exception(
				"can't create a NUMA buffer from a raw byte array on multiple nodes");
	}
	return NUMA_buffer::ptr(new NUMA_buffer(data, length, mapper));
}

class in_mem_byte_array: public page_byte_array
{
	off_t off;
	size_t size;
	const char *pages;

	void assign(in_mem_byte_array &arr) {
		this->off = arr.off;
		this->size = arr.size;
		this->pages = arr.pages;
	}

	in_mem_byte_array(in_mem_byte_array &arr) {
		assign(arr);
	}

	in_mem_byte_array &operator=(in_mem_byte_array &arr) {
		assign(arr);
		return *this;
	}
public:
	in_mem_byte_array(byte_array_allocator &alloc): page_byte_array(alloc) {
		off = 0;
		size = 0;
		pages = NULL;
	}

	in_mem_byte_array(const io_request &req, const char *pages,
			byte_array_allocator &alloc): page_byte_array(alloc) {
		this->off = req.get_offset();
		this->size = req.get_size();
		this->pages = pages;
	}

	virtual off_t get_offset() const {
		return off;
	}

	virtual off_t get_offset_in_first_page() const {
		return off % PAGE_SIZE;
	}

	virtual const char *get_page(int pg_idx) const {
		return pages + pg_idx * PAGE_SIZE;
	}

	virtual size_t get_size() const {
		return size;
	}

	void lock() {
		ABORT_MSG("lock isn't implemented");
	}

	void unlock() {
		ABORT_MSG("unlock isn't implemented");
	}

	page_byte_array *clone() {
		in_mem_byte_array *arr = (in_mem_byte_array *) get_allocator().alloc();
		*arr = *this;
		return arr;
	}
};

class in_mem_byte_array_allocator: public byte_array_allocator
{
	class array_initiator: public obj_initiator<in_mem_byte_array>
	{
		in_mem_byte_array_allocator *alloc;
	public:
		array_initiator(in_mem_byte_array_allocator *alloc) {
			this->alloc = alloc;
		}

		virtual void init(in_mem_byte_array *obj) {
			new (obj) in_mem_byte_array(*alloc);
		}
	};

	class array_destructor: public obj_destructor<in_mem_byte_array>
	{
	public:
		void destroy(in_mem_byte_array *obj) {
			obj->~in_mem_byte_array();
		}
	};

	obj_allocator<in_mem_byte_array> allocator;
public:
	in_mem_byte_array_allocator(thread *t): allocator(
			"byte-array-allocator", t->get_node_id(), false, 1024 * 1024,
			params.get_max_obj_alloc_size(),
			obj_initiator<in_mem_byte_array>::ptr(new array_initiator(this)),
			obj_destructor<in_mem_byte_array>::ptr(new array_destructor())) {
	}

	virtual page_byte_array *alloc() {
		return allocator.alloc_obj();
	}

	virtual void free(page_byte_array *arr) {
		allocator.free((in_mem_byte_array *) arr);
	}
};

void in_mem_io::process_req(const io_request &req)
{
	assert(req.get_req_type() == io_request::USER_COMPUTE);
	// The byte array assumes the data is stored in pages.
	off_t off = ROUND_PAGE(req.get_offset());
	size_t size = ROUNDUP_PAGE(req.get_offset() + req.get_size()) - off;
	NUMA_buffer::data_info info = data->get_data(off, size);
	assert(info.first);
	char *first_page = info.first;
	// If the data in the buffer isn't stored in contiguous memory,
	// we need to copy them to a piece of contiguous memory.
	// This should happen very rarely if the range size in the NUMA mapper
	// is very large. TODO but we can avoid the memory copy entirely.
	if (info.second < size) {
		first_page = (char *) malloc(size);
		data->copy_to(first_page, size, off);
	}

	in_mem_byte_array byte_arr(req, first_page, *array_allocator);
	user_compute *compute = req.get_compute();
	compute->run(byte_arr);
	comp_io_sched->post_comp_process(compute);

	// If we allocate memory to store the data, we need to delete it now.
	if (info.second < size)
		free(first_page);
}

void in_mem_io::process_computes()
{
	while (true) {
		comp_io_sched->get_requests(req_buf, req_buf.get_size());
		if (req_buf.is_empty())
			break;

		while (!req_buf.is_empty()) {
			io_request new_req = req_buf.pop_front();
			process_req(new_req);
		}
	}
}

io_status in_mem_io::access(char *buf, off_t off, ssize_t size, int access_method)
{
	if (access_method == READ)
		data->copy_to(buf, size, off);
	else
		data->copy_from(buf, size, off);
	return IO_OK;
}

void in_mem_io::access(io_request *requests, int num, io_status *)
{
	for (int i = 0; i < num; i++) {
		io_request &req = requests[i];
		if (req.get_req_type() == io_request::USER_COMPUTE) {
			// Let's possess a reference to the user compute first. process_req()
			// will release the reference when the user compute is completed.
			req.get_compute()->inc_ref();
			process_req(req);
		}
		else {
			assert(req.get_req_type() == io_request::BASIC_REQ);
			if (req.get_access_method() == READ)
				data->copy_to(req.get_buf(), req.get_size(), req.get_offset());
			else
				data->copy_from(req.get_buf(), req.get_size(), req.get_offset());
			io_request *reqs[1];
			reqs[0] = &req;
			if (this->have_callback())
				this->get_callback().invoke(reqs, 1);
		}
	}
	process_computes();
	comp_io_sched->gc_computes();
}

in_mem_io::in_mem_io(NUMA_buffer::ptr data, int file_id,
		thread *t): io_interface(t, safs_header()), req_buf(get_node_id(), 1024)
{
	this->data = data;
	this->file_id = file_id;
	array_allocator = std::unique_ptr<byte_array_allocator>(
			new in_mem_byte_array_allocator(t));
	comp_io_sched = comp_io_scheduler::ptr(
			new default_comp_io_scheduler(get_node_id()));
	comp_io_sched->set_io(this);
}

class in_mem_io_select: public io_select
{
public:
	virtual bool add_io(io_interface::ptr io) {
		return true;
	}
	virtual int num_pending_ios() const {
		return 0;
	}
	virtual int wait4complete(int num_to_complete) {
		return 0;
	}
};

io_select::ptr in_mem_io::create_io_select() const
{
	return io_select::ptr(new in_mem_io_select());
}

}
