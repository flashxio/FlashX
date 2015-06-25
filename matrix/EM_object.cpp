/*
 * Copyright 2015 Open Connectome Project (http://openconnecto.me)
 * Written by Da Zheng (zhengda1936@gmail.com)
 *
 * This file is part of FlashMatrix.
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

#include "EM_object.h"
#include "mem_worker_thread.h"

namespace fm
{

namespace detail
{

EM_object::file_holder::ptr EM_object::file_holder::create_temp(
		const std::string &name, size_t num_bytes)
{
	char *tmp = tempnam(".", name.c_str());
	std::string tmp_name = basename(tmp);
	safs::safs_file f(safs::get_sys_RAID_conf(), tmp_name);
	assert(!f.exist());
	bool ret = f.create_file(num_bytes);
	assert(ret);
	file_holder::ptr holder(new file_holder(tmp_name, false));
	free(tmp);
	return holder;
}

EM_object::file_holder::~file_holder()
{
	if (!persistent) {
		safs::safs_file f(safs::get_sys_RAID_conf(), file_name);
		assert(f.exist());
		f.delete_file();
	}
}

bool EM_object::file_holder::set_persistent(const std::string &new_name)
{
	safs::safs_file f(safs::get_sys_RAID_conf(), file_name);
	if (!f.rename(new_name))
		return false;
	persistent = true;
	this->file_name = new_name;
	return true;
}

safs::io_interface::ptr EM_object::io_set::create_io()
{
	thread *t = thread::get_curr_thread();
	assert(t);
	pthread_spin_lock(&io_lock);
	auto it = thread_ios.find(t);
	if (it == thread_ios.end()) {
		safs::io_interface::ptr io = safs::create_io(factory, t);
		io->set_callback(portion_callback::ptr(new portion_callback()));
		thread_ios.insert(std::pair<thread *, safs::io_interface::ptr>(t, io));
		pthread_setspecific(io_key, io.get());
		pthread_spin_unlock(&io_lock);
		return io;
	}
	else {
		safs::io_interface::ptr io = it->second;
		pthread_spin_unlock(&io_lock);
		return io;
	}
}

EM_object::io_set::io_set(safs::file_io_factory::shared_ptr factory)
{
	this->factory = factory;
	int ret = pthread_key_create(&io_key, NULL);
	assert(ret == 0);
	pthread_spin_init(&io_lock, PTHREAD_PROCESS_PRIVATE);
}

EM_object::io_set::~io_set()
{
	pthread_spin_destroy(&io_lock);
	pthread_key_delete(io_key);
	thread_ios.clear();
}

bool EM_object::io_set::has_io() const
{
	return pthread_getspecific(io_key) != NULL;
}

safs::io_interface &EM_object::io_set::get_curr_io() const
{
	void *io_addr = pthread_getspecific(io_key);
	if (io_addr)
		return *(safs::io_interface *) io_addr;
	else {
		safs::io_interface::ptr io = const_cast<io_set *>(this)->create_io();
		return *io;
	}
}

}

}
