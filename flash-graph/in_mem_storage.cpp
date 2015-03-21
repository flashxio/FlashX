/*
 * Copyright 2014 Open Connectome Project (http://openconnecto.me)
 * Written by Da Zheng (zhengda1936@gmail.com)
 *
 * This file is part of FlashGraph.
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

#include <stdlib.h>
#include <malloc.h>

#include <boost/format.hpp>

#include "log.h"
#include "safs_file.h"
#include "cache.h"
#include "slab_allocator.h"
#include "native_file.h"
#include "comp_io_scheduler.h"
#include "io_interface.h"
#include "in_mem_io.h"

#include "in_mem_storage.h"
#include "graph_file_header.h"

using namespace safs;

namespace fg
{

struct deleter
{
	void operator()(char *buf) const {
		free(buf);
	}
};

in_mem_graph::ptr in_mem_graph::load_graph(const std::string &file_name)
{
	native_file local_f(file_name);
	if (!local_f.exist())
		throw io_exception(file_name + std::string(" doesn't exist"));

	ssize_t size = local_f.get_size();
	assert(size > 0);

	in_mem_graph::ptr graph = in_mem_graph::ptr(new in_mem_graph());
	graph->graph_size = size;
	graph->graph_data = std::shared_ptr<char>((char *) malloc(size), deleter());
	assert(graph->graph_data);
	graph->graph_file_name = file_name;
	BOOST_LOG_TRIVIAL(info) << boost::format("load a graph of %1% bytes")
		% graph->graph_size;

	FILE *fd = fopen(file_name.c_str(), "r");
	if (fd == NULL)
		throw io_exception(std::string("can't open ") + file_name);
	if (fread(graph->graph_data.get(), size, 1, fd) != 1)
		throw io_exception(std::string("can't read from ") + file_name);
	fclose(fd);

	graph_header *header = (graph_header *) graph->graph_data.get();
	header->verify();

	return graph;
}

in_mem_graph::ptr in_mem_graph::load_safs_graph(const std::string &file_name)
{
	file_io_factory::shared_ptr io_factory = ::create_io_factory(file_name,
			REMOTE_ACCESS);

	in_mem_graph::ptr graph = in_mem_graph::ptr(new in_mem_graph());
	graph->graph_size = io_factory->get_file_size();
	size_t num_pages = ROUNDUP_PAGE(graph->graph_size) / PAGE_SIZE;
	char *graph_buf = NULL;
	int ret = posix_memalign((void **) &graph_buf, PAGE_SIZE,
			num_pages * PAGE_SIZE);
	graph->graph_data = std::shared_ptr<char>(graph_buf);
	BOOST_VERIFY(ret == 0);
	graph->graph_file_name = file_name;

	BOOST_LOG_TRIVIAL(info) << boost::format("load a graph of %1% bytes")
		% graph->graph_size;
#if 0
	graph->graph_file_id = io_factory->get_file_id();
#endif
	io_interface::ptr io = create_io(io_factory, thread::get_curr_thread());
	const size_t MAX_IO_SIZE = 256 * 1024 * 1024;
	for (off_t off = 0; (size_t) off < graph->graph_size; off += MAX_IO_SIZE) {
		data_loc_t loc(io_factory->get_file_id(), off);
		size_t req_size = min(MAX_IO_SIZE, graph->graph_size - off);
		io_request req(graph_buf + off, loc, req_size, READ);
		io->access(&req, 1);
		io->wait4complete(1);
	}

	graph_header *header = (graph_header *) graph_buf;
	header->verify();

	return graph;
}

file_io_factory::shared_ptr in_mem_graph::create_io_factory() const
{
	return file_io_factory::shared_ptr(new in_mem_io_factory(graph_data,
				graph_file_id, graph_file_name));
}

void in_mem_graph::dump(const std::string &file) const
{
	FILE *f = fopen(file.c_str(), "w");
	if (f == NULL) {
		perror("fopen");
		abort();
	}
	BOOST_VERIFY(fwrite(graph_data.get(), graph_size, 1, f));

	fclose(f);
}

}
