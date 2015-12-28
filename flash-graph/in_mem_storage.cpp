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
#include "graph_exception.h"

using namespace safs;

namespace fg
{

static const size_t GRAPH_CHUNK_SIZE_LOG = 30;

struct deleter
{
	void operator()(char *buf) const {
		free(buf);
	}
};

in_mem_graph::ptr in_mem_graph::create(const std::string &graph_name,
		std::shared_ptr<char> data, size_t size)
{
	// TODO It only supports one NUMA node in this setting right now.
	NUMA_mapper mapper(1, GRAPH_CHUNK_SIZE_LOG);
	safs::NUMA_buffer::ptr buf = safs::NUMA_buffer::create(data, size, mapper);
	return create(graph_name, buf, size);
}

in_mem_graph::ptr in_mem_graph::load_graph(const std::string &file_name)
{
	NUMA_mapper mapper(params.get_num_nodes(), GRAPH_CHUNK_SIZE_LOG);
	safs::NUMA_buffer::ptr numa_buf = safs::NUMA_buffer::load(file_name, mapper);
	assert(numa_buf);
	in_mem_graph::ptr graph = in_mem_graph::ptr(new in_mem_graph());
	graph->graph_size = numa_buf->get_length();
	graph->graph_data = numa_buf;
	graph->graph_file_name = file_name;
	BOOST_LOG_TRIVIAL(info) << boost::format("load a graph of %1% bytes")
		% graph->graph_size;

	safs::NUMA_buffer::cdata_info data = numa_buf->get_data(0, PAGE_SIZE);
	assert(data.first);
	graph_header *header = (graph_header *) data.first;
	if (!header->is_graph_file() || !header->is_right_version())
		throw wrong_format("wrong graph file or format version");
	return graph;
}

in_mem_graph::ptr in_mem_graph::load_safs_graph(const std::string &file_name)
{
	NUMA_mapper mapper(params.get_num_nodes(), GRAPH_CHUNK_SIZE_LOG);
	safs::NUMA_buffer::ptr numa_buf = safs::NUMA_buffer::load_safs(file_name,
			mapper);
	in_mem_graph::ptr graph = in_mem_graph::ptr(new in_mem_graph());
	graph->graph_size = numa_buf->get_length();
	graph->graph_data = numa_buf;
	graph->graph_file_name = file_name;

	BOOST_LOG_TRIVIAL(info) << boost::format("load a graph of %1% bytes")
		% graph->graph_size;
#if 0
	graph->graph_file_id = io_factory->get_file_id();
#endif

	safs::NUMA_buffer::cdata_info data = numa_buf->get_data(0, PAGE_SIZE);
	assert(data.first);
	graph_header *header = (graph_header *) data.first;
	if (!header->is_graph_file() || !header->is_right_version())
		throw wrong_format("wrong graph file or format version");
	return graph;
}

file_io_factory::shared_ptr in_mem_graph::create_io_factory() const
{
	return file_io_factory::shared_ptr(new in_mem_io_factory(graph_data,
				graph_file_id, graph_file_name));
}

void in_mem_graph::dump(const std::string &file) const
{
	graph_data->dump(file);
}

}
