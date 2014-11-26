#ifndef __IN_MEM_STORAGE_H__
#define __IN_MEM_STORAGE_H__

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

#include "io_interface.h"

class in_mem_io;
class thread_safe_page;

namespace fg
{

class in_mem_graph
{
	size_t graph_size;
	char *graph_data;
	int graph_file_id;
	std::string graph_file_name;

	in_mem_graph() {
		graph_data = NULL;
		graph_size = 0;
		graph_file_id = -1;
	}
public:
	typedef std::shared_ptr<in_mem_graph> ptr;

	/*
	 * in_mem_graph takes the memory buffer and is responsible to free it
	 * afterwards.
	 */
	static ptr create(const std::string &graph_name, char *buf, size_t size) {
		in_mem_graph *g = new in_mem_graph();
		g->graph_size = size;
		g->graph_data = buf;
		// TODO we don't init graph_file_id here.
		g->graph_file_name = graph_name;
		return ptr(g);
	}

	static ptr load_graph(const std::string &graph_file);
	static ptr load_safs_graph(const std::string &graph_file);

	~in_mem_graph() {
		free(graph_data);
	}

	safs::file_io_factory::shared_ptr create_io_factory() const;

	friend class in_mem_io;
};

}

#endif
